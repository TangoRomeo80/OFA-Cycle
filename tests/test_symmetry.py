"""Independent checks of Duval canonicalization and Burnside counting."""

import itertools
import os
import random
import unittest

from cycle_ofa.symmetry import bracelet_count, canonical, least_rotation


def brute_canonical(word):
    word = tuple(word)
    if not word:
        return ()
    reverse = word[::-1]
    return min(
        image[i:] + image[:i]
        for image in (word, reverse)
        for i in range(len(word))
    )


class SymmetryTests(unittest.TestCase):
    def test_exhaustive_ternary_through_length_six(self):
        for n in range(7):
            for word in itertools.product(range(3), repeat=n):
                self.assertEqual(canonical(word), brute_canonical(word), word)

    @unittest.skipUnless(os.environ.get("CYCLE_OFA_EXHAUSTIVE") == "1", "use validate --exhaustive")
    def test_exhaustive_ternary_through_length_eight(self):
        for n in (7, 8):
            for word in itertools.product(range(3), repeat=n):
                self.assertEqual(canonical(word), brute_canonical(word), word)

    def test_randomized_larger_words(self):
        rng = random.Random(20260822)
        for _ in range(300):
            word = tuple(rng.randrange(10) for _ in range(rng.randrange(9, 201)))
            self.assertEqual(canonical(word), brute_canonical(word))

    def test_least_rotation_orientation(self):
        word = (3, 1, 2, 1, 4)
        self.assertEqual(least_rotation(word), (1, 2, 1, 4, 3))
        self.assertEqual(least_rotation(()), ())
        self.assertEqual(least_rotation((2,) * 1000), (2,) * 1000)

    def test_bracelets_by_enumeration(self):
        for q in (1, 2, 3):
            for n in range(1, 8):
                count = len({brute_canonical(w) for w in itertools.product(range(q), repeat=n)})
                self.assertEqual(bracelet_count(n, q), count, (n, q))

    def test_paper_counts(self):
        for n, m, expected in (
            (8, 2, 498), (10, 2, 3210), (12, 2, 22913),
            (8, 3, 4435), (9, 3, 15084), (8, 4, 25395),
        ):
            self.assertEqual(bracelet_count(n, m + 1), expected)

    def test_invalid_counts(self):
        for n, q in ((0, 2), (3, 0), (-3, 2), (3, 1.5), (True, 2)):
            with self.assertRaises(ValueError):
                bracelet_count(n, q)


if __name__ == "__main__":
    unittest.main()
