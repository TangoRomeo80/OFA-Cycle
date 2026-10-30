"""Check the variance identity with exact small-instance moment values."""

from fractions import Fraction
from itertools import product
from random import Random
import statistics
import unittest

from cycle_ofa.diagnostics import conditional_diagnostic_trajectory, diagnose_variance, positional_variance
from cycle_ofa.simulation import conditional_trajectory


class Ranks:
    def __init__(self, values):
        self.values = iter(values)

    def randrange(self, stop):
        value = next(self.values)
        assert 0 <= value < stop
        return value


class DiagnosticTests(unittest.TestCase):
    def test_positional_second_moment(self):
        self.assertAlmostEqual(positional_variance(1, 1), 1 / 48)
        self.assertAlmostEqual(positional_variance(1, 2), 11 / 144)
        for n in (3, 8, 200):
            self.assertAlmostEqual(positional_variance(n, n), n * n / 48)
        for a, b in ((1, 7), (4, 9), (2, 2)):
            self.assertEqual(positional_variance(a, b), positional_variance(b, a))
            self.assertGreater(positional_variance(a, b), 0)
        with self.assertRaises(ValueError):
            positional_variance(0, 1)

    def test_seeded_cost_and_rng_consumption_match_existing_simulator(self):
        for n, m in ((3, 1), (3, 2), (8, 4), (20, 5), (200, 10)):
            for seed in range(7):
                rng_a, rng_b = Random(seed), Random(seed)
                y, removed = conditional_diagnostic_trajectory(n, m, rng_a)
                self.assertEqual(y, conditional_trajectory(n, m, rng_b))
                self.assertEqual(rng_a.getstate(), rng_b.getstate())
                self.assertGreater(removed, 0)

    def test_exact_total_variance_by_enumerating_rank_sequences(self):
        # Var(Y) and Var(C) independently obtained by exact moment recurrence.
        # All integer weighted-rank sequences are equally probable; unused
        # suffix ranks merely repeat a path that reached its singleton early.
        for n, m, mean_y, variance_y, variance_c in (
            (3, 1, Fraction(17, 12), Fraction(0), Fraction(41, 144)),
            (4, 1, Fraction(71, 32), Fraction(49, 3072), Fraction(1837, 3072)),
            (3, 2, Fraction(67, 27), Fraction(161, 2916), Fraction(2959, 5832)),
        ):
            count = (2 * n) ** (n * m - 1)
            ey = ey2 = es = 0.0
            for ranks in product(range(2 * n), repeat=n * m - 1):
                y, removed = conditional_diagnostic_trajectory(n, m, Ranks(ranks))
                ey += y / count
                ey2 += y * y / count
                es += removed / count
            self.assertAlmostEqual(ey, float(mean_y), places=10)
            self.assertAlmostEqual(ey2 - ey * ey, float(variance_y), places=10)
            self.assertAlmostEqual(es, float(variance_c - variance_y), places=10)
            self.assertAlmostEqual(ey2 - ey * ey + es, float(variance_c), places=10)

    def test_statistics_and_degenerate_ratio(self):
        summary = diagnose_variance(3, 1, runs=20, seed=9)
        self.assertEqual(summary.conditional_variance, 0)
        self.assertIsNone(summary.variance_ratio)
        self.assertAlmostEqual(summary.mean_removed_variance, 41 / 144)
        rng = Random(31)
        values = [conditional_diagnostic_trajectory(5, 2, rng) for _ in range(100)]
        summary = diagnose_variance(5, 2, runs=100, seed=31)
        self.assertAlmostEqual(summary.conditional_mean, statistics.mean(v[0] for v in values))
        self.assertAlmostEqual(summary.conditional_variance, statistics.variance(v[0] for v in values))
        self.assertAlmostEqual(summary.mean_removed_variance, statistics.mean(v[1] for v in values))
        self.assertAlmostEqual(summary.removed_variance_sample_variance, statistics.variance(v[1] for v in values))
        self.assertAlmostEqual(summary.implied_direct_variance, summary.conditional_variance + summary.mean_removed_variance)
        self.assertGreater(summary.variance_ratio, 1)
        self.assertEqual(summary.as_dict()["runs"], 100)

    def test_validation(self):
        for args in ((2, 1, 10, 1), (3, 0, 10, 1), (3, 1, 1, 1), (3, 1, 10, True)):
            with self.assertRaises(ValueError):
                diagnose_variance(*args)


if __name__ == "__main__":
    unittest.main()
