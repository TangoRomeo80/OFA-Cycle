"""Independent checks for sampling, wrap-around geometry, and simulation."""

from fractions import Fraction
from itertools import product
from math import sqrt
from random import Random
import random
import unittest

from cycle_ofa.simulation import (
    FenwickTree,
    conditional_trajectory,
    direct_trajectory,
    simulate,
)


class FenwickTests(unittest.TestCase):
    def test_ranks_with_zeros_and_boundaries(self):
        tree = FenwickTree([0, 2, 0, 3, 0])
        self.assertEqual([tree.find_by_order(i) for i in range(5)], [1, 1, 3, 3, 3])
        self.assertEqual([tree.prefix_sum(i) for i in range(6)], [0, 0, 2, 2, 5, 5])
        tree.set_weight(1, 0)
        tree.add(4, 2)
        self.assertEqual(tree.weights, (0, 0, 0, 3, 2))
        self.assertEqual(tree.total, 5)
        self.assertEqual([tree.find_by_order(i) for i in range(5)], [3, 3, 3, 4, 4])

    def test_random_updates_against_linear_selection(self):
        rng = Random(718)
        weights = [rng.randrange(8) for _ in range(19)]
        tree = FenwickTree(weights)
        for _ in range(150):
            index, weight = rng.randrange(19), rng.randrange(8)
            tree.set_weight(index, weight)
            weights[index] = weight
            expanded = [i for i, value in enumerate(weights) for _ in range(value)]
            self.assertEqual(tree.total, len(expanded))
            self.assertEqual([tree.find_by_order(i) for i in range(tree.total)], expanded)
            self.assertEqual([tree.prefix_sum(i) for i in range(20)], [sum(weights[:i]) for i in range(20)])

    def test_invalid_inputs(self):
        for weights in ([], [1, -1], [True], [1.0]):
            with self.assertRaises(ValueError):
                FenwickTree(weights)
        tree = FenwickTree([0, 0])
        for order in (-1, 0, 1, 0.5, True):
            with self.assertRaises(ValueError):
                tree.find_by_order(order)
        with self.assertRaises(ValueError):
            tree.add(0, -1)
        with self.assertRaises(ValueError):
            tree.set_weight(2, 1)
        with self.assertRaises(ValueError):
            tree.prefix_sum(3)


class ScriptedRandom:
    def __init__(self, points):
        self.points = iter(points)

    def random(self):
        return next(self.points)


class ScriptedRanks:
    def __init__(self, ranks):
        self.ranks = iter(ranks)

    def randrange(self, stop):
        rank = next(self.ranks)
        if not 0 <= rank < stop:
            raise AssertionError("scripted rank outside range")
        return rank


class TrajectoryTests(unittest.TestCase):
    def test_direct_exact_hits_and_tie_rule(self):
        self.assertEqual(direct_trajectory(3, 1, ScriptedRandom([0, 1 / 3, 2 / 3])), 0)
        self.assertEqual(direct_trajectory(3, 1, ScriptedRandom([1 / 6] * 3)), 2.5)
        # Arrivals near circumference endpoint must see vertex zero across wrap.
        self.assertAlmostEqual(direct_trajectory(3, 1, ScriptedRandom([0.99, 0, 0])), 2.03)

    def test_conditional_triangle_has_zero_variance(self):
        for seed in range(30):
            self.assertAlmostEqual(conditional_trajectory(3, 1, Random(seed)), 17 / 12)
        stats = simulate(3, 1, runs=200, seed=8)
        self.assertAlmostEqual(stats.mean, 17 / 12)
        self.assertEqual(stats.sample_variance, 0)

    def test_exhaustive_rank_paths_against_exact_moments(self):
        # These paths exercise every rank, every removal order, and k=2 -> 1.
        # n=3,m=2 also exercises repeated nonfinal capacity decrements.
        for n, m, expected_mean, expected_variance in (
            (3, 1, Fraction(17, 12), Fraction(0)),
            (4, 1, Fraction(71, 32), Fraction(49, 3072)),
            (3, 2, Fraction(67, 27), Fraction(161, 2916)),
        ):
            count = (2 * n) ** (n * m - 1)
            first = second = 0.0
            for ranks in product(range(2 * n), repeat=n * m - 1):
                value = conditional_trajectory(n, m, ScriptedRanks(ranks))
                first += value / count
                second += value * value / count
            self.assertAlmostEqual(first, float(expected_mean), places=11)
            self.assertAlmostEqual(second - first * first, float(expected_variance), places=11)

    def test_multicapacity_against_simple_independent_gap_sampler(self):
        # Same rank stream, but the oracle recomputes every active gap and
        # selects linearly; it has no linked list or Fenwick tree.
        def reference(n, m, rng):
            residual = [m] * n
            total = 0.0
            while sum(value > 0 for value in residual) > 1:
                active = [i for i, value in enumerate(residual) if value]
                gaps = [(active[(j + 1) % len(active)] - v) % n for j, v in enumerate(active)]
                rank = rng.randrange(2 * n)
                for j, v in enumerate(active):
                    a, b = gaps[j - 1], gaps[j]
                    if rank < a + b:
                        total += (a * a + b * b) / (4 * (a + b))
                        residual[v] -= 1
                        break
                    rank -= a + b
            return total + sum(residual) * n / 4

        for n, m in ((3, 2), (4, 5), (8, 2), (13, 3), (40, 2)):
            for seed in range(20):
                self.assertEqual(conditional_trajectory(n, m, Random(seed)), reference(n, m, Random(seed)))

    def test_direct_against_exhaustive_facility_search(self):
        def reference(n, m, points):
            residual = [m] * n
            total = 0.0
            for u in points:
                x = u * n
                distances = [(min(abs(x - v), n - abs(x - v)), v) for v in range(n) if residual[v]]
                cost, selected = min(distances)
                total += cost
                residual[selected] -= 1
            return total

        for n, m in ((3, 2), (5, 3), (20, 2)):
            points = [Random(551 + i).random() for i in range(n * m)]
            self.assertAlmostEqual(direct_trajectory(n, m, ScriptedRandom(points)), reference(n, m, points))

    def test_model_validation(self):
        for n, m in ((2, 1), (3, 0), (3.0, 1), (3, True)):
            for trajectory in (direct_trajectory, conditional_trajectory):
                with self.assertRaises(ValueError):
                    trajectory(n, m, Random(1))


class StatisticsTests(unittest.TestCase):
    def test_seed_reproducibility_and_global_rng_unchanged(self):
        state = random.getstate()
        first = simulate(5, 2, runs=50, seed=337, method="conditional")
        self.assertEqual(first, simulate(5, 2, runs=50, seed=337, method="conditional"))
        self.assertEqual(state, random.getstate())
        self.assertEqual(first.as_dict()["runs"], 50)

    def test_welford_against_batch_statistics(self):
        import statistics
        for method, trajectory in (("direct", direct_trajectory), ("conditional", conditional_trajectory)):
            rng = Random(112)
            values = [trajectory(4, 2, rng) for _ in range(100)]
            stats = simulate(4, 2, runs=100, seed=112, method=method)
            self.assertAlmostEqual(stats.mean, statistics.mean(values))
            self.assertAlmostEqual(stats.sample_variance, statistics.variance(values))
            self.assertAlmostEqual(stats.ci95_half_width, 1.96 * sqrt(statistics.variance(values) / 100))

    def test_invalid_inputs(self):
        for kwargs in ({"runs": 1}, {"runs": True}, {"seed": "1"}, {"method": "other"}):
            with self.assertRaises(ValueError):
                simulate(3, 1, **kwargs)


if __name__ == "__main__":
    unittest.main()
