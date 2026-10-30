"""Independent capacity-state assignment oracle and paired-statistic checks."""

from fractions import Fraction
from functools import lru_cache
from math import sqrt
from random import Random
import statistics
import unittest

from cycle_ofa.offline import compare_offline, direct_cost, offline_cost


def assignment_oracle(n, m, arrivals):
    """Exhaust all capacity-feasible assignments, using exact rational costs."""
    points = tuple(Fraction(point) for point in arrivals)

    @lru_cache(maxsize=None)
    def optimum(residual):
        index = n * m - sum(residual)
        if index == len(points):
            return Fraction(0)
        candidates = []
        for vertex, available in enumerate(residual):
            if available:
                distance = abs(points[index] - vertex)
                distance = min(distance, n - distance)
                successor = list(residual)
                successor[vertex] -= 1
                candidates.append(distance + optimum(tuple(successor)))
        return min(candidates)

    return optimum((m,) * n)


def greedy_oracle(n, m, arrivals):
    residual = [m] * n
    cost = 0.0
    for point in arrivals:
        options = []
        for vertex, available in enumerate(residual):
            if available:
                distance = abs(point - vertex)
                options.append((min(distance, n - distance), vertex))
        distance, selected = min(options)
        cost += distance
        residual[selected] -= 1
    return cost


class OfflineTests(unittest.TestCase):
    def test_small_independent_assignment_oracle(self):
        rng = Random(731)
        for n, m in ((3, 1), (3, 2), (4, 1), (4, 2), (5, 1)):
            for _ in range(30):
                # Dyadic coordinates are exactly representable as floats.
                points = [rng.randrange(n * 16) / 16 for _ in range(n * m)]
                self.assertAlmostEqual(offline_cost(n, m, points), float(assignment_oracle(n, m, points)), places=12)

    def test_repeated_zero_antipodal_and_wraparound(self):
        cases = (
            (3, 1, [0, 1, 2]),
            (3, 2, [0, 0, 1, 1, 2, 2]),
            (4, 1, [0, 0, 0, 0]),
            (4, 1, [2, 2, 2, 2]),
            (4, 1, [0.5, 1.5, 2.5, 3.5]),
            (3, 2, [0, 0, 1.5, 1.5, 2.9375, 2.9375]),
            (5, 1, [0.0625, 4.9375, 4.875, 0.125, 2.5]),
        )
        for n, m, points in cases:
            self.assertAlmostEqual(offline_cost(n, m, points), float(assignment_oracle(n, m, points)), places=12)
            self.assertGreaterEqual(direct_cost(n, m, points) + 1e-12, offline_cost(n, m, points))
        self.assertEqual(offline_cost(3, 1, [0, 1, 2]), 0)
        self.assertEqual(offline_cost(4, 1, [0, 0, 0, 0]), 4)

    def test_greedy_independent_nearest_facility_scan(self):
        rng = Random(2)
        for n, m in ((3, 1), (5, 2), (8, 4), (20, 5)):
            for _ in range(40):
                points = [rng.random() * n for _ in range(n * m)]
                self.assertAlmostEqual(direct_cost(n, m, points), greedy_oracle(n, m, points), places=11)
                self.assertGreaterEqual(direct_cost(n, m, points) + 1e-10, offline_cost(n, m, points))

    def test_order_and_rotation_properties(self):
        points = [0.125, 0.375, 1.25, 2.75, 3.5, 3.75, 0.75, 2.125]
        expected = offline_cost(4, 2, points)
        self.assertEqual(offline_cost(4, 2, list(reversed(points))), expected)
        for shift in range(4):
            shifted = [(point + shift) % 4 for point in points]
            self.assertEqual(offline_cost(4, 2, shifted), expected)
        self.assertEqual(offline_cost(4, 2, [(-point) % 4 for point in points]), expected)

    def test_paired_statistics_against_saved_sample_calculation(self):
        n, m, runs, seed = 4, 2, 120, 214
        stats = compare_offline(n, m, runs, seed)
        rng = Random(seed)
        pairs = []
        for _ in range(runs):
            points = [rng.random() * n for _ in range(n * m)]
            pairs.append((direct_cost(n, m, points), offline_cost(n, m, points)))
        greedy, offline = map(list, zip(*pairs))
        differences = [g - o for g, o in pairs]
        ratio = statistics.mean(greedy) / statistics.mean(offline)
        z = [g - ratio * o for g, o in pairs]
        self.assertAlmostEqual(stats.greedy_mean, statistics.mean(greedy), places=12)
        self.assertAlmostEqual(stats.offline_mean, statistics.mean(offline), places=12)
        self.assertAlmostEqual(stats.greedy_variance, statistics.variance(greedy), places=12)
        self.assertAlmostEqual(stats.offline_variance, statistics.variance(offline), places=12)
        covariance = sum((g - statistics.mean(greedy)) * (o - statistics.mean(offline)) for g, o in pairs) / (runs - 1)
        self.assertAlmostEqual(stats.covariance, covariance, places=12)
        self.assertAlmostEqual(stats.excess_mean, statistics.mean(differences), places=12)
        self.assertAlmostEqual(stats.excess_variance, statistics.variance(differences), places=12)
        self.assertAlmostEqual(stats.excess_ci95_half_width, 1.96 * statistics.stdev(differences) / sqrt(runs), places=12)
        self.assertAlmostEqual(stats.greedy_ci95_half_width, 1.96 * statistics.stdev(greedy) / sqrt(runs), places=12)
        self.assertAlmostEqual(stats.offline_ci95_half_width, 1.96 * statistics.stdev(offline) / sqrt(runs), places=12)
        self.assertAlmostEqual(stats.ratio_of_means, ratio, places=12)
        self.assertAlmostEqual(stats.ratio_ci95_half_width, 1.96 * statistics.stdev(z) / (sqrt(runs) * statistics.mean(offline)), places=12)
        self.assertEqual(stats, compare_offline(n, m, runs, seed))

    def test_private_random_generator(self):
        import random
        original = random.getstate()
        compare_offline(3, 1, 3, 1)
        self.assertEqual(random.getstate(), original)

    def test_larger_instance(self):
        rng = Random(18)
        points = [rng.random() * 200 for _ in range(2000)]
        self.assertGreaterEqual(direct_cost(200, 10, points) + 1e-9, offline_cost(200, 10, points))

    def test_invalid_inputs(self):
        for function in (offline_cost, direct_cost):
            for n, m, points in (
                (2, 1, [0, 1]), (3, 0, []), (True, 1, [0]),
                (3, 1, [0, 1]), (3, 1, [0, 1, 3]), (3, 1, [0, 1, -0.1]),
                (3, 1, [0, 1, float('nan')]), (3, 1, [0, 1, float('inf')]),
                (3, 1, [0, 1, True]), (3, 1, [0, 1, '2']), (3, 1, None),
            ):
                with self.assertRaises(ValueError):
                    function(n, m, points)
        for runs, seed in ((1, 1), (2.5, 1), (True, 1), (2, True), (2, 'seed')):
            with self.assertRaises(ValueError):
                compare_offline(3, 1, runs, seed)


if __name__ == '__main__':
    unittest.main()
