"""Gap DP checks against exact integration of the original distance model."""

from fractions import Fraction
from functools import lru_cache
import itertools
import unittest

from cycle_ofa.exact import ExactSolver, StateLimitError, gaps
from cycle_ofa.symmetry import bracelet_count, canonical


@lru_cache(maxsize=None)
def integrated_value(state):
    """Small-instance oracle integrating all edge-distance line arrangements.

    This uses d(e,t;v) from the model directly. It does not use gap lengths,
    selection probabilities, canonicalization, or the gap recurrence.
    """
    if not any(state):
        return Fraction(0)
    n = len(state)

    def delta(a, b):
        difference = abs(a - b)
        return min(difference, n - difference)

    result = Fraction(0)
    for edge in range(n):
        lines = []
        for vertex, residual in enumerate(state):
            if residual:
                lines.append((1, delta(edge, vertex), vertex))
                lines.append((-1, 1 + delta((edge + 1) % n, vertex), vertex))
        cuts = {Fraction(0), Fraction(1)}
        for slope, intercept, _ in lines:
            for other_slope, other_intercept, _ in lines:
                if slope != other_slope:
                    crossing = Fraction(other_intercept - intercept, slope - other_slope)
                    if 0 < crossing < 1:
                        cuts.add(crossing)
        cuts = sorted(cuts)
        for left, right in zip(cuts, cuts[1:]):
            midpoint = (left + right) / 2
            slope, intercept, vertex = min(
                lines, key=lambda line: (line[0] * midpoint + line[1], line[2])
            )
            successor = list(state)
            successor[vertex] -= 1
            cost_integral = Fraction(slope, 2) * (right**2 - left**2) + intercept * (right - left)
            result += (cost_integral + (right - left) * integrated_value(tuple(successor))) / n
    return result


class ExactTests(unittest.TestCase):
    def test_paper_regressions(self):
        for n, m, expected in (
            (3, 1, Fraction(17, 12)),
            (4, 1, Fraction(71, 32)),
            (5, 1, Fraction(1569, 500)),
            (3, 2, Fraction(67, 27)),
        ):
            solver = ExactSolver(n)
            self.assertEqual(solver.solve(m), expected)
            self.assertEqual(solver.state_count, bracelet_count(n, m + 1))

    def test_original_distance_integral_oracle(self):
        for n, q in ((3, 3), (4, 2), (5, 2)):
            solver = ExactSolver(n)
            for state in itertools.product(range(q), repeat=n):
                self.assertEqual(solver.value(state), integrated_value(state), state)

    def test_raw_and_canonical_match(self):
        for n, m in ((3, 3), (4, 2), (5, 1)):
            raw, reduced = ExactSolver(n, symmetry=False), ExactSolver(n)
            self.assertEqual(raw.solve(m), reduced.solve(m))
            self.assertEqual(raw.state_count, (m + 1) ** n)
            self.assertEqual(reduced.state_count, bracelet_count(n, m + 1))
            for state in itertools.product(range(m + 1), repeat=n):
                self.assertEqual(raw.value(state), reduced.value(state))

    def test_heterogeneous_invariance(self):
        state = (2, 0, 1, 3, 0)
        raw = ExactSolver(5, symmetry=False)
        self.assertEqual(raw.value(state), raw.value(canonical(state)))

    def test_gap_boundary_cases(self):
        self.assertEqual(gaps((0, 0, 0)), ())
        self.assertEqual(gaps((0, 2, 0, 0, 0)), (5,))
        self.assertEqual(gaps((0, 1, 0, 0, 2)), (3, 2))
        self.assertEqual(gaps((1, 1, 1)), (1, 1, 1))
        solver = ExactSolver(5)
        for residual in range(1, 6):
            self.assertEqual(solver.value((0, 0, residual, 0, 0)), Fraction(5 * residual, 4))

    def test_deep_state_does_not_use_python_recursion(self):
        self.assertEqual(ExactSolver(3).value((0, 0, 1200)), Fraction(900))

    def test_state_budget(self):
        solver = ExactSolver(4, max_states=3)
        with self.assertRaisesRegex(StateLimitError, "max_states=3"):
            solver.solve(2)
        self.assertLessEqual(solver.state_count, 3)
        self.assertEqual(ExactSolver(3, max_states=1).value((0, 0, 0)), 0)

    def test_invalid_inputs(self):
        for n in (2, True, 3.5):
            with self.assertRaises(ValueError):
                ExactSolver(n)
        for state in ((1, 1), (1, -1, 0), (True, 0, 1), (0.5, 0, 1), (1, 1, 1, 1)):
            with self.assertRaises(ValueError):
                ExactSolver(3).value(state)
        for capacity in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                ExactSolver(3).solve(capacity)
        for limit in (0, -1, 2.5, True):
            with self.assertRaises(ValueError):
                ExactSolver(3, max_states=limit)


if __name__ == "__main__":
    unittest.main()
