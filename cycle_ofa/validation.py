"""Mathematical self-checks available in both source and installed packages."""

from fractions import Fraction
from itertools import product
from random import Random

from .exact import ExactSolver
from .simulation import conditional_trajectory
from .symmetry import bracelet_count, canonical, least_rotation


def _definition(word: tuple[int, ...]) -> tuple[int, ...]:
    if not word:
        return ()
    reverse = word[::-1]
    return min(*(word[i:] + word[:i] for i in range(len(word))),
               *(reverse[i:] + reverse[:i] for i in range(len(word))))


def run_validation(exhaustive: bool = False) -> dict:
    """Raise RuntimeError on a failed check; return a compact result summary."""
    checks = 0

    def require(condition: bool, message: str) -> None:
        nonlocal checks
        if not condition:
            raise RuntimeError(message)
        checks += 1

    for n, m, expected in (
        (3, 1, Fraction(17, 12)), (4, 1, Fraction(71, 32)),
        (5, 1, Fraction(1569, 500)), (3, 2, Fraction(67, 27)),
    ):
        solver = ExactSolver(n)
        require(solver.solve(m) == expected, f"benchmark G_({n},{m}) failed")
        require(solver.state_count == bracelet_count(n, m + 1), "state count failed")
        require(ExactSolver(n, symmetry=False).solve(m) == expected, "raw DP failed")
    for n, q, expected in ((8, 3, 498), (10, 3, 3210), (12, 3, 22913),
                           (8, 4, 4435), (9, 4, 15084), (8, 5, 25395)):
        require(bracelet_count(n, q) == expected, "reported bracelet count failed")
    word_count = 0
    for length in range(0, 9 if exhaustive else 6):
        representatives = set()
        for word in product(range(3), repeat=length):
            actual = canonical(word)
            require(actual == _definition(word), f"canonicalization failed at {word}")
            if word:
                require(least_rotation(word) == min(word[i:] + word[:i] for i in range(length)),
                        f"least rotation failed at {word}")
            representatives.add(actual)
            word_count += 1
        if length:
            require(len(representatives) == bracelet_count(length, 3), "orbit enumeration failed")
    rng = Random(20260822)
    for length in (9, 16, 31, 64, 128):
        for _ in range(40 if exhaustive else 5):
            word = tuple(rng.randrange(8) for _ in range(length))
            require(canonical(word) == _definition(word), "randomized canonicalization failed")
    for seed in range(20):
        require(abs(conditional_trajectory(3, 1, Random(seed)) - 17 / 12) < 1e-12,
                "triangle conditional estimator must be constant")
    return dict(status="passed", checks=checks, exhaustive=exhaustive,
                exhaustive_ternary_max_length=8 if exhaustive else 5,
                enumerated_words=word_count)
