"""Decompose direct-cost variance into history and positional contributions.

Given the complete facility-selection sequence, arrival coordinates are
independent within their known Voronoi regions. Consequently, if ``Y`` is the
conditional cost and ``S`` is the sum of conditional one-step variances,
``Var(C) = Var(Y) + E[S]``. This identity explains measured variance ratios;
it does not assert an asymptotic limit as the number of facilities increases.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import sqrt
from random import Random

from .simulation import FenwickTree, _integer, _validate_model


def positional_variance(left_gap: int, right_gap: int) -> float:
    """Variance of arrival distance conditional on the selected facility.

    For adjacent gaps ``a,b``, the conditional second moment is
    ``(a**3+b**3)/(12*(a+b))`` and the mean is
    ``(a*a+b*b)/(4*(a+b))``. For a singleton use ``a=b=n``.
    """
    _integer(left_gap, "left_gap", 1)
    _integer(right_gap, "right_gap", 1)
    a, b = left_gap, right_gap
    mean = (a * a + b * b) / (4 * (a + b))
    return (a * a * a + b * b * b) / (12 * (a + b)) - mean * mean


def conditional_diagnostic_trajectory(n: int, m: int, rng: Random) -> tuple[float, float]:
    """Return conditional cost ``Y`` and removed positional variance ``S``.

    Sampling and arithmetic for ``Y`` follow ``conditional_trajectory`` in
    the same order, so an identically seeded RNG gives the same trajectory
    cost. Accumulating ``S`` consumes no additional randomness. The Fenwick
    tree and circular neighbor list retain O(mn log n) time/O(n) memory.
    """
    _validate_model(n, m)
    residual = [m] * n
    previous = [n - 1] + list(range(n - 1))
    following = list(range(1, n)) + [0]
    clockwise_gap = [1] * n
    tree = FenwickTree([2] * n)
    select = tree.find_by_order
    set_weight = tree.set_weight
    rank = rng.randrange
    circumference_twice = 2 * n
    active_count = n
    total_cost = 0.0
    removed_variance = 0.0
    while active_count > 1:
        selected = select(rank(circumference_twice))
        left = previous[selected]
        right = following[selected]
        a = clockwise_gap[left]
        b = clockwise_gap[selected]
        conditional_mean = (a * a + b * b) / (4 * (a + b))
        total_cost += conditional_mean
        removed_variance += (a * a * a + b * b * b) / (12 * (a + b)) - conditional_mean * conditional_mean
        residual[selected] -= 1
        if residual[selected]:
            continue
        set_weight(selected, 0)
        active_count -= 1
        if active_count == 1:
            # The remaining cost is deterministic after conditioning, but
            # the original direct process retains n^2/48 positional variance
            # per remaining arrival. Do not omit this singleton contribution.
            return (
                total_cost + residual[left] * n / 4,
                removed_variance + residual[left] * n * n / 48,
            )
        clockwise_gap[left] = a + b
        following[left] = right
        previous[right] = left
        set_weight(left, clockwise_gap[previous[left]] + clockwise_gap[left])
        set_weight(right, clockwise_gap[left] + clockwise_gap[right])
    raise AssertionError("a valid cycle must reach the singleton return")


@dataclass(frozen=True)
class VarianceDiagnostics:
    """Monte Carlo estimates for the law-of-total-variance decomposition.

    ``conditional_variance`` and ``removed_variance_sample_variance`` use
    the unbiased sample-variance divisor. ``variance_ratio`` estimates the
    population ratio by ``1 + mean(S)/sample_variance(Y)``; it is ``None``
    when the conditional sample variance vanishes. This diagnostic ratio
    differs from a ratio of variances from two independent simulation runs.
    """

    n: int
    m: int
    runs: int
    seed: int
    conditional_mean: float
    conditional_variance: float
    mean_removed_variance: float
    removed_variance_sample_variance: float
    removed_variance_ci95_half_width: float
    implied_direct_variance: float
    variance_ratio: float | None

    def as_dict(self) -> dict[str, int | float | None]:
        return asdict(self)


def diagnose_variance(
    n: int, m: int, runs: int = 10_000, seed: int = 20260822,
) -> VarianceDiagnostics:
    """Estimate history variance and mean positional noise on the same paths.

    Both moments are accumulated with Welford updates in constant memory.
    The uncertainty attached to ``mean_removed_variance`` uses the same
    normal approximation as other package Monte Carlo intervals. No
    confidence interval or asymptotic convergence claim is attached to the
    variance-ratio point estimate.
    """
    _validate_model(n, m)
    _integer(runs, "runs", 2)
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    rng = Random(seed)
    mean_y = mean_s = m2_y = m2_s = 0.0
    for count in range(1, runs + 1):
        y, s = conditional_diagnostic_trajectory(n, m, rng)
        delta_y, delta_s = y - mean_y, s - mean_s
        mean_y += delta_y / count
        mean_s += delta_s / count
        m2_y += delta_y * (y - mean_y)
        m2_s += delta_s * (s - mean_s)
    variance_y = max(0.0, m2_y / (runs - 1))
    variance_s = max(0.0, m2_s / (runs - 1))
    return VarianceDiagnostics(
        n=n, m=m, runs=runs, seed=seed, conditional_mean=mean_y,
        conditional_variance=variance_y, mean_removed_variance=mean_s,
        removed_variance_sample_variance=variance_s,
        removed_variance_ci95_half_width=1.96 * sqrt(variance_s / runs),
        implied_direct_variance=variance_y + mean_s,
        variance_ratio=1 + mean_s / variance_y if variance_y > 0 else None,
    )
