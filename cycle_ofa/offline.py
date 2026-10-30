"""Offline optimal transport and paired comparisons on the same arrivals.

The offline algorithm uses circular Earth Mover's Distance: the integral of
the absolute cumulative supply--demand imbalance after subtracting its
length-weighted median. See Rabin, Delon, and Gousseau, "Transportation
Distances on the Circle", Journal of Mathematical Imaging and Vision 41,
147--167 (2011), Section 2.2.2, equation (10), DOI 10.1007/s10851-011-0284-0.
https://perso.telecom-paristech.fr/gousseau/circle.pdf

Costs are evaluated in floating-point arithmetic. The algorithm solves the
offline optimum, rather than a heuristic assignment, but these numerical
values should not be described as exact rational results.
"""

from __future__ import annotations

from bisect import bisect_left
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from math import fsum, isfinite, sqrt
from numbers import Real
from random import Random


def _integer(value: int, name: str, minimum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _model(n: int, m: int) -> None:
    _integer(n, "n", 3)
    _integer(m, "m", 1)


def _arrivals(n: int, m: int, arrivals: Sequence[float]) -> list[float]:
    _model(n, m)
    try:
        values = list(arrivals)
    except TypeError as error:
        raise ValueError("arrivals must be a sequence of real coordinates") from error
    if len(values) != n * m:
        raise ValueError(f"expected exactly n*m={n * m} arrivals")
    converted = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, Real):
            raise ValueError("arrival coordinates must be finite real numbers in [0,n)")
        try:
            coordinate = float(value)
        except (OverflowError, ValueError) as error:
            raise ValueError("arrival coordinates must be finite real numbers in [0,n)") from error
        if not isfinite(coordinate) or not 0 <= coordinate < n:
            raise ValueError("arrival coordinates must be finite real numbers in [0,n)")
        converted.append(coordinate)
    return converted


def _offline_cost(n: int, m: int, arrivals: list[float]) -> float:
    events = [(coordinate, 1) for coordinate in arrivals]
    events.extend((float(vertex), -m) for vertex in range(n))
    events.sort()
    intervals: list[tuple[int, float]] = []
    weights: dict[int, float] = {}
    imbalance = 0
    index = 0
    while index < len(events):
        position = events[index][0]
        # Combine coincident supply and demand before assigning a gap's flow.
        while index < len(events) and events[index][0] == position:
            imbalance += events[index][1]
            index += 1
        next_position = events[index][0] if index < len(events) else float(n)
        length = next_position - position
        if length > 0:
            intervals.append((imbalance, length))
            weights[imbalance] = weights.get(imbalance, 0.0) + length
    # Vertex 0 is always an event, so intervals partition the whole [0,n).
    if imbalance != 0:
        raise AssertionError("offline supply and capacity must balance")
    cumulative_weight = 0.0
    median = 0
    for value, weight in sorted(weights.items()):
        median = value
        cumulative_weight += weight
        if cumulative_weight >= n / 2:
            break
    return fsum(length * abs(value - median) for value, length in intervals)


def offline_cost(n: int, m: int, arrivals: Sequence[float]) -> float:
    """Return minimum total assignment cost with advance arrival knowledge.

    All n*m arrivals are assigned, with exactly m assigned to each facility.
    Coordinates must belong to [0,n); coincident points are supported. Sorting
    events and weighted-median levels takes O((n*m+n) log(n*m+n)) time and
    O(n*m+n) memory. No coordinate binning or quadrature is performed.
    """
    return _offline_cost(n, m, _arrivals(n, m, arrivals))


def _direct_cost(n: int, m: int, arrivals: list[float]) -> float:
    active = list(range(n))
    residual = [m] * n
    costs = []
    for position in arrivals:
        right_index = bisect_left(active, position)
        if right_index == len(active):
            right_index = 0
        left_index = right_index - 1
        left, right = active[left_index], active[right_index]
        left_distance = (position - left) % n
        right_distance = (right - position) % n
        # Ties select the counterclockwise adjacent facility, as in simulation.
        if left_distance <= right_distance:
            selected, selected_index, cost = left, left_index, left_distance
        else:
            selected, selected_index, cost = right, right_index, right_distance
        costs.append(cost)
        residual[selected] -= 1
        if residual[selected] == 0:
            del active[selected_index]
    return fsum(costs)


def direct_cost(n: int, m: int, arrivals: Sequence[float]) -> float:
    """Run nearest-available online greedy on an explicit ordered arrival list."""
    return _direct_cost(n, m, _arrivals(n, m, arrivals))


@dataclass(frozen=True)
class OfflineComparison:
    """Paired Monte Carlo estimates of total costs, excess, and mean ratio.

    Variances/covariance are unbiased sample estimates using runs-1. Mean CIs
    use 1.96 standard errors. The ratio estimates E[greedy]/E[offline], with
    a paired delta-method CI; it is not a mean of trajectory ratios or a
    worst-case competitive ratio. A zero offline sample mean has no ratio.
    """

    n: int
    m: int
    runs: int
    seed: int
    greedy_mean: float
    offline_mean: float
    greedy_variance: float
    offline_variance: float
    covariance: float
    greedy_ci95_half_width: float
    offline_ci95_half_width: float
    excess_mean: float
    excess_variance: float
    excess_ci95_half_width: float
    ratio_of_means: float | None
    ratio_ci95_half_width: float | None

    def as_dict(self) -> dict:
        return asdict(self)


def compare_offline(
    n: int,
    m: int,
    runs: int = 1_000,
    seed: int = 20260822,
) -> OfflineComparison:
    """Compare greedy and offline costs on each of ``runs`` shared samples.

    ``seed`` is the actual private Random generator seed. Covariance is
    retained for uncertainty in the paired excess and the ratio of means.
    No conditional-estimator samples are used as realized greedy costs.
    """
    _model(n, m)
    _integer(runs, "runs", 2)
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    rng = Random(seed)
    greedy_mean = offline_mean = excess_mean = 0.0
    greedy_m2 = offline_m2 = excess_m2 = cross_m2 = 0.0
    for count in range(1, runs + 1):
        arrivals = [rng.random() * n for _ in range(n * m)]
        greedy = _direct_cost(n, m, arrivals)
        optimum = _offline_cost(n, m, arrivals)
        excess = greedy - optimum
        greedy_delta = greedy - greedy_mean
        offline_delta = optimum - offline_mean
        excess_delta = excess - excess_mean
        greedy_mean += greedy_delta / count
        offline_mean += offline_delta / count
        excess_mean += excess_delta / count
        greedy_m2 += greedy_delta * (greedy - greedy_mean)
        offline_m2 += offline_delta * (optimum - offline_mean)
        excess_m2 += excess_delta * (excess - excess_mean)
        cross_m2 += greedy_delta * (optimum - offline_mean)
    greedy_variance = max(0.0, greedy_m2 / (runs - 1))
    offline_variance = max(0.0, offline_m2 / (runs - 1))
    excess_variance = max(0.0, excess_m2 / (runs - 1))
    covariance = cross_m2 / (runs - 1)
    ratio = ratio_half_width = None
    if offline_mean > 0:
        ratio = greedy_mean / offline_mean
        centered_variance = max(
            0.0, greedy_variance + ratio * ratio * offline_variance - 2 * ratio * covariance,
        )
        ratio_half_width = 1.96 * sqrt(centered_variance / runs) / offline_mean
    return OfflineComparison(
        n=n, m=m, runs=runs, seed=seed,
        greedy_mean=greedy_mean, offline_mean=offline_mean,
        greedy_variance=greedy_variance, offline_variance=offline_variance,
        covariance=covariance,
        greedy_ci95_half_width=1.96 * sqrt(greedy_variance / runs),
        offline_ci95_half_width=1.96 * sqrt(offline_variance / runs),
        excess_mean=excess_mean, excess_variance=excess_variance,
        excess_ci95_half_width=1.96 * sqrt(excess_variance / runs),
        ratio_of_means=ratio, ratio_ci95_half_width=ratio_half_width,
    )
