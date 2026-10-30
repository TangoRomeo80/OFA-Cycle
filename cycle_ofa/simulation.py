"""Independent direct and Rao--Blackwell simulations on unit-edge cycles.

All randomness belongs to an explicitly supplied ``random.Random`` instance.
Conditional sampling uses integer weights, so categorical probabilities have
no cumulative floating-point rounding. Reported costs and Monte Carlo moments
use binary floating-point arithmetic; they are not exact rational results.
"""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import asdict, dataclass
from math import sqrt
from random import Random
from typing import Iterable


def _integer(value: int, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _validate_model(n: int, m: int) -> None:
    _integer(n, "n", 3)
    _integer(m, "m", 1)


class FenwickTree:
    """Nonnegative integer weights with logarithmic updates and selection.

    ``find_by_order(order)`` returns the index containing zero-based integer
    rank ``order`` in the expanded multiset of weights. Thus weights
    ``[2, 0, 3]`` map ranks ``0,1`` to index 0 and ``2,3,4`` to index 2.
    ``prefix_sum(stop)`` sums the half-open interval ``[0, stop)``.
    """

    def __init__(self, weights: Iterable[int]) -> None:
        values = list(weights)
        if not values:
            raise ValueError("weights must contain at least one entry")
        for value in values:
            _integer(value, "weight", 0)
        self._weights = values
        self._size = len(values)
        self._tree = [0] + values.copy()
        for index in range(1, self._size + 1):
            parent = index + (index & -index)
            if parent <= self._size:
                self._tree[parent] += self._tree[index]
        self._total = sum(values)
        self._top_bit = 1 << (self._size.bit_length() - 1)

    @property
    def weights(self) -> tuple[int, ...]:
        return tuple(self._weights)

    @property
    def total(self) -> int:
        return self._total

    def prefix_sum(self, stop: int) -> int:
        _integer(stop, "stop", 0)
        if stop > self._size:
            raise ValueError("stop exceeds the number of weights")
        value = 0
        while stop:
            value += self._tree[stop]
            stop -= stop & -stop
        return value

    def set_weight(self, index: int, weight: int) -> None:
        _integer(index, "index", 0)
        if index >= self._size:
            raise ValueError("index exceeds the number of weights")
        _integer(weight, "weight", 0)
        delta = weight - self._weights[index]
        self._weights[index] = weight
        self._total += delta
        index += 1
        while index <= self._size:
            self._tree[index] += delta
            index += index & -index

    def add(self, index: int, delta: int) -> None:
        _integer(index, "index", 0)
        if index >= self._size:
            raise ValueError("index exceeds the number of weights")
        if isinstance(delta, bool) or not isinstance(delta, int):
            raise ValueError("delta must be an integer")
        self.set_weight(index, self._weights[index] + delta)

    def find_by_order(self, order: int) -> int:
        _integer(order, "order", 0)
        if order >= self._total:
            raise ValueError("order must be smaller than the total weight")
        index = 0
        bit = self._top_bit
        tree = self._tree
        while bit:
            candidate = index + bit
            if candidate <= self._size and tree[candidate] <= order:
                order -= tree[candidate]
                index = candidate
            bit >>= 1
        return index


def direct_trajectory(n: int, m: int, rng: Random) -> float:
    """Return realized total cost for ``mn`` uniform continuous arrivals.

    Facilities are stored in sorted order. Binary search locates the two
    adjacent active facilities; the nearer is chosen by clockwise/counter-
    clockwise distance. Exact ties select the counterclockwise facility.
    List deletion costs O(n), yielding O(mn log n + n^2) time and O(n) space.
    No gap-selection probabilities or conditional costs are used here.
    """
    _validate_model(n, m)
    active = list(range(n))
    residual = [m] * n
    total_cost = 0.0
    uniform = rng.random
    for _ in range(n * m):
        position = uniform() * n
        right_index = bisect_left(active, position)
        if right_index == len(active):
            right_index = 0
        left_index = right_index - 1
        left = active[left_index]
        right = active[right_index]
        left_distance = (position - left) % n
        right_distance = (right - position) % n
        if left_distance <= right_distance:
            selected, selected_index = left, left_index
            total_cost += left_distance
        else:
            selected, selected_index = right, right_index
            total_cost += right_distance
        residual[selected] -= 1
        if residual[selected] == 0:
            del active[selected_index]
    return total_cost


def conditional_trajectory(n: int, m: int, rng: Random) -> float:
    """Return the conditional gap estimator in O(mn log n) time/O(n) space.

    A Fenwick tree samples facility ``i`` with integer weight equal to the
    sum of its two adjacent gaps. A circular linked list updates neighbors
    only when a facility exhausts its final slot. Once one facility remains,
    its residual contribution is deterministic and is added without RNG draws.
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
    while active_count > 1:
        selected = select(rank(circumference_twice))
        left = previous[selected]
        right = following[selected]
        left_gap = clockwise_gap[left]
        right_gap = clockwise_gap[selected]
        total_cost += (left_gap * left_gap + right_gap * right_gap) / (
            4 * (left_gap + right_gap)
        )
        residual[selected] -= 1
        if residual[selected]:
            continue
        set_weight(selected, 0)
        active_count -= 1
        if active_count == 1:
            # With two active sites, left and right are the same survivor.
            previous[left] = following[left] = left
            clockwise_gap[left] = n
            set_weight(left, circumference_twice)
            return total_cost + residual[left] * n / 4
        clockwise_gap[left] = left_gap + right_gap
        following[left] = right
        previous[right] = left
        set_weight(left, clockwise_gap[previous[left]] + clockwise_gap[left])
        set_weight(right, clockwise_gap[left] + clockwise_gap[right])
    raise AssertionError("a valid cycle must reach the singleton return")


@dataclass(frozen=True)
class SimulationStats:
    """Streaming Monte Carlo summary; variance uses the ``runs - 1`` divisor."""

    runs: int
    mean: float
    sample_variance: float
    ci95_half_width: float
    seed: int
    method: str
    n: int
    m: int

    def as_dict(self) -> dict[str, int | float | str]:
        return asdict(self)


def simulate(
    n: int,
    m: int,
    runs: int = 10_000,
    seed: int = 20260822,
    method: str = "conditional",
) -> SimulationStats:
    """Estimate total cost using a private RNG and Welford's online moments.

    The interval half-width is ``1.96*s/sqrt(runs)``, matching the paper's
    normal approximation. It is a Monte Carlo interval, not a rigorous bound.
    Changing the trajectory implementation or Python RNG implementation may
    change the exact finite sample even when the base seed is held fixed.
    """
    _validate_model(n, m)
    _integer(runs, "runs", 2)
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    if method == "direct":
        trajectory = direct_trajectory
    elif method == "conditional":
        trajectory = conditional_trajectory
    else:
        raise ValueError("method must be 'direct' or 'conditional'")
    rng = Random(seed)
    mean = 0.0
    squared_deviations = 0.0
    for count in range(1, runs + 1):
        value = trajectory(n, m, rng)
        difference = value - mean
        mean += difference / count
        squared_deviations += difference * (value - mean)
    variance = max(0.0, squared_deviations / (runs - 1))
    return SimulationStats(
        runs=runs,
        mean=mean,
        sample_variance=variance,
        ci95_half_width=1.96 * sqrt(variance / runs),
        seed=seed,
        method=method,
        n=n,
        m=m,
    )
