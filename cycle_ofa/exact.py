"""Exact rational gap recurrence for greedy facility assignment on cycles."""

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

from .symmetry import canonical


State = tuple[int, ...]


def _validate_state(state: Sequence[int]) -> State:
    values = tuple(state)
    if len(values) < 3:
        raise ValueError("a cycle state must have at least three facilities")
    if any(isinstance(r, bool) or not isinstance(r, int) or r < 0 for r in values):
        raise ValueError("residual capacities must be nonnegative integers")
    return values


def _geometry(state: State) -> tuple[tuple[int, ...], tuple[int, ...]]:
    active = tuple(v for v, residual in enumerate(state) if residual)
    if not active:
        return (), ()
    n = len(state)
    if len(active) == 1:
        return active, (n,)
    return active, tuple(
        (active[(i + 1) % len(active)] - vertex) % n
        for i, vertex in enumerate(active)
    )


def gaps(state: Sequence[int]) -> tuple[int, ...]:
    """Clockwise gaps starting at the first active vertex in coordinate order.

    The terminal state has no gaps; a sole active facility has self-gap n.
    """
    return _geometry(_validate_state(state))[1]


class StateLimitError(RuntimeError):
    """The configured exact-state budget has been reached."""


@dataclass
class _Frame:
    state: State
    transitions: tuple[tuple[State, int], ...]
    weighted_value: Fraction
    next_index: int = 0


class ExactSolver:
    """Memoized Fraction-valued DP, optionally keyed by dihedral orbits.

    Values are computed by an explicit depth-first stack, so a large residual
    sum does not modify or depend on Python's recursion limit. Each solver
    keeps its cache across calls. ``state_count`` includes the terminal state.
    The optional ``max_states`` limits cached states plus pending stack frames.
    Rational bit complexity and memory still grow exponentially with n.
    """

    def __init__(self, n: int, symmetry: bool = True, max_states: int | None = None):
        if isinstance(n, bool) or not isinstance(n, int) or n < 3:
            raise ValueError("n must be an integer of at least 3")
        if not isinstance(symmetry, bool):
            raise ValueError("symmetry must be a boolean")
        if max_states is not None and (
            isinstance(max_states, bool)
            or not isinstance(max_states, int)
            or max_states < 1
        ):
            raise ValueError("max_states must be a positive integer or None")
        self.n = n
        self.symmetry = symmetry
        self.max_states = max_states
        self._values: dict[State, Fraction] = {(0,) * n: Fraction(0)}

    @property
    def state_count(self) -> int:
        """Number of completed memo entries, including the terminal state."""
        return len(self._values)

    def _key(self, state: State) -> State:
        return canonical(state) if self.symmetry else state

    def _frame(self, state: State) -> _Frame:
        active, gap_lengths = _geometry(state)
        # Group successors in the same orbit, retaining their total weight.
        weights: dict[State, int] = {}
        for i, vertex in enumerate(active):
            successor = list(state)
            successor[vertex] -= 1
            key = self._key(tuple(successor))
            weight = gap_lengths[i - 1] + gap_lengths[i]
            weights[key] = weights.get(key, 0) + weight
        # Multiplying the recurrence by 2n avoids a Fraction per transition.
        immediate_scaled = Fraction(sum(g * g for g in gap_lengths), 2)
        return _Frame(state, tuple(weights.items()), immediate_scaled)

    def value(self, state: Sequence[int]) -> Fraction:
        """Return the exact expected remaining total cost from ``state``."""
        values = _validate_state(state)
        if len(values) != self.n:
            raise ValueError(f"state must have length {self.n}")
        root = self._key(values)
        if root in self._values:
            return self._values[root]
        stack: list[_Frame] = []

        def push(key: State) -> None:
            if self.max_states is not None and len(self._values) + len(stack) >= self.max_states:
                raise StateLimitError(
                    f"exact DP reached max_states={self.max_states}; "
                    "increase the limit or use simulation"
                )
            stack.append(self._frame(key))

        push(root)
        while stack:
            frame = stack[-1]
            if frame.next_index == len(frame.transitions):
                self._values[frame.state] = frame.weighted_value / (2 * self.n)
                stack.pop()
                continue
            successor, weight = frame.transitions[frame.next_index]
            if successor not in self._values:
                push(successor)
                continue
            frame.weighted_value += weight * self._values[successor]
            frame.next_index += 1
        return self._values[root]

    def solve(self, m: int) -> Fraction:
        """Return G_(n,m), the exact total cost with uniform capacity m."""
        if isinstance(m, bool) or not isinstance(m, int) or m < 1:
            raise ValueError("m must be a positive integer")
        return self.value((m,) * self.n)
