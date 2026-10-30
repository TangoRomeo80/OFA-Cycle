"""Linear-time dihedral canonicalization and exact bracelet counts."""

from collections.abc import Sequence


def least_rotation(word: Sequence[int]) -> tuple[int, ...]:
    """Return the lexicographically least rotation using Duval's algorithm.

    The factorization scan of the doubled word takes O(len(word)) symbol
    comparisons. Empty words are accepted and have the empty rotation.
    """
    values = tuple(word)
    n = len(values)
    if n == 0:
        return ()
    doubled = values + values
    i = 0
    start = 0
    while i < n:
        start = i
        j, k = i + 1, i
        while j < 2 * n and doubled[k] <= doubled[j]:
            k = i if doubled[k] < doubled[j] else k + 1
            j += 1
        period = j - k
        while i <= k:
            i += period
    return doubled[start : start + n]


def canonical(state: Sequence[int]) -> tuple[int, ...]:
    """Return the least rotation of a word or its reversal in O(n) time."""
    values = tuple(state)
    return min(least_rotation(values), least_rotation(values[::-1]))


def _positive_integer(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _totient(value: int) -> int:
    result = value
    prime = 2
    while prime * prime <= value:
        if value % prime == 0:
            while value % prime == 0:
                value //= prime
            result -= result // prime
        prime += 1
    if value > 1:
        result -= result // value
    return result


def bracelet_count(n: int, q: int) -> int:
    """Count length-n q-color words modulo rotations and reflections.

    This is the Burnside formula B(n,q) in the paper; zero-capacity states
    are included when q=m+1. The formula also applies to n=1 and n=2.
    """
    _positive_integer(n, "n")
    _positive_integer(q, "q")
    divisors = []
    divisor = 1
    while divisor * divisor <= n:
        if n % divisor == 0:
            divisors.append(divisor)
            if divisor * divisor != n:
                divisors.append(n // divisor)
        divisor += 1
    rotations = sum(_totient(d) * q ** (n // d) for d in divisors)
    if n % 2:
        reflections = n * q ** ((n + 1) // 2)
    else:
        reflections = (n // 2) * (q ** (n // 2) + q ** (n // 2 + 1))
    return (rotations + reflections) // (2 * n)
