"""A small, dependency-free interval index.

This is a behaviour-compatible replacement for the subset of the ``interlap``
package used by Liftoff. ``interlap`` is only distributed as a source archive,
which cannot be installed in Pyodide, and the functionality required here is
only a few lines of code.

Intervals are closed (both ends inclusive) and are kept sorted by
``(start, end, key)``. Queries return matches in that order, which matters
because the order of reported overlaps influences which features are selected
for re-mapping.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Iterable, Iterator
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Interval[T]:
    """A closed interval with an attached payload.

    Attributes
    ----------
    start : int
        Inclusive start coordinate.
    end : int
        Inclusive end coordinate.
    key : str
        Sort tie-breaker; for lifted features this is the copy ID.
    value : T
        Arbitrary payload returned with query results.
    """

    start: int
    end: int
    key: str
    value: T


class IntervalIndex[T]:
    """Sorted collection of intervals supporting overlap queries.

    Parameters
    ----------
    intervals : iterable of Interval, optional
        Initial intervals to index.
    """

    def __init__(self, intervals: Iterable[Interval[T]] = ()) -> None:
        self._intervals: list[Interval[T]] = sorted(
            intervals, key=lambda iv: (iv.start, iv.end, iv.key)
        )
        self._starts = [iv.start for iv in self._intervals]
        # Longest interval length bounds how far left of a query an
        # overlapping interval can start.
        self._max_length = max((iv.end - iv.start + 1 for iv in self._intervals), default=0)

    def __len__(self) -> int:
        return len(self._intervals)

    def __iter__(self) -> Iterator[Interval[T]]:
        return iter(self._intervals)

    def find(self, start: int, end: int) -> list[Interval[T]]:
        """Return all intervals overlapping ``[start, end]``.

        Parameters
        ----------
        start : int
            Inclusive query start.
        end : int
            Inclusive query end.

        Returns
        -------
        list of Interval
            Overlapping intervals in index order.
        """
        # Only intervals whose start lies in [start - max_length, end] can
        # possibly overlap the query, so restrict the scan to that slice.
        lo = bisect_left(self._starts, start - self._max_length)
        hi = bisect_right(self._starts, end)
        return [iv for iv in self._intervals[lo:hi] if iv.start <= end and iv.end >= start]
