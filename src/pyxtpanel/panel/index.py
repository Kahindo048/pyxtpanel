"""Internal index of a declared panel: rows sorted by (unit, period), units delimited by offsets.

The layout is the one of a compressed sparse row matrix: unit ``i`` owns the sorted rows
``offsets[i]:offsets[i + 1]``.  A kernel can loop over units without any Python object, and the
position of each sorted row in the user's frame is kept in ``source_row`` so that results can be
mapped back to the original order.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyxtpanel._types import IntArray


@dataclass(frozen=True, eq=False)
class PanelIndex:
    """Sorted index of the rows of a panel (arrays are read-only).

    Attributes
    ----------
    unit_labels
        Original identifier of each unit, in sorted order (length ``n_units``).
    offsets
        Unit ``i`` owns the sorted rows ``offsets[i]:offsets[i + 1]`` (length ``n_units + 1``).
    unit
        Unit code (``0..n_units - 1``) of each sorted row.
    rank
        Calendar rank of each sorted row, increasing within a unit.
    source_row
        Position, in the frame the panel was declared from, of each sorted row.
    """

    unit_labels: np.ndarray
    offsets: IntArray
    unit: IntArray
    rank: IntArray
    source_row: IntArray

    def __post_init__(self) -> None:
        for array in (self.unit_labels, self.offsets, self.unit, self.rank, self.source_row):
            array.flags.writeable = False

    @property
    def n_units(self) -> int:
        """Number of units."""
        return int(self.offsets.size - 1)

    @property
    def n_obs(self) -> int:
        """Number of rows."""
        return int(self.rank.size)

    @property
    def counts(self) -> IntArray:
        """Number of observations of each unit (``T_i``)."""
        return np.diff(self.offsets)


def subset_index(index: PanelIndex, keep: np.ndarray) -> PanelIndex:
    """Index restricted to the sorted rows where ``keep`` is true; empty units disappear.

    ``source_row`` still points into the frame the panel was declared from.
    """
    keep = np.asarray(keep, dtype=np.bool_)
    unit = index.unit[keep]
    present, counts = np.unique(unit, return_counts=True)
    renumber = np.zeros(index.n_units, dtype=np.int64)
    renumber[present] = np.arange(present.size)
    return PanelIndex(
        unit_labels=index.unit_labels[present],
        offsets=np.concatenate([[0], np.cumsum(counts)]).astype(np.int64),
        unit=renumber[unit],
        rank=np.ascontiguousarray(index.rank[keep]),
        source_row=np.ascontiguousarray(index.source_row[keep]),
    )


def build_index(
    unit_code: IntArray,
    rank: IntArray,
    positions: IntArray,
    unit_labels: np.ndarray,
    order: IntArray | None = None,
) -> PanelIndex:
    """Sort the usable rows by (unit, rank) and delimit the units.

    Parameters
    ----------
    unit_code
        Unit code of each usable row, in the order of ``positions``.
    rank
        Calendar rank of each usable row.
    positions
        Position of each usable row in the user's frame.
    unit_labels
        Original identifier of each unit code.
    order
        Positions that already put the rows in (unit, rank) order, when the caller has sorted
        them (``check_unique_pairs`` does); sorted here otherwise.

    Returns
    -------
    PanelIndex
        Units that have no usable row are dropped and the codes renumbered.
    """
    if order is None:
        order = np.lexsort((rank, unit_code))
    unit_sorted = unit_code[order]
    per_unit = np.bincount(unit_code, minlength=unit_labels.size)  # no second sort to count
    present = np.flatnonzero(per_unit)
    counts = per_unit[present]
    renumber = np.zeros(unit_labels.size, dtype=np.int64)
    renumber[present] = np.arange(present.size)
    offsets = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    return PanelIndex(
        unit_labels=np.asarray(unit_labels)[present],
        offsets=offsets,
        unit=renumber[unit_sorted],
        rank=np.ascontiguousarray(rank[order]),
        source_row=np.ascontiguousarray(positions[order]).astype(np.int64),
    )


def shift_positions(index: PanelIndex, k: int) -> IntArray:
    """Sorted-row position of the observation ``k`` periods earlier in the same unit.

    Parameters
    ----------
    index
        Sorted index of the panel.
    k
        Periods back (negative: forward).

    Returns
    -------
    numpy.ndarray
        One position per sorted row, ``-1`` where the unit has no observation in that period.
    """
    n = index.n_obs
    if k == 0 or n == 0:
        return np.arange(n, dtype=np.int64)
    rank, unit = index.rank, index.unit
    span = int(np.max(rank)) + 1
    target = rank - k
    inside = (target >= 0) & (target < span)
    key = unit * span + rank  # increasing: rows are sorted by (unit, rank)
    wanted = unit * span + target
    where = np.minimum(np.searchsorted(key, wanted), n - 1)
    found = inside & (key[where] == wanted)
    return np.where(found, where, -1).astype(np.int64)
