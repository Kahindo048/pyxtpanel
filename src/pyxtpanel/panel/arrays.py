"""Array contract for the estimation packages (layer 2, ARCHITECTURE.md section 3).

This is the narrow, versioned surface that ``pyxtabond2``, ``pyxthreg`` and ``pyxtdpthresh``
build on.  :meth:`PanelBase.arrays` returns a frozen :class:`PanelArrays`: plain NumPy arrays in
a fixed layout, with no pandas object and no compiled kernel (Numba stays in the estimators).

Guarantees (a breaking change raises :data:`CONTRACT_VERSION`; adding a field does not):

* ``values`` is a ``float64`` C-contiguous ``(n_obs, k)`` array that the caller owns (a copy);
  missing values and infinities are ``NaN``.
* Rows are sorted by unit then calendar period, with no repeated (unit, period) pair.
* ``offsets``, ``unit_of_row``, ``time_index`` and ``source_row`` are read-only ``int64`` arrays
  shared with the panel (not copied): unit ``i`` owns the rows ``offsets[i]:offsets[i + 1]``.
* ``time_index`` is the position on the calendar, ``(t − t0) / delta``: consecutive periods
  differ by exactly 1 whatever ``delta`` or the date frequency.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, PyxtpanelWarning
from pyxtpanel._types import BoolArray, FloatArray, IntArray
from pyxtpanel.panel.declare import _column
from pyxtpanel.panel.index import PanelIndex, shift_positions

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

#: Version of the array contract; an estimator checks it at import time.
CONTRACT_VERSION = 1

_i18n.register(
    "panel",
    {
        "panel.not_numeric": {
            "en": "variable {name!r} is not numeric (type {dtype}); operators and transformations "
            "need numbers.",
            "fr": "la variable {name!r} n'est pas numérique (type {dtype}) ; les opérateurs et "
            "les transformations exigent des nombres.",
        },
        "panel.infinite_to_missing": {
            "en": "{count} infinite value(s) in {names} were read as missing.",
            "fr": "{count} valeur(s) infinie(s) dans {names} ont été lues comme manquantes.",
        },
        "panel.not_strongly_balanced": {
            "en": "as_cube needs a strongly balanced panel (every unit observed in every period); "
            "this one is {kind}. Use balance() or tsfill() first.",
            "fr": "as_cube exige un panel cylindré (chaque unité observée à chaque période) ; "
            "celui-ci est {kind}. Utilisez d'abord PanelData.balance() ou tsfill().",
        },
        "panel.bad_reindex": {
            "en": "reindex_back expected {expected} rows (one per observation), got {got}.",
            "fr": "reindex_back attendait {expected} lignes (une par observation), reçu {got}.",
        },
    },
)


def sorted_columns(panel: PanelBase, names: Sequence[str]) -> FloatArray:
    """Return the columns ``names`` as a ``(n_obs, len(names))`` C-contiguous ``float64`` array.

    Rows follow the sorted index of the panel (unit, then period).  Missing values (``None``,
    ``pd.NA``, ``NaN``) and infinities become ``NaN``; infinities raise a
    :class:`~pyxtpanel.PyxtpanelWarning`.  The array is a copy: changing it never touches the
    panel.

    Raises
    ------
    DataError
        If a variable does not exist or is not numeric (booleans count as numbers).
    """
    out = np.empty((panel.n_obs, len(names)), dtype=np.float64, order="C")
    rows = panel.index.source_row
    for j, name in enumerate(names):
        series = _column(panel.frame, name)
        try:
            out[:, j] = series.to_numpy(dtype=np.float64, na_value=np.nan)[rows]
        except (TypeError, ValueError):
            raise DataError(
                _i18n.message("panel.not_numeric", name=name, dtype=str(series.dtype)),
                key="panel.not_numeric",
            ) from None
    infinite = np.isinf(out)
    if infinite.any():
        hit = [names[int(j)] for j in np.flatnonzero(infinite.any(axis=0))]
        warnings.warn(
            _i18n.message(
                "panel.infinite_to_missing", count=int(infinite.sum()), names=", ".join(hit)
            ),
            PyxtpanelWarning,
            stacklevel=3,
        )
        out[infinite] = np.nan
    return out


@dataclass(frozen=True, eq=False, kw_only=True)
class PanelArrays:
    """Frozen arrays of a panel in the layout of the contract (see the module documentation).

    Attributes
    ----------
    values
        ``(n_obs, k)`` ``float64`` C-contiguous values of the requested columns, ``NaN`` where
        missing.
    columns
        Names of the columns of ``values``.
    offsets
        ``(n_units + 1,)`` ``int64``: unit ``i`` owns the rows ``offsets[i]:offsets[i + 1]``.
    unit_labels
        Original identifier of each unit, in sorted order.
    unit_of_row
        ``(n_obs,)`` ``int64`` unit code (``0..n_units - 1``) of each row.
    time_index
        ``(n_obs,)`` ``int64`` calendar position of each row, increasing within a unit.
    source_row
        ``(n_obs,)`` ``int64`` position of each row in the frame the panel was declared from.
    valid
        ``(n_obs,)`` boolean: no requested column is missing on the row.
    n_source
        Number of rows of the declared frame (length of what :meth:`reindex_back` returns).
    delta
        Calendar units between two consecutive periods of the original time variable.
    contract_version
        :data:`CONTRACT_VERSION` this object follows.
    """

    values: FloatArray
    columns: tuple[str, ...]
    offsets: IntArray
    unit_labels: np.ndarray
    unit_of_row: IntArray
    time_index: IntArray
    source_row: IntArray
    valid: BoolArray
    n_source: int
    delta: int = 1
    contract_version: int = CONTRACT_VERSION

    @property
    def n_obs(self) -> int:
        """Number of rows."""
        return int(self.time_index.size)

    @property
    def n_units(self) -> int:
        """Number of units."""
        return int(self.offsets.size - 1)

    @property
    def counts(self) -> IntArray:
        """Rows of each unit (``T_i``)."""
        return np.diff(self.offsets)

    @property
    def n_time(self) -> int:
        """Length of the calendar span: ``max(time_index) + 1`` (0 without rows)."""
        return int(np.max(self.time_index)) + 1 if self.n_obs else 0

    def unit_slices(self) -> list[slice]:
        """Return one ``slice`` of rows per unit, in unit order."""
        starts, stops = self.offsets[:-1].tolist(), self.offsets[1:].tolist()
        return [slice(a, b) for a, b in zip(starts, stops, strict=True)]

    def lag_index(self, k: int = 1) -> IntArray:
        """For each row, the row of the same unit ``k`` periods earlier, or ``-1`` if none.

        Computed on the calendar (``time_index``), never from the previous row; ``k`` may be
        negative for a lead.
        """
        return shift_positions(
            PanelIndex(
                unit_labels=self.unit_labels,
                offsets=self.offsets,
                unit=self.unit_of_row,
                rank=self.time_index,
                source_row=self.source_row,
            ),
            int(k),
        )

    def reindex_back(self, array: np.ndarray, fill: float = np.nan) -> np.ndarray:
        """Put a result computed in sorted order back in the order of the declared frame.

        Parameters
        ----------
        array
            One entry per row of the arrays (``(n_obs,)`` or ``(n_obs, m)``).
        fill
            Value for the rows of the declared frame that are not in the panel.

        Returns
        -------
        numpy.ndarray
            ``(n_source,)`` or ``(n_source, m)`` array (``float64``).
        """
        data = np.asarray(array, dtype=np.float64)
        if data.shape[0] != self.n_obs:
            raise DataError(
                _i18n.message("panel.bad_reindex", expected=self.n_obs, got=data.shape[0]),
                key="panel.bad_reindex",
            )
        out = np.full((self.n_source, *data.shape[1:]), fill, dtype=np.float64)
        out[self.source_row] = data
        return out

    def as_cube(self) -> FloatArray:
        """Return ``values`` as a ``(n_units, n_time, k)`` view (no copy).

        Only a strongly balanced panel has this shape: every unit is observed once in each of the
        same periods, so the sorted rows are the cells of a rectangle.

        Raises
        ------
        DataError
            If the panel is not strongly balanced.
        """
        periods = np.unique(self.time_index)
        if not np.array_equal(self.time_index, np.tile(periods, self.n_units)):
            raise DataError(
                _i18n.message("panel.not_strongly_balanced", kind="not strongly balanced"),
                key="panel.not_strongly_balanced",
            )
        cube: FloatArray = self.values.reshape(self.n_units, periods.size, len(self.columns))
        return cube


def build_arrays(panel: PanelBase, columns: Sequence[str]) -> PanelArrays:
    """Build the :class:`PanelArrays` of ``columns`` (see :meth:`PanelBase.arrays`)."""
    names = tuple(columns)
    values = sorted_columns(panel, names)
    index = panel.index
    return PanelArrays(
        values=values,
        columns=names,
        offsets=index.offsets,
        unit_labels=index.unit_labels,
        unit_of_row=index.unit,
        time_index=index.rank,
        source_row=index.source_row,
        valid=~np.isnan(values).any(axis=1),
        n_source=len(panel.frame),
        delta=panel.axis.delta,
    )
