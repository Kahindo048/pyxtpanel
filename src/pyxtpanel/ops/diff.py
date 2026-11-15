"""Differences on the calendar (Stata's ``D.`` and ``S.``)."""

from __future__ import annotations

import numpy as np

from pyxtpanel._types import FloatArray, IntArray
from pyxtpanel.ops._spec import OpSpec, operator
from pyxtpanel.ops.lags import lag, shift_positions
from pyxtpanel.panel.index import PanelIndex


@operator(
    "difference",
    OpSpec(
        kind="difference",
        missing="either of the two periods is not observed for the unit, or its value is missing",
        stata="D.x, D2.x, ...",
        summary="Difference with the previous period, repeated `order` times.",
    ),
)
def difference(index: PanelIndex, values: FloatArray, order: int = 1) -> FloatArray:
    """``D{order}.x``: ``x_t − x_{t−1}`` repeated ``order`` times.

    The second difference is the difference of the first one.
    """
    out = values
    for _ in range(int(order)):
        out = out - lag(index, out, 1)
    return out


@operator(
    "seasonal",
    OpSpec(
        kind="difference",
        missing="either of the two periods is not observed for the unit, or its value is missing",
        stata="S12.x",
        summary="Difference with the value `period` periods earlier.",
    ),
)
def seasonal(index: PanelIndex, values: FloatArray, period: int = 1) -> FloatArray:
    """``S{period}.x``: ``x_t − x_{t−period}``."""
    return values - lag(index, values, int(period))


@operator(
    "fod",
    OpSpec(
        kind="difference",
        missing="the row's own value is missing, or no later observation of the unit is observed",
        stata="xtabond2, orthogonal  (forward orthogonal deviations)",
        summary="Deviation from the mean of the later observed values of the unit, scaled.",
    ),
)
def fod(index: PanelIndex, values: FloatArray, *, scale: bool = True) -> FloatArray:
    """Forward orthogonal deviations (Arellano and Bover).

    For a row with value ``x`` and ``F`` later observed values of the same unit, the result is
    ``√(F / (F + 1)) · (x − mean of the later values)``.  "Later" means the next observed rows of
    the unit, whatever the gaps between their dates (the definition of ``xtabond2``).  With
    ``scale=False`` the factor is omitted.  Missing where ``x`` is missing or ``F = 0``.
    """
    data = np.asarray(values, dtype=np.float64)
    flat = data.ndim == 1
    data = data.reshape(-1, 1) if flat else data
    ends = index.offsets[1:][index.unit] - 1  # last sorted row of the unit of each row
    out = np.full_like(data, np.nan)
    for j in range(data.shape[1]):
        x = data[:, j]
        valid = ~np.isnan(x)
        # translating a unit's values leaves its deviations unchanged: centre them first so that
        # the running sums stay small
        centre = np.bincount(index.unit, weights=np.where(valid, x, 0.0), minlength=index.n_units)
        seen = np.bincount(index.unit, weights=valid.astype(np.float64), minlength=index.n_units)
        with np.errstate(invalid="ignore", divide="ignore"):
            centre = np.where(seen > 0, centre / seen, 0.0)
        shifted = np.where(valid, x - centre[index.unit], 0.0)
        total, count = np.cumsum(shifted), np.cumsum(valid)
        later_sum = total[ends] - total
        later_n = (count[ends] - count).astype(np.float64)
        usable = valid & (later_n > 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            deviation = (shifted - later_sum / later_n) * (
                np.sqrt(later_n / (later_n + 1.0)) if scale else 1.0
            )
        out[usable, j] = deviation[usable]
    result: FloatArray = out[:, 0] if flat else out
    return result


def fd_equations(index: PanelIndex) -> tuple[IntArray, IntArray]:
    """Rows of the first-difference equations: ``(row, previous_row)`` for each ``D.x`` that exists.

    A first difference exists where the unit is observed in the period just before; both are
    sorted-row positions.  Estimators build their transformed equations from these pairs.
    """
    previous = shift_positions(index, 1)
    rows = np.flatnonzero(previous >= 0).astype(np.int64)
    return rows, previous[rows]
