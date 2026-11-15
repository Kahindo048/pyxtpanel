"""Means and counts by unit or by period, broadcast on the rows (D.2.2).

The Stata equivalents quoted in the registry (``egen``) are idioms, not commands: they are
checked against pandas and by hand.

All functions take the sorted index and a ``float64`` array in sorted order, ``(n_obs,)`` or
``(n_obs, k)``; missing values are ``NaN`` and are ignored by the means and the counts.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from pyxtpanel._types import FloatArray, IntArray
from pyxtpanel.ops._spec import OpSpec, operator
from pyxtpanel.panel.index import PanelIndex

By = Literal["unit", "time"]


def as_columns(values: FloatArray) -> tuple[FloatArray, bool]:
    """``values`` as a 2-D float array and whether it was 1-D (to give the shape back)."""
    array = np.asarray(values, dtype=np.float64)
    return (array.reshape(-1, 1), True) if array.ndim == 1 else (array, False)


def group_codes(index: PanelIndex, by: By) -> tuple[IntArray, int]:
    """Group code of each sorted row and the number of groups (units, or distinct periods)."""
    if by == "unit":
        return index.unit, index.n_units
    present, codes = np.unique(index.rank, return_inverse=True)
    return codes.astype(np.int64), int(present.size)


def _mean_1d(codes: IntArray, size: int, x: FloatArray, at_missing: bool) -> FloatArray:
    valid = ~np.isnan(x)
    total = np.bincount(codes, weights=np.where(valid, x, 0.0), minlength=size)
    count = np.bincount(codes, weights=valid.astype(np.float64), minlength=size)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(count > 0, total / count, np.nan)
    out: FloatArray = mean[codes]
    if not at_missing:
        out[~valid] = np.nan
    return out


@operator(
    "group_mean",
    OpSpec(
        kind="transform",
        missing="missing where the row's own value is missing (unless at_missing), or where "
        "the group has no observed value",
        stata="bysort id: egen m = mean(x)",
        summary="Mean of the group (unit or period) repeated on each row.",
    ),
)
def group_mean(
    index: PanelIndex, values: FloatArray, by: By = "unit", *, at_missing: bool = False
) -> FloatArray:
    """Mean of each unit (or period) over its observed values, repeated on its rows.

    With ``at_missing=True`` the mean is also given on the rows whose own value is missing.
    """
    codes, size = group_codes(index, by)
    data, flat = as_columns(values)
    out = np.column_stack(
        [_mean_1d(codes, size, data[:, j], at_missing) for j in range(data.shape[1])]
    )
    result: FloatArray = out[:, 0] if flat else out
    return result


@operator(
    "group_count",
    OpSpec(
        kind="transform",
        missing="missing where the row's own value is missing (unless at_missing)",
        stata="bysort id: egen n = count(x)",
        summary="Number of observed values of the group (unit or period) repeated on each row.",
    ),
)
def group_count(
    index: PanelIndex, values: FloatArray, by: By = "unit", *, at_missing: bool = False
) -> FloatArray:
    """Count the observed values of each unit (or period) and repeat it on the rows."""
    codes, size = group_codes(index, by)
    data, flat = as_columns(values)
    columns = []
    for j in range(data.shape[1]):
        valid = ~np.isnan(data[:, j])
        count = np.bincount(codes, weights=valid.astype(np.float64), minlength=size)[codes]
        columns.append(count if at_missing else np.where(valid, count, np.nan))
    out = np.column_stack(columns)
    result: FloatArray = out[:, 0] if flat else out
    return result


def group_size(index: PanelIndex, by: By = "unit") -> FloatArray:
    """Return the rows of the group (``T_i`` or ``N_t``) repeated on each of its rows."""
    codes, size = group_codes(index, by)
    counts = np.bincount(codes, minlength=size).astype(np.float64)
    result: FloatArray = counts[codes]
    return result
