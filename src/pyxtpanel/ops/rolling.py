"""Rolling and cumulative statistics per unit, and growth rates, on the calendar (D.2.7).

The Stata equivalents quoted in the registry (``tsegen``, ``sum()``, ``D.x / L.x``) are idioms,
not commands: these functions are checked by hand and against pandas.

A window of ``w`` periods ends at the row's period and covers the ``w`` calendar periods
``t−w+1 … t``: a period the unit skipped is missing from the window instead of being replaced by
an older row, exactly like a lag.
"""

from __future__ import annotations

import warnings
from typing import Literal

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._types import FloatArray
from pyxtpanel.ops._spec import OpSpec, operator
from pyxtpanel.ops.groups import as_columns
from pyxtpanel.ops.lags import lag
from pyxtpanel.panel.index import PanelIndex

Stat = Literal["mean", "sum", "std", "min", "max"]

_i18n.register(
    "ops",
    {
        "ops.bad_window": {
            "en": "window must be a positive integer and min_obs between 1 and the window, "
            "got window={window}, min_obs={min_obs}.",
            "fr": "window doit être un entier positif et min_obs compris entre 1 et la fenêtre, "
            "reçu window={window}, min_obs={min_obs}.",
        },
        "ops.unknown_stat": {
            "en": "unknown statistic {stat!r}; available: mean, sum, std, min, max.",
            "fr": "statistique inconnue {stat!r} ; disponibles : mean, sum, std, min, max.",
        },
    },
)


@operator(
    "rolling",
    OpSpec(
        kind="transform",
        missing="fewer than min_obs observed values among the window's calendar periods",
        stata="rangestat / tsegen (rowmean(L(0/w-1).x))",
        summary="Mean, sum, standard deviation, min or max over the last w calendar periods.",
    ),
)
def rolling(
    index: PanelIndex,
    values: FloatArray,
    window: int,
    stat: Stat = "mean",
    *,
    min_obs: int | None = None,
) -> FloatArray:
    """Statistic of the ``window`` calendar periods ending at each row (the row included).

    Parameters
    ----------
    window
        Number of calendar periods in the window.
    stat
        ``"mean"``, ``"sum"``, ``"std"`` (``n − 1`` divisor), ``"min"`` or ``"max"``.
    min_obs
        Observed values required in the window (default: all of them, ``window``).
    """
    needed = window if min_obs is None else min_obs
    if int(window) < 1 or not 1 <= int(needed) <= int(window):
        raise SpecificationError(
            _i18n.message("ops.bad_window", window=window, min_obs=min_obs), key="ops.bad_window"
        )
    if stat not in ("mean", "sum", "std", "min", "max"):
        raise SpecificationError(
            _i18n.message("ops.unknown_stat", stat=stat), key="ops.unknown_stat"
        )
    data, flat = as_columns(values)
    columns = []
    for j in range(data.shape[1]):
        stack = np.vstack([lag(index, data[:, j], k) for k in range(int(window))])
        seen = np.sum(~np.isnan(stack), axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-missing windows give NaN
            if stat == "std":
                column = np.nanstd(stack, axis=0, ddof=1)
            elif stat == "mean":
                column = np.nanmean(stack, axis=0)
            elif stat == "sum":
                column = np.nansum(stack, axis=0)
            elif stat == "min":
                column = np.nanmin(stack, axis=0)
            else:
                column = np.nanmax(stack, axis=0)
        columns.append(np.where(seen >= int(needed), column, np.nan))
    out = np.column_stack(columns)
    result: FloatArray = out[:, 0] if flat else out
    return result


@operator(
    "cumulative",
    OpSpec(
        kind="transform",
        missing="missing where the row's own value is missing",
        stata="bysort id (t): gen c = sum(x)",
        summary="Running sum or mean of the unit's observed values up to each row.",
    ),
)
def cumulative(
    index: PanelIndex, values: FloatArray, stat: Literal["sum", "mean"] = "sum"
) -> FloatArray:
    """Return the running sum (or mean) of the unit's observed values up to each row."""
    data, flat = as_columns(values)
    out = np.full_like(data, np.nan)
    for j in range(data.shape[1]):
        valid = ~np.isnan(data[:, j])
        column = np.where(valid, data[:, j], 0.0)
        seen = valid.astype(np.float64)
        for i in range(index.n_units):  # per unit, so that no sum carries another unit's total
            rows = slice(int(index.offsets[i]), int(index.offsets[i + 1]))
            running = np.cumsum(column[rows])
            if stat == "mean":
                with np.errstate(invalid="ignore", divide="ignore"):
                    running = running / np.cumsum(seen[rows])
            out[rows, j] = running
        out[~valid, j] = np.nan
    result: FloatArray = out[:, 0] if flat else out
    return result


@operator(
    "growth",
    OpSpec(
        kind="difference",
        missing="either period is not observed, or the earlier value is zero (or not positive "
        "for log)",
        stata="D.x / L.x  or  D.ln(x)",
        summary="Growth rate over `periods` calendar periods (relative or log difference).",
    ),
)
def growth(
    index: PanelIndex, values: FloatArray, periods: int = 1, *, log: bool = False
) -> FloatArray:
    """``x_t / x_{t−p} − 1`` (or ``ln x_t − ln x_{t−p}`` with ``log=True``)."""
    data = np.asarray(values, dtype=np.float64)
    before = lag(index, data, int(periods))
    with np.errstate(invalid="ignore", divide="ignore"):
        if log:
            result: FloatArray = np.log(data) - np.log(before)
        else:
            result = data / before - 1.0
    result[~np.isfinite(result)] = np.nan
    return result
