"""Sample selection: Stata's ``if`` and ``markout``, restriction to units and periods.

A selection never copies the data: the sub-panel keeps the frame, the declaration and the
calendar, and only its sorted index is restricted.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, DeclarationError
from pyxtpanel._types import BoolArray
from pyxtpanel.panel.checks import unknown_variable

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "select",
    {
        "select.where_not_boolean": {
            "en": "where= must give one True/False per row ({n} rows), got {got}.",
            "fr": "where= doit donner un Vrai/Faux par ligne ({n} lignes), reçu {got}.",
        },
    },
)


def resolve_where(frame: pd.DataFrame, where: str | Sequence[bool] | BoolArray | None) -> BoolArray:
    """Turn ``where`` into a boolean mask over the rows of ``frame``.

    Parameters
    ----------
    frame
        The data.
    where
        ``None`` keeps every row; a boolean array or Series (aligned by position); or a string
        evaluated with :meth:`pandas.DataFrame.eval`, the counterpart of Stata's ``if``.

    Raises
    ------
    DataError
        If ``where`` is not boolean or has the wrong length.
    """
    # STATA-QUIRK (inverted): a missing value in an ``if`` condition counts as TRUE in Stata;
    # here it counts as False, because keeping a row on an unknown condition is the riskier error.
    n = len(frame)
    if where is None:
        return np.ones(n, dtype=np.bool_)
    result = frame.eval(where) if isinstance(where, str) else where
    mask = np.asarray(result)
    if mask.dtype == object:
        mask = pd.array(mask, dtype="boolean").fillna(False).to_numpy(dtype=bool)
    if mask.dtype.kind != "b" or mask.shape != (n,):
        raise DataError(
            _i18n.message("select.where_not_boolean", n=n, got=f"{mask.dtype} {mask.shape}"),
            key="select.where_not_boolean",
        )
    return np.asarray(mask, dtype=np.bool_)


def complete_rows(frame: pd.DataFrame, columns: Sequence[str]) -> BoolArray:
    """Rows of ``frame`` where every listed column is present (Stata's ``markout``).

    Numeric columns must also be finite (an infinite value is not a number in Stata).
    """
    keep = np.ones(len(frame), dtype=np.bool_)
    for name in columns:
        if name not in frame.columns:
            raise unknown_variable(name, frame.columns)
        column = frame[name]
        if pd.api.types.is_numeric_dtype(column) and not pd.api.types.is_bool_dtype(column):
            values = column.to_numpy(dtype=np.float64, na_value=np.nan)
            keep &= np.isfinite(values)
        else:
            keep &= column.notna().to_numpy()
    return keep


def selection_mask(
    panel: PanelBase,
    *,
    where: str | Sequence[bool] | BoolArray | None = None,
    units: Sequence[Any] | None = None,
    periods: Sequence[Any] | slice | None = None,
    complete: Sequence[str] | None = None,
) -> BoolArray:
    """Mask over the sorted rows of ``panel`` that satisfy every given condition."""
    index = panel.index
    keep = np.ones(index.n_obs, dtype=np.bool_)
    if where is not None:
        keep &= resolve_where(panel.frame, where)[index.source_row]
    if complete is not None:
        keep &= complete_rows(panel.frame, complete)[index.source_row]
    if units is not None:
        keep &= np.isin(index.unit_labels, list(units))[index.unit]
    if periods is not None:
        axis = panel.axis
        if isinstance(periods, slice):
            low = -np.inf if periods.start is None else float(axis.rank_of([periods.start])[0])
            high = np.inf if periods.stop is None else float(axis.rank_of([periods.stop])[0])
            keep &= (index.rank >= low) & (index.rank <= high)
        else:
            keep &= np.isin(index.rank, axis.rank_of(list(periods)))
    return keep


def check_not_empty(keep: BoolArray) -> None:
    """Raise if a selection leaves no observation."""
    if not keep.any():
        raise DeclarationError(
            _i18n.message("panel.no_observations"), stata_code=2000, key="panel.no_observations"
        )
