"""Declaration of a panel (Stata's ``xtset``): the one place where a DataFrame is validated.

``declare`` reads the identifier and the time variable, encodes the time on a regular calendar,
checks the pairs, and builds the sorted index.  It never modifies the user's frame.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import DeclarationError, SpecificationError
from pyxtpanel._types import IntArray
from pyxtpanel.panel.checks import MAX_ROWS_LISTED, check_unique_pairs, unknown_variable
from pyxtpanel.panel.index import PanelIndex, build_index
from pyxtpanel.time.axis import TimeAxis, check_span, encode_time

MissingKeys = Literal["error", "drop"]


@dataclass(frozen=True)
class Declaration:
    """What ``declare`` found: names, calendar, sorted index and rows set aside."""

    id_name: str
    time_name: str
    axis: TimeAxis
    index: PanelIndex
    n_dropped: int


def _column(frame: pd.DataFrame, name: str) -> pd.Series:
    """Return a key as a positional Series, from a column or an index level."""
    if name in frame.columns:
        found = frame[name]
        if isinstance(found, pd.DataFrame):
            raise DeclarationError(
                _i18n.message("panel.duplicate_column", name=name), key="panel.duplicate_column"
            )
        flat: pd.Series = found.reset_index(drop=True)
        return flat
    if name in frame.index.names:
        level: pd.Series = pd.Series(frame.index.get_level_values(name))
        return level
    raise unknown_variable(name, [*frame.columns, *(n for n in frame.index.names if n is not None)])


def _keys(frame: pd.DataFrame, id: str | None, time: str | None) -> tuple[str, str]:  # noqa: A002
    """Names of the identifier and the time variable (the two index levels by default)."""
    if (id is None) != (time is None):
        raise SpecificationError(
            _i18n.message("panel.one_key_given"), stata_code=198, key="panel.one_key_given"
        )
    if id is not None and time is not None:
        return id, time
    if frame.index.nlevels != 2:
        raise SpecificationError(
            _i18n.message("panel.keys_not_found"), stata_code=198, key="panel.keys_not_found"
        )
    first, second = frame.index.names
    return str(first or "id"), str(second or "time")


def _levels(frame: pd.DataFrame, id: str | None, time: str | None) -> tuple[pd.Series, pd.Series]:  # noqa: A002
    if id is None:  # two-level index, possibly unnamed
        return (
            pd.Series(frame.index.get_level_values(0)),
            pd.Series(frame.index.get_level_values(1)),
        )
    assert time is not None
    return _column(frame, id), _column(frame, time)


def _refuse_missing(name: str, what: str, missing: np.ndarray) -> None:
    rows = np.flatnonzero(missing)
    raise DeclarationError(
        _i18n.message(
            "panel.key_missing",
            what=what,
            name=name,
            n=int(missing.sum()),
            rows=rows[:MAX_ROWS_LISTED].tolist(),
        ),
        stata_code=459,
        key="panel.key_missing",
    )


def _factorize(ids: pd.Series, name: str) -> tuple[IntArray, np.ndarray]:
    try:
        codes, labels = pd.factorize(ids, sort=True)
    except TypeError as error:
        raise DeclarationError(
            _i18n.message("panel.id_unsortable", name=name, error=error), key="panel.id_unsortable"
        ) from error
    return codes.astype(np.int64), np.asarray(labels)


def declare(
    frame: pd.DataFrame,
    id: str | None = None,  # noqa: A002 - the name Stata users expect
    time: str | None = None,
    *,
    delta: int | None = None,
    freq: str | None = None,
    format: str | None = None,  # noqa: A002
    missing_keys: MissingKeys = "error",
) -> Declaration:
    """Validate ``frame`` as a panel and build its sorted index.

    Parameters
    ----------
    frame
        Long data, one row per (unit, period).
    id, time
        Names of the identifier and the time variable, columns or index levels.  When both are
        omitted, the two levels of a MultiIndex are used, in that order.
    delta
        Calendar units between two consecutive periods (default 1).
    freq
        Calendar unit of dates, see :func:`pyxtpanel.time.axis.encode_time`.
    format
        Stata format of Stata-coded integer times, see the same function.
    missing_keys
        ``"error"`` (default) refuses a missing identifier or time and names the rows;
        ``"drop"`` ignores those observations, as Stata does, and counts them.

    Raises
    ------
    DeclarationError
        On missing keys (459), repeated (unit, period) pairs (451), duplicated key columns,
        unsortable identifiers, or no usable observation (2000).
    TimeAxisError
        If the time variable is not numeric, a date or a period, or is off the calendar grid.
    SpecificationError
        On invalid options, an unknown variable, or only one of ``id`` and ``time`` given.
    """
    if missing_keys not in ("error", "drop"):
        raise SpecificationError(
            _i18n.message("panel.missing_keys_option", value=missing_keys),
            stata_code=198,
            key="panel.missing_keys_option",
        )
    id_name, time_name = _keys(frame, id, time)
    id_series, time_series = _levels(frame, id, time)

    missing_id = id_series.isna().to_numpy()
    missing_time = time_series.isna().to_numpy()
    if missing_keys == "error":
        for name, what, flags in ((id_name, "id", missing_id), (time_name, "time", missing_time)):
            if flags.any():
                _refuse_missing(name, what, flags)
    keep = ~(missing_id | missing_time)
    positions = np.flatnonzero(keep).astype(np.int64)
    if positions.size == 0:
        raise DeclarationError(
            _i18n.message("panel.no_observations"), stata_code=2000, key="panel.no_observations"
        )

    kept_id = id_series[keep].reset_index(drop=True)
    kept_time = time_series[keep].reset_index(drop=True)
    unit_code, labels = _factorize(kept_id, id_name)
    encoded = encode_time(
        kept_time,
        name=time_name,
        delta=delta,
        freq=freq,
        format=format,
        units=kept_id.to_numpy(),
        unit_name=id_name,
        positions=positions,
    )
    check_span(encoded.axis, labels.size)
    check_unique_pairs(
        unit_code, encoded.rank, id_name, time_name, (kept_id.to_numpy(), kept_time.to_numpy())
    )
    index = build_index(unit_code, encoded.rank, positions, labels)
    return Declaration(
        id_name=id_name,
        time_name=time_name,
        axis=encoded.axis,
        index=index,
        n_dropped=int((~keep).sum()),
    )
