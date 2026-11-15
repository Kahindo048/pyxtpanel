"""Conversions between long, wide and array forms, and guided handling of duplicates (K.5).

Additions of the package built on the panel structure; ``to_wide`` and ``to_long`` follow
Stata's ``reshape wide`` and ``reshape long`` (column names ``stub`` + time, one row per unit in
wide form, all the combinations of unit and period in long form).  Validated by hand on small
panels, by round trips, by invariance tests and, for the wide form, against Stata.

* :func:`to_wide`: one row per unit, one column per variable and observed period;
* :func:`to_long`: the inverse, from a wide DataFrame (``reshape long``), returning a plain
  DataFrame to declare with :func:`pyxtpanel.xtset`;
* :func:`to_array3d`: an array ``variable x period x unit`` with the mask of observed cells;
* :func:`resolve_duplicates`: merge the repeated (unit, period) pairs that ``xtset`` refuses and
  report what was merged.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import BoolArray, FloatArray, Lang
from pyxtpanel.explore._spec import ToolSpec, tool
from pyxtpanel.explore.variation import numeric_names
from pyxtpanel.panel.arrays import sorted_columns
from pyxtpanel.panel.declare import _levels

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "explore",
    {
        "explore.reshape.no_columns": {
            "en": "no column of the wide frame starts with {stubs} followed by {sep!r} and a "
            "period; check stubs= and sep=.",
            "fr": "aucune colonne du tableau large ne commence par {stubs} suivi de {sep!r} et "
            "d'une période ; vérifiez stubs= et sep=.",
        },
        "explore.reshape.suffix": {
            "en": "the suffix {suffix!r} of column {column!r} is not an integer period; pass "
            "times={{suffix: period}}.",
            "fr": "le suffixe {suffix!r} de la colonne {column!r} n'est pas une période entière ; "
            "passez times={{suffixe: période}}.",
        },
        "explore.reshape.id": {
            "en": "the identifier {name!r} is not a column of the wide frame.",
            "fr": "l'identifiant {name!r} n'est pas une colonne du tableau large.",
        },
        "explore.dup.how": {
            "en": "how must be one of {allowed}, got {value!r}.",
            "fr": "how doit valoir l'une de {allowed}, reçu {value!r}.",
        },
        "explore.dup.title": {
            "en": "Duplicates resolved with {how!r}: {groups} repeated pairs, {removed} rows "
            "removed ({n_in} rows in, {n_out} out)",
            "fr": "Doublons résolus avec {how!r} : {groups} paires répétées, {removed} lignes "
            "supprimées ({n_in} lignes en entrée, {n_out} en sortie)",
        },
        "explore.dup.none": {"en": "No repeated pair.", "fr": "Aucune paire répétée."},
    },
)


# ---------------------------------------------------------------- wide


@tool(
    "to_wide",
    ToolSpec(
        summary="One row per unit, one column per variable and observed period.",
        validation="hand-computed panels, round trip with to_long, comparison with reshape wide",
        stata="reshape wide",
        similar="panelr::widen_panel, plm reshape",
    ),
)
def to_wide(
    panel: PanelBase, variables: Sequence[str] | None = None, *, sep: str = ""
) -> pd.DataFrame:
    """Reshape to one row per unit (Stata's ``reshape wide``).

    Parameters
    ----------
    panel
        A declared panel.
    variables
        Variables to spread; by default every numeric variable other than the keys.
    sep
        Text between the variable name and the period in the column names (Stata uses none:
        ``x1990``).

    Returns
    -------
    pandas.DataFrame
        Indexed by unit.  Columns are ordered by variable then period and named
        ``variable + sep + period``; a period a unit was not observed in is missing.  Only the
        periods observed by at least one unit get a column, as in Stata.
    """
    names = list(variables) if variables is not None else numeric_names(panel)
    index = panel.index
    periods = np.unique(index.rank)
    position = np.searchsorted(periods, index.rank)
    labels = panel.axis.labels(periods)
    data = sorted_columns(panel, names)
    blocks, columns = [], []
    for j, name in enumerate(names):
        block = np.full((index.n_units, periods.size), np.nan)
        block[index.unit, position] = data[:, j]
        blocks.append(block)
        columns += [f"{name}{sep}{label}" for label in labels]
    return pd.DataFrame(
        np.hstack(blocks),
        index=pd.Index(index.unit_labels, name=panel.id_name),
        columns=columns,
    )


@tool(
    "to_long",
    ToolSpec(
        summary="Back from the wide form: one row per unit and period, ready for xtset.",
        validation="hand-computed frames, round trip with to_wide, invariance tests",
        stata="reshape long",
        similar="panelr::long_panel",
        on_panel=False,
    ),
)
def to_long(
    wide: pd.DataFrame,
    stubs: Sequence[str],
    *,
    id: str | None = None,  # noqa: A002 - the name Stata users expect
    time: str = "t",
    sep: str = "",
    times: Mapping[str, Any] | None = None,
    drop_empty: bool = False,
) -> pd.DataFrame:
    """Reshape a wide frame to one row per unit and period (Stata's ``reshape long``).

    Parameters
    ----------
    wide
        One row per unit; columns ``stub + sep + period`` for each stub.
    stubs
        Names of the variables spread over columns.
    id
        Column that identifies the unit; the index of ``wide`` when omitted.
    time
        Name of the new time column.
    sep
        Text between the stub and the period in the column names.
    times
        Maps a column suffix to the time value it stands for, for suffixes that are not integers
        (for example ``{"q1": pd.Timestamp("2020-01-01")}``).
    drop_empty
        Drop the rows in which every stub is missing (Stata keeps them).

    Returns
    -------
    pandas.DataFrame
        Sorted by unit then period, with the other columns of ``wide`` repeated on every row.
        Declare it with :func:`pyxtpanel.xtset`.
    """
    stubs = list(stubs)
    if id is not None and id not in wide.columns:
        raise SpecificationError(
            _i18n.message("explore.reshape.id", name=id), stata_code=198, key="explore.reshape.id"
        )
    id_name = id if id is not None else str(wide.index.name or "id")
    found: dict[str, dict[str, str]] = {}  # suffix -> {stub: column}
    for column in wide.columns:
        for stub in stubs:
            prefix = f"{stub}{sep}"
            name = str(column)
            if name.startswith(prefix) and len(name) > len(prefix):
                suffix = name[len(prefix) :]
                if times is None and not re.fullmatch(r"-?\d+", suffix):
                    continue
                if times is not None and suffix not in times:
                    continue
                found.setdefault(suffix, {})[stub] = column
    if not found:
        raise SpecificationError(
            _i18n.message("explore.reshape.no_columns", stubs=", ".join(stubs), sep=sep),
            stata_code=198,
            key="explore.reshape.no_columns",
        )
    used = {column for columns in found.values() for column in columns.values()}
    carried = [c for c in wide.columns if c not in used and c != id]
    ids = wide[id] if id is not None else pd.Series(wide.index, index=wide.index)
    parts = []
    for suffix in sorted(found, key=lambda s: int(s) if times is None else list(times).index(s)):
        value = int(suffix) if times is None else times[suffix]
        part = pd.DataFrame({id_name: ids.to_numpy()})
        part[time] = value
        for stub in stubs:
            source = found[suffix].get(stub)
            part[stub] = wide[source].to_numpy() if source is not None else np.nan
        for name in carried:
            part[name] = wide[name].to_numpy()
        parts.append(part)
    long = pd.concat(parts, ignore_index=True)
    long = long.sort_values([id_name, time], kind="stable").reset_index(drop=True)
    if drop_empty:
        long = long[long[stubs].notna().any(axis=1)].reset_index(drop=True)
    return long


# ---------------------------------------------------------------- array


@dataclass(frozen=True, kw_only=True)
class Array3D:
    """A panel as an array ``variable x period x unit`` (result of :func:`to_array3d`).

    Attributes
    ----------
    values
        Array of shape ``(n_variables, n_periods, n_units)``; missing where the unit is not
        observed in the period or the value is missing.
    variables
        Names along the first axis.
    periods
        Time values along the second axis: every period of the calendar from the first to the
        last, in steps of ``delta``.
    units
        Unit identifiers along the third axis, in sorted order.
    observed
        Boolean ``(n_periods, n_units)``: the unit has a row in that period (its value may
        still be missing).
    """

    values: FloatArray
    variables: tuple[str, ...]
    periods: pd.Index
    units: np.ndarray
    observed: BoolArray


@tool(
    "to_array3d",
    ToolSpec(
        summary="Array variable x period x unit over the whole calendar, with the mask of rows.",
        validation="hand-computed panels, agreement with the array contract of a balanced panel",
        similar="plm values3d, Stata psmat (one variable)",
    ),
)
def to_array3d(panel: PanelBase, variables: Sequence[str] | None = None) -> Array3D:
    """Lay the panel out as an array over the full calendar.

    Unlike :meth:`PanelData.arrays`, which keeps only the rows that exist, this array has a
    cell for every (period, unit), so that lags are positions along the second axis.
    """
    names = list(variables) if variables is not None else numeric_names(panel)
    index = panel.index
    low, high = int(np.min(index.rank)), int(np.max(index.rank))
    data = sorted_columns(panel, names)
    cube = np.full((len(names), high - low + 1, index.n_units), np.nan)
    cube[:, index.rank - low, index.unit] = data.T
    observed = np.zeros((high - low + 1, index.n_units), dtype=np.bool_)
    observed[index.rank - low, index.unit] = True
    return Array3D(
        values=cube,
        variables=tuple(names),
        periods=panel.axis.values(np.arange(low, high + 1)),
        units=index.unit_labels,
        observed=observed,
    )


# ---------------------------------------------------------------- duplicates

#: Ways of merging the rows that share a (unit, period) pair.
HOW = ("first", "last", "mean", "sum", "min", "max", "median")


@dataclass(frozen=True, kw_only=True)
class DuplicateResolution(Result):
    """Result of :func:`resolve_duplicates`.

    Attributes
    ----------
    frame
        The data with one row per (unit, period) pair, in the original order of the rows kept.
    groups
        One row per repeated pair: the keys, ``n_rows`` and ``differing`` (the columns whose
        values are not identical across the repeated rows, separated by commas).
    how
        The rule applied.
    n_rows_in, n_rows_out
        Rows before and after.
    """

    frame: pd.DataFrame
    groups: pd.DataFrame
    how: str
    n_rows_in: int
    n_rows_out: int

    def to_frame(self) -> pd.DataFrame:
        """Return the table of repeated pairs."""
        return self.groups

    def summary(self, lang: Lang | None = None) -> str:
        """Return the totals and the repeated pairs as text."""
        head = _i18n.message(
            "explore.dup.title",
            lang,
            how=self.how,
            groups=len(self.groups),
            removed=self.n_rows_in - self.n_rows_out,
            n_in=self.n_rows_in,
            n_out=self.n_rows_out,
        )
        if not len(self.groups):
            return _i18n.message("explore.dup.none", lang)
        return head + "\n\n" + self.groups.head(30).to_string(index=False)


@tool(
    "resolve_duplicates",
    ToolSpec(
        summary="Merge repeated (unit, period) pairs by first, last, mean, sum... and report "
        "what was merged.",
        validation="hand-computed frames, independent loop implementation, invariance tests",
        similar="Stata duplicates drop (identical rows only)",
        on_panel=False,
    ),
)
def resolve_duplicates(
    frame: pd.DataFrame,
    id: str | None = None,  # noqa: A002 - the name Stata users expect
    time: str | None = None,
    *,
    how: Literal["first", "last", "mean", "sum", "min", "max", "median"] = "first",
) -> DuplicateResolution:
    """Merge the rows that share a (unit, period) pair so that ``xtset`` accepts the data.

    ``"first"`` and ``"last"`` keep one of the repeated rows whole.  The other rules replace the
    numeric columns (not the keys, not booleans) by their statistic over the repeated rows,
    ignoring missing values (``sum`` is missing when all are), and keep the first row for every
    other column; a numeric column of integers becomes float when it is averaged.  Rows with a
    missing key are left untouched.

    Parameters
    ----------
    frame
        The raw data.
    id, time
        Names of the key columns (or ``None`` when a two-level index holds them).
    how
        The merging rule.
    """
    if how not in HOW:
        raise SpecificationError(
            _i18n.message("explore.dup.how", allowed=", ".join(HOW), value=how),
            stata_code=198,
            key="explore.dup.how",
        )
    id_series, time_series = _levels(frame, id, time)
    keys = pd.DataFrame({"i": id_series.to_numpy(), "t": time_series.to_numpy()})
    usable = keys.notna().all(axis=1).to_numpy()
    repeated = keys.duplicated(keep=False).to_numpy() & usable
    n_in = len(frame)
    if not repeated.any():
        empty = pd.DataFrame(columns=[id or "id", time or "t", "n_rows", "differing"])
        return DuplicateResolution(
            frame=frame.copy(), groups=empty, how=how, n_rows_in=n_in, n_rows_out=n_in
        )
    positions = np.flatnonzero(repeated)  # rows of ``frame`` that belong to a repeated pair
    group = keys.iloc[positions].groupby(["i", "t"], sort=False).ngroup().to_numpy()
    sub = frame.iloc[positions]
    _, first_of = np.unique(group, return_index=True)  # first row of each pair, within sub
    id_name, time_name = id or "id", time or "t"
    distinct = (
        sub.assign(_g=group).groupby("_g").nunique(dropna=False).drop(columns="_g", errors="ignore")
    )
    groups = pd.DataFrame(
        {
            id_name: keys["i"].to_numpy()[positions][first_of],
            time_name: keys["t"].to_numpy()[positions][first_of],
            "n_rows": np.bincount(group),
            "differing": [
                ", ".join(str(c) for c in distinct.columns[row > 1]) for row in distinct.to_numpy()
            ],
        }
    )
    keep = ~keys.duplicated(keep="last" if how == "last" else "first").to_numpy() | ~usable
    out = frame.iloc[np.flatnonzero(keep)].copy()
    if how not in ("first", "last"):
        skip = {c for c in (id, time) if c is not None}
        numeric = [
            c
            for c in frame.columns
            if c not in skip
            and pd.api.types.is_numeric_dtype(frame[c])
            and not pd.api.types.is_bool_dtype(frame[c])
        ]
        merged = sub.assign(_g=group).groupby("_g", sort=True)[numeric]
        stat = merged.sum(min_count=1) if how == "sum" else getattr(merged, how)()
        at = np.searchsorted(np.flatnonzero(keep), positions[first_of])
        for column in numeric:
            values = out[column].to_numpy(dtype=stat[column].dtype, copy=True)
            values[at] = stat[column].to_numpy()
            out[column] = values
    return DuplicateResolution(
        frame=out, groups=groups, how=how, n_rows_in=n_in, n_rows_out=len(out)
    )
