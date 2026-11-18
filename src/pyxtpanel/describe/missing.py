"""Missing values by variable, by period and by unit, and complete cases (E.5).

Only the rows that belong to the panel are looked at (rows set aside for a missing key or left
out by a selection do not count).  A value is missing when pandas says so.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.describe._spec import DescribeSpec, description
from pyxtpanel.describe.structure import PERCENTILES, stata_quantile

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "describe",
    {
        "describe.ms.no_variable": {
            "en": "no variable to look at: give variables= or add columns besides the keys.",
            "fr": "aucune variable à examiner : donnez variables= ou ajoutez des colonnes "
            "en plus des clés.",
        },
    },
)


def _flags(panel: PanelBase, variables: Sequence[str] | None) -> tuple[list[str], pd.DataFrame]:
    """Names and the missing-value flags of the panel's rows, in sorted order."""
    frame = panel.frame.reset_index(drop=True)
    names = (
        list(variables)
        if variables is not None
        else [c for c in frame.columns if c not in (panel.id_name, panel.time_name)]
    )
    if not names:
        raise SpecificationError(
            _i18n.message("describe.ms.no_variable"), key="describe.ms.no_variable"
        )
    for name in names:
        if name not in frame.columns:
            from pyxtpanel.panel.checks import unknown_variable

            raise unknown_variable(name, frame.columns)
    flags = frame.loc[panel.index.source_row, names].isna()
    return names, flags.reset_index(drop=True)


@dataclass(frozen=True, kw_only=True)
class MissingReport(Result):
    """Result of :func:`missing_report`.

    Attributes
    ----------
    by_variable
        Per variable: observations, missing, percent missing, units entirely missing and units
        partly missing.
    by_period
        Per period (rows) and variable (columns): share of missing values; ``n_obs`` is the
        number of observations of the period.
    by_unit
        Per unit (rows) and variable (columns): number of missing values; ``T_i`` is the
        number of rows of the unit.
    patterns
        Joint patterns across the variables (``1`` observed, ``.`` missing, one symbol per
        variable), with the number and share of rows, most frequent first.
    """

    by_variable: pd.DataFrame
    by_period: pd.DataFrame
    by_unit: pd.DataFrame
    patterns: pd.DataFrame

    def to_frame(self) -> pd.DataFrame:
        """Return the table by variable."""
        return self.by_variable

    def summary(self, lang: Lang | None = None) -> str:
        """Return the tables by variable and the joint patterns as text."""
        return (
            self.by_variable.round(2).to_string()
            + "\n\n"
            + self.patterns.round(2).to_string(index=False)
        )


@description(
    "missing_report",
    DescribeSpec(
        stata=None,
        summary="Missing values by variable, period and unit, with joint patterns.",
    ),
)
def missing_report(panel: PanelBase, variables: Sequence[str] | None = None) -> MissingReport:
    """Report missing values by variable, by period, by unit and jointly.

    Parameters
    ----------
    panel
        The declared panel.
    variables
        Columns to examine (default: every column except the two keys).

    Returns
    -------
    MissingReport
        Result of :func:`missing_report`.
    """
    names, flags = _flags(panel, variables)
    index = panel.index
    missing = flags.to_numpy()
    n_rows = missing.shape[0]
    per_unit = np.column_stack(
        [
            np.bincount(index.unit, weights=missing[:, j], minlength=index.n_units)
            for j in range(len(names))
        ]
    )
    t_i = index.counts
    by_variable = pd.DataFrame(
        {
            "n_obs": n_rows,
            "n_missing": missing.sum(axis=0),
            "pct_missing": 100.0 * missing.mean(axis=0),
            "units_all_missing": (per_unit == t_i[:, None]).sum(axis=0),
            "units_some_missing": ((per_unit > 0) & (per_unit < t_i[:, None])).sum(axis=0),
        },
        index=pd.Index(names, name="variable"),
    )
    periods, code = np.unique(index.rank, return_inverse=True)
    per_period = np.column_stack(
        [
            np.bincount(code, weights=missing[:, j], minlength=periods.size)
            for j in range(len(names))
        ]
    )
    n_period = np.bincount(code, minlength=periods.size)
    labels = pd.Index(panel.axis.labels(periods), name=panel.time_name)
    by_period = pd.DataFrame(per_period / n_period[:, None], index=labels, columns=names)
    by_period.insert(0, "n_obs", n_period)
    by_unit = pd.DataFrame(
        per_unit.astype(np.int64),
        index=pd.Index(index.unit_labels, name=panel.id_name),
        columns=names,
    )
    by_unit.insert(0, "T_i", t_i)
    keys = ["".join("." if m else "1" for m in row) for row in missing]
    counts = pd.Series(keys).value_counts(sort=False)
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    patterns = pd.DataFrame(ordered, columns=["pattern", "freq"])
    patterns["percent"] = 100.0 * patterns["freq"] / n_rows
    patterns["cum_percent"] = patterns["percent"].cumsum()
    return MissingReport(
        by_variable=by_variable, by_period=by_period, by_unit=by_unit, patterns=patterns
    )


@dataclass(frozen=True, kw_only=True)
class CompleteCases(Result):
    """Result of :func:`complete_cases`.

    Attributes
    ----------
    mask
        Boolean Series on the index of the declared frame: the row is in the panel and all the
        variables are observed.
    n_rows, n_rows_kept
        Rows of the panel before and after.
    n_units, n_units_kept
        Units before and after.
    ti
        Distribution of ``T_i`` before and after (columns ``before`` and ``after``; rows are the
        statistics ``mean``, ``min``, the percentiles and ``max``), for the units kept.
    panel
        The sub-panel of the complete rows.
    """

    mask: pd.Series
    n_rows: int
    n_rows_kept: int
    n_units: int
    n_units_kept: int
    ti: pd.DataFrame
    panel: Any

    def to_frame(self) -> pd.DataFrame:
        """Return the distribution of ``T_i`` before and after."""
        return self.ti

    def summary(self, lang: Lang | None = None) -> str:
        """Return a short report."""
        return (
            f"rows: {self.n_rows} -> {self.n_rows_kept}   "
            f"units: {self.n_units} -> {self.n_units_kept}\n" + self.ti.round(2).to_string()
        )


@description(
    "complete_cases",
    DescribeSpec(
        stata=None,
        summary="Rows and units kept when all of some variables must be observed, and T_i effect.",
    ),
)
def complete_cases(panel: PanelBase, variables: Sequence[str]) -> CompleteCases:
    """Keep the rows where all ``variables`` are observed and show the effect on ``T_i``.

    Parameters
    ----------
    panel
        The declared panel.
    variables
        Variables that must all be observed for a row to be kept.

    Returns
    -------
    CompleteCases
        Result of :func:`complete_cases`.
    """
    _, flags = _flags(panel, variables)
    keep_sorted = ~flags.to_numpy().any(axis=1)
    index = panel.index
    sub = panel.select(complete=list(variables)) if keep_sorted.any() else None
    mask = pd.Series(False, index=panel.frame.index)
    mask.iloc[index.source_row[keep_sorted]] = True
    kept = np.bincount(index.unit[keep_sorted], minlength=index.n_units)
    before = index.counts.astype(float)
    after = kept[kept > 0].astype(float)

    def stats(values: np.ndarray) -> dict[str, float]:
        out = {"mean": float(values.mean()) if values.size else float("nan")}
        out["min"] = float(values.min()) if values.size else float("nan")
        for p in PERCENTILES:
            out[f"p{p}"] = stata_quantile(values, p) if values.size else float("nan")
        out["max"] = float(values.max()) if values.size else float("nan")
        return out

    return CompleteCases(
        mask=mask,
        n_rows=int(index.n_obs),
        n_rows_kept=int(keep_sorted.sum()),
        n_units=int(index.n_units),
        n_units_kept=int((kept > 0).sum()),
        # Series, not dicts: the rows keep the order mean, min, p..., max on every pandas version
        ti=pd.DataFrame({"before": pd.Series(stats(before)), "after": pd.Series(stats(after))}),
        panel=sub,
    )
