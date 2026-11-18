"""``xtsum``: overall, between and within summary statistics (E.2).

Conventions of Stata (ANALYSE_STATA.md, section 3), unbalanced panels included:

* overall: mean and standard deviation of all observations (divisor ``N − 1``);
* between: standard deviation of the unit means ``x̄_i``, **not weighted by ``T_i``** (divisor
  ``n − 1``), minimum and maximum of the ``x̄_i``;
* within: ``√(Σ (x_it − x̄_i)² / (N − 1))`` (divisor ``N − 1``, not ``N − n``); minimum and
  maximum of ``x_it − x̄_i + x̄``, the overall mean ``x̄`` being added back.
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
from pyxtpanel._types import FloatArray, IntArray, Lang
from pyxtpanel.describe._format import number
from pyxtpanel.describe._spec import DescribeSpec, description
from pyxtpanel.ops.frame import term_columns

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "describe",
    {
        "describe.xs.no_variable": {
            "en": "xtsum needs at least one variable.",
            "fr": "xtsum exige au moins une variable.",
        },
        "describe.xs.by_unknown": {
            "en": "by= must name a column of the data, got {name!r}.",
            "fr": "by= doit nommer une colonne des données, reçu {name!r}.",
        },
        "describe.xs.header": {
            "en": "Variable         |      Mean   Std. dev.       Min        Max |    Observations",
            "fr": "Variable         |   Moyenne  Écart type       Min        Max |    Observations",
        },
    },
)


@dataclass(frozen=True)
class SumStats:
    """The 13 numbers of ``xtsum`` for one variable (``r()`` of Stata)."""

    N: float
    n: float
    Tbar: float
    mean: float
    sd: float
    min: float
    max: float
    sd_b: float
    min_b: float
    max_b: float
    sd_w: float
    min_w: float
    max_w: float


def sum_stats(unit: IntArray, n_units: int, x: FloatArray) -> SumStats:
    """Overall, between and within statistics of one variable (missing values ignored)."""
    keep = ~np.isnan(x)
    values, codes = x[keep], unit[keep]
    nan = float("nan")
    count = float(values.size)
    if count == 0:
        return SumStats(0.0, 0.0, nan, *([nan] * 10))
    size = np.bincount(codes, minlength=n_units)
    present = size > 0
    total = np.bincount(codes, weights=values, minlength=n_units)
    unit_mean = total[present] / size[present]
    n = float(present.sum())
    mean = float(values.mean())
    sd = float(values.std(ddof=1)) if count > 1 else nan
    sd_b = float(unit_mean.std(ddof=1)) if n > 1 else nan
    # STATA-QUIRK: xtsum within uses the divisor N - 1 (not N - n); min, max add back the mean
    deviation = values - (total / np.maximum(size, 1))[codes]
    sd_w = float(np.sqrt(np.sum(deviation**2) / (count - 1))) if count > 1 else nan
    return SumStats(
        N=count,
        n=n,
        Tbar=count / n,
        mean=mean,
        sd=sd,
        min=float(np.min(values)),
        max=float(np.max(values)),
        sd_b=sd_b,
        min_b=float(np.min(unit_mean)),
        max_b=float(np.max(unit_mean)),
        sd_w=sd_w,
        min_w=float(np.min(deviation)) + mean,
        max_w=float(np.max(deviation)) + mean,
    )


@dataclass(frozen=True, kw_only=True)
class XtSum(Result):
    """Result of :func:`xtsum`.

    Attributes
    ----------
    table
        Rows ``(variable, overall | between | within)`` (preceded by the group with ``by=``),
        columns ``mean``, ``sd``, ``min``, ``max`` and ``obs`` (``N``, ``n`` or ``T-bar``).
    scalars
        One row per variable (and group) with the 13 values Stata returns in ``r()``, plus
        ``within_share = s_w² / (s_b² + s_w²)``, a diagnostic that Stata does not compute.
    by
        Name of the grouping column, if any.
    """

    table: pd.DataFrame
    scalars: pd.DataFrame
    by: str | None = None

    def to_frame(self) -> pd.DataFrame:
        """Return the table (variable × component by statistic)."""
        return self.table

    def summary(self, lang: Lang | None = None) -> str:
        """Return the table in the layout of Stata's ``xtsum``."""
        rule = "-" * 17 + "+" + "-" * 45 + "+" + "-" * 16
        lines = [_i18n.message("describe.xs.header", lang), rule]
        for key, row in self.table.iterrows():
            parts = key if isinstance(key, tuple) else (key,)
            name, kind = str(parts[-2]), str(parts[-1])
            head = f"{name[:8]:<8} {kind:<8}" if kind == "overall" else f"{'':<8} {kind:<8}"
            cells = ("mean", "sd", "min", "max")
            stats = "".join(
                f"{number(row[c]) if c != 'mean' or kind == 'overall' else '':>10}" for c in cells
            )
            tag = {"overall": "N =", "between": "n =", "within": "T-bar ="}[kind]
            obs = f"{row['obs']:.4g}" if kind == "within" else f"{row['obs']:.0f}"
            lines.append(f"{head:<17}|{stats} | {tag:>7} {obs:>8}")
        return "\n".join(lines)


def _block(
    name: str, s: SumStats, prefix: tuple[Any, ...]
) -> list[tuple[tuple[Any, ...], dict[str, float]]]:
    nan = float("nan")
    return [
        ((*prefix, name, "overall"), dict(mean=s.mean, sd=s.sd, min=s.min, max=s.max, obs=s.N)),
        ((*prefix, name, "between"), dict(mean=nan, sd=s.sd_b, min=s.min_b, max=s.max_b, obs=s.n)),
        (
            (*prefix, name, "within"),
            dict(mean=nan, sd=s.sd_w, min=s.min_w, max=s.max_w, obs=s.Tbar),
        ),
    ]


@description(
    "xtsum",
    DescribeSpec(
        stata="xtsum",
        summary="Overall, between and within mean, standard deviation, min and max.",
    ),
)
def xtsum(panel: PanelBase, variables: str | Sequence[str], *, by: str | None = None) -> XtSum:
    """Summarize variables overall, between and within units, as Stata's ``xtsum``.

    Parameters
    ----------
    panel
        The declared panel.
    variables
        Column names, or Stata-style terms (``"L.x"``, ``"D.x"``).
    by
        Column that defines groups; each group is summarized on its own rows.

    Returns
    -------
    XtSum
        Result of :func:`xtsum`.

    Raises
    ------
    SpecificationError
        If no variable is given or ``by`` is not a column.
    DataError
        If a variable is unknown or not numeric.
    """
    wanted = [variables] if isinstance(variables, str) else list(variables)
    if not wanted:
        raise SpecificationError(
            _i18n.message("describe.xs.no_variable"), key="describe.xs.no_variable"
        )
    groups: list[tuple[tuple[Any, ...], PanelBase]] = [((), panel)]
    if by is not None:
        if by not in panel.frame.columns:
            raise SpecificationError(
                _i18n.message("describe.xs.by_unknown", name=by), key="describe.xs.by_unknown"
            )
        column = panel.frame[by].reset_index(drop=True)
        groups = []
        for value in sorted(column.dropna().unique().tolist()):
            groups.append(((value,), panel.select((column == value).to_numpy())))
    rows: list[tuple[tuple[Any, ...], dict[str, float]]] = []
    scalar_rows: list[tuple[tuple[Any, ...], dict[str, float]]] = []
    for prefix, sub in groups:
        names, values = term_columns(sub, wanted)
        for j, name in enumerate(names):
            stats = sum_stats(sub.index.unit, sub.n_units, values[:, j])
            rows += _block(name, stats, prefix)
            spread = stats.sd_b**2 + stats.sd_w**2
            share = stats.sd_w**2 / spread if spread > 0 else float("nan")
            scalar_rows.append(((*prefix, name), {**vars(stats), "within_share": share}))
    level_names = ([by] if by else []) + ["variable"]
    table = pd.DataFrame(
        [r for _, r in rows],
        index=pd.MultiIndex.from_tuples([k for k, _ in rows], names=[*level_names, "component"]),
        columns=["mean", "sd", "min", "max", "obs"],
    )
    keys = [k for k, _ in scalar_rows]
    index = (
        pd.MultiIndex.from_tuples(keys, names=level_names)
        if by
        else pd.Index([k[0] for k in keys], name="variable")
    )
    scalars = pd.DataFrame([r for _, r in scalar_rows], index=index)
    return XtSum(table=table, scalars=scalars, by=by)
