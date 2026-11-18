"""``xttab``: overall, between and within frequency of the values of a variable (E.3.1).

Definitions of Stata (ANALYSE_STATA.md, section 4), for each value ``v``:

* overall: observations equal to ``v`` and their share of all observations;
* between: units that take ``v`` at least once and their share of the units (the shares add up
  to more than 100);
* within: average, over the units that take ``v``, of the share of their observations equal to
  ``v``.  The total row weights the within shares by the between frequencies.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.describe._format import number
from pyxtpanel.describe._spec import DescribeSpec, description
from pyxtpanel.ops.frame import term_columns

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

#: More distinct values than this is almost surely a continuous variable.
MAX_VALUES = 1000

_i18n.register(
    "describe",
    {
        "describe.xt.too_many": {
            "en": "{name!r} has {count} distinct values (limit {limit}); discretize it first or "
            "raise max_values.",
            "fr": "« {name} » a {count} valeurs distinctes (limite {limit}) ; discrétisez-la ou "
            "augmentez max_values.",
        },
        "describe.xt.one_variable": {
            "en": "{command} takes exactly one variable, got {count}.",
            "fr": "{command} exige exactement une variable, reçu {count}.",
        },
        "describe.xt.no_observations": {
            "en": "{name!r} has no observed value.",
            "fr": "« {name} » n'a aucune valeur observée.",
        },
        "describe.xt.header": {
            "en": "           |       Overall      |      Between      |  Within\n"
            "        {name:<3}|  Freq.   Percent   |  Freq.   Percent  |  Percent",
            "fr": "           |       Globale      |       Entre       |  Intra\n"
            "        {name:<3}|  Eff.    Pourc.    |  Eff.    Pourc.   |  Pourc.",
        },
    },
)


@dataclass(frozen=True, kw_only=True)
class XtTab(Result):
    """Result of :func:`xttab`.

    Attributes
    ----------
    variable
        Name of the tabulated variable.
    table
        One row per value plus a ``Total`` row: ``overall_freq``, ``overall_pct``,
        ``between_freq``, ``between_pct``, ``within_pct`` (Stata's ``r(results)``).
    n_units
        Units with at least one observed value (``r(n_panels)``).
    n_obs
        Observed values.
    """

    variable: str
    table: pd.DataFrame
    n_units: int
    n_obs: int

    def to_frame(self) -> pd.DataFrame:
        """Return the table."""
        return self.table

    def summary(self, lang: Lang | None = None) -> str:
        """Return the table in the layout of Stata's ``xttab``."""
        lines = _i18n.message("describe.xt.header", lang, name=self.variable[:3]).split("\n")
        lines.append("-----------+" + "-" * 20 + "+" + "-" * 19 + "+" + "-" * 10)
        for key, row in self.table.iterrows():
            label = number(float(str(key)), 6) if key != "Total" else "Total"
            lines.append(
                f"{label:>10} |{row['overall_freq']:>7.0f}{row['overall_pct']:>10.2f}   |"
                f"{row['between_freq']:>7.0f}{row['between_pct']:>10.2f}  |"
                f"{row['within_pct']:>8.2f}"
            )
        return "\n".join(lines)


def tabulate(
    unit: np.ndarray, n_units: int, x: np.ndarray, max_values: int = MAX_VALUES, name: str = "x"
) -> pd.DataFrame:
    """Overall, between and within table of the observed values of ``x``."""
    keep = ~np.isnan(x)
    if not keep.any():
        raise DataError(
            _i18n.message("describe.xt.no_observations", name=name),
            stata_code=2000,
            key="describe.xt.no_observations",
        )
    values, codes = x[keep], unit[keep]
    levels, which = np.unique(values, return_inverse=True)
    if levels.size > max_values:
        raise SpecificationError(
            _i18n.message("describe.xt.too_many", name=name, count=levels.size, limit=max_values),
            stata_code=134,
            key="describe.xt.too_many",
        )
    k = levels.size
    cells = np.bincount(codes * k + which, minlength=n_units * k).reshape(n_units, k)
    t_i = cells.sum(axis=1)
    present = t_i > 0
    cells, t_i = cells[present], t_i[present]
    n = float(cells.shape[0])
    overall = cells.sum(axis=0).astype(np.float64)
    between = (cells > 0).sum(axis=0).astype(np.float64)
    share = cells / t_i[:, None]
    within = np.where(between > 0, share.sum(axis=0) / np.maximum(between, 1), 0.0) * 100.0
    # Stata: xttab, r(results)
    # STATA-QUIRK: the total within percent weights the values by their between frequencies
    total_between = between.sum()
    table = pd.DataFrame(
        {
            "overall_freq": [*overall, overall.sum()],
            "overall_pct": [*(100 * overall / overall.sum()), 100.0],
            "between_freq": [*between, total_between],
            "between_pct": [*(100 * between / n), 100 * total_between / n],
            "within_pct": [*within, float(np.sum(between / total_between * within))],
        },
        index=pd.Index([*levels, "Total"], name="value"),
    )
    return table


@description(
    "xttab",
    DescribeSpec(
        stata="xttab",
        summary="Overall, between and within frequencies of the values of a variable.",
    ),
)
def xttab(panel: PanelBase, variable: str, *, max_values: int = MAX_VALUES) -> XtTab:
    """Tabulate a variable overall, between and within units, as Stata's ``xttab``.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        A column (or a Stata-style term such as ``"L.x"``); missing values are left out.
    max_values
        Refuse variables with more distinct values than this.

    Returns
    -------
    XtTab
        Result of :func:`xttab`.

    Raises
    ------
    SpecificationError
        If the variable has too many distinct values (Stata code 134).
    DataError
        If the variable is unknown or not numeric.
    """
    names, values = term_columns(panel, variable)
    if len(names) != 1:
        raise SpecificationError(
            _i18n.message("describe.xt.one_variable", command="xttab", count=len(names)),
            key="describe.xt.one_variable",
        )
    table = tabulate(panel.index.unit, panel.n_units, values[:, 0], max_values, names[0])
    n_obs = int(table["overall_freq"].iloc[-1])
    seen = np.unique(panel.index.unit[~np.isnan(values[:, 0])]).size
    return XtTab(variable=names[0], table=table, n_units=int(seen), n_obs=n_obs)
