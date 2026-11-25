"""Where the variation of each variable comes from (K.3.1) and what is left of it (K.3.3).

An addition of the package (Stata's ``xtsum`` gives standard deviations, not the count of units
or periods that vary; R's ``pvar`` and Python's ``are_varying`` give a yes/no per variable);
validated by hand on small panels, against an independent loop implementation and by invariance
tests.

:func:`variation` classifies each variable:

* ``constant``: one observed value in the whole panel;
* ``time_invariant``: constant within every unit, so only the differences *between* units
  exist (a fixed effect absorbs it);
* ``common``: constant across the units of every period, so only the changes *over time*
  exist (a period effect absorbs it);
* ``varying``: varies within units and across units.

:func:`residual_variation` applies the within, first-difference and forward-orthogonal-deviation
transformations of :mod:`pyxtpanel.ops` and reports what each leaves: observations, units with
a non-zero value, standard deviation and share of the variance of the variable, to choose the
estimator knowing what it can still use.  Values are compared exactly (float64), missing values
are ignored.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import FloatArray, Lang
from pyxtpanel.explore._spec import ToolSpec, tool
from pyxtpanel.ops.diff import difference, fod
from pyxtpanel.ops.within import within
from pyxtpanel.panel.arrays import sorted_columns

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

#: Transformations of :func:`residual_variation`, in the order reported.
TRANSFORMS = ("within", "fd", "fod")

_i18n.register(
    "explore",
    {
        "explore.var.title": {
            "en": "Variation of {n} variable(s): {units} units, {periods} periods",
            "fr": "Variation de {n} variable(s) : {units} unités, {periods} périodes",
        },
        "explore.var.resid": {
            "en": "Variation left by each transformation (share of the variance of the variable):",
            "fr": "Variation restante après chaque transformation (part de la variance de la "
            "variable) :",
        },
        "explore.var.no_numeric": {
            "en": "the panel has no numeric variable besides its keys.",
            "fr": "le panel n'a aucune variable numérique en dehors de ses clés.",
        },
        "explore.var.transform": {
            "en": "unknown transformation {value!r}; use within, fd or fod.",
            "fr": "transformation inconnue {value!r} ; utilisez within, fd ou fod.",
        },
    },
)


def numeric_names(panel: PanelBase) -> list[str]:
    """Numeric variables of the panel other than its two keys."""
    names = list(panel.numeric_names)
    if not names:
        raise SpecificationError(
            _i18n.message("explore.var.no_numeric"), key="explore.var.no_numeric"
        )
    return names


@dataclass(frozen=True, kw_only=True)
class Variation(Result):
    """Result of :func:`variation`.

    Attributes
    ----------
    table
        One row per variable: ``n_obs``, ``n_distinct``, ``units`` (with an observed value),
        ``units_varying`` (more than one value), ``units_single`` (exactly one observation, so
        no within information), ``share_units_varying``, ``periods`` (with an observed value),
        ``periods_varying`` (the units differ in that period) and ``kind``.
    n_units, n_periods
        Size of the panel.
    """

    table: pd.DataFrame
    n_units: int
    n_periods: int

    def to_frame(self) -> pd.DataFrame:
        """Return the table by variable."""
        return self.table

    def summary(self, lang: Lang | None = None) -> str:
        """Return the table as text."""
        head = _i18n.message(
            "explore.var.title",
            lang,
            n=len(self.table),
            units=self.n_units,
            periods=self.n_periods,
        )
        return head + "\n\n" + self.table.to_string(float_format=lambda v: f"{v:.3f}")


@tool(
    "variation",
    ToolSpec(
        summary="For each variable: units and periods that vary, and whether it is constant, "
        "time-invariant or common to all units.",
        validation="hand-computed panels, independent loop implementation, invariance tests",
        similar="plm::pvar, linearmodels are_varying",
    ),
)
def variation(panel: PanelBase, variables: Sequence[str] | None = None) -> Variation:
    """Tell where each variable varies: across units, over time, within units.

    Parameters
    ----------
    panel
        A declared panel.
    variables
        Variables to examine; by default every numeric variable other than the keys.

    Returns
    -------
    Variation
        Result of :func:`variation`.
    """
    names = list(variables) if variables is not None else numeric_names(panel)
    index = panel.index
    data = sorted_columns(panel, names)
    low = int(np.min(index.rank))
    period = index.rank - low
    size = int(np.max(index.rank)) - low + 1
    rows = []
    for j in range(len(names)):
        x = data[:, j]
        seen = ~np.isnan(x)
        lo_u = np.full(index.n_units, np.inf)
        hi_u = np.full(index.n_units, -np.inf)
        np.minimum.at(lo_u, index.unit[seen], x[seen])
        np.maximum.at(hi_u, index.unit[seen], x[seen])
        count_u = np.bincount(index.unit[seen], minlength=index.n_units)
        has = count_u > 0
        varies_u = has & (hi_u > lo_u)
        lo_t = np.full(size, np.inf)
        hi_t = np.full(size, -np.inf)
        np.minimum.at(lo_t, period[seen], x[seen])
        np.maximum.at(hi_t, period[seen], x[seen])
        has_t = np.bincount(period[seen], minlength=size) > 0
        varies_t = has_t & (hi_t > lo_t)
        n_distinct = int(np.unique(x[seen]).size)
        if n_distinct <= 1:
            kind = "constant"
        elif not varies_u.any():
            kind = "time_invariant"
        elif not varies_t.any():
            kind = "common"
        else:
            kind = "varying"
        rows.append(
            {
                "n_obs": int(seen.sum()),
                "n_distinct": n_distinct,
                "units": int(has.sum()),
                "units_varying": int(varies_u.sum()),
                "units_single": int((count_u == 1).sum()),
                "share_units_varying": float(np.sum(varies_u)) / float(np.sum(has))
                if has.any()
                else np.nan,
                "periods": int(has_t.sum()),
                "periods_varying": int(varies_t.sum()),
                "kind": kind,
            }
        )
    table = pd.DataFrame(rows, index=pd.Index(names, name="variable"))
    return Variation(table=table, n_units=index.n_units, n_periods=size)


@dataclass(frozen=True, kw_only=True)
class ResidualVariation(Result):
    """Result of :func:`residual_variation`.

    Attributes
    ----------
    table
        One row per (variable, transformation): ``n`` (observations with a value),
        ``units`` (with at least one non-zero value), ``sd``, ``variance_share`` (variance of
        the transformed values over the variance of the variable) and ``usable`` (``n`` and
        ``sd`` are positive).
    """

    table: pd.DataFrame

    def to_frame(self) -> pd.DataFrame:
        """Return the table by variable and transformation."""
        return self.table

    def summary(self, lang: Lang | None = None) -> str:
        """Return the table as text."""
        return (
            _i18n.message("explore.var.resid", lang)
            + "\n\n"
            + self.table.to_string(float_format=lambda v: f"{v:.4f}")
        )


def _transformed(panel: PanelBase, x: FloatArray, kind: str) -> FloatArray:
    if kind == "within":
        return within(panel.index, x)
    if kind == "fd":
        return difference(panel.index, x)
    return fod(panel.index, x)


@tool(
    "residual_variation",
    ToolSpec(
        summary="Observations, units and variance left by the within, first-difference and "
        "orthogonal-deviation transformations.",
        validation="hand-computed panels, independent loop implementation, invariance tests",
    ),
)
def residual_variation(
    panel: PanelBase,
    variables: Sequence[str],
    transforms: Sequence[str] = TRANSFORMS,
) -> ResidualVariation:
    """Measure what is left of each variable once the unit effect is removed.

    A value counts as non-zero when it exceeds ``1e-12`` times the largest absolute value of the
    variable, so that rounding error does not make a constant look varying.

    Parameters
    ----------
    panel
        The declared panel.
    variables
        Variables to examine.
    transforms
        Transformations whose residual variation is measured (``"within"``, ``"fd"``, ``"fod"``).

    Returns
    -------
    ResidualVariation
        Result of :func:`residual_variation`.

    Raises
    ------
    SpecificationError
        If a transformation is not one of those listed.
    """
    for kind in transforms:
        if kind not in TRANSFORMS:
            raise SpecificationError(
                _i18n.message("explore.var.transform", value=kind),
                stata_code=198,
                key="explore.var.transform",
            )
    index = panel.index
    data = sorted_columns(panel, list(variables))
    out = {}
    for j, name in enumerate(variables):
        x = data[:, j]
        seen = ~np.isnan(x)
        scale = float(np.max(np.abs(x[seen]))) if seen.any() else 0.0
        total = float(np.var(x[seen], ddof=1)) if seen.sum() > 1 else np.nan
        for kind in transforms:
            moved = _transformed(panel, x, kind)
            ok = ~np.isnan(moved)
            nonzero = ok & (np.abs(moved) > 1e-12 * max(scale, 1.0))
            sd = float(np.std(moved[ok], ddof=1)) if ok.sum() > 1 else np.nan
            share = sd**2 / total if total and total > 0 else np.nan
            out[(name, kind)] = {
                "n": int(ok.sum()),
                "units": int(np.unique(index.unit[nonzero]).size),
                "sd": sd,
                "variance_share": share,
                "usable": bool(ok.sum() > 0 and sd > 0),
            }
    table = pd.DataFrame.from_dict(out, orient="index")
    table.index = pd.MultiIndex.from_tuples(table.index, names=["variable", "transform"])
    return ResidualVariation(table=table)
