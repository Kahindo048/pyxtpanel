"""Correlations of panel variables: overall, within and between, and within autocorrelation (E.4).

Each pair of variables uses the rows where **both** are observed (pairwise deletion, as
``pwcorr``); the count of each pair is returned with the correlation.

* overall: Pearson correlation of the pair's observations;
* within: correlation of the deviations from the unit means, the means taken over the pair's
  rows (what ``xtdata, fe`` followed by ``correlate`` gives on the pair's sample);
* between: correlation of the unit means over the pair's rows, one point per unit, unweighted
  (``xtdata, be``).

``autocorr`` is an addition of the package (Stata has no equivalent).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import FloatArray, IntArray, Lang
from pyxtpanel.describe._spec import DescribeSpec, description
from pyxtpanel.ops.frame import term_columns
from pyxtpanel.ops.groups import group_mean
from pyxtpanel.ops.lags import lag

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

Kind = Literal["overall", "within", "between"]
KINDS: tuple[Kind, ...] = ("overall", "within", "between")

_i18n.register(
    "describe",
    {
        "describe.xc.few_variables": {
            "en": "xtcorr needs at least two variables.",
            "fr": "xtcorr exige au moins deux variables.",
        },
        "describe.xc.bad_kind": {
            "en": "kind must be one of overall, within, between or all, got {kind!r}.",
            "fr": "kind doit valoir overall, within, between ou all, reçu {kind!r}.",
        },
        "describe.xc.bad_lags": {
            "en": "lags must be positive integers, got {lags!r}.",
            "fr": "lags doit contenir des entiers positifs, reçu {lags!r}.",
        },
    },
)


#: A variable whose deviations are below this fraction of its own scale is constant up to rounding.
RELATIVE_ZERO = 1e-12


def pearson(a: FloatArray, b: FloatArray, scale: tuple[float, float] | None = None) -> float:
    """Pearson correlation of two samples.

    Missing when there are fewer than 2 points or when a sample is constant, up to rounding: a
    sample whose deviations are within ``RELATIVE_ZERO`` of its scale has no correlation (the
    subtraction of a mean from a constant leaves noise of about 1e-16, which must not be
    reported as a correlation).

    Parameters
    ----------
    scale
        Scale (largest absolute value) of each variable before its means were removed; by default
        that of the samples themselves.
    """
    if a.size < 2:
        return float("nan")
    scale_a, scale_b = scale or (float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    da, db = a - a.mean(), b - b.mean()
    ss_a, ss_b = float(np.sum(da**2)), float(np.sum(db**2))
    if (
        ss_a <= a.size * (RELATIVE_ZERO * scale_a) ** 2
        or ss_b <= a.size * (RELATIVE_ZERO * scale_b) ** 2
    ):
        return float("nan")
    return float(np.sum(da * db) / np.sqrt(ss_a * ss_b))


def unit_means(unit: IntArray, n_units: int, x: FloatArray) -> tuple[FloatArray, np.ndarray]:
    """Mean of each unit over the given rows and the mask of units that have rows."""
    size = np.bincount(unit, minlength=n_units)
    total = np.bincount(unit, weights=x, minlength=n_units)
    return total / np.maximum(size, 1), size > 0


def pair_correlations(
    unit: IntArray, n_units: int, a: FloatArray, b: FloatArray, kinds: Sequence[Kind]
) -> dict[Kind, tuple[float, int]]:
    """Correlations of one pair on its common rows, and the number of points each uses.

    The rows where both variables are observed are selected once and shared by the kinds.
    """
    both = ~np.isnan(a) & ~np.isnan(b)
    a, b, u = a[both], b[both], unit[both]
    raw = (float(np.max(np.abs(a))), float(np.max(np.abs(b)))) if a.size else (0.0, 0.0)
    out: dict[Kind, tuple[float, int]] = {}
    means: tuple[FloatArray, FloatArray, np.ndarray] | None = None
    for kind in kinds:
        if kind == "overall":
            out[kind] = (pearson(a, b), int(a.size))
            continue
        if means is None:
            mean_a, present = unit_means(u, n_units, a)
            means = (mean_a, unit_means(u, n_units, b)[0], present)
        mean_a, mean_b, present = means
        if kind == "between":
            out[kind] = (
                pearson(mean_a[present], mean_b[present], raw),
                int(present.sum()),
            )
        else:
            out[kind] = (pearson(a - mean_a[u], b - mean_b[u], raw), int(a.size))
    return out


@dataclass(frozen=True, kw_only=True)
class XtCorr(Result):
    """Result of :func:`xtcorr`.

    Attributes
    ----------
    variables
        Names of the variables.
    corr
        ``{kind: DataFrame}`` correlation matrices.
    n
        ``{kind: DataFrame}`` points behind each correlation (observations for overall and
        within, units for between).
    """

    variables: tuple[str, ...]
    corr: dict[str, pd.DataFrame]
    n: dict[str, pd.DataFrame]

    def to_frame(self) -> pd.DataFrame:
        """Return the correlations, one block of rows per kind."""
        return pd.concat(self.corr, names=["kind", "variable"])

    def summary(self, lang: Lang | None = None) -> str:
        """Return the matrices as text, one per kind."""
        blocks = []
        for kind, matrix in self.corr.items():
            blocks.append(f"{kind}\n{matrix.round(4).to_string(na_rep='.')}")
        return "\n\n".join(blocks)


@description(
    "xtcorr",
    DescribeSpec(
        stata=None,
        summary="Overall, within and between correlations with the number of points per pair.",
    ),
)
def xtcorr(
    panel: PanelBase,
    variables: Sequence[str],
    *,
    kind: Kind | Literal["all"] = "all",
) -> XtCorr:
    """Correlate variables overall, within units and between units.

    Parameters
    ----------
    panel
        The declared panel.
    variables
        At least two columns or Stata-style terms.
    kind
        ``"overall"``, ``"within"``, ``"between"`` or ``"all"`` (default).

    Raises
    ------
    SpecificationError
        If fewer than two variables or an invalid ``kind``.
    """
    if kind != "all" and kind not in KINDS:
        raise SpecificationError(
            _i18n.message("describe.xc.bad_kind", kind=kind), key="describe.xc.bad_kind"
        )
    names, values = term_columns(panel, list(variables))
    if len(names) < 2:
        raise SpecificationError(
            _i18n.message("describe.xc.few_variables"), key="describe.xc.few_variables"
        )
    wanted = KINDS if kind == "all" else (kind,)
    k = len(names)
    r = {one: np.full((k, k), np.nan) for one in wanted}
    n = {one: np.zeros((k, k)) for one in wanted}
    for i in range(k):
        for j in range(i, k):
            pair = pair_correlations(
                panel.index.unit, panel.n_units, values[:, i], values[:, j], wanted
            )
            for one in wanted:
                r[one][i, j] = r[one][j, i] = pair[one][0]
                n[one][i, j] = n[one][j, i] = pair[one][1]
    corr: dict[str, pd.DataFrame] = {
        one: pd.DataFrame(r[one], index=names, columns=names) for one in wanted
    }
    count: dict[str, pd.DataFrame] = {
        one: pd.DataFrame(n[one], index=names, columns=names) for one in wanted
    }
    return XtCorr(variables=tuple(names), corr=corr, n=count)


@dataclass(frozen=True, kw_only=True)
class AutoCorr(Result):
    """Result of :func:`autocorr`: within autocorrelation by lag (``rho`` and ``n`` pairs)."""

    variable: str
    table: pd.DataFrame

    def to_frame(self) -> pd.DataFrame:
        """Return the table (index: lag; columns ``rho`` and ``n``)."""
        return self.table

    def summary(self, lang: Lang | None = None) -> str:
        """Return the table as text."""
        return f"{self.variable}\n{self.table.to_string(na_rep='.')}"


@description(
    "autocorr",
    DescribeSpec(
        stata=None,
        summary="Within autocorrelation by calendar lag (addition of the package).",
    ),
)
def autocorr(panel: PanelBase, variable: str, lags: int | Sequence[int] = 1) -> AutoCorr:
    """Autocorrelation of the within deviations ``x_it − x̄_i`` at calendar lags.

    The correlation is Pearson's on the pooled pairs ``(w_it, w_i,t−k)`` where both periods are
    observed; a gap in the calendar removes the pair instead of pairing distant observations.
    """
    wanted = [lags] if isinstance(lags, (int, np.integer)) else list(lags)
    if not wanted or any(int(k) < 1 for k in wanted):
        raise SpecificationError(
            _i18n.message("describe.xc.bad_lags", lags=lags), key="describe.xc.bad_lags"
        )
    names, values = term_columns(panel, variable)
    x = values[:, 0]
    deviation = x - group_mean(panel.index, x, "unit", at_missing=True)
    rows = []
    for k in wanted:
        earlier = lag(panel.index, deviation, int(k))
        both = ~np.isnan(deviation) & ~np.isnan(earlier)
        rows.append((int(k), pearson(deviation[both], earlier[both]), int(both.sum())))
    table = pd.DataFrame(rows, columns=["lag", "rho", "n"]).set_index("lag")
    return AutoCorr(variable=names[0], table=table)
