"""Fisher-type panel unit-root test built from per-unit ADF p-values (category A, L.2.5).

Validated against Stata 17: ``xtunitroot fisher, dfuller`` with fixed lags, a constant, a trend
or a drift, with and without cross-sectional means removed, on a balanced panel, on panels with
entry, exit and holes and on a panel with units too short to be tested; the four statistics (P, Z,
L and Pm), their degrees of freedom and p-values are compared, and so are the per-unit p-values
that Stata returns through ``generate()``.  The cases are in
``validation/cases/diagnostics_unitroot_fisher.py``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.adf import AdfCase, AdfFit, adf_fit, adf_pvalue

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.fisher.title": {
            "en": "Fisher-type unit-root test (augmented Dickey-Fuller)",
            "fr": "Test de racine unitaire de type Fisher (Dickey-Fuller augmenté)",
        },
        "diag.fisher.null": {
            "en": "all panels contain unit roots",
            "fr": "toutes les séries du panel ont une racine unitaire",
        },
        "diag.fisher.alternative": {
            "en": "at least one panel is stationary",
            "fr": "au moins une série du panel est stationnaire",
        },
        "diag.fisher.skipped": {
            "en": "{count} unit(s) were left out: their ADF regression has no equation or no "
            "p-value (too few consecutive observations for {lags} lag(s)).",
            "fr": "{count} unité(s) ont été écartées : leur régression ADF n'a aucune équation "
            "ou aucune p-value (trop peu d'observations consécutives pour {lags} retard(s)).",
        },
        "diag.fisher.bad_method": {
            "en": "method must be 'inverse_chi2', 'inverse_normal', 'logit' or 'modified', "
            "got {value!r}.",
            "fr": "method doit valoir 'inverse_chi2', 'inverse_normal', 'logit' ou 'modified', "
            "reçu {value!r}.",
        },
        "diag.fisher.bad_deterministics": {
            "en": "deterministics must be 'constant', 'trend' or 'drift', got {value!r} (the "
            "p-values are not available without a constant).",
            "fr": "deterministics doit valoir 'constant', 'trend' ou 'drift', reçu {value!r} (les "
            "p-values n'existent pas sans constante).",
        },
        "diag.fisher.bad_lags": {
            "en": "lags must be a non-negative integer, got {value!r}.",
            "fr": "lags doit être un entier positif ou nul, reçu {value!r}.",
        },
        "diag.fisher.none_computed": {
            "en": "no unit has enough observations for an ADF regression with {lags} lag(s).",
            "fr": "aucune unité n'a assez d'observations pour une régression ADF à {lags} "
            "retard(s).",
        },
    },
)

#: ``p`` is floored here before taking logarithms (Stata's ``c(epsdouble)``).
_EPSILON = 2.220446049250313e-16

_METHODS = {
    "inverse_chi2": "P",
    "inverse_normal": "Z",
    "logit": "L",
    "modified": "Pm",
}


@register_test(
    "fisher",
    TestSpec(
        family="unitroot",
        category="A",
        wave=1,
        reference="Maddala, G. S. and Wu, S. (1999), Oxford Bull. Econ. Stat. 61(S1), 631-652; "
        "Choi, I. (2001), J. of International Money and Finance 20(2), 249-272",
        stata="xtunitroot fisher",
        validation=(
            "fi_lags1",
            "fi_trend_demean",
            "fi_drift",
            "fi_entry_exit",
            "fi_gaps",
            "fi_short_units",
        ),
        inputs="panel",
        balanced="both",
        r="plm::purtest(test='madwu')",
        summary="Combination of the p-values of one ADF test per unit.",
    ),
)
def fisher(
    panel: PanelBase,
    variable: str,
    *,
    lags: int,
    deterministics: Literal["constant", "trend", "drift"] = "constant",
    demean: bool = False,
    method: Literal["inverse_chi2", "inverse_normal", "logit", "modified"] = "inverse_chi2",
) -> TestResult:
    """Return a Fisher-type combination of the per-unit ADF p-values.

    One augmented Dickey-Fuller test with ``lags`` lagged differences is run per unit; the
    p-values ``p_i`` of the ``N`` units that could be tested are combined by ``method``:

    * ``"inverse_chi2"`` (Maddala-Wu, ``P``): ``-2 sum ln p_i`` on ``2N`` df, which needs a
      finite number of units;
    * ``"inverse_normal"`` (``Z``): ``sum Phi^{-1}(p_i) / sqrt(N)``, standard normal;
    * ``"logit"`` (``L``): ``sqrt(k) sum ln(p_i / (1 - p_i))`` on ``5N + 4`` df;
    * ``"modified"`` (``Pm``): ``sum (-2 ln p_i - 2) / (2 sqrt(N))``, standard normal, right tail.

    The null is a unit root in every unit.  The panel may be unbalanced and have holes: each
    unit is tested on its own calendar.  A unit whose regression cannot be run is left out and
    reported in a note (Stata's ``xtunitroot fisher`` does the same).  Inspect
    ``result.details`` for the per-unit statistic, p-value and sample.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of the numeric variable tested.
    lags
        Number of lagged differences in every ADF regression (Stata requires it).
    deterministics
        ``"constant"`` (default), ``"trend"`` (constant and linear trend, MacKinnon p-values) or
        ``"drift"`` (constant; p-values from Student's t, as ``dfuller, drift``).
    demean
        Remove the mean across units in each period first.
    method
        How to combine the p-values.

    Raises
    ------
    SpecificationError
        If ``lags``, ``deterministics`` or ``method`` is invalid.
    NumericalError
        If no unit can be tested.
    """
    # Stata: xtunitroot fisher, dfuller lags(#) [trend|drift] [demean] -> dfuller on every unit
    # (r(p)), p floored at epsdouble; P = -2 sum ln p (df 2N); Z = sum invnormal(p) / sqrt(N);
    # L = sqrt(3 (5N + 3) / (pi^2 N (5N + 2))) sum logit(p) (t, 5N + 4 df); Pm = sum(-2 ln p - 2)
    # / (2 sqrt(N)) (right tail).  Each statistic counts the units whose term is not missing: a
    # unit with p = 1 has no inverse normal or logit.
    # STATA-QUIRK: a unit with equations but no defined statistic (saturated regression) gets
    # p = 1 and counts in the degrees of freedom; Z and L drop it (invnormal(1) is missing).  A
    # unit without any equation is skipped, but Stata stops when that unit is the first one.
    if isinstance(lags, bool) or not isinstance(lags, (int, np.integer)) or lags < 0:
        raise SpecificationError(
            _i18n.message("diag.fisher.bad_lags", value=lags), key="diag.fisher.bad_lags"
        )
    if deterministics not in ("constant", "trend", "drift"):
        raise SpecificationError(
            _i18n.message("diag.fisher.bad_deterministics", value=deterministics),
            key="diag.fisher.bad_deterministics",
        )
    if method not in _METHODS:
        raise SpecificationError(
            _i18n.message("diag.fisher.bad_method", value=method), key="diag.fisher.bad_method"
        )
    arrays = panel_variables(panel, [variable], test="fisher", min_units=2)
    values = arrays.values[:, 0].copy()
    if demean:
        period = arrays.time_index
        totals = np.bincount(period, weights=values)
        counts = np.bincount(period)
        values -= (totals / np.where(counts > 0, counts, 1))[period]
    case: AdfCase = deterministics
    labels, sizes, taus, pvalues = [], [], [], []
    skipped = 0
    for unit, window in enumerate(arrays.unit_slices()):
        time = arrays.time_index[window]
        series = np.full(int(time[-1] - time[0]) + 1, np.nan)
        series[time - time[0]] = values[window]
        fit = adf_fit(series, int(lags), case)
        if fit is None:
            skipped += 1
            continue
        p_unit = adf_pvalue(fit, case)
        if math.isnan(p_unit):  # a drift test with an undefined statistic has no p-value
            skipped += 1
            continue
        labels.append(arrays.unit_labels[unit])
        sizes.append(fit)
        taus.append(fit.tau)
        pvalues.append(p_unit)
    if not pvalues:
        raise NumericalError(
            _i18n.message("diag.fisher.none_computed", lags=int(lags)),
            key="diag.fisher.none_computed",
        )
    p = np.maximum(np.array(pvalues), _EPSILON)
    statistic, pvalue, df = _combine(method, p)
    notes = (Note("diag.fisher.skipped", {"count": skipped, "lags": int(lags)}),) if skipped else ()
    details = _unit_table(labels, sizes, np.array(taus), np.array(pvalues))
    return TestResult(
        name="fisher",
        statistic=statistic,
        pvalue=pvalue,
        statistic_name=_METHODS[method],
        df=df,
        n_obs=arrays.n_obs,
        n_groups=arrays.n_units,
        variable=variable,
        options={
            "lags": int(lags),
            "deterministics": deterministics,
            "demean": demean,
            "method": method,
        },
        details=details,
        notes=notes,
    )


def _combine(method: str, p: np.ndarray) -> tuple[float, float, tuple[float, ...]]:
    """Return the statistic, p-value and degrees of freedom of the combination ``method``."""
    n = p.size
    if method == "inverse_chi2":
        statistic = float(-2.0 * np.sum(np.log(p)))
        return statistic, float(stats.chi2.sf(statistic, 2 * n)), (float(2 * n),)
    if method == "modified":
        statistic = float(np.sum(-2.0 * np.log(p) - 2.0) / (2.0 * math.sqrt(n)))
        return statistic, float(stats.norm.sf(statistic)), ()
    usable = p[p < 1.0]  # the inverse normal and the logit of p = 1 are missing
    m = usable.size
    if m == 0:
        return math.nan, math.nan, ()
    if method == "inverse_normal":
        statistic = float(np.sum(stats.norm.ppf(usable)) / math.sqrt(m))
        return statistic, float(stats.norm.cdf(statistic)), ()
    k = 3.0 * (5.0 * m + 3.0) / (math.pi**2 * m * (5.0 * m + 2.0))
    statistic = float(math.sqrt(k) * np.sum(np.log(usable / (1.0 - usable))))
    df = 5.0 * m + 4.0
    return statistic, float(stats.t.cdf(statistic, df)), (df,)


def _unit_table(
    labels: list[object], fits: list[AdfFit], taus: np.ndarray, pvalues: np.ndarray
) -> pd.DataFrame:
    """Return the per-unit table of the ADF regressions."""
    return pd.DataFrame(
        {
            "n": [fit.n for fit in fits],
            "lags": [fit.lags for fit in fits],
            "statistic": taus,
            "pvalue": pvalues,
        },
        index=pd.Index(labels, name="unit"),
    )
