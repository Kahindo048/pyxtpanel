"""Breitung unit-root test on a strongly balanced panel (category A, L.2.6).

Validated against Stata 17: ``xtunitroot breitung`` with a constant, without constant and with
panel trends, with and without prewhitening lags, with cross-sectional means removed and with the
statistic robust to cross-sectional correlation (Breitung and Das 2005).  The cases are in
``validation/cases/diagnostics_unitroot_breitung.py``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, panel_variables, register_test

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.breitung.title": {
            "en": "Breitung unit-root test",
            "fr": "Test de racine unitaire de Breitung",
        },
        "diag.breitung.null": {
            "en": "panels contain unit roots",
            "fr": "les séries du panel ont une racine unitaire",
        },
        "diag.breitung.alternative": {
            "en": "panels are stationary (common AR parameter)",
            "fr": "les séries du panel sont stationnaires (paramètre AR commun)",
        },
        "diag.breitung.bad_deterministics": {
            "en": "deterministics must be 'none', 'constant' or 'trend', got {value!r}.",
            "fr": "deterministics doit valoir 'none', 'constant' ou 'trend', reçu {value!r}.",
        },
        "diag.breitung.bad_lags": {
            "en": "lags must be a non-negative integer, got {value!r}.",
            "fr": "lags doit être un entier positif ou nul, reçu {value!r}.",
        },
        "diag.breitung.robust_needs_periods": {
            "en": "the robust statistic needs more periods, net of lags, than panels: "
            "T - lags - 1 = {available} is below N = {units}.",
            "fr": "la statistique robuste exige plus de périodes, retards déduits, que de séries : "
            "T - lags - 1 = {available} est inférieur à N = {units}.",
        },
        "diag.breitung.degenerate": {
            "en": "a unit has no residual variation after removing the lagged differences (or "
            "the transformed lagged level is zero everywhere): the statistic is not defined.",
            "fr": "une unité n'a plus de variation résiduelle après retrait des différences "
            "retardées (ou le niveau retardé transformé est nul partout) : la statistique n'est "
            "pas définie.",
        },
    },
)

_MODELS = ("none", "constant", "trend")


@register_test(
    "breitung",
    TestSpec(
        family="unitroot",
        category="A",
        wave=2,
        reference="Breitung, J. (2000), Advances in Econometrics 15, 161-178; Breitung, J. and "
        "Das, S. (2005), Statistica Neerlandica 59(4), 414-433",
        stata="xtunitroot breitung",
        validation=(
            "br_constant",
            "br_none",
            "br_trend",
            "br_lags2",
            "br_trend_lags",
            "br_none_lags",
            "br_demean",
            "br_robust",
            "br_robust_trend_lags",
        ),
        inputs="panel",
        balanced="required",
        summary="Unit-root test of a common AR parameter on forward-orthogonalised data.",
    ),
)
def breitung(
    panel: PanelBase,
    variable: str,
    *,
    deterministics: Literal["none", "constant", "trend"] = "constant",
    lags: int = 0,
    demean: bool = False,
    robust: bool = False,
) -> TestResult:
    """Return the Breitung ``lambda`` statistic, asymptotically standard normal.

    Each unit is purged of ``lags`` lagged differences (prewhitening); the lagged level is
    transformed to remove the deterministic terms (the first value for a constant, the
    first value and the trend for a trend) and the difference is, with a trend, put in forward
    orthogonal form (Helmert); the statistic is the pooled t ratio of their regression, with each
    unit scaled by its residual variance.  The test rejects for small values (left tail).  With
    ``robust`` the variance is replaced by an estimate that allows correlation across units
    (Breitung and Das 2005), which needs more periods than units.  The panel must be strongly
    balanced and without gaps.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of the numeric variable tested.
    deterministics
        ``"constant"`` (default), ``"none"`` or ``"trend"`` (panel means and linear trends).
    lags
        Number of lagged differences used to prewhiten each unit.
    demean
        Remove the mean across units in each period first.
    robust
        Report the statistic robust to cross-sectional correlation.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If ``deterministics`` or ``lags`` is invalid.
    DataError
        If the panel is not strongly balanced, has a gap, is too short, or ``robust`` is asked
        with no more periods (net of lags) than units.
    NumericalError
        If a unit has no residual variation.
    """
    # Stata: xtunitroot breitung -> per unit: constant: lagy = L.y - y[lags+1], dify = D.y, both
    # purged of L(1/lags)D.y (no constant); trend: du, u from the regression of D.y on L(1/lags)D.y
    # and a constant, lagy = u - u[lags+2] - (n - lags - 2) mean(du), dify = Helmert(du);
    # sig2 = sum(dify_or_du^2)/(T - lags - 2); lambda = sum(lagy dify/sig2) / sqrt(sum(lagy^2 /
    # sig2));
    # robust: lambda* = sum_t v1_t' dv_t / sqrt(sum_t v1_t' Omega v1_t), Omega[i, j] =
    # sum_t u_i u_j / (T - lags - 2) from u = dify - phi lagy, phi = sum s10 / sum s11.
    if deterministics not in _MODELS:
        raise SpecificationError(
            _i18n.message("diag.breitung.bad_deterministics", value=deterministics),
            key="diag.breitung.bad_deterministics",
        )
    if isinstance(lags, bool) or not isinstance(lags, (int, np.integer)) or lags < 0:
        raise SpecificationError(
            _i18n.message("diag.breitung.bad_lags", value=lags), key="diag.breitung.bad_lags"
        )
    lags = int(lags)
    trend = deterministics == "trend"
    arrays = panel_variables(
        panel,
        [variable],
        test="breitung",
        balanced="required",
        consecutive=True,
        min_obs=max(lags + 4, 2 * lags + 3) if trend else max(lags + 3, 2 * lags + 2),
        min_units=2,
    )
    y = arrays.as_cube()[:, :, 0].copy()
    n_units, periods = y.shape
    if robust and periods - lags - 1 < n_units:
        raise DataError(
            _i18n.message(
                "diag.breitung.robust_needs_periods", available=periods - lags - 1, units=n_units
            ),
            key="diag.breitung.robust_needs_periods",
        )
    if demean:
        y -= y.mean(axis=0)
    lagged = np.empty((n_units, periods))
    differenced = np.empty((n_units, periods))
    sigma2 = np.empty(n_units)
    for unit in range(n_units):
        lagged[unit], differenced[unit], sigma2[unit] = _unit(y[unit], lags, deterministics)
    if not np.all(np.isfinite(sigma2)) or np.any(sigma2 <= 0.0):
        raise NumericalError(
            _i18n.message("diag.breitung.degenerate"), key="diag.breitung.degenerate"
        )
    both = np.isfinite(lagged) & np.isfinite(differenced)
    product = np.where(both, lagged * differenced, 0.0)
    square = np.where(both, lagged**2, 0.0)
    sum_10 = float(np.sum(product.sum(axis=1) / sigma2))
    sum_11 = float(np.sum(square.sum(axis=1) / sigma2))
    if sum_11 <= 0.0:
        raise NumericalError(
            _i18n.message("diag.breitung.degenerate"), key="diag.breitung.degenerate"
        )
    lam = sum_10 / math.sqrt(sum_11)
    quantities = {"lambda": lam, "pvalue_lambda": float(stats.norm.cdf(lam))}
    statistic = lam
    if robust:
        statistic = _robust(lagged, differenced, sum_10 / sum_11, periods, lags, trend)
        quantities["lambda_robust"] = statistic
    details = pd.DataFrame(
        {"value": list(quantities.values())}, index=pd.Index(list(quantities), name="quantity")
    )
    return TestResult(
        name="breitung",
        statistic=float(statistic),
        pvalue=float(stats.norm.cdf(statistic)),
        statistic_name="lambda*" if robust else "lambda",
        df=(),
        n_obs=n_units * periods,
        n_groups=n_units,
        variable=variable,
        options={
            "deterministics": deterministics,
            "lags": lags,
            "demean": demean,
            "robust": robust,
        },
        details=details,
    )


def _residuals(target: FloatArray, regressors: FloatArray) -> FloatArray:
    """Return the least-squares residuals of ``target`` on the columns of ``regressors``."""
    coefficients, *_ = np.linalg.lstsq(regressors, target, rcond=None)
    return np.asarray(target - regressors @ coefficients)


def _unit(
    series: FloatArray, lags: int, deterministics: str
) -> tuple[FloatArray, FloatArray, float]:
    """Return the transformed lagged level, the transformed difference and ``sig2`` of one unit.

    Both series have ``NaN`` where they are not defined.
    """
    periods = series.size
    diff = np.full(periods, np.nan)
    diff[1:] = np.diff(series)
    level = np.full(periods, np.nan)
    level[1:] = series[:-1]
    rows = np.arange(lags + 1, periods)  # the equations: t = lags + 1, ..., T - 1 (0-based)
    lag_columns = [diff[rows - j] for j in range(1, lags + 1)]
    if deterministics != "trend":
        if deterministics == "constant":
            level = level - series[lags]
        out_dif, out_lag = np.full(periods, np.nan), np.full(periods, np.nan)
        if lags > 0:
            x = np.column_stack(lag_columns)
            out_dif[rows] = _residuals(diff[rows], x)
            out_lag[rows] = _residuals(level[rows], x)
        else:
            out_dif[1:], out_lag[1:] = diff[1:], level[1:]
        sig2 = float(np.nansum(out_dif**2)) / (periods - lags - 2)
        return out_lag, out_dif, sig2
    du = np.full(periods, np.nan)
    u = np.full(periods, np.nan)
    if lags > 0:
        x = np.column_stack([*lag_columns, np.ones(rows.size)])
        beta, *_ = np.linalg.lstsq(x, diff[rows], rcond=None)
        du[rows] = diff[rows] - sum(beta[j - 1] * diff[rows - j] for j in range(1, lags + 1))
        utmp = np.full(periods, np.nan)
        at = np.arange(lags, periods)
        utmp[at] = series[at] - sum(beta[j - 1] * series[at - j] for j in range(1, lags + 1))
        u[1:] = utmp[:-1]
    else:
        du[1:], u[1:] = diff[1:], level[1:]
    count = periods - lags - 1
    mean_du = float(np.nansum(du)) / count
    sig2 = float(np.nansum((du - mean_du) * du)) / (periods - lags - 2)
    first = lags + 1
    out_lag = np.full(periods, np.nan)
    out_lag[first:] = u[first:] - u[first] - (np.arange(first, periods) - lags - 1) * mean_du
    cumulative = np.cumsum(np.nan_to_num(du))
    t = np.arange(periods)
    out_dif = np.full(periods, np.nan)
    ahead = periods - 1 - t  # number of observations after t
    valid = (ahead > 0) & np.isfinite(du)
    forward = (cumulative[-1] - cumulative[valid]) / ahead[valid]
    out_dif[valid] = (du[valid] - forward) * np.sqrt(ahead[valid] / (ahead[valid] + 1.0))
    return out_lag, out_dif, sig2


def _robust(
    lagged: FloatArray,
    differenced: FloatArray,
    phi: float,
    periods: int,
    lags: int,
    trend: bool,
) -> float:
    """Return the statistic robust to cross-sectional correlation (Breitung and Das 2005).

    Stata sums over the periods ``lags + 1`` to ``T - 1`` (``lags`` to ``T - 2`` with a trend,
    where the last forward-orthogonal difference does not exist), a missing value counting as 0.
    """
    residual = np.where(
        np.isfinite(lagged) & np.isfinite(differenced), differenced - phi * lagged, 0.0
    )
    omega = residual @ residual.T / (periods - lags - 2)
    window = slice(lags, periods - 1) if trend else slice(lags + 1, periods)
    v1 = np.nan_to_num(lagged)[:, window]
    dv = np.nan_to_num(differenced)[:, window]
    a = float(np.sum(v1 * dv))
    b = float(np.sum((omega @ v1) * v1))
    if b <= 0.0:
        raise NumericalError(
            _i18n.message("diag.breitung.degenerate"), key="diag.breitung.degenerate"
        )
    return a / math.sqrt(b)
