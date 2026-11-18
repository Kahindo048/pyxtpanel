"""Kao residual-based panel cointegration tests (category A, L.4.1).

Validated against Stata 17: ``xtcointtest kao`` with the default Newey-West lag selection and
with fixed kernel lags, with the Bartlett and Parzen kernels, with fixed and selected augmented
lags, with cross-sectional means removed, with one and two regressors, on a cointegrated panel
and on a panel without cointegration.  The cases are in
``validation/cases/diagnostics_cointegration_kao.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.long_run import long_run_covariance
from pyxtpanel.diagnostics._shared.ols import is_perfect_fit

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.kao.title": {
            "en": "Kao test for cointegration",
            "fr": "Test de cointégration de Kao",
        },
        "diag.kao.null": {"en": "no cointegration", "fr": "pas de cointégration"},
        "diag.kao.alternative": {
            "en": "all panels are cointegrated (same cointegrating vector)",
            "fr": "toutes les séries du panel sont cointégrées (même vecteur de cointégration)",
        },
        "diag.kao.no_regressor": {
            "en": "kao needs at least one regressor besides the dependent variable.",
            "fr": "kao exige au moins un régresseur en plus de la variable dépendante.",
        },
        "diag.kao.bad_option": {
            "en": "{name} must be one of {allowed}, got {value!r}.",
            "fr": "{name} doit valoir l'une de ces valeurs : {allowed} ; reçu {value!r}.",
        },
        "diag.kao.bad_lags": {
            "en": "{name} must be a non-negative integer, got {value!r}.",
            "fr": "{name} doit être un entier positif ou nul, reçu {value!r}.",
        },
        "diag.kao.degenerate": {
            "en": "the residuals have no usable variation (a zero variance or a singular "
            "covariance of the regressors' innovations): the statistics are not defined.",
            "fr": "les résidus n'ont pas de variation exploitable (variance nulle ou covariance "
            "singulière des innovations des régresseurs) : les statistiques ne sont pas définies.",
        },
        "diag.kao.undefined_bandwidth": {
            "en": "the Newey-West lag selection is not defined for a unit (its long-run variance "
            "estimate is not positive); give a fixed number of kernel lags.",
            "fr": "la sélection de retards de Newey-West n'est pas définie pour une unité (son "
            "estimation de variance de long terme n'est pas positive) ; donnez un nombre fixe de "
            "retards du noyau.",
        },
        "diag.kao.quirk": {
            "en": "p-values are the left tail of the standard normal; Stata's xtcointtest reports "
            "1 - Phi(|z|), which is the same for a negative statistic and differs for a positive "
            "one.",
            "fr": "les p-values sont la queue gauche de la loi normale ; xtcointtest de Stata "
            "affiche 1 - Phi(|z|), identique pour une statistique négative et différente pour une "
            "positive.",
        },
    },
)

_STATISTICS = ("modified_df", "df", "adf", "unadjusted_modified_df", "unadjusted_df")
_KERNELS = ("bartlett", "parzen", "quadraticspectral")


@register_test(
    "kao",
    TestSpec(
        family="cointegration",
        category="A",
        wave=2,
        reference="Kao, C. (1999), J. of Econometrics 90(1), 1-44",
        stata="xtcointtest kao",
        validation=(
            "kao_default",
            "kao_default_t",
            "kao_nocoint",
            "kao_nocoint_t",
            "kao_bartlett3",
            "kao_bartlett3_t",
            "kao_parzen",
            "kao_parzen_t",
            "kao_parzen_auto",
            "kao_parzen_auto_t",
            "kao_adf2",
            "kao_adf2_t",
            "kao_adf_aic",
            "kao_adf_aic_t",
            "kao_adf_bic",
            "kao_adf_bic_t",
            "kao_adf_hqic",
            "kao_adf_hqic_t",
            "kao_demean",
            "kao_demean_t",
            "kao_one_regressor",
            "kao_one_regressor_t",
            "kao_qs",
            "kao_qs_t",
            "kao_qs_auto",
            "kao_qs_auto_t",
        ),
        inputs="panel",
        balanced="required",
        summary="Residual-based tests of no cointegration with a common cointegrating vector.",
    ),
)
def kao(
    panel: PanelBase,
    dependent: str,
    regressors: Sequence[str],
    *,
    statistic: Literal[
        "modified_df", "df", "adf", "unadjusted_modified_df", "unadjusted_df"
    ] = "adf",
    adf_lags: int = 1,
    adf_lag_selection: Literal["aic", "bic", "hqic"] | None = None,
    kernel: Literal["bartlett", "parzen", "quadraticspectral"] = "bartlett",
    kernel_lags: int | None = None,
    demean: bool = False,
) -> TestResult:
    """Return the Kao statistics of no cointegration, asymptotically standard normal.

    The model is ``y_it = a_i + x_it' b + e_it`` with one cointegrating vector ``b`` for all
    units.  It is fitted by within least squares, and ``e`` is tested for a unit root by a pooled
    Dickey-Fuller regression (common AR parameter ``rho``) and a pooled augmented one.  Five
    statistics are computed, all rejecting for small values:

    * ``"modified_df"`` and ``"df"`` correct the Dickey-Fuller statistics for the endogeneity of
      the regressors and for serial correlation with the long-run covariance of the innovations;
    * ``"adf"`` is the corrected statistic of the augmented regression;
    * ``"unadjusted_modified_df"`` and ``"unadjusted_df"`` assume strictly exogenous regressors
      and no serial correlation.

    ``statistic`` chooses the one reported as the result's statistic and p-value; all five, with
    their p-values, are in ``details``.  The panel must be strongly balanced and without gaps.

    Parameters
    ----------
    panel
        The declared panel.
    dependent
        Name of the dependent variable.
    regressors
        Names of the regressors (at least one), without a constant.
    statistic
        The statistic reported as ``result.statistic`` (default the augmented one).
    adf_lags
        Number of lagged differences of the augmented regression, or their maximum when
        ``adf_lag_selection`` is given (between 1 and ``adf_lags``).
    adf_lag_selection
        ``"aic"``, ``"bic"`` or ``"hqic"`` to choose the number of augmented lags.
    kernel
        Kernel of the long-run covariance.
    kernel_lags
        Number of lags of the kernel, or ``None`` (default) for the Newey-West selection of each
        unit.
    demean
        Remove the mean across units in each period first.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If an option is invalid or there is no regressor.
    DataError
        If the panel is not strongly balanced, has a gap or is too short.
    NumericalError
        If the residual covariances are singular, or the automatic selection is not defined.
    """
    # Stata: xtcointtest kao -> err_v = residual of the pooled `reg v L.v, noconstant` (v = y, x),
    # per-unit centred; Sigma_i, Omega_i (Newey-West, kernel lags per unit, bw = int(m), lags =
    # min(bw, int(4 (n/100)^(2/9))) from the trace of the autocovariances) averaged over units;
    # s2 = S_uu - S_ue S_ee^-1 S_eu, s02 likewise from Omega; e = within residual of y on x over
    # t >= 2; rho from the pooled DF regression of e on its lag over t >= 3 (T = periods left);
    # DF*_rho = (sqrt(N) T (rho - 1) + 3 sqrt(N) s2/s02) / sqrt(3 + 36 s2^2/(5 s02^2));
    # DF*_t = (t + sqrt(6 N s2/s02)/2) / sqrt(s02/(2 s2) + 3 s2/(10 s02)); ADF likewise from the
    # pooled augmented t; DF_rho = (sqrt(N) T (rho - 1) + 3 sqrt(N)) / sqrt(51/5);
    # DF_t = sqrt(5/4) t + sqrt(15 N/8); p = 1 - normal(|z|).
    names = list(regressors)
    if not names:
        raise SpecificationError(
            _i18n.message("diag.kao.no_regressor"), key="diag.kao.no_regressor"
        )
    _validate(statistic, kernel, kernel_lags, adf_lags, adf_lag_selection)
    arrays = panel_variables(
        panel,
        [dependent, *names],
        test="kao",
        balanced="required",
        consecutive=True,
        min_obs=adf_lags + len(names) + 5,
        min_units=2,
    )
    z = arrays.as_cube().copy()  # (units, periods, series)
    n_units, periods, _ = z.shape
    if demean:
        z -= z.mean(axis=0)
    long_run = long_run_covariance(_innovations(z), kernel, kernel_lags)
    if long_run is None:
        raise NumericalError(
            _i18n.message("diag.kao.undefined_bandwidth"), key="diag.kao.undefined_bandwidth"
        )
    ratio = _variance_ratio(long_run.sigma, long_run.omega)
    e = _within_residuals(z[:, 1:, :])
    lag, cur = e[:, :-1], e[:, 1:]
    sum_lag2 = float(np.sum(lag**2))
    rho = float(np.sum(lag * cur)) / sum_lag2
    t_periods = periods - 2
    resid = cur - rho * lag
    s2_df = float(np.sum(resid**2)) / (n_units * t_periods)
    t_rho = (rho - 1.0) * math.sqrt(sum_lag2 / s2_df)
    p_adf = _selected_lags(e, adf_lags, adf_lag_selection)
    t_adf = _adf_t(e, p_adf, t_periods)
    root_n = math.sqrt(n_units)
    r = ratio
    denominator_t = math.sqrt(1.0 / (2.0 * r) + 3.0 * r / 10.0)
    values = {
        "modified_df": (root_n * t_periods * (rho - 1.0) + 3.0 * root_n * r)
        / math.sqrt(3.0 + 36.0 * r**2 / 5.0),
        "df": (t_rho + math.sqrt(6.0 * n_units * r) / 2.0) / denominator_t,
        "adf": (t_adf + math.sqrt(6.0 * n_units * r) / 2.0) / denominator_t,
        "unadjusted_modified_df": (root_n * t_periods * (rho - 1.0) + 3.0 * root_n)
        / math.sqrt(51.0 / 5.0),
        "unadjusted_df": math.sqrt(5.0 / 4.0) * t_rho + math.sqrt(15.0 * n_units / 8.0),
    }
    details = pd.DataFrame(
        {
            "statistic": [values[k] for k in _STATISTICS],
            "pvalue": [float(stats.norm.cdf(values[k])) for k in _STATISTICS],
        },
        index=pd.Index(list(_STATISTICS), name="statistic"),
    )
    chosen = values[statistic]
    return TestResult(
        name="kao",
        statistic=float(chosen),
        pvalue=float(stats.norm.cdf(chosen)),
        statistic_name=statistic,
        df=(),
        n_obs=n_units * t_periods,
        n_groups=n_units,
        variable=dependent,
        options={
            "regressors": tuple(names),
            "adf_lags": p_adf,
            "kernel": kernel,
            "kernel_lags": long_run.lags,
            "demean": demean,
        },
        details=details,
    )


def _validate(
    statistic: str, kernel: str, kernel_lags: object, adf_lags: object, selection: object
) -> None:
    """Refuse an invalid option with the message of the catalogue."""
    for name, value, allowed in (
        ("statistic", statistic, _STATISTICS),
        ("kernel", kernel, _KERNELS),
        ("adf_lag_selection", selection, (None, "aic", "bic", "hqic")),
    ):
        if value not in allowed:
            raise SpecificationError(
                _i18n.message("diag.kao.bad_option", name=name, allowed=list(allowed), value=value),
                key="diag.kao.bad_option",
            )
    for name, value, minimum in (("adf_lags", adf_lags, 1), ("kernel_lags", kernel_lags, 0)):
        if name == "kernel_lags" and value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
            raise SpecificationError(
                _i18n.message("diag.kao.bad_lags", name=name, value=value), key="diag.kao.bad_lags"
            )


def _innovations(z: FloatArray) -> FloatArray:
    """Residuals of the pooled regression ``v_t = a v_{t-1}`` (no constant) of every series."""
    lag, cur = z[:, :-1, :], z[:, 1:, :]
    slope = np.sum(lag * cur, axis=(0, 1)) / np.sum(lag * lag, axis=(0, 1))
    return np.asarray(cur - slope * lag)


def _variance_ratio(sigma: FloatArray, omega: FloatArray) -> float:
    """Return ``s2 / s02``: conditional contemporaneous over conditional long-run variance."""

    def conditional(m: FloatArray) -> float:
        try:
            return float(m[0, 0] - m[0, 1:] @ np.linalg.solve(m[1:, 1:], m[1:, 0]))
        except np.linalg.LinAlgError:
            return math.nan

    s2, s02 = conditional(sigma), conditional(omega)
    if not (math.isfinite(s2) and math.isfinite(s02) and s2 > 0.0 and s02 > 0.0):
        raise NumericalError(_i18n.message("diag.kao.degenerate"), key="diag.kao.degenerate")
    return s2 / s02


def _within_residuals(z: FloatArray) -> FloatArray:
    """Residuals of the within regression of series 0 on the others, ``(units, periods)``."""
    centred = z - z.mean(axis=1, keepdims=True)
    y, x = centred[:, :, 0].ravel(), centred[:, :, 1:].reshape(-1, z.shape[2] - 1)
    slope, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ slope
    # a perfect fit (x collinear with y) leaves rounding noise, which no statistic should read
    if is_perfect_fit(float(residual @ residual), float(y @ y)):
        raise NumericalError(_i18n.message("diag.kao.degenerate"), key="diag.kao.degenerate")
    return np.asarray(residual.reshape(z.shape[0], z.shape[1]))


def _adf_design(e: FloatArray, p: int, first: int) -> tuple[FloatArray, FloatArray]:
    """Return the pooled response and regressors (lagged level, ``p`` lagged differences)."""
    diff = np.diff(e, axis=1)  # diff[:, s] = e[:, s + 1] - e[:, s]
    rows = np.arange(first, e.shape[1])
    columns = [e[:, rows - 1]]
    columns.extend(diff[:, rows - 1 - j] for j in range(1, p + 1))
    x = np.stack(columns, axis=2).reshape(-1, p + 1)
    return e[:, rows].ravel(), x


def _adf_t(e: FloatArray, p: int, t_periods: int) -> float:
    """Return Stata's ``t_ADF``: the t ratio of the lagged level in the pooled ADF regression.

    The coefficient and the residuals are those of the pooled regression; the variance of the
    residuals is divided by ``N T`` (``T`` the periods of the Dickey-Fuller regression) and the
    lagged level enters the standard error after each unit's own projection on its lagged
    differences, which differs from the ordinary least-squares standard error when ``N > 1``.
    """
    n_units = e.shape[0]
    y, x = _adf_design(e, p, p + 1)
    coefficients, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ coefficients
    s2 = float(residual @ residual) / (n_units * t_periods)
    level = x[:, 0].reshape(n_units, -1)
    lagged = x[:, 1:].reshape(n_units, level.shape[1], p)
    squared = 0.0
    for unit in range(n_units):
        z, v = lagged[unit], level[unit]
        v_tilde = v - z @ np.linalg.lstsq(z, v, rcond=None)[0]
        squared += float(v_tilde @ v_tilde)
    return float((coefficients[0] - 1.0) / math.sqrt(s2 / squared))


def _selected_lags(e: FloatArray, maximum: int, criterion: str | None) -> int:
    """Return the number of augmented lags: ``maximum``, or the one minimising ``criterion``."""
    if criterion is None:
        return maximum
    best, best_value = maximum, math.inf
    for p in range(1, maximum + 1):
        y, x = _adf_design(e, p, maximum + 1)
        coefficients, *_ = np.linalg.lstsq(x, y, rcond=None)
        n = y.size
        rss = float(np.sum((y - x @ coefficients) ** 2))
        log_likelihood = -0.5 * n * (math.log(2.0 * math.pi) + 1.0 + math.log(rss / n))
        penalty = {"aic": 2.0, "bic": math.log(n), "hqic": 2.0 * math.log(math.log(n))}[criterion]
        value = -2.0 * log_likelihood + penalty * (p + 1)
        if value < best_value:
            best, best_value = p, value
    return best
