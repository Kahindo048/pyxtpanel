"""Levin-Lin-Chu unit-root test on a strongly balanced panel (category A, L.2.2).

Validated against Stata 17: ``xtunitroot llc`` with a constant, without constant and with panel
trends, with fixed lags and with lags chosen by information criteria, with the Bartlett, Parzen
and quadratic-spectral kernels and the three bandwidth rules (fixed, ``llc``, ``nwest``), with
cross-sectional means removed and without interpolation of the adjustment table.  The cases are
in ``validation/cases/diagnostics_unitroot_levin_lin_chu.py``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import (
    TestResult,
    TestSpec,
    panel_variables,
    plain_label,
    register_test,
)
from pyxtpanel.diagnostics._shared.adf import select_lags
from pyxtpanel.diagnostics._shared.ols import is_perfect_fit, ols_residuals
from pyxtpanel.diagnostics._shared.tables import load_table

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.levin_lin_chu.title": {
            "en": "Levin-Lin-Chu unit-root test",
            "fr": "Test de racine unitaire de Levin-Lin-Chu",
        },
        "diag.levin_lin_chu.null": {
            "en": "panels contain unit roots",
            "fr": "les séries du panel ont une racine unitaire",
        },
        "diag.levin_lin_chu.alternative": {
            "en": "panels are stationary (common AR parameter)",
            "fr": "les séries du panel sont stationnaires (paramètre AR commun)",
        },
        "diag.levin_lin_chu.bad_option": {
            "en": "{name} must be one of {allowed}, got {value!r}.",
            "fr": "{name} doit valoir l'une de ces valeurs : {allowed} ; reçu {value!r}.",
        },
        "diag.levin_lin_chu.bad_lags": {
            "en": "lags must be a non-negative integer, got {value!r}.",
            "fr": "lags doit être un entier positif ou nul, reçu {value!r}.",
        },
        "diag.levin_lin_chu.bad_bandwidth": {
            "en": "bandwidth must be 'llc', 'nwest' or a non-negative number (an integer for the "
            "Bartlett and Parzen kernels), got {value!r}.",
            "fr": "bandwidth doit valoir 'llc', 'nwest' ou un nombre positif ou nul (entier pour "
            "les noyaux de Bartlett et de Parzen), reçu {value!r}.",
        },
        "diag.levin_lin_chu.undefined_bandwidth": {
            "en": "the Newey-West bandwidth is not defined for a unit (its long-run variance "
            "estimate is not positive); give a fixed bandwidth.",
            "fr": "la fenêtre de Newey-West n'est pas définie pour une unité (l'estimation de sa "
            "variance de long terme n'est pas positive) ; donnez une fenêtre fixe.",
        },
        "diag.levin_lin_chu.constant_unit": {
            "en": "unit {unit!r} has no variation left once the ADF terms are removed (it is "
            "constant, or its augmented regression fits perfectly): the standardised residuals "
            "are not defined.",
            "fr": "l'unité {unit!r} ne varie plus une fois les termes ADF retirés (elle est "
            "constante, ou sa régression augmentée est parfaite) : les résidus standardisés ne "
            "sont pas définis.",
        },
        "diag.levin_lin_chu.degenerate": {
            "en": "the pooled regression of the standardised residuals is degenerate (the lagged "
            "level has no variation left after the ADF terms).",
            "fr": "la régression empilée des résidus standardisés est dégénérée (le niveau retardé "
            "n'a plus de variation après les termes ADF).",
        },
    },
)

_MODELS = {"none": 1, "constant": 2, "trend": 3}
_KERNELS = ("bartlett", "parzen", "quadraticspectral")


@register_test(
    "levin_lin_chu",
    TestSpec(
        family="unitroot",
        category="A",
        wave=1,
        reference="Levin, A., Lin, C.-F. and Chu, C.-S. J. (2002), J. of Econometrics 108(1), 1-24",
        stata="xtunitroot llc",
        validation=(
            "llc_lags1",
            "llc_none",
            "llc_none_aic",
            "llc_trend",
            "llc_demean",
            "llc_aic",
            "llc_hqic_parzen",
            "llc_parzen_nwest",
            "llc_bartlett_nwest",
            "llc_qs_fixed",
            "llc_qs_nwest",
            "llc_nointerpolate",
        ),
        inputs="panel",
        aliases=("llc",),
        balanced="required",
        tables=("llc_adjustments.csv",),
        r="plm::purtest(test='levinlin')",
        summary="Unit-root test of a common AR parameter, with the bias-adjusted t statistic.",
    ),
)
def levin_lin_chu(
    panel: PanelBase,
    variable: str,
    *,
    deterministics: Literal["none", "constant", "trend"] = "constant",
    lags: int = 1,
    lag_selection: Literal["aic", "bic", "hqic"] | None = None,
    kernel: Literal["bartlett", "parzen", "quadraticspectral"] = "bartlett",
    bandwidth: int | float | Literal["llc", "nwest"] = "llc",
    demean: bool = False,
    interpolate: bool = True,
) -> TestResult:
    """Return the Levin-Lin-Chu adjusted t statistic ``t*``, asymptotically standard normal.

    For each unit the first difference and the lagged level are purged of ``lags`` lagged
    differences and of the deterministic terms; the residuals are standardised by the unit's
    residual variance and pooled, and ``delta`` is their common slope.  The unadjusted
    ``t_delta`` is corrected for the bias and variance that the deterministic terms and the
    estimated long-run variance cause (Levin, Lin and Chu 2002, Table 2), with the long-run
    variance of each unit estimated by ``kernel`` with ``bandwidth``.  The test rejects for small
    values (left tail).  The panel must be strongly balanced and without gaps.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of the numeric variable tested.
    deterministics
        ``"constant"`` (default), ``"none"`` or ``"trend"`` (unit means and linear trends).
    lags
        Number of lagged differences, or their maximum when ``lag_selection`` is given.
    lag_selection
        ``"aic"``, ``"bic"`` or ``"hqic"`` to choose the number of lags of each unit between 0
        and ``lags``.
    kernel
        Kernel of the long-run variance.
    bandwidth
        Number of lags of the kernel: a fixed number, ``"llc"`` (``int(3.21 T^(1/3))``) or
        ``"nwest"`` (Newey-West 1994 automatic selection).
    demean
        Remove the mean across units in each period first.
    interpolate
        Interpolate the adjustment table in the average number of equations (default); with
        ``False`` take the row of the first tabulated value at or above it.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If an option has an invalid value.
    DataError
        If the panel is not strongly balanced, has a gap or is too short.
    """
    # Stata: xtunitroot llc -> per unit: regress D.y and L.y on L(1/p).D.y and the deterministic
    # terms (_pregress), mrssdf = rss(e on nu) / (T - p - 1), standardise, s_i = sqrt(varhac /
    # mrssdf); pooled regression of e on nu without constant, se scaled by sqrt((N - 1) / N);
    # t* = (t - N~ sbar / sig_e2 * se * mu) / sigma, with (mu, sigma) from Table 2 at the average
    # number of equations.
    # STATA-QUIRK: with lag selection and no constant, Stata passes a mistyped macro to its lag
    # selection, so the information criteria are computed with a constant in the regressions.
    model = _validated(deterministics, lags, lag_selection, kernel, bandwidth)
    min_obs = lags + 4 + (1 if model == 3 else 0)
    arrays = panel_variables(
        panel,
        [variable],
        test="levin_lin_chu",
        balanced="required",
        consecutive=True,
        min_obs=min_obs,
        min_units=2,
    )
    y = arrays.as_cube()[:, :, 0].copy()
    n_units, periods = y.shape
    if demean:
        y -= y.mean(axis=0)
    e_all, nu_all, s_all, used, hac_used, equations = [], [], [], [], [], []
    for position, row in enumerate(y):
        try:
            e, nu, s_i, p_i, hac = _unit(row, model, lags, lag_selection, kernel, bandwidth)
        except _ConstantUnitError:
            raise NumericalError(
                _i18n.message(
                    "diag.levin_lin_chu.constant_unit",
                    unit=plain_label(arrays.unit_labels[position]),
                ),
                key="diag.levin_lin_chu.constant_unit",
            ) from None
        e_all.append(e)
        nu_all.append(nu)
        s_all.append(s_i)
        used.append(p_i)
        hac_used.append(hac)
        equations.append(e.size)
    e_pool, nu_pool = np.concatenate(e_all), np.concatenate(nu_all)
    sum_nu2 = float(nu_pool @ nu_pool)
    if sum_nu2 <= 0.0:
        raise NumericalError(
            _i18n.message("diag.levin_lin_chu.degenerate"), key="diag.levin_lin_chu.degenerate"
        )
    n = e_pool.size
    delta = float(nu_pool @ e_pool) / sum_nu2
    rss = float(np.sum((e_pool - delta * nu_pool) ** 2))
    se = math.sqrt(rss / n / sum_nu2)
    td = delta / se
    t_tilde = float(np.mean(equations))
    s_bar = float(np.mean(s_all))
    var_e = rss / n
    mu_adj, sigma_adj = _adjustments(t_tilde, model, interpolate)
    total = float(np.sum(equations))
    tds = (td - total * s_bar / var_e * se * mu_adj) / sigma_adj
    quantities = {
        "delta": delta,
        "se_delta": se,
        "t_delta": td,
        "mu_adj": mu_adj,
        "sigma_adj": sigma_adj,
        "t_tilde": t_tilde,
        "s_bar": s_bar,
        "var_e": var_e,
        "adf_lags": float(np.mean(used)),
        "hac_lags": float(np.mean(hac_used)),
    }
    if model == 1:
        quantities["p_t_delta"] = float(stats.norm.cdf(td))
    details = pd.DataFrame(
        {"value": list(quantities.values())}, index=pd.Index(list(quantities), name="quantity")
    )
    return TestResult(
        name="levin_lin_chu",
        statistic=float(tds),
        pvalue=float(stats.norm.cdf(tds)),
        statistic_name="t*",
        df=(),
        n_obs=n_units * periods,
        n_groups=n_units,
        variable=variable,
        options={
            "deterministics": deterministics,
            "lags": lags,
            "lag_selection": lag_selection,
            "kernel": kernel,
            "bandwidth": bandwidth,
            "demean": demean,
            "interpolate": interpolate,
        },
        details=details,
    )


def _validated(
    deterministics: object,
    lags: object,
    lag_selection: object,
    kernel: object,
    bandwidth: object,
) -> int:
    """Check the options and return Levin, Lin and Chu's model number (1, 2 or 3)."""
    for name, value, allowed in (
        ("deterministics", deterministics, tuple(_MODELS)),
        ("lag_selection", lag_selection, (None, "aic", "bic", "hqic")),
        ("kernel", kernel, _KERNELS),
    ):
        if value not in allowed:
            raise SpecificationError(
                _i18n.message(
                    "diag.levin_lin_chu.bad_option",
                    name=name,
                    allowed=", ".join(map(repr, allowed)),
                    value=value,
                ),
                key="diag.levin_lin_chu.bad_option",
            )
    if isinstance(lags, bool) or not isinstance(lags, (int, np.integer)) or lags < 0:
        raise SpecificationError(
            _i18n.message("diag.levin_lin_chu.bad_lags", value=lags),
            key="diag.levin_lin_chu.bad_lags",
        )
    valid: bool
    if isinstance(bandwidth, str):
        valid = bandwidth in ("llc", "nwest")
    elif isinstance(bandwidth, (int, float, np.integer, np.floating)) and not isinstance(
        bandwidth, bool
    ):
        valid = bool(bandwidth >= 0) and (
            kernel == "quadraticspectral" or float(bandwidth) == int(bandwidth)
        )
    else:
        valid = False
    if not valid:
        raise SpecificationError(
            _i18n.message("diag.levin_lin_chu.bad_bandwidth", value=bandwidth),
            key="diag.levin_lin_chu.bad_bandwidth",
        )
    return _MODELS[str(deterministics)]


class _ConstantUnitError(Exception):
    """A unit has no variation left after the ADF terms (turned into a refusal naming it)."""


def _unit(
    y: FloatArray,
    model: int,
    lags: int,
    lag_selection: str | None,
    kernel: str,
    bandwidth: int | float | str,
) -> tuple[FloatArray, FloatArray, float, int, float]:
    """Return a unit's standardised residuals, ``s_i``, lags used and kernel bandwidth used."""
    periods = y.size
    dy, level = np.diff(y), y[:-1]
    trend = np.arange(2.0, periods + 1.0)
    with_trend = model == 3
    if lag_selection is not None:
        p = select_lags(y, lags, lag_selection, trend=with_trend, constant=True)  # type: ignore[arg-type]
    else:
        p = lags
    columns = [dy[p - k : periods - 1 - k] for k in range(1, p + 1)]
    if with_trend:
        columns.append(trend[p:])
    if model >= 2:
        columns.append(np.ones(periods - 1 - p))
    z = np.column_stack(columns) if columns else np.empty((periods - 1 - p, 0))
    e = ols_residuals(dy[p:], z)
    nu = ols_residuals(level[p:], z)
    sum_nu2 = float(nu @ nu)
    # a constant series leaves rounding noise, not an exact zero: compare with the data's scale
    if is_perfect_fit(sum_nu2, float(level @ level)):
        raise _ConstantUnitError
    slope = float(nu @ e) / sum_nu2
    rss = float(np.sum((e - slope * nu) ** 2))
    if is_perfect_fit(rss, float(dy @ dy)):
        raise _ConstantUnitError
    sigma2 = rss / (periods - p - 1)
    if model == 1:
        w = dy
    elif model == 2:
        w = dy - dy.mean()
    else:
        w = ols_residuals(dy, np.column_stack([np.ones(periods - 1), trend]))
    if bandwidth == "llc":
        hac = float(int(3.21 * periods ** (1.0 / 3.0)))
    elif bandwidth == "nwest":
        hac = _newey_west_bandwidth(w, kernel)
    else:
        hac = float(bandwidth)
    long_run = _long_run_variance(w, kernel, hac)
    scale = math.sqrt(sigma2)
    return e / scale, nu / scale, math.sqrt(long_run / sigma2), p, hac


def _autocovariance_sum(w: FloatArray, lag: int) -> float:
    """Return ``sum_t w_t w_{t-lag}`` over the observations that have the lag."""
    return float(np.dot(w[lag:], w[:-lag])) if 0 < lag < w.size else 0.0


def _long_run_variance(w: FloatArray, kernel: str, lags: float) -> float:
    """Return the kernel estimate of the long-run variance of ``w`` (Stata's ``_gethacvar``)."""
    n = w.size
    total = float(w @ w)
    accumulated = 0.0
    if kernel == "bartlett":
        count = int(lags)
        for i in range(1, count + 1):
            accumulated += 2.0 * (1.0 - i / (count + 1.0)) * _autocovariance_sum(w, i)
    elif kernel == "parzen":
        count = int(lags)
        for i in range(1, count + 1):
            x = i / (count + 1.0)
            weight = 1 - 6 * x**2 + 6 * abs(x) ** 3 if x <= 0.5 else 2 * (1 - abs(x)) ** 3
            accumulated += 2.0 * weight * _autocovariance_sum(w, i)
    else:
        for i in range(1, n + 1):
            x = i / (lags + 1.0)
            argument = 6.0 * math.pi * x / 5.0
            weight = (
                25.0
                / (12.0 * (math.pi * x) ** 2)
                * (math.sin(argument) / argument - math.cos(argument))
            )
            accumulated += 2.0 * weight * _autocovariance_sum(w, i)
    return (total + accumulated) / n


def _newey_west_bandwidth(w: FloatArray, kernel: str) -> float:
    """Return Newey and West's (1994) automatic number of lags (Stata's ``nwest``)."""
    n = w.size
    if kernel == "bartlett":
        q, cy, integer, steps = 1, 1.1447, True, int(4 * (n / 100) ** (2 / 9))
    elif kernel == "parzen":
        q, cy, integer, steps = 2, 2.6614, True, int(4 * (n / 100) ** (4 / 25))
    else:
        q, cy, integer, steps = 2, 1.3221, False, int(4 * (n / 100) ** (2 / 25))
    s_q, s_0 = 0.0, float(w @ w)
    for j in range(1, steps + 1):
        sigma_j = _autocovariance_sum(w, j)
        s_q += 2.0 * j**q * sigma_j
        s_0 += 2.0 * sigma_j
    exponent = 1.0 / (2 * q + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = cy * ((s_q / s_0) ** 2) ** exponent if s_0 != 0 else math.nan
    if not math.isfinite(gamma):
        raise NumericalError(
            _i18n.message("diag.levin_lin_chu.undefined_bandwidth"),
            key="diag.levin_lin_chu.undefined_bandwidth",
        )
    if integer:
        band = int(gamma * n**exponent)
        return float(band - 1 if band < n else n - 2)
    band_real = gamma * n**exponent
    return band_real - 1.0 if band_real <= n else n - 1.0


def _adjustments(t: float, model: int, interpolate: bool) -> tuple[float, float]:
    """Return the mean and standard-deviation adjustments for ``t`` equations per unit."""
    table = load_table("llc_adjustments.csv")
    rows = table[table["model"] == model].reset_index(drop=True)
    keys = rows["t"].to_numpy()
    finite = int(np.sum(np.isfinite(keys)))
    position = int(np.searchsorted(keys[:finite], t, side="left"))  # first key >= t
    row = position if position < finite else finite  # past the last key: the T = infinity row
    if row in (0, finite) or not interpolate:
        return float(rows.loc[row, "mu"]), float(rows.loc[row, "sigma"])
    below, above = rows.loc[row - 1], rows.loc[row]
    weight = (t - below["t"]) / (above["t"] - below["t"])
    return (
        float((1 - weight) * below["mu"] + weight * above["mu"]),
        float((1 - weight) * below["sigma"] + weight * above["sigma"]),
    )
