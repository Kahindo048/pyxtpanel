"""Hadri stationarity test on a strongly balanced panel (category A, L.2.7).

Validated against Stata 17: ``xtunitroot hadri`` with a constant and with panel trends, with
cross-sectional means removed, with the heteroskedasticity-robust statistic and with the
Bartlett, Parzen and quadratic-spectral long-run variances.  The cases are in
``validation/cases/diagnostics_unitroot_hadri.py``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._core.inputs import plain_label

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.hadri.title": {
            "en": "Hadri LM stationarity test",
            "fr": "Test de stationnarité LM de Hadri",
        },
        "diag.hadri.null": {
            "en": "all panels are stationary",
            "fr": "toutes les séries du panel sont stationnaires",
        },
        "diag.hadri.alternative": {
            "en": "some panels contain unit roots",
            "fr": "certaines séries du panel ont une racine unitaire",
        },
        "diag.hadri.bad_deterministics": {
            "en": "deterministics must be 'constant' or 'trend', got {value!r}.",
            "fr": "deterministics doit valoir 'constant' ou 'trend', reçu {value!r}.",
        },
        "diag.hadri.bad_kernel": {
            "en": "kernel must be None, 'bartlett', 'parzen' or 'quadraticspectral', got "
            "{value!r}.",
            "fr": "kernel doit valoir None, 'bartlett', 'parzen' ou 'quadraticspectral', reçu "
            "{value!r}.",
        },
        "diag.hadri.bad_lags": {
            "en": "lags must be a non-negative number (an integer for the Bartlett and Parzen "
            "kernels) and needs a kernel, got {value!r}.",
            "fr": "lags doit être un nombre positif ou nul (entier pour les noyaux de Bartlett et "
            "de Parzen) et exige un noyau, reçu {value!r}.",
        },
        "diag.hadri.robust_and_kernel": {
            "en": "robust and kernel cannot be combined: robust scales each unit by its own "
            "variance, kernel estimates one long-run variance.",
            "fr": "robust et kernel ne se combinent pas : robust divise chaque unité par sa "
            "propre variance, kernel estime une seule variance de long terme.",
        },
        "diag.hadri.degenerate_unit": {
            "en": "the residual variance of unit {unit} is zero (it is constant after removing its "
            "mean or trend): the robust statistic is not defined.",
            "fr": "la variance résiduelle de l'unité {unit} est nulle (elle est constante une "
            "fois sa moyenne ou sa tendance retirée) : la statistique robuste n'est pas définie.",
        },
        "diag.hadri.degenerate": {
            "en": "the residual variance is zero (a unit, or every unit, is constant after "
            "removing its mean or trend): the statistic is not defined.",
            "fr": "la variance résiduelle est nulle (une unité, ou toutes, est constante une "
            "fois sa moyenne ou sa tendance retirée) : la statistique n'est pas définie.",
        },
    },
)

_KERNELS = ("bartlett", "parzen", "quadraticspectral")


@register_test(
    "hadri",
    TestSpec(
        family="unitroot",
        category="A",
        wave=2,
        reference="Hadri, K. (2000), Econometrics Journal 3(2), 148-161",
        stata="xtunitroot hadri",
        validation=(
            "ha_constant",
            "ha_trend",
            "ha_demean",
            "ha_robust",
            "ha_robust_trend",
            "ha_bartlett",
            "ha_parzen_trend",
            "ha_qs",
            "ha_kernel_default",
        ),
        inputs="panel",
        balanced="required",
        r="plm::purtest(test='hadri')",
        summary="LM test of the null that every series of the panel is stationary.",
    ),
)
def hadri(
    panel: PanelBase,
    variable: str,
    *,
    deterministics: Literal["constant", "trend"] = "constant",
    demean: bool = False,
    robust: bool = False,
    kernel: Literal["bartlett", "parzen", "quadraticspectral"] | None = None,
    lags: int | float | None = None,
) -> TestResult:
    """Return the Hadri ``z`` statistic of the stationarity LM test.

    Each unit's series is regressed on a constant (or a constant and a trend); with ``S_it`` the
    partial sums of the residuals, the LM statistic is the mean of ``S_it^2 / T^2`` over units
    and periods, scaled by an estimate of the residual variance, and
    ``z = sqrt(N) (LM - mean) / sqrt(variance)`` is standard normal.  The null is that every unit
    is stationary, so the test rejects for large values (right tail), the opposite of the other
    unit-root tests.

    The variance is the pooled residual variance, with the degrees-of-freedom correction
    ``T - 1`` (``T - 2`` with a trend); ``robust`` uses each unit's own variance (Hadri's
    heteroskedasticity-robust statistic); ``kernel`` uses one long-run variance, the average over
    units of the kernel estimate with ``lags`` lags (1 by default).  The panel must be strongly
    balanced and without gaps.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of the numeric variable tested.
    deterministics
        ``"constant"`` (default) or ``"trend"`` (a unit-specific linear trend).
    demean
        Remove the mean across units in each period first.
    robust
        Allow a different variance in each unit.
    kernel
        ``"bartlett"``, ``"parzen"`` or ``"quadraticspectral"`` to estimate a long-run variance
        that allows serial correlation.
    lags
        Number of lags of the kernel (an integer for Bartlett and Parzen; any non-negative number
        for the quadratic-spectral kernel, whose sum runs over its integer part).

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If an option has an invalid value, or ``robust`` is combined with ``kernel``.
    DataError
        If the panel is not strongly balanced, has a gap or is too short.
    NumericalError
        If the residual variance is zero.
    """
    # Stata: xtunitroot hadri -> e = residual of y on (1, t) per unit, S = cumulative sum of e;
    # LM = sum_i sum_t S^2 / (N T^2 sigma^2) with sigma^2 = sum e^2 / (N (T - 1 - trend)); robust:
    # sum_i (sum_t S^2 / sigma_i^2) / (N T^2), sigma_i^2 = sum_t e^2 / (T - 1 - trend); kernel:
    # sigma^2 replaced by the mean over units of e'e/T + 2 sum_l w(l/(L+1)) sum_t e_t e_{t-l}/T;
    # z = sqrt(N)(LM - mu)/sqrt(v), (mu, v) = (1/6, 1/45) or (1/15, 11/6300); p = 1 - Phi(z).
    if deterministics not in ("constant", "trend"):
        raise SpecificationError(
            _i18n.message("diag.hadri.bad_deterministics", value=deterministics),
            key="diag.hadri.bad_deterministics",
        )
    if kernel is not None and kernel not in _KERNELS:
        raise SpecificationError(
            _i18n.message("diag.hadri.bad_kernel", value=kernel), key="diag.hadri.bad_kernel"
        )
    if robust and kernel is not None:
        raise SpecificationError(
            _i18n.message("diag.hadri.robust_and_kernel"), key="diag.hadri.robust_and_kernel"
        )
    if lags is not None and not _valid_lags(lags, kernel):
        raise SpecificationError(
            _i18n.message("diag.hadri.bad_lags", value=lags), key="diag.hadri.bad_lags"
        )
    trend = deterministics == "trend"
    arrays = panel_variables(
        panel,
        [variable],
        test="hadri",
        balanced="required",
        consecutive=True,
        min_obs=3 if trend else 2,
        min_units=2,
    )
    y = arrays.as_cube()[:, :, 0].copy()
    n_units, periods = y.shape
    if demean:
        y -= y.mean(axis=0)
    e = _residuals(y, trend)
    # a unit that is constant (or exactly linear with a trend) leaves only rounding noise
    unit_noise = np.sum(e**2, axis=1) <= 1e-20 * np.maximum(np.sum(y**2, axis=1), 1e-300)
    partial = np.cumsum(e, axis=1)
    sum_squares = np.sum(partial**2, axis=1)  # per unit
    divisor = periods - 1 - (1 if trend else 0)
    unit_variance = np.sum(e**2, axis=1) / divisor
    bandwidth = None
    if kernel is not None:
        bandwidth = 1 if lags is None else lags
        long_run = float(np.mean(_long_run_variance(e, kernel, float(bandwidth))))
        scale = long_run
    elif robust:
        scale = float("nan")
    else:
        scale = float(np.sum(e**2)) / (n_units * divisor)
    if robust:
        if np.any(unit_noise):
            first = plain_label(arrays.unit_labels[int(np.argmax(unit_noise))])
            raise NumericalError(
                _i18n.message("diag.hadri.degenerate_unit", unit=first),
                key="diag.hadri.degenerate_unit",
            )
        lm = float(np.sum(sum_squares / unit_variance)) / (n_units * periods**2)
    else:
        if not scale > 0.0 or bool(np.all(unit_noise)):
            raise NumericalError(
                _i18n.message("diag.hadri.degenerate"), key="diag.hadri.degenerate"
            )
        lm = float(np.sum(sum_squares)) / (n_units * periods**2) / scale
    mean, variance = (1.0 / 15.0, 11.0 / 6300.0) if trend else (1.0 / 6.0, 1.0 / 45.0)
    z = math.sqrt(n_units) * (lm - mean) / math.sqrt(variance)
    details = pd.DataFrame(
        {"value": [lm, mean, variance]},
        index=pd.Index(["lm", "mean_lm", "variance_lm"], name="quantity"),
    )
    options: dict[str, object] = {
        "deterministics": deterministics,
        "demean": demean,
        "robust": robust,
        "kernel": kernel,
    }
    if bandwidth is not None:
        options["lags"] = bandwidth
    return TestResult(
        name="hadri",
        statistic=float(z),
        pvalue=float(stats.norm.sf(z)),
        statistic_name="z",
        df=(),
        n_obs=n_units * periods,
        n_groups=n_units,
        variable=variable,
        options=options,
        details=details,
    )


def _valid_lags(lags: object, kernel: str | None) -> bool:
    """Whether ``lags`` is a usable kernel length (it needs a kernel; an integer unless QS)."""
    if kernel is None or isinstance(lags, bool):
        return False
    if not isinstance(lags, (int, float, np.integer, np.floating)):
        return False
    value = float(lags)
    if not math.isfinite(value) or value < 0.0:
        return False
    return kernel == "quadraticspectral" or value == int(value)


def _residuals(y: FloatArray, trend: bool) -> FloatArray:
    """Return the residuals of each row of ``y`` on a constant (and a linear trend)."""
    periods = y.shape[1]
    columns = [np.ones(periods)]
    if trend:
        columns.append(np.arange(1.0, periods + 1))
    basis, _ = np.linalg.qr(np.column_stack(columns))
    return np.asarray(y - (y @ basis) @ basis.T)


def _kernel_weight(kernel: str, z: float) -> float:
    """Return the weight of ``kernel`` at ``z = l / (L + 1)``."""
    if kernel == "bartlett":
        return 1.0 - z
    if kernel == "parzen":
        return 1.0 - 6.0 * z**2 + 6.0 * z**3 if z <= 0.5 else 2.0 * (1.0 - z) ** 3
    if z == 1.0:
        return 1.0
    theta = 6.0 * math.pi * z / 5.0
    return 3.0 * (math.sin(theta) / theta - math.cos(theta)) / theta**2


def _long_run_variance(e: FloatArray, kernel: str, lags: float) -> FloatArray:
    """Return the kernel estimate of the long-run variance of each row of ``e``."""
    periods = e.shape[1]
    out = np.sum(e**2, axis=1) / periods
    for lag in range(1, int(math.floor(lags)) + 1):
        weight = _kernel_weight(kernel, lag / (lags + 1.0))
        out = out + 2.0 * weight * np.sum(e[:, lag:] * e[:, :-lag], axis=1) / periods
    return np.asarray(out)
