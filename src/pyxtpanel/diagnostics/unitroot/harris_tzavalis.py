"""Harris-Tzavalis unit-root test on a strongly balanced panel (category A, L.2.3).

Validated against Stata 17: ``xtunitroot ht`` with a constant, without constant, with panel
trends, with cross-sectional means removed (``demean``) and with the small-sample adjustment of
``T`` (``altt``), on panels of random walks, of stationary series and of series with a trend.  The
cases are in ``validation/cases/diagnostics_unitroot_harris_tzavalis.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, panel_variables, register_test

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.harris_tzavalis.title": {
            "en": "Harris-Tzavalis unit-root test",
            "fr": "Test de racine unitaire de Harris-Tzavalis",
        },
        "diag.harris_tzavalis.null": {
            "en": "panels contain unit roots",
            "fr": "les séries du panel ont une racine unitaire",
        },
        "diag.harris_tzavalis.alternative": {
            "en": "panels are stationary (common AR parameter, N large, T fixed)",
            "fr": "les séries du panel sont stationnaires (paramètre AR commun, N grand, T fixe)",
        },
        "diag.harris_tzavalis.bad_deterministics": {
            "en": "deterministics must be 'none', 'constant' or 'trend', got {value!r}.",
            "fr": "deterministics doit valoir 'none', 'constant' ou 'trend', reçu {value!r}.",
        },
    },
)

#: Smallest number of periods per unit for each specification (the variance has a factor
#: ``T - 1`` or ``T - 2``, and a trend fitted on ``T - 1`` points needs at least 3).
_MIN_PERIODS = {"none": 3, "constant": 3, "trend": 4}


@register_test(
    "harris_tzavalis",
    TestSpec(
        family="unitroot",
        category="A",
        wave=1,
        reference="Harris, R. D. F. and Tzavalis, E. (1999), J. of Econometrics 91(2), 201-226",
        stata="xtunitroot ht",
        validation=(
            "ht_walk",
            "ht_stationary",
            "ht_none",
            "ht_trend",
            "ht_demean",
            "ht_altt",
        ),
        inputs="panel",
        aliases=("ht",),
        balanced="required",
        summary="Unit-root test of a common AR parameter for panels with large N and fixed T.",
    ),
)
def harris_tzavalis(
    panel: PanelBase,
    variable: str,
    *,
    deterministics: Literal["none", "constant", "trend"] = "constant",
    demean: bool = False,
    small_sample_t: bool = False,
) -> TestResult:
    """Return the Harris-Tzavalis ``z`` statistic of the pooled AR parameter ``rho``.

    ``rho`` is the least-squares slope of ``y_it`` on ``y_i,t-1`` after removing the unit means
    (``"constant"``), the unit means and trends (``"trend"``) or nothing (``"none"``).  Its mean
    and variance under the unit root are exact functions of ``T`` (``N -> infinity`` with ``T``
    fixed), ``z = (rho - mean) / sqrt(variance)`` is standard normal, and the test rejects for
    small values: the p-value is the left tail.

    The panel must be strongly balanced and without gaps.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of the numeric variable tested.
    deterministics
        ``"constant"`` (unit means removed, the default), ``"none"`` or ``"trend"`` (unit means
        and linear trends removed).
    demean
        Remove the mean across units in each period first, which reduces the effect of
        cross-sectional dependence.
    small_sample_t
        Use ``T - 1`` instead of ``T`` in the moments (Stata's ``altt``: one observation is lost
        to the lag).

    Raises
    ------
    SpecificationError
        If ``deterministics`` is not one of the three values.
    DataError
        If the panel is not strongly balanced, has a gap, or has too few periods.
    """
    # Stata: xtunitroot ht -> rho from reg y L.y (noconstant), xtreg y L.y, fe, or the same on
    # series with the unit mean and trend removed over t = 2..T; mean and variance of rho from
    # the closed forms of Harris and Tzavalis (models 1, 2, 3); z = (rho - mean)/sd; p = Phi(z).
    if deterministics not in _MIN_PERIODS:
        raise SpecificationError(
            _i18n.message("diag.harris_tzavalis.bad_deterministics", value=deterministics),
            key="diag.harris_tzavalis.bad_deterministics",
        )
    arrays = panel_variables(
        panel,
        [variable],
        test="harris_tzavalis",
        balanced="required",
        consecutive=True,
        min_obs=_MIN_PERIODS[deterministics],
        min_units=2,
    )
    y = arrays.as_cube()[:, :, 0].copy()  # (units, periods)
    n, periods = y.shape
    if demean:
        y -= y.mean(axis=0)
    dependent, lagged = y[:, 1:], y[:, :-1]
    if deterministics != "none":
        dependent, lagged = _remove(dependent, deterministics), _remove(lagged, deterministics)
    rho = float(np.sum(lagged * dependent) / np.sum(lagged * lagged))
    t = periods - 1 if small_sample_t else periods
    mean, variance = _moments(deterministics, float(t), n)
    z = (rho - mean) / np.sqrt(variance)
    details = pd.DataFrame(
        {"value": [rho, mean, variance]},
        index=pd.Index(["rho", "mean_rho", "variance_rho"], name="quantity"),
    )
    return TestResult(
        name="harris_tzavalis",
        statistic=float(z),
        pvalue=float(stats.norm.cdf(z)),
        statistic_name="z",
        df=(),
        n_obs=n * periods,
        n_groups=n,
        variable=variable,
        options={
            "deterministics": deterministics,
            "demean": demean,
            "small_sample_t": small_sample_t,
        },
        details=details,
    )


def _remove(x: np.ndarray, deterministics: str) -> np.ndarray:
    """Remove each unit's mean (and linear trend) from the rows of ``x``."""
    if deterministics == "constant":
        return x - x.mean(axis=1, keepdims=True)
    count = x.shape[1]
    basis, _ = np.linalg.qr(np.column_stack([np.ones(count), np.arange(1.0, count + 1)]))
    return x - (x @ basis) @ basis.T


def _moments(deterministics: str, t: float, n: int) -> tuple[float, float]:
    """Mean and variance of ``rho`` under the unit root, for ``n`` units and ``t`` periods."""
    if deterministics == "none":
        return 1.0, 2.0 / (t * (t - 1.0)) / n
    if deterministics == "constant":
        mean = 1.0 - 3.0 / (t + 1.0)
        variance = 3.0 * (17.0 * t**2 - 20.0 * t + 17.0) / (5.0 * (t - 1.0) * (t + 1.0) ** 3)
        return mean, variance / n
    mean = 1.0 - 15.0 / (2.0 * (t + 2.0))
    variance = 15.0 * (193.0 * t**2 - 728.0 * t + 1147.0) / (112.0 * (t + 2.0) ** 3 * (t - 2.0))
    return mean, variance / n
