"""Augmented Dickey-Fuller regression of one series and its MacKinnon p-value (L.2.1).

The common brick of the panel unit-root tests that run one ADF regression per unit (Fisher,
Im-Pesaran-Shin).  It follows Stata's ``dfuller``: ``D.y`` is regressed on ``L.y``, ``lags``
lagged differences, a constant and, with a trend, the time since the first observation; the
statistic is the t ratio of ``L.y``.  Differences and lags follow the calendar, so a hole in the
series removes the equations it touches.

The p-value is MacKinnon's approximation (table ``mackinnon_adf.csv``) for a constant or a
constant and a trend, and the Student t distribution for the ``drift`` variant, exactly as
``dfuller`` does.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import stats

from pyxtpanel._types import FloatArray
from pyxtpanel.diagnostics._shared.tables import load_table

#: Deterministic terms of the regression: a constant, a constant and a trend, or a constant
#: with the Student t distribution under the null of a random walk with drift.
AdfCase = Literal["constant", "trend", "drift"]


@dataclass(frozen=True)
class AdfFit:
    """Result of one ADF regression.

    Attributes
    ----------
    tau
        t ratio of the lagged level (the Dickey-Fuller statistic).
    coefficient
        Coefficient of the lagged level.
    n
        Number of equations used.
    df_resid
        Residual degrees of freedom.
    lags
        Number of lagged differences.
    standard_error
        Standard error of the lagged level (``NaN`` when not defined).
    rmse
        Root mean squared error of the regression (``NaN`` when not defined).
    """

    tau: float
    coefficient: float
    n: int
    df_resid: int
    lags: int
    standard_error: float = math.nan
    rmse: float = math.nan


def adf_fit(series: FloatArray, lags: int, case: AdfCase = "constant") -> AdfFit | None:
    """Run the ADF regression of ``series`` and return ``None`` when it cannot be run.

    Parameters
    ----------
    series
        The values on the calendar, from the first to the last period of the unit, ``NaN`` where
        a period is missing (so that differences and lags follow the calendar).
    lags
        Number of lagged differences.
    case
        ``"constant"``, ``"trend"`` (constant and linear trend) or ``"drift"`` (constant).

    Returns
    -------
    AdfFit or None
        ``None`` when there are fewer equations than coefficients (``dfuller`` stops with
        "insufficient observations" or "no observations").  When the regression runs but the
        statistic is not defined (as many equations as coefficients, a collinear design, a
        perfect fit), ``tau`` is ``NaN``; :func:`adf_pvalue` then gives 1, as ``dfuller`` does
        with its missing statistic.
    """
    y = np.asarray(series, dtype=np.float64)
    count = y.size
    if count < lags + 2:
        return None
    first = lags + 1  # first position with a lagged level and `lags` lagged differences
    diff = np.full(count, np.nan)
    diff[1:] = y[1:] - y[:-1]
    columns = [y[first - 1 : count - 1]]
    columns.extend(diff[first - k : count - k] for k in range(1, lags + 1))
    if case == "trend":
        columns.append(np.arange(first, count, dtype=np.float64))
    columns.append(np.ones(count - first))
    design = np.column_stack(columns)
    response = diff[first:]
    usable = np.isfinite(design).all(axis=1) & np.isfinite(response)
    design, response = design[usable], response[usable]
    n, k = design.shape
    if n < k:
        return None  # Stata's regress stops with "insufficient observations"
    if n == k:
        return AdfFit(math.nan, math.nan, n, 0, lags)
    q, r = np.linalg.qr(design)
    if np.abs(np.diag(r)).min() <= 1e-10 * np.abs(np.diag(r)).max():
        return AdfFit(math.nan, math.nan, n, n - k, lags)
    coefficients = np.linalg.solve(r, q.T @ response)
    residual = response - design @ coefficients
    df_resid = n - k
    sigma2 = float(residual @ residual) / df_resid
    r_inverse = np.linalg.inv(r)
    error = float(np.sqrt(sigma2 * (r_inverse @ r_inverse.T)[0, 0]))
    if not np.isfinite(error) or error <= 0.0:
        return AdfFit(math.nan, float(coefficients[0]), n, df_resid, lags)
    return AdfFit(
        float(coefficients[0] / error),
        float(coefficients[0]),
        n,
        df_resid,
        lags,
        error,
        math.sqrt(sigma2),
    )


def adf_pvalue(fit: AdfFit, case: AdfCase = "constant") -> float:
    """Return the p-value of the Dickey-Fuller statistic as ``dfuller`` reports it.

    MacKinnon's normal approximation for ``"constant"`` and ``"trend"`` (0 below the lowest
    tabulated value and 1 above the highest), Student's t with the residual degrees of freedom
    for ``"drift"``.  A statistic that is not defined gives 1 with MacKinnon's approximation
    (Stata compares a missing value as larger than any number) and ``NaN`` with Student's t.
    """
    if math.isnan(fit.tau):
        return math.nan if case == "drift" else 1.0
    if case == "drift":
        return float(stats.t.cdf(fit.tau, fit.df_resid))
    row = load_table("mackinnon_adf.csv").set_index("case").loc["c" if case == "constant" else "ct"]
    tau = fit.tau
    if tau < row["tau_min"]:
        return 0.0
    if tau > row["tau_max"]:
        return 1.0
    part = "small" if tau <= row["tau_star"] else "large"
    h = sum(row[f"{part}_g{i}"] * tau**i for i in range(4))
    return float(stats.norm.cdf(h))


def select_lags(
    series: FloatArray,
    max_lags: int,
    criterion: Literal["aic", "bic", "hqic"],
    *,
    trend: bool = False,
    constant: bool = True,
) -> int:
    """Return the number of lags (0 to ``max_lags``) that minimises an information criterion.

    As in Stata's ``xtunitroot`` (``_xturadfreg``), every candidate is fitted on the equations of
    the regression with ``max_lags`` lags, so that the criteria compare the same sample; the
    criterion is ``(-2 ll + c (k + 1)) / n`` with ``c`` equal to 2 (AIC), ``ln n`` (BIC) or
    ``2 ln ln n`` (HQIC), ``ll`` the Gaussian log-likelihood and ``k`` the number of retained
    regressors without the constant.  Ties go to the larger number of lags.

    Parameters
    ----------
    series
        The unit's values on the calendar (``NaN`` where missing).
    max_lags
        Largest number of lags tried.
    criterion
        ``"aic"``, ``"bic"`` or ``"hqic"``.
    trend, constant
        Deterministic terms of the regressions.
    """
    y = np.asarray(series, dtype=np.float64)
    count = y.size
    first = max_lags + 1
    diff = np.full(count, np.nan)
    diff[1:] = y[1:] - y[:-1]
    response = diff[first:]
    level = y[first - 1 : count - 1]
    lagged = [diff[first - k : count - k] for k in range(1, max_lags + 1)]
    extra = [np.arange(first, count, dtype=np.float64)] if trend else []
    if constant:
        extra.append(np.ones(count - first))
    full = np.column_stack([level, *lagged, *extra]) if count > first else np.empty((0, 1))
    usable = np.isfinite(full).all(axis=1) & np.isfinite(response)
    n = int(usable.sum())
    best, best_value = max_lags, math.inf
    for p in range(max_lags, -1, -1):
        columns = [level, *lagged[:p], *extra]
        x = np.column_stack(columns)[usable]
        target = response[usable]
        coefficients, _, rank, _ = np.linalg.lstsq(x, target, rcond=1e-10)
        rss = max(float(np.sum((target - x @ coefficients) ** 2)), 1e-300)
        log_likelihood = -0.5 * n * (math.log(2.0 * math.pi) + 1.0 + math.log(rss / n))
        penalty = {"aic": 2.0, "bic": math.log(n), "hqic": 2.0 * math.log(math.log(n))}[criterion]
        value = (-2.0 * log_likelihood + penalty * (int(rank) - (1 if constant else 0) + 1)) / n
        if value < best_value:
            best, best_value = p, value
    return best
