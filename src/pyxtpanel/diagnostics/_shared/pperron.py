"""Phillips-Perron regression of one series and its MacKinnon p-value (L.2.10).

The brick of the Fisher-type panel unit-root test built on Phillips-Perron p-values.  It follows
Stata's ``pperron``: ``y`` is regressed on its lag, a constant and, with a trend, the time since the
first observation; ``Z(rho)`` and ``Z(t)`` correct the Dickey-Fuller statistics for serial
correlation and heteroskedasticity with a Bartlett estimate of the long-run variance of the
residuals (``lags`` Newey-West lags).  Lags follow the calendar, so a hole in the series removes the
equations it touches.  The p-value is MacKinnon's approximation for ``Z(t)``, the same table as the
ADF test (:func:`pyxtpanel.diagnostics._shared.adf.adf_pvalue`).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np

from pyxtpanel._types import FloatArray
from pyxtpanel.diagnostics._shared.adf import AdfFit, adf_pvalue


@dataclass(frozen=True)
class PpFit:
    """Result of one Phillips-Perron regression.

    Attributes
    ----------
    z_t
        The ``Z(t)`` statistic (``NaN`` when not defined).
    z_rho
        The ``Z(rho)`` statistic (``NaN`` when not defined).
    n
        Number of equations used.
    lags
        Number of Newey-West lags.
    """

    z_t: float
    z_rho: float
    n: int
    lags: int


def pp_fit(series: FloatArray, lags: int, case: Literal["constant", "trend"]) -> PpFit | None:
    """Run the Phillips-Perron regression of ``series``; ``None`` when Stata's ``pperron`` fails.

    Parameters
    ----------
    series
        The values on the calendar, from the first to the last period of the unit, ``NaN`` where a
        period is missing (so that lags follow the calendar).
    lags
        Number of Newey-West lags of the long-run variance.
    case
        ``"constant"`` or ``"trend"`` (constant and linear trend).

    Returns
    -------
    PpFit or None
        ``None`` with fewer than two equations (``pperron`` stops with "no observations" or
        "insufficient observations").  With as many equations as coefficients, a collinear design or
        a perfect fit the regression runs but the statistics are not defined: both are ``NaN`` and
        :func:`pp_pvalue` gives 1, as ``pperron`` does.
    """
    y = np.asarray(series, dtype=np.float64)
    count = y.size
    columns = [y[:-1]]
    if case == "trend":
        columns.append(np.arange(1, count, dtype=np.float64))
    columns.append(np.ones(count - 1))
    design = np.column_stack(columns) if count > 1 else np.empty((0, len(columns)))
    response = y[1:]
    usable = np.isfinite(design).all(axis=1) & np.isfinite(response)
    rows = np.flatnonzero(usable)
    n, k = int(rows.size), len(columns)
    if n < 2:
        return None
    undefined = PpFit(math.nan, math.nan, n, lags)
    if n <= k:
        return undefined
    x, target = design[usable], response[usable]
    q, r = np.linalg.qr(x)
    if np.abs(np.diag(r)).min() <= 1e-10 * np.abs(np.diag(r)).max():
        return undefined
    coefficients = np.linalg.solve(r, q.T @ target)
    residual = target - x @ coefficients
    rss = float(residual @ residual)
    s2 = rss / (n - k)
    r_inverse = np.linalg.inv(r)
    sig = math.sqrt(s2 * float((r_inverse @ r_inverse.T)[0, 0]))
    s = math.sqrt(s2)
    if not (math.isfinite(sig) and sig > 0.0 and s > 0.0):
        return undefined
    on_calendar = np.full(count, np.nan)
    on_calendar[rows + 1] = residual
    gamma0 = rss / n
    lam = gamma0
    for i in range(1, lags + 1):
        a, b = on_calendar[i:], on_calendar[:-i]
        both = np.isfinite(a) & np.isfinite(b)
        lam += 2.0 * (1.0 - i / (lags + 1.0)) * float(np.sum(a[both] * b[both])) / n
    rho = float(coefficients[0])
    z_rho = n * (rho - 1.0) - ((n * sig / s) ** 2) * (lam - gamma0) / 2.0
    if lam <= 0.0:
        return undefined
    z_t = (rho - 1.0) / sig * math.sqrt(gamma0 / lam) - 0.5 * (lam - gamma0) * (
        n * sig / s
    ) / math.sqrt(lam)
    return PpFit(float(z_t), float(z_rho), n, lags)


def pp_pvalue(fit: PpFit, case: Literal["constant", "trend"]) -> float:
    """Return MacKinnon's p-value of ``Z(t)`` as ``pperron`` reports it (1 when it is undefined)."""
    return adf_pvalue(AdfFit(fit.z_t, math.nan, fit.n, 0, fit.lags), case)
