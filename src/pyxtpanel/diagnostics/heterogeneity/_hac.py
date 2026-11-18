"""Blomquist-Westerlund HAC version of Swamy's slope heterogeneity statistic (L.8.13).

The computation of ``deltatesthac`` of ``xthst`` (Stata): from the data with the constant, the
partialled variables and the cross-sectional averages already removed, it estimates for every unit
a long-run covariance matrix of the scores ``x_it e_it`` with a kernel, weights the unit's slope
by it and returns the statistic ``S_HAC`` with the average bandwidth.  Several details of the Mata
code are reproduced as they are (marked ``STATA-QUIRK``).
"""

from __future__ import annotations

import math

import numpy as np
from scipy import linalg

from pyxtpanel._types import FloatArray, IntArray

#: Singular values of the whitening matrix are capped here (Andrews and Monahan, 1992, p. 957).
SINGULAR_CAP = 0.97

#: A pivot below this fraction of the original diagonal is a dependent row for ``invsym``.
TOLERANCE = 1e-13


def hac_swamy(
    y: FloatArray,
    x: FloatArray,
    offsets: IntArray,
    *,
    kernel: str,
    bandwidth: int | None,
    whitening: bool,
) -> tuple[float, float]:
    """Return ``S_HAC`` and the average bandwidth.

    Parameters
    ----------
    y, x
        ``(n,)`` and ``(n, K)``: the data after the constant, the partialled variables and the
        cross-sectional averages are removed unit by unit, sorted by unit and period.
    offsets
        ``(N + 1,)`` start of each unit.
    kernel
        ``"bartlett"``, ``"qs"`` (quadratic spectral) or ``"truncated"``.
    bandwidth
        Bandwidth of every unit, or ``None`` for the automatic choice of ``xthst`` (it differs by
        kernel and by unit; the average is returned).
    whitening
        Pre-whiten the scores with a first-order vector autoregression (Andrews and Monahan).
    """
    # Stata: xthst, hac -> deltatesthac(): FE slope on the partialled data, scores x_it e_it, per
    # unit V_i = G_0 + sum_j k(j) (G_j + c G_j') with the unit's bandwidth, Q_i = X_i'X_i / T_i,
    # beta = (sum T_i Q_i V_i^-1 Q_i)^-1 sum Q_i V_i^-1 X_i'Y_i, S = sum T_i (b_i - beta)' Q_i
    # V_i^-1 Q_i (b_i - beta); V_i is inverted with invsym (upper triangle only)
    k = x.shape[1]
    pooled = np.linalg.solve(x.T @ x, x.T @ y)
    scores = x * (y - x @ pooled)[:, None]
    units = offsets.size - 1
    beta_low, beta_up = np.zeros((k, k)), np.zeros(k)
    weights: list[FloatArray] = []
    bandwidth_total = 0.0
    for i in range(units):
        rows = slice(int(offsets[i]), int(offsets[i + 1]))
        series, whitener = scores[rows], None
        if whitening:
            series, whitener = _whiten(series)
        if bandwidth is None:
            used, lags = _automatic_bandwidth(series, kernel)
            bandwidth_total += lags
        else:
            used, lags = float(bandwidth), float(bandwidth)
            bandwidth_total = float(bandwidth) * units
        long_run = _long_run(series, kernel, used, lags)
        if whitener is not None:
            factor = _invsym(np.eye(k) - whitener)
            long_run = factor @ long_run @ factor.T
        xi, yi = x[rows], y[rows]
        q = xi.T @ xi / xi.shape[0]
        inverse = _invsym(long_run)
        qvq = q @ inverse @ q
        weights.append(qvq)
        beta_low += xi.shape[0] * qvq
        beta_up += q @ inverse @ (xi.T @ yi)  # STATA-QUIRK: X'Y is not divided by T here
    beta = _invsym(beta_low) @ beta_up
    total = 0.0
    for i, qvq in enumerate(weights):
        rows = slice(int(offsets[i]), int(offsets[i + 1]))
        xi, yi = x[rows], y[rows]
        own = np.linalg.solve(xi.T @ xi, xi.T @ yi)
        total += xi.shape[0] * float((own - beta) @ qvq @ (own - beta))
    return total, bandwidth_total / units


def _invsym(matrix: FloatArray) -> FloatArray:
    """Inverse of the symmetric matrix whose upper triangle is ``matrix`` (Mata's ``invsym``).

    The matrix is swept pivot by pivot, always on the largest remaining diagonal element, which is
    what Mata does on the indefinite matrices of the truncated kernel (the order decides which
    row is dropped).  A pivot that is not positive, or is negligible next to the original
    diagonal, marks a dependent row: its row and column of the result are zero.
    """
    upper = np.triu(matrix)
    a = upper + np.triu(matrix, 1).T
    n = a.shape[0]
    original = np.abs(np.diag(a)).copy()
    remaining = list(range(n))
    while remaining:
        k = max(remaining, key=lambda index: a[index, index])
        remaining.remove(k)
        pivot = a[k, k]
        if not pivot > TOLERANCE * original[k]:
            a[k, :] = 0.0
            a[:, k] = 0.0
            continue
        a[k, :] /= pivot
        for i in range(n):
            if i != k:
                factor = a[i, k]
                a[i, :] -= factor * a[k, :]
                a[i, k] = -factor / pivot
        a[k, k] = 1.0 / pivot
    return a


def _whiten(series: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Return the residuals of a VAR(1) of the scores and its (capped) coefficient matrix."""
    # STATA-QUIRK: the singular values are capped at 0.97 and the matrix is rebuilt as
    # U diag(s) V, not U diag(s) V' (the third factor of Mata's svd() is already V'), so the
    # coefficient matrix of a unit with several regressors is not the one that was estimated;
    # the result then depends on the signs of the singular vectors, which Mata and SciPy give
    # alike with one or two regressors and differently for about a quarter of the units with three
    lag, current = series[:-1], series[1:]
    coefficients = _invsym(lag.T @ lag) @ (lag.T @ current)
    u, s, vt = linalg.svd(coefficients, lapack_driver="gesvd")
    s = np.minimum(s, SINGULAR_CAP)
    rebuilt = u @ np.diag(s) @ vt.T
    return current - lag @ rebuilt, rebuilt


def _automatic_bandwidth(series: FloatArray, kernel: str) -> tuple[float, float]:
    """Return the bandwidth of a unit and the number of lags summed (Newey-West, Andrews)."""
    t = series.shape[0]
    if kernel == "truncated":
        # Newey and West (1994), p. 641
        lags = float(math.floor(4.0 * (t / 100.0) ** (1.0 / 5.0)))
        return lags, lags
    if kernel == "qs":
        up = low = 0.0
        for column in range(series.shape[1]):  # Andrews (1991), p. 830, AR(1) of each score
            previous, current = series[:-1, column], series[1:, column]
            rho = float(previous @ current) / float(previous @ previous)
            innovation = current - previous * rho
            sigma2 = float(innovation @ innovation) / (t - 1)
            up += (2.0 * rho * sigma2 / (1.0 - rho) ** 4) ** 2
            low += (sigma2 / (1.0 - rho) ** 2) ** 2
        # STATA-QUIRK: the ratio is squared once more than in Andrews' alpha(2)
        return 1.3221 * ((up / low) ** 2 * t) ** (1.0 / 5.0), float(t - 1)
    return _bartlett_bandwidth(series)


def _bartlett_bandwidth(series: FloatArray) -> tuple[float, float]:
    t = series.shape[0]
    top = int(4.0 * (t / 100.0) ** (2.0 / 9.0))  # mstar, truncated toward zero
    with np.errstate(invalid="ignore", divide="ignore"):
        sigma0 = np.sqrt(series.T @ series / t)  # elementwise: a negative entry is missing
        sigma1 = np.zeros_like(sigma0)
        for j in range(1, top + 1):
            cross = np.sqrt(series[: t - j].T @ series[j:] / (t - j))
            sigma0 = sigma0 + 2.0 * cross
            sigma1 = sigma1 + 2.0 * cross * j
        candidate = np.floor(1.1447 * ((sigma1 / sigma0) ** 2 * t) ** (1.0 / 3.0))
    finite = candidate[np.isfinite(candidate)]
    # Mata's min() of the matrix ignores missing values and is missing when all are; the outer
    # min() with mstar then gives mstar
    bandwidth = min(float(finite.min()), float(top)) if finite.size else float(top)
    lags = float(t - 1) if bandwidth > t - 1 else bandwidth
    return bandwidth, lags


def _long_run(series: FloatArray, kernel: str, bandwidth: float, lags: float) -> FloatArray:
    t = series.shape[0]
    total = series.T @ series / t
    j = 1
    while j <= lags:
        gamma = series[j:].T @ series[: t - j] / t
        sign = 1.0
        if kernel == "bartlett":
            weight = max(1.0 - j / (bandwidth + 1.0), 0.0)
            sign = -1.0  # STATA-QUIRK: Bartlett adds the difference G_j - G_j', not the sum
        elif kernel == "qs":
            argument = j / bandwidth
            angle = 6.0 * math.pi * argument / 5.0
            weight = (
                25.0
                / (12.0 * math.pi**2 * argument**2)
                * (math.sin(angle) / angle - math.cos(angle))
            )
        else:
            weight = 1.0
        total = total + weight * (gamma + sign * gamma.T)
        j += 1
    return total
