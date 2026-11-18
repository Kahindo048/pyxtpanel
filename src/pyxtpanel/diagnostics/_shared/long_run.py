"""Long-run covariance with the Newey-West (1994) automatic lag selection (L.4).

The common brick of the panel cointegration tests.  It follows Stata's ``xtcointtest``: each
unit's residual matrix ``w`` (observations x series) is centred, its contemporaneous covariance
``Sigma`` and its kernel long-run covariance ``Omega`` are computed, and the panel's matrices are
the averages over units.  With the automatic selection every unit picks its own number of lags
from the *trace* of its autocovariance matrices (Newey and West's weights ``h = (1, ..., 1)``
applied to a vector process, which Stata reduces to the trace of the autocovariance matrices).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np

from pyxtpanel._types import FloatArray

Kernel = Literal["bartlett", "parzen", "quadraticspectral"]

#: kernel -> (q, c, exponent of the preliminary lag count, exponent of ``n`` in the bandwidth)
_NEWEY_WEST = {
    "bartlett": (1, 1.1447, 2.0 / 9.0, 1.0 / 3.0),
    "parzen": (2, 2.6614, 4.0 / 25.0, 1.0 / 5.0),
    "quadraticspectral": (2, 1.3221, 2.0 / 25.0, 1.0 / 5.0),
}


@dataclass(frozen=True)
class LongRun:
    """Average contemporaneous and long-run covariance of the units' residuals.

    Attributes
    ----------
    sigma, omega
        Contemporaneous and long-run covariance matrices (``k x k``), averaged over units.
    lags
        Mean of the lag counts Stata reports (capped for the automatic selection).
    bandwidth
        Mean bandwidth of the kernel (the number of lags plus one when the lags are fixed).
    """

    sigma: FloatArray
    omega: FloatArray
    lags: float
    bandwidth: float


def kernel_weight(kernel: str, lag: int, bandwidth: float) -> float:
    """Return the weight of ``kernel`` at ``lag`` for a bandwidth (the number of lags plus one)."""
    z = lag / bandwidth
    if kernel == "bartlett":
        return 1.0 - z
    if kernel == "parzen":
        return 1.0 - 6.0 * z**2 + 6.0 * z**3 if z <= 0.5 else 2.0 * (1.0 - z) ** 3
    theta = 6.0 * math.pi * z / 5.0
    return 1.0 if z == 0.0 else 3.0 * (math.sin(theta) / theta - math.cos(theta)) / theta**2


def newey_west_lags(w: FloatArray, kernel: str) -> tuple[float, float]:
    """Return Newey and West's (1994) number of lags for ``w`` and the uncapped bandwidth.

    The preliminary autocovariances are the traces of ``w_t' w_{t-j} / n`` for ``j = 1`` to
    ``int(4 (n / 100)^a)``; the lag count is ``min(int(m), cap)`` for the Bartlett and Parzen
    kernels and ``min(m, cap)`` for the quadratic-spectral one.
    """
    n = w.shape[0]
    q, c, pre, exponent = _NEWEY_WEST[kernel]
    cap = int(4.0 * (n / 100.0) ** pre)
    s_0 = float(np.sum(w * w)) / n
    s_q = 0.0
    for j in range(1, cap + 1):
        trace = float(np.sum(w[j:] * w[:-j])) / n
        s_0 += 2.0 * trace
        s_q += 2.0 * j**q * trace
    if s_0 == 0.0:
        return math.nan, math.nan
    gamma = c * ((s_q / s_0) ** 2) ** (1.0 / (2 * q + 1))
    bandwidth = gamma * n**exponent
    if kernel == "quadraticspectral":
        return min(bandwidth, float(cap)), bandwidth
    return float(min(int(bandwidth), cap)), float(int(bandwidth))


def long_run_covariance(
    residuals: FloatArray, kernel: str, lags: int | float | None
) -> LongRun | None:
    """Average ``Sigma`` and ``Omega`` of the units' ``(units, observations, series)`` residuals.

    Parameters
    ----------
    residuals
        Residuals of every unit, ``(units, observations, series)``.
    kernel
        ``"bartlett"``, ``"parzen"`` or ``"quadraticspectral"``.
    lags
        A fixed number of lags, or ``None`` for the Newey-West selection of each unit.

    With a fixed number of lags ``m`` the bandwidth of the kernel is ``m + 1``; with the automatic
    selection it is the selected ``int(m_hat)`` itself, so the sum has ``int(m_hat) - 1`` terms
    (a unit with bandwidth 1 gets none).  The quadratic-spectral kernel always sums all ``n - 1``
    lags.  Stata reports the *capped* count ``min(int(m_hat), int(4 (n / 100)^a))`` as the lags
    used, but computes the kernel with the uncapped bandwidth: ``lags`` holds the reported mean
    and ``bandwidth`` the mean bandwidth.

    Returns
    -------
    LongRun or None
        ``None`` when the automatic selection is not defined for a unit (a zero long-run
        variance estimate of its preliminary sum).
    """
    n_units, n, k = residuals.shape
    sigma = np.zeros((k, k))
    omega = np.zeros((k, k))
    reported, raw = [], []
    for unit in range(n_units):
        w = residuals[unit] - residuals[unit].mean(axis=0)
        if lags is None:
            capped, bandwidth = newey_west_lags(w, kernel)
            if math.isnan(capped):
                return None
        else:
            capped, bandwidth = float(lags), float(lags) + 1.0
        s = w.T @ w / n
        o = s.copy()
        top = n - 1 if kernel == "quadraticspectral" else min(int(math.ceil(bandwidth)) - 1, n - 1)
        for j in range(1, top + 1):
            gamma = w[j:].T @ w[:-j] / n
            o += kernel_weight(kernel, j, bandwidth) * (gamma + gamma.T)
        sigma += s
        omega += o
        reported.append(capped)
        raw.append(bandwidth)
    return LongRun(sigma / n_units, omega / n_units, float(np.mean(reported)), float(np.mean(raw)))
