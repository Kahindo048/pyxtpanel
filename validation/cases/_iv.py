"""Linear IV fits in NumPy that feed the over-identification tests of the Stata cases (L.6.2).

This is validation code, not part of the package: the package tests take the moments and the
weight matrix from the estimator.  Here a small two-stage least squares and two-step efficient
GMM (cluster-robust weight) build them, exactly as Stata's ``ivregress`` does, so that the
statistic the package computes from them can be compared with ``estat overid``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyxtpanel.diagnostics._shared.overid import moment_covariance
from validation import panels


@dataclass(frozen=True)
class IVFit:
    """Moments, weight and scale at the one-step and two-step estimates."""

    n_parameters: int
    n_groups: int
    n_instruments: int
    sargan_moments: np.ndarray
    sargan_weight: np.ndarray
    sargan_scale: float
    hansen_moments: np.ndarray
    hansen_weight: np.ndarray


def fit(panel: str, instruments: tuple[str, ...]) -> IVFit:
    """Two-stage least squares then two-step GMM of ``y`` on ``x`` (endogenous) and a constant."""
    frame = panels.load(panel)
    n = len(frame)
    y = frame["y"].to_numpy(float)
    x_mat = np.column_stack([frame["x"].to_numpy(float), np.ones(n)])
    z_mat = np.column_stack([*(frame[v].to_numpy(float) for v in instruments), np.ones(n)])
    ids = frame["id"].to_numpy()
    starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])

    w0 = np.linalg.inv(z_mat.T @ z_mat)
    zx, zy = z_mat.T @ x_mat, z_mat.T @ y
    beta1 = np.linalg.solve(zx.T @ w0 @ zx, zx.T @ w0 @ zy)
    e1 = y - x_mat @ beta1

    unit_moments = np.add.reduceat(z_mat * e1[:, None], starts, axis=0)
    w1 = np.linalg.inv(moment_covariance(unit_moments))
    beta2 = np.linalg.solve(zx.T @ w1 @ zx, zx.T @ w1 @ zy)
    e2 = y - x_mat @ beta2
    return IVFit(
        n_parameters=x_mat.shape[1],
        n_groups=starts.size,
        n_instruments=z_mat.shape[1],
        sargan_moments=z_mat.T @ e1,
        sargan_weight=w0,
        sargan_scale=float(e1 @ e1) / n,
        hansen_moments=z_mat.T @ e2,
        hansen_weight=w1,
    )
