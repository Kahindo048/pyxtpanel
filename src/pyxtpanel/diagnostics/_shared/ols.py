"""Least squares for the auxiliary regressions inside the tests (L.1.4).

The heteroskedasticity, specification and serial-correlation tests regress a transformation of
the residuals on some variables and read the explained sum of squares, ``R^2`` or ``F``.  This is
not a model of the package: it returns the few sums of squares those statistics need.  Like
Stata's ``regress`` it drops variables that are constant or collinear with the others, so the
degrees of freedom are the rank of the retained variables, not their number.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyxtpanel._types import FloatArray

#: Relative tolerance below which a centred column is constant or collinear.
RANK_TOLERANCE = 1e-10

#: A residual sum of squares below this fraction of the scale of the data is rounding noise: the
#: model fits perfectly, and a statistic built on its residuals would be noise too.
PERFECT_FIT = 1e-24


def is_perfect_fit(rss: float, scale: float) -> bool:
    """Return whether ``rss`` is nothing but rounding noise next to ``scale``.

    ``scale`` is the sum of squares of the variable that was fitted (about its mean when the model
    has a constant).  A zero ``scale`` (a constant variable) always counts as a perfect fit.
    """
    return rss <= PERFECT_FIT * scale


@dataclass(frozen=True)
class OLSFit:
    """Sums of squares of a regression on a constant and some variables.

    Attributes
    ----------
    n
        Number of observations.
    df_model
        Rank of the retained variables (the constant is not counted).
    df_resid
        ``n - df_model - 1``.
    rss, mss, tss
        Residual, explained and total sums of squares about the mean.
    """

    n: int
    df_model: int
    df_resid: int
    rss: float
    mss: float
    tss: float

    @property
    def r2(self) -> float:
        """``mss / tss`` (``NaN`` when the dependent variable is constant)."""
        return self.mss / self.tss if self.tss > 0 else float("nan")

    @property
    def f_statistic(self) -> float:
        """``(mss / df_model) / (rss / df_resid)`` (``NaN`` when not defined)."""
        if self.df_model < 1 or self.df_resid < 1 or self.rss <= 0:
            return float("nan")
        return (self.mss / self.df_model) / (self.rss / self.df_resid)


def ols_fit(y: FloatArray, x: FloatArray) -> OLSFit:
    """Regress ``y`` on a constant and the columns of ``x`` and return its sums of squares.

    Columns that are constant or exactly linear in the others are dropped (the rank is read after
    centring and scaling each column, which is what makes the tolerance relative).

    Parameters
    ----------
    y
        ``(n,)`` dependent variable.
    x
        ``(n, m)`` regressors, without the constant.
    """
    n = y.size
    centred = x - x.mean(axis=0)
    norm = np.sqrt(np.sum(centred**2, axis=0))
    raw = np.sqrt(np.sum(x**2, axis=0))
    keep = norm > RANK_TOLERANCE * np.where(raw > 0, raw, 1.0)
    yc = y - y.mean()
    tss = float(yc @ yc)
    if not keep.any():
        return OLSFit(n, 0, n - 1, tss, 0.0, tss)
    scaled = centred[:, keep] / norm[keep]
    u, s, _ = np.linalg.svd(scaled, full_matrices=False)
    rank = int(np.sum(s > RANK_TOLERANCE * s[0]))
    basis = u[:, :rank]  # orthonormal basis of the column space of the retained variables
    fitted = basis @ (basis.T @ yc)
    mss = float(fitted @ fitted)
    rss = float((yc - fitted) @ (yc - fitted))
    return OLSFit(n, rank, n - rank - 1, rss, mss, tss)


def ols_residuals(y: FloatArray, x: FloatArray) -> FloatArray:
    """Return the residuals of ``y`` on the columns of ``x`` as they are (no constant is added).

    Collinear columns are harmless: the residuals are those of the projection on the column
    space, whatever the coefficients chosen.  An empty ``x`` returns ``y``.
    """
    if x.shape[1] == 0:
        return y
    coefficients, *_ = np.linalg.lstsq(x, y, rcond=RANK_TOLERANCE)
    residual: FloatArray = y - x @ coefficients
    return residual
