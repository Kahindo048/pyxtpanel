"""Ramsey RESET test of omitted variables, from bare arrays (category B, L.8.8).

Validated against Stata 17: ``estat ovtest`` after ``regress`` (powers of the fitted values and
``rhs``, with a dummy among the regressors, with a collinear regressor) on a correctly specified
and on a misspecified model, and against an independent double-precision computation in Stata
(``estat ovtest`` builds the powers as single-precision variables, so the two agree to about
1e-6 only).  The regressors and the dependent variable are exported from Stata and given to the
function.  The cases are in ``validation/cases/diagnostics_specification_reset.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, as_array, register_test
from pyxtpanel.diagnostics._shared.ols import RANK_TOLERANCE

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray

_i18n.register(
    "diag",
    {
        "diag.reset.title": {
            "en": "Ramsey RESET test for omitted variables",
            "fr": "Test RESET de Ramsey de variables omises",
        },
        "diag.reset.null": {
            "en": "the model has no omitted variables",
            "fr": "le modèle n'a pas de variable omise",
        },
        "diag.reset.alternative": {
            "en": "powers of the fitted values (or of the regressors) are omitted variables",
            "fr": "des puissances des valeurs ajustées (ou des régresseurs) sont des variables "
            "omises",
        },
        "diag.reset.rows": {
            "en": "x must have one row per observation of y ({n}), got {got}.",
            "fr": "x doit avoir une ligne par observation de y ({n}), reçu {got}.",
        },
        "diag.reset.collinear": {
            "en": "the powers are collinear with the explanatory variables (typically because "
            "they are all indicator variables): the test is not possible.",
            "fr": "les puissances sont colinéaires aux variables explicatives (typiquement parce "
            "qu'elles sont toutes des indicatrices) : le test n'est pas possible.",
        },
        "diag.reset.perfect_fit": {
            "en": "the augmented model fits perfectly: the F statistic is not defined.",
            "fr": "le modèle augmenté s'ajuste parfaitement : la statistique F n'est pas définie.",
        },
        "diag.reset.dropped": {
            "en": "{dropped} of the {added} powers are collinear with the others and left out; "
            "the test has {df} degrees of freedom.",
            "fr": "{dropped} des {added} puissances sont colinéaires aux autres et écartées ; "
            "le test a {df} degrés de liberté.",
        },
    },
)

#: Residual sum of squares, relative to the total, under which the augmented fit is perfect.
PERFECT_FIT = 1e-22


@register_test(
    "reset",
    TestSpec(
        family="specification",
        category="B",
        wave=2,
        reference="Ramsey, J. B. (1969), Journal of the Royal Statistical Society B 31(2), 350-371",
        stata="estat ovtest; test",
        validation=(
            "reset_fitted",
            "reset_misspecified",
            "reset_rhs",
            "reset_rhs_dummy",
            "reset_collinear",
            "reset_fitted_exact",
            "reset_rhs_exact",
        ),
        inputs="arrays",
        balanced="na",
        r="lmtest::resettest",
        aliases=("ovtest", "ramsey"),
        summary="F test that powers of the fitted values (or regressors) add explanatory power.",
    ),
)
def reset(y: object, x: object, *, rhs: bool = False) -> TestResult:
    """Return Ramsey's RESET, the test that powers 2 to 4 of the fitted values add information.

    The model ``y = [1, x] b + e`` is fitted by least squares, then re-fitted with the powers 2, 3
    and 4 of the fitted values (default, as ``estat ovtest``) or of every regressor (``rhs``)
    added; the statistic is the ``F`` test that the added coefficients are zero.  As Stata, the
    fitted values or the regressors are first rescaled to ``[0, 1]`` (this does not change the
    test) and indicator regressors are not raised to powers.  Powers collinear with the others
    are left out and not counted in the degrees of freedom.

    Parameters
    ----------
    y
        ``(n,)`` dependent variable.
    x
        ``(n,)`` or ``(n, k)`` regressors, without the constant (which is always included).
    rhs
        Use the powers of the regressors instead of those of the fitted values.

    Raises
    ------
    SpecificationError
        If the shapes do not match.
    NumericalError
        If the powers are collinear with the model (all-indicator regressors) or the augmented
        model fits perfectly.
    """
    # Stata: estat ovtest -> regress y x yh2 yh3 yh4 with yh = (xb - min) / (max - min); test
    # yh2 yh3 yh4 (F, df, df_r); with rhs, powers 2-4 of each non-binary regressor in [0, 1]
    dep = as_array("y", y, ndim=1)
    raw = np.asarray(x, dtype=np.float64)
    reg = as_array("x", raw[:, None] if raw.ndim == 1 else raw, ndim=2)
    n = dep.size
    if reg.shape[0] != n:
        raise SpecificationError(
            _i18n.message("diag.reset.rows", n=n, got=reg.shape[0]), key="diag.reset.rows"
        )
    if rhs:
        base, powers = _rhs_columns(reg)
    else:
        base = reg
        fitted = _fit(dep, np.column_stack([np.ones(n), reg]))[0]
        powers = _powers(fitted) if np.max(fitted) > np.min(fitted) else np.empty((n, 0))
    if powers.shape[1] == 0:
        raise NumericalError(_i18n.message("diag.reset.collinear"), key="diag.reset.collinear")
    restricted = np.column_stack([np.ones(n), base])
    _, rss_r, rank_r = _fit(dep, restricted)
    _, rss_u, rank_u = _fit(dep, np.column_stack([restricted, powers]))
    q = rank_u - rank_r
    df_resid = n - rank_u
    if q < 1:
        raise NumericalError(_i18n.message("diag.reset.collinear"), key="diag.reset.collinear")
    if df_resid < 1 or rss_u <= PERFECT_FIT * float(np.sum((dep - dep.mean()) ** 2)):
        raise NumericalError(_i18n.message("diag.reset.perfect_fit"), key="diag.reset.perfect_fit")
    statistic = ((rss_r - rss_u) / q) / (rss_u / df_resid)
    notes: tuple[Note, ...] = ()
    if q < powers.shape[1]:
        notes = (
            Note(
                "diag.reset.dropped",
                {"dropped": powers.shape[1] - q, "added": powers.shape[1], "df": q},
            ),
        )
    return TestResult(
        name="reset",
        statistic=float(statistic),
        pvalue=float(stats.f.sf(statistic, q, df_resid)),
        statistic_name="F",
        df=(float(q), float(df_resid)),
        n_obs=n,
        options={"rhs": rhs},
        notes=notes,
    )


def _fit(y: FloatArray, z: FloatArray) -> tuple[FloatArray, float, int]:
    """Return the fitted values, the residual sum of squares and the rank of ``z``."""
    norm = np.sqrt(np.sum(z**2, axis=0))
    scaled = z / np.where(norm > 0, norm, 1.0)
    coefficients, _, rank, _ = np.linalg.lstsq(scaled, y, rcond=RANK_TOLERANCE)
    fitted: FloatArray = scaled @ coefficients
    resid = y - fitted
    return fitted, float(resid @ resid), int(rank)


def _powers(v: FloatArray) -> FloatArray:
    """Return the powers 2, 3 and 4 of ``v`` rescaled to ``[0, 1]``."""
    low, high = np.min(v), np.max(v)
    u = (v - low) / (high - low)
    return np.column_stack([u**2, u**3, u**4])


def _rhs_columns(reg: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Return the retained regressors and the powers of the non-binary ones (``rhs`` option)."""
    kept: list[FloatArray] = []
    powers: list[FloatArray] = []
    for column in reg.T:
        if np.min(column) == np.max(column):  # constant: Stata drops it from the model
            continue
        kept.append(column)
        if not np.all((column == 0.0) | (column == 1.0)):
            powers.append(_powers(column))
    n = reg.shape[0]
    base = np.column_stack(kept) if kept else np.empty((n, 0))
    return base, np.column_stack(powers) if powers else np.empty((n, 0))
