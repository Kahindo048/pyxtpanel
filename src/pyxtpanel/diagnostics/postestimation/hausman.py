"""Hausman specification test from two sets of coefficients and covariances (category B, L.6.4).

Validated against Stata 17: ``hausman`` after ``xtreg, fe`` / ``xtreg, re`` (the coefficients and
covariances are exported from Stata and given to the function), including a case whose
difference of covariances is not positive definite (negative statistic, as Stata reports it) and
a case with ``df()``; the cases are in
``validation/cases/diagnostics_postestimation_hausman.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, as_array, register_test

_i18n.register(
    "diag",
    {
        "diag.hausman.title": {
            "en": "Hausman specification test",
            "fr": "Test de spécification de Hausman",
        },
        "diag.hausman.null": {
            "en": "the difference in coefficients is not systematic",
            "fr": "la différence entre les coefficients n'est pas systématique",
        },
        "diag.hausman.alternative": {
            "en": "the difference in coefficients is systematic",
            "fr": "la différence entre les coefficients est systématique",
        },
        "diag.hausman.not_positive_definite": {
            "en": "V_consistent - V_efficient is not positive definite; the statistic uses its "
            "generalised inverse as Stata does.",
            "fr": "V_cohérent - V_efficient n'est pas définie positive ; la statistique utilise "
            "son inverse généralisée, comme Stata.",
        },
        "diag.hausman.rank": {
            "en": "the rank of the differenced covariance matrix ({rank}) is below the {k} "
            "coefficients tested; the degrees of freedom are the rank.",
            "fr": "le rang de la différence des matrices de covariance ({rank}) est inférieur "
            "aux {k} coefficients testés ; les degrés de liberté sont le rang.",
        },
        "diag.hausman.negative": {
            "en": "the statistic is negative ({value:.4g}): the model assumptions are not met "
            "(the efficient estimator is not efficient); no p-value.",
            "fr": "la statistique est négative ({value:.4g}) : les hypothèses du modèle ne "
            "sont pas satisfaites (l'estimateur efficient ne l'est pas) ; pas de p-value.",
        },
        "diag.hausman.shape": {
            "en": "{name} must have {expected} entries, got {got}.",
            "fr": "{name} doit avoir {expected} éléments, reçu {got}.",
        },
    },
)


@register_test(
    "hausman",
    TestSpec(
        family="postestimation",
        category="B",
        wave=1,
        reference="Hausman, J. A. (1978), Econometrica 46(6), 1251-1271",
        stata="hausman",
        validation=(
            "hausman_fe_re",
            "hausman_fe_re_invalid",
            "hausman_swapped_negative",
            "hausman_df_option",
        ),
        inputs="arrays",
        balanced="na",
        r="plm::phtest",
        summary="Hausman test comparing a consistent and an efficient estimator.",
    ),
)
def hausman(
    b_consistent: object,
    b_efficient: object,
    v_consistent: object,
    v_efficient: object,
    *,
    keep: Sequence[int] | None = None,
    names: Sequence[str] | None = None,
    df: int | None = None,
) -> TestResult:
    """Return the Hausman statistic ``(b1 - b2)' (V1 - V2)^- (b1 - b2)``.

    ``b1`` and ``V1`` come from the estimator that is consistent under both hypotheses (fixed
    effects, instrumental variables), ``b2`` and ``V2`` from the one that is efficient under the
    null (random effects, least squares).  The generalised inverse is the one of Stata's
    ``hausman``: singular values below ``max * k * 1e-9`` are set to zero, and the degrees of
    freedom are the rank of ``V1 - V2``.  A negative statistic is returned as such, with no
    p-value, like Stata ("chi2 < 0").

    Parameters
    ----------
    b_consistent, b_efficient
        ``(k,)`` coefficient vectors of the two estimators, in the same order.
    v_consistent, v_efficient
        ``(k, k)`` covariance matrices of the coefficients.
    keep
        Positions of the coefficients to test (to leave the constant out, as Stata does by
        default); all by default.
    names
        Coefficient names, for the table of differences.
    df
        Degrees of freedom to use instead of the rank (Stata's ``df()``).

    Raises
    ------
    SpecificationError
        If the vectors and matrices do not have matching sizes.
    """
    b1 = as_array("b_consistent", b_consistent, ndim=1)
    k = b1.size
    b2 = as_array("b_efficient", b_efficient, ndim=1, shape=(k,))
    v1 = as_array("v_consistent", v_consistent, ndim=2, shape=(k, k))
    v2 = as_array("v_efficient", v_efficient, ndim=2, shape=(k, k))
    index = np.arange(k) if keep is None else np.asarray(list(keep), dtype=np.int64)
    if names is not None and len(names) != index.size:
        raise SpecificationError(
            _i18n.message("diag.hausman.shape", name="names", expected=index.size, got=len(names)),
            key="diag.hausman.shape",
        )
    diff = b1[index] - b2[index]
    # Stata: V2 = V1 - V2, forced symmetric when it is not (hausman.ado, SymPosDef)
    vd = v1[np.ix_(index, index)] - v2[np.ix_(index, index)]
    vd = (vd + vd.T) / 2.0
    positive = bool(np.min(np.linalg.eigvalsh(vd)) > 0)
    inverse, rank = _ginv(vd)
    value = float(diff @ inverse @ diff)
    degrees = rank if df is None else int(df)
    notes: list[Note] = []
    if not positive:
        notes.append(Note("diag.hausman.not_positive_definite"))
    if rank != index.size:
        notes.append(Note("diag.hausman.rank", {"rank": rank, "k": index.size}))
    if value < 0:
        notes.append(Note("diag.hausman.negative", {"value": value}))
        pvalue = math.nan
    else:
        pvalue = float(stats.chi2.sf(value, degrees)) if degrees > 0 else math.nan
    details = None
    if names is not None:
        sd = np.sqrt(np.where(np.diag(vd) >= 0, np.diag(vd), np.nan))
        details = pd.DataFrame(
            np.column_stack([b1[index], b2[index], diff, sd]),
            index=pd.Index(list(names)),
            columns=["consistent", "efficient", "difference", "std_err"],
        )
    return TestResult(
        name="hausman",
        statistic=value,
        pvalue=pvalue,
        statistic_name="chi2",
        df=(float(degrees),),
        options={"rank": rank, "k": int(index.size)},
        details=details,
        notes=tuple(notes),
    )


def _ginv(a: np.ndarray) -> tuple[np.ndarray, int]:
    """Moore-Penrose inverse with Stata's tolerance ``max(singular value) * k * 1e-9``."""
    k = a.shape[0]
    if k == 0:
        return a, 0
    u, s, vt = np.linalg.svd(a)
    tolerance = float(np.max(s)) * k * 1e-9
    keep = s > tolerance
    inverse = (vt.T[:, keep] / s[keep]) @ u.T[keep, :]
    return inverse, int(keep.sum())
