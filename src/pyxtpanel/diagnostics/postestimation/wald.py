"""Joint Wald test of linear restrictions from coefficients and covariance (category B, L.6.5).

Validated against Stata 17: ``test`` after ``regress`` and ``xtreg, fe`` with a clustered
covariance (F version) and after ``ivregress 2sls`` (chi-square version), for restrictions
equal to zero and to non-zero values; the coefficients and covariances are exported from Stata
and given to the function.  The cases are in
``validation/cases/diagnostics_postestimation_wald.py``.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, as_array, register_test

_i18n.register(
    "diag",
    {
        "diag.wald.title": {"en": "Joint Wald test", "fr": "Test de Wald conjoint"},
        "diag.wald.null": {
            "en": "the linear restrictions R b = r hold",
            "fr": "les restrictions linéaires R b = r sont vérifiées",
        },
        "diag.wald.alternative": {
            "en": "at least one restriction does not hold",
            "fr": "au moins une restriction n'est pas vérifiée",
        },
        "diag.wald.dependent": {
            "en": "the restrictions are not independent: R V R' is singular (rank {rank} for "
            "{q} restrictions). Remove the redundant restrictions.",
            "fr": "les restrictions ne sont pas indépendantes : R V R' est singulière (rang "
            "{rank} pour {q} restrictions). Retirez les restrictions redondantes.",
        },
        "diag.wald.bad_df": {
            "en": "df_resid must be a positive number, got {value}.",
            "fr": "df_resid doit être un nombre positif, reçu {value}.",
        },
    },
)


@register_test(
    "wald",
    TestSpec(
        family="postestimation",
        category="B",
        wave=1,
        reference="Wald, A. (1943), Transactions of the American Mathematical Society 54(3)",
        stata="test",
        validation=(
            "wald_regress_f",
            "wald_xtreg_fe_cluster_f",
            "wald_ivregress_chi2",
            "wald_nonzero_values",
        ),
        inputs="arrays",
        balanced="na",
        r="plm::pwaldtest",
        summary="Joint Wald chi-square or F test of linear restrictions R b = r.",
    ),
)
def wald(
    b: object,
    v: object,
    restrictions: object | None = None,
    values: object | None = None,
    *,
    df_resid: float | None = None,
) -> TestResult:
    """Return the Wald test of ``R b = r``: ``(Rb - r)' (R V R')^-1 (Rb - r)``.

    The statistic is a chi-square with ``q`` degrees of freedom (``q`` restrictions); when
    ``df_resid`` is given it is divided by ``q`` and compared with an ``F(q, df_resid)``, which is
    what Stata's ``test`` reports after ``regress`` and ``xtreg``.

    Parameters
    ----------
    b
        ``(k,)`` coefficient vector.
    v
        ``(k, k)`` covariance matrix of the coefficients (robust or not).
    restrictions
        ``(q, k)`` matrix ``R``; the identity (all coefficients are zero) by default.
    values
        ``(q,)`` vector ``r``; zeros by default.
    df_resid
        Residual degrees of freedom for the ``F`` version (``e(df_r)`` in Stata); ``None`` for the
        chi-square version.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    NumericalError
        If the restrictions are redundant (``R V R'`` is singular).
    SpecificationError
        If sizes do not match or ``df_resid`` is not positive.
    """
    coefficients = as_array("b", b, ndim=1)
    k = coefficients.size
    covariance = as_array("v", v, ndim=2, shape=(k, k))
    r_matrix = (
        np.eye(k)
        if restrictions is None
        else as_array("restrictions", restrictions, ndim=2, shape=(None, k))
    )
    q = r_matrix.shape[0]
    r_values = np.zeros(q) if values is None else as_array("values", values, ndim=1, shape=(q,))
    if df_resid is not None and not (df_resid > 0 and np.isfinite(df_resid)):
        raise SpecificationError(
            _i18n.message("diag.wald.bad_df", value=df_resid), key="diag.wald.bad_df"
        )
    # Stata: test -> chi2 = (Rb - r)'(R V R')^-1 (Rb - r), F = chi2 / q with df_r = e(df_r)
    gap = r_matrix @ coefficients - r_values
    middle = r_matrix @ ((covariance + covariance.T) / 2.0) @ r_matrix.T
    rank = int(np.linalg.matrix_rank(middle))
    if rank < q:
        raise NumericalError(
            _i18n.message("diag.wald.dependent", rank=rank, q=q), key="diag.wald.dependent"
        )
    chi2 = float(gap @ np.linalg.solve(middle, gap))
    if df_resid is None:
        return TestResult(
            name="wald",
            statistic=chi2,
            pvalue=float(stats.chi2.sf(chi2, q)),
            statistic_name="chi2",
            df=(float(q),),
            options={"restrictions": q},
        )
    f_value = chi2 / q
    return TestResult(
        name="wald",
        statistic=f_value,
        pvalue=float(stats.f.sf(f_value, q, df_resid)),
        statistic_name="F",
        df=(float(q), float(df_resid)),
        options={"restrictions": q},
    )
