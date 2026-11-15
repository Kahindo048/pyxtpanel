"""Breusch-Pagan LM test against an independent Mata implementation (domain ``diagnostics_crosssection``).

``xttest2`` (Baum, SSC) is not installed.  Stata fits ``xtreg y x, fe``, keeps the within
residuals, reshapes them wide (one column per unit) and a Mata loop computes, for every pair of
columns, the correlation over the rows where both are observed and accumulates
``LM = sum T_ij * rho_ij^2``; the p-value is ``chi2tail(N (N - 1) / 2, LM)``.  The code shares
nothing with the package.
"""

from __future__ import annotations

from pyxtpanel.diagnostics.crosssection.lm_cross_section import lm_cross_section
from validation.cases._residuals import register_residual_case

DOMAIN = "diagnostics_crosssection"

MATA_LM = [
    "preserve",
    "quietly keep if !missing(ehat)",
    "quietly keep id t ehat",
    "quietly reshape wide ehat, i(t) j(id)",
    "mata:",
    "X = st_data(., 2..st_nvar())",
    "n = cols(X)",
    "lm = 0",
    "for (i = 1; i < n; i++) {",
    "    for (j = i + 1; j <= n; j++) {",
    "        ok = (X[., i] :< .) :& (X[., j] :< .)",
    "        a = select(X[., i], ok)",
    "        b = select(X[., j], ok)",
    "        r = quadcorrelation((a, b))[1, 2]",
    "        lm = lm + rows(a) * r^2",
    "    }",
    "}",
    'st_numscalar("lmv", lm)',
    'st_numscalar("lmdf", n * (n - 1) / 2)',
    "end",
    "restore",
    "scalar lmp = chi2tail(lmdf, lmv)",
]


def _python(panel: object) -> list[float]:
    result = lm_cross_section(panel, "ehat")  # type: ignore[arg-type]
    return [result.statistic, result.df[0], result.pvalue]


for _tag, _dataset, _model in (
    ("lmcs_factor", "csd", "y x"),
    ("lmcs_unbalanced", "csd_unbalanced", "y x"),
    ("lmcs_independent", "balanced", "x"),
):
    register_residual_case(
        _tag,
        domain=DOMAIN,
        dataset=_dataset,
        model=_model,
        stata_lines=MATA_LM,
        scalars=["lmv", "lmdf", "lmp"],
        python=_python,
    )
