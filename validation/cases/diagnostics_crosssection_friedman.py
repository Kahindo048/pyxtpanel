"""Friedman's test against ``xtcsd, friedman`` in Stata 17 (domain ``diagnostics_crosssection``).

Stata fits ``xtreg y x, fe``, takes the within residuals and runs ``xtcsd, friedman``; the
statistic is ``r(fried)`` and the p-value is the chi-square tail with ``N_g - 1`` degrees of
freedom, computed in Stata.  The trimmed panel has units of different lengths, so Stata ranks the
residuals over the periods common to every unit.
"""

from __future__ import annotations

from pyxtpanel.diagnostics.crosssection.friedman import friedman
from validation.cases._residuals import register_residual_case

DOMAIN = "diagnostics_crosssection"

for _tag, _dataset in (("fri_balanced", "csd"), ("fri_trimmed", "csd_trimmed")):
    register_residual_case(
        _tag,
        domain=DOMAIN,
        dataset=_dataset,
        model="y x",
        stata_lines=[
            "quietly xtcsd, friedman",
            "scalar fr = r(fried)",
            "scalar frdf = e(N_g) - 1",
            "scalar frp = chi2tail(frdf, fr)",
        ],
        scalars=["fr", "frdf", "frp"],
        python=lambda panel: (lambda r: [r.statistic, r.df[0], r.pvalue])(friedman(panel, "ehat")),
    )

# The same statistic with the degrees of freedom of its distribution, T - 1 (``dof="periods"``):
# xtcsd does not offer it, so Stata computes the tail itself from r(fried) and e(g_min).
for _tag, _dataset in (("fri_periods_balanced", "csd"), ("fri_periods_trimmed", "csd_trimmed")):
    register_residual_case(
        _tag,
        domain=DOMAIN,
        dataset=_dataset,
        model="y x",
        stata_lines=[
            "quietly xtcsd, friedman",
            "scalar fr = r(fried)",
            "scalar frdf = e(g_min) - 1",
            "scalar frp = chi2tail(frdf, fr)",
        ],
        scalars=["fr", "frdf", "frp"],
        python=lambda panel: (lambda r: [r.statistic, r.df[0], r.pvalue])(
            friedman(panel, "ehat", dof="periods")
        ),
    )
