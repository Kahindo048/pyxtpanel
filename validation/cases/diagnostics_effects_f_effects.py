"""F test of the unit effects against ``xtreg, fe`` in Stata 17 (domain ``diagnostics_effects``).

Stata fits the pooled regression and ``xtreg, fe``, and stores ``e(F_f)`` with its numerator and
denominator degrees of freedom (``e(df_a)``, ``e(df_r)``).  The two residual sums of squares and
the sample sizes are exported and given to the function.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.effects.f_effects import f_effects
from validation.cases._frozen import stata_objects
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_effects"

#: tag -> (dataset, model)
SPECS = {
    "fe_balanced": ("sim_iv_valid", "y x z1 z2"),
    "fe_unbalanced": ("sim_entry_exit", "x z"),
    "fe_gaps": ("sim_gaps", "x z"),
}


def _stata(tag: str, model: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly regress {model}",
        "scalar rssp = e(rss)",
        f"quietly xtreg {model}, fe",
        "scalar rsse = e(rss)",
        "scalar nobs = e(N)",
        "scalar ngrp = e(N_g)",
        "scalar kreg = colsof(e(b)) - 1",
        "scalar ff = e(F_f)",
        "scalar fdf1 = e(df_a)",
        "scalar fdf2 = e(df_r)",
        "scalar fp = Ftail(e(df_a), e(df_r), e(F_f))",
        export_scalars(tag, "inputs", ["rssp", "rsse", "nobs", "ngrp", "kreg"]),
        export_scalars(tag, "f", ["ff", "fdf1", "fdf2", "fp"]),
    ]


def _register(tag: str, dataset: str, model: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        rss_pooled, rss_effects, n_obs, n_groups, k = s["inputs"].ravel()
        result = f_effects(rss_pooled, rss_effects, int(n_obs), int(n_groups), int(k))
        return {
            "f": np.array([result.statistic, *result.df, result.pvalue]),
            "inputs": s["inputs"],
        }

    case(tag, domain=DOMAIN, dataset=dataset, stata=_stata(tag, model), atol=1e-9, rtol=1e-9)(run)


for _tag, (_dataset, _model) in SPECS.items():
    _register(_tag, _dataset, _model)
