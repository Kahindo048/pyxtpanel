"""Modified Wald test against ``xttest3`` in Stata 17 (domain ``diagnostics_heteroskedasticity``).

Stata fits ``xtreg, fe``, takes the within residuals (``predict, e``) and runs ``xttest3``; the
residuals and the identifier are exported and given to the function, which finds the unit
boundaries from the identifier.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.heteroskedasticity.wald_groupwise import wald_groupwise
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_heteroskedasticity"

#: tag -> (dataset, model)
SPECS = {
    "wg_balanced": ("sim_groupwise", "y x z"),
    "wg_unbalanced": ("sim_entry_exit", "x z"),
    # T = 6 per unit: the test rejects although the errors are homoskedastic, as Baum warns for
    # short panels; the case checks that the function reproduces that behaviour.
    "wg_short_panel": ("sim_iv_valid", "y x z1 z2"),
}


def _stata(tag: str, model: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly xtreg {model}, fe",
        "quietly predict double ehat, e",
        "quietly xttest3",
        "scalar ww = r(wald)",
        "scalar wdf = r(df)",
        "scalar wp = r(p)",
        "preserve",
        "quietly keep if !missing(ehat)",
        export(tag, "e", 'st_data(., "ehat")'),
        export(tag, "id", 'st_data(., "id")'),
        "restore",
        export_scalars(tag, "w", ["ww", "wdf", "wp"]),
    ]


def _register(tag: str, dataset: str, model: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        ids = s["id"].ravel()
        starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])
        result = wald_groupwise(s["e"].ravel(), np.append(starts, ids.size))
        return {
            "w": np.array([result.statistic, result.df[0], result.pvalue]),
            "e": s["e"],
            "id": s["id"],
        }

    case(tag, domain=DOMAIN, dataset=dataset, stata=_stata(tag, model), atol=1e-9, rtol=1e-9)(run)


for _tag, (_dataset, _model) in SPECS.items():
    _register(_tag, _dataset, _model)
