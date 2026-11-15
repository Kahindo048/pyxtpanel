"""Breusch-Pagan LM test against ``xttest0`` in Stata 17 (domain ``diagnostics_effects``).

Stata fits the pooled regression (its residuals are exported), then ``xtreg, re`` and
``xttest0``.  The function receives the pooled residuals, the unit boundaries and the
random-effects ``sigma_u``; two of the panels have no unit effect (``sigma_u`` is zero and
``xttest0`` reports a statistic of 0 and a p-value of 1).
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.effects.lm_random_effects import lm_random_effects
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_effects"

#: tag -> (dataset, model)
SPECS = {
    "lm_re_effect": ("sim_iv_valid", "y x z1 z2"),
    "lm_re_unbalanced": ("sim_entry_exit", "x z"),
    "lm_re_no_effect": ("sim_balanced", "x z"),
    "lm_re_gaps": ("sim_gaps", "x z"),
}


def _stata(tag: str, model: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly regress {model}",
        "quietly predict double ehat, resid",
        f"quietly xtreg {model}, re",
        export(tag, "sigma_u", 'st_numscalar("e(sigma_u)")'),
        "quietly xttest0",
        "scalar lm = r(lm)",
        "scalar lmp = r(p)",
        "preserve",
        "quietly keep if !missing(ehat)",
        export(tag, "e", 'st_data(., "ehat")'),
        export(tag, "id", 'st_data(., "id")'),
        "restore",
        export_scalars(tag, "lm", ["lm", "lmp"]),
    ]


def _register(tag: str, dataset: str, model: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        ids, e = s["id"].ravel(), s["e"].ravel()
        starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])
        result = lm_random_effects(
            e, np.append(starts, ids.size), variance_u=float(s["sigma_u"].ravel()[0])
        )
        return {
            "lm": np.array([result.statistic, result.pvalue]),
            "e": s["e"],
            "id": s["id"],
            "sigma_u": s["sigma_u"],
        }

    case(tag, domain=DOMAIN, dataset=dataset, stata=_stata(tag, model), atol=1e-9, rtol=1e-9)(run)


for _tag, (_dataset, _model) in SPECS.items():
    _register(_tag, _dataset, _model)
