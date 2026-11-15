"""Difference-in-Hansen against ``xtabond2`` in Stata 17 (domain ``diagnostics_postestimation``).

Stata fits a two-step robust difference GMM with several instrument groups and stores, in
``e(diffsargan)``, for each group the Hansen statistic without the group, the difference, the
number of instruments dropped and the p-values.  The two Hansen statistics and the number of
instruments dropped are exported and given to the function, which must return Stata's
difference, degrees of freedom and p-value.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.postestimation.diff_hansen import diff_hansen
from validation.cases._frozen import stata_objects
from validation.registry import case, export

DOMAIN = "diagnostics_postestimation"
ROWS = "(2 \\ 3 \\ 5)"  # difference, instruments dropped, p-value of the difference


def _stata(tag: str, groups: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly xtabond2 x L.x z, {groups} noleveleq noconstant twostep robust",
        export(tag, "jfull", 'st_numscalar("e(hansen)")'),
        export(tag, "jexcl", 'st_matrix("e(diffsargan)")[1,.]'),
        export(tag, "stata", f'st_matrix("e(diffsargan)")[{ROWS},.]'),
    ]


def _register(tag: str, dataset: str, groups: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        j_full = float(s["jfull"].ravel()[0])
        j_excl, stata = s["jexcl"].ravel(), s["stata"]
        results = [diff_hansen(j_full, float(j), int(stata[1, g])) for g, j in enumerate(j_excl)]
        table = np.array([[r.statistic, r.df[0], r.pvalue] for r in results]).T
        return {"jfull": s["jfull"], "jexcl": s["jexcl"], "stata": table}

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=_stata(tag, groups),
        atol=1e-9,
        rtol=1e-9,
    )(run)


_register("diffhansen_balanced", "sim_balanced", "gmm(x, lag(2 .)) iv(z)")
_register("diffhansen_gaps_groups", "sim_gaps", "gmm(x, lag(2 .)) gmm(z, lag(1 2)) iv(z)")
