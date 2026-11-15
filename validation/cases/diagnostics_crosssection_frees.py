"""Frees' test against ``xtcsd, frees`` in Stata 17 (domain ``diagnostics_crosssection``).

``xtcsd, frees`` returns the statistic in ``r(frees)``.  Below 31 common periods it prints
critical values and no p-value, so only the statistic is compared (the critical values of the
table are checked in the unit tests against the numbers Stata printed).  On the long panel (40
periods) Stata prints a p-value that it does not return; the case computes it in Stata with the
same normal approximation, ``2 * (1 - normal(frees / seq))``, from the number of periods.
"""

from __future__ import annotations

import math

from pyxtpanel.diagnostics.crosssection.frees import frees
from validation.cases._residuals import register_residual_case

DOMAIN = "diagnostics_crosssection"

NORMAL_APPROXIMATION = [
    "quietly summarize t",
    "scalar mm = r(max) - r(min) + 1",
    "scalar seq1 = (32/25) * ((mm + 2)^2 / ((mm - 1)^3 * (mm + 1)^2))",
    "scalar seq2 = (4/25) * (((5*mm + 6)^2 * (mm - 3)) / ((mm * (mm - 1)^2) * (mm + 1)^2))",
    "scalar frp = 2 * (1 - normal(fr / sqrt(seq1 + seq2)))",
]


def _short(panel: object) -> list[float]:
    return [frees(panel, "ehat").statistic]  # type: ignore[arg-type]


def _long(panel: object) -> list[float]:
    result = frees(panel, "ehat")  # type: ignore[arg-type]
    return [result.statistic, result.pvalue]


for _tag, _dataset in (("frees_balanced", "csd"), ("frees_trimmed", "csd_trimmed")):
    register_residual_case(
        _tag,
        domain=DOMAIN,
        dataset=_dataset,
        model="y x",
        stata_lines=["quietly xtcsd, frees", "scalar fr = r(frees)"],
        scalars=["fr"],
        python=_short,
    )

register_residual_case(
    "frees_long",
    domain=DOMAIN,
    dataset="csd_long",
    model="y x",
    stata_lines=["quietly xtcsd, frees", "scalar fr = r(frees)", *NORMAL_APPROXIMATION],
    scalars=["fr", "frp"],
    python=_long,
)
del math
