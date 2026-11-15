"""Jarque-Bera against its formula computed in Stata 17 (domain ``diagnostics_normality``).

Stata has no official Jarque-Bera command.  The case computes the statistic in Stata from the
skewness, kurtosis and count of ``summarize, detail`` with the published formula
``n/6 * (S^2 + (K-3)^2/4)`` and its ``chi2tail`` p-value, independently of the package; the
residuals are exported and given to the function.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.normality.jarque_bera import jarque_bera
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_normality"

#: tag -> (dataset, rows kept, model)
SPECS = {
    "jb_heavy_tails": ("sim_het", "", "y x z1 z2"),
    "jb_normal": ("sim_iv_valid", "", "y x z1 z2"),
    "jb_small": ("sim_balanced", "in 1/15", "x z"),
}


def _stata(tag: str, rows: str, model: str) -> list[str]:
    keep = [f"quietly keep {rows}"] if rows else []
    return [
        *keep,
        f"quietly regress {model}",
        "quietly predict double ehat if e(sample), resid",
        "quietly keep if !missing(ehat)",
        export(tag, "e", 'st_data(., "ehat")'),
        "quietly summarize ehat, detail",
        "scalar jb = r(N) / 6 * (r(skewness)^2 + (r(kurtosis) - 3)^2 / 4)",
        "scalar jbp = chi2tail(2, jb)",
        export_scalars(tag, "jb", ["jb", "jbp"]),
    ]


def _register(tag: str, dataset: str, rows: str, model: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        result = jarque_bera(s["e"].ravel())
        return {"jb": np.array([result.statistic, result.pvalue]), "e": s["e"]}

    case(tag, domain=DOMAIN, dataset=dataset, stata=_stata(tag, rows, model), atol=1e-9, rtol=1e-9)(
        run
    )


for _tag, (_dataset, _rows, _model) in SPECS.items():
    _register(_tag, _dataset, _rows, _model)
