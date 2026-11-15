"""Skewness and kurtosis test against ``sktest`` in Stata 17 (domain ``diagnostics_normality``).

Stata fits ``regress``, keeps the residuals and runs ``sktest`` (with and without Royston's
adjustment).  The residuals are exported and given to the function, which must return the
joint chi-square, its p-value and the two separate p-values.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.normality.skewness_kurtosis import skewness_kurtosis
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_normality"

#: tag -> (dataset, rows kept, model, sktest options)
SPECS = {
    "sk_heavy_tails": ("sim_het", "", "y x z1 z2", ""),
    "sk_normal": ("sim_iv_valid", "", "y x z1 z2", ""),
    "sk_unadjusted": ("sim_het", "", "y x z1 z2", ", noadjust"),
    "sk_small": ("sim_balanced", "in 1/15", "x z", ""),
}


def _stata(tag: str, rows: str, model: str, options: str) -> list[str]:
    keep = [f"quietly keep {rows}"] if rows else []
    return [
        *keep,
        f"quietly regress {model}",
        "quietly predict double ehat if e(sample), resid",
        "quietly keep if !missing(ehat)",
        export(tag, "e", 'st_data(., "ehat")'),
        f"quietly sktest ehat{options}",
        "scalar k2 = r(chi2)",
        "scalar kp = r(p_chi2)",
        "scalar ps = r(p_skew)",
        "scalar pk = r(p_kurt)",
        export_scalars(tag, "sk", ["k2", "kp", "ps", "pk"]),
    ]


def _register(tag: str, dataset: str, rows: str, model: str, options: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        result = skewness_kurtosis(s["e"].ravel(), adjust="noadjust" not in options)
        assert result.details is not None
        p_skew, p_kurt = result.details["pvalue"]
        return {
            "sk": np.array([result.statistic, result.pvalue, p_skew, p_kurt]),
            "e": s["e"],
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=_stata(tag, rows, model, options),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _rows, _model, _options) in SPECS.items():
    _register(_tag, _dataset, _rows, _model, _options)
