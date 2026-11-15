"""Shapiro-Francia test against ``sfrancia`` in Stata 17 (domain ``diagnostics_normality``).

Stata takes the residuals of a ``regress`` (or a raw variable for the small samples), runs
``sfrancia`` and returns ``r(W)`` (the statistic W'), ``r(V)``, ``r(z)`` and ``r(p)``; the values
are exported and given to the function.  The cases cover the log transformation (n >= 10) and the
Box-Cox one (``boxcox``, n >= 5), heavy tails, normal data, a sample below the smallest size where
Stata returns missing values, and data with many tied values.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.normality.shapiro_francia import shapiro_francia
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_normality"

#: Cases run with Stata's default storage type (float): ``swilk`` and ``sfrancia`` create their
#: working variables with it, so the data and the coefficients are rounded to single precision.
#: The others run with ``set type double`` and compare the algorithm itself.
SINGLE = {"sf_single_precision"}

#: tag -> (dataset, rows kept, model or "" for the raw variable x, rounding, sfrancia options)
SPECS = {
    "sf_heavy_tails": ("sim_het", "", "y x z1 z2", None, ""),
    "sf_single_precision": ("sim_het", "", "y x z1 z2", None, ""),
    "sf_normal": ("sim_iv_valid", "", "y x z1 z2", None, ""),
    "sf_boxcox": ("sim_het", "", "y x z1 z2", None, ", boxcox"),
    "sf_boxcox_small": ("sim_balanced", "in 1/8", "", None, ", boxcox"),
    "sf_small": ("sim_balanced", "in 1/12", "", None, ""),
    "sf_ties": ("sim_het", "", "y x z1 z2", 0.5, ""),
    "sf_noties": ("sim_het", "", "y x z1 z2", 0.5, ", noties"),
}


def _stata(tag: str, rows: str, model: str, rounding: float | None, options: str) -> list[str]:
    keep = [f"quietly keep {rows}"] if rows else []
    if model:
        fit = [f"quietly regress {model}", "quietly predict double ehat if e(sample), resid"]
    else:
        fit = ["quietly generate double ehat = x"]
    tie = [f"quietly replace ehat = round(ehat, {rounding})"] if rounding else []
    return [
        *keep,
        *fit,
        *tie,
        "quietly keep if !missing(ehat)",
        export(tag, "e", 'st_data(., "ehat")'),
        "set type double" if tag not in SINGLE else "set type float",
        f"quietly sfrancia ehat{options}",
        "scalar fw = r(W)",
        "scalar fv = r(V)",
        "scalar fz = r(z)",
        "scalar fp = r(p)",
        "scalar fn = r(N)",
        export_scalars(tag, "sf", ["fw", "fv", "fz", "fp", "fn"]),
        "set type float",
    ]


def _register(
    tag: str, dataset: str, rows: str, model: str, rounding: float | None, options: str
) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        result = shapiro_francia(
            s["e"].ravel(),
            ties="noties" not in options,
            transform="boxcox" if "boxcox" in options else "log",
        )
        assert result.details is not None
        w, v, z = result.details["value"]
        return {
            "sf": np.array([w, v, z, result.pvalue, result.n_obs]),
            "e": s["e"],
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=_stata(tag, rows, model, rounding, options),
        atol=1e-8 if tag in SINGLE else 1e-9,
        rtol=1e-5 if tag in SINGLE else 1e-9,
    )(run)


for _tag, _spec in SPECS.items():
    _register(_tag, *_spec)
