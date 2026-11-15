"""Shapiro-Wilk test against ``swilk`` in Stata 17 (domain ``diagnostics_normality``).

Stata takes the residuals of a ``regress`` (or a raw variable for the very small samples), runs
``swilk`` and returns ``r(W)``, ``r(V)``, ``r(z)`` and ``r(p)``; the values are exported and given
to the function.  The cases cover each branch of Royston's approximation (n = 3, exact; n = 4 and
5, coefficients with one edge value; n = 7 and 11, the small-sample moments; n = 12 and more),
heavy tails, normal data, and data with many tied values with and without ``noties``.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.normality.shapiro_wilk import shapiro_wilk
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_normality"

#: Cases run with Stata's default storage type (float): ``swilk`` and ``sfrancia`` create their
#: working variables with it, so the data and the coefficients are rounded to single precision.
#: The others run with ``set type double`` and compare the algorithm itself.
SINGLE = {"sw_single_precision"}

#: tag -> (dataset, rows kept, model or "" for the raw variable x, round the values to, swilk options)
SPECS = {
    "sw_heavy_tails": ("sim_het", "", "y x z1 z2", None, ""),
    "sw_single_precision": ("sim_het", "", "y x z1 z2", None, ""),
    "sw_normal": ("sim_iv_valid", "", "y x z1 z2", None, ""),
    "sw_small": ("sim_balanced", "in 1/15", "x z", None, ""),
    "sw_n3": ("sim_balanced", "in 1/3", "", None, ""),
    "sw_n4": ("sim_balanced", "in 1/4", "", None, ""),
    "sw_n5": ("sim_balanced", "in 1/5", "", None, ""),
    "sw_n7": ("sim_balanced", "in 1/7", "", None, ""),
    "sw_n11": ("sim_balanced", "in 1/11", "", None, ""),
    "sw_n12": ("sim_balanced", "in 1/12", "", None, ""),
    "sw_ties": ("sim_het", "", "y x z1 z2", 0.5, ""),
    "sw_noties": ("sim_het", "", "y x z1 z2", 0.5, ", noties"),
    "sw_large": ("sim_wide", "", "", None, ""),
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
        f"quietly swilk ehat{options}",
        "scalar sw = r(W)",
        "scalar sv = r(V)",
        "scalar sz = r(z)",
        "scalar sp = r(p)",
        "scalar sn = r(N)",
        export_scalars(tag, "sw", ["sw", "sv", "sz", "sp", "sn"]),
        "set type float",
        "set type float",
    ]


def _register(
    tag: str, dataset: str, rows: str, model: str, rounding: float | None, options: str
) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        result = shapiro_wilk(s["e"].ravel(), ties="noties" not in options)
        assert result.details is not None
        w, v, z = result.details["value"]
        return {
            "sw": np.array([w, v, z, result.pvalue, result.n_obs]),
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
