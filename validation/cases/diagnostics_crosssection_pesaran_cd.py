"""Pesaran CD test against Stata 17 (domain ``diagnostics_crosssection``).

Two references, because Stata has two commands:

* ``xtcsd, pesaran`` after ``xtreg, fe`` takes the within residuals; every pair of units uses its
  common periods, with ``sqrt(T_ij)`` weights, which is the formula of Pesaran (2015) for an
  unbalanced panel.  The residuals, the identifier and the period are exported from Stata and
  rebuilt into a panel on the Python side, so that the function reads exactly Stata's residuals.
* ``xtcd2, noest`` on a variable of a strongly balanced panel (``xtcd2`` is the user-written
  command of Ditzen).  On an unbalanced panel ``xtcd2`` uses another estimator (it fills the panel
  and takes the correlation of unit means and of all non-missing values, Chudik and Pesaran
  2013), which gives other numbers than the formula above; no case compares them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pyxtpanel import xtset
from pyxtpanel.diagnostics.crosssection.pesaran_cd import pesaran_cd
from validation import panels
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_crosssection"

#: tag -> (dataset, xtreg model) for the residual cases
RESIDUALS = {
    "cd_balanced": ("csd", "y x"),
    "cd_unbalanced": ("csd_unbalanced", "y x"),
    "cd_gaps": ("gaps", "x"),
}

#: tag -> (dataset, variable) for the cases on a variable
VARIABLES = {
    "cd_variable": ("csd", "s"),
    "cd_independent": ("balanced", "x"),
}


def _residual_stata(tag: str, model: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly xtreg {model}, fe",
        "quietly predict double ehat, e",
        "quietly xtcsd, pesaran",
        "scalar cd = r(pesaran)",
        "scalar cdp = 2 * normal(-abs(cd))",
        "preserve",
        "quietly keep if !missing(ehat)",
        export(tag, "e", 'st_data(., "ehat")'),
        export(tag, "id", 'st_data(., "id")'),
        export(tag, "t", 'st_data(., "t")'),
        "restore",
        export_scalars(tag, "cd", ["cd", "cdp"]),
    ]


def _register_residuals(tag: str, dataset: str, model: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        frame = pd.DataFrame({"id": s["id"].ravel(), "t": s["t"].ravel(), "ehat": s["e"].ravel()})
        result = pesaran_cd(xtset(frame, "id", "t"), "ehat")
        return {
            "cd": np.array([result.statistic, result.pvalue]),
            "e": s["e"],
            "id": s["id"],
            "t": s["t"],
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_residual_stata(tag, model),
        atol=1e-9,
        rtol=1e-9,
    )(run)


def _variable_stata(tag: str, variable: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly xtcd2 {variable}, noest",
        "matrix cdm = r(CD)",
        "scalar cd = cdm[1, 1]",
        "scalar cdp = 2 * normal(-abs(cd))",
        export_scalars(tag, "cd", ["cd", "cdp"]),
    ]


def _register_variable(tag: str, dataset: str, variable: str) -> None:
    def run() -> dict[str, np.ndarray]:
        result = pesaran_cd(xtset(panels.load(dataset), "id", "t"), variable)
        return {"cd": np.array([result.statistic, result.pvalue])}

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_variable_stata(tag, variable),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _model) in RESIDUALS.items():
    _register_residuals(_tag, _dataset, _model)
for _tag, (_dataset, _variable) in VARIABLES.items():
    _register_variable(_tag, _dataset, _variable)
