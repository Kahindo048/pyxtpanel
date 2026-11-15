"""Breusch-Pagan / Cook-Weisberg against ``estat hettest`` in Stata 17 (domain
``diagnostics_heteroskedasticity``).

Stata runs ``regress y x z1 z2`` then ``estat hettest``; the residuals, fitted values and
regressors are exported and given to the function, which must return Stata's statistic,
degrees of freedom and p-value for the variables tested.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.heteroskedasticity.breusch_pagan import breusch_pagan
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_heteroskedasticity"
REGRESSORS = ("x", "z1", "z2")

#: tag -> (dataset, estat options, columns tested, variant); columns: "fit", "rhs" or names
SPECS = {
    "bp_fitted": ("sim_het", "", "fit", "normal"),
    "bp_rhs": ("sim_het", ", rhs", "rhs", "normal"),
    "bp_variables": ("sim_het", " z1 z2", ("z1", "z2"), "normal"),
    "bp_iid": ("sim_het", ", iid", "fit", "iid"),
    "bp_fstat": ("sim_het", ", fstat", "fit", "fstat"),
    "bp_homoskedastic": ("sim_iv_valid", "", "fit", "normal"),
}


def _stata(tag: str, options: str, variant: str) -> list[str]:
    fstat = variant == "fstat"
    names = ["r(F)", "r(df_m)", "r(df_r)", "r(p)"] if fstat else ["r(chi2)", "r(df)", "r(p)"]
    scalars = [f"bp{i}" for i in range(len(names))]
    return [
        "quietly xtset id t",
        "quietly regress y x z1 z2",
        "quietly predict double ehat, resid",
        "quietly predict double fit, xb",
        export(tag, "e", 'st_data(., "ehat")'),
        export(tag, "fit", 'st_data(., "fit")'),
        export(tag, "x", 'st_data(., ("x", "z1", "z2"))'),
        f"quietly estat hettest{options}",
        *(f"scalar {s} = {n}" for s, n in zip(scalars, names, strict=True)),
        export_scalars(tag, "bp", scalars),
    ]


def _register(tag: str, dataset: str, options: str, tested: object, variant: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        if tested == "fit":
            z = s["fit"]
        elif tested == "rhs":
            z = s["x"]
        else:
            z = s["x"][:, [REGRESSORS.index(v) for v in tested]]  # type: ignore[union-attr]
        result = breusch_pagan(s["e"].ravel(), z, variant=variant)  # type: ignore[arg-type]
        stat = [result.statistic, *result.df, result.pvalue]
        return {"bp": np.array(stat), "e": s["e"], "fit": s["fit"], "x": s["x"]}

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=_stata(tag, options, variant),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _options, _tested, _variant) in SPECS.items():
    _register(_tag, _dataset, _options, _tested, _variant)
