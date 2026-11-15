"""White's test against ``estat imtest, white`` in Stata 17 (domain
``diagnostics_heteroskedasticity``).

Stata runs ``regress`` then ``estat imtest, white double`` (the ``double`` option keeps the
products in double precision; without it they are floats and the statistic differs in the
seventh digit).  The residuals and regressors are exported and given to the function.  The
four-regressor model contains a dummy whose square is collinear with it, which Stata drops from
the auxiliary regression and the function must drop too.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.heteroskedasticity.white import white
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_heteroskedasticity"

#: tag -> (dataset, regressors of the model)
SPECS = {
    "white_one": ("sim_het", ("x",)),
    "white_three": ("sim_het", ("x", "z1", "z2")),
    "white_four": ("sim_het", ("x", "z1", "z2", "d")),
    "white_homoskedastic": ("sim_iv_valid", ("x", "z1", "z2")),
}


def _stata(tag: str, regressors: tuple[str, ...]) -> list[str]:
    columns = ", ".join(f'"{v}"' for v in regressors)
    names = ", ".join(regressors)
    return [
        "quietly xtset id t",
        f"quietly regress y {names.replace(',', '')}",
        "quietly predict double ehat, resid",
        export(tag, "e", 'st_data(., "ehat")'),
        export(tag, "x", f"st_data(., ({columns}))"),
        "quietly estat imtest, white double",
        "scalar w0 = r(chi2)",
        "scalar w1 = r(df)",
        "scalar w2 = r(p)",
        export_scalars(tag, "w", ["w0", "w1", "w2"]),
    ]


def _register(tag: str, dataset: str, regressors: tuple[str, ...]) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        result = white(s["e"].ravel(), s["x"])
        return {
            "w": np.array([result.statistic, result.df[0], result.pvalue]),
            "e": s["e"],
            "x": s["x"],
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=_stata(tag, regressors),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _regressors) in SPECS.items():
    _register(_tag, _dataset, _regressors)
