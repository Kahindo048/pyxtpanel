"""Variance inflation factors against ``estat vif`` in Stata 17 (domain
``diagnostics_specification``).

Stata fits ``regress`` and runs ``estat vif`` (``nosort`` keeps the order of the regressors, so
the vector compares position by position; the sorted order is a sort of the same numbers) and
``estat vif, uncentered``; the regressors are exported and given to the function, which must
return the factors of Stata.  ``estat vif`` keeps the factors in a ``float`` variable, so the
cases run with ``set type double`` to compare the computation rather than the storage.  The regressors are correlated by construction.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.specification.vif import vif
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_specification"

#: tag -> (dataset, regressors, options)
SPECS = {
    "vif_correlated": ("sim_het", "x w1 w2", ""),
    "vif_dummy": ("sim_het", "x w1 d", ""),
    "vif_uncentered": ("sim_het", "x w1 w2", ", uncentered"),
    "vif_independent": ("sim_iv_valid", "x z1 z2 z3", ""),
}


def _stata(tag: str, regressors: str, options: str) -> list[str]:
    names = regressors.split()
    count = len(names) + (1 if "uncentered" in options else 0)
    scalars = [f"v{i}" for i in range(1, count + 1)]
    quoted = ", ".join(f'"{v}"' for v in names)
    return [
        "quietly generate double w1 = 0.8 * x + 0.6 * z1" if "w1" in regressors else "",
        "quietly generate double w2 = 0.5 * x - 0.5 * w1 + 0.3 * z2" if "w2" in regressors else "",
        f"quietly regress y {regressors}",
        export(tag, "x", f"st_data(., ({quoted}))"),
        "set type double",
        f"quietly estat vif{options or ','} nosort",
        "set type float",
        *(f"scalar {s} = r(vif_{i})" for i, s in enumerate(scalars, start=1)),
        export_scalars(tag, "v", scalars),
    ]


def _register(tag: str, dataset: str, regressors: str, options: str) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        result = vif(s["x"], uncentered="uncentered" in options, sort=False)
        return {"v": result.vif.to_numpy(), "x": s["x"]}

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=[line for line in _stata(tag, regressors, options) if line],
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _regressors, _options) in SPECS.items():
    _register(_tag, _dataset, _regressors, _options)
