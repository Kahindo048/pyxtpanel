"""Correlation matrix between units against an independent Mata computation in Stata 17.

The variable is reshaped wide (one column per unit) and a Mata loop correlates every pair of
columns over the rows where both are observed (``quadcorrelation``), leaving the pair missing
when it has fewer than 3 common rows or no variation; the matrices of correlations and of common
periods and the two averages over the defined pairs are exported and compared.  The sparse panel
has several pairs too short to be defined.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.crosssection.cross_correlation import cross_correlation
from validation import panels
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_crosssection"

#: tag -> (dataset, variable)
SPECS = {
    "cc_balanced": ("csd", "s"),
    "cc_unbalanced": ("csd_unbalanced", "s"),
    "cc_sparse": ("sparse_t", "x"),
}


def _stata(tag: str, variable: str) -> list[str]:
    return [
        "quietly xtset id t",
        "preserve",
        f"quietly keep if !missing({variable})",
        f"quietly keep id t {variable}",
        f"quietly reshape wide {variable}, i(t) j(id)",
        "mata:",
        "X = st_data(., 2..st_nvar())",
        "n = cols(X)",
        "R = J(n, n, .)",
        "N = J(n, n, .)",
        "total = 0",
        "totalabs = 0",
        "count = 0",
        "for (i = 1; i <= n; i++) {",
        "    for (j = i; j <= n; j++) {",
        "        ok = (X[., i] :< .) :& (X[., j] :< .)",
        "        a = select(X[., i], ok)",
        "        b = select(X[., j], ok)",
        "        N[i, j] = N[j, i] = rows(a)",
        "        if (i == j) {",
        "            R[i, j] = 1",
        "        }",
        "        else if (rows(a) >= 3) {",
        "            r = quadcorrelation((a, b))[1, 2]",
        "            R[i, j] = R[j, i] = r",
        "            if (r < .) {",
        "                total = total + r",
        "                totalabs = totalabs + abs(r)",
        "                count = count + 1",
        "            }",
        "        }",
        "    }",
        "}",
        'st_numscalar("ccavg", total / count)',
        'st_numscalar("ccabs", totalabs / count)',
        'st_numscalar("ccn", count)',
        "end",
        export(tag, "rho", "R"),
        export(tag, "overlap", "N"),
        "restore",
        export_scalars(tag, "avg", ["ccavg", "ccabs", "ccn"]),
    ]


def _register(tag: str, dataset: str, variable: str) -> None:
    def run() -> dict[str, np.ndarray]:
        panel = xtset(panels.load(dataset), "id", "t", delta=2 if dataset == "sparse_t" else 1)
        result = cross_correlation(panel, variable)
        return {
            "rho": result.correlation.to_numpy(),
            "overlap": result.overlap.to_numpy().astype(float),
            "avg": np.array([result.average, result.average_abs, result.n_pairs]),
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, variable),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _variable) in SPECS.items():
    _register(_tag, _dataset, _variable)
