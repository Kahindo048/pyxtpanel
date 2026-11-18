"""Dumitrescu-Hurlin test against an independent computation in Stata 17.

``xtgcause`` (Lopez and Weber, SSC) is not installed.  For every unit Stata fits
``regress y L(1/K).y L(1/K).x`` and takes the Wald statistic of the lags of ``x`` from ``test``
(``K`` times the F statistic, which uses the unrestricted variance estimate with its degrees of
freedom); a Mata block then forms ``W-bar``, ``Z-bar`` and ``Z-tilde`` with their two-sided
p-values.  The package is also checked against ``plm::pgrangertest`` of R in its unit test.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.causality.dumitrescu_hurlin import dumitrescu_hurlin
from validation import panels
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_causality"

#: tag -> (dataset, lags)
SPECS = {
    "dh_lag1": ("granger_yes", 1),
    "dh_lag2": ("granger_yes", 2),
    "dh_lag3": ("granger_yes", 3),
    "dh_no_causality": ("granger_no", 1),
    "dh_no_causality_lag2": ("granger_no", 2),
}


def _stata(tag: str, k: int) -> list[str]:
    return [
        "quietly xtset id t",
        "quietly summarize t",
        "scalar tt = r(max) - r(min) + 1",
        "quietly levelsof id",
        "scalar nn = r(r)",
        "matrix W = J(nn, 1, .)",
        "local lx",
        f"forvalues j = 1/{k} {{",
        "    local lx `lx' L`j'.x",
        "}",
        "forvalues i = 1/`=nn' {",
        f"    quietly regress y L(1/{k}).y L(1/{k}).x if id == `i'",
        "    quietly test `lx'",
        f"    matrix W[`i', 1] = r(F) * {k}",
        "}",
        "mata:",
        'W = st_matrix("W")',
        "N = rows(W)",
        f"K = {k}",
        'T = st_numscalar("tt")',
        "wbar = mean(W)",
        "zbar = sqrt(N / (2 * K)) * (wbar - K)",
        "zt = sqrt(N / (2 * K) * (T - 3 * K - 5) / (T - 2 * K - 3)) "
        "* ((T - 3 * K - 3) / (T - 3 * K - 1) * wbar - K)",
        'st_numscalar("wbar", wbar)',
        'st_numscalar("zbar", zbar)',
        'st_numscalar("pzbar", 2 * normal(-abs(zbar)))',
        'st_numscalar("zt", zt)',
        'st_numscalar("pzt", 2 * normal(-abs(zt)))',
        "end",
        export(tag, "unit_w", 'st_matrix("W")'),
        export_scalars(tag, "agg", ["wbar", "zbar", "pzbar", "zt", "pzt", "nn", "tt"]),
    ]


def _register(tag: str, dataset: str, k: int) -> None:
    def run() -> dict[str, np.ndarray]:
        panel = xtset(panels.load(dataset), "id", "t")
        ztilde = dumitrescu_hurlin(panel, "y", "x", lags=k)
        zbar = dumitrescu_hurlin(panel, "y", "x", lags=k, statistic="zbar")
        assert ztilde.details is not None and ztilde.n_groups is not None
        return {
            "unit_w": ztilde.details[["wald"]].to_numpy(),
            "agg": np.array(
                [
                    ztilde.options["wbar"],
                    zbar.statistic,
                    zbar.pvalue,
                    ztilde.statistic,
                    ztilde.pvalue,
                    ztilde.n_groups,
                    ztilde.n_obs // ztilde.n_groups + k,  # type: ignore[operator]
                ],
                dtype=float,
            ),
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, k),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _k) in SPECS.items():
    _register(_tag, _dataset, _k)
