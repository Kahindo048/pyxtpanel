"""Pesaran-Yamagata slope homogeneity test against ``xthst`` in Stata 17 (domain
``diagnostics_heterogeneity``).

``xthst depvar indepvars`` returns the standard ``Delta`` and the adjusted ``Delta`` in
``r(delta)``.  Its two-sided p-values, ``2 * (1 - normal(|Delta|))``, lose every digit below 1e-15,
so the cases recompute them as ``2 * normal(-|Delta|)`` from Stata's ``Delta``.  ``xthst`` builds
the cross-sectional averages with ``egen``, whose default type is ``float``: the cases run with
``set type double`` to compare the computation rather than the storage.  Swamy's ``S`` is not
returned: it is recovered in Stata from ``Delta = sqrt(N) (S / N - K) / sqrt(2 K)``, so its value
rests on ``xthst``'s.  The cases cover a panel with common slopes, one with slopes drawn per unit,
an unbalanced panel with holes, variables partialled out, no constant, the ``ar`` correction and
cross-sectional averages with and without lags, and the HAC version (Bartlett, quadratic
spectral and truncated kernels, automatic and fixed bandwidth, whitening, partialled variables,
an unbalanced panel), whose average bandwidth is compared too.  Whitening is checked with one
and two regressors only: with three, the signs of the singular vectors that Mata's ``svd`` returns
differ from LAPACK's under SciPy for about a quarter of the units, and ``xthst`` rebuilds its
whitening matrix in a way that depends on them (see ``_hac._whiten``).
"""

from __future__ import annotations

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.heterogeneity.slope_homogeneity import slope_homogeneity
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_heterogeneity"

#: tag -> (panel, xthst options, python options)
SPECS = {
    "sh_homogeneous": ("slopes_homog", "", {}),
    "sh_heterogeneous": ("slopes_het", "", {}),
    "sh_unbalanced": ("slopes_unbalanced", "", {}),
    "sh_partial": ("slopes_het", ", partial(w)", {"partial": ["w"]}),
    "sh_noconstant": ("slopes_het", ", noconstant", {"constant": False}),
    "sh_ar": ("slopes_het", ", ar", {"ar": True}),
    "sh_csa": ("slopes_het", ", crosssectional(y x z)", {"crosssectional": ["y", "x", "z"]}),
    "sh_csa_lags": (
        "slopes_het",
        ", crosssectional(y x z, cr_lags(2))",
        {"crosssectional": ["y", "x", "z"], "cr_lags": 2},
    ),
    "sh_hac_bartlett": ("slopes_het", ", hac", {"hac": True}),
    "sh_hac_qs": ("slopes_het", ", hac kernel(qs)", {"hac": True, "kernel": "qs"}),
    "sh_hac_truncated": (
        "slopes_het",
        ", hac kernel(truncated)",
        {"hac": True, "kernel": "truncated"},
    ),
    "sh_hac_whitening": ("slopes_het", ", whitening", {"whitening": True}),
    "sh_hac_bandwidth": ("slopes_het", ", bw(3)", {"bandwidth": 3}),
    "sh_hac_qs_bandwidth": ("slopes_het", ", kernel(qs) bw(4)", {"kernel": "qs", "bandwidth": 4}),
    "sh_hac_truncated_whitening": (
        "slopes_het",
        ", whitening kernel(truncated)",
        {"whitening": True, "kernel": "truncated"},
    ),
    "sh_hac_partial": (
        "slopes_het",
        ", hac partial(w) crosssectional(y x z)",
        {"hac": True, "partial": ["w"], "crosssectional": ["y", "x", "z"]},
    ),
    "sh_hac_unbalanced": ("slopes_unbalanced", ", hac", {"hac": True}),
    "sh_hac_unbalanced_qs": (
        "slopes_unbalanced",
        ", hac kernel(qs)",
        {"hac": True, "kernel": "qs"},
    ),
    "sh_hac_unbalanced_whitening": (
        "slopes_unbalanced",
        ", whitening kernel(truncated)",
        {"whitening": True, "kernel": "truncated"},
    ),
    "sh_hac_bandwidth_zero": ("slopes_het", ", bw(0)", {"bandwidth": 0}),
    "sh_hac_csa_lags": (
        "slopes_het",
        ", hac crosssectional(y x z, cr_lags(1))",
        {"hac": True, "crosssectional": ["y", "x", "z"], "cr_lags": 1},
    ),
    "sh_hac_one_regressor": ("slopes_het", ", hac", {"hac": True}),
    "sh_hac_one_whitening": (
        "slopes_het",
        ", whitening kernel(qs)",
        {"whitening": True, "kernel": "qs"},
    ),
    "sh_hac_one_truncated": (
        "slopes_het",
        ", hac kernel(truncated)",
        {"hac": True, "kernel": "truncated"},
    ),
    "sh_hac_three_regressors": ("slopes_het", ", hac", {"hac": True}),
    "sh_hac_three_truncated": (
        "slopes_het",
        ", hac kernel(truncated)",
        {"hac": True, "kernel": "truncated"},
    ),
}
#: cases that test another set of regressors than ``x`` and ``z``
REGRESSOR_SETS = {
    "sh_hac_one_regressor": ("x",),
    "sh_hac_one_whitening": ("x",),
    "sh_hac_one_truncated": ("x",),
    "sh_hac_three_regressors": ("x", "z", "w"),
    "sh_hac_three_truncated": ("x", "z", "w"),
}
#: the cases that use the HAC version also compare the average bandwidth (``r(bw)``)
HAC = {tag for tag in SPECS if tag.startswith("sh_hac")}
REGRESSORS = ("x", "z")


def _stata(tag: str, options: str) -> list[str]:
    regressors = REGRESSOR_SETS.get(tag, REGRESSORS)
    k = len(regressors)
    return [
        "quietly xtset id t",
        "set type double",
        "quietly xtdescribe",
        "scalar nunits = r(N)",
        f"quietly xthst y {' '.join(regressors)}{options}, nooutput"
        if not options.startswith(",")
        else f"quietly xthst y {' '.join(regressors)}{options} nooutput",
        "matrix dd = r(delta)",
        "scalar d1 = dd[1, 1]",
        "scalar d2 = dd[2, 1]",
        "scalar p1 = 2 * normal(-abs(d1))",
        "scalar p2 = 2 * normal(-abs(d2))",
        f"scalar sw = nunits * ({k} + d1 * sqrt(2 * {k}) / sqrt(nunits))",
        f"scalar swp = chi2tail({k} * (nunits - 1), sw)",
        "capture scalar bw = r(bw)",
    ]


def _register(tag: str, dataset: str, options: str, python: dict[str, object]) -> None:
    def run() -> dict[str, np.ndarray]:
        panel = xtset(panels.load(dataset), "id", "t")
        result = slope_homogeneity(panel, "y", list(REGRESSOR_SETS.get(tag, REGRESSORS)), **python)  # type: ignore[arg-type]
        assert result.details is not None
        table = result.details
        values = [
            table.loc["delta", "statistic"],
            table.loc["delta_adjusted", "statistic"],
            table.loc["delta", "pvalue"],
            table.loc["delta_adjusted", "pvalue"],
            table.loc["swamy", "statistic"],
            table.loc["swamy", "pvalue"],
        ]
        if tag in HAC:
            values.append(float(result.options["bandwidth"]))  # type: ignore[arg-type]
        return {"sh": np.array(values)}

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=[
            *_stata(tag, options),
            export_scalars(
                tag, "sh", ["d1", "d2", "p1", "p2", "sw", "swp", *(["bw"] if tag in HAC else [])]
            ),
            "set type float",
        ],
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _options, _python) in SPECS.items():
    _register(_tag, _dataset, _options, _python)
