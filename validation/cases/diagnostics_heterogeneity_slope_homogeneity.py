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
cross-sectional averages.
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
}
REGRESSORS = ("x", "z")


def _stata(options: str) -> list[str]:
    return [
        "quietly xtset id t",
        "set type double",
        "quietly xtdescribe",
        "scalar nunits = r(N)",
        f"quietly xthst y {' '.join(REGRESSORS)}{options}, nooutput"
        if not options.startswith(",")
        else f"quietly xthst y {' '.join(REGRESSORS)}{options} nooutput",
        "matrix dd = r(delta)",
        "scalar d1 = dd[1, 1]",
        "scalar d2 = dd[2, 1]",
        "scalar p1 = 2 * normal(-abs(d1))",
        "scalar p2 = 2 * normal(-abs(d2))",
        f"scalar sw = nunits * ({len(REGRESSORS)} + d1 * sqrt(2 * {len(REGRESSORS)}) / sqrt(nunits))",
        f"scalar swp = chi2tail({len(REGRESSORS)} * (nunits - 1), sw)",
    ]


def _register(tag: str, dataset: str, options: str, python: dict[str, object]) -> None:
    def run() -> dict[str, np.ndarray]:
        panel = xtset(panels.load(dataset), "id", "t")
        result = slope_homogeneity(panel, "y", list(REGRESSORS), **python)  # type: ignore[arg-type]
        assert result.details is not None
        table = result.details
        return {
            "sh": np.array(
                [
                    table.loc["delta", "statistic"],
                    table.loc["delta_adjusted", "statistic"],
                    table.loc["delta", "pvalue"],
                    table.loc["delta_adjusted", "pvalue"],
                    table.loc["swamy", "statistic"],
                    table.loc["swamy", "pvalue"],
                ]
            )
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=[
            *_stata(options),
            export_scalars(tag, "sh", ["d1", "d2", "p1", "p2", "sw", "swp"]),
            "set type float",
        ],
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _options, _python) in SPECS.items():
    _register(_tag, _dataset, _options, _python)
