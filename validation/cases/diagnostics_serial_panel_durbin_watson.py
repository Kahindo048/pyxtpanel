"""Panel Durbin-Watson and Baltagi-Wu LBI against ``xtregar, fe lbi`` in Stata 17 (domain
``diagnostics_serial``).

``xtregar depvar indepvars, fe lbi`` stores the modified Bhargava et al. Durbin-Watson statistic in
``e(d1)`` and the Baltagi-Wu LBI statistic in ``e(LBI)``, computed on the within residuals with the
calendar of the ``xtset`` declaration.  The cases cover errors with first-order autocorrelation,
independent errors, holes, entry and exit, a panel observed every second period (``delta(2)``) and
a panel with a unit seen once.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.serial.panel_durbin_watson import panel_durbin_watson
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_serial"

#: tag -> (dataset, dependent, regressors)
SPECS = {
    "pd_autocorrelated": ("serial", "y", ("x", "z")),
    "pd_independent": ("balanced", "x", ("z",)),
    "pd_gaps": ("gaps", "x", ("z",)),
    "pd_entry_exit": ("entry_exit", "x", ("z",)),
    "pd_delta2": ("delta2", "x", ("z",)),
    "pd_singleton": ("one_obs_unit", "x", ("z",)),
}


def _stata(dataset: str, dependent: str, regressors: tuple[str, ...]) -> list[str]:
    delta = panels.SPECS[dataset].delta
    return [
        f"quietly xtset id t, delta({delta})",
        f"quietly xtregar {dependent} {' '.join(regressors)}, fe lbi",
        "scalar dw = e(d1)",
        "scalar lbi = e(LBI)",
        "scalar nlbi = e(N_LBI)",
    ]


def _register(tag: str, dataset: str, dependent: str, regressors: tuple[str, ...]) -> None:
    def run() -> dict[str, np.ndarray]:
        spec = panels.SPECS[dataset]
        panel = xtset(panels.load(dataset), "id", "t", delta=spec.delta)
        result = panel_durbin_watson(panel, dependent, regressors)
        assert result.details is not None
        return {"d": np.array([*result.details["value"], result.n_obs])}

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=[
            *_stata(dataset, dependent, regressors),
            export_scalars(tag, "d", ["dw", "lbi", "nlbi"]),
        ],
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _dependent, _regressors) in SPECS.items():
    _register(_tag, _dataset, _dependent, _regressors)
