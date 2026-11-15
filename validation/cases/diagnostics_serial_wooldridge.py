"""Wooldridge autocorrelation test against ``xtserial`` in Stata 17 (domain ``diagnostics_serial``).

``xtserial depvar indepvars`` regresses the first differences without a constant, then the
residuals on their lag, and tests the slope against -1/2 with a covariance clustered by unit.  The
cases cover errors with first-order autocorrelation, independent errors, and panels with holes and
with entry and exit (where Stata's ``D.`` and ``L.`` follow the calendar).
"""

from __future__ import annotations

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.serial.wooldridge import wooldridge
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_serial"

#: tag -> (dataset, dependent, regressors)
SPECS = {
    "wr_autocorrelated": ("serial", "y", ("x", "z")),
    "wr_independent": ("balanced", "x", ("z",)),
    "wr_gaps": ("gaps", "x", ("z",)),
    "wr_entry_exit": ("entry_exit", "x", ("z",)),
}


def _stata(tag: str, dependent: str, regressors: tuple[str, ...]) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly xtserial {dependent} {' '.join(regressors)}",
        "scalar ff = r(F)",
        "scalar fdf = r(df)",
        "scalar fdfr = r(df_r)",
        "scalar fp = r(p)",
        "scalar frho = r(corr)",
        export_scalars(tag, "w", ["ff", "fdf", "fdfr", "fp", "frho"]),
    ]


def _register(tag: str, dataset: str, dependent: str, regressors: tuple[str, ...]) -> None:
    def run() -> dict[str, np.ndarray]:
        result = wooldridge(xtset(panels.load(dataset), "id", "t"), dependent, regressors)
        assert result.details is not None
        return {
            "w": np.array(
                [
                    result.statistic,
                    result.df[0],
                    result.df[1],
                    result.pvalue,
                    result.details.loc["residual_autocorrelation", "value"],
                ]
            )
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, dependent, regressors),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _dependent, _regressors) in SPECS.items():
    _register(_tag, _dataset, _dependent, _regressors)
