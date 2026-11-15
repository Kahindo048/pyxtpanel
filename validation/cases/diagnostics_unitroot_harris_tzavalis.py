"""Harris-Tzavalis test against ``xtunitroot ht`` in Stata 17 (domain ``diagnostics_unitroot``).

The cases cover the three specifications (constant, no constant, panel trends), the removal of
cross-sectional means and the small-sample adjustment of ``T``, on panels of random walks,
stationary series, series with a trend and a mixture.  Stata returns ``rho``, its mean and
variance under the null, ``z`` and the p-value; all are compared.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.unitroot.harris_tzavalis import harris_tzavalis
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_unitroot"

#: tag -> (dataset, Stata options, Python keyword arguments)
SPECS: dict[str, tuple[str, str, dict[str, Any]]] = {
    "ht_walk": ("ur_walk", "", {}),
    "ht_stationary": ("ur_stationary", "", {}),
    "ht_none": ("ur_walk", "noconstant", {"deterministics": "none"}),
    "ht_trend": ("ur_trend", "trend", {"deterministics": "trend"}),
    "ht_demean": ("ur_stationary", "demean", {"demean": True}),
    "ht_altt": ("ur_mixed", "altt", {"small_sample_t": True}),
}


def _stata(tag: str, options: str) -> list[str]:
    after = f", {options}" if options else ""
    return [
        "quietly xtset id t",
        f"quietly xtunitroot ht y{after}",
        "scalar htrho = r(rho)",
        "scalar htmean = r(mean_rho)",
        "scalar htvar = r(Var_rho)",
        "scalar htz = r(z)",
        "scalar htp = r(p)",
        export_scalars(tag, "ht", ["htrho", "htmean", "htvar", "htz", "htp"]),
    ]


def _register(tag: str, dataset: str, options: str, kwargs: dict[str, Any]) -> None:
    def run() -> dict[str, np.ndarray]:
        result = harris_tzavalis(xtset(panels.load(dataset), "id", "t"), "y", **kwargs)
        assert result.details is not None
        rho, mean, variance = result.details["value"]
        return {"ht": np.array([rho, mean, variance, result.statistic, result.pvalue])}

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, options),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _options, _kwargs) in SPECS.items():
    _register(_tag, _dataset, _options, _kwargs)
