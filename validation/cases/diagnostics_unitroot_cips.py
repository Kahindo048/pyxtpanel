"""CIPS test against the user-written ``pescadf`` in Stata 17 (domain ``diagnostics_unitroot``).

``pescadf`` (Lewandowski, SSC) runs Pesaran's cross-sectionally augmented Dickey-Fuller
regressions and returns t-bar, its standardised form Z[t-bar] with the p-value, and for a balanced
panel the critical values of t-bar.  The cases cover the constant and the trend, the removal of the
initial cross-section mean (``nodemean``, which changes nothing: it is a constant shift of the
series absorbed by each regression's constant, so the Python call has no matching option), the
truncation of the unit statistics, one lag length and one per unit,
balanced and unbalanced panels, a panel driven by a common factor and one with fewer than ten
units (where the tables are read at their first column).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.unitroot.cips import cips
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_unitroot"

#: lags of each of the 40 units of ur_mixed (cycling through 0, 1, 2, 1)
PER_UNIT = [(0, 1, 2, 1)[i % 4] for i in range(40)]

#: tag -> (dataset, Stata options, Python keyword arguments, balanced)
SPECS: dict[str, tuple[str, str, dict[str, Any], bool]] = {
    "cips_lag0": ("ur_walk", "lags(0)", {"lags": 0}, True),
    "cips_lag1": ("ur_mixed", "lags(1)", {"lags": 1}, True),
    "cips_lag2_trend": (
        "ur_long",
        "lags(2) trend",
        {"lags": 2, "deterministics": "trend"},
        True,
    ),
    "cips_nodemean": ("ur_mixed", "lags(1) nodemean", {"lags": 1}, True),
    "cips_notrunc": ("ur_stationary", "lags(1) trunoff", {"lags": 1, "truncate": False}, True),
    "cips_lag3_trend_notrunc": (
        "ur_long",
        "lags(3) trend trunoff nodemean",
        {"lags": 3, "deterministics": "trend", "truncate": False},
        True,
    ),
    "cips_factor": ("cips_factor", "lags(1)", {"lags": 1}, True),
    "cips_factor_trend": (
        "cips_factor",
        "lags(2) trend",
        {"lags": 2, "deterministics": "trend"},
        True,
    ),
    "cips_unbalanced": ("ur_unbalanced", "lags(1)", {"lags": 1}, False),
    "cips_unbalanced_lags": (
        "ur_unbalanced",
        "lags(2) trend nodemean trunoff",
        {"lags": 2, "deterministics": "trend", "truncate": False},
        False,
    ),
    "cips_per_unit_lags": (
        "ur_mixed",
        "lags(" + " ".join(map(str, PER_UNIT)) + ")",
        {"lags": PER_UNIT},
        True,
    ),
    "cips_small_n": ("cips_small", "lags(1)", {"lags": 1}, True),
    "cips_short": ("ur_trend", "lags(1) trend", {"lags": 1, "deterministics": "trend"}, True),
}


def _stata(tag: str, options: str, balanced: bool) -> list[str]:
    scalars = ["tb", "zt", "pv", "nn"]
    lines = [
        "quietly xtset id t",
        f"quietly pescadf y, {options}",
        "scalar tb = r(tbar)",
        "scalar zt = r(Ztbar)",
        "scalar pv = r(pval)",
        "scalar nn = r(N)",
    ]
    if balanced:
        lines += ["scalar c1 = r(cv1)", "scalar c5 = r(cv5)", "scalar c10 = r(cv10)"]
        scalars += ["c1", "c5", "c10"]
    return [*lines, export_scalars(tag, "c", scalars)]


def _register(tag: str, dataset: str, options: str, kwargs: dict[str, Any], balanced: bool) -> None:
    def run() -> dict[str, np.ndarray]:
        result = cips(xtset(panels.load(dataset), "id", "t"), "y", **kwargs)
        out = [result.options["t_bar"], result.statistic, result.pvalue, result.n_groups]
        if balanced:
            out += [result.options["cv_1"], result.options["cv_5"], result.options["cv_10"]]
        return {"c": np.array(out, dtype=float)}

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, options, balanced),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _options, _kwargs, _balanced) in SPECS.items():
    _register(_tag, _dataset, _options, _kwargs, _balanced)
