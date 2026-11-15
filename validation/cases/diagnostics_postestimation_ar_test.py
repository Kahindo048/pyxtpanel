"""Arellano-Bond AR(1) and AR(2) against Stata 17 (domain ``diagnostics_postestimation``).

Stata fits a difference GMM (one lag of ``x``, instruments from ``x``, robust variance) with
``xtabond2`` (``e(ar1)``, ``e(ar2)``) and with ``xtabond`` followed by ``estat abond``
(``r(arm)``).  The residuals, regressors, instruments, weight matrix and robust covariance come
from the engine of ``pyxtabond2``, frozen in ``tests/golden/ar_pyxtabond2.json`` by
``python -m validation.cross_pyxtabond2_ar``; the package test takes them as input.  The weight
matrix of ``pyxtabond2`` is a pseudo-inverse, so the agreement is at 1e-10, not at machine
precision (tolerance 1e-8).  Both Stata commands give the same numbers on these fits.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import numpy as np

from pyxtpanel.diagnostics import Equations
from pyxtpanel.diagnostics.postestimation.ar_test import ar_test
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_postestimation"
FIXTURE = Path(__file__).resolve().parents[2] / "tests" / "golden" / "ar_pyxtabond2.json"

#: fixture tag -> (panel, two-step)
FITS = {
    "ab1_balanced": ("balanced", False),
    "ab2_balanced": ("balanced", True),
    "ab2_gaps": ("gaps", True),
    "ab2_entry_exit": ("entry_exit", True),
}


def _statistics(fit: str) -> list[tuple[float, float]]:
    data: dict[str, dict[str, object]] = json.loads(FIXTURE.read_text(encoding="utf-8"))
    e = data[fit]
    system = Equations(
        np.array(e["residuals"]),
        np.array(e["regressors"]),
        np.array(e["offsets"]),
        np.array(e["time_index"]),
        np.array(e["instruments"]),
    )
    weight, covariance = np.array(e["weight"]), np.array(e["covariance"])
    out = []
    for order in (1, 2):
        result = ar_test(order, system, weight=weight, covariance=covariance)
        out.append((result.statistic, result.pvalue))
    return out


def _xtabond2(fit: str) -> str:
    _, twostep = FITS[fit]
    options = ["noleveleq", "noconstant", "nodiffsargan", "robust"] + (
        ["twostep"] if twostep else []
    )
    return f"quietly xtabond2 x L.x, gmm(x, lag(2 .)) {' '.join(options)}"


def _xtabond(fit: str) -> str:
    _, twostep = FITS[fit]
    options = ["noconstant", "vce(robust)"] + (["twostep"] if twostep else [])
    return f"quietly xtabond x, lags(1) {' '.join(options)}"


def _from_xtabond2(fit: str) -> Callable[[], dict[str, np.ndarray]]:
    def run() -> dict[str, np.ndarray]:
        (z1, p1), (z2, p2) = _statistics(fit)
        return {"ar": np.array([z1, p1, z2, p2])}

    return run


def _from_abond(fit: str) -> Callable[[], dict[str, np.ndarray]]:
    def run() -> dict[str, np.ndarray]:
        (z1, p1), (z2, p2) = _statistics(fit)
        return {"arm": np.array([[1.0, z1, p1], [2.0, z2, p2]])}

    return run


for _fit, (_panel, _twostep) in FITS.items():
    case(
        f"ar_{_fit}",
        domain=DOMAIN,
        dataset=f"sim_{_panel}",
        stata=[
            "quietly xtset id t",
            _xtabond2(_fit),
            "scalar a1 = e(ar1)",
            "scalar p1 = e(ar1p)",
            "scalar a2 = e(ar2)",
            "scalar p2 = e(ar2p)",
            export_scalars(f"ar_{_fit}", "ar", ["a1", "p1", "a2", "p2"]),
        ],
        atol=1e-8,
        rtol=1e-8,
    )(_from_xtabond2(_fit))

# ``estat abond`` is checked on the balanced panel only: on a panel with gaps ``xtabond`` and
# ``xtabond2`` give different statistics in Stata itself (AR(1) -3.049 against -2.963 here), so
# the equations of ``pyxtabond2``, which follow ``xtabond2``, are not those of ``xtabond``.
for _fit in ("ab2_balanced", "ab1_balanced"):
    _tag = {"ab2_balanced": "ar_abond_balanced", "ab1_balanced": "ar_abond1_balanced"}[_fit]
    case(
        _tag,
        domain=DOMAIN,
        dataset=f"sim_{FITS[_fit][0]}",
        stata=[
            "quietly xtset id t",
            _xtabond(_fit),
            "quietly estat abond, artests(2)",
            export(_tag, "arm", 'st_matrix("r(arm)")'),
        ],
        atol=1e-8,
        rtol=1e-8,
    )(_from_abond(_fit))
