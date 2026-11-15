"""Hansen and Sargan against ``xtabond2`` in Stata 17 (domain ``diagnostics_postestimation``).

Stata fits a difference GMM (one lag of ``x``, instruments from ``x``) and stores
``e(sargan)``, ``e(hansen)`` and their p-values.  The moments, the weight matrices and the
scale of Sargan come from the engine of ``pyxtabond2``, frozen in
``tests/golden/overid_pyxtabond2.json`` by ``python -m validation.cross_pyxtabond2_overid``;
the package tests take them as input.  The weight matrix of ``pyxtabond2`` is a pseudo-inverse,
so the agreement is at 1e-10 and not at machine precision (tolerance 1e-8).  A disagreement with Stata would therefore show either a
fault of the tests or a difference between ``pyxtabond2`` and ``xtabond2``.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import numpy as np

from pyxtpanel.diagnostics.postestimation.hansen import hansen
from pyxtpanel.diagnostics.postestimation.sargan import sargan
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_postestimation"
FIXTURE = Path(__file__).resolve().parents[2] / "tests" / "golden" / "overid_pyxtabond2.json"

#: tag -> (panel, two-step, orthogonal deviations)
FITS = {
    "ab1_balanced": ("balanced", False, False),
    "ab2_balanced": ("balanced", True, False),
    "ab2_gaps": ("gaps", True, False),
    "ab2_entry_exit": ("entry_exit", True, False),
}


def _fixture(tag: str) -> dict[str, object]:
    data: dict[str, dict[str, object]] = json.loads(FIXTURE.read_text(encoding="utf-8"))
    return data[tag]


def _xtabond2(tag: str) -> str:
    _, twostep, orthogonal = FITS[tag]
    lags = "lag(1 .)" if orthogonal else "lag(2 .)"
    options = ["noleveleq", "noconstant", "nodiffsargan"]
    options += ["twostep"] if twostep else []
    options += ["orthogonal"] if orthogonal else []
    return f"quietly xtabond2 x L.x, gmm(x, {lags}) {' '.join(options)}"


def _stata(tag: str, kind: str) -> list[str]:
    name = f"{kind}_{tag}"
    if kind == "hansen":
        body = [
            "scalar st = e(hansen)",
            "scalar sdf = e(hansen_df)",
            "scalar sp = e(hansenp)",
        ]
    else:
        body = [
            "scalar st = e(sargan)",
            "scalar sdf = e(j) - colsof(e(b))",
            "scalar sp = e(sarganp)",
        ]
    return [
        "quietly xtset id t",
        _xtabond2(tag),
        *body,
        export_scalars(name, "s", ["st", "sdf", "sp"]),
    ]


def _hansen(tag: str) -> Callable[[], dict[str, np.ndarray]]:
    def run() -> dict[str, np.ndarray]:
        e = _fixture(tag)
        result = hansen(
            e["hansen_moments"],
            e["hansen_weight"],
            e["n_parameters"],  # type: ignore[arg-type]
            n_instruments=e["n_instruments"],  # type: ignore[arg-type]
            n_groups=e["n_groups"],  # type: ignore[arg-type]
        )
        return {"s": np.array([result.statistic, result.df[0], result.pvalue])}

    return run


def _sargan(tag: str) -> Callable[[], dict[str, np.ndarray]]:
    def run() -> dict[str, np.ndarray]:
        e = _fixture(tag)
        result = sargan(
            e["sargan_moments"],
            e["sargan_weight"],
            e["n_parameters"],  # type: ignore[arg-type]
            scale=e["sargan_scale"],  # type: ignore[arg-type]
            n_instruments=e["n_instruments"],  # type: ignore[arg-type]
            n_groups=e["n_groups"],  # type: ignore[arg-type]
        )
        return {"s": np.array([result.statistic, result.df[0], result.pvalue])}

    return run


for _tag, (_panel, _twostep, _orth) in FITS.items():
    case(
        f"sargan_{_tag}",
        domain=DOMAIN,
        dataset=f"sim_{_panel}",
        stata=_stata(_tag, "sargan"),
        atol=1e-8,
        rtol=1e-8,
    )(_sargan(_tag))
    if _twostep:
        case(
            f"hansen_{_tag}",
            domain=DOMAIN,
            dataset=f"sim_{_panel}",
            stata=_stata(_tag, "hansen"),
            atol=1e-8,
            rtol=1e-8,
        )(_hansen(_tag))
