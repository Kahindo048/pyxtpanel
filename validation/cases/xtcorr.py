"""Overall, within and between correlations against Stata (E.6.4).

For each pair Stata keeps the rows where both variables are observed, then ``correlate`` gives
the overall correlation, and ``correlate`` after ``xtdata, fe`` and ``xtdata, be`` the within and
between ones.
"""

from __future__ import annotations

import numpy as np

from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "xtcorr"
PANELS = (
    "balanced",
    "gaps",
    "entry_exit",
    "weak",
    "delta2",
    "unsorted_ids",
    "one_obs_unit",
    "re_effects",
)
PAIRS = (("x", "z"), ("x", "lx"), ("z", "lx"))
KINDS = ("overall", "within", "between")
TERMS = {"x": "x", "z": "z", "lx": "L.x"}


def _option(name: str) -> str:
    delta = panels.SPECS[name].delta
    return f", delta({delta})" if delta != 1 else ""


def _python(name: str) -> dict[str, np.ndarray]:
    from pyxtpanel import xtset

    panel = xtset(panels.load(name), "id", "t", delta=panels.SPECS[name].delta)
    result = panel.xtcorr([TERMS["x"], TERMS["z"], TERMS["lx"]])
    names = {v: k for k, v in TERMS.items()}
    rho, n = [], []
    for a, b in PAIRS:
        for kind in KINDS:
            rho.append(result.corr[kind].loc[TERMS[a], TERMS[b]])
            n.append(result.n[kind].loc[TERMS[a], TERMS[b]])
    del names
    return {"rho": np.array(rho).reshape(-1, 1), "n": np.array(n).reshape(-1, 1)}


def _register(name: str) -> None:
    tag = f"xtcorr_{name}"
    lines = [f"xtset id t{_option(name)}", "generate double lx = L.x"]
    scalars = []
    for a, b in PAIRS:
        for kind in KINDS:
            key = f"{a}_{b}_{kind}"
            scalars.append(key)
            lines += [
                f"scalar r_{key} = .",
                f"scalar n_{key} = .",
                "preserve",
                f"quietly keep if !missing({a}, {b})",
            ]
            if kind == "within":
                lines.append(f"quietly xtdata {a} {b}, fe clear")
            if kind == "between":
                lines.append(f"quietly xtdata {a} {b}, be clear")
            lines += [
                f"capture quietly correlate {a} {b}",
                f"if _rc == 0 scalar r_{key} = r(rho)",
                f"if _rc == 0 scalar n_{key} = r(N)",
                "restore",
            ]
    lines += [
        export_scalars(tag, "rho", [f"r_{s}" for s in scalars]),
        export_scalars(tag, "n", [f"n_{s}" for s in scalars]),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        return _python(name)


for _name in PANELS:
    _register(_name)
