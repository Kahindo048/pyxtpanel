"""``xtsum`` against the ``r()`` scalars of Stata (E.6.2).

Balanced and unbalanced panels, a unit with one observation, one unit, one period, a variable
that is constant within units and one that is constant everywhere, and time-series operators.
"""

from __future__ import annotations

import numpy as np

from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "xtsum"
SCALARS = [
    "N",
    "n",
    "Tbar",
    "mean",
    "sd",
    "min",
    "max",
    "sd_b",
    "min_b",
    "max_b",
    "sd_w",
    "min_w",
    "max_w",
]

CASES: list[tuple[str, tuple[str, ...]]] = [
    ("balanced", ("x", "z")),
    ("gaps", ("x", "z", "L.x", "D.z", "L2.z")),
    ("entry_exit", ("x", "z")),
    ("weak", ("x", "z")),
    ("delta2", ("x", "z", "L.x")),
    ("quarterly", ("x", "z", "D.x")),
    ("unsorted_ids", ("x", "z")),
    ("one_obs_unit", ("x", "z")),
    ("single_unit", ("x", "z")),
    ("single_period", ("x", "z")),
    ("constants", ("x", "cu", "ck")),
    ("re_effects", ("x", "z")),
]


def _option(name: str) -> str:
    delta = panels.SPECS[name].delta
    return f", delta({delta})" if delta != 1 else ""


def _safe(term: str) -> str:
    return term.replace(".", "_")


def _python(name: str, term: str) -> dict[str, np.ndarray]:
    from pyxtpanel import xtset

    panel = xtset(panels.load(name), "id", "t", delta=panels.SPECS[name].delta)
    row = panel.xtsum([term]).scalars.iloc[0]
    return {"r": row[SCALARS].to_numpy(float).reshape(-1, 1)}


def _register(name: str, term: str) -> None:
    tag = f"xtsum_{name}_{_safe(term)}"
    lines = [
        f"xtset id t{_option(name)}",
        f"xtsum {term}",
        *(f"scalar s_{s} = r({s})" for s in SCALARS),
        export_scalars(tag, "r", [f"s_{s}" for s in SCALARS]),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name, term: str = term) -> dict[str, np.ndarray]:
        return _python(name, term)


for _name, _terms in CASES:
    for _term in _terms:
        _register(_name, _term)
