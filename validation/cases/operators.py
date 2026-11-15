"""Time-series operators ``L. F. D. S.`` and their combinations against Stata (D.4.1).

Stata computes each term with ``generate`` after ``xtset``; the package evaluates the same
spelling with :meth:`PanelData.operate`.  Panels with holes, late entry, a step of 2, quarterly
and monthly dates, unsorted identifiers, a one-observation unit and ignored missing keys exercise
the rule that a lag is the value of the *period* before, not of the previous row.
"""

from __future__ import annotations

import numpy as np

from validation import panels
from validation.registry import case, export

DOMAIN = "operators"
PANELS = (
    "gaps",
    "entry_exit",
    "weak",
    "delta2",
    "quarterly",
    "monthly",
    "unsorted_ids",
    "one_obs_unit",
    "single_unit",
    "single_period",
    "all_missing_var",
)
TERMS_X = [
    "L.x",
    "L2.x",
    "L3.x",
    "F.x",
    "F2.x",
    "D.x",
    "D2.x",
    "S2.x",
    "S3.x",
    "LD.x",
    "L2D.x",
    "FD.x",
    "DS2.x",
    "LS2.x",
    "D2S2.x",
    "L(1/2).x",
]
TERMS_Z = ["L.z", "F.z", "D.z", "L2.z", "LD.z", "S2.z", "D2.z"]


def _expand(terms: list[str]) -> list[str]:
    """Terms as Stata writes them, ranges expanded (``L(1/2).x`` -> ``L.x L2.x``)."""
    from pyxtpanel.ops.parse import parse_terms

    return list(dict.fromkeys(t.name for t in parse_terms(terms)))


def _python(name: str, terms: list[str], missing_keys: str) -> dict[str, np.ndarray]:
    from pyxtpanel import xtset

    spec = panels.SPECS[name]
    frame = panels.load(name)
    panel = xtset(frame, "id", "t", delta=spec.delta, missing_keys=missing_keys)  # type: ignore[arg-type]
    values = panel.operate(*terms)
    rows = panel.index.source_row
    table = np.column_stack(
        [frame["id"].to_numpy(float)[rows], *(values[c].to_numpy()[rows] for c in values.columns)]
    )
    return {"data": table}


def _register(name: str, missing_keys: str = "error") -> None:
    spec = panels.SPECS[name]
    tag = f"operators_{name}"
    option = f", delta({spec.delta})" if spec.delta != 1 else ""
    terms = _expand([*TERMS_X, *TERMS_Z])
    lines = [f"xtset id t{option}"]
    lines += [f"quietly generate double c{j} = {term}" for j, term in enumerate(terms, start=1)]
    columns = ", ".join(['"id"', *(f'"c{j}"' for j in range(1, len(terms) + 1))])
    lines += [
        "drop if id >= . | t >= .",
        "sort id t",
        export(tag, "data", f"st_data(., ({columns}))"),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        return _python(name, terms, missing_keys)


for _name in PANELS:
    _register(_name)
_register("missing_keys", missing_keys="drop")
