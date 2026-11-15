"""``tsfill`` and ``tsfill, full`` against Stata (C.6.2): rows added and their content.

Stata fills the data, sorts them by (id, t) and exports ``id, t, x, z``; the package does the same
with :meth:`PanelData.tsfill`.  Added rows have missing ``x`` and ``z`` on both sides.
"""

from __future__ import annotations

import numpy as np

from validation import panels
from validation.registry import case, export

DOMAIN = "tsfill"
PANELS = ("gaps", "entry_exit", "weak", "delta2", "quarterly", "monthly", "unsorted_ids")


def _python(name: str, full: bool) -> dict[str, np.ndarray]:
    from pyxtpanel import xtset

    spec = panels.SPECS[name]
    panel = xtset(panels.load(name), "id", "t", delta=spec.delta).tsfill(full=full)
    data = panel.to_frame(order="sorted")
    time = panel.axis.stata_values(panel.index.rank).astype(float)
    table = np.column_stack([data["id"].to_numpy(float), time, data["x"], data["z"]])
    return {"data": table, "nrows": np.array([[float(len(data))]])}


def _register(name: str, full: bool) -> None:
    spec = panels.SPECS[name]
    tag = f"tsfill_{'full_' if full else ''}{name}"
    option = f", delta({spec.delta})" if spec.delta != 1 else ""
    lines = [
        f"xtset id t{option}",
        "tsfill, full" if full else "tsfill",
        "sort id t",
        export(tag, "data", 'st_data(., ("id", "t", "x", "z"))'),
        export(tag, "nrows", "st_nobs()"),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name, full: bool = full) -> dict[str, np.ndarray]:
        return _python(name, full)


for _name in PANELS:
    _register(_name, full=False)
    _register(_name, full=True)
