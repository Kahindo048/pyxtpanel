"""``xttab`` and ``xttrans`` against Stata (E.6.3).

Discrete variables are made from ``z`` (three states, missing kept) and ``x`` (four states) in
Stata and in Python with the same cut points.  ``xttab`` is compared with ``r(results)``;
``xttrans`` with its ``r()`` scalars and with the ``tabulate`` of ``was`` and ``is`` rebuilt as
the ado file does (the matrix is not returned).
"""

from __future__ import annotations

import numpy as np

from validation import panels
from validation.registry import case, export, export_scalars

DOMAIN = "xttab"
PANELS = (
    "balanced",
    "gaps",
    "entry_exit",
    "weak",
    "delta2",
    "quarterly",
    "unsorted_ids",
    "one_obs_unit",
    "single_unit",
    "re_effects",
)

MAKE = [
    "generate byte dz = 0 if z < -0.5",
    "replace dz = 1 if z >= -0.5 & z < 0.5",
    "replace dz = 2 if z >= 0.5 & !missing(z)",
    "generate byte dx = 0 if x < -1",
    "replace dx = 1 if x >= -1 & x < 0",
    "replace dx = 2 if x >= 0 & x < 1",
    "replace dx = 3 if x >= 1 & !missing(x)",
]


def _option(name: str) -> str:
    delta = panels.SPECS[name].delta
    return f", delta({delta})" if delta != 1 else ""


def _panel(name: str):  # type: ignore[no-untyped-def]
    from pyxtpanel import xtset

    frame = panels.load(name)
    z, x = frame["z"], frame["x"]
    frame["dz"] = np.select([z < -0.5, z < 0.5, z >= 0.5], [0.0, 1.0, 2.0], default=np.nan)
    frame["dx"] = np.select([x < -1, x < 0, x < 1, x >= 1], [0.0, 1.0, 2.0, 3.0], default=np.nan)
    return xtset(frame, "id", "t", delta=panels.SPECS[name].delta)


def _tab(name: str, var: str) -> dict[str, np.ndarray]:
    table = _panel(name).xttab(var).table
    levels = [float(v) for v in table.index[:-1]] + [float("nan")]
    matrix = np.column_stack([levels, table.to_numpy()])
    return {"results": matrix}


def _trans(name: str, var: str) -> dict[str, np.ndarray]:
    result = _panel(name).xttrans(var)
    scalars = np.array([[result.n_trans], [result.n_rows], [result.n_cols]], dtype=float)
    return {
        "r": scalars,
        "cells": result.counts.to_numpy(float),
        "rows": result.counts.index.to_numpy(float).reshape(-1, 1),
        "cols": result.counts.columns.to_numpy(float).reshape(-1, 1),
    }


def _register(name: str, var: str) -> None:
    head = [f"xtset id t{_option(name)}", *MAKE]
    tag = f"xttab_{name}_{var}"
    lines = [*head, f"xttab {var}", export(tag, "results", 'st_matrix("r(results)")')]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines, atol=1e-9, rtol=1e-9)
    def _run_tab(name: str = name, var: str = var) -> dict[str, np.ndarray]:
        return _tab(name, var)

    tag = f"xttrans_{name}_{var}"
    lines = [
        *head,
        f"xttrans {var}",
        "scalar s_ntrans = r(n_trans)",
        "scalar s_nrows = r(n_rows)",
        "scalar s_ncols = r(n_cols)",
        export_scalars(tag, "r", ["s_ntrans", "s_nrows", "s_ncols"]),
        "sort id t",
        f"quietly by id: generate float was = {var} if _n < _N",
        f"quietly by id: generate float is = {var}[_n + 1] if _n < _N",
        "quietly generate byte touse = !missing(id, t, " + var + ")",
        "quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N",
        "quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)",
        export(tag, "cells", 'st_matrix("M")'),
        export(tag, "rows", 'st_matrix("R")'),
        export(tag, "cols", 'st_matrix("C")\''),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run_trans(name: str = name, var: str = var) -> dict[str, np.ndarray]:
        return _trans(name, var)


for _name in PANELS:
    for _var in ("dz", "dx"):
        _register(_name, _var)
