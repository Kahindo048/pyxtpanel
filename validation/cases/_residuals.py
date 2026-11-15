"""Cases on the within residuals of ``xtreg, fe``, exported from Stata and tested in Python.

The cross-section dependence tests read a variable of a declared panel.  For residuals, Stata
fits ``xtreg <model>, fe`` and ``predict ehat, e``; the residuals with their identifier and
period are exported and rebuilt into a panel on the Python side, so the function reads exactly
Stata's residuals.  ``stata_lines`` compute the Stata scalars named in ``scalars`` (the reference),
``python`` returns the same numbers from the rebuilt panel.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd

from pyxtpanel import xtset
from pyxtpanel.panel.base import PanelBase
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars


def register_residual_case(
    tag: str,
    *,
    domain: str,
    dataset: str,
    model: str,
    stata_lines: Sequence[str],
    scalars: Sequence[str],
    python: Callable[[PanelBase], Sequence[float]],
) -> None:
    """Register the case ``tag``: Stata's scalars against ``python`` on Stata's residuals."""
    lines = [
        "quietly xtset id t",
        f"quietly xtreg {model}, fe",
        "quietly predict double ehat, e",
        *stata_lines,
        "preserve",
        "quietly keep if !missing(ehat)",
        export(tag, "e", 'st_data(., "ehat")'),
        export(tag, "id", 'st_data(., "id")'),
        export(tag, "t", 'st_data(., "t")'),
        "restore",
        export_scalars(tag, "r", list(scalars)),
    ]

    def run() -> dict[str, np.ndarray]:
        s = stata_objects(domain, tag)
        frame = pd.DataFrame({"id": s["id"].ravel(), "t": s["t"].ravel(), "ehat": s["e"].ravel()})
        values = python(xtset(frame, "id", "t"))
        return {"r": np.asarray(values, dtype=float), "e": s["e"], "id": s["id"], "t": s["t"]}

    case(tag, domain=domain, dataset=f"sim_{dataset}", stata=lines, atol=1e-9, rtol=1e-9)(run)
