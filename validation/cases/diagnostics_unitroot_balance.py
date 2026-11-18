"""The unit-root tests that need a balanced panel refuse an unbalanced one, as Stata does (L.7.1).

``xtunitroot llc``, ``ht``, ``breitung`` and ``hadri`` stop on a panel whose units have different
numbers of periods (r(498)); ``ips`` runs on it (its cases are in the file of the test) but not on
one with holes.  Each case runs the command on the unbalanced panel ``ur_unbalanced`` (40 units of
12 to 24 consecutive periods) or on ``gaps`` (interior holes), stores the return code, and the test
of this package must refuse the same panel with the same code.  ``pescadf`` (CIPS) does run on a
panel with holes, which ``cips`` refuses: a stated limit, not a case.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from pyxtpanel import PyxtpanelError, xtset
from pyxtpanel.diagnostics import (
    breitung,
    hadri,
    harris_tzavalis,
    im_pesaran_shin,
    levin_lin_chu,
)
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_unitroot"

#: tag -> (Stata command line, the test of this package)
REFUSALS: dict[str, tuple[str, Callable[..., object]]] = {
    "bal_llc": ("xtunitroot llc y, lags(1)", lambda p: levin_lin_chu(p, "y", lags=1)),
    "bal_ht": ("xtunitroot ht y", lambda p: harris_tzavalis(p, "y")),
    "bal_breitung": ("xtunitroot breitung y, lags(1)", lambda p: breitung(p, "y", lags=1)),
    "bal_hadri": ("xtunitroot hadri y", lambda p: hadri(p, "y")),
}


#: tag -> (Stata command line, the test of this package) on the panel ``gaps`` (interior holes)
GAPS: dict[str, tuple[str, Callable[..., object]]] = {
    "gap_ips": ("xtunitroot ips x, lags(0)", lambda p: im_pesaran_shin(p, "x", lags=0)),
}


def _register(
    tag: str, command: str, run: Callable[..., object], panel: str = "ur_unbalanced"
) -> None:
    @case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{panel}",
        stata=[
            "quietly xtset id t",
            f"capture quietly {command}",
            "scalar rc = _rc",
            export_scalars(tag, "rc", ["rc"]),
        ],
    )
    def _refused() -> dict[str, np.ndarray]:
        try:
            run(xtset(panels.load(panel), "id", "t"))
        except PyxtpanelError as error:
            return {"rc": np.array([[float(error.stata_code or 0)]])}
        return {"rc": np.array([[0.0]])}


for _tag, (_command, _run) in REFUSALS.items():
    _register(_tag, _command, _run)


for _tag, (_command, _run) in GAPS.items():
    _register(_tag, _command, _run, panel="gaps")
