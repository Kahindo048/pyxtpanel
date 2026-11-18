"""The cointegration tests refuse the panels Stata refuses, with the same code (L.7.1).

``xtcointtest kao`` needs a balanced panel, ``pedroni`` and ``westerlund`` accept a panel whose units
have different lengths but not one with holes.  Each case runs the command on the panel ``gaps``
(interior holes), stores the return code, and the test of this package must refuse the same panel
with the same code.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from pyxtpanel import PyxtpanelError, xtset
from pyxtpanel.diagnostics import kao, pedroni, westerlund
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_cointegration"

#: tag -> (Stata command line, the test of this package)
REFUSALS: dict[str, tuple[str, Callable[..., object]]] = {
    "gap_kao": ("xtcointtest kao x z", lambda p: kao(p, "x", ["z"])),
    "gap_pedroni": (
        "xtcointtest pedroni x z, lags(1)",
        lambda p: pedroni(p, "x", ["z"], kernel_lags=1),
    ),
    "gap_westerlund": ("xtcointtest westerlund x z", lambda p: westerlund(p, "x", ["z"])),
}


def _register(tag: str, command: str, run: Callable[..., object]) -> None:
    @case(
        tag,
        domain=DOMAIN,
        dataset="sim_gaps",
        stata=[
            "quietly xtset id t",
            f"capture quietly {command}",
            "scalar rc = _rc",
            export_scalars(tag, "rc", ["rc"]),
        ],
    )
    def _refused() -> dict[str, np.ndarray]:
        try:
            run(xtset(panels.load("gaps"), "id", "t"))
        except PyxtpanelError as error:
            return {"rc": np.array([[float(error.stata_code or 0)]])}
        return {"rc": np.array([[0.0]])}


for _tag, (_command, _run) in REFUSALS.items():
    _register(_tag, _command, _run)
