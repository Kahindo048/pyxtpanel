"""Lag and lead on the calendar (Stata's ``L.`` and ``F.``).

``L.x`` is the value of the *period* before, not of the previous row: a period a unit skipped
gives a missing value.  With ``delta(2)`` or dates, one period is one step of the calendar.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel._types import FloatArray, IntArray
from pyxtpanel.ops._spec import OpSpec, operator
from pyxtpanel.panel.index import PanelIndex, shift_positions

__all__ = ["lag", "lead", "shift_positions", "take"]


def take(values: FloatArray, positions: IntArray) -> FloatArray:
    """``values`` at ``positions`` (rows), ``NaN`` where the position is ``-1``."""
    out: FloatArray = np.asarray(values, dtype=np.float64)[np.maximum(positions, 0)]
    out[positions < 0] = np.nan
    return out


@operator(
    "lag",
    OpSpec(
        kind="shift",
        missing="the period k steps before is not observed for the unit, or its value is missing",
        stata="L.x, L2.x, ...",
        summary="Value k periods earlier on the calendar of the panel.",
    ),
)
def lag(index: PanelIndex, values: FloatArray, k: int = 1) -> FloatArray:
    """Value ``k`` periods earlier (``L{k}.x``); ``k`` may be negative (a lead)."""
    # Stata: L.x looks for the period before, not the previous row
    return take(values, shift_positions(index, int(k)))


@operator(
    "lead",
    OpSpec(
        kind="shift",
        missing="the period k steps after is not observed for the unit, or its value is missing",
        stata="F.x, F2.x, ...",
        summary="Value k periods later on the calendar of the panel.",
    ),
)
def lead(index: PanelIndex, values: FloatArray, k: int = 1) -> FloatArray:
    """Value ``k`` periods later (``F{k}.x``)."""
    return take(values, shift_positions(index, -int(k)))
