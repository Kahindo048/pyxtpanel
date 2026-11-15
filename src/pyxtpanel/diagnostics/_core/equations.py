"""``Equations``: a system of estimation or test equations stacked by unit (L.6.1).

The serial-correlation test of Arellano and Bond reads the residuals, the regressors and the
instruments of the equations of a dynamic-panel estimator.  The estimator builds them (first
differences, forward orthogonal deviations, levels...) and hands them over in this layout, the
one of the array contract: rows sorted by unit and, within a unit, by calendar position, with the
start of each unit in ``offsets``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel._types import FloatArray, IntArray
from pyxtpanel.diagnostics._core.inputs import as_array

_i18n.register(
    "diag",
    {
        "diag.equations.offsets": {
            "en": "offsets must start at 0, increase and end at the number of rows ({n}); "
            "got {got}.",
            "fr": "offsets doit commencer à 0, croître et finir au nombre de lignes ({n}) ; "
            "reçu {got}.",
        },
        "diag.equations.time": {
            "en": "time_index must increase strictly within each unit; unit {unit} breaks the "
            "order at row {row}.",
            "fr": "time_index doit croître strictement dans chaque unité ; l'unité {unit} "
            "rompt l'ordre à la ligne {row}.",
        },
        "diag.equations.rows": {
            "en": "{name} must have one entry per equation ({n}), got {got}.",
            "fr": "{name} doit avoir une entrée par équation ({n}), reçu {got}.",
        },
    },
)


@dataclass(frozen=True, eq=False)
class Equations:
    """Equations of one system, stacked by unit.

    Attributes
    ----------
    residuals
        ``(n,)`` residual of each equation.
    regressors
        ``(n, p)`` regressors of each equation (transformed like the residuals).
    offsets
        ``(n_units + 1,)`` start of each unit: unit ``i`` owns the rows ``offsets[i]:offsets[i+1]``.
    time_index
        ``(n,)`` calendar position of each equation, strictly increasing within a unit; a gap in
        the calendar is a gap here (the lag of a missing period does not exist).
    instruments
        ``(n, L)`` instrument matrix (row ``j`` holds the instruments of equation ``j``), or
        ``None`` for a system that only serves as the test sample.
    """

    residuals: FloatArray
    regressors: FloatArray
    offsets: IntArray
    time_index: IntArray
    instruments: FloatArray | None = None

    def __post_init__(self) -> None:
        resid = as_array("residuals", self.residuals, ndim=1)
        n = resid.size
        object.__setattr__(self, "residuals", resid)
        object.__setattr__(
            self, "regressors", as_array("regressors", self.regressors, ndim=2, shape=(n, None))
        )
        if self.instruments is not None:
            object.__setattr__(
                self,
                "instruments",
                as_array("instruments", self.instruments, ndim=2, shape=(n, None)),
            )
        offsets = np.asarray(self.offsets, dtype=np.int64)
        valid = offsets.ndim == 1 and offsets.size >= 2 and offsets[0] == 0 and offsets[-1] == n
        if not valid or np.any(np.diff(offsets) < 0):
            raise SpecificationError(
                _i18n.message("diag.equations.offsets", n=n, got=offsets.tolist()[:6]),
                key="diag.equations.offsets",
            )
        time = np.asarray(self.time_index, dtype=np.int64)
        if time.shape != (n,):
            raise SpecificationError(
                _i18n.message("diag.equations.rows", name="time_index", n=n, got=time.size),
                key="diag.equations.rows",
            )
        step = np.diff(time)
        inside = np.ones(step.size, dtype=bool)
        boundary = offsets[1:-1] - 1  # the step between two units is not a step
        inside[boundary[(boundary >= 0) & (boundary < step.size)]] = False
        broken = np.flatnonzero(inside & (step <= 0))
        if broken.size:
            row = int(broken[0]) + 1
            unit = int(np.searchsorted(offsets, row, side="right") - 1)
            raise DataError(
                _i18n.message("diag.equations.time", unit=unit, row=row),
                key="diag.equations.time",
            )
        object.__setattr__(self, "offsets", offsets)
        object.__setattr__(self, "time_index", time)

    @property
    def n_units(self) -> int:
        """Number of units."""
        return int(self.offsets.size - 1)

    @property
    def unit_of_row(self) -> IntArray:
        """``(n,)`` unit code of each equation."""
        out: IntArray = np.repeat(np.arange(self.n_units, dtype=np.int64), np.diff(self.offsets))
        return out
