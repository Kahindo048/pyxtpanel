"""Read-only view of the calendar of a panel, for the packages that build on ``pyxtpanel``.

``panel.calendar`` is what a companion package (a plotting package, an estimator) may rely on to
turn the calendar positions of :attr:`pyxtpanel.contract.PanelArrays.time_index` back into time
values, dates or printable labels.  It stands in front of :class:`~pyxtpanel.time.axis.TimeAxis`,
whose internals (ordinals, Stata formats, lookup keys) can change without notice.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import TimeAxisError
from pyxtpanel._types import IntArray
from pyxtpanel.time.axis import Kind, TimeAxis
from pyxtpanel.time.formats import Freq

_i18n.register(
    "time",
    {
        "time.no_dates": {
            "en": "This calendar counts plain integers and has no dates: use values() or labels(), "
            "or declare the panel with a date frequency or a Stata format.",
            "fr": "Ce calendrier compte des entiers simples et n'a pas de dates : utilisez "
            "values() ou labels(), ou déclarez le panel avec une fréquence de dates ou un format "
            "Stata.",
        },
    },
)


class Calendar:
    """The calendar shared by all the units of a panel (read only).

    Obtain it with :attr:`pyxtpanel.PanelData.calendar`.  A *position* is the calendar position
    of a period, ``(t - t0) / delta``: the values of
    :attr:`pyxtpanel.contract.PanelArrays.time_index`.  Consecutive periods differ by exactly
    one position, whatever the unit of the time variable.
    """

    __slots__ = ("_axis",)

    def __init__(self, axis: TimeAxis) -> None:
        self._axis = axis

    @property
    def kind(self) -> Kind:
        """``"integer"``, ``"datetime"``, ``"period"`` or ``"stata"``: how the time was given."""
        return self._axis.kind

    @property
    def freq(self) -> Freq | None:
        """Calendar unit (``"D"``, ``"W"``, ``"SW"``, ``"M"``, ``"Q"``, ``"H"``, ``"Y"``).

        ``None`` for plain integers.
        """
        return self._axis.freq

    @property
    def delta(self) -> int:
        """Calendar units between two consecutive positions (Stata's ``delta()``)."""
        return self._axis.delta

    @property
    def n_steps(self) -> int:
        """Number of positions from the first to the last period, inclusive."""
        return self._axis.n_steps

    def values(self, positions: IntArray | Sequence[int] | int) -> pd.Index:
        """Time values, in the type of the declared time variable, at ``positions``.

        Parameters
        ----------
        positions
            Calendar positions (for example ``PanelArrays.time_index``).

        Returns
        -------
        pandas.Index
            Integers, timestamps (first day of the period) or periods, as the data were given.
        """
        return self._axis.values(np.asarray(positions, dtype=np.int64))

    def datetimes(self, positions: IntArray | Sequence[int] | int) -> pd.DatetimeIndex:
        """First day of the periods at ``positions``, as timestamps.

        Parameters
        ----------
        positions
            Calendar positions.

        Returns
        -------
        pandas.DatetimeIndex
            One timestamp per position.

        Raises
        ------
        TimeAxisError
            If the time variable holds plain integers: there are no dates to return.
        """
        if self._axis.freq is None:
            raise TimeAxisError(_i18n.message("time.no_dates"), key="time.no_dates")
        return self._axis.datetimes(np.asarray(positions, dtype=np.int64))

    def labels(self, positions: IntArray | Sequence[int] | int) -> list[str]:
        """Printable labels (``"2020q1"``, ``"1987"``, ``"01jan2020"``) of the periods.

        Parameters
        ----------
        positions
            Calendar positions.

        Returns
        -------
        list of str
            One label per position.
        """
        return self._axis.labels(np.asarray(positions, dtype=np.int64))

    def positions_of(self, values: Sequence[Any] | np.ndarray | pd.Index | pd.Series) -> IntArray:
        """Calendar positions of time values given in the type of the declared time variable.

        Parameters
        ----------
        values
            Time values (integers, dates, periods or Stata-coded integers, as declared).

        Returns
        -------
        numpy.ndarray
            ``int64`` positions; they may be negative or beyond ``n_steps - 1`` for periods
            outside the observed range.

        Raises
        ------
        TimeAxisError
            If a value is missing, not numeric or a date, or not on the calendar grid.
        """
        return self._axis.rank_of(values)

    def __repr__(self) -> str:
        return (
            f"Calendar(kind={self.kind!r}, freq={self.freq!r}, delta={self.delta}, "
            f"n_steps={self.n_steps})"
        )
