"""The declared panel (layer 2): data, declaration and sorted index, immutable.

``PanelBase`` wraps the user's DataFrame, never modifies it, and caches what was computed once at
declaration.  The public class ``PanelData`` (API layer) inherits from it and adds one-line
convenience methods.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from functools import cached_property
from pathlib import Path
from typing import Any, Literal, TypeVar

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError
from pyxtpanel._types import BoolArray, IntArray
from pyxtpanel.panel.arrays import PanelArrays, build_arrays
from pyxtpanel.panel.bridges import frame_of_linearmodels, to_linearmodels
from pyxtpanel.panel.checks import BalanceSummary, balance_summary
from pyxtpanel.panel.declare import Declaration, MissingKeys, _column, declare
from pyxtpanel.panel.fill import balanced_mask, tsfill
from pyxtpanel.panel.gaps import gaps_table
from pyxtpanel.panel.index import PanelIndex, subset_index
from pyxtpanel.panel.io import panel_from_stata, write_stata
from pyxtpanel.panel.select import check_not_empty, selection_mask
from pyxtpanel.time.axis import TimeAxis

P = TypeVar("P", bound="PanelBase")

_i18n.register(
    "panel",
    {
        "panel.with_frame_mismatch": {
            "en": "with_frame needs a frame with the same rows and the same keys as the declared "
            "one ({reason}); declare a new panel with xtset instead.",
            "fr": "with_frame exige un DataFrame de mêmes lignes et de mêmes clés que celui "
            "déclaré ({reason}) ; déclarez plutôt un nouveau panel avec xtset.",
        },
    },
)


def _snapshot(frame: pd.DataFrame) -> pd.DataFrame:
    """Copy of the user's frame that later edits of the original cannot reach.

    Copy-on-write (pandas 3) makes a shallow copy safe; before it, in-place edits of the original
    would show through a shallow copy, so the data are copied.
    """
    return frame.copy(deep=int(pd.__version__.split(".")[0]) < 3)


def _scalar(value: Any) -> Any:
    """Python scalar for a NumPy scalar (stable display across NumPy versions)."""
    return value.item() if isinstance(value, np.generic) else value


@dataclass(frozen=True, eq=False, repr=False)
class PanelBase:
    """A DataFrame declared as a panel.

    Build one with :meth:`from_frame` (the ``xtset`` of the package), not directly.

    Attributes
    ----------
    frame
        The data, as the user gave them (a private copy: the user's frame is never modified and
        later edits of it do not reach the panel).
    id_name, time_name
        Names of the identifier and the time variable.
    axis
        The calendar shared by all units.
    index
        Rows sorted by (unit, period) with the units delimited by offsets.
    n_dropped
        Observations ignored because of a missing key (``missing_keys="drop"``).
    """

    frame: pd.DataFrame
    id_name: str
    time_name: str
    axis: TimeAxis
    index: PanelIndex
    n_dropped: int = 0

    # ---------------------------------------------------------------- construction

    @classmethod
    def from_frame(
        cls: type[P],
        frame: pd.DataFrame,
        id: str | None = None,  # noqa: A002
        time: str | None = None,
        *,
        delta: int | None = None,
        freq: str | None = None,
        format: str | None = None,  # noqa: A002
        missing_keys: MissingKeys = "error",
    ) -> P:
        """Declare ``frame`` as a panel (Stata's ``xtset id time, delta()``).

        Parameters
        ----------
        frame
            Long data, one row per (unit, period).
        id, time
            Identifier and time variable, columns or index levels.  Both omitted: the two
            levels of a MultiIndex.
        delta
            Calendar units between two consecutive periods (default 1).
        freq
            ``"D"``, ``"W"``, ``"M"``, ``"Q"``, ``"H"`` or ``"Y"`` for dates; inferred when
            omitted.
        format
            Stata ``%t`` format when the time holds Stata-coded integers.
        missing_keys
            ``"error"`` or ``"drop"`` (ignore observations with a missing key, as Stata does).

        Returns
        -------
        PanelBase
            The declared panel (an instance of the class this is called on).

        Raises
        ------
        DeclarationError
            Missing keys, repeated (unit, period) pairs, no usable observation.
        TimeAxisError
            Time that is not integer, a date or a period, or off the calendar grid.
        SpecificationError
            Invalid option, unknown variable.
        """
        found: Declaration = declare(
            frame,
            id,
            time,
            delta=delta,
            freq=freq,
            format=format,
            missing_keys=missing_keys,
        )
        return cls(
            frame=_snapshot(frame),
            id_name=found.id_name,
            time_name=found.time_name,
            axis=found.axis,
            index=found.index,
            n_dropped=found.n_dropped,
        )

    # ---------------------------------------------------------------- structure

    @property
    def n_units(self) -> int:
        """Number of units."""
        return self.index.n_units

    @property
    def n_obs(self) -> int:
        """Number of observations (rows with both keys present)."""
        return self.index.n_obs

    def __len__(self) -> int:
        return self.n_obs

    @property
    def delta(self) -> int:
        """Calendar units between two consecutive periods."""
        return self.axis.delta

    @property
    def counts(self) -> IntArray:
        """Observations of each unit (``T_i``), in unit order."""
        return self.index.counts

    @property
    def unit_labels(self) -> np.ndarray:
        """Identifier of each unit, in sorted order."""
        return self.index.unit_labels

    @cached_property
    def balance_summary(self) -> BalanceSummary:
        """Balance with the counts behind it (Stata's ``r(balanced)``)."""
        return balance_summary(self.index.offsets, self.index.unit, self.index.rank)

    @property
    def balanced(self) -> Literal["strongly balanced", "weakly balanced", "unbalanced"]:
        """``"strongly balanced"``, ``"weakly balanced"`` or ``"unbalanced"`` (``r(balanced)``)."""
        return self.balance_summary.kind

    @property
    def n_periods(self) -> int:
        """Distinct periods observed in the whole panel."""
        return self.balance_summary.n_periods

    @cached_property
    def n_gaps(self) -> int:
        """Observations that follow their unit's previous one by more than one period."""
        same_unit = self.index.unit[1:] == self.index.unit[:-1]
        step = np.diff(self.index.rank)
        return int(np.count_nonzero(same_unit & (step != 1)))

    @property
    def has_gaps(self) -> bool:
        """Whether some unit skips a period (Stata's ``r(gaps)``)."""
        return self.n_gaps > 0

    # ---------------------------------------------------------------- time

    @property
    def t_min(self) -> Any:
        """First period, in the type of the time variable."""
        return _scalar(self.axis.values(int(np.min(self.index.rank)))[0]) if self.n_obs else None

    @property
    def t_max(self) -> Any:
        """Last period, in the type of the time variable."""
        return _scalar(self.axis.values(int(np.max(self.index.rank)))[0]) if self.n_obs else None

    @property
    def t_mean(self) -> float:
        """Average number of observations per unit (``T-bar`` of ``xtsum``)."""
        return float(self.n_obs / self.n_units) if self.n_units else float("nan")

    # ---------------------------------------------------------------- display

    def __repr__(self) -> str:
        first = self.axis.labels(int(np.min(self.index.rank)))[0]
        last = self.axis.labels(int(np.max(self.index.rank)))[0]
        note = ""
        if self.has_gaps:
            note = ", but with a gap" if self.n_gaps == 1 else ", but with gaps"
        lines = [
            f"{type(self).__name__}: {self.n_units} units, {self.n_obs} observations",
            f"Panel variable: {self.id_name} ({self.balanced})",
            f" Time variable: {self.time_name}, {first} to {last}{note}",
            f"         Delta: {self.axis.delta_text}",
        ]
        if self.n_dropped:
            lines.append(f"       Ignored: {self.n_dropped} observation(s) with a missing key")
        return "\n".join(lines)

    def _repr_html_(self) -> str:
        """HTML for Jupyter: the declaration, then the first rows in sorted order."""
        head = self.to_frame(order="sorted").head(10)
        title = (
            f"<p><b>{type(self).__name__}</b>: {self.n_units} units, {self.n_obs} observations "
            f"({self.balanced}); panel variable <code>{self.id_name}</code>, "
            f"time variable <code>{self.time_name}</code>, delta {self.axis.delta_text}</p>"
        )
        return title + head.to_html(max_rows=10)

    # ---------------------------------------------------------------- back to a DataFrame

    def to_frame(
        self,
        *,
        order: Literal["original", "sorted"] = "original",
        index: Literal["original", "keys"] = "original",
    ) -> pd.DataFrame:
        """Return the data of the panel as a new DataFrame.

        Parameters
        ----------
        order
            ``"original"`` keeps the order of the rows of the frame that was declared;
            ``"sorted"`` sorts them by unit then period.
        index
            ``"original"`` keeps the index of the frame; ``"keys"`` makes (id, time) the
            two-level index.

        Returns
        -------
        pandas.DataFrame
            A copy: only observations with both keys present (all of them unless
            ``missing_keys="drop"`` set some aside).
        """
        rows = self.index.source_row if order == "sorted" else np.sort(self.index.source_row)
        out = self.frame.iloc[rows]
        if index == "keys":
            names = (self.id_name, self.time_name)
            if self.id_name not in out.columns and self.id_name not in out.index.names:
                out = out.rename_axis(list(names))  # declared from an unnamed two-level index
            else:
                moved = [n for n in names if n in out.index.names]
                if moved:
                    out = out.reset_index(level=moved)
                out = out.set_index(list(names))
        result: pd.DataFrame = out.copy()
        return result

    # ---------------------------------------------------------------- selection and structure

    def select(
        self: P,
        where: str | Sequence[bool] | BoolArray | None = None,
        *,
        units: Sequence[Any] | None = None,
        periods: Sequence[Any] | slice | None = None,
        complete: Sequence[str] | None = None,
    ) -> P:
        """Restrict the panel to the observations that satisfy every condition.

        The sub-panel keeps the declaration and the calendar of the panel; the data are not
        copied.

        Parameters
        ----------
        where
            Boolean mask over the rows of the declared frame, or a string evaluated with
            :meth:`pandas.DataFrame.eval` (Stata's ``if``).  A missing value counts as False.
        units
            Identifiers of the units to keep.
        periods
            Time values to keep, or ``slice(start, stop)`` with both ends included.
        complete
            Columns that must be present and finite (Stata's ``markout``).

        Returns
        -------
        PanelBase
            The sub-panel (same class as ``self``).  Units left empty disappear.

        Raises
        ------
        DeclarationError
            If nothing is left (Stata code 2000).
        DataError
            If ``where`` is not boolean or a column is unknown.
        TimeAxisError
            If a period is not on the calendar grid.
        """
        keep = selection_mask(self, where=where, units=units, periods=periods, complete=complete)
        check_not_empty(keep)
        return replace(self, index=subset_index(self.index, keep))

    def gaps(self, *, only_with_gaps: bool = False) -> pd.DataFrame:
        """Per unit: first and last period, observations, breaks and missing periods inside.

        Parameters
        ----------
        only_with_gaps
            Keep only the units that have a gap.

        Returns
        -------
        pandas.DataFrame
            One row per unit.
        """
        return gaps_table(self, only_with_gaps=only_with_gaps)

    def tsfill(self: P, *, full: bool = False, indicator: str | None = None) -> P:
        """Add the missing periods (Stata's ``tsfill``) and return the new panel.

        Parameters
        ----------
        full
            Fill every unit from the first to the last period of the whole panel
            (``tsfill, full``) instead of between its own first and last observation.
        indicator
            Name of a boolean column that marks the added rows (none by default).

        Returns
        -------
        PanelBase
            A new panel whose frame has the added rows (identifier and time filled, every other
            column missing), sorted by unit then period, with a fresh index.  Key columns that
            were index levels become columns.
        """
        filled = tsfill(self, full=full, indicator=indicator)
        return type(self).from_frame(
            filled,
            self.id_name,
            self.time_name,
            delta=self.axis.delta,
            freq=self.axis.freq,
            format=self.axis.stata_format,
        )

    def balance(
        self: P,
        method: Literal["drop", "largest_block"] = "drop",
        *,
        min_obs: int | None = None,
    ) -> P:
        """Make the panel strongly balanced by dropping observations, never by inventing them.

        Parameters
        ----------
        method
            ``"drop"`` keeps the units observed in every period.  ``"largest_block"`` keeps the
            window of consecutive periods and the units complete in it that hold the most
            observations.
        min_obs
            Remove first the units with fewer observations than this.

        Returns
        -------
        PanelBase
            The balanced sub-panel.

        Raises
        ------
        DeclarationError
            If no unit is left.
        """
        keep = balanced_mask(self, method, min_obs)
        check_not_empty(keep)
        return replace(self, index=subset_index(self.index, keep))

    # ---------------------------------------------------------------- bridges to pandas

    @classmethod
    def from_multiindex(cls: type[P], frame: pd.DataFrame, **options: Any) -> P:
        """Declare a frame indexed by (unit, period) as a panel; options as in ``from_frame``."""
        return cls.from_frame(frame, **options)

    @classmethod
    def from_stata(
        cls: type[P],
        path: str | Path,
        id: str | None = None,  # noqa: A002
        time: str | None = None,
        **options: Any,
    ) -> P:
        """Read a ``.dta`` file and declare the panel it records (``xtset``), or the given one.

        The panel variable, the time variable and the step ``delta`` come from the file when it
        was saved after ``xtset``; ``id``, ``time`` and the other ``options`` of
        :meth:`from_frame` override them.  Names are never guessed.

        Raises
        ------
        SpecificationError
            If the file carries no ``xtset`` and ``id`` and ``time`` are not given.
        """
        return panel_from_stata(cls, path, id=id, time=time, **options)  # type: ignore[return-value]

    def to_stata(self, path: str | Path) -> None:
        """Write a Stata ``.dta`` file that keeps the ``xtset`` declaration of the panel.

        Parameters
        ----------
        path
            Destination ``.dta`` file.
        """
        write_stata(self, path)

    @classmethod
    def from_linearmodels(cls: type[P], obj: Any, **options: Any) -> P:
        """Declare the data of a ``linearmodels`` panel container (duck-typed) as a panel."""
        return cls.from_frame(frame_of_linearmodels(obj), **options)

    def to_linearmodels(self) -> Any:
        """Return a ``linearmodels.panel.data.PanelData`` (needs ``linearmodels``).

        Returns
        -------
        Any
            Special type indicating an unconstrained type.
        """
        return to_linearmodels(self)

    def to_multiindex(self, *, order: Literal["original", "sorted"] = "sorted") -> pd.DataFrame:
        """Return the data indexed by (id, time), sorted by default.

        Parameters
        ----------
        order
            ``"sorted"`` (default) or ``"original"``: the order of the rows.

        Returns
        -------
        pandas.DataFrame
            The data indexed by (id, time).
        """
        return self.to_frame(order=order, index="keys")

    def arrays(
        self,
        columns: Sequence[str],
        *,
        where: str | Sequence[bool] | BoolArray | None = None,
    ) -> PanelArrays:
        """Return the array contract for the estimators (see :class:`PanelArrays`).

        Parameters
        ----------
        columns
            Numeric variables to extract, in the column order wanted.
        where
            Optional restriction, as in :meth:`select` (a mask over the declared frame or a
            string for :meth:`pandas.DataFrame.eval`).

        Returns
        -------
        PanelArrays
            Frozen ``float64`` arrays sorted by unit then period, the unit offsets, the calendar
            position of each row and the way back to the declared frame.
        """
        panel = self if where is None else self.select(where)
        return build_arrays(panel, columns)

    def with_frame(self: P, frame: pd.DataFrame) -> P:
        """Return a panel of the same structure holding new values, without sorting again.

        For loops where only the values change (for example iterated factor models): the sorted
        index is reused, so no sorting or validation is repeated.  The new frame must have the
        same number of rows and the same identifier and time values, row by row.

        Parameters
        ----------
        frame
            New values: same number of rows, same identifier and time values, row by row.

        Returns
        -------
        PanelData
            A panel of the same structure holding the new values.

        Raises
        ------
        DataError
            If the rows or the keys differ from the declared frame.
        """
        reason = None
        if len(frame) != len(self.frame):
            reason = f"{len(frame)} rows instead of {len(self.frame)}"
        elif not self._same_keys(frame):
            reason = "the identifier or time values differ"
        if reason is not None:
            raise DataError(
                _i18n.message("panel.with_frame_mismatch", reason=reason),
                key="panel.with_frame_mismatch",
            )
        return replace(self, frame=_snapshot(frame))

    def _same_keys(self, frame: pd.DataFrame) -> bool:
        for name, level in ((self.id_name, 0), (self.time_name, 1)):
            try:
                new, old = _column(frame, name), _column(self.frame, name)
            except DataError:  # declared from an unnamed two-level index
                new = pd.Series(frame.index.get_level_values(level))
                old = pd.Series(self.frame.index.get_level_values(level))
            if not new.equals(old):
                return False
        return True

    def matches(self, frame: pd.DataFrame) -> bool:
        """Whether ``frame`` still equals the data this panel was declared from.

        Parameters
        ----------
        frame
            The DataFrame to compare with the data the panel was declared from.

        Returns
        -------
        bool
            ``True`` when ``frame`` holds the data the panel was declared from.
        """
        return bool(frame.equals(self.frame))

    def align(self, values: Any, name: str | None = None) -> pd.Series | pd.DataFrame:
        """Put values computed in sorted order back on the index of the declared frame.

        Parameters
        ----------
        values
            One value per observation of the panel, in sorted order (length ``n_obs``), or a
            2-D array with one column per output.
        name
            Name of the resulting Series.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            On the index and in the row order of the declared frame.  Rows that are not in the
            panel (set aside for a missing key, or outside a selection) are missing.
        """
        array = np.asarray(values)
        if array.shape[0] != self.n_obs:
            raise ValueError(f"expected {self.n_obs} values, got {array.shape[0]}")
        table = array.reshape(self.n_obs, -1)
        if table.dtype.kind == "f":  # scatter into a missing-filled array: no reindex of a frame
            scattered = np.full((len(self.frame), table.shape[1]), np.nan, dtype=table.dtype)
            scattered[self.index.source_row] = table
            data = pd.DataFrame(scattered, index=self.frame.index)
        else:
            data = pd.DataFrame(table, index=self.index.source_row)
            data = data.reindex(range(len(self.frame)))
            data.index = self.frame.index
        if array.ndim > 1:
            return data
        series: pd.Series = data.iloc[:, 0]
        series.name = name
        return series
