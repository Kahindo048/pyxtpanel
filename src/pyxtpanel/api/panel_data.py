"""The public panel object (layer 6).

``PanelData`` is a :class:`~pyxtpanel.panel.base.PanelBase` with one-line convenience methods
that delegate to the functions of the lower layers.  The methods are added as the features that
back them are written (G.1.5).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd

from pyxtpanel.api._lazy import run_panel_test
from pyxtpanel.describe.attrition import Attrition, attrition
from pyxtpanel.describe.corr import AutoCorr, XtCorr, autocorr, xtcorr
from pyxtpanel.describe.missing import (
    CompleteCases,
    MissingReport,
    complete_cases,
    missing_report,
)
from pyxtpanel.describe.structure import XtDescribe, xtdescribe
from pyxtpanel.describe.summary import XtSum, xtsum
from pyxtpanel.describe.tab import XtTab, xttab
from pyxtpanel.describe.trans import XtTrans, xttrans
from pyxtpanel.explore.audit import AuditReport, audit
from pyxtpanel.explore.balance import BalanceCost, Unbalancedness, balance_cost, unbalancedness
from pyxtpanel.explore.equations import Transform, UsableEquations, usable_equations
from pyxtpanel.explore.event_time import EventTime, event_time
from pyxtpanel.explore.gaps import GapTypes, gap_types
from pyxtpanel.explore.outliers import LevelShifts, WithinOutliers, level_shifts, within_outliers
from pyxtpanel.explore.reshape import Array3D, to_array3d, to_wide
from pyxtpanel.explore.spells import Spells, spells
from pyxtpanel.explore.variation import ResidualVariation, Variation, residual_variation, variation
from pyxtpanel.ops.frame import apply_operator, eliminate_effects, operate, term_columns, transform
from pyxtpanel.ops.groups import group_size
from pyxtpanel.ops.within import theta_re
from pyxtpanel.panel.arrays import PanelArrays, arrays_from_values
from pyxtpanel.panel.base import PanelBase
from pyxtpanel.plot._spec import plot as _plot

if TYPE_CHECKING:
    from pyxtpanel.diagnostics import TestResult


@dataclass(frozen=True, eq=False, repr=False)
class PanelData(PanelBase):
    """A DataFrame declared as a panel, as Stata's ``xtset`` does.

    Create one with :func:`pyxtpanel.xtset`.  The object is immutable and never modifies the
    DataFrame it wraps.

    Examples
    --------
    >>> import pandas as pd
    >>> import pyxtpanel as xt
    >>> data = pd.DataFrame({"firm": [1, 1, 2, 2], "year": [2000, 2001, 2000, 2001]})
    >>> panel = xt.xtset(data, "firm", "year")
    >>> panel.balanced, panel.n_units, panel.t_min, panel.t_max
    ('strongly balanced', 2, 2000, 2001)
    """

    # ------------------------------------------------------------ time-series operators

    def operate(self, *terms: str) -> pd.DataFrame:
        """Evaluate Stata-style terms such as ``"L2.x"``, ``"D.y"``, ``"LD.z"`` or ``"L(1/3).x"``.

        Returns a DataFrame with one column per term, named in canonical form, on the index of
        the declared frame.  A lag is the value of the *period* before: a missing period gives a
        missing value, as in Stata.

        Parameters
        ----------
        terms
            Terms such as ``"L2.x"``, ``"D.y"``, ``"LD.z"`` or ``"L(1/3).x"``.

        Returns
        -------
        pandas.DataFrame
            One column per term, named in canonical form, on the index of the declared frame.
        """
        return operate(self, terms)

    def term_arrays(self, items: str | Sequence[str]) -> PanelArrays:
        """Return the array contract for variables or Stata-style terms (``"x"``, ``"L.x"``).

        The same object as :meth:`arrays`, for terms: a name that is a column of the frame is
        read as it is, anything else is read as a term (``"L2.x"``, ``"D.x"``).  The names in
        ``columns`` are canonical (``"L2.x"``).  There is no ``where``: restrict the panel with
        :meth:`select` first, since a restriction applied before or after a lag does not give the
        same values.

        Parameters
        ----------
        items
            Column names or terms, in the order wanted (one string or a sequence).

        Returns
        -------
        PanelArrays
            Frozen ``float64`` arrays sorted by unit then period, as :meth:`arrays` returns.

        Raises
        ------
        SpecificationError
            If a term is malformed.
        DataError
            If a variable is unknown or not numeric.
        """
        names, values = term_columns(self, items)
        return arrays_from_values(self, names, values)

    def lag(
        self, variables: str | Sequence[str], k: int | Sequence[int] = 1
    ) -> pd.Series | pd.DataFrame:
        """Value ``k`` periods earlier (``L.x``, ``L2.x``...).

        A Series for one variable and one ``k``, otherwise a DataFrame.

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        k
            Lag in periods, or several lags (a DataFrame comes back).

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return apply_operator(self, "lag", variables, k)

    def lead(
        self, variables: str | Sequence[str], k: int | Sequence[int] = 1
    ) -> pd.Series | pd.DataFrame:
        """Value ``k`` periods later (``F.x``, ``F2.x``...).

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        k
            Lead in periods, or several leads (a DataFrame comes back).

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return apply_operator(self, "lead", variables, k)

    def diff(
        self, variables: str | Sequence[str], order: int | Sequence[int] = 1
    ) -> pd.Series | pd.DataFrame:
        """Difference with the previous period (``D.x``); ``order=2`` is ``D2.x``.

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        order
            Order of the difference, or several orders (``order=2`` is ``D2.x``).

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return apply_operator(self, "diff", variables, order)

    def seasonal(
        self, variables: str | Sequence[str], period: int | Sequence[int] = 1
    ) -> pd.Series | pd.DataFrame:
        """Difference with the value ``period`` periods earlier (``S12.x``).

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        period
            Length of the season in periods (12 for monthly data), or several.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return apply_operator(self, "seasonal", variables, period)

    # ------------------------------------------------------------ transformations

    def transform(
        self, name: str, variables: str | Sequence[str], *, listwise: bool = False, **params: Any
    ) -> pd.Series | pd.DataFrame:
        """Apply a registered transformation (``"within"``, ``"fod"``, ``"rolling"``...) by name.

        ``variables`` is one name (a Series comes back) or several (a DataFrame with the same
        column names), always on the index of the declared frame.  ``listwise=True`` first drops
        the rows where any of the variables is missing, as Stata's ``xtdata`` does.

        Parameters
        ----------
        name
            A name of :data:`pyxtpanel.ops._spec.OPERATORS` (``"within"``, ``"between"``,
            ``"fod"``, ``"rolling"``...).
        variables
            One name (a Series is returned) or several (a DataFrame with the same column names).
        listwise
            Drop, for all the variables, the rows where any of them is missing, before computing
            (what Stata's ``xtdata`` does).  Default: each variable uses its own observed rows.
        params
            Parameters of the transformation.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(self, name, variables, listwise=listwise, **params)

    def within(
        self,
        variables: str | Sequence[str],
        *,
        effects: Literal["unit", "time", "twoway"] = "unit",
        add_mean: bool = False,
        listwise: bool = False,
        tol: float = 1e-12,
        maxiter: int = 10_000,
    ) -> pd.Series | pd.DataFrame:
        """Deviation from the mean of the unit, of the period, or of both (two-way).

        ``add_mean=True`` adds the overall mean back (``xtdata, fe``).  The two-way result of an
        unbalanced panel is the residual of the regression on unit and period dummies, found by
        alternating projections (``tol``, ``maxiter``); its convergence report is in
        ``result.attrs["convergence"]``.

        Parameters
        ----------
        variables
            One name (a Series is returned) or several (a DataFrame with the same column names).
        effects
            ``"unit"`` (default), ``"time"`` or ``"twoway"``: the effects removed.
        add_mean
            Add the overall mean back (what ``xtdata, fe`` does).
        listwise
            Drop, for all the variables, the rows where any of them is missing, before computing
            (what Stata's ``xtdata`` does).  Default: each variable uses its own observed rows.
        tol
            Convergence tolerance of the alternating projections (two-way, unbalanced panels).
        maxiter
            Maximum number of alternating projections.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return eliminate_effects(
            self,
            variables,
            effects,
            add_mean=add_mean,
            listwise=listwise,
            tol=tol,
            maxiter=maxiter,
        )

    def between(
        self, variables: str | Sequence[str], *, listwise: bool = False
    ) -> pd.Series | pd.DataFrame:
        """Mean of the unit over its observed values, repeated on the rows with a value.

        Parameters
        ----------
        variables
            One name (a Series is returned) or several (a DataFrame with the same column names).
        listwise
            Drop, for all the variables, the rows where any of them is missing, before computing
            (what Stata's ``xtdata`` does).  Default: each variable uses its own observed rows.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(self, "between", variables, listwise=listwise)

    def quasi_demean(
        self,
        variables: str | Sequence[str],
        theta: float | Sequence[float],
        *,
        listwise: bool = False,
    ) -> pd.Series | pd.DataFrame:
        """``x − θ x̄_i``, ``θ`` a scalar or one value per unit (order of ``unit_labels``).

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        theta
            Scalar, or one value per unit in the order of ``unit_labels``.
        listwise
            Drop, for all the variables, the rows where any of them is missing, before computing
            (what Stata's ``xtdata`` does).  Default: each variable uses its own observed rows.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(
            self, "quasi_demean", variables, listwise=listwise, theta=np.asarray(theta, float)
        )

    def theta_re(self, sigma_u: float, sigma_e: float) -> pd.Series:
        """Random-effects ``θ_i`` of each unit for the given standard deviations.

        Parameters
        ----------
        sigma_u
            Standard deviation of the unit effect.
        sigma_e
            Standard deviation of the idiosyncratic error.

        Returns
        -------
        pandas.Series
            ``theta_i`` of each unit, in the order of ``unit_labels``.
        """
        return pd.Series(
            theta_re(sigma_u, sigma_e, self.counts), index=self.unit_labels, name="theta"
        )

    def fod(
        self, variables: str | Sequence[str], *, scale: bool = True
    ) -> pd.Series | pd.DataFrame:
        """Forward orthogonal deviations: later observed rows of the unit, ``√(F/(F+1))`` scale.

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        scale
            Multiply by ``sqrt(F/(F+1))``, ``F`` being the number of later observed rows, so
            that uncorrelated errors stay uncorrelated and keep their variance (default);
            ``False`` leaves the plain deviation.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(self, "fod", variables, scale=scale)

    def group_mean(
        self,
        variables: str | Sequence[str],
        by: Literal["unit", "time"] = "unit",
        *,
        at_missing: bool = False,
    ) -> pd.Series | pd.DataFrame:
        """Mean of the unit (or period) repeated on the rows.

        ``at_missing=True`` also fills the rows that have no value of their own.

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        by
            ``"unit"`` (default) or ``"time"``: the group the mean is taken over.
        at_missing
            Also fill the rows that have no value of their own.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(self, "group_mean", variables, by=by, at_missing=at_missing)

    def group_count(
        self,
        variables: str | Sequence[str],
        by: Literal["unit", "time"] = "unit",
        *,
        at_missing: bool = False,
    ) -> pd.Series | pd.DataFrame:
        """Count the observed values of the unit (or period) and repeat it on the rows.

        ``at_missing=True`` also fills the rows that have no value of their own.

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        by
            ``"unit"`` (default) or ``"time"``: the group the values are counted in.
        at_missing
            Also fill the rows that have no value of their own.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(self, "group_count", variables, by=by, at_missing=at_missing)

    def group_size(self, by: Literal["unit", "time"] = "unit") -> pd.Series:
        """Rows of the unit (``T_i``) or of the period (``N_t``), repeated on each row.

        Parameters
        ----------
        by
            ``"unit"`` (default) for the rows of the unit ``T_i``, ``"time"`` for the rows of
            the period ``N_t``.

        Returns
        -------
        pandas.Series
            Rows of the unit or of the period, on the index of the declared frame.
        """
        series = self.align(group_size(self.index, by), f"T_{'i' if by == 'unit' else 't'}")
        assert isinstance(series, pd.Series)
        return series

    def rolling(
        self,
        variables: str | Sequence[str],
        window: int,
        stat: Literal["mean", "sum", "std", "min", "max"] = "mean",
        *,
        min_obs: int | None = None,
    ) -> pd.Series | pd.DataFrame:
        """Statistic over the ``window`` calendar periods ending at each row (gaps respected).

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        window
            Length of the window, in calendar periods.
        stat
            ``"mean"``, ``"sum"``, ``"std"``, ``"min"`` or ``"max"``.
        min_obs
            Observed values needed in the window for the result not to be missing; default: all
            of them.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(self, "rolling", variables, window=window, stat=stat, min_obs=min_obs)

    def cumulative(
        self, variables: str | Sequence[str], stat: Literal["sum", "mean"] = "sum"
    ) -> pd.Series | pd.DataFrame:
        """Return the running sum (or mean) of the unit's observed values.

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        stat
            ``"sum"`` or ``"mean"``.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(self, "cumulative", variables, stat=stat)

    def growth(
        self, variables: str | Sequence[str], periods: int = 1, *, log: bool = False
    ) -> pd.Series | pd.DataFrame:
        """Growth rate over ``periods`` calendar periods (relative, or log difference).

        Parameters
        ----------
        variables
            One column name (a Series is returned) or several (a DataFrame with the same column
            names).
        periods
            Number of calendar periods the growth is measured over.
        log
            Log difference instead of the relative change.

        Returns
        -------
        pandas.Series or pandas.DataFrame
            A Series for one variable, a DataFrame for several, on the index of the declared
            frame.
        """
        return transform(self, "growth", variables, periods=periods, log=log)

    # ------------------------------------------------------------ description

    def xtdescribe(
        self,
        *,
        variables: Sequence[str] | None = None,
        patterns: int = 9,
        width: int = 100,
    ) -> XtDescribe:
        """Structure of the panel and participation patterns (Stata's ``xtdescribe``).

        Parameters
        ----------
        variables
            Count only the observations where all these variables are observed (the
            participation of the variables, not of the rows).
        patterns
            Number of patterns shown (default 9); ``0`` shows none; the rest is grouped.
        width
            Maximum columns of a pattern; beyond it a column represents several periods.

        Returns
        -------
        XtDescribe
            Result of :func:`xtdescribe`.
        """
        return xtdescribe(self, variables=variables, patterns=patterns, width=width)

    def xtsum(self, variables: str | Sequence[str], *, by: str | None = None) -> XtSum:
        """Overall, between and within summary statistics (Stata's ``xtsum``).

        Parameters
        ----------
        variables
            Column names, or Stata-style terms (``"L.x"``, ``"D.x"``).
        by
            Column that defines groups; each group is summarized on its own rows.

        Returns
        -------
        XtSum
            Result of :func:`xtsum`.
        """
        return xtsum(self, variables, by=by)

    def xttab(self, variable: str, *, max_values: int = 1000) -> XtTab:
        """Overall, between and within frequencies of a variable (Stata's ``xttab``).

        Parameters
        ----------
        variable
            A column (or a Stata-style term such as ``"L.x"``); missing values are left out.
        max_values
            Refuse variables with more distinct values than this.

        Returns
        -------
        XtTab
            Result of :func:`xttab`.
        """
        return xttab(self, variable, max_values=max_values)

    def xttrans(
        self,
        variable: str,
        *,
        gaps: Literal["bridge", "strict"] = "bridge",
        bins: int | Sequence[float] | None = None,
        max_values: int = 1000,
    ) -> XtTrans:
        """Transition matrix between consecutive observations (Stata's ``xttrans``).

        Parameters
        ----------
        variable
            A column (or a Stata-style term).  Like Stata, values are compared as ``float32``.
        gaps
            ``"bridge"`` (default, Stata): pair each observation with the unit's next row,
            whatever the calendar gap.  ``"strict"``: only pairs of consecutive periods.
        bins
            Discretize a continuous variable first (an addition of the package, not in Stata):
            an integer ``k`` cuts it at its ``k``-quantiles, a list gives the cut points.
            Values are replaced by the number of the class, ``0`` for the lowest.
        max_values
            Refuse variables with more distinct values than this.

        Returns
        -------
        XtTrans
            Result of :func:`xttrans`.
        """
        return xttrans(self, variable, gaps=gaps, bins=bins, max_values=max_values)

    def xtcorr(
        self,
        variables: Sequence[str],
        *,
        kind: Literal["overall", "within", "between", "all"] = "all",
    ) -> XtCorr:
        """Overall, within and between correlations with the points behind each pair.

        Parameters
        ----------
        variables
            At least two columns or Stata-style terms.
        kind
            ``"overall"``, ``"within"``, ``"between"`` or ``"all"`` (default).

        Returns
        -------
        XtCorr
            Result of :func:`xtcorr`.
        """
        return xtcorr(self, variables, kind=kind)

    def autocorr(self, variable: str, lags: int | Sequence[int] = 1) -> AutoCorr:
        """Autocorrelation of the within deviations at calendar lags.

        Parameters
        ----------
        variable
            Name of a numeric variable; rows where it is missing are ignored.
        lags
            One lag or several, in calendar periods (each at least 1).

        Returns
        -------
        AutoCorr
            Result of :func:`autocorr`: within autocorrelation by lag (``rho`` and ``n`` pairs).
        """
        return autocorr(self, variable, lags)

    def attrition(self) -> Attrition:
        """Entries, exits and re-entries by period, and Kaplan-Meier survival of the units.

        Returns
        -------
        Attrition
            Result of :func:`attrition`.
        """
        return attrition(self)

    def missing_report(self, variables: Sequence[str] | None = None) -> MissingReport:
        """Missing values by variable, period and unit, with joint patterns.

        Parameters
        ----------
        variables
            Columns to examine (default: every column except the two keys).

        Returns
        -------
        MissingReport
            Result of :func:`missing_report`.
        """
        return missing_report(self, variables)

    def complete_cases(self, variables: Sequence[str]) -> CompleteCases:
        """Rows and units kept when all ``variables`` must be observed, and the effect on T_i.

        Parameters
        ----------
        variables
            Variables that must all be observed for a row to be kept.

        Returns
        -------
        CompleteCases
            Result of :func:`complete_cases`.
        """
        return complete_cases(self, variables)

    def gap_types(self) -> GapTypes:
        """Late entry, early exit, internal gap and re-entry of each unit, and gap lengths.

        Returns
        -------
        GapTypes
            Result of :func:`gap_types`.
        """
        return gap_types(self)

    def unbalancedness(self) -> Unbalancedness:
        """Ahrens-Pincus measures of imbalance, gamma and nu (as ``plm::punbalancedness``).

        Returns
        -------
        Unbalancedness
            Result of :func:`unbalancedness`.
        """
        return unbalancedness(self)

    def balance_cost(self) -> BalanceCost:
        """Units and observations kept by each minimum length T_i >= k and each window.

        Returns
        -------
        BalanceCost
            Result of :func:`balance_cost`.
        """
        return balance_cost(self)

    def usable_equations(
        self,
        variable: str,
        *,
        lags: int = 1,
        transform: Transform = "fd",
        exog: Sequence[str] = (),
        instrument_lag: int | None = None,
    ) -> UsableEquations:
        """Equations a dynamic-panel estimator can use, by unit and period (FD or FOD).

        Parameters
        ----------
        variable
            The dependent variable.
        lags
            Number of lags of ``variable`` among the regressors (0 for a static model).
        transform
            ``"fd"`` (first differences) or ``"fod"`` (forward orthogonal deviations).
        exog
            Other regressors that must exist after the transformation.
        instrument_lag
            Oldest lag allowed to serve as an instrument; default 2 for ``"fd"`` and 1 for
            ``"fod"``.

        Returns
        -------
        UsableEquations
            Result of :func:`usable_equations`.
        """
        return usable_equations(
            self,
            variable,
            lags=lags,
            transform=transform,
            exog=exog,
            instrument_lag=instrument_lag,
        )

    def variation(self, variables: Sequence[str] | None = None) -> Variation:
        """Where each variable varies: units and periods that vary, constant or time-invariant.

        Parameters
        ----------
        variables
            Variables to examine; by default every numeric variable other than the keys.

        Returns
        -------
        Variation
            Result of :func:`variation`.
        """
        return variation(self, variables)

    def residual_variation(
        self, variables: Sequence[str], transforms: Sequence[str] = ("within", "fd", "fod")
    ) -> ResidualVariation:
        """Observations, units and variance left by within, first differences and FOD.

        Parameters
        ----------
        variables
            Variables to examine.
        transforms
            Transformations whose residual variation is measured (``"within"``, ``"fd"``,
            ``"fod"``).

        Returns
        -------
        ResidualVariation
            Result of :func:`residual_variation`.
        """
        return residual_variation(self, variables, transforms)

    def spells(self, variable: str, *, max_states: int = 50) -> Spells:
        """Spells of a discrete state: durations, changes per unit, units that change.

        Parameters
        ----------
        variable
            A discrete variable (numeric, with at most ``max_states`` distinct values).
        max_states
            Refuse a variable with more distinct values, which is almost surely continuous.

        Returns
        -------
        Spells
            Result of :func:`spells`.
        """
        return spells(self, variable, max_states=max_states)

    def event_time(
        self, *, indicator: str | None = None, event_period: str | None = None
    ) -> EventTime:
        """Time relative to an event and the periods each unit is observed around it.

        Parameters
        ----------
        indicator
            Name of a 0/1 variable; the event is the first period in which it is non-zero.
        event_period
            Name of a column holding the period of the event (a value of the time variable,
            constant within the unit, missing if the unit is never treated).

        Returns
        -------
        EventTime
            Result of :func:`event_time`.
        """
        return event_time(self, indicator=indicator, event_period=event_period)

    def within_outliers(
        self, variable: str, *, threshold: float = 3.5, min_obs: int = 5
    ) -> WithinOutliers:
        """Within outliers by robust scores (median and MAD in the unit), flagged not removed.

        Parameters
        ----------
        variable
            The variable to examine.
        threshold
            A value is flagged when the absolute modified z-score exceeds it.
        min_obs
            Units with fewer observed values are not scored.

        Returns
        -------
        WithinOutliers
            Result of :func:`within_outliers`.
        """
        return within_outliers(self, variable, threshold=threshold, min_obs=min_obs)

    def level_shifts(
        self,
        variable: str,
        *,
        window: int = 3,
        threshold: float = 5.0,
        influence: float = 3.0,
    ) -> LevelShifts:
        """Level shifts of each unit and the units that weigh most in the within variation.

        Parameters
        ----------
        variable
            The variable to examine.
        window
            Number of observed values on each side of a candidate shift.
        threshold
            A unit is flagged when its largest absolute score exceeds it.
        influence
            A unit is influential when its share of the within sum of squares is above this
            multiple of the average share.

        Returns
        -------
        LevelShifts
            Result of :func:`level_shifts`.
        """
        return level_shifts(self, variable, window=window, threshold=threshold, influence=influence)

    def to_wide(self, variables: Sequence[str] | None = None, *, sep: str = "") -> pd.DataFrame:
        """One row per unit, one column per variable and observed period (``reshape wide``).

        Parameters
        ----------
        variables
            Variables to spread; by default every numeric variable other than the keys.
        sep
            Text between the variable name and the period in the column names (Stata uses none:
            ``x1990``).

        Returns
        -------
        pandas.DataFrame
            One row per unit, one column per variable and observed period.
        """
        return to_wide(self, variables, sep=sep)

    def to_array3d(self, variables: Sequence[str] | None = None) -> Array3D:
        """Array variable x period x unit over the whole calendar, with the mask of rows.

        Parameters
        ----------
        variables
            Variables to lay out; by default every numeric variable other than the keys.

        Returns
        -------
        Array3D
            A panel as an array ``variable x period x unit`` (result of :func:`to_array3d`).
        """
        return to_array3d(self, variables)

    def audit(
        self,
        variables: Sequence[str] | None = None,
        *,
        dependent: str | None = None,
        lags: int = 1,
        thresholds: Mapping[str, float] | None = None,
        checks: Sequence[str] | None = None,
    ) -> AuditReport:
        """Examine the panel in one call: alerts and advice in French and English.

        Parameters
        ----------
        variables
            Variables to examine; by default every numeric variable other than the keys.
        dependent
            Dependent variable of a dynamic model; when given, the audit counts the equations
            left by first differences and by orthogonal deviations with ``lags`` lags.
        lags
            Lags of ``dependent`` in the model.
        thresholds
            Overrides of :class:`~pyxtpanel.explore.checks._context.AuditThresholds`, as a
            mapping (``{"min_units": 50}``) or an instance.
        checks
            Names of the checks to run (default: all of them, see ``CHECKS.names()``).

        Returns
        -------
        AuditReport
            Result of :func:`audit`.
        """
        return audit(
            self,
            variables,
            dependent=dependent,
            lags=lags,
            thresholds=thresholds,
            checks=checks,
        )

    def test(self, name: str, *args: Any, **options: Any) -> TestResult:
        """Run a test on the data by name: ``panel.test("ips", "y", lags=1)``.

        Same as ``run_test(name, panel, ...)`` for the tests on a panel (category A); the tests
        on arrays (category B) are refused with a pointer to ``run_test``.

        Parameters
        ----------
        name
            Name or alias of the test (``"ips"``, ``"im_pesaran_shin"``); required.
        args
            Arguments of the test after the panel (the variable names).
        options
            The options of the test; an option the test does not have is refused with the list
            of those it has.

        Returns
        -------
        TestResult
            The test statistic, its p-value and the quantities behind them.
        """
        return run_panel_test(self, name, *args, **options)

    def unit_root(self, variable: str, *, test: str | None = None, **options: Any) -> TestResult:
        """Test ``variable`` for a unit root (see :func:`pyxtpanel.unit_root`).

        Parameters
        ----------
        variable
            Name of the numeric variable tested.
        test
            Name or alias of the test (required).
        options
            The options of the test (``deterministics``, ``lags``, ``demean``...); an option the
            test does not have is refused with the list of those it has.

        Returns
        -------
        TestResult
            The test statistic, its p-value and the quantities behind them.
        """
        return run_panel_test(self, test, variable, family="unitroot", **options)

    def cross_section_dependence(
        self, variable: str, *, test: str | None = None, **options: Any
    ) -> TestResult:
        """Test ``variable`` for dependence between units (see ``cross_section_dependence``).

        Parameters
        ----------
        variable
            Name of the numeric variable tested.
        test
            Name or alias of the test (required).
        options
            The options of the test.

        Returns
        -------
        TestResult
            The test statistic, its p-value and the quantities behind them.
        """
        return run_panel_test(self, test, variable, family="crosssection", **options)

    def serial_correlation(
        self,
        dependent: str,
        regressors: Sequence[str],
        *,
        test: str | None = None,
        **options: Any,
    ) -> TestResult:
        """Test the errors of a model for serial correlation (see ``serial_correlation``).

        Parameters
        ----------
        dependent
            Names of the dependent variable and of the regressors of the model.
        regressors
            Names of the dependent variable and of the regressors of the model.
        test
            Name or alias of the test (required).
        options
            The options of the test.

        Returns
        -------
        TestResult
            The test statistic, its p-value and the quantities behind them.
        """
        return run_panel_test(self, test, dependent, regressors, family="serial", **options)

    def causality(
        self,
        dependent: str,
        cause: str,
        *,
        test: str | None = None,
        **options: Any,
    ) -> TestResult:
        """Test whether ``cause`` Granger-causes ``dependent`` (see ``causality``).

        Parameters
        ----------
        dependent
            Names of the variable that is caused and of the variable whose lags are tested.
        cause
            Names of the variable that is caused and of the variable whose lags are tested.
        test
            Name or alias of the test (required).
        options
            The options of the test (``lags``, ``statistic``).

        Returns
        -------
        TestResult
            The test statistic, its p-value and the quantities behind them.
        """
        return run_panel_test(self, test, dependent, cause, family="causality", **options)

    def cointegration(
        self,
        dependent: str,
        regressors: Sequence[str],
        *,
        test: str | None = None,
        **options: Any,
    ) -> TestResult:
        """Test a model for cointegration (see ``cointegration``).

        Parameters
        ----------
        dependent
            Names of the dependent variable and of the regressors of the cointegrating relation.
        regressors
            Names of the dependent variable and of the regressors of the cointegrating relation.
        test
            Name or alias of the test (required).
        options
            The options of the test (``statistic``, ``adf_lags``, ``kernel``,
            ``kernel_lags``...).

        Returns
        -------
        TestResult
            The test statistic, its p-value and the quantities behind them.
        """
        return run_panel_test(self, test, dependent, regressors, family="cointegration", **options)

    def slope_heterogeneity(
        self,
        dependent: str,
        regressors: Sequence[str],
        *,
        test: str | None = None,
        **options: Any,
    ) -> TestResult:
        """Test whether the slopes of a model are equal across units (``slope_heterogeneity``).

        Parameters
        ----------
        dependent
            Names of the dependent variable and of the regressors whose slopes are tested.
        regressors
            Names of the dependent variable and of the regressors whose slopes are tested.
        test
            Name or alias of the test (required).
        options
            The options of the test (``partial``, ``constant``, ``ar``, ``crosssectional``,
            ``adjusted``).

        Returns
        -------
        TestResult
            The test statistic, its p-value and the quantities behind them.
        """
        return run_panel_test(self, test, dependent, regressors, family="heterogeneity", **options)

    def plot(self, kind: str, **options: Any) -> Any:
        """Draw a figure (needs the ``plot`` extra).

        ``kind`` is ``"lines"``, ``"scatter"``, ``"distribution"``, ``"participation"`` or
        ``"transition"``; the options are those of the kind.

        Returns the matplotlib ``Figure``; save it with :func:`pyxtpanel.plot.save_figure`.

        Parameters
        ----------
        kind
            ``"lines"``, ``"scatter"``, ``"distribution"``, ``"participation"`` or
            ``"transition"``.
        options
            Options of the kind (``variable``, ``x`` and ``y``, ``lang``...).

        Returns
        -------
        matplotlib.figure.Figure
            The figure; save it with :func:`pyxtpanel.plot.save_figure`.
        """
        return _plot(self, kind, **options)
