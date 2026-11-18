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
from pyxtpanel.ops.frame import apply_operator, eliminate_effects, operate, transform
from pyxtpanel.ops.groups import group_size
from pyxtpanel.ops.within import theta_re
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
        """
        return operate(self, terms)

    def lag(
        self, variables: str | Sequence[str], k: int | Sequence[int] = 1
    ) -> pd.Series | pd.DataFrame:
        """Value ``k`` periods earlier (``L.x``, ``L2.x``...).

        A Series for one variable and one ``k``, otherwise a DataFrame.
        """
        return apply_operator(self, "lag", variables, k)

    def lead(
        self, variables: str | Sequence[str], k: int | Sequence[int] = 1
    ) -> pd.Series | pd.DataFrame:
        """Value ``k`` periods later (``F.x``, ``F2.x``...)."""
        return apply_operator(self, "lead", variables, k)

    def diff(
        self, variables: str | Sequence[str], order: int | Sequence[int] = 1
    ) -> pd.Series | pd.DataFrame:
        """Difference with the previous period (``D.x``); ``order=2`` is ``D2.x``."""
        return apply_operator(self, "diff", variables, order)

    def seasonal(
        self, variables: str | Sequence[str], period: int | Sequence[int] = 1
    ) -> pd.Series | pd.DataFrame:
        """Difference with the value ``period`` periods earlier (``S12.x``)."""
        return apply_operator(self, "seasonal", variables, period)

    # ------------------------------------------------------------ transformations

    def transform(
        self, name: str, variables: str | Sequence[str], *, listwise: bool = False, **params: Any
    ) -> pd.Series | pd.DataFrame:
        """Apply a registered transformation (``"within"``, ``"fod"``, ``"rolling"``...) by name.

        ``variables`` is one name (a Series comes back) or several (a DataFrame with the same
        column names), always on the index of the declared frame.  ``listwise=True`` first drops
        the rows where any of the variables is missing, as Stata's ``xtdata`` does.
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
        """Mean of the unit over its observed values, repeated on the rows with a value."""
        return transform(self, "between", variables, listwise=listwise)

    def quasi_demean(
        self,
        variables: str | Sequence[str],
        theta: float | Sequence[float],
        *,
        listwise: bool = False,
    ) -> pd.Series | pd.DataFrame:
        """``x − θ x̄_i``, ``θ`` a scalar or one value per unit (order of ``unit_labels``)."""
        return transform(
            self, "quasi_demean", variables, listwise=listwise, theta=np.asarray(theta, float)
        )

    def theta_re(self, sigma_u: float, sigma_e: float) -> pd.Series:
        """Random-effects ``θ_i`` of each unit for the given standard deviations."""
        return pd.Series(
            theta_re(sigma_u, sigma_e, self.counts), index=self.unit_labels, name="theta"
        )

    def fod(
        self, variables: str | Sequence[str], *, scale: bool = True
    ) -> pd.Series | pd.DataFrame:
        """Forward orthogonal deviations: later observed rows of the unit, ``√(F/(F+1))`` scale."""
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
        """
        return transform(self, "group_count", variables, by=by, at_missing=at_missing)

    def group_size(self, by: Literal["unit", "time"] = "unit") -> pd.Series:
        """Rows of the unit (``T_i``) or of the period (``N_t``), repeated on each row."""
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
        """Statistic over the ``window`` calendar periods ending at each row (gaps respected)."""
        return transform(self, "rolling", variables, window=window, stat=stat, min_obs=min_obs)

    def cumulative(
        self, variables: str | Sequence[str], stat: Literal["sum", "mean"] = "sum"
    ) -> pd.Series | pd.DataFrame:
        """Return the running sum (or mean) of the unit's observed values."""
        return transform(self, "cumulative", variables, stat=stat)

    def growth(
        self, variables: str | Sequence[str], periods: int = 1, *, log: bool = False
    ) -> pd.Series | pd.DataFrame:
        """Growth rate over ``periods`` calendar periods (relative, or log difference)."""
        return transform(self, "growth", variables, periods=periods, log=log)

    # ------------------------------------------------------------ description

    def xtdescribe(
        self,
        *,
        variables: Sequence[str] | None = None,
        patterns: int = 9,
        width: int = 100,
    ) -> XtDescribe:
        """Structure of the panel and participation patterns (Stata's ``xtdescribe``)."""
        return xtdescribe(self, variables=variables, patterns=patterns, width=width)

    def xtsum(self, variables: str | Sequence[str], *, by: str | None = None) -> XtSum:
        """Overall, between and within summary statistics (Stata's ``xtsum``)."""
        return xtsum(self, variables, by=by)

    def xttab(self, variable: str, *, max_values: int = 1000) -> XtTab:
        """Overall, between and within frequencies of a variable (Stata's ``xttab``)."""
        return xttab(self, variable, max_values=max_values)

    def xttrans(
        self,
        variable: str,
        *,
        gaps: Literal["bridge", "strict"] = "bridge",
        bins: int | Sequence[float] | None = None,
        max_values: int = 1000,
    ) -> XtTrans:
        """Transition matrix between consecutive observations (Stata's ``xttrans``)."""
        return xttrans(self, variable, gaps=gaps, bins=bins, max_values=max_values)

    def xtcorr(
        self,
        variables: Sequence[str],
        *,
        kind: Literal["overall", "within", "between", "all"] = "all",
    ) -> XtCorr:
        """Overall, within and between correlations with the points behind each pair."""
        return xtcorr(self, variables, kind=kind)

    def autocorr(self, variable: str, lags: int | Sequence[int] = 1) -> AutoCorr:
        """Autocorrelation of the within deviations at calendar lags."""
        return autocorr(self, variable, lags)

    def attrition(self) -> Attrition:
        """Entries, exits and re-entries by period, and Kaplan-Meier survival of the units."""
        return attrition(self)

    def missing_report(self, variables: Sequence[str] | None = None) -> MissingReport:
        """Missing values by variable, period and unit, with joint patterns."""
        return missing_report(self, variables)

    def complete_cases(self, variables: Sequence[str]) -> CompleteCases:
        """Rows and units kept when all ``variables`` must be observed, and the effect on T_i."""
        return complete_cases(self, variables)

    def gap_types(self) -> GapTypes:
        """Late entry, early exit, internal gap and re-entry of each unit, and gap lengths."""
        return gap_types(self)

    def unbalancedness(self) -> Unbalancedness:
        """Ahrens-Pincus measures of imbalance, gamma and nu (as ``plm::punbalancedness``)."""
        return unbalancedness(self)

    def balance_cost(self) -> BalanceCost:
        """Units and observations kept by each minimum length T_i >= k and each window."""
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
        """Equations a dynamic-panel estimator can use, by unit and period (FD or FOD)."""
        return usable_equations(
            self,
            variable,
            lags=lags,
            transform=transform,
            exog=exog,
            instrument_lag=instrument_lag,
        )

    def variation(self, variables: Sequence[str] | None = None) -> Variation:
        """Where each variable varies: units and periods that vary, constant or time-invariant."""
        return variation(self, variables)

    def residual_variation(
        self, variables: Sequence[str], transforms: Sequence[str] = ("within", "fd", "fod")
    ) -> ResidualVariation:
        """Observations, units and variance left by within, first differences and FOD."""
        return residual_variation(self, variables, transforms)

    def spells(self, variable: str, *, max_states: int = 50) -> Spells:
        """Spells of a discrete state: durations, changes per unit, units that change."""
        return spells(self, variable, max_states=max_states)

    def event_time(
        self, *, indicator: str | None = None, event_period: str | None = None
    ) -> EventTime:
        """Time relative to an event and the periods each unit is observed around it."""
        return event_time(self, indicator=indicator, event_period=event_period)

    def within_outliers(
        self, variable: str, *, threshold: float = 3.5, min_obs: int = 5
    ) -> WithinOutliers:
        """Within outliers by robust scores (median and MAD in the unit), flagged not removed."""
        return within_outliers(self, variable, threshold=threshold, min_obs=min_obs)

    def level_shifts(
        self,
        variable: str,
        *,
        window: int = 3,
        threshold: float = 5.0,
        influence: float = 3.0,
    ) -> LevelShifts:
        """Level shifts of each unit and the units that weigh most in the within variation."""
        return level_shifts(self, variable, window=window, threshold=threshold, influence=influence)

    def to_wide(self, variables: Sequence[str] | None = None, *, sep: str = "") -> pd.DataFrame:
        """One row per unit, one column per variable and observed period (``reshape wide``)."""
        return to_wide(self, variables, sep=sep)

    def to_array3d(self, variables: Sequence[str] | None = None) -> Array3D:
        """Array variable x period x unit over the whole calendar, with the mask of rows."""
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
        """Examine the panel in one call: alerts and advice in French and English."""
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
        """
        return run_panel_test(self, name, *args, **options)

    def unit_root(self, variable: str, *, test: str | None = None, **options: Any) -> TestResult:
        """Test ``variable`` for a unit root (see :func:`pyxtpanel.unit_root`)."""
        return run_panel_test(self, test, variable, family="unitroot", **options)

    def cross_section_dependence(
        self, variable: str, *, test: str | None = None, **options: Any
    ) -> TestResult:
        """Test ``variable`` for dependence between units (see ``cross_section_dependence``)."""
        return run_panel_test(self, test, variable, family="crosssection", **options)

    def serial_correlation(
        self,
        dependent: str,
        regressors: Sequence[str],
        *,
        test: str | None = None,
        **options: Any,
    ) -> TestResult:
        """Test the errors of a model for serial correlation (see ``serial_correlation``)."""
        return run_panel_test(self, test, dependent, regressors, family="serial", **options)

    def causality(
        self,
        dependent: str,
        cause: str,
        *,
        test: str | None = None,
        **options: Any,
    ) -> TestResult:
        """Test whether ``cause`` Granger-causes ``dependent`` (see ``causality``)."""
        return run_panel_test(self, test, dependent, cause, family="causality", **options)

    def cointegration(
        self,
        dependent: str,
        regressors: Sequence[str],
        *,
        test: str | None = None,
        **options: Any,
    ) -> TestResult:
        """Test a model for cointegration (see ``cointegration``)."""
        return run_panel_test(self, test, dependent, regressors, family="cointegration", **options)

    def slope_heterogeneity(
        self,
        dependent: str,
        regressors: Sequence[str],
        *,
        test: str | None = None,
        **options: Any,
    ) -> TestResult:
        """Test whether the slopes of a model are equal across units (``slope_heterogeneity``)."""
        return run_panel_test(self, test, dependent, regressors, family="heterogeneity", **options)

    def plot(self, kind: str, **options: Any) -> Any:
        """Draw a figure (needs the ``plot`` extra).

        ``kind`` is ``"lines"``, ``"scatter"``, ``"distribution"``, ``"participation"`` or
        ``"transition"``; the options are those of the kind.

        Returns the matplotlib ``Figure``; save it with :func:`pyxtpanel.plot.save_figure`.
        """
        return _plot(self, kind, **options)
