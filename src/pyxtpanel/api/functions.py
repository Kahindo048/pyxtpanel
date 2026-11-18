"""Top-level functions of the package (layer 6)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.api._lazy import run_panel_test
from pyxtpanel.api.panel_data import PanelData
from pyxtpanel.describe.attrition import attrition
from pyxtpanel.describe.corr import autocorr, xtcorr
from pyxtpanel.describe.missing import complete_cases, missing_report
from pyxtpanel.describe.structure import xtdescribe
from pyxtpanel.describe.summary import xtsum
from pyxtpanel.describe.tab import xttab
from pyxtpanel.describe.trans import xttrans
from pyxtpanel.explore.audit import AuditReport, duplicates_report
from pyxtpanel.explore.audit import audit as audit_panel
from pyxtpanel.explore.balance import balance_cost, unbalancedness
from pyxtpanel.explore.equations import usable_equations
from pyxtpanel.explore.event_time import event_time
from pyxtpanel.explore.gaps import gap_types
from pyxtpanel.explore.outliers import level_shifts, within_outliers
from pyxtpanel.explore.reshape import resolve_duplicates, to_array3d, to_long, to_wide
from pyxtpanel.explore.spells import spells
from pyxtpanel.explore.variation import residual_variation, variation
from pyxtpanel.ops.frame import operate, transform
from pyxtpanel.panel.base import PanelBase
from pyxtpanel.panel.declare import MissingKeys
from pyxtpanel.panel.gaps import find_duplicates as duplicates
from pyxtpanel.panel.io import read_stata, read_table

if TYPE_CHECKING:
    from pyxtpanel.diagnostics import TestResult

__all__ = [
    "PanelLike",
    "as_panel",
    "attrition",
    "audit",
    "autocorr",
    "balance_cost",
    "complete_cases",
    "cross_section_dependence",
    "duplicates",
    "event_time",
    "gap_types",
    "level_shifts",
    "missing_report",
    "operate",
    "read_stata",
    "read_table",
    "residual_variation",
    "resolve_duplicates",
    "serial_correlation",
    "causality",
    "cointegration",
    "slope_heterogeneity",
    "spells",
    "to_array3d",
    "to_long",
    "to_wide",
    "transform",
    "unbalancedness",
    "unit_root",
    "usable_equations",
    "variation",
    "within_outliers",
    "xtcorr",
    "xtdescribe",
    "xtset",
    "xtsum",
    "xttab",
    "xttrans",
]

_i18n.register(
    "api",
    {
        "api.as_panel_needs_keys": {
            "en": "as_panel received a DataFrame without a two-level index: pass id= and time=.",
            "fr": "as_panel a reçu un DataFrame sans index à deux niveaux : passez id= et time=.",
        },
        "api.as_panel_keys": {
            "en": "the panel is declared by ({id!r}, {time!r}), not ({got_id!r}, {got_time!r}); "
            "declare a new panel with xtset.",
            "fr": "le panel est déclaré par ({id!r}, {time!r}), pas par "
            "({got_id!r}, {got_time!r}) ; déclarez un nouveau panel avec xtset.",
        },
        "api.as_panel_type": {
            "en": "as_panel expects a PanelData or a DataFrame, got {kind}.",
            "fr": "as_panel attend un PanelData ou un DataFrame, reçu {kind}.",
        },
    },
)


@runtime_checkable
class PanelLike(Protocol):
    """What an estimation package may rely on in a declared panel.

    Any :class:`PanelData` satisfies it.  Estimators accept a ``PanelLike`` (or call
    :func:`as_panel` on their input) and read the data through :meth:`arrays`.
    """

    id_name: str
    time_name: str

    @property
    def n_units(self) -> int:
        """Number of units."""
        ...

    @property
    def n_obs(self) -> int:
        """Number of observations."""
        ...

    def arrays(self, columns: Any, *, where: Any = None) -> Any:
        """Return the array contract (see :class:`pyxtpanel.contract.PanelArrays`)."""
        ...


def as_panel(
    obj: PanelBase | pd.DataFrame,
    id: str | None = None,  # noqa: A002
    time: str | None = None,
    **options: Any,
) -> PanelData:
    """Return ``obj`` as a :class:`PanelData`: the entry point of the estimation packages.

    Parameters
    ----------
    obj
        A declared panel (returned as it is, as a ``PanelData``), or a DataFrame.
    id, time
        Names of the identifier and the time variable for a DataFrame; leave both out when the
        frame has a two-level index (unit, period).  For a declared panel they must be omitted
        or equal to its own.
    **options
        Options of :func:`xtset` (``delta``, ``freq``, ``format``, ``missing_keys``).

    Returns
    -------
    PanelData
        A DataFrame declared as a panel, as Stata's ``xtset`` does.

    Raises
    ------
    SpecificationError
        If a DataFrame has no way to find its keys, or the keys given contradict the panel's,
        or ``obj`` is neither a panel nor a DataFrame.
    """
    if isinstance(obj, PanelBase):
        if (id is not None and id != obj.id_name) or (time is not None and time != obj.time_name):
            raise SpecificationError(
                _i18n.message(
                    "api.as_panel_keys",
                    id=obj.id_name,
                    time=obj.time_name,
                    got_id=id,
                    got_time=time,
                ),
                key="api.as_panel_keys",
            )
        if isinstance(obj, PanelData):
            return obj
        return PanelData(
            frame=obj.frame,
            id_name=obj.id_name,
            time_name=obj.time_name,
            axis=obj.axis,
            index=obj.index,
            n_dropped=obj.n_dropped,
        )
    if not isinstance(obj, pd.DataFrame):
        raise SpecificationError(
            _i18n.message("api.as_panel_type", kind=type(obj).__name__), key="api.as_panel_type"
        )
    if id is None and time is None and obj.index.nlevels != 2:
        raise SpecificationError(
            _i18n.message("api.as_panel_needs_keys"), key="api.as_panel_needs_keys"
        )
    return xtset(obj, id, time, **options)


def xtset(
    frame: pd.DataFrame,
    id: str | None = None,  # noqa: A002 - the name Stata users expect
    time: str | None = None,
    *,
    delta: int | None = None,
    freq: str | None = None,
    format: str | None = None,  # noqa: A002
    missing_keys: MissingKeys = "error",
) -> PanelData:
    """Declare a DataFrame as a panel (Stata's ``xtset id time, delta()``).

    The declaration is validated once: repeated (unit, period) pairs, missing keys and a time
    variable that is not on a regular calendar are refused with a message that names the rows.
    The DataFrame is never modified.

    Parameters
    ----------
    frame
        Long data, one row per (unit, period).
    id, time
        Identifier and time variable, columns or index levels.  Leave both out to use the two
        levels of a MultiIndex, in that order.
    delta
        Calendar units between two consecutive periods (default 1, as in Stata).
    freq
        ``"D"``, ``"W"``, ``"M"``, ``"Q"``, ``"H"`` or ``"Y"`` for dates; inferred when omitted.
    format
        Stata ``%t`` format (``"%tm"``, ``"%tq"``...) when the time holds Stata-coded integers.
    missing_keys
        ``"error"`` (default) refuses a missing identifier or time; ``"drop"`` ignores those
        observations as Stata does and counts them.

    Returns
    -------
    PanelData
        The declared panel.

    Raises
    ------
    DeclarationError
        Missing keys (Stata code 459), repeated pairs (451), no usable observation (2000).
    TimeAxisError
        Time that is not an integer, a date or a period, or is off the calendar grid.
    SpecificationError
        Invalid option, unknown variable, or only one of ``id`` and ``time`` given.

    Examples
    --------
    >>> import pandas as pd
    >>> import pyxtpanel as xt
    >>> data = pd.DataFrame({"firm": [1, 1, 2], "year": [2000, 2001, 2001], "y": [1.0, 2, 3]})
    >>> print(xt.xtset(data, "firm", "year"))
    PanelData: 2 units, 3 observations
    Panel variable: firm (unbalanced)
     Time variable: year, 2000 to 2001
             Delta: 1 unit
    """
    return PanelData.from_frame(
        frame, id, time, delta=delta, freq=freq, format=format, missing_keys=missing_keys
    )


def audit(
    data: PanelBase | pd.DataFrame,
    id: str | None = None,  # noqa: A002
    time: str | None = None,
    *,
    variables: Sequence[str] | None = None,
    dependent: str | None = None,
    lags: int = 1,
    thresholds: Mapping[str, float] | None = None,
    checks: Sequence[str] | None = None,
    **options: Any,
) -> AuditReport:
    """Examine a panel (or raw data) in one call and report what to look at before estimating.

    Raw data whose (unit, period) pairs repeat cannot be declared: the report then holds a
    critical alert saying how many pairs repeat and how to merge them.

    Parameters
    ----------
    data
        A declared panel, or a DataFrame (then ``id`` and ``time`` name its keys, as in
        :func:`xtset`, and ``**options`` are the options of :func:`xtset`).
    variables, dependent, lags, thresholds, checks
        See :func:`pyxtpanel.explore.audit.audit`.
    id
        Column that identifies the unit, when ``data`` is a DataFrame.
    time
        Column that holds the period, when ``data`` is a DataFrame.
    options
        Options of :func:`xtset` (``delta``, ``freq``, ``format``, ``missing_keys``), used when
        ``data`` is a DataFrame.

    Returns
    -------
    AuditReport
        Result of :func:`audit`.
    """
    if isinstance(data, pd.DataFrame):
        repeated = duplicates(data, id, time)
        if len(repeated):
            pairs = float((1.0 / repeated["n_repeats"]).sum())
            return duplicates_report(len(repeated), round(pairs))
    return audit_panel(
        as_panel(data, id, time, **options),
        variables,
        dependent=dependent,
        lags=lags,
        thresholds=thresholds,
        checks=checks,
    )


def unit_root(
    panel: PanelBase, variable: str, *, test: str | None = None, **options: Any
) -> TestResult:
    """Test a variable for a unit root in every unit of a declared panel.

    The tests are Harris-Tzavalis (``"harris_tzavalis"`` or ``"ht"``), Levin-Lin-Chu
    (``"levin_lin_chu"`` or ``"llc"``), Breitung (``"breitung"``), Im-Pesaran-Shin
    (``"im_pesaran_shin"`` or ``"ips"``), Fisher (``"fisher"``, with ADF or Phillips-Perron
    p-values), Pesaran's CIPS (``"cips"``, valid when the units are correlated through common
    factors) and Hadri (``"hadri"``, whose null is stationarity, the opposite of the others).
    There is no default: they answer different questions (a common autoregressive parameter or
    one per unit, balanced panels or not, a null of unit root or of stationarity).

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of the numeric variable tested.
    test
        Name or alias of the test (required).
    **options
        The options of the test (``deterministics``, ``lags``, ``demean``...); an option the test
        does not have is refused with the list of those it has.

    Returns
    -------
    TestResult
        The test statistic, its p-value and the quantities behind them.
    """
    return run_panel_test(panel, test, variable, family="unitroot", **options)


def cross_section_dependence(
    panel: PanelBase, variable: str, *, test: str | None = None, **options: Any
) -> TestResult:
    """Test a variable (or residuals added to the frame) for dependence between units.

    The tests are Pesaran's CD (``"pesaran_cd"`` or ``"cd"``), the Breusch-Pagan LM
    (``"lm_cross_section"``), Friedman's (``"friedman"``) and Frees' (``"frees"``).  They answer
    different questions (many units or few, normal errors or not), so there is no default test.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of the numeric variable tested.
    test
        Name or alias of the test (required).
    **options
        The options of the test.

    Returns
    -------
    TestResult
        The test statistic, its p-value and the quantities behind them.
    """
    return run_panel_test(panel, test, variable, family="crosssection", **options)


def serial_correlation(
    panel: PanelBase,
    dependent: str,
    regressors: Sequence[str],
    *,
    test: str | None = None,
    **options: Any,
) -> TestResult:
    """Test the errors of a model for first-order serial correlation, from the data.

    The tests are Wooldridge's (``"wooldridge"`` or ``"xtserial"``) and the panel Durbin-Watson
    and Baltagi-Wu LBI statistics (``"panel_durbin_watson"`` or ``"xtregar"``, with no p-value).
    There is no default test.

    Parameters
    ----------
    panel
        The declared panel.
    dependent, regressors
        Names of the dependent variable and of the regressors of the model.
    test
        Name or alias of the test (required).
    **options
        The options of the test.

    Returns
    -------
    TestResult
        The test statistic, its p-value and the quantities behind them.
    """
    return run_panel_test(panel, test, dependent, regressors, family="serial", **options)


def causality(
    panel: PanelBase,
    dependent: str,
    cause: str,
    *,
    test: str | None = None,
    **options: Any,
) -> TestResult:
    """Test whether ``cause`` Granger-causes ``dependent`` in the units of a panel.

    The test is Dumitrescu and Hurlin's (``"dumitrescu_hurlin"``, ``"granger"`` or
    ``"xtgcause"``).  There is no default test.

    Parameters
    ----------
    panel
        The declared panel.
    dependent, cause
        Names of the variable that is caused and of the variable whose lags are tested.
    test
        Name or alias of the test (required).
    **options
        The options of the test (``lags``, ``statistic``).

    Returns
    -------
    TestResult
        The test statistic, its p-value and the quantities behind them.
    """
    return run_panel_test(panel, test, dependent, cause, family="causality", **options)


def cointegration(
    panel: PanelBase,
    dependent: str,
    regressors: Sequence[str],
    *,
    test: str | None = None,
    **options: Any,
) -> TestResult:
    """Test a model for cointegration (no cointegration under the null), from the data.

    The test is Kao's (``"kao"``): residual-based, with the same cointegrating vector in every
    unit.  There is no default test.

    Parameters
    ----------
    panel
        The declared panel.
    dependent, regressors
        Names of the dependent variable and of the regressors of the cointegrating relation.
    test
        Name or alias of the test (required).
    **options
        The options of the test (``statistic``, ``adf_lags``, ``kernel``, ``kernel_lags``...).

    Returns
    -------
    TestResult
        The test statistic, its p-value and the quantities behind them.
    """
    return run_panel_test(panel, test, dependent, regressors, family="cointegration", **options)


def slope_heterogeneity(
    panel: PanelBase,
    dependent: str,
    regressors: Sequence[str],
    *,
    test: str | None = None,
    **options: Any,
) -> TestResult:
    """Test whether the slopes of a model are the same in every unit, from the data.

    The test is Pesaran and Yamagata's (``"slope_homogeneity"`` or ``"xthst"``), the standardised
    form of Swamy's.  There is no default test.

    Parameters
    ----------
    panel
        The declared panel.
    dependent, regressors
        Names of the dependent variable and of the regressors whose slopes are tested.
    test
        Name or alias of the test (required).
    **options
        The options of the test (``partial``, ``constant``, ``ar``, ``crosssectional``,
        ``adjusted``).

    Returns
    -------
    TestResult
        The test statistic, its p-value and the quantities behind them.
    """
    return run_panel_test(panel, test, dependent, regressors, family="heterogeneity", **options)
