"""Degenerate panels meet every test on a panel with a result or a clean refusal (plan L.7.6).

For each category A test and each degenerate panel (short T, small N, constant series, collinear
regressors, a unit entirely missing, holes, infinite values) the test must either return a result
with a finite statistic or refuse with a ``PyxtpanelError`` that carries a catalogue key and says
what is wrong.  It must never leak a NumPy or ``ZeroDivisionError`` exception, a numerical warning
or a NaN statistic.  Some refusals are pinned: a panel that carries no information (one period,
nothing observed, a constant series) never gets a number, with the one exception Stata itself
gives for the Fisher test.  Numerical warnings are turned into errors here; the package's own
warnings (an infinite value read as missing) are the documented way to say what was done.
"""

from __future__ import annotations

import contextlib
import math
import warnings
from collections.abc import Callable

import numpy as np
import pandas as pd
import pytest
from tests.contract.test_diagnostics_panels import CATEGORY_A, RECIPES

import pyxtpanel as xt
from pyxtpanel import _i18n
from pyxtpanel._errors import PyxtpanelError, PyxtpanelWarning
from pyxtpanel.diagnostics import TestResult, run_test

#: Tests that read the regressor ``x`` of the recipe (the others look at ``y`` alone).
USES_X = ("dumitrescu_hurlin", "kao", "panel_durbin_watson", "pedroni", "slope_homogeneity")
USES_X += ("westerlund", "wooldridge")


def frame(units: int = 12, periods: int = 20) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, units + 1), periods),
            "t": np.tile(np.arange(1, periods + 1), units),
            "y": rng.normal(size=(units, periods)).cumsum(axis=1).ravel(),
            "x": rng.normal(size=units * periods),
        }
    )


def _hole(f: pd.DataFrame) -> pd.DataFrame:
    return f[~((f["id"] == 1) & f["t"].between(8, 10))]


def _one_observation_unit(f: pd.DataFrame) -> pd.DataFrame:
    return f[~((f["id"] == 1) & (f["t"] > 1))]


def _assign(column: str, rows: str | None, value: float) -> Callable[[pd.DataFrame], pd.DataFrame]:
    def apply(f: pd.DataFrame) -> pd.DataFrame:
        mask = f["id"] == 1 if rows == "unit1" else f.index == 3 if rows == "row3" else f.index >= 0
        f.loc[mask, column] = value
        return f

    return apply


def _collinear(f: pd.DataFrame) -> pd.DataFrame:
    f["x"] = 2.0 * f["y"] + 1.0
    return f


def _scattered_missing(f: pd.DataFrame) -> pd.DataFrame:
    f.loc[f.sample(frac=0.05, random_state=1).index, "y"] = np.nan
    return f


#: Scenarios that edit the 12 x 20 panel; the shape scenarios (``T=``, ``N=``) rebuild it.
EDITS: dict[str, Callable[[pd.DataFrame], pd.DataFrame]] = {
    "constant_y": _assign("y", None, 3.0),
    "constant_unit": _assign("y", "unit1", 5.0),
    "constant_x": _assign("x", None, 1.0),
    "collinear_x": _collinear,
    "missing_unit_y": _assign("y", "unit1", np.nan),
    "all_missing_y": _assign("y", None, np.nan),
    "hole": _hole,
    "scattered_missing": _scattered_missing,
    "one_observation_unit": _one_observation_unit,
    "infinite": _assign("y", "row3", np.inf),
}


def degenerate(kind: str) -> pd.DataFrame:
    """Return the panel of one scenario (12 units x 20 periods unless the scenario says so)."""
    if kind.startswith("T="):
        return frame(periods=int(kind[2:]))
    if kind.startswith("N="):
        return frame(units=int(kind[2:]))
    return EDITS[kind](frame()).reset_index(drop=True)


SCENARIOS = (
    ["T=1", "T=2", "T=3", "T=5", "N=1", "N=2", "N=3", "constant_y", "constant_unit"]
    + ["constant_x", "collinear_x", "missing_unit_y", "all_missing_y", "hole"]
    + ["scattered_missing", "one_observation_unit", "infinite"]
)


def attempt(name: str, f: pd.DataFrame) -> TestResult | PyxtpanelError:
    """Run a test on a frame: its result, or the refusal it raises; anything else propagates."""
    args, options = RECIPES[name]
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # a numerical warning is a defect...
        warnings.simplefilter("ignore", PyxtpanelWarning)  # ...the package's own are not
        try:
            return run_test(name, xt.xtset(f, "id", "t"), *args, **options)
        except PyxtpanelError as refusal:
            return refusal


@pytest.mark.parametrize("scenario", SCENARIOS)
@pytest.mark.parametrize("name", CATEGORY_A)
def test_a_degenerate_panel_gives_a_result_or_a_clean_refusal(name: str, scenario: str) -> None:
    outcome = attempt(name, degenerate(scenario))
    if isinstance(outcome, PyxtpanelError):
        assert outcome.key is not None and _i18n.has_language(outcome.key, "fr"), outcome
        assert str(outcome), "a refusal must say what is wrong"
    else:
        assert math.isfinite(outcome.statistic), "a NaN statistic must be a refusal instead"
        assert outcome.n_groups is None or outcome.n_groups >= 1


@pytest.mark.parametrize("scenario", ["T=1", "all_missing_y"])
@pytest.mark.parametrize("name", CATEGORY_A)
def test_a_panel_with_nothing_to_estimate_is_always_refused(name: str, scenario: str) -> None:
    assert isinstance(attempt(name, degenerate(scenario)), PyxtpanelError)


@pytest.mark.parametrize("name", [n for n in CATEGORY_A if n != "fisher"])
def test_a_constant_series_never_gets_a_number(name: str) -> None:
    # Stata prints a p-value of 1 for the Fisher test of a constant series, and so does the test
    assert isinstance(attempt(name, degenerate("constant_y")), PyxtpanelError)


@pytest.mark.parametrize("name", USES_X)
def test_a_regressor_collinear_with_the_dependent_variable_is_refused(name: str) -> None:
    # a perfect fit leaves rounding noise as residuals, which no statistic should be read from
    assert isinstance(attempt(name, degenerate("collinear_x")), PyxtpanelError)


@pytest.mark.parametrize("name", [n for n in CATEGORY_A if n not in ("panel_durbin_watson",)])
def test_a_single_unit_is_refused(name: str) -> None:
    assert isinstance(attempt(name, degenerate("N=1")), PyxtpanelError)


def test_the_durbin_watson_statistic_of_a_single_unit_is_a_time_series_statistic() -> None:
    outcome = attempt("panel_durbin_watson", degenerate("N=1"))
    assert isinstance(outcome, TestResult) and outcome.n_groups == 1


@pytest.mark.parametrize("name", CATEGORY_A)
def test_a_unit_entirely_missing_is_left_out_as_if_it_were_not_there(name: str) -> None:
    full, without = degenerate("missing_unit_y"), degenerate("missing_unit_y")
    without = without[without["id"] != 1].reset_index(drop=True)
    with_missing, left_out = attempt(name, full), attempt(name, without)
    assert isinstance(with_missing, TestResult) and isinstance(left_out, TestResult)
    assert with_missing.n_groups == left_out.n_groups == 11
    assert with_missing.n_obs == left_out.n_obs
    assert with_missing.statistic == pytest.approx(left_out.statistic, rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("name", USES_X)
def test_a_unit_missing_its_regressor_is_left_out_too(name: str) -> None:
    f = frame()
    f.loc[f["id"] == 1, "x"] = np.nan
    left_out = f[f["id"] != 1].reset_index(drop=True)
    with_missing, without = attempt(name, f), attempt(name, left_out)
    assert isinstance(with_missing, TestResult) and isinstance(without, TestResult)
    assert with_missing.n_groups == without.n_groups == 11
    assert with_missing.statistic == pytest.approx(without.statistic, rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("name", CATEGORY_A)
def test_an_infinite_value_is_read_as_missing_and_said_so(name: str) -> None:
    args, options = RECIPES[name]
    panel = xt.xtset(degenerate("infinite"), "id", "t")
    # a balanced-only test refuses the panel the missing value leaves unbalanced
    with pytest.warns(PyxtpanelWarning, match="infinite"), contextlib.suppress(PyxtpanelError):
        run_test(name, panel, *args, **options)


# --- category B: arrays that carry no information are refused by every test that reads them ---

_NAN, _INF = float("nan"), float("inf")
_M, _W = np.array([1.0, 2.0, 3.0, 4.0]), np.eye(4) * 0.5
_B, _V = np.array([1.0, 2.0]), np.array([[1.0, 0.2], [0.2, 2.0]])
_U = np.random.default_rng(0).normal(size=40)
_OFFSETS = np.array([0, 10, 20, 30, 40])


def _call(name: str, *args: object, **options: object) -> Callable[[], TestResult]:
    return lambda: run_test(name, *args, **options)


#: label -> the call that must be refused with a catalogue message
B_REFUSED: dict[str, Callable[[], TestResult]] = {
    "hansen: missing moment": _call("hansen", np.array([_NAN, 1, 1, 1]), _W, 2),
    "hansen: infinite moment": _call("hansen", np.array([_INF, 1, 1, 1]), _W, 2),
    "hansen: no moment": _call("hansen", np.array([]), np.zeros((0, 0)), 0),
    "hansen: zero weight": _call("hansen", _M, np.zeros((4, 4)), 2),
    "hansen: negative weight": _call("hansen", _M, -np.eye(4), 2),
    "hansen: weight of another size": _call("hansen", _M, np.eye(3), 2),
    "hansen: missing weight": _call("hansen", _M, np.full((4, 4), _NAN), 2),
    "hansen: no group": _call("hansen", _M, _W, 2, n_groups=0),
    "sargan: no scale": _call("sargan", _M, _W, 2, scale=0.0),
    "sargan: missing scale": _call("sargan", _M, _W, 2, scale=_NAN),
    "diff_hansen: missing J": _call("diff_hansen", _NAN, 6.0, 2),
    "diff_hansen: infinite J": _call("diff_hansen", _INF, 6.0, 2),
    "diff_hansen: nothing dropped": _call("diff_hansen", 10.0, 6.0, 0),
    "hausman: missing coefficient": _call("hausman", np.array([_NAN, 1]), _B, _V, _V),
    "hausman: coefficients of another size": _call("hausman", _B, np.array([1.0]), _V, _V),
    "hausman: missing covariance": _call("hausman", _B, _B, np.full((2, 2), _NAN), _V),
    "hausman: kept coefficient out of range": _call("hausman", _B, _B, _V, _V, keep=[5]),
    "hausman: kept coefficient repeated": _call("hausman", _B, _B, _V, _V, keep=[0, 0]),
    "hausman: nothing kept": _call("hausman", _B, _B, _V, _V, keep=[]),
    "hausman: no coefficient": _call("hausman", np.array([]), np.array([]), _V[:0, :0], _V[:0, :0]),
    "wald: zero covariance": _call("wald", _B, np.zeros((2, 2))),
    "wald: missing coefficient": _call("wald", np.array([_NAN, 1]), _V),
    "wald: missing covariance": _call("wald", _B, np.full((2, 2), _NAN)),
    "wald: restriction of another size": _call("wald", _B, _V, np.ones((1, 3))),
    "wald: null restriction": _call("wald", _B, _V, np.zeros((1, 2))),
    "wald: repeated restriction": _call("wald", _B, _V, np.array([[1.0, 0], [1.0, 0]])),
    "wald: no residual degrees of freedom": _call("wald", _B, _V, df_resid=0.0),
    "wald: no coefficient": _call("wald", np.array([]), np.zeros((0, 0))),
    "f_effects: perfect fit": _call("f_effects", 10.0, 0.0, 100, 10, 3),
    "f_effects: pooled fit better": _call("f_effects", 8.0, 10.0, 100, 10, 3),
    "f_effects: one group": _call("f_effects", 10.0, 8.0, 100, 1, 3),
    "f_effects: no degree of freedom": _call("f_effects", 10.0, 8.0, 13, 10, 3),
    "f_effects: missing sum of squares": _call("f_effects", _NAN, 8.0, 100, 10, 3),
    "lm_random_effects: offsets too short": _call(
        "lm_random_effects", _U, np.array([0, 10, 40])[:2]
    ),
    "lm_random_effects: one observation per unit": _call("lm_random_effects", _U[:4], np.arange(5)),
    "baltagi_li: all residuals zero": _call("baltagi_li", np.zeros(40), _OFFSETS),
    "baltagi_li: offsets of another length": _call("baltagi_li", _U, _OFFSETS[:-1]),
    "baltagi_li: unit of two observations": _call(
        "baltagi_li", np.array([1.0, 2.0, -1.0, -2.0]), np.array([0, 2, 4])
    ),
    "baltagi_li: one observation per unit": _call(
        "baltagi_li", np.array([1.0, 2.0, -1.0]), np.arange(4), statistic="serial"
    ),
    "baltagi_li: unit with a gap": _call(
        "baltagi_li", _U[:7], np.array([0, 4, 7]), time_index=np.array([1, 2, 4, 5, 1, 2, 3])
    ),
    "wald_groupwise: offsets not ending at n": _call("wald_groupwise", _U, _OFFSETS[:-1]),
}


@pytest.mark.parametrize("label", list(B_REFUSED))
def test_degenerate_arrays_are_refused_with_a_catalogue_message(label: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(PyxtpanelError) as refusal:
            B_REFUSED[label]()
    assert refusal.value.key is not None and _i18n.has_language(refusal.value.key, "fr")


def test_a_test_that_is_not_identified_says_so_instead_of_giving_a_number() -> None:
    for name, args, options in (
        ("hansen", (_M, _W, 4), {}),
        ("sargan", (_M, _W, 4), {"scale": 1.0}),
    ):
        result = run_test(name, *args, **options)
        assert math.isnan(result.statistic) and math.isnan(result.pvalue)
        assert [note.key for note in result.notes] == ["diag.overid.not_overidentified"]
