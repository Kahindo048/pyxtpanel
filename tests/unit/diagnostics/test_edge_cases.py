"""Degenerate inputs of the residual and specification tests (plan L.8.12).

Whatever the input, a test must either return a result or refuse with a package error that says
what is wrong; it must never leak a NumPy exception, a ``ZeroDivisionError`` or a
``RuntimeWarning``, and it must never return a number for an input that has none (empty, missing
or infinite values, residuals that are all zero).  Warnings are turned into errors here, so a
numerical warning inside a test fails the case.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, PyxtpanelError, SpecificationError
from pyxtpanel.diagnostics import (
    TestResult,
    breusch_pagan,
    cips,
    dumitrescu_hurlin,
    jarque_bera,
    kao,
    lm_random_effects,
    panel_durbin_watson,
    pedroni,
    reset,
    shapiro_francia,
    shapiro_wilk,
    skewness_kurtosis,
    slope_homogeneity,
    vif,
    wald_groupwise,
    westerlund,
    white,
)

N = 60
OFFSETS = np.arange(0, N + 1, 6)

#: test -> call from (residuals or dependent variable, regressors)
CALLS: dict[str, Callable[[np.ndarray, np.ndarray], TestResult]] = {
    "breusch_pagan": lambda e, x: breusch_pagan(e, x),
    "white": lambda e, x: white(e, x),
    "wald_groupwise": lambda e, x: wald_groupwise(e, OFFSETS),
    "lm_random_effects": lambda e, x: lm_random_effects(e, OFFSETS),
    "skewness_kurtosis": lambda e, x: skewness_kurtosis(e),
    "jarque_bera": lambda e, x: jarque_bera(e),
    "shapiro_wilk": lambda e, x: shapiro_wilk(e),
    "shapiro_francia": lambda e, x: shapiro_francia(e),
    "reset": lambda e, x: reset(e, x),
}


def good() -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(0)
    return rng.normal(size=N), rng.normal(size=(N, 2))


def outcome(call: Callable[[], TestResult]) -> TestResult | PyxtpanelError:
    """Run ``call`` with warnings as errors; return its result or the package error it raised."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        try:
            return call()
        except PyxtpanelError as error:
            return error


def is_refusal_or_empty_result(value: TestResult | PyxtpanelError) -> bool:
    """A package error, or a result with no number and a remark saying why."""
    if isinstance(value, PyxtpanelError):
        return True
    return math.isnan(value.statistic) and bool(value.notes)


@pytest.mark.parametrize("name", CALLS)
def test_a_valid_input_gives_a_finite_statistic(name: str) -> None:
    e, x = good()
    result = outcome(lambda: CALLS[name](e, x))
    assert isinstance(result, TestResult) and math.isfinite(result.statistic)


@pytest.mark.parametrize("name", CALLS)
@pytest.mark.parametrize("kind", ["nan", "inf", "empty", "zeros"])
def test_an_input_without_numbers_is_refused_and_never_gets_a_statistic(
    name: str, kind: str
) -> None:
    e, x = good()
    if kind == "nan":
        e = e.copy()
        e[3] = np.nan
    elif kind == "inf":
        e = e.copy()
        e[3] = np.inf
    elif kind == "empty":
        e, x = e[:0], x[:0]
    else:
        e = np.zeros(N)
    result = outcome(lambda: CALLS[name](e, x))
    assert is_refusal_or_empty_result(result), (name, kind, result)


@pytest.mark.parametrize("name", CALLS)
def test_constant_residuals_do_not_leak_a_numerical_warning_or_a_raw_error(name: str) -> None:
    e, x = good()
    outcome(lambda: CALLS[name](np.ones(N), x))  # any package outcome is fine, a leak raises


@pytest.mark.parametrize("name", CALLS)
def test_a_single_observation_is_refused_or_gives_an_empty_result(name: str) -> None:
    e, x = good()
    if name in ("wald_groupwise", "lm_random_effects"):  # one unit with one observation
        test = {"wald_groupwise": wald_groupwise, "lm_random_effects": lm_random_effects}[name]
        result = outcome(lambda: test(e[:1], np.array([0, 1])))
    else:
        result = outcome(lambda: CALLS[name](e[:1], x[:1]))
    assert is_refusal_or_empty_result(result), name


@pytest.mark.parametrize("name", ["breusch_pagan", "white", "reset"])
def test_regressors_with_a_different_length_are_refused(name: str) -> None:
    e, x = good()
    result = outcome(lambda: CALLS[name](e, x[:-1]))
    assert isinstance(result, SpecificationError)


@pytest.mark.parametrize("name", ["breusch_pagan", "white", "reset"])
def test_constant_or_missing_regressors_are_refused(name: str) -> None:
    e, _ = good()
    result = outcome(lambda: CALLS[name](e, np.ones((N, 2))))
    assert isinstance(result, NumericalError), name


def test_collinear_regressors_are_handled_as_stata_does() -> None:
    e, x = good()
    twin = np.column_stack([x, x[:, 0] + x[:, 1]])
    assert outcome(lambda: breusch_pagan(e, twin)).statistic == pytest.approx(  # type: ignore[union-attr]
        breusch_pagan(e, x).statistic, rel=1e-8
    )
    assert math.isfinite(reset(e, twin).statistic)  # the collinear column adds no parameter
    assert vif(x).mean < 1.5
    with pytest.raises(NumericalError, match="linear combination"):
        vif(twin)  # Stata omits the variable; the function names it and asks to drop it


def test_squared_residuals_without_variation_are_refused_by_the_variance_tests() -> None:
    for call in (
        lambda: breusch_pagan(np.ones(N), good()[1]),
        lambda: white(np.ones(N), good()[1]),
    ):
        assert isinstance(outcome(call), NumericalError)


def test_every_unit_without_a_usable_variance_is_refused_by_the_groupwise_test() -> None:
    assert isinstance(outcome(lambda: wald_groupwise(np.ones(N), OFFSETS)), NumericalError)


def test_a_constant_dependent_variable_is_refused_by_reset() -> None:
    _, x = good()
    assert isinstance(outcome(lambda: reset(np.ones(N), x)), NumericalError)


# ---------------------------------------------------------------- tests on a panel


def frame(units: int = 10, periods: int = 8) -> pd.DataFrame:
    rng = np.random.default_rng(1)
    out = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, units + 1), periods),
            "t": np.tile(np.arange(1, periods + 1), units),
        }
    )
    out["x"] = rng.normal(size=len(out))
    out["z"] = rng.normal(size=len(out))
    out["y"] = out["x"] + rng.normal(size=len(out))
    return out


PANEL_CALLS: dict[str, Callable[[xt.PanelData], TestResult]] = {
    "slope_homogeneity": lambda p: slope_homogeneity(p, "y", ["x", "z"]),
    "panel_durbin_watson": lambda p: panel_durbin_watson(p, "y", ["x", "z"]),
    "cips": lambda p: cips(p, "y", lags=0),
    "kao": lambda p: kao(p, "y", ["x", "z"]),
    "westerlund": lambda p: westerlund(p, "y", ["x", "z"]),
    "pedroni": lambda p: pedroni(p, "y", ["x", "z"], kernel_lags=1),
    "dumitrescu_hurlin": lambda p: dumitrescu_hurlin(p, "y", "x", statistic="zbar"),
}


@pytest.mark.parametrize("name", PANEL_CALLS)
def test_a_valid_panel_gives_a_finite_statistic(name: str) -> None:
    result = outcome(lambda: PANEL_CALLS[name](xt.xtset(frame(), "id", "t")))
    assert isinstance(result, TestResult) and math.isfinite(result.statistic)


@pytest.mark.parametrize("name", PANEL_CALLS)
@pytest.mark.parametrize(
    "kind", ["constant_y", "all_missing", "one_period", "one_obs_unit", "collinear", "tiny_t"]
)
def test_a_degenerate_panel_never_leaks_a_raw_error_or_a_warning(name: str, kind: str) -> None:
    f = frame()
    if kind == "constant_y":
        f["y"] = 1.0
    elif kind == "all_missing":
        f["y"] = np.nan
    elif kind == "one_period":
        f = f[f["t"] == 1]
    elif kind == "one_obs_unit":
        f = f[~((f["id"] == 1) & (f["t"] > 1))]
    elif kind == "collinear":
        f["z"] = 2 * f["x"]
    else:
        f = frame(periods=2)
    result = outcome(lambda: PANEL_CALLS[name](xt.xtset(f.reset_index(drop=True), "id", "t")))
    if kind in ("constant_y", "all_missing", "one_period"):
        assert isinstance(result, PyxtpanelError), (name, kind, result)


def test_a_unit_entirely_missing_is_dropped_and_a_single_unit_is_refused_by_the_slope_test() -> (
    None
):
    f = frame()
    f.loc[f["id"] == 1, "y"] = np.nan
    result = slope_homogeneity(xt.xtset(f, "id", "t"), "y", ["x", "z"])
    assert result.n_groups == 9
    one = frame()
    one = one[one["id"] == 1]
    assert isinstance(
        outcome(lambda: slope_homogeneity(xt.xtset(one, "id", "t"), "y", ["x", "z"])), DataError
    )
