"""Fisher-type unit-root test and the ADF brick: Stata's dfuller values, oracle and refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import fisher, run_test
from pyxtpanel.diagnostics._shared.adf import adf_fit, adf_pvalue
from pyxtpanel.diagnostics.unitroot.fisher import _combine

PANEL = pd.read_csv(
    __import__("pathlib").Path(__file__).resolve().parents[2] / "golden" / "data" / "ur_mixed.csv"
)


def unit(i: int) -> np.ndarray:
    return PANEL.loc[PANEL["id"] == i, "y"].to_numpy()


def test_the_adf_statistic_and_p_value_are_those_of_dfuller() -> None:
    # Values printed by Stata 17: dfuller y if id == 1, lags(1) [trend | drift] on sim_ur_mixed.
    fit = adf_fit(unit(1), 1, "constant")
    assert fit is not None and fit.n == 18 and fit.lags == 1
    assert fit.tau == pytest.approx(-2.130130106729879, rel=1e-12)
    assert adf_pvalue(fit, "constant") == pytest.approx(0.232572959312111, rel=1e-12)
    trend = adf_fit(unit(1), 1, "trend")
    assert trend is not None and trend.tau == pytest.approx(-4.651, abs=5e-4)
    assert adf_pvalue(trend, "trend") == pytest.approx(0.0008, abs=5e-5)
    drift = adf_fit(unit(1), 1, "drift")
    assert drift is not None and adf_pvalue(drift, "drift") == pytest.approx(0.0251, abs=5e-5)


def test_the_adf_regression_matches_an_explicit_least_squares() -> None:
    y = unit(2)
    dy = np.diff(y)
    x = np.column_stack([y[1:-1], dy[:-1], np.ones(dy.size - 1)])
    target = dy[1:]
    beta = np.linalg.lstsq(x, target, rcond=None)[0]
    resid = target - x @ beta
    se = math.sqrt(resid @ resid / (len(target) - 3) * np.linalg.inv(x.T @ x)[0, 0])
    fit = adf_fit(y, 1, "constant")
    assert fit is not None and fit.tau == pytest.approx(beta[0] / se, rel=1e-10)
    assert fit.coefficient == pytest.approx(beta[0], rel=1e-10) and fit.df_resid == len(target) - 3


def test_a_hole_removes_the_equations_it_touches() -> None:
    y = unit(3)
    holed = y.copy()
    holed[7] = np.nan
    fit, base = adf_fit(holed, 0, "constant"), adf_fit(y, 0, "constant")
    assert fit is not None and base is not None and fit.n == base.n - 2


def test_degenerate_regressions() -> None:
    assert adf_fit(np.array([1.0, 2.0]), 1, "constant") is None  # no equation for one lag
    assert adf_fit(np.array([1.0, 2.0, 4.0]), 1, "constant") is None  # one equation, 3 coefficients
    saturated = adf_fit(np.array([1.0, 2.0, 4.0]), 0, "constant")  # 2 equations, 2 coefficients
    assert saturated is not None and math.isnan(saturated.tau)
    assert adf_pvalue(saturated, "constant") == 1.0 and math.isnan(adf_pvalue(saturated, "drift"))
    constant = adf_fit(np.full(10, 3.0), 0, "constant")
    assert constant is not None and math.isnan(constant.tau)


def test_the_four_combinations_follow_their_formulas() -> None:
    p = np.array([0.01, 0.2, 0.5, 0.9])
    n = 4
    stat, pval, df = _combine("inverse_chi2", p)
    assert stat == pytest.approx(-2 * np.log(p).sum()) and df == (8.0,)
    assert pval == pytest.approx(stats.chi2.sf(stat, 8))
    stat, pval, df = _combine("inverse_normal", p)
    assert stat == pytest.approx(stats.norm.ppf(p).sum() / 2) and df == ()
    assert pval == pytest.approx(stats.norm.cdf(stat))
    stat, pval, df = _combine("logit", p)
    k = 3 * (5 * n + 3) / (math.pi**2 * n * (5 * n + 2))
    assert stat == pytest.approx(math.sqrt(k) * np.log(p / (1 - p)).sum()) and df == (24.0,)
    assert pval == pytest.approx(stats.t.cdf(stat, 24))
    stat, pval, df = _combine("modified", p)
    assert stat == pytest.approx((-2 * np.log(p) - 2).sum() / 4) and pval == pytest.approx(
        stats.norm.sf(stat)
    )


def test_a_unit_with_p_one_counts_in_chi2_and_modified_but_not_in_z_and_logit() -> None:
    p = np.array([0.01, 0.2, 1.0])
    assert _combine("inverse_chi2", p)[2] == (6.0,)
    assert _combine("logit", p)[2] == (14.0,)  # 5 * 2 + 4: two units enter
    assert _combine("inverse_normal", p)[0] == pytest.approx(
        stats.norm.ppf(p[:2]).sum() / math.sqrt(2)
    )


def declare(frame: pd.DataFrame) -> xt.PanelData:
    return xt.xtset(frame, "id", "t")


def test_the_panel_test_combines_the_per_unit_p_values() -> None:
    result = fisher(declare(PANEL), "y", lags=1)
    assert result.details is not None and len(result.details) == 40
    expected = -2 * np.log(np.maximum(result.details["pvalue"], 2.220446049250313e-16)).sum()
    assert result.statistic == pytest.approx(expected, rel=1e-13) and result.df == (80.0,)
    assert result.details.loc[1, "pvalue"] == pytest.approx(0.232572959312111, rel=1e-12)
    assert result.rejects(0.001) and result.n_obs == 800 and result.n_groups == 40


def test_demean_equals_subtracting_the_period_means_by_hand() -> None:
    frame = PANEL.copy()
    frame["y"] -= frame.groupby("t")["y"].transform("mean")
    assert fisher(declare(PANEL), "y", lags=1, demean=True).statistic == pytest.approx(
        fisher(declare(frame), "y", lags=1).statistic, rel=1e-12
    )


def test_units_without_enough_equations_are_left_out_with_a_note() -> None:
    frame = PANEL[~((PANEL["id"] == 5) & (PANEL["t"] > 2))]  # unit 5 keeps two observations
    result = fisher(declare(frame), "y", lags=1)
    assert result.df == (78.0,) and result.n_groups == 40
    assert "1 unit(s) were left out" in result.summary("en")


def test_refusals() -> None:
    panel = declare(PANEL)
    for bad in (-1, 1.5, True):
        with pytest.raises(SpecificationError, match="lags must be a non-negative integer"):
            fisher(panel, "y", lags=bad)  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="deterministics must be"):
        fisher(panel, "y", lags=1, deterministics="none")  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="method must be"):
        fisher(panel, "y", lags=1, method="pm")  # type: ignore[arg-type]
    with pytest.raises(NumericalError, match="no unit has enough observations"):
        fisher(panel, "y", lags=18)


def test_by_name_in_both_languages() -> None:
    result = run_test("fisher", declare(PANEL), "y", lags=1)
    assert "Fisher-type unit-root test" in result.summary("en")
    assert "type Fisher" in result.summary("fr")
