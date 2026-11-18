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


# ---------------------------------------------------------------- Phillips-Perron (L.2.10)


def naive_pp(y: np.ndarray, lags: int, trend: bool) -> tuple[float, float]:
    """Z(t) and Z(rho) of ``pperron`` from the normal equations and explicit autocovariances."""
    t_index = np.arange(1.0, y.size)
    columns = [y[:-1], *([t_index] if trend else []), np.ones(y.size - 1)]
    x = np.column_stack(columns)
    target = y[1:]
    xtx_inverse = np.linalg.inv(x.T @ x)
    beta = xtx_inverse @ x.T @ target
    u = target - x @ beta
    n, k = x.shape
    s2 = float(u @ u) / (n - k)
    sig, s = math.sqrt(s2 * xtx_inverse[0, 0]), math.sqrt(s2)
    gamma0 = float(u @ u) / n
    lam = gamma0
    for i in range(1, lags + 1):
        lam += 2 * (1 - i / (lags + 1)) * float(sum(u[j] * u[j - i] for j in range(i, n))) / n
    z_rho = n * (beta[0] - 1) - (n * sig / s) ** 2 * (lam - gamma0) / 2
    z_t = (beta[0] - 1) / sig * math.sqrt(gamma0 / lam) - 0.5 * (lam - gamma0) * (
        n * sig / s
    ) / math.sqrt(lam)
    return z_t, z_rho


@pytest.mark.parametrize("trend", [False, True])
@pytest.mark.parametrize("lags", [0, 1, 3])
def test_the_phillips_perron_statistics_match_a_direct_computation(trend: bool, lags: int) -> None:
    from pyxtpanel.diagnostics._shared.pperron import pp_fit, pp_pvalue

    y = unit(4)
    case = "trend" if trend else "constant"
    fit = pp_fit(y, lags, case)
    z_t, z_rho = naive_pp(y, lags, trend)
    assert fit is not None and fit.n == y.size - 1
    assert fit.z_t == pytest.approx(z_t, rel=1e-10) and fit.z_rho == pytest.approx(z_rho, rel=1e-10)
    assert 0.0 <= pp_pvalue(fit, case) <= 1.0


def test_without_lags_the_phillips_perron_statistic_is_the_dickey_fuller_one() -> None:
    from pyxtpanel.diagnostics._shared.pperron import pp_fit

    y = unit(6)
    pp, adf = pp_fit(y, 0, "constant"), adf_fit(y, 0, "constant")
    assert pp is not None and adf is not None
    assert pp.z_t == pytest.approx(adf.tau, rel=1e-10)
    assert pp.z_rho == pytest.approx(adf.n * (adf.coefficient), rel=1e-10)  # n (rho - 1), rho - 1


def test_values_printed_by_stata_pperron() -> None:
    from pyxtpanel.diagnostics._shared.pperron import pp_fit, pp_pvalue

    # Stata 17: pperron y if id == 1, lags(2) on sim_ur_mixed (see the validation case fp_lags2).
    fit = pp_fit(unit(1), 2, "constant")
    assert fit is not None and fit.lags == 2
    assert pp_pvalue(fit, "constant") == pytest.approx(
        fisher(declare(PANEL), "y", lags=2, test="pperron").details.loc[1, "pvalue"],  # type: ignore[union-attr]
        rel=1e-13,
    )


def test_degenerate_phillips_perron_regressions_as_pperron_does() -> None:
    from pyxtpanel.diagnostics._shared.pperron import pp_fit, pp_pvalue

    assert pp_fit(np.array([1.0]), 0, "constant") is None  # no equation
    assert pp_fit(np.array([1.0, 2.0]), 0, "constant") is None  # one equation
    saturated = pp_fit(np.array([1.0, 2.0, 4.0]), 0, "constant")  # two equations, two coefficients
    assert (
        saturated is not None
        and math.isnan(saturated.z_t)
        and pp_pvalue(saturated, "constant") == 1.0
    )
    trend = pp_fit(
        np.array([1.0, 2.0, 4.0, 3.0]), 0, "trend"
    )  # three equations, three coefficients
    assert trend is not None and math.isnan(trend.z_t) and pp_pvalue(trend, "trend") == 1.0
    constant = pp_fit(np.full(10, 3.0), 1, "constant")
    assert constant is not None and math.isnan(constant.z_t)
    holed = unit(3).copy()
    holed[7] = np.nan
    assert pp_fit(holed, 1, "constant").n == unit(3).size - 3  # type: ignore[union-attr]


def test_the_fisher_test_with_phillips_perron_p_values() -> None:
    from pyxtpanel.diagnostics._shared.pperron import pp_fit, pp_pvalue

    result = fisher(declare(PANEL), "y", lags=2, test="pperron")
    assert result.options["test"] == "pperron"
    assert result.details is not None and {"statistic", "pvalue", "z_rho"} <= set(result.details)
    expected = [pp_pvalue(pp_fit(unit(i), 2, "constant"), "constant") for i in range(1, 41)]  # type: ignore[arg-type]
    assert result.details["pvalue"].to_numpy() == pytest.approx(expected, rel=1e-13)
    assert result.statistic == pytest.approx(
        -2 * np.log(np.maximum(expected, 2.220446049250313e-16)).sum(), rel=1e-12
    )
    assert result.df == (80.0,) and result.n_groups == 40
    adf = fisher(declare(PANEL), "y", lags=2)
    assert adf.statistic != result.statistic and "z_rho" not in (adf.details.columns)  # type: ignore[union-attr]


def test_phillips_perron_demean_and_trend_and_the_unit_note() -> None:
    frame = PANEL.copy()
    frame["y"] -= frame.groupby("t")["y"].transform("mean")
    assert fisher(
        declare(PANEL), "y", lags=1, demean=True, test="pperron"
    ).statistic == pytest.approx(
        fisher(declare(frame), "y", lags=1, test="pperron").statistic, rel=1e-12
    )
    assert math.isfinite(
        fisher(declare(PANEL), "y", lags=1, deterministics="trend", test="pperron").statistic
    )
    short = PANEL[~((PANEL["id"] == 5) & (PANEL["t"] > 1))]  # unit 5 keeps a single observation
    result = fisher(declare(short), "y", lags=1, test="pperron")
    assert result.df == (78.0,) and "Phillips-Perron regression" in result.summary("en")


def test_phillips_perron_refusals() -> None:
    panel = declare(PANEL)
    with pytest.raises(SpecificationError, match="only available with test='dfuller'"):
        fisher(panel, "y", lags=1, deterministics="drift", test="pperron")
    with pytest.raises(SpecificationError, match="test must be 'dfuller' or 'pperron'"):
        fisher(panel, "y", lags=1, test="kpss")  # type: ignore[arg-type]
