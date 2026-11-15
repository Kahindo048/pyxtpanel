"""Shapiro-Wilk test: scipy's W as oracle, Royston's branches, ties, edge cases and refusals."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from pyxtpanel._errors import NumericalError
from pyxtpanel.diagnostics import run_test, shapiro_wilk


@pytest.mark.parametrize("n", [3, 4, 5, 6, 7, 11, 12, 40, 400])
def test_w_agrees_with_scipy_to_the_precision_of_royston_s_coefficients(n: int) -> None:
    x = np.random.default_rng(n).standard_t(df=5, size=n)
    result = shapiro_wilk(x)
    w, p = stats.shapiro(x)
    # scipy implements AS R94 (Royston 1995), Stata the older approximations: close, not equal.
    assert result.statistic == pytest.approx(w, abs=5e-4)
    assert result.pvalue == pytest.approx(p, abs=2e-2) or result.pvalue < 1e-3


def test_the_details_give_w_v_and_z_and_the_p_value_is_the_normal_tail_of_z() -> None:
    x = np.random.default_rng(3).standard_t(df=4, size=200)
    result = shapiro_wilk(x)
    assert result.details is not None
    w, v, z = result.details["value"]
    assert w == result.statistic
    assert result.pvalue == pytest.approx(float(stats.norm.sf(z)), rel=1e-13)
    assert v > 1  # heavy tails: the index of departure from normality exceeds 1


def test_the_statistic_does_not_depend_on_the_order_nor_on_a_location_scale_change() -> None:
    x = np.random.default_rng(8).normal(size=60)
    base = shapiro_wilk(x)
    shuffled = shapiro_wilk(np.random.default_rng(9).permutation(x))
    scaled = shapiro_wilk(10 * x + 3)
    assert shuffled.statistic == pytest.approx(base.statistic, rel=1e-13)
    assert scaled.statistic == pytest.approx(base.statistic, rel=1e-12)
    assert scaled.pvalue == pytest.approx(base.pvalue, rel=1e-9)


def test_rejects_heavy_tails_and_not_normal_data() -> None:
    rng = np.random.default_rng(4)
    assert shapiro_wilk(rng.standard_t(df=3, size=800)).rejects(0.001)
    assert shapiro_wilk(rng.exponential(size=100)).rejects(0.001)
    assert not shapiro_wilk(rng.normal(size=800)).rejects(0.001)


def test_ties_change_the_scores_and_noties_ranks_them_in_sorted_order() -> None:
    x = np.round(np.random.default_rng(5).normal(size=300), 1)  # many tied values
    with_ties, without = shapiro_wilk(x), shapiro_wilk(x, ties=False)
    assert with_ties.statistic != without.statistic
    assert with_ties.options == {"ties": True} and without.options == {"ties": False}
    unique = np.random.default_rng(5).normal(size=300)  # no ties: both rank the same way
    assert shapiro_wilk(unique).statistic == shapiro_wilk(unique, ties=False).statistic


def test_three_values_use_the_exact_p_value_and_a_remark_about_the_range() -> None:
    result = shapiro_wilk(np.array([0.946041, -0.692712, 0.154832]))
    assert result.statistic == pytest.approx(0.9996062348, rel=1e-8)  # swilk, Stata 17
    assert result.pvalue == pytest.approx(0.9621001270, rel=1e-8)
    assert "4 <= n <= 2000" in result.summary("en")


def test_n3_at_the_lower_bound_has_a_zero_p_value_as_stata_shows_it_negative() -> None:
    # Stata 17 prints W = .75, V = 3.7320633, z missing and Prob>z = -0.00005 for (0, 0, 1)
    result = shapiro_wilk(np.array([0.0, 0.0, 1.0]))
    assert result.details is not None
    assert result.statistic == pytest.approx(0.75, rel=1e-13)
    assert result.details.loc["V", "value"] == pytest.approx(3.7320633, rel=1e-7)
    assert math.isnan(result.details.loc["z", "value"])
    assert result.pvalue == 0.0
    assert "lower bound" in result.summary("en") and "borne inférieure" in result.summary("fr")


def test_a_large_sample_gets_a_remark_and_a_value() -> None:
    result = shapiro_wilk(np.random.default_rng(6).normal(size=2500))
    assert math.isfinite(result.statistic) and "valid for 4 <= n <= 2000" in result.summary("en")


def test_fewer_than_three_observations_give_nan_and_a_note() -> None:
    result = shapiro_wilk(np.array([1.0, 2.0]))
    assert math.isnan(result.statistic) and math.isnan(result.pvalue)
    assert "at least 3" in result.summary("en") and "au moins 3" in result.summary("fr")


def test_a_constant_sample_is_refused() -> None:
    with pytest.raises(NumericalError, match="constant"):
        shapiro_wilk(np.ones(20))


def test_by_name_and_alias_in_both_languages() -> None:
    x = np.random.default_rng(7).normal(size=50)
    result = run_test("swilk", x)
    assert result.name == "shapiro_wilk"
    assert "Shapiro-Wilk W test for normality" in result.summary("en")
    assert "Test W de normalité de Shapiro-Wilk" in result.summary("fr")
