"""Shapiro-Francia test: the squared correlation with the normal scores as oracle and refusals."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import run_test, shapiro_francia


def sample(seed: int = 0, n: int = 300) -> np.ndarray:
    return np.random.default_rng(seed).standard_t(df=5, size=n)


def test_the_statistic_is_the_squared_correlation_with_the_blom_scores() -> None:
    x = sample()
    scores = stats.norm.ppf((stats.rankdata(x) - 0.375) / (x.size + 0.25))
    result = shapiro_francia(x)
    assert result.statistic == pytest.approx(float(np.corrcoef(x, scores)[0, 1] ** 2), rel=1e-12)


def test_the_details_give_w_v_and_z_and_the_p_value_is_the_normal_tail_of_z() -> None:
    result = shapiro_francia(sample())
    assert result.details is not None
    w, v, z = result.details["value"]
    assert w == result.statistic
    assert result.pvalue == pytest.approx(float(stats.norm.sf(z)), rel=1e-13)
    assert v > 1


def test_the_boxcox_transformation_gives_a_different_z_from_the_same_statistic() -> None:
    x = sample()
    log, boxcox = shapiro_francia(x), shapiro_francia(x, transform="boxcox")
    assert log.statistic == pytest.approx(boxcox.statistic, rel=1e-14)
    assert log.details is not None and boxcox.details is not None
    assert log.details.loc["z", "value"] != boxcox.details.loc["z", "value"]
    assert boxcox.options["transform"] == "boxcox"


def test_the_statistic_does_not_depend_on_the_order_nor_on_a_location_scale_change() -> None:
    x = np.random.default_rng(8).normal(size=60)
    base = shapiro_francia(x)
    assert shapiro_francia(np.random.default_rng(9).permutation(x)).statistic == pytest.approx(
        base.statistic, rel=1e-13
    )
    assert shapiro_francia(10 * x + 3).statistic == pytest.approx(base.statistic, rel=1e-12)


def test_rejects_heavy_tails_and_not_normal_data() -> None:
    rng = np.random.default_rng(4)
    assert shapiro_francia(rng.standard_t(df=3, size=800)).rejects(0.001)
    assert not shapiro_francia(rng.normal(size=800)).rejects(0.001)


def test_ties_change_the_scores_and_noties_ranks_them_in_sorted_order() -> None:
    x = np.round(np.random.default_rng(5).normal(size=300), 1)
    assert shapiro_francia(x).statistic != shapiro_francia(x, ties=False).statistic


def test_below_the_smallest_size_of_the_transformation_the_result_is_nan_with_a_note() -> None:
    x = sample(n=9)
    log = shapiro_francia(x)
    assert math.isnan(log.statistic) and math.isnan(log.pvalue)
    assert "at least 10" in log.summary("en") and "au moins 10" in log.summary("fr")
    assert math.isfinite(shapiro_francia(x, transform="boxcox").statistic)  # 5 is enough
    assert math.isnan(shapiro_francia(x[:4], transform="boxcox").statistic)


def test_a_large_sample_gets_a_remark_and_a_value() -> None:
    result = shapiro_francia(np.random.default_rng(6).normal(size=5200))
    assert math.isfinite(result.statistic) and "10 <= n <= 5000" in result.summary("en")


def test_an_unknown_transformation_and_a_constant_sample_are_refused() -> None:
    with pytest.raises(SpecificationError, match="'log' or 'boxcox'"):
        shapiro_francia(sample(), transform="sqrt")  # type: ignore[arg-type]
    with pytest.raises(NumericalError, match="constant"):
        shapiro_francia(np.ones(20))


def test_by_name_and_alias_in_both_languages() -> None:
    result = run_test("sfrancia", sample())
    assert result.name == "shapiro_francia"
    assert "Shapiro-Francia W' test for normality" in result.summary("en")
    assert "Test W' de normalité de Shapiro-Francia" in result.summary("fr")
