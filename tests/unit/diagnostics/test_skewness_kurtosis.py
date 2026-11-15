"""Skewness and kurtosis test: scipy moments as oracle, Royston's adjustment and refusals."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from pyxtpanel._errors import NumericalError
from pyxtpanel.diagnostics import run_test, skewness_kurtosis


def sample(seed: int = 0, n: int = 400) -> np.ndarray:
    return np.random.default_rng(seed).standard_t(df=5, size=n)


def test_the_moments_are_the_biased_ones_of_scipy() -> None:
    x = sample()
    result = skewness_kurtosis(x)
    assert result.details is not None
    assert result.details.loc["skewness", "value"] == pytest.approx(stats.skew(x), rel=1e-12)
    assert result.details.loc["kurtosis", "value"] == pytest.approx(
        stats.kurtosis(x, fisher=False), rel=1e-12
    )


def test_the_unadjusted_statistic_is_the_sum_of_two_squared_normal_scores() -> None:
    x = sample()
    result = skewness_kurtosis(x, adjust=False)
    assert result.details is not None
    z = result.details["z"]
    assert result.statistic == pytest.approx(float(z @ z), rel=1e-13)
    assert result.pvalue == pytest.approx(float(stats.chi2.sf(result.statistic, 2)), rel=1e-12)
    for name in ("skewness", "kurtosis"):
        assert result.details.loc[name, "pvalue"] == pytest.approx(
            2 * float(stats.norm.sf(abs(z[name]))), rel=1e-12
        )


def test_royston_adjustment_changes_the_statistic_and_matches_its_p_value() -> None:
    x = sample(n=60)
    adjusted, plain = skewness_kurtosis(x), skewness_kurtosis(x, adjust=False)
    assert adjusted.statistic != plain.statistic
    assert adjusted.pvalue == pytest.approx(float(stats.chi2.sf(adjusted.statistic, 2)), rel=1e-9)


def test_a_symmetric_sample_has_zero_skewness_and_p_value_one() -> None:
    x = np.concatenate([-np.arange(1.0, 21.0), np.arange(1.0, 21.0)])
    result = skewness_kurtosis(x)
    assert result.details is not None
    assert result.details.loc["skewness", "z"] == pytest.approx(0.0, abs=1e-12)
    assert result.details.loc["skewness", "pvalue"] == pytest.approx(1.0, abs=1e-12)


def test_rejects_heavy_tails_and_not_normal_data() -> None:
    rng = np.random.default_rng(4)
    assert skewness_kurtosis(rng.standard_t(df=3, size=2000)).rejects(0.001)
    assert not skewness_kurtosis(rng.normal(size=2000)).rejects(0.001)


def test_fewer_than_eight_observations_give_nan_and_a_note() -> None:
    result = skewness_kurtosis(np.arange(7.0) ** 2)
    assert math.isnan(result.statistic) and math.isnan(result.pvalue)
    assert "at least 8" in result.summary("en") and "au moins 8" in result.summary("fr")


def test_a_constant_sample_is_refused_and_an_extreme_kurtosis_is_undefined() -> None:
    with pytest.raises(NumericalError, match="constant"):
        skewness_kurtosis(np.ones(20))
    two_points = np.array([0.0] * 10 + [1.0] * 10)  # kurtosis 1: far below the normal value
    assert math.isfinite(skewness_kurtosis(two_points).statistic)


def test_an_undefined_kurtosis_transformation_gives_nan_and_a_note(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Mata's negative ^ (1/3) is missing; no simple sample reaches it, so the branch is forced.
    from pyxtpanel.diagnostics.normality import skewness_kurtosis as module

    monkeypatch.setattr(module, "_z_kurtosis", lambda kurtosis, n: math.nan)
    result = skewness_kurtosis(sample())
    assert math.isnan(result.statistic) and "not defined" in result.summary("en")
    assert "n'est pas définie" in result.summary("fr")


def test_by_name_in_both_languages() -> None:
    result = run_test("skewness_kurtosis", sample())
    assert "Skewness and kurtosis test for normality" in result.summary("en")
    assert "asymétrie et d'aplatissement" in result.summary("fr")
