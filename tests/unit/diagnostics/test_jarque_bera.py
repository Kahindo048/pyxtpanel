"""Jarque-Bera test: scipy as the independent oracle, and refusals."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import jarque_bera, run_test


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_it_equals_scipy(seed: int) -> None:
    x = np.random.default_rng(seed).standard_t(df=4, size=300)
    ours, theirs = jarque_bera(x), stats.jarque_bera(x)
    assert ours.statistic == pytest.approx(theirs.statistic, rel=1e-12)
    assert ours.pvalue == pytest.approx(theirs.pvalue, rel=1e-10)
    assert ours.df == (2.0,) and ours.n_obs == 300


def test_hand_computed_small_sample() -> None:
    # x = (-2, -1, 0, 1, 2) is symmetric: S = 0; m2 = 2, m4 = 34/5: K = 1.7, so
    # JB = 5 / 6 * (1.7 - 3)^2 / 4.
    result = jarque_bera([-2.0, -1.0, 0.0, 1.0, 2.0])
    assert result.statistic == pytest.approx(5 / 6 * 1.3**2 / 4, rel=1e-13)
    assert result.options["skewness"] == pytest.approx(0.0, abs=1e-15)
    assert result.options["kurtosis"] == pytest.approx(1.7, rel=1e-14)


def test_rejects_heavy_tails_and_not_normal_data() -> None:
    rng = np.random.default_rng(3)
    assert jarque_bera(rng.standard_t(df=3, size=2000)).rejects(0.001)
    assert not jarque_bera(rng.normal(size=2000)).rejects(0.001)


def test_refusals() -> None:
    with pytest.raises(NumericalError, match="constant"):
        jarque_bera(np.ones(10))
    with pytest.raises(SpecificationError, match="n must be at least 2"):
        jarque_bera([1.0])


def test_by_name_in_both_languages() -> None:
    result = run_test("jarque_bera", np.random.default_rng(0).normal(size=50))
    assert "Jarque-Bera test for normality" in result.summary("en")
    assert "Test de normalité de Jarque-Bera" in result.summary("fr")
