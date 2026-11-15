"""Breusch-Pagan LM test for random effects: hand-computed values and Stata's zero rule."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import lm_random_effects, run_test


def test_hand_computed_two_units() -> None:
    # e = (1, 1 | -1, -1): unit sums 2 and -2, sum of squares 4: A = 1 - 8 / 4 = -1.
    # N = 4, sum T^2 = 8: LM = 16 / 2 * 1 / (8 - 4) = 2, p = chi2(1) tail / 2 = erfc(1) / 2.
    result = lm_random_effects([1.0, 1.0, -1.0, -1.0], [0, 2, 4])
    assert result.statistic == pytest.approx(2.0, rel=1e-14) and result.df == (1.0,)
    assert result.pvalue == pytest.approx(math.erfc(1.0) / 2.0, rel=1e-12)
    assert result.n_obs == 4 and result.n_groups == 2


def test_unequal_unit_sizes_enter_through_the_sum_of_squared_sizes() -> None:
    # Units of 1 and 3 observations: sum T^2 = 10, N = 4.  e = (2 | 1, 1, -3): sums 2 and -1.
    # A = 1 - (4 + 1) / 15 with sum e^2 = 15 ... computed below from the definition.
    e = np.array([2.0, 1.0, 1.0, -3.0])
    a = 1 - (2.0**2 + (-1.0) ** 2) / float(e @ e)
    result = lm_random_effects(e, [0, 1, 4])
    assert result.statistic == pytest.approx(16 / 2 * a**2 / (10 - 4), rel=1e-13)


def test_a_zero_random_effects_variance_gives_zero_and_one_like_xttest0() -> None:
    result = lm_random_effects([1.0, 1.0, -1.0, -1.0], [0, 2, 4], variance_u=0.0)
    assert result.statistic == 0.0 and result.pvalue == 1.0
    assert "is zero" in result.summary("en") and "est nulle" in result.summary("fr")
    assert lm_random_effects([1.0, 1.0, -1.0, -1.0], [0, 2, 4], variance_u=0.3).statistic == 2.0


def test_it_detects_a_unit_effect() -> None:
    rng = np.random.default_rng(1)
    effect = np.repeat(rng.normal(size=50), 6)
    offsets = np.arange(0, 301, 6)
    assert lm_random_effects(effect + rng.normal(size=300), offsets).rejects(0.001)
    assert not lm_random_effects(rng.normal(size=300), offsets).rejects(0.001)


def test_refusals() -> None:
    with pytest.raises(SpecificationError, match="offsets must start at 0"):
        lm_random_effects([1.0, 2.0], [0, 1])
    with pytest.raises(NumericalError, match="single observation"):
        lm_random_effects([1.0, 2.0, 3.0], [0, 1, 2, 3])
    with pytest.raises(NumericalError, match="all residuals are zero"):
        lm_random_effects(np.zeros(4), [0, 2, 4])


def test_by_name_in_both_languages() -> None:
    result = run_test("lm_random_effects", [1.0, 1.0, -1.0, -1.0], [0, 2, 4])
    assert "Breusch and Pagan Lagrangian multiplier test" in result.summary("en")
    assert "multiplicateur de Lagrange" in result.summary("fr")
