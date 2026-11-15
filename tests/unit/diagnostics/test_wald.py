"""Joint Wald test: values computed by hand and the chi-square / F versions."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import run_test
from pyxtpanel.diagnostics.postestimation.wald import wald

B = np.array([1.0, 2.0, 3.0])
V = np.diag([1.0, 4.0, 9.0])
FIRST_TWO = np.array([[1.0, 0, 0], [0, 1, 0]])


def test_chi_square_by_hand() -> None:
    # Rb = (1, 2), R V R' = diag(1, 4): chi2 = 1 + 1 = 2 on 2 df, p = exp(-1).
    result = wald(B, V, FIRST_TWO)
    assert result.statistic == pytest.approx(2.0, rel=1e-14) and result.statistic_name == "chi2"
    assert result.df == (2.0,) and result.pvalue == pytest.approx(math.exp(-1.0), rel=1e-12)


def test_the_identity_tests_every_coefficient_against_zero() -> None:
    result = wald(B, V)
    assert result.statistic == pytest.approx(3.0) and result.df == (3.0,)


def test_f_version_divides_by_the_number_of_restrictions() -> None:
    # F = 2 / 2 = 1 on (2, 10); for two numerator df the tail is (1 + 2 F / d)^(-d / 2).
    result = wald(B, V, FIRST_TWO, df_resid=10)
    assert result.statistic == pytest.approx(1.0) and result.statistic_name == "F"
    assert result.df == (2.0, 10.0)
    assert result.pvalue == pytest.approx((1 + 2 * 1.0 / 10) ** (-5), rel=1e-12)


def test_non_zero_values_and_a_single_restriction() -> None:
    zero = wald(B, V, FIRST_TWO, [1.0, 2.0])
    assert zero.statistic == 0.0 and zero.pvalue == 1.0
    one = wald(B, V, [[0.0, 0.0, 1.0]], [0.0])  # b3 = 0: 9 / 9
    assert one.statistic == pytest.approx(1.0) and one.df == (1.0,)
    assert one.pvalue == pytest.approx(math.erfc(math.sqrt(0.5)), rel=1e-12)


def test_a_linear_combination_uses_the_covariance_terms() -> None:
    v = np.array([[1.0, 0.5, 0.0], [0.5, 4.0, 0.0], [0.0, 0.0, 9.0]])
    # b1 - b2 = -1; var = 1 + 4 - 2 * 0.5 = 4: chi2 = 1 / 4.
    result = wald(B, v, [[1.0, -1.0, 0.0]])
    assert result.statistic == pytest.approx(0.25)


def test_redundant_restrictions_are_refused() -> None:
    twice = np.array([[1.0, 0, 0], [2.0, 0, 0]])
    with pytest.raises(NumericalError, match="not independent"):
        wald(B, V, twice)


def test_bad_inputs_are_refused() -> None:
    with pytest.raises(SpecificationError, match="restrictions must be n x 3"):
        wald(B, V, np.ones((2, 2)))
    with pytest.raises(SpecificationError, match="values must be 2"):
        wald(B, V, FIRST_TWO, [1.0])
    for bad in (0.0, -3.0, math.nan):
        with pytest.raises(SpecificationError, match="df_resid must be a positive number"):
            wald(B, V, FIRST_TWO, df_resid=bad)


def test_by_name_and_in_both_languages() -> None:
    result = run_test("wald", B, V, FIRST_TWO)
    assert "Joint Wald test" in result.summary("en") and "Wald conjoint" in result.summary("fr")
    assert result.to_frame().loc["wald", "df"] == 2.0
