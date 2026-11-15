"""Hausman test: values computed by hand and Stata's conventions."""

from __future__ import annotations

import math

import numpy as np
import pytest

import pyxtpanel as xt
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics import run_test
from pyxtpanel.diagnostics.postestimation.hausman import hausman


def test_hand_computed_with_diagonal_covariances() -> None:
    # b1 - b2 = (1, 2); V1 - V2 = diag(2, 4): H = 1/2 + 4/4 = 1.5 on 2 df, p = exp(-H / 2).
    result = hausman([1.0, 2.0], [0.0, 0.0], np.diag([3.0, 5.0]), np.diag([1.0, 1.0]))
    assert result.statistic == pytest.approx(1.5, rel=1e-14) and result.df == (2.0,)
    assert result.pvalue == pytest.approx(math.exp(-0.75), rel=1e-12)
    assert result.name == "hausman" and not result.notes


def test_a_general_covariance_difference() -> None:
    d = np.array([1.0, -1.0])
    vd = np.array([[2.0, 1.0], [1.0, 3.0]])  # inverse = [[3, -1], [-1, 2]] / 5
    expected = float(d @ np.linalg.inv(vd) @ d)
    result = hausman(d, np.zeros(2), vd + np.eye(2), np.eye(2))
    assert result.statistic == pytest.approx(expected, rel=1e-13) == pytest.approx(7 / 5)


def test_a_singular_difference_uses_the_rank_as_degrees_of_freedom_like_stata() -> None:
    result = hausman([1.0, 1.0], [0.0, 0.0], np.diag([3.0, 1.0]), np.diag([1.0, 1.0]))
    # V1 - V2 = diag(2, 0): generalised inverse diag(1/2, 0), rank 1.
    assert result.statistic == pytest.approx(0.5) and result.df == (1.0,)
    assert "rank of the differenced covariance matrix (1)" in result.summary("en")
    assert "rang de la différence" in result.summary("fr")


def test_a_negative_statistic_is_returned_with_no_p_value_and_a_note() -> None:
    # Efficient estimator with the larger variance: V1 - V2 = diag(-2, -4).
    result = hausman([1.0, 2.0], [0.0, 0.0], np.diag([1.0, 1.0]), np.diag([3.0, 5.0]))
    assert result.statistic == pytest.approx(-1.5) and math.isnan(result.pvalue)
    text = result.summary("en")
    assert "not positive definite" in text and "statistic is negative" in text
    assert not result.rejects()


def test_keep_leaves_out_coefficients_such_as_the_constant() -> None:
    b1, b2 = [1.0, 2.0, 99.0], [0.0, 0.0, 0.0]
    v1, v2 = np.diag([3.0, 5.0, 7.0]), np.eye(3)
    kept = hausman(b1, b2, v1, v2, keep=[0, 1])
    assert kept.statistic == pytest.approx(1.5) and kept.df == (2.0,)
    assert hausman(b1[:2], b2[:2], v1[:2, :2], v2[:2, :2]).statistic == kept.statistic


def test_df_option_replaces_the_rank() -> None:
    result = hausman([1.0, 2.0], [0.0, 0.0], np.diag([3.0, 5.0]), np.diag([1.0, 1.0]), df=1)
    assert result.df == (1.0,) and result.pvalue == pytest.approx(math.erfc(math.sqrt(0.75)))


def test_the_asymmetric_difference_is_symmetrised() -> None:
    v1 = np.array([[3.0, 0.5], [0.5 + 1e-9, 4.0]])
    result = hausman([1.0, 2.0], [0.0, 0.0], v1, np.eye(2))
    sym = hausman([1.0, 2.0], [0.0, 0.0], (v1 + v1.T) / 2, np.eye(2))
    assert result.statistic == pytest.approx(sym.statistic, rel=1e-12)


def test_names_give_a_table_of_the_differences() -> None:
    result = hausman([1.0, 2.0], [0.0, 0.5], np.diag([3.0, 5.0]), np.eye(2), names=["x", "z"])
    table = result.details
    assert table is not None and list(table.index) == ["x", "z"]
    assert table.loc["z", "difference"] == 1.5 and table.loc["x", "std_err"] == pytest.approx(
        math.sqrt(2.0)
    )
    assert "consistent" in result.summary("en")


def test_bad_sizes_are_refused() -> None:
    with pytest.raises(SpecificationError, match="b_efficient must be 2, got an array of shape"):
        hausman([1.0, 2.0], [0.0], np.eye(2), np.eye(2))
    with pytest.raises(SpecificationError, match=r"v_consistent must be 2 x 2"):
        hausman([1.0, 2.0], [0.0, 0.0], np.eye(3), np.eye(2))
    with pytest.raises(SpecificationError, match="names must have 1 entries, got 2"):
        hausman([1.0, 2.0], [0.0, 0.0], np.eye(2), np.eye(2), keep=[0], names=["a", "b"])
    with pytest.raises(xt.DataError, match="missing or infinite"):
        hausman([1.0, np.nan], [0.0, 0.0], np.eye(2), np.eye(2))


def test_it_runs_by_name_and_summarises_in_both_languages() -> None:
    result = run_test("hausman", [1.0, 2.0], [0.0, 0.0], np.diag([3.0, 5.0]), np.eye(2))
    assert "Hausman specification test" in result.summary("en")
    assert "Test de spécification de Hausman" in result.summary("fr")
    assert "systematic" in result.alternative("en") and "systématique" in result.alternative("fr")
