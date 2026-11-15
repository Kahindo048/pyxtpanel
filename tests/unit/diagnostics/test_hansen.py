"""Hansen J test: values computed by hand, refusals, notes and the shared quadratic form."""

from __future__ import annotations

import math

import numpy as np
import pytest

import pyxtpanel as xt
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import list_tests, run_test
from pyxtpanel.diagnostics._shared.overid import moment_covariance
from pyxtpanel.diagnostics.postestimation.hansen import hansen


def test_hand_computed_one_degree_of_freedom() -> None:
    # m = (1, 2), W = I: J = 5; two instruments, one parameter: chi2(1), p = erfc(sqrt(J / 2)).
    result = hansen([1.0, 2.0], np.eye(2), 1)
    assert result.statistic == 5.0 and result.df == (1.0,)
    assert result.pvalue == pytest.approx(math.erfc(math.sqrt(2.5)), rel=1e-12)
    assert result.statistic_name == "chi2" and result.name == "hansen"


def test_hand_computed_two_degrees_of_freedom_and_a_general_weight() -> None:
    # m'Wm with W = [[2, 1], [1, 3]] and m = (1, -1): 2 - 2 + 3 = 3; three instruments are not
    # needed: pass a third moment of zero with a zero row and column.
    weight = np.array([[2.0, 1, 0], [1, 3, 0], [0, 0, 1]])
    result = hansen([1.0, -1.0, 0.0], weight, 1)
    assert result.statistic == 3.0 and result.df == (2.0,)
    assert result.pvalue == pytest.approx(math.exp(-1.5), rel=1e-12)  # chi2(2) tail


def test_moments_of_each_unit_are_summed() -> None:
    per_unit = np.array([[0.5, 1.0], [0.5, 1.0]])
    assert hansen(per_unit, np.eye(2), 1).statistic == hansen([1.0, 2.0], np.eye(2), 1).statistic
    assert hansen(per_unit, np.eye(2), 1).n_groups == 2


def test_the_weight_is_symmetrised_but_a_clearly_asymmetric_one_is_refused() -> None:
    almost = np.array([[2.0, 1.0], [1.0 + 1e-12, 3.0]])
    assert hansen([1.0, -1.0], almost, 1).statistic == pytest.approx(3.0)
    with pytest.raises(SpecificationError, match="weight must be symmetric"):
        hansen([1.0, -1.0], np.array([[2.0, 1.0], [0.0, 3.0]]), 1)


def test_a_weight_that_is_not_positive_semidefinite_is_a_numerical_error() -> None:
    with pytest.raises(NumericalError, match="negative"):
        hansen([1.0, 0.0], np.diag([-1.0, 1.0]), 1)


def test_not_over_identified_gives_nan_and_a_note_in_both_languages() -> None:
    result = hansen([1.0, 2.0], np.eye(2), 2)
    assert math.isnan(result.statistic) and math.isnan(result.pvalue) and result.df == ()
    assert "not over-identified" in result.summary("en")
    assert "pas sur-identifié" in result.summary("fr")
    assert not result.rejects()


def test_rank_replaces_the_number_of_moments_in_the_degrees_of_freedom() -> None:
    result = hansen([1.0, 2.0, 0.0], np.eye(3), 1, n_instruments=2)
    assert result.df == (1.0,) and result.options["n_instruments"] == 2
    with pytest.raises(SpecificationError, match="n_instruments must be an integer >= 1 and <= 3"):
        hansen([1.0, 2.0, 0.0], np.eye(3), 1, n_instruments=4)


def test_more_instruments_than_groups_adds_a_proliferation_note() -> None:
    result = hansen([1.0, 2.0, 0.5], np.eye(3), 1, n_groups=2)
    assert "exceed the 2 groups" in result.summary("en")
    assert "dépassent les 2 groupes" in result.summary("fr")
    assert not hansen([1.0, 2.0, 0.5], np.eye(3), 1, n_groups=5).notes


def test_bad_shapes_counts_and_values_are_refused() -> None:
    with pytest.raises(SpecificationError, match=r"weight must be 2 x 2, got an array of shape"):
        hansen([1.0, 2.0], np.eye(3), 1)
    with pytest.raises(SpecificationError, match="n_parameters must be an integer >= 0"):
        hansen([1.0, 2.0], np.eye(2), -1)
    with pytest.raises(SpecificationError, match="n_parameters must be an integer"):
        hansen([1.0, 2.0], np.eye(2), 1.5)  # type: ignore[arg-type]
    with pytest.raises(xt.DataError, match="missing or infinite"):
        hansen([1.0, np.nan], np.eye(2), 1)


def test_the_statistic_is_the_quadratic_form_for_random_inputs() -> None:
    rng = np.random.default_rng(3)
    q = 6
    a = rng.normal(size=(q, q))
    weight, m = a @ a.T, rng.normal(size=q)
    result = hansen(m, weight, 2)
    assert result.statistic == pytest.approx(float(np.sum(np.outer(m, m) * weight)), rel=1e-12)
    assert result.df == (4.0,)


def test_moment_covariance_by_hand() -> None:
    g = np.array([[1.0, 2.0], [3.0, 4.0]])
    np.testing.assert_array_equal(moment_covariance(g), [[10.0, 14.0], [14.0, 20.0]])
    np.testing.assert_array_equal(moment_covariance(g, center=True), [[2.0, 2.0], [2.0, 2.0]])


def test_the_test_is_listed_and_run_by_name_and_exported_at_the_root() -> None:
    table = list_tests()
    assert table.loc["hansen", "category"] == "B" and "xtabond2" in table.loc["hansen", "stata"]
    assert run_test("hansen", [1.0, 2.0], np.eye(2), 1).statistic == 5.0
    assert xt.run_test is run_test


def test_the_summary_and_frame_carry_the_numbers() -> None:
    result = hansen([1.0, 2.0], np.eye(2), 1, n_groups=10)
    text = result.summary("en")
    assert "Hansen J test" in text and "chi2 = 5.0000" in text and "degrees of freedom = 1" in text
    assert "H0: les restrictions" in result.summary("fr")
    frame = result.to_frame()
    assert frame.loc["hansen", "statistic"] == 5.0 and frame.loc["hansen", "n_groups"] == 10
