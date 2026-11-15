"""Difference-in-Hansen: values computed by hand and refusals."""

from __future__ import annotations

import math

import pytest

from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics import diff_hansen, hansen, run_test


def test_hand_computed_one_degree_of_freedom() -> None:
    # 10 - 7.5 = 2.5 on 1 df: p = erfc(sqrt(2.5 / 2)).
    result = diff_hansen(10.0, 7.5, 1)
    assert result.statistic == 2.5 and result.df == (1.0,) and result.name == "diff_hansen"
    assert result.pvalue == pytest.approx(math.erfc(math.sqrt(1.25)), rel=1e-12)


def test_hand_computed_two_degrees_of_freedom() -> None:
    result = diff_hansen(6.0, 2.0, 2)  # chi2(2) tail: exp(-C / 2)
    assert result.statistic == 4.0 and result.pvalue == pytest.approx(math.exp(-2.0), rel=1e-12)


def test_equal_statistics_give_a_zero_difference_and_a_p_value_of_one() -> None:
    result = diff_hansen(5.0, 5.0, 3)
    assert result.statistic == 0.0 and result.pvalue == 1.0 and not result.rejects()


def test_a_negative_difference_is_reported_without_p_value() -> None:
    result = diff_hansen(3.0, 5.0, 2)
    assert result.statistic == -2.0 and math.isnan(result.pvalue)
    assert "difference is negative" in result.summary("en")
    assert "différence est négative" in result.summary("fr")


def test_it_composes_with_two_hansen_results() -> None:
    from numpy import eye

    full = hansen([1.0, 2.0, 1.5], eye(3), 1)  # J = 7.25 on 2 df
    restricted = hansen([1.0, 2.0], eye(2), 1)  # J = 5 on 1 df
    result = diff_hansen(full.statistic, restricted.statistic, 1, df_full=int(full.df[0]))
    assert result.statistic == pytest.approx(2.25) and result.df == (1.0,)


def test_bad_arguments_are_refused() -> None:
    for j in (-1.0, math.nan, math.inf):
        with pytest.raises(SpecificationError, match="finite non-negative statistic"):
            diff_hansen(j, 1.0, 1)
        with pytest.raises(SpecificationError, match="j_excluding must be a finite"):
            diff_hansen(1.0, j, 1)
    for n in (0, -1, 1.5):
        with pytest.raises(SpecificationError, match="n_dropped must be a positive integer"):
            diff_hansen(2.0, 1.0, n)  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="at most the 3 degrees of freedom"):
        diff_hansen(2.0, 1.0, 4, df_full=3)


def test_it_runs_by_name_and_reads_in_both_languages() -> None:
    result = run_test("diff_hansen", 10.0, 7.5, 1)
    assert "Difference-in-Hansen" in result.summary("en")
    assert "Hansen en différence" in result.summary("fr")
