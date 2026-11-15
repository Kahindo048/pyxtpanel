"""F test of the unit effects: hand-computed values and refusals."""

from __future__ import annotations

import pytest

from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics import f_effects, run_test


def test_hand_computed() -> None:
    # (20 - 10) / 2 over 10 / 10 = 5 on F(2, 10); for two numerator df the tail is
    # (1 + 2 F / d2)^(-d2 / 2) = 2^-5.
    result = f_effects(20.0, 10.0, n_obs=14, n_groups=3, n_regressors=1)
    assert result.statistic == pytest.approx(5.0, rel=1e-14) and result.df == (2.0, 10.0)
    assert result.pvalue == pytest.approx(1 / 32, rel=1e-12)
    assert result.statistic_name == "F" and result.n_groups == 3 and result.n_obs == 14


def test_equal_sums_of_squares_give_zero_and_p_one() -> None:
    result = f_effects(10.0, 10.0, 14, 3, 1)
    assert result.statistic == 0.0 and result.pvalue == 1.0


def test_refusals() -> None:
    with pytest.raises(SpecificationError, match="rss_pooled must be at least rss_effects"):
        f_effects(5.0, 10.0, 14, 3, 1)
    for bad in (0.0, -1.0, float("nan")):
        with pytest.raises(SpecificationError, match="rss_effects must be a positive"):
            f_effects(10.0, bad, 14, 3, 1)
    with pytest.raises(SpecificationError, match="n_groups must be at least 2"):
        f_effects(20.0, 10.0, 14, 1, 1)
    with pytest.raises(SpecificationError, match="n_obs must be larger than"):
        f_effects(20.0, 10.0, 4, 3, 1)


def test_by_name_in_both_languages() -> None:
    result = run_test("f_effects", 20.0, 10.0, 14, 3, 1)
    assert "F test that all unit effects are zero" in result.summary("en")
    assert "nullité de tous les effets" in result.summary("fr")
