"""Modified Wald test for groupwise heteroskedasticity: hand-computed values and refusals."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics import run_test, wald_groupwise

# Unit A: e^2 = (1, 1, 4, 4); unit B: e^2 = (9, 9, 0, 0).
E = np.array([1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.0, 0.0])
OFFSETS = [0, 4, 8]


def test_hand_computed_two_units() -> None:
    # s_A = 2.5, V_A = mean((e^2 - s_A)^2) / (T - 1) = 2.25 / 3 = 0.75; s_B = 4.5,
    # V_B = 20.25 / 3 = 6.75; overall variance = 28 / 8 = 3.5.
    # W = (2.5 - 3.5)^2 / 0.75 + (4.5 - 3.5)^2 / 6.75 on 2 df, p = exp(-W / 2).
    result = wald_groupwise(E, OFFSETS)
    expected = 1 / 0.75 + 1 / 6.75
    assert result.statistic == pytest.approx(expected, rel=1e-14) and result.df == (2.0,)
    assert result.pvalue == pytest.approx(math.exp(-expected / 2), rel=1e-12)
    assert result.n_obs == 8 and result.n_groups == 2 and not result.notes


def test_the_details_give_the_size_and_variance_of_each_unit() -> None:
    result = wald_groupwise(E, OFFSETS)
    assert result.details is not None and list(result.details["T"]) == [4, 4]
    assert result.details["variance"].tolist() == pytest.approx([2.5, 4.5])


def test_a_unit_without_variance_of_its_squares_is_skipped_but_counts_in_the_df() -> None:
    e = np.array([1.0, -1.0, 1.0, -1.0, 3.0, -3.0, 0.5, -0.5])  # unit A: e^2 is constant
    result = wald_groupwise(e, OFFSETS)
    assert result.df == (2.0,) and "1 units have no usable variance" in result.summary("en")
    single = wald_groupwise(np.array([1.0, 2.0, -3.0, 1.0, 0.5]), [0, 1, 5])  # T = 1 then T = 4
    assert single.df == (2.0,) and single.notes


def test_invariant_to_the_scale_of_the_residuals() -> None:
    assert wald_groupwise(10 * E, OFFSETS).statistic == pytest.approx(
        wald_groupwise(E, OFFSETS).statistic, rel=1e-12
    )


def within(rng: np.random.Generator, t: int) -> np.ndarray:
    e = rng.normal(size=(40, t))
    return (e - e.mean(axis=1, keepdims=True)).ravel()  # residuals of a fixed-effects fit


def test_a_common_variance_does_not_reject_and_different_ones_do_with_long_panels() -> None:
    rng = np.random.default_rng(0)
    same = within(rng, 400)
    sd = np.repeat(np.exp(rng.normal(scale=0.7, size=40)), 400)
    offsets = np.arange(0, 40 * 400 + 1, 400)
    assert not wald_groupwise(same, offsets).rejects(0.001)
    assert wald_groupwise(same * sd, offsets).rejects(0.001)


def test_the_test_over_rejects_in_short_panels_as_its_author_warns() -> None:
    # 200 simulated panels of 40 units with a common variance: at T = 8 the 5% test rejects
    # every time, at T = 200 about one time in ten (Baum's warning about the size).
    rng = np.random.default_rng(1)
    critical = 55.758  # 95th percentile of chi2(40)
    rates = {}
    for t in (8, 200):
        offsets = np.arange(0, 40 * t + 1, t)
        stats_ = [wald_groupwise(within(rng, t), offsets).statistic for _ in range(60)]
        rates[t] = float(np.mean(np.array(stats_) > critical))
    assert rates[8] > 0.9 and 0.0 < rates[200] < 0.4


def test_refusals() -> None:
    for bad in ([0, 4, 7], [1, 4, 8], [0, 8, 8], [0, 4, 4, 8], [[0, 8]]):
        with pytest.raises(SpecificationError, match="offsets must start at 0"):
            wald_groupwise(E, bad)


def test_by_name_in_both_languages() -> None:
    result = run_test("wald_groupwise", E, OFFSETS)
    assert "Modified Wald test for groupwise heteroskedasticity" in result.summary("en")
    assert "hétéroscédasticité par groupe" in result.summary("fr")
