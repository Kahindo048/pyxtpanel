"""Arellano-Bond AR(m) test: hand-computed values, calendar lags, refusals and notes."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.diagnostics import Equations, run_test
from pyxtpanel.diagnostics.postestimation.ar_test import ar_test


def system(
    residuals: list[float],
    offsets: list[int],
    time: list[int],
    *,
    p: int = 1,
    instruments: bool = True,
) -> Equations:
    n = len(residuals)
    regressors = np.ones((n, p))
    return Equations(
        np.array(residuals),
        regressors,
        np.array(offsets),
        np.array(time),
        np.ones((n, 1)) if instruments else None,
    )


# Two units of three periods; residuals chosen so that the products are easy to add up.
RESIDUALS = [1.0, 2.0, 3.0, 2.0, -1.0, 4.0]
OFFSETS = [0, 3, 6]
TIME = [0, 1, 2, 0, 1, 2]


def test_numerator_and_first_variance_term_by_hand() -> None:
    # Lag 1 pairs: unit 0: 2*1 + 3*2 = 8; unit 1: -1*2 + 4*(-1) = -6.
    # d0 = 2, T1 = 8^2 + 6^2 = 100.  Make the estimation-error terms vanish: zero covariance and
    # instruments orthogonal to the regressors (a zero weight matrix).
    result = ar_test(1, system(RESIDUALS, OFFSETS, TIME), weight=[[0.0]], covariance=[[0.0]])
    assert result.details is not None
    assert result.details.loc["d0", "value"] == 2.0 and result.details.loc["T1", "value"] == 100.0
    assert result.statistic == pytest.approx(2.0 / 10.0, rel=1e-14)
    assert result.pvalue == pytest.approx(math.erfc(0.2 / math.sqrt(2.0)), rel=1e-12)
    assert result.n_obs == 4 and result.n_groups == 2 and result.statistic_name == "m1"


def test_the_second_lag_pairs_periods_two_apart() -> None:
    # Lag 2 pairs: unit 0: 3*1 = 3; unit 1: 4*2 = 8; d0 = 11, T1 = 9 + 64 = 73.
    result = ar_test(2, system(RESIDUALS, OFFSETS, TIME), weight=[[0.0]], covariance=[[0.0]])
    assert result.statistic == pytest.approx(11.0 / math.sqrt(73.0), rel=1e-14)
    assert result.n_obs == 2


def test_the_estimation_error_terms_by_hand() -> None:
    # One regressor and one instrument, all ones, so that M = X'Z = 6, W = 1, V = 0.5:
    # g = sum of the lagged residuals = 1 + 2 + 2 + (-1) = 4  (one lag, found rows only),
    # s = sum_i d_i * sum of e = 8 * 6 + (-6) * 5 = 18, bread = 1 / (6 * 1 * 6) = 1/36,
    # T2 = -2 * 4 * (1/36) * (6 * 1 * 18) = -24, T3 = 4 * 0.5 * 4 = 8; v = 100 - 24 + 8 = 84.
    result = ar_test(1, system(RESIDUALS, OFFSETS, TIME), weight=[[1.0]], covariance=[[0.5]])
    assert result.details is not None
    assert result.details.loc["T2", "value"] == pytest.approx(-24.0, rel=1e-12)
    assert result.details.loc["T3", "value"] == pytest.approx(8.0, rel=1e-12)
    assert result.statistic == pytest.approx(2.0 / math.sqrt(84.0), rel=1e-12)


def test_the_lag_is_read_on_the_calendar_so_a_gap_gives_no_pair() -> None:
    # Unit 0 skips period 1: its rows are at 0 and 2, so lag 1 has no pair there.
    result = ar_test(
        1,
        system([1.0, 3.0, 2.0, -1.0, 4.0], [0, 2, 5], [0, 2, 0, 1, 2]),
        weight=[[0.0]],
        covariance=[[0.0]],
    )
    # Only unit 1 pairs: -1*2 + 4*(-1) = -6; d0 = -6, T1 = 36, one unit contributes.
    assert result.n_obs == 2 and result.n_groups == 1
    assert result.statistic == pytest.approx(-1.0, rel=1e-14)


def test_the_statistic_does_not_depend_on_how_units_are_labelled_or_the_calendar_origin() -> None:
    base = ar_test(1, system(RESIDUALS, OFFSETS, TIME), weight=[[1.0]], covariance=[[0.5]])
    shifted = ar_test(
        1,
        system(RESIDUALS, OFFSETS, [t + 17 for t in TIME]),
        weight=[[1.0]],
        covariance=[[0.5]],
    )
    assert shifted.statistic == base.statistic


def test_a_different_estimation_system_is_used_for_the_instruments() -> None:
    est = Equations(
        np.array([2.0, 1.0, 1.0, 3.0]),  # e.g. forward orthogonal deviations: other rows
        np.ones((4, 1)),
        np.array([0, 2, 4]),
        np.array([0, 1, 0, 1]),
        np.ones((4, 1)),
    )
    both = ar_test(
        1, system(RESIDUALS, OFFSETS, TIME), estimation=est, weight=[[1.0]], covariance=[[0.5]]
    )
    alone = ar_test(1, system(RESIDUALS, OFFSETS, TIME), weight=[[1.0]], covariance=[[0.5]])
    assert both.statistic != alone.statistic and math.isfinite(both.statistic)


def test_undefined_cases_give_nan_and_a_note_in_both_languages() -> None:
    nothing = ar_test(3, system(RESIDUALS, OFFSETS, TIME), weight=[[1.0]], covariance=[[0.5]])
    assert math.isnan(nothing.statistic) and math.isnan(nothing.pvalue) and nothing.n_obs == 0
    assert "no pair of residuals 3 periods apart" in nothing.summary("en")
    assert "aucune paire de résidus distants de 3 périodes" in nothing.summary("fr")
    few = ar_test(
        1, system(RESIDUALS, OFFSETS, TIME), weight=[[1.0]], covariance=[[0.5]], min_pairs=5
    )
    assert math.isnan(few.statistic) and "fewer than the 5 required" in few.summary("en")
    negative = ar_test(1, system(RESIDUALS, OFFSETS, TIME), weight=[[1.0]], covariance=[[-100.0]])
    assert math.isnan(negative.statistic) and "not positive" in negative.summary("en")
    assert not nothing.rejects()


def test_bad_arguments_and_sizes_are_refused() -> None:
    ok = system(RESIDUALS, OFFSETS, TIME)
    with pytest.raises(SpecificationError, match="order must be an integer >= 1"):
        ar_test(0, ok, weight=[[1.0]], covariance=[[1.0]])
    with pytest.raises(SpecificationError, match="min_pairs must be an integer >= 1"):
        ar_test(1, ok, weight=[[1.0]], covariance=[[1.0]], min_pairs=0)
    with pytest.raises(SpecificationError, match="instruments must be given"):
        ar_test(
            1,
            system(RESIDUALS, OFFSETS, TIME, instruments=False),
            weight=[[1.0]],
            covariance=[[1.0]],
        )
    with pytest.raises(SpecificationError, match=r"weight must be 1 x 1"):
        ar_test(1, ok, weight=np.eye(2), covariance=[[1.0]])
    with pytest.raises(SpecificationError, match=r"covariance must be 1 x 1"):
        ar_test(1, ok, weight=[[1.0]], covariance=np.eye(2))
    other = system(RESIDUALS[:3], [0, 3], TIME[:3])
    with pytest.raises(SpecificationError, match="the same units: 2 against 1"):
        ar_test(1, ok, estimation=other, weight=[[1.0]], covariance=[[1.0]])


def test_equations_refuse_inconsistent_layouts() -> None:
    with pytest.raises(SpecificationError, match="offsets must start at 0"):
        system(RESIDUALS, [0, 3, 5], TIME)
    with pytest.raises(DataError, match="time_index must increase strictly within each unit"):
        system(RESIDUALS, OFFSETS, [0, 0, 2, 0, 1, 2])
    with pytest.raises(SpecificationError, match="time_index must have one entry per equation"):
        Equations(np.ones(3), np.ones((3, 1)), np.array([0, 3]), np.array([0, 1]))
    with pytest.raises(DataError, match="missing or infinite"):
        Equations(np.array([1.0, np.nan]), np.ones((2, 1)), np.array([0, 2]), np.array([0, 1]))
    assert system(RESIDUALS, OFFSETS, TIME).n_units == 2
    assert system(RESIDUALS, OFFSETS, TIME).unit_of_row.tolist() == [0, 0, 0, 1, 1, 1]


def test_it_runs_by_name_and_reads_in_both_languages() -> None:
    result = run_test(
        "ar_test", 1, system(RESIDUALS, OFFSETS, TIME), weight=[[1.0]], covariance=[[0.5]]
    )
    assert "Arellano-Bond test for serial correlation" in result.summary("en")
    assert "Test d'Arellano-Bond" in result.summary("fr")
    assert result.to_frame().loc["ar_test", "n_obs"] == 4
