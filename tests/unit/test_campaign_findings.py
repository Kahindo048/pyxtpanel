"""What the random campaign found (J.4.2): the defect that was fixed and the divergence kept.

Each test fails if the behaviour it documents changes; ``validation/campaign.py`` carries the Stata
side (``tests/golden/campaign_golden.json``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt


def panel(rows: list[tuple[int, int, float, float]]) -> xt.PanelData:
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "x", "z"]), "id", "t")


def test_a_constant_variable_has_no_correlation_not_a_correlation_of_rounding_noise() -> None:
    """Defect: x - mean(x) of a constant left 1e-16, and xtcorr reported -2e-18, not missing."""
    rng = np.random.default_rng(1)
    rows = [(u, t, 0.123457, float(rng.normal())) for u in range(1, 6) for t in range(1, 8)]
    result = panel(rows).xtcorr(["x", "z"])
    for kind in ("overall", "within", "between"):
        assert np.isnan(result.corr[kind].loc["x", "z"]), kind


def test_a_variable_that_is_constant_in_each_unit_has_no_within_correlation() -> None:
    rng = np.random.default_rng(2)
    rows = [(u, t, 0.1 * u, float(rng.normal())) for u in range(1, 6) for t in range(1, 8)]
    result = panel(rows).xtcorr(["x", "z"])
    assert np.isnan(result.corr["within"].loc["x", "z"])
    assert np.isfinite(result.corr["overall"].loc["x", "z"])
    assert np.isfinite(result.corr["between"].loc["x", "z"])


def test_a_small_but_real_variation_is_not_taken_for_noise() -> None:
    rng = np.random.default_rng(3)
    rows = [
        (u, t, 1e6 + 1e-2 * float(rng.normal()), float(rng.normal()))
        for u in range(1, 6)
        for t in range(1, 30)
    ]
    assert np.isfinite(panel(rows).xtcorr(["x", "z"]).corr["overall"].loc["x", "z"])


def test_xtdescribe_refuses_a_panel_observed_in_a_single_period_as_stata_does() -> None:
    """Stata: r(459).  The campaign found that this refusal was already in the package."""
    one_period = panel([(1, 5, 1.0, 2.0), (2, 5, 3.0, 4.0), (3, 5, 5.0, 6.0)])
    with pytest.raises(xt.DataError) as caught:
        one_period.xtdescribe()
    assert caught.value.stata_code == 459


def test_divergence_xtcorr_reports_the_points_of_a_pair_with_fewer_than_two() -> None:
    """Stata's correlate stores nothing below 2 observations; the package gives n, r missing."""
    result = panel([(1, 1, 1.0, 2.0), (1, 2, np.nan, 3.0), (2, 1, np.nan, 4.0)]).xtcorr(["x", "z"])
    assert result.n["overall"].loc["x", "z"] == 1
    assert np.isnan(result.corr["overall"].loc["x", "z"])
