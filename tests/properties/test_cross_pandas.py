"""Within, between, first difference and lag, written with pandas alone (B.3.5, M.4.2).

Stata has no command that returns these transformations row by row, and the campaign compares
their summaries (sums, moments).  Here the transformed values themselves, row by row, are compared
with an independent implementation written with ``groupby`` on the 160 frozen scenarios of the
campaign (sizes from one row to 60 units, holes, missing values, constant variables, steps of 1 to
7): the same panels, a different code.

A lag is the value of the *calendar period* before: a hole is not skipped, and a unit with a hole
gets a missing value after it, as Stata does.  pandas' ``shift`` and ``diff`` pair neighbours in
the rows, so the reference code first checks that the neighbour is one step away.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from validation import campaign

import pyxtpanel as xt

#: Scenarios whose time is an integer (a Stata-coded time is converted to dates by the package).
SCENARIOS = [s for s in campaign.load() if not s.fmt]

pytestmark = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")


def declared(scenario: campaign.Scenario) -> xt.PanelData | None:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return xt.xtset(scenario.frame, "id", "t", delta=scenario.delta)
    except xt.DeclarationError:
        return None  # the campaign also holds panels Stata refuses (repeated pairs, missing keys)


def reference_lag(frame: pd.DataFrame, delta: int, name: str, k: int = 1) -> pd.Series:
    """``x`` of the row ``k`` steps earlier in the same unit, missing when no such row exists."""
    ordered = frame.sort_values(["id", "t"])
    grouped = ordered.groupby("id", sort=False)
    previous = grouped[name].shift(k)
    gap = ordered["t"] - grouped["t"].shift(k)
    return previous.where(gap == k * delta).reindex(frame.index)


@pytest.mark.parametrize("scenario", SCENARIOS, ids=lambda s: str(s.k))
def test_within_and_between_equal_the_group_means_of_pandas(scenario: campaign.Scenario) -> None:
    panel = declared(scenario)
    if panel is None:
        pytest.skip("refused by xtset")
    frame = scenario.frame
    mean = frame.groupby("id")["x"].transform("mean")  # skips missing values
    between = mean.where(frame["x"].notna())  # only the rows that have a value of their own
    np.testing.assert_allclose(
        panel.between("x").to_numpy(float),
        between.to_numpy(float),
        rtol=1e-12,
        atol=1e-12,
        equal_nan=True,
    )
    within = frame["x"] - between
    np.testing.assert_allclose(
        panel.within("x").to_numpy(float),
        within.to_numpy(float),
        rtol=1e-12,
        atol=1e-12,
        equal_nan=True,
    )


@pytest.mark.parametrize("scenario", SCENARIOS, ids=lambda s: str(s.k))
def test_lag_and_difference_equal_the_shifted_values_of_pandas(
    scenario: campaign.Scenario,
) -> None:
    panel = declared(scenario)
    if panel is None:
        pytest.skip("refused by xtset")
    frame = scenario.frame
    lag = reference_lag(frame, scenario.delta, "x")
    np.testing.assert_allclose(
        panel.lag("x").to_numpy(float), lag.to_numpy(float), rtol=0, atol=0, equal_nan=True
    )
    np.testing.assert_allclose(
        panel.diff("x").to_numpy(float),
        (frame["x"] - lag).to_numpy(float),
        rtol=1e-12,
        atol=1e-12,
        equal_nan=True,
    )
    lead = reference_lag(frame, scenario.delta, "x", k=-1)
    np.testing.assert_allclose(
        panel.lead("x").to_numpy(float), lead.to_numpy(float), rtol=0, atol=0, equal_nan=True
    )


def test_the_campaign_gives_enough_panels_to_compare() -> None:
    compared = [s for s in SCENARIOS if declared(s) is not None]
    assert len(compared) >= 80
    assert any(declared(s).has_gaps for s in compared)  # type: ignore[union-attr]
