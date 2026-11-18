"""Within, between and first difference against ``linearmodels`` (B.3.5, M.4.2).

``linearmodels`` has its own panel container with ``demean("entity")``, ``mean("entity")`` and
``first_difference()``.  The same frozen scenarios as in ``test_cross_pandas.py`` go through both
libraries, restricted to what ``linearmodels`` defines the same way:

* no missing value in ``x`` (``linearmodels`` drops incomplete rows before it transforms);
* no hole in the calendar of a unit for the first difference (``linearmodels`` subtracts the
  previous *row*, ``pyxtpanel`` the previous *period*, as Stata does: the two agree exactly when
  no period is skipped).

The test is skipped where ``linearmodels`` is not installed; the CI has a job that installs it.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from validation import campaign

import pyxtpanel as xt

pytestmark = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")

linearmodels_data = pytest.importorskip("linearmodels.panel.data")


def usable(scenario: campaign.Scenario) -> xt.PanelData | None:
    """The declared panel of a scenario ``linearmodels`` can be compared on, else ``None``."""
    frame = scenario.frame
    if scenario.fmt or frame["x"].isna().any() or len(frame) < 2:
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return xt.xtset(frame, "id", "t", delta=scenario.delta)
    except xt.DeclarationError:
        return None


SCENARIOS = [(s, p) for s in campaign.load() if (p := usable(s)) is not None]


def by_key(values: np.ndarray, frame: pd.DataFrame) -> pd.Series:
    """Values of the rows of ``frame`` as a Series indexed by (unit, period)."""
    return pd.Series(values, index=pd.MultiIndex.from_frame(frame[["id", "t"]]), dtype=float)


def container(frame: pd.DataFrame) -> object:
    return linearmodels_data.PanelData(frame.set_index(["id", "t"])[["x"]])


def test_the_campaign_gives_enough_panels_to_compare() -> None:
    assert len(SCENARIOS) >= 30


@pytest.mark.parametrize(("scenario", "panel"), SCENARIOS, ids=lambda v: str(getattr(v, "k", "")))
def test_within_and_between_equal_linearmodels(
    scenario: campaign.Scenario, panel: xt.PanelData
) -> None:
    frame = scenario.frame
    lm = container(frame)
    theirs = lm.demean("entity").dataframe["x"]  # type: ignore[attr-defined]
    ours = by_key(panel.within("x").to_numpy(float), frame).reindex(theirs.index)
    np.testing.assert_allclose(ours.to_numpy(), theirs.to_numpy(), rtol=1e-10, atol=1e-12)

    means = lm.mean("entity")["x"]  # type: ignore[attr-defined]
    between = pd.Series(panel.between("x").to_numpy(float), index=frame["id"].to_numpy())
    np.testing.assert_allclose(
        between.groupby(level=0).first().reindex(means.index).to_numpy(),
        means.to_numpy(),
        rtol=1e-10,
        atol=1e-12,
    )


@pytest.mark.parametrize(("scenario", "panel"), SCENARIOS, ids=lambda v: str(getattr(v, "k", "")))
def test_first_difference_equals_linearmodels_where_no_period_is_skipped(
    scenario: campaign.Scenario, panel: xt.PanelData
) -> None:
    if panel.has_gaps:
        pytest.skip("linearmodels subtracts the previous row, not the previous period")
    frame = scenario.frame
    theirs = container(frame).first_difference().dataframe["x"]  # type: ignore[attr-defined]
    ours = by_key(panel.diff("x").to_numpy(float), frame).reindex(theirs.index)
    np.testing.assert_allclose(ours.to_numpy(), theirs.to_numpy(), rtol=1e-10, atol=1e-12)
