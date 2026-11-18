"""Invariances of the core on the frozen campaign scenarios (J.4.1), without Stata.

The 160 scenarios of ``validation/campaign.py`` (sizes from one row to 60 units, holes, entry and
exit, missing and constant variables, steps of 1 to 7, Stata quarters and months) go through the
battery of the package after transformations that must change nothing:

* shuffling the rows of the frame;
* relabelling the identifiers in the same order (text, large or negative integers);
* relabelling them in any order (the statistics that weigh a position then change, the others not);
* shifting the origin of the time variable;
* multiplying the time by a factor and the step by the same factor.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from validation import campaign

SCENARIOS = campaign.load()
INTEGER_TIME = [s for s in SCENARIOS if not s.fmt]
POSITION = [i for i, name in enumerate(campaign.COLUMNS) if name.endswith("wsum")]
VALUE_OF_TIME = [campaign.COLUMNS.index(c) for c in ("tmin", "tmax")]
#: Columns unaffected by the order of the units or by the values of the time variable.
ORDER_FREE = [i for i in range(len(campaign.COLUMNS)) if i not in POSITION + VALUE_OF_TIME]
LABEL_FREE = [i for i in range(len(campaign.COLUMNS)) if i not in VALUE_OF_TIME]

pytestmark = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")


def row(
    frame: pd.DataFrame, scenario: campaign.Scenario, *, delta: int | None = None
) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return campaign.python_row(frame, scenario.delta if delta is None else delta, scenario.fmt)


BASE = {s.k: row(s.frame, s) for s in SCENARIOS}


def same(a: np.ndarray, b: np.ndarray, columns: list[int]) -> None:
    np.testing.assert_allclose(a[columns], b[columns], rtol=1e-12, atol=1e-12, equal_nan=True)


@pytest.mark.parametrize("scenario", SCENARIOS, ids=lambda s: f"{s.k}")
def test_shuffling_the_rows_changes_nothing(scenario: campaign.Scenario) -> None:
    frame = scenario.frame.sample(frac=1, random_state=scenario.k).reset_index(drop=True)
    same(BASE[scenario.k], row(frame, scenario), list(range(len(campaign.COLUMNS))))


@pytest.mark.parametrize("scenario", SCENARIOS, ids=lambda s: f"{s.k}")
def test_relabelling_the_identifiers_in_the_same_order_changes_nothing(
    scenario: campaign.Scenario,
) -> None:
    ids = scenario.frame["id"]
    for labels in (
        ids.map(lambda u: f"firm_{u + 1000:07d}"),  # text that sorts like the integers
        ids * 1000 + 7,  # large integers
        ids - 100_000,  # negative integers
    ):
        same(
            BASE[scenario.k],
            row(scenario.frame.assign(id=labels), scenario),
            list(range(len(campaign.COLUMNS))),
        )


@pytest.mark.parametrize("scenario", SCENARIOS, ids=lambda s: f"{s.k}")
def test_relabelling_the_identifiers_in_any_order_keeps_what_does_not_weigh_a_position(
    scenario: campaign.Scenario,
) -> None:
    units = scenario.frame["id"].unique()
    shuffled = np.random.default_rng(scenario.k).permutation(len(units))
    names = {u: f"unit_{shuffled[i]:05d}" for i, u in enumerate(units)}
    same(
        BASE[scenario.k],
        row(scenario.frame.assign(id=scenario.frame["id"].map(names)), scenario),
        ORDER_FREE,
    )


@pytest.mark.parametrize("scenario", SCENARIOS, ids=lambda s: f"{s.k}")
def test_shifting_the_origin_of_time_moves_only_its_own_values(scenario: campaign.Scenario) -> None:
    shift = 37 * scenario.delta
    shifted = row(scenario.frame.assign(t=scenario.frame["t"] + shift), scenario)
    base = BASE[scenario.k]
    same(base, shifted, LABEL_FREE)
    if not np.isnan(base[VALUE_OF_TIME]).any():
        np.testing.assert_allclose(shifted[VALUE_OF_TIME], base[VALUE_OF_TIME] + shift)


@pytest.mark.parametrize("scenario", INTEGER_TIME, ids=lambda s: f"{s.k}")
def test_scaling_the_time_and_its_step_together_changes_nothing(
    scenario: campaign.Scenario,
) -> None:
    factor = 3
    scaled = row(
        scenario.frame.assign(t=scenario.frame["t"] * factor),
        scenario,
        delta=scenario.delta * factor,
    )
    columns = [i for i in LABEL_FREE if campaign.COLUMNS[i] != "tdelta"]
    same(BASE[scenario.k], scaled, columns)
