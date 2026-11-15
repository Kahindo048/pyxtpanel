"""Cost of balancing (K.2.3): by hand, window enumeration oracle, agreement with ``balance``."""

from __future__ import annotations

import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pyxtpanel as xt

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")

#: Presence over the calendar 1..4: A 1111, B 0011, C 1100, D 1101, E 0101.
HAND = {"A": [1, 2, 3, 4], "B": [3, 4], "C": [1, 2], "D": [1, 2, 4], "E": [2, 4]}


def declare(units: dict[object, list[int]]) -> xt.PanelData:
    rows = [(u, t) for u, ts in units.items() for t in ts]
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t"]), "id", "t")


def test_balance_cost_by_hand() -> None:
    cost = declare(HAND).balance_cost()
    assert cost.thresholds["units"].to_dict() == {1: 5, 2: 5, 3: 2, 4: 1}
    assert cost.thresholds["obs"].to_dict() == {1: 13, 2: 13, 3: 7, 4: 4}
    assert cost.windows["units"].to_dict() == {1: 4, 2: 3, 3: 1, 4: 1}
    assert cost.windows["start"].to_dict() == {1: 2, 2: 1, 3: 1, 4: 1}  # earliest best start
    assert cost.windows["end"].to_dict() == {1: 2, 2: 2, 3: 3, 4: 4}
    assert cost.windows["obs"].to_dict() == {1: 4, 2: 6, 3: 3, 4: 4}
    big = cost.largest
    assert (big.start, big.end, big.n_units, big.n_periods, big.n_obs) == (1, 2, 3, 2, 6)
    assert cost.windows["obs_share"].iloc[1] == pytest.approx(6 / 13)
    assert "6 observations" in cost.summary("en")
    assert "6 observations" in cost.summary("fr")


def test_the_largest_block_is_the_one_balance_keeps() -> None:
    panel = declare(HAND)
    kept = panel.balance("largest_block")
    big = panel.balance_cost().largest
    assert (kept.n_units, kept.n_obs) == (big.n_units, big.n_obs)
    assert (kept.t_min, kept.t_max) == (big.start, big.end)


def test_a_balanced_panel_costs_nothing() -> None:
    cost = declare({"A": [1, 2, 3], "B": [1, 2, 3]}).balance_cost()
    assert cost.windows["units"].tolist() == [2, 2, 2] and cost.largest.n_obs == 6
    assert cost.thresholds["obs_share"].tolist() == [1.0, 1.0, 1.0]


def window_oracle(units: dict[int, list[int]]) -> tuple[dict[int, int], tuple[int, int, int]]:
    """Best units by window length and the largest block, by enumerating every window."""
    low = min(min(ts) for ts in units.values())
    high = max(max(ts) for ts in units.values())
    best: dict[int, int] = {}
    block = (0, 0, 0)  # (observations, length, -start): the order balance() maximises
    for length in range(1, high - low + 2):
        for start in range(low, high - length + 2):
            window = range(start, start + length)
            n = sum(all(t in set(ts) for t in window) for ts in units.values())
            best[length] = max(best.get(length, 0), n)
            if n and (n * length, length, -start) > block:
                block = (n * length, length, -start)
    return best, block


panels = st.dictionaries(
    st.integers(0, 8),
    st.sets(st.integers(0, 14), min_size=1, max_size=15).map(sorted),
    min_size=1,
    max_size=9,
)


@quiet
@settings(max_examples=100, deadline=None)
@given(units=panels)
def test_balance_cost_equals_window_enumeration_and_balance(units) -> None:  # type: ignore[no-untyped-def]
    panel = declare(units)
    cost = panel.balance_cost()
    best, block = window_oracle(units)
    assert cost.windows["units"].to_dict() == best
    big = cost.largest
    assert (big.n_obs, big.n_periods, -big.start) == (block[0], block[1], block[2])
    kept = panel.balance("largest_block")
    assert (kept.n_units, kept.n_obs) == (big.n_units, big.n_obs)
    assert (kept.t_min, kept.t_max) == (big.start, big.end)
    counts = [len(ts) for ts in units.values()]
    for k, row in cost.thresholds.iterrows():
        assert row["units"] == sum(c >= k for c in counts)
        assert row["obs"] == sum(c for c in counts if c >= k)
