"""Gap typology and imbalance measures (K.2.1, K.2.2): hand-computed, loop oracle, invariances."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pyxtpanel as xt
from pyxtpanel.explore._spec import TOOLS
from pyxtpanel.explore.balance import unbalancedness
from pyxtpanel.explore.gaps import gap_types

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")

#: Calendar 1..4.  A complete; B enters late; C leaves early; D has one gap; E enters late + gap.
HAND = {
    "A": [1, 2, 3, 4],
    "B": [3, 4],
    "C": [1, 2],
    "D": [1, 2, 4],
    "E": [2, 4],
}


def declare(units: dict[object, list[int]], **kw: object) -> xt.PanelData:
    rows = [(u, t) for u, ts in units.items() for t in ts]
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t"]), "id", "t", **kw)  # type: ignore[arg-type]


# ----------------------------------------------------------------- gap_types, by hand


def test_gap_types_by_hand() -> None:
    result = declare(HAND).gap_types()
    table = result.by_unit
    assert table["type"].to_dict() == {
        "A": "complete",
        "B": "late_entry",
        "C": "early_exit",
        "D": "internal_gap",
        "E": "late_entry+internal_gap",
    }
    assert table["lead_missing"].to_dict() == {"A": 0, "B": 2, "C": 0, "D": 0, "E": 1}
    assert table["trail_missing"].to_dict() == {"A": 0, "B": 0, "C": 2, "D": 0, "E": 0}
    assert table["n_gaps"].to_dict() == {"A": 0, "B": 0, "C": 0, "D": 1, "E": 1}
    assert result.n_periods == 4 and (result.first, result.last) == (1, 4)
    assert result.counts["units"].to_dict() == {
        "complete": 1,
        "early_exit": 1,
        "internal_gap": 1,
        "late_entry": 1,
        "late_entry+internal_gap": 1,
    }
    assert result.counts["share"].sum() == pytest.approx(1.0)


def test_the_gap_table_gives_the_period_of_reentry_and_the_length() -> None:
    gaps = declare({"A": [1, 2, 6], "B": [1, 2, 3, 4, 5, 6], "C": [1, 4, 6]}).gap_types().gaps
    assert gaps.to_dict("records") == [
        {"unit": "A", "after": 2, "resume": 6, "length": 3},
        {"unit": "C", "after": 1, "resume": 4, "length": 2},
        {"unit": "C", "after": 4, "resume": 6, "length": 1},
    ]


def test_length_distribution_counts_each_length_once() -> None:
    result = declare({"A": [1, 2, 6], "C": [1, 4, 6], "D": [1, 3, 5]}).gap_types()
    assert result.length_distribution["gaps"].to_dict() == {1: 3, 2: 1, 3: 1}


def test_a_balanced_panel_is_all_complete_with_no_gap() -> None:
    result = declare({"A": [1, 2, 3], "B": [1, 2, 3]}).gap_types()
    assert result.counts["units"].to_dict() == {"complete": 2}
    assert result.gaps.empty and "No internal gap" in result.summary("en")
    assert "Aucun trou interne" in result.summary("fr")


def test_gap_types_uses_the_calendar_step_not_the_integers() -> None:
    quarters = declare({"A": [0, 4, 8, 12], "B": [4, 12]}, delta=4)
    result = quarters.gap_types()
    assert result.by_unit.loc["B", "type"] == "late_entry+internal_gap"
    assert result.gaps.to_dict("records") == [{"unit": "B", "after": 4, "resume": 12, "length": 1}]
    assert result.n_periods == 4


def test_gap_types_on_dated_panels_gives_dates() -> None:
    frame = pd.DataFrame(
        {"id": [1, 1, 1], "t": pd.to_datetime(["2020-01-01", "2020-02-01", "2020-05-01"])}
    )
    result = xt.xtset(frame, "id", "t", freq="M").gap_types()
    assert result.gaps.loc[0, "after"] == pd.Timestamp("2020-02-01")
    assert result.gaps.loc[0, "length"] == 2


# ----------------------------------------------------------------- unbalancedness, by hand


def test_unbalancedness_by_hand_and_against_the_published_formulas() -> None:
    result = declare(HAND).unbalancedness()
    # T_i = 4, 2, 2, 3, 2; N = 5; n = 13; Tbar = 2.6
    assert result.gamma == pytest.approx(12 / 13, abs=1e-15)
    assert result.nu == pytest.approx(169 / 185, abs=1e-15)
    assert (result.n_units, result.n_obs, result.n_periods) == (5, 13, 4)
    assert result.t_mean == pytest.approx(2.6) and (result.t_min, result.t_max) == (2, 4)
    assert result.t_harmonic == pytest.approx(5 / (0.25 + 0.5 + 0.5 + 1 / 3 + 0.5))
    assert result.fill_rate == pytest.approx(13 / 20)


def test_unbalancedness_of_a_balanced_panel_is_one() -> None:
    result = declare({"A": [1, 2, 3], "B": [1, 2, 3]}).unbalancedness()
    assert (result.gamma, result.nu, result.fill_rate) == (1.0, 1.0, 1.0)
    assert "balanced" in result.summary("en") and "cylindré" in result.summary("fr")


def test_unbalancedness_stays_in_its_bounds() -> None:
    result = declare({"A": list(range(1, 11)), "B": [1], "C": [1], "D": [1]}).unbalancedness()
    assert 0 < result.gamma < 1 and 1 / 4 < result.nu < 1


def test_the_results_have_a_table_form() -> None:
    panel = declare(HAND)
    assert list(panel.unbalancedness().to_frame().index[:2]) == ["gamma", "nu"]
    assert "units" in panel.gap_types().to_markdown()
    assert "tabular" in panel.gap_types().to_latex()


# ----------------------------------------------------------------- independent loop oracle


def oracle(units: dict[int, list[int]]) -> dict[str, object]:
    """Gap typology and imbalance written with plain loops and sets (no numpy)."""
    low = min(min(ts) for ts in units.values())
    high = max(max(ts) for ts in units.values())
    types, gaps = {}, []
    for unit, ts in units.items():
        seen = set(ts)
        flags = []
        if min(ts) > low:
            flags.append("late_entry")
        if max(ts) < high:
            flags.append("early_exit")
        missing = [t for t in range(min(ts), max(ts) + 1) if t not in seen]
        if missing:
            flags.append("internal_gap")
        types[unit] = "+".join(flags) or "complete"
        run = 0
        for t in range(min(ts), max(ts) + 1):
            if t in seen:
                if run:
                    gaps.append((unit, t - run - 1, t, run))
                run = 0
            else:
                run += 1
    counts = [len(ts) for ts in units.values()]
    n, total = len(counts), sum(counts)
    tbar = total / n
    gamma = n / (tbar * sum(1 / c for c in counts))
    nu = 1 / (n * sum((c / total) ** 2 for c in counts))
    return {"types": types, "gaps": sorted(gaps), "gamma": gamma, "nu": nu}


panels = st.dictionaries(
    st.integers(0, 8),
    st.sets(st.integers(0, 14), min_size=1, max_size=15).map(sorted),
    min_size=1,
    max_size=9,
)


@quiet
@settings(max_examples=120, deadline=None)
@given(units=panels)
def test_the_vectorised_versions_equal_the_loop_implementation(units) -> None:  # type: ignore[no-untyped-def]
    panel = declare(units)
    expected = oracle(units)
    result = gap_types(panel)
    assert result.by_unit["type"].to_dict() == expected["types"]
    got = sorted(map(tuple, result.gaps[["unit", "after", "resume", "length"]].to_numpy().tolist()))
    assert got == expected["gaps"]
    measures = unbalancedness(panel)
    assert measures.gamma == pytest.approx(expected["gamma"], rel=1e-13)  # type: ignore[arg-type]
    assert measures.nu == pytest.approx(expected["nu"], rel=1e-13)  # type: ignore[arg-type]


@quiet
@settings(max_examples=60, deadline=None)
@given(units=panels, shift=st.integers(-500, 500), seed=st.integers(0, 10_000))
def test_results_do_not_depend_on_row_order_unit_names_or_the_origin_of_time(  # type: ignore[no-untyped-def]
    units, shift, seed
) -> None:
    base = declare(units)
    rows = [(f"unit_{u:02d}", t + shift) for u, ts in units.items() for t in ts]
    frame = pd.DataFrame(rows, columns=["id", "t"])
    frame = frame.iloc[np.random.default_rng(seed).permutation(len(frame))]
    other = xt.xtset(frame, "id", "t")
    a, b = base.gap_types(), other.gap_types()
    assert a.by_unit["type"].tolist() == b.by_unit["type"].tolist()
    assert a.gaps["length"].tolist() == b.gaps["length"].tolist()
    assert (a.gaps["after"] + shift).tolist() == b.gaps["after"].tolist()
    assert a.counts["units"].tolist() == b.counts["units"].tolist()
    u, v = base.unbalancedness(), other.unbalancedness()
    assert (u.gamma, u.nu, u.fill_rate) == pytest.approx((v.gamma, v.nu, v.fill_rate))


# ----------------------------------------------------------------- registry


def test_the_tools_state_their_validation_and_equivalents() -> None:
    assert {"gap_types", "unbalancedness"} <= set(TOOLS.names())
    assert TOOLS.spec("unbalancedness").similar == "plm::punbalancedness"
    for _name, spec in TOOLS.specs():
        assert spec.validation and spec.summary


def test_clean_panels_raise_no_warning() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        declare(HAND).gap_types()
        declare(HAND).unbalancedness()
