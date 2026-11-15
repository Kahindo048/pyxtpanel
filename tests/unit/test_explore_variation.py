"""Variation and residual variation (K.3.1, K.3.3): by hand, loop oracle, invariances."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pyxtpanel as xt
from pyxtpanel._errors import SpecificationError

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")


def hand() -> xt.PanelData:
    # k constant; f time-invariant; c common to all units; v varies; unit D has one observation
    rows = []
    for unit, f, v in (("A", 1, [1, 2, 3]), ("B", 2, [2, 2, 5]), ("C", 3, [1, 1, 1])):
        for t in (1, 2, 3):
            rows.append((unit, t, 4.0, float(f), 10.0 * t, float(v[t - 1])))
    rows.append(("D", 1, 4.0, 4.0, 10.0, 9.0))
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "k", "f", "c", "v"]), "id", "t")


def test_variation_by_hand() -> None:
    table = hand().variation().table
    assert table["kind"].to_dict() == {
        "k": "constant",
        "f": "time_invariant",
        "c": "common",
        "v": "varying",
    }
    assert table["n_distinct"].to_dict() == {"k": 1, "f": 4, "c": 3, "v": 5}
    assert table["units_varying"].to_dict() == {"k": 0, "f": 0, "c": 3, "v": 2}
    assert table["units_single"].to_dict() == {"k": 1, "f": 1, "c": 1, "v": 1}
    assert table["periods_varying"].to_dict() == {"k": 0, "f": 3, "c": 0, "v": 3}
    assert table.loc["v", "share_units_varying"] == pytest.approx(2 / 4)
    assert table["n_obs"].unique().tolist() == [10]


def test_variation_takes_a_list_and_skips_the_keys_and_text() -> None:
    frame = pd.DataFrame({"id": [1, 1, 2], "t": [1, 2, 1], "x": [1.0, 2, 3], "s": ["a", "b", "c"]})
    panel = xt.xtset(frame, "id", "t")
    assert panel.variation().table.index.tolist() == ["x"]
    assert panel.variation(["x"]).summary("fr").startswith("Variation de 1 variable(s)")
    with pytest.raises(SpecificationError):
        xt.xtset(frame[["id", "t", "s"]], "id", "t").variation()


def test_residual_variation_by_hand() -> None:
    frame = pd.DataFrame(
        {
            "id": ["A"] * 3 + ["B"] * 3,
            "t": [1, 2, 3] * 2,
            "x": [1.0, 2, 3, 10, 10, 10],
        }
    )
    table = xt.xtset(frame, "id", "t").residual_variation(["x"]).table
    total = np.var([1, 2, 3, 10, 10, 10], ddof=1)
    within = np.var([-1, 0, 1, 0, 0, 0], ddof=1)
    assert table.loc[("x", "within"), "sd"] == pytest.approx(math.sqrt(within))
    assert table.loc[("x", "within"), "variance_share"] == pytest.approx(within / total)
    assert (table.loc[("x", "within"), "n"], table.loc[("x", "within"), "units"]) == (6, 1)
    fd = np.var([1, 1, 0, 0], ddof=1)
    assert table.loc[("x", "fd"), "variance_share"] == pytest.approx(fd / total)
    assert table.loc[("x", "fd"), "n"] == 4
    a = [math.sqrt(2 / 3) * (1 - 2.5), math.sqrt(1 / 2) * (2 - 3)]
    assert table.loc[("x", "fod"), "sd"] == pytest.approx(np.std([*a, 0, 0], ddof=1))
    assert table["usable"].all()


def test_a_time_invariant_variable_has_nothing_left() -> None:
    table = hand().residual_variation(["f", "v"]).table
    for kind in ("within", "fd", "fod"):
        assert not table.loc[("f", kind), "usable"] and table.loc[("f", kind), "units"] == 0
        assert table.loc[("v", kind), "usable"]


def test_residual_variation_refuses_an_unknown_transformation() -> None:
    with pytest.raises(SpecificationError) as bad:
        hand().residual_variation(["v"], ["twoway"])
    assert bad.value.stata_code == 198


# ----------------------------------------------------------------- loop oracle


def loop_variation(data: dict[int, dict[int, float]]) -> dict[str, object]:
    values = [y for rows in data.values() for y in rows.values() if not math.isnan(y)]
    per_unit = {u: [y for y in rows.values() if not math.isnan(y)] for u, rows in data.items()}
    per_unit = {u: ys for u, ys in per_unit.items() if ys}
    periods = sorted({t for rows in data.values() for t in rows})
    per_period = {
        t: [data[u][t] for u in data if t in data[u] and not math.isnan(data[u][t])]
        for t in periods
    }
    per_period = {t: ys for t, ys in per_period.items() if ys}
    units_varying = sum(len(set(ys)) > 1 for ys in per_unit.values())
    periods_varying = sum(len(set(ys)) > 1 for ys in per_period.values())
    if len(set(values)) <= 1:
        kind = "constant"
    elif units_varying == 0:
        kind = "time_invariant"
    elif periods_varying == 0:
        kind = "common"
    else:
        kind = "varying"
    return {
        "n_distinct": len(set(values)),
        "units_varying": units_varying,
        "units_single": sum(len(ys) == 1 for ys in per_unit.values()),
        "periods_varying": periods_varying,
        "kind": kind,
    }


def loop_residual(data: dict[int, dict[int, float]], kind: str) -> tuple[int, int, float]:
    """(n, units with a non-zero value, sum of squares and count for the sd) by loops."""
    out: list[tuple[int, float]] = []
    for u, rows in data.items():
        seen = [(t, y) for t, y in sorted(rows.items()) if not math.isnan(y)]
        mean = sum(y for _, y in seen) / len(seen) if seen else math.nan
        for t, y in seen:
            if kind == "within":
                out.append((u, y - mean))
            elif kind == "fd":
                if t - 1 in rows and not math.isnan(rows[t - 1]):
                    out.append((u, y - rows[t - 1]))
            else:
                later = [z for s, z in seen if s > t]
                if later:
                    f = len(later)
                    out.append((u, math.sqrt(f / (f + 1)) * (y - sum(later) / f)))
    scale = max(
        [abs(y) for rows in data.values() for y in rows.values() if not math.isnan(y)] + [1.0]
    )
    units = len({u for u, v in out if abs(v) > 1e-12 * scale})
    vals = [v for _, v in out]
    sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else math.nan
    return len(vals), units, sd


cells = st.one_of(st.integers(-3, 3).map(float), st.just(math.nan))
panels = st.dictionaries(
    st.integers(0, 6),
    st.dictionaries(st.integers(0, 8), cells, min_size=1, max_size=9),
    min_size=1,
    max_size=7,
)


def build(data: dict[int, dict[int, float]], name=lambda u: u, shift: int = 0) -> xt.PanelData:  # type: ignore[no-untyped-def]
    rows = [(name(u), t + shift, y) for u, rs in data.items() for t, y in rs.items()]
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "x"]), "id", "t")


@quiet
@settings(max_examples=120, deadline=None)
@given(data=panels)
def test_variation_equals_the_loop_implementation(data) -> None:  # type: ignore[no-untyped-def]
    got = build(data).variation(["x"]).table.loc["x"]
    expected = loop_variation(data)
    for key, value in expected.items():
        assert got[key] == value, key


@quiet
@settings(max_examples=120, deadline=None)
@given(data=panels, kind=st.sampled_from(["within", "fd", "fod"]))
def test_residual_variation_equals_the_loop_implementation(data, kind) -> None:  # type: ignore[no-untyped-def]
    row = build(data).residual_variation(["x"], [kind]).table.loc[("x", kind)]
    n, units, sd = loop_residual(data, kind)
    assert (row["n"], row["units"]) == (n, units)
    if math.isnan(sd):
        assert math.isnan(row["sd"])
    else:
        assert row["sd"] == pytest.approx(sd, abs=1e-9)


@quiet
@settings(max_examples=50, deadline=None)
@given(data=panels, shift=st.integers(-200, 200))
def test_results_do_not_depend_on_unit_names_or_the_origin_of_time(data, shift) -> None:  # type: ignore[no-untyped-def]
    a, b = build(data), build(data, name=lambda u: f"u{u:02d}", shift=shift)
    pd.testing.assert_frame_equal(a.variation(["x"]).table, b.variation(["x"]).table)
    pd.testing.assert_frame_equal(
        a.residual_variation(["x"]).table, b.residual_variation(["x"]).table
    )
