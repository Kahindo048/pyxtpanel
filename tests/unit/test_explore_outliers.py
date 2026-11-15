"""Within outliers and level shifts (K.4): by hand, loop oracle, invariances."""

from __future__ import annotations

import math
import statistics

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pyxtpanel as xt
from pyxtpanel._errors import SpecificationError

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")
nan = math.nan


def frame_of(units: dict[str, list[float]]) -> pd.DataFrame:
    rows = [(u, t + 1, x) for u, xs in units.items() for t, x in enumerate(xs)]
    return pd.DataFrame(rows, columns=["id", "t", "x"])


# ----------------------------------------------------------------- within_outliers, by hand


def test_within_outliers_by_hand() -> None:
    panel = xt.xtset(
        frame_of(
            {
                "A": [10, 11, 10, 12, 11, 50],  # median 11, MAD 1
                "B": [5, 5, 5, 5, 5, 9],  # MAD 0: the mean absolute deviation is used
                "C": [3, 3, 3, 3, 3, 3],  # all equal: z = 0
                "D": [1, 100, 2],  # fewer than min_obs values: not scored
            }
        ),
        "id",
        "t",
    )
    result = panel.within_outliers("x")
    assert (result.n_scored, result.n_flagged) == (18, 2)
    assert result.rows["unit"].tolist() == ["A", "B"]  # largest |score| first
    scores = result.scores.to_numpy()
    assert scores[5] == pytest.approx(0.6745 * 39)
    assert scores[:5] == pytest.approx([-0.6745, 0, -0.6745, 0.6745, 0], abs=1e-12)
    assert scores[11] == pytest.approx(4 / (1.253314 * 4 / 6))
    assert (scores[12:18] == 0).all() and np.isnan(scores[18:]).all()
    assert result.flagged.sum() == 2 and result.flagged.dtype == bool
    row = result.rows.set_index("unit")
    assert row.loc["A", "value"] == 50 and row.loc["A", "median"] == 11
    assert result.by_unit["n_flagged"].to_dict() == {"A": 1, "B": 1, "C": 0, "D": 0}
    assert "2 of 18" in result.summary("en") and "2 valeurs sur 18" in result.summary("fr")


def test_within_outliers_scores_the_values_in_the_order_of_the_declared_frame() -> None:
    frame = frame_of({"A": [10, 11, 10, 12, 11, 50]}).iloc[::-1]
    result = xt.xtset(frame, "id", "t").within_outliers("x")
    assert result.scores.index.equals(frame.index)
    assert result.flagged.to_numpy()[0] and result.flagged.sum() == 1


def test_within_outliers_refuse_bad_arguments() -> None:
    panel = xt.xtset(frame_of({"A": [1, 2, 3, 4, 5, 6]}), "id", "t")
    for kwargs in ({"threshold": 0}, {"min_obs": 1}):
        with pytest.raises(SpecificationError) as bad:
            panel.within_outliers("x", **kwargs)
        assert bad.value.stata_code == 198


def loop_outliers(
    units: dict[int, list[float]], threshold: float, min_obs: int
) -> dict[tuple[int, int], float]:
    out = {}
    for u, xs in units.items():
        seen = [x for x in xs if not math.isnan(x)]
        for t, x in enumerate(xs):
            if math.isnan(x) or len(seen) < min_obs:
                out[(u, t)] = nan
                continue
            med = statistics.median(seen)
            dev = [abs(y - med) for y in seen]
            mad, mean_ad = statistics.median(dev), sum(dev) / len(dev)
            if mad > 0:
                out[(u, t)] = 0.6745 * (x - med) / mad
            elif mean_ad > 0:
                out[(u, t)] = (x - med) / (1.253314 * mean_ad)
            else:
                out[(u, t)] = 0.0
    return out


cells = st.one_of(st.integers(-5, 5).map(float), st.just(nan))
panels = st.dictionaries(
    st.integers(0, 5), st.lists(cells, min_size=1, max_size=10), min_size=1, max_size=6
)


def build(units: dict[int, list[float]], name=lambda u: u, shift: int = 0) -> xt.PanelData:  # type: ignore[no-untyped-def]
    rows = [(name(u), t + shift, x) for u, xs in units.items() for t, x in enumerate(xs)]
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "x"]), "id", "t")


@quiet
@settings(max_examples=150, deadline=None)
@given(units=panels, threshold=st.sampled_from([1.0, 3.5]), min_obs=st.integers(2, 5))
def test_outliers_equal_the_loop_implementation(units, threshold, min_obs) -> None:  # type: ignore[no-untyped-def]
    frame = build(units).to_frame().reset_index(drop=True)
    result = xt.xtset(frame, "id", "t").within_outliers("x", threshold=threshold, min_obs=min_obs)
    expected = loop_outliers(units, threshold, min_obs)
    got = dict(
        zip(zip(frame["id"], frame["t"], strict=True), result.scores.to_numpy(), strict=True)
    )
    for key, value in expected.items():
        np.testing.assert_allclose(got[key], value, rtol=1e-12, atol=1e-12)
    finite = [v for v in expected.values() if not math.isnan(v)]
    assert result.n_scored == len(finite)
    assert result.n_flagged == sum(abs(v) > threshold for v in finite)


@quiet
@settings(max_examples=50, deadline=None)
@given(units=panels, shift=st.integers(-100, 100), scale=st.sampled_from([0.5, 2.0, 10.0]))
def test_scores_ignore_names_origin_of_time_and_scale(units, shift, scale) -> None:  # type: ignore[no-untyped-def]
    a = build(units).within_outliers("x", min_obs=2)
    scaled = {u: [x * scale + 7 for x in xs] for u, xs in units.items()}
    b = build(scaled, name=lambda u: f"u{u:02d}", shift=shift).within_outliers("x", min_obs=2)
    np.testing.assert_allclose(
        a.scores.to_numpy(), b.scores.to_numpy(), rtol=1e-9, atol=1e-9, equal_nan=True
    )
    assert a.n_flagged == b.n_flagged


# ----------------------------------------------------------------- level_shifts, by hand


def test_level_shifts_by_hand() -> None:
    panel = xt.xtset(
        frame_of(
            {
                "A": [0, 0, 0, 0, 10, 10, 10, 10],  # one lasting jump; scale 0: fallback to B's
                "B": [0, 1, 3, 2, 4, 3, 5, 4],  # diffs 1 2 -1 2 -1 2 -1: MAD 1, scale 1.4826
            }
        ),
        "id",
        "t",
    )
    result = panel.level_shifts("x", window=1)
    table = result.by_unit
    scale = 1.4826
    assert table.loc["A", "max_score"] == pytest.approx(10 / scale)
    assert table.loc["A", "shift"] == 10 and table.loc["A", "period"] == 5
    assert table.loc["B", "max_score"] == pytest.approx(2 / scale)
    assert table["flagged"].to_dict() == {"A": True, "B": False}
    assert (result.n_flagged, result.n_scored) == (1, 2)
    assert table["within_share"].sum() == pytest.approx(1.0)
    assert table.loc["A", "within_share"] == pytest.approx(200 / 219.5)
    assert not table["influential"].any()  # 0.91 is below 3 times the average 0.5
    assert panel.level_shifts("x", window=1, influence=1.5).by_unit["influential"].tolist() == [
        True,
        False,
    ]
    assert "1 of 2 units" in result.summary("en") and "1 unités sur 2" in result.summary("fr")


def test_a_single_spike_is_not_a_level_shift() -> None:
    values = [1, 2, 1, 2, 1, 40, 1, 2, 1, 2, 1, 2]  # a spike at the 6th value
    shifted = [1, 2, 1, 2, 1, 2, 40, 41, 40, 41, 40, 41]  # the level moves for good
    table = (
        xt.xtset(frame_of({"spike": values, "shift": shifted}), "id", "t")
        .level_shifts("x", window=3)
        .by_unit
    )
    assert table.loc["spike", "max_score"] < 1.0 < table.loc["shift", "max_score"]


def test_short_units_are_not_scored() -> None:
    table = (
        xt.xtset(frame_of({"A": [1, 2, 3, 4, 5], "B": [1, 2]}), "id", "t")
        .level_shifts("x", window=3)
        .by_unit
    )
    assert table["max_score"].isna().all() and not table["flagged"].any()


def loop_shifts(units: dict[int, list[float]], window: int) -> dict[int, float]:
    scales: dict[int, float] = {}
    for u, xs in units.items():
        v = [x for x in xs if not math.isnan(x)]
        if len(v) >= max(2 * window, 2):
            d = [b - a for a, b in zip(v, v[1:], strict=False)]
            m = statistics.median(d)
            scales[u] = 1.4826 * statistics.median([abs(y - m) for y in d])
    positive = [s for s in scales.values() if s > 0]
    fallback = statistics.median(positive) if positive else nan
    out: dict[int, float] = {}
    for u, xs in units.items():
        v = [x for x in xs if not math.isnan(x)]
        if u not in scales:
            out[u] = nan
            continue
        used = scales[u] if scales[u] > 0 else fallback
        if not used > 0:
            out[u] = nan
            continue
        best = 0.0
        for j in range(window, len(v) - window + 1):
            shift = statistics.median(v[j : j + window]) - statistics.median(v[j - window : j])
            best = max(best, abs(shift) / used)
        out[u] = best
    return out


@quiet
@settings(max_examples=120, deadline=None)
@given(units=panels, window=st.integers(1, 3))
def test_level_shifts_equal_the_loop_implementation(units, window) -> None:  # type: ignore[no-untyped-def]
    table = build(units).level_shifts("x", window=window).by_unit
    expected = loop_shifts(units, window)
    got = dict(zip(table.index.tolist(), table["max_score"].to_numpy(), strict=True))
    for u, value in expected.items():
        np.testing.assert_allclose(got[u], value, rtol=1e-10, atol=1e-10, equal_nan=True)


@quiet
@settings(max_examples=50, deadline=None)
@given(units=panels, shift=st.integers(-100, 100))
def test_level_shifts_ignore_names_and_the_origin_of_time(units, shift) -> None:  # type: ignore[no-untyped-def]
    a = build(units).level_shifts("x", window=2).by_unit
    b = build(units, name=lambda u: f"u{u:02d}", shift=shift).level_shifts("x", window=2).by_unit
    np.testing.assert_allclose(a["max_score"].to_numpy(), b["max_score"].to_numpy(), equal_nan=True)
    assert a["flagged"].tolist() == b["flagged"].tolist()
