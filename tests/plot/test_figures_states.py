"""States map, attrition and balancing-cost figures: drawn data equals computed (F.2.7, F.2.8)."""

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")

import pyxtpanel as xt  # noqa: E402
from pyxtpanel._errors import SpecificationError  # noqa: E402
from pyxtpanel.plot.kinds.states import states_data  # noqa: E402

nan = np.nan


def hand() -> xt.PanelData:
    rows = []
    for unit, periods, values in (
        ("A", [1, 2, 3, 4], [0, 0, 1, 1]),
        ("B", [2, 3, 4], [0, 0, 0]),
        ("C", [1, 2, 3, 4], [0, 0, nan, 1]),
        ("D", [1, 2], [0, 0]),
    ):
        rows += [(unit, t, v, float(t)) for t, v in zip(periods, values, strict=True)]
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "d", "x"]), "id", "t")


def texts(fig) -> list[str]:  # type: ignore[no-untyped-def]
    return [t.get_text() for ax in fig.axes for t in ax.texts] + [
        ax.get_title(loc) for ax in fig.axes for loc in ("center", "left")
    ]


# ---------------------------------------------------------------- states


def test_the_state_map_by_hand() -> None:
    data = states_data(hand(), "d")
    # -2 not observed, -1 value missing, 0 and 1 the states; A and C change, earliest first
    np.testing.assert_array_equal(
        data.grid,
        [[0, 0, 1, 1], [0, 0, -1, 1], [-2, 0, 0, 0], [0, 0, -2, -2]],
    )
    assert data.order.tolist() == [0, 2, 1, 3] and data.states.tolist() == [0.0, 1.0]
    assert data.n_changing == 2 and data.ranks.tolist() == [0, 1, 2, 3]
    by_unit = states_data(hand(), "d", sort="unit")
    assert by_unit.order.tolist() == [0, 1, 2, 3]


def test_the_state_figure_draws_the_grid_and_names_missing_and_absent_cells() -> None:
    panel = hand()
    fig = panel.plot("states", variable="d")
    ax = fig.axes[0]
    np.testing.assert_array_equal(
        np.asarray(ax.images[0].get_array()), states_data(panel, "d").grid + 2
    )
    assert [t.get_text() for t in ax.get_yticklabels()] == ["A", "C", "B", "D"]
    legend = [t.get_text() for t in ax.get_legend().get_texts()]
    assert legend == ["0", "1", "value missing", "not observed"]
    assert "d: 2 of 4 units change state" in texts(fig)
    french = " ".join(texts(panel.plot("states", variable="d", lang="fr")))
    assert "2 unités sur 4 changent d'état" in french
    assert [
        t.get_text()
        for t in panel.plot("states", variable="d", lang="fr").axes[0].get_legend().get_texts()
    ][2:] == [
        "valeur manquante",
        "non observée",
    ]


def test_the_legend_only_lists_what_is_drawn() -> None:
    balanced = xt.xtset(
        pd.DataFrame({"id": [1, 1, 2, 2], "t": [1, 2, 1, 2], "d": [0.0, 1, 0, 0]}), "id", "t"
    )
    legend = [
        t.get_text() for t in balanced.plot("states", variable="d").axes[0].get_legend().get_texts()
    ]
    assert legend == ["0", "1"]


def test_state_figure_errors() -> None:
    panel = hand()
    with pytest.raises(SpecificationError) as bad:
        panel.plot("states", variable="x", max_states=3)
    assert bad.value.stata_code == 134
    with pytest.raises(SpecificationError) as wrong:
        panel.plot("states", variable="d", sort="random")
    assert wrong.value.stata_code == 198


# ---------------------------------------------------------------- attrition


def test_attrition_figure_draws_the_computed_tables() -> None:
    panel = hand()
    result = panel.attrition()
    fig = panel.plot("attrition")
    top, bottom, side = fig.axes
    np.testing.assert_array_equal(top.lines[0].get_ydata(), result.by_period["active"].to_numpy())
    heights = [np.asarray(p.get_paths()[0].vertices)[:, 1] for p in bottom.collections]
    assert max(heights[0]) == result.by_period["entries"].max()
    assert min(heights[1]) == -result.by_period["exits"].max()
    survival = result.survival
    x, y = side.lines[0].get_data()
    np.testing.assert_array_equal(np.asarray(x), np.concatenate([[0], survival.index.to_numpy()]))
    np.testing.assert_allclose(
        np.asarray(y), np.concatenate([[1.0], survival["survival"].to_numpy()])
    )
    assert (
        f"{int(survival['censored'].sum())} units still observed in the last period (censored)"
        in texts(fig)
    )
    assert "entrées" in " ".join(texts(panel.plot("attrition", lang="fr")))


# ---------------------------------------------------------------- balance cost


def test_balance_cost_figure_draws_the_computed_curves_and_marks_the_largest_block() -> None:
    panel = hand()
    cost = panel.balance_cost()
    fig = panel.plot("balance_cost")
    left, right = fig.axes
    np.testing.assert_allclose(left.lines[0].get_ydata(), cost.thresholds["units_share"])
    np.testing.assert_allclose(left.lines[1].get_ydata(), cost.thresholds["obs_share"])
    np.testing.assert_allclose(right.lines[0].get_ydata(), cost.windows["units_share"])
    np.testing.assert_allclose(right.lines[1].get_ydata(), cost.windows["obs_share"])
    point = right.collections[0].get_offsets().data[0]
    assert point[0] == cost.largest.n_periods
    assert point[1] == pytest.approx(cost.largest.n_obs / cost.n_obs)
    big = cost.largest
    assert f"largest block: {big.n_units} units x {big.n_periods} periods" in texts(fig)
    french = texts(panel.plot("balance_cost", lang="fr"))
    assert f"plus grand bloc : {big.n_units} unités x {big.n_periods} périodes" in french


def test_the_three_kinds_are_registered_and_exported(tmp_path) -> None:  # type: ignore[no-untyped-def]
    from pyxtpanel.plot import KINDS, save_figure

    assert {"states", "attrition", "balance_cost"} <= set(KINDS.names())
    for kind, options in (("states", {"variable": "d"}), ("attrition", {}), ("balance_cost", {})):
        written = save_figure(
            hand().plot(kind, **options), tmp_path / f"{kind}.png", formats=(".png", ".pdf", ".svg")
        )
        assert all(p.stat().st_size > 500 for p in written)
