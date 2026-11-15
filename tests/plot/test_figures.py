"""Figures: what is drawn equals what is computed, texts, export, no style leak (F)."""

import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")

import pyxtpanel as xt  # noqa: E402
from pyxtpanel._errors import SpecificationError  # noqa: E402
from pyxtpanel.plot import KINDS, plot, save_figure  # noqa: E402
from pyxtpanel.plot._data import period_mean, sample_units, unit_lines  # noqa: E402
from pyxtpanel.plot.kinds.distribution import distribution_data  # noqa: E402
from pyxtpanel.plot.kinds.participation import participation_data  # noqa: E402
from pyxtpanel.plot.kinds.scatter import scatter_data  # noqa: E402
from pyxtpanel.plot.style import BLUE, GREY, ORANGE  # noqa: E402


@pytest.fixture
def panel() -> xt.PanelData:
    rng = np.random.default_rng(5)
    rows = []
    for unit in range(1, 15):
        effect = rng.normal()
        for year in range(2000, 2010):
            if rng.random() < 0.15:
                continue
            z = rng.normal()
            rows.append((unit, year, effect + 0.5 * z + rng.normal(scale=0.5), z))
    frame = pd.DataFrame(rows, columns=["id", "t", "x", "z"])
    frame["d"] = np.select([frame.z < -0.4, frame.z < 0.4], [0.0, 1.0], default=2.0)
    return xt.xtset(frame, "id", "t")


def texts(fig) -> list[str]:  # type: ignore[no-untyped-def]
    return [t.get_text() for ax in fig.axes for t in ax.texts] + [
        ax.get_title(loc) for ax in fig.axes for loc in ("center", "left")
    ]


# ---------------------------------------------------------------- lines


def test_lines_draw_the_computed_trajectories_and_the_period_mean(panel) -> None:  # type: ignore[no-untyped-def]
    fig = panel.plot("lines", variable="x", n_units=5, seed=3, highlight=[1])
    expected_codes = sample_units(panel, None, 5, 3)
    values = panel.arrays(["x"]).values[:, 0]
    expected = unit_lines(panel, values, expected_codes)
    drawn = [line for line in fig.axes[0].lines]
    assert len(drawn) == 5 + 1  # five units and the mean
    for line, unit in zip(drawn[:5], expected, strict=True):
        np.testing.assert_array_equal(line.get_ydata(), unit.y)
    periods, mean = period_mean(panel, values)
    np.testing.assert_allclose(drawn[-1].get_ydata(), mean)
    assert drawn[-1].get_linestyle() == "--"
    assert "5 of 14 units" in texts(fig)


def test_the_unit_sample_is_reproducible_and_limited(panel) -> None:  # type: ignore[no-untyped-def]
    first = sample_units(panel, None, 6, seed=11)
    assert first == sample_units(panel, None, 6, seed=11) and len(first) == 6
    assert first != sample_units(panel, None, 6, seed=12)
    assert sample_units(panel, None, 100, seed=0) == list(range(14))
    assert sample_units(panel, [3, 5], 1, 0) == [2, 4]


def test_highlighted_units_are_coloured_and_named(panel) -> None:  # type: ignore[no-untyped-def]
    fig = panel.plot("lines", variable="x", n_units=14, highlight=[2, 4], mean=False)
    colours = {line.get_color() for line in fig.axes[0].lines}
    assert colours == {GREY, BLUE, ORANGE}
    names = [t.get_text() for t in fig.axes[0].texts if t.get_text() in ("2", "4")]
    assert sorted(names) == ["2", "4"]


def test_calendar_gaps_break_the_line() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1, 1], "t": [1, 2, 4, 5], "x": [1.0, 2, 3, 4]})
    panel = xt.xtset(frame, "id", "t")
    (line,) = unit_lines(panel, panel.arrays(["x"]).values[:, 0], [0])
    assert np.isnan(line.y).tolist() == [False, False, True, False, False]
    assert line.x.tolist() == [1, 2, 2, 4, 5]


def test_small_multiples_one_panel_per_unit(panel) -> None:  # type: ignore[no-untyped-def]
    fig = panel.plot("lines", variable="x", n_units=6, panels=True)
    visible = [ax for ax in fig.axes if ax.get_visible()]
    assert len(visible) == 6 and all(len(ax.lines) == 1 for ax in visible)


def test_dated_panels_draw_timestamps() -> None:
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1],
            "t": pd.to_datetime(["2020-01-01", "2020-04-01", "2020-07-01"]),
            "x": [1.0, 2, 3],
        }
    )
    panel = xt.xtset(frame, "id", "t")
    fig = panel.plot("lines", variable="x")
    assert np.issubdtype(fig.axes[0].lines[0].get_xdata().dtype, np.datetime64)


# ---------------------------------------------------------------- scatter


def test_scatter_points_and_line_equal_the_computed_ones(panel) -> None:  # type: ignore[no-untyped-def]
    arrays = panel.arrays(["z", "x"])
    z, x = arrays.values[:, 0], arrays.values[:, 1]
    fig = panel.plot("scatter", x="z", y="x")
    assert len(fig.axes) == 2
    for ax, kind in zip(fig.axes, ("between", "within"), strict=True):
        data = scatter_data(panel, z, x, kind)  # type: ignore[arg-type]
        offsets = ax.collections[0].get_offsets()
        np.testing.assert_allclose(offsets[:, 0], data.x)
        np.testing.assert_allclose(offsets[:, 1], data.y)
        line = ax.lines[0]
        np.testing.assert_allclose(
            line.get_ydata(), data.intercept + data.slope * np.asarray(line.get_xdata())
        )
    frame = panel.to_frame()
    means = frame.groupby("id")[["z", "x"]].mean()
    np.testing.assert_allclose(
        scatter_data(panel, z, x, "between").slope, np.polyfit(means.z, means.x, 1)[0]
    )
    assert f"n = {len(frame)}" in " ".join(texts(fig))


def test_scatter_kinds_and_errors(panel) -> None:  # type: ignore[no-untyped-def]
    assert len(panel.plot("scatter", x="z", y="x", part="within").axes) == 1
    with pytest.raises(SpecificationError, match="part must be"):
        panel.plot("scatter", x="z", y="x", part="total")


# ---------------------------------------------------------------- distribution


def test_distribution_histograms_count_units_and_observations(panel) -> None:  # type: ignore[no-untyped-def]
    fig = panel.plot("distribution", variable="x", bins=10)
    between, within = fig.axes
    assert sum(p.get_height() for p in between.patches) == panel.n_units
    assert sum(p.get_height() for p in within.patches) == panel.n_obs
    data = distribution_data(panel, panel.arrays(["x"]).values[:, 0])
    assert data.between.size == panel.n_units and data.within.size == panel.n_obs
    overall = panel.arrays(["x"]).values[:, 0].mean()
    np.testing.assert_allclose(data.within.mean(), overall)  # within values keep the mean
    sd = panel.xtsum("x").scalars.loc["x"]
    np.testing.assert_allclose(np.std(data.between, ddof=1), sd["sd_b"])
    assert len(panel.plot("distribution", variable="x", style="box").axes) == 2
    with pytest.raises(SpecificationError, match="hist"):
        panel.plot("distribution", variable="x", style="violin")


# ---------------------------------------------------------------- participation


def test_participation_image_is_the_presence_matrix_sorted_by_pattern(panel) -> None:  # type: ignore[no-untyped-def]
    fig = panel.plot("participation")
    image = fig.axes[0].images[0].get_array()
    data = participation_data(panel)
    np.testing.assert_array_equal(np.asarray(image), data.present.astype(int))
    assert data.present.shape == (panel.n_units, 10)
    assert data.present.sum() == panel.n_obs
    assert sorted(data.order.tolist()) == list(range(panel.n_units))
    # identical patterns are adjacent: the number of changes of pattern equals patterns - 1
    changes = int(np.any(data.present[1:] != data.present[:-1], axis=1).sum())
    assert changes == data.n_patterns - 1
    # the most frequent pattern comes first
    rows, counts = np.unique(data.present, axis=0, return_counts=True)
    first = tuple(data.present[0])
    assert counts[[tuple(r) == first for r in rows].index(True)] == counts.max()


# ---------------------------------------------------------------- transition


def test_transition_heat_map_shows_the_xttrans_percentages(panel) -> None:  # type: ignore[no-untyped-def]
    fig = panel.plot("transition", variable="d")
    expected = panel.xttrans("d").percent.to_numpy()
    np.testing.assert_allclose(np.asarray(fig.axes[0].images[0].get_array()), expected)
    cells = [t.get_text() for t in fig.axes[0].texts]
    assert cells == [f"{v:.1f}" for v in expected.ravel()]
    strict = panel.plot("transition", variable="d", gaps="strict")
    assert "left out" in " ".join(texts(strict)) and "bridge" not in " ".join(texts(strict))
    assert "bridge a calendar gap" in " ".join(texts(fig))


# ---------------------------------------------------------------- language, export, isolation


def test_french_texts_and_decimal_comma(panel) -> None:  # type: ignore[no-untyped-def]
    fig = panel.plot("transition", variable="d", lang="fr")
    cells = [t.get_text() for t in fig.axes[0].texts]
    assert all("," in c and "." not in c for c in cells)
    assert "Valeur suivante" in fig.axes[0].get_xlabel()
    assert "moyenne de toutes les unités" in texts(panel.plot("lines", variable="x", lang="fr"))
    xt.set_options(lang="fr")
    try:
        assert "Valeur suivante" in panel.plot("transition", variable="d").axes[0].get_xlabel()
    finally:
        xt.set_options(lang="en")
    with pytest.raises(SpecificationError, match="lang"):
        panel.plot("lines", variable="x", lang="de")


def test_export_png_pdf_svg(panel, tmp_path) -> None:  # type: ignore[no-untyped-def]
    fig = panel.plot("lines", variable="x")
    written = save_figure(fig, tmp_path / "f.png", formats=(".png", ".pdf", ".svg"))
    assert [p.suffix for p in written] == [".png", ".pdf", ".svg"]
    assert all(p.stat().st_size > 500 for p in written)
    assert written[0].read_bytes()[:4] == b"\x89PNG" and written[1].read_bytes()[:4] == b"%PDF"
    assert save_figure(fig, tmp_path / "g.svg")[0].suffix == ".svg"
    with pytest.raises(SpecificationError, match="unknown figure format"):
        save_figure(fig, tmp_path / "f.gif")


def test_no_style_leaks_into_matplotlib(panel, tmp_path) -> None:  # type: ignore[no-untyped-def]
    before = dict(matplotlib.rcParams)
    fig = panel.plot("scatter", x="z", y="x")
    save_figure(fig, tmp_path / "a.png")
    assert dict(matplotlib.rcParams) == before


def test_unknown_kind_lists_the_available_ones(panel) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(SpecificationError, match="participation"):
        plot(panel, "pie")
    assert {"lines", "scatter", "distribution", "participation", "transition"} <= set(KINDS.names())


def test_importing_the_package_does_not_import_matplotlib_and_the_error_says_how_to_install() -> (
    None
):
    code = (
        "import sys\n"
        "sys.modules['matplotlib'] = None\n"  # makes 'import matplotlib' fail
        "import pyxtpanel as xt, pandas as pd\n"
        "import pyxtpanel.plot\n"
        "p = xt.xtset(pd.DataFrame({'id': [1, 1], 't': [1, 2], 'x': [1.0, 2]}), 'id', 't')\n"
        "p.xtsum('x')\n"
        "try:\n"
        "    p.plot('lines', variable='x')\n"
        "except ImportError as error:\n"
        "    sys.exit(0 if 'pip install pyxtpanel[plot]' in str(error) else 2)\n"
        "sys.exit(3)\n"
    )
    assert subprocess.run([sys.executable, "-c", code], check=False).returncode == 0
    clean = (
        "import sys, pyxtpanel, pyxtpanel.plot, pyxtpanel.api.panel_data\n"
        "sys.exit(int('matplotlib' in sys.modules))\n"
    )
    assert subprocess.run([sys.executable, "-c", clean], check=False).returncode == 0


def test_the_gallery_is_reproducible(tmp_path) -> None:  # type: ignore[no-untyped-def]
    from tools import galerie

    first = galerie.gallery_panel()
    pd.testing.assert_frame_equal(first, galerie.gallery_panel())
    written = galerie.render(tmp_path)
    assert len(written) == 10 and all(
        p.suffix == ".png" and p.stat().st_size > 1000 for p in written
    )
