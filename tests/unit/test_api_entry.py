"""as_panel, PanelLike, display and export of results (G.1.3, G.3)."""

import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel._errors import SpecificationError
from pyxtpanel.panel.base import PanelBase


def frame() -> pd.DataFrame:
    return pd.DataFrame({"id": [1, 1, 2], "t": [1, 2, 2], "x": [1.0, 2, 3]})


def test_as_panel_accepts_panels_frames_and_indexed_frames() -> None:
    panel = xt.xtset(frame(), "id", "t")
    assert xt.as_panel(panel) is panel
    assert xt.as_panel(panel, "id", "t") is panel
    built = xt.as_panel(frame(), "id", "t", delta=1)
    assert isinstance(built, xt.PanelData) and built.n_obs == 3
    indexed = frame().set_index(["id", "t"])
    assert xt.as_panel(indexed).n_units == 2


def test_as_panel_upgrades_a_base_panel_and_refuses_contradictions() -> None:
    base = PanelBase.from_frame(frame(), "id", "t")
    upgraded = xt.as_panel(base)
    assert isinstance(upgraded, xt.PanelData) and upgraded.index is base.index
    with pytest.raises(SpecificationError, match="declared by"):
        xt.as_panel(upgraded, "t", "id")
    with pytest.raises(SpecificationError, match="two-level index"):
        xt.as_panel(frame())
    with pytest.raises(SpecificationError, match="expects a PanelData"):
        xt.as_panel([1, 2, 3])  # type: ignore[arg-type]
    xt.set_options(lang="fr")
    try:
        with pytest.raises(SpecificationError, match="attend un PanelData"):
            xt.as_panel("x")  # type: ignore[arg-type]
    finally:
        xt.set_options(lang="en")


def test_panel_data_satisfies_the_panel_like_protocol() -> None:
    panel = xt.xtset(frame(), "id", "t")
    assert isinstance(panel, xt.PanelLike)
    assert not isinstance(frame(), xt.PanelLike)


def test_html_repr_of_a_panel_and_of_results() -> None:
    panel = xt.xtset(frame(), "id", "t")
    html = panel._repr_html_()
    assert "PanelData" in html and "<table" in html and "2 units" in html
    assert "<table" in panel.xtsum("x")._repr_html_()


def test_results_export_to_markdown_and_latex() -> None:
    panel = xt.xtset(frame().assign(d=[0.0, 1, 1]), "id", "t")
    result = panel.xttab("d")
    markdown = result.to_markdown()
    assert markdown.splitlines()[0].startswith("| value") and "Total" in markdown
    latex = result.to_latex()
    assert latex.startswith(r"\begin{tabular}") and latex.rstrip().endswith(r"\end{tabular}")
    assert r"overall\_freq" in latex  # underscores are escaped
    assert panel.attrition().to_markdown().count("\n") >= 2
