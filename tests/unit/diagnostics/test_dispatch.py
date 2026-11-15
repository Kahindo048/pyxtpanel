"""Family functions, PanelData methods and the generic test runner: same results, clear refusals."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel import _config
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics import (
    harris_tzavalis,
    im_pesaran_shin,
    panel_tests,
    pesaran_cd,
    resolve_test,
    run_test,
    wooldridge,
)


@pytest.fixture(scope="module")
def panel() -> xt.PanelData:
    rng = np.random.default_rng(0)
    n, t = 30, 16
    walk = rng.normal(size=(n, t)).cumsum(axis=1)
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n + 1), t),
            "t": np.tile(np.arange(1, t + 1), n),
            "y": walk.ravel(),
            "x": rng.normal(size=n * t),
        }
    )
    return xt.xtset(frame, "id", "t")


def test_the_family_functions_return_what_the_tests_return(panel: xt.PanelData) -> None:
    expected = im_pesaran_shin(panel, "y", lags=1)
    assert xt.unit_root(panel, "y", test="ips", lags=1).statistic == expected.statistic
    assert panel.unit_root("y", test="im_pesaran_shin", lags=1).statistic == expected.statistic
    assert panel.test("ips", "y", lags=1).statistic == expected.statistic
    assert run_test("ips", panel, "y", lags=1).statistic == expected.statistic
    assert (
        xt.cross_section_dependence(panel, "y", test="cd").statistic
        == pesaran_cd(panel, "y").statistic
    )
    assert (
        panel.cross_section_dependence("y", test="pesaran_cd").pvalue
        == pesaran_cd(panel, "y").pvalue
    )
    assert (
        xt.serial_correlation(panel, "y", ["x"], test="xtserial").statistic
        == wooldridge(panel, "y", ["x"]).statistic
    )
    assert panel.serial_correlation("y", ["x"], test="wooldridge").pvalue == pytest.approx(
        wooldridge(panel, "y", ["x"]).pvalue
    )
    assert xt.unit_root(panel, "y", test="HT").statistic == harris_tzavalis(panel, "y").statistic


def test_aliases_resolve_to_names_and_unknown_names_list_the_choices() -> None:
    assert [resolve_test(a) for a in ("ht", "LLC", "ips", "cd", "xtserial", "fisher")] == [
        "harris_tzavalis",
        "levin_lin_chu",
        "im_pesaran_shin",
        "pesaran_cd",
        "wooldridge",
        "fisher",
    ]
    with pytest.raises(SpecificationError, match=r"unknown test 'nope'; available: .*llc"):
        resolve_test("nope")


def test_panel_tests_lists_the_category_a_tests_by_family() -> None:
    assert panel_tests("unitroot") == [
        "fisher",
        "harris_tzavalis",
        "im_pesaran_shin",
        "levin_lin_chu",
    ]
    assert panel_tests("serial") == ["wooldridge"]
    assert panel_tests("heterogeneity") == ["slope_homogeneity"]
    assert "hansen" not in panel_tests()


def test_no_default_test(panel: xt.PanelData) -> None:
    with pytest.raises(SpecificationError, match="choose a test with test=.*fisher, harris"):
        xt.unit_root(panel, "y")
    with pytest.raises(
        SpecificationError,
        match="for cross-section dependence: frees, friedman, lm_cross_section, pesaran_cd",
    ):
        panel.cross_section_dependence("y")


def test_a_test_of_another_family_or_on_arrays_is_refused(panel: xt.PanelData) -> None:
    with pytest.raises(SpecificationError, match="'pesaran_cd' is not a unit-root test"):
        xt.unit_root(panel, "y", test="cd")
    with pytest.raises(SpecificationError, match="test on arrays \(category B\)"):
        panel.test("hansen")


def test_unknown_options_and_missing_arguments_name_what_the_test_accepts(
    panel: xt.PanelData,
) -> None:
    with pytest.raises(
        SpecificationError, match="harris_tzavalis: got an unexpected keyword argument 'lags'"
    ) as info:
        xt.unit_root(panel, "y", test="ht", lags=1)
    assert "deterministics='constant'" in str(info.value) and "small_sample_t" in str(info.value)
    with pytest.raises(SpecificationError, match="fisher: missing a required keyword-only"):
        xt.unit_root(panel, "y", test="fisher")
    with pytest.raises(SpecificationError, match="wooldridge: missing a required argument"):
        panel.test("wooldridge", "y")


def test_a_frame_is_refused_with_the_way_to_declare_it(panel: xt.PanelData) -> None:
    with pytest.raises(SpecificationError, match="needs a declared panel.*got DataFrame"):
        xt.unit_root(panel.frame, "y", test="ips")  # type: ignore[arg-type]


def test_the_refusals_are_translated(panel: xt.PanelData) -> None:
    with _config.option_context(lang="fr"):
        with pytest.raises(SpecificationError, match="choisissez un test avec test="):
            xt.unit_root(panel, "y")
        with pytest.raises(SpecificationError, match="n'est pas un test de la famille"):
            xt.unit_root(panel, "y", test="cd")
        with pytest.raises(SpecificationError, match="test inconnu 'zzz'"):
            run_test("zzz")
