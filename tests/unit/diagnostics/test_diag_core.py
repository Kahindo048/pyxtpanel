"""Infrastructure of the statistical tests: result, registry and normalised inputs (L.1)."""

from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel._registry import Registry
from pyxtpanel.diagnostics import TestResult, TestSpec, list_tests, register_test, run_test
from pyxtpanel.diagnostics._core import as_array, panel_variables, spec

_i18n.register(
    "diag",
    {
        "diag.probe.title": {"en": "Probe test", "fr": "Test sonde"},
        "diag.probe.null": {"en": "nothing happens", "fr": "il ne se passe rien"},
        "diag.probe.alternative": {"en": "something happens", "fr": "il se passe quelque chose"},
    },
)


def result(**changes: object) -> TestResult:
    base: dict[str, object] = {
        "name": "probe",
        "statistic": 2.5,
        "pvalue": 0.0124,
        "statistic_name": "z",
        "n_obs": 100,
        "n_groups": 10,
        "variable": "x",
    }
    return TestResult(**{**base, **changes})  # type: ignore[arg-type]


# ---------------------------------------------------------------- TestResult


def test_result_to_frame_has_one_row_with_the_degrees_of_freedom() -> None:
    frame = result(df=(3,)).to_frame()
    assert list(frame.columns) == ["statistic", "pvalue", "df", "n_obs", "n_groups"]
    assert frame.index.tolist() == ["probe"]
    two = result(df=(2, 40)).to_frame()
    assert {"df1", "df2"} <= set(two.columns)
    assert two.loc["probe", "df2"] == 40.0


def test_result_summary_is_bilingual_and_carries_the_numbers() -> None:
    en, fr = result(df=(3,)).summary("en"), result(df=(3,)).summary("fr")
    assert "Probe test (x)" in en and "Test sonde (x)" in fr
    assert "Ho: nothing happens" in en and "H0: il ne se passe rien" in fr
    assert "z = 2.5000" in en and "p-value = 0.0124" in en
    assert "100 observations, 10 units" in en and "100 observations, 10 unités" in fr


def test_result_summary_lists_options_notes_and_details() -> None:
    details = pd.DataFrame({"t": [1.5, 2.0]}, index=pd.Index(["a", "b"], name="unit"))
    text = result(options={"lags": 2}, notes=("one unit dropped",), details=details).summary("en")
    assert "Options: lags=2" in text and "Notes: one unit dropped" in text
    assert "1.5000" in text
    assert result(details=details)._repr_html_().count("<table") == 2


def test_result_is_immutable_and_options_are_read_only() -> None:
    res = result(options={"lags": 1})
    with pytest.raises(dataclasses.FrozenInstanceError):
        res.statistic = 0.0  # type: ignore[misc]
    with pytest.raises(TypeError):
        res.options["lags"] = 3  # type: ignore[index]


def test_result_rejects_at_a_level_and_a_nan_pvalue_never_rejects() -> None:
    assert result(pvalue=0.04).rejects(0.05) and not result(pvalue=0.04).rejects(0.01)
    assert not result(pvalue=np.nan).rejects(0.05)
    assert result(statistic=3, pvalue=0.5).statistic == 3.0  # coerced to float


def test_result_is_not_collected_as_a_test_class() -> None:
    assert TestResult.__test__ is False and TestSpec.__test__ is False


# ---------------------------------------------------------------- registry


SPEC = TestSpec(
    family="serial",
    category="B",
    wave=1,
    reference="Author (2000)",
    stata="mata",
    validation=("probe",),
    inputs="arrays",
    summary="a probe",
)


def test_register_list_and_run_use_the_registry(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(spec, "TESTS", Registry("test"))

    @register_test("probe", SPEC)
    def probe(x: object) -> TestResult:
        return result(statistic=float(np.sum(np.asarray(x))))

    table = list_tests()
    assert table.index.tolist() == ["probe"]
    assert table.loc["probe", "category"] == "B" and table.loc["probe", "stata"] == "mata"
    assert run_test("probe", [1, 2]).statistic == 3.0
    with pytest.raises(SpecificationError, match="available: probe"):
        run_test("nothing")


def test_an_empty_registry_still_lists_its_columns(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(spec, "TESTS", Registry("test"))
    table = list_tests()
    assert len(table) == 0 and {"family", "category", "wave", "stata"} <= set(table.columns)


def test_the_root_package_exports_the_test_api() -> None:
    assert xt.TestResult is TestResult
    assert xt.list_tests is list_tests and xt.run_test is run_test


# ---------------------------------------------------------------- panel inputs


def panel(rows: list[tuple[int, int, float]]) -> xt.PanelData:
    frame = pd.DataFrame(rows, columns=["id", "t", "x"])
    return xt.xtset(frame, "id", "t")


BALANCED = [(i, t, float(10 * i + t)) for i in (1, 2, 3) for t in range(1, 5)]


def test_a_balanced_panel_passes_every_requirement() -> None:
    arrays = panel_variables(
        panel(BALANCED), ["x"], test="T", balanced="required", consecutive=True, min_obs=4
    )
    assert arrays.n_units == 3 and arrays.counts.tolist() == [4, 4, 4]
    assert arrays.values[:, 0].tolist() == [11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34]


def test_an_unbalanced_panel_is_refused_with_the_first_incomplete_unit() -> None:
    rows = [r for r in BALANCED if r != (2, 3, 23.0)]
    with pytest.raises(DataError, match=r"unit 2 is observed 3 times out of 4") as info:
        panel_variables(panel(rows), ["x"], test="llc", balanced="required")
    assert info.value.stata_code == 498 and "llc needs" in str(info.value)
    assert panel_variables(panel(rows), ["x"], test="ips", balanced="both").n_obs == 11


def test_a_gap_is_refused_with_the_unit_and_the_period() -> None:
    rows = [r for r in BALANCED if r != (2, 3, 23.0)]
    with pytest.raises(DataError, match=r"unit 2 has a gap after period 2") as info:
        panel_variables(panel(rows), ["x"], test="xtserial", consecutive=True)
    assert info.value.stata_code == 498


def test_a_missing_value_removes_its_row_and_can_make_a_gap() -> None:
    rows = [(i, t, np.nan if (i, t) == (1, 2) else x) for i, t, x in BALANCED]
    arrays = panel_variables(panel(rows), ["x"], test="T")
    assert arrays.n_obs == 11 and arrays.counts.tolist() == [3, 4, 4]
    assert arrays.offsets.tolist() == [0, 3, 7, 11] and arrays.valid.all()
    assert arrays.unit_of_row.tolist() == [0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]
    with pytest.raises(DataError, match="unit 1 has a gap after period 1"):
        panel_variables(panel(rows), ["x"], test="T", consecutive=True)


def test_a_unit_entirely_missing_disappears() -> None:
    rows = [(i, t, np.nan if i == 2 else x) for i, t, x in BALANCED]
    arrays = panel_variables(panel(rows), ["x"], test="T", balanced="required")
    assert arrays.unit_labels.tolist() == [1, 3] and arrays.n_units == 2
    assert arrays.unit_of_row.tolist() == [0] * 4 + [1] * 4


def test_short_units_and_too_few_units_are_refused() -> None:
    with pytest.raises(DataError, match="at least 5 observations per unit: unit 1 has 4"):
        panel_variables(panel(BALANCED), ["x"], test="adf", min_obs=5)
    with pytest.raises(DataError, match="at least 4 units, got 3"):
        panel_variables(panel(BALANCED), ["x"], test="cd", min_units=4)


def test_no_observation_left_is_refused() -> None:
    rows = [(i, t, np.nan) for i, t, _ in BALANCED]
    with pytest.raises(DataError, match="no observation is left"):
        panel_variables(panel(rows), ["x"], test="T")


def test_the_balance_is_judged_on_the_periods_that_remain() -> None:
    # Period 4 is missing for every unit: still strongly balanced on periods 1-3.
    rows = [(i, t, x) for i, t, x in BALANCED if t != 4]
    assert panel_variables(panel(rows), ["x"], test="T", balanced="required").n_obs == 9


def test_two_variables_drop_the_row_when_either_is_missing() -> None:
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 2, 2, 2],
            "t": [1, 2, 3] * 2,
            "x": [1, 2, 3, 4, 5, 6.0],
            "z": [1, np.nan, 3, 4, 5, 6.0],
        }
    )
    arrays = panel_variables(xt.xtset(frame, "id", "t"), ["x", "z"], test="T")
    assert arrays.n_obs == 5 and arrays.values.shape == (5, 2)


# ---------------------------------------------------------------- array inputs


def test_as_array_accepts_lists_series_and_checks_shape() -> None:
    assert as_array("b", [1, 2, 3], ndim=1).dtype == np.float64
    assert as_array("V", pd.DataFrame(np.eye(2)), ndim=2, shape=(2, 2)).shape == (2, 2)
    with pytest.raises(
        SpecificationError, match=r"V must be 2 x 2, got an array of shape \(3, 2\)"
    ):
        as_array("V", np.ones((3, 2)), ndim=2, shape=(2, 2))
    with pytest.raises(SpecificationError, match="b must be n, got an array of shape"):
        as_array("b", np.ones((2, 2)), ndim=1)
    with pytest.raises(SpecificationError, match="must be numeric"):
        as_array("b", ["a", "b"], ndim=1)


def test_as_array_refuses_missing_values_unless_allowed() -> None:
    with pytest.raises(DataError, match="b contains 1 missing or infinite"):
        as_array("b", [1.0, np.nan], ndim=1)
    assert np.isnan(as_array("b", [1.0, np.nan], ndim=1, finite=False)[1])
    with pytest.raises(DataError, match="2 missing or infinite"):
        as_array("b", [np.inf, -np.inf], ndim=1)


def test_the_category_b_tests_are_importable_from_the_diagnostics_package() -> None:
    # The estimation packages import them from here: this path is part of the contract.
    from pyxtpanel import diagnostics

    for name in (
        "hansen",
        "sargan",
        "hausman",
        "wald",
        "ar_test",
        "moment_covariance",
        "Equations",
    ):
        assert hasattr(diagnostics, name) and name in diagnostics.__all__
    for name in ("hansen", "sargan", "hausman", "wald", "ar_test", "diff_hansen"):
        assert getattr(diagnostics, name) is spec.TESTS.get(name)
