"""The panel audit (K.1): each alert on a panel built to raise it, exports, invariances."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from validation import panels

import pyxtpanel as xt
from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.explore.checks._spec import CHECKS

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")


def make(units: dict[int, list[int]], seed: int = 0, **columns: object) -> xt.PanelData:
    """A panel with the given periods per unit and a normal variable ``x`` (plus ``columns``)."""
    rng = np.random.default_rng(seed)
    rows = [(u, t) for u, ts in units.items() for t in ts]
    frame = pd.DataFrame(rows, columns=["id", "t"])
    frame["x"] = rng.normal(size=len(frame))
    for name, build in columns.items():
        frame[name] = build(frame, rng) if callable(build) else build  # type: ignore[operator]
    return xt.xtset(frame, "id", "t")


def balanced(n_units: int = 40, n_periods: int = 8, **columns: object) -> xt.PanelData:
    return make({u: list(range(n_periods)) for u in range(n_units)}, **columns)


def alerts(report: xt.audit.__class__ | object) -> dict[str, str]:  # noqa: ANN401
    return {a.code: a.level for a in report.alerts}  # type: ignore[attr-defined]


# ----------------------------------------------------------------- a clean panel


def test_a_clean_panel_raises_no_warning() -> None:
    panel = xt.xtset(panels.load("balanced"), "id", "t")
    report = panel.audit()
    assert report.n_critical == 0 and report.n_warning == 0
    assert "Nothing to report" in report.summary("en") or report.n_info > 0
    assert report.facts.loc["n_units", "value"] == 40 and report.facts.loc["gamma", "value"] == 1.0


# ----------------------------------------------------------------- structure


def test_size_alerts() -> None:
    assert alerts(make({1: list(range(8))}).audit()).get("single_unit") == "critical"
    assert alerts(balanced(n_units=5).audit()).get("few_units") == "warning"
    assert alerts(balanced(n_periods=2).audit()).get("very_short") == "critical"
    short = alerts(balanced(n_periods=3).audit())
    assert short.get("short_panel") == "warning" and "very_short" not in short
    units = {u: list(range(8)) for u in range(40)}
    units[40], units[41] = [3], [5]
    assert alerts(make(units).audit()).get("single_obs_units") == "info"
    many = {u: [0] for u in range(10)} | {u + 10: list(range(8)) for u in range(30)}
    assert alerts(make(many).audit()).get("single_obs_units") == "warning"


def test_thresholds_can_be_overridden_and_unknown_ones_are_refused() -> None:
    panel = balanced(n_units=40)
    assert "few_units" not in alerts(panel.audit())
    assert alerts(panel.audit(thresholds={"min_units": 100})).get("few_units") == "warning"
    with pytest.raises(SpecificationError) as bad:
        panel.audit(thresholds={"min_unit": 5})
    assert bad.value.stata_code == 198 and "min_units" in str(bad.value)


# ----------------------------------------------------------------- balance, gaps, attrition


def uneven() -> xt.PanelData:
    units = {}
    for u in range(40):
        periods = list(range(10))
        if u < 8:
            periods.remove(4)  # an internal gap: 20 % of the units
        elif u < 20:
            periods = periods[:2]  # leaves early: 30 %
        elif u < 32:
            periods = periods[6:]  # enters late: 30 %
        units[u] = periods
    return make(units)


def test_balance_alerts() -> None:
    found = alerts(uneven().audit())
    assert found["internal_gaps"] == "warning"
    assert found["heavy_attrition"] == "warning"
    assert found["late_entry"] == "info"
    assert found["unbalanced"] == "warning"


def test_a_few_gaps_are_only_information() -> None:
    units = {u: list(range(8)) for u in range(40)}
    units[0].remove(3)
    found = alerts(make(units).audit())
    assert found["internal_gaps"] == "info" and "unbalanced" not in found


# ----------------------------------------------------------------- missing values


def test_missing_value_alerts() -> None:
    def gen(share: float):  # type: ignore[no-untyped-def]
        return lambda frame, rng: np.where(
            rng.random(len(frame)) < share, np.nan, rng.normal(size=len(frame))
        )

    panel = balanced(none=gen(1.0), most=gen(0.7), some=gen(0.3))
    found = alerts(panel.audit())
    assert found["all_missing"] == "critical"
    assert found["mostly_missing"] == "critical"
    assert found["high_missing"] == "warning"
    assert found["listwise_loss"] == "warning"
    text = panel.audit().summary("en")
    assert "none" in text and "most" in text and "some" in text
    # an all-missing variable is not also reported as constant
    assert "constant_variable" not in found


# ----------------------------------------------------------------- variation


def test_variation_alerts() -> None:
    panel = balanced(
        const=lambda f, r: 1.0,
        fixed=lambda f, r: f["id"].astype(float),
        common=lambda f, r: f["t"].astype(float),
        sticky=lambda f, r: f["id"] * 1000.0 + r.normal(size=len(f)) * 0.01,
    )
    report = panel.audit()
    found = alerts(report)
    assert found["constant_variable"] == "critical"
    assert found["time_invariant"] == "warning"
    assert found["common_variable"] == "info"
    assert found["little_within"] == "warning"
    by_name = {a.code: a.fields["names"] for a in report.alerts}
    assert by_name["constant_variable"] == "const"
    assert by_name["time_invariant"] == "fixed" and by_name["common_variable"] == "common"
    assert by_name["little_within"] == "sticky"
    assert report.variables.loc["const", "kind"] == "constant"
    assert report.variables.loc["x", "kind"] == "varying"


# ----------------------------------------------------------------- outliers


def test_outlier_alert_counts_the_flagged_values() -> None:
    def spikes(every: int):  # type: ignore[no-untyped-def]
        def build(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
            values = rng.normal(size=len(frame))
            values[::every] = 50.0
            return values

        return build

    many = balanced(spike=spikes(8)).audit()  # one spike in every unit: 12.5 % of the values
    alert = next(a for a in many.alerts if a.code == "outliers")
    assert alert.level == "warning" and "spike" in str(alert.fields["names"])
    assert many.variables.loc["spike", "outliers_flagged"] >= 40
    few = balanced(spike=spikes(40)).audit()  # 8 spikes in 320 values: below the 5 % limit
    assert next(a for a in few.alerts if a.code == "outliers").level == "info"


# ----------------------------------------------------------------- equations


@quiet
def test_equation_alerts_need_a_dependent_variable_and_compare_fd_and_fod() -> None:
    units = {u: [t for t in range(10) if t % 3 != 1] for u in range(40)}  # a gap every 3 periods
    panel = make(units)
    assert not {"no_equations", "few_equations", "prefer_fod"} & set(alerts(panel.audit()))
    found = alerts(panel.audit(dependent="x", lags=1))
    assert found["few_equations"] == "warning" and found["prefer_fod"] == "info"
    none = make({u: [0, 2, 4, 6] for u in range(40)}).audit(dependent="x", lags=1)
    assert alerts(none)["no_equations"] == "critical"


# ----------------------------------------------------------------- raw data and arguments


def test_raw_data_with_repeated_pairs_gets_a_critical_alert() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1, 2, 2], "t": [1, 1, 2, 1, 1], "x": [1.0, 2, 3, 4, 5]})
    report = xt.audit(frame, "id", "t")
    assert report.codes() == ["duplicates"] and report.n_critical == 1
    assert report.alerts[0].fields == {"rows": 4, "pairs": 2}  # the rows of the repeated pairs
    assert "resolve_duplicates" in report.alerts[0].advice("fr")
    assert "duplicates" not in xt.audit(frame.drop(index=[1, 4]), "id", "t").codes()


def test_a_frame_is_declared_and_a_panel_is_audited_as_it_is() -> None:
    panel = balanced()
    frame = panel.to_frame()
    assert xt.audit(frame, "id", "t").codes() == panel.audit().codes()
    assert xt.audit(panel).codes() == panel.audit().codes()


def test_checks_can_be_chosen_and_unknown_ones_are_refused() -> None:
    panel = balanced(n_units=5)
    assert alerts(panel.audit(checks=["missing"])) == {}
    assert "few_units" in alerts(panel.audit(checks=["structure"]))
    with pytest.raises(SpecificationError) as bad:
        panel.audit(checks=["nothing"])
    assert bad.value.stata_code == 198


def test_alerts_are_sorted_by_level() -> None:
    report = uneven().audit()
    order = [a.level for a in report.alerts]
    assert order == sorted(order, key=("critical", "warning", "info").index)


def test_the_audit_does_not_modify_the_data() -> None:
    panel = uneven()
    before = panel.frame.copy()
    panel.audit(dependent="x")
    pd.testing.assert_frame_equal(panel.frame, before)


# ----------------------------------------------------------------- exports


def test_exports_in_both_languages() -> None:
    report = uneven().audit()
    assert "Audit du panel" in report.summary("fr") and "Panel audit" in report.summary("en")
    assert "ATTENTION" in report.summary("fr") and "WARNING" in report.summary("en")
    markdown = report.to_markdown(lang="fr")
    assert markdown.startswith("## Audit du panel") and "### Chiffres clés" in markdown
    html = report.to_html("en")
    assert html.startswith("<div class='pyxtpanel-audit'>") and "border-left:4px solid" in html
    assert report._repr_html_() == report.to_html()
    latex = report.to_latex(lang="en")
    assert latex.startswith(r"\begin{tabular}{lp{6.5cm}p{6.5cm}}") and r"\%" in latex
    frame = report.to_frame("fr")
    assert list(frame.columns) == ["level", "check", "code", "message", "advice"]
    assert len(frame) == len(report.alerts)


def test_html_escapes_the_text_of_the_messages() -> None:
    frame = pd.DataFrame(
        {"id": np.repeat(np.arange(40), 8), "t": np.tile(np.arange(8), 40), "<b>x</b>": 1.0}
    )
    html = xt.audit(frame, "id", "t").to_html("en")
    assert "<b>x</b>" not in html and "&lt;b&gt;x&lt;/b&gt;" in html


# ----------------------------------------------------------------- catalogue and registry


def test_every_code_of_every_check_has_a_message_and_an_advice_in_both_languages() -> None:
    for name, spec in CHECKS.specs():
        assert spec.codes, name
        for code in spec.codes:
            for suffix in ("", ".advice"):
                key = f"explore.audit.{code}{suffix}"
                for lang in ("en", "fr"):
                    assert _i18n._CATALOGUE[key][lang], key


def test_the_checks_run_in_the_order_they_declare() -> None:
    orders = [spec.order for _, spec in CHECKS.specs()]
    assert len(set(orders)) == len(orders)
    assert {"structure", "balance", "missing", "variation", "outliers", "equations"} <= set(
        CHECKS.names()
    )


# ----------------------------------------------------------------- invariances


@quiet
@settings(max_examples=40, deadline=None)
@given(
    units=st.dictionaries(
        st.integers(0, 8),
        st.sets(st.integers(0, 9), min_size=1, max_size=10).map(sorted),
        min_size=1,
        max_size=9,
    ),
    shift=st.integers(-100, 100),
    seed=st.integers(0, 1000),
)
def test_the_alerts_do_not_depend_on_row_order_unit_names_or_the_origin_of_time(  # type: ignore[no-untyped-def]
    units, shift, seed
) -> None:
    a = make(units, seed=1)
    frame = a.to_frame().copy()
    frame["id"] = frame["id"].map(lambda u: f"unit_{u:02d}")
    frame["t"] = frame["t"] + shift
    frame = frame.iloc[np.random.default_rng(seed).permutation(len(frame))]
    b = xt.xtset(frame, "id", "t")
    assert a.audit(dependent="x").codes() == b.audit(dependent="x").codes()
    pd.testing.assert_series_equal(
        a.audit().facts["value"], b.audit().facts["value"], check_names=False
    )
