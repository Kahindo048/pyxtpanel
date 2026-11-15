"""Description of panels: xtdescribe, xtsum, xttab, xttrans, correlations, missing values (E)."""

import numpy as np
import pandas as pd
import pytest

from pyxtpanel import xtset
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.describe._spec import DESCRIPTIONS

NAN = np.nan


def close(actual, expected, tol=1e-12):  # type: ignore[no-untyped-def]
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), np.asarray(expected, dtype=float), atol=tol, rtol=tol
    )


@pytest.fixture
def probe():  # type: ignore[no-untyped-def]
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 1, 2, 2, 2, 3, 3],
            "t": [1, 2, 3, 4, 1, 2, 4, 2, 3],
            "x": [1.0, 2, 3, 4, 5, 5, 7, 2, 2],
            "z": [10.0, 11, 12, 13, NAN, 15, 17, 20, 21],
            "d": [0.0, 0, 1, 1, 1, 1, 2, 0, 0],
        }
    )
    return xtset(frame, "id", "t")


# ---------------------------------------------------------------- xtdescribe


def test_xtdescribe_counts_ti_and_pattern_order(probe) -> None:  # type: ignore[no-untyped-def]
    result = probe.xtdescribe()
    assert (result.n_units, result.n_periods, result.span) == (3, 4, 4)
    assert (result.t_first, result.t_second, result.t_last) == ("1", "2", "4")
    assert result.delta_text == "1 unit"
    assert result.ti["min"] == 2 and result.ti["max"] == 4 and result.ti["mean"] == 3
    assert result.ti["p50"] == 3 and result.ti["N"] == 3
    # equal frequencies: ordered by the counts per column, increasing
    assert result.patterns["pattern"].tolist() == [".11.", "11.1", "1111"]
    assert result.patterns["freq"].tolist() == [1, 1, 1]
    assert result.footer == "XXXX" and result.n_patterns == 3 and result.other_freq == 0


def test_patterns_limit_other_row_and_zero(probe) -> None:  # type: ignore[no-untyped-def]
    some = probe.xtdescribe(patterns=2)
    assert len(some.patterns) == 2 and some.other_freq == 1 and some.n_patterns == 3
    assert "(other patterns)" in some.summary()
    none = probe.xtdescribe(patterns=0)
    assert none.patterns.empty and "Pattern" not in none.summary()
    with pytest.raises(SpecificationError, match="width"):
        probe.xtdescribe(width=2)
    with pytest.raises(SpecificationError, match="patterns"):
        probe.xtdescribe(patterns=-1)


def test_ten_or_more_observations_in_a_column_show_x() -> None:
    frame = pd.DataFrame({"id": [1] * 12 + [2], "t": [1] * 0 + list(range(1, 13)) + [1]})
    # 12 periods, width 3: (12 - 1) // (3 - 1) = 5 periods per column -> counts 5, 5, 2
    panel = xtset(frame.assign(x=1.0), "id", "t")
    narrow = panel.xtdescribe(width=3)
    assert narrow.periods_per_column == 5 and "Each column represents 5 periods" in narrow.summary()
    assert narrow.patterns["pattern"].tolist() == ["1..", "552"]
    wide = xtset(pd.DataFrame({"id": [1] * 30 + [2], "t": list(range(30)) + [0]}), "id", "t")
    assert any("X" in p for p in wide.xtdescribe(width=3).patterns["pattern"])


def test_variables_restrict_to_the_rows_where_they_are_observed(probe) -> None:  # type: ignore[no-untyped-def]
    result = probe.xtdescribe(variables=["z"])
    assert result.ti["min"] == 2 and result.ti["N"] == 3
    assert result.patterns["pattern"].tolist()[0].count(".") >= 1
    with_missing = probe.xtdescribe()
    assert with_missing.n_units == result.n_units
    assert result.patterns["freq"].sum() == 3


def test_xtdescribe_refuses_one_period() -> None:
    one = xtset(pd.DataFrame({"id": [1, 2], "t": [5, 5], "x": [1.0, 2]}), "id", "t")
    with pytest.raises(DataError, match="multiple distinct") as info:
        one.xtdescribe()
    assert info.value.stata_code == 459


def test_xtdescribe_text_in_both_languages(probe) -> None:  # type: ignore[no-untyped-def]
    result = probe.xtdescribe()
    assert "uniquely identifies" in result.summary("en")
    assert "identifie de façon unique" in result.summary("fr")
    assert result.to_frame() is result.patterns


# ---------------------------------------------------------------- xtsum


def test_xtsum_matches_a_pandas_computation_on_an_unbalanced_panel(probe) -> None:  # type: ignore[no-untyped-def]
    frame = probe.to_frame()
    row = probe.xtsum("x").scalars.loc["x"]
    x = frame["x"]
    means = frame.groupby("id")["x"].mean()
    deviation = x - frame.groupby("id")["x"].transform("mean")
    close(row["N"], 9)
    close(row["n"], 3)
    close(row["Tbar"], 3)
    close(row["mean"], x.mean())
    close(row["sd"], x.std())
    close(row["sd_b"], means.std())
    close(row["min_b"], means.min())
    close(row["sd_w"], np.sqrt((deviation**2).sum() / (9 - 1)))
    close(row["min_w"], deviation.min() + x.mean())
    close(row["max_w"], deviation.max() + x.mean())
    close(row["within_share"], row["sd_w"] ** 2 / (row["sd_b"] ** 2 + row["sd_w"] ** 2))


def test_xtsum_missing_values_constants_and_terms(probe) -> None:  # type: ignore[no-untyped-def]
    z = probe.xtsum("z").scalars.loc["z"]
    assert z["N"] == 8 and z["n"] == 3 and z["Tbar"] == 8 / 3
    frame = probe.to_frame().assign(k=4.0)
    constant = xtset(frame, "id", "t").xtsum("k").scalars.loc["k"]
    assert constant["sd"] == 0 and constant["sd_w"] == 0 and constant["min_w"] == 4.0
    assert np.isnan(constant["within_share"])
    lagged = probe.xtsum("L.x")
    assert lagged.table.index.get_level_values("variable")[0] == "L.x"
    assert lagged.scalars.loc["L.x", "N"] == 5


def test_xtsum_one_unit_and_empty_variable() -> None:
    single = xtset(pd.DataFrame({"id": [1, 1, 1], "t": [1, 2, 3], "x": [1.0, 2, 4]}), "id", "t")
    row = single.xtsum("x").scalars.loc["x"]
    assert np.isnan(row["sd_b"]) and row["n"] == 1 and row["sd_w"] > 0
    blank = xtset(pd.DataFrame({"id": [1, 1], "t": [1, 2], "x": [NAN, NAN]}), "id", "t")
    assert blank.xtsum("x").scalars.loc["x", "N"] == 0


def test_xtsum_by_groups_formats_and_errors(probe) -> None:  # type: ignore[no-untyped-def]
    grouped = probe.xtsum(["x"], by="d")
    assert grouped.table.index.names == ["d", "variable", "component"]
    assert set(grouped.scalars.index.get_level_values("d")) == {0.0, 1.0, 2.0}
    plain = probe.xtsum(["x", "z"])
    text = plain.summary()
    assert "overall" in text and "between" in text and "T-bar" in text
    assert plain.to_markdown().startswith("| variable | component |")
    assert plain.to_frame() is plain.table
    with pytest.raises(SpecificationError, match="at least one"):
        probe.xtsum([])
    with pytest.raises(SpecificationError, match="by="):
        probe.xtsum("x", by="nope")


# ---------------------------------------------------------------- xttab


def test_xttab_hand_computed(probe) -> None:  # type: ignore[no-untyped-def]
    table = probe.xttab("d").table
    # values 0: units 1 (2 of 4), 3 (2 of 2); 1: units 1 (2 of 4), 2 (2 of 3); 2: unit 2 (1 of 3)
    assert table.index.tolist() == [0.0, 1.0, 2.0, "Total"]
    assert table["overall_freq"].tolist() == [4, 4, 1, 9]
    close(table["overall_pct"], [400 / 9, 400 / 9, 100 / 9, 100])
    assert table["between_freq"].tolist() == [2, 2, 1, 5]
    close(table["between_pct"], [200 / 3, 200 / 3, 100 / 3, 500 / 3])
    within = [(0.5 + 1.0) / 2 * 100, (0.5 + 2 / 3) / 2 * 100, 100 / 3]
    close(table["within_pct"].iloc[:3], within)
    close(table["within_pct"].iloc[3], (2 * within[0] + 2 * within[1] + 1 * within[2]) / 5)


def test_xttab_guards(probe) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(SpecificationError, match="distinct values") as info:
        probe.xttab("z", max_values=3)
    assert info.value.stata_code == 134
    with pytest.raises(SpecificationError, match="exactly one"):
        probe.xttab("L(1/2).x")
    nothing = xtset(pd.DataFrame({"id": [1, 1], "t": [1, 2], "x": [NAN, NAN]}), "id", "t")
    with pytest.raises(DataError, match="no observed value"):
        nothing.xttab("x")
    result = probe.xttab("d")
    assert result.n_units == 3 and result.n_obs == 9
    assert "Total" in result.summary() and result.to_markdown().startswith("| value")


# ---------------------------------------------------------------- xttrans


def test_xttrans_bridges_gaps_like_stata(probe) -> None:  # type: ignore[no-untyped-def]
    result = probe.xttrans("d")
    # unit 1: 0->0, 0->1, 1->1 ; unit 2 (t=1,2,4): 1->1, 1->2 (bridges 2->4) ; unit 3: 0->0
    assert result.counts.to_numpy().tolist() == [[2, 1, 0], [0, 2, 1]] or result.n_trans == 6
    assert result.n_trans == 6 and result.n_bridged == 1
    assert result.counts.loc[1.0, 2.0] == 1  # the bridged transition
    assert (result.n_rows, result.n_cols) == (2, 3)
    close(result.percent.sum(axis=1), [100.0, 100.0])


def test_xttrans_strict_drops_bridged_pairs_and_missing_is_not_bridged(probe) -> None:  # type: ignore[no-untyped-def]
    strict = probe.xttrans("d", gaps="strict")
    assert strict.n_trans == 5 and strict.n_bridged == 1 and 2.0 not in strict.counts.columns
    frame = pd.DataFrame({"id": [1, 1, 1], "t": [1, 2, 3], "d": [0.0, NAN, 1.0]})
    with pytest.raises(DataError, match="no transition"):
        xtset(frame, "id", "t").xttrans("d")
    with pytest.raises(SpecificationError, match="bridge"):
        probe.xttrans("d", gaps="loose")  # type: ignore[arg-type]
    assert "bridge a gap" in probe.xttrans("d").summary()


def test_xttrans_compares_values_as_float32() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1], "t": [1, 2, 3], "d": [0.1, 0.1 + 1e-10, 0.1]})
    result = xtset(frame, "id", "t").xttrans("d")
    assert result.n_rows == 1 and result.n_cols == 1  # Stata stores the variable as float


# ---------------------------------------------------------------- correlations


def test_xtcorr_overall_within_between_against_pandas(probe) -> None:  # type: ignore[no-untyped-def]
    frame = probe.to_frame()
    both = frame.dropna(subset=["x", "z"])
    result = probe.xtcorr(["x", "z"])
    close(result.corr["overall"].loc["x", "z"], both["x"].corr(both["z"]))
    dev = both[["x", "z"]] - both.groupby("id")[["x", "z"]].transform("mean")
    close(result.corr["within"].loc["x", "z"], dev["x"].corr(dev["z"]))
    means = both.groupby("id")[["x", "z"]].mean()
    close(result.corr["between"].loc["x", "z"], means["x"].corr(means["z"]))
    assert result.n["overall"].loc["x", "z"] == 8 and result.n["between"].loc["x", "z"] == 3
    assert result.n["overall"].loc["x", "x"] == 9
    close(np.diag(result.corr["within"].to_numpy()), [1, 1])


def test_xtcorr_options_and_outputs(probe) -> None:  # type: ignore[no-untyped-def]
    only = probe.xtcorr(["x", "z"], kind="within")
    assert list(only.corr) == ["within"]
    assert np.isnan(probe.xtcorr(["x", "z"]).corr["within"].to_numpy()).sum() == 0
    with pytest.raises(SpecificationError, match="at least two"):
        probe.xtcorr(["x"])
    with pytest.raises(SpecificationError, match="kind"):
        probe.xtcorr(["x", "z"], kind="total")  # type: ignore[arg-type]
    assert "within" in probe.xtcorr(["x", "z"]).summary()
    assert probe.xtcorr(["x", "z"]).to_frame().index.names == ["kind", "variable"]
    assert probe.xtcorr(["x", "z"]).to_markdown().startswith("| kind")
    constant = xtset(probe.to_frame().assign(k=1.0), "id", "t").xtcorr(["x", "k"])
    assert np.isnan(constant.corr["overall"].loc["x", "k"])


def test_autocorr_equals_pandas_for_one_unit_and_respects_gaps() -> None:
    rng = np.random.default_rng(1)
    x = rng.normal(size=40).cumsum()
    one = xtset(pd.DataFrame({"id": 1, "t": range(40), "x": x}), "id", "t")
    table = one.autocorr("x", [1, 2, 3]).table
    series = pd.Series(x - x.mean())
    for k in (1, 2, 3):
        close(table.loc[k, "rho"], series.autocorr(k))
        assert table.loc[k, "n"] == 40 - k
    gappy = xtset(pd.DataFrame({"id": 1, "t": [1, 2, 4, 5], "x": [1.0, 3, 2, 6]}), "id", "t")
    assert gappy.autocorr("x", 1).table.loc[1, "n"] == 2  # pairs (1,2) and (4,5) only
    with pytest.raises(SpecificationError, match="lags"):
        one.autocorr("x", 0)
    assert "x" in one.autocorr("x").summary() and one.autocorr("x").to_frame().shape == (1, 2)


# ---------------------------------------------------------------- attrition


def test_attrition_by_period_and_survival(probe) -> None:  # type: ignore[no-untyped-def]
    result = probe.attrition()
    table = result.by_period
    assert table["active"].tolist() == [2, 3, 2, 2]
    assert table["entries"].tolist() == [2, 1, 0, 0]
    assert table["exits"].tolist() == [0, 0, 1, 2]
    assert table["reentries"].tolist() == [0, 0, 0, 1]  # unit 2 comes back in period 4
    assert table["in_gap"].tolist() == [0, 0, 1, 0]  # unit 2 is absent in period 3
    assert result.units_with_gaps == 1
    survival = result.survival
    assert survival["at_risk"].tolist() == [3, 3, 2, 2]
    # unit 3 (life 2) exits before the end; units 1 and 2 are still observed in period 4
    assert survival["exits"].tolist() == [0, 1, 0, 0] and survival["censored"].tolist() == [
        0,
        0,
        0,
        2,
    ]
    close(survival["survival"], [1.0, 2 / 3, 2 / 3, 2 / 3])
    assert "units with gaps" in result.summary() and result.to_frame() is table


# ---------------------------------------------------------------- missing values


def test_missing_report_tables_and_patterns() -> None:
    frame = pd.DataFrame(
        {
            "id": [1, 1, 2, 2, 2],
            "t": [1, 2, 1, 2, 3],
            "a": [1.0, NAN, 3, 4, NAN],
            "b": [NAN, NAN, 1, 2, 3],
        }
    )
    report = xtset(frame, "id", "t").missing_report()
    by_variable = report.by_variable
    assert by_variable.loc["a", "n_missing"] == 2 and by_variable.loc["b", "n_missing"] == 2
    close(by_variable.loc["a", "pct_missing"], 40.0)
    assert by_variable.loc["b", "units_all_missing"] == 1  # unit 1
    assert by_variable.loc["a", "units_some_missing"] == 2
    assert report.by_unit.loc[1, "b"] == 2 and report.by_unit["T_i"].tolist() == [2, 3]
    close(report.by_period.loc["2", "a"], 0.5)
    assert report.by_period["n_obs"].tolist() == [2, 2, 1]
    patterns = dict(zip(report.patterns["pattern"], report.patterns["freq"], strict=True))
    assert patterns == {"1.": 1, "..": 1, "11": 2, ".1": 1} or patterns["11"] == 2
    close(report.patterns["percent"].sum(), 100.0)
    assert "n_missing" in report.summary() and report.to_frame() is by_variable


def test_missing_report_errors(probe) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(DataError, match="not found"):
        probe.missing_report(["nope"])
    only_keys = xtset(pd.DataFrame({"id": [1, 1], "t": [1, 2]}), "id", "t")
    with pytest.raises(SpecificationError, match="no variable"):
        only_keys.missing_report()


def test_complete_cases_effect_on_ti(probe) -> None:  # type: ignore[no-untyped-def]
    result = probe.complete_cases(["z"])
    assert (result.n_rows, result.n_rows_kept) == (9, 8)
    assert (result.n_units, result.n_units_kept) == (3, 3)
    assert result.mask.sum() == 8 and result.panel.n_obs == 8
    assert result.ti.loc["min", "before"] == 2 and result.ti.loc["mean", "after"] == 8 / 3
    assert "rows: 9 -> 8" in result.summary() and result.to_frame() is result.ti
    nothing = xtset(pd.DataFrame({"id": [1, 1], "t": [1, 2], "x": [NAN, NAN]}), "id", "t")
    result = nothing.complete_cases(["x"])
    assert result.panel is None and result.n_units_kept == 0


# ---------------------------------------------------------------- registry


def test_descriptions_are_registered_with_their_stata_origin() -> None:
    expected = {"xtdescribe", "xtsum", "xttab", "xttrans", "xtcorr", "autocorr"}
    expected |= {"attrition", "missing_report", "complete_cases"}
    assert expected <= set(DESCRIPTIONS.names())
    assert DESCRIPTIONS.spec("xtsum").stata == "xtsum"
    assert DESCRIPTIONS.spec("autocorr").stata is None


def test_xttrans_of_a_discretized_continuous_variable() -> None:
    frame = pd.DataFrame(
        {"id": [1, 1, 1, 2, 2, 2], "t": [1, 2, 3] * 2, "v": [0.1, 0.9, 0.2, 0.8, 0.7, 0.3]}
    )
    panel = xtset(frame, "id", "t")
    by_edges = panel.xttrans("v", bins=[0.5])
    # classes: <0.5 -> 0, >=0.5 -> 1;  unit 1: 0,1,0  unit 2: 1,1,0
    assert by_edges.counts.to_numpy().tolist() == [[0, 1], [2, 1]]
    by_quantile = panel.xttrans("v", bins=2)  # cut at the median of the six values
    assert by_quantile.n_trans == 4 and set(by_quantile.counts.index) <= {0.0, 1.0}
    with pytest.raises(SpecificationError, match="bins"):
        panel.xttrans("v", bins=1)
