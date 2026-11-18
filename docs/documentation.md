# pyxtpanel documentation

`pyxtpanel` is the Python foundation for panel data: it declares a panel the way Stata's `xtset`
does, manages and explores it, and runs the statistical tests of panel data that do not depend on an
estimator. It estimates no model. Everything that has a Stata 17 equivalent is compared with Stata
on shared datasets; the frozen reference results run in the test suite without Stata
(annex E).

This is the single reference document. Parts 0 to 9 are written by hand and **every Python block in
parts 0 to 8 is executed by the test suite**, so the examples stay true. The annexes A to E are
generated from the registries of the package (`python -m tools.doc --write`) and cannot drift from the
code.

| Part | Content |
|---|---|
| [0](#0-getting-started) | Install, a first panel in ten lines, map of the package |
| [1](#1-concepts) | `PanelData`, the calendar, balance, missing values, results, errors, language |
| [2](#2-declaring-and-reading-a-panel) | `xtset`, reading and writing `.dta`, selecting, filling, balancing |
| [3](#3-operators-and-transformations) | `L F D S`, within, between, forward orthogonal deviations, windows |
| [4](#4-describing-a-panel) | `xtdescribe`, `xtsum`, `xttab`, `xttrans`, correlations, missing values |
| [5](#5-exploration-tools) | Gaps, balance, variation, spells, events, outliers, conversions, audit |
| [6](#6-figures) | The optional `plot` extra |
| [7](#7-statistical-tests-guide) | The tests of panel data: which one, how to call it, how to read it |
| [8](#8-authors-of-estimation-packages) | The array contract and what an estimator may rely on |
| [9](#9-method-and-validation) | Formulas, Stata conventions reproduced, validation |
| Annexes | A functions, B test cards, C errors, D Stata behaviours reproduced, E validation journal |

## 0. Getting started

Install the package (the figures need the optional extra):

```text
pip install pyxtpanel            # numpy, pandas, scipy
pip install "pyxtpanel[plot]"    # adds matplotlib for PanelData.plot
```

The examples of this document use one simulated country-by-year panel. `full` is complete (24
countries, 25 years); `df` has 5 % of its rows removed at random, so it has holes, as real data do.
`reform` is a 0/1 policy state that switches on in some countries and never in others.

```python
import numpy as np
import pandas as pd
import pyxtpanel as xt

rng = np.random.default_rng(2026)
rows = []
for i in range(24):
    country, level = f"C{i:02d}", rng.normal(8, 1)
    reform_year = int(rng.choice([1998, 2002, 2006, 3000]))  # 3000: never reformed
    for year in range(1990, 2015):
        level += 0.02 + 0.03 * rng.normal()
        rows.append((country, year, level, float(year >= reform_year)))
full = pd.DataFrame(rows, columns=["country", "year", "lgdp", "reform"])
full["invest"] = 0.3 * full["lgdp"] + rng.normal(size=len(full))
df = full[rng.random(len(full)) > 0.05].reset_index(drop=True)
```

Declare the panel once, then read it:

```python
panel = xt.xtset(df, "country", "year")  # like `xtset country year`
print(panel)  # structure, as xtset prints it
print(panel.xtsum(["lgdp", "invest"]).summary())  # overall, between, within
growth = panel.operate("D.lgdp")  # first difference on the calendar
print(growth.head())
```

### Map of the package

`pyxtpanel` is built in layers; a layer imports only from the ones below it, and a test
(import-linter) enforces it. From the top:

| Layer | Folders | Role |
|---|---|---|
| 1 | `api/`, `contract` | `PanelData`, the top-level functions, the array contract |
| 2 | `plot/`, `report/`, `datasets/` | figures (the only place that imports matplotlib), exports, bundled data |
| 3 | `explore/` | exploration tools (never import the statistical tests) |
| 4 | `describe/`, `diagnostics/` | `xtdescribe`, `xtsum`, ...; the statistical tests |
| 5 | `ops/` | operators and transformations, a registry of named operators |
| 6 | `panel/` | declaration, sorted index, balance, gaps, `tsfill`, `.dta` input and output, arrays |
| 7 | `time/` | the calendar axis |
| 8 | `_results`, `_registry`, `_i18n`, `_config`, `_errors`, `_types` | result base classes, registries, French and English messages, options, typed errors |

## 1. Concepts

**`PanelData`.** `xt.xtset(frame, id, time)` validates a long DataFrame and returns an immutable
`PanelData`. It never modifies the user's frame; every function reads the panel and returns a new
object. The panel keeps the frame and a sorted index (unit, then period) that every operation reuses.
Most functions exist both as methods (`panel.xtsum(...)`) and as top-level functions
(`xt.xtsum(panel, ...)`); the two give the same result. The tools that work on a plain DataFrame
(`to_long`, `resolve_duplicates`) are top-level functions only.

**The calendar.** The time variable is an integer (with `delta` for a step other than 1), a date
(`freq`) or a Stata-coded integer (`format`: `%td`, `%tw`, `%tm`, `%tq`, `%th` or `%ty`). The panel
works on the *calendar*, not on rows: the lag of 2005 in a unit that has no 2004 row is missing, not
the 2003 value. Periods are compared by their position `(t - t0) / delta`.

Stata's weeks (`%tw`) are not calendar weeks: a year has 52 of them, week `w` starts on 1 January plus
`7 * (w - 1)` days, and the 52nd lasts 8 or 9 days, so that 2000w52 is followed by 2001w1. They are
counted from 1960w1 and printed `2000w3`; with dates, `freq="SW"` puts each date in its Stata week.
`freq="W"` is the other kind, the calendar week of seven days from a Monday, which has no Stata format
and cannot be written to a `.dta`.

**Balance.** A panel is *strongly balanced* (every unit observed in every period), *weakly balanced*
(same number of periods in each unit, not the same periods) or *unbalanced*. `panel.balanced` gives
the class; `panel.has_gaps` and `panel.n_gaps` count holes inside a unit's observed range. Entry and
exit (a unit that starts late or stops early) is unbalance without a hole. Many tests accept entry and
exit but refuse holes: part 7 gives the table.

**Missing values.** Missing values (`None`, `pd.NA`, `NaN`) are `NaN` in every computation; infinities
are read as missing with a `PyxtpanelWarning`. A missing identifier or time is an error (459) unless
`missing_keys="drop"`, which drops the rows as Stata does and counts them in `panel.n_dropped`.

**Results.** Every function returns a typed object (`XtSum`, `Spells`, `TestResult`, ...) with
`.to_frame()`, `.summary()` (text), `.to_markdown()` and `.to_latex()`. A statistical test returns a
`TestResult`: `statistic`, `pvalue`, `df`, `n_obs`, `n_groups`, `options`, `details` (a table of the
quantities behind the statistic) and `notes` (remarks such as "the panel is not balanced").

**Errors.** Every error derives from `PyxtpanelError` and carries the Stata return code it
corresponds to in `stata_code` (annex C): 451 repeated (unit, period), 459 missing key or too few
distinct values, 498 a panel that must be balanced or free of holes is not, 2000 no observations, 198
invalid argument. Messages are in French or English; each has a catalogue key in `key`.

**Language.** English by default. `xt.set_options(lang="fr")` switches the messages and the texts of
the results; the `summary` of a result that carries text (a test, the audit) also takes `lang="fr"` for
one call.

```python
print(panel.n_units, panel.n_obs, panel.balanced, panel.has_gaps)
result = panel.xtsum("lgdp")
print(type(result).__name__, result.to_frame().shape)

try:
    xt.xtset(pd.concat([df, df.head(1)]), "country", "year")  # one repeated (unit, period)
except xt.DeclarationError as error:
    print(error.stata_code, "-", str(error).splitlines()[0])

previous = xt.set_options(lang="fr")
print(panel.audit().summary().splitlines()[0])  # the audit speaks French now
_ = xt.set_options(lang=previous.lang)
```

## 2. Declaring and reading a panel

**`xtset(frame, id, time, delta=, freq=, format=, missing_keys=)`** declares the panel. The identifier
may be text or numbers; the time variable may be integers, dates or Stata-coded. A frame whose index
is (unit, period) can be declared without names, and `xt.as_panel(obj)` accepts either a frame or a
panel (the entry point for estimation packages).

**Reading and writing Stata files.** `xt.read_stata(path)` returns the DataFrame and keeps what the file
says about its `xtset`; `xt.PanelData.from_stata(path)` declares it directly; `panel.to_stata(path)`
writes a `.dta` with the declaration as characteristics. `xt.read_table(path)` reads CSV or similar.

**Selecting.** `panel.select(where=, units=, periods=, complete=)` keeps rows by a condition, by units,
by periods or by the rows where listed variables are all observed (Stata's `markout`).

**Filling and balancing.** `panel.tsfill()` adds the missing rows of each unit's range (Stata's
`tsfill`); `tsfill(full=True)` fills to the range of the whole panel and gives a strongly balanced
panel; `indicator="added"` marks the new rows. `panel.balance("drop")` keeps the units observed in
every period; `balance("largest_block")` keeps the largest block of units and consecutive periods.
`xt.duplicates(frame, id, time)` reports repeated pairs before declaring.

```python
sub = panel.select(where="year >= 2000", units=["C00", "C01", "C02"])
print(sub.n_units, sub.n_periods)
complete = panel.select(complete=["lgdp", "invest"])

filled = panel.tsfill(indicator="added")  # holes become rows with missing values
print(panel.n_obs, "->", filled.n_obs, int(filled.to_frame()["added"].sum()))
strong = panel.tsfill(full=True)
print(strong.balanced)
kept = panel.balance("largest_block")
print(kept.balanced, kept.n_units, kept.n_periods)

report = xt.duplicates(pd.concat([df, df.head(2)]), "country", "year")
print(len(report))
```

## 3. Operators and transformations

**Stata notation.** `panel.operate("L.x", "L2.x", "D.x", "F.x", "S12.x", "LD.x", "L(1/3).x")` evaluates
the time-series operators `L` (lag), `F` (lead), `D` (difference) and `S` (seasonal difference) on the
calendar, with the same notation as Stata; it returns a DataFrame with one column per term. A lag
over a hole is missing. The methods `lag`, `lead`, `diff` and `seasonal` do the same for several
lags at once.

**Registered transformations.** `panel.transform("name", variables, **params)` applies a named
transformation, and each also has a method. The registry is annex A; the main ones:

| Name | Method | What it computes |
|---|---|---|
| `lag`, `lead`, `difference`, `seasonal` | `lag`, `lead`, `diff`, `seasonal` | calendar lags, leads, differences, seasonal differences |
| `growth` | `growth(periods=, log=)` | relative or log growth over `periods` calendar periods |
| `within` | `within(effects=, add_mean=)` | deviation from the mean of the unit (`effects="unit"`), of the period (`"time"`) or from both (`"twoway"`) |
| `between` | `between` | the mean of the unit, repeated on its rows |
| `twoway` | `within(effects="twoway")` | residual of x on unit and period effects, exact if balanced and iterated if not |
| `fod` | `fod(scale=)` | forward orthogonal deviations (Arellano-Bover), computed over a hole, as `xtabond2` |
| `quasi_demean` | `quasi_demean(theta)` | `x - theta * unit mean`, for random-effects transformations; `theta_re(sigma_u, sigma_e)` gives the unit-specific `theta` |
| `group_mean`, `group_count` | `group_mean(by=)` | mean or count of the group, `by="unit"` or `"time"`, repeated on each row |
| `rolling` | `rolling(window, stat=)` | mean, sum, standard deviation, min or max over the last `window` calendar periods |
| `cumulative` | `cumulative(stat=)` | running sum or mean within the unit |

Every transformation returns a Series or DataFrame on the index of the frame you declared, with
missing values where the row's own value is missing; `listwise=True` first drops, for all the
variables, the rows where any is missing (Stata's `xtdata` behaviour).

```python
terms = panel.operate("L.lgdp", "D.lgdp", "L2D.lgdp", "S5.lgdp")
print(list(terms.columns))
print(panel.lag("lgdp", [1, 2]).shape)

w = panel.within(["lgdp", "invest"])  # unit within
b = panel.between(["lgdp", "invest"])
two = panel.within("lgdp", effects="twoway")
fe = panel.within("lgdp", add_mean=True)  # like `xtdata, fe`
fod = panel.fod("lgdp")
theta = panel.theta_re(sigma_u=0.5, sigma_e=0.3)  # one theta per unit
re = panel.quasi_demean("lgdp", theta)
print(w.shape, b.shape, two.shape, theta.shape)

m3 = panel.rolling("lgdp", 3, "mean")
growth = panel.growth("lgdp", periods=1, log=True)
mean_by_year = panel.group_mean("lgdp", by="time")
by_unit = panel.transform("group_count", "lgdp")
running = panel.cumulative("invest")
print(m3.name, growth.name)
```

## 4. Describing a panel

All the descriptions return objects that print like Stata's output (`summary()`) and give a DataFrame
(`to_frame()`).

* **`xtdescribe()`** (`xtdescribe`): units, periods, delta, the distribution of the number of periods
  per unit `T_i`, and the most frequent participation patterns (`patterns=`, `width=`).
* **`xtsum(variables, by=)`** (`xtsum`): overall, between and within mean, standard deviation, min and
  max, with the numbers of observations, units and the mean `T`. The within deviation adds the
  overall mean back, as Stata does.
* **`xttab(variable)`** (`xttab`): overall, between and within frequencies of the values of a variable.
* **`xttrans(variable, gaps=, bins=)`** (`xttrans`): transition matrix between a unit's consecutive
  observations; `gaps="bridge"` (default, as Stata) pairs a row with the unit's next row even across a
  gap, `gaps="strict"` only consecutive periods; `bins=` bins a continuous variable.
* **`xtcorr(variables, kind=)`**: overall, within and between correlations, with the number of points
  behind each pair; it has no Stata command and is compared with `correlate` (after `xtdata, fe` and
  `xtdata, be` for the within and between ones) on the rows where both variables are observed. A
  correlation with a variable that is constant (up to rounding) is missing, not a number made of
  rounding noise; with fewer than two points the correlation is missing and `n` says how many
  points there were (Stata's `correlate` stores nothing then).
* **`autocorr(variable, lags=)`**: within autocorrelation at calendar lags.
* **`attrition()`**: entries, exits and re-entries by period, and unit survival.
* **`missing_report(variables)`** and **`complete_cases(variables)`**: missing values by variable,
  period and unit with joint patterns, and the rows and units kept if all of some variables must be
  observed.

```python
print(panel.xtdescribe().summary().splitlines()[0])
table = panel.xtsum(["lgdp", "invest"], by="reform")
print(table.to_frame().head())

tab = panel.xttab("reform")
trans = panel.xttrans("reform")
binned = panel.xttrans("lgdp", bins=3)
print(trans.to_frame())

corr = panel.xtcorr(["lgdp", "invest"])
print(corr.corr["within"].round(2))
ac = panel.autocorr("lgdp", lags=[1, 2])
at = panel.attrition()
miss = panel.missing_report(["lgdp", "invest"])
cc = panel.complete_cases(["lgdp", "invest"])
```

## 5. Exploration tools

The exploration tools have no Stata command. Each states in its docstring and in the registry (annex A)
what the nearest equivalent is and how it is validated: hand-computed small panels, an independent
implementation in plain Python loops compared on random panels, invariance tests (shuffled rows,
renamed units, shifted time origin) and, where a Stata computation of the same quantity exists,
a comparison with Stata.

| Question | Tool |
|---|---|
| Why is my panel unbalanced? | `gap_types()`, `unbalancedness()`, `balance_cost()` |
| What does an estimator have left? | `usable_equations()`, `residual_variation()` |
| Which variables carry information? | `variation()`, `spells()`, `event_time()` |
| Are there entry errors? | `within_outliers()`, `level_shifts()` |
| Reshaping | `to_wide()`, `to_long()`, `to_array3d()`, `resolve_duplicates()` |
| Everything at once | `audit()` |

* **`gap_types()`** classifies each unit as complete, late entry, early exit, internal gap or
  re-entry, lists the gaps with their length and the period of re-entry, and counts units by type.
* **`unbalancedness()`** gives the Ahrens-Pincus measures `gamma` and `nu` (as `plm::punbalancedness`),
  the mean and harmonic mean of `T_i` and the fill rate.
* **`balance_cost()`** tabulates what balancing costs: for each threshold `T_i >= k` and each window of
  consecutive periods, the units and observations kept, and the largest balanced block. The block is
  the one `balance("largest_block")` keeps.
* **`usable_equations(variable, lags=, transform="fd" | "fod", exog=)`** counts, by unit and by period,
  the equations a dynamic panel GMM can use after the transformation and the lags; the counts match
  `pyxtabond2` on five panels.
* **`residual_variation(variables, transforms=("within", "fd", "fod"))`** gives the standard deviation
  and variance share left after each transformation, to choose an estimator.
* **`variation(variables)`** classifies each variable as constant, time-invariant, unit-invariant or
  varying, and counts the units and periods in which it varies.
* **`spells(variable)`** measures the episodes of a discrete state, cut by a change, a gap or a missing
  value, with open ends marked, changes per unit and the duration of complete episodes.
* **`event_time(indicator=)`** or **`event_time(event_period=)`** gives the time relative to an event, the coverage by
  relative time and the bounds common to all units.
* **`within_outliers(variable)`** flags values by a modified z-score (median and MAD inside the unit)
  without deleting anything; **`level_shifts(variable, window=)`** finds jumps of the median between
  windows and the units that dominate the within variation.
* **`to_wide(variables)`** (`reshape wide`) gives one row per unit; **`xt.to_long(wide, stubs, id=)`**
  (`reshape long`) goes back; **`to_array3d(variables)`** gives a variable x period x unit array with
  the mask of observed cells; **`xt.resolve_duplicates(frame, id, time, how=)`** merges repeated
  (unit, period) pairs (`first`, `last`, `mean`, `sum`, `min`, `max`, `median`) and reports what was merged.
* **`audit(dependent=)`** runs all the checks (structure, balance, missing values, variation,
  outliers, equations) and returns an `AuditReport` of alerts with advice in French and English,
  exportable with `summary`, `to_markdown`, `to_html` (shown by Jupyter) and `to_latex`.

```python
gaps = panel.gap_types()
print(gaps.counts)
print(panel.unbalancedness().gamma)
cost = panel.balance_cost()
print(cost.to_frame().head(3))

eq = panel.usable_equations("lgdp", lags=1, transform="fod")
rv = panel.residual_variation(["lgdp"])
var = panel.variation(["lgdp", "reform"])
print(var.to_frame()[["kind", "units_varying"]])

sp = panel.spells("reform")
ev = panel.event_time(indicator="reform")
print(ev.to_frame().head(3))

out = panel.within_outliers("lgdp")
shifts = panel.level_shifts("lgdp", window=3)

wide = panel.to_wide(["lgdp"])
back = xt.to_long(wide.reset_index(), ["lgdp"], id="country", time="year", drop_empty=True)
cube = panel.to_array3d(["lgdp"])
print(wide.shape, len(back), cube.values.shape)

dup = pd.concat([df, df.head(3)])
resolved = xt.resolve_duplicates(dup, "country", "year", how="mean")
print(resolved.n_rows_in - resolved.n_rows_out)

report = panel.audit(dependent="lgdp")
print(report.summary().splitlines()[0])
```

## 6. Figures

The figures need the extra `pyxtpanel[plot]` (matplotlib); nothing else in the package imports it.
`panel.plot(kind, **options)` returns a matplotlib `Figure` in one common style, in English or French
(`lang=`, decimal comma in French), and the figure saves to PNG, PDF or SVG with `fig.savefig`.

| Kind | Shows | Main options |
|---|---|---|
| `"lines"` | trajectories of a sample of units, with the mean by period (like `xtline`) | `variable`, `n_units`, `units`, `highlight`, `panels` |
| `"scatter"` | between (unit means) and within (deviations) scatter of x and y | `x`, `y`, `part` |
| `"distribution"` | distribution of a variable between and within units | `variable`, `style` (`"hist"` or `"box"`) |
| `"participation"` | map of the periods each unit is observed in, sorted by pattern | |
| `"states"` | units by periods coloured by a discrete state (like `panelView`) | `variable`, `sort` |
| `"transition"` | heat map of the transition matrix | `variable`, `gaps`, `bins` |
| `"attrition"` | units observed by period, entries, exits, Kaplan-Meier survival | |
| `"balance_cost"` | share kept by minimum length and by window, largest block marked | |

```python
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

figures = [
    panel.plot("lines", variable="lgdp", n_units=6),
    panel.plot("scatter", x="invest", y="lgdp"),
    panel.plot("distribution", variable="lgdp"),
    panel.plot("participation"),
    panel.plot("states", variable="reform"),
    panel.plot("transition", variable="reform"),
    panel.plot("attrition"),
    panel.plot("balance_cost"),
    panel.plot("lines", variable="lgdp", n_units=4, lang="fr"),
]
print(len(figures))
plt.close("all")
```

```python skip
figures[0].savefig("lgdp_lines.png", dpi=150)
```

A figure called wrongly is refused before anything is drawn, with a `SpecificationError` that
names the figure, says what is wrong and lists the arguments the figure accepts (in the language
given by `lang=`, or the global one): `panel.plot("states")` answers
`figure 'states': missing a required argument: 'variable'. Arguments of the figure: variable,
sort='first_change', max_states=8, lang=None.` instead of a bare `TypeError`.

## 7. Statistical tests guide

The package contains the tests of panel data that do not depend on an estimator. Every one of them
has been run in Stata 17 on shared datasets, the results frozen in `tests/golden` and compared to
1e-9 or better in the test suite: **a test without a Stata case cannot be registered** (decision of
2026-11-07). The list with its Stata command, its R equivalent and its reference is annex B.

### 7.1 Three categories

| Category | Input | Where it lives | Example |
|---|---|---|---|
| **A** | a declared panel and variable names | here | `im_pesaran_shin(panel, "x", lags=1)` |
| **B** | bare arrays (residuals, coefficients, moments) | here | `hansen(moments, weight, n_parameters)` |
| **C** | the internals of an estimator (its transformed regressors, its instruments, its weight matrix) | in the estimation packages (`pyxtabond2`, `pyxthreg`, `pyxtdpthresh`) | the specification tests of one estimator |

A test goes in category A when its input is a panel and a variable (it needs the unit and the period),
in category B when it needs only numbers, and in category C when it cannot be computed without the
estimator's own conventions. The estimation packages call the category B functions with the arrays
they already hold (`from pyxtpanel.diagnostics import hansen, ar_test`); a user calls the same
functions with arrays from any estimator.

### 7.2 Calling a test

```python
fp = xt.xtset(full, "country", "year")  # the complete panel: no holes

# category A: through the method of its family ...
ips = fp.unit_root("lgdp", test="im_pesaran_shin", lags=1)
cd = fp.cross_section_dependence("invest", test="pesaran_cd")
wool = fp.serial_correlation("lgdp", ["invest"], test="wooldridge")
kao = fp.cointegration("lgdp", ["invest"], test="kao")
slopes = fp.slope_heterogeneity("lgdp", ["invest"], test="slope_homogeneity")
granger = fp.causality("lgdp", "invest", test="dumitrescu_hurlin")

# ... or by name, which also takes an alias ("ips" is "im_pesaran_shin")
same = xt.run_test("ips", fp, "lgdp", lags=1)
assert same.statistic == ips.statistic
pedroni = xt.run_test("pedroni", fp, "lgdp", ["invest"], kernel_lags=2)

# the registry
catalogue = xt.list_tests()
print(catalogue[["family", "category", "balanced"]].head(3))
print(len(catalogue), "tests")
```

The families and their methods: `unit_root` (`levin_lin_chu`, `harris_tzavalis`, `breitung`,
`im_pesaran_shin`, `fisher`, `cips`, `hadri`), `cross_section_dependence` (`pesaran_cd`,
`lm_cross_section`, `friedman`, `frees`), `serial_correlation` (`wooldridge`, `panel_durbin_watson`),
`cointegration` (`kao`, `pedroni`, `westerlund`), `slope_heterogeneity` (`slope_homogeneity`) and
`causality` (`dumitrescu_hurlin`). A call without its `test=` is refused with the list of the ones the
family has: there is no default, because the tests answer different questions.

### 7.3 Which test for which question

| You want to know | Use | Null hypothesis | Notes |
|---|---|---|---|
| Does a series have a **common** unit root? | `levin_lin_chu`, `harris_tzavalis`, `breitung` | unit root | `harris_tzavalis` is for short `T` and large `N`; `breitung` needs no bias correction, `robust=True` for cross-sectional dependence |
| Do **some** units have a unit root, each with its own root? | `im_pesaran_shin`, `fisher` | all units have a unit root | `fisher` accepts holes; four combinations by `method` |
| The same, with a **common factor** across units | `cips` | unit root | the t-bar against the critical values of the paper; the p-value of Z-t-bar is oversized under a strong factor |
| Is a series **stationary** (null reversed)? | `hadri` | every unit is stationary | rejecting means a unit root; `robust=True` for heteroskedasticity |
| Are the **units correlated** in the cross-section? | `pesaran_cd`, `lm_cross_section`, `friedman`, `frees` | no cross-sectional dependence | `pesaran_cd` for many units and small `T`; `lm_cross_section` when `N` is small; `friedman` and `frees` are rank-based |
| Do two series **cointegrate**? | `kao`, `pedroni`, `westerlund` | no cointegration | `kao` and `pedroni` have panel and group forms; `westerlund` is the error-correction test |
| Does `x` **Granger-cause** `y` in the panel? | `dumitrescu_hurlin` | no unit has causality | allows heterogeneous coefficients |
| Is there **serial correlation** in a panel regression? | `wooldridge`, `panel_durbin_watson` | no serial correlation | `panel_durbin_watson` gives a statistic without a p-value |
| Are residuals **heteroskedastic**? | `breusch_pagan`, `white`, `wald_groupwise` | homoskedasticity | `wald_groupwise` (`xttest3`) tests group-wise variances and is oversized at moderate `T` |
| Are there **unit effects**? | `f_effects`, `lm_random_effects`, `baltagi_li` | no effects | `lm_random_effects` is one-sided and needs `variance_u`; `baltagi_li` is joint with serial correlation |
| Are the **slopes** the same in every unit? | `slope_homogeneity` | slope homogeneity | `hac=True` for serial correlation; `adjusted=False` for the unadjusted statistic |
| Are residuals **normal**? | `jarque_bera`, `shapiro_wilk`, `shapiro_francia`, `skewness_kurtosis` | normality | |
| Is the **specification** right? | `reset`, and `vif` for collinearity | no omitted nonlinearity | `vif` returns the factors, not a test |
| Do the **instruments** hold after a GMM? | `hansen`, `sargan`, `diff_hansen`, `ar_test` | valid instruments; no `AR(k)` | `hansen` after an efficient GMM, `sargan` after two-stage least squares |
| **Fixed or random effects**? | `hausman` | the efficient estimator is consistent | takes the two coefficient vectors and covariances |
| A **joint restriction** on coefficients? | `wald` | `R b = r` | chi-square, or `F` with `df_resid` |

`cross_correlation(panel, variable)` gives the correlation of a variable between every pair of units
(the matrix behind the cross-sectional tests); unlike the tests, it does not stop on a pair that is too
short.

### 7.4 Reading a result

```python
print(ips.summary())  # the Stata-like text: hypotheses, statistic, p-value
print(ips.statistic, ips.pvalue, ips.df, ips.n_obs, ips.n_groups)
print(ips.rejects(0.05))  # True when the null is rejected at 5 %
print(ips.to_frame())  # one row, to stack the results of several tests
print(ips.details)  # the quantities behind the statistic
print(ips.notes)  # remarks: unbalanced panel, instruments outnumber units, ...
print(ips.summary(lang="fr"))
```

`statistic` and `pvalue` are floats (`NaN` when the test is not defined, for example a model that is not
over-identified: a remark in `notes` says why). `null` and `alternative` are the hypotheses in words;
`options` is what was actually used, defaults included, so that a result can be reproduced from its own
print.

### 7.5 What a test needs

**Balance and holes.** A panel with **entry and exit** (units that start late or stop early) is
*unbalanced*; a panel with **holes** is missing periods inside a unit's range. The table below comes
from running each category A test on a balanced panel, on a panel with late entry and on a panel with
two missing periods; a test of this document runs it again, so it cannot go stale.

<!-- balance-table:start -->
| Test | Balanced | Entry and exit | Holes |
|---|---|---|---|
| `breitung` | yes | refused (498) | refused (498) |
| `hadri` | yes | refused (498) | refused (498) |
| `harris_tzavalis` | yes | refused (498) | refused (498) |
| `levin_lin_chu` | yes | refused (498) | refused (498) |
| `kao` | yes | refused (498) | refused (498) |
| `dumitrescu_hurlin` | yes | refused (498) | refused (498) |
| `im_pesaran_shin` | yes | yes | refused (498) |
| `cips` | yes | yes | refused (498) |
| `pedroni` | yes | yes | refused (498) |
| `westerlund` | yes | yes | refused (498) |
| `fisher` | yes | yes | yes |
| `frees` | yes | yes | yes |
| `friedman` | yes | yes | yes |
| `lm_cross_section` | yes | yes | yes |
| `pesaran_cd` | yes | yes | yes |
| `slope_homogeneity` | yes | yes | yes |
| `panel_durbin_watson` | yes | yes | yes |
| `wooldridge` | yes | yes | yes |
<!-- balance-table:end -->

A refusal is the Stata return code 498 with a message that says what to do (`balance()`, `tsfill()`,
`select(...)`). Category B tests take no panel, so the question does not arise, but the arrays they
receive must be consistent: `offsets` starts at 0, increases, and ends at the number of rows.

**Length and size.** A panel too short for the test is refused with 459 and a message that gives the
minimum (fewer equations than coefficients in a unit's regression, fewer common periods than
`min_overlap` in the cross-sectional tests). Where the choice matters there is no hidden default for
the lags: `cips`, `fisher` and `pedroni` ask for `lags` or `kernel_lags`.

**Missing values.** Rows with a missing value are ignored for that variable; a unit left with too few
observations is dropped with a note, or the call is refused, as Stata does (annex D says where).

**Degenerate data.** A constant series, a perfect fit, a zero weight or an empty panel is refused with a
message from the catalogue, in French and in English, never a raw numerical warning or a `NaN`
statistic. The test suite runs every category A test on a set of degenerate panels.

### 7.6 Known behaviours: where the size is not the nominal one

The rejection rates under the null and under a close alternative are simulated in `docs/size_power.md`
(2000 replications, 5 % level). Some tests reject too often at a moderate `T` or `N`. **This is the
behaviour of Stata's command, which the package reproduces to 1e-9, not a defect of the port**: the same
design simulated in Stata gives the same rate.

| Test | Size at 5 % | Design | Advice |
|---|---|---|---|
| `levin_lin_chu` | about 27 % | 30 units, 50 periods, one lag | prefer `im_pesaran_shin` or `breitung` |
| `harris_tzavalis` | about 9 % | 30 units, 50 periods | it is made for fixed `T` and large `N` |
| `fisher` | about 8 % | 30 units, 50 periods | acceptable; the Phillips-Perron version about 11 % |
| `cips` (Z-t-bar p-value) | about 10 % | strong common factor | use the t-bar against the critical values, as the paper does (size 4 %) |
| `im_pesaran_shin` with a common factor | about 24 % | strong common factor | use `cips` |
| `pedroni` | 12 to 15 % | 20 units, 40 periods | oversized at short `T`, as Pedroni (2004) reports |
| `wald_groupwise` | about 46 % | 20 units, 30 periods | reproduces `xttest3`, which is oversized |
| `friedman` | right only if `N = T` | | the p-value uses `N - 1` degrees of freedom as `xtcsd` does; `dof="periods"` uses `T - 1` |
| `baltagi_li` (adjusted versions) | about 20 % with a large unit effect | | the adjustments are robust only locally |
| `lm_random_effects` without `variance_u` | about 10 % | | pass the random-effects variance so that the sign rule applies |
| `dumitrescu_hurlin` | about 7 to 8 % | 20 units, 30 periods | the standardisation assumes normal errors |

### 7.7 Category B: from any estimator

The inputs are plain arrays. `fp.arrays([...])` (part 8) gives them in the layout the tests expect,
with the `offsets` that delimit each unit.

```python
a = fp.arrays(["lgdp", "invest"])
y, x = a.values[:, 0], a.values[:, 1:]
X = np.column_stack([np.ones(len(y)), x])
beta = np.linalg.lstsq(X, y, rcond=None)[0]
e = y - X @ beta  # pooled residuals

effects = xt.run_test("lm_random_effects", e, a.offsets)
joint = xt.run_test("baltagi_li", e, a.offsets)
het = xt.run_test("breusch_pagan", e, x)
white = xt.run_test("white", e, x)
groupwise = xt.run_test("wald_groupwise", e, a.offsets)
normal = [
    xt.run_test(name, e)
    for name in ("jarque_bera", "shapiro_wilk", "shapiro_francia", "skewness_kurtosis")
]
spec = xt.run_test("reset", y, x)

# an F test of the unit effects from two residual sums of squares
within = fp.within(["lgdp", "invest"]).to_numpy()
bw = np.linalg.lstsq(within[:, 1:], within[:, 0], rcond=None)[0]
rss_fe = float(np.sum((within[:, 0] - within[:, 1:] @ bw) ** 2))
fe = xt.run_test("f_effects", float(e @ e), rss_fe, len(y), fp.n_units, 1)

# a Wald test and a Hausman test from coefficients and covariances
V = np.linalg.inv(X.T @ X) * float(e @ e) / (len(y) - 2)
joint_zero = xt.run_test("wald", beta, V)
hausman = xt.run_test("hausman", bw, beta[1:], np.eye(1) * 0.02, np.eye(1) * 0.01)
print(effects.pvalue < 0.05, joint_zero.n_obs, len(normal))
```

`hansen` and `sargan` take the moments `sum_i Z_i' e_i` (or the unit-by-unit moments), the weight matrix
of the criterion and the number of parameters; `diff_hansen` takes two `J` statistics; `ar_test` takes
`Equations` (residuals, regressors, `offsets`, `time_index`, instruments of a system) with the weight
and covariance matrices of the estimate, and gives the `AR(k)` of `xtabond2` on balanced panels, with
holes and with entry and exit. They need the internals of an estimator, so only `hansen` is run here
with simulated moments.

```python
from pyxtpanel.diagnostics import moment_covariance

moments = rng.normal(size=(24, 4))  # unit-by-unit moments of 4 instruments
weight = np.linalg.inv(moment_covariance(moments))
j = xt.run_test("hansen", moments, weight, 2)
print(j.df, j.statistic >= 0)
```

## 8. Authors of estimation packages

*Provisional version.* The plan finalises this part after the migrations of `pyxtabond2`, `pyxthreg`
and `pyxtdpthresh` to this package (item H.4.1); it describes the contract as it stands.

An estimation package depends on `pyxtpanel` for the panel, never for the estimation. What it may
rely on is deliberately narrow, versioned, and made of plain NumPy arrays.

### 8.1 Receiving a panel

Accept either a declared panel or a DataFrame, and read the data through `arrays`:

```python
from pyxtpanel import PanelLike, as_panel
from pyxtpanel.contract import CONTRACT_VERSION, PanelArrays


def fit(data, dependent, regressors):
    """What an estimator does with a panel (a sketch: no estimation here)."""
    assert CONTRACT_VERSION == 1  # checked once, at import time
    panel = as_panel(data, "country", "year") if not isinstance(data, PanelLike) else data
    arr: PanelArrays = panel.arrays([dependent, *regressors])
    return arr


arr = fit(df, "lgdp", ["invest"])
print(arr.values.shape, arr.offsets[:3], arr.contract_version)
```

`as_panel(obj, id, time)` returns a `PanelData` from a frame (using `xtset`) or the panel unchanged;
`PanelLike` is the protocol an estimator may test (`id_name`, `time_name`, `n_units`, `n_obs`,
`arrays`). Every `PanelData` satisfies it.

### 8.2 The array contract

`panel.arrays(columns, where=)` returns a frozen `PanelArrays` with these guarantees; a breaking change
raises `CONTRACT_VERSION` (adding a field does not):

| Field | Meaning |
|---|---|
| `values` | `float64`, C-contiguous `(n_obs, k)`, a copy the caller owns; missing values and infinities are `NaN` |
| `columns` | names of the columns of `values` |
| `offsets` | `int64` `(n_units + 1,)`: unit `i` owns the rows `offsets[i]:offsets[i + 1]` |
| `unit_labels`, `unit_of_row` | original identifier of each unit; unit code of each row |
| `time_index` | calendar position `(t - t0) / delta` of each row, increasing within a unit; consecutive periods differ by exactly 1 whatever `delta` or the date frequency |
| `source_row` | position of each row in the declared frame |
| `valid` | no requested column is missing on the row |
| `n_source`, `delta`, `contract_version` | length of the declared frame, calendar step, version |

Rows are sorted by unit then period, with no repeated (unit, period) pair. The read-only index arrays
are shared with the panel, not copied. The methods help an estimator without pandas:

* `lag_index(k)`: for each row, the row of the same unit `k` calendar periods earlier, or `-1`;
* `unit_slices()`: one `slice` of rows per unit;
* `reindex_back(array, fill=)`: put a result computed in sorted order back in the order of the frame;
* `as_cube()`: `values` as a `(n_units, n_time, k)` view, for a strongly balanced panel only.

```python
lag = arr.lag_index(1)
has_lag = lag >= 0
print(has_lag.sum(), len(arr.unit_slices()))

fitted = np.full(len(arr.values), 1.0)  # any result in sorted order ...
back = arr.reindex_back(fitted)  # ... returned in the order of the user's frame
print(back.shape)
cube = fp.arrays(["lgdp"]).as_cube()
print(cube.shape)
```

### 8.3 Equations and tests

An estimator holds *equations* (the transformed observations it estimates, after first differences or
forward orthogonal deviations, with their instruments). `pyxtpanel.diagnostics.Equations` carries them:
`residuals`, `regressors`, `offsets`, `time_index` and optional `instruments`. A calendar gap is a gap
in `time_index`, so the lag of a missing period does not exist. The category B functions take these
objects or bare arrays (part 7): `hansen`, `sargan`, `diff_hansen`, `hausman`, `wald`, `ar_test` and the
residual tests.

```python
from pyxtpanel.diagnostics import Equations

system = Equations(
    residuals=arr.values[:, 0] - arr.values[:, 0].mean(),
    regressors=arr.values[:, 1:],
    offsets=arr.offsets,
    time_index=arr.time_index,
)
print(system.residuals.shape, system.offsets.shape)
```

### 8.4 Errors, language and versioning

* **Errors.** Raise or let through `PyxtpanelError` subclasses; they carry `stata_code` and a catalogue
  `key`. Messages exist in French and English, selected by `xt.set_options(lang=...)`. Annex C lists
  the codes.
* **Version.** `CONTRACT_VERSION` is an integer that changes only with a breaking change of the
  array layout. The package follows semantic versioning; its public API is a snapshot tested
  (`tests/contract/api_snapshot.json`) and frozen at 1.0. Pin `pyxtpanel>=0.3,<2` once released.
* **Dependencies.** `numpy`, `pandas`, `scipy`; Numba stays in the estimators (the contract has no
  compiled kernel).

## 9. Method and validation

This part gives the formulas where Stata's convention is not the textbook, the way the work was
validated and how to read the evidence.

### 9.1 Conventions

* **Calendar, not rows.** Operators and transformations work on the calendar of the panel: `L.x` is
  the value one `delta` earlier, missing across a hole. Stata's `L.` after `xtset` does the same; a
  row shift (`groupby().shift`) does not.
* **Missing values.** An observation with a missing value in a variable is out of that variable's
  computation; unless `listwise=True`, each variable uses its own observed rows.
* **Unbalanced panels.** Every statistic that depends on `T` uses the unit's own `T_i`, with the
  documented exception of the tests whose reference command uses the longest unit (annex D).

### 9.2 Formulas

**`xtsum`.** For a variable with $N$ observed values in $n$ units, $\bar x$ the overall mean and
$\bar x_i$ the mean of unit $i$: the overall standard deviation is
$\sqrt{\sum (x_{it}-\bar x)^2/(N-1)}$; the between standard deviation is
$\sqrt{\sum_i (\bar x_i-\bar x_B)^2/(n-1)}$ with $\bar x_B$ the mean of the unit means; the within
standard deviation is $\sqrt{\sum (x_{it}-\bar x_i)^2/(N-1)}$, and the within minimum and maximum are
those of $x_{it}-\bar x_i+\bar x$. The divisor $N-1$, not $N-n$, is Stata's and makes the three pieces
comparable. In an unbalanced panel $\bar T = N/n$.

**Within and two-way within.** `within` is $x_{it}-\bar x_i$; with `add_mean=True` the overall mean of the
observed rows is added (`xtdata, fe`). The two-way within of a complete, strongly balanced panel is
$x_{it}-\bar x_i-\bar x_t+\bar x$; otherwise the unit means and period means of the observed rows are
removed in turn until the largest correction is below `tol`, which converges to the residual of the
regression on unit and period dummies (a non-convergence is reported, or raised).

**Forward orthogonal deviations.** For a row with $F$ later observed values of the same unit,
$$x^*_{it}=\sqrt{\frac{F}{F+1}}\Bigl(x_{it}-\frac1F\sum_{s\ \mathrm{later}}x_{is}\Bigr),$$
missing where $F=0$. "Later" means the next observed rows, whatever the gap between their dates: this is
`xtabond2`'s definition, and the transformed errors stay uncorrelated when the originals are.

**Random-effects transformation.** `quasi_demean` is $x_{it}-\theta_i\bar x_i$; `theta_re(sigma_u,
sigma_e)` gives $\theta_i = 1-\sigma_e/\sqrt{\sigma_e^2+T_i\sigma_u^2}$ with the unit's own $T_i$.

**`xttrans`.** Pairs of a unit's consecutive observed rows: the row of period $t$ with the unit's next row
(gaps bridged, as Stata); `gaps="strict"` keeps pairs one period apart only.

**Unbalancedness (`unbalancedness`).** With $n$ units, $T_i$ the periods of unit $i$, $\bar T$ their mean
and $N=\sum_i T_i$ the number of observations, as `plm::punbalancedness` (Ahrens and Pincus):
$$\gamma=\frac{n}{\bar T\sum_i 1/T_i},\qquad \nu=\frac{1}{n\sum_i (T_i/N)^2}.$$
Both equal 1 for a balanced panel and decrease as unbalance grows ($0<\gamma\le 1$, $1/n<\nu\le 1$).

**The tests.** Each test follows the formulas of its reference (annex B) with Stata's finite-sample
choices where they differ; the places where Stata departs from the paper are listed in annex D and each
has a test that fixes the behaviour.

### 9.3 How the package is validated

1. **Against Stata 17 MP.** Every function with a Stata equivalent is run in Stata (batch, through
   PowerShell) on shared simulated datasets, and the outputs are frozen in `tests/golden/*.json`; the
   test suite compares the package with the frozen values without Stata. The rule is an absolute
   error of at most 5e-5 *or* a relative error of at most 1e-6, the floor taken from `pyxtdpthresh`;
   most comparisons reach 1e-9 or better. Annex E lists each frozen domain, the Stata build and the
   number of cases. The cases include the awkward ones: unbalanced panels, entry and
   exit, holes, one-observation units, constant series, `delta` greater than 1.
2. **Where Stata has no command.** The exploration tools (part 5) are checked by hand-computed small
   panels, by an independent implementation in plain Python loops compared with the vectorised code on
   random panels, by invariance (shuffled rows, renamed units, shifted origin of time), and by a
   direct Stata computation of the same quantity where one exists (`egen`, `bysort`, `reshape`,
   `tsreport`, lags on the calendar).
3. **Against other software.** `usable_equations` counts agree with `pyxtabond2` on five panels; the
   over-identification and `AR(k)` tests agree with `xtabond2` and `pyxtdpthresh`; `baltagi_li`,
   `dumitrescu_hurlin` and `cips` are also compared with R's `plm` 2.6.7 (`pbsytest`, `pgrangertest`,
   `cipstest`), and the MacKinnon critical values with `statsmodels`.
4. **A random campaign.** `validation/campaign.py` draws 160 panels at random where a hand-written
   list is thin: from one row to 60 units, holes, entry and exit, units with a single row, missing
   values from none to all, constant variables, a repeated pair, any identifiers in any order, any
   origin and step of time, Stata days, weeks, months, quarters, half-years and years, shuffled rows.
   Stata and the package each
   compute 73 quantities per panel (declaration, `xtsum` of two variables, six operators, within and
   between, `xtdescribe`, `xttrans`, `xtcorr`, `tsfill`) and every cell agrees to 1e-9. A wider run of
   3 700 fresh panels (not frozen, three seeds) agrees too. The same panels, with identifiers relabelled, rows
   shuffled and time shifted or rescaled, give the same statistics
   (`tests/properties/test_campaign_invariance.py`, no Stata needed). The campaign found one defect
   (a constant variable had a correlation of rounding noise in `xtcorr`), fixed, and one divergence
   (9.4).

   **How close can an agreement be?** Stata does not agree with itself to the last digit, so a
   tolerance tighter than its own noise would ask the impossible. `validation/self_agreement.py` ran
   the 160 panels again in Stata with the units renamed in another order (which changes the order of
   its floating-point additions) and measured, cell by cell, the error `|a - b| / max(|b|, 1)`
   (`tests/golden/stata_self_agreement.json`): Stata against itself, at most **5e-14**; the package
   against Stata, at most **4e-13**, on 11 680 cells, none above 1e-12. Closed forms (counts, moments,
   transformations) therefore agree to about 1e-12, three orders of magnitude inside the 1e-9 the
   campaign checks and the 1e-10 the project aims at. Reading the value from the screen instead of
   from `r()` is far worse: written with `%12.0g` (the default of `string()`) the error reaches 2e-6,
   with `%9.0g` (the width of Stata's tables) 5e-3 on large sums. This is why the frozen references
   are exported from Mata at 15 digits and never read from a log. The 5e-5 / 1e-6 floor stays for the
   quantities that are not closed forms (iterative or simulated).
5. **Edge cases.** Every category A test is run on degenerate panels (a single unit, a constant, a
   unit entirely missing, one period); each call returns a finite statistic or a catalogue error in
   both languages, never a raw exception, warning or `NaN`.
6. **Size and power.** The tests are simulated under the null and under close alternatives
   (`docs/size_power.md`); a rate far from nominal is documented, with the Stata rate beside it where
   the cause is Stata's formula (part 7.6).
7. **Performance.** `docs/performance.md`: one million rows are declared and operated on in well under
   a second, and a guard (`pytest -m bench`) keeps it so.
8. **Version matrix.** The whole suite runs on pandas 2 and on pandas 3.

### 9.4 Limits

* `pedroni`'s automatic kernel bandwidth is not implemented: `kernel_lags` is required.
* `kao` and `dumitrescu_hurlin` require a balanced panel; `cips` refuses holes (Stata's `pescadf`
  accepts them and the package does not reproduce that).
* `slope_homogeneity` with `hac=True`, `whitening=True` and three or more regressors can differ from
  `xthst` because the signs of the singular vectors of Mata's `svd` differ from LAPACK's for about a
  quarter of the units; it is documented and not validated.
* `xtcorr` reports the number of points of a pair that has fewer than two, with a missing
  correlation; Stata's `correlate` stores nothing (the campaign reads neither).
* Nested panels (a firm in a sector in a country) are not built: decided on 2026-11-17.
* The package estimates no model: use the estimation packages for that.


# Annexes

## Annex A. Functions

<!-- annex:a:start -->

### Declaring, selecting, reading

| Call | Stata | R or Python equivalent | What it does |
|---|---|---|---|
| `xtset(frame, id, time, delta=)` | xtset id t, delta() | plm::pdata.frame; linearmodels PanelData | — |
| `PanelData.tsfill()` | tsset / tsfill | plm::make.pbalanced, make.pconsecutive; panelr balance_panel | — |
| `PanelData.tsfill(full=True)` | tsfill, full | plm::make.pbalanced | — |
| `duplicates(frame, id, time)` | duplicates report / list | plm (refuses duplicates) | — |
| `read_stata(path)` / `PanelData.to_stata(path)` | use / save (xtset kept) | — | — |
| `PanelData.select(where=, units=, periods=, complete=)` | if / in | — | — |
| `PanelData.select(complete=[...])` | markout | — | — |
| `PanelData.operate(*terms)`, `lag`, `lead`, `diff`, `seasonal` | L. F. D. S. LD. L(1/3). | — | — |
| `PanelData.within(add_mean=True, listwise=True)` / `between(listwise=True)` | xtdata, fe / be | — | — |

### Describing

| Call | Stata | R or Python equivalent | What it does |
|---|---|---|---|
| `PanelData.attrition()` | — | — | Entries, exits and re-entries of units by period, and unit survival. |
| `PanelData.autocorr()` | — | — | Within autocorrelation by calendar lag (addition of the package). |
| `PanelData.complete_cases()` | — | — | Rows and units kept when all of some variables must be observed, and T_i effect. |
| `PanelData.missing_report()` | — | collapse::fnobs | Missing values by variable, period and unit, with joint patterns. |
| `PanelData.xtcorr()` | — | plm::cortab; collapse::pwcor, pwnobs | Overall, within and between correlations with the number of points per pair. |
| `PanelData.xtdescribe()` | xtdescribe | plm::pdim (size only) | Panel structure: units, periods, distribution of T_i, participation patterns. |
| `PanelData.xtsum()` | xtsum | collapse::qsu(by = pid) (equivalence to xtsum not confirmed) | Overall, between and within mean, standard deviation, min and max. |
| `PanelData.xttab()` | xttab | — | Overall, between and within frequencies of the values of a variable. |
| `PanelData.xttrans()` | xttrans | — | Transition matrix between a unit's consecutive observations. |

### Operators and transformations

| Call | Stata | R or Python equivalent | What it does |
|---|---|---|---|
| `PanelData.transform("between")` | xtdata, be (rows kept: those with a value) | plm::Between; collapse::fbetween; linearmodels mean | Mean of the unit over its observed values, repeated on each row. |
| `PanelData.transform("cumulative")` | bysort id (t): gen c = sum(x) | — | Running sum or mean of the unit's observed values up to each row. |
| `PanelData.transform("difference")` | D.x, D2.x, ... | plm::diff; collapse::fdiff | Difference with the previous period, repeated `order` times. |
| `PanelData.transform("fod")` | xtabond2, orthogonal  (forward orthogonal deviations) | — | Deviation from the mean of the later observed values of the unit, scaled. |
| `PanelData.transform("group_count")` | bysort id: egen n = count(x) | collapse::fnobs; linearmodels count | Number of observed values of the group (unit or period) repeated on each row. |
| `PanelData.transform("group_mean")` | bysort id: egen m = mean(x) | collapse::fbetween | Mean of the group (unit or period) repeated on each row. |
| `PanelData.transform("growth")` | D.x / L.x  or  D.ln(x) | collapse::fgrowth | Growth rate over `periods` calendar periods (relative or log difference). |
| `PanelData.transform("lag")` | L.x, L2.x, ... | plm::lag; collapse::flag | Value k periods earlier on the calendar of the panel. |
| `PanelData.transform("lead")` | F.x, F2.x, ... | plm::lead; collapse::flag | Value k periods later on the calendar of the panel. |
| `PanelData.transform("quasi_demean")` | xtreg, re  (x - theta * unit mean) | — | Quasi-demeaning x − θ x̄_i with θ given by the user. |
| `PanelData.transform("rolling")` | rangestat / tsegen (rowmean(L(0/w-1).x)) | — | Mean, sum, standard deviation, min or max over the last w calendar periods. |
| `PanelData.transform("seasonal")` | S12.x | plm::diff (lag argument) | Difference with the value `period` periods earlier. |
| `PanelData.transform("twoway")` | — | linearmodels general_demean | Residual of x on unit and period effects (exact if balanced, iterated if not). |
| `PanelData.transform("within")` | xtdata, fe  (x - unit mean + overall mean with add_mean=True) | plm::Within; collapse::fwithin; linearmodels demean | Deviation from the mean of the unit. |

### Exploration tools

| Call | Stata | R or Python equivalent | What it does |
|---|---|---|---|
| `PanelData.audit()` | — | expdpy explore_panel_structure; pyxtdpthresh diagnose() | The whole panel examined in one call: structure, gaps, missing values, variation, outliers and equations, with alerts and advice in French and English. |
| `PanelData.balance_cost()` | — | Stata xtbalfront, xtbalance2 (SSC); plm::make.pbalanced applies one rule | Units and observations kept by each minimum length T_i >= k and each window. |
| `PanelData.event_time()` | — | — | Time relative to an event for each row, and the periods each unit is observed around it. |
| `PanelData.gap_types()` | tsreport, panel detail (lists the gaps, does not classify them) | Stata xtpattern (SSC), moderndid.scan_gaps | Late entry, early exit, internal gap and re-entry of each unit, gap lengths. |
| `PanelData.level_shifts()` | — | — | Level shifts of each unit (change of median) and units that weigh most in the within variation. |
| `PanelData.residual_variation()` | — | — | Observations, units and variance left by the within, first-difference and orthogonal-deviation transformations. |
| `xt.resolve_duplicates()` | — | Stata duplicates drop (identical rows only) | Merge repeated (unit, period) pairs by first, last, mean, sum... and report what was merged. |
| `PanelData.spells()` | — | Stata tsspell (SSC) | Spells of a discrete state: duration, changes per unit, units that change. |
| `PanelData.to_array3d()` | — | plm values3d, Stata psmat (one variable) | Array variable x period x unit over the whole calendar, with the mask of rows. |
| `xt.to_long()` | reshape long | panelr::long_panel | Back from the wide form: one row per unit and period, ready for xtset. |
| `PanelData.to_wide()` | reshape wide | panelr::widen_panel, plm reshape | One row per unit, one column per variable and observed period. |
| `PanelData.unbalancedness()` | — | plm::punbalancedness | Ahrens-Pincus measures of the imbalance of a panel (gamma and nu). |
| `PanelData.usable_equations()` | — | — | Equations left for a dynamic-panel estimator, by unit and period (FD or FOD). |
| `PanelData.variation()` | — | plm::pvar, linearmodels are_varying | For each variable: units and periods that vary, and whether it is constant, time-invariant or common to all units. |
| `PanelData.within_outliers()` | — | — | Within-unit outliers by robust scores (median and MAD in the unit), flagged and never removed. |

### Figures

| Call | Stata | R or Python equivalent | What it does |
|---|---|---|---|
| `PanelData.plot("attrition")` | — | — | Units observed by period with entries and exits, and Kaplan-Meier survival. |
| `PanelData.plot("balance_cost")` | — | — | Share of units and observations kept by each minimum length and window. |
| `PanelData.plot("distribution")` | — | — | Distribution of a variable between units and within units. |
| `PanelData.plot("lines")` | xtline | — | Trajectories of a sample of units, with the mean by period. |
| `PanelData.plot("participation")` | — | — | Map of the periods in which each unit is observed, sorted by pattern. |
| `PanelData.plot("scatter")` | — | — | Between (unit means) and within (deviations) scatter plots of x and y. |
| `PanelData.plot("states")` | — | — | Map of units by periods coloured by the state of a discrete variable (treatment status), missing cells distinct. |
| `PanelData.plot("transition")` | — | — | Heat map of the transition matrix between consecutive observations. |

### Statistical tests

| Call | Stata | R or Python equivalent | What it does |
|---|---|---|---|
| `xt.run_test("ar_test")` | xtabond2; estat abond | plm::mtest | Arellano-Bond m_k test of serial correlation of order k in differenced residuals. |
| `xt.run_test("baltagi_li")` | xttest1; mata | plm::pbsytest | Joint, marginal and robust LM tests of random effects and AR(1) errors. |
| `xt.run_test("breitung")` | xtunitroot breitung | — | Unit-root test of a common AR parameter on forward-orthogonalised data. |
| `xt.run_test("breusch_pagan")` | estat hettest | lmtest::bptest | Test that the variance of the errors depends on given variables. |
| `xt.run_test("cips")` | pescadf | plm::cipstest | Unit-root test averaging the cross-sectionally augmented Dickey-Fuller t statistics. |
| `xt.run_test("diff_hansen")` | xtabond2 | — | Difference of two Hansen statistics: exogeneity of a subset of instruments. |
| `xt.run_test("dumitrescu_hurlin")` | mata | plm::pgrangertest | Average of the unit Wald statistics of Granger non-causality. |
| `xt.run_test("f_effects")` | xtreg | plm::pFtest | F test of the fixed effects from the pooled and within residual sums of squares. |
| `xt.run_test("fisher")` | xtunitroot fisher | plm::purtest(test='madwu') | Combination of the p-values of one ADF test per unit. |
| `xt.run_test("frees")` | xtcsd | plm::pcdtest(test='rho') (related) | Test that the average squared Spearman correlation between units is 1 / (T - 1). |
| `xt.run_test("friedman")` | xtcsd | plm::pcdtest(test='rho') (related) | Test that the average Spearman correlation between units is zero. |
| `xt.run_test("hadri")` | xtunitroot hadri | plm::purtest(test='hadri') | LM test of the null that every series of the panel is stationary. |
| `xt.run_test("hansen")` | ivregress gmm; xtabond2 | — | J test of the overidentifying restrictions of an efficient GMM fit. |
| `xt.run_test("harris_tzavalis")` | xtunitroot ht | — | Unit-root test of a common AR parameter for panels with large N and fixed T. |
| `xt.run_test("hausman")` | hausman | plm::phtest | Hausman test comparing a consistent and an efficient estimator. |
| `xt.run_test("im_pesaran_shin")` | xtunitroot ips | plm::purtest(test='ips') | Unit-root test averaging the Dickey-Fuller t statistics of the units. |
| `xt.run_test("jarque_bera")` | summarize | tseries::jarque.bera.test | Asymptotic joint test of zero skewness and kurtosis three. |
| `xt.run_test("kao")` | xtcointtest kao | — | Residual-based tests of no cointegration with a common cointegrating vector. |
| `xt.run_test("levin_lin_chu")` | xtunitroot llc | plm::purtest(test='levinlin') | Unit-root test of a common AR parameter, with the bias-adjusted t statistic. |
| `xt.run_test("lm_cross_section")` | mata | plm::pcdtest(test='lm') | Sum of the squared correlations between units, weighted by the common periods. |
| `xt.run_test("lm_random_effects")` | xttest0 | plm::plmtest(type='bp') | Test that the variance of the unit effects is zero, from pooled residuals. |
| `xt.run_test("panel_durbin_watson")` | xtregar | plm::pbnftest | Modified Durbin-Watson and Baltagi-Wu LBI statistics of the within residuals. |
| `xt.run_test("pedroni")` | xtcointtest pedroni | — | Residual-based tests of no cointegration with panel-specific cointegrating vectors. |
| `xt.run_test("pesaran_cd")` | xtcsd; xtcd2 | plm::pcdtest(test='cd') | Test that the correlations between units are zero on average. |
| `xt.run_test("reset")` | estat ovtest; test | lmtest::resettest | F test that powers of the fitted values (or regressors) add explanatory power. |
| `xt.run_test("sargan")` | ivregress 2sls; xtabond2 | plm::sargan | Sargan test of the overidentifying restrictions of a one-step fit. |
| `xt.run_test("shapiro_francia")` | sfrancia | nortest::sf.test | Squared correlation between the data and their normal scores. |
| `xt.run_test("shapiro_wilk")` | swilk | stats::shapiro.test (Royston 1995: close, not identical) | Squared correlation between the ordered data and the normal order statistics. |
| `xt.run_test("skewness_kurtosis")` | sktest | moments::agostino.test / anscombe.test | Joint test of normality from the skewness and the kurtosis of the residuals. |
| `xt.run_test("slope_homogeneity")` | xthst | — | Test that the slope coefficients are the same in every unit. |
| `xt.run_test("wald")` | test | plm::pwaldtest | Joint Wald chi-square or F test of linear restrictions R b = r. |
| `xt.run_test("wald_groupwise")` | xttest3 | plm::pwtest (related); lmtest | Test that the error variance is the same in every unit of a fixed-effects fit. |
| `xt.run_test("westerlund")` | xtcointtest westerlund | — | Variance ratio of the cumulated residuals of the cointegrating regressions. |
| `xt.run_test("white")` | estat imtest | lmtest::bptest(studentize=TRUE) with squares and products | White's general test that the variance of the errors is constant. |
| `xt.run_test("wooldridge")` | xtserial | plm::pwfdtest | Test of first-order serial correlation of the errors from first differences. |

<!-- annex:a:end -->

## Annex B. Statistical test cards

<!-- annex:b:start -->

### Causality

#### `dumitrescu_hurlin`

Average of the unit Wald statistics of Granger non-causality.

```python
dumitrescu_hurlin(panel: PanelBase, dependent: str, cause: str, *, lags: int = 1,
    statistic: Literal['ztilde', 'zbar', 'wbar'] = 'ztilde')
```

- **Input**: panel (category A); balanced panel required
- **Family**: causality
- **Stata**: mata
- **R**: plm::pgrangertest
- **Reference**: Dumitrescu, E.-I. and Hurlin, C. (2012), Economic Modelling 29(4), 1450-1460; Lopez, L. and Weber, S. (2017), Stata Journal 17(4), 972-984
- **Validation**: 5 cases run in Stata 17 and frozen
- **Aliases**: `granger`, `xtgcause`, `pgrangertest`

### Cointegration

#### `kao`

Residual-based tests of no cointegration with a common cointegrating vector.

```python
kao(panel: PanelBase, dependent: str, regressors: Sequence[str], *, statistic:
    Literal['modified_df', 'df', 'adf', 'unadjusted_modified_df', 'unadjusted_df'] =
    'adf', adf_lags: int = 1, adf_lag_selection: Literal['aic', 'bic', 'hqic'] | None =
    None, kernel: Literal['bartlett', 'parzen', 'quadraticspectral'] = 'bartlett',
    kernel_lags: int | None = None, demean: bool = False)
```

- **Input**: panel (category A); balanced panel required
- **Family**: cointegration
- **Stata**: xtcointtest kao
- **R**: —
- **Reference**: Kao, C. (1999), J. of Econometrics 90(1), 1-44
- **Validation**: 26 cases run in Stata 17 and frozen

#### `pedroni`

Residual-based tests of no cointegration with panel-specific cointegrating vectors.

```python
pedroni(panel: PanelBase, dependent: str, regressors: Sequence[str], *, kernel_lags:
    int, ar: Literal['panelspecific', 'same'] = 'panelspecific', statistic:
    Literal['nu', 'rho', 'pp', 'adf'] = 'adf', deterministics: Literal['constant',
    'none', 'trend'] = 'constant', kernel: Literal['bartlett', 'parzen',
    'quadraticspectral'] = 'bartlett', adf_lags: int = 1, adf_lag_selection:
    Literal['aic', 'bic', 'hqic'] | None = None, demean: bool = False)
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: cointegration
- **Stata**: xtcointtest pedroni
- **R**: —
- **Reference**: Pedroni, P. (1999), Oxford Bulletin of Economics and Statistics 61, 653-670; Pedroni, P. (2004), Econometric Theory 20(3), 597-625
- **Validation**: 20 cases run in Stata 17 and frozen

#### `westerlund`

Variance ratio of the cumulated residuals of the cointegrating regressions.

```python
westerlund(panel: PanelBase, dependent: str, regressors: Sequence[str], *, alternative:
    Literal['some', 'all'] = 'some', deterministics: Literal['constant', 'trend'] =
    'constant', demean: bool = False)
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: cointegration
- **Stata**: xtcointtest westerlund
- **R**: —
- **Reference**: Westerlund, J. (2005), Econometric Reviews 24(3), 297-316
- **Validation**: 11 cases run in Stata 17 and frozen

### Crosssection

#### `frees`

Test that the average squared Spearman correlation between units is 1 / (T - 1).

```python
frees(panel: PanelBase, variable: str)
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: crosssection
- **Stata**: xtcsd
- **R**: plm::pcdtest(test='rho') (related)
- **Reference**: Frees, E. W. (1995), J. of Econometrics 69(2), 393-414
- **Validation**: 3 cases run in Stata 17 and frozen
- **Aliases**: `xtcsd_frees`

#### `friedman`

Test that the average Spearman correlation between units is zero.

```python
friedman(panel: PanelBase, variable: str, *, dof: Literal['stata', 'periods'] = 'stata')
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: crosssection
- **Stata**: xtcsd
- **R**: plm::pcdtest(test='rho') (related)
- **Reference**: Friedman, M. (1937), J. of the American Statistical Association 32(200), 675-701; Frees, E. W. (1995), J. of Econometrics 69(2), 393-414
- **Validation**: 4 cases run in Stata 17 and frozen
- **Aliases**: `xtcsd_friedman`

#### `lm_cross_section`

Sum of the squared correlations between units, weighted by the common periods.

```python
lm_cross_section(panel: PanelBase, variable: str, *, min_overlap: int = 3)
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: crosssection
- **Stata**: mata
- **R**: plm::pcdtest(test='lm')
- **Reference**: Breusch, T. S. and Pagan, A. R. (1980), Rev. Econ. Studies 47(1), 239-253
- **Validation**: 3 cases run in Stata 17 and frozen
- **Aliases**: `bp_lm`, `xttest2`

#### `pesaran_cd`

Test that the correlations between units are zero on average.

```python
pesaran_cd(panel: PanelBase, variable: str, *, min_overlap: int = 3)
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: crosssection
- **Stata**: xtcsd; xtcd2
- **R**: plm::pcdtest(test='cd')
- **Reference**: Pesaran, M. H. (2015), Econometric Reviews 34(6-10), 1089-1117
- **Validation**: 5 cases run in Stata 17 and frozen
- **Aliases**: `cd`

### Effects

#### `baltagi_li`

Joint, marginal and robust LM tests of random effects and AR(1) errors.

```python
baltagi_li(residuals: object, offsets: object, *, statistic: str = 'joint', time_index:
    object | None = None)
```

- **Input**: arrays (category B); no panel structure
- **Family**: effects
- **Stata**: xttest1; mata
- **R**: plm::pbsytest
- **Reference**: Baltagi, B. H. and Li, Q. (1991), Economics Letters 37(3), 277-280; Bera, A. K., Sosa-Escudero, W. and Yoon, M. (2001), J. of Econometrics 101(1), 123-137; Sosa-Escudero, W. and Bera, A. K. (2008), Stata Journal 8(1), 134-135
- **Validation**: 9 cases run in Stata 17 and frozen
- **Aliases**: `xttest1`, `pbsytest`

#### `f_effects`

F test of the fixed effects from the pooled and within residual sums of squares.

```python
f_effects(rss_pooled: float, rss_effects: float, n_obs: int, n_groups: int,
    n_regressors: int)
```

- **Input**: arrays (category B); no panel structure
- **Family**: effects
- **Stata**: xtreg
- **R**: plm::pFtest
- **Reference**: Baltagi, B. H. (2021), Econometric Analysis of Panel Data, 6th ed., ch. 4
- **Validation**: 3 cases run in Stata 17 and frozen

#### `lm_random_effects`

Test that the variance of the unit effects is zero, from pooled residuals.

```python
lm_random_effects(residuals: object, offsets: object, *, variance_u: float | None =
    None)
```

- **Input**: arrays (category B); no panel structure
- **Family**: effects
- **Stata**: xttest0
- **R**: plm::plmtest(type='bp')
- **Reference**: Breusch, T. S. and Pagan, A. R. (1980), Rev. Econ. Studies 47(1), 239-253
- **Validation**: 4 cases run in Stata 17 and frozen

### Heterogeneity

#### `slope_homogeneity`

Test that the slope coefficients are the same in every unit.

```python
slope_homogeneity(panel: PanelBase, dependent: str, regressors: Sequence[str], *,
    partial: Sequence[str] = (), constant: bool = True, ar: bool = False,
    crosssectional: Sequence[str] = (), cr_lags: int | Sequence[int] = 0, adjusted: bool
    = True, hac: bool = False, kernel: Kernel = 'bartlett', bandwidth: int | None =
    None, whitening: bool = False)
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: heterogeneity
- **Stata**: xthst
- **R**: —
- **Reference**: Pesaran, M. H. and Yamagata, T. (2008), Journal of Econometrics 142(1), 50-93; Swamy, P. A. V. B. (1970), Econometrica 38(2), 311-323
- **Validation**: 26 cases run in Stata 17 and frozen
- **Aliases**: `xthst`, `pesaran_yamagata`, `swamy`

### Heteroskedasticity

#### `breusch_pagan`

Test that the variance of the errors depends on given variables.

```python
breusch_pagan(residuals: object, variables: object, *, variant: Literal['normal', 'iid',
    'fstat'] = 'normal')
```

- **Input**: arrays (category B); no panel structure
- **Family**: heteroskedasticity
- **Stata**: estat hettest
- **R**: lmtest::bptest
- **Reference**: Breusch, T. S. and Pagan, A. R. (1979), Econometrica 47(5), 1287-1294; Cook, R. D. and Weisberg, S. (1983), Biometrika 70(1), 1-10
- **Validation**: 6 cases run in Stata 17 and frozen

#### `wald_groupwise`

Test that the error variance is the same in every unit of a fixed-effects fit.

```python
wald_groupwise(residuals: object, offsets: object)
```

- **Input**: arrays (category B); no panel structure
- **Family**: heteroskedasticity
- **Stata**: xttest3
- **R**: plm::pwtest (related); lmtest
- **Reference**: Greene, W. H. (2000), Econometric Analysis, 4th ed., p. 598
- **Validation**: 3 cases run in Stata 17 and frozen

#### `white`

White's general test that the variance of the errors is constant.

```python
white(residuals: object, regressors: object)
```

- **Input**: arrays (category B); no panel structure
- **Family**: heteroskedasticity
- **Stata**: estat imtest
- **R**: lmtest::bptest(studentize=TRUE) with squares and products
- **Reference**: White, H. (1980), Econometrica 48(4), 817-838
- **Validation**: 4 cases run in Stata 17 and frozen

### Normality

#### `jarque_bera`

Asymptotic joint test of zero skewness and kurtosis three.

```python
jarque_bera(residuals: object)
```

- **Input**: arrays (category B); no panel structure
- **Family**: normality
- **Stata**: summarize
- **R**: tseries::jarque.bera.test
- **Reference**: Jarque, C. M. and Bera, A. K. (1987), International Statistical Review 55(2)
- **Validation**: 3 cases run in Stata 17 and frozen

#### `shapiro_francia`

Squared correlation between the data and their normal scores.

```python
shapiro_francia(residuals: object, *, ties: bool = True, transform: Literal['log',
    'boxcox'] = 'log')
```

- **Input**: arrays (category B); no panel structure
- **Family**: normality
- **Stata**: sfrancia
- **R**: nortest::sf.test
- **Reference**: Shapiro, S. S. and Francia, R. S. (1972), JASA 67(337), 215-216; Royston, P. (1983), Applied Statistics 32(2), 121-133; Royston, P. (1993), Statistician 42(1), 37-43
- **Validation**: 8 cases run in Stata 17 and frozen
- **Aliases**: `sfrancia`

#### `shapiro_wilk`

Squared correlation between the ordered data and the normal order statistics.

```python
shapiro_wilk(residuals: object, *, ties: bool = True)
```

- **Input**: arrays (category B); no panel structure
- **Family**: normality
- **Stata**: swilk
- **R**: stats::shapiro.test (Royston 1995: close, not identical)
- **Reference**: Shapiro, S. S. and Wilk, M. B. (1965), Biometrika 52(3-4), 591-611; Royston, P. (1982), Applied Statistics 31(2), 115-124 and 31(2), 176-180
- **Validation**: 13 cases run in Stata 17 and frozen
- **Aliases**: `swilk`

#### `skewness_kurtosis`

Joint test of normality from the skewness and the kurtosis of the residuals.

```python
skewness_kurtosis(residuals: object, *, adjust: bool = True)
```

- **Input**: arrays (category B); no panel structure
- **Family**: normality
- **Stata**: sktest
- **R**: moments::agostino.test / anscombe.test
- **Reference**: D'Agostino, R. B., Belanger, A. and D'Agostino Jr., R. B. (1990), The American Statistician 44(4), 316-321; Royston, P. (1991), Stata Technical Bulletin 3
- **Validation**: 4 cases run in Stata 17 and frozen

### Postestimation

#### `ar_test`

Arellano-Bond m_k test of serial correlation of order k in differenced residuals.

```python
ar_test(order: int, test: Equations, *, weight: object, covariance: object, estimation:
    Equations | None = None, min_pairs: int = 1)
```

- **Input**: arrays (category B); no panel structure
- **Family**: postestimation
- **Stata**: xtabond2; estat abond
- **R**: plm::mtest
- **Reference**: Arellano, M. and Bond, S. (1991), Review of Economic Studies 58(2), 277-297
- **Validation**: 6 cases run in Stata 17 and frozen

#### `diff_hansen`

Difference of two Hansen statistics: exogeneity of a subset of instruments.

```python
diff_hansen(j_full: float, j_excluding: float, n_dropped: int, *, df_full: int | None =
    None)
```

- **Input**: arrays (category B); no panel structure
- **Family**: postestimation
- **Stata**: xtabond2
- **R**: —
- **Reference**: Eichenbaum, M. S., Hansen, L. P. and Singleton, K. J. (1988), QJE 103(1), 51-78
- **Validation**: 2 cases run in Stata 17 and frozen

#### `hansen`

J test of the overidentifying restrictions of an efficient GMM fit.

```python
hansen(moments: object, weight: object, n_parameters: int, *, n_instruments: int | None
    = None, n_groups: int | None = None)
```

- **Input**: arrays (category B); no panel structure
- **Family**: postestimation
- **Stata**: ivregress gmm; xtabond2
- **R**: —
- **Reference**: Hansen, L. P. (1982), Econometrica 50(4), 1029-1054
- **Validation**: 6 cases run in Stata 17 and frozen

#### `hausman`

Hausman test comparing a consistent and an efficient estimator.

```python
hausman(b_consistent: object, b_efficient: object, v_consistent: object, v_efficient:
    object, *, keep: Sequence[int] | None = None, names: Sequence[str] | None = None,
    df: int | None = None)
```

- **Input**: arrays (category B); no panel structure
- **Family**: postestimation
- **Stata**: hausman
- **R**: plm::phtest
- **Reference**: Hausman, J. A. (1978), Econometrica 46(6), 1251-1271
- **Validation**: 4 cases run in Stata 17 and frozen

#### `sargan`

Sargan test of the overidentifying restrictions of a one-step fit.

```python
sargan(moments: object, weight: object, n_parameters: int, *, scale: float,
    n_instruments: int | None = None, n_groups: int | None = None)
```

- **Input**: arrays (category B); no panel structure
- **Family**: postestimation
- **Stata**: ivregress 2sls; xtabond2
- **R**: plm::sargan
- **Reference**: Sargan, J. D. (1958), Econometrica 26(3), 393-415
- **Validation**: 7 cases run in Stata 17 and frozen

#### `wald`

Joint Wald chi-square or F test of linear restrictions R b = r.

```python
wald(b: object, v: object, restrictions: object | None = None, values: object | None =
    None, *, df_resid: float | None = None)
```

- **Input**: arrays (category B); no panel structure
- **Family**: postestimation
- **Stata**: test
- **R**: plm::pwaldtest
- **Reference**: Wald, A. (1943), Transactions of the American Mathematical Society 54(3)
- **Validation**: 4 cases run in Stata 17 and frozen

### Serial

#### `panel_durbin_watson`

Modified Durbin-Watson and Baltagi-Wu LBI statistics of the within residuals.

```python
panel_durbin_watson(panel: PanelBase, dependent: str, regressors: Sequence[str])
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: serial
- **Stata**: xtregar
- **R**: plm::pbnftest
- **Reference**: Bhargava, A., Franzini, L. and Narendranathan, W. (1982), Review of Economic Studies 49(4), 533-549; Baltagi, B. H. and Wu, P. X. (1999), Econometric Theory 15(6), 814-823
- **Validation**: 6 cases run in Stata 17 and frozen
- **Aliases**: `xtregar`, `lbi`, `bhargava`

#### `wooldridge`

Test of first-order serial correlation of the errors from first differences.

```python
wooldridge(panel: PanelBase, dependent: str, regressors: Sequence[str])
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: serial
- **Stata**: xtserial
- **R**: plm::pwfdtest
- **Reference**: Wooldridge, J. M. (2002), Econometric Analysis of Cross Section and Panel Data, MIT Press, sec. 10.4.1; Drukker, D. M. (2003), Stata Journal 3(2), 168-177
- **Validation**: 4 cases run in Stata 17 and frozen
- **Aliases**: `xtserial`

### Specification

#### `reset`

F test that powers of the fitted values (or regressors) add explanatory power.

```python
reset(y: object, x: object, *, rhs: bool = False)
```

- **Input**: arrays (category B); no panel structure
- **Family**: specification
- **Stata**: estat ovtest; test
- **R**: lmtest::resettest
- **Reference**: Ramsey, J. B. (1969), Journal of the Royal Statistical Society B 31(2), 350-371
- **Validation**: 7 cases run in Stata 17 and frozen
- **Aliases**: `ovtest`, `ramsey`

### Unitroot

#### `breitung`

Unit-root test of a common AR parameter on forward-orthogonalised data.

```python
breitung(panel: PanelBase, variable: str, *, deterministics: Literal['none', 'constant',
    'trend'] = 'constant', lags: int = 0, demean: bool = False, robust: bool = False)
```

- **Input**: panel (category A); balanced panel required
- **Family**: unitroot
- **Stata**: xtunitroot breitung
- **R**: —
- **Reference**: Breitung, J. (2000), Advances in Econometrics 15, 161-178; Breitung, J. and Das, S. (2005), Statistica Neerlandica 59(4), 414-433
- **Validation**: 9 cases run in Stata 17 and frozen

#### `cips`

Unit-root test averaging the cross-sectionally augmented Dickey-Fuller t statistics.

```python
cips(panel: PanelBase, variable: str, *, lags: int | Sequence[int], deterministics:
    Literal['constant', 'trend'] = 'constant', truncate: bool = True)
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: unitroot
- **Stata**: pescadf
- **R**: plm::cipstest
- **Reference**: Pesaran, M. H. (2007), Journal of Applied Econometrics 22(2), 265-312
- **Validation**: 13 cases run in Stata 17 and frozen
- **Aliases**: `pescadf`, `pesaran_cips`

#### `fisher`

Combination of the p-values of one ADF test per unit.

```python
fisher(panel: PanelBase, variable: str, *, lags: int, deterministics:
    Literal['constant', 'trend', 'drift'] = 'constant', demean: bool = False, method:
    Literal['inverse_chi2', 'inverse_normal', 'logit', 'modified'] = 'inverse_chi2',
    test: Literal['dfuller', 'pperron'] = 'dfuller')
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: unitroot
- **Stata**: xtunitroot fisher
- **R**: plm::purtest(test='madwu')
- **Reference**: Maddala, G. S. and Wu, S. (1999), Oxford Bull. Econ. Stat. 61(S1), 631-652; Choi, I. (2001), J. of International Money and Finance 20(2), 249-272
- **Validation**: 6 cases run in Stata 17 and frozen

#### `hadri`

LM test of the null that every series of the panel is stationary.

```python
hadri(panel: PanelBase, variable: str, *, deterministics: Literal['constant', 'trend'] =
    'constant', demean: bool = False, robust: bool = False, kernel: Literal['bartlett',
    'parzen', 'quadraticspectral'] | None = None, lags: int | float | None = None)
```

- **Input**: panel (category A); balanced panel required
- **Family**: unitroot
- **Stata**: xtunitroot hadri
- **R**: plm::purtest(test='hadri')
- **Reference**: Hadri, K. (2000), Econometrics Journal 3(2), 148-161
- **Validation**: 9 cases run in Stata 17 and frozen

#### `harris_tzavalis`

Unit-root test of a common AR parameter for panels with large N and fixed T.

```python
harris_tzavalis(panel: PanelBase, variable: str, *, deterministics: Literal['none',
    'constant', 'trend'] = 'constant', demean: bool = False, small_sample_t: bool =
    False)
```

- **Input**: panel (category A); balanced panel required
- **Family**: unitroot
- **Stata**: xtunitroot ht
- **R**: —
- **Reference**: Harris, R. D. F. and Tzavalis, E. (1999), J. of Econometrics 91(2), 201-226
- **Validation**: 6 cases run in Stata 17 and frozen
- **Aliases**: `ht`

#### `im_pesaran_shin`

Unit-root test averaging the Dickey-Fuller t statistics of the units.

```python
im_pesaran_shin(panel: PanelBase, variable: str, *, lags: int | None = None,
    lag_selection: Literal['aic', 'bic', 'hqic'] | None = None, deterministics:
    Literal['constant', 'trend'] = 'constant', demean: bool = False, interpolate: bool =
    True)
```

- **Input**: panel (category A); balanced or unbalanced
- **Family**: unitroot
- **Stata**: xtunitroot ips
- **R**: plm::purtest(test='ips')
- **Reference**: Im, K. S., Pesaran, M. H. and Shin, Y. (2003), J. of Econometrics 115(1), 53-74
- **Validation**: 12 cases run in Stata 17 and frozen
- **Aliases**: `ips`

#### `levin_lin_chu`

Unit-root test of a common AR parameter, with the bias-adjusted t statistic.

```python
levin_lin_chu(panel: PanelBase, variable: str, *, deterministics: Literal['none',
    'constant', 'trend'] = 'constant', lags: int = 1, lag_selection: Literal['aic',
    'bic', 'hqic'] | None = None, kernel: Literal['bartlett', 'parzen',
    'quadraticspectral'] = 'bartlett', bandwidth: int | float | Literal['llc', 'nwest']
    = 'llc', demean: bool = False, interpolate: bool = True)
```

- **Input**: panel (category A); balanced panel required
- **Family**: unitroot
- **Stata**: xtunitroot llc
- **R**: plm::purtest(test='levinlin')
- **Reference**: Levin, A., Lin, C.-F. and Chu, C.-S. J. (2002), J. of Econometrics 108(1), 1-24
- **Validation**: 12 cases run in Stata 17 and frozen
- **Aliases**: `llc`

<!-- annex:b:end -->

## Annex C. Return codes and errors

<!-- annex:c:start -->

### Stata return codes kept

| Code | Meaning |
|---:|---|
| 110 | a variable that must be new already exists |
| 134 | too many distinct values to tabulate |
| 198 | invalid syntax, option or argument |
| 451 | repeated time values within a panel |
| 459 | something that should not be missing is missing (or too few distinct values) |
| 498 | the panel is not balanced or has gaps, for a command that needs neither |
| 2000 | no observations |

### Error classes

| Class | Raised when |
|---|---|
| `DataError` | The data do not satisfy a requirement (balance, duplicates, missing values...). |
| `DeclarationError` | The panel cannot be declared (``xtset``): repeated pairs, missing keys, bad time. |
| `NumericalError` | A numerical step failed (singular matrix, no convergence, empty sample). |
| `PyxtpanelError` | Base class of every error raised by pyxtpanel. |
| `PyxtpanelWarning` | Non-fatal condition worth reporting (infinite value read as missing, fallback used). |
| `SpecificationError` | An option or an argument is invalid or inconsistent. |
| `TimeAxisError` | The time variable does not define a regular axis (frequency, step, irregular dates). |

<!-- annex:c:end -->

## Annex D. Behaviours of Stata reproduced on purpose

<!-- annex:d:start -->

Each item is a place where the package follows what Stata does rather than the textbook, and where a test compares the result with Stata 17.

**`pyxtpanel.describe.patterns`**

- xtdescribe sizes a column from the smallest observed step, not from delta

**`pyxtpanel.describe.summary`**

- xtsum within uses the divisor N - 1 (not N - n); min, max add back the mean

**`pyxtpanel.describe.tab`**

- the total within percent weights the values by their between frequencies

**`pyxtpanel.describe.trans`**

- xttrans pairs a row with the unit's next row, so calendar gaps are bridged
- xttrans stores the variable as float (single precision) before tabulating

**`pyxtpanel.diagnostics.cointegration.pedroni`**

- T is the longest unit on an unbalanced panel; every p-value is reported as 1 - normal(|z|), the left tail normal(z) (right tail for nu) is used here; the table of adjustment terms has +25.889 (not -25.889) for the panel rho, no constant, 7 regressors.

**`pyxtpanel.diagnostics.cointegration.westerlund`**

- T is the largest unit length on an unbalanced panel (not the mean), and the reported p-value is 1 - normal(|z|); the left tail normal(z) is used here.

**`pyxtpanel.diagnostics.crosssection.frees`**

- that p-value doubles the normal tail although Q is only large under the alternative, and exceeds 1 for a negative Q; it is reproduced, capped at 1.

**`pyxtpanel.diagnostics.crosssection.pesaran_cd`**

- on an unbalanced panel xtcd2 does not use this formula (it fills the panel and correlates the unit means with all non-missing values); its CD differs from xtcsd's, e.g. 16.08 against 2.52 on the 'gaps' panel.  This function follows Pesaran (2015) and xtcsd.

**`pyxtpanel.diagnostics.effects.baltagi_li`**

- xttest1 stores the lagged residual as a float, so its statistics agree with the double-precision ones to about 2e-6 (relative), and returns the one-sided p-values as 1 - normprob(z), which loses digits from z = 5 and is 0 beyond z = 8; these are exact here

**`pyxtpanel.diagnostics.heterogeneity._hac`**

- the singular values are capped at 0.97 and the matrix is rebuilt as U diag(s) V, not U diag(s) V' (the third factor of Mata's svd() is already V'), so the coefficient matrix of a unit with several regressors is not the one that was estimated; the result then depends on the signs of the singular vectors, which Mata and SciPy give alike with one or two regressors and differently for about a quarter of the units with three
- the ratio is squared once more than in Andrews' alpha(2)

**`pyxtpanel.diagnostics.postestimation.ar_test`**

- on a panel with gaps xtabond and xtabond2 give different AR statistics (here AR(1) -3.049 and -2.963); the equations handed to this function decide which one is reproduced

**`pyxtpanel.diagnostics.unitroot.fisher`**

- a unit with equations but no defined statistic (saturated regression) gets p = 1 and counts in the degrees of freedom; Z and L drop it (invnormal(1) is missing).  A unit without any equation is skipped, but Stata stops when that unit is the first one.

**`pyxtpanel.diagnostics.unitroot.levin_lin_chu`**

- with lag selection and no constant, Stata passes a mistyped macro to its lag selection, so the information criteria are computed with a constant in the regressions.

<!-- annex:d:end -->

## Annex E. Validation journal

<!-- annex:e:start -->

| Domain | Stata | Frozen on | Cases |
|---|---|---|---:|
| `campaign` | Stata 17 MP (batch, PowerShell) | 2026-11-18 | 1 |
| `diagnostics_causality` | Stata 17 MP (batch, PowerShell) | 2026-11-17 | 5 |
| `diagnostics_cointegration` | Stata 17 MP (batch, PowerShell) | 2026-11-17 | 60 |
| `diagnostics_crosssection` | Stata 17 MP (batch, PowerShell) | 2026-11-17 | 18 |
| `diagnostics_effects` | Stata 17 MP (batch, PowerShell) | 2026-11-17 | 16 |
| `diagnostics_heterogeneity` | Stata 17 MP (batch, PowerShell) | 2026-11-17 | 26 |
| `diagnostics_heteroskedasticity` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 13 |
| `diagnostics_normality` | Stata 17 MP (batch, PowerShell) | 2026-11-15 | 28 |
| `diagnostics_postestimation` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 29 |
| `diagnostics_serial` | Stata 17 MP (batch, PowerShell) | 2026-11-15 | 10 |
| `diagnostics_specification` | Stata 17 MP (batch, PowerShell) | 2026-11-15 | 11 |
| `diagnostics_unitroot` | Stata 17 MP (batch, PowerShell) | 2026-11-17 | 80 |
| `explore` | Stata 17 MP (batch, PowerShell) | 2026-11-17 | 56 |
| `operators` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 12 |
| `transforms` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 38 |
| `tsfill` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 14 |
| `xtcorr` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 8 |
| `xtdescribe` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 58 |
| `xtset` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 16 |
| `xtsum` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 30 |
| `xttab` | Stata 17 MP (batch, PowerShell) | 2026-11-07 | 40 |
| **Total** | | | **569** |

<!-- annex:e:end -->
