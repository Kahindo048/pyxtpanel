# pyxtpanel

**Panel-data management, exploration and statistical tests for Python, validated against Stata.**

`pyxtpanel` brings the `xt` layer of Stata to Python: declare a panel once (`xtset`), describe it
(`xtdescribe`, `xtsum`, `xttab`, `xttrans`), take lags, differences and within transformations on the
calendar, and run the statistical tests of panel data that need no estimator. It estimates no model:
it is the foundation under the estimation packages
[`pyxtabond2`](https://pypi.org/project/pyxtabond2/), [`pyxthreg`](https://pypi.org/project/pyxthreg/)
and [`pyxtdpthresh`](https://pypi.org/project/pyxtdpthresh/), and a tool of its own for anyone who
holds panel data in pandas.

```text
pip install pyxtpanel            # numpy, pandas, scipy
pip install "pyxtpanel[plot]"    # adds matplotlib for the figures
```

Python 3.10 or later. Alpha (0.x): the API may change until 1.0.

## A first panel

```python
import numpy as np
import pandas as pd
import pyxtpanel as xt

rng = np.random.default_rng(0)
df = pd.DataFrame(
    {
        "country": np.repeat([f"C{i:02d}" for i in range(20)], 15),
        "year": np.tile(np.arange(2000, 2015), 20),
        "gdp": rng.normal(8, 1, size=300),
    }
)
df = df.drop(index=[5, 6, 47, 130])  # real data have holes
panel = xt.xtset(df, "country", "year")  # xtset country year
print(panel.xtsum("gdp").summary())  # xtsum gdp: overall, between, within
print(panel.operate("L.gdp", "D.gdp").head(3))  # L.gdp D.gdp, on the calendar
print(panel.audit().summary())  # what to know before modelling
print(panel.cross_section_dependence("gdp", test="pesaran_cd").pvalue)  # xtcsd
```

The panel is an immutable object: your DataFrame is never modified, every function returns a new
object, and every result has `.to_frame()`, `.summary()`, `.to_markdown()` and `.to_latex()`.
Messages and texts are in English or French (`xt.set_options(lang="fr")`).

## What is inside

| You want to | Look at |
|---|---|
| Declare and manage a panel: dates, `delta`, holes, `tsfill`, selection, `.dta` input and output | `xtset`, `PanelData.tsfill`, `xt.read_stata` |
| Lag, lead, difference, seasonal difference, with Stata's notation (`"L2D.x"`, `"L(1/3).x"`) | `PanelData.operate`, `lag`, `diff` |
| Remove unit means, take between means, forward orthogonal deviations, rolling windows | `within`, `between`, `fod`, `rolling` |
| Describe a panel as Stata does | `xtdescribe`, `xtsum`, `xttab`, `xttrans`, `xtcorr` |
| Understand an unbalanced panel (tools that Stata and R do not have) | `gap_types`, `balance_cost`, `unbalancedness`, `usable_equations`, `spells`, `audit` |
| Test unit roots, cross-sectional dependence, serial correlation, cointegration, causality, slope heterogeneity | 35 tests: `xt.list_tests()`, `PanelData.unit_root`, `xt.run_test` |
| Test the residuals of any estimator from bare arrays | `hansen`, `sargan`, `hausman`, `wald`, `ar_test`, `breusch_pagan`, `white`, ... |
| Draw the panel (extra `plot`) | `PanelData.plot("lines")`, `"states"`, `"attrition"`, ... |
| Feed an estimation package | `PanelData.arrays`, the frozen array contract |

## Why trust it

Every function that has a Stata equivalent is compared with **Stata 17 MP** on shared datasets, and
the reference results are frozen in the test suite, so the 3 500 tests run without Stata. The
statistical tests agree with Stata to 1e-6 or better, in most cases to 1e-9. A random campaign draws 160 panels (holes, entry
and exit, any identifiers, any time origin, Stata date formats) and compares 73 quantities per panel
with Stata; fresh runs of 3 700 more panels agree too. Where Stata does something odd, the package
either reproduces it on purpose and says so, or documents the difference
([documentation](https://github.com/Kahindo048/pyxtpanel/blob/main/docs/documentation.md), annex D and
part 9).

## Learn more

* [Documentation](https://github.com/Kahindo048/pyxtpanel/blob/main/docs/documentation.md): one
  reference document, from a first panel to the formulas; every Python example in it is run by the
  test suite.
* [Examples](https://github.com/Kahindo048/pyxtpanel/tree/main/examples): four scripts, a country-by-year
  panel, the Hansen (1999) firm investment data, an unbalanced panel, and the residual tests of an
  estimator.
* [CHANGELOG](https://github.com/Kahindo048/pyxtpanel/blob/main/CHANGELOG.md),
  [CONTRIBUTING](https://github.com/Kahindo048/pyxtpanel/blob/main/CONTRIBUTING.md) and
  [CITATION.cff](https://github.com/Kahindo048/pyxtpanel/blob/main/CITATION.cff).

## Licence

MIT, see [LICENSE](https://github.com/Kahindo048/pyxtpanel/blob/main/LICENSE). Where the formulas, the
embedded numerical tables and the reused code come from, and what is not copied (no Stata code), is in
[NOTICE.md](https://github.com/Kahindo048/pyxtpanel/blob/main/NOTICE.md). Stata is a registered
trademark of StataCorp LLC; this package is not affiliated with it.
