# Examples

Four scripts you can run as they are (`python examples/01_country_year.py`). They need only
`pyxtpanel` and its three dependencies; the test suite runs each of them, so they stay true.

| Script | What it shows | Data |
|---|---|---|
| [01_country_year.py](01_country_year.py) | The first hour with a panel: declare it, describe it, lags and differences on the calendar, within, the audit, one test, Markdown output | Simulated, 24 countries, 1990 to 2014, with holes |
| [02_hansen_investment.py](02_hansen_investment.py) | A real balanced panel: where the variation lives (between or within firms), unit root, cross-sectional dependence and serial correlation tests, LaTeX table | Hansen (1999), 565 firms, 1973 to 1987 (`data/hansen_1999_invest.dta`) |
| [03_unbalanced_panel.py](03_unbalanced_panel.py) | Why a panel is unbalanced (late entry, early exit, holes), what balancing costs, how many equations a dynamic-panel estimator keeps with first differences and with forward orthogonal deviations | Simulated, 40 units with the three causes mixed |
| [04_after_an_estimator.py](04_after_an_estimator.py) | The two directions with an estimator: the lag, the orthogonal deviations and the array contract it receives; the panel tests it runs on its residuals | Hansen (1999) |

The documentation ([docs/documentation.md](../docs/documentation.md)) explains every function used here.

**Data.** `data/hansen_1999_invest.dta` is the firm investment panel of Hansen (B. E., 1999,
*Threshold effects in non-dynamic panels*, Journal of Econometrics 93(2)); see
[NOTICE.md](../NOTICE.md). The unbalanced real-world panel `nlswork` is a Stata sample dataset and is
not redistributed here; example 03 builds an unbalanced panel with known causes instead, so its answers
can be checked.
