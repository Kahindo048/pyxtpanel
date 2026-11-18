# pyxtpanel

Panel-data management and exploration for Python, modelled on the `xt` layer of Stata (`xtset`, `xtdescribe`, `xtsum`, `xttab`, `xttrans`, `tsfill`, time-series operators), with new exploration tools and the statistical tests of panel data that do not depend on an estimator. It estimates no model: it is the foundation under the estimation packages `pyxtabond2`, `pyxthreg` and `pyxtdpthresh`.

**Status: under construction (0.0.1).** Nothing is released yet. The design is in [ARCHITECTURE.md](../ARCHITECTURE.md) and the work plan in [PLAN.md](../PLAN.md) (both in the parent folder while the package is being built).

## Design in one paragraph

A panel is declared once with `xtset(df, id, time, delta=None)` and becomes an immutable `PanelData` object that wraps the DataFrame, validates it (duplicates, missing keys, time axis) and caches a sorted index. Every function reads that object and returns new objects: the user's DataFrame is never modified. Results are typed objects with `.to_frame()`; messages are in French and English. Dependencies: `numpy`, `pandas`, `scipy`; plotting is the optional extra `pyxtpanel[plot]`.

## Documentation

One reference document, [docs/documentation.md](docs/documentation.md): getting started, concepts, every function with its Stata equivalent, the guide to the statistical tests, the contract for estimation packages, the method and the validation journal. Its Python examples are run by the test suite. Performance: [docs/performance.md](docs/performance.md); size and power of the tests: [docs/size_power.md](docs/size_power.md).

## Validation

Every function that has a Stata equivalent is compared with Stata 17 MP on shared datasets; the frozen reference results run in the test suite without Stata. See `validation/`.

## Licence

MIT, see [LICENSE](LICENSE).
