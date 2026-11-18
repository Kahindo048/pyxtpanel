# Performance of the panel foundation

Measured on one Windows machine (Python 3.12, pandas 3.1, numpy 2.5, Stata 17 MP), on simulated panels
of 100 000 and 1 000 000 rows: 50 periods per unit, 5 % of the rows removed at random, three numeric
columns, **rows shuffled** (the worst case for a declaration that has to sort).  Times are the best of
three runs from a freshly declared panel.  The tables are generated, so they can be refreshed after a
change:

```
python -m tools.perf --rows 100000 1000000 --memory --pandas --write docs/performance_tables.md
python -m validation.perf_stata --write docs/performance_stata.md      # needs Stata, through PowerShell
python -m tools.perf --profile xtset tsfill --rows 1000000             # cProfile of chosen operations
python -m tools.perf --baseline tests/benchmarks/baseline.json         # reference of the speed guard
```

## Targets (I.3.3)

| Target | Measured at 1 000 000 rows | Met |
|---|---|---|
| `xtset` plus one operator in under one second | `xtset` 0.27 s + `lag` 0.10 s = 0.37 s (0.60 s before the work below) | yes |
| Every listed operation under one second | slowest: `within` two-way 0.55 s, `xtcorr` 0.48 s | yes |
| Peak allocation of an operation | at most 111 MB (two-way `within`, `tsfill`), the frame itself being about 48 MB | yes |

No corrective plan is needed.

## Where the time went (I.1.1, I.1.2)

`cProfile` of each operation at 1 000 000 rows, before the optimisations.  The choice column is the one
of I.1.2: an algorithm, numpy, or an optional accelerator.

| Operation | Before (s) | Bottleneck found by the profile | Choice | After (s) |
|---|---:|---|---|---:|
| `xtset` | 0.435 | the (unit, period) sort ran twice (duplicate check, then index) and `np.unique` sorted a third time to count rows per unit | algorithm: sort once and share the order; count with `bincount` | 0.27 |
| `lag`, `diff`, `L2D.x` | 0.158 / 0.168 / 0.247 | putting the result back on the frame's row order went through a `DataFrame.reindex` (0.1 s) | numpy: scatter into a missing-filled array | 0.10 / 0.11 / 0.20 |
| `within`, `between`, `group_mean` | 0.169 / 0.151 / 0.105 | the same reindex | numpy: same scatter | 0.13 / 0.11 / 0.05 |
| `xtcorr` | 0.721 | the common-row mask and the fancy indexing were redone for each of the three kinds | algorithm: one mask per pair, shared by the kinds | 0.48 |
| `tsfill` | 0.422 | declares the filled frame again (inherits the `xtset` gain) | same as `xtset` | 0.34 |
| `within` two-way | 0.506 | alternating projections (`argsort` and `unique` to compact the effects, 0.15 s per sweep) | none: intrinsic to the method, under target | 0.55 (same, noise) |
| `xtdescribe` | 0.276 | `np.unique` on the participation patterns (12 sorts) | none: under target | 0.29 |
| `xtsum`, `gaps`, `gap_types`, `to_wide`, `arrays` | at most 0.09 | nothing above the noise | none | unchanged |

The operators already used a sorted-key search (`searchsorted` on `unit * span + rank`) and the per-unit
aggregates already used `bincount` (I.2.2): the profile shows no `groupby`, no Python loop over units
and no per-row Python call in any listed operation.

## Memory (I.2.3)

The index holds three integer arrays of one value per row (unit, rank, source row; the unit labels
and the offsets are one value per unit), 8 bytes each, so 24 MB at one million rows; one float column
is 8 MB.  The peaks of the generated table are the allocation of the operation on top of the declared
panel, measured with `tracemalloc`: from 16 MB (`gaps`) to 111 MB (two-way `within`, `tsfill`), that is
a few copies of a column (sort buffers, masks, the result).  No operation approaches the memory of the
frame times ten, so no change was needed; `with_frame` shares the sorted index when only the values
change.

## Optional accelerator (I.2.4)

**Numba is not adopted.**  The decision rule was to retain it only if an operation stayed above its
target after the algorithmic work; none did.  Every operation is below 0.6 s at one million rows and
the heaviest are sort-bound, where compiled loops would not help.  Only the Python 3.13 installation of the
development machine has `numba`, so an extra would also be hard to test.  The refusal is reconsidered if a future
operation (for example a bootstrap over units) cannot reach its target with numpy.

## Against pandas (I.3.2)

Same result, same inputs, pandas written as an expert would (`groupby().shift`, `groupby().transform`,
`describe`, `reindex` on the full product).  Where pandas has a one-line equivalent it is as fast or a
little faster for the bare computation; the package returns more (calendar lags that respect gaps, a
result aligned on the user's rows, the declared structure).  The `lag` pandas line is not the same
computation: `shift` moves by row, so the comparison adds the `diff == 1` test that makes it a calendar
lag.

Tables: [performance_tables.md](performance_tables.md) (columns "pandas (s)").

## Against Stata (I.3.1)

Stata 17 MP with its own `timer`, on the same panel written to a `.dta` (`xtset id t`,
`gen double lx = L.x`, `gen double dx = D.x`, within and between deviations by `bysort id: egen`,
`quietly xtsum`, `quietly xtdescribe`, `tsfill`).

Table: [performance_stata.md](performance_stata.md).  At one million rows the package is 4 times faster
than Stata on `xtset`, about 6 times on within and between, more than 11 times on `xtsum`, 6 times on `tsfill`
and 48 times on `xtdescribe` (the participation patterns).  The operators `L.` and `D.` are the one case
where Stata is faster (0.08 s against 0.11 s): Stata writes the column in place, the package also
puts the result back on the row order of the user's frame.

## Guards (I.4)

`pytest -m bench` (left out of the default run: a timing needs a quiet machine) runs two suites:

* `tests/benchmarks/test_speed.py`: each operation on a 200 000-row panel with `pytest-benchmark`
  (best of five rounds), divided by a fixed numpy workload measured in the same process
  (`tools.perf.calibration`); the ratio must stay within three times the one recorded in
  `tests/benchmarks/baseline.json`.  A slower machine moves both numbers, so the ratio travels.
* `tests/benchmarks/test_import_time.py`: `import pyxtpanel` loads neither numpy, pandas, scipy nor
  matplotlib and takes under half the time of `import numpy, pandas`; declaring and describing a
  panel (`xtset`, `xtsum`, `lag`, `within`, `audit`) loads neither matplotlib nor scipy; first use adds
  at most 60 % to the import time of numpy and pandas.

Scipy used to be loaded at the first use of `xt.xtset` (about one second), because the API module
imported the statistical tests at the top.  The tests are now reached through `api/_lazy.py`, which
imports them at the first call of a test.
