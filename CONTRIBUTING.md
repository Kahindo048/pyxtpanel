# Contributing

The design is in `ARCHITECTURE.md` (structure, extension points, dependency rules, coding rules) and the work plan in `PLAN.md`. The essentials:

## Set up

```bash
py -3.12 -m venv .venv && .venv\Scripts\activate
pip install -e ".[dev]"
pre-commit install
```

## Before every commit

```bash
ruff check src tests && ruff format src tests
mypy
lint-imports
pytest -m "not stata"
python -m tools.doc --write   # after adding a function, operator, figure, tool or test: refresh the annexes
python -m tools.docstrings    # every public docstring: summary, each parameter, Returns, Raises
```

A new feature must also be named in its part of `docs/documentation.md` (the check `python -m tools.doc` and the contract test fail otherwise); the Python blocks of that document are run by the test suite. `python -m tools.doc --pdf` rebuilds `docs/documentation.tex` and `docs/documentation.pdf` (generated files, ignored by git; the Markdown stays the source).

`pre-commit install` runs `ruff`, `mypy` and `lint-imports` before each commit (`.pre-commit-config.yaml`; the tools of the active environment are used).

## Continuous integration and releases

GitHub Actions (`.github/workflows/ci.yml`) runs on every push to `main` and every pull request: ruff, mypy and import-linter, the generated files, the docstrings of the public API, the whole test suite on Linux and Windows with Python 3.10 to 3.13, with pandas 2.x and 3.x and with the oldest numpy, pandas and scipy the package promises (without matplotlib), the coverage gate (90 %), and the build of the wheel and sdist installed alone in a fresh environment. No Stata is needed: its results are frozen in `tests/golden/`. A release is a tag: see [docs/RELEASING.md](docs/RELEASING.md).

## Adding a feature (one file, no existing file modified)

```bash
python tools/nouveau.py test --famille unitroot --nom fisher --categorie A
```

The generator creates the module, its unit test, its Stata validation case and a changelog line. The meta-tests in `tests/contract/` fail until the new item has its reference, `# Stata:` provenance, French and English messages, test and validation case (or the explicit note "no Stata reference: validated by simulation"). The generator comes with A.5.9; until then follow the layout of the existing modules.

## Rules in one screen

1. **Layers import downward only**, and sibling features do not import each other (`lint-imports`, and the contract tests).
2. **Exactness first**: `float64`; tolerances in one module with their origin.
3. **Provenance**: `# Stata: command (detail)`; deliberate quirks are marked `# STATA-QUIRK: reason`.
4. **No Numba, no compiled extension**: vectorised NumPy.
5. **Randomness** only through an explicit `numpy.random.Generator`.
6. **Errors** derive from `PyxtpanelError`; no `print` outside the explicit rendering of a result.
7. **Messages** live in `_i18n/` in French and English, with domain-prefixed keys.
8. **Results are immutable**, and no function modifies the user's DataFrame.
9. Commits follow the conventional format: `feat:`, `fix:`, `perf:`, `test:`, `docs:`, `refactor:`, `chore:`.

## Validation against Stata

`validation/run_stata.ps1` runs Stata 17 MP in batch **from PowerShell only**. Reference results are frozen in `tests/golden/`, so the test suite runs without Stata (marker `stata` is for local runs).
