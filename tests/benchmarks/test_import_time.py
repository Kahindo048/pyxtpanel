"""Import time and what a first use loads (I.4.2).

Each measure runs a fresh interpreter, so that nothing is cached.  ``import pyxtpanel`` only
registers names (the rest loads at first use); declaring and describing a panel must not load
matplotlib (an optional extra) nor ``scipy`` (only the statistical tests need it).
"""

from __future__ import annotations

import subprocess
import sys
import time

import pytest

pytestmark = pytest.mark.bench

SAMPLE = (
    "import pandas as pd, pyxtpanel as xt;"
    "f = pd.DataFrame({'id': [1, 1, 2, 2], 't': [1, 2, 1, 2], 'x': [1., 2., 3., 5.]});"
    "p = xt.xtset(f, 'id', 't'); p.xtsum('x'); p.lag('x'); p.within('x'); p.audit();"
)


def seconds(code: str, repeat: int = 3) -> float:
    best = float("inf")
    for _ in range(repeat):
        start = time.perf_counter()
        subprocess.run([sys.executable, "-c", code], check=True)
        best = min(best, time.perf_counter() - start)
    return best


def loaded(code: str, *modules: str) -> list[str]:
    probe = code + f"import sys; print([m for m in {modules!r} if m in sys.modules])"
    out = subprocess.run([sys.executable, "-c", probe], check=True, capture_output=True, text=True)
    return eval(out.stdout.strip().splitlines()[-1])  # noqa: S307 - our own list of names


def test_importing_the_package_loads_none_of_the_heavy_libraries() -> None:
    assert loaded("import pyxtpanel;", "numpy", "pandas", "scipy", "matplotlib") == []


def test_describing_a_panel_loads_neither_matplotlib_nor_scipy() -> None:
    assert loaded(SAMPLE, "matplotlib", "scipy") == []


def test_the_package_adds_little_to_the_import_of_its_dependencies() -> None:
    base = seconds("import numpy, pandas")
    package = seconds("import pandas, pyxtpanel as xt; xt.xtset")
    assert package <= 1.6 * base, f"{package:.2f} s against {base:.2f} s for numpy and pandas"


def test_the_import_of_the_package_itself_is_instant() -> None:
    base = seconds("import numpy, pandas")
    assert seconds("import pyxtpanel") <= 0.5 * base
