"""The documentation stays true (J.2): its examples run and its generated parts are current.

* every ``python`` block of parts 0 to 8 of ``docs/documentation.md`` is executed, in order, in one
  namespace (a block marked ``python skip`` is shown but not run);
* the annexes equal what ``python -m tools.doc --write`` would write;
* every registered operator, description, exploration tool, figure and test is named in the part
  that teaches it.
"""

from __future__ import annotations

import importlib.util
import re
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from tools import doc

import pyxtpanel as xt
from pyxtpanel._errors import DataError
from pyxtpanel.diagnostics._core import TESTS

TEXT = doc.TARGET.read_text(encoding="utf-8")


def test_the_annexes_are_up_to_date() -> None:
    assert doc.rendered(TEXT) == TEXT, "run: python -m tools.doc --write"


def test_every_registered_feature_is_named_in_its_part() -> None:
    assert doc.coverage(TEXT) == []


def test_the_document_has_ten_parts_and_five_annexes() -> None:
    assert sorted(doc.parts(TEXT)) == list(range(10))
    for key in doc.ANNEXES:
        assert doc.MARK.format(key=key, edge="start") in TEXT


def test_the_python_blocks_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    has_matplotlib = importlib.util.find_spec("matplotlib") is not None
    namespace: dict[str, object] = {"__name__": "documentation"}
    blocks = doc.python_blocks(TEXT)
    assert len(blocks) >= 12
    for number, (part, code) in enumerate(blocks, start=1):
        if part == 6 and not has_matplotlib:
            continue  # figures need the optional extra
        try:
            exec(compile(code, f"documentation part {part}, block {number}", "exec"), namespace)  # noqa: S102
        except Exception:
            pytest.fail(f"part {part}, block {number} failed:\n{code}\n{traceback.format_exc()}")


# ---------------------------------------------------------------- the balance table of part 7

CALLS: dict[str, tuple[tuple[object, ...], dict[str, object]]] = {
    "dumitrescu_hurlin": (("y", "x"), {}),
    "kao": (("y", ["x"]), {}),
    "pedroni": (("y", ["x"]), {"kernel_lags": 2}),
    "westerlund": (("y", ["x"]), {}),
    "frees": (("x",), {}),
    "friedman": (("x",), {}),
    "lm_cross_section": (("x",), {}),
    "pesaran_cd": (("x",), {}),
    "slope_homogeneity": (("y", ["x"]), {}),
    "panel_durbin_watson": (("y", ["x"]), {}),
    "wooldridge": (("y", ["x"]), {}),
    "breitung": (("x",), {}),
    "cips": (("x",), {"lags": 1}),
    "fisher": (("x",), {"lags": 1}),
    "hadri": (("x",), {}),
    "harris_tzavalis": (("x",), {}),
    "im_pesaran_shin": (("x",), {"lags": 1}),
    "levin_lin_chu": (("x",), {}),
}


def shaped_panel(kind: str) -> xt.PanelData:
    """12 units, 30 periods; one unit in three enters late or loses two periods."""
    rng = np.random.default_rng(3)
    rows = []
    for i in range(12):
        x = np.cumsum(rng.normal(size=30))
        y = 0.5 * x + 0.3 * np.cumsum(rng.normal(size=30))
        for t in range(30):
            if kind == "entry" and i % 3 == 0 and t < 4:
                continue
            if kind == "holes" and i % 3 == 0 and t in (10, 11):
                continue
            rows.append((i, t, x[t], y[t]))
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "x", "y"]), "id", "t")


def outcome(name: str, kind: str) -> str:
    args, options = CALLS[name]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            TESTS.get(name)(shaped_panel(kind), *args, **options)  # type: ignore[operator]
        except DataError as error:
            return f"refused ({error.stata_code})"
    return "yes"


def test_the_balance_table_of_part_7_is_what_the_tests_do() -> None:
    body = TEXT.split("<!-- balance-table:start -->")[1].split("<!-- balance-table:end -->")[0]
    rows = [line for line in body.splitlines() if line.startswith("| `")]
    documented = {
        re.match(r"\| `(\w+)`", line).group(1): [c.strip() for c in line.split("|")[2:5]]  # type: ignore[union-attr]
        for line in rows
    }
    panel_tests = {n for n, s in TESTS.specs() if s.category == "A"}
    assert set(documented) == panel_tests, "every category A test is in the table"
    for name, cells in documented.items():
        found = [outcome(name, kind) for kind in ("balanced", "entry", "holes")]
        assert cells == found, f"{name}: documented {cells}, the test does {found}"


# ---------------------------------------------------------------- the LaTeX version


def test_the_latex_preparation_keeps_tables_and_fixes_what_pandoc_cannot_read() -> None:
    prepared = doc.latex_markdown(TEXT)
    assert not prepared.startswith("# pyxtpanel documentation")
    assert "```python skip" not in prepared and "](#0-getting-started)" not in prepared
    assert "](#getting-started)" in prepared
    before = [line.count("|") for line in TEXT.splitlines() if line.startswith("|")]
    after = [line.count("|") for line in prepared.splitlines() if line.startswith("|")]
    assert before == after, "the widths of the columns change the separator rows only"
    separators = [line for line in prepared.splitlines() if re.fullmatch(r"\|[-:| ]+\|", line)]
    assert separators and all(max(map(len, s.strip("|").split("|"))) >= 6 for s in separators)
