"""The README and the examples stay true (J.2.1, J.2.5).

* every ``python`` block of ``README.md`` runs;
* every script of ``examples/`` runs and says what its comments claim;
* ``examples/README.md`` lists every script;
* the figures quoted in the README (number of tests) are the registry's.
"""

from __future__ import annotations

import re
import runpy
import traceback
from pathlib import Path

import pytest

import pyxtpanel as xt

ROOT = Path(__file__).resolve().parents[2]
README = (ROOT / "README.md").read_text(encoding="utf-8")
EXAMPLES = ROOT / "examples"
SCRIPTS = sorted(EXAMPLES.glob("[0-9][0-9]_*.py"))


def test_the_readme_python_blocks_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    blocks = re.findall(r"```python\n(.*?)```", README, flags=re.DOTALL)
    assert blocks, "the README has a first example"
    namespace: dict[str, object] = {"__name__": "readme"}
    for number, code in enumerate(blocks, start=1):
        try:
            exec(compile(code, f"README block {number}", "exec"), namespace)  # noqa: S102
        except Exception:
            pytest.fail(f"README block {number} failed:\n{code}\n{traceback.format_exc()}")


def test_the_readme_quotes_the_real_number_of_tests() -> None:
    quoted = re.search(r"(\d+) tests: `xt\.list_tests\(\)`", README)
    assert quoted, "the table of the README names the number of statistical tests"
    assert int(quoted.group(1)) == len(xt.list_tests())


def test_the_examples_folder_is_what_the_index_says() -> None:
    index = (EXAMPLES / "README.md").read_text(encoding="utf-8")
    assert len(SCRIPTS) == 4
    for script in SCRIPTS:
        assert f"[{script.name}]({script.name})" in index
    for data in re.findall(r"`(data/[\w.]+)`", index):
        assert (EXAMPLES / data).exists(), data


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.stem)
def test_each_example_runs(script: Path, capsys: pytest.CaptureFixture[str]) -> None:
    runpy.run_path(str(script), run_name="__main__")
    out = capsys.readouterr().out
    assert len(out.splitlines()) >= 5, "the example prints what it computes"


def test_the_unbalanced_example_finds_the_causes_it_planted(
    capsys: pytest.CaptureFixture[str],
) -> None:
    runpy.run_path(str(EXAMPLES / "03_unbalanced_panel.py"), run_name="__main__")
    out = capsys.readouterr().out
    for cause, count in (
        ("complete", 12),
        ("internal_gap", 10),
        ("late_entry", 10),
        ("early_exit", 8),
    ):
        assert re.search(rf"{cause}\s+{count}\b", out), f"{cause}: {count} units"
    assert "holes: 0" in out, "tsfill removes the holes"


def test_the_hansen_example_reads_the_whole_panel(capsys: pytest.CaptureFixture[str]) -> None:
    runpy.run_path(str(EXAMPLES / "02_hansen_investment.py"), run_name="__main__")
    out = capsys.readouterr().out
    assert "565 units, 8475 observations" in out and "strongly balanced" in out
    assert "missing after one lag: 565 of 8475" in out
