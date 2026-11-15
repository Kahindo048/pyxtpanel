"""Stata reads the .dta files that pyxtpanel writes and finds the panel declared (needs Stata)."""

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from pyxtpanel.panel.base import PanelBase

ROOT = Path(__file__).resolve().parents[2]
STATA = Path(r"C:\Program Files\Stata17\StataMP-64.exe")
RUNNER = ROOT / "validation" / "run_stata.ps1"

pytestmark = [
    pytest.mark.stata,
    pytest.mark.skipif(
        not (STATA.exists() and shutil.which("powershell")), reason="Stata 17 is not available"
    ),
]


def test_stata_finds_the_xtset_written_by_pyxtpanel(tmp_path: Path) -> None:
    quarters = pd.DataFrame(
        {
            "id": [1, 1, 1, 2, 2],
            "q": pd.to_datetime(
                ["2020-01-01", "2020-04-01", "2020-10-01", "2020-01-01", "2020-04-01"]
            ),
            "x": [1.0, 2.0, 3.0, 4.0, 5.0],
        }
    )
    PanelBase.from_frame(quarters, "id", "q").to_stata(tmp_path / "quarters.dta")
    stepped = pd.DataFrame({"id": [1, 1, 1], "t": [1, 3, 7], "x": [1.0, 2.0, 3.0]})
    PanelBase.from_frame(stepped, "id", "t", delta=2).to_stata(tmp_path / "stepped.dta")
    do_file = tmp_path / "check.do"
    do_file.write_text(
        'use "quarters.dta", clear\n'
        "xtset\n"
        'display "Q:" r(tdeltas) "|" r(tmins) "|" r(tmaxs) "|" r(gaps) "|" r(balanced)\n'
        'use "stepped.dta", clear\n'
        "xtset\n"
        'display "S:" r(tdelta) "|" r(tmin) "|" r(tmax) "|" r(gaps) "|" r(balanced)\n',
        encoding="ascii",
    )
    subprocess.run(
        ["powershell", "-NoProfile", "-File", str(RUNNER), "-Do", str(do_file), "-Quiet"],
        check=True,
        timeout=600,
    )
    log = (tmp_path / "_run_check.log").read_text(encoding="utf-8", errors="replace")
    assert "Q:1 quarter|2020q1|2020q4|1|unbalanced" in log
    assert "S:2|1|7|1|strongly balanced" in log
