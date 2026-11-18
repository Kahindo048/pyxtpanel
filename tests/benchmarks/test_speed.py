"""Speed guard (I.4.1): each operation against a reference measured in the same process.

Run with ``pytest -m bench`` (the benchmarks are left out of the default run: a timing is only
meaningful on a quiet machine).  The time of an operation is divided by the time of a fixed numpy
workload measured just before (``tools.perf.calibration``), so a slower machine moves both; the
ratio is compared with ``baseline.json`` and may grow by ``MARGIN`` before the test fails.
Refresh the reference after a deliberate change with ``python -m tools.perf --baseline
tests/benchmarks/baseline.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from tools.perf import calibration, guard_rows, operations, simulate

import pyxtpanel as xt

pytest.importorskip("pytest_benchmark")  # a development extra: absent, the guard is skipped
pytestmark = pytest.mark.bench

BASELINE = json.loads((Path(__file__).parent / "baseline.json").read_text(encoding="utf-8"))
#: A ratio may grow to this multiple of the recorded one before the guard fails.
MARGIN = 3.0
#: Operations shorter than this are dominated by noise: they must stay under it, whatever the ratio.
FLOOR_SECONDS = 0.003


@pytest.fixture(scope="module")
def frame():  # type: ignore[no-untyped-def]
    assert BASELINE["rows"] == guard_rows()
    return simulate(guard_rows())


@pytest.fixture(scope="module")
def reference() -> float:
    return calibration()


def check(seconds: float, name: str, reference: float) -> None:
    ratio = seconds / reference
    allowed = BASELINE["ratios"][name] * MARGIN
    assert seconds <= FLOOR_SECONDS or ratio <= allowed, (
        f"{name}: {seconds * 1000:.1f} ms = {ratio:.2f} x the reference, "
        f"recorded {BASELINE['ratios'][name]:.2f} x (allowed {allowed:.2f} x)"
    )


def test_every_operation_has_a_recorded_ratio() -> None:
    assert set(BASELINE["ratios"]) == {"xtset (shuffled rows)", *operations()}


def test_declaration(benchmark, frame, reference) -> None:  # type: ignore[no-untyped-def]
    benchmark.pedantic(lambda: xt.xtset(frame, "id", "t"), rounds=5, warmup_rounds=1)
    if benchmark.stats:
        check(benchmark.stats.stats.min, "xtset (shuffled rows)", reference)


@pytest.mark.parametrize("name", list(operations()))
def test_operation(benchmark, frame, reference, name) -> None:  # type: ignore[no-untyped-def]
    panel = xt.xtset(frame, "id", "t")
    function = operations()[name]
    function(panel)  # parts that belong to the panel (cached on first use) are not measured
    benchmark.pedantic(lambda: function(panel), rounds=5, warmup_rounds=1)
    if benchmark.stats:
        check(benchmark.stats.stats.min, name, reference)
