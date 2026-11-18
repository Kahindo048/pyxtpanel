"""Hansen and Arellano-Bond tests against the frozen Stata 17 runs of ``xtdpthresh`` (L.7.4).

``tests/golden/overid_ar_pyxtdpthresh.json`` holds, for 11 golden cases of ``pyxtdpthresh``
(``python validation/cross_pyxtdpthresh.py``), the inputs its estimator hands to its Hansen and
AR tests (moments, weight, residual systems, instruments, covariance) and the statistics that
``xtdpthresh`` printed in Stata 17 for the same case.  The functions of this package must give
Stata's numbers from those inputs: the J statistic with its degrees of freedom (the threshold
counts as a parameter) and p-value, and AR(1) and AR(2) with their pieces ``b0``, ``T1`` and
``T2 + T3``, for the covariance that treats the threshold as known and for the joint one.  The
weight matrix is stored divided by the number of equations, the scale at which ``J = m' W m``.
"""

from __future__ import annotations

import base64
import json
import math
import zlib
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from pyxtpanel.diagnostics import Equations, ar_test, hansen

pytestmark = pytest.mark.golden

FIXTURE = Path(__file__).resolve().parent / "overid_ar_pyxtdpthresh.json"
MASTER: dict[str, Any] = json.loads(FIXTURE.read_text(encoding="utf-8"))
CASES: dict[str, Any] = MASTER["cases"]
ARRAYS: dict[str, Any] = MASTER["arrays"]  # matrices shared by the cases, under their digest
#: the estimator and Stata agree to rounding; the stored inputs are exact doubles
RTOL = 1e-8


def unpack(entry: dict[str, Any]) -> np.ndarray:
    raw = zlib.decompress(base64.b64decode(entry["data"]))
    out: np.ndarray = np.frombuffer(raw, dtype="<f8").reshape(entry["shape"]).copy()
    return out


def system(entry: dict[str, Any], instruments: np.ndarray | None = None) -> Equations:
    return Equations(
        unpack(entry["residuals"]),
        unpack(entry["regressors"]),
        np.array(entry["offsets"]),
        np.array(entry["time_index"]),
        instruments,
    )


def close(got: float, want: float, what: str) -> None:
    if math.isnan(want):
        assert math.isnan(got), f"{what}: got {got}, Stata has a missing value"
    else:
        assert got == pytest.approx(want, rel=RTOL, abs=1e-12), what


@pytest.mark.parametrize("name", list(CASES))
def test_hansen_is_the_statistic_stata_prints(name: str) -> None:
    case = CASES[name]
    stata, inputs = case["stata"], case["hansen"]
    if not inputs["twostep"]:  # the weight of a one-step fit is not efficient: no J in Stata
        assert math.isnan(stata["hansen"])
        return
    result = hansen(
        np.array(inputs["moments"]),
        unpack(ARRAYS[inputs["weight"]]),
        inputs["n_parameters"],
        n_instruments=inputs["n_instruments"],
        n_groups=inputs["n_groups"],
    )
    close(result.statistic, stata["hansen"], "e(hansen)")
    close(result.pvalue, stata["hansen_p"], "e(hansen_p)")
    assert tuple(float(d) for d in result.df) == (stata["hansen_df"],)


def _statistics(name: str, label: str) -> dict[int, Any]:
    systems = CASES[name]["ar"][label]
    instruments = unpack(ARRAYS[systems["instruments"]])
    test = system(systems["test"])
    estimation = system(systems["estimation"], instruments)
    weight, covariance = unpack(ARRAYS[systems["weight"]]), unpack(ARRAYS[systems["covariance"]])
    return {
        order: ar_test(
            order,
            test,
            weight=weight,
            covariance=covariance,
            estimation=estimation,
            min_pairs=5,  # the threshold of xtdpthresh
        )
        for order in systems["orders"]
    }


@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("name", list(CASES))
def test_ar_reported_by_stata_is_reproduced(name: str, order: int) -> None:
    stata = CASES[name]["stata"]
    assert stata["ar_joint"] == 1.0  # every chosen case reports the joint variance
    result = _statistics(name, "reported")[order]
    close(result.statistic, stata[f"ar{order}"], f"e(ar{order})")
    close(result.pvalue, stata[f"ar{order}_p"], f"e(ar{order}_p)")
    pieces = result.details["value"]
    close(float(pieces["d0"]), stata[f"ar{order}_b0"], f"e(ar{order}_b0)")
    close(float(pieces["T1"]), stata[f"ar{order}_T1"], f"e(ar{order}_T1)")
    close(float(pieces["T2"] + pieces["T3"]), stata[f"ar{order}_TT"], f"e(ar{order}_TT)")


@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("name", list(CASES))
def test_ar_with_the_threshold_taken_as_known_is_reproduced(name: str, order: int) -> None:
    stata = CASES[name]["stata"]
    result = _statistics(name, "conditional")[order]
    close(result.statistic, stata[f"ar{order}_cond"], f"e(ar{order}_cond)")


def test_the_frozen_cases_cover_both_transformations_and_an_unbalanced_panel() -> None:
    names = set(CASES)
    assert {"fd_jump", "fod_jump_wind", "fd_gaps_wind", "fod_gaps_wind"} <= names
    assert any(not c["hansen"]["twostep"] for c in CASES.values())
