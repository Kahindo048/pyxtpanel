"""Freeze the Hansen and AR inputs of ``pyxtdpthresh`` next to Stata's own statistics (L.7.4).

    python validation/cross_pyxtdpthresh.py [path to the pyxtdpthresh repository]

``pyxtdpthresh`` reproduces the Hansen J and the Arellano-Bond AR(1) and AR(2) statistics of
``xtdpthresh`` in Stata 17 (``tests/golden/variance_golden.json`` of that repository, to 1e-10).
This script runs its estimator on the golden cases chosen in :data:`CASES` and writes to
``tests/golden/overid_ar_pyxtdpthresh.json``:

* the inputs of ``hansen`` (the moments ``Z'e``, the weight matrix, the numbers of instruments
  and parameters) and of ``ar_test`` (the residuals, regressors, instruments, weight and
  covariance of the estimation and test systems), exactly as the estimator hands them to its own
  tests;
* Stata's statistics for the same case, copied from the golden file (``e(hansen)``,
  ``e(ar1)``, ``e(ar1_p)``, ``e(ar1_cond)``, ``e(ar1_b0)``, ``e(ar1_T1)``, ``e(ar1_TT)`` and the
  same for AR(2)).

``tests/golden/test_against_pyxtdpthresh.py`` feeds ``hansen`` and ``ar_test`` of this package
with the inputs and compares the result with Stata's numbers; the package never imports
``pyxtdpthresh`` at test time.  The script needs an interpreter that has ``numba`` (the
estimator's kernels) and must not run from a directory that puts this package's ``validation``
before the one of ``pyxtdpthresh``; it does not import ``pyxtpanel``.
"""

from __future__ import annotations

import base64
import hashlib
import json
import sys
import warnings
import zlib
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "tests" / "golden" / "overid_ar_pyxtdpthresh.json"
DEFAULT_REPOSITORY = Path(r"C:\Users\Fred-derick\Desktop\Dynamic threshold\pyxtdpthresh")

#: Golden cases of ``pyxtdpthresh`` (variance block): first differences and forward orthogonal
#: deviations (tested in first differences), jump and kink, conditional and Windmeijer variances,
#: collapsed and static instruments, a panel with gaps, one where Hansen is not defined, and the
#: 565-firm panel of the help examples.
CASES = (
    "fd_jump",
    "fd_jump_wind",
    "fd_kink_wind",
    "fd_nocenter",
    "fd_static",
    "fd_collapse_wind",
    "fd_gaps_wind",
    "fd_alternate_wind",
    "fod_jump_wind",
    "fod_gaps_wind",
    "full_unbalanced_gaps_wind",
)

STATA_SCALARS = (
    "hansen", "hansen_df", "hansen_p", "n", "n_units", "ar_joint",
    "ar1", "ar1_p", "ar1_cond", "ar1_b0", "ar1_T1", "ar1_TT",
    "ar2", "ar2_p", "ar2_cond", "ar2_b0", "ar2_T1", "ar2_TT",
)  # fmt: skip


def pack(array: np.ndarray) -> dict[str, Any]:
    """Array as zlib-compressed little-endian doubles (the layout of the other golden files)."""
    flat = np.ascontiguousarray(array, dtype="<f8")
    return {
        "shape": list(flat.shape),
        "data": base64.b64encode(zlib.compress(flat.tobytes(), 9)).decode(),
    }


def share(store: dict[str, Any], array: np.ndarray) -> str:
    """Store an array once under the digest of its bytes (cases share their instruments)."""
    flat = np.ascontiguousarray(array, dtype="<f8")
    key = hashlib.sha256(flat.tobytes() + str(flat.shape).encode()).hexdigest()[:16]
    store.setdefault(key, pack(flat))
    return key


def offsets_of(unit: np.ndarray) -> np.ndarray:
    """Start of each unit in rows sorted by unit, with the number of rows at the end."""
    unit = np.asarray(unit)
    return np.concatenate([[0], np.flatnonzero(np.diff(unit)) + 1, [unit.size]]).astype(np.int64)


def within_unit_position(offsets: np.ndarray) -> np.ndarray:
    """Position of each row within its unit (the estimation rows carry no calendar)."""
    return np.concatenate(
        [np.arange(b - a) for a, b in zip(offsets[:-1], offsets[1:], strict=True)]
    )


def equations(resid: Any, regressors: Any, unit: Any, time: Any | None) -> dict[str, Any]:
    offsets = offsets_of(unit)
    position = within_unit_position(offsets) if time is None else np.asarray(time)
    return {
        "residuals": pack(np.asarray(resid, dtype=float)),
        "regressors": pack(np.asarray(regressors, dtype=float)),
        "offsets": offsets.tolist(),
        "time_index": np.asarray(position, dtype=np.int64).tolist(),
    }


def main(argv: list[str]) -> int:
    repository = Path(argv[0]) if argv else DEFAULT_REPOSITORY
    sys.path.insert(0, str(repository / "src"))
    sys.path.insert(0, str(repository))
    warnings.simplefilter("ignore")
    from pyxtdpthresh.models.xtdpthresh import diagnostics
    from pyxtdpthresh.models.xtdpthresh.fit import fit_point

    from validation.variance_cases import CASES as ALL_CASES
    from validation.variance_cases import GOLDEN, frame

    meta = json.loads(GOLDEN.read_text(encoding="utf-8"))
    stata = meta["cases"]
    known = {c.name: c for c in ALL_CASES}
    frozen: dict[str, Any] = {
        "meta": {
            "stata": meta["meta"]["stata"],
            "xtdpthresh": meta["meta"]["xtdpthresh"],
            "source": "tests/golden/variance_golden.json of pyxtdpthresh",
            "scale": "weight = A / n_rows, so that the Hansen J is m' W m",
        },
        "arrays": {},
        "cases": {},
    }
    arrays: dict[str, Any] = frozen["arrays"]
    original = diagnostics.ar_test
    for name in CASES:
        case = known[name]
        args = case.full_args()
        y, q = args.pop("y"), args.pop("q")
        x = args.pop("x", ())
        calls: list[tuple[int, Any, Any, Any]] = []

        def recorder(order, test, est, cov, _original=original, _calls=calls):  # noqa: ANN001, ANN202
            _calls.append((order, test, est, cov))
            return _original(order, test, est, cov)

        diagnostics.ar_test = recorder
        try:
            fit = fit_point(frame(case), y, x, q=q, id="firm", time="year", **args)
        finally:
            diagnostics.ar_test = original
        stack, search = fit.stack, fit.search.fit
        residual = stack.dy - stack.regressors(search.gamma) @ search.theta
        moments = stack.z.dense.T @ residual
        scalars = {k: float(stata[name]["scalars"][k]) for k in STATA_SCALARS}
        # the calls are the conditional pair, then (when it exists) the joint pair
        pairs = {"cond": calls[:2], "joint": calls[2:4]}
        reported = pairs["joint" if fit.serial.joint else "cond"]
        entry: dict[str, Any] = {
            "stata": scalars,
            "hansen": {
                "moments": moments.tolist(),
                "weight": share(arrays, search.weight / stack.n_rows),
                "n_parameters": int(search.theta.size) + 1,  # the threshold is a parameter too
                "n_instruments": int(stack.z.dense.shape[1]),
                "n_groups": int(stack.n_units),
                # a one-step weight is not efficient: the estimator and Stata report no J for it
                "twostep": bool(search.twostep),
            },
            "ar": {},
        }
        for label, chosen in (("reported", reported), ("conditional", pairs["cond"])):
            order, test, est, cov = chosen[0]
            systems = {
                "test": equations(test.resid, test.regressors, test.unit, test.time),
                "estimation": equations(est.resid, est.regressors, est.unit, None),
                "instruments": share(arrays, est.z.dense),
                "weight": share(arrays, est.weight),
                "covariance": share(arrays, cov),
                "orders": [int(c[0]) for c in chosen],
            }
            # the second order shares the systems; only the order changes
            assert all(c[1] is test and c[2] is est for c in chosen)
            entry["ar"][label] = systems
        frozen["cases"][name] = entry
        print(
            f"{name}: hansen {fit.hansen.statistic:.6g} (Stata {scalars['hansen']:.6g}), "
            f"joint AR {fit.serial.joint}"
        )
    TARGET.write_text(json.dumps(frozen) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {TARGET} ({TARGET.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
