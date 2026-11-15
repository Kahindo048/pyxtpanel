"""Freeze the equations, weight and covariance of ``pyxtabond2`` for the AR(m) test (L.6.1).

    python -m validation.cross_pyxtabond2_ar [path to the pyxtabond2 repository]

The engine of ``pyxtabond2`` builds the difference-GMM system (differenced residuals and
regressors, instruments, weight matrix, robust covariance).  They are written to
``tests/golden/ar_pyxtabond2.json`` and given to ``ar_test`` by the cases of
``validation/cases/diagnostics_postestimation_ar_test.py``, whose reference is Stata's
``xtabond2`` and ``estat abond``; the package never needs ``pyxtabond2`` at test time.  Needs
``statsmodels`` in the interpreter that runs it.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "tests" / "golden" / "ar_pyxtabond2.json"
DEFAULT_REPOSITORY = Path(r"C:\Users\Fred-derick\Desktop\pyxtabond2")

#: (tag, panel, two-step); robust variances, one lag of x, instruments ``gmm(x, lag(2 .))``.
SPECS = (
    ("ab1_balanced", "balanced", False),
    ("ab2_balanced", "balanced", True),
    ("ab2_gaps", "gaps", True),
    ("ab2_entry_exit", "entry_exit", True),
)


def main(argv: list[str]) -> int:
    repository = Path(argv[0]) if argv else DEFAULT_REPOSITORY
    sys.path.insert(0, str(repository))
    sys.path.insert(0, str(ROOT))
    warnings.simplefilter("ignore")
    from pyxtabond2 import PanelData, PyXtabond2

    from validation import panels

    frozen: dict[str, dict[str, object]] = {}
    for tag, name, twostep in SPECS:
        frame = panels.load(name)[["id", "t", "x"]].copy()
        old = PanelData(frame.copy(), "id", "t")
        ready = old.data.reset_index()
        ready["L1_x"] = old.get_lag("x", 1).to_numpy()
        model = PyXtabond2(
            ready,
            id_col="id",
            time_col="t",
            dep_var="x",
            x_vars=["L1_x"],
            gmm_vars=["x"],
            model_type="difference",
            lag_limits_diff=(2, None),
            twostep=twostep,
            robust=True,
            artests=2,
        )
        result = model.fit()
        eng = result.engine
        residuals = (eng.e2 if twostep else eng.e1).ravel()
        weight = eng.W2 if twostep else eng.W1
        covariance = eng.V2_robust  # also the robust covariance of a one-step fit
        groups = np.asarray(eng.group_ids)
        time = np.asarray(eng.t_index, dtype=np.int64)
        order = np.lexsort((time, groups))
        _, first = np.unique(groups[order], return_index=True)
        frozen[tag] = {
            "panel": name,
            "residuals": residuals[order].tolist(),
            "regressors": eng.X[order].tolist(),
            "instruments": eng.Z[order].tolist(),
            "offsets": np.append(first, order.size).tolist(),
            "time_index": time[order].tolist(),
            "weight": weight.tolist(),
            "covariance": np.asarray(covariance).tolist(),
            "pyxtabond2_ar": {
                str(k): [float(v[0]), float(v[1])] for k, v in result.diag["ar"].items()
            },
        }
    TARGET.write_text(json.dumps(frozen) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
