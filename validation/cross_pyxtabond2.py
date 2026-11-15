"""Freeze the observations ``pyxtabond2`` uses, to validate ``usable_equations`` (K.2.4, K.7).

    python -m validation.cross_pyxtabond2 [path to the pyxtabond2 repository]

``pyxtabond2`` is not a dependency of this package: the numbers it produces (the number of
equations of a difference or orthogonal-deviations GMM with one lag of the dependent variable)
are written to ``tests/golden/equations_pyxtabond2.json`` and compared by a plain test.  Needs
``statsmodels`` in the interpreter that runs it.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "tests" / "golden" / "equations_pyxtabond2.json"
PANELS = ("balanced", "gaps", "entry_exit", "weak", "unsorted_ids")
DEFAULT_REPOSITORY = Path(r"C:\Users\Fred-derick\Desktop\pyxtabond2")


def main(argv: list[str]) -> int:
    repository = Path(argv[0]) if argv else DEFAULT_REPOSITORY
    sys.path.insert(0, str(repository))
    sys.path.insert(0, str(ROOT))
    warnings.simplefilter("ignore")
    from pyxtabond2 import PanelData, PyXtabond2

    from validation import panels

    frozen: dict[str, dict[str, int]] = {}
    for name in PANELS:
        frame = panels.load(name)[["id", "t", "x"]].copy()
        old = PanelData(frame.copy(), "id", "t")
        ready = old.data.reset_index()
        ready["L1_x"] = old.get_lag("x", 1).to_numpy()
        for transform, orthogonal, limits in (("fd", False, (2, None)), ("fod", True, (1, None))):
            model = PyXtabond2(
                ready,
                id_col="id",
                time_col="t",
                dep_var="x",
                x_vars=["L1_x"],
                gmm_vars=["x"],
                model_type="difference",
                orthogonal=orthogonal,
                lag_limits_diff=limits,
            )
            frozen.setdefault(name, {})[transform] = int(model.fit().engine.n_obs)
    TARGET.write_text(json.dumps(frozen, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
