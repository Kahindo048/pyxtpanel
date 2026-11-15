"""Freeze the moments and weight matrices of ``pyxtabond2`` for the over-identification tests.

    python -m validation.cross_pyxtabond2_overid [path to the pyxtabond2 repository]

``pyxtabond2`` is not a dependency of this package.  Its engine builds the instrument matrix,
the residuals and the weight matrices of a difference GMM; the numbers the tests take as input
(``m = Z'e``, ``W`` and the scale of Sargan) are written to
``tests/golden/overid_pyxtabond2.json``.  The cases of
``validation/cases/diagnostics_postestimation_xtabond2.py`` give them to ``hansen`` and
``sargan`` and compare the result with ``e(hansen)`` and ``e(sargan)`` of Stata's ``xtabond2``,
so the package never needs ``pyxtabond2`` at test time.  Needs ``statsmodels`` in the
interpreter that runs it.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "tests" / "golden" / "overid_pyxtabond2.json"
DEFAULT_REPOSITORY = Path(r"C:\Users\Fred-derick\Desktop\pyxtabond2")

#: (tag, panel, two-step, forward orthogonal deviations); one lag of x, instruments from x.
#: The orthogonal-deviations fit is not a case: with ``lag(1 .)`` pyxtabond2 builds 21 instruments
#: on the balanced panel where ``xtabond2 ..., orthogonal`` builds 27, so the two estimators do not
#: share inputs and the comparison would not test the over-identification tests.
SPECS = (
    ("ab1_balanced", "balanced", False, False),
    ("ab2_balanced", "balanced", True, False),
    ("ab2_gaps", "gaps", True, False),
    ("ab2_entry_exit", "entry_exit", True, False),
)


def main(argv: list[str]) -> int:
    repository = Path(argv[0]) if argv else DEFAULT_REPOSITORY
    sys.path.insert(0, str(repository))
    sys.path.insert(0, str(ROOT))
    warnings.simplefilter("ignore")
    from pyxtabond2 import PanelData, PyXtabond2

    from validation import panels

    frozen: dict[str, dict[str, object]] = {}
    for tag, name, twostep, orthogonal in SPECS:
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
            orthogonal=orthogonal,
            lag_limits_diff=(1, None) if orthogonal else (2, None),
            twostep=twostep,
        )
        result = model.fit()
        eng = result.engine
        sargan, hansen = result.diag["sargan"], result.diag["hansen"]
        m1 = (eng.Z.T @ eng.e1).ravel()
        quad1 = float(m1 @ eng.W1 @ m1)
        entry: dict[str, object] = {
            "panel": name,
            "n_parameters": int(eng.k_vars),
            "n_instruments": int(eng.n_instruments),
            "n_groups": int(eng.N_groups),
            "sargan_moments": m1.tolist(),
            "sargan_weight": eng.W1.tolist(),
            "sargan_scale": quad1 / float(sargan[0]),
        }
        if twostep:
            m2 = (eng.Z.T @ eng.e2).ravel()
            entry["hansen_moments"] = m2.tolist()
            entry["hansen_weight"] = eng.W2.tolist()
            entry["pyxtabond2_hansen"] = float(hansen[0])
        entry["pyxtabond2_sargan"] = float(sargan[0])
        frozen[tag] = entry
    TARGET.write_text(json.dumps(frozen, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
