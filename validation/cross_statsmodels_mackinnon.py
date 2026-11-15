"""Freeze statsmodels' MacKinnon (1994) coefficients, to check the embedded table (L.1.3).

    python -m validation.cross_statsmodels_mackinnon

``statsmodels`` implements the same paper independently; its tabulated coefficients and a few
p-values from ``mackinnonp`` are written to ``tests/golden/mackinnon_statsmodels.json``.
``tests/unit/diagnostics/test_tables.py`` compares the embedded ``mackinnon_adf.csv`` with them.
The package never needs ``statsmodels`` at test time.  Needs ``statsmodels`` in the interpreter
that runs it.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from statsmodels.tsa import adfvalues
from statsmodels.tsa.stattools import mackinnonp

TARGET = Path(__file__).resolve().parents[1] / "tests" / "golden" / "mackinnon_statsmodels.json"
TAUS = [-25.0, -17.0, -10.0, -4.0, -3.5, -3.0, -2.89, -2.5, -1.9, -1.61, -1.0, 0.0, 0.5, 2.0, 3.0]


def main() -> None:
    out: dict[str, object] = {"version": __import__("statsmodels").__version__, "cases": {}}
    cases = out["cases"]
    assert isinstance(cases, dict)
    for case in ("c", "ct"):
        star, low, high = (getattr(adfvalues, f"tau_{k}_{case}") for k in ("star", "min", "max"))
        small = getattr(adfvalues, f"tau_{case}_smallp")
        large = getattr(adfvalues, f"tau_{case}_largep")
        cases[case] = {  # the first row/entry of each table is the one for a single series (N = 1)
            "tau_star": float(star[0]),
            "tau_min": float(low[0]),
            "tau_max": float(high[0]),
            "small": [float(v) for v in small[0]] + [0.0],
            "large": [float(v) for v in large[0]],
            "p": {str(t): float(mackinnonp(t, regression=case, N=1)) for t in TAUS},
        }
    TARGET.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print("wrote", TARGET, np.__version__)


if __name__ == "__main__":
    main()
