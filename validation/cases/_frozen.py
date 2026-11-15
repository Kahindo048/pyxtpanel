"""Read back the objects Stata exported and froze, to use them as inputs of a test (L.6).

The Hausman and Wald tests take coefficients and covariances.  The cases export ``e(b)`` and
``e(V)`` from Stata's own estimators next to Stata's test result, and the Python side feeds the
function with those very numbers, so that no estimator of this package or of NumPy stands between
the function and Stata.  The inputs are read from the frozen file of the domain; they are also
returned by the case, so that they are compared (trivially) like any exported object.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

GOLDEN = Path(__file__).resolve().parents[2] / "tests" / "golden"


def stata_objects(domain: str, tag: str) -> dict[str, np.ndarray]:
    """Return the frozen objects of the case ``tag`` as arrays (matrices keep their shape)."""
    master = json.loads((GOLDEN / f"{domain}_golden.json").read_text(encoding="utf-8"))
    return {name: np.asarray(value, dtype=float) for name, value in master["cases"][tag].items()}
