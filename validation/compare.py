"""Comparison of Python results with the objects exported from Stata.

Tolerance rule: an element agrees when its absolute difference is at most ``atol`` **or** its
relative difference is at most ``rtol``.  The defaults (5e-5 and 1e-6) are the floor of the
project; closed-form quantities are compared with much tighter values (1e-10 relative is the
aim; measured in B.3.4 by ``validation/self_agreement.py``: Stata differs from itself by at most
5e-14 and the package from Stata by at most 4e-13 on the campaign).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

ATOL = 5e-5
RTOL = 1e-6


@dataclass(frozen=True)
class Diff:
    """Comparison of one named object."""

    name: str
    shape_python: tuple[int, ...]
    shape_stata: tuple[int, ...]
    max_abs: float
    max_rel: float
    ok: bool

    def line(self) -> str:
        status = "OK  " if self.ok else "FAIL"
        if self.shape_python != self.shape_stata:
            return f"{status} {self.name:<16} shape {self.shape_python} vs {self.shape_stata}"
        return (
            f"{status} {self.name:<16} max|abs| = {self.max_abs:.3e}  max|rel| = {self.max_rel:.3e}"
        )


def read_matrix(path: Path) -> np.ndarray:
    """Read a CSV written by the Stata ``pyx_export`` helper (missing values become NaN)."""
    rows = [
        [np.nan if f.strip() in ("", ".") else float(f) for f in raw.split(",")]
        for raw in path.read_text(encoding="ascii").splitlines()
        if raw.strip()
    ]
    return np.array(rows, dtype=float)


def load_stata(out_dir: Path, tag: str) -> dict[str, np.ndarray]:
    """All exported objects of one case: ``{name: array}`` from ``<tag>__<name>.csv``."""
    prefix = f"{tag}__"
    found = {p.stem[len(prefix) :]: read_matrix(p) for p in sorted(out_dir.glob(f"{prefix}*.csv"))}
    if not found:
        raise FileNotFoundError(f"No exported file for case {tag!r} in {out_dir}")
    return found


def compare_arrays(
    name: str, python: np.ndarray, stata: np.ndarray, atol: float = ATOL, rtol: float = RTOL
) -> Diff:
    """Compare two arrays element by element (NaN must be NaN on both sides)."""
    a = np.atleast_2d(np.asarray(python, dtype=float))
    b = np.atleast_2d(np.asarray(stata, dtype=float))
    if a.shape != b.shape and a.T.shape == b.shape and 1 in a.shape:
        a = a.T  # a vector and the same vector as a column are the same object here
    if a.shape != b.shape or (np.isnan(a) ^ np.isnan(b)).any():
        return Diff(name, a.shape, b.shape, np.inf, np.inf, False)
    both_nan = np.isnan(a) & np.isnan(b)
    diff = np.abs(a - b)[~both_nan]
    rel = diff / np.maximum(np.abs(b)[~both_nan], 1e-300)
    max_abs = float(diff.max()) if diff.size else 0.0
    max_rel = float(rel.max()) if rel.size else 0.0
    return Diff(
        name, a.shape, b.shape, max_abs, max_rel, bool(np.all((diff <= atol) | (rel <= rtol)))
    )


def compare_case(
    python: dict[str, np.ndarray],
    stata: dict[str, np.ndarray],
    atol: float = ATOL,
    rtol: float = RTOL,
) -> list[Diff]:
    """Compare every object present on the Stata side; a missing Python object fails."""
    diffs = []
    for name, expected in stata.items():
        if name not in python:
            diffs.append(Diff(name, (), expected.shape, np.inf, np.inf, False))
        else:
            diffs.append(compare_arrays(name, python[name], expected, atol, rtol))
    return diffs


def report(tag: str, diffs: list[Diff]) -> str:
    """Readable report of one case."""
    verdict = "PASS" if all(d.ok for d in diffs) else "FAIL"
    return "\n".join([f"[{verdict}] case {tag}", *(f"  {d.line()}" for d in diffs)])
