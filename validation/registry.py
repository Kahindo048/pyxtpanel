"""Registry of validation cases: one file per domain in ``validation/cases/``.

A case says which dataset Stata loads, which commands it runs and what it exports, and how to
compute the same objects in Python.  Adding a case is writing one decorated function in a file
of ``validation/cases/``; the do-file, the Stata run, the freezing and the golden tests discover
it by themselves.
"""

from __future__ import annotations

import importlib
import pkgutil
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

#: Mata helper that writes a matrix to a CSV with 15 significant digits (B.1.2).
MATA_EXPORT = """\
mata:
void pyx_export(string scalar fname, real matrix M)
{
    real scalar fh, i, j
    string scalar line
    if (fileexists(fname)) unlink(fname)
    fh = fopen(fname, "w")
    for (i = 1; i <= rows(M); i++) {
        line = strofreal(M[i,1], "%21.15g")
        for (j = 2; j <= cols(M); j++) line = line + "," + strofreal(M[i,j], "%21.15g")
        fput(fh, line)
    }
    fclose(fh)
}
end
"""


def export(tag: str, name: str, matrix_expr: str) -> str:
    """Stata line that exports the Mata expression ``matrix_expr`` as ``<tag>__<name>.csv``."""
    return f'mata: pyx_export("../out/{tag}__{name}.csv", {matrix_expr})'


def export_scalars(tag: str, name: str, scalars: list[str]) -> str:
    """Stata line that exports the named Stata scalars as a column vector."""
    stack = " \\ ".join(f'st_numscalar("{s}")' for s in scalars)
    return export(tag, name, f"({stack})")


@dataclass(frozen=True)
class Case:
    """One validation case.

    Attributes
    ----------
    tag
        Short unique name, also the prefix of the exported files.
    domain
        Name of the file in ``validation/cases/`` (and of the golden file).
    dataset
        ``.dta`` file (without extension) the do-file loads from ``validation/data``.
    stata
        Do-file lines run after loading the data; they must write the exports.
    python
        Function that returns the same objects, as ``{name: array}``.
    atol, rtol
        Tolerances of the comparison.
    """

    tag: str
    domain: str
    dataset: str
    stata: tuple[str, ...]
    python: Callable[[], dict[str, np.ndarray]]
    atol: float = 1e-12
    rtol: float = 1e-12


CASES: dict[str, Case] = {}

#: Mata code a domain needs before its cases (``register_prelude``).
PRELUDES: dict[str, str] = {}


def register_prelude(domain: str, code: str) -> None:
    """Add Mata (or Stata) code that the do-file of ``domain`` runs before its cases."""
    PRELUDES[domain] = code


def case(
    tag: str,
    *,
    domain: str,
    dataset: str,
    stata: list[str],
    atol: float = 1e-12,
    rtol: float = 1e-12,
) -> Callable[[Callable[[], dict[str, np.ndarray]]], Callable[[], dict[str, np.ndarray]]]:
    """Register the decorated function as the Python side of a case."""

    def wrap(func: Callable[[], dict[str, np.ndarray]]) -> Callable[[], dict[str, np.ndarray]]:
        if tag in CASES:
            raise ValueError(f"case {tag!r} is already registered")
        CASES[tag] = Case(tag, domain, dataset, tuple(stata), func, atol, rtol)
        return func

    return wrap


def discover() -> dict[str, Case]:
    """Import every module of ``validation.cases`` so that their decorators run."""
    import validation.cases as package

    for info in pkgutil.iter_modules(package.__path__, package.__name__ + "."):
        importlib.import_module(info.name)
    return CASES


def of_domain(domain: str) -> list[Case]:
    """Cases of one domain, in registration order."""
    return [c for c in discover().values() if c.domain == domain]
