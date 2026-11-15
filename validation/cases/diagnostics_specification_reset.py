"""Ramsey RESET against ``estat ovtest`` in Stata 17 (domain ``diagnostics_specification``).

Stata fits ``regress`` then ``estat ovtest`` (powers of the fitted values) or
``estat ovtest, rhs`` (powers of the regressors) and returns ``r(F)``, ``r(df)``, ``r(df_r)`` and
``r(p)``; the dependent variable and the regressors are exported and given to the function.

``estat ovtest`` builds the powers as ``float`` variables, so the F statistic carries a
single-precision error and these cases are compared at a relative 1e-5 (the observed differences are below 3e-6).  The ``*_exact`` cases
compute the same test in Stata in double precision (rescale, powers, ``regress``, ``test``) and
are compared at 1e-9: they show that the difference is Stata's storage, not the formula.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.specification.reset import reset
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_specification"

#: Cases compared with the single-precision ``estat ovtest``.
OFFICIAL = {
    "reset_fitted": ("y", "x z1 z2", False),
    "reset_misspecified": ("yq", "x z1 z2", False),
    "reset_rhs": ("yq", "x z1 z2", True),
    "reset_rhs_dummy": ("yq", "x z1 d", True),
    "reset_collinear": ("yq", "x z1 xc", False),
}

#: Cases compared with the double-precision computation of the same test in Stata.
EXACT = {
    "reset_fitted_exact": ("yq", "x z1 z2", False),
    "reset_rhs_exact": ("yq", "x z1 z2", True),
}


def _prelude(dep: str, regressors: str) -> list[str]:
    return [
        "quietly generate double yq = y + 1.5 * x * x",
        "quietly generate double xc = 2 * x + z1",
        f"quietly regress {dep} {regressors}",
    ]


def _exports(tag: str, dep: str, regressors: str) -> list[str]:
    names = ", ".join(f'"{v}"' for v in regressors.split())
    return [
        export(tag, "y", f'st_data(., "{dep}")'),
        export(tag, "x", f"st_data(., ({names}))"),
    ]


def _official(tag: str, dep: str, regressors: str, rhs: bool) -> list[str]:
    return [
        *_prelude(dep, regressors),
        *_exports(tag, dep, regressors),
        f"quietly estat ovtest{', rhs' if rhs else ''}",
        "scalar rf = r(F)",
        "scalar rdf = r(df)",
        "scalar rdr = r(df_r)",
        "scalar rp = r(p)",
        export_scalars(tag, "r", ["rf", "rdf", "rdr", "rp"]),
    ]


def _exact(tag: str, dep: str, regressors: str, rhs: bool) -> list[str]:
    if rhs:
        build = [
            f"foreach v in {regressors} {{",
            "    quietly summarize `v'",
            "    quietly generate double s_`v' = (`v' - r(min)) / (r(max) - r(min))",
            "    quietly generate double s2_`v' = s_`v'^2",
            "    quietly generate double s3_`v' = s_`v'^3",
            "    quietly generate double s4_`v' = s_`v'^4",
            "}",
        ]
        powers = " ".join(f"s{k}_{v}" for v in regressors.split() for k in (2, 3, 4))
    else:
        build = [
            "quietly predict double yh, xb",
            "quietly summarize yh",
            "quietly replace yh = (yh - r(min)) / (r(max) - r(min))",
            "quietly generate double yh2 = yh^2",
            "quietly generate double yh3 = yh^3",
            "quietly generate double yh4 = yh^4",
        ]
        powers = "yh2 yh3 yh4"
    return [
        *_prelude(dep, regressors),
        *_exports(tag, dep, regressors),
        *build,
        f"quietly regress {dep} {regressors} {powers}",
        f"quietly test {powers}",
        "scalar rf = r(F)",
        "scalar rdf = r(df)",
        "scalar rdr = r(df_r)",
        "scalar rp = r(p)",
        export_scalars(tag, "r", ["rf", "rdf", "rdr", "rp"]),
    ]


def _register(tag: str, spec: tuple[str, str, bool], *, exact: bool) -> None:
    dep, regressors, rhs = spec

    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        result = reset(s["y"].ravel(), s["x"], rhs=rhs)
        return {
            "r": np.array([result.statistic, *result.df, result.pvalue]),
            "y": s["y"],
            "x": s["x"],
        }

    build = _exact if exact else _official
    case(
        tag,
        domain=DOMAIN,
        dataset="sim_het",
        stata=build(tag, dep, regressors, rhs),
        atol=1e-9 if exact else 1e-7,
        rtol=1e-9 if exact else 1e-5,
    )(run)


for _tag, _spec in OFFICIAL.items():
    _register(_tag, _spec, exact=False)
for _tag, _spec in EXACT.items():
    _register(_tag, _spec, exact=True)
