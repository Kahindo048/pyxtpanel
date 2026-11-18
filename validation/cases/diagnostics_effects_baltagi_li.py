"""Baltagi-Li tests against ``xttest1`` and an independent Mata computation (domain ``diagnostics_effects``).

``xttest1`` (Sosa-Escudero and Bera, update ``sg164_1`` of the Stata Journal 8(1), installed for
this validation) runs after ``xtreg, re``; ``xttest1, unadjusted`` prints and returns the seven
statistics.  Stata fits the pooled regression (its residuals are exported) and the function
receives the residuals, the unit boundaries and the periods.  Each model is registered twice:
against ``xttest1``, and against a Mata block that rebuilds the statistics from the residuals
unit by unit (``panelsetup`` loops, the definitions of the papers, nothing shared with the ado).
``xttest1`` refuses a panel with gaps (r(498)), and so does the function given the periods.
The one-sided p-values are the right tail of the standard normal.

Two defects of ``xttest1`` limit the comparison with it.  It stores the lagged residual as a
``float``, so its statistics agree with the double-precision ones to about 2e-6 (relative
tolerance 1e-5 for its cases, 1e-9 for the Mata ones).  It returns the one-sided p-values as
``1 - normprob(z)``, which loses digits from ``z = 5`` and is 0 beyond ``z = 8``: they are compared
with the Mata ones only.  The models have no missing regressor, because a missing value inside a
unit leaves a hole in the residuals that ``xttest1`` reads through the calendar of the data.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel import PyxtpanelError
from pyxtpanel.diagnostics.effects.baltagi_li import baltagi_li
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_effects"

#: statistic key -> (name of the Stata scalar of the value, of the p-value; ``None`` = not returned)
KEYS = {
    "joint": ("jt", "jtp"),
    "re": ("re", None),  # xttest1 returns no p-value for the marginal random-effects test
    "re_adjusted": ("rea", "reap"),
    "serial": ("se", "sep"),
    "serial_adjusted": ("sea", "seap"),
    "re_onesided": ("reo", "reop"),
    "re_adjusted_onesided": ("reao", "reaop"),
}
#: tag -> (dataset, pooled model)
SPECS = {
    "bl_effect": ("sim_iv_valid", "y x z1 z2"),
    "bl_unbalanced": ("sim_entry_exit", "x t"),
    "bl_no_effect": ("sim_balanced", "x t"),
    "bl_serial": ("sim_ur_stationary", "y t"),
}

_XTTEST1 = {
    "jt": "r(jt_lm)",
    "jtp": "r(jt_lm_p)",
    "re": "r(rets_lm)",
    "rea": "r(rets_alm)",
    "reap": "r(rets_alm_p)",
    "se": "r(sc_lm)",
    "sep": "r(sc_lm_p)",
    "sea": "r(sc_alm)",
    "seap": "r(sc_alm_p)",
    "reo": "r(reos_lm)",
    "reop": "r(reos_lm_p)",
    "reao": "r(reos_alm)",
    "reaop": "r(reos_alm_p)",
}

#: the same statistics from the residuals, unit by unit, from the definitions
_MATA = [
    "mata:",
    'e = st_data(., "ehat")',
    'id = st_data(., "id")',
    "info = panelsetup(id, 1)",
    "G = rows(info)",
    "m = rows(e)",
    "S = 0",
    "a = 0",
    "num = 0",
    "den = 0",
    "for (g = 1; g <= G; g++) {",
    "    ei = panelsubmatrix(e, g, info)",
    "    Ti = rows(ei)",
    "    S = S + sum(ei)^2",
    "    a = a + Ti^2",
    "    if (Ti > 1) {",
    "        num = num + sum(ei[|2 \\ Ti|] :* ei[|1 \\ Ti - 1|])",
    "        den = den + sum(ei[|2 \\ Ti|] :^ 2)",
    "    }",
    "}",
    "A = 1 - S / sum(e :^ 2)",
    "B = num / den",
    "lre = m^2 * A^2 / (2 * (a - m))",
    "lse = m^2 * B^2 / (m - G)",
    "lrea = m^2 * (A + 2 * B)^2 / (2 * (a - 3 * m + 2 * G))",
    "lsea = lse - lre + lrea",
    'st_numscalar("jt", lre + lsea)',
    'st_numscalar("re", lre)',
    'st_numscalar("rea", lrea)',
    'st_numscalar("se", lse)',
    'st_numscalar("sea", lsea)',
    'st_numscalar("reo", -sqrt(m^2 / (2 * (a - m))) * A)',
    'st_numscalar("reao", -sqrt(m^2 / (2 * (a - 3 * m + 2 * G))) * (A + 2 * B))',
    "end",
    "scalar jtp = chi2tail(2, jt)",
    "scalar reap = chi2tail(1, rea)",
    "scalar sep = chi2tail(1, se)",
    "scalar seap = chi2tail(1, sea)",
    "scalar reop = normal(-reo)",
    "scalar reaop = normal(-reao)",
]


def _stata(tag: str, model: str, reference: str) -> list[str]:
    lines = [
        "quietly xtset id t",
        f"quietly regress {model}",
        "quietly predict double ehat, resid",
    ]
    if reference == "xttest1":
        lines += [f"quietly xtreg {model}, re", "quietly xttest1, unadjusted"]
        lines += [f"scalar {name} = {source}" for name, source in _XTTEST1.items()]
    else:
        lines += _MATA
    return [
        *lines,
        "preserve",
        "quietly keep if !missing(ehat)",
        export(tag, "e", 'st_data(., "ehat")'),
        export(tag, "id", 'st_data(., "id")'),
        export(tag, "t", 'st_data(., "t")'),
        "restore",
        export_scalars(tag, "r", _scalars(reference)),
    ]


def _arrays(tag: str) -> tuple[dict[str, np.ndarray], np.ndarray]:
    s = stata_objects(DOMAIN, tag)
    ids = s["id"].ravel()
    starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])
    return s, np.append(starts, ids.size)


def _scalars(reference: str) -> list[str]:
    """Stata scalars compared: ``xttest1``'s one-sided p-values are numerically poor, left out."""
    skip = {"reop", "reaop"} if reference == "xttest1" else set()
    return [s for pair in KEYS.values() for s in pair if s and s not in skip]


def _register(tag: str, dataset: str, model: str, reference: str) -> None:
    names = _scalars(reference)

    def run() -> dict[str, np.ndarray]:
        s, offsets = _arrays(tag)
        e, t = s["e"].ravel(), s["t"].ravel()
        values: list[float] = []
        for key, (name, p_name) in KEYS.items():
            result = baltagi_li(e, offsets, statistic=key, time_index=t)
            if name in names:
                values.append(result.statistic)
            if p_name in names:
                values.append(result.pvalue)
        return {"r": np.array(values), "e": s["e"], "id": s["id"], "t": s["t"]}

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=_stata(tag, model, reference),
        atol=1e-9,
        rtol=1e-5 if reference == "xttest1" else 1e-9,
    )(run)


for _tag, (_dataset, _model) in SPECS.items():
    _register(_tag, _dataset, _model, "xttest1")
    _register(_tag.replace("bl_", "bl_mata_", 1), _dataset, _model, "mata")


@case(
    "bl_gap_refused",
    domain=DOMAIN,
    dataset="sim_gaps",
    stata=[
        "quietly xtset id t",
        "quietly regress x z",
        "quietly predict double ehat, resid",
        "quietly xtreg x z, re",
        "capture quietly xttest1, unadjusted",
        "scalar rc = _rc",
        "preserve",
        "quietly keep if !missing(ehat)",
        export("bl_gap_refused", "e", 'st_data(., "ehat")'),
        export("bl_gap_refused", "id", 'st_data(., "id")'),
        export("bl_gap_refused", "t", 'st_data(., "t")'),
        "restore",
        export_scalars("bl_gap_refused", "rc", ["rc"]),
    ],
)
def _gap_refused() -> dict[str, np.ndarray]:
    s, offsets = _arrays("bl_gap_refused")
    try:
        baltagi_li(s["e"].ravel(), offsets, time_index=s["t"].ravel())
    except PyxtpanelError as error:
        code = float(error.stata_code or 0)
    else:
        code = 0.0
    return {"rc": np.array([[code]]), "e": s["e"], "id": s["id"], "t": s["t"]}
