"""Baltagi-Li tests of serial correlation and random effects, from bare arrays (L.8.14).

The Lagrange multiplier tests of Baltagi and Li (1991) and their robust versions of Bera,
Sosa-Escudero and Yoon (2001), for unbalanced panels (Sosa-Escudero and Bera, 2008), all built
from the residuals of the pooled least-squares regression.  Validated against Stata 17:
``xttest1, unadjusted`` after ``xtreg, re`` (the update ``sg164_1`` of the Stata Journal, installed
for the purpose) and an independent Mata computation; the cases are in
``validation/cases/diagnostics_effects_baltagi_li.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, as_array, register_test

_i18n.register(
    "diag",
    {
        "diag.baltagi_li.title": {
            "en": "Baltagi-Li tests of serial correlation and random effects",
            "fr": "Tests de Baltagi-Li d'autocorrélation et d'effets aléatoires",
        },
        "diag.baltagi_li.null": {
            "en": "Var(u) = 0 and lambda = 0 (joint test); each other statistic tests the part "
            "named by its key",
            "fr": "Var(u) = 0 et lambda = 0 (test joint) ; chacune des autres statistiques "
            "teste la partie que nomme sa clé",
        },
        "diag.baltagi_li.alternative": {
            "en": "random effects, first-order serial correlation of the errors, or both",
            "fr": "effets aléatoires, autocorrélation d'ordre un des erreurs, ou les deux",
        },
        "diag.baltagi_li.bad_statistic": {
            "en": "statistic must be one of {allowed}, got {value!r}.",
            "fr": "statistic doit valoir l'une de {allowed}, reçu {value!r}.",
        },
        "diag.baltagi_li.offsets": {
            "en": "offsets must start at 0, increase and end at the number of residuals ({n}).",
            "fr": "offsets doit commencer à 0, croître et finir au nombre de résidus ({n}).",
        },
        "diag.baltagi_li.time": {
            "en": "time_index must have one entry per residual ({n}), got {got}.",
            "fr": "time_index doit avoir une entrée par résidu ({n}), reçu {got}.",
        },
        "diag.baltagi_li.gap": {
            "en": "the panel has a gap (unit {unit} skips a period after row {row}): the test "
            "needs consecutive observations, as xttest1 (r(498)).",
            "fr": "le panel a un trou (l'unité {unit} saute une période après la ligne {row}) : "
            "le test exige des observations consécutives, comme xttest1 (r(498)).",
        },
        "diag.baltagi_li.zero_residuals": {
            "en": "all residuals are zero: the tests are not defined.",
            "fr": "tous les résidus sont nuls : les tests ne sont pas définis.",
        },
        "diag.baltagi_li.undefined": {
            "en": "the statistic {statistic!r} is not defined for these unit sizes: {need}.",
            "fr": "la statistique {statistic!r} n'est pas définie pour ces tailles d'unités : "
            "{need}.",
        },
        "diag.baltagi_li.need_repeated": {
            "en": "some unit needs two observations",
            "fr": "une unité au moins doit avoir deux observations",
        },
        "diag.baltagi_li.need_three": {
            "en": "some unit needs three observations",
            "fr": "une unité au moins doit avoir trois observations",
        },
    },
)

#: key -> (statistic name, number of degrees of freedom or ``None`` for a standard normal)
_STATISTICS: dict[str, tuple[str, int | None]] = {
    "joint": ("chi2", 2),
    "re": ("chi2", 1),
    "re_adjusted": ("chi2", 1),
    "serial": ("chi2", 1),
    "serial_adjusted": ("chi2", 1),
    "re_onesided": ("z", None),
    "re_adjusted_onesided": ("z", None),
}


@register_test(
    "baltagi_li",
    TestSpec(
        family="effects",
        category="B",
        wave=2,
        reference="Baltagi, B. H. and Li, Q. (1991), Economics Letters 37(3), 277-280; "
        "Bera, A. K., Sosa-Escudero, W. and Yoon, M. (2001), J. of Econometrics 101(1), 123-137; "
        "Sosa-Escudero, W. and Bera, A. K. (2008), Stata Journal 8(1), 134-135",
        stata="xttest1; mata",
        validation=(
            "bl_effect",
            "bl_unbalanced",
            "bl_no_effect",
            "bl_serial",
            "bl_mata_effect",
            "bl_mata_unbalanced",
            "bl_mata_no_effect",
            "bl_mata_serial",
            "bl_gap_refused",
        ),
        inputs="arrays",
        balanced="na",
        aliases=("xttest1", "pbsytest"),
        r="plm::pbsytest",
        summary="Joint, marginal and robust LM tests of random effects and AR(1) errors.",
    ),
)
def baltagi_li(
    residuals: object,
    offsets: object,
    *,
    statistic: str = "joint",
    time_index: object | None = None,
) -> TestResult:
    """Return one of the Baltagi-Li statistics of serial correlation and random effects.

    With ``e`` the residuals of the pooled least-squares regression (a constant and the
    regressors), ``m`` their number, ``G`` the units, ``T_i`` the size of unit ``i`` and
    ``S = sum_i T_i^2``::

        A = 1 - sum_i (sum_t e_it)^2 / sum e^2         B = sum e_it e_i,t-1 / sum_{t>1} e_it^2
        LM_re = m^2 A^2 / (2 (S - m))                  LM_serial = m^2 B^2 / (m - G)
        ALM_re = m^2 (A + 2B)^2 / (2 (S - 3m + 2G))    ALM_serial = LM_serial - LM_re + ALM_re
        LM_joint = LM_re + ALM_serial

    ``LM_re`` and ``LM_serial`` are the marginal tests, valid when the other component is absent;
    the ``ALM`` versions are robust to it (they tell random effects from serial correlation).  The
    one-sided tests are ``-sqrt(m^2 / (2 (S - m))) A`` and ``-sqrt(m^2 / (2 (S - 3m + 2G))) (A +
    2B)``, standard normal with the p-value in the right tail.  All seven are in ``details``.

    Parameters
    ----------
    residuals
        ``(m,)`` residuals of the pooled regression, sorted by unit and in time order within it,
        consecutive in time (the lag is the previous row).
    offsets
        ``(G + 1,)`` start of each unit in ``residuals``.
    statistic
        Statistic to return: ``"joint"`` (default, chi-square with 2 degrees of freedom), ``"re"``,
        ``"re_adjusted"``, ``"serial"``, ``"serial_adjusted"`` (chi-square with 1) or
        ``"re_onesided"``, ``"re_adjusted_onesided"`` (standard normal).
    time_index
        Optional ``(m,)`` integer period of each residual: a unit that skips a period is then
        refused, as ``xttest1`` does (r(498)).

    Raises
    ------
    SpecificationError
        If ``statistic`` is unknown or ``offsets`` or ``time_index`` do not describe the residuals.
    DataError
        If ``time_index`` shows a gap.
    NumericalError
        If all residuals are zero or the unit sizes leave the statistic undefined (no unit with
        two observations for the serial tests, none with three for the adjusted ones).
    """
    # STATA-QUIRK: xttest1 stores the lagged residual as a float, so its statistics agree with the
    # double-precision ones to about 2e-6 (relative), and returns the one-sided p-values as
    # 1 - normprob(z), which loses digits from z = 5 and is 0 beyond z = 8; these are exact here
    # Stata: xttest1, unadjusted (sg164_1) -> from the residuals of `regress` over e(sample):
    # A1 = 1 - sum_i(sum_t e)^2 / sum(e^2), B1 = sum(e e_lag) / sum_{t>1}(e^2), a = sum(T_i^2),
    # m = N, nn = G; LM1 = .5 m^2 A1^2 / (a - m), LM2 = m^2 B1^2 / (m - nn),
    # LM3 = m^2 (A1 + 2 B1)^2 / (2 (a - 3m + 2nn)), LM4 = LM2 - LM1 + LM3, LM5 = LM1 + LM4,
    # LM6 = -sqrt(.5 m^2 / (a - m)) A1, LM7 = -sqrt(m^2 / (2 (a - 3m + 2nn))) (A1 + 2 B1);
    # chi2 tails, and 1 - normal(LM) for the two one-sided ones; panel with gaps: r(498)
    if statistic not in _STATISTICS:
        raise SpecificationError(
            _i18n.message(
                "diag.baltagi_li.bad_statistic", allowed=list(_STATISTICS), value=statistic
            ),
            key="diag.baltagi_li.bad_statistic",
        )
    e = as_array("residuals", residuals, ndim=1)
    bounds = _check_offsets(e.size, offsets)
    if time_index is not None:
        _check_consecutive(e.size, bounds, time_index)
    sizes = np.diff(bounds).astype(np.float64)
    m, groups = np.float64(e.size), np.float64(sizes.size)
    sum_sq = float(e @ e)
    if sum_sq == 0.0:
        raise NumericalError(
            _i18n.message("diag.baltagi_li.zero_residuals"), key="diag.baltagi_li.zero_residuals"
        )
    s = np.float64(np.sum(sizes**2))
    a1 = 1.0 - np.sum(np.add.reduceat(e, bounds[:-1]) ** 2) / sum_sq
    first = np.zeros(e.size, dtype=bool)
    first[bounds[:-1]] = True
    pairs = ~first
    lagged = np.concatenate([[0.0], e[:-1]])
    denominator = float(e[pairs] @ e[pairs])
    b1 = (e[pairs] @ lagged[pairs]) / denominator if denominator > 0.0 else np.float64(math.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        lm_re = 0.5 * m**2 * a1**2 / (s - m)
        lm_serial = m**2 * b1**2 / (m - groups)
        alm_re = m**2 * (a1 + 2.0 * b1) ** 2 / (2.0 * (s - 3.0 * m + 2.0 * groups))
        alm_serial = lm_serial - lm_re + alm_re
        values = {
            "joint": lm_re + alm_serial,
            "re": lm_re,
            "re_adjusted": alm_re,
            "serial": lm_serial,
            "serial_adjusted": alm_serial,
            "re_onesided": -math.sqrt(0.5 * m**2 / (s - m)) * a1 if s > m else math.nan,
            "re_adjusted_onesided": (
                -math.sqrt(m**2 / (2.0 * (s - 3.0 * m + 2.0 * groups))) * (a1 + 2.0 * b1)
                if s - 3.0 * m + 2.0 * groups > 0
                else math.nan
            ),
        }
    chosen = values[statistic]
    if not math.isfinite(chosen):
        need = (
            "diag.baltagi_li.need_three"
            if "adjusted" in statistic or statistic == "joint"
            else ("diag.baltagi_li.need_repeated")
        )
        raise NumericalError(
            _i18n.message(
                "diag.baltagi_li.undefined",
                statistic=statistic,
                need=_i18n.message(need),
            ),
            key="diag.baltagi_li.undefined",
        )
    name, df = _STATISTICS[statistic]
    pvalues = {
        key: float(stats.norm.sf(v))
        if _STATISTICS[key][1] is None
        else float(stats.chi2.sf(v, _STATISTICS[key][1]))
        for key, v in values.items()
    }
    details = pd.DataFrame(
        {
            "statistic": list(values.values()),
            "pvalue": [pvalues[k] for k in values],
            "df": [float(_STATISTICS[k][1] or math.nan) for k in values],
        },
        index=pd.Index(list(values), name="test"),
    )
    return TestResult(
        name="baltagi_li",
        statistic=float(chosen),
        pvalue=pvalues[statistic],
        statistic_name=name,
        df=() if df is None else (float(df),),
        n_obs=e.size,
        n_groups=int(sizes.size),
        options={"statistic": statistic},
        details=details,
    )


def _check_offsets(n: int, offsets: object) -> np.ndarray:
    bounds = np.asarray(offsets, dtype=np.int64)
    if (
        bounds.ndim != 1
        or bounds.size < 2
        or bounds[0] != 0
        or bounds[-1] != n
        or np.any(np.diff(bounds) < 1)
    ):
        raise SpecificationError(
            _i18n.message("diag.baltagi_li.offsets", n=n), key="diag.baltagi_li.offsets"
        )
    return bounds


def _check_consecutive(n: int, bounds: np.ndarray, time_index: object) -> None:
    time = np.asarray(time_index, dtype=np.int64)
    if time.shape != (n,):
        raise SpecificationError(
            _i18n.message("diag.baltagi_li.time", n=n, got=time.size), key="diag.baltagi_li.time"
        )
    step = np.diff(time)
    inside = np.ones(step.size, dtype=bool)
    inside[bounds[1:-1] - 1] = False  # the step between two units is not a step
    bad = np.flatnonzero(inside & (step != 1))
    if bad.size:
        row = int(bad[0])
        unit = int(np.searchsorted(bounds, row, side="right") - 1)
        raise DataError(
            _i18n.message("diag.baltagi_li.gap", unit=unit, row=row),
            stata_code=498,
            key="diag.baltagi_li.gap",
        )
