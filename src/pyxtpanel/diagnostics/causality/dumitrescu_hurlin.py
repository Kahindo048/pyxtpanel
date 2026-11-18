"""Dumitrescu-Hurlin panel Granger non-causality test (category A, L.5.2).

Validated against Stata 17 through an independent computation in Stata (``xtgcause`` is SSC and
is not installed): for every unit Stata fits ``regress y L(1/K).y L(1/K).x`` and the Wald test of
the lags of ``x`` with ``test``, then Mata forms the aggregate statistics and their p-values; the
statistics are also checked against ``plm::pgrangertest`` of R.  The cases are in
``validation/cases/diagnostics_causality_dumitrescu_hurlin.py``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.dumitrescu_hurlin.title": {
            "en": "Dumitrescu-Hurlin panel Granger non-causality test",
            "fr": "Test de non-causalité de Granger en panel de Dumitrescu-Hurlin",
        },
        "diag.dumitrescu_hurlin.null": {
            "en": "the cause does not Granger-cause the variable in any unit",
            "fr": "la cause ne cause pas la variable au sens de Granger dans aucune unité",
        },
        "diag.dumitrescu_hurlin.alternative": {
            "en": "the cause Granger-causes the variable in at least one unit",
            "fr": "la cause cause la variable au sens de Granger dans au moins une unité",
        },
        "diag.dumitrescu_hurlin.bad_lags": {
            "en": "lags must be a positive integer, got {value!r}.",
            "fr": "lags doit être un entier strictement positif, reçu {value!r}.",
        },
        "diag.dumitrescu_hurlin.bad_statistic": {
            "en": "statistic must be 'ztilde', 'zbar' or 'wbar', got {value!r}.",
            "fr": "statistic doit valoir 'ztilde', 'zbar' ou 'wbar', reçu {value!r}.",
        },
        "diag.dumitrescu_hurlin.same_variable": {
            "en": "the dependent variable and the cause must be different variables.",
            "fr": "la variable dépendante et la cause doivent être des variables différentes.",
        },
        "diag.dumitrescu_hurlin.perfect_fit": {
            "en": "unit {unit} fits perfectly with the lags of both variables: its Wald "
            "statistic is not defined.",
            "fr": "l'unité {unit} est ajustée parfaitement par les retards des deux variables : "
            "sa statistique de Wald n'est pas définie.",
        },
        "diag.dumitrescu_hurlin.collinear": {
            "en": "the regressors of unit {unit} are collinear: its Wald statistic is not defined.",
            "fr": "les régresseurs de l'unité {unit} sont colinéaires : sa statistique de Wald "
            "n'est pas définie.",
        },
        "diag.dumitrescu_hurlin.no_pvalue": {
            "en": "no p-value for W-bar: its exact null distribution is not tabulated; use "
            "statistic='ztilde' or 'zbar'.",
            "fr": "pas de p-value pour W-bar : sa loi nulle exacte n'est pas tabulée ; utilisez "
            "statistic='ztilde' ou 'zbar'.",
        },
    },
)

#: Relative tolerance on the residual sum of squares below which a unit fits perfectly.
PERFECT_FIT = 1e-24


@register_test(
    "dumitrescu_hurlin",
    TestSpec(
        family="causality",
        category="A",
        wave=2,
        reference="Dumitrescu, E.-I. and Hurlin, C. (2012), Economic Modelling 29(4), 1450-1460; "
        "Lopez, L. and Weber, S. (2017), Stata Journal 17(4), 972-984",
        stata="mata",
        validation=(
            "dh_lag1",
            "dh_lag2",
            "dh_lag3",
            "dh_no_causality",
            "dh_no_causality_lag2",
        ),
        inputs="panel",
        balanced="required",
        aliases=("granger", "xtgcause", "pgrangertest"),
        r="plm::pgrangertest",
        summary="Average of the unit Wald statistics of Granger non-causality.",
    ),
)
def dumitrescu_hurlin(
    panel: PanelBase,
    dependent: str,
    cause: str,
    *,
    lags: int = 1,
    statistic: Literal["ztilde", "zbar", "wbar"] = "ztilde",
) -> TestResult:
    """Return the Dumitrescu-Hurlin statistic of ``cause`` not Granger-causing ``dependent``.

    For each unit, ``dependent`` is regressed on a constant and ``lags`` lags of itself and of
    ``cause``; ``W_i`` is the Wald statistic of the joint nullity of the lags of ``cause`` (the
    increase of the residual sum of squares over the unrestricted variance estimate, a chi-square
    with ``lags`` degrees of freedom for a large ``T``).  With ``N`` units, ``T`` periods and ``K``
    lags, the statistics built on the average ``W-bar`` are:

    * ``"wbar"``: ``W-bar`` itself (no p-value);
    * ``"zbar"``: ``sqrt(N / 2K) (W-bar - K)``, standard normal for ``T`` then ``N`` large;
    * ``"ztilde"`` (default): the same standardised with the exact finite-sample mean and variance
      of ``W_i``, ``sqrt(N/2K * (T-3K-5)/(T-2K-3)) ((T-3K-3)/(T-3K-1) W-bar - K)``, standard
      normal for a fixed ``T > 3K + 5`` and ``N`` large.

    The null is that ``cause`` does not Granger-cause ``dependent`` in any unit; the alternative
    is causality in at least one unit.  p-values are two-sided.  The panel must be strongly
    balanced and without gaps.  ``details`` gives each unit's Wald statistic and its chi-square
    p-value.

    Parameters
    ----------
    panel
        The declared panel.
    dependent
        Name of the variable that is caused.
    cause
        Name of the variable whose lags are tested.
    lags
        Number of lags ``K`` of both variables (the same in every unit).
    statistic
        The statistic reported as the result's statistic and p-value.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If ``lags`` or ``statistic`` is invalid or the two variables are the same.
    DataError
        If the panel is not strongly balanced, has a gap or is too short.
    NumericalError
        If a unit's regression is degenerate.
    """
    # Stata: independent computation, per unit `regress y L(1/K).y L(1/K).x`, W_i = K * F of
    # `test L(1/K).x`; Wbar = mean W_i; Zbar = sqrt(N/(2K)) (Wbar - K); Ztilde = sqrt(N/(2K)
    # (T-3K-5)/(T-2K-3)) ((T-3K-3)/(T-3K-1) Wbar - K); p = 2 normal(-|z|).
    if isinstance(lags, bool) or not isinstance(lags, (int, np.integer)) or lags < 1:
        raise SpecificationError(
            _i18n.message("diag.dumitrescu_hurlin.bad_lags", value=lags),
            key="diag.dumitrescu_hurlin.bad_lags",
        )
    if statistic not in ("ztilde", "zbar", "wbar"):
        raise SpecificationError(
            _i18n.message("diag.dumitrescu_hurlin.bad_statistic", value=statistic),
            key="diag.dumitrescu_hurlin.bad_statistic",
        )
    if dependent == cause:
        raise SpecificationError(
            _i18n.message("diag.dumitrescu_hurlin.same_variable"),
            key="diag.dumitrescu_hurlin.same_variable",
        )
    k = int(lags)
    arrays = panel_variables(
        panel,
        [dependent, cause],
        test="dumitrescu_hurlin",
        balanced="required",
        consecutive=True,
        min_obs=3 * k + (6 if statistic == "ztilde" else 2),
        min_units=2,
    )
    cube = arrays.as_cube()  # (units, periods, 2)
    n_units, periods, _ = cube.shape
    walds = np.array(
        [
            _unit_wald(cube[unit, :, 0], cube[unit, :, 1], k, arrays.unit_labels[unit])
            for unit in range(n_units)
        ]
    )
    wbar = float(walds.mean())
    value = _aggregate(statistic, wbar, n_units, periods, k)
    pvalue = math.nan if statistic == "wbar" else float(2.0 * stats.norm.sf(abs(value)))
    details = pd.DataFrame(
        {"wald": walds, "pvalue": stats.chi2.sf(walds, k), "n_obs": periods - k},
        index=pd.Index(arrays.unit_labels, name="unit"),
    )
    notes = (Note("diag.dumitrescu_hurlin.no_pvalue"),) if statistic == "wbar" else ()
    return TestResult(
        name="dumitrescu_hurlin",
        statistic=value,
        pvalue=pvalue,
        statistic_name={"ztilde": "Z-tilde", "zbar": "Z-bar", "wbar": "W-bar"}[statistic],
        df=(),
        n_obs=n_units * (periods - k),
        n_groups=n_units,
        variable=dependent,
        options={"cause": cause, "lags": k, "statistic": statistic, "wbar": wbar},
        details=details,
        notes=notes,
    )


def _aggregate(statistic: str, wbar: float, n_units: int, periods: int, k: int) -> float:
    """Return ``W-bar``, ``Z-bar`` or ``Z-tilde`` for ``N`` units, ``T`` periods and ``K`` lags."""
    if statistic == "wbar":
        return wbar
    if statistic == "zbar":
        return math.sqrt(n_units / (2.0 * k)) * (wbar - k)
    t = periods
    scale = n_units / (2.0 * k) * (t - 3.0 * k - 5.0) / (t - 2.0 * k - 3.0)
    return math.sqrt(scale) * ((t - 3.0 * k - 3.0) / (t - 3.0 * k - 1.0) * wbar - k)


def _unit_wald(y: FloatArray, x: FloatArray, k: int, label: object) -> float:
    """Return the Wald statistic ``W_i`` of the ``k`` lags of ``x`` in the regression of ``y``."""
    n = y.size - k
    response = y[k:]
    own = [y[k - j : y.size - j] for j in range(1, k + 1)]
    other = [x[k - j : x.size - j] for j in range(1, k + 1)]
    restricted = np.column_stack([np.ones(n), *own])
    full = np.column_stack([restricted, *other])
    rss_u, rank_u = _rss(response, full)
    rss_r, _ = _rss(response, restricted)
    if rank_u < full.shape[1]:
        raise NumericalError(
            _i18n.message("diag.dumitrescu_hurlin.collinear", unit=label),
            key="diag.dumitrescu_hurlin.collinear",
        )
    if rss_u <= PERFECT_FIT * float(response @ response + 1.0):
        raise NumericalError(
            _i18n.message("diag.dumitrescu_hurlin.perfect_fit", unit=label),
            key="diag.dumitrescu_hurlin.perfect_fit",
        )
    return (rss_r - rss_u) / (rss_u / (n - 2 * k - 1))


def _rss(y: FloatArray, x: FloatArray) -> tuple[float, int]:
    """Return the residual sum of squares of the regression of ``y`` on ``x`` and its rank."""
    coefficients, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ coefficients
    return float(residual @ residual), int(rank)
