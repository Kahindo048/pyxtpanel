"""Pedroni residual-based panel cointegration tests (category A, L.4.2).

Validated against Stata 17: ``xtcointtest pedroni`` with panel-specific (group statistics) and
common (panel statistics) autoregressive parameters, with no constant, a constant and a trend,
with one to seven regressors, with the three kernels at a fixed number of lags, with one and
several augmented lags and with lags chosen by information criteria, with and without the
cross-sectional means removed, on balanced and unbalanced panels.  The cases are in
``validation/cases/diagnostics_cointegration_pedroni.py``.  The adjustment terms that standardise
the statistics (means and variances of Pedroni, 1999) are compiled in Stata's Mata library and
are not published by it: they were deduced from its output, see
``pyxtpanel/_data/tables/MANIFEST.json``.  The automatic Newey-West bandwidth of Stata (its
default) is not reproduced: ``kernel_lags`` is required.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.long_run import long_run_covariance
from pyxtpanel.diagnostics._shared.tables import load_table

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.pedroni.title": {
            "en": "Pedroni panel cointegration test",
            "fr": "Test de cointégration en panel de Pedroni",
        },
        "diag.pedroni.null": {
            "en": "no cointegration",
            "fr": "pas de cointégration",
        },
        "diag.pedroni.alternative": {
            "en": "all panels are cointegrated",
            "fr": "toutes les unités sont cointégrées",
        },
        "diag.pedroni.bad_option": {
            "en": "{name} must be one of {allowed}, got {value!r}.",
            "fr": "{name} doit valoir l'une de ces valeurs : {allowed} ; reçu {value!r}.",
        },
        "diag.pedroni.regressors": {
            "en": "give between 1 and 7 distinct regressors different from the dependent "
            "variable, got {value!r}.",
            "fr": "donnez de 1 à 7 régresseurs distincts et différents de la variable "
            "dépendante, reçu {value!r}.",
        },
        "diag.pedroni.bad_lags": {
            "en": "{name} must be {rule}, got {value!r}.",
            "fr": "{name} doit être {rule}, reçu {value!r}.",
        },
        "diag.pedroni.statistic": {
            "en": "with ar={ar!r} the statistic must be one of {allowed}, got {value!r}.",
            "fr": "avec ar={ar!r} la statistique doit valoir l'une de ces valeurs : {allowed} ; "
            "reçu {value!r}.",
        },
        "diag.pedroni.degenerate": {
            "en": "the residuals of unit {unit!r} are degenerate (no variation, a perfect fit or "
            "a singular long-run covariance): the statistics are not defined.",
            "fr": "les résidus de l'unité {unit!r} sont dégénérés (aucune variation, ajustement "
            "parfait ou covariance de long terme singulière) : les statistiques ne sont pas "
            "définies.",
        },
        "diag.pedroni.sign_quirk": {
            "en": "Stata's adjustment term for the panel modified Phillips-Perron t with no "
            "deterministic term and seven regressors is +25.889, while the series of the other "
            "numbers of regressors and Pedroni's table (1999) give -25.889: the sign looks wrong "
            "in Stata, and this result reproduces it.",
            "fr": "le terme d'ajustement de Stata pour le t de Phillips-Perron modifié en panel, "
            "sans terme déterministe et avec sept régresseurs, vaut +25,889, alors que la suite "
            "des autres nombres de régresseurs et la table de Pedroni (1999) donnent -25,889 : le "
            "signe semble faux dans Stata, et ce résultat le reproduit.",
        },
        "diag.pedroni.selection_needs_lags": {
            "en": "adf_lag_selection chooses among 1 to adf_lags lags: give the largest number in "
            "adf_lags.",
            "fr": "adf_lag_selection choisit entre 1 et adf_lags retards : donnez le maximum dans "
            "adf_lags.",
        },
    },
)

_KERNELS = ("bartlett", "parzen", "quadraticspectral")
_FAMILY = {
    "same": ("nu", "rho", "pp", "adf"),
    "panelspecific": ("rho", "pp", "adf"),
}
_NAMES = {
    ("same", "nu"): ("panel_nu", "Panel variance ratio"),
    ("same", "rho"): ("panel_rho", "Panel modified Phillips-Perron t"),
    ("same", "pp"): ("panel_pp", "Panel Phillips-Perron t"),
    ("same", "adf"): ("panel_adf", "Panel ADF t"),
    ("panelspecific", "rho"): ("group_rho", "Group modified Phillips-Perron t"),
    ("panelspecific", "pp"): ("group_pp", "Group Phillips-Perron t"),
    ("panelspecific", "adf"): ("group_adf", "Group ADF t"),
}


@register_test(
    "pedroni",
    TestSpec(
        family="cointegration",
        category="A",
        wave=2,
        reference="Pedroni, P. (1999), Oxford Bulletin of Economics and Statistics 61, 653-670; "
        "Pedroni, P. (2004), Econometric Theory 20(3), 597-625",
        stata="xtcointtest pedroni",
        validation=(
            "pe_group_m1",
            "pe_panel_m1",
            "pe_group_m3_trend",
            "pe_panel_m3_trend",
            "pe_group_m2_none",
            "pe_panel_m2_none",
            "pe_group_m7",
            "pe_panel_m7_trend",
            "pe_parzen",
            "pe_quadraticspectral",
            "pe_kernel0",
            "pe_adf2",
            "pe_adf_aic",
            "pe_adf_bic",
            "pe_adf_hqic",
            "pe_demean",
            "pe_unbalanced",
            "pe_unbalanced_panel_trend",
            "pe_unbalanced_none",
            "pe_coint",
        ),
        inputs="panel",
        balanced="both",
        tables=("pedroni_moments.csv",),
        summary="Residual-based tests of no cointegration with panel-specific cointegrating "
        "vectors.",
    ),
)
def pedroni(
    panel: PanelBase,
    dependent: str,
    regressors: Sequence[str],
    *,
    kernel_lags: int,
    ar: Literal["panelspecific", "same"] = "panelspecific",
    statistic: Literal["nu", "rho", "pp", "adf"] = "adf",
    deterministics: Literal["constant", "none", "trend"] = "constant",
    kernel: Literal["bartlett", "parzen", "quadraticspectral"] = "bartlett",
    adf_lags: int = 1,
    adf_lag_selection: Literal["aic", "bic", "hqic"] | None = None,
    demean: bool = False,
) -> TestResult:
    """Return one of Pedroni's statistics for no cointegration (the others are in ``details``).

    Each unit's dependent variable is regressed on its own constant (and trend) and regressors; the
    residuals ``e`` are tested for a unit root with panel-specific (``ar="panelspecific"``, the
    group statistics ``rho``, ``pp`` and ``adf``) or a common autoregressive parameter
    (``ar="same"``, the panel statistics ``nu``, ``rho``, ``pp`` and ``adf``), the latter weighting
    each unit by the inverse of its long-run conditional variance ``L11`` of the differenced
    series.  The non-parametric statistics correct the numerator by ``lambda = (sigma^2 - s^2) / 2``
    of the autoregression of ``e`` without constant (long-run variance with ``kernel`` and
    ``kernel_lags``); the ADF statistics use ``adf_lags`` lagged differences.  Each statistic is
    standardised as ``(Z - mean sqrt(N)) / sqrt(variance)`` with the adjustment terms of Pedroni
    (1999) for the number of regressors and the deterministic terms, and is standard normal under
    the null.  The ``nu`` statistic rejects in the right tail, the six others in the left tail.

    The panel may be unbalanced but must have no gaps; ``T`` in the formulas is the longest unit.
    At most seven regressors.  The automatic bandwidth of Stata is not available: give
    ``kernel_lags``.

    Parameters
    ----------
    panel
        The declared panel.
    dependent
        Name of the dependent variable.
    regressors
        Names of the regressors (1 to 7).
    kernel_lags
        Number of lags of the kernel in the long-run variances (``0`` or more).
    ar
        ``"panelspecific"`` (group statistics) or ``"same"`` (panel statistics).
    statistic
        The statistic returned as ``result.statistic``: ``"rho"``, ``"pp"`` or ``"adf"``, and
        ``"nu"`` with ``ar="same"``.
    deterministics
        ``"constant"`` (default), ``"none"`` or ``"trend"`` (constant and trend).
    kernel
        ``"bartlett"`` (default), ``"parzen"`` or ``"quadraticspectral"``.
    adf_lags
        Number of lagged differences of the ADF regressions (at least 1), or their maximum with
        ``adf_lag_selection``.
    adf_lag_selection
        ``"aic"``, ``"bic"`` or ``"hqic"``: choose the number of lags between 1 and ``adf_lags``
        on the pooled regression.
    demean
        Remove the mean across units in each period from every variable first.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If an option is invalid or the regressors are not 1 to 7 distinct names.
    DataError
        If a unit has a gap or too few observations.
    NumericalError
        If a unit's residuals are degenerate.
    """
    # Stata: xtcointtest pedroni y x [, noconstant trend kernel(k #) lags(#) ar(same) demean] ->
    # e = per-unit residuals of y on 1 [t] x; L11 = long-run conditional variance of (dy, dx)
    # (centred, Sigma/Omega of xtlrcov, bandwidth lags + 1); DF on e without constant: u,
    # s2 = centred variance, sigma2 = centred long-run variance, lambda = (sigma2 - s2) / 2,
    # numerator sum_t (e_{t-1} de_t - lambda); ADF: pooled lag choice, per unit regression of de
    # on e_{t-1}, de lags without constant, s2 centred; raw statistics (T = longest unit):
    # panel nu T^2 N^1.5 / sum e2/L11, rho T sqrt(N) sum num/L11 / sum e2/L11, pp, adf; group
    # rho T / sqrt(N) sum num/e2, pp, adf; z = (Z - mean sqrt(N)) / sqrt(variance).
    # STATA-QUIRK: T is the longest unit on an unbalanced panel; every p-value is reported as
    # 1 - normal(|z|), the left tail normal(z) (right tail for nu) is used here; the table of
    # adjustment terms has +25.889 (not -25.889) for the panel rho, no constant, 7 regressors.
    names = _names(dependent, regressors)
    _validate(kernel_lags, ar, statistic, deterministics, kernel, adf_lags, adf_lag_selection)
    trend = deterministics == "trend"
    constant = deterministics != "none"
    arrays = panel_variables(
        panel,
        [dependent, *names],
        test="pedroni",
        consecutive=True,
        min_obs=_min_periods(
            len(names), int(constant) + int(trend), int(adf_lags), int(kernel_lags)
        ),
        min_units=2,
    )
    values = arrays.values.copy()
    if demean:
        counts = np.bincount(arrays.time_index)
        for column in range(values.shape[1]):
            totals = np.bincount(arrays.time_index, weights=values[:, column])
            values[:, column] -= (totals / np.where(counts > 0, counts, 1))[arrays.time_index]
    pieces = _unit_pieces(
        values, arrays.unit_slices(), arrays.unit_labels, deterministics, kernel, int(kernel_lags)
    )
    p = _selected_lags([unit.e for unit in pieces], int(adf_lags), adf_lag_selection)
    adf = [_adf_unit(unit, p, arrays.unit_labels[i]) for i, unit in enumerate(pieces)]
    raw = _raw_statistics(pieces, adf, float(np.max(arrays.counts)))
    n_units = len(pieces)
    table = load_table("pedroni_moments.csv")
    rows = {}
    for key in _FAMILY[ar]:
        label, title = _NAMES[(ar, key)]
        found = table[
            (table["deterministics"] == deterministics)
            & (table["m"] == len(names))
            & (table["statistic"] == label)
        ].iloc[0]
        z = (raw[label] - float(found["mean"]) * math.sqrt(n_units)) / math.sqrt(
            float(found["variance"])
        )
        pvalue = float(stats.norm.sf(z)) if key == "nu" else float(stats.norm.cdf(z))
        rows[title] = {"statistic": z, "pvalue": pvalue, "raw": raw[label]}
    label, title = _NAMES[(ar, statistic)]
    details = pd.DataFrame.from_dict(rows, orient="index")
    details.index.name = "statistic"
    return TestResult(
        name="pedroni",
        statistic=float(details.loc[title, "statistic"]),
        pvalue=float(details.loc[title, "pvalue"]),
        statistic_name=title,
        df=(),
        n_obs=arrays.n_obs - n_units,  # the rows left after differencing, as r(N)
        n_groups=n_units,
        variable=dependent,
        options={
            "regressors": list(names),
            "ar": ar,
            "statistic": statistic,
            "deterministics": deterministics,
            "kernel": kernel,
            "kernel_lags": int(kernel_lags),
            "adf_lags": p,
            "adf_lag_selection": adf_lag_selection,
            "demean": demean,
        },
        details=details,
        notes=(
            (Note("diag.pedroni.sign_quirk"),)
            if (ar, deterministics, len(names)) == ("same", "none", 7)
            else ()
        ),
    )


class _Unit:
    """The quantities of one unit that the statistics are built from."""

    __slots__ = ("e", "l11", "n", "num", "sigma2", "sxx")

    def __init__(
        self, e: FloatArray, l11: float, num: float, sxx: float, sigma2: float, n: int
    ) -> None:
        self.e = e
        self.l11 = l11
        self.num = num
        self.sxx = sxx
        self.sigma2 = sigma2
        self.n = n


def _names(dependent: str, regressors: Sequence[str]) -> list[str]:
    names = [regressors] if isinstance(regressors, str) else list(regressors)
    if not 1 <= len(names) <= 7 or len(set(names)) != len(names) or dependent in names:
        raise SpecificationError(
            _i18n.message("diag.pedroni.regressors", value=regressors),
            key="diag.pedroni.regressors",
        )
    return names


def _min_periods(m: int, deterministic: int, adf_lags: int, kernel_lags: int) -> int:
    """Fewest periods a unit needs: the cointegrating, differenced, ADF and kernel regressions."""
    return max(m + deterministic + 1, m + 2, 2 * adf_lags + 3, kernel_lags + 2)


def _is_count(value: object, minimum: int) -> bool:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        return False
    return int(value) >= minimum


def _validate(
    kernel_lags: object,
    ar: object,
    statistic: object,
    deterministics: object,
    kernel: object,
    adf_lags: object,
    adf_lag_selection: object,
) -> None:
    for name, value, allowed in (
        ("ar", ar, ("panelspecific", "same")),
        ("deterministics", deterministics, ("constant", "none", "trend")),
        ("kernel", kernel, _KERNELS),
        ("adf_lag_selection", adf_lag_selection, (None, "aic", "bic", "hqic")),
    ):
        if value not in allowed:
            raise SpecificationError(
                _i18n.message(
                    "diag.pedroni.bad_option",
                    name=name,
                    allowed=", ".join(map(repr, allowed)),
                    value=value,
                ),
                key="diag.pedroni.bad_option",
            )
    if statistic not in _FAMILY[str(ar)]:
        raise SpecificationError(
            _i18n.message(
                "diag.pedroni.statistic",
                ar=ar,
                allowed=", ".join(map(repr, _FAMILY[str(ar)])),
                value=statistic,
            ),
            key="diag.pedroni.statistic",
        )
    if not _is_count(kernel_lags, 0):
        raise SpecificationError(
            _i18n.message(
                "diag.pedroni.bad_lags",
                name="kernel_lags",
                rule="a non-negative integer (the automatic bandwidth of Stata is not available)",
                value=kernel_lags,
            ),
            key="diag.pedroni.bad_lags",
        )
    if not _is_count(adf_lags, 1):
        raise SpecificationError(
            _i18n.message(
                "diag.pedroni.bad_lags",
                name="adf_lags",
                rule="an integer of at least 1",
                value=adf_lags,
            ),
            key="diag.pedroni.bad_lags",
        )


def _plain(label: object) -> object:
    return label.item() if isinstance(label, np.generic) else label


def _degenerate(label: object) -> NumericalError:
    return NumericalError(
        _i18n.message("diag.pedroni.degenerate", unit=_plain(label)), key="diag.pedroni.degenerate"
    )


def _unit_pieces(
    values: FloatArray,
    windows: list[slice],
    labels: Any,
    deterministics: str,
    kernel: str,
    lags: int,
) -> list[_Unit]:
    """Residuals, ``L11`` and the Dickey-Fuller quantities of each unit."""
    pieces = []
    for unit, window in enumerate(windows):
        block = values[window]
        y, x = block[:, 0], block[:, 1:]
        t = y.size
        columns = []
        if deterministics != "none":
            columns.append(np.ones(t))
        if deterministics == "trend":
            columns.append(np.arange(1.0, t + 1.0))
        design = np.column_stack([*columns, x])
        e = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
        w = np.column_stack([np.diff(y), np.diff(x, axis=0)])
        long_run = long_run_covariance(w[None, :, :], kernel, lags)
        lag, cur = e[:-1], e[1:]
        sxx = float(lag @ lag)
        if long_run is None or sxx <= 0.0:
            raise _degenerate(labels[unit])
        omega = long_run.omega
        try:
            l11 = float(omega[0, 0] - omega[0, 1:] @ np.linalg.solve(omega[1:, 1:], omega[1:, 0]))
        except np.linalg.LinAlgError as error:
            raise _degenerate(labels[unit]) from error
        u = cur - (float(lag @ cur) / sxx) * lag
        scalar = long_run_covariance(u[None, :, None], kernel, lags)
        if scalar is None or not l11 > 0.0:
            raise _degenerate(labels[unit])
        s2, sigma2 = float(scalar.sigma[0, 0]), float(scalar.omega[0, 0])
        correction = 0.5 * (sigma2 - s2)
        num = float(lag @ (cur - lag)) - lag.size * correction
        pieces.append(_Unit(e, l11, num, sxx, sigma2, lag.size))
    return pieces


def _adf_rows(e: FloatArray, p: int, first: int) -> tuple[FloatArray, FloatArray]:
    """Response ``e_t`` and regressors (lagged level, ``p`` lagged differences), from ``first``."""
    diff = np.diff(e)
    rows = np.arange(first, e.size)
    columns = [e[rows - 1]]
    columns.extend(diff[rows - 1 - j] for j in range(1, p + 1))
    return e[rows], np.column_stack(columns)


def _selected_lags(units: list[FloatArray], maximum: int, criterion: str | None) -> int:
    """Return ``maximum`` or the lag count that minimises ``criterion`` on the pooled regression."""
    if criterion is None:
        return maximum
    best, best_value = 1, math.inf
    for p in range(1, maximum + 1):
        responses, designs = zip(*(_adf_rows(e, p, maximum + 1) for e in units), strict=True)
        y, x = np.concatenate(responses), np.vstack(designs)
        coefficients, *_ = np.linalg.lstsq(x, y, rcond=None)
        n = y.size
        rss = float(np.sum((y - x @ coefficients) ** 2))
        log_likelihood = -0.5 * n * (math.log(2.0 * math.pi) + 1.0 + math.log(rss / n))
        penalty = {"aic": 2.0, "bic": math.log(n), "hqic": 2.0 * math.log(math.log(n))}[criterion]
        value = -2.0 * log_likelihood + penalty * (p + 1)
        if value < best_value:
            best, best_value = p, value
    return best


def _adf_unit(unit: _Unit, p: int, label: object) -> tuple[float, float, float]:
    """Return ``sum e_{t-1} de_t``, ``sum e_{t-1}^2`` and the centred residual variance."""
    e = unit.e
    diff = np.diff(e)
    rows = np.arange(p + 1, e.size)
    response = diff[rows - 1]
    level = e[rows - 1]
    design = np.column_stack([level, *(diff[rows - 1 - j] for j in range(1, p + 1))])
    residual = response - design @ np.linalg.lstsq(design, response, rcond=None)[0]
    centred = residual - residual.mean()
    s2 = float(centred @ centred) / residual.size
    if not s2 > 0.0:
        raise _degenerate(label)
    return float(level @ response), float(level @ level), s2


def _raw_statistics(
    units: list[_Unit], adf: list[tuple[float, float, float]], longest: float
) -> dict[str, float]:
    n_units = len(units)
    sqrt_n = math.sqrt(n_units)
    l11 = np.array([u.l11 for u in units])
    sxx = np.array([u.sxx for u in units])
    num = np.array([u.num for u in units])
    sigma2 = np.array([u.sigma2 for u in units])
    a_num, a_sxx, a_s2 = (np.array(column) for column in zip(*adf, strict=True))
    t = longest - 1.0
    weighted_sxx = float(np.sum(sxx / l11))
    weighted_num = float(np.sum(num / l11))
    return {
        "panel_nu": t * t * n_units * sqrt_n / weighted_sxx,
        "panel_rho": t * sqrt_n * weighted_num / weighted_sxx,
        "panel_pp": weighted_num / math.sqrt(float(np.mean(sigma2 / l11)) * weighted_sxx),
        "panel_adf": float(np.sum(a_num / l11))
        / math.sqrt(float(np.mean(a_s2 / l11)) * float(np.sum(a_sxx / l11))),
        "group_rho": t / sqrt_n * float(np.sum(num / sxx)),
        "group_pp": float(np.sum(num / np.sqrt(sigma2 * sxx))) / sqrt_n,
        "group_adf": float(np.sum(a_num / np.sqrt(a_s2 * a_sxx))) / sqrt_n,
    }
