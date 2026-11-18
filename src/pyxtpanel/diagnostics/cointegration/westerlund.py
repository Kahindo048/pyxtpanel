"""Westerlund variance-ratio test for cointegration (category A, L.4.3).

Validated against Stata 17: ``xtcointtest westerlund`` with panel-specific (``somepanels``) and
common (``allpanels``) autoregressive parameters, with a constant and with a trend, with one to
seven regressors, with and without the cross-sectional means removed, on balanced and unbalanced
panels.  The cases are in ``validation/cases/diagnostics_cointegration_westerlund.py``.  The
moments that standardise the statistic are not published by Stata (they are compiled in its Mata
library): they were deduced from its output on simulated panels, see
``pyxtpanel/_data/tables/MANIFEST.json``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.tables import load_table

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.westerlund.title": {
            "en": "Westerlund variance-ratio test for cointegration",
            "fr": "Test de cointégration de Westerlund (rapport de variances)",
        },
        "diag.westerlund.null": {
            "en": "no cointegration",
            "fr": "pas de cointégration",
        },
        "diag.westerlund.alternative": {
            "en": "some or all panels are cointegrated (see the option alternative)",
            "fr": "certaines ou toutes les unités sont cointégrées (voir l'option alternative)",
        },
        "diag.westerlund.bad_option": {
            "en": "{name} must be one of {allowed}, got {value!r}.",
            "fr": "{name} doit valoir l'une de ces valeurs : {allowed} ; reçu {value!r}.",
        },
        "diag.westerlund.regressors": {
            "en": "give between 1 and 7 distinct regressors different from the dependent "
            "variable, got {value!r}.",
            "fr": "donnez de 1 à 7 régresseurs distincts et différents de la variable "
            "dépendante, reçu {value!r}.",
        },
        "diag.westerlund.degenerate": {
            "en": "the residuals of unit {unit!r} are zero: the regression fits perfectly.",
            "fr": "les résidus de l'unité {unit!r} sont nuls : la régression est parfaite.",
        },
    },
)


@register_test(
    "westerlund",
    TestSpec(
        family="cointegration",
        category="A",
        wave=2,
        reference="Westerlund, J. (2005), Econometric Reviews 24(3), 297-316",
        stata="xtcointtest westerlund",
        validation=(
            "we_some_m1",
            "we_all_m1",
            "we_some_m3_trend",
            "we_all_m3_trend",
            "we_some_m7",
            "we_all_m7_trend",
            "we_demean",
            "we_coint",
            "we_unbalanced",
            "we_unbalanced_all_trend",
            "we_unbalanced_demean_all",
        ),
        inputs="panel",
        balanced="both",
        tables=("westerlund_moments.csv",),
        summary="Variance ratio of the cumulated residuals of the cointegrating regressions.",
    ),
)
def westerlund(
    panel: PanelBase,
    dependent: str,
    regressors: Sequence[str],
    *,
    alternative: Literal["some", "all"] = "some",
    deterministics: Literal["constant", "trend"] = "constant",
    demean: bool = False,
) -> TestResult:
    """Return Westerlund's variance-ratio statistic for no cointegration.

    For each unit the dependent variable is regressed on a constant (and a trend with
    ``deterministics="trend"``) and the regressors, with its own coefficients.  With the residuals
    ``e`` of unit ``i``, ``E_t`` their cumulated sum and ``R_i`` the sum of their squares, the
    unit's variance ratio is ``sum_t E_t^2 / (T^2 R_i)``.  ``alternative="some"`` (the default)
    uses its average over the units, with a null of no cointegration against cointegration in
    some units; ``alternative="all"`` divides the sum of the ``E_t^2`` over all units by the sum of
    the ``R_i`` (a common autoregressive parameter), against cointegration in all units.  The
    statistic is ``sqrt(N) (VR - mean) / sqrt(variance)`` with the mean and variance of Westerlund
    (2005), by number of regressors and deterministic terms, and is standard normal under the
    null; it rejects for small values.  ``T`` is the longest unit (the number of periods of a
    balanced panel).

    The panel may be unbalanced but must have no gaps.  At most seven regressors.

    Parameters
    ----------
    panel
        The declared panel.
    dependent
        Name of the dependent variable.
    regressors
        Names of the regressors (1 to 7).
    alternative
        ``"some"`` or ``"all"`` (see above).
    deterministics
        ``"constant"`` (default) or ``"trend"``.
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
        If a unit's regression fits perfectly.
    """
    # Stata: xtcointtest westerlund y x [, trend allpanels demean] -> per unit regress y on 1 [t]
    # x; E = cumulated residuals, R_i = sum e^2; somepanels: VR = mean_i sum E^2 / R_i; allpanels:
    # VR = sum_i sum E^2 / sum_i R_i; z = sqrt(N) (VR / T^2 - mean) / sqrt(variance), mean and
    # variance deduced from Stata's output for 1 to 7 regressors; p-value 1 - normal(|z|).
    # STATA-QUIRK: T is the largest unit length on an unbalanced panel (not the mean), and the
    # reported p-value is 1 - normal(|z|); the left tail normal(z) is used here.
    names = _names(dependent, regressors)
    for name, value, allowed in (
        ("alternative", alternative, ("some", "all")),
        ("deterministics", deterministics, ("constant", "trend")),
    ):
        if value not in allowed:
            raise SpecificationError(
                _i18n.message(
                    "diag.westerlund.bad_option",
                    name=name,
                    allowed=", ".join(map(repr, allowed)),
                    value=value,
                ),
                key="diag.westerlund.bad_option",
            )
    trend = deterministics == "trend"
    arrays = panel_variables(
        panel,
        [dependent, *names],
        test="westerlund",
        consecutive=True,
        min_obs=len(names) + int(trend) + 3,
        min_units=2,
    )
    values = arrays.values.copy()
    if demean:
        counts = np.bincount(arrays.time_index)
        for column in range(values.shape[1]):
            totals = np.bincount(arrays.time_index, weights=values[:, column])
            values[:, column] -= (totals / np.where(counts > 0, counts, 1))[arrays.time_index]
    numerators: list[float] = []
    squares: list[float] = []
    for unit, window in enumerate(arrays.unit_slices()):
        e = _residuals(values[window], trend)
        total = float(e @ e)
        if total <= 1e-24 * max(float(np.sum(values[window][:, 0] ** 2)), 1.0):
            raise NumericalError(
                _i18n.message("diag.westerlund.degenerate", unit=_plain(arrays.unit_labels[unit])),
                key="diag.westerlund.degenerate",
            )
        numerators.append(float(np.sum(np.cumsum(e) ** 2)))
        squares.append(total)
    num = np.array(numerators, dtype=np.float64)
    sq = np.array(squares, dtype=np.float64)
    longest = float(np.max(arrays.counts))
    n_units = arrays.n_units
    ratio = float(np.mean(num / sq)) if alternative == "some" else float(num.sum() / sq.sum())
    variance_ratio = ratio / (longest * longest)
    mean, variance = _moments(len(names), deterministics, alternative)
    statistic = math.sqrt(n_units) * (variance_ratio - mean) / math.sqrt(variance)
    details = pd.DataFrame(
        {"variance_ratio": num / sq / (longest * longest), "periods": arrays.counts.astype(float)},
        index=pd.Index(arrays.unit_labels, name="unit"),
    )
    return TestResult(
        name="westerlund",
        statistic=float(statistic),
        pvalue=float(stats.norm.cdf(statistic)),
        statistic_name="VR",
        df=(),
        n_obs=arrays.n_obs,
        n_groups=n_units,
        variable=dependent,
        options={
            "regressors": list(names),
            "alternative": alternative,
            "deterministics": deterministics,
            "demean": demean,
            "variance_ratio": variance_ratio,
        },
        details=details,
    )


def _plain(label: object) -> object:
    return label.item() if isinstance(label, np.generic) else label


def _names(dependent: str, regressors: Sequence[str]) -> list[str]:
    names = [regressors] if isinstance(regressors, str) else list(regressors)
    if not 1 <= len(names) <= 7 or len(set(names)) != len(names) or dependent in names:
        raise SpecificationError(
            _i18n.message("diag.westerlund.regressors", value=regressors),
            key="diag.westerlund.regressors",
        )
    return names


def _residuals(block: FloatArray, trend: bool) -> FloatArray:
    """Residuals of the first column on a constant, a trend if asked, and the other columns."""
    t = block.shape[0]
    columns = [np.ones(t)]
    if trend:
        columns.append(np.arange(1.0, t + 1.0))
    design = np.column_stack([*columns, block[:, 1:]])
    coefficients, *_ = np.linalg.lstsq(design, block[:, 0], rcond=None)
    residual: FloatArray = block[:, 0] - design @ coefficients
    return residual


def _moments(m: int, deterministics: str, alternative: str) -> tuple[float, float]:
    table = load_table("westerlund_moments.csv")
    row = table[
        (table["m"] == m)
        & (table["deterministics"] == deterministics)
        & (table["alternative"] == alternative)
    ].iloc[0]
    return float(row["mean"]), float(row["variance"])
