"""Pesaran-Yamagata test of slope homogeneity (category A, L.8.11).

Validated against Stata 17: ``xthst`` (Bersvendsen and Ditzen, Stata Journal 21(1)) with the
standard delta statistics, with variables partialled out, without a constant, with ``ar`` and
with cross-sectional averages, on balanced and unbalanced panels with homogeneous and
heterogeneous slopes.  The cases are in
``validation/cases/diagnostics_heterogeneity_slope_homogeneity.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, panel_variables, register_test

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray, IntArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.slope_homogeneity.title": {
            "en": "Pesaran-Yamagata test of slope homogeneity",
            "fr": "Test d'homogénéité des pentes de Pesaran-Yamagata",
        },
        "diag.slope_homogeneity.null": {
            "en": "the slope coefficients are the same in every unit",
            "fr": "les coefficients de pente sont identiques dans toutes les unités",
        },
        "diag.slope_homogeneity.alternative": {
            "en": "the slope coefficients differ across units",
            "fr": "les coefficients de pente diffèrent selon les unités",
        },
        "diag.slope_homogeneity.no_regressor": {
            "en": "slope_homogeneity needs at least one regressor besides the dependent variable.",
            "fr": "slope_homogeneity exige au moins un régresseur en plus de la variable "
            "dépendante.",
        },
        "diag.slope_homogeneity.overlap": {
            "en": "{names} cannot be both a tested variable and a partialled one.",
            "fr": "{names} ne peut pas être à la fois une variable testée et une variable "
            "partialisée.",
        },
        "diag.slope_homogeneity.short_unit": {
            "en": "slope_homogeneity needs more observations than parameters in every unit: "
            "unit {unit!r} has {have}, the {need} parameters (slopes and partialled variables) "
            "need at least {more}.",
            "fr": "slope_homogeneity exige plus d'observations que de paramètres dans chaque "
            "unité : l'unité {unit!r} en a {have}, les {need} paramètres (pentes et variables "
            "partialisées) en exigent au moins {more}.",
        },
        "diag.slope_homogeneity.zero_variance": {
            "en": "unit {unit!r} fits its regression perfectly: its error variance is zero and "
            "the statistic is not defined.",
            "fr": "l'unité {unit!r} s'ajuste parfaitement : sa variance d'erreur est nulle et la "
            "statistique n'est pas définie.",
        },
        "diag.slope_homogeneity.singular": {
            "en": "the regressors are collinear after partialling out: the slopes are not "
            "identified.",
            "fr": "les régresseurs sont colinéaires une fois les variables partialisées : les "
            "pentes ne sont pas identifiées.",
        },
    },
)


@register_test(
    "slope_homogeneity",
    TestSpec(
        family="heterogeneity",
        category="A",
        wave=2,
        reference="Pesaran, M. H. and Yamagata, T. (2008), Journal of Econometrics 142(1), 50-93; "
        "Swamy, P. A. V. B. (1970), Econometrica 38(2), 311-323",
        stata="xthst",
        validation=(
            "sh_homogeneous",
            "sh_heterogeneous",
            "sh_unbalanced",
            "sh_partial",
            "sh_noconstant",
            "sh_ar",
            "sh_csa",
        ),
        inputs="panel",
        balanced="both",
        aliases=("xthst", "pesaran_yamagata", "swamy"),
        summary="Test that the slope coefficients are the same in every unit.",
    ),
)
def slope_homogeneity(
    panel: PanelBase,
    dependent: str,
    regressors: Sequence[str],
    *,
    partial: Sequence[str] = (),
    constant: bool = True,
    ar: bool = False,
    crosssectional: Sequence[str] = (),
    adjusted: bool = True,
) -> TestResult:
    """Return the Pesaran-Yamagata ``Delta`` test that the slopes of ``regressors`` are equal.

    Every unit is regressed on its own (OLS) after the constant, the ``partial`` variables and the
    cross-sectional averages have been partialled out of it, which lets those vary freely across
    units.  The slopes of that regression are compared with the weighted fixed-effects slope with
    Swamy's statistic ``S = sum_i (b_i - b_wfe)' (X_i'X_i / s_i^2) (b_i - b_wfe)``, where
    ``s_i^2`` is the error variance of unit ``i`` from the pooled fit.  ``Delta`` standardises it,
    ``sqrt(N) (S / N - K) / sqrt(2 K)``, and is standard normal for large ``N`` and ``T``; the
    adjusted form, which has the exact mean and variance under normal errors and is better in
    small samples, is the default statistic.  ``details`` gives the two forms with their
    two-sided p-values (the lines ``Delta`` and ``adj.`` of ``xthst``) and Swamy's ``S`` with its
    chi-square p-value for fixed ``N``.

    Parameters
    ----------
    panel
        The declared panel (balanced or not).
    dependent
        Name of the dependent variable.
    regressors
        Names of the regressors whose slopes are tested.
    partial
        Variables partialled out in each unit before testing (their slopes may differ).
    constant
        Partial out a constant (the unit effect); ``False`` is ``xthst``'s ``noconstant``.
    ar
        Autoregressive model (``xthst``'s ``ar``): the degrees-of-freedom corrections do not
        count the partialled variables.
    crosssectional
        Variables whose cross-sectional mean (over the units, period by period) is partialled
        out in each unit, which lets the test cope with cross-section dependence.
    adjusted
        Use the adjusted statistic as the result (``True``) or the standard one.

    Raises
    ------
    SpecificationError
        If there is no regressor or a variable is both tested and partialled out.
    DataError
        If a unit has too few observations for its parameters.
    NumericalError
        If a unit fits perfectly or the regressors are collinear after partialling out.
    """
    # Stata: xthst -> Mata deltatest(): partial Z out of y and X unit by unit, pooled FE slope,
    # sigma2_i = e_i'e_i / (T_i - Kpartial) (Kpartial = 0 with ar), weighted FE slope, S,
    # Delta = sqrt(N) (S/N - K) / sqrt(2K), adjusted: variance 2K (Tavg - K - Kpartial - 1) /
    # (Tavg - Kpartial + 1); two-sided normal p-values
    names, extra, averaged = list(regressors), list(partial), list(crosssectional)
    if not names:
        raise SpecificationError(
            _i18n.message("diag.slope_homogeneity.no_regressor"),
            key="diag.slope_homogeneity.no_regressor",
        )
    clash = sorted({dependent, *names} & set(extra))
    if clash:
        raise SpecificationError(
            _i18n.message("diag.slope_homogeneity.overlap", names=", ".join(clash)),
            key="diag.slope_homogeneity.overlap",
        )
    wanted = list(dict.fromkeys([dependent, *names, *extra, *averaged]))
    arrays = panel_variables(panel, wanted, test="slope_homogeneity")
    column = {name: arrays.values[:, i] for i, name in enumerate(wanted)}
    y = column[dependent]
    x = np.column_stack([column[name] for name in names])
    z_columns = [column[name] for name in extra]
    z_columns += [_period_mean(column[name], arrays.time_index) for name in averaged]
    if constant:
        z_columns.insert(0, np.ones(arrays.n_obs))
    z = np.column_stack(z_columns) if z_columns else np.empty((arrays.n_obs, 0))
    k_regressors, k_partial = x.shape[1], z.shape[1]
    # Partial Z out of y and X, unit by unit.
    y_t, x_t = y.copy(), x.copy()
    labels = arrays.unit_labels
    for i in range(arrays.n_units):
        rows = slice(int(arrays.offsets[i]), int(arrays.offsets[i + 1]))
        have = rows.stop - rows.start
        if have <= k_regressors + k_partial:
            raise DataError(
                _i18n.message(
                    "diag.slope_homogeneity.short_unit",
                    unit=_plain(labels[i]),
                    have=have,
                    need=k_regressors + k_partial,
                    more=k_regressors + k_partial + 1,
                ),
                key="diag.slope_homogeneity.short_unit",
            )
        if k_partial:
            y_t[rows], x_t[rows] = _residualise(y[rows], x[rows], z[rows])
    # Pooled fixed-effects slope on the partialled data and the residual variance of each unit.
    xx_total = x_t.T @ x_t
    _check_full_rank(xx_total)
    pooled = np.linalg.solve(xx_total, x_t.T @ y_t)
    resid = y_t - x_t @ pooled
    dof_partial = 0 if ar else k_partial
    n_units = arrays.n_units
    sigma2 = np.empty(n_units)
    unit_xx: list[FloatArray] = []
    unit_b: list[FloatArray] = []
    for i in range(n_units):
        rows = slice(int(arrays.offsets[i]), int(arrays.offsets[i + 1]))
        e, xi = resid[rows], x_t[rows]
        sigma2[i] = float(e @ e) / (rows.stop - rows.start - dof_partial)
        if sigma2[i] <= 0.0:
            raise NumericalError(
                _i18n.message("diag.slope_homogeneity.zero_variance", unit=_plain(labels[i])),
                key="diag.slope_homogeneity.zero_variance",
            )
        xxi = xi.T @ xi
        _check_full_rank(xxi)
        unit_xx.append(xxi)
        unit_b.append(np.linalg.solve(xxi, xi.T @ y_t[rows]))
    weighted_xx = sum(xx / s2 for xx, s2 in zip(unit_xx, sigma2, strict=True))
    weighted_xy = sum((xx @ b) / s2 for xx, b, s2 in zip(unit_xx, unit_b, sigma2, strict=True))
    b_wfe = np.linalg.solve(weighted_xx, weighted_xy)
    swamy = float(
        sum(
            (b - b_wfe) @ (xx / s2) @ (b - b_wfe)
            for xx, b, s2 in zip(unit_xx, unit_b, sigma2, strict=True)
        )
    )
    k = float(k_regressors)
    t_average = arrays.n_obs / n_units
    delta = math.sqrt(n_units) * (swamy / n_units - k) / math.sqrt(2.0 * k)
    variance = 2.0 * k * (t_average - k - dof_partial - 1.0) / (t_average - dof_partial + 1.0)
    delta_adjusted = math.sqrt(n_units) * (swamy / n_units - k) / math.sqrt(variance)
    p_delta = 2.0 * float(stats.norm.sf(abs(delta)))
    p_adjusted = 2.0 * float(stats.norm.sf(abs(delta_adjusted)))
    df_swamy = k * (n_units - 1.0)
    details = pd.DataFrame(
        {
            "statistic": [delta, delta_adjusted, swamy],
            "pvalue": [p_delta, p_adjusted, float(stats.chi2.sf(swamy, df_swamy))],
        },
        index=pd.Index(["delta", "delta_adjusted", "swamy"], name="statistic_name"),
    )
    return TestResult(
        name="slope_homogeneity",
        statistic=delta_adjusted if adjusted else delta,
        pvalue=p_adjusted if adjusted else p_delta,
        statistic_name="Delta_adj" if adjusted else "Delta",
        n_obs=arrays.n_obs,
        n_groups=n_units,
        variable=dependent,
        options={
            "regressors": tuple(names),
            "partial": tuple(extra),
            "constant": constant,
            "ar": ar,
            "crosssectional": tuple(averaged),
            "adjusted": adjusted,
        },
        details=details,
    )


def _period_mean(values: FloatArray, time_index: IntArray) -> FloatArray:
    """Return, for every row, the mean of ``values`` over the units observed in its period."""
    _, code = np.unique(time_index, return_inverse=True)
    totals = np.bincount(code, weights=values)
    counts = np.bincount(code)
    mean: FloatArray = (totals / counts)[code]
    return mean


def _residualise(y: FloatArray, x: FloatArray, z: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Return ``y`` and ``x`` with the projection on the columns of ``z`` removed."""
    target = np.column_stack([y, x])
    coefficients, *_ = np.linalg.lstsq(z, target, rcond=1e-12)
    cleaned = target - z @ coefficients
    return cleaned[:, 0], cleaned[:, 1:]


def _check_full_rank(matrix: FloatArray) -> None:
    if np.linalg.matrix_rank(matrix) < matrix.shape[0]:
        raise NumericalError(
            _i18n.message("diag.slope_homogeneity.singular"),
            key="diag.slope_homogeneity.singular",
        )


def _plain(label: object) -> object:
    return label.item() if isinstance(label, np.generic) else label
