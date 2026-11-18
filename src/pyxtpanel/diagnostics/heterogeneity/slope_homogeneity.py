"""Pesaran-Yamagata test of slope homogeneity (category A, L.8.11).

Validated against Stata 17: ``xthst`` (Bersvendsen and Ditzen, Stata Journal 21(1)) with the
standard delta statistics, with variables partialled out, without a constant, with ``ar`` and
with cross-sectional averages and their lags, on balanced and unbalanced panels with homogeneous
and heterogeneous slopes, and with the heteroskedasticity and autocorrelation consistent version
of Blomquist and Westerlund (2013) (``hac``: Bartlett, quadratic spectral and truncated kernels,
automatic or fixed bandwidth, whitening).  The cases are in
``validation/cases/diagnostics_heterogeneity_slope_homogeneity.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.ols import PERFECT_FIT
from pyxtpanel.diagnostics.heterogeneity._hac import hac_swamy

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray, IntArray
    from pyxtpanel.panel.arrays import PanelArrays
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
        "diag.slope_homogeneity.too_short": {
            "en": "slope_homogeneity needs more periods: with {k} slopes and {partial} partialled "
            "variables the adjusted statistic needs an average of more than {need} observations "
            "per unit, the panel has {have:.1f}.",
            "fr": "slope_homogeneity exige plus de périodes : avec {k} pentes et {partial} "
            "variables partialisées, la statistique ajustée exige plus de {need} observations "
            "par unité en moyenne, le panel en a {have:.1f}.",
        },
        "diag.slope_homogeneity.zero_variance": {
            "en": "unit {unit!r} fits its regression perfectly: its error variance is zero and "
            "the statistic is not defined.",
            "fr": "l'unité {unit!r} s'ajuste parfaitement : sa variance d'erreur est nulle et la "
            "statistique n'est pas définie.",
        },
        "diag.slope_homogeneity.hac_ar": {
            "en": "hac and ar cannot be combined (xthst, r(184)).",
            "fr": "hac et ar ne peuvent pas être combinés (xthst, r(184)).",
        },
        "diag.slope_homogeneity.bad_option": {
            "en": "{name} must be {expected}, got {value!r}.",
            "fr": "{name} doit être {expected}, reçu {value!r}.",
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
            "sh_csa_lags",
            "sh_hac_bartlett",
            "sh_hac_qs",
            "sh_hac_truncated",
            "sh_hac_whitening",
            "sh_hac_bandwidth",
            "sh_hac_partial",
            "sh_hac_unbalanced",
            "sh_hac_unbalanced_qs",
            "sh_hac_unbalanced_whitening",
            "sh_hac_bandwidth_zero",
            "sh_hac_csa_lags",
            "sh_hac_one_regressor",
            "sh_hac_one_whitening",
            "sh_hac_one_truncated",
            "sh_hac_three_regressors",
            "sh_hac_three_truncated",
            "sh_hac_qs_bandwidth",
            "sh_hac_truncated_whitening",
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
    cr_lags: int | Sequence[int] = 0,
    adjusted: bool = True,
    hac: bool = False,
    kernel: Kernel = "bartlett",
    bandwidth: int | None = None,
    whitening: bool = False,
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

    With ``hac`` the weights ``s_i^2`` are replaced by a kernel estimate of the long-run covariance
    of the scores ``x_it e_it`` of each unit (Blomquist and Westerlund, 2013), which makes the
    test robust to heteroskedasticity and to autocorrelation of the errors; ``S`` is then the
    HAC statistic.

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
        count the partialled variables.  Cannot be combined with ``hac``.
    crosssectional
        Variables whose cross-sectional mean (over the units, period by period) is partialled
        out in each unit, which lets the test cope with cross-section dependence.
    cr_lags
        Lags of the cross-sectional averages that are partialled out too: one number for all the
        variables, or one per variable (the last one applies to the rest, as in ``xthst``).  The
        lag is read on the calendar and a row whose lagged average does not exist is dropped.
    adjusted
        Use the adjusted statistic as the result (``True``) or the standard one.
    hac
        Use the HAC version of the statistic.  Implied by ``whitening`` and by ``bandwidth``.
    kernel
        Kernel of the long-run covariance: ``"bartlett"`` (default), ``"qs"`` (quadratic
        spectral) or ``"truncated"``.
    bandwidth
        Bandwidth of every unit; ``None`` is ``xthst``'s automatic choice, which depends on the
        kernel and on the unit (the average is in ``options``).
    whitening
        Pre-whiten the scores with a first-order vector autoregression (Andrews and Monahan).
        With three regressors or more the result can differ from ``xthst``'s, which rebuilds the
        whitening matrix from singular vectors whose signs depend on the LAPACK in use.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If there is no regressor, a variable is both tested and partialled out, an option is
        invalid, or ``hac`` is combined with ``ar``.
    DataError
        If a unit has too few observations for its parameters, the average number of periods is
        too small for the adjusted statistic, or there are fewer than 2 units.
    NumericalError
        If a unit fits perfectly or the regressors are collinear after partialling out.
    """
    # Stata: xthst -> Mata deltatest(): partial Z out of y and X unit by unit, pooled FE slope,
    # sigma2_i = e_i'e_i / (T_i - Kpartial) (Kpartial = 0 with ar), weighted FE slope, S,
    # Delta = sqrt(N) (S/N - K) / sqrt(2K), adjusted: variance 2K (Tavg - K - Kpartial - 1) /
    # (Tavg - Kpartial + 1); two-sided normal p-values; with hac (also implied by whitening and
    # bw()) deltatesthac() in _hac.py gives S; cross-sectional averages come from xtdcce2_csa
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
    hac = hac or whitening or bandwidth is not None
    lags = _check_hac_options(ar, hac, kernel, bandwidth, cr_lags, len(averaged))
    wanted = list(dict.fromkeys([dependent, *names, *extra, *averaged]))
    arrays = panel_variables(panel, wanted, test="slope_homogeneity", min_units=2)
    design = _design(arrays, wanted, dependent, names, extra, averaged, lags, constant)
    k_regressors, k_partial = design.x.shape[1], design.z.shape[1]
    y_t, x_t = _partial_out(design)
    dof_partial = 0 if ar else k_partial
    n_units = design.offsets.size - 1
    t_average = design.y.size / n_units
    need = k_regressors + dof_partial + 1
    if t_average <= need:  # the variance of the adjusted statistic is not positive
        raise DataError(
            _i18n.message(
                "diag.slope_homogeneity.too_short",
                k=k_regressors,
                partial=k_partial,
                need=need,
                have=t_average,
            ),
            key="diag.slope_homogeneity.too_short",
        )
    average_bandwidth: float | None = None
    if hac:
        _check_slopes_identified(design, x_t)
        swamy, average_bandwidth = hac_swamy(
            y_t, x_t, design.offsets, kernel=kernel, bandwidth=bandwidth, whitening=whitening
        )
    else:
        swamy = _swamy(design, y_t, x_t, dof_partial)
    k = float(k_regressors)
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
    options: dict[str, object] = {
        "regressors": tuple(names),
        "partial": tuple(extra),
        "constant": constant,
        "ar": ar,
        "crosssectional": tuple(averaged),
        "adjusted": adjusted,
    }
    if averaged:
        options["cr_lags"] = tuple(lags)
    if hac:
        options.update(
            kernel=kernel,
            bandwidth=average_bandwidth,
            whitening=whitening,
            bandwidth_fixed=bandwidth is not None,
        )
    return TestResult(
        name="slope_homogeneity",
        statistic=delta_adjusted if adjusted else delta,
        pvalue=p_adjusted if adjusted else p_delta,
        statistic_name="Delta_adj" if adjusted else "Delta",
        n_obs=design.y.size,
        n_groups=n_units,
        variable=dependent,
        options=options,
        details=details,
    )


Kernel = Literal["bartlett", "qs", "truncated"]


@dataclass(frozen=True)
class _Design:
    """Data of the test after the cross-sectional lags have dropped their rows."""

    y: FloatArray
    x: FloatArray
    z: FloatArray
    offsets: IntArray
    labels: list[object]


def _check_hac_options(
    ar: bool,
    hac: bool,
    kernel: str,
    bandwidth: int | None,
    cr_lags: int | Sequence[int],
    n_averaged: int,
) -> list[int]:
    """Refuse an invalid option and return the lags of each cross-sectional average."""
    if hac and ar:
        raise SpecificationError(
            _i18n.message("diag.slope_homogeneity.hac_ar"), key="diag.slope_homogeneity.hac_ar"
        )
    if kernel not in ("bartlett", "qs", "truncated"):
        _bad_option("kernel", "'bartlett', 'qs' or 'truncated'", kernel)
    if bandwidth is not None and (
        isinstance(bandwidth, bool) or not isinstance(bandwidth, (int, np.integer)) or bandwidth < 0
    ):
        _bad_option("bandwidth", "None or an integer >= 0", bandwidth)
    given = [cr_lags] if isinstance(cr_lags, (int, np.integer)) else list(cr_lags)
    for value in given:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 0:
            _bad_option("cr_lags", "an integer >= 0 or a sequence of them", cr_lags)
    lags = [int(v) for v in given]
    # xthst: the lag of a variable is the last one given up to it
    return [lags[min(i, len(lags) - 1)] if lags else 0 for i in range(n_averaged)]


def _bad_option(name: str, expected: str, value: object) -> None:
    raise SpecificationError(
        _i18n.message(
            "diag.slope_homogeneity.bad_option", name=name, expected=expected, value=value
        ),
        key="diag.slope_homogeneity.bad_option",
    )


def _design(
    arrays: PanelArrays,
    wanted: list[str],
    dependent: str,
    names: list[str],
    extra: list[str],
    averaged: list[str],
    lags: list[int],
    constant: bool,
) -> _Design:
    """Assemble the variables, the partialled columns and the cross-sectional averages."""
    column = {name: arrays.values[:, i] for i, name in enumerate(wanted)}
    z_columns = [column[name] for name in extra]
    for name, lag in zip(averaged, lags, strict=True):
        mean = _period_mean(column[name], arrays.time_index)
        for k in range(lag + 1):
            if k == 0:
                z_columns.append(mean)
            else:  # a row whose lagged average does not exist is dropped (xthst: markout)
                source = arrays.lag_index(k)
                z_columns.append(np.where(source >= 0, mean[np.maximum(source, 0)], np.nan))
    if constant:
        z_columns.insert(0, np.ones(arrays.n_obs))
    y = column[dependent]
    x = np.column_stack([column[name] for name in names])
    z = np.column_stack(z_columns) if z_columns else np.empty((arrays.n_obs, 0))
    unit = arrays.unit_of_row
    keep = np.isfinite(z).all(axis=1)
    labels = list(arrays.unit_labels)
    if not keep.all():
        y, x, z, unit = y[keep], x[keep], z[keep], unit[keep]
        present, counts = np.unique(unit, return_counts=True)
        labels = [labels[int(u)] for u in present]
        return _Design(y, x, z, np.concatenate([[0], np.cumsum(counts)]).astype(np.int64), labels)
    return _Design(y, x, z, arrays.offsets, labels)


def _partial_out(design: _Design) -> tuple[FloatArray, FloatArray]:
    """Partial the constant, the variables and the averages out of y and X, unit by unit."""
    k_regressors, k_partial = design.x.shape[1], design.z.shape[1]
    y_t, x_t = design.y.copy(), design.x.copy()
    for i, label in enumerate(design.labels):
        rows = slice(int(design.offsets[i]), int(design.offsets[i + 1]))
        have = rows.stop - rows.start
        if have <= k_regressors + k_partial:
            raise DataError(
                _i18n.message(
                    "diag.slope_homogeneity.short_unit",
                    unit=_plain(label),
                    have=have,
                    need=k_regressors + k_partial,
                    more=k_regressors + k_partial + 1,
                ),
                key="diag.slope_homogeneity.short_unit",
            )
        if k_partial:
            y_t[rows], x_t[rows] = _residualise(design.y[rows], design.x[rows], design.z[rows])
    return y_t, x_t


def _check_slopes_identified(design: _Design, x_t: FloatArray) -> None:
    _check_full_rank(x_t.T @ x_t)
    for i in range(len(design.labels)):
        rows = slice(int(design.offsets[i]), int(design.offsets[i + 1]))
        _check_full_rank(x_t[rows].T @ x_t[rows])


def _swamy(design: _Design, y_t: FloatArray, x_t: FloatArray, dof_partial: int) -> float:
    """Swamy's statistic: unit slopes against the weighted fixed-effects slope."""
    _check_slopes_identified(design, x_t)
    resid = y_t - x_t @ np.linalg.solve(x_t.T @ x_t, x_t.T @ y_t)
    n_units = len(design.labels)
    sigma2 = np.empty(n_units)
    unit_xx: list[FloatArray] = []
    unit_b: list[FloatArray] = []
    for i, label in enumerate(design.labels):
        rows = slice(int(design.offsets[i]), int(design.offsets[i + 1]))
        e, xi = resid[rows], x_t[rows]
        sigma2[i] = float(e @ e) / (rows.stop - rows.start - dof_partial)
        # a perfect fit leaves rounding noise, not an exact zero: compare with the data's scale
        if sigma2[i] <= PERFECT_FIT * float(design.y[rows] @ design.y[rows]):
            raise NumericalError(
                _i18n.message("diag.slope_homogeneity.zero_variance", unit=_plain(label)),
                key="diag.slope_homogeneity.zero_variance",
            )
        xxi = xi.T @ xi
        unit_xx.append(xxi)
        unit_b.append(np.linalg.solve(xxi, xi.T @ y_t[rows]))
    weighted_xx = sum(xx / s2 for xx, s2 in zip(unit_xx, sigma2, strict=True))
    weighted_xy = sum((xx @ b) / s2 for xx, b, s2 in zip(unit_xx, unit_b, sigma2, strict=True))
    b_wfe = np.linalg.solve(weighted_xx, weighted_xy)
    return float(
        sum(
            (b - b_wfe) @ (xx / s2) @ (b - b_wfe)
            for xx, b, s2 in zip(unit_xx, unit_b, sigma2, strict=True)
        )
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
