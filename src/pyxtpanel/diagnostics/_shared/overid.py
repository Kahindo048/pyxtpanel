"""Quadratic form of the moments behind the over-identification tests (L.6.2, L.6.3).

Hansen's J and Sargan's statistic are both ``m' W m`` (divided by a scale for Sargan), where
``m = sum_i Z_i' e_i`` is the vector of moments summed over units and ``W`` the weight matrix of
the estimator.  What makes them differ is only the weight matrix and the scale, and those are
conventions of the estimator (which residuals, which instruments, which variance), so the
caller passes them.  This module holds the arithmetic, the degrees of freedom and the checks;
the difference-in-Hansen test of L.6.3 reuses it.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, as_array

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray

_i18n.register(
    "diag",
    {
        "diag.overid.weight_not_symmetric": {
            "en": "weight must be symmetric (largest asymmetry {gap:.3g}); pass the weight matrix "
            "of the estimator, not a transposed or partial one.",
            "fr": "weight doit être symétrique (plus grand écart {gap:.3g}) ; passez la matrice "
            "de poids de l'estimateur, non une version transposée ou partielle.",
        },
        "diag.overid.negative": {
            "en": "the statistic is negative ({value:.6g}): the weight matrix is not positive "
            "semi-definite.",
            "fr": "la statistique est négative ({value:.6g}) : la matrice de poids n'est pas "
            "semi-définie positive.",
        },
        "diag.overid.bad_count": {
            "en": "{name} must be {expected}, got {value}.",
            "fr": "{name} doit être {expected}, reçu {value}.",
        },
        "diag.overid.not_overidentified": {
            "en": "{n_instruments} instruments for {n_parameters} parameters: the model is not "
            "over-identified, so there is nothing to test.",
            "fr": "{n_instruments} instruments pour {n_parameters} paramètres : le modèle n'est "
            "pas sur-identifié, il n'y a rien à tester.",
        },
        "diag.overid.proliferation": {
            "en": "{n_instruments} instruments exceed the {n_groups} groups: the test is "
            "weakened and its p-value cannot be trusted (instrument proliferation).",
            "fr": "{n_instruments} instruments dépassent les {n_groups} groupes : le test est "
            "affaibli et sa p-value n'est pas fiable (prolifération d'instruments).",
        },
    },
)


def moment_covariance(unit_moments: object, *, center: bool = False) -> FloatArray:
    """Return ``sum_i g_i g_i'`` for the moments ``g_i = Z_i' e_i`` of each unit.

    Its inverse (or generalised inverse) is the optimal weight matrix of two-step GMM.

    Parameters
    ----------
    unit_moments
        ``(n_units, q)`` array, one row of moments per unit.
    center
        Subtract the mean of the moments first (Stata's ``center`` option); default ``False``.
    """
    g = as_array("unit_moments", unit_moments, ndim=2)
    if center:
        g = g - g.mean(axis=0)
    out: FloatArray = g.T @ g
    return out


def overidentification(
    name: str,
    moments: object,
    weight: object,
    n_parameters: int,
    *,
    scale: float = 1.0,
    n_instruments: int | None = None,
    n_groups: int | None = None,
) -> TestResult:
    """Return the test ``m' W m / scale`` against a chi-square with the excess instruments.

    Parameters
    ----------
    name
        Registered name of the test, key of its texts.
    moments
        ``m = sum_i Z_i' e_i`` as a ``(q,)`` vector, or the ``(n_units, q)`` moments of each
        unit (summed here, and ``n_groups`` defaults to the number of rows).
    weight
        ``(q, q)`` weight matrix of the estimator; symmetrised.
    n_parameters
        Number of estimated parameters; the degrees of freedom are
        ``n_instruments - n_parameters``.
    scale
        Positive divisor of the quadratic form (the residual variance for Sargan).
    n_instruments
        Number of independent instruments, ``q`` by default (pass the rank when instruments are
        collinear).
    n_groups
        Number of units; adds a remark when instruments outnumber them.

    Raises
    ------
    SpecificationError
        If a shape, a count or the scale is invalid, or the weight matrix is not symmetric.
    NumericalError
        If the weight matrix is not positive semi-definite (negative statistic).
    """
    raw = np.asarray(moments, dtype=np.float64)
    if raw.ndim == 2:
        g = as_array("moments", raw, ndim=2)
        m = g.sum(axis=0)
        n_groups = g.shape[0] if n_groups is None else n_groups
    else:
        m = as_array("moments", raw, ndim=1)
    q = m.size
    w = as_array("weight", weight, ndim=2, shape=(q, q))
    gap = float(np.max(np.abs(w - w.T))) if q else 0.0
    if gap > 1e-8 * max(float(np.max(np.abs(w))) if q else 0.0, np.finfo(float).tiny):
        raise SpecificationError(
            _i18n.message("diag.overid.weight_not_symmetric", gap=gap),
            key="diag.overid.weight_not_symmetric",
        )
    w = (w + w.T) / 2.0
    instruments = q if n_instruments is None else n_instruments
    _positive_int("n_parameters", n_parameters, minimum=0)
    _positive_int("n_instruments", instruments, minimum=1, maximum=q)
    if not (math.isfinite(scale) and scale > 0):
        raise SpecificationError(
            _i18n.message("diag.overid.bad_count", name="scale", expected="positive", value=scale),
            key="diag.overid.bad_count",
        )
    df = instruments - n_parameters
    value = float(m @ w @ m) / scale
    if value < -1e-9 * max(abs(value), 1.0) - 1e-12:
        raise NumericalError(
            _i18n.message("diag.overid.negative", value=value), key="diag.overid.negative"
        )
    value = max(value, 0.0)
    notes: list[Note] = []
    if df <= 0:
        notes.append(
            Note(
                "diag.overid.not_overidentified",
                {"n_instruments": instruments, "n_parameters": n_parameters},
            )
        )
        statistic, pvalue = math.nan, math.nan
    else:
        statistic, pvalue = value, float(stats.chi2.sf(value, df))
    if n_groups is not None and instruments > n_groups:
        notes.append(
            Note("diag.overid.proliferation", {"n_instruments": instruments, "n_groups": n_groups})
        )
    options: dict[str, object] = {"n_instruments": instruments, "n_parameters": n_parameters}
    if scale != 1.0:
        options["scale"] = scale
    return TestResult(
        name=name,
        statistic=statistic,
        pvalue=pvalue,
        statistic_name="chi2",
        df=(float(df),) if df > 0 else (),
        n_groups=n_groups,
        options=options,
        notes=tuple(notes),
    )


def _positive_int(name: str, value: int, *, minimum: int, maximum: int | None = None) -> None:
    ok = float(value) == int(value) and value >= minimum and (maximum is None or value <= maximum)
    if not ok:
        expected = f"an integer >= {minimum}" + (f" and <= {maximum}" if maximum else "")
        raise SpecificationError(
            _i18n.message("diag.overid.bad_count", name=name, expected=expected, value=value),
            key="diag.overid.bad_count",
        )
