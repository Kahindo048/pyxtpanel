"""Within, between, quasi-demeaning and two-way elimination (D.2.1, D.2.3, D.2.4).

Every function takes the sorted index and a ``float64`` array in sorted order (``(n_obs,)`` or
``(n_obs, k)``) and ignores missing values: a variable is transformed on its own observed rows
(the ``listwise`` option of the frame level reproduces Stata's ``xtdata``, which drops a row as
soon as one variable is missing).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel._types import FloatArray, IntArray
from pyxtpanel.ops._spec import OpSpec, operator
from pyxtpanel.ops.groups import as_columns, group_codes, group_mean
from pyxtpanel.panel.index import PanelIndex

_i18n.register(
    "ops",
    {
        "ops.no_convergence": {
            "en": "two-way elimination did not converge in {maxiter} iterations (last change "
            "{change:.3g}, tolerance {tol:.3g}); raise maxiter or relax tol.",
            "fr": "l'élimination à deux effets n'a pas convergé en {maxiter} itérations "
            "(dernier écart {change:.3g}, tolérance {tol:.3g}) ; augmentez maxiter ou "
            "relâchez tol.",
        },
        "ops.bad_theta": {
            "en": "theta must be a scalar or one value per unit ({n_units} units), got {got}.",
            "fr": "theta doit être un scalaire ou une valeur par unité ({n_units} unités), "
            "reçu {got}.",
        },
    },
)


@operator(
    "between",
    OpSpec(
        kind="transform",
        missing="missing where the row's own value is missing",
        stata="xtdata, be (rows kept: those with a value)",
        summary="Mean of the unit over its observed values, repeated on each row.",
    ),
)
def between(index: PanelIndex, values: FloatArray) -> FloatArray:
    """Mean of each unit over its observed values, repeated on the rows where ``x`` is observed."""
    return group_mean(index, values, "unit")


@operator(
    "within",
    OpSpec(
        kind="transform",
        missing="missing where the row's own value is missing",
        stata="xtdata, fe  (x - unit mean + overall mean with add_mean=True)",
        summary="Deviation from the mean of the unit.",
    ),
)
def within(index: PanelIndex, values: FloatArray, *, add_mean: bool = False) -> FloatArray:
    """``x_it − x̄_i``; with ``add_mean=True`` the overall mean is added back (``xtdata, fe``).

    The overall mean is the mean of the observed rows of the variable, not the mean of the unit
    means (the two differ in an unbalanced panel); it is the one ``xtdata, fe`` uses.
    """
    data, flat = as_columns(values)
    # Stata: xtdata, fe  (x - unit mean + overall mean)
    out = data - group_mean(index, data, "unit")
    if add_mean:
        out = out + np.nanmean(data, axis=0, keepdims=True)
    result: FloatArray = out[:, 0] if flat else out
    return result


def theta_re(sigma_u: float, sigma_e: float, counts: FloatArray | IntArray) -> FloatArray:
    """Random-effects ``θ_i = 1 − σ_e / √(σ_e² + T_i σ_u²)`` for each unit.

    Parameters
    ----------
    sigma_u, sigma_e
        Standard deviations of the unit effect and of the idiosyncratic error.
    counts
        ``T_i`` of each unit (``panel.counts``).
    """
    if sigma_e < 0 or sigma_u < 0:
        raise SpecificationError(
            _i18n.message("ops.bad_theta", n_units=len(counts), got="negative sigma"),
            key="ops.bad_theta",
        )
    t = np.asarray(counts, dtype=np.float64)
    result: FloatArray = 1.0 - sigma_e / np.sqrt(sigma_e**2 + t * sigma_u**2)
    return result


@operator(
    "quasi_demean",
    OpSpec(
        kind="transform",
        missing="missing where the row's own value is missing",
        stata="xtreg, re  (x - theta * unit mean)",
        summary="Quasi-demeaning x − θ x̄_i with θ given by the user.",
    ),
)
def quasi_demean(
    index: PanelIndex, values: FloatArray, theta: float | FloatArray = 1.0
) -> FloatArray:
    """``x_it − θ_i x̄_i`` with a scalar ``θ`` or one value per unit (``θ = 1`` is the within)."""
    scale = np.asarray(theta, dtype=np.float64)
    if scale.ndim == 0:
        per_row: float | FloatArray = float(scale)
    elif scale.shape == (index.n_units,):
        per_row = scale[index.unit]
    else:
        raise SpecificationError(
            _i18n.message("ops.bad_theta", n_units=index.n_units, got=str(scale.shape)),
            key="ops.bad_theta",
        )
    data, flat = as_columns(values)
    weight = per_row if np.isscalar(per_row) else np.asarray(per_row).reshape(-1, 1)
    out = data - weight * group_mean(index, data, "unit")
    result: FloatArray = out[:, 0] if flat else out
    return result


@dataclass(frozen=True)
class Convergence:
    """How a two-way elimination ended (one entry per column).

    Attributes
    ----------
    exact
        True where the closed form was used (complete, strongly balanced column): no iteration.
    iterations
        Rounds of alternating projections used (0 for the closed form).
    max_change
        Largest correction applied by the last round.
    converged
        Whether ``max_change`` fell under the tolerance (always true when ``exact``).
    """

    exact: tuple[bool, ...]
    iterations: tuple[int, ...]
    max_change: tuple[float, ...]
    converged: tuple[bool, ...]

    def as_dict(self) -> dict[str, tuple[object, ...]]:
        """Return the report as a plain dictionary (for ``DataFrame.attrs``)."""
        return {
            "exact": self.exact,
            "iterations": self.iterations,
            "max_change": self.max_change,
            "converged": self.converged,
        }


def _closed_form(
    unit: IntArray, time: IntArray, n_units: int, n_times: int, x: FloatArray
) -> FloatArray:
    """``x − x̄_i − x̄_t + x̄`` (exact for a complete, strongly balanced column)."""
    mean_i = np.bincount(unit, weights=x, minlength=n_units) / n_times
    mean_t = np.bincount(time, weights=x, minlength=n_times) / n_units
    result: FloatArray = x - mean_i[unit] - mean_t[time] + x.mean()
    return result


def twoway_demean(
    index: PanelIndex,
    values: FloatArray,
    *,
    tol: float = 1e-12,
    maxiter: int = 10_000,
    strict: bool = True,
) -> tuple[FloatArray, Convergence]:
    """Remove unit and period effects from each column.

    A complete column of a strongly balanced panel uses the closed form
    ``x − x̄_i − x̄_t + x̄``.  Otherwise the unit and period means of the observed rows are
    removed in turn until the largest correction is below ``tol`` (alternating projections),
    which converges to the residual of the regression on unit and period dummies.

    Parameters
    ----------
    index
        Sorted index of the panel.
    values
        Columns to transform, ``NaN`` for missing; only observed rows take part.
    tol
        Stopping tolerance on the largest correction, relative to the largest absolute value.
    maxiter
        Maximum number of rounds per column.
    strict
        Raise :class:`~pyxtpanel.NumericalError` when a column does not converge; otherwise
        return the last iterate and say so in the report.

    Returns
    -------
    (numpy.ndarray, Convergence)
        The transformed array (same shape as ``values``) and the convergence report.
    """
    data, flat = as_columns(values)
    time = group_codes(index, "time")[0]
    out = np.full_like(data, np.nan)
    exact: list[bool] = []
    iterations: list[int] = []
    change: list[float] = []
    converged: list[bool] = []
    for j in range(data.shape[1]):
        usable = ~np.isnan(data[:, j])
        x = data[usable, j]
        if x.size == 0:
            _log(exact, iterations, change, converged, True, 0, 0.0, True)
            continue
        u, n_units = _compact(index.unit[usable])
        t, n_times = _compact(time[usable])
        if x.size == n_units * n_times and np.unique(u * n_times + t).size == x.size:
            out[usable, j] = _closed_form(u, t, n_units, n_times, x)
            _log(exact, iterations, change, converged, True, 0, 0.0, True)
            continue
        residual, rounds, last = _alternate(x, u, t, n_units, n_times, tol, maxiter)
        out[usable, j] = residual
        ok = last <= tol * max(1.0, float(np.max(np.abs(x))))
        if not ok and strict:
            raise NumericalError(
                _i18n.message("ops.no_convergence", maxiter=maxiter, change=last, tol=tol),
                key="ops.no_convergence",
            )
        _log(exact, iterations, change, converged, False, rounds, last, ok)
    report = Convergence(tuple(exact), tuple(iterations), tuple(change), tuple(converged))
    result: FloatArray = out[:, 0] if flat else out
    return result, report


def _compact(codes: IntArray) -> tuple[IntArray, int]:
    """Renumber the codes present as ``0..k-1``."""
    present, compact = np.unique(codes, return_inverse=True)
    return compact.astype(np.int64), int(present.size)


def _log(
    exact: list[bool],
    iterations: list[int],
    change: list[float],
    converged: list[bool],
    was_exact: bool,
    rounds: int,
    last: float,
    ok: bool,
) -> None:
    exact.append(was_exact)
    iterations.append(rounds)
    change.append(last)
    converged.append(ok)


def _alternate(
    x: FloatArray, u: IntArray, t: IntArray, n_units: int, n_times: int, tol: float, maxiter: int
) -> tuple[FloatArray, int, float]:
    """Alternate unit and period means on one column; return (residual, rounds, last change)."""
    count_u = np.bincount(u, minlength=n_units).astype(np.float64)
    count_t = np.bincount(t, minlength=n_times).astype(np.float64)
    stop = tol * max(1.0, float(np.max(np.abs(x))))
    r = x.copy()
    biggest = np.inf
    rounds = 0
    while rounds < maxiter and biggest > stop:
        rounds += 1
        step_u = np.bincount(u, weights=r, minlength=n_units) / count_u
        r -= step_u[u]
        step_t = np.bincount(t, weights=r, minlength=n_times) / count_t
        r -= step_t[t]
        biggest = float(np.max(np.abs(step_u))) if rounds > 1 else np.inf
    return r, rounds, biggest


@operator(
    "twoway",
    OpSpec(
        kind="transform",
        missing="missing where the row's own value is missing",
        stata=None,
        summary="Residual of x on unit and period effects (exact if balanced, iterated if not).",
    ),
)
def twoway(
    index: PanelIndex,
    values: FloatArray,
    *,
    tol: float = 1e-12,
    maxiter: int = 10_000,
) -> FloatArray:
    """Two-way demeaned values (see :func:`twoway_demean` for the report)."""
    return twoway_demean(index, values, tol=tol, maxiter=maxiter)[0]
