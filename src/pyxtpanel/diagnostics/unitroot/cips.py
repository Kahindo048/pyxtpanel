"""Pesaran's CIPS unit-root test with cross-section dependence (category A, L.2.8).

Validated against Stata 17 through the user-written ``pescadf`` (Lewandowski, SSC): the average
of the cross-sectionally augmented Dickey-Fuller t statistics, its standardised form with its
p-value and its critical values, with a constant and with a trend, with and without the
truncation of the unit statistics, with one lag length or one per unit, on balanced and
unbalanced panels.  The unit statistics are also checked against ``plm::cipstest`` of R.  The
cases are in ``validation/cases/diagnostics_unitroot_cips.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.tables import load_table

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.cips.title": {
            "en": "Pesaran CIPS panel unit-root test (cross-sectionally augmented)",
            "fr": "Test de racine unitaire en panel CIPS de Pesaran (augmenté en coupe)",
        },
        "diag.cips.null": {
            "en": "all panels contain unit roots",
            "fr": "toutes les séries du panel ont une racine unitaire",
        },
        "diag.cips.alternative": {
            "en": "some panels are stationary",
            "fr": "certaines séries du panel sont stationnaires",
        },
        "diag.cips.bad_option": {
            "en": "{name} must be one of {allowed}, got {value!r}.",
            "fr": "{name} doit valoir l'une de ces valeurs : {allowed} ; reçu {value!r}.",
        },
        "diag.cips.bad_lags": {
            "en": "lags must be a non-negative integer or one such integer per unit, got "
            "{value!r}.",
            "fr": "lags doit être un entier positif ou nul ou un tel entier par unité, reçu "
            "{value!r}.",
        },
        "diag.cips.lags_length": {
            "en": "lags gives {given} values for {units} units: give one value or one per unit.",
            "fr": "lags donne {given} valeurs pour {units} unités : donnez une valeur ou une par "
            "unité.",
        },
        "diag.cips.no_regression": {
            "en": "the augmented regression of unit {unit!r} has no degree of freedom with "
            "{lags} lag(s) and {periods} periods; use fewer lags or longer series.",
            "fr": "la régression augmentée de l'unité {unit!r} n'a aucun degré de liberté avec "
            "{lags} retard(s) et {periods} périodes ; utilisez moins de retards ou des séries "
            "plus longues.",
        },
        "diag.cips.collinear": {
            "en": "the regressors of unit {unit!r} are collinear (for instance the cross-section "
            "mean equals the series): its t statistic is not defined.",
            "fr": "les régresseurs de l'unité {unit!r} sont colinéaires (par exemple la moyenne "
            "transversale égale la série) : sa statistique t n'est pas définie.",
        },
        "diag.cips.not_tabulated": {
            "en": "the moments of the statistic are tabulated up to 999 units and 999 periods; "
            "got {units} units and {periods} periods.",
            "fr": "les moments de la statistique sont tabulés jusqu'à 999 unités et 999 périodes ; "
            "reçu {units} unités et {periods} périodes.",
        },
        "diag.cips.dependence": {
            "en": "Z-t-bar standardises as if the unit statistics were independent, so with a "
            "strong common factor its p-value rejects too often (about 10 % at 5 % in "
            "simulations, as Stata's pescadf); for a balanced panel decide with t-bar against its "
            "critical values (options cv_1, cv_5, cv_10).",
            "fr": "Z-t-bar standardise comme si les statistiques des unités étaient "
            "indépendantes : avec un facteur commun fort sa p-value rejette trop souvent (environ "
            "10 % à 5 % en simulation, comme pescadf de Stata) ; pour un panel cylindré, décidez "
            "avec t-bar et ses valeurs critiques (options cv_1, cv_5, cv_10).",
        },
        "diag.cips.no_critical": {
            "en": "the critical values of t-bar are given for balanced panels only; compare "
            "Z-t-bar with the normal distribution (its p-value).",
            "fr": "les valeurs critiques de t-bar ne sont données que pour un panel cylindré ; "
            "comparez Z-t-bar à la loi normale (sa p-value).",
        },
    },
)

#: Values of N and T the tables of moments and critical values are tabulated for.
KEYS = (10, 15, 20, 30, 50, 70, 100, 200, 999)

#: Truncation bounds (K1, K2) of the unit statistics for a constant and for a constant and trend.
TRUNCATION = {"constant": (-6.19, 2.61), "trend": (-6.42, 1.70)}


@register_test(
    "cips",
    TestSpec(
        family="unitroot",
        category="A",
        wave=2,
        reference="Pesaran, M. H. (2007), Journal of Applied Econometrics 22(2), 265-312",
        stata="pescadf",
        validation=(
            "cips_lag0",
            "cips_lag1",
            "cips_lag2_trend",
            "cips_nodemean",
            "cips_notrunc",
            "cips_lag3_trend_notrunc",
            "cips_factor",
            "cips_factor_trend",
            "cips_unbalanced",
            "cips_unbalanced_lags",
            "cips_per_unit_lags",
            "cips_small_n",
            "cips_short",
        ),
        inputs="panel",
        aliases=("pescadf", "pesaran_cips"),
        balanced="both",
        tables=("cips_moments.csv", "cips_critical.csv"),
        r="plm::cipstest",
        summary="Unit-root test averaging the cross-sectionally augmented Dickey-Fuller t "
        "statistics.",
    ),
)
def cips(
    panel: PanelBase,
    variable: str,
    *,
    lags: int | Sequence[int],
    deterministics: Literal["constant", "trend"] = "constant",
    truncate: bool = True,
) -> TestResult:
    """Return Pesaran's ``Z[t-bar]`` of the cross-sectionally augmented Dickey-Fuller regressions.

    For each unit the first difference of the series is regressed on a constant, the lagged
    level, the lagged cross-section mean of the series and its first difference, ``lags`` lagged
    differences of the series and of the cross-section mean, and a trend if requested; the unit
    statistic is the t ratio of the lagged level.  The cross-section mean (over the units
    observed in each period) absorbs a common factor, which is what makes the test valid when the
    units are correlated.  By default the series are first taken in deviation from their
    cross-section mean in the first period (the paper's advice for a fixed ``T``, which
    ``pescadf`` does unless told ``nodemean``): a constant shift of the series, which the constant
    of each regression absorbs, so it changes nothing and has no option here.  Each unit
    statistic is truncated to ``[K1, K2]`` (Pesaran's bounds) so that its moments are finite.

    ``t-bar`` is the mean of the unit statistics (the CIPS statistic); ``Z-t-bar`` is the sum of
    their deviations from the tabulated moments over the root of the sum of the tabulated
    variances, standard normal under the null and rejecting for small values; it is the statistic
    returned, with its p-value.  The moments and the critical values are the tables of the paper
    read at the first tabulated ``N`` and ``T`` at or above the panel's, as ``pescadf`` does.  The
    critical values need a balanced panel.  ``details`` has one row per unit; ``options`` carries
    ``t_bar`` and the 1 %, 5 % and 10 % critical values (``NaN`` for an unbalanced panel).

    Parameters
    ----------
    panel
        The declared panel (may be unbalanced, without gaps: ``pescadf`` also runs on a panel with
        holes, this test refuses one).
    variable
        Name of the numeric variable tested.
    lags
        Number of lagged differences: one integer for all units, or one per unit in the order of
        the units.
    deterministics
        ``"constant"`` (default) or ``"trend"`` (constant and linear trend).
    truncate
        Truncate the unit statistics (default).

    Raises
    ------
    SpecificationError
        If an option has an invalid value or ``lags`` does not match the units.
    DataError
        If a unit has a gap or too few observations, or the panel exceeds the tables.
    NumericalError
        If a unit's regression is degenerate.
    """
    # Stata: pescadf y, lags(p) [trend] [trunoff] (nodemean changes nothing) -> per unit
    # regress D.y L.y L.ybar D.ybar L(1/p)D.y L(1/p)D.ybar [trend], t = b / se of L.y, truncated
    # to [K1, K2] = (-6.19, 2.61) (constant) or (-6.42, 1.70) (trend); Z = sum(t_i - mu_i) /
    # sqrt(sum sd_i^2) with mu_i, sd_i of Pesaran's Table II at the first tabulated N and T_i at
    # or above; p = normal(Z); t-bar critical values of Table I at (N, mean T), balanced panel.
    _validate(deterministics)
    trend = deterministics == "trend"
    arrays = panel_variables(panel, [variable], test="cips", consecutive=True, min_units=2)
    per_unit = _unit_lags(lags, arrays.n_units)
    values = arrays.values[:, 0].copy()
    calendar = arrays.time_index
    counts_by_period = np.bincount(calendar)
    totals = np.bincount(calendar, weights=values)
    cross_mean = totals / np.where(counts_by_period > 0, counts_by_period, 1)
    lower, upper = TRUNCATION[deterministics]
    rows: list[dict[str, float]] = []
    for unit, window in enumerate(arrays.unit_slices()):
        label = arrays.unit_labels[unit]
        series = values[window]
        if series.size < 6 + 3 * int(per_unit[unit]) + int(trend):
            raise NumericalError(
                _i18n.message(
                    "diag.cips.no_regression",
                    unit=_plain(label),
                    lags=int(per_unit[unit]),
                    periods=series.size,
                ),
                key="diag.cips.no_regression",
            )
        raw = _unit_statistic(series, cross_mean[calendar[window]], int(per_unit[unit]), trend)
        if raw is None:
            raise NumericalError(
                _i18n.message("diag.cips.collinear", unit=_plain(label)), key="diag.cips.collinear"
            )
        used = min(max(raw, lower), upper) if truncate else raw
        rows.append(
            {
                "periods": float(series.size),
                "lags": float(per_unit[unit]),
                "n_obs": float(series.size - 1 - per_unit[unit]),
                "t_raw": raw,
                "t": used,
            }
        )
    frame = pd.DataFrame(rows, index=pd.Index(arrays.unit_labels, name="unit"))
    n_units = len(frame)
    moments = _moments(deterministics, n_units, frame["periods"].to_numpy())
    frame["mean"], frame["sd"] = moments[:, 0], moments[:, 1]
    statistic = float((frame["t"] - frame["mean"]).sum() / math.sqrt((frame["sd"] ** 2).sum()))
    t_bar = float(frame["t"].mean())
    balanced = bool(np.all(arrays.counts == arrays.counts[0]))
    critical = (
        _critical_values(deterministics, truncate, n_units, float(arrays.counts.mean()))
        if balanced
        else (math.nan,) * 3
    )
    return TestResult(
        name="cips",
        statistic=statistic,
        pvalue=float(stats.norm.cdf(statistic)),
        statistic_name="Z-t-bar",
        df=(),
        n_obs=int(frame["n_obs"].sum()),
        n_groups=n_units,
        variable=variable,
        options={
            "lags": lags if isinstance(lags, (int, np.integer)) else list(lags),
            "deterministics": deterministics,
            "truncate": truncate,
            "t_bar": t_bar,
            "cv_1": critical[0],
            "cv_5": critical[1],
            "cv_10": critical[2],
        },
        details=frame,
        notes=(
            (Note("diag.cips.dependence"),)
            if balanced
            else (Note("diag.cips.dependence"), Note("diag.cips.no_critical"))
        ),
    )


def _plain(label: object) -> object:
    return label.item() if isinstance(label, np.generic) else label


def _validate(deterministics: object) -> None:
    if deterministics not in ("constant", "trend"):
        raise SpecificationError(
            _i18n.message(
                "diag.cips.bad_option",
                name="deterministics",
                allowed="'constant', 'trend'",
                value=deterministics,
            ),
            key="diag.cips.bad_option",
        )


def _unit_lags(lags: object, n_units: int) -> np.ndarray:
    """Return the lag order of each unit (one value for all, or one per unit)."""
    values: list[Any]
    if isinstance(lags, (int, np.integer)) and not isinstance(lags, bool):
        values = [int(lags)]
    elif isinstance(lags, Sequence) and not isinstance(lags, str):
        values = [v for v in lags]
    else:
        values = [lags]
    if not values or any(
        isinstance(v, bool) or not isinstance(v, (int, np.integer)) or v < 0 for v in values
    ):
        raise SpecificationError(
            _i18n.message("diag.cips.bad_lags", value=lags), key="diag.cips.bad_lags"
        )
    if len(values) not in (1, n_units):
        raise SpecificationError(
            _i18n.message("diag.cips.lags_length", given=len(values), units=n_units),
            key="diag.cips.lags_length",
        )
    return np.full(n_units, values[0], dtype=np.int64) if len(values) == 1 else np.array(values)


def _unit_statistic(
    series: FloatArray, cross_mean: FloatArray, lags: int, trend: bool
) -> float | None:
    """Return the t ratio of the lagged level in the CADF regression, ``None`` if degenerate."""
    p = lags
    d_y, d_m = np.diff(series), np.diff(cross_mean)
    sample = np.arange(p + 1, series.size)  # periods with the lagged level and p lagged differences
    response = d_y[sample - 1]
    columns = [series[sample - 1], cross_mean[sample - 1], d_m[sample - 1]]
    for j in range(1, p + 1):
        columns += [d_y[sample - 1 - j], d_m[sample - 1 - j]]
    if trend:
        columns.append(sample.astype(float))
    design = np.column_stack([np.ones(sample.size), *columns])
    n, k = design.shape
    if n <= k or np.linalg.matrix_rank(design) < k:
        return None
    # t ratio of the first regressor by partialling out the others (Frisch-Waugh)
    others = np.delete(design, 1, axis=1)
    x1 = design[:, 1]
    basis, _ = np.linalg.qr(others)
    x_tilde = x1 - basis @ (basis.T @ x1)
    y_tilde = response - basis @ (basis.T @ response)
    sxx = float(x_tilde @ x_tilde)
    coefficient = float(x_tilde @ y_tilde) / sxx
    residual = y_tilde - coefficient * x_tilde
    sigma2 = float(residual @ residual) / (n - k)
    if sigma2 <= 0.0 or not math.isfinite(sigma2):
        return None
    return coefficient / math.sqrt(sigma2 / sxx)


def _first_key(value: float) -> int | None:
    """Return the first tabulated key at or above ``value`` (``None`` above 999)."""
    return next((key for key in KEYS if value <= key), None)


def _moments(case: str, n_units: int, periods: np.ndarray) -> np.ndarray:
    """Return the tabulated mean and standard deviation of each unit's statistic."""
    table = load_table("cips_moments.csv")
    table = table[table["case"] == case]
    lookup = {
        (int(t), int(n)): (m, s)
        for t, n, m, s in zip(table["t"], table["n"], table["mean"], table["sd"], strict=True)
    }
    column = _first_key(n_units)
    out = np.empty((periods.size, 2))
    for i, t in enumerate(periods):
        row = _first_key(float(t))
        if column is None or row is None:
            raise DataError(
                _i18n.message("diag.cips.not_tabulated", units=n_units, periods=int(t)),
                key="diag.cips.not_tabulated",
            )
        out[i] = lookup[(row, column)]
    return out


def _critical_values(
    case: str, truncated: bool, n_units: int, mean_periods: float
) -> tuple[float, float, float]:
    """Return the 1 %, 5 % and 10 % critical values of t-bar at ``(N, T)`` for a balanced panel."""
    table = load_table("cips_critical.csv")
    table = table[(table["case"] == case) & (table["truncated"] == int(truncated))]
    row, column = _first_key(mean_periods), _first_key(n_units)
    if row is None or column is None:
        raise DataError(
            _i18n.message(
                "diag.cips.not_tabulated", units=n_units, periods=int(math.ceil(mean_periods))
            ),
            key="diag.cips.not_tabulated",
        )
    found = table[(table["t"] == row) & (table["n"] == column)].iloc[0]
    return float(found["cv_1"]), float(found["cv_5"]), float(found["cv_10"])
