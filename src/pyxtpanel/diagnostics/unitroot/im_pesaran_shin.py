"""Im-Pesaran-Shin unit-root test for heterogeneous panels (category A, L.2.4).

Validated against Stata 17: ``xtunitroot ips`` without lags (t-bar with its exact critical
values, t-tilde-bar and Z-t-tilde-bar, with and without trend and cross-sectional means, on a
balanced and on an unbalanced panel, short and long) and with lags (W-t-bar with fixed lags and
with lags chosen by information criteria, with interpolation of the moment table or without).
The cases are in ``validation/cases/diagnostics_unitroot_im_pesaran_shin.py``.
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
from pyxtpanel.diagnostics._shared.adf import AdfCase, adf_fit, select_lags
from pyxtpanel.diagnostics._shared.tables import load_table

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.im_pesaran_shin.title": {
            "en": "Im-Pesaran-Shin unit-root test",
            "fr": "Test de racine unitaire d'Im-Pesaran-Shin",
        },
        "diag.im_pesaran_shin.null": {
            "en": "all panels contain unit roots",
            "fr": "toutes les séries du panel ont une racine unitaire",
        },
        "diag.im_pesaran_shin.alternative": {
            "en": "some panels are stationary",
            "fr": "certaines séries du panel sont stationnaires",
        },
        "diag.im_pesaran_shin.no_pvalue": {
            "en": "the normality of Z-t-tilde-bar needs at least {need} observations per unit; "
            "no p-value is given (compare t-bar with its exact critical values when the panel is "
            "balanced).",
            "fr": "la normalité de Z-t-tilde-bar exige au moins {need} observations par unité ; "
            "aucune p-value n'est donnée (comparez t-bar à ses valeurs critiques exactes si le "
            "panel est cylindré).",
        },
        "diag.im_pesaran_shin.bad_option": {
            "en": "{name} must be one of {allowed}, got {value!r}.",
            "fr": "{name} doit valoir l'une de ces valeurs : {allowed} ; reçu {value!r}.",
        },
        "diag.im_pesaran_shin.bad_lags": {
            "en": "lags must be None or a non-negative integer, got {value!r}.",
            "fr": "lags doit valoir None ou un entier positif ou nul, reçu {value!r}.",
        },
        "diag.im_pesaran_shin.selection_needs_lags": {
            "en": "lag_selection chooses among 0 to `lags` lags: give the largest number in lags.",
            "fr": "lag_selection choisit entre 0 et `lags` retards : donnez le maximum dans lags.",
        },
        "diag.im_pesaran_shin.no_moments": {
            "en": "the mean and variance of the t statistic are not tabulated for T = {t} with "
            "{lags} lag(s) (unit {unit!r}); use fewer lags or longer series.",
            "fr": "la moyenne et la variance de la statistique t ne sont pas tabulées pour T = "
            "{t} avec {lags} retard(s) (unité {unit!r}) ; utilisez moins de retards ou des séries "
            "plus longues.",
        },
        "diag.im_pesaran_shin.no_regression": {
            "en": "the Dickey-Fuller regression of unit {unit!r} cannot be run (too few "
            "observations for {lags} lag(s)).",
            "fr": "la régression de Dickey-Fuller de l'unité {unit!r} est impossible (trop peu "
            "d'observations pour {lags} retard(s)).",
        },
    },
)

_ROWS_TABLE1 = (6, 7, 8, 9, 10, 15, 20, 25, 30, 40, 50, 100, 500, 1000)
_N_KEYS = (5, 7, 10, 15, 20, 25, 50)
_T_KEYS = (5, 10, 15, 20, 25, 30, 40, 50, 60, 70)
_T_KEYS3 = (10, 15, 20, 25, 30, 40, 50, 60, 70)


@register_test(
    "im_pesaran_shin",
    TestSpec(
        family="unitroot",
        category="A",
        wave=1,
        reference="Im, K. S., Pesaran, M. H. and Shin, Y. (2003), J. of Econometrics 115(1), 53-74",
        stata="xtunitroot ips",
        validation=(
            "ips_balanced",
            "ips_short",
            "ips_trend",
            "ips_demean",
            "ips_long",
            "ips_unbalanced",
            "ips_unbalanced_short",
            "ips_lags1",
            "ips_lags_aic",
            "ips_lags_trend_unbalanced",
            "ips_lags_hqic_long",
            "ips_nointerpolate",
        ),
        inputs="panel",
        aliases=("ips",),
        balanced="both",
        tables=("ips_table1.csv", "ips_table2.csv", "ips_table3.csv"),
        r="plm::purtest(test='ips')",
        summary="Unit-root test averaging the Dickey-Fuller t statistics of the units.",
    ),
)
def im_pesaran_shin(
    panel: PanelBase,
    variable: str,
    *,
    lags: int | None = None,
    lag_selection: Literal["aic", "bic", "hqic"] | None = None,
    deterministics: Literal["constant", "trend"] = "constant",
    demean: bool = False,
    interpolate: bool = True,
) -> TestResult:
    """Return the Im-Pesaran-Shin statistic: ``Z-t-tilde-bar`` without lags, ``W-t-bar`` with lags.

    One Dickey-Fuller regression is run per unit (own autoregressive parameter).  Without
    ``lags`` the regressions have no lagged differences and the statistics are ``t-bar`` (the
    mean of the unit t ratios, with the exact critical values of the paper for a balanced
    panel), ``t-tilde-bar`` and its standardised form ``Z-t-tilde-bar``, which is the statistic
    returned; they use the moments of the paper's Table 1 at each unit's ``T``.  With ``lags``
    (a number, or the maximum when ``lag_selection`` is given) the statistic is ``W-t-bar``,
    standardised with the moments of Table 3 for that many lags.  Both are standard normal under
    the null and reject for small values.

    The panel may be unbalanced but must have no gaps.  ``result.details`` holds ``t-bar`` and,
    without lags, ``t-tilde-bar``, ``Z-t-bar`` and the 1%, 5% and 10% critical values (``NaN``
    for an unbalanced panel).

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of the numeric variable tested.
    lags
        ``None`` for regressions without lagged differences; otherwise the number of lagged
        differences (0 to 8), or their maximum with ``lag_selection``.
    lag_selection
        ``"aic"``, ``"bic"`` or ``"hqic"``: choose each unit's number of lags between 0 and
        ``lags``.
    deterministics
        ``"constant"`` (default) or ``"trend"`` (constant and linear trend).
    demean
        Remove the mean across units in each period first.
    interpolate
        Interpolate the tables of moments in ``T`` (default); with ``False`` take the entry of
        the first tabulated ``T`` at or above it.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If an option has an invalid value.
    DataError
        If a unit has a gap or too few observations.
    NumericalError
        If a regression cannot be run or the moments are not tabulated for a unit's ``T`` and
        number of lags.
    """
    # Stata: xtunitroot ips -> per unit regress D.y L.y [lags] [trend] (constant); without lags
    # t_i, ttilde_i = b / (se / rmse * sd(D.y)), moments from Table 1 at T_i (interpolated up to
    # T = 1000), Z = sqrt(N) (mean - mean of moments) / sqrt(mean of variances); with lags
    # W = sqrt(N) (tbar - mean of mu_i) / sqrt(mean of var_i) with Table 3 at (T_i, p_i).
    # Z-t-tilde-bar has a p-value when every unit has more than 5 (6 with trend) observations
    # in a balanced panel, more than 9 in an unbalanced one.
    _validate(lags, lag_selection, deterministics)
    max_lags = 0 if lags is None else int(lags)
    trend = deterministics == "trend"
    arrays = panel_variables(
        panel,
        [variable],
        test="im_pesaran_shin",
        consecutive=True,
        min_obs=2 * max_lags + 4 + int(trend),
        min_units=2,
    )
    values = arrays.values[:, 0].copy()
    if demean:
        totals = np.bincount(arrays.time_index, weights=values)
        counts = np.bincount(arrays.time_index)
        values -= (totals / np.where(counts > 0, counts, 1))[arrays.time_index]
    case: AdfCase = "trend" if trend else "constant"
    balanced = bool(np.all(arrays.counts == arrays.counts[0])) and (
        np.unique(arrays.time_index).size == arrays.counts[0]
    )
    rows: list[dict[str, float]] = []
    for unit, window in enumerate(arrays.unit_slices()):
        series = values[window]
        label = arrays.unit_labels[unit]
        p_i = max_lags
        if lags is not None and lag_selection is not None:
            p_i = select_lags(series, max_lags, lag_selection, trend=trend, constant=True)
        fit = adf_fit(series, p_i, case)
        if fit is None or not math.isfinite(fit.tau):
            raise NumericalError(
                _i18n.message("diag.im_pesaran_shin.no_regression", unit=_plain(label), lags=p_i),
                key="diag.im_pesaran_shin.no_regression",
            )
        row = {"t": float(series.size), "tau": fit.tau, "n": float(fit.n), "lags": float(p_i)}
        if lags is None:
            sd = float(np.std(np.diff(series), ddof=1))
            row["tau_tilde"] = fit.coefficient / (fit.standard_error / fit.rmse * sd)
            row.update(_table1(series.size, interpolate))
        else:
            mu, var = _table3(p_i, series.size, trend, interpolate)
            if math.isnan(mu) or math.isnan(var):
                raise NumericalError(
                    _i18n.message(
                        "diag.im_pesaran_shin.no_moments",
                        t=series.size,
                        lags=p_i,
                        unit=_plain(label),
                    ),
                    key="diag.im_pesaran_shin.no_moments",
                )
            row.update({"mu": mu, "var": var})
        rows.append(row)
    frame = pd.DataFrame(rows)
    n_units = len(frame)
    sqrt_n = math.sqrt(n_units)
    t_bar = float(frame["tau"].mean())
    quantities: dict[str, float] = {"t_bar": t_bar}
    notes: tuple[Note, ...] = ()
    if lags is None:
        statistic_name = "Z-t-tilde-bar"
        t_tilde_bar = float(frame["tau_tilde"].mean())
        z_t = sqrt_n * (t_bar - frame["mu_t"].mean()) / math.sqrt(frame["var_t"].mean())
        statistic = (
            sqrt_n * (t_tilde_bar - frame["mu_tilde"].mean()) / math.sqrt(frame["var_tilde"].mean())
        )
        minimum = int(np.min(arrays.counts))
        need = (7 if trend else 6) if balanced else 10
        pvalue = float(stats.norm.cdf(statistic)) if minimum >= need else math.nan
        if math.isnan(pvalue):
            notes = (Note("diag.im_pesaran_shin.no_pvalue", {"need": need}),)
        cv = _critical_values(n_units, float(arrays.counts.mean()), trend) if balanced else None
        quantities.update(
            {
                "t_tilde_bar": t_tilde_bar,
                "z_t_bar": float(z_t),
                "cv_1": cv[0] if cv else math.nan,
                "cv_5": cv[1] if cv else math.nan,
                "cv_10": cv[2] if cv else math.nan,
            }
        )
    else:
        statistic_name = "W-t-bar"
        statistic = sqrt_n * (t_bar - frame["mu"].mean()) / math.sqrt(frame["var"].mean())
        pvalue = float(stats.norm.cdf(statistic))
        quantities.update(
            {"mu_bar": float(frame["mu"].mean()), "var_bar": float(frame["var"].mean())}
        )
    quantities["lags"] = float(frame["lags"].mean())
    quantities["n_equations"] = float(frame["n"].sum())
    quantities["average_periods"] = float(arrays.counts.mean())
    details = pd.DataFrame(
        {"value": list(quantities.values())}, index=pd.Index(list(quantities), name="quantity")
    )
    return TestResult(
        name="im_pesaran_shin",
        statistic=float(statistic),
        pvalue=pvalue,
        statistic_name=statistic_name,
        df=(),
        n_obs=arrays.n_obs,
        n_groups=arrays.n_units,
        variable=variable,
        options={
            "lags": lags,
            "lag_selection": lag_selection,
            "deterministics": deterministics,
            "demean": demean,
            "interpolate": interpolate,
        },
        details=details,
        notes=notes,
    )


def _plain(label: object) -> object:
    return label.item() if isinstance(label, np.generic) else label


def _validate(lags: object, lag_selection: object, deterministics: object) -> None:
    """Check the options."""
    for name, value, allowed in (
        ("lag_selection", lag_selection, (None, "aic", "bic", "hqic")),
        ("deterministics", deterministics, ("constant", "trend")),
    ):
        if value not in allowed:
            raise SpecificationError(
                _i18n.message(
                    "diag.im_pesaran_shin.bad_option",
                    name=name,
                    allowed=", ".join(map(repr, allowed)),
                    value=value,
                ),
                key="diag.im_pesaran_shin.bad_option",
            )
    if lags is not None and (
        isinstance(lags, bool) or not isinstance(lags, (int, np.integer)) or lags < 0
    ):
        raise SpecificationError(
            _i18n.message("diag.im_pesaran_shin.bad_lags", value=lags),
            key="diag.im_pesaran_shin.bad_lags",
        )
    if lag_selection is not None and lags is None:
        raise SpecificationError(
            _i18n.message("diag.im_pesaran_shin.selection_needs_lags"),
            key="diag.im_pesaran_shin.selection_needs_lags",
        )


def _table1(t: int, interpolate: bool) -> dict[str, float]:
    """Return the moments of Table 1 at ``t`` periods (interpolated between 10 and 1000)."""
    table = load_table("ips_table1.csv")
    row = next((i for i, key in enumerate(_ROWS_TABLE1) if t <= key), len(_ROWS_TABLE1))
    names = {
        "mu_t": "mean_t",
        "var_t": "var_t",
        "mu_tilde": "mean_ttilde",
        "var_tilde": "var_ttilde",
    }
    columns = {k: table[v].to_numpy(dtype=float) for k, v in names.items()}
    if t <= 10 or t > 1000 or not interpolate:
        return {k: float(v[row]) for k, v in columns.items()}
    bottom, top = _ROWS_TABLE1[row - 1], _ROWS_TABLE1[row]
    below = 1.0 - (t - bottom) / (top - bottom)
    above = 1.0 - (top - t) / (top - bottom)
    return {k: float(below * v[row - 1] + above * v[row]) for k, v in columns.items()}


def _critical_values(n_units: int, periods: float, trend: bool) -> tuple[float, float, float]:
    """Return the 1%, 5% and 10% exact critical values of t-bar (Table 2)."""
    table = load_table("ips_table2.csv")
    n_key = next((k for k in _N_KEYS if n_units <= k), 100)
    t_key = next((k for k in _T_KEYS if periods <= k), 100)
    chosen = table[(table["trend"] == int(trend)) & (table["n"] == n_key) & (table["t"] == t_key)]
    values = dict(
        zip(chosen["level"].to_numpy(dtype=int), chosen["value"].to_numpy(dtype=float), strict=True)
    )
    return float(values[1]), float(values[5]), float(values[10])


def _table3(lags: int, t: int, trend: bool, interpolate: bool) -> tuple[float, float]:
    """Return the mean and variance of the ADF t statistic (Table 3) for ``lags`` and ``t``."""
    if lags < 0 or lags > 8:
        return math.nan, math.nan
    table = load_table("ips_table3.csv")
    part = table[(table["trend"] == int(trend)) & (table["lags"] == lags)]
    keys = part["t"].to_numpy(dtype=int)
    means = dict(zip(keys, part["mean"].to_numpy(dtype=float), strict=True))
    variances = dict(zip(keys, part["var"].to_numpy(dtype=float), strict=True))
    heads = (*_T_KEYS3, 100)
    column = next((i for i, key in enumerate(_T_KEYS3) if t <= key), len(heads) - 1)
    if not interpolate or column == 0 or t >= 100 or t in heads:
        key = heads[column]
        return float(means[key]), float(variances[key])
    bottom, top = heads[column - 1], heads[column]
    below = 1.0 - (t - bottom) / (top - bottom)
    above = 1.0 - (top - t) / (top - bottom)
    return (
        float(below * means[bottom] + above * means[top]),
        float(below * variances[bottom] + above * variances[top]),
    )
