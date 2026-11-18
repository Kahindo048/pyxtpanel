"""Size and power of the panel and residual tests by simulation (plans L.7.5 and L.8.12).

Each scenario draws data from a known process, runs one test and records whether it rejects at
5 %.  Under the null the rejection rate (the *size*) must be close to the level; under the
alternative (the *power*) it must be high.  Seeds are fixed, so a run is reproducible.

    python -m validation.simulate --reps 2000 --write docs/size_power.md

writes the table of the documentation; ``tests/integration/test_size_power.py`` runs the same
scenarios with fewer replications and checks the documented bounds.
"""

from __future__ import annotations

import argparse
import warnings
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import pyxtpanel as xt
from pyxtpanel.diagnostics import (
    baltagi_li,
    breitung,
    breusch_pagan,
    cips,
    dumitrescu_hurlin,
    f_effects,
    fisher,
    frees,
    friedman,
    hadri,
    harris_tzavalis,
    im_pesaran_shin,
    jarque_bera,
    kao,
    levin_lin_chu,
    lm_cross_section,
    lm_random_effects,
    panel_durbin_watson,
    pedroni,
    pesaran_cd,
    reset,
    shapiro_francia,
    shapiro_wilk,
    skewness_kurtosis,
    slope_homogeneity,
    wald_groupwise,
    westerlund,
    white,
    wooldridge,
)

LEVEL = 0.05

#: A scenario returns the p-value (or the statistic for tests without one) of one replication.
Draw = Callable[[np.random.Generator], float]


@dataclass(frozen=True)
class Study:
    """One test: how it is simulated under the null and under an alternative.

    Attributes
    ----------
    null, alternative
        Functions that draw one replication and return its p-value.
    size_bounds
        Interval the rejection rate under the null must fall in (wide enough for the Monte Carlo
        error of the replications used by the automatic test, narrow enough to catch a broken
        statistic).  A test known to be oversized has bounds that say so.
    min_power
        Smallest rejection rate under the alternative.
    note
        What the scenario is.
    """

    null: Draw
    alternative: Draw
    size_bounds: tuple[float, float]
    min_power: float
    note: str


# ---------------------------------------------------------------- data generators


def _ols_residuals(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    z = np.column_stack([np.ones(y.size), x])
    fitted = z @ np.linalg.lstsq(z, y, rcond=None)[0]
    return y - fitted, fitted


def _regression(rng: np.random.Generator, n: int, *, scale: str) -> tuple[np.ndarray, ...]:
    """Residuals, fitted values and regressors of a linear model, with a chosen error scale."""
    x = rng.normal(size=(n, 2))
    sd = np.exp(0.5 * x[:, 0]) if scale == "heteroskedastic" else 1.0
    y = 1.0 + x @ np.array([1.0, -0.5]) + sd * rng.normal(size=n)
    e, fitted = _ols_residuals(y, x)
    return e, fitted, x


def _breusch_pagan(scale: str) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        e, fitted, _ = _regression(rng, 300, scale=scale)
        return breusch_pagan(e, fitted).pvalue

    return draw


def _white(scale: str) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        e, _, x = _regression(rng, 300, scale=scale)
        return white(e, x).pvalue

    return draw


def _within_residuals(
    rng: np.random.Generator, n_units: int, periods: int, sd: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    e = rng.normal(size=(n_units, periods)) * sd[:, None]
    e = e - e.mean(axis=1, keepdims=True)
    return e.ravel(), np.arange(n_units + 1) * periods


def _wald_groupwise(*, heteroskedastic: bool) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        sd = np.linspace(0.5, 2.0, 20) if heteroskedastic else np.ones(20)
        e, offsets = _within_residuals(rng, 20, 30, sd)
        return wald_groupwise(e, offsets).pvalue

    return draw


def _pooled_residuals(
    rng: np.random.Generator, n_units: int, periods: int, effect_sd: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    unit = np.repeat(np.arange(n_units), periods)
    x = rng.normal(size=n_units * periods)
    y = 1.0 + x + effect_sd * rng.normal(size=n_units)[unit] + rng.normal(size=n_units * periods)
    e, _ = _ols_residuals(y, x[:, None])
    return y, x, e


def _lm_random_effects(effect_sd: float) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        _, _, e = _pooled_residuals(rng, 40, 5, effect_sd)
        # the random-effects estimator gives Var(u); its sign is that of this moment estimate
        variance_u = float(np.sum(e.reshape(40, 5).sum(axis=1) ** 2) - e @ e)
        return lm_random_effects(e, np.arange(41) * 5, variance_u=variance_u).pvalue

    return draw


def _baltagi_li(statistic: str, effect_sd: float, rho: float) -> Draw:
    """p-value of a Baltagi-Li statistic: 40 units x 8 periods, unit effect and AR(1) errors."""

    def draw(rng: np.random.Generator) -> float:
        n_units, periods = 40, 8
        u = np.repeat(rng.normal(size=n_units) * effect_sd, periods)
        v = np.zeros((n_units, periods))
        for t in range(periods):
            v[:, t] = (rho * v[:, t - 1] if t else 0.0) + rng.normal(size=n_units)
        x = rng.normal(size=n_units * periods)
        y = 1.0 + x + u + v.ravel()
        design = np.column_stack([np.ones(y.size), x])
        e = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
        return baltagi_li(e, np.arange(n_units + 1) * periods, statistic=statistic).pvalue

    return draw


def _f_effects(effect_sd: float) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        n_units, periods = 40, 5
        y, x, pooled = _pooled_residuals(rng, n_units, periods, effect_sd)
        yc = y.reshape(n_units, periods) - y.reshape(n_units, periods).mean(axis=1, keepdims=True)
        xc = x.reshape(n_units, periods) - x.reshape(n_units, periods).mean(axis=1, keepdims=True)
        slope = float(xc.ravel() @ yc.ravel()) / float(xc.ravel() @ xc.ravel())
        within = yc.ravel() - slope * xc.ravel()
        return f_effects(
            float(pooled @ pooled), float(within @ within), n_units * periods, n_units, 1
        ).pvalue

    return draw


def _normality(test: Callable[..., object], *, heavy: bool) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        x = rng.standard_t(df=4, size=200) if heavy else rng.normal(size=200)
        return float(test(x).pvalue)  # type: ignore[attr-defined]

    return draw


def _reset(*, curved: bool) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        n = 200
        x = rng.normal(size=(n, 2))
        y = 1.0 + x[:, 0] - x[:, 1] + (0.5 * x[:, 0] ** 2 if curved else 0.0) + rng.normal(size=n)
        return reset(y, x).pvalue

    return draw


def _slope_panel(rng: np.random.Generator, *, heterogeneous: bool) -> xt.PanelData:
    n_units, periods = 30, 20
    unit = np.repeat(np.arange(1, n_units + 1), periods)
    spread = 0.6 * rng.normal(size=n_units) if heterogeneous else np.zeros(n_units)
    slope = (1.0 + spread)[unit - 1]
    x = rng.normal(size=n_units * periods)
    y = rng.normal(size=n_units)[unit - 1] + slope * x + rng.normal(size=n_units * periods)
    frame = pd.DataFrame({"id": unit, "t": np.tile(np.arange(1, periods + 1), n_units)})
    frame["x"], frame["y"] = x, y
    return xt.xtset(frame, "id", "t")


def _slope_homogeneity(*, heterogeneous: bool) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        return slope_homogeneity(_slope_panel(rng, heterogeneous=heterogeneous), "y", ["x"]).pvalue

    return draw


def _slope_homogeneity_hac(kernel: str, *, heterogeneous: bool) -> Draw:
    """HAC version: 40 units x 30 periods, unit-specific error variances and AR(0.4) errors."""

    def draw(rng: np.random.Generator) -> float:
        n_units, periods = 40, 30
        unit = np.repeat(np.arange(n_units), periods)
        scale = np.exp(0.6 * rng.normal(size=n_units))[unit]
        e = np.zeros((n_units, periods))
        for k in range(periods):
            e[:, k] = (0.4 * e[:, k - 1] if k else 0.0) + rng.normal(size=n_units)
        x = rng.normal(size=n_units * periods)
        slope = 1.0 + (0.6 * rng.normal(size=n_units)[unit] if heterogeneous else 0.0)
        y = rng.normal(size=n_units)[unit] + slope * x + scale * e.ravel()
        frame = pd.DataFrame({"id": unit + 1, "t": np.tile(np.arange(1, periods + 1), n_units)})
        frame["x"], frame["y"] = x, y
        panel = xt.xtset(frame, "id", "t")
        return slope_homogeneity(panel, "y", ["x"], hac=True, kernel=kernel).pvalue  # type: ignore[arg-type]

    return draw


def _durbin_watson(rho: float) -> Draw:
    """The statistic itself (it has no p-value): its mean over replications is reported."""

    def draw(rng: np.random.Generator) -> float:
        n_units, periods = 50, 10
        e = np.empty((n_units, periods))
        e[:, 0] = rng.normal(size=n_units)
        for k in range(1, periods):
            e[:, k] = rho * e[:, k - 1] + np.sqrt(1 - rho**2) * rng.normal(size=n_units)
        x = rng.normal(size=n_units * periods)
        unit = np.repeat(np.arange(1, n_units + 1), periods)
        y = rng.normal(size=n_units)[unit - 1] + x + e.ravel()
        frame = pd.DataFrame({"id": unit, "t": np.tile(np.arange(1, periods + 1), n_units)})
        frame["x"], frame["y"] = x, y
        return panel_durbin_watson(xt.xtset(frame, "id", "t"), "y", ["x"]).statistic

    return draw


# ---------------------------------------------------------------- panel tests (plan L.7.5)


def _declared(columns: dict[str, np.ndarray]) -> xt.PanelData:
    """Declare a panel from ``(units, periods)`` arrays (one per variable), units 1..N."""
    n_units, periods = next(iter(columns.values())).shape
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n_units + 1), periods),
            "t": np.tile(np.arange(1, periods + 1), n_units),
        }
    )
    for name, values in columns.items():
        frame[name] = values.ravel()
    return xt.xtset(frame, "id", "t")


def _series(
    rng: np.random.Generator, n_units: int, periods: int, *, rho: float, scale: float = 1.0
) -> np.ndarray:
    """AR(1) series with a unit mean (``rho = 1``: a random walk started at that level).

    For ``rho < 1`` the series starts from its stationary distribution.
    """
    shocks = scale * rng.normal(size=(n_units, periods))
    level = rng.normal(scale=2.0, size=n_units)
    y = np.empty((n_units, periods))
    y[:, 0] = shocks[:, 0] / np.sqrt(1.0 - rho**2) if rho < 1.0 else shocks[:, 0]
    for k in range(1, periods):
        y[:, k] = rho * y[:, k - 1] + shocks[:, k]
    return y + level[:, None]


def _unit_root(procedure: Callable[..., object], *, rho: float, **options: object) -> Draw:
    """A unit-root test on 30 series of 50 periods with autoregressive root ``rho``."""

    def draw(rng: np.random.Generator) -> float:
        y = _series(rng, 30, 50, rho=rho)
        return float(procedure(_declared({"y": y}), "y", **options).pvalue)  # type: ignore[attr-defined]

    return draw


def _factor_series(rng: np.random.Generator, rho: float) -> np.ndarray:
    """30 series of 50 periods whose innovations share a common factor (loadings around 1).

    The units are strongly correlated, which is what the tests that ignore it cannot handle.
    """
    n_units, periods = 30, 50
    common = rng.normal(size=periods)
    loading = 1.0 + 0.5 * rng.normal(size=(n_units, 1))
    shocks = rng.normal(size=(n_units, periods)) + loading * common
    y = np.empty((n_units, periods))
    y[:, 0] = shocks[:, 0]
    for k in range(1, periods):
        y[:, k] = rho * y[:, k - 1] + shocks[:, k]
    return y + rng.normal(scale=2.0, size=(n_units, 1))


def _factor_unit_root(procedure: Callable[..., object], *, rho: float, **options: object) -> Draw:
    """A unit-root test applied to the series of :func:`_factor_series`."""

    def draw(rng: np.random.Generator) -> float:
        y = _factor_series(rng, rho)
        return float(procedure(_declared({"y": y}), "y", **options).pvalue)  # type: ignore[attr-defined]

    return draw


def _cips_critical(*, rho: float) -> Draw:
    """CIPS decided as the paper does: t-bar against its 5 % critical value (0 = reject)."""

    def draw(rng: np.random.Generator) -> float:
        result = cips(_declared({"y": _factor_series(rng, rho)}), "y", lags=1)
        return 0.0 if result.options["t_bar"] < result.options["cv_5"] else 1.0  # type: ignore[operator]

    return draw


def _hadri(*, share: float) -> Draw:
    """Hadri's null is stationarity: white noise around unit means.

    Under the alternative a share of the units are random walks (the others stay white noise).
    """

    def draw(rng: np.random.Generator) -> float:
        y = _series(rng, 30, 50, rho=0.0)
        walks = int(round(share * 30))
        if walks:
            y[:walks] = _series(rng, walks, 50, rho=1.0)
        return hadri(_declared({"y": y}), "y").pvalue

    return draw


def _kao(*, rho: float) -> Draw:
    """``y = x + u`` with ``x`` a random walk; ``u`` has root ``rho`` (1: no cointegration)."""

    def draw(rng: np.random.Generator) -> float:
        n_units, periods = 30, 40
        shocks = rng.normal(size=(n_units, periods, 2))
        x = np.cumsum(shocks[:, :, 0] + 0.5 * shocks[:, :, 1], axis=1)
        u = np.empty((n_units, periods))
        u[:, 0] = shocks[:, 0, 1]
        for k in range(1, periods):
            u[:, k] = rho * u[:, k - 1] + shocks[:, k, 1]
        y = rng.normal(scale=2.0, size=(n_units, 1)) + x + u
        return kao(_declared({"y": y, "x": x}), "y", ["x"]).pvalue

    return draw


def _cointegration_panel(rng: np.random.Generator, rho: float) -> xt.PanelData:
    """20 units x 40 periods: ``y = x + u``, ``x`` a random walk, ``u`` with root ``rho``."""
    n_units, periods = 20, 40
    shocks = rng.normal(size=(n_units, periods, 2))
    x = np.cumsum(shocks[:, :, 0] + 0.5 * shocks[:, :, 1], axis=1)
    u = np.empty((n_units, periods))
    u[:, 0] = shocks[:, 0, 1]
    for k in range(1, periods):
        u[:, k] = rho * u[:, k - 1] + shocks[:, k, 1]
    y = rng.normal(scale=2.0, size=(n_units, 1)) + x + u
    return _declared({"y": y, "x": x})


def _westerlund(*, rho: float) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        return westerlund(_cointegration_panel(rng, rho), "y", ["x"]).pvalue

    return draw


def _pedroni(*, rho: float, ar: str, statistic: str) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        panel = _cointegration_panel(rng, rho)
        return pedroni(
            panel,
            "y",
            ["x"],
            kernel_lags=2,
            ar=ar,
            statistic=statistic,  # type: ignore[arg-type]
        ).pvalue

    return draw


def _granger(*, slope: float) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        n_units, periods = 20, 30
        x = np.empty((n_units, periods))
        y = np.empty((n_units, periods))
        x[:, 0], y[:, 0] = rng.normal(size=n_units), rng.normal(size=n_units)
        for k in range(1, periods):
            x[:, k] = 0.5 * x[:, k - 1] + rng.normal(size=n_units)
            y[:, k] = 0.4 * y[:, k - 1] + slope * x[:, k - 1] + rng.normal(size=n_units)
        return dumitrescu_hurlin(_declared({"y": y, "x": x}), "y", "x").pvalue

    return draw


def _cross_section(
    test: Callable[..., object],
    *,
    loading: float,
    shape: tuple[int, int] = (10, 30),
    frees_rule: bool = False,
    **options: object,
) -> Draw:
    """Residuals of ``shape`` = (units, periods) with a common factor of mean loading ``loading``.

    Frees' test has critical values instead of a p-value: the draw is 0 when the statistic
    exceeds the 5 % value and 1 otherwise, so that its rejection rate is counted like the others.
    """

    def draw(rng: np.random.Generator) -> float:
        e = rng.normal(size=shape)
        if loading:
            e += (loading + 0.5 * loading * rng.normal(size=(shape[0], 1))) * rng.normal(
                size=(1, shape[1])
            )
        result = test(_declared({"e": e}), "e", **options)
        if frees_rule:
            critical = result.details.loc["cv_5", "value"]  # type: ignore[attr-defined]
            return 0.0 if result.statistic > critical else 1.0  # type: ignore[attr-defined]
        return float(result.pvalue)  # type: ignore[attr-defined]

    return draw


def _wooldridge(*, rho: float) -> Draw:
    def draw(rng: np.random.Generator) -> float:
        n_units, periods = 40, 10
        e = _series(rng, n_units, periods, rho=rho)
        x = rng.normal(size=(n_units, periods))
        y = rng.normal(size=(n_units, 1)) + x + e
        return wooldridge(_declared({"y": y, "x": x}), "y", ["x"]).pvalue

    return draw


STUDIES: dict[str, Study] = {
    "breusch_pagan": Study(
        _breusch_pagan("homoskedastic"),
        _breusch_pagan("heteroskedastic"),
        (0.01, 0.10),
        0.8,
        "n = 300, fitted values; error sd proportional to exp(x1 / 2) under the alternative",
    ),
    "white": Study(
        _white("homoskedastic"),
        _white("heteroskedastic"),
        (0.01, 0.10),
        0.8,
        "n = 300, two regressors; error sd proportional to exp(x1 / 2) under the alternative",
    ),
    "wald_groupwise": Study(
        _wald_groupwise(heteroskedastic=False),
        _wald_groupwise(heteroskedastic=True),
        (0.15, 1.0),
        0.9,
        "20 units x 30 periods; oversized at this T, as xttest3 (documented in plan 8.3); "
        "unit sd from 0.5 to 2 under the alternative",
    ),
    "lm_random_effects": Study(
        _lm_random_effects(0.0),
        _lm_random_effects(1.0),
        (0.01, 0.10),
        0.9,
        "40 units x 5 periods; unit effect with sd 1 under the alternative; the sign rule of "
        "xttest0 applied through variance_u (without it the size is 10 %)",
    ),
    "baltagi_li": Study(
        _baltagi_li("joint", 0.0, 0.0),
        _baltagi_li("joint", 1.0, 0.0),
        (0.01, 0.10),
        0.9,
        "40 units x 8 periods; independent errors under the null, a unit effect with sd 1 under "
        "the alternative; the joint statistic (chi-square with 2 degrees of freedom)",
    ),
    "baltagi_li_serial_adjusted": Study(
        _baltagi_li("serial_adjusted", 0.5, 0.0),
        _baltagi_li("serial_adjusted", 0.5, 0.3),
        (0.01, 0.12),
        0.9,
        "40 units x 8 periods with a unit effect of sd 0.5 in both: AR(1) errors with root 0.3 "
        "under the alternative; robust to the effect only when it is small (local), the size "
        "is about 21 % with an effect of sd 1",
    ),
    "baltagi_li_re_adjusted": Study(
        _baltagi_li("re_adjusted_onesided", 0.0, 0.1),
        _baltagi_li("re_adjusted_onesided", 1.0, 0.1),
        (0.01, 0.10),
        0.9,
        "40 units x 8 periods with AR(1) errors of root 0.1 in both: a unit effect of sd 1 under "
        "the alternative; one-sided robust statistic, robust to serial correlation only when it "
        "is small (local), the size is about 20 % with root 0.3 and 67 % with root 0.5",
    ),
    "f_effects": Study(
        _f_effects(0.0),
        _f_effects(1.0),
        (0.01, 0.10),
        0.9,
        "40 units x 5 periods; unit effect with sd 1 under the alternative",
    ),
    "skewness_kurtosis": Study(
        _normality(skewness_kurtosis, heavy=False),
        _normality(skewness_kurtosis, heavy=True),
        (0.01, 0.10),
        0.7,
        "n = 200; Student t with 4 degrees of freedom under the alternative",
    ),
    "jarque_bera": Study(
        _normality(jarque_bera, heavy=False),
        _normality(jarque_bera, heavy=True),
        (0.005, 0.10),
        0.6,
        "n = 200; Student t with 4 degrees of freedom under the alternative",
    ),
    "shapiro_wilk": Study(
        _normality(shapiro_wilk, heavy=False),
        _normality(shapiro_wilk, heavy=True),
        (0.01, 0.10),
        0.8,
        "n = 200; Student t with 4 degrees of freedom under the alternative",
    ),
    "shapiro_francia": Study(
        _normality(shapiro_francia, heavy=False),
        _normality(shapiro_francia, heavy=True),
        (0.01, 0.10),
        0.8,
        "n = 200; Student t with 4 degrees of freedom under the alternative",
    ),
    "reset": Study(
        _reset(curved=False),
        _reset(curved=True),
        (0.01, 0.10),
        0.6,
        "n = 200, two regressors; an omitted 0.5 * x1^2 under the alternative",
    ),
    "slope_homogeneity": Study(
        _slope_homogeneity(heterogeneous=False),
        _slope_homogeneity(heterogeneous=True),
        (0.005, 0.12),
        0.9,
        "30 units x 20 periods; slopes with sd 0.6 around 1 under the alternative",
    ),
    "slope_homogeneity_hac_bartlett": Study(
        _slope_homogeneity_hac("bartlett", heterogeneous=False),
        _slope_homogeneity_hac("bartlett", heterogeneous=True),
        (0.01, 0.12),
        0.9,
        "HAC version, Bartlett kernel with the automatic bandwidth; 40 units x 30 periods, "
        "unit-specific error variances and AR(1) errors with root 0.4; slopes with sd 0.6 around "
        "1 under the alternative",
    ),
    "slope_homogeneity_hac_qs": Study(
        _slope_homogeneity_hac("qs", heterogeneous=False),
        _slope_homogeneity_hac("qs", heterogeneous=True),
        (0.01, 0.12),
        0.9,
        "as the Bartlett version with the quadratic spectral kernel",
    ),
    "breitung": Study(
        _unit_root(breitung, rho=1.0),
        _unit_root(breitung, rho=0.95),
        (0.005, 0.10),
        0.9,
        "30 units x 50 periods; random walks under the null, AR(1) with root 0.95 under the "
        "alternative; unit means",
    ),
    "harris_tzavalis": Study(
        _unit_root(harris_tzavalis, rho=1.0),
        _unit_root(harris_tzavalis, rho=0.95),
        (0.03, 0.15),
        0.9,
        "30 units x 50 periods, root 0.95 under the alternative; oversized (about 9 %) at this "
        "T, as Stata's xtunitroot ht: the moments of the statistic are those of a fixed T",
    ),
    "levin_lin_chu": Study(
        _unit_root(levin_lin_chu, rho=1.0),
        _unit_root(levin_lin_chu, rho=0.95),
        (0.15, 0.35),
        0.3,
        "30 units x 50 periods, one ADF lag, root 0.95 under the alternative; oversized (about "
        "25 %) with the default kernel, as Stata's xtunitroot llc (24.7 % in 300 replications "
        "of the same design)",
    ),
    "im_pesaran_shin": Study(
        _unit_root(im_pesaran_shin, rho=1.0, lags=1),
        _unit_root(im_pesaran_shin, rho=0.95, lags=1),
        (0.01, 0.10),
        0.4,
        "30 units x 50 periods, one ADF lag, root 0.95 under the alternative",
    ),
    "fisher": Study(
        _unit_root(fisher, rho=1.0, lags=1),
        _unit_root(fisher, rho=0.95, lags=1),
        (0.03, 0.14),
        0.3,
        "30 units x 50 periods, Dickey-Fuller with one lag, inverse chi-square combination, "
        "root 0.95 under the alternative; a little oversized (about 8 %) as Stata's xtunitroot "
        "fisher (7.8 % in 1000 replications of the same design)",
    ),
    "fisher_phillips_perron": Study(
        _unit_root(fisher, rho=1.0, lags=2, test="pperron"),
        _unit_root(fisher, rho=0.95, lags=2, test="pperron"),
        (0.03, 0.17),
        0.3,
        "as fisher with the Phillips-Perron unit test (two Newey-West lags)",
    ),
    "cips": Study(
        _cips_critical(rho=1.0),
        _cips_critical(rho=0.85),
        (0.01, 0.11),
        0.7,
        "30 units x 50 periods, innovations with a strong common factor (loadings of mean 1), "
        "one lag; rejection when t-bar is below its 5 % critical value (the paper's rule); root "
        "0.85 under the alternative",
    ),
    "cips_z": Study(
        _factor_unit_root(cips, rho=1.0, lags=1),
        _factor_unit_root(cips, rho=0.85, lags=1),
        (0.07, 0.20),
        0.7,
        "as cips with the p-value of Z-t-bar: it standardises as if the unit statistics were "
        "independent, so it is oversized (about 10 to 12 %) under a common factor, as Stata's pescadf "
        "(9.9 % in 1000 replications of the same design)",
    ),
    "im_pesaran_shin_common_factor": Study(
        _factor_unit_root(im_pesaran_shin, rho=1.0, lags=1),
        _factor_unit_root(im_pesaran_shin, rho=0.85, lags=1),
        (0.15, 0.40),
        0.7,
        "the design of cips: the test that ignores the common factor is oversized (about 27 %)",
    ),
    "hadri": Study(
        _hadri(share=0.0),
        _hadri(share=1 / 30),
        (0.01, 0.10),
        0.7,
        "30 units x 50 periods; white noise around unit means under the null, one unit out of "
        "30 a random walk under the alternative; the null is stationarity (rejection = a unit "
        "root)",
    ),
    "kao": Study(
        _kao(rho=1.0),
        _kao(rho=0.9),
        (0.01, 0.10),
        0.7,
        "30 units x 40 periods, ADF statistic with one lag; y = x + u with x a random walk, u a "
        "random walk under the null and AR(1) with root 0.9 under the alternative",
    ),
    "westerlund": Study(
        _westerlund(rho=1.0),
        _westerlund(rho=0.9),
        (0.01, 0.10),
        0.5,
        "20 units x 40 periods, variance ratio (some panels); y = x + u with x a random walk, u "
        "a random walk under the null and AR(1) with root 0.9 under the alternative",
    ),
    "pedroni": Study(
        _pedroni(rho=1.0, ar="panelspecific", statistic="adf"),
        _pedroni(rho=0.9, ar="panelspecific", statistic="adf"),
        (0.06, 0.22),
        0.45,
        "as westerlund, group ADF statistic with two kernel lags and one augmented lag; "
        "oversized (about 13 %) at this T, which Pedroni (2004) reports too",
    ),
    "pedroni_panel": Study(
        _pedroni(rho=1.0, ar="same", statistic="adf"),
        _pedroni(rho=0.9, ar="same", statistic="adf"),
        (0.07, 0.24),
        0.55,
        "as pedroni with a common autoregressive parameter (panel ADF statistic); oversized "
        "(about 15 %) at this T",
    ),
    "dumitrescu_hurlin": Study(
        _granger(slope=0.0),
        _granger(slope=0.15),
        (0.02, 0.12),
        0.4,
        "20 units x 30 periods, one lag, Z-tilde; x enters the equation of y with slope 0.15 "
        "under the alternative",
    ),
    "pesaran_cd": Study(
        _cross_section(pesaran_cd, loading=0.0),
        _cross_section(pesaran_cd, loading=0.3),
        (0.01, 0.10),
        0.5,
        "10 units x 30 periods; one common factor, loadings of mean 0.3 under the alternative",
    ),
    "lm_cross_section": Study(
        _cross_section(lm_cross_section, loading=0.0),
        _cross_section(lm_cross_section, loading=0.5),
        (0.01, 0.10),
        0.75,
        "10 units x 30 periods; one common factor, loadings of mean 0.5 under the alternative",
    ),
    "friedman": Study(
        _cross_section(friedman, loading=0.0, shape=(15, 15)),
        _cross_section(friedman, loading=0.5, shape=(15, 15)),
        (0.01, 0.12),
        0.85,
        "15 units x 15 periods (N = T: the p-value of xtcsd uses N - 1 degrees of freedom, "
        "which is right only then); one common factor, loadings of mean 0.5 under the "
        "alternative",
    ),
    "friedman_periods": Study(
        _cross_section(friedman, loading=0.0, shape=(10, 30), dof="periods"),
        _cross_section(friedman, loading=0.5, shape=(10, 30), dof="periods"),
        (0.01, 0.11),
        0.8,
        "10 units x 30 periods with dof='periods' (T - 1 degrees of freedom): the right size "
        "where the default of xtcsd rejects almost always under the null",
    ),
    "frees": Study(
        _cross_section(frees, loading=0.0, frees_rule=True),
        _cross_section(frees, loading=0.5, frees_rule=True),
        (0.01, 0.11),
        0.7,
        "10 units x 30 periods; rejection when the statistic exceeds the 5 % critical value; "
        "one common factor, loadings of mean 0.5 under the alternative",
    ),
    "wooldridge": Study(
        _wooldridge(rho=0.0),
        _wooldridge(rho=0.2),
        (0.01, 0.10),
        0.5,
        "40 units x 10 periods, one regressor and unit effects; errors AR(1) with root 0.2 "
        "under the alternative",
    ),
}

#: The statistic of a test without p-value: (draw under the null, under rho = 0.5, bounds of the
#: mean under the null, upper bound of the mean under the alternative).
DURBIN_WATSON = (_durbin_watson(0.0), _durbin_watson(0.5), (1.85, 2.15), 1.6)


def rejection_rate(draw: Draw, reps: int, seed: int) -> float:
    """Return the share of ``reps`` replications whose p-value is below the level."""
    rng = np.random.default_rng(seed)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # a numerical warning inside a test is a defect
        hits = sum(draw(rng) < LEVEL for _ in range(reps))
    return hits / reps


def mean_statistic(draw: Draw, reps: int, seed: int) -> float:
    """Return the mean of the statistic over ``reps`` replications."""
    rng = np.random.default_rng(seed)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        return float(np.mean([draw(rng) for _ in range(reps)]))


def run_all(reps: int) -> pd.DataFrame:
    """Return the size and power of every study (rows) at ``reps`` replications."""
    rows = {}
    for index, (name, study) in enumerate(STUDIES.items()):
        seed = 1000 + index
        rows[name] = {
            "size": rejection_rate(study.null, reps, seed),
            "power": rejection_rate(study.alternative, reps, seed + 500),
            "scenario": study.note,
        }
    return pd.DataFrame.from_dict(rows, orient="index")


def markdown(table: pd.DataFrame, reps: int, dw: tuple[float, float]) -> str:
    """Render the table of the documentation."""
    lines = [
        "# Size and power of the panel and residual tests",
        "",
        f"Rejection rates at the 5 % level over {reps} replications per cell (seeds fixed in "
        "`validation/simulate.py`; regenerate with `python -m validation.simulate --write`).  "
        "The Monte Carlo standard error of a rate near 5 % is "
        f"{100 * np.sqrt(0.05 * 0.95 / reps):.1f} points.",
        "",
        "| Test | Size (H0 true) | Power (H1 true) | Scenario |",
        "|---|---|---|---|",
    ]
    for name, row in table.iterrows():
        lines.append(f"| `{name}` | {row['size']:.3f} | {row['power']:.3f} | {row['scenario']} |")
    lines += [
        "",
        "`panel_durbin_watson` has no p-value, so its mean statistic is given instead: "
        f"{dw[0]:.3f} without autocorrelation (about 2) and {dw[1]:.3f} with rho = 0.5 "
        "(50 units x 10 periods).",
        "",
        "`wald_groupwise` rejects far too often at this T: this is the behaviour of Stata's "
        "`xttest3`, which the function reproduces (plan item 8.3), not a defect.",
        "",
        "`lm_random_effects` is one-sided: it needs `variance_u` (the random-effects estimate) to "
        "apply the sign rule of `xttest0`; without it a negative within-unit correlation also "
        "rejects and the size is about 10 %.  Pooled residuals alone cannot tell the two cases "
        "apart, which is why the function adds a remark when the correlation is negative.",
        "",
        "`shapiro_francia` rejects slightly too often at n = 200 under the normal (about 6 %): "
        "the normal approximation of Royston (1993) that `sfrancia` uses is coarser than "
        "Shapiro-Wilk's.",
        "",
        "## Panel tests of category A",
        "",
        "Panels are simulated by `validation/simulate.py`; the alternatives are close ones, so "
        "that the power is informative (a rate of 1.000 would only say the scenario is easy).  "
        "A rejection is a rejection of the null of the test: a unit root for `breitung`, "
        "`harris_tzavalis`, `levin_lin_chu`, `im_pesaran_shin` and `fisher`, stationarity for "
        "`hadri`, no cointegration for `kao`, no Granger causality for `dumitrescu_hurlin`, "
        "no cross-section dependence for `pesaran_cd`, `lm_cross_section`, `friedman` and "
        "`frees`, no first-order autocorrelation for `wooldridge`.",
        "",
        "Three unit-root tests are oversized at 30 units x 50 periods, and this is the "
        "behaviour of Stata, not of the package: the same design simulated in Stata 17 "
        "(`xtunitroot`) gives 24.7 % (300 replications) for `llc`, 8.9 % (1000) for `ht` and "
        "7.8 % (1000) for `fisher`, against about 25 %, 9 % and 8 to 9 % here.  `llc` is the "
        "worst: its adjusted t statistic has a mean of about -1 under the null, so that "
        "its rate of rejection at 5 % is five times the level; treat its p-value as "
        "indicative at this T, or prefer `im_pesaran_shin` or `breitung`.",
        "",
        "`friedman` is simulated with as many periods as units.  The statistic is "
        "`(T - 1) ((N - 1) r + 1)`, a chi-square with `T - 1` degrees of freedom under the "
        "null, but `xtcsd, friedman` (and so this function) takes `N - 1`: the p-value is "
        "right only when `N = T`.  With `T > N` it rejects almost always under the null (97 % "
        "at 10 units x 30 periods) and with `T < N` it almost never does (0 % at 30 x 10).  "
        '`dof="periods"` gives the p-value of the distribution (`T - 1` degrees of freedom): '
        "the `friedman_periods` row, 10 units x 30 periods, has the right size.  The default "
        "stays that of Stata, and the result carries a remark when `N` and `T` differ.",
        "",
        "`cips` is simulated with a strong common factor in the innovations, the case it is made "
        "for.  The p-value of its `Z-t-bar` standardises as if the unit statistics were "
        "independent: it is oversized (about 10 to 12 %, 9.9 % in `pescadf` of Stata 17 over 1000 "
        "replications of the same design), whereas the decision of the paper, t-bar against its "
        "critical value, has about the right size (4 to 5 %, the `cips` row; 3.8 % in Stata).  `im_pesaran_shin`, which ignores "
        "the factor, rejects about a quarter of the time (the `im_pesaran_shin_common_factor` "
        "row).  The critical values exist for balanced panels only.",
        "",
        "`pedroni` and `pedroni_panel` are oversized at 20 units x 40 periods (13 to 15 % at "
        "5 %), as Pedroni (2004) finds for his tests at a short T: the statistics are those of "
        "Stata's `xtcointtest pedroni`, to 1e-9 on the validation cases, so the behaviour is "
        "theirs.  The variance-ratio test of `westerlund` has about the right size in the same "
        "design.",
        "",
        "`dumitrescu_hurlin` is slightly oversized (8 %) with 20 units and 30 periods: its "
        "standardisation by the exact moments of the Wald statistics assumes normal errors "
        "and a large number of units.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=int, default=2000)
    parser.add_argument("--write", type=Path, default=None, help="write the Markdown table here")
    args = parser.parse_args(argv)
    table = run_all(args.reps)
    null, alternative, _, _ = DURBIN_WATSON
    dw = (mean_statistic(null, args.reps // 4, 7), mean_statistic(alternative, args.reps // 4, 8))
    text = markdown(table, args.reps, dw)
    if args.write is None:
        print(text)
    else:
        args.write.write_text(text, encoding="utf-8")
        print(f"wrote {args.write}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
