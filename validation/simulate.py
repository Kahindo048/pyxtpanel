"""Size and power of the residual and specification tests by simulation (plan L.8.12).

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
    breusch_pagan,
    f_effects,
    jarque_bera,
    lm_random_effects,
    panel_durbin_watson,
    reset,
    shapiro_francia,
    shapiro_wilk,
    skewness_kurtosis,
    slope_homogeneity,
    wald_groupwise,
    white,
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
        "# Size and power of the residual and specification tests",
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
