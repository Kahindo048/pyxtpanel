"""Simulated and hostile panels shared by Stata and Python (blocks B.2.3 and B.2.4).

Every panel is generated once from a fixed seed, then **frozen** as a CSV under
``tests/golden/data`` so that the golden tests never depend on a random generator.  The Stata
side reads the ``.dta`` version written by :func:`write_dta`.  The time variable is always called
``t``; the identifier ``id``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "tests" / "golden" / "data"
DTA_DIR = Path(__file__).resolve().parent / "data"


@dataclass(frozen=True)
class PanelSpec:
    """How a panel is built and how it must be declared."""

    name: str
    build: Callable[[np.random.Generator], pd.DataFrame]
    time_kind: str = "int"  # "int", "quarter" or "month"
    delta: int = 1
    note: str = ""


SPECS: dict[str, PanelSpec] = {}


def _register(
    name: str, *, time_kind: str = "int", delta: int = 1, note: str = ""
) -> Callable[[Callable[[np.random.Generator], pd.DataFrame]], Callable[..., pd.DataFrame]]:
    def wrap(build: Callable[[np.random.Generator], pd.DataFrame]) -> Callable[..., pd.DataFrame]:
        SPECS[name] = PanelSpec(name, build, time_kind, delta, note)
        return build

    return wrap


# ---------------------------------------------------------------- building blocks


def _cells(units: list[int], ranks_of: Callable[[int], np.ndarray]) -> pd.DataFrame:
    rows = [(u, int(r)) for u in units for r in ranks_of(u)]
    return pd.DataFrame(rows, columns=["id", "rank"])


def _holes(rng: np.random.Generator, n_units: int, n_periods: int, p: float) -> pd.DataFrame:
    """Every unit spans all periods but loses interior cells with probability ``p``."""
    kept = []
    for u in range(1, n_units + 1):
        ranks = np.arange(n_periods)
        drop = rng.random(n_periods) < p
        drop[0] = drop[-1] = False
        kept.append(pd.DataFrame({"id": u, "rank": ranks[~drop]}))
    return pd.concat(kept, ignore_index=True)


def _fill(
    rng: np.random.Generator, cells: pd.DataFrame, time: np.ndarray | pd.Series
) -> pd.DataFrame:
    n = len(cells)
    x = rng.normal(size=n).round(6)
    z = rng.normal(size=n).round(6)
    z[rng.random(n) < 0.1] = np.nan
    out = pd.DataFrame({"id": cells["id"].to_numpy(), "t": time, "x": x, "z": z})
    return out


def _int_time(cells: pd.DataFrame, start: int = 1, step: int = 1) -> np.ndarray:
    return start + step * cells["rank"].to_numpy()


def _quarters(cells: pd.DataFrame) -> pd.Series:
    start = pd.Period("2015Q1", freq="Q")
    return pd.Series([(start + int(r)).start_time for r in cells["rank"]])


def _months(cells: pd.DataFrame) -> pd.Series:
    start = pd.Period("2016-01", freq="M")
    return pd.Series([(start + int(r)).start_time for r in cells["rank"]])


def _shuffle(rng: np.random.Generator, frame: pd.DataFrame) -> pd.DataFrame:
    return frame.iloc[rng.permutation(len(frame))].reset_index(drop=True)


def _iv_panel(rng: np.random.Generator, *, bad_instrument: bool) -> pd.DataFrame:
    """One endogenous regressor, four instruments, errors correlated within a unit (L.6.2).

    ``z4`` is correlated with the error when ``bad_instrument`` is set, so the J tests reject.
    """
    n_units, n_periods = 150, 6
    unit = np.repeat(np.arange(1, n_units + 1), n_periods)
    n = unit.size
    effect = np.repeat(rng.normal(scale=0.7, size=n_units), n_periods)
    u = effect + rng.normal(size=n)
    z = rng.normal(size=(n, 4))
    if bad_instrument:
        z[:, 3] += 0.6 * u
    x = z @ np.array([0.5, 0.4, 0.3, 0.3]) + 0.6 * u + rng.normal(size=n)
    y = 1.0 + 0.8 * x + u
    out = pd.DataFrame({"id": unit, "t": np.tile(np.arange(1, n_periods + 1), n_units)})
    out["y"], out["x"] = y.round(6), x.round(6)
    for j in range(4):
        out[f"z{j + 1}"] = z[:, j].round(6)
    return out


@_register("iv_valid", note="150 units x 6 periods, valid instruments (L.6.2)")
def _iv_valid(rng: np.random.Generator) -> pd.DataFrame:
    return _iv_panel(rng, bad_instrument=False)


@_register("iv_invalid", note="same, one instrument correlated with the error (L.6.2)")
def _iv_invalid(rng: np.random.Generator) -> pd.DataFrame:
    return _iv_panel(rng, bad_instrument=True)


@_register("het", note="300 units x 5 periods, error variance depends on x and z1 (L.8)")
def _het(rng: np.random.Generator) -> pd.DataFrame:
    n_units, n_periods = 300, 5
    n = n_units * n_periods
    x, z1, z2 = rng.normal(size=n), rng.normal(size=n), rng.normal(size=n)
    d = (rng.random(n) < 0.4).astype(float)
    u = rng.normal(size=n) * np.exp(0.5 * x + 0.3 * z1)
    y = 1.0 + x + 0.5 * z1 - 0.3 * z2 + 0.4 * d + u
    out = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n_units + 1), n_periods),
            "t": np.tile(np.arange(1, n_periods + 1), n_units),
        }
    )
    for name, values in (("y", y), ("x", x), ("z1", z1), ("z2", z2), ("d", d)):
        out[name] = values.round(6)
    return out


@_register("groupwise", note="60 units x 8 periods, error variance differs by unit (L.8.3)")
def _groupwise(rng: np.random.Generator) -> pd.DataFrame:
    n_units, n_periods = 60, 8
    n = n_units * n_periods
    sd = np.repeat(np.exp(rng.normal(scale=0.6, size=n_units)), n_periods)
    effect = np.repeat(rng.normal(scale=1.0, size=n_units), n_periods)
    x, z = rng.normal(size=n), rng.normal(size=n)
    y = effect + 0.8 * x - 0.4 * z + sd * rng.normal(size=n)
    out = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n_units + 1), n_periods),
            "t": np.tile(np.arange(1, n_periods + 1), n_units),
        }
    )
    for name, values in (("y", y), ("x", x), ("z", z)):
        out[name] = values.round(6)
    return out


def _factor_panel(
    rng: np.random.Generator, cells: pd.DataFrame, time: np.ndarray | pd.Series
) -> pd.DataFrame:
    """Variables driven by a common factor, for the cross-section dependence tests (L.3.1).

    ``s`` loads on the factor (strong dependence); ``x`` does not; ``y`` has a unit effect, an
    effect of ``x`` and the factor in its error, so the residuals of ``y`` on ``x`` stay dependent.
    """
    n_periods = int(cells["rank"].max()) + 1
    factor = rng.normal(size=n_periods)
    unit = cells["id"].to_numpy()
    loading = rng.normal(loc=1.0, scale=0.5, size=int(unit.max()) + 1)[unit]
    effect = rng.normal(size=int(unit.max()) + 1)[unit]
    common = loading * factor[cells["rank"].to_numpy()]
    n = len(cells)
    x = rng.normal(size=n)
    s = common + rng.normal(size=n)
    y = effect + 0.7 * x + 0.8 * common + rng.normal(size=n)
    out = pd.DataFrame({"id": unit, "t": time})
    for name, values in (("y", y), ("x", x), ("s", s)):
        out[name] = values.round(6)
    return out


@_register("csd", note="25 units x 12 periods balanced, common factor (L.3.1)")
def _csd(rng: np.random.Generator) -> pd.DataFrame:
    cells = _cells(list(range(1, 26)), lambda u: np.arange(12))
    return _factor_panel(rng, cells, _int_time(cells, 1))


@_register("csd_unbalanced", note="30 units, entry, exit and holes, common factor (L.3.1)")
def _csd_unbalanced(rng: np.random.Generator) -> pd.DataFrame:
    starts = rng.integers(0, 3, size=30)
    stops = rng.integers(9, 12, size=30)
    cells = _cells(list(range(1, 31)), lambda u: np.arange(starts[u - 1], stops[u - 1] + 1))
    keep = rng.random(len(cells)) > 0.08
    return _factor_panel(rng, cells[keep].reset_index(drop=True), _int_time(cells[keep], 1))


@_register("csd_trimmed", note="25 units x 12 periods, up to 2 periods cut at each end (L.3)")
def _csd_trimmed(rng: np.random.Generator) -> pd.DataFrame:
    head = rng.integers(0, 3, size=25)
    tail = rng.integers(0, 3, size=25)
    cells = _cells(list(range(1, 26)), lambda u: np.arange(head[u - 1], 12 - tail[u - 1]))
    return _factor_panel(rng, cells, _int_time(cells, 1))


@_register("csd_long", note="15 units x 40 periods balanced, common factor (L.3)")
def _csd_long(rng: np.random.Generator) -> pd.DataFrame:
    cells = _cells(list(range(1, 16)), lambda u: np.arange(40))
    return _factor_panel(rng, cells, _int_time(cells, 1))


@_register("serial", note="80 units x 8 periods, AR(1) errors with rho 0.6 (L.5.1)")
def _serial(rng: np.random.Generator) -> pd.DataFrame:
    n_units, n_periods = 80, 8
    u = np.empty((n_units, n_periods))
    u[:, 0] = rng.normal(size=n_units)
    for t in range(1, n_periods):
        u[:, t] = 0.6 * u[:, t - 1] + 0.8 * rng.normal(size=n_units)
    n = n_units * n_periods
    effect = np.repeat(rng.normal(size=n_units), n_periods)
    x, z = rng.normal(size=n), rng.normal(size=n)
    y = effect + 0.8 * x - 0.4 * z + u.ravel()
    out = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n_units + 1), n_periods),
            "t": np.tile(np.arange(1, n_periods + 1), n_units),
        }
    )
    for name, values in (("y", y), ("x", x), ("z", z)):
        out[name] = values.round(6)
    return out


def _slopes_panel(
    rng: np.random.Generator,
    cells: pd.DataFrame,
    time: np.ndarray | pd.Series,
    *,
    heterogeneous: bool,
) -> pd.DataFrame:
    """Panel for the slope homogeneity test (L.8.11): ``y`` on ``x`` and ``z``, plus ``w``.

    The slopes of ``x`` and ``z`` are the same in every unit, or drawn per unit when
    ``heterogeneous``; ``w`` enters with a common slope and is meant to be partialled out; a
    common factor drives ``x`` and the error of ``y`` so that cross-sectional averages matter.
    """
    n_periods = int(cells["rank"].max()) + 1
    factor = rng.normal(size=n_periods)[cells["rank"].to_numpy()]
    unit = cells["id"].to_numpy()
    n_units = int(unit.max()) + 1
    spread = 1.0 if heterogeneous else 0.0
    slope_x = (1.0 + 0.5 * spread * rng.normal(size=n_units))[unit]
    slope_z = (0.5 + 0.3 * spread * rng.normal(size=n_units))[unit]
    effect = rng.normal(size=n_units)[unit]
    n = len(cells)
    x = 0.6 * factor + rng.normal(size=n)
    z, w = rng.normal(size=n), rng.normal(size=n)
    y = effect + slope_x * x + slope_z * z + 0.4 * w + 0.7 * factor + rng.normal(size=n)
    out = pd.DataFrame({"id": unit, "t": time})
    for name, values in (("y", y), ("x", x), ("z", z), ("w", w)):
        out[name] = values.round(6)
    return out


@_register("slopes_homog", note="40 units x 20 periods, the same slopes in every unit (L.8.11)")
def _slopes_homog(rng: np.random.Generator) -> pd.DataFrame:
    cells = _cells(list(range(1, 41)), lambda u: np.arange(20))
    return _slopes_panel(rng, cells, _int_time(cells, 1), heterogeneous=False)


@_register("slopes_het", note="40 units x 20 periods, slopes drawn per unit (L.8.11)")
def _slopes_het(rng: np.random.Generator) -> pd.DataFrame:
    cells = _cells(list(range(1, 41)), lambda u: np.arange(20))
    return _slopes_panel(rng, cells, _int_time(cells, 1), heterogeneous=True)


@_register("slopes_unbalanced", note="40 units, entry, exit and holes, per-unit slopes (L.8.11)")
def _slopes_unbalanced(rng: np.random.Generator) -> pd.DataFrame:
    starts = rng.integers(0, 4, size=40)
    stops = rng.integers(15, 20, size=40)
    cells = _cells(list(range(1, 41)), lambda u: np.arange(starts[u - 1], stops[u - 1] + 1))
    keep = rng.random(len(cells)) > 0.06
    return _slopes_panel(
        rng, cells[keep].reset_index(drop=True), _int_time(cells[keep], 1), heterogeneous=True
    )


def _unit_root_panel(
    rng: np.random.Generator, *, n_units: int, n_periods: int, kind: str
) -> pd.DataFrame:
    """Series of ``kind`` ("walk", "stationary", "trend", "mixed") with a common innovation shock."""
    shock = 0.8 * rng.normal(size=n_periods)
    level = rng.normal(scale=2.0, size=(n_units, 1))
    e = rng.normal(size=(n_units, n_periods)) + shock  # the shock enters every unit's innovation
    y = np.empty((n_units, n_periods))
    walk = np.ones(n_units, dtype=bool)
    ar = np.zeros(n_units)
    slope = np.zeros(n_units)
    if kind == "stationary":
        walk[:], ar[:] = False, 0.5
    elif kind == "trend":
        walk[:], ar[:] = False, 0.4
        slope = rng.normal(loc=0.3, scale=0.1, size=n_units)
    elif kind == "mixed":
        walk = np.arange(n_units) % 2 == 0
        ar = np.where(walk, 1.0, 0.6)
    else:
        ar[:] = 1.0
    state = e[:, 0]
    y[:, 0] = state
    for t in range(1, n_periods):
        state = ar * state + e[:, t]
        y[:, t] = state
    y = y + level + slope[:, None] * np.arange(n_periods)
    return pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n_units + 1), n_periods),
            "t": np.tile(np.arange(1, n_periods + 1), n_units),
            "y": y.ravel().round(6),
        }
    )


@_register("ur_walk", note="50 units x 12 periods, random walks with a period shock (L.2)")
def _ur_walk(rng: np.random.Generator) -> pd.DataFrame:
    return _unit_root_panel(rng, n_units=50, n_periods=12, kind="walk")


@_register("ur_stationary", note="50 units x 12 periods, AR(0.5) with a period shock (L.2)")
def _ur_stationary(rng: np.random.Generator) -> pd.DataFrame:
    return _unit_root_panel(rng, n_units=50, n_periods=12, kind="stationary")


@_register("ur_trend", note="50 units x 12 periods, trend plus AR(0.4) (L.2)")
def _ur_trend(rng: np.random.Generator) -> pd.DataFrame:
    return _unit_root_panel(rng, n_units=50, n_periods=12, kind="trend")


@_register("ur_mixed", note="40 units x 20 periods, half random walks, half AR(0.6) (L.2)")
def _ur_mixed(rng: np.random.Generator) -> pd.DataFrame:
    return _unit_root_panel(rng, n_units=40, n_periods=20, kind="mixed")


@_register("ur_long", note="30 units x 60 periods, half random walks, half AR(0.6) (L.2)")
def _ur_long(rng: np.random.Generator) -> pd.DataFrame:
    return _unit_root_panel(rng, n_units=30, n_periods=60, kind="mixed")


@_register("ur_unbalanced", note="40 units, 12 to 24 consecutive periods, no holes (L.2)")
def _ur_unbalanced(rng: np.random.Generator) -> pd.DataFrame:
    full = _unit_root_panel(rng, n_units=40, n_periods=24, kind="mixed")
    starts = rng.integers(0, 7, size=40)
    stops = rng.integers(18, 24, size=40)
    keep = (full["t"].to_numpy() - 1 >= starts[full["id"].to_numpy() - 1]) & (
        full["t"].to_numpy() - 1 <= stops[full["id"].to_numpy() - 1]
    )
    return full[keep].reset_index(drop=True)


# ---------------------------------------------------------------- simulated panels (B.2.3)


@_register("balanced", note="40 units x 8 periods")
def _balanced(rng: np.random.Generator) -> pd.DataFrame:
    cells = _cells(list(range(1, 41)), lambda u: np.arange(8))
    return _fill(rng, cells, _int_time(cells, 1990))


@_register("gaps", note="interior holes")
def _gaps(rng: np.random.Generator) -> pd.DataFrame:
    cells = _holes(rng, 40, 10, 0.2)
    return _fill(rng, cells, _int_time(cells, 2000))


@_register("entry_exit", note="late entry and early exit, no holes")
def _entry_exit(rng: np.random.Generator) -> pd.DataFrame:
    starts = rng.integers(0, 4, size=40)
    stops = rng.integers(6, 10, size=40)
    cells = _cells(list(range(1, 41)), lambda u: np.arange(starts[u - 1], stops[u - 1] + 1))
    return _fill(rng, cells, _int_time(cells, 1))


@_register("weak", note="same number of periods per unit, different windows")
def _weak(rng: np.random.Generator) -> pd.DataFrame:
    starts = rng.integers(0, 4, size=30)
    cells = _cells(list(range(1, 31)), lambda u: np.arange(starts[u - 1], starts[u - 1] + 5))
    return _fill(rng, cells, _int_time(cells, 1))


@_register("delta2", delta=2, note="time step 2 with holes")
def _delta2(rng: np.random.Generator) -> pd.DataFrame:
    cells = _holes(rng, 30, 8, 0.25)
    return _fill(rng, cells, _int_time(cells, 1, 2))


@_register("quarterly", time_kind="quarter", note="quarters with holes")
def _quarterly(rng: np.random.Generator) -> pd.DataFrame:
    cells = _holes(rng, 30, 12, 0.2)
    return _fill(rng, cells, _quarters(cells))


@_register("monthly", time_kind="month", note="months with holes")
def _monthly(rng: np.random.Generator) -> pd.DataFrame:
    cells = _holes(rng, 25, 24, 0.15)
    return _fill(rng, cells, _months(cells))


@_register("unsorted_ids", note="arbitrary identifiers, shuffled rows, holes")
def _unsorted_ids(rng: np.random.Generator) -> pd.DataFrame:
    cells = _holes(rng, 35, 9, 0.2)
    labels = rng.choice(np.arange(1, 5000), size=35, replace=False)
    cells["id"] = labels[cells["id"].to_numpy() - 1]
    return _shuffle(rng, _fill(rng, cells, _int_time(cells, 1980)))


@_register("constants", note="a variable constant within each unit and one constant everywhere")
def _constants(rng: np.random.Generator) -> pd.DataFrame:
    frame = _balanced(rng)
    level = rng.normal(size=40).round(6)
    frame["cu"] = level[frame["id"].to_numpy() - 1]
    frame["ck"] = 3.5
    return frame


@_register("re_effects", note="a unit effect in x, unbalanced with holes (random-effects theta)")
def _re_effects(rng: np.random.Generator) -> pd.DataFrame:
    cells = _holes(rng, 40, 8, 0.2)
    frame = _fill(rng, cells, _int_time(cells, 2000))
    effect = rng.normal(size=40).round(6)
    z = rng.normal(size=len(frame)).round(6)
    noise = rng.normal(scale=0.8, size=len(frame))
    frame["z"] = z
    frame["x"] = (0.5 * z + effect[frame["id"].to_numpy() - 1] + noise).round(6)
    return frame


@_register("wide", note="130 periods, scattered absences: more columns than xtdescribe's width")
def _wide(rng: np.random.Generator) -> pd.DataFrame:
    cells = _holes(rng, 30, 130, 0.1)
    return _fill(rng, cells, _int_time(cells, 1900))


@_register("sparse_t", note="observed periods 2000, 2002, 2004, 2006 only (smallest step is 2)")
def _sparse_t(rng: np.random.Generator) -> pd.DataFrame:
    cells = _holes(rng, 25, 4, 0.3)
    return _fill(rng, cells, _int_time(cells, 2000, 2))


# ---------------------------------------------------------------- hostile panels (B.2.4)


@_register("one_obs_unit", note="a unit with a single observation")
def _one_obs_unit(rng: np.random.Generator) -> pd.DataFrame:
    frame = _gaps(rng)
    extra = pd.DataFrame({"id": [999], "t": [2003], "x": [0.5], "z": [np.nan]})
    return pd.concat([frame, extra], ignore_index=True)


@_register("single_unit", note="one unit only")
def _single_unit(rng: np.random.Generator) -> pd.DataFrame:
    cells = _cells([1], lambda u: np.arange(6))
    return _fill(rng, cells, _int_time(cells, 1))


@_register("single_period", note="one period only")
def _single_period(rng: np.random.Generator) -> pd.DataFrame:
    cells = _cells(list(range(1, 21)), lambda u: np.arange(1))
    return _fill(rng, cells, _int_time(cells, 2010))


@_register("all_missing_var", note="a variable that is missing everywhere")
def _all_missing_var(rng: np.random.Generator) -> pd.DataFrame:
    frame = _balanced(rng)
    frame["z"] = np.nan
    return frame


@_register("missing_keys", note="rows with a missing identifier or time (Stata ignores them)")
def _missing_keys(rng: np.random.Generator) -> pd.DataFrame:
    frame = _gaps(rng)
    frame["id"] = frame["id"].astype(float)
    frame["t"] = frame["t"].astype(float)
    frame.loc[[3, 50, 120], "id"] = np.nan
    frame.loc[[7, 80], "t"] = np.nan
    return frame


#: Refused by ``xtset`` in Stata and in the package: the case checks the return code.
HOSTILE_REFUSED: dict[str, pd.DataFrame] = {
    "duplicates": pd.DataFrame({"id": [1, 1, 1, 2], "t": [1, 2, 2, 1], "x": [1.0, 2.0, 3.0, 4.0]}),
    "non_integer_time": pd.DataFrame({"id": [1, 1, 2], "t": [1.0, 2.5, 1.0], "x": [1.0, 2.0, 3.0]}),
}


# ---------------------------------------------------------------- freezing and loading


def generate(name: str) -> pd.DataFrame:
    """Build a panel from its seed (use :func:`load` in tests)."""
    seed = sum(ord(c) for c in name) * 7919
    return SPECS[name].build(np.random.default_rng(seed))


def freeze_panels() -> None:
    """Write every panel to ``tests/golden/data/<name>.csv`` (ISO dates, 17 digits)."""
    FROZEN.mkdir(parents=True, exist_ok=True)
    for name in SPECS:
        generate(name).to_csv(FROZEN / f"{name}.csv", index=False, float_format="%.17g")


def load(name: str) -> pd.DataFrame:
    """Read a frozen panel; time columns are parsed according to the panel's kind."""
    frame = pd.read_csv(FROZEN / f"{name}.csv")
    if SPECS[name].time_kind != "int":
        frame["t"] = pd.to_datetime(frame["t"])
    return frame


def write_dta(name: str) -> Path:
    """Write the frozen panel as a Stata file ``validation/data/sim_<name>.dta``."""
    DTA_DIR.mkdir(parents=True, exist_ok=True)
    frame = load(name)
    kind = SPECS[name].time_kind
    dates = {"t": {"quarter": "tq", "month": "tm"}[kind]} if kind != "int" else None
    path = DTA_DIR / f"sim_{name}.dta"
    frame.to_stata(path, write_index=False, convert_dates=dates, version=118)
    return path
