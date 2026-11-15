"""Render the gallery of figures from a reproducible simulated panel (F.3.3).

    python -m tools.galerie [output_directory]     # default: docs/gallery

The panel is generated from a fixed seed (countries x years with entries, exits, holes and a
unit effect), so the same command gives the same figures.  Needs the ``plot`` extra.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def gallery_panel(seed: int = 2024) -> pd.DataFrame:
    """A panel of 30 'countries' over 1995-2014 with a unit effect, holes and late entries."""
    rng = np.random.default_rng(seed)
    rows = []
    for unit in range(1, 31):
        start = int(rng.integers(0, 4))
        stop = 20 - int(rng.integers(0, 4))
        effect = rng.normal(scale=1.2)
        level = effect
        for step in range(start, stop):
            level = 0.7 * level + effect * 0.3 + rng.normal(scale=0.5)
            if 0 < step < stop - 1 and rng.random() < 0.08:
                continue
            gdp = effect + 0.4 * level + rng.normal(scale=0.4)
            rows.append((f"C{unit:02d}", 1995 + step, level, gdp))
    frame = pd.DataFrame(rows, columns=["country", "year", "growth", "income"])
    frame["regime"] = np.select(
        [frame["growth"] < -0.5, frame["growth"] < 0.5], [0.0, 1.0], default=2.0
    )
    # a staggered reform: 60 % of the countries adopt it in a year between 2000 and 2008
    # (a separate generator, so the other columns do not change)
    later = np.random.default_rng(seed + 1)
    adopts = {
        f"C{unit:02d}": int(later.integers(2000, 2009)) if later.random() < 0.6 else 9999
        for unit in range(1, 31)
    }
    frame["reform"] = (frame["year"] >= frame["country"].map(adopts)).astype(float)
    return frame


def render(directory: Path | None = None) -> list[Path]:
    """Write the figures of the gallery (PNG) and return their paths."""
    import pyxtpanel as xt
    from pyxtpanel.plot import save_figure

    out = Path(directory) if directory else ROOT / "docs" / "gallery"
    out.mkdir(parents=True, exist_ok=True)
    panel = xt.xtset(gallery_panel(), "country", "year")
    figures = {
        "lines": panel.plot("lines", variable="growth", n_units=14, highlight=["C03", "C17"]),
        "lines_panels": panel.plot("lines", variable="growth", n_units=8, panels=True),
        "scatter": panel.plot("scatter", x="income", y="growth"),
        "distribution": panel.plot("distribution", variable="growth"),
        "participation": panel.plot("participation"),
        "transition": panel.plot("transition", variable="regime"),
        "transition_fr": panel.plot("transition", variable="regime", lang="fr"),
        "states": panel.plot("states", variable="reform"),
        "attrition": panel.plot("attrition"),
        "balance_cost": panel.plot("balance_cost"),
    }
    return [path for name, fig in figures.items() for path in save_figure(fig, out / f"{name}.png")]


if __name__ == "__main__":
    for path in render(Path(sys.argv[1]) if len(sys.argv) > 1 else None):
        print("wrote", path)
