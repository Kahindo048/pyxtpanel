"""An unbalanced panel: why it is unbalanced, what balancing costs, what a GMM estimator keeps.

Real panels lose units and periods.  This example builds one with the three causes mixed (late
entry, early exit, internal holes), then asks the exploration tools to say which cause hits which
unit, what balancing would cost, and how many equations a dynamic-panel estimator has left.  The
data are simulated, so the answers can be checked against the recipe below.

Run:  python examples/03_unbalanced_panel.py
"""

import numpy as np
import pandas as pd

import pyxtpanel as xt


def simulate(seed: int = 7) -> pd.DataFrame:
    """40 units over 2000 to 2019: 10 enter late, 8 exit early, 10 have holes, 12 are complete."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(40):
        first, last, holes = 2000, 2019, ()
        if i < 10:
            first = 2000 + int(rng.integers(2, 7))  # late entry
        elif i < 18:
            last = 2019 - int(rng.integers(2, 7))  # early exit
        elif i < 28:
            holes = tuple(rng.choice(np.arange(2003, 2017), size=2, replace=False))  # holes
        for year in range(first, last + 1):
            if year in holes:
                continue
            rows.append((f"U{i:02d}", year, rng.normal() + 0.1 * (year - 2000)))
    return pd.DataFrame(rows, columns=["unit", "year", "y"])


def main() -> None:
    panel = xt.xtset(simulate(), "unit", "year")
    print(panel)

    # 1. Which cause hits which unit?
    gaps = panel.gap_types()
    print(gaps.counts)  # complete, late entry, early exit, internal gap, re-entry

    # 2. How unbalanced is it, in one number? (Ahrens-Pincus gamma: 1 is balanced)
    print(f"gamma = {panel.unbalancedness().gamma:.3f}")

    # 3. What would balancing cost?  Units and rows kept by each minimum length T_i >= k.
    print(panel.balance_cost().to_frame().head(8))

    # 4. Attrition: entries, exits and the survival of units.
    print(panel.attrition().to_frame().head(5))

    # 5. What does a dynamic-panel estimator keep with one lag?  First differences lose an
    #    equation after every hole; forward orthogonal deviations do not.
    for transform in ("fd", "fod"):
        eq = panel.usable_equations("y", lags=1, transform=transform)
        print(f"{transform}: {eq.n_usable} usable equations of {eq.n_observed} observed rows")

    # 6. Fill the calendar (like `tsfill`): the holes become rows with missing values.
    filled = panel.tsfill()
    print(f"after tsfill: {filled.n_obs} rows (was {panel.n_obs}), holes: {filled.n_gaps}")


if __name__ == "__main__":
    main()
