"""A country-by-year panel, from the data frame to a first reading.

Declares a panel, describes it, takes a lag and a first difference on the calendar, removes the unit
means (within), audits the data and runs one test.  The data are simulated: 24 countries, 1990 to
2014, 5 % of the rows missing at random.

Run:  python examples/01_country_year.py
"""

import numpy as np
import pandas as pd

import pyxtpanel as xt


def simulate(seed: int = 2026) -> pd.DataFrame:
    """A long frame (one row per country and year) with holes, as real data have."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(24):
        level = rng.normal(8, 1)
        for year in range(1990, 2015):
            level += 0.02 + 0.03 * rng.normal()
            rows.append((f"C{i:02d}", year, level))
    frame = pd.DataFrame(rows, columns=["country", "year", "lgdp"])
    frame["invest"] = 0.3 * frame["lgdp"] + rng.normal(size=len(frame))
    return frame[rng.random(len(frame)) > 0.05].reset_index(drop=True)


def main() -> None:
    df = simulate()

    # 1. Declare the panel once (like `xtset country year`); the frame is never modified.
    panel = xt.xtset(df, "country", "year")
    print(panel)
    print(
        f"{panel.n_units} countries, {panel.n_obs} rows, {panel.balanced}, holes: {panel.n_gaps}\n"
    )

    # 2. Describe it (like `xtdescribe`, `xtsum`).
    print(panel.xtdescribe().summary())
    print(panel.xtsum(["lgdp", "invest"]).summary())

    # 3. Operators work on the calendar: the lag of a year after a hole is missing, not the
    #    value of two years before.
    terms = panel.operate("L.lgdp", "D.lgdp")
    print(terms.head(4))
    within = panel.within("invest")  # invest minus the mean of its country
    print(f"within variance share: {within.var() / df['invest'].var():.2f}\n")

    # 4. The audit lists what to know before modelling, in English or French.
    report = panel.audit()
    print(report.summary())

    # 5. A statistical test straight from the panel, with no estimator: Pesaran's test of
    #    cross-sectional dependence accepts holes (each pair of countries uses its common years).
    result = panel.cross_section_dependence("invest", test="pesaran_cd")
    print(f"Pesaran CD: statistic {result.statistic:.3f}, p-value {result.pvalue:.3f}")

    # 6. Every result converts to a table, Markdown or LaTeX.
    print(panel.xtsum("lgdp").to_markdown())


if __name__ == "__main__":
    main()
