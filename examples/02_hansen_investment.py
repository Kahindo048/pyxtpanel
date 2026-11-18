"""A real balanced panel: the firm investment data of Hansen (1999), 565 firms, 1973 to 1987.

Reads the Stata file, declares the panel, shows where the variation of each variable lives
(between firms or within a firm), and runs the tests that need no estimator.

The data are in ``examples/data/hansen_1999_invest.dta`` (Hansen, B. E. (1999), Threshold effects in
non-dynamic panels, Journal of Econometrics 93(2); see NOTICE.md).

Run:  python examples/02_hansen_investment.py
"""

from pathlib import Path

import pyxtpanel as xt

DATA = Path(__file__).parent / "data" / "hansen_1999_invest.dta"


def main() -> None:
    frame = xt.read_stata(DATA)
    panel = xt.xtset(frame, "firm", "year")
    print(panel)

    # Where does the variation live?  xtsum splits it as Stata does.
    sums = panel.xtsum(["invest", "tobin_q", "cashflow", "debt"])
    print(sums.summary())
    sd = sums.to_frame()["sd"].unstack("component")
    print("share of the variance that is within a firm (within sd squared / overall sd squared):")
    print(((sd["within"] / sd["overall"]) ** 2).round(2).to_string())

    # The calendar makes lags and differences safe; the panel is balanced, so nothing is lost
    # but the first year of each firm.
    terms = panel.operate("L.invest", "D.cashflow")
    print(
        f"missing after one lag: {int(terms['L.invest'].isna().sum())} of {len(terms)} rows "
        f"(one per firm: {panel.n_units})"
    )

    # Tests that need no estimator, on this balanced panel.
    unit_root = panel.unit_root("invest", test="im_pesaran_shin", lags=1)
    dependence = panel.cross_section_dependence("invest", test="pesaran_cd")
    serial = panel.serial_correlation("invest", ["tobin_q", "cashflow"], test="wooldridge")
    for label, result in (
        ("Im-Pesaran-Shin, unit root in invest", unit_root),
        ("Pesaran CD, cross-sectional dependence", dependence),
        ("Wooldridge, serial correlation in first differences", serial),
    ):
        print(f"{label}: statistic {result.statistic:9.3f}, p-value {result.pvalue:.4f}")

    # A table for a paper, in LaTeX.
    print(sums.to_latex()[:300])


if __name__ == "__main__":
    main()
