"""Working with an estimator: prepare the data for it, then test its residuals.

``pyxtpanel`` estimates nothing.  An estimation package (or your own code) asks it for the arrays
and the calendar, and sends back the residuals for the panel tests that every serious regression
should pass.  Here a fixed-effects regression is written by hand in numpy, on the Hansen (1999)
firm investment data, to show both directions.

Run:  python examples/04_after_an_estimator.py
"""

from pathlib import Path

import numpy as np

import pyxtpanel as xt

DATA = Path(__file__).parent / "data" / "hansen_1999_invest.dta"


def main() -> None:
    panel = xt.xtset(xt.read_stata(DATA), "firm", "year")

    # 1. What a dynamic-panel estimator needs: the lag and the forward orthogonal deviations
    #    (they keep every equation of a unit except its last), both computed on the calendar.
    lagged = panel.lag("invest", 1)
    fod = panel.fod("invest")
    usable = panel.usable_equations("invest", lags=1, transform="fod")
    print(
        f"rows: {panel.n_obs}, with a lag: {int(lagged.notna().sum())}, "
        f"with a FOD value: {int(fod.notna().sum())}, FOD equations usable: {usable.n_usable}"
    )

    # 2. The array contract: plain float arrays, one slice of rows per unit, the calendar
    #    position of each row.  This is all an estimator may rely on.
    arrays = panel.arrays(["invest", "tobin_q", "cashflow"])
    print(
        f"values {arrays.values.shape}, offsets {arrays.offsets.shape}, "
        f"contract version {arrays.contract_version}"
    )
    print("first firm owns rows", arrays.unit_slices()[0])

    # 3. A fixed-effects regression by hand: remove the firm means, then least squares.
    within = panel.within(["invest", "tobin_q", "cashflow"]).to_numpy()
    y, x = within[:, 0], within[:, 1:]
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residuals = y - x @ beta
    print(f"within coefficients: tobin_q {beta[0]:.4f}, cashflow {beta[1]:.4f}")

    # 4. Category B tests take bare arrays, from any estimator: the residuals, the regressors
    #    and the offsets that delimit each firm.  The test of the firm effects and the
    #    heteroskedasticity test read the residuals of the pooled regression; the others, those of
    #    the fixed-effects regression above.
    offsets, regressors = arrays.offsets, arrays.values[:, 1:]
    design = np.column_stack([np.ones(len(regressors)), regressors])
    pooled = (
        arrays.values[:, 0] - design @ np.linalg.lstsq(design, arrays.values[:, 0], rcond=None)[0]
    )
    checks = {
        "lm_random_effects (variance of the firm effects)": xt.run_test(
            "lm_random_effects", pooled, offsets
        ),
        "breusch_pagan (heteroskedasticity)": xt.run_test("breusch_pagan", pooled, regressors),
        "wald_groupwise (firm-by-firm variances)": xt.run_test(
            "wald_groupwise", residuals, offsets
        ),
        "jarque_bera (normality of the residuals)": xt.run_test("jarque_bera", residuals),
    }
    for label, result in checks.items():
        print(f"{label}: statistic {result.statistic:12.2f}, p-value {result.pvalue:.4f}")


if __name__ == "__main__":
    main()
