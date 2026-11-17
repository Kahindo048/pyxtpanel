# Size and power of the residual and specification tests

Rejection rates at the 5 % level over 2000 replications per cell (seeds fixed in `validation/simulate.py`; regenerate with `python -m validation.simulate --write`).  The Monte Carlo standard error of a rate near 5 % is 0.5 points.

| Test | Size (H0 true) | Power (H1 true) | Scenario |
|---|---|---|---|
| `breusch_pagan` | 0.045 | 1.000 | n = 300, fitted values; error sd proportional to exp(x1 / 2) under the alternative |
| `white` | 0.054 | 1.000 | n = 300, two regressors; error sd proportional to exp(x1 / 2) under the alternative |
| `wald_groupwise` | 0.461 | 1.000 | 20 units x 30 periods; oversized at this T, as xttest3 (documented in plan 8.3); unit sd from 0.5 to 2 under the alternative |
| `lm_random_effects` | 0.043 | 1.000 | 40 units x 5 periods; unit effect with sd 1 under the alternative; the sign rule of xttest0 applied through variance_u (without it the size is 10 %) |
| `f_effects` | 0.049 | 1.000 | 40 units x 5 periods; unit effect with sd 1 under the alternative |
| `skewness_kurtosis` | 0.057 | 0.921 | n = 200; Student t with 4 degrees of freedom under the alternative |
| `jarque_bera` | 0.051 | 0.944 | n = 200; Student t with 4 degrees of freedom under the alternative |
| `shapiro_wilk` | 0.055 | 0.926 | n = 200; Student t with 4 degrees of freedom under the alternative |
| `shapiro_francia` | 0.065 | 0.947 | n = 200; Student t with 4 degrees of freedom under the alternative |
| `reset` | 0.053 | 0.894 | n = 200, two regressors; an omitted 0.5 * x1^2 under the alternative |
| `slope_homogeneity` | 0.053 | 1.000 | 30 units x 20 periods; slopes with sd 0.6 around 1 under the alternative |

`panel_durbin_watson` has no p-value, so its mean statistic is given instead: 1.996 without autocorrelation (about 2) and 1.218 with rho = 0.5 (50 units x 10 periods).

`wald_groupwise` rejects far too often at this T: this is the behaviour of Stata's `xttest3`, which the function reproduces (plan item 8.3), not a defect.

`lm_random_effects` is one-sided: it needs `variance_u` (the random-effects estimate) to apply the sign rule of `xttest0`; without it a negative within-unit correlation also rejects and the size is about 10 %.  Pooled residuals alone cannot tell the two cases apart, which is why the function adds a remark when the correlation is negative.

`shapiro_francia` rejects slightly too often at n = 200 under the normal (about 6 %): the normal approximation of Royston (1993) that `sfrancia` uses is coarser than Shapiro-Wilk's.
