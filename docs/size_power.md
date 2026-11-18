# Size and power of the panel and residual tests

Rejection rates at the 5 % level over 2000 replications per cell (seeds fixed in `validation/simulate.py`; regenerate with `python -m validation.simulate --write`).  The Monte Carlo standard error of a rate near 5 % is 0.5 points.

| Test | Size (H0 true) | Power (H1 true) | Scenario |
|---|---|---|---|
| `breusch_pagan` | 0.045 | 1.000 | n = 300, fitted values; error sd proportional to exp(x1 / 2) under the alternative |
| `white` | 0.054 | 1.000 | n = 300, two regressors; error sd proportional to exp(x1 / 2) under the alternative |
| `wald_groupwise` | 0.461 | 1.000 | 20 units x 30 periods; oversized at this T, as xttest3 (documented in plan 8.3); unit sd from 0.5 to 2 under the alternative |
| `lm_random_effects` | 0.043 | 1.000 | 40 units x 5 periods; unit effect with sd 1 under the alternative; the sign rule of xttest0 applied through variance_u (without it the size is 10 %) |
| `baltagi_li` | 0.070 | 1.000 | 40 units x 8 periods; independent errors under the null, a unit effect with sd 1 under the alternative; the joint statistic (chi-square with 2 degrees of freedom) |
| `baltagi_li_serial_adjusted` | 0.070 | 0.992 | 40 units x 8 periods with a unit effect of sd 0.5 in both: AR(1) errors with root 0.3 under the alternative; robust to the effect only when it is small (local), the size is about 21 % with an effect of sd 1 |
| `baltagi_li_re_adjusted` | 0.059 | 1.000 | 40 units x 8 periods with AR(1) errors of root 0.1 in both: a unit effect of sd 1 under the alternative; one-sided robust statistic, robust to serial correlation only when it is small (local), the size is about 20 % with root 0.3 and 67 % with root 0.5 |
| `f_effects` | 0.048 | 1.000 | 40 units x 5 periods; unit effect with sd 1 under the alternative |
| `skewness_kurtosis` | 0.059 | 0.927 | n = 200; Student t with 4 degrees of freedom under the alternative |
| `jarque_bera` | 0.042 | 0.944 | n = 200; Student t with 4 degrees of freedom under the alternative |
| `shapiro_wilk` | 0.044 | 0.920 | n = 200; Student t with 4 degrees of freedom under the alternative |
| `shapiro_francia` | 0.048 | 0.950 | n = 200; Student t with 4 degrees of freedom under the alternative |
| `reset` | 0.051 | 0.895 | n = 200, two regressors; an omitted 0.5 * x1^2 under the alternative |
| `slope_homogeneity` | 0.046 | 1.000 | 30 units x 20 periods; slopes with sd 0.6 around 1 under the alternative |
| `slope_homogeneity_hac_bartlett` | 0.044 | 1.000 | HAC version, Bartlett kernel with the automatic bandwidth; 40 units x 30 periods, unit-specific error variances and AR(1) errors with root 0.4; slopes with sd 0.6 around 1 under the alternative |
| `slope_homogeneity_hac_qs` | 0.021 | 1.000 | as the Bartlett version with the quadratic spectral kernel |
| `breitung` | 0.054 | 0.992 | 30 units x 50 periods; random walks under the null, AR(1) with root 0.95 under the alternative; unit means |
| `harris_tzavalis` | 0.090 | 0.967 | 30 units x 50 periods, root 0.95 under the alternative; oversized (about 9 %) at this T, as Stata's xtunitroot ht: the moments of the statistic are those of a fixed T |
| `levin_lin_chu` | 0.269 | 0.529 | 30 units x 50 periods, one ADF lag, root 0.95 under the alternative; oversized (about 25 %) with the default kernel, as Stata's xtunitroot llc (24.7 % in 300 replications of the same design) |
| `im_pesaran_shin` | 0.044 | 0.583 | 30 units x 50 periods, one ADF lag, root 0.95 under the alternative |
| `fisher` | 0.080 | 0.434 | 30 units x 50 periods, Dickey-Fuller with one lag, inverse chi-square combination, root 0.95 under the alternative; a little oversized (about 8 %) as Stata's xtunitroot fisher (7.8 % in 1000 replications of the same design) |
| `fisher_phillips_perron` | 0.107 | 0.570 | as fisher with the Phillips-Perron unit test (two Newey-West lags) |
| `cips` | 0.042 | 0.915 | 30 units x 50 periods, innovations with a strong common factor (loadings of mean 1), one lag; rejection when t-bar is below its 5 % critical value (the paper's rule); root 0.85 under the alternative |
| `cips_z` | 0.100 | 0.983 | as cips with the p-value of Z-t-bar: it standardises as if the unit statistics were independent, so it is oversized (about 10 to 12 %) under a common factor, as Stata's pescadf (9.9 % in 1000 replications of the same design) |
| `im_pesaran_shin_common_factor` | 0.237 | 0.948 | the design of cips: the test that ignores the common factor is oversized (about 27 %) |
| `hadri` | 0.054 | 0.905 | 30 units x 50 periods; white noise around unit means under the null, one unit out of 30 a random walk under the alternative; the null is stationarity (rejection = a unit root) |
| `kao` | 0.051 | 0.856 | 30 units x 40 periods, ADF statistic with one lag; y = x + u with x a random walk, u a random walk under the null and AR(1) with root 0.9 under the alternative |
| `westerlund` | 0.039 | 0.640 | 20 units x 40 periods, variance ratio (some panels); y = x + u with x a random walk, u a random walk under the null and AR(1) with root 0.9 under the alternative |
| `pedroni` | 0.116 | 0.605 | as westerlund, group ADF statistic with two kernel lags and one augmented lag; oversized (about 13 %) at this T, which Pedroni (2004) reports too |
| `pedroni_panel` | 0.153 | 0.769 | as pedroni with a common autoregressive parameter (panel ADF statistic); oversized (about 15 %) at this T |
| `dumitrescu_hurlin` | 0.068 | 0.622 | 20 units x 30 periods, one lag, Z-tilde; x enters the equation of y with slope 0.15 under the alternative |
| `pesaran_cd` | 0.053 | 0.673 | 10 units x 30 periods; one common factor, loadings of mean 0.3 under the alternative |
| `lm_cross_section` | 0.064 | 0.864 | 10 units x 30 periods; one common factor, loadings of mean 0.5 under the alternative |
| `friedman` | 0.044 | 0.949 | 15 units x 15 periods (N = T: the p-value of xtcsd uses N - 1 degrees of freedom, which is right only then); one common factor, loadings of mean 0.5 under the alternative |
| `friedman_periods` | 0.035 | 0.965 | 10 units x 30 periods with dof='periods' (T - 1 degrees of freedom): the right size where the default of xtcsd rejects almost always under the null |
| `frees` | 0.067 | 0.851 | 10 units x 30 periods; rejection when the statistic exceeds the 5 % critical value; one common factor, loadings of mean 0.5 under the alternative |
| `wooldridge` | 0.053 | 0.635 | 40 units x 10 periods, one regressor and unit effects; errors AR(1) with root 0.2 under the alternative |

`panel_durbin_watson` has no p-value, so its mean statistic is given instead: 1.996 without autocorrelation (about 2) and 1.218 with rho = 0.5 (50 units x 10 periods).

`wald_groupwise` rejects far too often at this T: this is the behaviour of Stata's `xttest3`, which the function reproduces (plan item 8.3), not a defect.

`lm_random_effects` is one-sided: it needs `variance_u` (the random-effects estimate) to apply the sign rule of `xttest0`; without it a negative within-unit correlation also rejects and the size is about 10 %.  Pooled residuals alone cannot tell the two cases apart, which is why the function adds a remark when the correlation is negative.

`shapiro_francia` rejects slightly too often at n = 200 under the normal (about 6 %): the normal approximation of Royston (1993) that `sfrancia` uses is coarser than Shapiro-Wilk's.

## Panel tests of category A

Panels are simulated by `validation/simulate.py`; the alternatives are close ones, so that the power is informative (a rate of 1.000 would only say the scenario is easy).  A rejection is a rejection of the null of the test: a unit root for `breitung`, `harris_tzavalis`, `levin_lin_chu`, `im_pesaran_shin` and `fisher`, stationarity for `hadri`, no cointegration for `kao`, no Granger causality for `dumitrescu_hurlin`, no cross-section dependence for `pesaran_cd`, `lm_cross_section`, `friedman` and `frees`, no first-order autocorrelation for `wooldridge`.

Three unit-root tests are oversized at 30 units x 50 periods, and this is the behaviour of Stata, not of the package: the same design simulated in Stata 17 (`xtunitroot`) gives 24.7 % (300 replications) for `llc`, 8.9 % (1000) for `ht` and 7.8 % (1000) for `fisher`, against about 25 %, 9 % and 8 to 9 % here.  `llc` is the worst: its adjusted t statistic has a mean of about -1 under the null, so that its rate of rejection at 5 % is five times the level; treat its p-value as indicative at this T, or prefer `im_pesaran_shin` or `breitung`.

`friedman` is simulated with as many periods as units.  The statistic is `(T - 1) ((N - 1) r + 1)`, a chi-square with `T - 1` degrees of freedom under the null, but `xtcsd, friedman` (and so this function) takes `N - 1`: the p-value is right only when `N = T`.  With `T > N` it rejects almost always under the null (97 % at 10 units x 30 periods) and with `T < N` it almost never does (0 % at 30 x 10).  `dof="periods"` gives the p-value of the distribution (`T - 1` degrees of freedom): the `friedman_periods` row, 10 units x 30 periods, has the right size.  The default stays that of Stata, and the result carries a remark when `N` and `T` differ.

`cips` is simulated with a strong common factor in the innovations, the case it is made for.  The p-value of its `Z-t-bar` standardises as if the unit statistics were independent: it is oversized (about 10 to 12 %, 9.9 % in `pescadf` of Stata 17 over 1000 replications of the same design), whereas the decision of the paper, t-bar against its critical value, has about the right size (4 to 5 %, the `cips` row; 3.8 % in Stata).  `im_pesaran_shin`, which ignores the factor, rejects about a quarter of the time (the `im_pesaran_shin_common_factor` row).  The critical values exist for balanced panels only.

`pedroni` and `pedroni_panel` are oversized at 20 units x 40 periods (13 to 15 % at 5 %), as Pedroni (2004) finds for his tests at a short T: the statistics are those of Stata's `xtcointtest pedroni`, to 1e-9 on the validation cases, so the behaviour is theirs.  The variance-ratio test of `westerlund` has about the right size in the same design.

`dumitrescu_hurlin` is slightly oversized (8 %) with 20 units and 30 periods: its standardisation by the exact moments of the Wald statistics assumes normal errors and a large number of units.
