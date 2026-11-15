version 17
set more off
clear all

mata:
void pyx_export(string scalar fname, real matrix M)
{
    real scalar fh, i, j
    string scalar line
    if (fileexists(fname)) unlink(fname)
    fh = fopen(fname, "w")
    for (i = 1; i <= rows(M); i++) {
        line = strofreal(M[i,1], "%21.15g")
        for (j = 2; j <= cols(M); j++) line = line + "," + strofreal(M[i,j], "%21.15g")
        fput(fh, line)
    }
    fclose(fh)
}
end

mata: pyx_export("../out/_meta__stata.csv", (st_numscalar("c(stata_version)") \ st_numscalar("c(MP)")))

* ---- case fi_lags1
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot fisher y, dfuller lags(1)  generate(pu)
scalar fp = r(P)
scalar fdfp = r(df_P)
scalar fpp = r(p_P)
scalar fz = r(Z)
scalar fpz = r(p_Z)
scalar fl = r(L)
scalar fdfl = r(df_L)
scalar fpl = r(p_L)
scalar fpm = r(Pm)
scalar fppm = r(p_Pm)
scalar fng = r(N_g)
scalar fn = r(N)
preserve
quietly keep if !missing(pu)
mata: pyx_export("../out/fi_lags1__unit_p.csv", st_data(., "pu"))
restore
mata: pyx_export("../out/fi_lags1__f.csv", (st_numscalar("fp") \ st_numscalar("fdfp") \ st_numscalar("fpp") \ st_numscalar("fz") \ st_numscalar("fpz") \ st_numscalar("fl") \ st_numscalar("fdfl") \ st_numscalar("fpl") \ st_numscalar("fpm") \ st_numscalar("fppm") \ st_numscalar("fng") \ st_numscalar("fn")))

* ---- case fi_trend_demean
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot fisher y, dfuller lags(1) trend demean generate(pu)
scalar fp = r(P)
scalar fdfp = r(df_P)
scalar fpp = r(p_P)
scalar fz = r(Z)
scalar fpz = r(p_Z)
scalar fl = r(L)
scalar fdfl = r(df_L)
scalar fpl = r(p_L)
scalar fpm = r(Pm)
scalar fppm = r(p_Pm)
scalar fng = r(N_g)
scalar fn = r(N)
preserve
quietly keep if !missing(pu)
mata: pyx_export("../out/fi_trend_demean__unit_p.csv", st_data(., "pu"))
restore
mata: pyx_export("../out/fi_trend_demean__f.csv", (st_numscalar("fp") \ st_numscalar("fdfp") \ st_numscalar("fpp") \ st_numscalar("fz") \ st_numscalar("fpz") \ st_numscalar("fl") \ st_numscalar("fdfl") \ st_numscalar("fpl") \ st_numscalar("fpm") \ st_numscalar("fppm") \ st_numscalar("fng") \ st_numscalar("fn")))

* ---- case fi_drift
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot fisher y, dfuller lags(2) drift generate(pu)
scalar fp = r(P)
scalar fdfp = r(df_P)
scalar fpp = r(p_P)
scalar fz = r(Z)
scalar fpz = r(p_Z)
scalar fl = r(L)
scalar fdfl = r(df_L)
scalar fpl = r(p_L)
scalar fpm = r(Pm)
scalar fppm = r(p_Pm)
scalar fng = r(N_g)
scalar fn = r(N)
preserve
quietly keep if !missing(pu)
mata: pyx_export("../out/fi_drift__unit_p.csv", st_data(., "pu"))
restore
mata: pyx_export("../out/fi_drift__f.csv", (st_numscalar("fp") \ st_numscalar("fdfp") \ st_numscalar("fpp") \ st_numscalar("fz") \ st_numscalar("fpz") \ st_numscalar("fl") \ st_numscalar("fdfl") \ st_numscalar("fpl") \ st_numscalar("fpm") \ st_numscalar("fppm") \ st_numscalar("fng") \ st_numscalar("fn")))

* ---- case fi_entry_exit
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly xtunitroot fisher x, dfuller lags(0) demean generate(pu)
scalar fp = r(P)
scalar fdfp = r(df_P)
scalar fpp = r(p_P)
scalar fz = r(Z)
scalar fpz = r(p_Z)
scalar fl = r(L)
scalar fdfl = r(df_L)
scalar fpl = r(p_L)
scalar fpm = r(Pm)
scalar fppm = r(p_Pm)
scalar fng = r(N_g)
scalar fn = r(N)
preserve
quietly keep if !missing(pu)
mata: pyx_export("../out/fi_entry_exit__unit_p.csv", st_data(., "pu"))
restore
mata: pyx_export("../out/fi_entry_exit__f.csv", (st_numscalar("fp") \ st_numscalar("fdfp") \ st_numscalar("fpp") \ st_numscalar("fz") \ st_numscalar("fpz") \ st_numscalar("fl") \ st_numscalar("fdfl") \ st_numscalar("fpl") \ st_numscalar("fpm") \ st_numscalar("fppm") \ st_numscalar("fng") \ st_numscalar("fn")))

* ---- case fi_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly xtunitroot fisher x, dfuller lags(0)  generate(pu)
scalar fp = r(P)
scalar fdfp = r(df_P)
scalar fpp = r(p_P)
scalar fz = r(Z)
scalar fpz = r(p_Z)
scalar fl = r(L)
scalar fdfl = r(df_L)
scalar fpl = r(p_L)
scalar fpm = r(Pm)
scalar fppm = r(p_Pm)
scalar fng = r(N_g)
scalar fn = r(N)
preserve
quietly keep if !missing(pu)
mata: pyx_export("../out/fi_gaps__unit_p.csv", st_data(., "pu"))
restore
mata: pyx_export("../out/fi_gaps__f.csv", (st_numscalar("fp") \ st_numscalar("fdfp") \ st_numscalar("fpp") \ st_numscalar("fz") \ st_numscalar("fpz") \ st_numscalar("fl") \ st_numscalar("fdfl") \ st_numscalar("fpl") \ st_numscalar("fpm") \ st_numscalar("fppm") \ st_numscalar("fng") \ st_numscalar("fn")))

* ---- case fi_short_units
use "../data/sim_one_obs_unit.dta", clear
quietly xtset id t
quietly xtunitroot fisher x, dfuller lags(0)  generate(pu)
scalar fp = r(P)
scalar fdfp = r(df_P)
scalar fpp = r(p_P)
scalar fz = r(Z)
scalar fpz = r(p_Z)
scalar fl = r(L)
scalar fdfl = r(df_L)
scalar fpl = r(p_L)
scalar fpm = r(Pm)
scalar fppm = r(p_Pm)
scalar fng = r(N_g)
scalar fn = r(N)
preserve
quietly keep if !missing(pu)
mata: pyx_export("../out/fi_short_units__unit_p.csv", st_data(., "pu"))
restore
mata: pyx_export("../out/fi_short_units__f.csv", (st_numscalar("fp") \ st_numscalar("fdfp") \ st_numscalar("fpp") \ st_numscalar("fz") \ st_numscalar("fpz") \ st_numscalar("fl") \ st_numscalar("fdfl") \ st_numscalar("fpl") \ st_numscalar("fpm") \ st_numscalar("fppm") \ st_numscalar("fng") \ st_numscalar("fn")))

* ---- case ht_walk
use "../data/sim_ur_walk.dta", clear
quietly xtset id t
quietly xtunitroot ht y
scalar htrho = r(rho)
scalar htmean = r(mean_rho)
scalar htvar = r(Var_rho)
scalar htz = r(z)
scalar htp = r(p)
mata: pyx_export("../out/ht_walk__ht.csv", (st_numscalar("htrho") \ st_numscalar("htmean") \ st_numscalar("htvar") \ st_numscalar("htz") \ st_numscalar("htp")))

* ---- case ht_stationary
use "../data/sim_ur_stationary.dta", clear
quietly xtset id t
quietly xtunitroot ht y
scalar htrho = r(rho)
scalar htmean = r(mean_rho)
scalar htvar = r(Var_rho)
scalar htz = r(z)
scalar htp = r(p)
mata: pyx_export("../out/ht_stationary__ht.csv", (st_numscalar("htrho") \ st_numscalar("htmean") \ st_numscalar("htvar") \ st_numscalar("htz") \ st_numscalar("htp")))

* ---- case ht_none
use "../data/sim_ur_walk.dta", clear
quietly xtset id t
quietly xtunitroot ht y, noconstant
scalar htrho = r(rho)
scalar htmean = r(mean_rho)
scalar htvar = r(Var_rho)
scalar htz = r(z)
scalar htp = r(p)
mata: pyx_export("../out/ht_none__ht.csv", (st_numscalar("htrho") \ st_numscalar("htmean") \ st_numscalar("htvar") \ st_numscalar("htz") \ st_numscalar("htp")))

* ---- case ht_trend
use "../data/sim_ur_trend.dta", clear
quietly xtset id t
quietly xtunitroot ht y, trend
scalar htrho = r(rho)
scalar htmean = r(mean_rho)
scalar htvar = r(Var_rho)
scalar htz = r(z)
scalar htp = r(p)
mata: pyx_export("../out/ht_trend__ht.csv", (st_numscalar("htrho") \ st_numscalar("htmean") \ st_numscalar("htvar") \ st_numscalar("htz") \ st_numscalar("htp")))

* ---- case ht_demean
use "../data/sim_ur_stationary.dta", clear
quietly xtset id t
quietly xtunitroot ht y, demean
scalar htrho = r(rho)
scalar htmean = r(mean_rho)
scalar htvar = r(Var_rho)
scalar htz = r(z)
scalar htp = r(p)
mata: pyx_export("../out/ht_demean__ht.csv", (st_numscalar("htrho") \ st_numscalar("htmean") \ st_numscalar("htvar") \ st_numscalar("htz") \ st_numscalar("htp")))

* ---- case ht_altt
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot ht y, altt
scalar htrho = r(rho)
scalar htmean = r(mean_rho)
scalar htvar = r(Var_rho)
scalar htz = r(z)
scalar htp = r(p)
mata: pyx_export("../out/ht_altt__ht.csv", (st_numscalar("htrho") \ st_numscalar("htmean") \ st_numscalar("htvar") \ st_numscalar("htz") \ st_numscalar("htp")))

* ---- case ips_balanced
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot ips y
scalar i_zttildebar = r(zttildebar)
scalar i_p_zttildebar = r(p_zttildebar)
scalar i_tbar = r(tbar)
scalar i_ttildebar = r(ttildebar)
scalar i_zt = r(zt)
scalar i_cv_1 = r(cv_1)
scalar i_cv_5 = r(cv_5)
scalar i_cv_10 = r(cv_10)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_balanced__ips.csv", (st_numscalar("i_zttildebar") \ st_numscalar("i_p_zttildebar") \ st_numscalar("i_tbar") \ st_numscalar("i_ttildebar") \ st_numscalar("i_zt") \ st_numscalar("i_cv_1") \ st_numscalar("i_cv_5") \ st_numscalar("i_cv_10") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_short
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtunitroot ips x
scalar i_zttildebar = r(zttildebar)
scalar i_p_zttildebar = r(p_zttildebar)
scalar i_tbar = r(tbar)
scalar i_ttildebar = r(ttildebar)
scalar i_zt = r(zt)
scalar i_cv_1 = r(cv_1)
scalar i_cv_5 = r(cv_5)
scalar i_cv_10 = r(cv_10)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_short__ips.csv", (st_numscalar("i_zttildebar") \ st_numscalar("i_p_zttildebar") \ st_numscalar("i_tbar") \ st_numscalar("i_ttildebar") \ st_numscalar("i_zt") \ st_numscalar("i_cv_1") \ st_numscalar("i_cv_5") \ st_numscalar("i_cv_10") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_trend
use "../data/sim_ur_trend.dta", clear
quietly xtset id t
quietly xtunitroot ips y, trend
scalar i_zttildebar = r(zttildebar)
scalar i_p_zttildebar = r(p_zttildebar)
scalar i_tbar = r(tbar)
scalar i_ttildebar = r(ttildebar)
scalar i_zt = r(zt)
scalar i_cv_1 = r(cv_1)
scalar i_cv_5 = r(cv_5)
scalar i_cv_10 = r(cv_10)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_trend__ips.csv", (st_numscalar("i_zttildebar") \ st_numscalar("i_p_zttildebar") \ st_numscalar("i_tbar") \ st_numscalar("i_ttildebar") \ st_numscalar("i_zt") \ st_numscalar("i_cv_1") \ st_numscalar("i_cv_5") \ st_numscalar("i_cv_10") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_demean
use "../data/sim_ur_walk.dta", clear
quietly xtset id t
quietly xtunitroot ips y, demean
scalar i_zttildebar = r(zttildebar)
scalar i_p_zttildebar = r(p_zttildebar)
scalar i_tbar = r(tbar)
scalar i_ttildebar = r(ttildebar)
scalar i_zt = r(zt)
scalar i_cv_1 = r(cv_1)
scalar i_cv_5 = r(cv_5)
scalar i_cv_10 = r(cv_10)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_demean__ips.csv", (st_numscalar("i_zttildebar") \ st_numscalar("i_p_zttildebar") \ st_numscalar("i_tbar") \ st_numscalar("i_ttildebar") \ st_numscalar("i_zt") \ st_numscalar("i_cv_1") \ st_numscalar("i_cv_5") \ st_numscalar("i_cv_10") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_long
use "../data/sim_ur_long.dta", clear
quietly xtset id t
quietly xtunitroot ips y
scalar i_zttildebar = r(zttildebar)
scalar i_p_zttildebar = r(p_zttildebar)
scalar i_tbar = r(tbar)
scalar i_ttildebar = r(ttildebar)
scalar i_zt = r(zt)
scalar i_cv_1 = r(cv_1)
scalar i_cv_5 = r(cv_5)
scalar i_cv_10 = r(cv_10)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_long__ips.csv", (st_numscalar("i_zttildebar") \ st_numscalar("i_p_zttildebar") \ st_numscalar("i_tbar") \ st_numscalar("i_ttildebar") \ st_numscalar("i_zt") \ st_numscalar("i_cv_1") \ st_numscalar("i_cv_5") \ st_numscalar("i_cv_10") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_unbalanced
use "../data/sim_ur_unbalanced.dta", clear
quietly xtset id t
quietly xtunitroot ips y
scalar i_zttildebar = r(zttildebar)
scalar i_p_zttildebar = r(p_zttildebar)
scalar i_tbar = r(tbar)
scalar i_ttildebar = r(ttildebar)
scalar i_zt = r(zt)
scalar i_cv_1 = r(cv_1)
scalar i_cv_5 = r(cv_5)
scalar i_cv_10 = r(cv_10)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_unbalanced__ips.csv", (st_numscalar("i_zttildebar") \ st_numscalar("i_p_zttildebar") \ st_numscalar("i_tbar") \ st_numscalar("i_ttildebar") \ st_numscalar("i_zt") \ st_numscalar("i_cv_1") \ st_numscalar("i_cv_5") \ st_numscalar("i_cv_10") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_unbalanced_short
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly xtunitroot ips x
scalar i_zttildebar = r(zttildebar)
scalar i_p_zttildebar = r(p_zttildebar)
scalar i_tbar = r(tbar)
scalar i_ttildebar = r(ttildebar)
scalar i_zt = r(zt)
scalar i_cv_1 = r(cv_1)
scalar i_cv_5 = r(cv_5)
scalar i_cv_10 = r(cv_10)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_unbalanced_short__ips.csv", (st_numscalar("i_zttildebar") \ st_numscalar("i_p_zttildebar") \ st_numscalar("i_tbar") \ st_numscalar("i_ttildebar") \ st_numscalar("i_zt") \ st_numscalar("i_cv_1") \ st_numscalar("i_cv_5") \ st_numscalar("i_cv_10") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_lags1
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot ips y, lags(1)
scalar i_wtbar = r(wtbar)
scalar i_p_wtbar = r(p_wtbar)
scalar i_lags = r(lags)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_lags1__ips.csv", (st_numscalar("i_wtbar") \ st_numscalar("i_p_wtbar") \ st_numscalar("i_lags") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_lags_aic
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot ips y, lags(aic 3)
scalar i_wtbar = r(wtbar)
scalar i_p_wtbar = r(p_wtbar)
scalar i_lags = r(lagm)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_lags_aic__ips.csv", (st_numscalar("i_wtbar") \ st_numscalar("i_p_wtbar") \ st_numscalar("i_lags") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_lags_trend_unbalanced
use "../data/sim_ur_unbalanced.dta", clear
quietly xtset id t
quietly xtunitroot ips y, lags(2) trend
scalar i_wtbar = r(wtbar)
scalar i_p_wtbar = r(p_wtbar)
scalar i_lags = r(lags)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_lags_trend_unbalanced__ips.csv", (st_numscalar("i_wtbar") \ st_numscalar("i_p_wtbar") \ st_numscalar("i_lags") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_lags_hqic_long
use "../data/sim_ur_long.dta", clear
quietly xtset id t
quietly xtunitroot ips y, lags(hqic 4) trend
scalar i_wtbar = r(wtbar)
scalar i_p_wtbar = r(p_wtbar)
scalar i_lags = r(lagm)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_lags_hqic_long__ips.csv", (st_numscalar("i_wtbar") \ st_numscalar("i_p_wtbar") \ st_numscalar("i_lags") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case ips_nointerpolate
use "../data/sim_ur_unbalanced.dta", clear
quietly xtset id t
quietly xtunitroot ips y, lags(1) nointerpolate
scalar i_wtbar = r(wtbar)
scalar i_p_wtbar = r(p_wtbar)
scalar i_lags = r(lags)
scalar i_N = r(N)
scalar i_N_g = r(N_g)
scalar i_N_t = r(N_t)
mata: pyx_export("../out/ips_nointerpolate__ips.csv", (st_numscalar("i_wtbar") \ st_numscalar("i_p_wtbar") \ st_numscalar("i_lags") \ st_numscalar("i_N") \ st_numscalar("i_N_g") \ st_numscalar("i_N_t")))

* ---- case llc_lags1
use "../data/sim_ur_long.dta", clear
quietly xtset id t
quietly xtunitroot llc y, lags(1)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_lags1__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_none
use "../data/sim_ur_long.dta", clear
quietly xtset id t
quietly xtunitroot llc y, noconstant lags(1)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_none__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_none_aic
use "../data/sim_ur_walk.dta", clear
quietly xtset id t
quietly xtunitroot llc y, noconstant lags(aic 3)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lagm)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_none_aic__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_trend
use "../data/sim_ur_long.dta", clear
quietly xtset id t
quietly xtunitroot llc y, trend lags(2)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_trend__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_demean
use "../data/sim_ur_stationary.dta", clear
quietly xtset id t
quietly xtunitroot llc y, demean lags(1)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_demean__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_aic
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot llc y, lags(aic 3)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lagm)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_aic__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_hqic_parzen
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot llc y, lags(hqic 3) kernel(parzen 3)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lagm)
scalar l_hac = r(hac_lags)
mata: pyx_export("../out/llc_hqic_parzen__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_parzen_nwest
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot llc y, lags(2) kernel(parzen nwest)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_parzen_nwest__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_bartlett_nwest
use "../data/sim_ur_long.dta", clear
quietly xtset id t
quietly xtunitroot llc y, lags(1) kernel(bartlett nwest)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_bartlett_nwest__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_qs_fixed
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot llc y, lags(1) kernel(andrews 4.5)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lags)
mata: pyx_export("../out/llc_qs_fixed__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_qs_nwest
use "../data/sim_ur_mixed.dta", clear
quietly xtset id t
quietly xtunitroot llc y, lags(1) kernel(andrews nwest)
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_qs_nwest__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

* ---- case llc_nointerpolate
use "../data/sim_ur_long.dta", clear
quietly xtset id t
quietly xtunitroot llc y, lags(1) nointerpolate
scalar l_tds = r(tds)
scalar l_td = r(td)
scalar l_delta = r(delta)
scalar l_se_delta = r(se_delta)
scalar l_mu_adj = r(mu_adj)
scalar l_sig_adj = r(sig_adj)
scalar l_ttilde = r(ttilde)
scalar l_sbar = r(sbar)
scalar l_Var_ep = r(Var_ep)
scalar l_p_tds = r(p_tds)
scalar l_N = r(N)
scalar l_N_g = r(N_g)
scalar l_N_t = r(N_t)
scalar l_adf = r(adf_lags)
scalar l_hac = r(hac_lagm)
mata: pyx_export("../out/llc_nointerpolate__llc.csv", (st_numscalar("l_tds") \ st_numscalar("l_td") \ st_numscalar("l_delta") \ st_numscalar("l_se_delta") \ st_numscalar("l_mu_adj") \ st_numscalar("l_sig_adj") \ st_numscalar("l_ttilde") \ st_numscalar("l_sbar") \ st_numscalar("l_Var_ep") \ st_numscalar("l_p_tds") \ st_numscalar("l_N") \ st_numscalar("l_N_g") \ st_numscalar("l_N_t") \ st_numscalar("l_adf") \ st_numscalar("l_hac")))

exit, clear
