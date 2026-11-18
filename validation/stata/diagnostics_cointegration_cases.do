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

* ---- case gap_kao
use "../data/sim_gaps.dta", clear
quietly xtset id t
capture quietly xtcointtest kao x z
scalar rc = _rc
mata: pyx_export("../out/gap_kao__rc.csv", (st_numscalar("rc")))

* ---- case gap_pedroni
use "../data/sim_gaps.dta", clear
quietly xtset id t
capture quietly xtcointtest pedroni x z, lags(1)
scalar rc = _rc
mata: pyx_export("../out/gap_pedroni__rc.csv", (st_numscalar("rc")))

* ---- case gap_westerlund
use "../data/sim_gaps.dta", clear
quietly xtset id t
capture quietly xtcointtest westerlund x z
scalar rc = _rc
mata: pyx_export("../out/gap_westerlund__rc.csv", (st_numscalar("rc")))

* ---- case kao_default
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_default__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_default_t
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_default_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_nocoint
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_nocoint__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_nocoint_t
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_nocoint_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_bartlett3
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(bartlett 3)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_bartlett3__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_bartlett3_t
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(bartlett 3)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_bartlett3_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_parzen
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(parzen 2) lags(2)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_parzen__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_parzen_t
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(parzen 2) lags(2)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_parzen_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_parzen_auto
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(parzen nwest)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_parzen_auto__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_parzen_auto_t
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(parzen nwest)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_parzen_auto_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_adf2
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, lags(2)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_adf2__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_adf2_t
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, lags(2)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_adf2_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_adf_aic
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, lags(aic 3)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_adf_aic__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_adf_aic_t
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, lags(aic 3)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_adf_aic_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_adf_bic
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, lags(bic 3)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_adf_bic__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_adf_bic_t
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, lags(bic 3)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_adf_bic_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_adf_hqic
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, lags(hqic 4)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_adf_hqic__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_adf_hqic_t
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, lags(hqic 4)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_adf_hqic_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_demean
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, demean
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_demean__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_demean_t
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, demean
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_demean_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_one_regressor
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_one_regressor__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_one_regressor_t
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_one_regressor_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_qs
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(quadraticspectral 3)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_qs__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_qs_t
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(quadraticspectral 3)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_qs_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case kao_qs_auto
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(quadraticspectral nwest)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_qs_auto__exact.csv", (st_numscalar("s_modified_df") \ st_numscalar("s_unadjusted_modified_df") \ st_numscalar("p_modified_df") \ st_numscalar("p_unadjusted_modified_df") \ st_numscalar("k_lags") \ st_numscalar("k_adf") \ st_numscalar("k_ng") \ st_numscalar("k_nt") \ st_numscalar("k_n")))

* ---- case kao_qs_auto_t
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest kao y x z, kernel(quadraticspectral nwest)
matrix S = r(stats)
scalar s_modified_df = S[1, 1]
scalar p_modified_df = normal(S[1, 1])
scalar s_df = S[1, 2]
scalar p_df = normal(S[1, 2])
scalar s_adf = S[1, 3]
scalar p_adf = normal(S[1, 3])
scalar s_unadjusted_modified_df = S[1, 4]
scalar p_unadjusted_modified_df = normal(S[1, 4])
scalar s_unadjusted_df = S[1, 5]
scalar p_unadjusted_df = normal(S[1, 5])
scalar k_lags = r(hac_lagm)
scalar k_adf = r(adf_lags)
scalar k_ng = r(N_g)
scalar k_nt = r(N_t)
scalar k_n = r(N)
set type float
mata: pyx_export("../out/kao_qs_auto_t__rounded.csv", (st_numscalar("s_df") \ st_numscalar("s_adf") \ st_numscalar("s_unadjusted_df") \ st_numscalar("p_df") \ st_numscalar("p_adf") \ st_numscalar("p_unadjusted_df")))

* ---- case pe_group_m1
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x, kernel(bartlett 2) 
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_group_m1__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_panel_m1
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x, kernel(bartlett 2) ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_panel_m1__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_group_m3_trend
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x1 x2 x3, kernel(bartlett 2) trend
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_group_m3_trend__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_panel_m3_trend
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x1 x2 x3, kernel(bartlett 2) trend ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_panel_m3_trend__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_group_m2_none
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) noconstant
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_group_m2_none__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_panel_m2_none
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) noconstant ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_panel_m2_none__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_group_m7
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x1 x2 x3 x4 x5 x6 x7, kernel(bartlett 2) 
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_group_m7__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_panel_m7_trend
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x1 x2 x3 x4 x5 x6 x7, kernel(bartlett 2) trend ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_panel_m7_trend__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_parzen
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(parzen 3) ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_parzen__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_quadraticspectral
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(quadraticspectral 3)
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_quadraticspectral__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_kernel0
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 0) ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_kernel0__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_adf2
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) lags(2)
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_adf2__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_adf_aic
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) lags(aic 3)
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_adf_aic__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_adf_bic
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x1 x2, kernel(bartlett 2) lags(bic 3) ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_adf_bic__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_adf_hqic
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) lags(hqic 4)
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_adf_hqic__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_demean
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) demean
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_demean__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_unbalanced
use "../data/sim_coint_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_unbalanced__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_unbalanced_panel_trend
use "../data/sim_coint_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) trend
matrix S = r(stats)
scalar z_grho = S[1,1]
scalar p_grho = normal(z_grho)
scalar z_gtpp = S[1,2]
scalar p_gtpp = normal(z_gtpp)
scalar z_gadf = S[1,3]
scalar p_gadf = normal(z_gadf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_unbalanced_panel_trend__p.csv", (st_numscalar("z_grho") \ st_numscalar("p_grho") \ st_numscalar("z_gtpp") \ st_numscalar("p_gtpp") \ st_numscalar("z_gadf") \ st_numscalar("p_gadf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_unbalanced_none
use "../data/sim_coint_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x, kernel(bartlett 2) noconstant ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_unbalanced_none__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case pe_coint
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest pedroni y x z, kernel(bartlett 2) ar(same)
matrix S = r(stats)
scalar z_nu = S[1,1]
scalar p_nu = 1 - normal(z_nu)
scalar z_prho = S[1,2]
scalar p_prho = normal(z_prho)
scalar z_ptpp = S[1,3]
scalar p_ptpp = normal(z_ptpp)
scalar z_padf = S[1,4]
scalar p_padf = normal(z_padf)
scalar nl = r(adf_lags)
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/pe_coint__p.csv", (st_numscalar("z_nu") \ st_numscalar("p_nu") \ st_numscalar("z_prho") \ st_numscalar("p_prho") \ st_numscalar("z_ptpp") \ st_numscalar("p_ptpp") \ st_numscalar("z_padf") \ st_numscalar("p_padf") \ st_numscalar("nl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_some_m1
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_some_m1__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_all_m1
use "../data/sim_coint_no.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x, allpanels
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_all_m1__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_some_m3_trend
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x1 x2 x3, trend
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_some_m3_trend__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_all_m3_trend
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x1 x2 x3, trend allpanels
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_all_m3_trend__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_some_m7
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x1 x2 x3 x4 x5 x6 x7
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_some_m7__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_all_m7_trend
use "../data/sim_coint_wide.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x1 x2 x3 x4 x5 x6 x7, trend allpanels
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_all_m7_trend__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_demean
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x z, demean
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_demean__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_coint
use "../data/sim_coint_yes.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x z
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_coint__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_unbalanced
use "../data/sim_coint_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x z
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_unbalanced__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_unbalanced_all_trend
use "../data/sim_coint_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x z, trend allpanels
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_unbalanced_all_trend__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

* ---- case we_unbalanced_demean_all
use "../data/sim_coint_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtcointtest westerlund y x, demean allpanels
scalar st = r(stat)
scalar pl = normal(r(stat))
scalar ng = r(N_g)
scalar nn = r(N)
set type float
mata: pyx_export("../out/we_unbalanced_demean_all__w.csv", (st_numscalar("st") \ st_numscalar("pl") \ st_numscalar("ng") \ st_numscalar("nn")))

exit, clear
