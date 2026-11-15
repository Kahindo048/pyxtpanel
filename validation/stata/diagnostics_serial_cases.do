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

* ---- case pd_autocorrelated
use "../data/sim_serial.dta", clear
quietly xtset id t, delta(1)
quietly xtregar y x z, fe lbi
scalar dw = e(d1)
scalar lbi = e(LBI)
scalar nlbi = e(N_LBI)
mata: pyx_export("../out/pd_autocorrelated__d.csv", (st_numscalar("dw") \ st_numscalar("lbi") \ st_numscalar("nlbi")))

* ---- case pd_independent
use "../data/sim_balanced.dta", clear
quietly xtset id t, delta(1)
quietly xtregar x z, fe lbi
scalar dw = e(d1)
scalar lbi = e(LBI)
scalar nlbi = e(N_LBI)
mata: pyx_export("../out/pd_independent__d.csv", (st_numscalar("dw") \ st_numscalar("lbi") \ st_numscalar("nlbi")))

* ---- case pd_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t, delta(1)
quietly xtregar x z, fe lbi
scalar dw = e(d1)
scalar lbi = e(LBI)
scalar nlbi = e(N_LBI)
mata: pyx_export("../out/pd_gaps__d.csv", (st_numscalar("dw") \ st_numscalar("lbi") \ st_numscalar("nlbi")))

* ---- case pd_entry_exit
use "../data/sim_entry_exit.dta", clear
quietly xtset id t, delta(1)
quietly xtregar x z, fe lbi
scalar dw = e(d1)
scalar lbi = e(LBI)
scalar nlbi = e(N_LBI)
mata: pyx_export("../out/pd_entry_exit__d.csv", (st_numscalar("dw") \ st_numscalar("lbi") \ st_numscalar("nlbi")))

* ---- case pd_delta2
use "../data/sim_delta2.dta", clear
quietly xtset id t, delta(2)
quietly xtregar x z, fe lbi
scalar dw = e(d1)
scalar lbi = e(LBI)
scalar nlbi = e(N_LBI)
mata: pyx_export("../out/pd_delta2__d.csv", (st_numscalar("dw") \ st_numscalar("lbi") \ st_numscalar("nlbi")))

* ---- case pd_singleton
use "../data/sim_one_obs_unit.dta", clear
quietly xtset id t, delta(1)
quietly xtregar x z, fe lbi
scalar dw = e(d1)
scalar lbi = e(LBI)
scalar nlbi = e(N_LBI)
mata: pyx_export("../out/pd_singleton__d.csv", (st_numscalar("dw") \ st_numscalar("lbi") \ st_numscalar("nlbi")))

* ---- case wr_autocorrelated
use "../data/sim_serial.dta", clear
quietly xtset id t
quietly xtserial y x z
scalar ff = r(F)
scalar fdf = r(df)
scalar fdfr = r(df_r)
scalar fp = r(p)
scalar frho = r(corr)
mata: pyx_export("../out/wr_autocorrelated__w.csv", (st_numscalar("ff") \ st_numscalar("fdf") \ st_numscalar("fdfr") \ st_numscalar("fp") \ st_numscalar("frho")))

* ---- case wr_independent
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtserial x z
scalar ff = r(F)
scalar fdf = r(df)
scalar fdfr = r(df_r)
scalar fp = r(p)
scalar frho = r(corr)
mata: pyx_export("../out/wr_independent__w.csv", (st_numscalar("ff") \ st_numscalar("fdf") \ st_numscalar("fdfr") \ st_numscalar("fp") \ st_numscalar("frho")))

* ---- case wr_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly xtserial x z
scalar ff = r(F)
scalar fdf = r(df)
scalar fdfr = r(df_r)
scalar fp = r(p)
scalar frho = r(corr)
mata: pyx_export("../out/wr_gaps__w.csv", (st_numscalar("ff") \ st_numscalar("fdf") \ st_numscalar("fdfr") \ st_numscalar("fp") \ st_numscalar("frho")))

* ---- case wr_entry_exit
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly xtserial x z
scalar ff = r(F)
scalar fdf = r(df)
scalar fdfr = r(df_r)
scalar fp = r(p)
scalar frho = r(corr)
mata: pyx_export("../out/wr_entry_exit__w.csv", (st_numscalar("ff") \ st_numscalar("fdf") \ st_numscalar("fdfr") \ st_numscalar("fp") \ st_numscalar("frho")))

exit, clear
