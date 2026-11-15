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
