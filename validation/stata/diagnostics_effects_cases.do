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

* ---- case fe_balanced
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly regress y x z1 z2
scalar rssp = e(rss)
quietly xtreg y x z1 z2, fe
scalar rsse = e(rss)
scalar nobs = e(N)
scalar ngrp = e(N_g)
scalar kreg = colsof(e(b)) - 1
scalar ff = e(F_f)
scalar fdf1 = e(df_a)
scalar fdf2 = e(df_r)
scalar fp = Ftail(e(df_a), e(df_r), e(F_f))
mata: pyx_export("../out/fe_balanced__inputs.csv", (st_numscalar("rssp") \ st_numscalar("rsse") \ st_numscalar("nobs") \ st_numscalar("ngrp") \ st_numscalar("kreg")))
mata: pyx_export("../out/fe_balanced__f.csv", (st_numscalar("ff") \ st_numscalar("fdf1") \ st_numscalar("fdf2") \ st_numscalar("fp")))

* ---- case fe_unbalanced
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly regress x z
scalar rssp = e(rss)
quietly xtreg x z, fe
scalar rsse = e(rss)
scalar nobs = e(N)
scalar ngrp = e(N_g)
scalar kreg = colsof(e(b)) - 1
scalar ff = e(F_f)
scalar fdf1 = e(df_a)
scalar fdf2 = e(df_r)
scalar fp = Ftail(e(df_a), e(df_r), e(F_f))
mata: pyx_export("../out/fe_unbalanced__inputs.csv", (st_numscalar("rssp") \ st_numscalar("rsse") \ st_numscalar("nobs") \ st_numscalar("ngrp") \ st_numscalar("kreg")))
mata: pyx_export("../out/fe_unbalanced__f.csv", (st_numscalar("ff") \ st_numscalar("fdf1") \ st_numscalar("fdf2") \ st_numscalar("fp")))

* ---- case fe_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly regress x z
scalar rssp = e(rss)
quietly xtreg x z, fe
scalar rsse = e(rss)
scalar nobs = e(N)
scalar ngrp = e(N_g)
scalar kreg = colsof(e(b)) - 1
scalar ff = e(F_f)
scalar fdf1 = e(df_a)
scalar fdf2 = e(df_r)
scalar fp = Ftail(e(df_a), e(df_r), e(F_f))
mata: pyx_export("../out/fe_gaps__inputs.csv", (st_numscalar("rssp") \ st_numscalar("rsse") \ st_numscalar("nobs") \ st_numscalar("ngrp") \ st_numscalar("kreg")))
mata: pyx_export("../out/fe_gaps__f.csv", (st_numscalar("ff") \ st_numscalar("fdf1") \ st_numscalar("fdf2") \ st_numscalar("fp")))

* ---- case lm_re_effect
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
quietly xtreg y x z1 z2, re
mata: pyx_export("../out/lm_re_effect__sigma_u.csv", st_numscalar("e(sigma_u)"))
quietly xttest0
scalar lm = r(lm)
scalar lmp = r(p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/lm_re_effect__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/lm_re_effect__id.csv", st_data(., "id"))
restore
mata: pyx_export("../out/lm_re_effect__lm.csv", (st_numscalar("lm") \ st_numscalar("lmp")))

* ---- case lm_re_unbalanced
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly regress x z
quietly predict double ehat, resid
quietly xtreg x z, re
mata: pyx_export("../out/lm_re_unbalanced__sigma_u.csv", st_numscalar("e(sigma_u)"))
quietly xttest0
scalar lm = r(lm)
scalar lmp = r(p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/lm_re_unbalanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/lm_re_unbalanced__id.csv", st_data(., "id"))
restore
mata: pyx_export("../out/lm_re_unbalanced__lm.csv", (st_numscalar("lm") \ st_numscalar("lmp")))

* ---- case lm_re_no_effect
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly regress x z
quietly predict double ehat, resid
quietly xtreg x z, re
mata: pyx_export("../out/lm_re_no_effect__sigma_u.csv", st_numscalar("e(sigma_u)"))
quietly xttest0
scalar lm = r(lm)
scalar lmp = r(p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/lm_re_no_effect__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/lm_re_no_effect__id.csv", st_data(., "id"))
restore
mata: pyx_export("../out/lm_re_no_effect__lm.csv", (st_numscalar("lm") \ st_numscalar("lmp")))

* ---- case lm_re_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly regress x z
quietly predict double ehat, resid
quietly xtreg x z, re
mata: pyx_export("../out/lm_re_gaps__sigma_u.csv", st_numscalar("e(sigma_u)"))
quietly xttest0
scalar lm = r(lm)
scalar lmp = r(p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/lm_re_gaps__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/lm_re_gaps__id.csv", st_data(., "id"))
restore
mata: pyx_export("../out/lm_re_gaps__lm.csv", (st_numscalar("lm") \ st_numscalar("lmp")))

exit, clear
