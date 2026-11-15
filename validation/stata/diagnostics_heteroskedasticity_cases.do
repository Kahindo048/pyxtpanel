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

* ---- case bp_fitted
use "../data/sim_het.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
quietly predict double fit, xb
mata: pyx_export("../out/bp_fitted__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bp_fitted__fit.csv", st_data(., "fit"))
mata: pyx_export("../out/bp_fitted__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat hettest
scalar bp0 = r(chi2)
scalar bp1 = r(df)
scalar bp2 = r(p)
mata: pyx_export("../out/bp_fitted__bp.csv", (st_numscalar("bp0") \ st_numscalar("bp1") \ st_numscalar("bp2")))

* ---- case bp_rhs
use "../data/sim_het.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
quietly predict double fit, xb
mata: pyx_export("../out/bp_rhs__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bp_rhs__fit.csv", st_data(., "fit"))
mata: pyx_export("../out/bp_rhs__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat hettest, rhs
scalar bp0 = r(chi2)
scalar bp1 = r(df)
scalar bp2 = r(p)
mata: pyx_export("../out/bp_rhs__bp.csv", (st_numscalar("bp0") \ st_numscalar("bp1") \ st_numscalar("bp2")))

* ---- case bp_variables
use "../data/sim_het.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
quietly predict double fit, xb
mata: pyx_export("../out/bp_variables__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bp_variables__fit.csv", st_data(., "fit"))
mata: pyx_export("../out/bp_variables__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat hettest z1 z2
scalar bp0 = r(chi2)
scalar bp1 = r(df)
scalar bp2 = r(p)
mata: pyx_export("../out/bp_variables__bp.csv", (st_numscalar("bp0") \ st_numscalar("bp1") \ st_numscalar("bp2")))

* ---- case bp_iid
use "../data/sim_het.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
quietly predict double fit, xb
mata: pyx_export("../out/bp_iid__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bp_iid__fit.csv", st_data(., "fit"))
mata: pyx_export("../out/bp_iid__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat hettest, iid
scalar bp0 = r(chi2)
scalar bp1 = r(df)
scalar bp2 = r(p)
mata: pyx_export("../out/bp_iid__bp.csv", (st_numscalar("bp0") \ st_numscalar("bp1") \ st_numscalar("bp2")))

* ---- case bp_fstat
use "../data/sim_het.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
quietly predict double fit, xb
mata: pyx_export("../out/bp_fstat__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bp_fstat__fit.csv", st_data(., "fit"))
mata: pyx_export("../out/bp_fstat__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat hettest, fstat
scalar bp0 = r(F)
scalar bp1 = r(df_m)
scalar bp2 = r(df_r)
scalar bp3 = r(p)
mata: pyx_export("../out/bp_fstat__bp.csv", (st_numscalar("bp0") \ st_numscalar("bp1") \ st_numscalar("bp2") \ st_numscalar("bp3")))

* ---- case bp_homoskedastic
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
quietly predict double fit, xb
mata: pyx_export("../out/bp_homoskedastic__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bp_homoskedastic__fit.csv", st_data(., "fit"))
mata: pyx_export("../out/bp_homoskedastic__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat hettest
scalar bp0 = r(chi2)
scalar bp1 = r(df)
scalar bp2 = r(p)
mata: pyx_export("../out/bp_homoskedastic__bp.csv", (st_numscalar("bp0") \ st_numscalar("bp1") \ st_numscalar("bp2")))

* ---- case wg_balanced
use "../data/sim_groupwise.dta", clear
quietly xtset id t
quietly xtreg y x z, fe
quietly predict double ehat, e
quietly xttest3
scalar ww = r(wald)
scalar wdf = r(df)
scalar wp = r(p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/wg_balanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/wg_balanced__id.csv", st_data(., "id"))
restore
mata: pyx_export("../out/wg_balanced__w.csv", (st_numscalar("ww") \ st_numscalar("wdf") \ st_numscalar("wp")))

* ---- case wg_unbalanced
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly xtreg x z, fe
quietly predict double ehat, e
quietly xttest3
scalar ww = r(wald)
scalar wdf = r(df)
scalar wp = r(p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/wg_unbalanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/wg_unbalanced__id.csv", st_data(., "id"))
restore
mata: pyx_export("../out/wg_unbalanced__w.csv", (st_numscalar("ww") \ st_numscalar("wdf") \ st_numscalar("wp")))

* ---- case wg_short_panel
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly xtreg y x z1 z2, fe
quietly predict double ehat, e
quietly xttest3
scalar ww = r(wald)
scalar wdf = r(df)
scalar wp = r(p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/wg_short_panel__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/wg_short_panel__id.csv", st_data(., "id"))
restore
mata: pyx_export("../out/wg_short_panel__w.csv", (st_numscalar("ww") \ st_numscalar("wdf") \ st_numscalar("wp")))

* ---- case white_one
use "../data/sim_het.dta", clear
quietly xtset id t
quietly regress y x
quietly predict double ehat, resid
mata: pyx_export("../out/white_one__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/white_one__x.csv", st_data(., ("x")))
quietly estat imtest, white double
scalar w0 = r(chi2)
scalar w1 = r(df)
scalar w2 = r(p)
mata: pyx_export("../out/white_one__w.csv", (st_numscalar("w0") \ st_numscalar("w1") \ st_numscalar("w2")))

* ---- case white_three
use "../data/sim_het.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
mata: pyx_export("../out/white_three__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/white_three__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat imtest, white double
scalar w0 = r(chi2)
scalar w1 = r(df)
scalar w2 = r(p)
mata: pyx_export("../out/white_three__w.csv", (st_numscalar("w0") \ st_numscalar("w1") \ st_numscalar("w2")))

* ---- case white_four
use "../data/sim_het.dta", clear
quietly xtset id t
quietly regress y x z1 z2 d
quietly predict double ehat, resid
mata: pyx_export("../out/white_four__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/white_four__x.csv", st_data(., ("x", "z1", "z2", "d")))
quietly estat imtest, white double
scalar w0 = r(chi2)
scalar w1 = r(df)
scalar w2 = r(p)
mata: pyx_export("../out/white_four__w.csv", (st_numscalar("w0") \ st_numscalar("w1") \ st_numscalar("w2")))

* ---- case white_homoskedastic
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
mata: pyx_export("../out/white_homoskedastic__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/white_homoskedastic__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat imtest, white double
scalar w0 = r(chi2)
scalar w1 = r(df)
scalar w2 = r(p)
mata: pyx_export("../out/white_homoskedastic__w.csv", (st_numscalar("w0") \ st_numscalar("w1") \ st_numscalar("w2")))

exit, clear
