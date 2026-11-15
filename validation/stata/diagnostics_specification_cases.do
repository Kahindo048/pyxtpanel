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

* ---- case reset_fitted
use "../data/sim_het.dta", clear
quietly generate double yq = y + 1.5 * x * x
quietly generate double xc = 2 * x + z1
quietly regress y x z1 z2
mata: pyx_export("../out/reset_fitted__y.csv", st_data(., "y"))
mata: pyx_export("../out/reset_fitted__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat ovtest
scalar rf = r(F)
scalar rdf = r(df)
scalar rdr = r(df_r)
scalar rp = r(p)
mata: pyx_export("../out/reset_fitted__r.csv", (st_numscalar("rf") \ st_numscalar("rdf") \ st_numscalar("rdr") \ st_numscalar("rp")))

* ---- case reset_misspecified
use "../data/sim_het.dta", clear
quietly generate double yq = y + 1.5 * x * x
quietly generate double xc = 2 * x + z1
quietly regress yq x z1 z2
mata: pyx_export("../out/reset_misspecified__y.csv", st_data(., "yq"))
mata: pyx_export("../out/reset_misspecified__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat ovtest
scalar rf = r(F)
scalar rdf = r(df)
scalar rdr = r(df_r)
scalar rp = r(p)
mata: pyx_export("../out/reset_misspecified__r.csv", (st_numscalar("rf") \ st_numscalar("rdf") \ st_numscalar("rdr") \ st_numscalar("rp")))

* ---- case reset_rhs
use "../data/sim_het.dta", clear
quietly generate double yq = y + 1.5 * x * x
quietly generate double xc = 2 * x + z1
quietly regress yq x z1 z2
mata: pyx_export("../out/reset_rhs__y.csv", st_data(., "yq"))
mata: pyx_export("../out/reset_rhs__x.csv", st_data(., ("x", "z1", "z2")))
quietly estat ovtest, rhs
scalar rf = r(F)
scalar rdf = r(df)
scalar rdr = r(df_r)
scalar rp = r(p)
mata: pyx_export("../out/reset_rhs__r.csv", (st_numscalar("rf") \ st_numscalar("rdf") \ st_numscalar("rdr") \ st_numscalar("rp")))

* ---- case reset_rhs_dummy
use "../data/sim_het.dta", clear
quietly generate double yq = y + 1.5 * x * x
quietly generate double xc = 2 * x + z1
quietly regress yq x z1 d
mata: pyx_export("../out/reset_rhs_dummy__y.csv", st_data(., "yq"))
mata: pyx_export("../out/reset_rhs_dummy__x.csv", st_data(., ("x", "z1", "d")))
quietly estat ovtest, rhs
scalar rf = r(F)
scalar rdf = r(df)
scalar rdr = r(df_r)
scalar rp = r(p)
mata: pyx_export("../out/reset_rhs_dummy__r.csv", (st_numscalar("rf") \ st_numscalar("rdf") \ st_numscalar("rdr") \ st_numscalar("rp")))

* ---- case reset_collinear
use "../data/sim_het.dta", clear
quietly generate double yq = y + 1.5 * x * x
quietly generate double xc = 2 * x + z1
quietly regress yq x z1 xc
mata: pyx_export("../out/reset_collinear__y.csv", st_data(., "yq"))
mata: pyx_export("../out/reset_collinear__x.csv", st_data(., ("x", "z1", "xc")))
quietly estat ovtest
scalar rf = r(F)
scalar rdf = r(df)
scalar rdr = r(df_r)
scalar rp = r(p)
mata: pyx_export("../out/reset_collinear__r.csv", (st_numscalar("rf") \ st_numscalar("rdf") \ st_numscalar("rdr") \ st_numscalar("rp")))

* ---- case reset_fitted_exact
use "../data/sim_het.dta", clear
quietly generate double yq = y + 1.5 * x * x
quietly generate double xc = 2 * x + z1
quietly regress yq x z1 z2
mata: pyx_export("../out/reset_fitted_exact__y.csv", st_data(., "yq"))
mata: pyx_export("../out/reset_fitted_exact__x.csv", st_data(., ("x", "z1", "z2")))
quietly predict double yh, xb
quietly summarize yh
quietly replace yh = (yh - r(min)) / (r(max) - r(min))
quietly generate double yh2 = yh^2
quietly generate double yh3 = yh^3
quietly generate double yh4 = yh^4
quietly regress yq x z1 z2 yh2 yh3 yh4
quietly test yh2 yh3 yh4
scalar rf = r(F)
scalar rdf = r(df)
scalar rdr = r(df_r)
scalar rp = r(p)
mata: pyx_export("../out/reset_fitted_exact__r.csv", (st_numscalar("rf") \ st_numscalar("rdf") \ st_numscalar("rdr") \ st_numscalar("rp")))

* ---- case reset_rhs_exact
use "../data/sim_het.dta", clear
quietly generate double yq = y + 1.5 * x * x
quietly generate double xc = 2 * x + z1
quietly regress yq x z1 z2
mata: pyx_export("../out/reset_rhs_exact__y.csv", st_data(., "yq"))
mata: pyx_export("../out/reset_rhs_exact__x.csv", st_data(., ("x", "z1", "z2")))
foreach v in x z1 z2 {
    quietly summarize `v'
    quietly generate double s_`v' = (`v' - r(min)) / (r(max) - r(min))
    quietly generate double s2_`v' = s_`v'^2
    quietly generate double s3_`v' = s_`v'^3
    quietly generate double s4_`v' = s_`v'^4
}
quietly regress yq x z1 z2 s2_x s3_x s4_x s2_z1 s3_z1 s4_z1 s2_z2 s3_z2 s4_z2
quietly test s2_x s3_x s4_x s2_z1 s3_z1 s4_z1 s2_z2 s3_z2 s4_z2
scalar rf = r(F)
scalar rdf = r(df)
scalar rdr = r(df_r)
scalar rp = r(p)
mata: pyx_export("../out/reset_rhs_exact__r.csv", (st_numscalar("rf") \ st_numscalar("rdf") \ st_numscalar("rdr") \ st_numscalar("rp")))

* ---- case vif_correlated
use "../data/sim_het.dta", clear
quietly generate double w1 = 0.8 * x + 0.6 * z1
quietly generate double w2 = 0.5 * x - 0.5 * w1 + 0.3 * z2
quietly regress y x w1 w2
mata: pyx_export("../out/vif_correlated__x.csv", st_data(., ("x", "w1", "w2")))
set type double
quietly estat vif, nosort
set type float
scalar v1 = r(vif_1)
scalar v2 = r(vif_2)
scalar v3 = r(vif_3)
mata: pyx_export("../out/vif_correlated__v.csv", (st_numscalar("v1") \ st_numscalar("v2") \ st_numscalar("v3")))

* ---- case vif_dummy
use "../data/sim_het.dta", clear
quietly generate double w1 = 0.8 * x + 0.6 * z1
quietly regress y x w1 d
mata: pyx_export("../out/vif_dummy__x.csv", st_data(., ("x", "w1", "d")))
set type double
quietly estat vif, nosort
set type float
scalar v1 = r(vif_1)
scalar v2 = r(vif_2)
scalar v3 = r(vif_3)
mata: pyx_export("../out/vif_dummy__v.csv", (st_numscalar("v1") \ st_numscalar("v2") \ st_numscalar("v3")))

* ---- case vif_uncentered
use "../data/sim_het.dta", clear
quietly generate double w1 = 0.8 * x + 0.6 * z1
quietly generate double w2 = 0.5 * x - 0.5 * w1 + 0.3 * z2
quietly regress y x w1 w2
mata: pyx_export("../out/vif_uncentered__x.csv", st_data(., ("x", "w1", "w2")))
set type double
quietly estat vif, uncentered nosort
set type float
scalar v1 = r(vif_1)
scalar v2 = r(vif_2)
scalar v3 = r(vif_3)
scalar v4 = r(vif_4)
mata: pyx_export("../out/vif_uncentered__v.csv", (st_numscalar("v1") \ st_numscalar("v2") \ st_numscalar("v3") \ st_numscalar("v4")))

* ---- case vif_independent
use "../data/sim_iv_valid.dta", clear
quietly regress y x z1 z2 z3
mata: pyx_export("../out/vif_independent__x.csv", st_data(., ("x", "z1", "z2", "z3")))
set type double
quietly estat vif, nosort
set type float
scalar v1 = r(vif_1)
scalar v2 = r(vif_2)
scalar v3 = r(vif_3)
scalar v4 = r(vif_4)
mata: pyx_export("../out/vif_independent__v.csv", (st_numscalar("v1") \ st_numscalar("v2") \ st_numscalar("v3") \ st_numscalar("v4")))

exit, clear
