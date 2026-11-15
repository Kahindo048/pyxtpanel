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

* ---- case jb_heavy_tails
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/jb_heavy_tails__e.csv", st_data(., "ehat"))
quietly summarize ehat, detail
scalar jb = r(N) / 6 * (r(skewness)^2 + (r(kurtosis) - 3)^2 / 4)
scalar jbp = chi2tail(2, jb)
mata: pyx_export("../out/jb_heavy_tails__jb.csv", (st_numscalar("jb") \ st_numscalar("jbp")))

* ---- case jb_normal
use "../data/sim_iv_valid.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/jb_normal__e.csv", st_data(., "ehat"))
quietly summarize ehat, detail
scalar jb = r(N) / 6 * (r(skewness)^2 + (r(kurtosis) - 3)^2 / 4)
scalar jbp = chi2tail(2, jb)
mata: pyx_export("../out/jb_normal__jb.csv", (st_numscalar("jb") \ st_numscalar("jbp")))

* ---- case jb_small
use "../data/sim_balanced.dta", clear
quietly keep in 1/15
quietly regress x z
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/jb_small__e.csv", st_data(., "ehat"))
quietly summarize ehat, detail
scalar jb = r(N) / 6 * (r(skewness)^2 + (r(kurtosis) - 3)^2 / 4)
scalar jbp = chi2tail(2, jb)
mata: pyx_export("../out/jb_small__jb.csv", (st_numscalar("jb") \ st_numscalar("jbp")))

* ---- case sf_heavy_tails
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sf_heavy_tails__e.csv", st_data(., "ehat"))
set type double
quietly sfrancia ehat
scalar fw = r(W)
scalar fv = r(V)
scalar fz = r(z)
scalar fp = r(p)
scalar fn = r(N)
mata: pyx_export("../out/sf_heavy_tails__sf.csv", (st_numscalar("fw") \ st_numscalar("fv") \ st_numscalar("fz") \ st_numscalar("fp") \ st_numscalar("fn")))
set type float

* ---- case sf_single_precision
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sf_single_precision__e.csv", st_data(., "ehat"))
set type float
quietly sfrancia ehat
scalar fw = r(W)
scalar fv = r(V)
scalar fz = r(z)
scalar fp = r(p)
scalar fn = r(N)
mata: pyx_export("../out/sf_single_precision__sf.csv", (st_numscalar("fw") \ st_numscalar("fv") \ st_numscalar("fz") \ st_numscalar("fp") \ st_numscalar("fn")))
set type float

* ---- case sf_normal
use "../data/sim_iv_valid.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sf_normal__e.csv", st_data(., "ehat"))
set type double
quietly sfrancia ehat
scalar fw = r(W)
scalar fv = r(V)
scalar fz = r(z)
scalar fp = r(p)
scalar fn = r(N)
mata: pyx_export("../out/sf_normal__sf.csv", (st_numscalar("fw") \ st_numscalar("fv") \ st_numscalar("fz") \ st_numscalar("fp") \ st_numscalar("fn")))
set type float

* ---- case sf_boxcox
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sf_boxcox__e.csv", st_data(., "ehat"))
set type double
quietly sfrancia ehat, boxcox
scalar fw = r(W)
scalar fv = r(V)
scalar fz = r(z)
scalar fp = r(p)
scalar fn = r(N)
mata: pyx_export("../out/sf_boxcox__sf.csv", (st_numscalar("fw") \ st_numscalar("fv") \ st_numscalar("fz") \ st_numscalar("fp") \ st_numscalar("fn")))
set type float

* ---- case sf_boxcox_small
use "../data/sim_balanced.dta", clear
quietly keep in 1/8
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sf_boxcox_small__e.csv", st_data(., "ehat"))
set type double
quietly sfrancia ehat, boxcox
scalar fw = r(W)
scalar fv = r(V)
scalar fz = r(z)
scalar fp = r(p)
scalar fn = r(N)
mata: pyx_export("../out/sf_boxcox_small__sf.csv", (st_numscalar("fw") \ st_numscalar("fv") \ st_numscalar("fz") \ st_numscalar("fp") \ st_numscalar("fn")))
set type float

* ---- case sf_small
use "../data/sim_balanced.dta", clear
quietly keep in 1/12
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sf_small__e.csv", st_data(., "ehat"))
set type double
quietly sfrancia ehat
scalar fw = r(W)
scalar fv = r(V)
scalar fz = r(z)
scalar fp = r(p)
scalar fn = r(N)
mata: pyx_export("../out/sf_small__sf.csv", (st_numscalar("fw") \ st_numscalar("fv") \ st_numscalar("fz") \ st_numscalar("fp") \ st_numscalar("fn")))
set type float

* ---- case sf_ties
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly replace ehat = round(ehat, 0.5)
quietly keep if !missing(ehat)
mata: pyx_export("../out/sf_ties__e.csv", st_data(., "ehat"))
set type double
quietly sfrancia ehat
scalar fw = r(W)
scalar fv = r(V)
scalar fz = r(z)
scalar fp = r(p)
scalar fn = r(N)
mata: pyx_export("../out/sf_ties__sf.csv", (st_numscalar("fw") \ st_numscalar("fv") \ st_numscalar("fz") \ st_numscalar("fp") \ st_numscalar("fn")))
set type float

* ---- case sf_noties
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly replace ehat = round(ehat, 0.5)
quietly keep if !missing(ehat)
mata: pyx_export("../out/sf_noties__e.csv", st_data(., "ehat"))
set type double
quietly sfrancia ehat, noties
scalar fw = r(W)
scalar fv = r(V)
scalar fz = r(z)
scalar fp = r(p)
scalar fn = r(N)
mata: pyx_export("../out/sf_noties__sf.csv", (st_numscalar("fw") \ st_numscalar("fv") \ st_numscalar("fz") \ st_numscalar("fp") \ st_numscalar("fn")))
set type float

* ---- case sw_heavy_tails
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_heavy_tails__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_heavy_tails__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_single_precision
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_single_precision__e.csv", st_data(., "ehat"))
set type float
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_single_precision__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_normal
use "../data/sim_iv_valid.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_normal__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_normal__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_small
use "../data/sim_balanced.dta", clear
quietly keep in 1/15
quietly regress x z
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_small__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_small__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_n3
use "../data/sim_balanced.dta", clear
quietly keep in 1/3
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_n3__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_n3__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_n4
use "../data/sim_balanced.dta", clear
quietly keep in 1/4
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_n4__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_n4__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_n5
use "../data/sim_balanced.dta", clear
quietly keep in 1/5
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_n5__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_n5__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_n7
use "../data/sim_balanced.dta", clear
quietly keep in 1/7
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_n7__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_n7__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_n11
use "../data/sim_balanced.dta", clear
quietly keep in 1/11
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_n11__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_n11__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_n12
use "../data/sim_balanced.dta", clear
quietly keep in 1/12
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_n12__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_n12__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_ties
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly replace ehat = round(ehat, 0.5)
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_ties__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_ties__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_noties
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly replace ehat = round(ehat, 0.5)
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_noties__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat, noties
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_noties__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sw_large
use "../data/sim_wide.dta", clear
quietly generate double ehat = x
quietly keep if !missing(ehat)
mata: pyx_export("../out/sw_large__e.csv", st_data(., "ehat"))
set type double
quietly swilk ehat
scalar sw = r(W)
scalar sv = r(V)
scalar sz = r(z)
scalar sp = r(p)
scalar sn = r(N)
mata: pyx_export("../out/sw_large__sw.csv", (st_numscalar("sw") \ st_numscalar("sv") \ st_numscalar("sz") \ st_numscalar("sp") \ st_numscalar("sn")))
set type float
set type float

* ---- case sk_heavy_tails
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sk_heavy_tails__e.csv", st_data(., "ehat"))
quietly sktest ehat
scalar k2 = r(chi2)
scalar kp = r(p_chi2)
scalar ps = r(p_skew)
scalar pk = r(p_kurt)
mata: pyx_export("../out/sk_heavy_tails__sk.csv", (st_numscalar("k2") \ st_numscalar("kp") \ st_numscalar("ps") \ st_numscalar("pk")))

* ---- case sk_normal
use "../data/sim_iv_valid.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sk_normal__e.csv", st_data(., "ehat"))
quietly sktest ehat
scalar k2 = r(chi2)
scalar kp = r(p_chi2)
scalar ps = r(p_skew)
scalar pk = r(p_kurt)
mata: pyx_export("../out/sk_normal__sk.csv", (st_numscalar("k2") \ st_numscalar("kp") \ st_numscalar("ps") \ st_numscalar("pk")))

* ---- case sk_unadjusted
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sk_unadjusted__e.csv", st_data(., "ehat"))
quietly sktest ehat, noadjust
scalar k2 = r(chi2)
scalar kp = r(p_chi2)
scalar ps = r(p_skew)
scalar pk = r(p_kurt)
mata: pyx_export("../out/sk_unadjusted__sk.csv", (st_numscalar("k2") \ st_numscalar("kp") \ st_numscalar("ps") \ st_numscalar("pk")))

* ---- case sk_small
use "../data/sim_balanced.dta", clear
quietly keep in 1/15
quietly regress x z
quietly predict double ehat if e(sample), resid
quietly keep if !missing(ehat)
mata: pyx_export("../out/sk_small__e.csv", st_data(., "ehat"))
quietly sktest ehat
scalar k2 = r(chi2)
scalar kp = r(p_chi2)
scalar ps = r(p_skew)
scalar pk = r(p_kurt)
mata: pyx_export("../out/sk_small__sk.csv", (st_numscalar("k2") \ st_numscalar("kp") \ st_numscalar("ps") \ st_numscalar("pk")))

exit, clear
