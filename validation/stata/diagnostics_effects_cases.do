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

* ---- case bl_effect
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
quietly xtreg y x z1 z2, re
quietly xttest1, unadjusted
scalar jt = r(jt_lm)
scalar jtp = r(jt_lm_p)
scalar re = r(rets_lm)
scalar rea = r(rets_alm)
scalar reap = r(rets_alm_p)
scalar se = r(sc_lm)
scalar sep = r(sc_lm_p)
scalar sea = r(sc_alm)
scalar seap = r(sc_alm_p)
scalar reo = r(reos_lm)
scalar reop = r(reos_lm_p)
scalar reao = r(reos_alm)
scalar reaop = r(reos_alm_p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_effect__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_effect__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_effect__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_effect__r.csv", (st_numscalar("jt") \ st_numscalar("jtp") \ st_numscalar("re") \ st_numscalar("rea") \ st_numscalar("reap") \ st_numscalar("se") \ st_numscalar("sep") \ st_numscalar("sea") \ st_numscalar("seap") \ st_numscalar("reo") \ st_numscalar("reao")))

* ---- case bl_mata_effect
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly regress y x z1 z2
quietly predict double ehat, resid
mata:
e = st_data(., "ehat")
id = st_data(., "id")
info = panelsetup(id, 1)
G = rows(info)
m = rows(e)
S = 0
a = 0
num = 0
den = 0
for (g = 1; g <= G; g++) {
    ei = panelsubmatrix(e, g, info)
    Ti = rows(ei)
    S = S + sum(ei)^2
    a = a + Ti^2
    if (Ti > 1) {
        num = num + sum(ei[|2 \ Ti|] :* ei[|1 \ Ti - 1|])
        den = den + sum(ei[|2 \ Ti|] :^ 2)
    }
}
A = 1 - S / sum(e :^ 2)
B = num / den
lre = m^2 * A^2 / (2 * (a - m))
lse = m^2 * B^2 / (m - G)
lrea = m^2 * (A + 2 * B)^2 / (2 * (a - 3 * m + 2 * G))
lsea = lse - lre + lrea
st_numscalar("jt", lre + lsea)
st_numscalar("re", lre)
st_numscalar("rea", lrea)
st_numscalar("se", lse)
st_numscalar("sea", lsea)
st_numscalar("reo", -sqrt(m^2 / (2 * (a - m))) * A)
st_numscalar("reao", -sqrt(m^2 / (2 * (a - 3 * m + 2 * G))) * (A + 2 * B))
end
scalar jtp = chi2tail(2, jt)
scalar reap = chi2tail(1, rea)
scalar sep = chi2tail(1, se)
scalar seap = chi2tail(1, sea)
scalar reop = normal(-reo)
scalar reaop = normal(-reao)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_mata_effect__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_mata_effect__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_mata_effect__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_mata_effect__r.csv", (st_numscalar("jt") \ st_numscalar("jtp") \ st_numscalar("re") \ st_numscalar("rea") \ st_numscalar("reap") \ st_numscalar("se") \ st_numscalar("sep") \ st_numscalar("sea") \ st_numscalar("seap") \ st_numscalar("reo") \ st_numscalar("reop") \ st_numscalar("reao") \ st_numscalar("reaop")))

* ---- case bl_unbalanced
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly regress x t
quietly predict double ehat, resid
quietly xtreg x t, re
quietly xttest1, unadjusted
scalar jt = r(jt_lm)
scalar jtp = r(jt_lm_p)
scalar re = r(rets_lm)
scalar rea = r(rets_alm)
scalar reap = r(rets_alm_p)
scalar se = r(sc_lm)
scalar sep = r(sc_lm_p)
scalar sea = r(sc_alm)
scalar seap = r(sc_alm_p)
scalar reo = r(reos_lm)
scalar reop = r(reos_lm_p)
scalar reao = r(reos_alm)
scalar reaop = r(reos_alm_p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_unbalanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_unbalanced__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_unbalanced__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_unbalanced__r.csv", (st_numscalar("jt") \ st_numscalar("jtp") \ st_numscalar("re") \ st_numscalar("rea") \ st_numscalar("reap") \ st_numscalar("se") \ st_numscalar("sep") \ st_numscalar("sea") \ st_numscalar("seap") \ st_numscalar("reo") \ st_numscalar("reao")))

* ---- case bl_mata_unbalanced
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly regress x t
quietly predict double ehat, resid
mata:
e = st_data(., "ehat")
id = st_data(., "id")
info = panelsetup(id, 1)
G = rows(info)
m = rows(e)
S = 0
a = 0
num = 0
den = 0
for (g = 1; g <= G; g++) {
    ei = panelsubmatrix(e, g, info)
    Ti = rows(ei)
    S = S + sum(ei)^2
    a = a + Ti^2
    if (Ti > 1) {
        num = num + sum(ei[|2 \ Ti|] :* ei[|1 \ Ti - 1|])
        den = den + sum(ei[|2 \ Ti|] :^ 2)
    }
}
A = 1 - S / sum(e :^ 2)
B = num / den
lre = m^2 * A^2 / (2 * (a - m))
lse = m^2 * B^2 / (m - G)
lrea = m^2 * (A + 2 * B)^2 / (2 * (a - 3 * m + 2 * G))
lsea = lse - lre + lrea
st_numscalar("jt", lre + lsea)
st_numscalar("re", lre)
st_numscalar("rea", lrea)
st_numscalar("se", lse)
st_numscalar("sea", lsea)
st_numscalar("reo", -sqrt(m^2 / (2 * (a - m))) * A)
st_numscalar("reao", -sqrt(m^2 / (2 * (a - 3 * m + 2 * G))) * (A + 2 * B))
end
scalar jtp = chi2tail(2, jt)
scalar reap = chi2tail(1, rea)
scalar sep = chi2tail(1, se)
scalar seap = chi2tail(1, sea)
scalar reop = normal(-reo)
scalar reaop = normal(-reao)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_mata_unbalanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_mata_unbalanced__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_mata_unbalanced__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_mata_unbalanced__r.csv", (st_numscalar("jt") \ st_numscalar("jtp") \ st_numscalar("re") \ st_numscalar("rea") \ st_numscalar("reap") \ st_numscalar("se") \ st_numscalar("sep") \ st_numscalar("sea") \ st_numscalar("seap") \ st_numscalar("reo") \ st_numscalar("reop") \ st_numscalar("reao") \ st_numscalar("reaop")))

* ---- case bl_no_effect
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly regress x t
quietly predict double ehat, resid
quietly xtreg x t, re
quietly xttest1, unadjusted
scalar jt = r(jt_lm)
scalar jtp = r(jt_lm_p)
scalar re = r(rets_lm)
scalar rea = r(rets_alm)
scalar reap = r(rets_alm_p)
scalar se = r(sc_lm)
scalar sep = r(sc_lm_p)
scalar sea = r(sc_alm)
scalar seap = r(sc_alm_p)
scalar reo = r(reos_lm)
scalar reop = r(reos_lm_p)
scalar reao = r(reos_alm)
scalar reaop = r(reos_alm_p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_no_effect__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_no_effect__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_no_effect__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_no_effect__r.csv", (st_numscalar("jt") \ st_numscalar("jtp") \ st_numscalar("re") \ st_numscalar("rea") \ st_numscalar("reap") \ st_numscalar("se") \ st_numscalar("sep") \ st_numscalar("sea") \ st_numscalar("seap") \ st_numscalar("reo") \ st_numscalar("reao")))

* ---- case bl_mata_no_effect
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly regress x t
quietly predict double ehat, resid
mata:
e = st_data(., "ehat")
id = st_data(., "id")
info = panelsetup(id, 1)
G = rows(info)
m = rows(e)
S = 0
a = 0
num = 0
den = 0
for (g = 1; g <= G; g++) {
    ei = panelsubmatrix(e, g, info)
    Ti = rows(ei)
    S = S + sum(ei)^2
    a = a + Ti^2
    if (Ti > 1) {
        num = num + sum(ei[|2 \ Ti|] :* ei[|1 \ Ti - 1|])
        den = den + sum(ei[|2 \ Ti|] :^ 2)
    }
}
A = 1 - S / sum(e :^ 2)
B = num / den
lre = m^2 * A^2 / (2 * (a - m))
lse = m^2 * B^2 / (m - G)
lrea = m^2 * (A + 2 * B)^2 / (2 * (a - 3 * m + 2 * G))
lsea = lse - lre + lrea
st_numscalar("jt", lre + lsea)
st_numscalar("re", lre)
st_numscalar("rea", lrea)
st_numscalar("se", lse)
st_numscalar("sea", lsea)
st_numscalar("reo", -sqrt(m^2 / (2 * (a - m))) * A)
st_numscalar("reao", -sqrt(m^2 / (2 * (a - 3 * m + 2 * G))) * (A + 2 * B))
end
scalar jtp = chi2tail(2, jt)
scalar reap = chi2tail(1, rea)
scalar sep = chi2tail(1, se)
scalar seap = chi2tail(1, sea)
scalar reop = normal(-reo)
scalar reaop = normal(-reao)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_mata_no_effect__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_mata_no_effect__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_mata_no_effect__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_mata_no_effect__r.csv", (st_numscalar("jt") \ st_numscalar("jtp") \ st_numscalar("re") \ st_numscalar("rea") \ st_numscalar("reap") \ st_numscalar("se") \ st_numscalar("sep") \ st_numscalar("sea") \ st_numscalar("seap") \ st_numscalar("reo") \ st_numscalar("reop") \ st_numscalar("reao") \ st_numscalar("reaop")))

* ---- case bl_serial
use "../data/sim_ur_stationary.dta", clear
quietly xtset id t
quietly regress y t
quietly predict double ehat, resid
quietly xtreg y t, re
quietly xttest1, unadjusted
scalar jt = r(jt_lm)
scalar jtp = r(jt_lm_p)
scalar re = r(rets_lm)
scalar rea = r(rets_alm)
scalar reap = r(rets_alm_p)
scalar se = r(sc_lm)
scalar sep = r(sc_lm_p)
scalar sea = r(sc_alm)
scalar seap = r(sc_alm_p)
scalar reo = r(reos_lm)
scalar reop = r(reos_lm_p)
scalar reao = r(reos_alm)
scalar reaop = r(reos_alm_p)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_serial__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_serial__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_serial__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_serial__r.csv", (st_numscalar("jt") \ st_numscalar("jtp") \ st_numscalar("re") \ st_numscalar("rea") \ st_numscalar("reap") \ st_numscalar("se") \ st_numscalar("sep") \ st_numscalar("sea") \ st_numscalar("seap") \ st_numscalar("reo") \ st_numscalar("reao")))

* ---- case bl_mata_serial
use "../data/sim_ur_stationary.dta", clear
quietly xtset id t
quietly regress y t
quietly predict double ehat, resid
mata:
e = st_data(., "ehat")
id = st_data(., "id")
info = panelsetup(id, 1)
G = rows(info)
m = rows(e)
S = 0
a = 0
num = 0
den = 0
for (g = 1; g <= G; g++) {
    ei = panelsubmatrix(e, g, info)
    Ti = rows(ei)
    S = S + sum(ei)^2
    a = a + Ti^2
    if (Ti > 1) {
        num = num + sum(ei[|2 \ Ti|] :* ei[|1 \ Ti - 1|])
        den = den + sum(ei[|2 \ Ti|] :^ 2)
    }
}
A = 1 - S / sum(e :^ 2)
B = num / den
lre = m^2 * A^2 / (2 * (a - m))
lse = m^2 * B^2 / (m - G)
lrea = m^2 * (A + 2 * B)^2 / (2 * (a - 3 * m + 2 * G))
lsea = lse - lre + lrea
st_numscalar("jt", lre + lsea)
st_numscalar("re", lre)
st_numscalar("rea", lrea)
st_numscalar("se", lse)
st_numscalar("sea", lsea)
st_numscalar("reo", -sqrt(m^2 / (2 * (a - m))) * A)
st_numscalar("reao", -sqrt(m^2 / (2 * (a - 3 * m + 2 * G))) * (A + 2 * B))
end
scalar jtp = chi2tail(2, jt)
scalar reap = chi2tail(1, rea)
scalar sep = chi2tail(1, se)
scalar seap = chi2tail(1, sea)
scalar reop = normal(-reo)
scalar reaop = normal(-reao)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_mata_serial__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_mata_serial__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_mata_serial__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_mata_serial__r.csv", (st_numscalar("jt") \ st_numscalar("jtp") \ st_numscalar("re") \ st_numscalar("rea") \ st_numscalar("reap") \ st_numscalar("se") \ st_numscalar("sep") \ st_numscalar("sea") \ st_numscalar("seap") \ st_numscalar("reo") \ st_numscalar("reop") \ st_numscalar("reao") \ st_numscalar("reaop")))

* ---- case bl_gap_refused
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly regress x z
quietly predict double ehat, resid
quietly xtreg x z, re
capture quietly xttest1, unadjusted
scalar rc = _rc
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/bl_gap_refused__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/bl_gap_refused__id.csv", st_data(., "id"))
mata: pyx_export("../out/bl_gap_refused__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/bl_gap_refused__rc.csv", (st_numscalar("rc")))

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
