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

* ---- case cc_balanced
use "../data/sim_csd.dta", clear
quietly xtset id t
preserve
quietly keep if !missing(s)
quietly keep id t s
quietly reshape wide s, i(t) j(id)
mata:
X = st_data(., 2..st_nvar())
n = cols(X)
R = J(n, n, .)
N = J(n, n, .)
total = 0
totalabs = 0
count = 0
for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
        ok = (X[., i] :< .) :& (X[., j] :< .)
        a = select(X[., i], ok)
        b = select(X[., j], ok)
        N[i, j] = N[j, i] = rows(a)
        if (i == j) {
            R[i, j] = 1
        }
        else if (rows(a) >= 3) {
            r = quadcorrelation((a, b))[1, 2]
            R[i, j] = R[j, i] = r
            if (r < .) {
                total = total + r
                totalabs = totalabs + abs(r)
                count = count + 1
            }
        }
    }
}
st_numscalar("ccavg", total / count)
st_numscalar("ccabs", totalabs / count)
st_numscalar("ccn", count)
end
mata: pyx_export("../out/cc_balanced__rho.csv", R)
mata: pyx_export("../out/cc_balanced__overlap.csv", N)
restore
mata: pyx_export("../out/cc_balanced__avg.csv", (st_numscalar("ccavg") \ st_numscalar("ccabs") \ st_numscalar("ccn")))

* ---- case cc_unbalanced
use "../data/sim_csd_unbalanced.dta", clear
quietly xtset id t
preserve
quietly keep if !missing(s)
quietly keep id t s
quietly reshape wide s, i(t) j(id)
mata:
X = st_data(., 2..st_nvar())
n = cols(X)
R = J(n, n, .)
N = J(n, n, .)
total = 0
totalabs = 0
count = 0
for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
        ok = (X[., i] :< .) :& (X[., j] :< .)
        a = select(X[., i], ok)
        b = select(X[., j], ok)
        N[i, j] = N[j, i] = rows(a)
        if (i == j) {
            R[i, j] = 1
        }
        else if (rows(a) >= 3) {
            r = quadcorrelation((a, b))[1, 2]
            R[i, j] = R[j, i] = r
            if (r < .) {
                total = total + r
                totalabs = totalabs + abs(r)
                count = count + 1
            }
        }
    }
}
st_numscalar("ccavg", total / count)
st_numscalar("ccabs", totalabs / count)
st_numscalar("ccn", count)
end
mata: pyx_export("../out/cc_unbalanced__rho.csv", R)
mata: pyx_export("../out/cc_unbalanced__overlap.csv", N)
restore
mata: pyx_export("../out/cc_unbalanced__avg.csv", (st_numscalar("ccavg") \ st_numscalar("ccabs") \ st_numscalar("ccn")))

* ---- case cc_sparse
use "../data/sim_sparse_t.dta", clear
quietly xtset id t
preserve
quietly keep if !missing(x)
quietly keep id t x
quietly reshape wide x, i(t) j(id)
mata:
X = st_data(., 2..st_nvar())
n = cols(X)
R = J(n, n, .)
N = J(n, n, .)
total = 0
totalabs = 0
count = 0
for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
        ok = (X[., i] :< .) :& (X[., j] :< .)
        a = select(X[., i], ok)
        b = select(X[., j], ok)
        N[i, j] = N[j, i] = rows(a)
        if (i == j) {
            R[i, j] = 1
        }
        else if (rows(a) >= 3) {
            r = quadcorrelation((a, b))[1, 2]
            R[i, j] = R[j, i] = r
            if (r < .) {
                total = total + r
                totalabs = totalabs + abs(r)
                count = count + 1
            }
        }
    }
}
st_numscalar("ccavg", total / count)
st_numscalar("ccabs", totalabs / count)
st_numscalar("ccn", count)
end
mata: pyx_export("../out/cc_sparse__rho.csv", R)
mata: pyx_export("../out/cc_sparse__overlap.csv", N)
restore
mata: pyx_export("../out/cc_sparse__avg.csv", (st_numscalar("ccavg") \ st_numscalar("ccabs") \ st_numscalar("ccn")))

* ---- case frees_balanced
use "../data/sim_csd.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, frees
scalar fr = r(frees)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/frees_balanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/frees_balanced__id.csv", st_data(., "id"))
mata: pyx_export("../out/frees_balanced__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/frees_balanced__r.csv", (st_numscalar("fr")))

* ---- case frees_trimmed
use "../data/sim_csd_trimmed.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, frees
scalar fr = r(frees)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/frees_trimmed__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/frees_trimmed__id.csv", st_data(., "id"))
mata: pyx_export("../out/frees_trimmed__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/frees_trimmed__r.csv", (st_numscalar("fr")))

* ---- case frees_long
use "../data/sim_csd_long.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, frees
scalar fr = r(frees)
quietly summarize t
scalar mm = r(max) - r(min) + 1
scalar seq1 = (32/25) * ((mm + 2)^2 / ((mm - 1)^3 * (mm + 1)^2))
scalar seq2 = (4/25) * (((5*mm + 6)^2 * (mm - 3)) / ((mm * (mm - 1)^2) * (mm + 1)^2))
scalar frp = 2 * (1 - normal(fr / sqrt(seq1 + seq2)))
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/frees_long__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/frees_long__id.csv", st_data(., "id"))
mata: pyx_export("../out/frees_long__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/frees_long__r.csv", (st_numscalar("fr") \ st_numscalar("frp")))

* ---- case fri_balanced
use "../data/sim_csd.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, friedman
scalar fr = r(fried)
scalar frdf = e(N_g) - 1
scalar frp = chi2tail(frdf, fr)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/fri_balanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/fri_balanced__id.csv", st_data(., "id"))
mata: pyx_export("../out/fri_balanced__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/fri_balanced__r.csv", (st_numscalar("fr") \ st_numscalar("frdf") \ st_numscalar("frp")))

* ---- case fri_trimmed
use "../data/sim_csd_trimmed.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, friedman
scalar fr = r(fried)
scalar frdf = e(N_g) - 1
scalar frp = chi2tail(frdf, fr)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/fri_trimmed__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/fri_trimmed__id.csv", st_data(., "id"))
mata: pyx_export("../out/fri_trimmed__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/fri_trimmed__r.csv", (st_numscalar("fr") \ st_numscalar("frdf") \ st_numscalar("frp")))

* ---- case fri_periods_balanced
use "../data/sim_csd.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, friedman
scalar fr = r(fried)
scalar frdf = e(g_min) - 1
scalar frp = chi2tail(frdf, fr)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/fri_periods_balanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/fri_periods_balanced__id.csv", st_data(., "id"))
mata: pyx_export("../out/fri_periods_balanced__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/fri_periods_balanced__r.csv", (st_numscalar("fr") \ st_numscalar("frdf") \ st_numscalar("frp")))

* ---- case fri_periods_trimmed
use "../data/sim_csd_trimmed.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, friedman
scalar fr = r(fried)
scalar frdf = e(g_min) - 1
scalar frp = chi2tail(frdf, fr)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/fri_periods_trimmed__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/fri_periods_trimmed__id.csv", st_data(., "id"))
mata: pyx_export("../out/fri_periods_trimmed__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/fri_periods_trimmed__r.csv", (st_numscalar("fr") \ st_numscalar("frdf") \ st_numscalar("frp")))

* ---- case lmcs_factor
use "../data/sim_csd.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
preserve
quietly keep if !missing(ehat)
quietly keep id t ehat
quietly reshape wide ehat, i(t) j(id)
mata:
X = st_data(., 2..st_nvar())
n = cols(X)
lm = 0
for (i = 1; i < n; i++) {
    for (j = i + 1; j <= n; j++) {
        ok = (X[., i] :< .) :& (X[., j] :< .)
        a = select(X[., i], ok)
        b = select(X[., j], ok)
        r = quadcorrelation((a, b))[1, 2]
        lm = lm + rows(a) * r^2
    }
}
st_numscalar("lmv", lm)
st_numscalar("lmdf", n * (n - 1) / 2)
end
restore
scalar lmp = chi2tail(lmdf, lmv)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/lmcs_factor__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/lmcs_factor__id.csv", st_data(., "id"))
mata: pyx_export("../out/lmcs_factor__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/lmcs_factor__r.csv", (st_numscalar("lmv") \ st_numscalar("lmdf") \ st_numscalar("lmp")))

* ---- case lmcs_unbalanced
use "../data/sim_csd_unbalanced.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
preserve
quietly keep if !missing(ehat)
quietly keep id t ehat
quietly reshape wide ehat, i(t) j(id)
mata:
X = st_data(., 2..st_nvar())
n = cols(X)
lm = 0
for (i = 1; i < n; i++) {
    for (j = i + 1; j <= n; j++) {
        ok = (X[., i] :< .) :& (X[., j] :< .)
        a = select(X[., i], ok)
        b = select(X[., j], ok)
        r = quadcorrelation((a, b))[1, 2]
        lm = lm + rows(a) * r^2
    }
}
st_numscalar("lmv", lm)
st_numscalar("lmdf", n * (n - 1) / 2)
end
restore
scalar lmp = chi2tail(lmdf, lmv)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/lmcs_unbalanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/lmcs_unbalanced__id.csv", st_data(., "id"))
mata: pyx_export("../out/lmcs_unbalanced__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/lmcs_unbalanced__r.csv", (st_numscalar("lmv") \ st_numscalar("lmdf") \ st_numscalar("lmp")))

* ---- case lmcs_independent
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtreg x, fe
quietly predict double ehat, e
preserve
quietly keep if !missing(ehat)
quietly keep id t ehat
quietly reshape wide ehat, i(t) j(id)
mata:
X = st_data(., 2..st_nvar())
n = cols(X)
lm = 0
for (i = 1; i < n; i++) {
    for (j = i + 1; j <= n; j++) {
        ok = (X[., i] :< .) :& (X[., j] :< .)
        a = select(X[., i], ok)
        b = select(X[., j], ok)
        r = quadcorrelation((a, b))[1, 2]
        lm = lm + rows(a) * r^2
    }
}
st_numscalar("lmv", lm)
st_numscalar("lmdf", n * (n - 1) / 2)
end
restore
scalar lmp = chi2tail(lmdf, lmv)
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/lmcs_independent__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/lmcs_independent__id.csv", st_data(., "id"))
mata: pyx_export("../out/lmcs_independent__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/lmcs_independent__r.csv", (st_numscalar("lmv") \ st_numscalar("lmdf") \ st_numscalar("lmp")))

* ---- case cd_balanced
use "../data/sim_csd.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, pesaran
scalar cd = r(pesaran)
scalar cdp = 2 * normal(-abs(cd))
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/cd_balanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/cd_balanced__id.csv", st_data(., "id"))
mata: pyx_export("../out/cd_balanced__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/cd_balanced__cd.csv", (st_numscalar("cd") \ st_numscalar("cdp")))

* ---- case cd_unbalanced
use "../data/sim_csd_unbalanced.dta", clear
quietly xtset id t
quietly xtreg y x, fe
quietly predict double ehat, e
quietly xtcsd, pesaran
scalar cd = r(pesaran)
scalar cdp = 2 * normal(-abs(cd))
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/cd_unbalanced__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/cd_unbalanced__id.csv", st_data(., "id"))
mata: pyx_export("../out/cd_unbalanced__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/cd_unbalanced__cd.csv", (st_numscalar("cd") \ st_numscalar("cdp")))

* ---- case cd_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly xtreg x, fe
quietly predict double ehat, e
quietly xtcsd, pesaran
scalar cd = r(pesaran)
scalar cdp = 2 * normal(-abs(cd))
preserve
quietly keep if !missing(ehat)
mata: pyx_export("../out/cd_gaps__e.csv", st_data(., "ehat"))
mata: pyx_export("../out/cd_gaps__id.csv", st_data(., "id"))
mata: pyx_export("../out/cd_gaps__t.csv", st_data(., "t"))
restore
mata: pyx_export("../out/cd_gaps__cd.csv", (st_numscalar("cd") \ st_numscalar("cdp")))

* ---- case cd_variable
use "../data/sim_csd.dta", clear
quietly xtset id t
quietly xtcd2 s, noest
matrix cdm = r(CD)
scalar cd = cdm[1, 1]
scalar cdp = 2 * normal(-abs(cd))
mata: pyx_export("../out/cd_variable__cd.csv", (st_numscalar("cd") \ st_numscalar("cdp")))

* ---- case cd_independent
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtcd2 x, noest
matrix cdm = r(CD)
scalar cd = cdm[1, 1]
scalar cdp = 2 * normal(-abs(cd))
mata: pyx_export("../out/cd_independent__cd.csv", (st_numscalar("cd") \ st_numscalar("cdp")))

exit, clear
