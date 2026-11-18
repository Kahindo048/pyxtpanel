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

* ---- case dh_lag1
use "../data/sim_granger_yes.dta", clear
quietly xtset id t
quietly summarize t
scalar tt = r(max) - r(min) + 1
quietly levelsof id
scalar nn = r(r)
matrix W = J(nn, 1, .)
local lx
forvalues j = 1/1 {
    local lx `lx' L`j'.x
}
forvalues i = 1/`=nn' {
    quietly regress y L(1/1).y L(1/1).x if id == `i'
    quietly test `lx'
    matrix W[`i', 1] = r(F) * 1
}
mata:
W = st_matrix("W")
N = rows(W)
K = 1
T = st_numscalar("tt")
wbar = mean(W)
zbar = sqrt(N / (2 * K)) * (wbar - K)
zt = sqrt(N / (2 * K) * (T - 3 * K - 5) / (T - 2 * K - 3)) * ((T - 3 * K - 3) / (T - 3 * K - 1) * wbar - K)
st_numscalar("wbar", wbar)
st_numscalar("zbar", zbar)
st_numscalar("pzbar", 2 * normal(-abs(zbar)))
st_numscalar("zt", zt)
st_numscalar("pzt", 2 * normal(-abs(zt)))
end
mata: pyx_export("../out/dh_lag1__unit_w.csv", st_matrix("W"))
mata: pyx_export("../out/dh_lag1__agg.csv", (st_numscalar("wbar") \ st_numscalar("zbar") \ st_numscalar("pzbar") \ st_numscalar("zt") \ st_numscalar("pzt") \ st_numscalar("nn") \ st_numscalar("tt")))

* ---- case dh_lag2
use "../data/sim_granger_yes.dta", clear
quietly xtset id t
quietly summarize t
scalar tt = r(max) - r(min) + 1
quietly levelsof id
scalar nn = r(r)
matrix W = J(nn, 1, .)
local lx
forvalues j = 1/2 {
    local lx `lx' L`j'.x
}
forvalues i = 1/`=nn' {
    quietly regress y L(1/2).y L(1/2).x if id == `i'
    quietly test `lx'
    matrix W[`i', 1] = r(F) * 2
}
mata:
W = st_matrix("W")
N = rows(W)
K = 2
T = st_numscalar("tt")
wbar = mean(W)
zbar = sqrt(N / (2 * K)) * (wbar - K)
zt = sqrt(N / (2 * K) * (T - 3 * K - 5) / (T - 2 * K - 3)) * ((T - 3 * K - 3) / (T - 3 * K - 1) * wbar - K)
st_numscalar("wbar", wbar)
st_numscalar("zbar", zbar)
st_numscalar("pzbar", 2 * normal(-abs(zbar)))
st_numscalar("zt", zt)
st_numscalar("pzt", 2 * normal(-abs(zt)))
end
mata: pyx_export("../out/dh_lag2__unit_w.csv", st_matrix("W"))
mata: pyx_export("../out/dh_lag2__agg.csv", (st_numscalar("wbar") \ st_numscalar("zbar") \ st_numscalar("pzbar") \ st_numscalar("zt") \ st_numscalar("pzt") \ st_numscalar("nn") \ st_numscalar("tt")))

* ---- case dh_lag3
use "../data/sim_granger_yes.dta", clear
quietly xtset id t
quietly summarize t
scalar tt = r(max) - r(min) + 1
quietly levelsof id
scalar nn = r(r)
matrix W = J(nn, 1, .)
local lx
forvalues j = 1/3 {
    local lx `lx' L`j'.x
}
forvalues i = 1/`=nn' {
    quietly regress y L(1/3).y L(1/3).x if id == `i'
    quietly test `lx'
    matrix W[`i', 1] = r(F) * 3
}
mata:
W = st_matrix("W")
N = rows(W)
K = 3
T = st_numscalar("tt")
wbar = mean(W)
zbar = sqrt(N / (2 * K)) * (wbar - K)
zt = sqrt(N / (2 * K) * (T - 3 * K - 5) / (T - 2 * K - 3)) * ((T - 3 * K - 3) / (T - 3 * K - 1) * wbar - K)
st_numscalar("wbar", wbar)
st_numscalar("zbar", zbar)
st_numscalar("pzbar", 2 * normal(-abs(zbar)))
st_numscalar("zt", zt)
st_numscalar("pzt", 2 * normal(-abs(zt)))
end
mata: pyx_export("../out/dh_lag3__unit_w.csv", st_matrix("W"))
mata: pyx_export("../out/dh_lag3__agg.csv", (st_numscalar("wbar") \ st_numscalar("zbar") \ st_numscalar("pzbar") \ st_numscalar("zt") \ st_numscalar("pzt") \ st_numscalar("nn") \ st_numscalar("tt")))

* ---- case dh_no_causality
use "../data/sim_granger_no.dta", clear
quietly xtset id t
quietly summarize t
scalar tt = r(max) - r(min) + 1
quietly levelsof id
scalar nn = r(r)
matrix W = J(nn, 1, .)
local lx
forvalues j = 1/1 {
    local lx `lx' L`j'.x
}
forvalues i = 1/`=nn' {
    quietly regress y L(1/1).y L(1/1).x if id == `i'
    quietly test `lx'
    matrix W[`i', 1] = r(F) * 1
}
mata:
W = st_matrix("W")
N = rows(W)
K = 1
T = st_numscalar("tt")
wbar = mean(W)
zbar = sqrt(N / (2 * K)) * (wbar - K)
zt = sqrt(N / (2 * K) * (T - 3 * K - 5) / (T - 2 * K - 3)) * ((T - 3 * K - 3) / (T - 3 * K - 1) * wbar - K)
st_numscalar("wbar", wbar)
st_numscalar("zbar", zbar)
st_numscalar("pzbar", 2 * normal(-abs(zbar)))
st_numscalar("zt", zt)
st_numscalar("pzt", 2 * normal(-abs(zt)))
end
mata: pyx_export("../out/dh_no_causality__unit_w.csv", st_matrix("W"))
mata: pyx_export("../out/dh_no_causality__agg.csv", (st_numscalar("wbar") \ st_numscalar("zbar") \ st_numscalar("pzbar") \ st_numscalar("zt") \ st_numscalar("pzt") \ st_numscalar("nn") \ st_numscalar("tt")))

* ---- case dh_no_causality_lag2
use "../data/sim_granger_no.dta", clear
quietly xtset id t
quietly summarize t
scalar tt = r(max) - r(min) + 1
quietly levelsof id
scalar nn = r(r)
matrix W = J(nn, 1, .)
local lx
forvalues j = 1/2 {
    local lx `lx' L`j'.x
}
forvalues i = 1/`=nn' {
    quietly regress y L(1/2).y L(1/2).x if id == `i'
    quietly test `lx'
    matrix W[`i', 1] = r(F) * 2
}
mata:
W = st_matrix("W")
N = rows(W)
K = 2
T = st_numscalar("tt")
wbar = mean(W)
zbar = sqrt(N / (2 * K)) * (wbar - K)
zt = sqrt(N / (2 * K) * (T - 3 * K - 5) / (T - 2 * K - 3)) * ((T - 3 * K - 3) / (T - 3 * K - 1) * wbar - K)
st_numscalar("wbar", wbar)
st_numscalar("zbar", zbar)
st_numscalar("pzbar", 2 * normal(-abs(zbar)))
st_numscalar("zt", zt)
st_numscalar("pzt", 2 * normal(-abs(zt)))
end
mata: pyx_export("../out/dh_no_causality_lag2__unit_w.csv", st_matrix("W"))
mata: pyx_export("../out/dh_no_causality_lag2__agg.csv", (st_numscalar("wbar") \ st_numscalar("zbar") \ st_numscalar("pzbar") \ st_numscalar("zt") \ st_numscalar("pzt") \ st_numscalar("nn") \ st_numscalar("tt")))

exit, clear
