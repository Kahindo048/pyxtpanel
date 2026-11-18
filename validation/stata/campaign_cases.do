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

mata:
real rowvector pyx_triple(string scalar v)
{
    real colvector x, k, keep
    x = st_data(., v); k = (1::rows(x)); keep = (x :< .)
    return((sum(keep), sum(select(x, keep)), sum(select(x :* k, keep))))
}
real rowvector pyx_pair(string scalar v)
{
    real colvector x, k, keep
    x = st_data(., v); k = (1::rows(x)); keep = (x :< .)
    return((sum(select(x :^ 2, keep)), sum(select(x :* k, keep)), sum(select(x, keep))))
}
end

capture program drop pyx_battery
program define pyx_battery
    args k delta fmt
    preserve
    quietly keep if scenario == `k'
    matrix row = J(1, 73, .)
    if "`fmt'" != "" format t `fmt'
    capture xtset id t, delta(`delta')
    matrix row[1, 1] = _rc
    if _rc == 0 {
    matrix row[1, 2] = r(tmin)
    matrix row[1, 3] = r(tmax)
    matrix row[1, 4] = r(tdelta)
    matrix row[1, 5] = r(gaps)
    matrix row[1, 6] = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
        scalar nobs_ = _N
        tempvar f1 f2
        bysort id: gen byte `f1' = (_n == 1)
        quietly summarize `f1'
        scalar nunits_ = r(sum)
        bysort t: gen byte `f2' = (_n == 1)
        quietly summarize `f2'
        scalar nper_ = r(sum)
        sort id t
    matrix row[1, 7] = nobs_
    matrix row[1, 8] = nunits_
    matrix row[1, 9] = nper_
        capture quietly xtsum x
        if _rc == 0 {
        matrix row[1, 10] = r(N)
        matrix row[1, 11] = r(n)
        matrix row[1, 12] = r(Tbar)
        matrix row[1, 13] = r(mean)
        matrix row[1, 14] = r(sd)
        matrix row[1, 15] = r(min)
        matrix row[1, 16] = r(max)
        matrix row[1, 17] = r(sd_b)
        matrix row[1, 18] = r(min_b)
        matrix row[1, 19] = r(max_b)
        matrix row[1, 20] = r(sd_w)
        matrix row[1, 21] = r(min_w)
        matrix row[1, 22] = r(max_w)
        }
        capture quietly xtsum z
        if _rc == 0 {
        matrix row[1, 23] = r(N)
        matrix row[1, 24] = r(n)
        matrix row[1, 25] = r(Tbar)
        matrix row[1, 26] = r(mean)
        matrix row[1, 27] = r(sd)
        matrix row[1, 28] = r(min)
        matrix row[1, 29] = r(max)
        matrix row[1, 30] = r(sd_b)
        matrix row[1, 31] = r(min_b)
        matrix row[1, 32] = r(max_b)
        matrix row[1, 33] = r(sd_w)
        matrix row[1, 34] = r(min_w)
        matrix row[1, 35] = r(max_w)
        }
        tempvar o_L
        quietly gen double `o_L' = L.x
        mata: st_matrix("t3", pyx_triple("`o_L'"))
        matrix row[1, 36] = t3[1, 1]
        matrix row[1, 37] = t3[1, 2]
        matrix row[1, 38] = t3[1, 3]
        tempvar o_L2
        quietly gen double `o_L2' = L2.x
        mata: st_matrix("t3", pyx_triple("`o_L2'"))
        matrix row[1, 39] = t3[1, 1]
        matrix row[1, 40] = t3[1, 2]
        matrix row[1, 41] = t3[1, 3]
        tempvar o_F
        quietly gen double `o_F' = F.x
        mata: st_matrix("t3", pyx_triple("`o_F'"))
        matrix row[1, 42] = t3[1, 1]
        matrix row[1, 43] = t3[1, 2]
        matrix row[1, 44] = t3[1, 3]
        tempvar o_D
        quietly gen double `o_D' = D.x
        mata: st_matrix("t3", pyx_triple("`o_D'"))
        matrix row[1, 45] = t3[1, 1]
        matrix row[1, 46] = t3[1, 2]
        matrix row[1, 47] = t3[1, 3]
        tempvar o_D2
        quietly gen double `o_D2' = D2.x
        mata: st_matrix("t3", pyx_triple("`o_D2'"))
        matrix row[1, 48] = t3[1, 1]
        matrix row[1, 49] = t3[1, 2]
        matrix row[1, 50] = t3[1, 3]
        tempvar o_S2
        quietly gen double `o_S2' = S2.x
        mata: st_matrix("t3", pyx_triple("`o_S2'"))
        matrix row[1, 51] = t3[1, 1]
        matrix row[1, 52] = t3[1, 2]
        matrix row[1, 53] = t3[1, 3]
        tempvar mx dev bx
        quietly bysort id: egen double `mx' = mean(x)
        quietly gen double `dev' = x - `mx'
        quietly gen double `bx' = `mx' if !missing(x)
        mata: st_matrix("w2", pyx_pair("`dev'"))
        mata: st_matrix("b2", pyx_pair("`bx'"))
        matrix row[1, 54] = w2[1, 1]
        matrix row[1, 55] = w2[1, 2]
        matrix row[1, 56] = b2[1, 3]
        matrix row[1, 57] = b2[1, 2]
        capture quietly xtdescribe
        matrix row[1, 58] = _rc
        if _rc == 0 {
        matrix row[1, 59] = r(N)
        matrix row[1, 60] = r(min)
        matrix row[1, 61] = r(max)
        }
        tempfile base
        quietly save `base'
        capture tsfill
        if _rc == 0 {
        matrix row[1, 62] = _N
            quietly count if !missing(x)
        matrix row[1, 63] = r(N)
        }
        quietly use `base', clear
        quietly xtset id t, delta(`delta')
        capture tsfill, full
        if _rc == 0 {
        matrix row[1, 64] = _N
        }
        quietly use `base', clear
        quietly xtset id t, delta(`delta')
        quietly generate byte dz = 0 if z < -0.5
        quietly replace dz = 1 if z >= -0.5 & z < 0.5
        quietly replace dz = 2 if z >= 0.5 & !missing(z)
        capture quietly xttrans dz
        if _rc == 0 {
        matrix row[1, 65] = r(n_trans)
        matrix row[1, 66] = r(n_rows)
        matrix row[1, 67] = r(n_cols)
        }
        quietly use `base', clear
        quietly xtset id t, delta(`delta')
        quietly keep if !missing(x, z)
        local rcdata = 0
        // a variable constant up to rounding has no correlation: Stata's xtdata leaves noise
        // of 1e-16 in it and correlate reports a number made of that noise, so the battery
        // does not read it (the package returns a missing value, see describe/corr.py)
        if `rcdata' == 0 {
            quietly summarize x
            local degenerate = missing(r(sd)) | r(sd) <= 1e-12 * max(abs(r(min)), abs(r(max)))
            quietly summarize z
            local degenerate = `degenerate' | missing(r(sd)) | r(sd) <= 1e-12 * max(abs(r(min)), abs(r(max)))
            capture quietly correlate x z
            if _rc == 0 & !`degenerate' {
        matrix row[1, 68] = r(rho)
        matrix row[1, 69] = r(N)
            }
        }
        quietly use `base', clear
        quietly xtset id t, delta(`delta')
        quietly keep if !missing(x, z)
        capture quietly xtdata x z, fe clear
        local rcdata = _rc
        // a variable constant up to rounding has no correlation: Stata's xtdata leaves noise
        // of 1e-16 in it and correlate reports a number made of that noise, so the battery
        // does not read it (the package returns a missing value, see describe/corr.py)
        if `rcdata' == 0 {
            quietly summarize x
            local degenerate = missing(r(sd)) | r(sd) <= 1e-12 * max(abs(r(min)), abs(r(max)))
            quietly summarize z
            local degenerate = `degenerate' | missing(r(sd)) | r(sd) <= 1e-12 * max(abs(r(min)), abs(r(max)))
            capture quietly correlate x z
            if _rc == 0 & !`degenerate' {
        matrix row[1, 70] = r(rho)
        matrix row[1, 71] = r(N)
            }
        }
        quietly use `base', clear
        quietly xtset id t, delta(`delta')
        quietly keep if !missing(x, z)
        capture quietly xtdata x z, be clear
        local rcdata = _rc
        // a variable constant up to rounding has no correlation: Stata's xtdata leaves noise
        // of 1e-16 in it and correlate reports a number made of that noise, so the battery
        // does not read it (the package returns a missing value, see describe/corr.py)
        if `rcdata' == 0 {
            quietly summarize x
            local degenerate = missing(r(sd)) | r(sd) <= 1e-12 * max(abs(r(min)), abs(r(max)))
            quietly summarize z
            local degenerate = `degenerate' | missing(r(sd)) | r(sd) <= 1e-12 * max(abs(r(min)), abs(r(max)))
            capture quietly correlate x z
            if _rc == 0 & !`degenerate' {
        matrix row[1, 72] = r(rho)
        matrix row[1, 73] = r(N)
            }
        }
    }
    mata: R = R \ st_matrix("row")
    restore
end

* ---- case campaign
use "../data/sim_campaign.dta", clear
mata: R = J(0, 73, .)
pyx_battery 0 1 ""
pyx_battery 1 1 ""
pyx_battery 2 1 ""
pyx_battery 3 1 ""
pyx_battery 4 1 "%td"
pyx_battery 5 7 "%td"
pyx_battery 6 1 ""
pyx_battery 7 3 "%tm"
pyx_battery 8 3 ""
pyx_battery 9 1 ""
pyx_battery 10 1 ""
pyx_battery 11 1 "%th"
pyx_battery 12 7 "%td"
pyx_battery 13 1 "%td"
pyx_battery 14 1 ""
pyx_battery 15 1 "%tw"
pyx_battery 16 1 "%tm"
pyx_battery 17 1 ""
pyx_battery 18 7 ""
pyx_battery 19 3 "%tw"
pyx_battery 20 1 ""
pyx_battery 21 1 ""
pyx_battery 22 1 ""
pyx_battery 23 3 "%tq"
pyx_battery 24 7 "%tq"
pyx_battery 25 1 ""
pyx_battery 26 5 "%ty"
pyx_battery 27 1 ""
pyx_battery 28 5 "%tq"
pyx_battery 29 3 "%ty"
pyx_battery 30 5 "%tm"
pyx_battery 31 5 "%th"
pyx_battery 32 1 "%tq"
pyx_battery 33 7 ""
pyx_battery 34 1 ""
pyx_battery 35 2 "%tq"
pyx_battery 36 2 "%tq"
pyx_battery 37 1 "%tw"
pyx_battery 38 3 ""
pyx_battery 39 1 "%ty"
pyx_battery 40 5 "%td"
pyx_battery 41 3 ""
pyx_battery 42 1 "%tw"
pyx_battery 43 1 ""
pyx_battery 44 1 "%tm"
pyx_battery 45 1 ""
pyx_battery 46 3 "%th"
pyx_battery 47 1 ""
pyx_battery 48 1 ""
pyx_battery 49 2 "%ty"
pyx_battery 50 1 "%td"
pyx_battery 51 5 ""
pyx_battery 52 1 "%tm"
pyx_battery 53 3 "%ty"
pyx_battery 54 1 "%tw"
pyx_battery 55 3 "%td"
pyx_battery 56 1 ""
pyx_battery 57 1 ""
pyx_battery 58 1 "%th"
pyx_battery 59 2 ""
pyx_battery 60 1 ""
pyx_battery 61 2 ""
pyx_battery 62 1 "%tm"
pyx_battery 63 3 "%tm"
pyx_battery 64 5 "%td"
pyx_battery 65 7 ""
pyx_battery 66 1 "%td"
pyx_battery 67 3 ""
pyx_battery 68 3 "%tm"
pyx_battery 69 1 ""
pyx_battery 70 1 ""
pyx_battery 71 7 ""
pyx_battery 72 7 "%th"
pyx_battery 73 1 "%ty"
pyx_battery 74 3 ""
pyx_battery 75 7 ""
pyx_battery 76 1 ""
pyx_battery 77 1 "%tw"
pyx_battery 78 1 "%tq"
pyx_battery 79 1 ""
pyx_battery 80 5 ""
pyx_battery 81 5 ""
pyx_battery 82 2 ""
pyx_battery 83 1 ""
pyx_battery 84 1 ""
pyx_battery 85 1 ""
pyx_battery 86 2 ""
pyx_battery 87 5 ""
pyx_battery 88 7 ""
pyx_battery 89 1 "%tq"
pyx_battery 90 1 ""
pyx_battery 91 1 ""
pyx_battery 92 2 ""
pyx_battery 93 2 "%th"
pyx_battery 94 1 ""
pyx_battery 95 1 ""
pyx_battery 96 1 "%ty"
pyx_battery 97 1 ""
pyx_battery 98 1 ""
pyx_battery 99 2 ""
pyx_battery 100 1 ""
pyx_battery 101 3 "%tm"
pyx_battery 102 1 "%th"
pyx_battery 103 3 ""
pyx_battery 104 1 "%ty"
pyx_battery 105 1 ""
pyx_battery 106 7 "%tw"
pyx_battery 107 1 ""
pyx_battery 108 5 ""
pyx_battery 109 1 ""
pyx_battery 110 1 "%tm"
pyx_battery 111 1 "%tq"
pyx_battery 112 2 "%tq"
pyx_battery 113 1 "%tw"
pyx_battery 114 1 "%tq"
pyx_battery 115 5 ""
pyx_battery 116 5 "%td"
pyx_battery 117 1 "%td"
pyx_battery 118 1 ""
pyx_battery 119 7 "%tm"
pyx_battery 120 1 "%td"
pyx_battery 121 5 "%tq"
pyx_battery 122 7 "%tm"
pyx_battery 123 1 ""
pyx_battery 124 1 ""
pyx_battery 125 1 ""
pyx_battery 126 1 ""
pyx_battery 127 7 ""
pyx_battery 128 1 ""
pyx_battery 129 1 ""
pyx_battery 130 1 "%ty"
pyx_battery 131 1 ""
pyx_battery 132 3 "%tm"
pyx_battery 133 7 "%td"
pyx_battery 134 1 ""
pyx_battery 135 3 "%tw"
pyx_battery 136 1 "%th"
pyx_battery 137 2 "%tw"
pyx_battery 138 1 "%tm"
pyx_battery 139 2 "%tq"
pyx_battery 140 1 "%ty"
pyx_battery 141 1 "%th"
pyx_battery 142 3 "%tq"
pyx_battery 143 1 ""
pyx_battery 144 1 "%tm"
pyx_battery 145 1 "%td"
pyx_battery 146 5 ""
pyx_battery 147 2 "%tw"
pyx_battery 148 1 "%tq"
pyx_battery 149 7 ""
pyx_battery 150 1 ""
pyx_battery 151 1 "%td"
pyx_battery 152 5 ""
pyx_battery 153 3 "%tq"
pyx_battery 154 5 "%tw"
pyx_battery 155 3 ""
pyx_battery 156 1 ""
pyx_battery 157 2 ""
pyx_battery 158 2 ""
pyx_battery 159 7 "%td"
mata: pyx_export("../out/campaign__battery.csv", R)

exit, clear
