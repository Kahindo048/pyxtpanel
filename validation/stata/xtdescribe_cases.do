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
real matrix pyx_patterns(string scalar fname)
{
    real scalar fh, p, k, j, nT
    string scalar line, left, right, c
    string rowvector tok
    real matrix M
    real rowvector row
    fh = fopen(fname, "r")
    nT = 0
    while ((line = fget(fh)) != J(0,0,"")) {
        p = strpos(line, "|")
        if (p == 0) continue
        right = strtrim(substr(line, p + 1, .))
        if (substr(right, 1, 1) != "(" & right != "Pattern" & substr(right, 1, 7) != "Pattern") nT = max((nT, strlen(right)))
    }
    fclose(fh)
    fh = fopen(fname, "r")
    M = J(0, 3 + nT, .)
    while ((line = fget(fh)) != J(0,0,"")) {
        p = strpos(line, "|")
        if (p == 0) continue
        left = strtrim(substr(line, 1, p - 1))
        right = strtrim(substr(line, p + 1, .))
        tok = tokens(left)
        if (cols(tok) < 2 | cols(tok) > 3) continue
        if (missing(strtoreal(tok[1]))) continue
        row = J(1, 3 + nT, -1)
        for (k = 1; k <= cols(tok); k++) row[k] = strtoreal(tok[k])
        if (substr(right, 1, 1) != "(") {
            for (j = 1; j <= strlen(right); j++) {
                c = substr(right, j, 1)
                row[3 + j] = (c == "." ? 0 : (c == "X" ? 10 : strtoreal(c)))
            }
        }
        M = M \ row
    }
    fclose(fh)
    return(M)
}
end


* ---- case xtdescribe_balanced_p9_w100
use "../data/sim_balanced.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_balanced_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_balanced_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_balanced_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_balanced_p9_w100.log"))

* ---- case xtdescribe_balanced_p3_w100
use "../data/sim_balanced.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_balanced_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_balanced_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_balanced_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_balanced_p3_w100.log"))

* ---- case xtdescribe_balanced_p1000_w100
use "../data/sim_balanced.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_balanced_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_balanced_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_balanced_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_balanced_p1000_w100.log"))

* ---- case xtdescribe_balanced_p0_w100
use "../data/sim_balanced.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_balanced_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_balanced_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_balanced_p9_w40
use "../data/sim_balanced.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_balanced_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_balanced_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_balanced_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_balanced_p9_w40.log"))

* ---- case xtdescribe_gaps_p9_w100
use "../data/sim_gaps.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_gaps_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_gaps_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_gaps_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_gaps_p9_w100.log"))

* ---- case xtdescribe_gaps_p3_w100
use "../data/sim_gaps.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_gaps_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_gaps_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_gaps_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_gaps_p3_w100.log"))

* ---- case xtdescribe_gaps_p1000_w100
use "../data/sim_gaps.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_gaps_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_gaps_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_gaps_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_gaps_p1000_w100.log"))

* ---- case xtdescribe_gaps_p0_w100
use "../data/sim_gaps.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_gaps_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_gaps_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_gaps_p9_w40
use "../data/sim_gaps.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_gaps_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_gaps_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_gaps_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_gaps_p9_w40.log"))

* ---- case xtdescribe_entry_exit_p9_w100
use "../data/sim_entry_exit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_entry_exit_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_entry_exit_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_entry_exit_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_entry_exit_p9_w100.log"))

* ---- case xtdescribe_entry_exit_p3_w100
use "../data/sim_entry_exit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_entry_exit_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_entry_exit_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_entry_exit_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_entry_exit_p3_w100.log"))

* ---- case xtdescribe_entry_exit_p1000_w100
use "../data/sim_entry_exit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_entry_exit_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_entry_exit_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_entry_exit_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_entry_exit_p1000_w100.log"))

* ---- case xtdescribe_entry_exit_p0_w100
use "../data/sim_entry_exit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_entry_exit_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_entry_exit_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_entry_exit_p9_w40
use "../data/sim_entry_exit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_entry_exit_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_entry_exit_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_entry_exit_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_entry_exit_p9_w40.log"))

* ---- case xtdescribe_weak_p9_w100
use "../data/sim_weak.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_weak_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_weak_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_weak_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_weak_p9_w100.log"))

* ---- case xtdescribe_weak_p3_w100
use "../data/sim_weak.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_weak_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_weak_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_weak_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_weak_p3_w100.log"))

* ---- case xtdescribe_weak_p1000_w100
use "../data/sim_weak.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_weak_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_weak_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_weak_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_weak_p1000_w100.log"))

* ---- case xtdescribe_weak_p0_w100
use "../data/sim_weak.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_weak_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_weak_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_weak_p9_w40
use "../data/sim_weak.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_weak_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_weak_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_weak_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_weak_p9_w40.log"))

* ---- case xtdescribe_delta2_p9_w100
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
set linesize 255
quietly log using "../out/xtdescribe_delta2_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_delta2_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_delta2_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_delta2_p9_w100.log"))

* ---- case xtdescribe_delta2_p3_w100
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
set linesize 255
quietly log using "../out/xtdescribe_delta2_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_delta2_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_delta2_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_delta2_p3_w100.log"))

* ---- case xtdescribe_delta2_p1000_w100
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
set linesize 255
quietly log using "../out/xtdescribe_delta2_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_delta2_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_delta2_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_delta2_p1000_w100.log"))

* ---- case xtdescribe_delta2_p0_w100
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
set linesize 255
quietly log using "../out/xtdescribe_delta2_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_delta2_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_delta2_p9_w40
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
set linesize 255
quietly log using "../out/xtdescribe_delta2_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_delta2_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_delta2_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_delta2_p9_w40.log"))

* ---- case xtdescribe_quarterly_p9_w100
use "../data/sim_quarterly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_quarterly_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_quarterly_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_quarterly_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_quarterly_p9_w100.log"))

* ---- case xtdescribe_quarterly_p3_w100
use "../data/sim_quarterly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_quarterly_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_quarterly_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_quarterly_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_quarterly_p3_w100.log"))

* ---- case xtdescribe_quarterly_p1000_w100
use "../data/sim_quarterly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_quarterly_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_quarterly_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_quarterly_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_quarterly_p1000_w100.log"))

* ---- case xtdescribe_quarterly_p0_w100
use "../data/sim_quarterly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_quarterly_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_quarterly_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_quarterly_p9_w40
use "../data/sim_quarterly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_quarterly_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_quarterly_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_quarterly_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_quarterly_p9_w40.log"))

* ---- case xtdescribe_monthly_p9_w100
use "../data/sim_monthly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_monthly_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_monthly_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_monthly_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_monthly_p9_w100.log"))

* ---- case xtdescribe_monthly_p3_w100
use "../data/sim_monthly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_monthly_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_monthly_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_monthly_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_monthly_p3_w100.log"))

* ---- case xtdescribe_monthly_p1000_w100
use "../data/sim_monthly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_monthly_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_monthly_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_monthly_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_monthly_p1000_w100.log"))

* ---- case xtdescribe_monthly_p0_w100
use "../data/sim_monthly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_monthly_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_monthly_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_monthly_p9_w40
use "../data/sim_monthly.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_monthly_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_monthly_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_monthly_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_monthly_p9_w40.log"))

* ---- case xtdescribe_unsorted_ids_p9_w100
use "../data/sim_unsorted_ids.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_unsorted_ids_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_unsorted_ids_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_unsorted_ids_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_unsorted_ids_p9_w100.log"))

* ---- case xtdescribe_unsorted_ids_p3_w100
use "../data/sim_unsorted_ids.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_unsorted_ids_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_unsorted_ids_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_unsorted_ids_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_unsorted_ids_p3_w100.log"))

* ---- case xtdescribe_unsorted_ids_p1000_w100
use "../data/sim_unsorted_ids.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_unsorted_ids_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_unsorted_ids_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_unsorted_ids_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_unsorted_ids_p1000_w100.log"))

* ---- case xtdescribe_unsorted_ids_p0_w100
use "../data/sim_unsorted_ids.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_unsorted_ids_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_unsorted_ids_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_unsorted_ids_p9_w40
use "../data/sim_unsorted_ids.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_unsorted_ids_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_unsorted_ids_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_unsorted_ids_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_unsorted_ids_p9_w40.log"))

* ---- case xtdescribe_one_obs_unit_p9_w100
use "../data/sim_one_obs_unit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_one_obs_unit_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_one_obs_unit_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_one_obs_unit_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_one_obs_unit_p9_w100.log"))

* ---- case xtdescribe_one_obs_unit_p3_w100
use "../data/sim_one_obs_unit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_one_obs_unit_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_one_obs_unit_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_one_obs_unit_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_one_obs_unit_p3_w100.log"))

* ---- case xtdescribe_one_obs_unit_p1000_w100
use "../data/sim_one_obs_unit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_one_obs_unit_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_one_obs_unit_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_one_obs_unit_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_one_obs_unit_p1000_w100.log"))

* ---- case xtdescribe_one_obs_unit_p0_w100
use "../data/sim_one_obs_unit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_one_obs_unit_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_one_obs_unit_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_one_obs_unit_p9_w40
use "../data/sim_one_obs_unit.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_one_obs_unit_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_one_obs_unit_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_one_obs_unit_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_one_obs_unit_p9_w40.log"))

* ---- case xtdescribe_re_effects_p9_w100
use "../data/sim_re_effects.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_re_effects_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_re_effects_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_re_effects_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_re_effects_p9_w100.log"))

* ---- case xtdescribe_re_effects_p3_w100
use "../data/sim_re_effects.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_re_effects_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_re_effects_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_re_effects_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_re_effects_p3_w100.log"))

* ---- case xtdescribe_re_effects_p1000_w100
use "../data/sim_re_effects.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_re_effects_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_re_effects_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_re_effects_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_re_effects_p1000_w100.log"))

* ---- case xtdescribe_re_effects_p0_w100
use "../data/sim_re_effects.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_re_effects_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_re_effects_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_re_effects_p9_w40
use "../data/sim_re_effects.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_re_effects_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_re_effects_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_re_effects_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_re_effects_p9_w40.log"))

* ---- case xtdescribe_sparse_t_p9_w100
use "../data/sim_sparse_t.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_sparse_t_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_sparse_t_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_sparse_t_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_sparse_t_p9_w100.log"))

* ---- case xtdescribe_sparse_t_p3_w100
use "../data/sim_sparse_t.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_sparse_t_p3_w100.log", text replace name(pyxlog)
xtdescribe, patterns(3) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_sparse_t_p3_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_sparse_t_p3_w100__patterns.csv", pyx_patterns("../out/xtdescribe_sparse_t_p3_w100.log"))

* ---- case xtdescribe_sparse_t_p1000_w100
use "../data/sim_sparse_t.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_sparse_t_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_sparse_t_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_sparse_t_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_sparse_t_p1000_w100.log"))

* ---- case xtdescribe_sparse_t_p0_w100
use "../data/sim_sparse_t.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_sparse_t_p0_w100.log", text replace name(pyxlog)
xtdescribe, patterns(0) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_sparse_t_p0_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))

* ---- case xtdescribe_sparse_t_p9_w40
use "../data/sim_sparse_t.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_sparse_t_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_sparse_t_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_sparse_t_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_sparse_t_p9_w40.log"))

* ---- case xtdescribe_wide_p9_w100
use "../data/sim_wide.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_wide_p9_w100.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_wide_p9_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_wide_p9_w100__patterns.csv", pyx_patterns("../out/xtdescribe_wide_p9_w100.log"))

* ---- case xtdescribe_wide_p1000_w100
use "../data/sim_wide.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_wide_p1000_w100.log", text replace name(pyxlog)
xtdescribe, patterns(1000) width(100)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_wide_p1000_w100__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_wide_p1000_w100__patterns.csv", pyx_patterns("../out/xtdescribe_wide_p1000_w100.log"))

* ---- case xtdescribe_wide_p9_w40
use "../data/sim_wide.dta", clear
xtset id t
set linesize 255
quietly log using "../out/xtdescribe_wide_p9_w40.log", text replace name(pyxlog)
xtdescribe, patterns(9) width(40)
scalar s_N = r(N)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_p1 = r(p1)
scalar s_p5 = r(p5)
scalar s_p10 = r(p10)
scalar s_p25 = r(p25)
scalar s_p50 = r(p50)
scalar s_p75 = r(p75)
scalar s_p90 = r(p90)
scalar s_p95 = r(p95)
scalar s_p99 = r(p99)
scalar s_max = r(max)
quietly log close pyxlog
mata: pyx_export("../out/xtdescribe_wide_p9_w40__ti.csv", (st_numscalar("s_N") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_p1") \ st_numscalar("s_p5") \ st_numscalar("s_p10") \ st_numscalar("s_p25") \ st_numscalar("s_p50") \ st_numscalar("s_p75") \ st_numscalar("s_p90") \ st_numscalar("s_p95") \ st_numscalar("s_p99") \ st_numscalar("s_max")))
mata: pyx_export("../out/xtdescribe_wide_p9_w40__patterns.csv", pyx_patterns("../out/xtdescribe_wide_p9_w40.log"))

exit, clear
