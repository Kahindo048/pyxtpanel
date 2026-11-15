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

* ---- case variation_balanced
use "../data/sim_balanced.dta", clear
preserve
quietly keep if !missing(x)
bysort id: egen double mn = min(x)
bysort id: egen double mx = max(x)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_x = r(N)
quietly count if u1 & cnt == 1
scalar us_x = r(N)
bysort t: egen double tmn = min(x)
bysort t: egen double tmx = max(x)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_x = r(N)
bysort x: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_x = r(N)
restore
preserve
quietly keep if !missing(z)
bysort id: egen double mn = min(z)
bysort id: egen double mx = max(z)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_z = r(N)
quietly count if u1 & cnt == 1
scalar us_z = r(N)
bysort t: egen double tmn = min(z)
bysort t: egen double tmx = max(z)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_z = r(N)
bysort z: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_z = r(N)
restore
mata: pyx_export("../out/variation_balanced__counts.csv", (st_numscalar("uv_x") \ st_numscalar("us_x") \ st_numscalar("pv_x") \ st_numscalar("nd_x") \ st_numscalar("uv_z") \ st_numscalar("us_z") \ st_numscalar("pv_z") \ st_numscalar("nd_z")))

* ---- case spells_balanced
use "../data/sim_balanced.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate byte ld = L.d
generate byte start = !missing(d) & (missing(ld) | ld != d)
generate byte change = !missing(d) & !missing(ld) & ld != d
quietly count if start
scalar n_spells = r(N)
quietly count if change
scalar n_changes = r(N)
bysort id: egen byte any_change = max(change)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & any_change == 1
scalar n_units_changing = r(N)
mata: pyx_export("../out/spells_balanced__counts.csv", (st_numscalar("n_spells") \ st_numscalar("n_changes") \ st_numscalar("n_units_changing")))

* ---- case event_balanced
use "../data/sim_balanced.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate double on = cond(d == 1, t, .)
bysort id: egen double ev = min(on)
generate double rel = (t - ev) / 1
quietly summarize rel
scalar n_rel = r(N)
scalar sum_rel = r(sum)
generate double rel2 = rel^2
quietly summarize rel2
scalar sumsq_rel = r(sum)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & !missing(ev)
scalar n_treated = r(N)
mata: pyx_export("../out/event_balanced__counts.csv", (st_numscalar("n_rel") \ st_numscalar("sum_rel") \ st_numscalar("sumsq_rel") \ st_numscalar("n_treated")))

* ---- case outliers_balanced
use "../data/sim_balanced.dta", clear
bysort id: egen long n_x = count(x)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_x = pctile(x), p(50)
generate double dev_x = abs(x - med_x)
bysort id: egen double mad_x = pctile(dev_x), p(50)
generate double z_x = 0.6745 * (x - med_x) / mad_x if n_x >= 3 & mad_x > 0 & !missing(x)
generate double az_x = abs(z_x)
quietly summarize az_x
scalar ns_x = r(N)
scalar sa_x = r(sum)
quietly count if az_x > 2 & !missing(az_x)
scalar nf_x = r(N)
bysort id: egen long n_z = count(z)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_z = pctile(z), p(50)
generate double dev_z = abs(z - med_z)
bysort id: egen double mad_z = pctile(dev_z), p(50)
generate double z_z = 0.6745 * (z - med_z) / mad_z if n_z >= 3 & mad_z > 0 & !missing(z)
generate double az_z = abs(z_z)
quietly summarize az_z
scalar ns_z = r(N)
scalar sa_z = r(sum)
quietly count if az_z > 2 & !missing(az_z)
scalar nf_z = r(N)
mata: pyx_export("../out/outliers_balanced__counts.csv", (st_numscalar("ns_x") \ st_numscalar("nf_x") \ st_numscalar("sa_x") \ st_numscalar("ns_z") \ st_numscalar("nf_z") \ st_numscalar("sa_z")))

* ---- case wide_balanced
use "../data/sim_balanced.dta", clear
keep id t x z
reshape wide x z, i(id) j(t)
scalar ncols = c(k)
scalar nunits = _N
egen double tot = rowtotal(x* z*)
quietly summarize tot
scalar total = r(sum)
egen long nm = rownonmiss(x* z*)
quietly summarize nm
scalar nonmiss = r(sum)
mata: M = editmissing(st_data(., "x* z*"), 0)
mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))
mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))
mata: pyx_export("../out/wide_balanced__counts.csv", (st_numscalar("ncols") \ st_numscalar("nunits") \ st_numscalar("total") \ st_numscalar("nonmiss") \ st_numscalar("colw") \ st_numscalar("roww")))

* ---- case gaps_balanced
use "../data/sim_balanced.dta", clear
xtset id t
quietly tsreport, panel
scalar n_gaps = r(N_gaps)
scalar n_obs = r(N_obs)
mata: pyx_export("../out/gaps_balanced__counts.csv", (st_numscalar("n_gaps") \ st_numscalar("n_obs")))

* ---- case variation_gaps
use "../data/sim_gaps.dta", clear
preserve
quietly keep if !missing(x)
bysort id: egen double mn = min(x)
bysort id: egen double mx = max(x)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_x = r(N)
quietly count if u1 & cnt == 1
scalar us_x = r(N)
bysort t: egen double tmn = min(x)
bysort t: egen double tmx = max(x)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_x = r(N)
bysort x: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_x = r(N)
restore
preserve
quietly keep if !missing(z)
bysort id: egen double mn = min(z)
bysort id: egen double mx = max(z)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_z = r(N)
quietly count if u1 & cnt == 1
scalar us_z = r(N)
bysort t: egen double tmn = min(z)
bysort t: egen double tmx = max(z)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_z = r(N)
bysort z: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_z = r(N)
restore
mata: pyx_export("../out/variation_gaps__counts.csv", (st_numscalar("uv_x") \ st_numscalar("us_x") \ st_numscalar("pv_x") \ st_numscalar("nd_x") \ st_numscalar("uv_z") \ st_numscalar("us_z") \ st_numscalar("pv_z") \ st_numscalar("nd_z")))

* ---- case spells_gaps
use "../data/sim_gaps.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate byte ld = L.d
generate byte start = !missing(d) & (missing(ld) | ld != d)
generate byte change = !missing(d) & !missing(ld) & ld != d
quietly count if start
scalar n_spells = r(N)
quietly count if change
scalar n_changes = r(N)
bysort id: egen byte any_change = max(change)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & any_change == 1
scalar n_units_changing = r(N)
mata: pyx_export("../out/spells_gaps__counts.csv", (st_numscalar("n_spells") \ st_numscalar("n_changes") \ st_numscalar("n_units_changing")))

* ---- case event_gaps
use "../data/sim_gaps.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate double on = cond(d == 1, t, .)
bysort id: egen double ev = min(on)
generate double rel = (t - ev) / 1
quietly summarize rel
scalar n_rel = r(N)
scalar sum_rel = r(sum)
generate double rel2 = rel^2
quietly summarize rel2
scalar sumsq_rel = r(sum)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & !missing(ev)
scalar n_treated = r(N)
mata: pyx_export("../out/event_gaps__counts.csv", (st_numscalar("n_rel") \ st_numscalar("sum_rel") \ st_numscalar("sumsq_rel") \ st_numscalar("n_treated")))

* ---- case outliers_gaps
use "../data/sim_gaps.dta", clear
bysort id: egen long n_x = count(x)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_x = pctile(x), p(50)
generate double dev_x = abs(x - med_x)
bysort id: egen double mad_x = pctile(dev_x), p(50)
generate double z_x = 0.6745 * (x - med_x) / mad_x if n_x >= 3 & mad_x > 0 & !missing(x)
generate double az_x = abs(z_x)
quietly summarize az_x
scalar ns_x = r(N)
scalar sa_x = r(sum)
quietly count if az_x > 2 & !missing(az_x)
scalar nf_x = r(N)
bysort id: egen long n_z = count(z)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_z = pctile(z), p(50)
generate double dev_z = abs(z - med_z)
bysort id: egen double mad_z = pctile(dev_z), p(50)
generate double z_z = 0.6745 * (z - med_z) / mad_z if n_z >= 3 & mad_z > 0 & !missing(z)
generate double az_z = abs(z_z)
quietly summarize az_z
scalar ns_z = r(N)
scalar sa_z = r(sum)
quietly count if az_z > 2 & !missing(az_z)
scalar nf_z = r(N)
mata: pyx_export("../out/outliers_gaps__counts.csv", (st_numscalar("ns_x") \ st_numscalar("nf_x") \ st_numscalar("sa_x") \ st_numscalar("ns_z") \ st_numscalar("nf_z") \ st_numscalar("sa_z")))

* ---- case wide_gaps
use "../data/sim_gaps.dta", clear
keep id t x z
reshape wide x z, i(id) j(t)
scalar ncols = c(k)
scalar nunits = _N
egen double tot = rowtotal(x* z*)
quietly summarize tot
scalar total = r(sum)
egen long nm = rownonmiss(x* z*)
quietly summarize nm
scalar nonmiss = r(sum)
mata: M = editmissing(st_data(., "x* z*"), 0)
mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))
mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))
mata: pyx_export("../out/wide_gaps__counts.csv", (st_numscalar("ncols") \ st_numscalar("nunits") \ st_numscalar("total") \ st_numscalar("nonmiss") \ st_numscalar("colw") \ st_numscalar("roww")))

* ---- case gaps_gaps
use "../data/sim_gaps.dta", clear
xtset id t
quietly tsreport, panel
scalar n_gaps = r(N_gaps)
scalar n_obs = r(N_obs)
mata: pyx_export("../out/gaps_gaps__counts.csv", (st_numscalar("n_gaps") \ st_numscalar("n_obs")))

* ---- case variation_entry_exit
use "../data/sim_entry_exit.dta", clear
preserve
quietly keep if !missing(x)
bysort id: egen double mn = min(x)
bysort id: egen double mx = max(x)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_x = r(N)
quietly count if u1 & cnt == 1
scalar us_x = r(N)
bysort t: egen double tmn = min(x)
bysort t: egen double tmx = max(x)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_x = r(N)
bysort x: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_x = r(N)
restore
preserve
quietly keep if !missing(z)
bysort id: egen double mn = min(z)
bysort id: egen double mx = max(z)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_z = r(N)
quietly count if u1 & cnt == 1
scalar us_z = r(N)
bysort t: egen double tmn = min(z)
bysort t: egen double tmx = max(z)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_z = r(N)
bysort z: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_z = r(N)
restore
mata: pyx_export("../out/variation_entry_exit__counts.csv", (st_numscalar("uv_x") \ st_numscalar("us_x") \ st_numscalar("pv_x") \ st_numscalar("nd_x") \ st_numscalar("uv_z") \ st_numscalar("us_z") \ st_numscalar("pv_z") \ st_numscalar("nd_z")))

* ---- case spells_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate byte ld = L.d
generate byte start = !missing(d) & (missing(ld) | ld != d)
generate byte change = !missing(d) & !missing(ld) & ld != d
quietly count if start
scalar n_spells = r(N)
quietly count if change
scalar n_changes = r(N)
bysort id: egen byte any_change = max(change)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & any_change == 1
scalar n_units_changing = r(N)
mata: pyx_export("../out/spells_entry_exit__counts.csv", (st_numscalar("n_spells") \ st_numscalar("n_changes") \ st_numscalar("n_units_changing")))

* ---- case event_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate double on = cond(d == 1, t, .)
bysort id: egen double ev = min(on)
generate double rel = (t - ev) / 1
quietly summarize rel
scalar n_rel = r(N)
scalar sum_rel = r(sum)
generate double rel2 = rel^2
quietly summarize rel2
scalar sumsq_rel = r(sum)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & !missing(ev)
scalar n_treated = r(N)
mata: pyx_export("../out/event_entry_exit__counts.csv", (st_numscalar("n_rel") \ st_numscalar("sum_rel") \ st_numscalar("sumsq_rel") \ st_numscalar("n_treated")))

* ---- case outliers_entry_exit
use "../data/sim_entry_exit.dta", clear
bysort id: egen long n_x = count(x)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_x = pctile(x), p(50)
generate double dev_x = abs(x - med_x)
bysort id: egen double mad_x = pctile(dev_x), p(50)
generate double z_x = 0.6745 * (x - med_x) / mad_x if n_x >= 3 & mad_x > 0 & !missing(x)
generate double az_x = abs(z_x)
quietly summarize az_x
scalar ns_x = r(N)
scalar sa_x = r(sum)
quietly count if az_x > 2 & !missing(az_x)
scalar nf_x = r(N)
bysort id: egen long n_z = count(z)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_z = pctile(z), p(50)
generate double dev_z = abs(z - med_z)
bysort id: egen double mad_z = pctile(dev_z), p(50)
generate double z_z = 0.6745 * (z - med_z) / mad_z if n_z >= 3 & mad_z > 0 & !missing(z)
generate double az_z = abs(z_z)
quietly summarize az_z
scalar ns_z = r(N)
scalar sa_z = r(sum)
quietly count if az_z > 2 & !missing(az_z)
scalar nf_z = r(N)
mata: pyx_export("../out/outliers_entry_exit__counts.csv", (st_numscalar("ns_x") \ st_numscalar("nf_x") \ st_numscalar("sa_x") \ st_numscalar("ns_z") \ st_numscalar("nf_z") \ st_numscalar("sa_z")))

* ---- case wide_entry_exit
use "../data/sim_entry_exit.dta", clear
keep id t x z
reshape wide x z, i(id) j(t)
scalar ncols = c(k)
scalar nunits = _N
egen double tot = rowtotal(x* z*)
quietly summarize tot
scalar total = r(sum)
egen long nm = rownonmiss(x* z*)
quietly summarize nm
scalar nonmiss = r(sum)
mata: M = editmissing(st_data(., "x* z*"), 0)
mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))
mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))
mata: pyx_export("../out/wide_entry_exit__counts.csv", (st_numscalar("ncols") \ st_numscalar("nunits") \ st_numscalar("total") \ st_numscalar("nonmiss") \ st_numscalar("colw") \ st_numscalar("roww")))

* ---- case gaps_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
quietly tsreport, panel
scalar n_gaps = r(N_gaps)
scalar n_obs = r(N_obs)
mata: pyx_export("../out/gaps_entry_exit__counts.csv", (st_numscalar("n_gaps") \ st_numscalar("n_obs")))

* ---- case variation_weak
use "../data/sim_weak.dta", clear
preserve
quietly keep if !missing(x)
bysort id: egen double mn = min(x)
bysort id: egen double mx = max(x)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_x = r(N)
quietly count if u1 & cnt == 1
scalar us_x = r(N)
bysort t: egen double tmn = min(x)
bysort t: egen double tmx = max(x)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_x = r(N)
bysort x: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_x = r(N)
restore
preserve
quietly keep if !missing(z)
bysort id: egen double mn = min(z)
bysort id: egen double mx = max(z)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_z = r(N)
quietly count if u1 & cnt == 1
scalar us_z = r(N)
bysort t: egen double tmn = min(z)
bysort t: egen double tmx = max(z)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_z = r(N)
bysort z: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_z = r(N)
restore
mata: pyx_export("../out/variation_weak__counts.csv", (st_numscalar("uv_x") \ st_numscalar("us_x") \ st_numscalar("pv_x") \ st_numscalar("nd_x") \ st_numscalar("uv_z") \ st_numscalar("us_z") \ st_numscalar("pv_z") \ st_numscalar("nd_z")))

* ---- case spells_weak
use "../data/sim_weak.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate byte ld = L.d
generate byte start = !missing(d) & (missing(ld) | ld != d)
generate byte change = !missing(d) & !missing(ld) & ld != d
quietly count if start
scalar n_spells = r(N)
quietly count if change
scalar n_changes = r(N)
bysort id: egen byte any_change = max(change)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & any_change == 1
scalar n_units_changing = r(N)
mata: pyx_export("../out/spells_weak__counts.csv", (st_numscalar("n_spells") \ st_numscalar("n_changes") \ st_numscalar("n_units_changing")))

* ---- case event_weak
use "../data/sim_weak.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate double on = cond(d == 1, t, .)
bysort id: egen double ev = min(on)
generate double rel = (t - ev) / 1
quietly summarize rel
scalar n_rel = r(N)
scalar sum_rel = r(sum)
generate double rel2 = rel^2
quietly summarize rel2
scalar sumsq_rel = r(sum)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & !missing(ev)
scalar n_treated = r(N)
mata: pyx_export("../out/event_weak__counts.csv", (st_numscalar("n_rel") \ st_numscalar("sum_rel") \ st_numscalar("sumsq_rel") \ st_numscalar("n_treated")))

* ---- case outliers_weak
use "../data/sim_weak.dta", clear
bysort id: egen long n_x = count(x)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_x = pctile(x), p(50)
generate double dev_x = abs(x - med_x)
bysort id: egen double mad_x = pctile(dev_x), p(50)
generate double z_x = 0.6745 * (x - med_x) / mad_x if n_x >= 3 & mad_x > 0 & !missing(x)
generate double az_x = abs(z_x)
quietly summarize az_x
scalar ns_x = r(N)
scalar sa_x = r(sum)
quietly count if az_x > 2 & !missing(az_x)
scalar nf_x = r(N)
bysort id: egen long n_z = count(z)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_z = pctile(z), p(50)
generate double dev_z = abs(z - med_z)
bysort id: egen double mad_z = pctile(dev_z), p(50)
generate double z_z = 0.6745 * (z - med_z) / mad_z if n_z >= 3 & mad_z > 0 & !missing(z)
generate double az_z = abs(z_z)
quietly summarize az_z
scalar ns_z = r(N)
scalar sa_z = r(sum)
quietly count if az_z > 2 & !missing(az_z)
scalar nf_z = r(N)
mata: pyx_export("../out/outliers_weak__counts.csv", (st_numscalar("ns_x") \ st_numscalar("nf_x") \ st_numscalar("sa_x") \ st_numscalar("ns_z") \ st_numscalar("nf_z") \ st_numscalar("sa_z")))

* ---- case wide_weak
use "../data/sim_weak.dta", clear
keep id t x z
reshape wide x z, i(id) j(t)
scalar ncols = c(k)
scalar nunits = _N
egen double tot = rowtotal(x* z*)
quietly summarize tot
scalar total = r(sum)
egen long nm = rownonmiss(x* z*)
quietly summarize nm
scalar nonmiss = r(sum)
mata: M = editmissing(st_data(., "x* z*"), 0)
mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))
mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))
mata: pyx_export("../out/wide_weak__counts.csv", (st_numscalar("ncols") \ st_numscalar("nunits") \ st_numscalar("total") \ st_numscalar("nonmiss") \ st_numscalar("colw") \ st_numscalar("roww")))

* ---- case gaps_weak
use "../data/sim_weak.dta", clear
xtset id t
quietly tsreport, panel
scalar n_gaps = r(N_gaps)
scalar n_obs = r(N_obs)
mata: pyx_export("../out/gaps_weak__counts.csv", (st_numscalar("n_gaps") \ st_numscalar("n_obs")))

* ---- case variation_delta2
use "../data/sim_delta2.dta", clear
preserve
quietly keep if !missing(x)
bysort id: egen double mn = min(x)
bysort id: egen double mx = max(x)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_x = r(N)
quietly count if u1 & cnt == 1
scalar us_x = r(N)
bysort t: egen double tmn = min(x)
bysort t: egen double tmx = max(x)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_x = r(N)
bysort x: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_x = r(N)
restore
preserve
quietly keep if !missing(z)
bysort id: egen double mn = min(z)
bysort id: egen double mx = max(z)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_z = r(N)
quietly count if u1 & cnt == 1
scalar us_z = r(N)
bysort t: egen double tmn = min(z)
bysort t: egen double tmx = max(z)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_z = r(N)
bysort z: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_z = r(N)
restore
mata: pyx_export("../out/variation_delta2__counts.csv", (st_numscalar("uv_x") \ st_numscalar("us_x") \ st_numscalar("pv_x") \ st_numscalar("nd_x") \ st_numscalar("uv_z") \ st_numscalar("us_z") \ st_numscalar("pv_z") \ st_numscalar("nd_z")))

* ---- case spells_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
generate byte d = cond(missing(x), ., x > 0)
generate byte ld = L.d
generate byte start = !missing(d) & (missing(ld) | ld != d)
generate byte change = !missing(d) & !missing(ld) & ld != d
quietly count if start
scalar n_spells = r(N)
quietly count if change
scalar n_changes = r(N)
bysort id: egen byte any_change = max(change)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & any_change == 1
scalar n_units_changing = r(N)
mata: pyx_export("../out/spells_delta2__counts.csv", (st_numscalar("n_spells") \ st_numscalar("n_changes") \ st_numscalar("n_units_changing")))

* ---- case event_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
generate byte d = cond(missing(x), ., x > 0)
generate double on = cond(d == 1, t, .)
bysort id: egen double ev = min(on)
generate double rel = (t - ev) / 2
quietly summarize rel
scalar n_rel = r(N)
scalar sum_rel = r(sum)
generate double rel2 = rel^2
quietly summarize rel2
scalar sumsq_rel = r(sum)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & !missing(ev)
scalar n_treated = r(N)
mata: pyx_export("../out/event_delta2__counts.csv", (st_numscalar("n_rel") \ st_numscalar("sum_rel") \ st_numscalar("sumsq_rel") \ st_numscalar("n_treated")))

* ---- case outliers_delta2
use "../data/sim_delta2.dta", clear
bysort id: egen long n_x = count(x)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_x = pctile(x), p(50)
generate double dev_x = abs(x - med_x)
bysort id: egen double mad_x = pctile(dev_x), p(50)
generate double z_x = 0.6745 * (x - med_x) / mad_x if n_x >= 3 & mad_x > 0 & !missing(x)
generate double az_x = abs(z_x)
quietly summarize az_x
scalar ns_x = r(N)
scalar sa_x = r(sum)
quietly count if az_x > 2 & !missing(az_x)
scalar nf_x = r(N)
bysort id: egen long n_z = count(z)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_z = pctile(z), p(50)
generate double dev_z = abs(z - med_z)
bysort id: egen double mad_z = pctile(dev_z), p(50)
generate double z_z = 0.6745 * (z - med_z) / mad_z if n_z >= 3 & mad_z > 0 & !missing(z)
generate double az_z = abs(z_z)
quietly summarize az_z
scalar ns_z = r(N)
scalar sa_z = r(sum)
quietly count if az_z > 2 & !missing(az_z)
scalar nf_z = r(N)
mata: pyx_export("../out/outliers_delta2__counts.csv", (st_numscalar("ns_x") \ st_numscalar("nf_x") \ st_numscalar("sa_x") \ st_numscalar("ns_z") \ st_numscalar("nf_z") \ st_numscalar("sa_z")))

* ---- case wide_delta2
use "../data/sim_delta2.dta", clear
keep id t x z
reshape wide x z, i(id) j(t)
scalar ncols = c(k)
scalar nunits = _N
egen double tot = rowtotal(x* z*)
quietly summarize tot
scalar total = r(sum)
egen long nm = rownonmiss(x* z*)
quietly summarize nm
scalar nonmiss = r(sum)
mata: M = editmissing(st_data(., "x* z*"), 0)
mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))
mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))
mata: pyx_export("../out/wide_delta2__counts.csv", (st_numscalar("ncols") \ st_numscalar("nunits") \ st_numscalar("total") \ st_numscalar("nonmiss") \ st_numscalar("colw") \ st_numscalar("roww")))

* ---- case gaps_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
quietly tsreport, panel
scalar n_gaps = r(N_gaps)
scalar n_obs = r(N_obs)
mata: pyx_export("../out/gaps_delta2__counts.csv", (st_numscalar("n_gaps") \ st_numscalar("n_obs")))

* ---- case variation_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
preserve
quietly keep if !missing(x)
bysort id: egen double mn = min(x)
bysort id: egen double mx = max(x)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_x = r(N)
quietly count if u1 & cnt == 1
scalar us_x = r(N)
bysort t: egen double tmn = min(x)
bysort t: egen double tmx = max(x)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_x = r(N)
bysort x: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_x = r(N)
restore
preserve
quietly keep if !missing(z)
bysort id: egen double mn = min(z)
bysort id: egen double mx = max(z)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_z = r(N)
quietly count if u1 & cnt == 1
scalar us_z = r(N)
bysort t: egen double tmn = min(z)
bysort t: egen double tmx = max(z)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_z = r(N)
bysort z: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_z = r(N)
restore
mata: pyx_export("../out/variation_unsorted_ids__counts.csv", (st_numscalar("uv_x") \ st_numscalar("us_x") \ st_numscalar("pv_x") \ st_numscalar("nd_x") \ st_numscalar("uv_z") \ st_numscalar("us_z") \ st_numscalar("pv_z") \ st_numscalar("nd_z")))

* ---- case spells_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate byte ld = L.d
generate byte start = !missing(d) & (missing(ld) | ld != d)
generate byte change = !missing(d) & !missing(ld) & ld != d
quietly count if start
scalar n_spells = r(N)
quietly count if change
scalar n_changes = r(N)
bysort id: egen byte any_change = max(change)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & any_change == 1
scalar n_units_changing = r(N)
mata: pyx_export("../out/spells_unsorted_ids__counts.csv", (st_numscalar("n_spells") \ st_numscalar("n_changes") \ st_numscalar("n_units_changing")))

* ---- case event_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate double on = cond(d == 1, t, .)
bysort id: egen double ev = min(on)
generate double rel = (t - ev) / 1
quietly summarize rel
scalar n_rel = r(N)
scalar sum_rel = r(sum)
generate double rel2 = rel^2
quietly summarize rel2
scalar sumsq_rel = r(sum)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & !missing(ev)
scalar n_treated = r(N)
mata: pyx_export("../out/event_unsorted_ids__counts.csv", (st_numscalar("n_rel") \ st_numscalar("sum_rel") \ st_numscalar("sumsq_rel") \ st_numscalar("n_treated")))

* ---- case outliers_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
bysort id: egen long n_x = count(x)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_x = pctile(x), p(50)
generate double dev_x = abs(x - med_x)
bysort id: egen double mad_x = pctile(dev_x), p(50)
generate double z_x = 0.6745 * (x - med_x) / mad_x if n_x >= 3 & mad_x > 0 & !missing(x)
generate double az_x = abs(z_x)
quietly summarize az_x
scalar ns_x = r(N)
scalar sa_x = r(sum)
quietly count if az_x > 2 & !missing(az_x)
scalar nf_x = r(N)
bysort id: egen long n_z = count(z)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_z = pctile(z), p(50)
generate double dev_z = abs(z - med_z)
bysort id: egen double mad_z = pctile(dev_z), p(50)
generate double z_z = 0.6745 * (z - med_z) / mad_z if n_z >= 3 & mad_z > 0 & !missing(z)
generate double az_z = abs(z_z)
quietly summarize az_z
scalar ns_z = r(N)
scalar sa_z = r(sum)
quietly count if az_z > 2 & !missing(az_z)
scalar nf_z = r(N)
mata: pyx_export("../out/outliers_unsorted_ids__counts.csv", (st_numscalar("ns_x") \ st_numscalar("nf_x") \ st_numscalar("sa_x") \ st_numscalar("ns_z") \ st_numscalar("nf_z") \ st_numscalar("sa_z")))

* ---- case wide_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
keep id t x z
reshape wide x z, i(id) j(t)
scalar ncols = c(k)
scalar nunits = _N
egen double tot = rowtotal(x* z*)
quietly summarize tot
scalar total = r(sum)
egen long nm = rownonmiss(x* z*)
quietly summarize nm
scalar nonmiss = r(sum)
mata: M = editmissing(st_data(., "x* z*"), 0)
mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))
mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))
mata: pyx_export("../out/wide_unsorted_ids__counts.csv", (st_numscalar("ncols") \ st_numscalar("nunits") \ st_numscalar("total") \ st_numscalar("nonmiss") \ st_numscalar("colw") \ st_numscalar("roww")))

* ---- case gaps_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
quietly tsreport, panel
scalar n_gaps = r(N_gaps)
scalar n_obs = r(N_obs)
mata: pyx_export("../out/gaps_unsorted_ids__counts.csv", (st_numscalar("n_gaps") \ st_numscalar("n_obs")))

* ---- case variation_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
preserve
quietly keep if !missing(x)
bysort id: egen double mn = min(x)
bysort id: egen double mx = max(x)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_x = r(N)
quietly count if u1 & cnt == 1
scalar us_x = r(N)
bysort t: egen double tmn = min(x)
bysort t: egen double tmx = max(x)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_x = r(N)
bysort x: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_x = r(N)
restore
preserve
quietly keep if !missing(z)
bysort id: egen double mn = min(z)
bysort id: egen double mx = max(z)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_z = r(N)
quietly count if u1 & cnt == 1
scalar us_z = r(N)
bysort t: egen double tmn = min(z)
bysort t: egen double tmx = max(z)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_z = r(N)
bysort z: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_z = r(N)
restore
mata: pyx_export("../out/variation_one_obs_unit__counts.csv", (st_numscalar("uv_x") \ st_numscalar("us_x") \ st_numscalar("pv_x") \ st_numscalar("nd_x") \ st_numscalar("uv_z") \ st_numscalar("us_z") \ st_numscalar("pv_z") \ st_numscalar("nd_z")))

* ---- case spells_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate byte ld = L.d
generate byte start = !missing(d) & (missing(ld) | ld != d)
generate byte change = !missing(d) & !missing(ld) & ld != d
quietly count if start
scalar n_spells = r(N)
quietly count if change
scalar n_changes = r(N)
bysort id: egen byte any_change = max(change)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & any_change == 1
scalar n_units_changing = r(N)
mata: pyx_export("../out/spells_one_obs_unit__counts.csv", (st_numscalar("n_spells") \ st_numscalar("n_changes") \ st_numscalar("n_units_changing")))

* ---- case event_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate double on = cond(d == 1, t, .)
bysort id: egen double ev = min(on)
generate double rel = (t - ev) / 1
quietly summarize rel
scalar n_rel = r(N)
scalar sum_rel = r(sum)
generate double rel2 = rel^2
quietly summarize rel2
scalar sumsq_rel = r(sum)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & !missing(ev)
scalar n_treated = r(N)
mata: pyx_export("../out/event_one_obs_unit__counts.csv", (st_numscalar("n_rel") \ st_numscalar("sum_rel") \ st_numscalar("sumsq_rel") \ st_numscalar("n_treated")))

* ---- case outliers_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
bysort id: egen long n_x = count(x)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_x = pctile(x), p(50)
generate double dev_x = abs(x - med_x)
bysort id: egen double mad_x = pctile(dev_x), p(50)
generate double z_x = 0.6745 * (x - med_x) / mad_x if n_x >= 3 & mad_x > 0 & !missing(x)
generate double az_x = abs(z_x)
quietly summarize az_x
scalar ns_x = r(N)
scalar sa_x = r(sum)
quietly count if az_x > 2 & !missing(az_x)
scalar nf_x = r(N)
bysort id: egen long n_z = count(z)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_z = pctile(z), p(50)
generate double dev_z = abs(z - med_z)
bysort id: egen double mad_z = pctile(dev_z), p(50)
generate double z_z = 0.6745 * (z - med_z) / mad_z if n_z >= 3 & mad_z > 0 & !missing(z)
generate double az_z = abs(z_z)
quietly summarize az_z
scalar ns_z = r(N)
scalar sa_z = r(sum)
quietly count if az_z > 2 & !missing(az_z)
scalar nf_z = r(N)
mata: pyx_export("../out/outliers_one_obs_unit__counts.csv", (st_numscalar("ns_x") \ st_numscalar("nf_x") \ st_numscalar("sa_x") \ st_numscalar("ns_z") \ st_numscalar("nf_z") \ st_numscalar("sa_z")))

* ---- case wide_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
keep id t x z
reshape wide x z, i(id) j(t)
scalar ncols = c(k)
scalar nunits = _N
egen double tot = rowtotal(x* z*)
quietly summarize tot
scalar total = r(sum)
egen long nm = rownonmiss(x* z*)
quietly summarize nm
scalar nonmiss = r(sum)
mata: M = editmissing(st_data(., "x* z*"), 0)
mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))
mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))
mata: pyx_export("../out/wide_one_obs_unit__counts.csv", (st_numscalar("ncols") \ st_numscalar("nunits") \ st_numscalar("total") \ st_numscalar("nonmiss") \ st_numscalar("colw") \ st_numscalar("roww")))

* ---- case gaps_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
quietly tsreport, panel
scalar n_gaps = r(N_gaps)
scalar n_obs = r(N_obs)
mata: pyx_export("../out/gaps_one_obs_unit__counts.csv", (st_numscalar("n_gaps") \ st_numscalar("n_obs")))

* ---- case variation_re_effects
use "../data/sim_re_effects.dta", clear
preserve
quietly keep if !missing(x)
bysort id: egen double mn = min(x)
bysort id: egen double mx = max(x)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_x = r(N)
quietly count if u1 & cnt == 1
scalar us_x = r(N)
bysort t: egen double tmn = min(x)
bysort t: egen double tmx = max(x)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_x = r(N)
bysort x: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_x = r(N)
restore
preserve
quietly keep if !missing(z)
bysort id: egen double mn = min(z)
bysort id: egen double mx = max(z)
bysort id: gen byte u1 = (_n == 1)
bysort id: gen long cnt = _N
quietly count if u1 & mx > mn
scalar uv_z = r(N)
quietly count if u1 & cnt == 1
scalar us_z = r(N)
bysort t: egen double tmn = min(z)
bysort t: egen double tmx = max(z)
bysort t: gen byte p1 = (_n == 1)
quietly count if p1 & tmx > tmn
scalar pv_z = r(N)
bysort z: gen byte d1 = (_n == 1)
quietly count if d1
scalar nd_z = r(N)
restore
mata: pyx_export("../out/variation_re_effects__counts.csv", (st_numscalar("uv_x") \ st_numscalar("us_x") \ st_numscalar("pv_x") \ st_numscalar("nd_x") \ st_numscalar("uv_z") \ st_numscalar("us_z") \ st_numscalar("pv_z") \ st_numscalar("nd_z")))

* ---- case spells_re_effects
use "../data/sim_re_effects.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate byte ld = L.d
generate byte start = !missing(d) & (missing(ld) | ld != d)
generate byte change = !missing(d) & !missing(ld) & ld != d
quietly count if start
scalar n_spells = r(N)
quietly count if change
scalar n_changes = r(N)
bysort id: egen byte any_change = max(change)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & any_change == 1
scalar n_units_changing = r(N)
mata: pyx_export("../out/spells_re_effects__counts.csv", (st_numscalar("n_spells") \ st_numscalar("n_changes") \ st_numscalar("n_units_changing")))

* ---- case event_re_effects
use "../data/sim_re_effects.dta", clear
xtset id t
generate byte d = cond(missing(x), ., x > 0)
generate double on = cond(d == 1, t, .)
bysort id: egen double ev = min(on)
generate double rel = (t - ev) / 1
quietly summarize rel
scalar n_rel = r(N)
scalar sum_rel = r(sum)
generate double rel2 = rel^2
quietly summarize rel2
scalar sumsq_rel = r(sum)
bysort id: gen byte u1 = (_n == 1)
quietly count if u1 & !missing(ev)
scalar n_treated = r(N)
mata: pyx_export("../out/event_re_effects__counts.csv", (st_numscalar("n_rel") \ st_numscalar("sum_rel") \ st_numscalar("sumsq_rel") \ st_numscalar("n_treated")))

* ---- case outliers_re_effects
use "../data/sim_re_effects.dta", clear
bysort id: egen long n_x = count(x)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_x = pctile(x), p(50)
generate double dev_x = abs(x - med_x)
bysort id: egen double mad_x = pctile(dev_x), p(50)
generate double z_x = 0.6745 * (x - med_x) / mad_x if n_x >= 3 & mad_x > 0 & !missing(x)
generate double az_x = abs(z_x)
quietly summarize az_x
scalar ns_x = r(N)
scalar sa_x = r(sum)
quietly count if az_x > 2 & !missing(az_x)
scalar nf_x = r(N)
bysort id: egen long n_z = count(z)
* egen median() averages the two middle values in float precision; pctile() does not
bysort id: egen double med_z = pctile(z), p(50)
generate double dev_z = abs(z - med_z)
bysort id: egen double mad_z = pctile(dev_z), p(50)
generate double z_z = 0.6745 * (z - med_z) / mad_z if n_z >= 3 & mad_z > 0 & !missing(z)
generate double az_z = abs(z_z)
quietly summarize az_z
scalar ns_z = r(N)
scalar sa_z = r(sum)
quietly count if az_z > 2 & !missing(az_z)
scalar nf_z = r(N)
mata: pyx_export("../out/outliers_re_effects__counts.csv", (st_numscalar("ns_x") \ st_numscalar("nf_x") \ st_numscalar("sa_x") \ st_numscalar("ns_z") \ st_numscalar("nf_z") \ st_numscalar("sa_z")))

* ---- case wide_re_effects
use "../data/sim_re_effects.dta", clear
keep id t x z
reshape wide x z, i(id) j(t)
scalar ncols = c(k)
scalar nunits = _N
egen double tot = rowtotal(x* z*)
quietly summarize tot
scalar total = r(sum)
egen long nm = rownonmiss(x* z*)
quietly summarize nm
scalar nonmiss = r(sum)
mata: M = editmissing(st_data(., "x* z*"), 0)
mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))
mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))
mata: pyx_export("../out/wide_re_effects__counts.csv", (st_numscalar("ncols") \ st_numscalar("nunits") \ st_numscalar("total") \ st_numscalar("nonmiss") \ st_numscalar("colw") \ st_numscalar("roww")))

* ---- case gaps_re_effects
use "../data/sim_re_effects.dta", clear
xtset id t
quietly tsreport, panel
scalar n_gaps = r(N_gaps)
scalar n_obs = r(N_obs)
mata: pyx_export("../out/gaps_re_effects__counts.csv", (st_numscalar("n_gaps") \ st_numscalar("n_obs")))

exit, clear
