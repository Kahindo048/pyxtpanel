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

* ---- case xttab_balanced_dz
use "../data/sim_balanced.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_balanced_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_balanced_dz
use "../data/sim_balanced.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_balanced_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_balanced_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_balanced_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_balanced_dz__cols.csv", st_matrix("C")')

* ---- case xttab_balanced_dx
use "../data/sim_balanced.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_balanced_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_balanced_dx
use "../data/sim_balanced.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_balanced_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_balanced_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_balanced_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_balanced_dx__cols.csv", st_matrix("C")')

* ---- case xttab_gaps_dz
use "../data/sim_gaps.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_gaps_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_gaps_dz
use "../data/sim_gaps.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_gaps_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_gaps_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_gaps_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_gaps_dz__cols.csv", st_matrix("C")')

* ---- case xttab_gaps_dx
use "../data/sim_gaps.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_gaps_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_gaps_dx
use "../data/sim_gaps.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_gaps_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_gaps_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_gaps_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_gaps_dx__cols.csv", st_matrix("C")')

* ---- case xttab_entry_exit_dz
use "../data/sim_entry_exit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_entry_exit_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_entry_exit_dz
use "../data/sim_entry_exit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_entry_exit_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_entry_exit_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_entry_exit_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_entry_exit_dz__cols.csv", st_matrix("C")')

* ---- case xttab_entry_exit_dx
use "../data/sim_entry_exit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_entry_exit_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_entry_exit_dx
use "../data/sim_entry_exit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_entry_exit_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_entry_exit_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_entry_exit_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_entry_exit_dx__cols.csv", st_matrix("C")')

* ---- case xttab_weak_dz
use "../data/sim_weak.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_weak_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_weak_dz
use "../data/sim_weak.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_weak_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_weak_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_weak_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_weak_dz__cols.csv", st_matrix("C")')

* ---- case xttab_weak_dx
use "../data/sim_weak.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_weak_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_weak_dx
use "../data/sim_weak.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_weak_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_weak_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_weak_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_weak_dx__cols.csv", st_matrix("C")')

* ---- case xttab_delta2_dz
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_delta2_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_delta2_dz
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_delta2_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_delta2_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_delta2_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_delta2_dz__cols.csv", st_matrix("C")')

* ---- case xttab_delta2_dx
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_delta2_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_delta2_dx
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_delta2_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_delta2_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_delta2_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_delta2_dx__cols.csv", st_matrix("C")')

* ---- case xttab_quarterly_dz
use "../data/sim_quarterly.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_quarterly_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_quarterly_dz
use "../data/sim_quarterly.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_quarterly_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_quarterly_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_quarterly_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_quarterly_dz__cols.csv", st_matrix("C")')

* ---- case xttab_quarterly_dx
use "../data/sim_quarterly.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_quarterly_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_quarterly_dx
use "../data/sim_quarterly.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_quarterly_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_quarterly_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_quarterly_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_quarterly_dx__cols.csv", st_matrix("C")')

* ---- case xttab_unsorted_ids_dz
use "../data/sim_unsorted_ids.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_unsorted_ids_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_unsorted_ids_dz
use "../data/sim_unsorted_ids.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_unsorted_ids_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_unsorted_ids_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_unsorted_ids_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_unsorted_ids_dz__cols.csv", st_matrix("C")')

* ---- case xttab_unsorted_ids_dx
use "../data/sim_unsorted_ids.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_unsorted_ids_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_unsorted_ids_dx
use "../data/sim_unsorted_ids.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_unsorted_ids_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_unsorted_ids_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_unsorted_ids_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_unsorted_ids_dx__cols.csv", st_matrix("C")')

* ---- case xttab_one_obs_unit_dz
use "../data/sim_one_obs_unit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_one_obs_unit_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_one_obs_unit_dz
use "../data/sim_one_obs_unit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_one_obs_unit_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_one_obs_unit_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_one_obs_unit_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_one_obs_unit_dz__cols.csv", st_matrix("C")')

* ---- case xttab_one_obs_unit_dx
use "../data/sim_one_obs_unit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_one_obs_unit_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_one_obs_unit_dx
use "../data/sim_one_obs_unit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_one_obs_unit_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_one_obs_unit_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_one_obs_unit_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_one_obs_unit_dx__cols.csv", st_matrix("C")')

* ---- case xttab_single_unit_dz
use "../data/sim_single_unit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_single_unit_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_single_unit_dz
use "../data/sim_single_unit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_single_unit_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_single_unit_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_single_unit_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_single_unit_dz__cols.csv", st_matrix("C")')

* ---- case xttab_single_unit_dx
use "../data/sim_single_unit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_single_unit_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_single_unit_dx
use "../data/sim_single_unit.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_single_unit_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_single_unit_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_single_unit_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_single_unit_dx__cols.csv", st_matrix("C")')

* ---- case xttab_re_effects_dz
use "../data/sim_re_effects.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dz
mata: pyx_export("../out/xttab_re_effects_dz__results.csv", st_matrix("r(results)"))

* ---- case xttrans_re_effects_dz
use "../data/sim_re_effects.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dz
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_re_effects_dz__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dz if _n < _N
quietly by id: generate float is = dz[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dz)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_re_effects_dz__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_re_effects_dz__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_re_effects_dz__cols.csv", st_matrix("C")')

* ---- case xttab_re_effects_dx
use "../data/sim_re_effects.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttab dx
mata: pyx_export("../out/xttab_re_effects_dx__results.csv", st_matrix("r(results)"))

* ---- case xttrans_re_effects_dx
use "../data/sim_re_effects.dta", clear
xtset id t
generate byte dz = 0 if z < -0.5
replace dz = 1 if z >= -0.5 & z < 0.5
replace dz = 2 if z >= 0.5 & !missing(z)
generate byte dx = 0 if x < -1
replace dx = 1 if x >= -1 & x < 0
replace dx = 2 if x >= 0 & x < 1
replace dx = 3 if x >= 1 & !missing(x)
xttrans dx
scalar s_ntrans = r(n_trans)
scalar s_nrows = r(n_rows)
scalar s_ncols = r(n_cols)
mata: pyx_export("../out/xttrans_re_effects_dx__r.csv", (st_numscalar("s_ntrans") \ st_numscalar("s_nrows") \ st_numscalar("s_ncols")))
sort id t
quietly by id: generate float was = dx if _n < _N
quietly by id: generate float is = dx[_n + 1] if _n < _N
quietly generate byte touse = !missing(id, t, dx)
quietly by id: replace touse = 0 if touse[_n + 1] == 0 & _n < _N
quietly tabulate was is if touse, matcell(M) matrow(R) matcol(C)
mata: pyx_export("../out/xttrans_re_effects_dx__cells.csv", st_matrix("M"))
mata: pyx_export("../out/xttrans_re_effects_dx__rows.csv", st_matrix("R"))
mata: pyx_export("../out/xttrans_re_effects_dx__cols.csv", st_matrix("C")')

exit, clear
