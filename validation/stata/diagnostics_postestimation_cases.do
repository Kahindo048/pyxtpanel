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

* ---- case ar_ab1_balanced
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan robust
scalar a1 = e(ar1)
scalar p1 = e(ar1p)
scalar a2 = e(ar2)
scalar p2 = e(ar2p)
mata: pyx_export("../out/ar_ab1_balanced__ar.csv", (st_numscalar("a1") \ st_numscalar("p1") \ st_numscalar("a2") \ st_numscalar("p2")))

* ---- case ar_ab2_balanced
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan robust twostep
scalar a1 = e(ar1)
scalar p1 = e(ar1p)
scalar a2 = e(ar2)
scalar p2 = e(ar2p)
mata: pyx_export("../out/ar_ab2_balanced__ar.csv", (st_numscalar("a1") \ st_numscalar("p1") \ st_numscalar("a2") \ st_numscalar("p2")))

* ---- case ar_ab2_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan robust twostep
scalar a1 = e(ar1)
scalar p1 = e(ar1p)
scalar a2 = e(ar2)
scalar p2 = e(ar2p)
mata: pyx_export("../out/ar_ab2_gaps__ar.csv", (st_numscalar("a1") \ st_numscalar("p1") \ st_numscalar("a2") \ st_numscalar("p2")))

* ---- case ar_ab2_entry_exit
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan robust twostep
scalar a1 = e(ar1)
scalar p1 = e(ar1p)
scalar a2 = e(ar2)
scalar p2 = e(ar2p)
mata: pyx_export("../out/ar_ab2_entry_exit__ar.csv", (st_numscalar("a1") \ st_numscalar("p1") \ st_numscalar("a2") \ st_numscalar("p2")))

* ---- case ar_abond_balanced
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtabond x, lags(1) noconstant vce(robust) twostep
quietly estat abond, artests(2)
mata: pyx_export("../out/ar_abond_balanced__arm.csv", st_matrix("r(arm)"))

* ---- case ar_abond1_balanced
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtabond x, lags(1) noconstant vce(robust)
quietly estat abond, artests(2)
mata: pyx_export("../out/ar_abond1_balanced__arm.csv", st_matrix("r(arm)"))

* ---- case diffhansen_balanced
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtabond2 x L.x z, gmm(x, lag(2 .)) iv(z) noleveleq noconstant twostep robust
mata: pyx_export("../out/diffhansen_balanced__jfull.csv", st_numscalar("e(hansen)"))
mata: pyx_export("../out/diffhansen_balanced__jexcl.csv", st_matrix("e(diffsargan)")[1,.])
mata: pyx_export("../out/diffhansen_balanced__stata.csv", st_matrix("e(diffsargan)")[(2 \ 3 \ 5),.])

* ---- case diffhansen_gaps_groups
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly xtabond2 x L.x z, gmm(x, lag(2 .)) gmm(z, lag(1 2)) iv(z) noleveleq noconstant twostep robust
mata: pyx_export("../out/diffhansen_gaps_groups__jfull.csv", st_numscalar("e(hansen)"))
mata: pyx_export("../out/diffhansen_gaps_groups__jexcl.csv", st_matrix("e(diffsargan)")[1,.])
mata: pyx_export("../out/diffhansen_gaps_groups__stata.csv", st_matrix("e(diffsargan)")[(2 \ 3 \ 5),.])

* ---- case hansen_ivregress
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly ivregress gmm y (x = z1 z2 z3 z4), wmatrix(cluster id)
quietly estat overid
scalar hj = r(HansenJ)
scalar hdf = r(J_df)
scalar hp = r(p_HansenJ)
mata: pyx_export("../out/hansen_ivregress__j.csv", (st_numscalar("hj") \ st_numscalar("hdf") \ st_numscalar("hp")))

* ---- case hansen_ivregress_invalid
use "../data/sim_iv_invalid.dta", clear
quietly xtset id t
quietly ivregress gmm y (x = z1 z2 z3 z4), wmatrix(cluster id)
quietly estat overid
scalar hj = r(HansenJ)
scalar hdf = r(J_df)
scalar hp = r(p_HansenJ)
mata: pyx_export("../out/hansen_ivregress_invalid__j.csv", (st_numscalar("hj") \ st_numscalar("hdf") \ st_numscalar("hp")))

* ---- case hansen_just_identified
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly ivregress gmm y (x = z1), wmatrix(cluster id)
capture estat overid
scalar refused = (_rc != 0)
mata: pyx_export("../out/hansen_just_identified__refused.csv", (st_numscalar("refused")))

* ---- case hausman_fe_re
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly xtreg y x z1 z2, fe
estimates store fe
mata: pyx_export("../out/hausman_fe_re__b_fe.csv", st_matrix("e(b)"))
mata: pyx_export("../out/hausman_fe_re__v_fe.csv", st_matrix("e(V)"))
quietly xtreg y x z1 z2, re
estimates store re
mata: pyx_export("../out/hausman_fe_re__b_re.csv", st_matrix("e(b)"))
mata: pyx_export("../out/hausman_fe_re__v_re.csv", st_matrix("e(V)"))
quietly hausman fe re
scalar hc = r(chi2)
scalar hdf = r(df)
scalar hp = r(p)
mata: pyx_export("../out/hausman_fe_re__h.csv", (st_numscalar("hc") \ st_numscalar("hdf") \ st_numscalar("hp")))

* ---- case hausman_fe_re_invalid
use "../data/sim_iv_invalid.dta", clear
quietly xtset id t
quietly xtreg y x z1 z2 z3, fe
estimates store fe
mata: pyx_export("../out/hausman_fe_re_invalid__b_fe.csv", st_matrix("e(b)"))
mata: pyx_export("../out/hausman_fe_re_invalid__v_fe.csv", st_matrix("e(V)"))
quietly xtreg y x z1 z2 z3, re
estimates store re
mata: pyx_export("../out/hausman_fe_re_invalid__b_re.csv", st_matrix("e(b)"))
mata: pyx_export("../out/hausman_fe_re_invalid__v_re.csv", st_matrix("e(V)"))
quietly hausman fe re
scalar hc = r(chi2)
scalar hdf = r(df)
scalar hp = r(p)
mata: pyx_export("../out/hausman_fe_re_invalid__h.csv", (st_numscalar("hc") \ st_numscalar("hdf") \ st_numscalar("hp")))

* ---- case hausman_swapped_negative
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly xtreg y x z1 z2, fe
estimates store fe
mata: pyx_export("../out/hausman_swapped_negative__b_fe.csv", st_matrix("e(b)"))
mata: pyx_export("../out/hausman_swapped_negative__v_fe.csv", st_matrix("e(V)"))
quietly xtreg y x z1 z2, re
estimates store re
mata: pyx_export("../out/hausman_swapped_negative__b_re.csv", st_matrix("e(b)"))
mata: pyx_export("../out/hausman_swapped_negative__v_re.csv", st_matrix("e(V)"))
quietly hausman re fe
scalar hc = r(chi2)
scalar hdf = r(df)
scalar hp = r(p)
mata: pyx_export("../out/hausman_swapped_negative__h.csv", (st_numscalar("hc") \ st_numscalar("hdf")))

* ---- case hausman_df_option
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly xtreg y x z1 z2, fe
estimates store fe
mata: pyx_export("../out/hausman_df_option__b_fe.csv", st_matrix("e(b)"))
mata: pyx_export("../out/hausman_df_option__v_fe.csv", st_matrix("e(V)"))
quietly xtreg y x z1 z2, re
estimates store re
mata: pyx_export("../out/hausman_df_option__b_re.csv", st_matrix("e(b)"))
mata: pyx_export("../out/hausman_df_option__v_re.csv", st_matrix("e(V)"))
quietly hausman fe re, df(2)
scalar hc = r(chi2)
scalar hdf = r(df)
scalar hp = r(p)
mata: pyx_export("../out/hausman_df_option__h.csv", (st_numscalar("hc") \ st_numscalar("hdf") \ st_numscalar("hp")))

* ---- case sargan_ivregress
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly ivregress 2sls y (x = z1 z2 z3 z4)
quietly estat overid
scalar ss = r(sargan)
scalar sdf = r(df)
scalar sp = r(p_sargan)
mata: pyx_export("../out/sargan_ivregress__s.csv", (st_numscalar("ss") \ st_numscalar("sdf") \ st_numscalar("sp")))

* ---- case sargan_ivregress_invalid
use "../data/sim_iv_invalid.dta", clear
quietly xtset id t
quietly ivregress 2sls y (x = z1 z2 z3 z4)
quietly estat overid
scalar ss = r(sargan)
scalar sdf = r(df)
scalar sp = r(p_sargan)
mata: pyx_export("../out/sargan_ivregress_invalid__s.csv", (st_numscalar("ss") \ st_numscalar("sdf") \ st_numscalar("sp")))

* ---- case sargan_just_identified
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly ivregress 2sls y (x = z1)
capture estat overid
scalar refused = (_rc != 0)
mata: pyx_export("../out/sargan_just_identified__refused.csv", (st_numscalar("refused")))

* ---- case wald_regress_f
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly regress y x z1 z2 z3, vce(cluster id)
mata: pyx_export("../out/wald_regress_f__b.csv", st_matrix("e(b)"))
mata: pyx_export("../out/wald_regress_f__v.csv", st_matrix("e(V)"))
quietly test x z1
scalar ws = r(F)
scalar wdf = r(df)
scalar wp = r(p)
scalar dfr = r(df_r)
mata: pyx_export("../out/wald_regress_f__dfr.csv", (st_numscalar("dfr")))
mata: pyx_export("../out/wald_regress_f__w.csv", (st_numscalar("ws") \ st_numscalar("wdf") \ st_numscalar("wp")))

* ---- case wald_xtreg_fe_cluster_f
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly xtreg y x z1 z2, fe vce(cluster id)
mata: pyx_export("../out/wald_xtreg_fe_cluster_f__b.csv", st_matrix("e(b)"))
mata: pyx_export("../out/wald_xtreg_fe_cluster_f__v.csv", st_matrix("e(V)"))
quietly test x z2
scalar ws = r(F)
scalar wdf = r(df)
scalar wp = r(p)
scalar dfr = r(df_r)
mata: pyx_export("../out/wald_xtreg_fe_cluster_f__dfr.csv", (st_numscalar("dfr")))
mata: pyx_export("../out/wald_xtreg_fe_cluster_f__w.csv", (st_numscalar("ws") \ st_numscalar("wdf") \ st_numscalar("wp")))

* ---- case wald_ivregress_chi2
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly ivregress 2sls y (x = z1 z2 z3) z4, vce(cluster id)
mata: pyx_export("../out/wald_ivregress_chi2__b.csv", st_matrix("e(b)"))
mata: pyx_export("../out/wald_ivregress_chi2__v.csv", st_matrix("e(V)"))
quietly test x z4
scalar ws = r(chi2)
scalar wdf = r(df)
scalar wp = r(p)
mata: pyx_export("../out/wald_ivregress_chi2__w.csv", (st_numscalar("ws") \ st_numscalar("wdf") \ st_numscalar("wp")))

* ---- case wald_nonzero_values
use "../data/sim_iv_valid.dta", clear
quietly xtset id t
quietly regress y x z1 z2 z3
mata: pyx_export("../out/wald_nonzero_values__b.csv", st_matrix("e(b)"))
mata: pyx_export("../out/wald_nonzero_values__v.csv", st_matrix("e(V)"))
quietly test (x = 0.8) (z1 = 0.1)
scalar ws = r(F)
scalar wdf = r(df)
scalar wp = r(p)
scalar dfr = r(df_r)
mata: pyx_export("../out/wald_nonzero_values__dfr.csv", (st_numscalar("dfr")))
mata: pyx_export("../out/wald_nonzero_values__w.csv", (st_numscalar("ws") \ st_numscalar("wdf") \ st_numscalar("wp")))

* ---- case sargan_ab1_balanced
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan
scalar st = e(sargan)
scalar sdf = e(j) - colsof(e(b))
scalar sp = e(sarganp)
mata: pyx_export("../out/sargan_ab1_balanced__s.csv", (st_numscalar("st") \ st_numscalar("sdf") \ st_numscalar("sp")))

* ---- case sargan_ab2_balanced
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan twostep
scalar st = e(sargan)
scalar sdf = e(j) - colsof(e(b))
scalar sp = e(sarganp)
mata: pyx_export("../out/sargan_ab2_balanced__s.csv", (st_numscalar("st") \ st_numscalar("sdf") \ st_numscalar("sp")))

* ---- case hansen_ab2_balanced
use "../data/sim_balanced.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan twostep
scalar st = e(hansen)
scalar sdf = e(hansen_df)
scalar sp = e(hansenp)
mata: pyx_export("../out/hansen_ab2_balanced__s.csv", (st_numscalar("st") \ st_numscalar("sdf") \ st_numscalar("sp")))

* ---- case sargan_ab2_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan twostep
scalar st = e(sargan)
scalar sdf = e(j) - colsof(e(b))
scalar sp = e(sarganp)
mata: pyx_export("../out/sargan_ab2_gaps__s.csv", (st_numscalar("st") \ st_numscalar("sdf") \ st_numscalar("sp")))

* ---- case hansen_ab2_gaps
use "../data/sim_gaps.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan twostep
scalar st = e(hansen)
scalar sdf = e(hansen_df)
scalar sp = e(hansenp)
mata: pyx_export("../out/hansen_ab2_gaps__s.csv", (st_numscalar("st") \ st_numscalar("sdf") \ st_numscalar("sp")))

* ---- case sargan_ab2_entry_exit
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan twostep
scalar st = e(sargan)
scalar sdf = e(j) - colsof(e(b))
scalar sp = e(sarganp)
mata: pyx_export("../out/sargan_ab2_entry_exit__s.csv", (st_numscalar("st") \ st_numscalar("sdf") \ st_numscalar("sp")))

* ---- case hansen_ab2_entry_exit
use "../data/sim_entry_exit.dta", clear
quietly xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan twostep
scalar st = e(hansen)
scalar sdf = e(hansen_df)
scalar sp = e(hansenp)
mata: pyx_export("../out/hansen_ab2_entry_exit__s.csv", (st_numscalar("st") \ st_numscalar("sdf") \ st_numscalar("sp")))

exit, clear
