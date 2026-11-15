* L.8.3-8.6 feasibility
use "../data/sim_iv_valid.dta", clear
xtset id t
quietly xtreg y x z1 z2, fe
display "F_f=" e(F_f) " df_a=" e(df_a) " df_r=" e(df_r) " rss=" e(rss) " N=" e(N) " Ng=" e(N_g) " p=" Ftail(e(N_g)-1, e(df_r), e(F_f))
ereturn list
quietly regress y x z1 z2
display "pooled rss=" e(rss)
foreach d in sim_balanced sim_entry_exit sim_gaps {
    use "../data/`d'.dta", clear
    xtset id t
    quietly xtreg x z, re
    display "`d' sigma_u=" e(sigma_u) " sigma_e=" e(sigma_e)
    quietly xttest0
    display "   LM=" r(lm) " p=" r(p)
}
use "../data/sim_het.dta", clear
quietly regress y x z1 z2
quietly predict double e, resid
sktest e
return list
quietly summarize e, detail
display "skew=" r(skewness) " kurt=" r(kurtosis) " N=" r(N)
