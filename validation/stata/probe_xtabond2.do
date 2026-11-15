* Probe: sigma_u and sigma_e of xtreg, re on the re_effects panel
version 17
clear all
use "../data/sim_re_effects.dta", clear
xtset id t
quietly xtreg x z, re
display %21.15g e(sigma_u) "|" %21.15g e(sigma_e)
ereturn list
