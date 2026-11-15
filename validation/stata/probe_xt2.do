set more off
display "=== H1 missing id"
clear
input id t x
1 1 1
1 2 2
. 1 3
2 1 4
2 2 5
end
capture noisily xtset id t
display "rc=" _rc
return list
display "=== H2 missing t"
clear
input id t x
1 1 1
1 2 2
1 . 3
2 1 4
2 2 5
end
capture noisily xtset id t
display "rc=" _rc
xtdescribe
display "=== H3 non integer t"
clear
input id t x
1 1 1
1 2.5 2
end
capture noisily xtset id t
display "rc=" _rc
display "=== H4 string id"
clear
input str3 id t x
a 1 1
a 2 2
b 1 3
end
capture noisily xtset id t
display "rc=" _rc
display "=== H5 duplicates"
clear
input id t x
1 1 1
1 1 2
end
capture noisily xtset id t
display "rc=" _rc
display "=== H6 delta(2) + off-grid"
clear
input id t x
1 1 1
1 3 2
1 7 3
2 1 4
2 3 5
end
xtset id t, delta(2)
return list
xtdescribe
capture noisily xtset id t, delta(3)
display "rc=" _rc
display "=== H7 quarterly dates"
clear
input id q x
1 1 1
1 2 2
1 4 3
2 1 4
end
gen tq = yq(2020,1) + q - 1
format tq %tq
xtset id tq
return list
xtdescribe
display "=== H8 time-series operators with gap"
xtset id tq
gen lx = L.x
gen dx = D.x
gen l2x = L2.x
gen sx = S2.x
gen fx = F.x
list id tq x lx dx l2x sx fx
display "=== H9 one-id panel / single period"
clear
input id t x
1 1 1
1 2 2
end
xtset id t
xtdescribe
xtsum x
display "=== H10 id one obs unit in xtsum"
clear
input id t x
1 1 1
2 1 4
2 2 5
end
xtset id t
xtsum x
