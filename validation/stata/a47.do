* A.4.7 (2): where the available commands come from, and their version
foreach c in xtunitroot xtcointtest xtserial xtcsd xtcd2 xttest0 xttest3 xthst pescadf ivreg2 ranktest xtabond2 {
    display "---- `c'"
    which `c'
}
* help of xtunitroot / xtcointtest subcommands as installed
capture ssc describe xtgcause
display "ssc rc = " _rc
