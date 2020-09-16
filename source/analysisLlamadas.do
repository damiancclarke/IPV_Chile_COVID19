/*analysisDenuncias.do v1.00     damiancclarke             yyyy-mm-dd:2020-08-27
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This code generates summary statistics, two-way FE models, and panel event study
analyses of the impact of quarantine imposition on phone calls to the national 
domestic violence hotline of Chile run by the Ministry for Women (MMEG). It req-
uires the following ados installed from the SSC.
 
 > eventdd
 > did_multiplegt
 > estout

All globals are set in section 0 to import data and export results and log files.
Replicating results should simply require altering these globals to addresses on 
your own computer, and installing the above mentioned ados.

*/

vers 15
clear all
set more off
cap log close

*-------------------------------------------------------------------------------
*--- (0) Globals, set-up
*-------------------------------------------------------------------------------
global DAT "..\C-19-Llamadas\llamadasExcel"
global OUT "..\Results\llamadas1455"
global LOG "..\Results\log"

log using "$LOG\analysisLlamadas.txt", text replace
cap mkdir $OUT

*-------------------------------------------------------------------------------
*--- (1) Open data
*-------------------------------------------------------------------------------
use "$DAT\llamadas1455", clear
drop if anio==2020&mes==1
gen callRate = Llamados/Population*100000

gen quarantine = Cuarentena>0 if Cuarentena!=.
drop if Region==.
bys Region (t): gen n = _n
bys Region (t): egen minn = min(n) if quarantine==1
bys Region: egen qstart = min(minn)
gen timeToQ = n-qstart
tab timeToQ
**DANIEL: Lo mismo que con SPD...  Colocamos promiedios de la 'linea base'
replace mobility_externo=0 if mobility_externo==. //REVISAR
replace mobility_interno=0 if mobility_interno==. //REVISAR

*-------------------------------------------------------------------------------
*--- (2) Event Study
*-------------------------------------------------------------------------------
foreach wt in no yes {
    if "`wt'"=="no"  {
    local opt
    local gn 
    }
	if "`wt'"=="yes" {
    local opt [aw=Population]
    local gn _Wt 
    }
    #delimit ;
    eventdd callRate i.t i.Region `opt', timevar(timeToQ) ci(rcap) lags(16) leads(3)
    wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
    graph_op(xlabel(-16 "{&le} -16" -12 "-12" -8 "-8" -4 "-4" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Calls to #1455 per 100,000 people"));
    graph export "$OUT\eventdd\event1`gn'_callRate.eps", replace;
	
    eventdd callRate i.t i.Region Population `opt', timevar(timeToQ) ci(rcap) lags(16) leads(3)
    wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
    graph_op(xlabel(-16 "{&le} -16" -12 "-12" -8 "-8" -4 "-4" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Calls to #1455 per 100,000 people"));
    graph export "$OUT\eventdd\event2`gn'_callRate.eps", replace;

    eventdd callRate i.t i.Region mobility_externo mobility_interno `opt',
    timevar(timeToQ) ci(rcap) lags(16) leads(3)
    wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
    graph_op(xlabel(-16 "{&le} -16" -12 "-12" -8 "-8" -4 "-4" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Calls to #1455 per 100,000 people"));
    graph export "$OUT\eventdd\event3`gn'_callRate.eps", replace;
	
    eventdd callRate i.t i.Region Population mobility_externo mobility_interno `opt',
    timevar(timeToQ) ci(rcap) lags(16) leads(3)
    wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
    graph_op(xlabel(-16 "{&le} -16" -12 "-12" -8 "-8" -4 "-4" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Calls to #1455 per 100,000 people"));
    graph export "$OUT\eventdd\event4`gn'_callRate.eps", replace;
    #delimit cr
}
graph drop _all

*-------------------------------------------------------------------------------
*--- (3) Two way FEs
*-------------------------------------------------------------------------------
xtset Region t
local ses abs(Region) cluster(Region) 

foreach q of varlist quarantine PropPopQuar {
    *Sin Pesos de Poblacion
    eststo: areg callRate i.t i.Region `q', `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate i.t i.Region Population `q', `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
    
    eststo: areg callRate i.t i.Region mobility_externo mobility_interno `q', `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate i.t i.Region mobility_externo mobility_interno Population `q', `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    *Con Pesos de Poblacion
    eststo: areg callRate i.t i.Region `q' [aw=Population], `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate i.t i.Region Population `q' [aw=Population], `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate i.t i.Region mobility_externo mobility_interno `q' [aw=Population], `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate i.t i.Region mobility_externo mobility_interno Population `q' [aw=Population], `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'

    #delimit ;
    esttab est1 est2 est3 est4 est5 est6 est7 est8 using "$OUT\areg\DD_`q'_`var'.tex",
    b(%-9.3f) se(%-9.3f) noobs keep(`q') nonotes nogaps
    mlabels(, none) nonumbers style(tex) fragment replace noline label
    starlevel ("*" 0.10 "**" 0.05 "***" 0.01);
    #delimit cr
    estimates clear
}
graph drop _all

*-------------------------------------------------------------------------------
*--- (4) Sharp Difference-in-Difference
*-------------------------------------------------------------------------------
foreach wt in no yes {
    if "`wt'"=="no"  {
        local opt
        local gn 
    }
    if "`wt'"=="yes" {
        local opt weight(Population)
        local gn _Wt 
    }
    foreach tv in quarantine PropPopQuar{
        if "`tv'"=="quarantine"  local v _quar
        if "`tv'"=="PropPopQuar" local v _pop
        did_multiplegt callRate Region t `tv', placebo(5) dynamic(3) breps(10) cluster(Region) `opt'
        ereturn list
        graph export "$OUT\did_multiplegt\did1`v'`gn'_callRate.eps", replace
        
        did_multiplegt callRate Region t `tv', placebo(5) dynamic(3) breps(10) cluster(Region) controls(Population) `opt'
        ereturn list
        graph export "$OUT\did_multiplegt\did2`v'`gn'_callRate.eps", replace
        
        did_multiplegt callRate Region t `tv', placebo(5) dynamic(3) breps(10) cluster(Region) controls(mobility_externo mobility_interno) `opt'
        ereturn list
        graph export "$OUT\did_multiplegt\did3`v'`gn'_callRate.eps", replace
        
        did_multiplegt callRate Region t `tv', placebo(5) dynamic(3) breps(10) cluster(Region) controls(mobility_externo mobility_interno Population) `opt'
        ereturn list
        graph export "$OUT\did_multiplegt\did4`v'`gn'_callRate.eps", replace
    }
    graph drop _all
}


cap log close


