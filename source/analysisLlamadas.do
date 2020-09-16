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
global ROOT "C:/Users/danie/Desktop/Proyectos/Asistente/replication/"

global DAT "$ROOT/data/calls/dta"
global OUT "$ROOT/results/calls"
global LOG "$ROOT/log"


cap mkdir "$OUT"
cap mkdir "$LOG"
log using "$LOG/analysisLlamadas.txt", text replace

foreach folder in eventdd areg did_multiplegt descriptives {
    cap mkdir "$OUT/`folder'"
}
foreach ado in eventdd did_multiplegt estout {
    cap which `ado'
    if _rc!=0 ssc install `ado'
}
*-------------------------------------------------------------------------------
*--- (1) Open data
*-------------------------------------------------------------------------------
use "$DAT/llamadas1455", clear
drop if anio==2020&mes==1
gen callRate = Llamados/Population*100000

gen quarantine = Cuarentena>0 if Cuarentena!=.
drop if Region==.
bys Region (t): gen n = _n
bys Region (t): egen minn = min(n) if quarantine==1
bys Region: egen qstart = min(minn)
gen timeToQ = n-qstart
tab timeToQ

levelsof Region, local(region)
foreach tp in externo interno {
    foreach c of local region {
        qui sum mobility_`tp' if t<=721&Region==`c' //solo se repite el dato de febrero
        replace mobility_`tp'=`r(mean)' if t<721&Region==`c'
    }
}
*-------------------------------------------------------------------------------
*--- (2) Event Study
*-------------------------------------------------------------------------------
local fes i.t i.Region

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
    eventdd callRate `fes' `opt', timevar(timeToQ) ci(rcap) lags(16) leads(3)
    wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
    graph_op(xlabel(-16 "{&le} -16" -12 "-12" -8 "-8" -4 "-4" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Calls to #1455 per 100,000 people"));
    graph export "$OUT/eventdd/event1`gn'_callRate.eps", replace;
	
    eventdd callRate `fes' Population `opt', timevar(timeToQ) ci(rcap) lags(16) leads(3)
    wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
    graph_op(xlabel(-16 "{&le} -16" -12 "-12" -8 "-8" -4 "-4" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Calls to #1455 per 100,000 people"));
    graph export "$OUT/eventdd/event2`gn'_callRate.eps", replace;

    eventdd callRate `fes' mobility_externo mobility_interno `opt',
    timevar(timeToQ) ci(rcap) lags(16) leads(3)
    wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
    graph_op(xlabel(-16 "{&le} -16" -12 "-12" -8 "-8" -4 "-4" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Calls to #1455 per 100,000 people"));
    graph export "$OUT/eventdd/event3`gn'_callRate.eps", replace;
	
    eventdd callRate `fes' Population mobility_externo mobility_interno `opt',
    timevar(timeToQ) ci(rcap) lags(16) leads(3)
    wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
    graph_op(xlabel(-16 "{&le} -16" -12 "-12" -8 "-8" -4 "-4" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Calls to #1455 per 100,000 people"));
    graph export "$OUT/eventdd/event4`gn'_callRate.eps", replace;
	#delimit cr
}
graph drop _all

*-------------------------------------------------------------------------------
*--- (3) Two way FEs
*-------------------------------------------------------------------------------
xtset Region t
local fes i.t i.Region
local ses abs(Region) cluster(Region) 

foreach q of varlist quarantine PropPopQuar {
    *Sin Pesos de Poblacion
    eststo: areg callRate `fes' `q', `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate `fes' Population `q', `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
    
    eststo: areg callRate `fes' mobility_externo mobility_interno `q', `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate `fes' mobility_externo mobility_interno Population `q', `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    *Con Pesos de Poblacion
    eststo: areg callRate `fes' `q' [aw=Population], `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate `fes' Population `q' [aw=Population], `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate `fes' mobility_externo mobility_interno `q' [aw=Population], `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
    eststo: areg callRate `fes' mobility_externo mobility_interno Population `q' [aw=Population], `ses'
    sum callRate if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
 
    local ests est1 est2 est3 est4 est5 est6 est7 est8
    #delimit ;
    esttab `ests' using "$OUT/areg/DD_`q'_callRate.tex", b(%-9.3f) se(%-9.3f) 
	noobs keep(mobility_externo mobility_interno Population `q') nonotes nogaps
    mlabels(, none) nonumbers style(tex) fragment replace noline label
    starlevel ("*" 0.10 "**" 0.05 "***" 0.01);
    #delimit cr
    estimates clear
}

*-------------------------------------------------------------------------------
*--- (4) Sharp Difference-in-Difference
*-------------------------------------------------------------------------------
local fes Region t

foreach wt in no yes {
    foreach tv in quarantine PropPopQuar {
	if "`wt'"=="no"  {
            if "`tv'"=="quarantine" {  
                local v _quar
	            local opt placebo(5) dynamic(3) breps(10) cluster(Region)
                local gn 
	    }
            if "`tv'"=="PropPopQuar" {
                local v _pop
                local opt placebo(5) dynamic(1) breps(10) cluster(Region)
                local gn 
            }
        }
        if "`wt'"=="yes" {
            if "`tv'"=="quarantine" {  
                local v _quar
	            local opt placebo(5) dynamic(3) breps(10) cluster(Region) weight(Population)
                local gn _Wt
            }
            if "`tv'"=="PropPopQuar" {
                local v _pop
                local opt placebo(5) dynamic(1) breps(10) cluster(Region) weight(Population)
                local gn _Wt
            }			
        }
		
	did_multiplegt callRate `fes' `tv', `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did1`v'`gn'_callRate.eps", replace
        
        did_multiplegt callRate `fes' `tv', controls(Population) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did2`v'`gn'_callRate.eps", replace
        
        did_multiplegt callRate `fes' `tv', controls(mobility_externo mobility_interno) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did3`v'`gn'_callRate.eps", replace
        
        did_multiplegt callRate `fes' `tv', controls(mobility_externo mobility_interno Population) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did4`v'`gn'_callRate.eps", replace
    }
    graph drop _all
}

cap log close
