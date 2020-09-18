/*analysisDenuncias.do v1.00     damiancclarke             yyyy-mm-dd:2020-08-31
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This code generates summary statistics, two-way FE models, and panel event study
analyses of the impact of quarantine imposition on "Casa de Acogida". It req-
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
set matsize 3000

*-------------------------------------------------------------------------------
*--- (0) Globals, set-up
*-------------------------------------------------------------------------------
global ROOT "/home/damian/investigacion/2020/IPV_COVID19/replication"

global DAT "$ROOT/data/CDA/dta"
global OUT "$ROOT/results/cda"
global LOG "$ROOT/log"

log using "$LOG/analysisCDA.txt", text replace
cap mkdir $OUT

*-------------------------------------------------------------------------------
*--- (1) Open data and creating new variables
*-------------------------------------------------------------------------------
use "$DAT/CDA.dta", clear

drop if t==21996 //21mar
drop if t==21995 //22mar

gen ingresspc = (ingreso_m+ingreso_n)/population*100000
gen occupancypc  = (ingreso_m/(ingreso_m+cupos_m))/population*100000

local outcomespc ingresspc occupancypc

bys Region (t): gen n = _n
bys Region (t): egen minn = min(n) if quarantine==1
bys Region: egen qstart = min(minn)
gen timeToQ = n-qstart
tab timeToQ

**DANIEL; AQUI LA MISMA SUGERENCIA DE SIEMPRE-  GRACIAS!!!
replace mobility_externo=0 if mobility_externo==. //REVISAR
replace mobility_interno=0 if mobility_interno==. //REVISAR

*-------------------------------------------------------------------------------
*--- (2) Event Study
*-------------------------------------------------------------------------------
***DANIEL: ¿Podemos tener 2 versiones de event studies aquí: 1 a nivel diario, y
*** otro a nivel semanal?  En la versión diario de repente agregamos más lags y
*** leads (50 lags, 50 leads, si es que da), y en la versión semanal, 10 lags y
*** 8 leads (si es que da).
                     
foreach wt in no yes {
    if "`wt'"=="no"  {
        local opt
        local gn 
    }
    if "`wt'"=="yes" {
        local opt [aw=population]
        local gn _Wt 
    }
    foreach var of varlist `outcomespc' {
	if "`var'"=="ingresspc"   local et="Ingress of Residents per 100,000 people"
	if "`var'"=="occupancypc" local et="Occupancy of Residents per 100,000 people"
        
        #delimit ;
        eventdd `var' i.t i.Region `opt', timevar(timeToQ) ci(rcap) lags(20) leads(3)
        wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-20 "{&le} -20" -15 "-15" -10 "-10" -5 "-5" -1 "-1"
                        0 "0" 1 "+1" 2 "+2" 3 "+3")
                 scheme(s1mono) xtitle("Days Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event1`gn'_`var'.eps", replace;
    	
        eventdd `var' i.t i.Region population `opt', timevar(timeToQ) ci(rcap) lags(20) leads(3)
        wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-20 "{&le} -20" -15 "-15" -10 "-10" -5 "-5" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
                 scheme(s1mono) xtitle("Days Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event2`gn'_`var'.eps", replace;
    
    	eventdd `var' i.t i.Region mobility_externo mobility_interno `opt',
        timevar(timeToQ) ci(rcap) lags(20) leads(3)
        wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-20 "{&le} -20" -15 "-15" -10 "-10" -5 "-5" -1 "-1"
                        0 "0" 1 "+1" 2 "+2" 3 "+3")
                 scheme(s1mono) xtitle("Days Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event3`gn'_`var'.eps", replace;
    	
    	eventdd `var' i.t i.Region population mobility_externo mobility_interno `opt', timevar(timeToQ) ci(rcap) lags(20) leads(3)
        wboot baseline(-1) coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-20 "{&le} -20" -15 "-15" -10 "-10" -5 "-5" -1 "-1" 0 "0" 1 "+1" 2 "+2" 3 "+3")
                 scheme(s1mono) xtitle("Days Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event4`gn'_`var'.eps", replace;
    	#delimit cr
    }
    graph drop _all
}

*-------------------------------------------------------------------------------
*--- (3) Two way FEs
*-------------------------------------------------------------------------------
**DANIEL, ¿Aquí hay alguna manera de automatizar la inclusión de los intervalos
** de confianza de boottest en la tabla de regresión?

xtset Region t
local se abs(Region) cluster(Region)

foreach q of varlist quarantine PropPopQuar {
    foreach var of varlist `outcomespc' {
        *Sin Pesos de Poblacion
        eststo: areg `var' i.t i.Region `q', `se' 
	sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
	eststo: areg `var' i.t i.Region population `q', `se'
	sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'

        
	eststo: areg `var' i.t i.Region mobility_externo mobility_interno `q', `se'
	sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
	eststo: areg `var' i.t i.Region mobility_externo mobility_interno population `q', `se'
	sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
        *Con Pesos de Poblacion
	eststo: areg `var' i.t i.Region `q' [aw=population], `se'
	sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
	eststo: areg `var' i.t i.Region population `q' [aw=population], `se'
	sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
	eststo: areg `var' i.t i.Region mobility_externo mobility_interno `q' [aw=population], `se'
	sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
	eststo: areg `var' i.t i.Region mobility_externo mobility_interno population `q' [aw=population], `se'
	sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
        
	#delimit ;
	esttab est1 est2 est3 est4 est5 est6 est7 est8 using "$OUT/areg/DD_`q'_`var'.tex",
	b(%-9.3f) se(%-9.3f) noobs keep(`q') nonotes nogaps
	mlabels(, none) nonumbers style(tex) fragment replace noline label
	starlevel ("*" 0.10 "**" 0.05 "***" 0.01);
	#delimit cr
	estimates clear
    }
    graph drop _all
}

*-------------------------------------------------------------------------------
*--- (4) Sharp Difference-in-Difference
*-------------------------------------------------------------------------------
***DANIEL ¿Aquí podemos utilizar datos semanales en vez de diarios?
foreach wt in no yes {
    if "`wt'"=="no"  {
        local opt
        local gn 
    }
    if "`wt'"=="yes" {
        local opt weight(population)
        local gn _Wt 
    }
    foreach tv in quarantine PropPopQuar{
        if "`tv'"=="quarantine"  local v _quar
        if "`tv'"=="PropPopQuar" local v _pop
        foreach var of varlist `outcomespc' {
            did_multiplegt `var' Region t `tv', placebo(5) dynamic(3) breps(50) cluster(Region) `opt'
            ereturn list
            graph export "$OUT/did_multiplegt/did1`v'`gn'_`var'.eps", replace
            
            did_multiplegt `var' Region t `tv', placebo(5) dynamic(3) breps(50) cluster(Region) controls(population) `opt'
            ereturn list
            graph export "$OUT/did_multiplegt/did2`v'`gn'_`var'.eps", replace
            
            did_multiplegt `var' Region t `tv', placebo(3) dynamic(3) breps(50) cluster(Region) controls(mobility_externo mobility_interno) `opt'
            ereturn list
            graph export "$OUT/did_multiplegt/did3`v'`gn'_`var'.eps", replace
            
            did_multiplegt `var' Region t `tv', placebo(3) dynamic(3) breps(50) cluster(Region) controls(mobility_externo mobility_interno population) `opt'
            ereturn list
            graph export "$OUT/did_multiplegt/did4`v'`gn'_`var'.eps", replace
        }
    }
    graph drop _all
}
exit

