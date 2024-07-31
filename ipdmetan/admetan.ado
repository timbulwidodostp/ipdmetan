* "Aggregate-data only" functionality, in separate program for simplicity
* originally written by David Fisher, June 2013

*! version 1.0  David Fisher  31jan2014

prog define admetan, rclass

	version 10
	syntax varlist(min=2 max=3 numeric) [if] [in] [, STUDY(varname) NPTS(varname) *]
	
	marksample touse
	if `"`study'"'==`""' {
		tempvar obs study
		gen int `obs' = _n
		qui bysort `touse' (`obs'): gen int `study' = _n if `touse'
		label var `study' "Study"
	}
	
	* Now run ipdmetan with "_admetan" after colon as a dummy command
	ipdmetan, study(`study') ad(, vars(`varlist') npts(`npts')) adonly `options' : admetan if `touse'
	return add

end


