* "over" functionality for ipdmetan
* In separate program on advice of Patrick Royston
* created by David Fisher, February 2013

*! version 1.0  David Fisher  31jan2014

prog define ipdover, rclass

	version 11
	* NOTE: mata requires v9 (??)
	* factor variable syntax requires 11.0
	* but neither is vital to ipdover

	// <ipdmetan_stuff> : <command>
	_on_colon_parse `0'
	local command `"`s(after)'"'
	local 0 `"`s(before)'"'

	* For ipdover, MUST specify: (1) an estimation command and (2) at least one over() option
	if `"`command'"'==`""' {
		disp as error `"Must specify an estimation command with ipdover"'
		exit 198
	}
	
	syntax [anything(name=exp_list equalok)], OVER(string) [STUDY(string) BY(string) AD(string) * ]
	if `"`study'"' != `""' {
		disp as err "cannot specify study() with ipdover; please use over() or the ipdmetan program"
		exit 198
	}
	if `"`by'"' != `""' {
		disp as err "cannot specify by() with ipdover; please use over() or the ipdmetan program"
		exit 198
	}
	if `"`ad'"' != `""' {
		disp as err "cannot use aggregate data with ipdover"
		exit 198
	}
	local Options `"`options'"'		// note capital "O" to distinguish
	local 0 `"`over'"'
	syntax varlist [, Missing]
	local over1 `"`varlist'"'
	local missing1 `"`missing'"'

	* Parse "over" options and map to study() and by() as appropriate
	* (N.B. "over" as an option is NOT passed to ipdmetan!)
	local 0 `", `Options'"'
	syntax [, OVER(string) * ]
	local Options `"`options'"'		// note capital "O" to distinguish
	local 0 `"`over'"'
	syntax [varlist(default=none)] [, Missing]
	local over2 `"`varlist'"'
	local missing2 `"`missing'"'	

	local 0 `", `Options'"'
	syntax [, OVER(string) * ]
	if `"`over'"'!=`""' {
		di as error "may not specify more than 2 over() options"
		exit 198
	}
	
	local nv2 : word count `over2'		// no. of vars in `over2'
	if `nv2'==0 {						// only `over1' supplied
		local before `"study(`over1')"'
		local smissing = cond(`"`missing1'"'!=`""', "smissing", "")
	}
	else if `nv2'==1 {		// `over2' supplied correctly (i.e. single var)
		local before `"study(`over1') by(`over2')"'
		local smissing = cond(`"`missing1'"'!=`""', "smissing", "")
		local bymissing = cond(`"`missing2'"'!=`""', "bymissing", "")
	}
	else {
		disp as err "cannot specify multiple vars to second 'over' option"
		exit 198
	}
	
	* Now run ipdmetan with "ipdover" option
	* This is a marker that data should not be pooled (i.e. is not a MA)
	ipdmetan `exp_list', `before' `options' ipdover `smissing' `bymissing' : `command'
	return add

end
	
	
	
