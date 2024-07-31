* Program to generate forestplots -- used by ipdmetan etc. but can also be run by itself
* April 2013
*   Forked from main ipdmetan code
* September 2013
*   Following UK Stata Users Meeting, reworked the plotid() option as recommended by Vince Wiggins

* version 1.0  David Fisher  31jan2014

* version 1.01  David Fisher  07feb2014
* Reason: fixed bug - random-effects note being overlaid on x-axis labels

* version 1.02  David Fisher  20feb2014
* Reason: allow user to affect null line options

*! version 1.03  David Fisher  23jul2014
* Reason: implented a couple of suggestions from Phil Jones
* Weighting is now consistent across plotid groups
* Tidying up some code that unnecessarily restricted where user-defined lcols/rcols could be plotted
* Minor bug fixes and code simplification
* New (improved?) textsize and aspect ratio algorithm

* Known issues, July 2014: xlabels when effect sizes are far from `h0', see line 323

* Coding of _USE:
* _USE == 0  subgroup labels (headings)
* _USE == 1  successfully estimated trial-level effects
* _USE == 2  unsuccessfully estimated trial-level effects ("Insufficient data")
* _USE == 3  subgroup effects
* _USE == 4  between-subgroup heterogeneity info
* _USE == 5  overall effect
* _USE == 6  blank lines
* _USE == 9  titles


program define forestplot, sortpreserve /*rclass*/

version 10		// metan is v9 and this doesn't use any more recent commands/syntaxes; v10 used only for sake of help file extension

syntax [namelist(min=3 max=5)] [if] [in] [, ///
	/// /* Sub-plot identifier for applying different appearance options, and dataset identifier to separate plots */
	PLOTID(string) DATAID(varname) ///
	/// /* General -forestplot- options (including any "passed through" from another program, e.g. ipdmetan) */
	BY(name) Classic DP(integer 2) EFORM EFFect(string) INTERaction IPDOVER LABels(name) LCOLs(namelist) NULLOFF RCOLs(namelist) ///
	noNAmes noNUll noOVerall noPRESERVE noSTATs noSUbgroup noWT ///
	/// /* x-axis options */
	XLABel(string) XTICk(string) XTItle(string asis) Range(numlist min=2 max=2) FAVours(string asis) FP(real 999) ///
	/// /* other "fine-tuning" options */
	noADJust ASText(integer 50) ///
	* ]

	* If forestplot is being run "stand-alone" (i.e. not called from ipdmetan etc.), parse eform option
	if "`preserve'" == "" {
		_get_eformopts, soptions eformopts(`options') allowed(__all__)
		local options `"`s(options)'"'
		local eform = cond(`"`s(eform)'"'=="", "", "eform")
		if `"`effect'"'==`""' local effect = cond(`"`s(str)'"'=="", "Effect", `"`s(str)'"')
		if `"`interaction'"'!=`""' local effect `"Interact. `effect'"'
	}
	marksample touse	// do this immediately, so that -syntax- can be used again
	
	if "`nulloff'"!="" local null "nonull"	// allow "nulloff" as alternative to "nonull" for compatability with -metan-
	local graphopts `"`options'"'			// "graph region" options (also includes plotopts for now)

	* Set up variable names
	if `"`namelist'"'==`""' {		// if not specified, assume "standard" varnames			
		local _ES "_ES"
		local _LCI "_LCI"
		local _UCI "_UCI"
		if `"`ipdover'"'==`""' local _WT "_WT"
		else local _NN "_NN"
		local _USE "_USE"
	}
	else {							// else check syntax of user-specified varnames
		local 0 `namelist'
		syntax varlist(min=3 max=5 numeric)
		tokenize `varlist'
		local _ES `1'
		local _LCI `2'
		local _UCI `3'
		if `"`ipdover'"'==`""' local _WT = cond(`"`4'"'!=`""', `"`4'"', "_WT")
		else local _NN = cond(`"`4'"'!=`""', `"`4'"', "_NN")
		local _USE = cond(`"`5'"'!=`""', `"`5'"', "_USE")
	}

	quietly {	// ends on line 778
	
		*** Set up data to use
		capture confirm numeric var `_USE'
		if _rc {
			if _rc!=7 {			// `_USE' does not exist
				tempvar _USE
				gen `_USE' = cond(missing(`_ES'*`_LCI'*`_UCI'), 2, 1)
			}
			else {
				nois disp as err `"_USE: variable `_USE' exists but is not numeric"'
				exit 198
			}
		}
		markout `touse' `_USE'
		if `"`overall'"'!=`""' {
			replace `touse'=0 if inlist(`_USE', 4, 5)
		}
		replace `touse'=0 if `"`subgroup'"'!=`""' & `_USE' == 3
		count if `touse'
		if !r(N) {
			nois disp as err "no observations"
			exit 2000
		}
		tempvar touse2 allobs obs
		gen long `allobs'=_n
		bysort `touse' (`allobs') : gen long `obs' = _n if `touse'
		drop `allobs'
		
		* Check existence of `_ES', `_LCI', `_UCI' (all required)
		foreach x in _ES _LCI _UCI {
			confirm numeric var ``x''
		}
		
		* Check existence of `_WT' (or `_NN' if ipdover)
		if "`ipdover'"!="" {
			capture confirm numeric var `_NN'
			if _rc {
				disp as err "variable `_NN' not found or not numeric"
				disp as err "(required by ipdover)"
				exit 198
			}
			local _WT "_NN"		// from now on, `_WT' contains whatever data is to be used as awweights (whether _WT, _NN or other)
		}
		else {
			capture confirm numeric var `_WT'
			if _rc {
				if _rc!=7 {
					tempvar _WT
					gen `_WT' = 1 if `touse'	// generate as constant if doesn't exist
				}
				else {
					nois disp as err `"_WT: variable `_WT' exists but is not numeric"'
					exit 198
				}
			}
		}
		local awweight "[aw= `_WT']"
		summ `_WT' if `touse' & inlist(`_USE', 1, 2), meanonly	// July 2014: global min & max weights, to maintain consistency across subgroups
		local minwt = r(min)
		local maxwt = r(max)
		
		* Check validity of `_USE' (already sorted out existence)
		capture assert !missing(`_ES'*`_LCI'*`_UCI') if `touse' & `_USE'==1
		local rctemp = _rc
		capture assert missing(`_ES'*`_LCI'*`_UCI') if `touse' & `_USE'==2
		if `rctemp' | _rc {
			nois disp as err `"effect sizes do not match with value of _USE"'
			exit 198
		}

		* Generate ordering variable (reverse sequential, since y axis runs bottom to top)
		assert inrange(`_USE', 0, 6)
		tempvar id
		bysort `touse' (`obs') : gen int `id' = _N - _n + 1 if `touse'
		
		* Check existence of `labels' and `by'
		foreach x in labels by {
			local X = upper("`x'")
			if `"``x''"'!=`""' confirm var ``x''
			else cap confirm var _`X'
			if !_rc local `x' "_`X'"			// use default varnames if they exist and option not explicitly given
		}
		
		* Sort out `plotid'
		if `"`plotid'"'==`""' {
			tempvar plotid
			gen byte `plotid' = 1 if `touse'	// create plotid as constant if not specified
			local np = 1
		}
		else {
			if "`preserve'" != "" disp _n _c		// spacing if following on from ipdmetan (etc.)

			capture confirm var _BY
			local _by = cond(_rc, "", "_BY")
			capture confirm var _OVER
			local _over = cond(_rc, "", "_OVER")
			
			local 0 `plotid'
			syntax name(name=plname id="plotid") [, List noGRaph]

			if "`plname'"!="_n" {
				confirm var `plname'
				tab `plname' if `touse', m
				if r(r)>20 {
					nois disp as err "plotid: variable takes too many values"
					exit 198
				}
				if `"`_over'"'==`""' {
					count if `touse' & inlist(`_USE', 1, 2) & missing(`plname')
					if r(N) {
						nois disp as err "Warning: plotid contains missing values"
						nois disp as err "plotid groups and/or allocated ordinal numbers may not be as expected"
						if "`list'"=="" nois disp as err "This may be checked using the 'list' suboption to 'plotid'"
					}
				}
			}
			
			* Create ordinal version of plotid...
			gen `touse2' = `touse' * inlist(`_USE', 1, 2, 3, 5)
			local plvar `plname'

			* ...extra tweaking if passed through from ipdmetan/ipdover (i.e. _STUDY, and possibly _OVER, exists)
			if inlist("`plname'", "_STUDY", "_n", "_LEVEL", "_OVER") {
				capture confirm var _STUDY
				local _study = cond(_rc, "_LEVEL", "_STUDY")
				tempvar smiss
				gen `smiss' = missing(`_study')
				
				if inlist("`plname'", "_STUDY", "_n") {
					tempvar plvar
					bysort `touse2' `smiss' (`_over' `_study') : gen `plvar' = _n if `touse2' & !`smiss'
				}
				else if "`plname'"=="_LEVEL" {
					tempvar plvar
					bysort `touse2' `smiss' `_by' (`_over' `_study') : gen `plvar' = _n if `touse2' & !`smiss'
				}
			}
			tempvar plobs plotid
			bysort `touse2' `smiss' `plvar' (`obs') : gen long `plobs' = `obs'[1] if `touse2'
			bysort `touse2' `smiss' `plobs' : gen long `plotid' = (_n==1) if `touse2'
			replace `plotid' = sum(`plotid')
			local np = `plotid'[_N]					// number of `plotid' levels
			label var `plotid' "plotid"
			sort `obs'
			
			* Optionally list observations contained within each plotid group
			if "`list'" != "" {
				nois disp as text _n "plotid: observations marked by " as res "`plname'" as text ":"
				forvalues p=1/`np' {
					nois disp as text _n "-> plotid = " as res `p' as text ":"
					nois list `dataid' `_USE' `_by' `_over' `labels' if `touse2' & `plotid'==`p', table noobs sep(0)
				}
				if `"`graph'"'!=`""' exit
			}
			drop `touse2' `plobs' `smiss'
		}
		

		// SORT OUT LCOLS AND RCOLS
		// Default "lcol1" (if using ipdmetan) is list of study names, headed "Study ID"
		// If "_LABELS" exists, check whether labels exist for use==1 | use==2
		// N.B. macro `names' (noNAmes) is "optionally off"; macro `snames' contains study names!
		tempvar snames
		if `"`labels'"'==`""' local names `"nonames"'		// turn on option noNAmes
		else {
			capture assert missing(`labels') if `touse' & inlist(`_USE', 1, 2)
			if !_rc local names `"nonames"'				// turn on option noNAmes
		}
		if `"`names'"'==`""' clonevar `snames' = `labels' if `touse' & inlist(`_USE', 1, 2)
		else gen `snames' = ""
		local lcols = trim(`"`snames' `lcols'"')	// `snames' is always the first element of `lcols', even if blank.

	
		// EFFECT SIZE AND WEIGHT COLUMNS
		// by default, rcols 1 and 2 are effect sizes and weights
		// `stats' and `wt' turn these optionally off.
		// (N.B. if "wtn" specified, _WT is replaced with _NN, and `nowt' turns of _NN instead)
		if "`stats'" == "" {
			tempvar estText
			if `"`eform'"'!=`""' local xexp "exp"
			gen str `estText' = string(`xexp'(`_ES'), "%10.`dp'f") ///
				+ " (" + string(`xexp'(`_LCI'), "%10.`dp'f") + ", " + string(`xexp'(`_UCI'), "%10.`dp'f") + ")" ///
				if `touse' & inlist(`_USE', 1, 3, 5) & !missing(`_ES')
			replace `estText' = "(Insufficient data)" if `touse' & `_USE' == 2
			if `"`eform'"'==`""' {
				replace `estText' = " " + `estText' if `touse' & `_ES'>=0 & !inlist(`_USE', 4, 6)	// indent by one character if non-negative, to line up
			}
			
			* Variable label
			if "`effect'" == "" {
				if "`interaction'"!="" local effect "Interaction effect"
				else local effect `"Effect"'
			}
			* else if "`effect'"=="Haz. Ratio" local effect "Haz.`=char(160)'Ratio"	// insert non-breaking space
			* else if "`effect'"=="Odds Ratio" local effect "Odds`=char(160)'Ratio"	// insert non-breaking space
			label var `estText' `"`effect' (`c(level)'% CI)"'
		}
		if "`wt'" == "" {
			if "`ipdover'" == "" {
				tempvar weightpc
				gen `weightpc' = 100*`_WT' if `touse' & inlist(`_USE', 1, 3, 5) & !missing(`_ES')
				format `weightpc' %4.2f
				label var `weightpc' "% Weight"
			}
			else local lcols = trim(`"`lcols' `_NN'"')			// by default, place _NN to left of plot if ipdover...
		}
		local rcols = trim(`"`estText' `weightpc' `rcols'"')	// ...but effect-size text and I-V weights to right of plot.
		
		// Test validity of lcols and rcols
		foreach x in lcols {					// "lcols" has to exist
			cap confirm var ``x'' 
			if _rc {
				nois disp as err "variable `x' not defined in option lcols"
				exit _rc
			}
		}
		if `"`rcols'"' != `""' {
			foreach x in rcols {
				cap confirm var ``x'' 
				if _rc {
					nois disp as err "variable `x' not defined in option rcols"
					exit _rc
				}
			}
		}
		local lcolsN : word count `lcols'
		local rcolsN : word count `rcols'

		
		// GET MIN AND MAX DISPLAY
		// SORT OUT TICKS- CODE PINCHED FROM MIKE AND FIDDLED. TURNS OUT I'VE BEEN USING SIMILAR NAMES...
		// AS SUGGESTED BY JS JUST ACCEPT ANYTHING AS TICKS AND RESPONSIBILITY IS TO USER!
		summ `_LCI' if `touse', meanonly
		local DXmin = r(min)			// minimum confidence limit
		summ `_UCI' if `touse', meanonly
		local DXmax = r(max)			// maximum confidence limit
		// DXmin & DXmax ARE THE LEFT AND RIGHT COORDS OF THE GRAPH PART

		local h0 = 0
		
		* July 2014: further work needed for cases where DXmin and DXmax are both far from zero
		* i.e. no need to have `h0' included, or at least not whichever of `ii' and -`ii' is more extreme relative to the data
		* Also has relevance to "xtitle" (since xmlabel is placed at `h0')
		* if `h0'<`DXmin' & `DXmin'-`h0' > `DXmax'-`DXmin' local h0
		* if `h0'>`DXmax' & `h0'-`DXmax' > `DXmax'-`DXmin' local h0
		* local h0 = cond(`"`null'"'==`""', 0, (`DXmax'-`DXmin')/2)
		
		* xlabel not supplied by user: choose sensible values
		* default is for symmetrical limits, with 3 labelled values including null
		* if `nonull', find symmetrical limits around midpoint between DXmin and DXmax
		if "`xlabel'" == "" {
			local Gmodxhi=max(abs(float(`DXmin')), abs(float(`DXmax')))
			if `Gmodxhi'==. {
				local Gmodxhi=2
			}
			local DXmin=-`Gmodxhi'
			local DXmax=`Gmodxhi'
			
			* DF added March 2013: choose "sensible" label values for x-axis
			if `"`eform'"'==`""' {		// linear scale
				local mag = ceil(abs(log10(abs(float(`Gmodxhi')))))*sign(log10(abs(float(`Gmodxhi'))))	// order of magnitude
				local xdiff = abs(float(`Gmodxhi')-`mag')
				local xlab = `"`h0'"'
				foreach i of numlist 1 2 5 10 {
					local ii = (`i'^`mag')
					if abs(float(`Gmodxhi') - `ii') <= float(`xdiff') {
						local xdiff = abs(float(`Gmodxhi') - `ii')
						local xlab = `"`xlab' `ii' -`ii'"'
					}
				}
			}
			else {						// log scale
				local mag = round(`Gmodxhi'/ln(2))
				local xdiff = abs(float(`Gmodxhi') - float(ln(2)))
				local xlab `"`h0'"'
				forvalues i=1/`mag' {
					local ii = ln(2^`i')
					local xlab = `"`xlab' `ii' -`ii'"'	// display all powers of 2
				}
				
				* If effect is small, use 1.5, 1.33, 1.25 or 1.11 instead, as appropriate
				foreach i of numlist 1.5 `=1/0.75' 1.25 `=1/0.9' {
					local ii = ln(`i')
					if abs(float(`Gmodxhi') - `ii') <= float(`xdiff') {
						local xdiff = abs(float(`Gmodxhi') - `ii')
						local xlab = `"`xlab' `ii' -`ii'"'
					}
				}					
			}
			numlist `"`xlab'"'
			local xlablist=r(numlist)
		}
		
		* xlabel supplied by user: parse and apply
		else {
			local 0 `"`xlabel'"'
			syntax anything(name=xlablist) [, FORCE *]
			local xlabopts `"`options'"'

			if `"`eform'"'!=`""' {					// assume given on exponentiated scale if "eform" specified, so need to take logs
				numlist "`xlablist'", range(>0)		// in which case, all values must be greater than zero
				local n : word count `r(numlist)'
				forvalues i=1/`n' {
					local xi : word `i' of `r(numlist)'
					local xlablist2 `"`xlablist2' `=ln(`xi')'"'
				}
				local xlablist "`xlablist2'"
			}
			if "`force'" == "" {
				numlist "`xlablist' `DXmin' `DXmax'", sort
				local n : word count `r(numlist)' 
				local DXmin2 : word 1 of `r(numlist)'
				local DXmax2 : word `n' of `r(numlist)'
			
				local Gmodxhi=max(abs(`DXmin'), abs(`DXmax'), abs(`DXmin2'), abs(`DXmax2'))	
				if `Gmodxhi'==.  local Gmodxhi=2
				local DXmin=-`Gmodxhi'
				local DXmax=`Gmodxhi'
			}										// "force" option only changes things if user supplies xlabel
			else {
				numlist "`xlablist'", sort
				local n : word count `r(numlist)' 
				local DXmin : word 1 of `r(numlist)'
				local DXmax : word `n' of `r(numlist)'
			}
		}
		
		* Ticks
		if "`xtick'" == "" {
			local xticklist `xlablist'		// if not specified, default to same as labels
		}
		else {
			gettoken xticklist : xtick, parse(",")
			if `"`eform'"'!=`""' {					// assume given on exponentiated scale if "eform" specified, so need to take logs
				numlist "`xticklist'", range(>0)		// in which case, all values must be greater than zero
				local n : word count `r(numlist)'
				forvalues i=1/`n' {
					local xi : word `i' of `r(numlist)'
					local xticklist2 `"`xticklist2' `=ln(`xi')'"'
				}
				local xticklist "`xticklist2'"
			}
			else {
				numlist "`xticklist'"
				local xticklist=r(numlist)
			}
		}
		
		* Range
		if "`range'" != `""' {
			if `"`eform'"'!=`""' {
				numlist "`range'", range(>0)
				tokenize "`range'"
				local range `"`=ln(`1')' `=ln(`2')'"'
			}
			else {
				numlist "`range'"
				local range=r(numlist)
			}
		}

		* Final calculation of DXmin and DXmax
		if "`range'" == `""' {
			numlist "`xlablist' `xticklist' `DXmin' `DXmax'", sort
			local n : word count `r(numlist)' 
			local DXmin : word 1 of `r(numlist)'
			local DXmax : word `n' of `r(numlist)'
		}
		else {
			numlist "`range'", sort
			local n : word count `r(numlist)' 
			local DXmin : word 1 of `r(numlist)'
			local DXmax : word `n' of `r(numlist)'
		}
			
		* If on exponentiated scale, re-label x-axis with exponentiated values (nothing else should need changing)
		if "`eform'" != "" {
			local xlblcmd
			foreach i of numlist `xlablist' {
				local lbl = string(`=exp(`i')',"%7.3g")
				local xlblcmd `"`xlblcmd' `i' "`lbl'""'
			}
		}
		else local xlblcmd `"`xlablist'"'
			
		local DXwidth = `DXmax'-`DXmin'
		* if `DXmin' > 0 local h0 = 1				// July 2014: don't think this is needed

		// END OF TICKS AND LABELS

		
		* Need to make changes to pre-existing data now, plus adding new obs to the dataset now
		* so use -preserve- before continuing
		* (if not already preserved by ipdmetan etc.) 
		
		if "`preserve'" == "" preserve
		
		// MAKE OFF-SCALE ARROWS -- fairly straightforward
		tempvar offscaleL offscaleR offLeftX offLeftX2 offRightX offRightX2 offYlo offYhi
		gen `touse2' = `touse' * (`_USE' == 1)
		gen `offscaleL' = `touse2' * (`_LCI' < `DXmin')
		gen `offscaleR' = `touse2' * (`_UCI' > `DXmax')

		replace `_LCI' = `DXmin' if `touse2' & `_LCI' < `DXmin'
		replace `_UCI' = `DXmax' if `touse2' & `_UCI' > `DXmax'
		replace `_LCI' = . if `touse2' & `_UCI' < `DXmin'
		replace `_UCI' = . if `touse2' & `_LCI' > `DXmax'
		replace `_ES' = . if `touse2' & `_ES' < `DXmin'
		replace `_ES' = . if `touse2' & `_ES' > `DXmax'
		drop `touse2'

		
		// DIAMONDS TAKE FOREVER...I DON'T THINK THIS IS WHAT MIKE DID
		tempvar DiamLeftX DiamRightX DiamBottomX DiamTopX DiamLeftY1 DiamRightY1 DiamLeftY2 DiamRightY2 DiamBottomY DiamTopY
		gen `touse2' = `touse' * inlist(`_USE', 3, 5)

		gen `DiamLeftX' = `_LCI' if `touse2'
		replace `DiamLeftX' = `DXmin' if `touse2' & `_LCI' < `DXmin'
		replace `DiamLeftX' = . if `touse2' & `_ES' < `DXmin'

		gen `DiamLeftY1' = `id' if `touse2'
		replace `DiamLeftY1' = `id' + 0.4*( abs((`DXmin'-`_LCI')/(`_ES'-`_LCI')) ) if `touse2' & `_LCI' < `DXmin'
		replace `DiamLeftY1' = . if `touse2' & `_ES' < `DXmin'
		
		gen `DiamLeftY2' = `id' if `touse2'
		replace `DiamLeftY2' = `id' - 0.4*( abs((`DXmin'-`_LCI')/(`_ES'-`_LCI')) ) if `touse2' & `_LCI' < `DXmin'
		replace `DiamLeftY2' = . if `touse2' & `_ES' < `DXmin'

		gen `DiamRightX' = `_UCI' if `touse2'
		replace `DiamRightX' = `DXmax' if `touse2' & `_UCI' > `DXmax'
		replace `DiamRightX' = . if `touse2' & `_ES' > `DXmax'
		
		gen `DiamRightY1' = `id' if `touse2'
		replace `DiamRightY1' = `id' + 0.4*( abs((`_UCI'-`DXmax')/(`_UCI'-`_ES')) ) if `touse2' & `_UCI' > `DXmax'
		replace `DiamRightY1' = . if `touse2' & `_ES' > `DXmax'
		
		gen `DiamRightY2' = `id' if `touse2'
		replace `DiamRightY2' = `id' - 0.4*( abs((`_UCI'-`DXmax')/(`_UCI'-`_ES')) ) if `touse2' & `_UCI' > `DXmax'
		replace `DiamRightY2' = . if `touse2' & `_ES' > `DXmax'
		
		gen `DiamBottomY' = `id' - 0.4 if `touse2'
		replace `DiamBottomY' = `id' - 0.4*( abs((`_UCI'-`DXmin')/(`_UCI'-`_ES')) ) if `touse2' & `_ES' < `DXmin'
		replace `DiamBottomY' = `id' - 0.4*( abs((`DXmax'-`_LCI')/(`_ES'-`_LCI')) ) if `touse2' & `_ES' > `DXmax'
		
		gen `DiamTopY' = `id' + 0.4 if `touse2'
		replace `DiamTopY' = `id' + 0.4*( abs((`_UCI'-`DXmin')/(`_UCI'-`_ES')) ) if `touse2' & `_ES' < `DXmin'
		replace `DiamTopY' = `id' + 0.4*( abs((`DXmax'-`_LCI')/(`_ES'-`_LCI')) ) if `touse2' & `_ES' > `DXmax'

		gen `DiamTopX' = `_ES' if `touse2'
		replace `DiamTopX' = `DXmin' if `touse2' & `_ES' < `DXmin'
		replace `DiamTopX' = `DXmax' if `touse2' & `_ES' > `DXmax'
		replace `DiamTopX' = . if `touse2' & (`_UCI' < `DXmin' | `_LCI' > `DXmax')
		gen `DiamBottomX' = `DiamTopX'
		
		drop `touse2'
		
		* Create dummy obs with global min & max weights, to maintain correct weighting
		if `"`plotid'"'!=`""' {
			local oldN = _N
			set obs `=`oldN' + 2*`np''
			forvalues i=1/`np' {
				replace `plotid' = `i' in `=`oldN' + (2*`i'-1)' / `=`oldN' + (2*`i')'
				replace `_WT' = `minwt' in `=`oldN' + (2*`i'-1)'
				replace `_WT' = `maxwt' in `=`oldN' + (2*`i')'
			}
			replace `_USE' = 1 in `=`oldN' + 1' / `=`oldN' + (2*`np')'
			replace `touse' = 1 in `=`oldN' + 1' / `=`oldN' + (2*`np')'
		}

		* Modify `touse' to take into account dummy obs
		tempvar toused
		gen `toused' = `touse'							// "touse + dummy obs", for scatter plots only
		replace `touse' = 0 if `touse' & missing(`id')	// general-purpose `touse'
		

		*** Left & right columns 
		// OPTIONS FOR L-R JUSTIFY?
		// HAVE ONE MORE COL POSITION THAN NECESSARY, COULD THEN R-JUSTIFY
		// BY ADDING 1 TO LOOP, ALSO HAVE MAX DIST FOR OUTER EDGE
		// HAVE USER SPECIFY % OF GRAPH USED FOR TEXT?

		// TITLES
		summ `id'
		local max = r(max)
		local oldN = _N
		set obs `=`oldN'+4'							// create four new observations
		forvalues i = 1/4 {							//   to store up to four lines for titles
			local Nnew`i' = `=`oldN' + `i''
			replace `id' = `max' + `i' + 1 in `Nnew`i''		// "+1" leaves a one-line gap between titles & main data
		}
		local borderline = `max' + 1 - 0.25
		replace `touse' = 1 in `=`oldN' + 1' / `=`oldN' + 4'
		replace `toused' = 1 in `=`oldN' + 1' / `=`oldN' + 4'	// mark these as "dummy obs" too, so they can be removed later
		
		tempvar strlen
		
		// LEFT COLUMNS
		forvalues i=1/`lcolsN' {
			tempvar left`i'			// for later use, to store x-axis positions of columns
			
			local lcoli : word `i' of `lcols'
			capture confirm string var `lcoli'
			if !_rc local leftLB`i' : copy local lcoli
			else {
				tempvar leftLB`i'
				capture decode `lcoli', gen(`leftLB`i'')
				if _rc {
					local f: format `lcoli'
					gen str `leftLB`i'' = string(`lcoli', "`f'")
					replace `leftLB`i'' = "" if `leftLB`i'' == "."
					
					gen long `strlen' = length(`leftLB`i'')	// July 2014: add leading spaces to right-justify
					summ `strlen', meanonly
					local maxlen=r(max)
					forvalues j=1/`=_N' {
						local diff = 2*(`maxlen' - `strlen'[`j'])	// x2 because numbers have width equal to two blank spaces
						if `diff' & `strlen'[`j'] {					// use non-breaking space (twoway doesn't seem to honour normal leading spaces)
							local blank : di _dup(`diff') `"`=char(160)'"'
							replace `leftLB`i'' = `"`blank'"' + `leftLB`i'' in `j'
						}
					}
					drop `strlen'
				}
			}
				
			// WORK OUT IF TITLE IS BIGGER THAN THE VARIABLE
			// SPREAD OVER UP TO FOUR LINES IF NECESSARY
			tempvar tmplen
			gen `tmplen' = length(`leftLB`i'')
			
			if `lcolsN'>1 local notuse0 `" & `_USE'!=0"'	// if more than one lcol, don't count width of subgroup cols
			summ `tmplen' if `touse' `notuse0', meanonly	// so they can span multiple lcols if necessary
			local otherlen = r(max)
			drop `tmplen'
				
			local colName: variable label `lcoli'
			if `"`colName'"' == "" & `"`lcoli'"' !=`"`snames'"' local colName = `"`lcoli'"'
			local target = max(`otherlen', int(length(`"`colName'"')/3))
			
			SpreadTitle `"`colName'"', target(`target') maxline(4)
			local nlines = r(nlines)
			forvalues j = `nlines'(-1)1 {
				if `"`r(title`j')'"'!=`""' {
					local k=`nlines'-`j'+1
					replace `leftLB`i'' = `"`r(title`j')'"' in `Nnew`k''
					replace `_USE' = 9 in `Nnew`k''
				}
			}
		}		// end of forvalues i=1/`lcolsN'

		
		* Now copy across previously generated titles (overall, sub est etc.)
		if `"`labels'"'!=`""' replace `leftLB1' = `labels' if `touse' & inlist(`_USE', 0, 3, 4, 5)

		// RIGHT COLUMNS
		forvalues i=1/`rcolsN' {		// if `rcols'==0, loop will be skipped
			tempvar right`i'			// for later use, to store x-axis positions of columns
			
			local rcoli : word `i' of `rcols'
			cap confirm string var `rcoli'
			if !_rc local rightLB`i' : copy local rcoli
			else {
				tempvar rightLB`i'
				capture decode `rcoli', gen(`rightLB`i'')
				if _rc {					// July 2014: pure numeric; add spaces to LHS to right-justify
					local f: format `rcoli'
					gen str `rightLB`i'' = string(`rcoli', "`f'")
					replace `rightLB`i'' = "" if `rightLB`i'' == "."
					
					gen long `strlen' = length(`rightLB`i'')
					summ `strlen', meanonly
					local maxlen=r(max)
					forvalues j=1/`=_N' {
						local diff = 2*(`maxlen' - `strlen'[`j'])	// x2 because numbers have width equal to two blank spaces
						if `diff' & `strlen'[`j'] {					// use non-breaking space (twoway doesn't seem to honour normal leading spaces)
							local blank : di _dup(`diff') `"`=char(160)'"'
							replace `rightLB`i'' = `"`blank'"' + `rightLB`i'' in `j'
						}
					}
					drop `strlen'
				}
			}
				
			// WORK OUT IF TITLE IS BIGGER THAN THE VARIABLE
			// SPREAD OVER UP TO FOUR LINES IF NECESSARY
			tempvar tmplen
			gen `tmplen' = length(`rightLB`i'')
			summ `tmplen' if `touse', meanonly
			local otherlen = r(max)
			drop `tmplen'

			local colName: variable label `rcoli'
			if `"`colName'"' == "" & `"`rcoli'"' !=`"`snames'"' local colName = `"`rcoli'"'
			local target = max(`otherlen', int(length(`"`colName'"')/3))
			
			SpreadTitle `"`colName'"', target(`target') maxline(4)
			local nlines = r(nlines)
			forvalues j = `nlines'(-1)1 {
				if `"`r(title`j')'"'!=`""' {
					local k=`nlines'-`j'+1
					replace `rightLB`i'' = `"`r(title`j')'"' in `Nnew`k''
					replace `_USE' = 9 in `Nnew`k''
				}
			}
		}		// end of forvalues i=1/`rcols'
		
		// now get rid of extra title rows if they weren't used
		drop if `toused' & missing(`_USE')
		
		// Calculate "leftWDtot" and "rightWDtot" -- the total widths to left and right of graph area
		// Unless multiple lcols/rcols, don't use titles or overall stats, just trial stats.
		local leftWDtot = 0
		local rightWDtot = 0
		local leftWDtotNoTi = 0

		tempvar tempWD
		forvalues i = 1/`rcolsN' {
			getWidth `rightLB`i'' `tempWD'
			summ `tempWD' if `touse', meanonly					// July 2014: use ALL rows (including titles)
			local rightWD`i' = cond(r(N), r(max), 0)
			local rightWDtot = `rightWDtot' + `rightWD`i''
			drop `tempWD'
		}
		local rightWDtot = `rightWDtot' + (31/24)			// add buffer

		* July 2014: if lcols>1, check for data when _USE==3,5 -- if found, include these obs in width calculations
		local include							// clear macro
		forvalues i = 2/`lcolsN' {				// if `lcolsN'==1, loop will be skipped
			getWidth `leftLB`i'' `tempWD'
			summ `tempWD' if `touse' & inlist(`_USE', 3, 5), meanonly
			drop `tempWD'
			if r(max)>1 {
				local include `" | inlist(`_USE', 3, 5)"'
				continue, break
			}
		}		
		forvalues i = 1/`lcolsN' {
			getWidth `leftLB`i'' `tempWD'
			summ `tempWD' if `touse' & (inlist(`_USE', 1, 2, 9) `include'), meanonly	// July 2014: include titles
			local leftWD`i' = cond(r(N), r(max), 0)
			local leftWDtot = `leftWDtot' + `leftWD`i''
			drop `tempWD'
		}
		
		* Compare width of text in subgroup/overall and header rows to current `leftWDtot'
		tempvar maxLeft
		getWidth `leftLB1' `maxLeft'
		count if `touse' & !inlist(`_USE', 1, 2, 9)
		if r(N) {
			summ `maxLeft' if `touse' & !inlist(`_USE', 1, 2, 9), meanonly
			local max = r(max)
		
			if `max' > `leftWDtot' {
				if "`adjust'" != "" local leftWDtot = `max'
				
				else {
					// WORK OUT HOW FAR INTO PLOT CAN EXTEND
					// Allow _USE=0,3,4,5 to extend into main plot by (lcimin-DXmin)/DXwidth
					// i.e. 1 + ((`lcimin'-`DXmin')/`DXwidth') * ((100-`astext')/`astext')) is the percentage increase
					// to apply to (`maxLeft'+`rightWDtot')/(`newleftWDtot'+`rightWDtot').
					// Then rearrange to find `newleftWDtot'.
					tempvar lci2
					gen `lci2' = cond(`_LCI'>0, 0, `_LCI')
					summ `lci2' if `touse' & !inlist(`_USE', 1, 2, 9), meanonly
					local lcimin = r(min)
					
					// BUT don't make it any less than before, unless there are no obs with inlist(`_USE', 1, 2)
					local newleftWDtot = ((`max'+`rightWDtot') / ( ( ((`lcimin'-`DXmin')/`DXwidth') * ((100-`astext')/`astext') ) + 1)) - `rightWDtot'
					count if `touse' & inlist(`_USE', 1, 2)
					if r(N) local leftWDtot = max(`leftWDtot', `newleftWDtot')
					else local leftWDtot = `newleftWDtot'

					drop `lci2'
				}
			}		// end if `max' > `leftWDtot'
		}		// end count if `touse' & !inlist(`_USE', 1, 2, 9) if r(N)
		
		// Generate position of lcols, using user-specified `astext'
		// (% of graph width taken by text)
		local textWD = (`DXwidth'/(1-`astext'/100) - `DXwidth') / (`leftWDtot'+`rightWDtot')

		// Now, carry on as before
		// N.B. although these are constants, they need to be stored variables for use with -twoway-
		local leftWDtot2 = `leftWDtot'
		forvalues i = 1/`lcolsN' {
			gen `left`i'' = `DXmin' - `leftWDtot2'*`textWD'
			local leftWDtot2 = `leftWDtot2' - `leftWD`i''
		}
		if `rcolsN' gen `right1' = `DXmax'		// July 2014
		forvalues i = 2/`rcolsN' {				// (P.S. if `rcolsN'=0 then loop will be skipped)
			gen `right`i'' = `right`=`i'-1'' + `rightWD`=`i'-1''*`textWD'
		}
		
		// AXmin AXmax ARE THE OVERALL LEFT AND RIGHT COORDS
		local AXmin = `left1'
		local AXmax = `DXmax' + `rightWDtot'*`textWD'
		
	}	// END QUIETLY


	// FIND OPTIMAL TEXT SIZE AND ASPECT RATIOS (given user input)
	// Notes:  (David Fisher, July 2014)
	
	// Let X, Y be dimensions of graphregion; x, y be dimensions of plotregion.
	// `approxChars' is the approximate width of the plot, in "character units" (i.e. width of text divided by `astext')
	// `height' is the approximate height of the plot, in terms of rows of text
	// If Y/X = `graphAspect'<1, `textSize' is the height of a row of text relative to Y; otherwise it is height relative to X.
	// We then let `approxChars' = x, and manipulate to find the optimum text size for the plot layout.
	// Finally, maximum text size is 100/y.

	// If y/x < Y/X < 1 then X = kx (i.e. plot takes up full width of "wide" graph, with an extra margin if overall title specified)
	//   then `textSize' = `textscale'/Y = `textscale'/(X * `graphAspect') -- but X = kx * `approxChars'
	//   ==> `textSize' = `textscale'/(k * `approxChars' * `graphAspect')
	
	// If Y/X < 1 and y/x > Y/X (i.e. plot is less wide than graph) then Y=ky where k = (`height'+delta)/`height'
	//   then `textSize' = `textscale'/ky = `textscale'/(x * k * `plotAspect') = `textscale'/(`approxChars' * k * `plotAspect')
	
	// If y/x > Y/X > 1 then Y = ky
	//   then `textSize' = `textscale'/X = `textscale' * `graphAspect'/ky = (`textscale' * `graphAspect')/(`approxChars' * k * `plotAspect')
	
	// If Y/X > 1 and y/x < Y/X (i.e. plot is less tall than graph) then assume X = kx
	//   then `textSize' = `textscale'/X = `textscale'/(k * `approxChars')

	//  - Note that this code has been changed considerably from the original -metan- code.

	local 0 `", `graphopts'"'
	syntax [, ASPECT(real 0) BOXscale(real 100.0) CAPTION(string asis) NOTE(string asis) RENOTE(string) SUBtitle(string asis) ///
		TEXTscale(real 100.0) TItle(string asis) XSIZe(real 5.5) YSIZe(real 4) * ]
	numlist "`xsize' `ysize'", range(>0)
	numlist "`aspect' `boxscale' `textscale'", range(>=0)

	local height=0
	qui count if `touse' & `_USE'==9
	if r(N) local height = r(N) + 1		// add 1 to height if titles (to take account of gap)
	qui count if `touse' & `_USE'!=9
	local height = `height' + r(N)
	local condtitle = 2*(`"`title'"'!=`""') + (`"`subtitle'"'!=`""') + (`"`caption'"'!=`""') + .5*(`"`note'"'!=`""')
	local yk = (`height' + (`"`xlblcmd'"'!=`""') + (`"`favours'"'!=`""') + `condtitle')/`height'
	local xk = (`height' + `condtitle')/`height'
	
	local approxChars = (`leftWDtot' + `rightWDtot')/(`astext'/100)
	if `aspect'==0 {
		local plotAspect = `height' / `approxChars'				// "natural aspect" in terms of text
		if `plotAspect'<1 local plotAspect = 2*`plotAspect'		// if "wide", use double spacing
		else local plotAspect = 1.5*`plotAspect'				// else just use 1.5-spacing
	}	
	else local plotAspect : copy local aspect					// user-specified aspect of plotregion	
	local graphAspect = `ysize'/`xsize'							// aspect of graphregion (assume 4/5.5 unless specified)
	
	if `graphAspect' < 1 & `plotAspect' < `graphAspect' {
		local textSize = `textscale'/(`xk' * `approxChars' * `graphAspect')
		local maxtextSize = 100 * `plotAspect' / (`xk' * `graphAspect' * `height')
	}
	if `graphAspect' < 1 & `plotAspect' > `graphAspect' {
		local textSize = `textscale'/(`yk' * `approxChars' * `plotAspect')
		local maxtextSize = 100 / (`yk' * `height')
	}
	if `graphAspect' > 1 & `plotAspect' > `graphAspect' {
		local textSize = (`textscale' * `graphAspect')/(`yk' * `approxChars' * `plotAspect')
		local maxtextSize = 100 * `graphAspect' / (`yk' * `height')
	}
	if `graphAspect' > 1 & `plotAspect' < `graphAspect' {
		local textSize = `textscale' / (`xk' * `approxChars')
		local maxtextSize = 100 * `plotAspect' / (`xk' * `height')
	}

	* Unless user has specified `textScale', cap the textsize so that lines of text to do not overlap
	if `textscale'==100 local textSize = min(`textSize', `maxtextSize')
	
	* Notes: for random-effects analyses, sample-size weights, or user-defined (will overwrite the first two)
	if "`wtn'"!="" & `"`renote'"'==`""' local renote "NOTE: Point estimates are weighted by sample size"
	if `"`renote'"'!=`""' & `"`note'"'==`""' {
		local note `""`renote'", size(`=`textSize'*.75')"'	// use 75% of text size used for rest of plot
	}
	
	local graphopts `"`options' aspect(`plotAspect') caption(`caption') note(`note') subtitle(`subtitle') title(`title') xsize(`xsize') ysize(`ysize')"'
	
	/*
	// Graphing quantities
	return scalar aspect = `plotAspect'
	return scalar approx_chars = `approxChars'
	return scalar height = `height'
	return scalar ysize = `ysize'
	return scalar xsize = `xsize'
	return scalar textsize = `textSize'
	*/

	// PLOT COLUMNS OF TEXT (lcols/rcols)
	forvalues i = 1/`lcolsN' {
		local lcolCommands `"`macval(lcolCommands)' || scatter `id' `left`i'' if `toused', msymbol(none) mlabel(`leftLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`textSize')"'
	}
	forvalues i = 1/`rcolsN' {
		local rcolCommands `"`macval(rcolCommands)' || scatter `id' `right`i'' if `toused', msymbol(none) mlabel(`rightLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`textSize')"'
	}	// July 2014: "& `_USE' != 4" removed from scatter `id' `right`i'' -- test for consequences


	// FAVOURS
	if `"`favours'"' != `""' {

		* DF added Jan 2013: allow multiple lines (cf twoway title option)
		gettoken leftfav rest : favours, parse("#") quotes
		if `"`leftfav'"'!=`"#"' {
			while `"`rest'"'!=`""' {
				gettoken next rest : rest, parse("#") quotes
				if `"`next'"'==`"#"' continue, break
				local leftfav `"`leftfav' `next'"'
			}
		}
		else local leftfav `""'
		local rightfav `"`rest'"'

		if `fp'>0 & `fp'<999 {					// 999 is a dummy "default"
			local leftfp = -`fp'
			local rightfp = `fp'
		}
		else if "`h0'" != "" & "`null'" == "" {
			local leftfp = `DXmin' + (`h0'-`DXmin')/2
			local rightfp = `h0' + (`DXmax'-`h0')/2
		}
		else {
			local leftfp = `DXmin'
			local rightfp = `DXmax'
		}
		local favopt `"xmlabel(`leftfp' `"`leftfav'"' `rightfp' `"`rightfav'"', noticks labels labsize(`textSize') labgap(5))"'
	}

	* xtitle - uses 'xmlabel' options, not 'title' options!  Parse all 'title' options to give suitable error message
	else if `"`xtitle'"' != `""' {
		local 0 `"`xtitle'"'
		syntax [anything] [, TSTYle(string) ORIENTation(string) SIze(string) Color(string) Justification(string) ///
			ALignment(string) Margin(string) LINEGAP(string) WIDTH(string) HEIGHT(string) BOX NOBOX ///
			BColor(string) FColor(string) LStyle(string) LPattern(string) LWidth(string) LColor(string) ///
			BMargin(string) BEXpand(string) PLACEment(string) *]
		if `"`size'"'!=`""' local labsizeopt `"labsize(`size')"'
		if `"`labsize'"'!=`""' local labsizeopt `"labsize(`labsize')"'
		else local labsizeopt `"labsize(`textSize')"'
		if `"`color'"'!=`""' local labcoloropt `"labcolor(`color')"'
		if `"`labcolor'"'!=`""' local labcoloropt `"labcolor(`labcolor')"'
		if !(`"`tstyle'"'==`""' & `"`orientation'"'==`""' & `"`justification'"'==`""' & `"`alignment'"'==`""' ///
			& `"`margin'"'==`""' & `"`linegap'"'==`""' & `"`width'"'==`""' & `"`height'"'==`""' & `"`box'"'==`""' & `"`nobox'"'==`""' ///
			& `"`bcolor'"'==`""' & `"`fcolor'"'==`""' & `"`lstyle'"'==`""' & `"`lpattern'"'==`""' & `"`lwidth'"'==`""' & `"`lcolor'"'==`""' ///
			& `"`bmargin'"'==`""' & `"`bexpand'"'==`""' & `"`placement'"'==`""') {
			disp as err `"option xtitle uses xmlabel options, not title options!  see help axis_label_options"'
			exit 198
		}
		local xtitleopt `"xmlabel(`h0' `"`anything'"', noticks labels `labsizeopt' `labcoloropt' labgap(5) `options')"'
	}


	// GRAPH APPEARANCE OPTIONS
	local boxSize = `boxscale'/150

	summ `id', meanonly
	local DYmin = r(min)-1
	local DYmax = r(max)+1

	tempvar useno
	qui gen byte `useno' = `_USE' * inlist(`_USE', 3, 5) if `touse'

	cap confirm var `dataid'
	if _rc {
		tempvar dataid
		qui gen byte `dataid'=1 if `touse'
	}
	sort `touse' `dataid' `id'
	qui replace `useno' = `useno'[_n-1] if `useno'<=`useno'[_n-1] & `dataid'==`dataid'[_n-1]	// find the largest value (from 3 & 5) "so far"

	* Flag obs through which the line should be drawn
	tempvar ovLine
	qui gen `ovLine'=.		// need this var to exist regardless

	summ `useno' if `touse', meanonly
	if r(max) {
		tempvar olinegroup check ovMin ovMax
		qui gen int `olinegroup' = (`_USE'==`useno') * (`useno'>0)
		qui bysort `touse' `dataid' (`id') : replace `olinegroup' = sum(`olinegroup') if inlist(`_USE', 1, 2, 3, 5)	// study obs & pooled results

		* Store min and max values for later plotting
		qui gen byte `check' = inlist(`_USE', 1, 2)
		qui bysort `touse' `dataid' `olinegroup' (`check') : replace `check' = `check'[_N]	// only draw oline if there are study obs in the same olinegroup
		qui replace `ovLine' = `_ES' if `touse' & `check' & `_USE'==`useno' & `useno'>0 & !(`_ES' > `DXmax' | `_ES' < `DXmin')

		sort `touse' `dataid' `olinegroup' `id'
		qui by `touse' `dataid' `olinegroup' : gen `ovMin' = `id'[1]-0.5 if `touse' & `_USE'==`useno' & `useno'>0 & !missing(`ovLine')
		qui by `touse' `dataid' `olinegroup' : gen `ovMax' = `id'[_N]+0.5 if `touse' & `_USE'==`useno' & `useno'>0 & !missing(`ovLine')
		drop `useno' `olinegroup' `check' `dataid'
	}


	*** Get options and store plot commands

	** "Global" options (includes null line)
	local 0 `", `graphopts'"'
	syntax [, ///
		/// /* standard options */
		BOXOPts(string asis) DIAMOPts(string asis) POINTOPts(string asis) CIOPts(string asis) OLINEOPts(string asis) ///
		/// /* non-diamond options */
		PPOINTOPts(string asis) PCIOPts(string asis) NLINEOPts(string asis) * ]

	local rest `"`options'"'

	* Global CI style (bare lines or capped lines)
	* (Also test for disallowed options during same parse)
	local 0 `", `ciopts'"'
	syntax [, HORizontal VERTical Connect(string asis) RCAP * ]
	if `"`horizontal'"'!=`""' | `"`vertical'"'!=`""' {
		disp as err "ciopts: options horizontal/vertical not allowed"
		exit 198
	}			
	if `"`connect'"' != `""' {
		disp as err "ciopts: option connect() not allowed"
		exit 198
	}
	local ciopts `"`options'"'
	local CIPlotType = cond("`rcap'"=="", "rspike", "rcap")
	local pCIPlotType `CIPlotType'

	* "Default" options
	local dispShape = cond("`interaction'"!="", "circle", "square")
	local DefColor = cond("`classic'"!="", "black", "180 180 180")
	local DefBoxopts = `"mcolor("`DefColor'") msymbol(`dispShape') msize(`boxSize')"'
	local DefCIopts `"lcolor(black) mcolor(black)"'		// includes "mcolor" for arrows (doesn't affect rspike/rcap)
	local DefPointopts `"msymbol(diamond) mcolor(black) msize(vsmall)"'
	local DefOlineopts `"lwidth(thin) lcolor(maroon) lpattern(shortdash)"'
	local DefDiamopts `"lcolor("0 0 100")"'
	local DefPPointopts `"msymbol("`dispShape'") mlcolor("0 0 100") mfcolor("none")"'
	local DefPCIopts `"lcolor("0 0 100")"'

	* Loop over possible values of `plotid' and test for plot#opts relating specifically to each value
	numlist "1/`np'"
	local plvals=r(numlist)			// need both of these as numlists,
	local pplvals `plvals'			//    for later macro manipulations
	forvalues p = 1/`np' {

		local 0 `", `rest'"'
		syntax [, ///
			/// /* standard options */
			BOX`p'opts(string asis) DIAM`p'opts(string asis) POINT`p'opts(string asis) CI`p'opts(string asis) OLINE`p'opts(string asis) ///
			/// /* non-diamond options
			PPOINT`p'opts(string asis) PCI`p'opts(string asis) * ]

		local rest `"`options'"'

		* Check if any options were found specifically for this value of `p'
		if `"`box`p'opts'`diam`p'opts'`point`p'opts'`ci`p'opts'`oline`p'opts'`ppoint`p'opts'`pci`p'opts'"' != `""' {
			
			local pplvals : list pplvals - p			// remove from list of "default" plotids
			
			* WEIGHTED SCATTER PLOT
			local 0 `", `box`p'opts'"'
			syntax [, MLABEL(string asis) MSIZe(string asis) * ]			// check for disallowed options
			if `"`mlabel'"' != `""' {
				disp as error "box`p'opts: option mlabel() not allowed"
				exit 198
			}
			if `"`msize'"' != `""' {
				disp as error "box`p'opts: option msize() not allowed"
				exit 198
			}
			qui count if `touse' & `_USE'==1 & `plotid'==`p'
			if r(N) {
				summ `_WT' if `touse' & `_USE'==1 & `plotid'==`p', meanonly
				if !r(N) disp as err `"No weights found for plotid `p'"'
				else local scPlot `"`macval(scPlot)' || scatter `id' `_ES' `awweight' if `toused' & `_USE'==1 & `plotid'==`p', `DefBoxopts' `boxopts' `box`p'opts'"'
			}		// N.B. scatter if `toused' <-- "dummy obs" for consistent weighting
			
			* CONFIDENCE INTERVAL PLOT
			local 0 `", `ci`p'opts'"'
			syntax [, HORizontal VERTical Connect(string asis) RCAP LColor(string asis) MColor(string asis) * ]	// check for disallowed options + rcap
			if `"`horizontal'"'!=`""' | `"`vertical'"'!=`""' {
				disp as error "ci`p'opts: options horizontal/vertical not allowed"
				exit 198
			}			
			if `"`connect'"' != `""' {
				di as error "ci`p'opts: option connect() not allowed"
				exit 198
			}
			if `"`lcolor'"'!=`""' & `"`mcolor'"'==`""' {
				local mcolor `lcolor'						// for pc(b)arrow
			}
			local ci`p'opts `"mcolor(`mcolor') lcolor(`lcolor') `options'"'
			local CIPlot`p'Type `CIPlotType'										// global status
			local CIPlot`p'Type = cond("`rcap'"=="", "`CIPlot`p'Type'", "rcap")		// overwrite global status if appropriate
			local CIPlot `"`macval(CIPlot)' || `CIPlot`p'Type' `_LCI' `_UCI' `id' if `touse' & `_USE'==1 & `plotid'==`p' & !`offscaleL' & !`offscaleR', hor `DefCIopts' `ciopts' `ci`p'opts'"'		

			qui count if `plotid'==`p' & `offscaleL' & `offscaleR'
			if r(N) {													// both ends off scale
				local CIPlot `"`macval(CIPlot)' || pcbarrow `id' `_LCI' `id' `_UCI' if `touse' & `plotid'==`p' & `offscaleL' & `offscaleR', `DefCIopts' `ciopts' `ci`p'opts'"'
			}
			qui count if `plotid'==`p' & `offscaleL' & !`offscaleR'
			if r(N) {													// only left off scale
				local CIPlot `"`macval(CIPlot)' || pcarrow `id' `_UCI' `id' `_LCI' if `touse' & `plotid'==`p' & `offscaleL' & !`offscaleR', `DefCIopts' `ciopts' `ci`p'opts'"'
				if "`CIPlot`p'Type'" == "rcap" {			// add cap to other end if appropriate
					local CIPlot `"`macval(CIPlot)' || rcap `_UCI' `_UCI' `id' if `touse' & `plotid'==`p' & `offscaleL' & !`offscaleR', hor `DefCIopts' `ciopts' `ci`p'opts'"'
				}
			}
			qui count if `plotid'==`p' & !`offscaleL' & `offscaleR'
			if r(N) {													// only right off scale
				local CIPlot `"`macval(CIPlot)' || pcarrow `id' `_LCI' `id' `_UCI' if `touse' & `plotid'==`p' & !`offscaleL' & `offscaleR', `DefCIopts' `ciopts' `ci`p'opts'"'
				if "`CIPlot`p'Type'" == "rcap" {			// add cap to other end if appropriate
					local CIPlot `"`macval(CIPlot)' || rcap `_LCI' `_LCI' `id' if `touse' & `plotid'==`p' & !`offscaleL' & `offscaleR', hor `DefCIopts' `ciopts' `ci`p'opts'"'
				}
			}

			* POINT PLOT (point estimates -- except if "classic")
			if "`classic'" == "" {
				local pointPlot `"`macval(pointPlot)' || scatter `id' `_ES' if `touse' & `_USE'==1 & `plotid'==`p', `DefPointopts' `pointopts' `point`p'opts'"'
			}
			
			* OVERALL LINE(S) (if appropriate)
			summ `ovLine' if `plotid'==`p', meanonly
			if r(N) {
				local olinePlot `"`macval(olinePlot)' || rspike `ovMin' `ovMax' `ovLine' if `touse' & `plotid'==`p', `DefOlineopts' `olineopts' `oline`p'opts'"'
			}		

			* POOLED EFFECT - DIAMOND
			* Assume diamond if no "pooled point/CI" options, and no "interaction" option
			if `"`ppointopts'`ppoint`p'opts'`pciopts'`pci`p'opts'`interaction'"' == `""' {
				local diamPlot `"`macval(diamPlot)' || pcspike `DiamLeftY1' `DiamLeftX' `DiamTopY' `DiamTopX' if `touse' & `plotid'==`p', `DefDiamopts' `diamopts' `diam`p'opts'"'
				local diamPlot `"`macval(diamPlot)' || pcspike `DiamTopY' `DiamTopX' `DiamRightY1' `DiamRightX' if `touse' & `plotid'==`p', `DefDiamopts' `diamopts' `diam`p'opts'"'
				local diamPlot `"`macval(diamPlot)' || pcspike `DiamRightY2' `DiamRightX' `DiamBottomY' `DiamBottomX' if `touse' & `plotid'==`p', `DefDiamopts' `diamopts' `diam`p'opts'"'
				local diamPlot `"`macval(diamPlot)' || pcspike `DiamBottomY' `DiamBottomX' `DiamLeftY2' `DiamLeftX' if `touse' & `plotid'==`p', `DefDiamopts' `diamopts' `diam`p'opts'"'
			}
			
			* POOLED EFFECT - PPOINT/PCI
			else {
				if `"`diam`p'opts'"'!=`""' {
					disp as err `"plotid `p': cannot specify options for both diamond and pooled point/CI"'
					disp as err `"diamond options will be ignored"'
				}	
			
				local 0 `", `pci`p'opts'"'
				syntax [, HORizontal VERTical Connect(string asis) RCAP *]
				if `"`horizontal'"'!=`""' | `"`vertical'"'!=`""' {
					disp as error "pci`p'opts: options horizontal/vertical not allowed"
					exit 198
				}			
				if `"`connect'"' != `""' {
					di as error "pci`p'opts: option connect() not allowed"
					exit 198
				}
				local pCIPlot`p'Type `pCIPlotType'											// global status
				local pCIPlot`p'Type = cond("`rcap'"=="", "`pCIPlot`p'Type'", "rcap")		// overwrite global status if appropriate
				local pCIPlot `"`macval(pCIPlot)' || `pCIPlotType' `_LCI' `_UCI' `id' if `touse' & inlist(`_USE', 3, 5) & `plotid'==`p', hor `DefPCIopts' `pciopts' `pci`p'opts'"'
				local ppointPlot `"`macval(ppointPlot)' || scatter `id' `_ES' if `touse' & inlist(`_USE', 3, 5) & `plotid'==`p', `DefPPointopts' `ppointopts' `ppoint`p'opts'"'
			}
		}
	}

	* Find invalid/repeated options
	* any such options would generate a suitable error message at the plotting stage
	* so just exit here with error, to save the user's time
	if regexm(`"`rest'"', "(box|diam|point|ci|oline|ppoint|pci)([0-9]+)") {
		local badopt = regexs(1)
		local badp = regexs(2)
		
		disp as err `"`badopt'`badp'opts: "' _c
		if `: list badp in plvals' disp as err "option supplied multiple times; should only be supplied once"
		else disp as err `"`badp' is not a valid plotid value"'
		exit 198
	}

	local graphopts `rest'		// this is now *just* the standard "twoway" options
								// i.e. the specialist "forestplot" options have been filtered out

					
	* FORM "DEFAULT" TWOWAY PLOT COMMAND (if appropriate)
	* Changed so that FOR WEIGHTED SCATTER each pplval is plotted separately (otherwise weights are messed up)
	* Other (nonweighted) plots can continue to be plotted as before
	if `"`pplvals'"'!=`""' {

		local pplvals2 : copy local pplvals						// copy, just for use in line 1171
		local pplvals : subinstr local pplvals " " ",", all		// so that "inlist" may be used

		* WEIGHTED SCATTER PLOT
		local 0 `", `boxopts'"'
		syntax [, MLABEL(string asis) MSIZe(string asis) * ]	// check for disallowed options
		if `"`mlabel'"' != `""' {
			disp as err "boxopts: option mlabel() not allowed"
			exit 198
		}
		if `"`msize'"' != `""' {
			disp as err "boxopts: option msize() not allowed"
			exit 198
		}
		if `"`pplvals'"'==`"`plvals'"' {		// if no plot#opts specified, can plot all plotid groups at once
			qui summ `_WT' if `_USE'==1 & inlist(`plotid', `pplvals')
			if r(N) local scPlot `"`macval(scPlot)' || scatter `id' `_ES' `awweight' if `toused' & `_USE'==1 & inlist(`plotid', `pplvals'), `DefBoxopts' `boxopts'"'
		}
		else {		// else, need to plot each group separately to maintain correct weighting (July 2014)
			foreach p of local pplvals2 {
				qui summ `_WT' if `_USE'==1 & `plotid'==`p'
				if r(N) local scPlot `"`macval(scPlot)' || scatter `id' `_ES' `awweight' if `toused' & `_USE'==1 & `plotid'==`p', `DefBoxopts' `boxopts'"'
			}
		}		// N.B. scatter if `toused' <-- "dummy obs" for consistent weighting
		
		* CONFIDENCE INTERVAL PLOT
		local 0 `", `ciopts'"'
		syntax [, HORizontal VERTical Connect(string asis) RCAP LColor(string asis) MColor(string asis) * ]	// check for disallowed options + rcap
		if `"`horizontal'"'!=`""' | `"`vertical'"'!=`""' {
			disp as error "ciopts: options horizontal/vertical not allowed"
			exit 198
		}			
		if `"`connect'"' != `""' {
			di as error "ciopts: option connect() not allowed"
			exit 198
		}
		if `"`lcolor'"'!=`""' & `"`mcolor'"'==`""' {
			local mcolor `lcolor'						// for pc(b)arrow
		}
		local ciopts `"mcolor(`mcolor') lcolor(`lcolor') `options'"'
		local CIPlotType = cond("`rcap'"=="", "`CIPlotType'", "rcap")		// overwrite global status if appropriate
		local CIPlot `"`macval(CIPlot)' || `CIPlotType' `_LCI' `_UCI' `id' if `touse' & `_USE'==1 & inlist(`plotid', `pplvals') & !`offscaleL' & !`offscaleR', hor `DefCIopts' `ciopts'"'

		qui count if inlist(`plotid', `pplvals') & `offscaleL' & `offscaleR'
		if r(N) {													// both ends off scale
			local CIPlot `"`macval(CIPlot)' || pcbarrow `id' `_LCI' `id' `_UCI' if `touse' & inlist(`plotid', `pplvals') & `offscaleL' & `offscaleR', `DefCIopts' `ciopts'"'
		}
		qui count if inlist(`plotid', `pplvals') & `offscaleL' & !`offscaleR'
		if r(N) {													// only left off scale
			local CIPlot `"`macval(CIPlot)' || pcarrow `id' `_UCI' `id' `_LCI' if `touse' & inlist(`plotid', `pplvals') & `offscaleL' & !`offscaleR', `DefCIopts' `ciopts'"'
			if "`CIPlotType'" == "rcap" {			// add cap to other end if appropriate
				local CIPlot `"`macval(CIPlot)' || rcap `_UCI' `_UCI' `id' if `touse' & inlist(`plotid', `pplvals') & `offscaleL' & !`offscaleR', hor `DefCIopts' `ciopts'"'
			}
		}
		qui count if inlist(`plotid', `pplvals') & !`offscaleL' & `offscaleR'
		if r(N) {													// only right off scale
			local CIPlot `"`macval(CIPlot)' || pcarrow `id' `_LCI' `id' `_UCI' if `touse' & inlist(`plotid', `pplvals') & !`offscaleL' & `offscaleR', `DefCIopts' `ciopts'"'
			if "`CIPlotType'" == "rcap" {			// add cap to other end if appropriate
				local CIPlot `"`macval(CIPlot)' || rcap `_LCI' `_LCI' `id' if `touse' & inlist(`plotid', `pplvals') & !`offscaleL' & `offscaleR', hor `DefCIopts' `ciopts'"'
			}
		}

		* POINT PLOT
		local pointPlot `"`macval(pointPlot)' || scatter `id' `_ES' if `touse' & `_USE'==1 & inlist(`plotid', `pplvals'), `DefPointopts' `pointopts'"'

		* OVERALL LINE(S) (if appropriate)
		summ `ovLine' if inlist(`plotid', `pplvals'), meanonly
		if r(N) {
			local olinePlot `"`macval(olinePlot)' || rspike `ovMin' `ovMax' `ovLine' if `touse' & inlist(`plotid', `pplvals'), `DefOlineopts' `olineopts'"'
		}

		* POOLED EFFECT - DIAMOND
		* Assume diamond if no "pooled point/CI" options, and no "interaction" option
		if `"`ppointopts'`pciopts'`interaction'"' == `""' {
			local diamPlot `"`macval(diamPlot)' || pcspike `DiamLeftY1' `DiamLeftX' `DiamTopY' `DiamTopX' if `touse' & inlist(`plotid', `pplvals'), `DefDiamopts' `diamopts'"'
			local diamPlot `"`macval(diamPlot)' || pcspike `DiamTopY' `DiamTopX' `DiamRightY1' `DiamRightX' if `touse' & inlist(`plotid', `pplvals'), `DefDiamopts' `diamopts'"'
			local diamPlot `"`macval(diamPlot)' || pcspike `DiamRightY2' `DiamRightX' `DiamBottomY' `DiamBottomX' if `touse' & inlist(`plotid', `pplvals'), `DefDiamopts' `diamopts'"'
			local diamPlot `"`macval(diamPlot)' || pcspike `DiamBottomY' `DiamBottomX' `DiamLeftY2' `DiamLeftX' if `touse' & inlist(`plotid', `pplvals'), `DefDiamopts' `diamopts'"'
		}
		
		* POOLED EFFECT - PPOINT/PCI
		else {
			if `"`diamopts'"'!=`""' {
				disp as err _n `"plotid: cannot specify options for both diamond and pooled point/CI"'
				disp as err `"diamond options will be ignored"'
			}
			
			local 0 `", `pciopts'"'
			syntax [, HORizontal VERTical Connect(string asis) RCAP *]		// check for disallowed options + rcap
			if `"`horizontal'"'!=`""' | `"`vertical'"'!=`""' {
				disp as error "pciopts: options horizontal/vertical not allowed"
				exit 198
			}			
			if `"`connect'"' != `""' {
				di as error "pciopts: option connect() not allowed"
				exit 198
			}
			local pCIPlotType = cond("`rcap'"=="", "`pCIPlotType'", "rcap")		// overwrite global status if appropriate
			local pCIPlot `"`macval(pCIPlot)' || `pCIPlotType' `_LCI' `_UCI' `id' if `touse' & inlist(`_USE', 3, 5) & inlist(`plotid', `pplvals'), hor `DefPCIopts' `pciopts'"'
			local ppointPlot `"`macval(ppointPlot)' || scatter `id' `_ES' if `touse' & inlist(`_USE', 3, 5) & inlist(`plotid', `pplvals'), `DefPPointopts' `ppointopts'"'
		}
	}
		
	// END GRAPH OPTS

	// DF: modified to use added line approach instead of pcspike (less complex & poss. more efficient as fewer vars)
	// null line (unless switched off)
	if "`null'" == "" {
		local 0 `", `nlineopts'"'
		syntax [, HORizontal VERTical Connect(string asis) * ]
		if `"`horizontal'"'!=`""' | `"`vertical'"'!=`""' {
			disp as err "nlineopts: options horizontal/vertical not allowed"
			exit 198
		}			
		if `"`connect'"' != `""' {
			disp as err "nlineopts: option connect() not allowed"
			exit 198
		}
		local nullCommand `"|| function y=`h0', horiz range(`DYmin' `borderline') n(2) lwidth(thin) lcolor(black) `options'"'
	}



	***************************
	***     DRAW GRAPH      ***
	***************************

	#delimit ;

	twoway
	/* OVERALL AND NULL LINES FIRST */ 
		`olinePlot' `nullCommand'
	/* PLOT BOXES AND PUT ALL THE GRAPH OPTIONS IN THERE, PLUS NOTE FOR RANDOM-EFFECTS */ 
		`scPlot' `notecmd'
			yscale(range(`DYmin' `DYmax') noline) ylabel(none) ytitle("")
			xscale(range(`AXmin' `AXmax')) xlabel(`xlblcmd', labsize(`textSize'))
			yline(`borderline', lwidth(thin) lcolor(gs12))
	/* FAVOURS OR XTITLE */
			`favopt' `xtitleopt'
	/* PUT LABELS UNDER xticks? Yes as labels now extended */
			xtitle("") legend(off) xtick("`xticklist'")
	/* NEXT, CONFIDENCE INTERVALS (plus offscale if necessary) */
		`CIPlot'
	/* DIAMONDS (or markers+CIs if appropriate) FOR SUMMARY ESTIMATES */
		`diamPlot' `ppointPlot' `pCIPlot'
	/* COLUMN VARIBLES (including effect sizes and weights on RHS by default) */
		`lcolCommands' `rcolCommands'
	/* LAST OF ALL PLOT EFFECT MARKERS TO CLARIFY */
		`pointPlot'
	/* Other options */
		|| , `graphopts' /* RMH added */ plotregion(margin(zero)) ;

	#delimit cr


end





program define getWidth
version 9.0

//	ROSS HARRIS, 13TH JULY 2006
//	TEXT SIZES VARY DEPENDING ON CHARACTER
//	THIS PROGRAM GENERATES APPROXIMATE DISPLAY WIDTH OF A STRING
//	FIRST ARG IS STRING TO MEASURE, SECOND THE NEW VARIABLE

//	PREVIOUS CODE DROPPED COMPLETELY AND REPLACED WITH SUGGESTION
//	FROM Jeff Pitblado

qui {
	gen `2' = 0
	count
	local N = r(N)
	forvalues i = 1/`N'{
		local this = `1'[`i']
		local width: _length `"`this'"'
		replace `2' =  `width' +1 in `i'
	}
} // end qui

end



* exit

//	METAN UPDATE
//	ROSS HARRIS, DEC 2006
//	MAIN UPDATE IS GRAPHICS IN THE _dispgby PROGRAM
//	ADDITIONAL OPTIONS ARE lcols AND rcols
//	THESE AFFECT DISPLAY ONLY AND ALLOW USER TO SPECIFY
//	VARIABLES AS A FORM OF TABLE. THIS EXTENDS THE label(namevar yearvar)
//	SYNTAX, ALLOWING AS MANY LEFT COLUMNS AS REQUIRED (WELL, LIMIT IS 10)
//	IF rcols IS OMMITTED DEFAULT IS THE STUDY EFFECT (95% CI) AND WEIGHT
//	AS BEFORE- THESE ARE ALWAYS IN UNLESS OMITTED USING OPTIONS
//	ANYTHING ADDED TO rcols COMES AFTER THIS.


********************
** May 2007 fixes **
********************

//	"nostandard" had disappeared from help file- back in
//	I sq. in return list
//	sorted out the extra top line that appears in column labels
//	fixed when using aspect ratio using xsize and ysize so inner bit matches graph area- i.e., get rid of spaces for long/wide graphs
//	variable display format preserved for lcols and rcols
//	abbreviated varlist now allowed
//	between groups het. only available with fixed
//	warnings if any heterogeneity with fixed (for between group het if any sub group has het, overall est if any het)
// 	nulloff option to get rid of line



* Subroutine to "spread" titles out over multiple lines if appropriate
* (copied from ipdmetan)
* Updated July 2014
program SpreadTitle, rclass

	syntax anything(name=title id="title string"), [TArget(integer 0) MAXWidth(integer 0) MAXLines(integer 0)]
	* Target = aim for this width, but allow expansion if alternative is wrapping "too early" (i.e before line is adequately filled)
	* Maxwidth = absolute maximum width.
	
	if !`target' {
		if !`maxwidth' {
			disp as err "must specify at least one of target or maxwidth"
			exit 198
		}
		local target = `maxwidth'
	}
	
	if `maxwidth' {
		cap assert `maxwidth'>=`target'
		if _rc {
			disp as err "maxwidth value must be greater than or equal to target value"
			exit 198
		}
	}
	
	local titlelen = length(`title')
	local spread = int(`titlelen'/`target') + 1
	* if `maxlines' & `spread'>`maxlines' local spread = `maxlines'
	
	local line = 1
	local end = 0
	local count = 2

	local title1 = word(`title', 1)
	local newwidth = length(`"`title1'"')
	local next = word(`title', 2)
	
	while `"`next'"' != "" {
		local check = trim(`"`title`line''"' + " " +`"`next'"')			// (potential) next iteration of `title`line''
		if length(`"`check'"') > `titlelen'/`spread' {					// if too long
																		// and further from target than before, or greater than maxwidth
			if abs(length(`"`check'"')-(`titlelen'/`spread')) > abs(length(`"`title`line''"')-(`titlelen'/`spread')) ///
				| (`maxwidth' & length(`"`check'"') > `maxwidth') {
				if `maxlines' & `line'==`maxlines' {					// if reached max no. of lines
					local title`line' `"`check'"'						//   - use next iteration anyway (to be truncated)
					continue, break										//   - break loop
				}
				else {													// otherwise:
					local ++line										//  - new line
					local title`line' `"`next'"'						//  - begin new line with next word
				}
			}
			else local title`line' `"`check'"'		// else use next iteration
			
		}
		else local title`line' `"`check'"'		// else use next iteration

		local ++count
		local next = word(`title', `count')
		local newwidth = max(`newwidth', length(`"`title`line''"'))		// update `newwidth'
	}

	* If last string is too long (including in above case), truncate
	if `newwidth' <= `target' local maxwidth = `target'
	else local maxwidth = cond(`maxwidth', min(`newwidth', `maxwidth'), `newwidth')
	if length(`"`title`line''"') > `maxwidth' local title`line' = substr(`"`title`line''"', 1, `maxwidth')
	
	* Return strings
	forvalues i=1/`line' {
		return local title`i' = trim(`"`title`i''"')
	}
	return scalar nlines = `line'
	return scalar maxwidth = min(`newwidth', `maxwidth')
	
end
