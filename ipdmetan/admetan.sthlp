{smcl}
{* *! version 1.01  David Fisher  15apr2014}{...}
{cmd:help admetan}{right: ({browse "http://www.stata-journal.com/article.html?article=st0384":SJ15-2: st0384})}
{hline}

{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{cmd:admetan} {hline 2}}Perform two-stage inverse-variance meta-analysis of aggregate data{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:admetan} {varlist} {ifin} [{cmd:,} {it:options}]

{pstd}
where {it:varlist} is 
{p_end}

{phang2}
{it:ES} {it:seES}{p_end}

{pstd}
or
{p_end}

{phang2}
{it:ES} {it:lci} {it:uci}

{synoptset 30}{...}
{synopthdr}
{synoptline}
{synopt :{opt npts(varname)}}specify variable containing participant numbers{p_end}
{synopt :{opt study(varname)}}specify study (trial) identifier{p_end}
{p2col  :{help ipdmetan##options:{it:ipdmetan_options}}}any {cmd:ipdmetan} options except those listed below{p_end}
{synopt	:{cmdab:forest:plot(}{help forestplot##options:{it:forestplot_options}}{cmd:)}}forestplot options{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:admetan} represents the special case of {helpb ipdmetan} where no
individual participant data are available, only aggregate data.  As such,
it may be considered a direct alternative to {helpb metan}.  {it:varlist}
must be supplied and must contain either two or three of the following variables:  the
effect size (on the normal scale), followed by either the standard error
of the effect size or the lower and upper 95% confidence limits.  Note
that {cmd:admetan} assumes that confidence intervals are symmetric and
calculates variances as confidence interval width or 2z.  Hence, supplied
confidence limits must be based on a Normal distribution, or the pooled
result will not be accurate.

{pstd}
Any {helpb ipdmetan##options:ipdmetan} options may be used 
except {opt ad()}, {opt nototal}, and {cmd:plotid(_BYAD)}.


{marker options}{...}
{title:Options}

{phang}
{opt npts(varname)} specifies a variable containing numbers of
participants in each study to display in tables and forest plots.

{phang}
{opt study(varname)} specifies the variable containing the study
identifier, which must be either integer valued or string.  The default
is for the studies to be labeled sequentially as "1", "2", etc.

{phang}
{it:ipdmetan_options} (given the caveat in the 
{help admetan##description:Description}); see 
{helpb ipdmetan##options:ipdmetan}.

{phang}
{opt forestplot(forestplot_options)}; see
{helpb forestplot##options:forestplot}.


{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:admetan} stores the same results in {cmd:r()} as {helpb ipdmetan}.


{title:Author}

{pstd}David Fisher{p_end}
{pstd}MRC Clinical Trials Unit at University College London{p_end}
{pstd}London, UK{p_end}
{pstd}{browse "mailto:d.fisher@ucl.ac.uk":d.fisher@ucl.ac.uk}{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 15, number 2: {browse "http://www.stata-journal.com/article.html?article=st0384":st0384}{p_end}
