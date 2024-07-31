{smcl}
{* *! version 1.0  David Fisher  30jan2014}{...}
{cmd:help petometan}{right: ({browse "http://www.stata-journal.com/article.html?article=st0384":SJ15-2: st0384})}
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{cmd:petometan} {hline 2}}Perform meta-analysis using the Peto log-rank method{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:petometan} {it:trt_var} [{cmd:,} {it:options}]

{synoptset 20}{...}
{synopthdr}
{synoptline}
{synopt :{opt stu:dy(varname)}}specify the study identifier to be displayed on screen and in the forest plot{p_end}
{synopt :{opt by(varname)}}specify an optional study-level subgroup identifier{p_end}
{synopt :{opt mat:save(name)}}store data in matrix {it:name}{p_end}
{synopt :{opt nogr:aph}}suppress forest plot{p_end}
{synopt :{opt nohet}}suppress heterogeneity statistics in forest plot{p_end}
{synopt :{opt noov:erall}}suppress overall pooling{p_end}
{synopt :{opt nosu:bgroup}}suppress pooling within subgroups{p_end}
{synopt :{cmd:ovstat(q)}}display Q statistics in forest plot instead of I-squared{p_end}
{synopt :{opt str:ata(varlist)}}specify further variables by which to stratify the log-rank calculations{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:petometan} performs aggregate or individual participant data
meta-analysis using the Peto log-rank method.  It can be seen as a
meta-analytic extension (or alternative) to {cmd:sts test} but only for
two-arm, single-failure data.  Data must be {helpb stset}.


{marker options}{...}
{title:Options}

{phang}
{opt study(varname)} specifies the meta-analysis study identifier to be
presented on screen and in the forest plot.  In terms of log-rank
calculation of the overall effect statistics, it is equivalent to
specifying {opt strata()}.

{phang}
{opt by(varname)} specifies a trial-level subgroup identifier.  It does
not affect the calculation of the overall pooled effect, but it will be
presented in the output along with relevant subgroup effects and
heterogeneity statistics.

{phang}
{opt matsave(name)} requests that study-level data be stored in
{it:name}.  The stored data consist of study or subgroup identifiers,
numbers of participants and observed events by treatment arm, expected
events, and O-E and V.

{phang}
{opt nograph}, {opt nohet}, {opt nooverall}, and {opt nosubgroup}
suppress, respectively, the forest plot, display of heterogeneity
statistics in the plot, overall pooling, and pooling within subgroups.

{phang}
{cmd:ovstat(q)} requests that Q statistics be used in the forest plot instead
of I-squared statistics.

{phang}
{opt strata(varlist)} specifies that further variables be used in log-rank
calculations but not presented in the output.


{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:petometan} stores the following in {cmd:r()}:

{synoptset 13 tabbed}{...}
{p2col 5 25 29 2: Scalars}{p_end}
{synopt:{cmd:r(OE)}}total difference (O-E) of observed and expected numbers of events{p_end}
{synopt:{cmd:r(V)}}total variance of the observed number of events (hypergeometric variance){p_end}
{synopt:{cmd:r(N)}}total number of participants{p_end}
{synopt:{cmd:r(o)}}total observed number of events{p_end}
{synopt:{cmd:r(chi2)}}overall chi-squared statistic{p_end}
{synopt:{cmd:r(lnHR)}}overall hazard ratio{p_end}
{synopt:{cmd:r(selnHR)}}standard error of overall hazard ratio{p_end}
{synopt:{cmd:r(Q)}}overall Q statistic for heterogeneity{p_end}
{synopt:{cmd:r(k)}}number of studies{p_end}


{title:Acknowledgment}

{pstd}
The code is based heavily on the built-in Stata command {helpb sts test}.


{title:Author}

{pstd}David Fisher{p_end}
{pstd}MRC Clinical Trials Unit at University College London{p_end}
{pstd}London, UK{p_end}
{pstd}{browse "mailto:d.fisher@ucl.ac.uk":d.fisher@ucl.ac.uk}{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 15, number 2: {browse "http://www.stata-journal.com/article.html?article=st0384":st0384}{p_end}
