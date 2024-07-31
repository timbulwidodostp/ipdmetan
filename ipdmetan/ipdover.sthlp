{smcl}
{* *! version 1.02  David Fisher  23jul2014}{...}
{cmd:help ipdover}{right: ({browse "http://www.stata-journal.com/article.html?article=st0384":SJ15-2: st0384})}
{hline}

{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{cmd:ipdover} {hline 2}}Generate data for forest plots outside the context of meta-analysis{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:ipdover}
	[{it:{help exp_list}}]{cmd:,}
	{cmd:over(}{it:varlist} [{cmd:, {ul:m}issing}]{cmd:)} 
	[{cmd:over(}{it:varname} [{cmd:, {ul:m}issing}]{cmd:)}]
	[{cmd:plotid(_BY}|{cmd:_OVER}|{cmd:_LEVEL}|{cmd:_n} 
        [{cmd:,} {cmdab:l:ist} {cmdab:nogr:aph}]{cmd:)} {it:ipdmetan_options}
	{opt forest:plot(forestplot_options)}]{cmd::} {it:command}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ipdover} extends the functionality of {helpb ipdmetan} outside the
context of meta-analysis.  It does not perform any pooling or heterogeneity
calculations; instead, it creates forest plots of subgroup
analyses within one trial dataset.  Basic syntax is

{phang2}
{cmd:. ipdover}{cmd:,} {opt over(varlist)}{cmd::} {it:command}

{pstd}
which fits the model {it:command} once within each level of each variable in
{it:varlist} and saves effect sizes and standard errors for screen output and
display of a forest plot.  Any e-class regression command (whether built-in or
user-defined) should be compatible with this basic syntax;
see {cmd:ipdmetan} for help with using non e-class commands.

{pstd}
{cmd:ipdover} functions similarly to {cmd:ipdmetan}, with the following
main differences in syntax:  the {cmd:over()} options replace
{cmd:study()} and {cmd:by()}; the {cmd:ad()} and {cmd:re()} options are not
permitted; and the {cmd:plotid()} option has a slightly different syntax (see
below).

{pstd}
Note that forest plots produced by {cmd:ipdover} are weighted by sample size rather than by the inverse of the variance.

{pstd}
Saved datasets (see {cmd:ipdmetan}) include the following identifier variables:{p_end}
{p2colset 8 24 24 8}
{p2col:{cmd:_BY}}subset of data (c.f. {helpb by}){p_end}
{p2col:{cmd:_OVER}}{cmd:over()} variable{p_end}
{p2col:{cmd:_LEVEL}}level of {cmd:over()} variable{p_end}


{marker options}{...}
{title:Options}

{phang}
{cmd:over(}{it:varlist} [{cmd:, missing}]{cmd:)}
[{cmd:over(}{it:varname} [{cmd:, missing}]{cmd:)}] specifies the variables
within whose levels {it:command} is to be fit.  The option can be repeated
at most once, in which case the second option must contain one
{it:varname} defining subsets of the data (c.f. {helpb by}).  {cmd:over()} is
required.

{pmore}
All variables must be either integer valued or string.  Variable and value
labels will appear in output where appropriate.

{phang2}
{opt missing} requests that missing values be treated as potential subgroups
or subsets (the default is to exclude them).

{phang}
{cmd:plotid(_BY}|{cmd:_OVER}|{cmd:_LEVEL}|{cmd:_n} 
[{cmd:, list nograph}]{cmd:)} functions similarly as in 
{helpb ipdmetan}; but instead of a {it:varname}, it accepts one of the
following values, which correspond to variables created in saved datasets created
by {cmd:ipdover}:{p_end}
{p2colset 9 24 24 8}
{p2col:{cmd:_BY}}group observations by levels of {cmd:_BY}{p_end}
{p2col:{cmd:_OVER}}group observations by levels of {cmd:_OVER}{p_end}
{p2col:{cmd:_LEVEL}}group observations by levels of {cmd:_LEVEL}{p_end}
{p2col:{cmd:_n}}allow each observation to be its own group{p_end}

{phang}
{it:ipdmetan_options} (given the caveat in the {help ipdover##description:Description});
see {helpb ipdmetan##options:ipdmetan}.

{phang}
{opt forestplot(forestplot_options)}; see
{helpb forestplot##options:forestplot}.


{title:Examples}

{pstd}
Example 1: Using the Hosmer and Lemeshow low birthweight data from {helpb logistic}{p_end}
{phang2}
{cmd:. webuse lbw}{p_end}
{phang2}
{cmd:. ipdover, over(race smoke ht) or forestplot(favours("Odds of LBW decrease" "as age increases" # "Odds of LBW increase" "as age increases") fp(0.5)): logistic low age}{p_end}

{pstd}
Example 2: Treatment effect by covariate subgroup by trial, using the example individual participant data meta-analysis dataset from {helpb ipdmetan}{p_end}
{phang2}
{cmd:. use ipdmetan_example}{p_end}
{phang2}
{cmd:. stset tcens, fail(fail)}{p_end}
{phang2}
{cmd:. ipdover, over(stage) over(trialid) hr nosubgroup nooverall forestplot(favours(Favours treatment # Favours control)): stcox trt}


{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:ipdover} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{p2col 5 25 29 2: Scalars}{p_end}
{synopt:{cmd:r(k)}}number of included trials k{p_end}
{synopt:{cmd:r(n)}}number of included patients{p_end}
{synopt:{cmd:r(mu_hat)}}overall effect size{p_end}
{synopt:{cmd:r(se_mu_hat)}}standard error of overall effect size{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 25 29 2: Macros}{p_end}
{synopt:{cmd:r(estvar)}}name of stored coefficient{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 25 29 2: Matrices}{p_end}
{synopt:{cmd:r(coeffs)}}matrix of study and subgroup identifiers, effect
coefficients, and numbers of participants{p_end}


{title:Acknowledgment}

{pstd}
I thank Phil Jones at the University of Western Ontario, Canada, for
suggesting improvements in functionality.


{title:Author} 

{pstd}David Fisher{p_end}
{pstd}MRC Clinical Trials Unit at University College London{p_end}
{pstd}London, UK{p_end}
{pstd}{browse "mailto:d.fisher@ucl.ac.uk":d.fisher@ucl.ac.uk}{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 15, number 2: {browse "http://www.stata-journal.com/article.html?article=st0384":st0384}{p_end}
