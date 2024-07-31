{smcl}
{* *! version 1.03  David Fisher  23jul2014}{...}
{cmd:help ipdmetan}{right: ({browse "http://www.stata-journal.com/article.html?article=st0384":SJ15-2: st0384})}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col:{cmd:ipdmetan} {hline 2}}Perform two-stage inverse-variance individual participant data meta-analysis{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:ipdmetan}
	[{it:{help exp_list}}]{cmd:, {ul:s}tudy(}{it:study_id}[{cmd:, {ul:m}issing}]{cmd:)} [{it:options}]{cmd::} {it:estimation_command}

{synoptset 40 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Main}
{p2coldent:* {cmd:{ul:s}tudy(}{it:study_id}[{cmd:, {ul:m}issing}]{cmd:)}}specify the study identifier variable{p_end}
{synopt :{cmd:by(}{it:subgroup_id}[{cmd:, {ul:m}issing}]{cmd:)}}group the studies in the output{p_end}
{synopt :{it:{help eform_option}}}exponentiate effect sizes and confidence limits{p_end}
{synopt :{opt eff:ect(string)}}title for "effect-size" column in the output{p_end}
{synopt :{opt inter:action}}automatically identify and pool a treatment-covariate interaction{p_end}
{synopt :{opt keepall}}display all studies in the output, even those for which no effect could be estimated{p_end}
{synopt :{opt me:ssages}}print messages relating to success of model fits{p_end}
{synopt :{opt nogr:aph}}suppress the forest plot{p_end}
{synopt :{opt notab:le}}suppress printing the table of effect sizes to screen{p_end}
{synopt :{opt nohet}}suppress all heterogeneity statistics{p_end}
{synopt :{opt nooverall}}suppress overall pooling{p_end}
{synopt :{opt nosu:bgroup}}suppress pooling within subgroups{p_end}
{synopt :{opt notot:al}}suppress fitting of {it:estimation_command} to the entire dataset{p_end}
{synopt :{opt ovwt sgwt}}override default choice of whether to display overall weights or within-subgroup weights{p_end}
{synopt :{opt pool:var(model_coefficient)}}specify the coefficient to pool{p_end}
{synopt :{opt random}}specify the DerSimonian-Laird random-effects model{p_end}
{synopt :{opt re}}specify the DerSimonian-Laird random-effects model{p_end}
{synopt :{cmd:re(}{help ipdmetan##re_model:{it:re_model}}{cmd:)}}specify alternative random-effects models{p_end}
{synopt :{cmd:sortby(}{it:varname}|{cmd:_n)}}specify ordering of studies in table and forest plot{p_end}

{syntab :Forest plots}
{synopt :{cmdab:lcol:s(}{help ipdmetan##cols_info:{it:cols_info}}{cmd:)} {cmdab:rcol:s(}{help ipdmetan##cols_info:{it:cols_info}}{cmd:)}}display (or save) columns of additional data{p_end}
{synopt :{cmd:plotid(}{it:varname}|{cmd:_BYAD}[{cmd:, {ul:l}ist {ul:nogr}aph}]{cmd:)}}define groups of observations in which to apply specific plot rendition options{p_end}
{synopt :{cmd:ovstat(q)}}display Q statistics instead of I-squared{p_end}
{synopt :{cmdab:sa:ving(}{it:{help filename}}[{cmd:, replace} {cmdab:stack:label}]{cmd:)}}save results in the form of a "forest plot dataset" to {it:filename}{p_end}
{synopt :{cmdab:forest:plot(}{help forestplot##options:{it:forestplot_options}}{cmd:)}}specify other options to pass to {helpb forestplot}{p_end}

{syntab :Combined IPD/aggregate data analysis}
{synopt :{cmd:ad(}{it:{help filename}} {ifin}{cmd:,} {help ipdmetan##ad_options:{it:ad_options}}{cmd:)}}combine IPD with aggregate data stored in {it:filename}{p_end}
{synoptline}
{pstd}* {cmd:study(}{it:study_id}[{cmd:, missing}]{cmd:)} is required,{p_end}

{pstd}
where {it:model_coefficient} is a variable name, a level indicator, an
interaction indicator, or an interaction involving continuous variables
(for example, syntax of {helpb test});

{marker cols_info}{...}
{pstd}
and where {it:cols_info} has the following syntax, which is based on that of {helpb collapse}:

{pmore}
[{opt (stat)}] [{it:newname}=]{it:item} [{it:%fmt} {cmd:"}{it:label}{cmd:"}] [[{it:newname}=]{it:item} [{it:%fmt} {cmd:"}{it:label}{cmd:"}] ] {it:...} [ [{opt (stat)}] {it:...}]

{pmore}
where {it:stat} is as defined in {helpb collapse};
{it:newname} is an optional user-specified variable name;
{it:item} is the name of either a numeric returned quantity from {it:estimation_command} (in parentheses, see {help exp_list})
or a variable currently in memory; {it:%fmt} is an optional {help format}; and {cmd:"}{it:label}{cmd:"} is an optional variable label.

{marker re_model}{...}
{synoptset 24}{...}
{synopthdr :re_model}
{synoptline}
{synopt :{opt dl}}DerSimonian-Laird estimator (equivalent to specifying {cmd:re}
alone, with no suboption); default{p_end}
{synopt :{opt dlt} or {opt hk}}DerSimonian-Laird with Hartung-Knapp t-based variance estimator{p_end}
{synopt :{opt bdl} or {opt dlb}}bootstrap DerSimonian-Laird estimator{p_end}
{synopt :{opt ca}, {opt he}, or {opt vc}}Hedges variance-component estimator, also known as Cochran ANOVA-like estimator{p_end}
{synopt :{opt sj}}Sidik-Jonkman two-step estimator{p_end}
{synopt :{opt b0}}Rukhin B0 estimator{p_end}
{synopt :{opt bp}}Rukhin BP estimator{p_end}
{synopt :{opt bs}, {opt bt}, or {opt gamma}}Biggerstaff-Tweedie approximate gamma model{p_end}
{synopt :{opt eb}, {opt gq}, {opt genq}, {cmd:mp}, or {opt q}}empirical Bayes estimator, also known as generalized Q and Mandel-Paule estimator{p_end}
{synopt :{opt ml}}maximum likelihood estimator{p_end}
{synopt :{opt pl}}profile likelihood model{p_end}
{synopt :{opt reml}}restricted maximum-likelihood estimator{p_end}
{synopt :{opt sa}[{cmd:, isq(}{it:real}{cmd:)}]}sensitivity analysis with user-defined I-squared; default is {cmd:isq(0.8)}{p_end}
{synoptline}

{synoptset 15}{...}
{marker ad_options}{...}
{synopthdr :ad_options}
{synoptline}
{synopt :{opt byad}}IPD and aggregate data are to be treated as subgroups (rather than as a single set of estimates){p_end}
{synopt :{opt npts(varname)}}specify variable containing participant numbers{p_end}
{synopt :{opt vars(varlist)}}specify variables containing effect size and either standard error or 95% confidence limits on the normal scale{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ipdmetan} performs two-stage individual participant data (IPD)
meta-analysis using the inverse-variance method.  Basic syntax is

{phang2}
{cmd:. ipdmetan}{cmd:,} {opt study(study_id)}{cmd::} {it:estimation_command}

{pstd}
which fits the model {it:estimation_command} once within each level of
{it:study_id} and saves effect sizes and standard errors for pooling and displaying in a forest plot.  Any e-class regression command
(whether built-in or user-defined) should be compatible with this basic
syntax of {cmd:ipdmetan}.

{pstd}
In the case of non e-class commands -- those which do not change the
contents of {cmd:e(b)} -- the effect size and standard-error statistics to
be collected from the execution of {it:estimation_command} must be
specified manually by supplying {it:{help exp_list}}.  If
{it:estimation_command} changes the contents in {cmd:e(b)},
{it:exp_list} defaults to {cmd:_b[}{it:varname}{cmd:]}
{cmd:_se[}{it:varname}{cmd:]}, where {it:varname} is the first
independent variable within {it:estimation_command}.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{cmd:study(}{it:study_id}[{cmd:, missing}]{cmd:)} specifies the
variable that contains the study identifier, which must be either
integer valued or string.  {cmd:study()} is required.

{phang2}
{opt missing} requests that missing values be treated as potential study
identifiers (the default is to exclude them).

{phang}
{cmd:by(}{it:subgroup_id}[{cmd:, missing}]{cmd:)} specifies a variable
identifying subgroups of studies (and must therefore be constant within
studies), which must be either integer valued or string.  If an aggregate data
dataset is specified and contains a variable named {it:subgroup_id},
subgrouping will be extended to the aggregate data observations.

{phang2}
{opt missing} requests that missing values be treated as potential subgroup identifiers (the default is to exclude them).

{phang}
{it:{help eform_option}} specifies that effect sizes and confidence limits
be exponentiated in the output (table and forest plot) and generates a
heading for the effect-size column.

{pmore}
Note that {cmd:ipdmetan} does not check the validity of the particular
{it:{help eform_option}}; for example, it does not check whether
{it:estimation_command} is a survival model if {opt hr} is supplied.

{phang}
{opt effect(string)} specifies a heading for the effect-size column,
overriding any heading generated by {it:eform_option}.

{phang}
{opt interaction} specifies that {it:estimation_command} contain one or more
interaction effects expressed using factor-variable syntax (see 
{help fvvarlist}) and that the first valid interaction effect be
pooled across studies.  This is a helpful shortcut for
simple interaction analyses, but it is not foolproof or comprehensive.
The output should be checked and, if
necessary, the analysis rerun with the {cmd:poolvar()} option.

{phang}
{opt keepall} specifies that all values of {it:study_id} be
visible in the output (table and forest plot), even if no effect could
be estimated (for example, because of insufficient observations or missing
data).  For such studies, {cmd:(Insufficient data)} will appear in place of
effect estimates and weights.

{phang}
{opt messages} requests that information be printed to screen 
whether the effect size and standard-error statistics have been successfully
obtained from each study and, if applicable, whether the iterative
random-effects calculations converged successfully.

{phang}
{opt nograph} and {opt notable} suppress construction of the forest plot
and the table of effect sizes, respectively.

{phang}
{opt nohet} suppresses all heterogeneity statistics.

{phang} 
{opt nooverall} suppresses the overall pooled effect so that, for
instance, subgroups are considered entirely independently.
It also suppresses between-subgroup heterogeneity statistics (if applicable).

{phang}
{opt nosubgroup} suppresses the within-subgroup pooled effects (if
applicable) so that subgroups are displayed separately but with a
single overall pooled effect with associated heterogeneity statistics.

{phang}
{opt nototal} requests that {it:estimation_command} not be fit within
the entire dataset, for example, for time-saving reasons.  By default,
such fitting is done to check for problems in convergence and in the
validity of requested coefficients and returned expressions.  If
{cmd:nototal} is specified, either {opt poolvar()} or {it:exp_list} must
be supplied, and a message appears above the table of results warning
that estimates should be double checked by the user.

{phang}
{opt ovwt} and {opt sgwt} override the default choice of whether to
display overall weights or within-subgroup weights in the screen output
and forest plot.  Note that because weights are normalized, these options
do not affect estimation of pooled effects or heterogeneity statistics.

{phang}
{opt poolvar(model_coefficient)} specifies the
coefficient from {cmd:e(b)} that is to be pooled when
the default behavior of {cmd:ipdmetan} fails or is incorrect.
{it:model_coefficient} should be a variable name, a level indicator, an
interaction indicator, or an interaction involving continuous variables
(c.f. syntax of {helpb test}).  Equation names can be specified using
the format {cmd:poolvar(}{it:eqname}{cmd::}{it:varname}{cmd:)}.  This option is
 appropriate only when {it:estimation_command}
changes the contents of {cmd:e(b)}; otherwise, see {it:exp_list}.

{phang}
{cmd:random} or {opt re} specifies DerSimonian and Laird random effects.

{phang}
{cmd:re(}{help ipdmetan##re_model:{it:re_model}}{cmd:)} specifies other
possible random-effects models.  The default is the
DerSimonian and Laird random-effects model.  Other currently supported models
are the following:

{phang2}
{cmd:dl} specifies the DerSimonian-Laird estimator (equivalent to
specifying {cmd:re} alone, with no suboption); this is the default.

{phang2}
{cmd:dlt} or {cmd:hk} specifies the DerSimonian-Laird with Hartung-Knapp
t-based variance estimator.

{phang2}
{cmd:bdl} or {cmd:dlb} specifies the bootstrapped DerSimonian-Laird
estimator.

{phang2}
{cmd:ca}, {cmd:he}, or {cmd:vc} specifies the Hedges variance-component
estimator, also known as the Cochran estimator.

{phang2}
{cmd:sj} specifies the Sidik-Jonkman two-step estimator.

{phang2}
{cmd:b0} and {cmd:bp} specify Rukhin's B0 and BP estimators, respectively.

{phang2}
{cmd:bs}, {cmd:bt}, or {cmd:gamma} specifies the Biggerstaff-Tweedie
approximate gamma model.

{phang2}
{cmd:eb}, {cmd:gq}, {cmd:genq}, {cmd:mp}, or {cmd:q} specifies the
empirical Bayes estimator, also known as the generalized Q estimator and
the Mandel-Paule estimator.

{phang2}
{cmd:ml} specifies the maximum likelihood estimator.

{phang2}
{cmd:pl} specifies the profile likelihood model.

{phang2}
{cmd:reml} specifies the restricted maximum-likelihood estimator.

{phang2}
{cmd:sa}[{cmd:, isq(}{it:real}{cmd:)}] specifies the
sensitivity analysis with a fixed user-specified value of I^2 (between
0 and 1, with default {cmd:isq(0.8)}).

{pmore}
Note that the approximate gamma, generalized Q, maximum likelihood,
profile likelihood, and restricted maximum-likelihood models require the
{cmd:mm_root()} function; and the bootstrapped DerSimonian-Laird model
requires
the {cmd:mm_bs()} and {cmd:mm_jk()} functions from the {cmd:moremata}
package ({cmd:ssc install moremata}).  The approximate gamma model also
requires the {cmd:integrate} command ({cmd:ssc install integrate}).

{phang}
{cmd:sortby(}{it:varname}|{cmd:_n)} allows user-specified ordering of
studies in the table and forest plot.  The default ordering is by
{it:study_id}.  Note that {opt sortby()} does not alter the data in memory.

{pmore}
Specify {cmd:sortby(_n)} to order the studies by their first appearance in the
data by using the current sort order.

{dlgtab:Forest plots}

{phang}
{cmd:lcols(}{help ipdmetan##cols_info:{it:cols_info}}{cmd:)} and
{cmd:rcols(}{help ipdmetan##cols_info:{it:cols_info}}{cmd:)} define
columns of additional data to be presented to the left or right of the
forest plot.  These options are carried over from {helpb metan}; however, for
IPD, they must first be generated from the existing dataset and thus
require more complex syntax.

{pmore}
Specifying {it:newname} is necessary only when the
name of the variable in the {cmd:forestplot} dataset is important.  For
example, users may have an aggregate dataset with a variable containing
data equivalent to an {it:item} and want all such data (whether IPD
or aggregate) to appear in one column in the forest plot.  To
achieve this, the user must specify {it:newname} as the name of the relevant variable
in the aggregate dataset.  The variables in the IPD and
aggregate datasets cannot have conflicting formats (for example, string
and numeric); if they do, a {helpb merge} error will be returned.

{pmore}
Note that {it:item} can be an existing string variable, in which case
the first nonempty observation for each study will be used, and the
{it:item} will not be displayed alongside overall or subgroup pooled
estimates.  To force this behavior for a numeric variable, the user must first convert item to string format by using {helpb recode} or {helpb tostring}.

{pmore}
{cmd:lcols()} and {cmd:rcols()} can also be supplied directly to
{cmd:forestplot}, but only as a list of existing variable names.

{phang}
{cmd:plotid(}{it:varname}|{cmd:_BYAD}[{cmd:, list nograph}]{cmd:)}
specifies one or more categorical variables to form a series of groups
of observations in which specific aspects of plot rendition may be
affected using {it:plot}[{it:#}]{cmd:opts}.  The groups of observations
will automatically be assigned ordinal labels (1, 2, ...) on the basis of
the ordering of {it:varlist}.  Note that {cmd:plotid()} does not alter the
data in memory.

{pmore}
For further details of this option and the {opt list} and {opt nograph}
suboptions, see {cmd:forestplot}.

{phang}
{cmd:ovstat(q)} displays Q statistics instead of I-squared statistics.

{phang}
{cmd:saving(}{it:{help filename}}[{cmd:, replace}
{cmd:stacklabel}]{cmd:)} saves the forest plot "results set" created by
{cmd:ipdmetan} in a dataset for further use or manipulation.  See
{cmd:forestplot} for further details.

{phang2}
{opt replace} overwrites {it:filename}.

{phang2}
{opt stacklabel} takes the {help label:variable label}
from the leftmost column variable (usually {it:study_id}), which would
usually appear outside the plot region as the column heading, and copies
it into a new first row in {it:filename}.  This allows multiple such
datasets to be {helpb append}ed without this information being
overwritten.

{phang}
{opt forestplot(forestplot_options)} specifies other options to pass to 
{helpb forestplot}.

{dlgtab:Combined IPD/aggregate data analysis}

{phang}
{cmd:ad(}{it:{help filename}} {ifin}{cmd:,} {opt vars(namelist)} 
[{help ipdmetan##ad_options:{it:ad_options}}]{cmd:)}
allows aggregate (summary) data to be included in the analysis alongside IPD,
for example, if some studies do not have IPD available.
If {cmd:ad()} is specified, {it:filename} and {opt vars(namelist)} are
required.

{phang2}
{opt vars(namelist)} contains the names of variables (within {it:filename})
containing the effect size and either a standard error or lower and upper 95%
confidence limits on the linear scale.
If confidence limits are supplied, they must be derived from a normal
distribution, or the pooled result will be incorrect (see {helpb admetan}).

{phang2}
{it:ad_options} can be the following:

{phang2}
{opt byad} specifies that IPD and aggregate data be treated as
subgroups rather than as a single set of estimates.

{phang2}
{opt npts(varname)} allows participant numbers (stored in {it:varname} within
{it:filename}) to be displayed in tables and forest plots.

{pmore}
Note that subgroups can be analyzed in the same way as IPD -- that
is, with the {opt by(varname)} option to {cmd:ipdmetan}.  {it:varname}
can be found in either the data in memory (IPD) or the aggregate
dataset or both.


{title:Examples}

{pstd}
Setup{p_end}
{phang2}
{cmd:. use ipdmetan_example}{p_end}
{phang2}
{cmd:. stset tcens, fail(fail)}

{pstd}
Basic use{p_end}
{phang2}
{cmd:. ipdmetan, study(trialid) hr by(region) nograph: stcox trt, strata(sex)}{p_end}

{pstd}
Use of {cmd:plotid()}{p_end}
{phang2}
{cmd:. ipdmetan, study(trialid) hr by(region) plotid(region) forestplot(favours(Favours treatment # Favours control) box1opts(mcolor(red)) ci1opts(lcolor(red) rcap) box2opts(mcolor(blue)) ci2opts(lcolor(blue)))}{cmd:: stcox trt, strata(sex)}

{pstd}
Treatment-covariate interactions{p_end}
{phang2}
{cmd:. ipdmetan, study(trialid) interaction hr keepall forestplot(favours("Favours greater treatment effect" "with higher disease stage" # "Favours greater treatment effect" "with lower disease stage") boxscale(200) fp(1)): stcox trt##c.stage}

{pstd}
Random effects{p_end}
{phang2}
{cmd:. ipdmetan, study(trialid) hr nograph re: stcox trt, strata(sex)}{p_end}
{phang2}
{cmd:. ipdmetan, study(trialid) hr nograph re(q): stcox trt, strata(sex)}

{pstd}
Aggregate data setup: Create aggregate dataset from IPD dataset{p_end}
{phang2}
{cmd:. quietly ipdmetan, study(trialid) hr nograph saving(region2.dta): stcox trt if region==2, strata(sex)}{p_end}
{phang2}
{cmd:. clonevar _STUDY = trialid}{p_end}

{pstd}
Include aggregate data in the analysis{p_end}
{phang2}
{cmd:. ipdmetan, study(_STUDY) hr ad(region2.dta if _USE==1, vars(_ES _seES) npts(_NN) byad) nooverall}{cmd:: stcox trt if region==1, strata(sex)}

{pstd}
Use of non e-class commands and {cmd:rcols()}: Peto log-rank analysis{p_end}
{phang2}
{cmd:. ipdmetan (u[1,1]/V[1,1]) (1/sqrt(V[1,1])), study(trialid) rcols((u[1,1]) %5.2f "o-E(o)" (V[1,1]) %5.1f "V(o)") by(region) plotid(region) hr forestplot(nostats nowt favours(Favours treatment # Favours control)): sts test trt, mat(u V)}


{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:ipdmetan} stores the following in {cmd:r()} (with some variation):

{synoptset 18 tabbed}{...}
{p2col 5 25 29 2: Scalars}{p_end}
{synopt:{cmd:r(k)}}number of included studies k{p_end}
{synopt:{cmd:r(n)}}number of included participants{p_end}
{synopt:{cmd:r(eff)}}overall pooled effect size{p_end}
{synopt:{cmd:r(se_eff)}}standard error of pooled effect size{p_end}
{synopt:{cmd:r(Q)}}Cochran Q heterogeneity statistic (on k-1 degrees of
freedom){p_end}
{synopt:{cmd:r(tausq)}}between-study variance tau-squared{p_end}
{synopt:{cmd:r(sigmasq)}}average within-study variance{p_end}
{synopt:{cmd:r(Isq)}}heterogeneity measure I-squared{p_end}
{synopt:{cmd:r(HsqM)}}heterogeneity measure H-squared (Mittlböck
modification){p_end}

{synoptset 18 tabbed}{...}
{p2col 5 25 29 2: Macros}{p_end}
{synopt:{cmd:r(command)}}full estimation command line {p_end}
{synopt:{cmd:r(cmdname)}}estimation command name{p_end}
{synopt:{cmd:r(estvar)}}name of pooled coefficient{p_end}
{synopt:{cmd:r(re_model)}}random-effects model used{p_end}

{synoptset 18 tabbed}{...}
{p2col 5 25 29 2: Matrices}{p_end}
{synopt:{cmd:r(coeffs)}}matrix of study and subgroup identifiers, effect
coefficients, numbers of participants, and weights{p_end}

{pstd}
Certain iterative random-effects models can save the following additional
results (see {helpb mf_mm_root} for interpretations of convergence success
values).

{synoptset 18 tabbed}{...}
{p2col 5 25 29 2: Scalars}{p_end}
{synopt:{cmd:r(tsq_var)}}estimated variance of tau-squared{p_end}
{synopt:{cmd:r(tsq_lci)}}lower confidence limit for tau-squared{p_end}
{synopt:{cmd:r(tsq_uci)}}upper confidence limit for tau-squared{p_end}
{synopt:{cmd:r(rc_tausq)}}convergence of tau-squared point estimate{p_end}
{synopt:{cmd:r(rc_tsq_lci)}}convergence of tau-squared lower confidence
limit{p_end}
{synopt:{cmd:r(rc_tsq_uci)}}convergence of tau-squared upper confidence
limit{p_end}
{synopt:{cmd:r(rc_eff_lci)}}convergence of effect-size lower confidence
limit{p_end}
{synopt:{cmd:r(rc_eff_uci)}}convergence of effect-size upper confidence
limit{p_end}


{title:Acknowledgments}

{pstd}
I thank the authors of {helpb metan}, the command upon which this code is based,
particularly Ross Harris for his comments and good wishes.


{title:Author}

{pstd}David Fisher{p_end}
{pstd}MRC Clinical Trials Unit at University College London{p_end}
{pstd}London, UK{p_end}
{pstd}{browse "mailto:d.fisher@ucl.ac.uk":d.fisher@ucl.ac.uk}{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 15, number 2: {browse "http://www.stata-journal.com/article.html?article=st0384":st0384}{p_end}
