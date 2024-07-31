{smcl}
{* *! version 1.1  David Fisher  23jul2014}{...}
{cmd:help forestplot}{right: ({browse "http://www.stata-journal.com/article.html?article=st0384":SJ15-2: st0384})}
{hline}

{title:Title}

{p2colset 5 19 21 2}{...}
{p2col:{cmd:forestplot} {hline 2}}Create forest plots from data currently in memory


{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:forestplot} [{varlist}] {ifin} [{cmd:, }
{it:options}]

{pstd}
where {it:varlist} is

{p 16 24 2}
{it:ES} {it:lci} {it:uci} [{it:wt}] [{it:use}]

{synoptset 35 tabbed}{...}
{synopthdr}
{synoptline}
{syntab: Main}
{synopt :{opt dataid(varname)}}define groups of observations forming a complete forest plot{p_end}
{synopt :{opt dp(#)}}set number of decimal places to display{p_end}
{synopt :{it:{help eform_option}}}exponentiate effect sizes and confidence limits{p_end}
{synopt :{opt eff:ect(string)}}set title for "effect size" column in the output{p_end}
{synopt :{opt fav:ours(string)}}set x-axis labeling specific to forest plots{p_end}
{synopt :{opt lab:els(varname)}}specify variable containing labels, for example, subgroup headings, heterogeneity info, study names{p_end}
{synopt :{opt lcol:s(varlist)} {opt rcol:s(varlist)}}display columns of additional data{p_end}
{synopt :{opt nona:me}}suppress display of study names in left-hand column{p_end}
{synopt :{opt nonu:ll}}suppress null hypothesis line{p_end}
{synopt :{opt nulloff}}remove null hypothesis line{p_end}
{synopt :{opt noov:erall}, {opt nosu:bgroup}}suppress display of overall or subgroup pooled estimates{p_end}
{synopt :{opt nostat:s}, {opt nowt}}suppress display of effect estimates or weights in right-hand columns{p_end}
{synopt :{cmd:plotid(}{it:varname} [{cmd:, {ul:l}ist {ul:nogr}aph}]{cmd:)}}define groups of observations in which
to apply specific plot rendition options{p_end}
{synopt :{it:{help twoway_options}}}specify other Stata {cmd:twoway} graph options{p_end}

{syntab: Fine-tuning}
{synopt :{opt noadj:ust}}suppress "space-saving" adjustment to text or data placement{p_end}
{synopt :{opt ast:ext(#)}}specify percentage of the graph to be occupied by text{p_end}
{synopt :{opt box:scale(#)}}box-size scaling{p_end}
{synopt :{opt r:ange(numlist)}}set x-axis limits of data-plotting area independently of axis ticks or labels{p_end}
{synopt :{opt text:scale(#)}}font-size scaling for text display on graph{p_end}
{synopt :{cmdab:xlab:el(}{it:numlist}{cmd:, force)}}force x-axis limits of data-plotting area together with ticks or labels{p_end}

{syntab: Plot rendition}
{synopt :{opt classic}}use "classic" set of plot options, as in {helpb metan}{p_end}
{synopt :{opt interaction}}use "interaction" set of plot options{p_end}
{synopt :{ul:{it:plot}}{cmd:{ul:op}ts(}{it:plot_options}{cmd:)}}affect rendition of all observations{p_end}
{synopt :{ul:{it:plot#}}{cmd:opts(}{it:plot_options}{cmd:)}}affect rendition of observations in {it:#}th {cmd:plotid()} group{p_end}
{synoptline}
{pstd}
{it:plot} may be {cmd:box}, {cmd:ci}, {cmd:diam}, {cmd:oline}, {cmd:point}, {cmd:pci}, or {cmd:ppoint}.{p_end}
{pstd}
{it:plot_options} are either {it:{help marker_options}} or {it:{help line_options}}, as appropriate.

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:forestplot} creates forest plots from variables containing point estimates and lower or upper confidence limits.
The meta-analysis program {helpb ipdmetan} and its sister programs {helpb admetan} and {helpb ipdover}
call {cmd:forestplot} by default, and they pass on any relevant options.
They can also save datasets with additional information that {cmd:forestplot} recognizes.
Alternatively, {cmd:forestplot} can use only data currently in memory.

{pstd}
{cmd:forestplot} requires three variables corresponding to {it:ES},
{it:lci}, and {it:uci} representing the effect size (point estimate)
and lower or upper confidence limits on the normal scale (that is, log odds
or log hazards, not exponentiated).  These can be supplied manually (in
the above order) using {it:varlist}.  Otherwise, {cmd:forestplot} will
expect to find variables in memory named {bf:_ES}, {bf:_LCI}, and
{bf:_UCI}.

{pstd}
{cmd:forestplot} will also check for variables corresponding to {it:wt}
and {it:use}, which represent, respectively, the weight (relative marker
size) and an indicator of the contents of each observation (study
effects, titles, spacing, description of heterogeneity, etc).  The
default names for these variables are {bf:_WT} and {bf:_USE},
respectively (although these can be overridden); {cmd:forestplot} will
assume they are constant if not found.

{pstd}
The values of the variable {it:use} are interpreted by {cmd:forestplot}
in the following way:{p_end}

{phang2}{cmd:0} = subgroup labels (headings){p_end}
{phang2}{cmd:1} = nonmissing study estimates{p_end}
{phang2}{cmd:2} = missing study estimates{p_end}
{phang2}{cmd:3} = subgroup pooled effects{p_end}
{phang2}{cmd:4} = description of between-subgroup heterogeneity{p_end}
{phang2}{cmd:5} = overall pooled effect{p_end}
{phang2}{cmd:6} = blank line{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt dataid(varname)} defines groups of observations forming a complete
forest plot.  The data in memory might come from multiple
separate meta-analyses whose forest plots must be plotted within
the same region (see {it:{help region_options}}).  Specifying
{cmd:dataid()} tells {cmd:forestplot} where the data from one
meta-analysis end and data from the next begin, and it results in correct placement
of the overall effect lines.  This option is unnecessary in
most circumstances.

{phang}
{opt dp(#)} specifies the number of decimal places at which to format the effect sizes.
This option is carried over from {helpb metan}, but it was previously undocumented.

{phang}
{it:{help eform_option}} specifies that the effect sizes and confidence limits
 be exponentiated.  The option also generates a heading for the effect-size column in the output (table and forest plot).

{pmore}
Note that {cmd:forestplot} expects effect sizes to be beta coefficients (that
is, on the linear scale).

{phang}
{opt effect(string)} specifies a heading for the effect-size column in the
output.  This overrides any heading generated by {it:{help eform_option}}.

{phang}
{opt favours(string)} applies a label saying something about the treatment
effect to either side of the graph (strings are separated by the # symbol).
This option is carried over from {helpb metan}.

{pmore}
Note that {opt favours()} and {opt xtitle()} use 
{it:{help axis_label_options}} rather than the usual 
{it:{help axis_title_options}}.  In the case of {opt favours()}, suboptions
must be passed to a separate {opt xmlabel()} option.  For example, 
{opt favours(string)} {cmd:xmlabel(, labgap(*1.5))}.

{phang}
{opt labels(varname)} specifies a text variable containing labels for
the left-hand side of the graph.  Examples include subgroup titles,
heterogeneity details, and study names.  You don't need to
specify this option if the variable {bf:_LABELS} exists and contains the
appropriate information.

{phang}
{opt lcols(varlist)} and {opt rcols(varlist)} define columns of additional
data to the left or right of the graph.  These options work in the same way as in {helpb metan} when specified directly to
{cmd:forestplot}.  (They have a different
syntax when specified to 
{helpb ipdmetan##forestplot_options:ipdmetan}.)

{pmore}
The first two columns on the right are automatically set to effect size and
weight, and the first on the left is set to study name or identifier (unless suppressed
using the options {opt nostats}, {opt nowt}, or {opt noname}, respectively).
{opt textsize()} can be used to fine-tune the size of the text to
achieve a satisfactory appearance.  The columns are labeled with the name of
the variable or macro.

{phang}
{opt noname} suppresses display of study names in left-hand column.

{phang}
{opt nonull} suppresses null hypothesis line.

{phang}
{opt nulloff} removes the null hypothesis line from the graph.  This option is
carried over from {helpb metan}.

{phang}
{cmd:nooverall} and {cmd:nosubgroup} suppress display of overall or
subgroup pooled estimates.

{phang}
{cmd:nostats} and {cmd:nowt} suppress display of effect estimates or
weights in right-hand columns.

{phang}
{cmd:plotid(}{it:varlist} [{cmd:, list nograph}]{cmd:)} specifies one or more
categorical variables to form a series of groups of observations in which
specific aspects of plot rendition may be affected using
{bf:box}{it:#}{bf:opts}, {bf:ci}{it:#}{bf:opts}, etc.  The groups of
observations will automatically be assigned ordinal labels (1, 2, ...) based
on the ordering of {it:varlist}.

{pmore}
The contents of each group can be inspected with the {bf:list} option.  In
complex situations, it may be helpful to view this list without creating the
plot itself; this can be achieved using the {bf:nograph} option.

{pmore}
Note that {opt plotid()} does not alter the placement or ordering of data
within the plot.

{phang}
{it:twoway_options} are other {cmd:twoway} graph options, as appropriate.

{dlgtab:Fine-tuning}

{pstd}
These options allow fine-tuning of the plot construction and text or data
placement.  Forest plots are nonstandard twoway plots and have features that
Stata graphs were not designed to accommodate (for example, columns of text
and a wide range of potential aspect ratios).  Occasionally,
{cmd:forestplot} will produce unsatisfactory results, which may be improved by
using one or more of these options.

{phang}
{opt noadjust} suppresses a calculation within {cmd:forestplot} that may
result in labeling text overlapping a pooled-estimate diamond.

{pmore}
This calculation, included in {cmd:metan}, takes advantage of
pooled-estimate diamonds having less width than individual study
estimates (although their labeling text is often longer) to "compress" the plot
and make it more aesthetic.

{phang}
{opt astext(#)} specifies the percentage of the graph to be occupied by text.
The default is {cmd:astext(50)}, and the percentage must be between
10-90.  This option is carried over from {cmd:metan}.

{phang}
{opt boxscale(#)} controls box scaling.  The default is {cmd:boxscale(100)}
(as in a percentage) and can be increased or decreased as such (for example,
80 or 120 for 20% smaller or larger, respectively).  This option is carried
over from {cmd:metan}.

{phang}
{opt range(numlist)} specifies the range of the x axis containing data,
independently of {it:{help axis_options}}, for the purposes of text placement.
For instance, a large blank space between the data and either the left or
right stats columns can be reduced or created.  Effect sizes or confidence
limits outside the range will be represented by off-scale arrows.

{phang}
{opt textscale(#)} specifies font size for text display on the graph.  The default
is {cmd:textscale(100)} (as in a percentage) and can be increased or decreased
as such (for example, 80 or 120 for 20% smaller or larger, respectively).
This option is carried over from {cmd:metan}.

{phang}
{cmd:xlabel(}{it:numlist}{cmd:, force)} with the {opt force} option operates
similarly to {opt range()} -- the smallest and largest values in {it:numlist}
become the range (unless {opt range()} itself is also specified).
{cmd:xlabel()} otherwise functions in the standard way.  This option is
carried over from {cmd:metan}, but with modifications.

{dlgtab:Plot rendition}

{phang}
{cmd:classic} uses the "classic" set of plot options, as in {cmd:metan}.

{phang}
{cmd:interaction} uses the "interaction" set of plot options.

{pstd}
These options specify suboptions for individual components forming the
forest plot, each of which is drawn using a separate {helpb twoway} command.
Any suboptions associated with a particular {cmd:twoway} command can be used
unless they would conflict with the appearance of the graph as a whole.  For
example, diamonds are plotted using the {helpb twoway_pcspike:twoway pcspike}
command, so {it:{help line_options}} can be used but not {opt horizontal} or
{opt vertical}.  These options are carried over from {cmd:metan} but with
modifications.

{pmore}
{cmd:boxopts(}{it:{help marker_options}}{cmd:)} and
{cmd:box}{it:#}{cmd:opts(}{it:marker_options}{cmd:)} affect the rendition of
weighted boxes representing point estimates and use options for a weighted
marker (for example, shape and color but not size).

{pmore}
{cmd:ciopts(}{it:{help line_options}} [{cmd:rcap}]{cmd:)} and
{cmd:ci}{it:#}{cmd:opts(}{it:line_options} [{cmd:rcap}]{cmd:)} affect the
rendition of confidence intervals. The additional option {cmd:rcap} requests
capped spikes.

{pmore}
{cmd:diamopts(}{it:{help line_options}}{cmd:)} and
{cmd:diam}{it:#}{cmd:opts(}{it:line_options}{cmd:)} affect the rendition of
diamonds representing pooled estimates.

{pmore}
{cmd:olineopts(}{it:{help line_options}}{cmd:)} and
{cmd:oline}{it:#}{cmd:opts(}{it:line_options}{cmd:)} affect the rendition of
overall effect lines.

{pmore}
{cmd:pointopts(}{it:{help marker_options}}{cmd:)} and
{cmd:point}{it:#}{cmd:opts(}{it:marker_options}{cmd:)} affect the rendition of
unweighted point-estimate markers (for example, to clarify the precise point
within a larger weighted box).

{pstd}
To represent pooled estimates by point estimates plus confidence intervals
(albeit with different {it:plot_options}) rather than by diamonds,
one can use the following options as replacements for {cmd:diamopts()} or
{cmd:diam}{it:#}{cmd:opts()}:

{pmore}
{cmd:pciopts(}{it:{help line_options}}{cmd:)} and
{cmd:pci}{it:#}{cmd:opts(}{it:line_options}{cmd:)} affect the rendition of
confidence intervals for pooled estimates.

{pmore}
{cmd:ppointopts(}{it:{help marker_options}}{cmd:)} and
{cmd:ppoint}{it:#}{cmd:opts(}{it:marker_options}{cmd:)} affect the rendition
of unweighted pooled estimate markers.


{title:Examples}

{pstd}
Setup{p_end}
{phang2}
{cmd:. use ipdmetan_example}{p_end}
{phang2}
{cmd:. stset tcens, fail(fail)}{p_end}
{phang2}
{cmd:. quietly ipdmetan, study(trialid) hr by(region) nograph saving(results.dta): stcox trt, strata(sex)}{p_end}
{phang2}
{cmd:. use results}

{pstd}
Use of {it:plot#}{cmd:opts()}.  (Note the use of {cmd:plotid(_BY)} instead of
{cmd:plotid(region)}: the reason is that the {cmd:by()} variable is given the generic name
{bf:_BY} in the results set.){p_end}
{phang2}
{cmd:. forestplot, hr favours(Favours treatment # Favours control) plotid(_BY) box1opts(mcolor(red)) ci1opts(lcolor(red) rcap) box2opts(mcolor(blue)) ci2opts(lcolor(blue))}

{pstd}
Replace weights with numbers of patients{p_end}
{phang2}
{cmd:. forestplot, hr favours(Favours treatment # Favours control) nowt rcols(_NN)}


{title:Acknowledgments}

{pstd}
I thank the authors of {helpb metan}, the command upon which this code is based;
I particularly thank Ross Harris for his comments and good wishes.  I also thank Vince Wiggins at StataCorp for advice.


{title:Author} 

{pstd}David Fisher{p_end}
{pstd}MRC Clinical Trials Unit at University College London{p_end}
{pstd}London, UK{p_end}
{pstd}{browse "mailto:d.fisher@ucl.ac.uk":d.fisher@ucl.ac.uk}{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 15, number 2: {browse "http://www.stata-journal.com/article.html?article=st0384":st0384}{p_end}

