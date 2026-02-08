{smcl}
{* *! version 1.0.0}{...}
{viewerjumpto "Syntax" "dlm##syntax"}{...}
{viewerjumpto "Description" "dlm##description"}{...}
{viewerjumpto "Options" "dlm##options"}{...}
{viewerjumpto "Stored results" "dlm##results"}{...}
{viewerjumpto "Examples" "dlm##examples"}{...}
{viewerjumpto "References" "dlm##references"}{...}
{title:Title}

{p2colset 5 15 17 2}{...}
{p2col:{cmd:dlm} {hline 2}}Distributed lag model estimation (Schmidheiny & Siegloch 2023){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:dlm} {depvar} {ifin}{cmd:,}
{opt exp:osure(varname)}
{opt unit(varname)}
{opt time(varname)}
{opt from(#)}
{opt to(#)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {opt exp:osure(varname)}}exposure/treatment variable (typically binary){p_end}
{p2coldent:* {opt unit(varname)}}panel unit identifier{p_end}
{p2coldent:* {opt time(varname)}}time period variable{p_end}
{p2coldent:* {opt from(#)}}starting relative period (negative integer, e.g. -3){p_end}
{p2coldent:* {opt to(#)}}ending relative period (positive integer, e.g. 3){p_end}
{synopt:{opt ref(#)}}reference period; default is {cmd:ref(-1)}{p_end}
{synopt:{opt cov:ariates(varlist)}}additional covariates{p_end}
{synopt:{opt addl_fes(varlist)}}additional fixed effects beyond unit and time{p_end}
{synopt:{opt ver:bose}}display detailed progress{p_end}
{synoptline}
{p 4 6 2}* required options{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dlm} estimates distributed lag models as described in
Schmidheiny and Siegloch (2023, {it:Journal of Applied Econometrics}).

{pstd}
The DLM approach regresses the outcome on leads and lags of the exposure
variable (the {it:gamma} coefficients), then transforms these into
cumulative treatment effects (the {it:beta} coefficients) via cumulative
summation. The betas are mathematically equivalent to event-study coefficients
from a binned-endpoint TWFE regression.

{pstd}
Unit and time fixed effects are absorbed via {cmd:reghdfe}, and standard
errors are clustered at the unit level. The gamma-to-beta transformation
properly propagates uncertainty through the variance-covariance matrix.


{marker options}{...}
{title:Options}

{phang}
{opt exposure(varname)} specifies the treatment/exposure variable. This is
typically a binary indicator (0/1) that turns on at treatment onset.

{phang}
{opt unit(varname)} specifies the panel unit identifier (e.g., individual,
firm, county).

{phang}
{opt time(varname)} specifies the time period variable. Must be numeric and
evenly spaced.

{phang}
{opt from(#)} specifies the earliest relative time period to estimate.
Must be a negative integer (e.g., {cmd:from(-3)} estimates effects from
3 periods before treatment).

{phang}
{opt to(#)} specifies the latest relative time period to estimate.
Must be a positive integer (e.g., {cmd:to(3)} estimates effects up to
3 periods after treatment).

{phang}
{opt ref(#)} specifies the reference (omitted) period. Default is -1
(one period before treatment). Must be between {cmd:from()} and {cmd:to()}.

{phang}
{opt covariates(varlist)} specifies additional control variables to include
in the regression.

{phang}
{opt addl_fes(varlist)} specifies additional fixed effects beyond the default
unit and time FEs. These are passed to {cmd:reghdfe}'s {cmd:absorb()} option.

{phang}
{opt verbose} displays detailed progress information during estimation.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:dlm} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters{p_end}
{synopt:{cmd:e(from)}}starting relative period{p_end}
{synopt:{cmd:e(to)}}ending relative period{p_end}
{synopt:{cmd:e(ref_period)}}reference period{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:dlm}{p_end}
{synopt:{cmd:e(outcome)}}outcome variable name{p_end}
{synopt:{cmd:e(exposure)}}exposure variable name{p_end}
{synopt:{cmd:e(unit)}}unit variable name{p_end}
{synopt:{cmd:e(time)}}time variable name{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(betas)}}K x 5 matrix with columns: time_to_event, coef, se, ci_lo, ci_hi{p_end}
{synopt:{cmd:e(gamma)}}1 x J vector of raw DLM (gamma) coefficients{p_end}
{synopt:{cmd:e(gamma_V)}}J x J variance-covariance matrix of gamma coefficients{p_end}

{pstd}
The {cmd:e(betas)} matrix includes a row for the reference period with
coef=0 and se=0, making it convenient for plotting.


{marker examples}{...}
{title:Examples}

{pstd}Basic estimation:{p_end}
{phang2}{cmd:. dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)}{p_end}

{pstd}With covariates and additional fixed effects:{p_end}
{phang2}{cmd:. dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) covariates(x1 x2) addl_fes(state)}{p_end}

{pstd}Change reference period:{p_end}
{phang2}{cmd:. dlm outcome, exposure(post) unit(unit) time(time) from(-5) to(5) ref(-2)}{p_end}

{pstd}Access results after estimation:{p_end}
{phang2}{cmd:. matrix list e(betas)}{p_end}
{phang2}{cmd:. display e(N)}{p_end}

{pstd}Generate test data and run:{p_end}
{phang2}{cmd:. dlm_gen_data, n_groups(500) n_times(20) seed(42) clear}{p_end}
{phang2}{cmd:. dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)}{p_end}


{marker references}{...}
{title:References}

{phang}
Schmidheiny, K. and S. Siegloch. 2023.
On event studies and distributed-lag models: Equivalence, generalization
and practical implications.
{it:Journal of Applied Econometrics} 38(5): 695-713.
{p_end}


{title:Also see}

{psee}
{space 2}Help: {helpb reghdfe}
{p_end}
