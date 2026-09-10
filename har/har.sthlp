{smcl}
{* *! version 13.1  11Apr2017}{...}
{viewerjumpto "Title" "har##title"}{...}
{viewerjumpto "Syntax" "har##syntax"}{...}
{viewerjumpto "Description" "har##description"}{...}
{viewerjumpto "Options" "har##options"}{...}
{viewerjumpto "Examples" "har##examples"}{...}
{viewerjumpto "Stored results" "har##results"}{...}
{viewerjumpto "References" "har##references"}{...}
{viewerjumpto "Authors" "har##authors"}{...}
{viewerjumpto "Also see" "har##alsosee"}{...}
{cmd:help har}{right: ({browse "http://www.stata-journal.com/article.html?article=st0548":SJ18-4: st0548})}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 12 14 2}{...}
{p2col:{cmd:har} {hline 2}}Regression with heteroskedasticity- and
autocorrelation-robust (HAR) standard errors{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 11 2}
{cmd:har} {depvar} [{it:{help varlist:varlist1}}]
{cmd:(}{it:{help varlist:varlist2}} {cmd:=}
        {it:{help varlist:varlist_iv}}{cmd:)} {ifin}{cmd:,} {opt kernel(string)} [{it:options}]

{synoptset 18 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {opt kernel(string)}}set the type of kernels:
Bartlett, Parzen, quadratic spectral (QS), or orthonormal series{p_end}
{synopt:{opt nocon:stant}}suppress the constant term{p_end}
{synopt:{opt l:evel(#)}}set the confidence level; default is {cmd:level(95)}{p_end}
{synoptline}
{p2colreset}{...}
{pstd}
* {cmd:kernel()} is required.{p_end}
{p 4 6 2}
You must {cmd:tsset} your data before using {opt har}; see
{helpb tsset:[TS] tsset}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
The {cmd:har} command estimates an instrumental-variable (IV) regression in the presence
of heteroskedasticity and autocorrelation.  Inferences are based on the
fixed-smoothing asymptotics.  The command is based on Sun (2013, 2014).  Sun
(2013) develops heteroskedasticity- and autocorrelation-robust F and t tests
using an orthonormal-series long-run variance matrix estimator.  The number of
orthonormal bases is selected by minimizing the type II error of the
associated test while controlling for its type I error.  Sun (2014) introduces
a new and easy-to-use asymptotic F test (and t test) based on kernel long-run
variance matrix estimators.  The proposed bandwidth selection rule is testing
optimal and similar to Sun (2013).  {cmd:har} allows for three types of
kernels: the Bartlett, Parzen, and QS kernels.

{pstd}
For the sake of simplicity, we refer to the orthonormal-series long-run
variance estimator in Sun (2013) as the kernel long-run variance estimator
with kernel "orthonormal series".  In total, {cmd:har} allows for four types
of kernels: Bartlett, Parzen, QS, and orthonormal series.


{marker options}{...}
{title:Options}

{phang}
{opt kernel(string)} sets the type of kernel.  For the Bartlett kernel, any of
the four usages -- {cmd:kernel(bartlett)}, {cmd:kernel(BARTLETT)},
{cmd:kernel(B)}, or {cmd:kernel(b)} -- produce the same results.  Similarly,
for the Parzen, QS, and orthonormal-series long-run variance estimators, we
can use any of the respective choices: ({cmd:PARZEN}, {cmd:parzen}, {cmd:P},
{cmd:p}), ({cmd:QUADRATIC}, {cmd:quadratic}, {cmd:Q}, {cmd:q}), and
({cmd:ORTHOSERIES}, {cmd:orthoseries}, {cmd:O}, {cmd:o}).  {cmd:kernel()} is
required.

{phang}
{opt noconstant} suppresses the constant term.

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for
confidence intervals.  The default is {cmd:level(95)}.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse idle2}{p_end}
{phang2}{cmd:. tsset time}

{pstd}Regression with HAR standard errors using the Bartlett kernel{p_end}
{phang2}{cmd:. har usr idle, kernel(bartlett)}

{pstd}Regression with HAR standard errors using the Parzen kernel{p_end}
{phang2}{cmd:. har usr idle, kernel(parzen)}

{pstd}Regression with HAR standard errors using the QS kernel{p_end}
{phang2}{cmd:. har usr idle, kernel(quadratic)}

{pstd}Regression with HAR standard errors using "orthonormal series" {p_end}
{phang2}{cmd:. har usr idle, kernel(orthoseries)}

{pstd}Regression with HAR standard errors using the Bartlett kernel (99% confidence 
interval with no constant included in the regression){p_end}
{phang2}{cmd:. har usr idle, kernel(bartlett) noconstant level(99)} 


{marker results}{...}
{title:Stored results}

{pstd}
In addition to the standard stored results from {helpb ivregress}, {cmd:har}
stores the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sF)}}adjusted F statistic (only for orthonormal series){p_end}
{synopt:{cmd:e(ssdf)}}second degrees of freedom (only for orthonormal series){p_end}
{synopt:{cmd:e(kopt)}}data-driven optimal number K of orthonormal bases (only for orthonormal series){p_end}
{synopt:{cmd:e(kF)}}adjusted F statistic (for Bartlett, Parzen, and QS){p_end}
{synopt:{cmd:e(ksdf)}}second degrees of freedom (for Bartlett, Parzen, and QS){p_end}
{synopt:{cmd:e(lag)}}data-driven optimal truncation lag (for Bartlett, Parzen, and QS){p_end}
{synopt:{cmd:e(fdf)}}first degrees of freedom{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:har}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(carg)}}{cmd:nocons} or {cmd:""} if specified{p_end}
{synopt:{cmd:e(varline)}}variable line as typed{p_end}
{synopt:{cmd:e(kerneltype)}}kernel in the estimation{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(sstderr)}}adjusted standard error for each individual coefficient (only for orthonormal series){p_end}
{synopt:{cmd:e(sdf)}}degrees of freedom of the approximating t distribution (only for orthonormal series){p_end}
{synopt:{cmd:e(st)}}t statistic (only for orthonormal series){p_end}
{synopt:{cmd:e(sbetahat)}}IV coefficient vector (only for orthonormal series){p_end}
{synopt:{cmd:e(kbetahat)}}IV coefficient vector (for Bartlett, Parzen, and QS){p_end}
{synopt:{cmd:e(kstderr)}}adjusted standard error for each individual
coefficient (for Bartlett, Parzen, and QS){p_end}
{synopt:{cmd:e(kdf)}}degrees of freedom of the approximating t
distribution (for Bartlett, Parzen, and QS){p_end}
{synopt:{cmd:e(kt)}}t statistic (for Bartlett, Parzen, and QS){p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks the estimation sample{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker S2013}{...}
{phang}
Sun, Y. 2013. A heteroskedasticity and autocorrelation robust F test using an orthonormal series variance estimator.
{it:Econometrics Journal} 16: 1-26.

{marker SK2012}{...}
{phang}
------. 2014. Let's fix it: Fixed-b asymptotics versus small-b asymptotics in heteroskedasticity and autocorrelation robust inference.
{it:Journal of Econometrics} 178: 659-677.


{marker authors}{...}
{title:Authors}

{pstd}
Xiaoqing Ye{break}
School of Mathematics and Statistics{break}
South-Central University for Nationalities{break}
Wuhan, China{break}
yshtim@126.com

{pstd}
Yixiao Sun{break}
Department of Economics{break}
University of California, San Diego{break}
La Jolla, CA{break}
yisun@ucsd.edu


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 18, number 4: {browse "http://www.stata-journal.com/article.html?article=st0548":st0548}{p_end}
