\name{survfit}
\alias{survfit}
\alias{[.survfit}
\alias{survfit.km}
\alias{survfit.coxph.null}
\alias{survfit.coxph}
\alias{basehaz}
\title{
Compute a Survival Curve for Censored Data
}
\description{
Computes an estimate of a survival curve for censored data
using either the Kaplan-Meier or the Fleming-Harrington method
or computes the predicted survivor function for a Cox proportional
hazards model.
}
\synopsis{survfit(formula, data, weights, subset, na.action, ...)}
\usage{
survfit(formula, data, weights, subset, na.action, 
        newdata, individual=F, conf.int=.95, se.fit=T, 
        type=c("kaplan-meier","fleming-harrington", "fh2"),
        error=c("greenwood","tsiatis"),
        conf.type=c("log","log-log","plain","none"),
        conf.lower=c("usual", "peto", "modified"))
\method{[}{survfit}(x,...,drop=FALSE)
basehaz(fit,centered=TRUE)
}
\arguments{
\item{formula}{
A formula object or a \code{coxph} object.
If a formula object is supplied it must have a \code{Surv} object as the 
response on the left of the \code{~} operator and, if desired, terms 
separated by + operators on the right.
One of the terms may be a \code{strata} object.  For a single survival curve
the \code{"~ 1"} part of the formula is not required.
}
\item{data}{
a data frame in which to interpret the variables named in the formula,
or in the \code{subset} and the \code{weights} argument.
}
\item{weights}{
The weights must be nonnegative and it is strongly recommended that 
they be strictly positive, since zero weights are ambiguous, compared
to use of the \code{subset} argument.
}
\item{subset}{
expression saying that only a subset of the rows of the data
should be used in the fit.
}
\item{na.action}{
a missing-data filter function, applied to the model frame, after any
\code{subset} argument has been used.
Default is \code{options()$na.action}.
}
\item{newdata}{
a data frame with the same variable names as those that appear
in the \code{coxph} formula.  Only applicable when \code{formula} is a \code{coxph} object.
The curve(s) produced will be representative of a cohort who's
covariates correspond to the values in \code{newdata}.
Default is the mean of the covariates used in the \code{coxph} fit.
}
\item{individual}{
a logical value indicating whether the data frame represents different
time epochs for only one individual (T), or whether multiple rows indicate
multiple individuals (F, the default).  If the former only one curve
will be produced; if the latter there will be one curve per row in \code{newdata}.
}
\item{conf.int}{
the level for a two-sided confidence interval on  the survival curve(s).
Default is 0.95.
}
\item{se.fit}{
a logical value indicating whether standard errors should be
computed.  Default is \code{TRUE}.
}
\item{type}{
a character string specifying the type of survival curve.
Possible values are \code{"kaplan-meier"}, \code{"fleming-harrington"} or \code{"fh2"}
if a formula is given
and \code{"aalen"} or \code{"kaplan-meier"} if the first argument is a \code{coxph} object,
(only the first two characters are necessary).
The default is \code{"aalen"} when a \code{coxph} object is given,
and it is \code{"kaplan-meier"} otherwise.
}
\item{error}{
either the string \code{"greenwood"} for the Greenwood formula or
\code{"tsiatis"} for the Tsiatis formula, (only the first character is
necessary).  The default is \code{"tsiatis"} when a \code{coxph} object is
given, and it is \code{"greenwood"} otherwise.
}
\item{conf.type}{
One of \code{"none"}, \code{"plain"}, \code{"log"} (the default), or \code{"log-log"}.  Only
enough of the string to uniquely identify it is necessary.
The first option causes confidence intervals not to be
generated.  The second causes the standard intervals
\code{curve +- k *se(curve)}, where k is determined from
\code{conf.int}.  The log option calculates intervals based on the
cumulative hazard or log(survival). The last option bases
intervals on the log hazard or log(-log(survival)).  These
last will never extend past 0 or 1.
}
\item{conf.lower}{
controls modified lower limits to the curve,
the upper limit remains unchanged.  The modified lower limit
is based on an 'effective n' argument.  The confidence
bands will agree with the usual calculation at each death time, but unlike
the usual bands the confidence interval becomes wider at each censored
observation.  The extra width is obtained by multiplying the usual
variance by a factor m/n, where n is the number currently at risk and
m is the number at risk at the last death time.  (The bands thus agree
with the un-modified bands at each death time.)
This is especially useful for survival curves with a long flat tail.


The Peto lower limit is based on the same 'effective n' argument as the
modified limit, but also replaces the usual Greenwood variance term with
a simple approximation.  It is known to be conservative.
}
\item{x}{a \code{survfit} object}
\item{fit}{a \code{coxph} object}
\item{centered}{Compute the baseline hazard at the covariate mean rather
than at zero?}
\item{drop}{Only \code{FALSE} is supported}
\item{...}{Other arguments for future expansion}
}
\value{
a \code{survfit} object; see the help on \code{survfit.object} for
details. Methods defined for \code{survfit} objects are provided for
\code{print}, \code{plot}, \code{lines}, and \code{points}.

For \code{basehaz}, a dataframe with the baseline hazard, times, and
strata.

The \code{"["} method returns a \code{survfit} object giving survival
for the selected groups.
}
\details{
Actually, the estimates used are the Kalbfleisch-Prentice
(Kalbfleisch and Prentice, 1980, p.86) and the Tsiatis/Link/Breslow,
which reduce to the Kaplan-Meier and Fleming-Harrington estimates,
respectively, when the weights are unity.  When curves are fit for a
Cox model, subject weights of \code{exp(sum(coef*(x-center)))} are used, 
ignoring any value for \code{weights} input by the user.  There is also an extra
term in the variance of the curve, due to the variance ofthe coefficients and
hence variance in the computed weights.


The Greenwood formula for the variance is a sum of terms
d/(n*(n-m)), where d is the number of deaths at a given time point, n
is the sum of \code{weights} for all individuals still at risk at that time, and
m is the sum of \code{weights} for the deaths at that time.  The
justification is based on a binomial argument when weights are all
equal to one; extension to the weighted case is ad hoc.  Tsiatis
(1981) proposes a sum of terms d/(n*n), based on a counting process
argument which includes the weighted case.


The two variants of the F-H estimate have to do with how ties are handled.
If there were 3 deaths out of 10 at risk, then the first would increment
the hazard by 3/10 and the second by 1/10 + 1/9 + 1/8.  For curves created
after a Cox model these correspond to the Breslow and Efron estimates,
respectively, and the proper choice is made automatically.
The \code{fh2} method will give results closer to the Kaplan-Meier.


Based on the work of Link (1984), the log transform is expected to produce
the most accurate confidence intervals.  If there is heavy censoring, then
based on the work of Dorey and Korn (1987) the modified estimate will give
a more reliable confidence band for the tails of the curve.
}
\references{
Dorey, F. J. and Korn, E. L. (1987).  Effective sample sizes for confidence
intervals for survival probabilities.  \emph{Statistics in Medicine} 6, 679-87.


Fleming, T. H. and Harrington, D.P. (1984).  Nonparametric estimation of the
survival distribution in censored data.  \emph{Comm. in Statistics} 13, 2469-86.


Kalbfleisch, J. D. and Prentice, R. L. (1980). 
\emph{The Statistical Analysis of Failure Time Data.}
Wiley, New York.


Link, C. L. (1984). Confidence intervals for the survival
function using Cox's proportional hazards model with 
covariates.  \emph{Biometrics} 40, 601-610.


Tsiatis, A. (1981). A large sample study of the estimate
for the integrated hazard function in Cox's regression
model for survival data. \emph{Annals of Statistics} 9, 93-108.
}
\seealso{
\code{\link{print.survfit}}, \code{\link{plot.survfit}},
\code{\link{lines.survfit}}, \code{\link{summary.survfit}}, \code{\link{survfit.object}}
\code{\link{coxph}}, \code{\link{Surv}}, \code{\link{strata}}.}
\examples{
#fit a Kaplan-Meier and plot it
fit <- survfit(Surv(time, status) ~ x, data=aml)
plot(fit)
# plot only 1 of the 2 curves from above
plot(fit[2])

#fit a cox proportional hazards model and plot the 
#predicted survival curve
fit <- coxph( Surv(futime,fustat)~resid.ds+rx+ecog.ps,data=ovarian)
plot( survfit( fit))
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
