\name{print.survfit}
\alias{print.survfit}
\title{
Print a Short Summary of a Survival Curve
}
\description{
  Print number of observations, number of events, the restricted
  mean survival and its standard error, and the median survival with confidence limits for the median.
}
\usage{
\method{print}{survfit}(x, scale=1,digits = max(options()$digits - 4,
3), print.n=getOption("survfit.print.n"), show.rmean=getOption("survfit.print.mean"),...)
}
\arguments{
\item{x}{
the result of a call to the \code{survfit} function.
}\item{print.n}{What to use for number of subjects (see below)}
\item{digits}{Number of digits to print}
\item{scale}{
a numeric value to rescale the survival time, e.g.,
if the input data to survfit were in days,
\code{scale=365} would scale the printout to years.
}\item{show.rmean}{Show the restricted mean survival?}\item{...}{other unused arguments}}
\value{
\code{x}, with the invisible flag set.
}
\section{Side Effects}{
The number of observations (see Details), the number of events, the
restricted mean survival (\code{rmean}) and its standard error, and the median survival
with its confidence interval are printed.  If there are multiple curves, there is one line of output for each.
}
\details{
The restricted mean (\code{rmean}) and its standard error \code{se(rmean)} are based on a truncated estimator.  If the
last observation(s) is not a death, then the survival curve estimate does not
go to zero and the mean survival time cannot be estimated. Instead, the
quantity reported is the mean of survival restricted to the time before
the last censoring. When the last censoring time is not random this
quantity is occasionally of interest.

Any randomness in the last censoring time is not taken into account in
computing the standard error of the restricted mean.  The restricted
mean is shown mainly for compatibility with S, and may not be shown by
default in future versions.


The median and its confidence interval are defined by drawing a horizontal
line at 0.5 on the plot of the survival curve and its confidence bands.
The intersection of the line with the lower CI band defines the lower limit
for the median's interval, and similarly for the upper band.  If any of the
intersections is not a point, then we use the smallest point of intersection,
e.g., if the survival curve were exactly equal to 0.5 over an
interval. 

The "number of observations" is not well-defined for counting process
data. Previous versions of this code used the number at risk at the
first time point. This is misleading if many individuals enter late or
change strata. The original S code for the current version uses the
number of records, which is misleading when the counting process data
actually represent a fixed cohort with time-dependent covariates.

Four possibilities are provided, controlled by \code{print.n} or by
\code{options(survfit.print.n)}: \code{"none"} prints \code{NA},
\code{"records"} prints the number of records, \code{"start"} prints the
number at the first time point and \code{"max"} prints the maximum
number at risk. The initial default is \code{"start"}.

}

\seealso{
\code{\link{summary.survfit}}, \code{\link{survfit.object}}, \code{\link{survfit}}
}
\examples{
##effect of print.n and show.rmean
data(heart)
a<-coxph(Surv(start,stop,event)~age+strata(transplant),data=heart)
b<-survfit(a)
print(b,print.n="none")
print(b,print.n="records")
print(b,print.n="start")
print(b,print.n="max")
print(b,show.rmean=FALSE)
}

\keyword{survival}
