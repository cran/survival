\name{print.survfit}
\alias{print.survfit}
\title{
Print a Short Summary of a Survival Curve
}
\description{
Print number of observations, number of events, the mean survival and its 
standard error, and the median survival with confidence limits for the 
median.
}
\usage{
\method{print}{survfit}(x, scale=1,digits = max(options()$digits - 4, 3),print.n=getOption("survfit.print.n"),...)
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
}\item{...}{other unused arguments}}
\value{
\code{x}, with the invisible flag set.
}
\section{Side Effects}{
the number of observations (see below), the number of events, the mean survival and its
standard error, and the median survival with its confidence interval are
printed.  If there are multiple curves, there is one line of output for each.
}
\details{
The mean and its variance are based on a truncated estimator.  That is, if the
last observation(s) is not a death, then the survival curve estimate does not
go to zero and the mean is undefined.  In such a case, the estimator is based
on an assumption that the true curve goes to zero just beyond the last
observed follow up time; it will systematically underestimate the true mean.


The median and its confidence interval are defined by drawing a horizontal
line at 0.5 on the plot of the survival curve and its confidence bands.
The intersection of the line with the lower CI band defines the lower limit
for the median's interval, and similarly for the upper band.  If any of the
intersections is not a point, then we use the smallest point of intersection,
e.g., if the survival curve were exactly equal to 0.5 over an interval.

The "number of observations" is not well-defined for counting process
data. Previous versions of this code used the number at risk at the
first time point. This is misleading if many individuals enter late or
change strata. The original S code for the current version uses the
number of records, which is misleading with time-dependent
covariates. Four possibilities are provided, controlled by
\code{print.n} or by \code{options(survfit.print.n)}: \code{"none"}
prints \code{NA}, \code{"records"} prints the number of records,
\code{"start"} prints the number at the first time point and
\code{"max"} prints the maximum number at risk. The initial default is
\code{"start"}.

}
\references{
Miller, Rupert G., Jr. (1981) 
\emph{Survival Analysis, }
Wiley, New York, p 71.
}
\seealso{
\code{\link{summary.survfit}}
}
\examples{
##effect of print.n
data(heart)
a<-coxph(Surv(start,stop,event)~age+strata(transplant),data=heart)
b<-survfit(a)
print(b,print.n="none")
print(b,print.n="records")
print(b,print.n="start")
print(b,print.n="max")
}

\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
