\name{summary.survfit}
\alias{summary.survfit}
\alias{print.summary.survfit}
\title{
Summary of a Survival Curve
}
\description{
Returns a list containing the survival curve, confidence limits for the
curve, and other information.
}
\usage{
summary.survfit(object, times, censored=F, scale=1, ...)
}
\arguments{
\item{object}{
output from a call to \code{survfit}.
}
\item{times}{
vector of times;
the returned matrix will contain 1 row for each time.
This must be in increasing order and missing values are not allowed.
If \code{censored=T}, the default \code{times} vector contains all the unique times in
\code{fit},
otherwise the default \code{times} vector uses only the event (death) times.
}
\item{censored}{
logical flag: should the censoring times be included in the output?
This is ignored if the \code{times} argument is present.
}
\item{scale}{
rescale the survival time, e.g., if the input data to \code{survfit} were in
days, \code{scale=365.25} would scale the output to years.
}\item{extend}{Estimate survival beyond the last observed time?}\item{...}{other unused arguments}}
\value{
a list with the following components
}
\item{time}{
the timepoint on the curve.
}
\item{surv}{
the value of the survival curve at time t+0.
}
\item{n.risk}{
the number of subjects at risk at time t-0
(but see the comments on weights in the \code{survfit} help file).
}
\item{n.event}{
if the \code{times} argument is missing, then this column is the number of
events that occurred at time t.
Otherwise, it is the cumulative number of events that have occurred
since the last time listed until time t+0.
}
\item{std.err}{
the standard error of the survival value.
}
\item{conf.int}{
level of confidence for the confidence intervals of survival.
}
\item{lower}{
lower confidence limits for the curve.
}
\item{upper}{
upper confidence limits for the curve.
}
\item{strata}{
indicates stratification of curve estimation.  If \code{strata} is not \code{NULL},
there are multiple curves in the result and the \code{surv}, \code{time}, \code{n.risk}, etc. 
vectors will contain multiple curves, pasted end to end. 
The levels of \code{strata} (a factor) are the labels for the curves.
}
\item{call}{
the statement used to create the \code{fit} object.
}
\item{na.action}{
passed through from \code{fit}, if present.
}}
\seealso{
\code{\link{survfit}}, \code{\link{print.summary.survfit}}.}
\examples{
data(ovarian)
summary( survfit( Surv(futime, fustat),data=ovarian))
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
