\name{survsum}
\alias{survsum}
\title{
Kaplan-Meier Survival Curve and Percentages at Selected Times
}
\description{
Calculates Kaplan-Meier survival percentages, standard error, and number 
at risk at specified times for defined groups.  Test for difference between
survival curves.
}
\usage{
survsum(formula, data=sys.parent(), sptms=NULL, xlim, 
tlines=T, log=F, xscale=1,yscale=100, mark.time=F, mark=3, 
cex=1, xlab="Time", ylab="Survival (\%)", lgd="bl",
ttl="K-M Survival",\dots)
}
\arguments{
\item{formula}{
a formula expression as for other survival models, of the  form
\code{Surv(time,  status) ~ predictors}.  Same formula expression as used 
in \code{survfit}. See help file for survfit.  Maximum of 6 groups.
}
\item{data}{
a data.frame in which to interpret the variables named in the formula.
}
\item{sptms}{
a specified vector of positive times at which to compute the survival 
percentages, standard errors, and numbers at risk.  A maximum of four 
times can be specified.
}
\item{xlim}{
a vector of the form: \code{c(x1,x2)}.  The approximate minimum and maximum values 
to  be  put  on  x-axis.  Default sets x1=0 and x2=maximum time value.
}
\item{tlines}{
a logical value indicating whether vertical lines and labels should be drawn
on the plot at the specified times.
}
\item{log}{
logical value: should the y axis be on a log scale?
}
\item{xscale}{
a scalar to be used to divide the x axis.  A value of 365, for instance, would
be used to convert from days to years.
}
\item{yscale}{
a scalar to be used to multiply the y axis.  The default value of 100 is used
to get a percent scale.  \code{yscale=1} would set the y axis from 0 to 1.
}
\item{mark.time}{
controls the labeling of the curves.  If set to True then curves are 
marked at each censoring time.  If \code{mark.time} is a numeric vector, then curves
are marked at these specified time points.
}
\item{mark}{
vector of mark parameters, which will be used to  label the  curves.  The
\code{lines} help file contains examples of the possible marks.  The 
vector is reused cyclically if it  is shorter than the number of curves.
}
\item{cex}{
parameter available to change the size of \code{mark}.   Not a vector; all marks
have the same size.
}
\item{xlab}{
character string label for the x axis.
}
\item{ylab}{
character string label for the y axis.
}
\item{lgd}{
legend placement.  \code{"tr"}=top right corner of the plot,
\code{"under"}=under the plot, \code{"n"} omits the legend.  The default is:
\code{"bl"}=bottom left corner of the plot.
}
\item{ttl}{
title to be printed in the center top.
}
\item{...}{
In addition, the high-level graphics arguments described under par and the
arguments to title may be supplied to this function.
}}
\value{
A list of class \code{"htest"}, containing the following components:



\item{no.pts}{
the total number of data points in each group.
}
\item{no.events}{
the total number of events in each group.
}
\item{chisq}{
the chisquare statistic for the test of a difference between survival curves.
}
\item{p}{
the p-value for the above test.
}
\item{t1}{
a matrix containing the survival percentages, standard errors, and numbers at
risk for all groups at time t1.
}
\item{t2,t3,t4}{
see above.
}}
\section{Side Effects}{
A plot with multiple survival curves is drawn (one for each group). This plot
includes: overall group statistics and group statistics at each specified 
time point and a test for a difference between survival curves.
}
\details{
The total number of points and events is reported for each group.  For each 
specified time point, group survival percentages (followed by standard error
and number left at risk) are computed.


The test for a difference between survival curves uses the chisquare
statistic from the \code{survdiff} function with \code{rho=0}.  This is the
log-rank test.
}
\author{
Mark Dietrich, Mayo Clinic Section of Medical Research Statistics
summer student 1992.
}
\seealso{
\code{\link{Surv}}, \code{\link{survdiff}}, \code{\link{survexp}}, \code{\link{survfit}}
\examples{
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
