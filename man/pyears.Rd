\name{pyears}
\alias{pyears}
\title{
Person Years
}
\description{
This function computes the person-years of follow-up time contributed by a
cohort of subjects, stratified into subgroups.
It also computes the number of subjects who contribute to each cell of the
output table, and optionally the number of events and/or expected number of
events in each cell.
}
\usage{
pyears(formula, data, weights, subset, na.action, ratetable=survexp.us,
scale=365.25, expect=c('event', 'pyears'), model=F, x=F, y=F)
}
\arguments{
\item{formula}{
a formula object.  The response variable will be a vector of follow-up times
for each subject, or a \code{Surv} object containing the follow-up time and an
event indicator.
The predictors consist of optional grouping variables
separated by + operators (exactly as in \code{survfit}), time-dependent grouping
variables such as age (specified with \code{tcut}), and optionally a
\code{ratetable()} term.  This latter matches each subject to his/her expected
cohort.
}
\item{data}{
a data frame in which to interpret the variables named in
the \code{formula}, or in the \code{subset} and the \code{weights} argument.
}
\item{weights}{
case weights.
}
\item{subset}{
expression saying that only a subset of the rows of the data
should be used in the fit.
}
\item{na.action}{
a missing-data filter function, applied to the model.frame, after any
\code{subset} argument has been used.  Default is \code{options()$na.action}.
}
\item{ratetable}{
a table of event rates, such as \code{survexp.uswhite}.
}
\item{scale}{
a scaling for the results.  As most rate tables are in units/day, the
default value of 365.25 causes the output to be reported in years.
}
\item{expect}{
should the output table include the expected number of events, or the
expected number of person-years of observation.  This is only valid with
a rate table.
}
\item{model,}{
If any of these is true, then the
model frame, the model matrix, and/or the vector of response times will be
returned as components of the final result.
}}
\value{
a list with components

\item{pyears}{
an array containing the person-years of exposure. (Or other units, depending
on the rate table and the scale).
}
\item{n}{
an array containing the number of subjects who contribute time to each cell
of the \code{pyears} array.
}
\item{event}{
an array containing the observed number of events.  This will be present only
if the response variable is a \code{Surv} object.
}
\item{expected}{
an array containing the expected number of events (or person years).
This will be present only
if there was a ratetable term.
}
\item{offtable}{
the number of person-years of exposure in the cohort that was not part of
any cell in the \code{pyears} array.  This is often useful as an error check; if
there is a mismatch of units between two variables, nearly all the person
years may be off table.
}
\item{summary}{
a summary of the rate-table matching.  This is also useful as an error
check.
}
\item{call}{
an image of the call to the function.
}
\item{na.action}{
the \code{na.action} attribute contributed by an \code{na.action} routine, if any.
}}
\details{
Because \code{pyears} may have several time variables, it is necessary that all
of them be in the same units.  For instance, in the call
 py <- pyears(futime ~ rx + ratetable(age=age, sex=sex, year=entry.dt))
with a \code{ratetable} whose natural unit is days, it is important that \code{futime},
\code{age} and \code{entry.dt} all be in days.  Given the wide range of possible inputs,
it is difficult for the routine to do sanity checks of this aspect.


A special function \code{tcut} is needed to specify time-dependent cutpoints.
For instance, assume that age is in years, and that the desired final
arrays have as one of their margins the age groups 0-2, 2-10, 10-25, and 25+.
A subject who enters the study at age 4 and remains under observation for
10 years will contribute follow-up time to both the 2-10 and 10-25
subsets.  If \code{cut(age, c(0,2,10,25,100))} were used in the formula, the
subject would be classified according to his starting age only.
The \code{tcut} function has the same arguments as \code{cut},
but produces a different output object which allows the \code{pyears} function
to correctly track the subject.


The results of \code{pyears} are normally used as input to further calculations.
The example below is from a study of hip fracture rates from 1930 - 1990
in Rochester, Minnesota.  Survival post hip fracture has increased over
that time, but so has the survival of elderly subjects in the population at
large.  A model of relative survival helps to clarify what has happened:
Poisson regression is used, but replacing exposure time with expected
exposure (for an age and sex matched control).
Death rates change with age, of course, so the result is carved into
1 year increments of time.  Males and females were done separately.
}
\seealso{
\code{\link{ratetable}}, \code{\link{survexp}}, \code{\link{Surv}}
\examples{\dontrun{
library(date)
data(ratetables)
# 
# Simple case: a single male subject, born 6/6/36 and entered on study 6/6/55.
#

temp1 <- mdy.date(6,6,36)
temp2 <- mdy.date(6,6,55)# Now compare the results from person-years
#
temp.age <- tcut(temp2-temp1, floor(c(-1, (18:31 * 365.24))),
        labels=c('0-18', paste(18:30, 19:31, sep='-')))
temp.yr  <- tcut(temp2, mdy.date(1,1,1954:1965), labels=1954:1964)
temp.time <- 3700   #total days of fu
py1 <- pyears(temp.time ~ temp.age + temp.yr, scale=1) #output in days
}}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
