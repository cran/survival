\name{is.ratetable}
\alias{is.ratetable}
\alias{Math.ratetable}
\alias{Ops.ratetable}
\title{
Verify that an Object is of Class `ratetable'
}
\description{
The function verifies not only the \code{class} attribute, but the
structure of the object.
}
\usage{
is.ratetable(x, verbose=F)
}
\arguments{
\item{x}{
the object to be verified.
}
\item{verbose}{
if \code{TRUE} and the object is not a ratetable, then return a character string
describing the way(s) in which \code{x} fails to be a proper ratetable object.
}}
\value{
returns \code{TRUE} if \code{x} is a ratetable, and \code{FALSE} or a description if it is not.
}
\details{
Rate tables are used by the \code{pyears} and \code{survexp} functions, and normally
contain death rates for some population, categorized by age, sex, or other
variables.  They have a fairly rigid structure, and the \code{verbose} option
can help in creating a new rate table.
}
\seealso{
\code{\link{pyears}}, \code{\link{survexp}}
}
\examples{\dontrun{
library(ratetables)
is.ratetable(survexp.us)  ##Yes
data(cancer)
is.ratetable(cancer)  ##No
}}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
