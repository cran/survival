\name{as.date}
\alias{as.date}
\alias{as.data.POSIXt}
\title{Coerce Data to Dates}
\description{
  Converts any of the following character forms to a Julian date:
  8/31/56, 8-31-1956, 31 8 56, 083156, 31Aug56, or August 31 1956.
}
\usage{
as.date(x, order = "mdy", ...)
}
\arguments{
  \item{x}{input data vector.}
  \item{order}{if x is character, defines the order in which the terms
    are assumed to appear in a xx/xx/xx date.  The default is
    month/day/year; any permutation of mdy is legal.}
  \item{...}{if x is character, then any other arguments from
    \code{mdy.date()} can be used as well.}
}
\value{
  For each date, the number of days between it and January 1, 1960.
  The date will be missing if the string is not interpretable.
}
\details{
  If \code{x} is numeric, then \code{floor(x)} is returned, e.g.,
  \code{as.date(35)} is the same as \code{as.date(35.2)} and gives
  \code{Feb 2, 1960}.  If \code{x} is character, the program attempts
  to parse it.
}
\seealso{
  \code{\link{mdy.date}},
  \code{\link{date.mmddyy}}, 
  \code{\link{date.ddmmmyy}}
}
\examples{
as.date(c("1jan1960", "2jan1960", "31mar1960", "30jul1960"))
}
\keyword{chron}
