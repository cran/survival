\name{date.object}
\alias{is.date}
\alias{Math.date}
\alias{Ops.date}
\alias{Summary.date}
\alias{[.date}
\alias{[[.date}
\alias{as.character.date}
\alias{as.data.frame.date}
\alias{as.vector.date}
\alias{is.na.date}
\alias{plot.date}
\alias{print.date}
\alias{summary.date}
\title{Date Objects}
\description{
  Creates and tests for objects of class \code{"date"}.
}
\usage{
as.date(x)
is.date(x)
}
\arguments{
  \item{x}{any R object.}
}
\details{
  Dates are stored as the number of days since 1/1/1960, and are kept
  in integer format.  (This is the same baseline value as is used by
  SAS).  The numerical methods for dates treat \code{date - date} as a
  numeric, and \code{date +- numeric} as a date.

  \code{is.date} returns \code{TRUE} if \code{x} has class
  \code{"date"}, and \code{FALSE} otherwise.  Its behavior is unaffected
  by any attributes of \code{x}; for example, \code{x} could be a date
  array (in contrast to the behavior of \code{is.vector}).

  \code{as.date} returns \code{x} if \code{x} is a simple object of
  class \code{"date"}, and otherwise a date vector of the same length as
  \code{x} and with data resulting from coercing the elements of
  \code{x} to class \code{"date"}.  See the manual page for
  \code{as.date()} for details.

  Logical operations as well as the numeric functions \code{exp()},
  \code{log()}, and so on are invalid.

  Other methods exist for missing value, \code{as.character()},
  printing, and summarizing. 
}
\seealso{
  \code{\link{date.mdy}},
  \code{\link{mdy.date}},
  \code{\link{date.ddmmmyy}}, 
  \code{\link{as.date}}.
}
\keyword{chron}
