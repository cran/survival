\name{date.mmddyy}
\alias{date.mmddyy}
\title{Format a Julian date}
\description{
  Given a vector of Julian dates, this returns them in the form
  ``10/11/89'', ``28/7/54'', etc.
}
\usage{
date.mmddyy(sdate, sep = "/")
}
\arguments{
  \item{sdate}{A vector of Julian dates, e.g., as returned by
    \code{mdy.date()}.}
  \item{sep}{Character string used to separate the month, day, and
    year portions of the returned string.}
}
\value{
  A vector of character strings containing the formatted dates.
}
\seealso{
  \code{\link{date.mdy}},
  \code{\link{mdy.date}}, 
  \code{\link{date.ddmmmyy}}
}
\examples{
date.mmddyy(as.date(10))
}
\keyword{chron}
