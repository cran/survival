\name{date.ddmmmyy}
\alias{date.ddmmmyy}
\title{Format a Julian date}
\description{
  Given a vector of Julian dates, this returns them in the form
  ``10Nov89'', ``28Jul54'', etc.
}
\usage{
date.ddmmmyy(sdate)
}
\arguments{
  \item{sdate}{A vector of Julian dates, e.g., as returned by
    \code{mdy.date()}.
}
\value{
  A vector of character strings containing the formatted dates.
}
\seealso{
  \code{\link{mdy.date}},
  \code{\link{date.mdy}}
}
\examples{
date.ddmmmyy(1:10)
}
\keyword{chron}
