\name{mdy.date}
\alias{mdy.date}
\title{Convert to Julian Dates}
\description{
  Given a month, day, and year, returns the number of days since
  January 1, 1960.
}
\usage{
mdy.date(month, day, year, nineteen = TRUE, fillday = FALSE,
         fillmonth = FALSE)
}
\arguments{
  \item{month}{vector of months.}
  \item{day}{vector of days.}
  \item{year}{vector of years.}
  \item{nineteen}{if \code{TRUE}, year values between 0 and 99 are
    assumed to be in the 20th century A.D.; otherwise, if \code{FALSE},
    they are assumed to be in the 1st century A.D.}
  \item{fillday}{if \code{TRUE}, then missing days are replaced with
    15.}
  \item{fillmonth}{if \code{TRUE}, then a missing month causes the
    month and day to be set to 7/1.}
}
\value{
  a vector of Julian dates.
}
\details{
  The date functions are particularly useful in computing time spans,
  such as number of days on test, and similar functions can be found in
  other statistical packages.  The baseline date of Jan 1, 1960 is, of
  course, completely arbitrary (it is the same one used by SAS).

  The \code{fillday} and \code{fillmonth} options are perhaps useful
  only to the author and a very few others: we sometimes deal with
  patients whose birth date was in the 1800's, and only the month or
  even only the year is known.  When the interval is greater than 80
  years, a filler seems defensible.
}
\references{
  Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery,
  B. P. (1992).
  \emph{Numerical Recipes: The Art of Scientific Computing (Second
    Edition)}.
  Cambridge University Press.

}
\seealso{
  \code{\link{date.mmddyy}}, 
  \code{\link{date.ddmmmyy}},
  \code{\link{date.mmddyyyy}}
}
\examples{
mdy.date(3, 10, 53)
xzt <-1:10
xzy <- as.date(xzt)
test <- data.frame(x = xzt, date = xzy)
summary(test)
}
\keyword{chron}
