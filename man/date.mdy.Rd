\name{date.mdy}
\alias{date.mdy}
\title{Convert from Julian Dates to Month, Day, and Year}
\description{
  Convert a vector of Julian dates to a list of vectors with the
  corresponding values of month, day and year, and optionally weekday.
}
\usage{
date.mdy(sdate, weekday = FALSE)
}
\arguments{
  \item{sdate}{a Julian date value, as returned by \code{mdy.date()},
    number of days since 1/1/1960.}
  \item{weekday}{if \code{TRUE}, then the returned list also will
    contain the day of the week (Sunday=1, Saturday=7).}
}
\value{
  A list with components \code{month}, \code{day}, and \code{year}.
}
\examples{
day <- 7
temp <- date.mdy(mdy.date(month = 7, day = day, year = 1960))
## Check for illegal dates, such as 29 Feb in a non leap year
if (temp$day != day) {
  cat("Some illegal dates\n")
} else {
  cat("All days are legal\n")
}
}
\references{
  Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery,
  B. P. (1992).
  \emph{Numerical Recipes: The Art of Scientific Computing (Second
    Edition)}.
  Cambridge University Press.
}
\keyword{chron}
