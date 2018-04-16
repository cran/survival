\name{ratetable}
\alias{ratetable}
\alias{[.ratetable}
\alias{print.ratetable}
\alias{summary.ratetable}
\title{
 Rate table structure
}
\description{
Description of the rate tables used by expected survival routines.
}
\details{
  A rate table contains event rates per unit time for some particular
  endpoint.  Death rates are the most common use, the \code{survexp.us}
  table, for instance, contains death rates for the United States by
  year of age, sex, and calendar year.

  A rate table is structured as a multi-way array with the following
  attributes:

  \describe{
    \item{dim}{the dimensions of the array}

    \item{dimnames}{a named list of dimnames. The names are used to
    match user data to the dimensions, e.g., see the \code{rmap}
    argument in the \code{pyears} example.  If a dimension is
    categorical, such as \code{sex} in \code{survexp.us}, then the
    dimname itself is matched against user's data values.  The matching
    ignored case and allows abbreviations, e.g., "M", "Male", and "male"
    all successfully match the \code{survexp.us} dimname of
    \code{sex=c("male", "female")}.}

    \item{type}{a vector giving the type of each dimension, which will
      be 1= categorical, 2= continuous, 3= date, 4= US calendar year.
      If \code{type} is 3 or 4, then the corresponding cutpoints must be
      one of the calendar date types: Date, POSIXt, date, or chron.
      This allows the code to properly match user data to the ratetable.}

    \item{cutpoints}{a list with one elment per dimension.  If
      \code{type=1} then the corresponding list element should be NULL,
      otherwise it should be a vector of length \code{dim[i]} containing
      the starting point of the interval to which the corresponding
      row/col of the array applies.  Cutpoints must be in the same units
      as the underlying table, e.g., the \code{survexp.us} table
      contains death rates per day, so the \code{age} cutpoint vector contains
      age in days while \code{year} contains a vector of Dates.
      Cutpoints do not need to be evenly
      spaced: the \code{survexp.us} table originally had age divided up
      as 0-1 days, 1-7 days, 7-28 days, 28 days - 1 year, 2, 3, \ldots
      119 years.  (Changes in the source of the tables made it difficult
      to continue splitting out the first year.)}

    \item{summary}{an optional summarization function.  If present, it
    will be called with a numeric matrix that has one column per
    dimension and one row per observation; and the function returns a
    character string giving a summary of the data.}

    \item{dimid}{optional attribute containing the names of the
    dimnames. In this case the dimnames list itself does not have
    names.}
}
}
\seealso{\code{\link{survexp}}, \code{\link{pyears}},
  \code{\link{survexp.us}}}
\keyword{survival}
