\name{ratetable}
\alias{ratetable}
\alias{[.ratetable}
\alias{[.ratetable2}
\alias{print.ratetable}
\alias{is.na.ratetable}
\alias{summary.ratetable}
\title{Ratetable reference in formula}
\description{
This function matches variable names in data to those in a ratetable for
\code{\link{survexp}}
}
\usage{
ratetable(...)
}

\arguments{
  \item{\dots}{tags matching dimensions of the ratetable and variables
    in the data frame (see example)}
}

\value{
  A data frame
}

\seealso{\code{\link{survexp}},\code{\link{survexp.us}},\code{\link{is.ratetable}}}

\examples{

fit <- survfit(Surv(time, status) ~ sex, pbc,subset=1:312)
if (require(date)){
tdata <- data.frame(agedays=pbc$age*365.24, sex=pbc$sex +1,
		    year= rep(mdy.date(1,1,1976), nrow(pbc)))
tdata<-tdata[1:312,] ## only the randomised people, with no missing data

efit <- survexp(~ sex+ratetable(sex=sex,age=agedays,year=year), data=tdata, ratetable=survexp.us, times=(0:24)*182)
}

}
\keyword{survival}%-- one or more ...
