\name{stanford2}
\alias{stanford2}
\title{More Stanford Heart Transplant data}
\usage{data(stanford2)}

\format{
    \tabular{ll}{
        id:  \tab ID number\cr
	time:\tab survival or censoring time\cr
	status:\tab censoring status\cr
	age: \tab in years\cr
        t5: \tab T5 mismatch score\cr
     }
 }
\description{This contains the Stanford Heart Transplant data in
a different format. The main data set is in \code{\link{heart}}}

\seealso{\code{\link{predict.survreg}},\code{\link{heart}}}
\source{Escobar & Meeker, Biometrics 1992 p519}
\keyword{datasets}
\keyword{survival}