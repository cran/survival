\name{bladder}
\alias{bladder}
\alias{bladder2}
\title{Bladder Cancer Recurrences}
\usage{data(bladder)}
\description{Data on recurrences of bladder cancer, used by many people
  to demonstrate methodology for recurrent event modelling.}
  
\format{
 bladder
    \tabular{ll}{
	id:\tab Patient id\cr
	rx:\tab Treatment 1=placebo  2=thiotepa\cr
	number:\tab Initial number of tumours (8=8 or more)\cr
	size:\tab size (cm) of largest initial tumour\cr
	stop:\tab recurrence or censoring time\cr
	enum:\tab which recurrence (up to 4)\cr
     }
 bladder2 
    \tabular{ll}{
	id:\tab Patient id\cr
	rx:\tab Treatment 1=placebo  2=thiotepa\cr
	number:\tab Initial number of tumours (8=8 or more)\cr
	size:\tab size (cm) of largest initial tumour\cr
	start:\tab start of interval (0 or previous recurrence time)\cr
	stop:\tab recurrence or censoring time\cr
	enum:\tab which recurrence (up to 4)\cr
     }

 }
\source{Wei, Lin, Weisfeld, JASA 1989}
\keyword{datasets}
\keyword{survival}
