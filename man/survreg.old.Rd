\name{survreg.old}
\alias{survreg.old}
\title{
Old (survival4) Regression for a parametric survival model
}
\usage{
survreg.old(formula, data=sys.frame(sys.parent()), \dots, link=c("log",
"identity"),dist=c("extreme", "logistic", "gaussian",
"exponential","rayleigh","weibull"), fixed=list())
}
\description{
  This routine is a backwards-compatible interface to the improved \link{survreg} function, which is better.
}
\arguments{
\item{formula}{
a formula expression as for other regression models.
See the documentation for \code{lm} and \code{formula} for details.
}
\item{data}{
optional data frame in which to interpret the variables occuring in the
formula.
}
\item{\dots}{other arguments to \code{\link{survreg}}}
\item{link}{
transformation to be used on the y variable.
}
\item{dist}{
assumed distribution for the transformed y variable.
}
\item{fixed}{
a list of fixed parameters, most often just the scale.
}
}
\value{
an object of class \code{survreg} is returned, which inherits from class \code{glm}.
}

\examples{
survreg.old(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='extreme',
		link='log', fixed=list(scale=1))   #Fit an exponential
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
