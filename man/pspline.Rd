\name{pspline}
\alias{pspline}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Penalised smoothing splines}
\usage{
pspline(x, df=4, theta, nterm=2.5 * df, degree=3, eps=0.1, method, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{predictor }
 \item{df}{ approximate degrees of freedom. \code{df=0} means use AIC }
 \item{theta}{roughness penalty }
 \item{nterm}{ number of splines in the basis }
 \item{degree}{ degree of splines }
 \item{eps}{accuracy for \code{df} }
 \item{method}{Method for automatic choice of \code{theta}}
 \item{\dots}{I don't know what this does }
}
\description{
 Specifies a penalised spline basis for the predictor.  This is done by fitting a comparatively small set of splines and penalising the integrated second derivative. Results are similar to smoothing splines with a knot at each data point but computationally simpler.
}

\value{
  Object of class \code{coxph.penalty} containing the spline basis with 
attributes specifying control functions.
}

\seealso{\code{\link{coxph}},\code{\link{survreg}},\code{\link{ridge}},\code{\link{frailty}}   }

\examples{
lfit6 <- survreg(Surv(time, status)~pspline(age, df=2), cancer)
plot(cancer$age, predict(lfit6), xlab='Age', ylab="Spline prediction")
title("Cancer Data")
fit0 <- coxph(Surv(time, status) ~ ph.ecog + age, cancer)
fit1 <- coxph(Surv(time, status) ~ ph.ecog + pspline(age,3), cancer)
fit3 <- coxph(Surv(time, status) ~ ph.ecog + pspline(age,8), cancer)
fit0
fit1
fit3
}
\keyword{ survival}%-- one or more ...
