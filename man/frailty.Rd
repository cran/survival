\name{frailty}
\alias{frailty}
\alias{frailty.gamma}
\alias{frailty.gaussian}
\alias{frailty.t}
\alias{frailty.brent}
\alias{frailty.controlaic}
\alias{frailty.controldf}
\alias{frailty.controlgam}
\alias{frailty.controlgauss}
\alias{frailty.gammacon}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{(Approximate) Frailty models}
\usage{
frailty(x, distribution="gamma", ...)
frailty.gamma(x, sparse = (nclass > 5), theta, df, eps = 1e-05, method = c("em","aic", "df", "fixed"), ...) 
frailty.gaussian(x, sparse = (nclass > 5), theta, df, method = c("reml","aic", "df", "fixed"), ...)
frailty.t(x, sparse = (nclass > 5), theta, df, eps = 1e-05, tdf = 5,method = c("aic", "df", "fixed"), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ group indicator}
 \item{distribution}{ frailty distribution }
 \item{\dots}{Arguments for specific distribution, including (but not
   limited to) }
 \item{sparse}{Use sparse Newton-Raphson algorithm}
 \item{df}{Approximate degrees of freedom}
 \item{theta}{Penalty}
 \item{eps}{Accuracy of \code{df}}
 \item{method}{maximisation algorithm}
 \item{tdf}{df of t-distribution}
}

\description{ 
When included in a \link{coxph} or \link{survreg}, fits
by penalised likelihood a random effects (frailty) model. \code{frailty} is generic, with methods for t, Gaussian and Gamma distributions.
}

\details{
The penalised likelihood method is equivalent to maximum (partial)
likelihood for the gamma frailty but not for the others.

The sparse algorithm uses the diagonal of the information matrix for
the random effects, which saves a lot of space. 

The frailty distributions are really the log-t and lognormal: t and
Gaussian are random effects on the scale of the linear predictor.
}
\value{
  An object of class \code{coxph.penalty} containing a factor with attributes specifying the control functions.
}

\seealso{\code{\link{coxph}},\code{\link{survreg}},\code{\link{ridge}},\code{\link{pspline}}  }
\references{Therneau TM, Grambsch PM, Pankratz VS (2003) "Penalized
  survival models and frailty" Journal of Computational and Graphical
  Statistics 12, 1: 156-175}
\examples{
kfit <- coxph(Surv(time, status)~ age + sex + disease + frailty(id), kidney)
kfit0 <- coxph(Surv(time, status)~ age + sex + disease, kidney)
kfitm1 <- coxph(Surv(time,status) ~ age + sex + disease + 
		frailty(id, dist='gauss'), kidney)
coxph(Surv(time, status) ~ age + sex + frailty(id, dist='gauss', method='aic',caic=TRUE), kidney)
# uncorrected aic
coxph(Surv(time, status) ~ age + sex + frailty(id, method='aic', caic=FALSE), kidney)

rfit2a <- survreg(Surv(time, status) ~ rx +
		  frailty.gaussian(litter, df=13, sparse=FALSE), rats )
rfit2b <- survreg(Surv(time, status) ~ rx +
		  frailty.gaussian(litter, df=13, sparse=TRUE), rats )
rfit2a
rfit2b
}
\keyword{survival}%-- one or more ...
