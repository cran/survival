\name{survreg.distributions}
\alias{survreg.distributions}
\title{Parametric Survival Distributions}
\usage{
survreg.distributions
}
\description{ List of distributions for accelerated failure models. These are
  location-scale families for some transformation of time. The entry
  describes  the cdf \eqn{F} and density \eqn{f} of a canonical member of
  the family.
}

\format{
 
  There are three basic formats; only the first two are used in the
  built-in distributions
  \tabular{ll}{
    name:\tab name of distribution\cr
    variance:\tab Variance \cr
    init(x,weights,...):\tab Function returning an initial\cr mean and
    \tab variance\cr
    deviance(y,scale,parms):\tab Function returning the deviance\cr
    density(x,parms):\tab Function returning \eqn{F},\cr
    \tab \eqn{1-F},\eqn{f},\eqn{f'/f},\eqn{f''/f}\cr
    quantile(p,parms):\tab Quantile function\cr
    scale:\tab Optional fixed value for scale parameter\cr
  }
  and for transformations of the time variable
  \tabular{ll}{
    name:\tab name of distribution\cr
    dist:\tab name of transformed distribution\cr
    trans:\tab transformation (eg log)\cr
    dtrans:\tab derivative of transformation\cr
    itrans:\tab inverse of transformation\cr
    scale:\tab Optional fixed value for scale parameter\cr
  }
  For transformations of user-defined families use
  \tabular{ll}{
    name:\tab name of distribution\cr
    dist:\tab transformed distribution in first format\cr
    trans:\tab transformation (eg log)\cr
    dtrans:\tab derivative of transformation\cr
    itrans:\tab inverse of transformation\cr
    scale:\tab Optional fixed value for scale parameter\cr
  }
  
}

\details{
There are four basic distributions:\code{extreme}, \code{gaussian},
\code{logistic} and \code{t}. The last three
are parametrised in the same way as the distributions already present in
\R. The extreme value cdf is
\deqn{F=1-e^{-e^t}.}

When the logarithm of survival time has one of the first three distributions
we obtain respectively \code{weibull}, \code{lognormal}, and
\code{loglogistic}. The Weibull distribution is not parameterised the
same way as in \code{\link{rweibull}}.

The other predefined distributions are defined in terms of these. The
\code{exponential} and \code{rayleigh} distributions are Weibull
distributions with fixed \code{scale} of 1 and 0.5 respectively, and
\code{loggaussian} is a synonym for \code{lognormal}.

Parts of the built-in distributions are hardcoded in C, so the elements
of \code{survreg.distributions} in the first format above must not be
changed and new ones must not be added.  The examples show how to
specify user-defined distributions to \code{survreg}.


}

\seealso{\code{\link{survreg}},
  \code{\link{pnorm}},\code{\link{plogis}}, \code{\link{pt}}
}
\examples{
## not a good fit, but a useful example
survreg(Surv(futime,fustat)~ecog.ps+rx,data=ovarian,dist='extreme')
## 
my.extreme<-survreg.distributions$extreme
my.extreme$name<-"Xtreme"
survreg(Surv(futime,fustat)~ecog.ps+rx,data=ovarian,dist=my.extreme)

## time transformation
survreg(Surv(futime,fustat)~ecog.ps+rx,data=ovarian,dist='weibull',scale=1)
my.weibull<-survreg.distributions$weibull
my.weibull$dist<-my.extreme
survreg(Surv(futime,fustat)~ecog.ps+rx,data=ovarian,dist=my.weibull,scale=1)

## change the transformation to work in years
## intercept changes by log(365), other coefficients stay the same
my.weibull$trans<-function(y) log(y/365)
my.weibull$itrans<-function(y) exp(365*y)
survreg(Surv(futime,fustat)~ecog.ps+rx,data=ovarian,dist=my.weibull,scale=1)

## Weibull parametrisation
y<-rweibull(1000, shape=2, scale=5)
survreg(Surv(y)~1, dist="weibull")
## survreg reports scale=1/2, intercept=log(5)
}
\keyword{survival}
