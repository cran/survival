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
 
  There are two basic formats
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
  
}

\details{
There are four basic distributions:\code{extreme}, \code{gaussian},
\code{logistic} and \code{t}. The last three
are parametrised in the same way as the distributions already present in
\R. The extreme value cdf is
\deqn{F=1-e^{-e^t}.}

When the logarithm of survival time has one of the first three distributions
we obtain respectively \code{lognormal},\code{loglogistic} and
\code{weibull}.

The other predefined distributions are defined in terms of these. The
\code{exponential} and \code{rayleigh} distributions are Weibull
distributions with fixed \code{scale} of 1 and 0.5 respectively, and
\code{loggaussian} is a synonym for \code{lognormal}.

} \seealso{\code{\link{survreg}},
\code{\link{pnorm}},\code{\link{plogis}}, \code{\link{pt}}}}
\keyword{survival}
