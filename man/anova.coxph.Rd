\name{anova.coxph}
\alias{anova.coxph}
\alias{anova.coxphlist}
\title{Analysis of Deviance for Cox model.}
\usage{
\method{anova}{coxph}(object, \dots,  test = NULL)
}
\description{
  Compute an analysis of deviance table for one Cox model fit.
}
\arguments{
  \item{object}{An object of class \code{coxph}}
  \item{\dots}{Further \code{coxph} objects}
  \item{test}{a character string, (partially) matching one of \code{"Chisq"},
    \code{"F"} or \code{"Cp"}. See \code{\link{stat.anova}}.}
}
\details{
  Specifying a single object gives a sequential analysis of deviance
  table for that fit.  That is, the reductions in the residual deviance
  as each term of the formula is added in turn are given in as
  the rows of a table, plus the residual deviances themselves.


  If more than one object is specified, the table has a row for the
  residual degrees of freedom and deviance for each model. For all
  but the first model, the change in degrees of freedom and deviance
  is also given. (This only make statistical sense if the models are
  nested.)  It is conventional to list the models from smallest to
  largest, but this is up to the user.
  
  The table will optionally contain test statistics (and P values)
  comparing the reduction in deviance for the row to the residuals.
    Mallows' \eqn{C_p}{Cp} statistic is the residual
    deviance plus twice the  residual degrees of freedom, which is closely
  related to AIC. Don't use \code{test=F}, it doesn't make sense.
}
\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\section{Warning}{
  The comparison between two or more models by \code{anova} or
   will only be valid if they
  are fitted to the same dataset. This may be a problem if there are
  missing values.}
\seealso{
  \code{\link{coxph}}, \code{\link{anova}}.
}
\examples{
data(ovarian)
fit <- coxph(Surv(futime, fustat) ~ resid.ds *rx + ecog.ps, data = ovarian) 
anova(fit)
anova(fit,test="Chisq")
fit2 <- coxph(Surv(futime, fustat) ~ resid.ds +rx + ecog.ps, data
=ovarian)
anova(fit2,fit)
 }
\keyword{models}
\keyword{regression}
\keyword{survival}

