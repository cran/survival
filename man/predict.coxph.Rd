\name{predict.coxph}
\alias{predict.coxph}
\alias{predict.coxph.penal}

\title{Predictions for Cox model}
\description{
Compute fitted values and regression terms for a model fitted by \code{\link{coxph}}
}
\usage{
\method{predict}{coxph}(object, newdata, type=c("lp", "risk", "expected", "terms"), se.fit=FALSE, terms=names(object$assign), collapse, safe=FALSE,...)
}

\arguments{
  \item{object}{A \code{coxph} object}
  \item{newdata}{Optional new data to predict at}
  \item{type}{Type of prediction wanted}
  \item{se.fit}{Return standard errors as well?}
  \item{terms}{If \code{type="terms"}, which terms to return.}
  \item{collapse}{identifier for which rows correspond to different individuals}
  \item{safe}{Use a more intensive `safe' prediction method}
  \item{...}{other unused arguments}
}

\value{
  A vector or matrix of fitted values. If \code{se.fit=TRUE} a list
  whose first component is the fitted values and second is the standard errors.
}

\seealso{\code{\link{predict}},\code{\link{coxph}},\code{\link{termplot}} }

\examples{
fit<-coxph(Surv(time,status)~x,data=aml)
predict(fit,type="lp")
predict(fit,type="risk")
predict(fit,type="expected")
predict(fit,type="terms")
predict(fit,type="lp",se.fit=TRUE)
predict(fit,type="risk",se.fit=TRUE)
predict(fit,type="expected",se.fit=TRUE)
predict(fit,type="terms",se.fit=TRUE)
}
\keyword{survival}%-- one or more ...
