\name{predict.survreg}
\alias{predict.survreg}
\alias{predict.survreg.penal}
\title{
Predicted Values for a `survreg' Object
}
\description{
Predicted values for a \code{survreg} object
}
\usage{
\method{predict}{survreg}(object, newdata, 
type=c("response", "link", "lp", "linear",  "terms", "quantile", 
	"uquantile"), 
se.fit=FALSE, terms=NULL, p=c(0.1, 0.9),...)
}
\arguments{
\item{object}{
result of a model fit using the \code{survreg} function.
}
\item{newdata}{
data for prediction.  If absent, predictions are for the
subjects used in the original fit.
}
\item{type}{
the type of predicted value. 
This can be on the original scale of the data (response),
the linear predictor (\code{"linear"}, with \code{"lp"} as an allowed abbreviation),
a predicted quantile on the original scale of the data (\code{"quantile"}),
a quantile on the linear predictor scale (\code{"uquantile"}),
or the matrix of terms for the linear predictor (\code{"terms"}).
At this time \code{"link"} and linear predictor (\code{"lp"}) are identical.
}
\item{se.fit}{
if TRUE, include the standard errors of the prediction in the result.
}
\item{terms}{
subset of terms.  The default for residual type \code{"terms"} is a matrix with
one column for every term (excluding the intercept) in the model.
}
\item{p}{
vector of percentiles.  This is used only for quantile predictions.
}\item{...}{other arguments}
}
\value{
a vector or matrix of predicted values.
}
\references{
Escobar and Meeker (1992). Assessing influence in regression analysis with
censored data. \emph{Biometrics,} 48, 507-528.
}
\seealso{
\code{\link{survreg}}, \code{\link{residuals.survreg}}}
\examples{
# Draw figure 1 from Escobar and Meeker
data(stanford2)
fit <- survreg(Surv(time,status) ~ age + age^2, data=stanford2,
	dist='lognormal')
plot(stanford2$age, stanford2$time, xlab='Age', ylab='Days',
	xlim=c(0,65), ylim=c(.01, 10^6), log='y')
pred <- predict(fit, newdata=list(age=1:65), type='quantile',
	         p=c(.1, .5, .9))
matlines(1:65, pred, lty=c(2,1,2), col=1)
}
\keyword{survival}

