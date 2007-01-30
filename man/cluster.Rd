\name{cluster}
\alias{cluster}
\title{
Identify Clusters
}
\description{
This is a special function used in the context of survival models.  It
identifies correlated groups of observations, and is used on the right hand
side of a formula.  Using \code{cluster()} in a formula implies that
robust sandwich variance estimators are desired.
}
\usage{
cluster(x)
}
\arguments{
\item{x}{
A character, factor, or numeric variable.
}}
\value{
\code{x}
}
\seealso{
\code{\link{coxph}}, \code{\link{Surv}}}
\examples{
frailty.model <- coxph(Surv(time, status)~ age + sex + disease + frailty(id), kidney)
marginal.model <- coxph(Surv(time, status)~ age + sex + disease+cluster(id), kidney)
summary(frailty.model)
summary(marginal.model)

frailty.model <- survreg(Surv(time, status) ~ rx +  frailty(litter), rats )
marginal.model <- survreg(Surv(time, status) ~ rx +  cluster(litter), rats )
summary(frailty.model)
summary(marginal.model)


}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
