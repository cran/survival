\name{survreg}
\alias{survreg}
\alias{model.frame.survreg}
\alias{labels.survreg}
\alias{print.survreg.penal}
\alias{print.summary.survreg}
\alias{survReg}
\alias{anova.survreg}
\alias{vcov.survreg}
\alias{anova.survreglist}
\alias{censorReg}
\title{
Regression for a Parametric Survival Model
}
\description{
Regression for a parametric survival model. These are all
time-transformed location models, with the most useful case being the
accelerated failure models that use a log transformation.
}
\usage{
survreg(formula=formula(data), data=parent.frame(), weights, 
subset,na.action,dist="weibull",  init=NULL, scale=0,
 control=survreg.control(),parms=NULL,model=FALSE, x=FALSE,
 y=TRUE, \dots)
}
\arguments{
\item{formula}{
a formula expression as for other regression models.
See the documentation for \code{lm} and \code{formula} for details.
}
\item{data}{
optional data frame in which to interpret the variables occurring in the
formula.
}
\item{weights}{Optional observation weights}
\item{subset}{
subset of the observations to be used in the fit.
}
\item{na.action}{
function to be used to handle any NAs in the data.
}
\item{dist}{
assumed distribution for y variable. 
If the argument is a character string, then it is assumed to name an
element from \code{\link{survreg.distributions}}. These include
\code{"weibull"}, \code{"exponential"}, \code{"gaussian"}, \code{"logistic"}, \code{"lognormal"} and \code{"loglogistic"}.
Otherwise, it is assumed to be a user defined list conforming to the
format described in \code{\link{survreg.distributions}}.
}
\item{parms}{
a list of fixed parameters.  For the t-distribution for instance this is
the degrees of freedom; most of the distributions have no parameters.
}
\item{init}{
optional vector of initial values for the parameters.
}
\item{scale}{
optional fixed value for the scale.  If set to <=0 then the scale is
estimated.
}
\item{control}{
a list of control values, in the format producted by \code{survreg.control}.
}
\item{model}{
if TRUE, the model frame is returned.
}
\item{x}{
if TRUE, then the X matrix is returned.
}
\item{y}{
if TRUE, then the y vector (or survival times) is returned.
}
\item{...}{
other arguments which will be passed to \code{survreg.control}.
}}
\value{
an object of class \code{survreg} is returned.
}

\section{Compatibility note}{
This routine underwent significant changes from survival4 to
survival5. The \link{survreg.old} function gives a backwards-compatible
interface.  In S-PLUS the new function is called \code{survReg} and
the old one \code{survreg}.
}

\seealso{
\code{\link{survreg.object}}, \code{\link{survreg.distributions}},
\code{\link{pspline}}, \code{\link{frailty}}, \code{\link{ridge}},
\code{\link{survreg.old}}
}
\examples{
data(ovarian)
## These are all the same
survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull',scale=1)
survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian,
dist="exponential")
survreg.old(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='extreme',fixed=list(scale=1),link="log")

}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
