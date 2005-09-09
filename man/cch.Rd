\alias{cch}
\name{cch}
\title{Fits proportional hazards regression model to case-cohort data}

\description{
Returns estimates and standard errors from relative risk regression fit to
data from case-cohort studies. A choice is available among the Prentice,
Self-Prentice and Lin-Ying methods for estimation of regression coefficients 
and standard errors.}
\usage{
cch(formula, data = sys.parent(), subcoh, id, cohort.size,
    method = c("Prentice", "SelfPrentice", "LinYing"))
}
\arguments{
\item{formula}{
A formula object that must have a \code{\link{Surv}} object as the response. 
The Surv object must be of type \code{"right"}, or of type \code{"counting"}.
}
\item{subcoh}{
  Vector of indicatorsfor subjects sampled as part of the
  sub-cohort. Code \code{1} or \code{TRUE} for members of the
  sub-cohort, \code{0} or \code{FALSE} for others. If \code{data} is a
  data frame then \code{subcoh} may be a one-sided formula.
}
\item{id}{
Vector of unique identifiers, or formula specifying such a vector. If \code{data} is a
  data frame then \code{subcoh} may be a one-sided formula.
}
\item{cohort.size}{
Scalar with size of original cohort from which subcohort was sampled
}
\item{data}{
An optional data frame in which to interpret the variables 
occurring in the formula.
}
\item{method}{
Three procedures are available. The default method is "Prentice", with 
options for "SelfPrentice" or "LinYing". 
}}
\value{
A whole bunch of stuff including list of estimated regression coefficients and two estimates of their 
asymptotic variance-covariance matrix.
\item{coef}{
regression coefficients.

}
\item{naive.var}{
Self-Prentice model based variance-covariance matrix.

}
\item{var}{
Lin-Ying empirical variance-covariance matrix. 

}}
\details{
Implements methods for case-cohort data analysis described by Therneau and
Li (1999). The three methods differ in the choice of "risk sets" used to
compare the covariate values of the failure with those of others at risk at
the time of failure. "Prentice" uses the sub-cohort members "at risk" plus
the failure if that occurs outside the sub-cohort and is score unbiased.
"SelfPren" (Self-Prentice) uses just the sub-cohort members "at risk". These
two have the same asymptotic variance-covariance matrix. "LinYing" (Lin-Ying)
uses the all members of the sub-cohort and all failures outside the sub-cohort
who are "at risk". The methods also differ in the weights given to different
score contributions.          
}
\author{Norman Breslow, modified by Thomas Lumley}
\references{
Prentice, RL (1986). A case-cohort design for epidemiologic cohort studies and
disease prevention trials. Biometrika 73: 1--11.


Self, S and Prentice, RL (1988). Asymptotic distribution theory and efficiency
results for case-cohort studies. Annals of Statistics 16: 64--81.


Lin, DY and Wei, LJ (1989). The robust inference for the Cox proportional hazards
model. Journal of the American Statistical Association 84: 1074--1078.


Lin, DY and Ying, Z (1993). Cox regression with incomplete covariate measurements.
Journal of the American Statistical Association 88: 1341--1349.


Barlow, WE (1994). Robust variance estimation for the case-cohort design. Biometrics
50: 1064--1072


Therneau, TM and Li, H (1999). Computing the Cox model for case-cohort designs.
Lifetime Data Analysis 5: 99--112.


Borgan, O et al. (1999). Exposure stratified case-cohort designs. Lifetime Data
Analysis 

}
\examples{
## The complete Wilms Tumor Data 
## (Breslow and Chatterjee, Applied Statistics, 1999)
## subcohort selected by simple random sampling.
##

data(nwtco)

subcoh <- nwtco$in.subcohort
selccoh <- with(nwtco, rel==1|subcoh==1)
ccoh.data <- nwtco[selccoh,]
ccoh.data$subcohort <- subcoh[selccoh]
ccoh.data$histol <- factor(ccoh.data$histol,labels=c("FH","UH"))  # Central histology 
ccoh.data$stage <- factor(ccoh.data$stage,labels=c("I","II","III","IV")) # Stage
ccoh.data$age <- ccoh.data$age/12 # Age in years

##
## Standard case-cohort analysis: simple random subcohort 
##

fit.ccP <- cch(Surv(edrel, rel) ~ stage + histol + age, data =ccoh.data,
   subcoh = ~subcohort, id=~seqno, cohort.size=4028)


fit.ccP

fit.ccSP <- cch(Surv(edrel, rel) ~ stage + histol + age, data =ccoh.data,
   subcoh = ~subcohort, id=~seqno, cohort.size=4028, method="SelfPren")

summary(fit.ccSP)
}
\keyword{survival}
