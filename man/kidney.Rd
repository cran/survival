\name{kidney}
\alias{kidney}
\title{Kidney data from survival5  }
\format{
\tabular{ll}{
patient:\tab id\cr
time:\tab time\cr
status:\tab event status\cr
age:\tab in years\cr
sex:\tab 1=male, 2=female\cr
disease:\tab  disease type (0=GN, 1=AN, 2=PKD, 3=Other)\cr
frail:\tab frailty estimate from original paper\cr
}}
\description{
Data on the recurrence times to infection, at the point of insertion of
  the catheter, for kidney patients using portable dialysis equipment.
  Catheters may be removed for reasons other than infection, in which case
  the observation is censored.  Each patient has exactly 2 observations.
}
\section{Note}{
  The original analysis had incorrect handling of ties and so is not
  exactly reproduced by survival.
}
\examples{
data(kidney)
kfit <- coxph(Surv(time, status)~ age + sex + disease + frailty(id), kidney)
kfit0 <- coxph(Surv(time, status)~ age + sex + disease, kidney)
kfitm1 <- coxph(Surv(time,status) ~ age + sex + disease + 
		frailty(id, dist='gauss'), kidney)

}
\source{McGilchrist and Aisbett, Biometrics 47, 461-66, 1991}
\keyword{survival}
