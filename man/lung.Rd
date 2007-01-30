\name{lung}
\docType{data}
\alias{cancer}
\alias{lung}
\title{Mayo Clinic Lung Cancer Data}
\description{
  Survival in patients with lung cancer at Mayo Clinic.  Performance
  scores rate how well the patient can perform usual daily activities.
}
\usage{
lung
cancer
}
\format{
  \tabular{ll}{
    inst:\tab Institution code\cr
    time:\tab Survival time in days\cr
    status:\tab censoring status 1=censored, 2=dead\cr
    age:\tab Age in years\cr
    sex:\tab  Male=1 Female=2\cr
    ph.ecog:\tab ECOG performance score (0=good 5=dead)\cr
    ph.karno:\tab Karnofsky performance score (bad=0-good=100) rated by physician\cr
    pat.karno:\tab Karnofsky performance score  rated by patient\cr
    meal.cal:\tab Calories consumed at meals\cr
    wt.loss:\tab Weight loss in last six months\cr
  }
}
\source{Terry Therneau}
\keyword{datasets}
