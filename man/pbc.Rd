\name{pbc}
\alias{pbc}
\docType{data}
\title{Mayo Clinic Primary Biliary Cirrhosis Data}
\description{Followup of 312 randomised and 108 unrandomised patients
  with primary biliary cirrhosis, a rare autoimmune liver disease, at
  Mayo Clinic. Most variables have some missing data, coded as -9.
}
\usage{pbc}
\format{
  \tabular{ll}{
    age:\tab in years\cr
    alb:\tab serum albumin\cr
    alkphos:\tab alkaline phosphotase\cr
    ascites:\tab presence of ascites\cr
    bili:\tab serum bilirunbin\cr
    chol:\tab serum cholesterol\cr
    edema:\tab presence of edema\cr
    edtrt:\tab 0 no edema, 0.5 untreated or successfully treated\cr
    \tab 1 unsuccessfully treated edema\cr
    hepmeg:\tab enlarged liver\cr
    time: \tab survival time\cr
    platelet:\tab platelet count\cr
    protime:\tab standardised blood clotting time\cr
    sex:\tab 1=male\cr
    sgot:\tab liver enzyme (now called AST)\cr
    spiders:\tab blood vessel malformations in the skin\cr
    stage:\tab histologic stage of disease (needs biopsy)\cr
    status:\tab censoring status\cr
    trt:\tab 1/2/-9 for control, treatment, not randomised\cr
    trig:\tab triglycerides\cr
    copper:\tab urine copper\cr
  }
}
\source{
  TR Fleming and DP Harrington (1991),
  \emph{Counting Processes and Survival Analysis}.
  Wiley, New York.
  ISBN: 0-471-52218-X.
}
\keyword{datasets}
