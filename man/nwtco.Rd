\name{nwtco}
\alias{nwtco}
\docType{data}
\title{Data from the National Wilm's Tumor Study}
\description{
  Missing data/masurement error example. Tumor histology predicts
  survival, but prediction is stronger with central lab histology than
  with the local institution determination.
}
\usage{data(nwtco)}
\format{
  A data frame with 4028 observations on the following 9 variables.
  \describe{
    \item{\code{seqno}}{id number}
    \item{\code{instit}}{Histology from local institution}
    \item{\code{histol}}{Histology from central lab}
    \item{\code{stage}}{Disease stage}
    \item{\code{study}}{study}
    \item{\code{rel}}{indicator for relapse}
    \item{\code{edrel}}{time to relapse}
    \item{\code{age}}{age in months}
    \item{\code{in.subcohort}}{Included in the subcohort for the example in the
    paper}
  }
}
\source{
 \url{http://faculty.washington.edu/norm/software.html}
}
\references{
 N.E. Breslow and N. Chatterjee, "Design and analysis of two-phase studies with binary outcome applied to Wilms tumour prognosis"  Applied Statistics  48:457-68, 1999
}
\examples{
data(nwtco)
with(nwtco, table(instit,histol))
anova(coxph(Surv(edrel,rel)~histol+instit,data=nwtco))
anova(coxph(Surv(edrel,rel)~instit+histol,data=nwtco))
  }
\keyword{datasets}
