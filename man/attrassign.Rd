\name{attrassign}
\alias{attrassign.default}
\alias{attrassign}
\alias{attrassign.lm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Create new-style "assign" attribute}
\description{
The \code{"assign"} attribute on model matrices describes which columns
come from which terms in the model formula. It has two version. R uses
the original version, but the newer version is sometimes useful.
}
\usage{
\method{attrassign}{default}(object, tt,...)
\method{attrassign}{lm}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{model matrix or linear model object}
  \item{tt}{terms object}
  \item{...}{ignored}
}

\value{
 A list with names corresponding to the term names and elements
 that are vectors indicating which columns come from which terms
}
\seealso{\code{\link{terms}},\code{\link{model.matrix}}}

\examples{
formula<-Surv(time,status)~factor(edtrt)
tt<-terms(formula)
mf<-model.frame(tt,data=pbc)
mm<-model.matrix(tt,mf)
## a few rows of data
mm[1:3,]
## old-style assign attribute
attr(mm,"assign")
## new-style assign attribute
attrassign(mm,tt)
}
\keyword{models}%-- one or more ...
