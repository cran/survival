\name{cluster}
\alias{cluster}
\title{
Identify Clusters
}
\description{
This is a special function used in the context of the Cox model.  It
identifies correlated groups of observations, and is used on the right hand
side of a formula.
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
\details{
The function's only action is semantic, to mark a variable as the
cluster indicator.
}
\seealso{
\code{\link{coxph}}, \code{\link{Surv}}}
\examples{\dontrun{
coxph(Surv(futime, fustat) ~ age + cluster(group), data = ovarian)
}}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
