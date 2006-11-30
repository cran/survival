\name{survival-internal}
\alias{survival-internal}
\alias{survreg.fit}
\alias{survpenal.fit}
\alias{survdiff.fit}
\alias{is.category}
\alias{is.na.ratetable2}
\alias{is.na.coxph.penalty}
\alias{match.ratetable}
\alias{as.matrix.ratetable}
\alias{model.newframe}
\alias{as.data.frame.difftime}
\alias{coxpenal.df}
\alias{coxpenal.fit}
\alias{coxph.fit}
\alias{coxph.wtest}
\alias{agreg.fit}
\alias{agexact.fit}
\title{Internal survival functions}
\description{Internal survival functions}
\usage{
survreg.fit(x, y, weights, offset, init, controlvals, dist, scale = 0,
    nstrat = 1, strata, parms = NULL)
survpenal.fit(x, y, weights, offset, init, controlvals, dist, scale = 0,
    nstrat = 1, strata, pcols, pattr, assign, parms = NULL)
survdiff.fit(y, x, strat, rho = 0)
is.category(x)
match.ratetable(R, ratetable)
\method{as.matrix}{ratetable}(x)
is.na.ratetable2(x)
is.na.coxph.penalty(x)
model.newframe(object, newdata, safe = FALSE, response = FALSE, ...)
\method{as.data.frame}{difftime}(x, row.names = NULL, optional = FALSE,...)
coxpenal.df(hmat, hinv, fdiag, assign.list, ptype, nvar, pen1,
    pen2, sparse)
coxpenal.fit(x, y, strata, offset, init, control, weights, method,
    rownames, pcols, pattr, assign)
coxph.fit(x, y, strata, offset, init, control, weights, method,
    rownames)
coxph.wtest(var, b, toler.chol = 1e-09)
agreg.fit(x, y, strata, offset, init, control, weights, method,
    rownames)
agexact.fit(x, y, strata, offset, init, control, weights, method,
    rownames) 
}
\details{These are not to be called by the user}
\keyword{survival}
\keyword{internal}
