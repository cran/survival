# $Id: points.survfit.S 11166 2008-11-24 22:10:34Z therneau $
points.survfit <- function(x, ...) {
    if (!is.matrix(x$surv))
	    points(x$time, x$surv, ...)
    else
	    matpoints(x$time, x$surv, ...)
    }
