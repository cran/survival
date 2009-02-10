# $Id: points.survfit.S 11059 2008-10-23 12:32:50Z therneau $
points.survfit <- function(x, ...) {
    if (!is.matrix(x$surv))
	    points(x$time, x$surv, ...)
    else
	    matpoints(x$time, x$surv, ...)
    }
