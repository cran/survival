#SCCS 04/14/92 @(#)points.survfit.s	4.2
points.survfit <- function(x, ...) {
    if (!is.matrix(x$surv))
	    points(x$time, x$surv, ...)
    else
	    matpoints(x$time, x$surv, ...)
    }
