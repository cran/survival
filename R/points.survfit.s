#SCCS 04/14/92 @(#)points.survfit.s	4.2
points.survfit <- function(object, ...) {
    if (!is.matrix(object$surv))
	    points(object$time, object$surv, ...)
    else
	    matpoints(object$time, object$surv, ...)
    }
