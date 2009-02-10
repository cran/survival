# $Id: lines.survfit.coxph.S 10788 2008-09-18 00:48:23Z therneau $
lines.survfit.coxph <- function(x, mark.time=FALSE, ...) {
    if (is.logical(mark.time) & mark.time)
	    stop("Invalid value for mark.time")
    invisible(NextMethod('lines', mark.time=mark.time))
    }
