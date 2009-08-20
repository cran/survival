# $Id: print.coxph.null.S 11059 2008-10-23 12:32:50Z therneau $
print.coxph.null <-
 function(x, digits=max(options()$digits - 4, 3), ...)
    {
    if (!is.null(cl<- x$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    cat("Null model\n  log likelihood=", format(x$loglik), "\n")
    omit <- x$na.action
    if (length(omit))
	cat("  n=", x$n, " (", naprint(omit), ")\n",
				sep="")
    else cat("  n=", x$n, "\n")
    }
