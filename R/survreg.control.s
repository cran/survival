#SCCS  @(#)survreg.control.s	4.4 02/21/99
survreg.control <- function(maxiter=30, rel.tolerance=1e-5, failure=1,
			    toler.chol=1e-10, iter.max, debug=0,
			    outer.max = 10) {

    if (missing(iter.max)) {
	iter.max <- maxiter
	}
    else  maxiter <- iter.max
    list(iter.max = iter.max, rel.tolerance = rel.tolerance, 
	 failure =failure, toler.chol= toler.chol, debug=debug,
	 maxiter=maxiter, outer.max=outer.max)
    }
