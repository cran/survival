# SCCS @(#)summary.survfit.s	5.1 08/30/98
summary.survfit <- function(object, times, censored=FALSE, scale=1, ...) {
    fit<-object
    if (!inherits(fit, 'survfit'))
	    stop("Invalid data")
    missing.times<-missing(times)
    n <- length(fit$time)
    stime <- fit$time/scale
    if (is.null(fit$strata)) {
	stemp <- rep(1,n)
	nstrat <- 1
	}
    else {
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat,fit$ntimes.strata)
	}

    surv <- as.matrix(fit$surv)
    if (is.null(fit$std.err)) std.err <- NULL
    else                      std.err <- fit$std.err * surv

    if (!is.null(fit$lower)) {
	lower <- as.matrix(fit$lower)
	upper <- as.matrix(fit$upper)
	}

    if (missing.times) {
	if (censored) {
	    times <- stime
	    n.risk<- fit$n.risk
	    n.event <- fit$n.event
	    }
	else {
	    who    <- (fit$n.event >0)
	    times  <-  stime[who]
	    n.risk <-  fit$n.risk[who]
	    n.event <- fit$n.event[who]
	    stemp <- stemp[who]
	    surv <- surv[who,,drop=FALSE]
	    if (!is.null(std.err)) std.err <- std.err[who,,drop=FALSE]
	    if (!is.null(fit$lower)) {
		lower <- lower[who,,drop=FALSE]
		upper <- upper[who,,drop=FALSE]
		}
	    }
	}

    else {  #this case is much harder
	if (any(times<0)) stop("Invalid time point requested")
        if(max(fit$time) < min(times))
            stop("Requested times are all beyond the end of the survival curve")
	if (length(times) >1 )
	    if (any(diff(times)<0)) stop("Times must be in increasing order")

	temp <- .C("survindex2", as.integer(n),
				  as.double(stime),
				  as.integer(stemp),
				  as.integer(length(times)),
				  as.double(times),
				  as.integer(nstrat),
				  indx = integer(nstrat*length(times)),
				  indx2= integer(nstrat*length(times)),
                   PACKAGE="survival")
	keep <- temp$indx >=0
	indx <- temp$indx[keep]
	ones <- (temp$indx2==1)[keep]
	ties <- (temp$indx2==2)[keep]  #data set time === requested time

	times <- rep(times, nstrat)[keep]
	n.risk <- fit$n.risk[indx+1 - (ties+ones)]
	surv   <- surv[indx,,drop=FALSE];   surv[ones,] <- 1
	if (!is.null(std.err)) {
	    std.err<- std.err[indx,,drop=FALSE]
	    std.err[ones,] <-0
	    }
	fit$n.event[stime>max(times)] <- 0
	n.event <- (cumsum(c(0,fit$n.event)))[ifelse(ones, indx, indx+1)]
	n.event<-  diff(c(0, n.event))

	if (!is.null(fit$lower)) {
	    lower <- lower[indx,,drop=FALSE];  lower[ones,] <- 1;
	    upper <- upper[indx,,drop=FALSE];  upper[ones,] <- 1;
	    }

	stemp <- stemp[indx]
	}

    ncurve <- ncol(surv)
    temp <- list(surv=surv, time=times, n.risk=n.risk, n.event=n.event,
			conf.int=fit$conf.int)
    if (ncurve==1) {
	temp$surv <- drop(temp$surv)
	if (!is.null(std.err)) temp$std.err <- drop(std.err)
	if (!is.null(fit$lower)) {
	    temp$lower <- drop(lower)
	    temp$upper <- drop(upper)
	    }
	}
    else {
	if (!is.null(std.err)) temp$std.err <- std.err
	if (!is.null(fit$lower)) {
	    temp$lower <- lower
	    temp$upper <- upper
	    }
	}
    if (!is.null(fit$strata))
	temp$strata <- factor(stemp,
	    labels = names(fit$strata)[sort(unique(stemp))])
    temp$call <- fit$call
    if (!is.null(fit$na.action)) temp$na.action <- fit$na.action
    class(temp) <- 'summary.survfit'
    temp
    }
