# SCCS @(#)lines.survfit.s	4.16  01/14/99
lines.survfit <- function(x, type='s', mark=3, col=1, lty=1, lwd=1,
			  mark.time =T, xscale=1, 
			  firstx=0, firsty=1, xmax, fun,
			  conf.int=F, ...) {

    if (inherits(x, 'survexp')) {
	if (missing(type)) type <- 'l'
	if (!is.numeric(mark.time)) mark.time <- F
	}
    if (inherits(x, 'survfit.coxph')) {
	if (!is.numeric(mark.time)) mark.time <- F
	}

    if (is.character(conf.int)) {
	if (conf.int=='only') {
	    conf.int <- T
	    plot.surv<- F
	    }
	else stop("Unrecognized option for conf.int")
	}
    else plot.surv <- T

    if (is.numeric(mark.time)) mark.time<- sort(unique(mark.time[mark.time>0]))

    if (is.matrix(x$surv)) {
	ncol.per.strat <- ncol(x$surv)
	ncurve <- ncol(x$surv)
	coffset <- nrow(x$surv)*(1:ncurve -1)     #within matrix offset
        }
    else {
	ncol.per.strat <- 1
	ncurve <- 1
	coffset <- 0
        }

    if (is.null(x$strata)) {
	nstrat <- 1
	stemp <- rep(1, length(x$time))
	}
    else {
	nstrat <- length(x$strata)
	ncurve <- ncurve * nstrat
	stemp <- rep(1:nstrat, x$strata)
	}

    ssurv <- x$surv
    stime <- x$time
    supper <- x$upper
    slower <- x$lower
    if (!missing(xmax) && any(x$time>xmax)) {
	# prune back the survival curves
	# I need to replace x's over the limit with xmax, and y's over the
	#  limit with either the prior y value or firsty
	keepx <- keepy <- NULL  # lines to keep
	yzero <- NULL           # if all points on a curve are < xmax
	tempn <- table(stemp)
	offset <- cumsum(c(0, tempn))
	for (i in 1:nstrat) {
	    ttime <-stime[stemp==i]
	    if (all(ttime <= xmax)) {
		keepx <- c(keepx, 1:tempn[i] + offset[i])
		keepy <- c(keepy, 1:tempn[i] + offset[i])
		}
	    else {
		bad <- min((1:tempn[i])[ttime>xmax])
		if (bad==1)  {
		    keepy <- c(keepy, 1+offset[i])
		    yzero <- c(yzero, 1+offset[i])
		    }
		else  keepy<- c(keepy, c(1:(bad-1), bad-1) + offset[i])
		keepx <- c(keepx, (1:bad)+offset[i])
		stime[bad+offset[i]] <- xmax
		x$n.event[bad+offset[i]] <- 1   #don't plot a tick mark
		}
	    }

	# ok, now actually prune it
	stime <- stime[keepx]
	stemp <- stemp[keepx]
	x$n.event <- x$n.event[keepx]
	if (is.matrix(ssurv)) {
	    if (length(yzero)) ssurv[yzero,] <- firsty
	    ssurv <- ssurv[keepy,,drop=F]
	    if (!is.null(supper)) {
		if (length(yzero)) supper[yzero,] <- slower[yzero,] <- firsty
		supper <- supper[keepy,,drop=F]
		slower <- slower[keepy,,drop=F]
		}
	    }
	else {
	    if (length(yzero)) ssurv[yzero] <- firsty
	    ssurv <- ssurv[keepy]
	    if (!is.null(supper)) {
		if (length(yzero)) supper[yzero] <- slower[yzero] <- firsty
		supper <- supper[keepy]
		slower <- slower[keepy]
		}
	    }
	}
	stime <- stime/xscale
    	
    if (!missing(fun)) {
	if (is.character(fun)) {
	    tfun <- switch(fun,
		            'log' = function(x) x,
			    'event'=function(x) 1-x,
			    'cumhaz'=function(x) -log(x),
			    'cloglog'=function(x) log(-log(x)),
			    'pct' = function(x) x*100,
			    'logpct'= function(x) 100*x,
			    stop("Unrecognized function argument")
			    )
	    }
	else if (is.function(fun)) tfun <- fun
	else stop("Invalid 'fun' argument")
	
	ssurv <- tfun(ssurv)
	if (!is.null(supper)) {
	    supper <- tfun(supper)
	    slower <- tfun(slower)
	    }
	firsty <- tfun(firsty)
        }
    else {
	firsty  <- firsty
	}

    strata <- table(stemp)
    soffset<- ncol.per.strat * c(0, cumsum(strata))
    mark <- rep(mark, length=ncurve)
    col  <- rep(col , length=ncurve)
    lty  <- rep(lty , length=ncurve)
    lwd  <- rep(lwd , length=ncurve)
    time <- rep(stime, ncol.per.strat)


    if (type=='s') {
	type=='l'
	dostep <- function(x,y) {
	    n <- length(x)
	    if (n >2) {
		# replace verbose horizonal sequences like
		# (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
		# with (1, .2), (3, .1).  They are slow, and can smear the 
		# looks of the line type.
		dupy <- c(T, diff(y[-n]) !=0, T)
		n2 <- sum(dupy)
		
		#create a step function
		xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
		yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
		list(x=xrep, y=yrep)
		}
	    else if (n==1) list(x=x, y=y)
	    else  list(x=x[c(1,2,2)], y=y[c(1,1,2)])
	    }	
	}
    else dostep <- function(x,y) list(x=x, y=y)

    k <- 0
    xend <- yend <- NULL
    for (i in 1:nstrat) {
      for (j in 1:ncol.per.strat) {
	k <- k +1  
	who <- seq(soffset[i]+ coffset[j]+1, length=strata[i])  
	if (is.finite(firstx) && is.finite(firsty)) {
	    xx <- c(firstx, time[who])
	    yy <- c(firsty, ssurv[who])
	    yyu<- c(firsty, supper[who])
	    yyl<- c(firsty, slower[who])
	    deaths <- c(-1, ssurv$n.event[who])
	    }
	else {
	    xx <- time[who]
	    yy <- ssurv[who]
	    yyu<- supper[who]
	    yyl<- slower[who]
	    deaths <- ssurv$n.event[who]
	    }
	nn <- length(xx)

	if (conf.int) {
	    lines(dostep(xx,yyl), type=type, col=col[k], 
		  lty=lty[k], lwd=lwd[k], ...)
	    lines(dostep(xx, yyu), type=type, col=col[k], 
		  lty=lty[k], lwd=lwd[k], ...)
	    }

	xend _ c(xend,max(xx))
	yend _ c(yend,min(yy))
	if (plot.surv) { 
	    lines(dostep(xx, yy), type=type, col=col[k], 
		      lty=lty[k], lwd=lwd[k], ...)
	    if (is.numeric(mark.time)) {
		indx <- mark.time
		for (k in seq(along=mark.time))
			indx[k] <- sum(mark.time[k] > xx)
		points(mark.time[indx<nn], yy[indx[indx<nn]],
		       pch=mark[k],col=col[k], ...)
		}
	    else if (mark.time==T) {
		if ( any(deaths==0))
			points(xx[deaths==0], yy[deaths==0],
				   pch=mark[k],col=col[k], ...)
		}
	    }
	}
      }
    invisible(list(x=xend, y=yend))
    }
