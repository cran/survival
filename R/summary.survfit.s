#SCCS @(#)summary.survfit.s	5.7 07/09/00
summary.survfit <- function(fit, times, censored=F, scale=1, extend=F, ...) {
    if (!inherits(fit, 'survfit'))
	    stop("Invalid data")

    ### In R, missing(times) stops being true when we assign to times
    missing.times<-missing(times)
    ###

    pfun <- function(nused, stime, min.time, surv, n.risk, n.event, lower,
		     upper) {
        #compute the mean, median, se(mean), and ci(median)
	minmin <- function(y, xx) {
            ww<-getOption("warn")
            options(warn=-1)
            on.exit(options(warn=ww))
	    if (any(!is.na(y) & y==.5)) {	
		if (any(!is.na(y) & y <.5))
			.5*(min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
		else .5*(min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
	        }
	    else min(xx[!is.na(y) & y<=.5])
	    }

	n <- length(stime)
	hh <- c(ifelse((n.risk[-n]-n.event[-n])==0, 0, 
		       n.event[-n]/(n.risk[-n]*(n.risk[-n]-n.event[-n]))),0)
	ndead<- sum(n.event)
	dif.time <- c(diff(c(min.time, stime)), 0)
	if (is.matrix(surv)) {
	    n <- nrow(surv)
	    mean <- dif.time * rbind(1, surv)
	    temp <- (apply(mean[(n+1):2,,drop=F], 2, cumsum))[n:1,,drop=F]
	    varmean <- c(hh %*% temp^2)
	    med <- apply(surv, 2, minmin, stime)
	    ##nused <- as.list(nused)
	    names(nused)<-NULL
	    if (!is.null(upper)) {
		upper <- apply(upper, 2, minmin, stime)
		lower <- apply(lower, 2, minmin, stime)
		cbind(nused, ndead, apply(mean, 2, sum),
		      sqrt(varmean), med, lower, upper)
	        }
	    else cbind(nused, ndead, apply(mean, 2, sum), sqrt(varmean), med)
	    }
	else {
	    mean <- dif.time*c(1, surv)
	    varmean <- sum(rev(cumsum(rev(mean))^2)[-1] * hh)
	    med <- minmin(surv, stime)
	    if (!is.null(upper)) {
		upper <- minmin(upper, stime)
		lower <- minmin(lower, stime)
		c(nused, ndead, sum(mean), sqrt(varmean), med, lower, upper)
	        }
	    else c(nused, ndead, sum(mean), sqrt(varmean), med)
	    }
        }

    n <- length(fit$time)
    stime <- fit$time/scale

    min.stime <- min(stime)
    min.time <- min(0, min.stime)
    surv <- fit$surv
    plab <- c("n", "events", "mean", "se(mean)", "median")
    if (!is.null(fit$conf.int))
	    plab2<- paste(fit$conf.int, c("LCL", "UCL"), sep='')

    #Four cases: strata Y/N  by  ncol(surv)>1 Y/N
    #  Repeat the code, with minor variations, for each one
    if (is.null(fit$strata)) {
	stemp <- rep(1,n)
	nstrat <- 1
	table <- pfun(fit$n, stime, min.time, surv, fit$n.risk, fit$n.event, 
		      fit$lower, fit$upper)
	if (is.matrix(table)) {
	    if (is.null(fit$lower)) dimnames(table) <- list(NULL, plab)
	    else dimnames(table) <- list(NULL, c(plab, plab2))
	    }
	else {
	    if (is.null(fit$lower)) names(table) <- plab
	    else names(table) <- c(plab, plab2)
  	    }
        }
    else {  #strata case
	nstrat <- length(fit$strata)
        if (is.null(fit$ntimes.strata))
            stemp <- rep(1:nstrat,fit$strata)
	else stemp <- rep(1:nstrat,fit$ntimes.strata)
	table <- NULL
	if (is.null(fit$strata.all)) strata.var <- fit$strata
	else strata.var <- fit$strata.all
	for (i in unique(stemp)) {
	    who <- (stemp==i)
	    if (is.matrix(surv)) {
		temp <- pfun(strata.var[i], stime[who], min.time, 
			     surv[who,,drop=F],
                             fit$n.risk[who], fit$n.event[who],
                             fit$lower[who,,drop=F], fit$upper[who,,drop=F])
		table <- rbind(table, temp)
            }
	    else  {
		temp <- pfun(strata.var[i], stime[who], min.time, 
			     surv[who], fit$n.risk[who], fit$n.event[who], 
			     fit$lower[who], fit$upper[who])
		table <- rbind(table, temp)
            }
        }
        
	temp <- names(fit$strata)
	if (nrow(table) > length(temp)) {
	    nrep <- nrow(table)/length(temp)
	    temp <- rep(temp, rep(nrep, length(temp)))
        }
        
	if (is.null(fit$lower))	dimnames(table) <- list(temp, plab)
	else dimnames(table) <- list(temp, c(plab, plab2))
    }
    
    surv <- as.matrix(fit$surv)
    if (is.null(fit$std.err))
        std.err <- NULL
    else
        std.err <- fit$std.err * surv

    if (!is.null(fit$lower)) {
	lower <- as.matrix(fit$lower)
	upper <- as.matrix(fit$upper)
        }

    if (missing.times) {
	if (censored) {
	    times <- stime
	    n.risk<- fit$n.risk
	    n.event <- fit$n.event
	    n.entered <- fit$enter
	    n.exit.censored <- fit$exit.censored
	    }
	else {
	    who <- (fit$n.event > 0)
	    times <- stime[who]
	    n.risk <- fit$n.risk[who]
	    n.event <- fit$n.event[who]
	    n.entered <- fit$enter[who]
	    n.exit.censored <- fit$exit.censored[who]
	    stemp <- stemp[who]
	    surv <- surv[who,,drop=F]
	    if (!is.null(std.err)) std.err <- std.err[who,,drop=F]
	    if (!is.null(fit$lower)) {
		lower <- lower[who,,drop=F]
		upper <- upper[who,,drop=F]
            }
        }
    }
    else {  #this case is much harder
	if (any(times<min.time)) stop("Invalid time point requested")
        if(any(times > max(stime)) && extend) 
            warning(paste("Time", list(times[times > max(stime)]), 
                          "beyond the end of the survival curve - using survival values of last time point", 
                          max(stime)))
	if (length(times) > 1) {
	    if (any(diff(times)<0)) stop("Times must be in increasing order")
        }
	ntime <- length(times)
	times.obs <- ntime * nstrat
        
	if (!is.null(fit$new.start)) new.start <- fit$new.start
	else new.start <- stime[1] - 1
        
	# changing NAs to -1 for C program
	upper <- ifelse(is.na(fit$upper), -1, fit$upper)
	lower <- ifelse(is.na(fit$lower), -1, fit$lower)
	std.err <- ifelse(is.na(std.err), -1, std.err)

	if (extend) num.extend <- 1
	else num.extend <- 0

	if (extend) {
	    temp.strata <- rep(ntime,nstrat)
	    times.strata <- rep(1:nstrat, temp.strata)
        }
        
	if (fit$type == 'right' || inherits(fit, 'survfit.cox')) {
	    temp <- .C("survindex2", 
		       as.integer(n),
		       as.double(stime),
		       as.integer(stemp),
		       as.integer(ntime),
		       as.double(times),
		       as.integer(nstrat),
		       as.integer(fit$n.risk),
		       as.integer(fit$n.event),
		       as.double(fit$surv),
		       as.double(std.err),
		       as.double(upper),
		       as.double(lower),
		       n.risk=integer(times.obs),
		       n.event=integer(times.obs),
		       surv=double(times.obs),
		       std.err=double(times.obs),
		       upper=double(times.obs),
		       lower=double(times.obs),
		       as.double(new.start),
		       as.integer(num.extend),
		       temp.strata=integer(nstrat),
		       temp.times=double(times.obs),PACKAGE="survival")
	    }
	if (fit$type == 'counting') {
	    temp <- .C("survindex3",
		       as.integer(n),
		       as.double(stime),
		       as.integer(stemp),
		       as.integer(ntime),
		       as.double(times),
		       as.integer(nstrat),
		       as.integer(fit$n.risk),
		       as.integer(fit$enter),
		       as.integer(fit$exit.censored),
		       as.integer(fit$n.event),
		       as.double(fit$surv),
		       as.double(std.err),
		       as.double(upper),
		       as.double(lower),
		       n.risk=integer(times.obs),
		       n.entered=integer(times.obs),
		       n.censored=integer(times.obs),
		       n.event=integer(times.obs),
		       surv=double(times.obs),
		       std.err=double(times.obs),
		       upper=double(times.obs),
		       lower=double(times.obs),
		       as.double(new.start),
		       as.integer(num.extend),
		       temp.strata=integer(nstrat),
		       temp.times=double(times.obs), PACKAGE="survival")
	    uncumm.times.strata <- rep(1:nstrat, temp$temp.strata)
	    if (!extend) 
		    temp$n.entered <- unlist(tapply(temp$n.entered[1:length(uncumm.times.strata)],
						    uncumm.times.strata, function(x) diff(c(0,x))))
	    else 
		    temp$n.entered <- unlist(tapply(temp$n.entered, times.strata, 
						    function(x) diff(c(0,x))))

	    if (!extend)
		    temp$n.censored <- unlist(tapply(temp$n.censored[1:length(uncumm.times.strata)],
						     uncumm.times.strata, function(x) diff(c(0,x))))
	    else
		    temp$n.censored <- unlist(tapply(temp$n.censored, times.strata, 
						     function(x) diff(c(0,x))))
	    }
	uncumm.times.strata <- rep(1:nstrat, temp$temp.strata)
	if (!extend)
            temp$n.event <- unlist(tapply(temp$n.event[1:length(uncumm.times.strata)], 
                                          uncumm.times.strata, function(x) diff(c(0,x))))
	else
            temp$n.event <- unlist(tapply(temp$n.event, times.strata, 
                                          function(x) diff(c(0,x))))
        ## if don't want to extend past time points, subset data
	if (!extend)
            times.strata <- rep(1:nstrat, temp$temp.strata)

        subs<-function(x,i1,i2){
            n<-length(x)
            if (i2<=n)
                return(x[-i1:-i2])
            return(x[(1:max(n,l2))[-l1:-l2]])
        }
 
	l1 <- length(times.strata)+1
	if (!extend)
            l2 <- length(temp$temp.times)
	else
            l2 <- length(times.strata)
	if (l1 > l2) l2 <- l2 + 1
	if (!extend)
            times <- subs(temp$temp.times,l1,l2)
	n.risk <- subs(temp$n.risk,l1,l2)
        n.event<-subs(temp$n.event,l1,l2)
	if (!is.null(temp$n.entered))
            n.entered <- subs(temp$n.entered,l1,l2)
	if (!is.null(temp$n.censored))
            n.exit.censored <- subs(temp$n.censored,l1,l2)
	surv <- subs(temp$surv,l1,l2)
	surv <- as.matrix(surv)

	# changing -1s back to NAs after C program 
	upper <- subs(temp$upper,l1,l2)
	upper <- ifelse(upper == -1, NA, upper)
	lower <- subs(temp$lower,l1,l2)
	lower <- ifelse(lower == -1, NA, lower)
	std.err <- subs(temp$std.err,l1,l2)
	std.err <- ifelse(std.err == -1, NA, std.err)

	std.err <- std.err
	upper <- as.matrix(upper)
	lower <- as.matrix(lower)
    }
    ncurve <- ncol(surv)
    if (fit$type == 'right' || inherits(fit, 'survfit.cox')) 
	    temp <- list(surv=surv, time=times, n.risk=n.risk, n.event=n.event,
			 conf.int=fit$conf.int, type=fit$type, table=table)
    
    if (fit$type == 'counting')
	    temp <- list(surv=surv, time=times, n.risk=n.risk, n.event=n.event,
			 conf.int=fit$conf.int, n.entered=n.entered,
			 n.exit.censored=n.exit.censored, type=fit$type, 
			 table=table)

    if (!is.null(fit$new.start)) temp$new.start <- fit$new.start
    if (!missing.times && (!is.null(fit$strata)))
        temp$times.strata <- factor(times.strata,
                                    labels=names(fit$strata)[sort(unique(times.strata))])
    if (!is.null(ncurve)) {
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
				      labels = 
				      names(fit$strata)[sort(unique(stemp))])

	temp$call <- fit$call
	if (!is.null(fit$na.action)) temp$na.action <- fit$na.action
        }
    class(temp) <- 'summary.survfit'
    temp
    }









