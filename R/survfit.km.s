#SCCS @(#)survfit.km.s	4.16 07/09/00
survfit.km <- function(x, y, casewt=rep(1,n),
		       type=c('kaplan-meier', 'fleming-harrington', 'fh2'),
		       error=c('greenwood', "tsiatis"), se.fit=TRUE,
		       conf.int= .95,
		       conf.type=c('log',  'log-log',  'plain', 'none'),
		       conf.lower=c('usual', 'peto', 'modified'),
		       new.start) {
    type <- match.arg(type)
    method <- match(type, c("kaplan-meier", "fleming-harrington", "fh2"))

    error <- match.arg(error)
    error.int <- match(error, c("greenwood", "tsiatis"))
    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)

    if (!is.Surv(y)) stop("y must be a Surv object")
    if (!is.factor(x)) stop("x must be a factor")
    if (attr(y, 'type') != 'right' && attr(y, 'type') != 'counting')
	    stop("Can only handle right censored or counting data")

    ny.all <- ncol(y) # getting total number of observations and strata
    n.all <- nrow(y)  # before possible subsetting
    sorted <- (1:n.all)[order(x, y[,ny.all-1])]
    strata.all.temp <- as.numeric(x[sorted])
    strata.all <- as.vector(table(strata.all.temp))

    if (!missing(new.start)) { # checking if not starting at first time
	keep <- (y[,ny.all-1] >= new.start)
	if (all(keep==FALSE))
		stop(paste("new.start =", new.start,
			   "is greater than all time points."))
	x <- x[keep]
	y <- y[keep,,drop=FALSE]
        }

    ny <- ncol(y) # getting working total observations and strata
    n <- nrow(y)
    sorted <- (1:n)[order(x, y[,ny-1])]
    y <- y[sorted,]
    strata.temp <- as.numeric(x[sorted])
    newstrat <- 1
    if (n > 1) {
	newstrat <- as.integer(c(1*(diff(strata.temp)!=0), 1))
	if (sum(newstrat) > n/2)
		warning("Number of strata > number of observations/2")
        }
    if (method==3 && any(floor(casewt) != casewt))
	    stop("The fh2 method is not valid for fractional case weights")

    nstrat <- length(unique(strata.temp))

    if (attr(y, 'type') == 'right') # list of times is created differently 
	    times <- y[,1, drop=FALSE]  # depending on censoring
    else if (attr(y, 'type') == 'counting')
	    times <- cbind(y[,1],y[,2])

    sort.times <- tapply(times, strata.temp[row(times)], 
			 function(x) sort(unique(x)))
    ntimes.strata <- sapply(sort.times, length)
    times <- unlist(sort.times)

    if (attr(y, 'type') == 'right') # call correct surv program based 
	    surv <- .C("survfit2",  # on censoring 
		       as.integer(n),
		       y = as.double(y),
		       as.double(casewt[sorted]),
		       strata= as.integer(newstrat),
		       as.integer(method),
		       as.integer(error.int),
		       mark=double(n),
		       surv=double(n),
		       varhaz=double(n),
		       risksum=double(n), PACKAGE="survival")
    else if (attr(y, 'type') == 'counting')
	    surv <- .C("survfit3",
		       as.integer(n),
		       as.double(y),
		       as.double(casewt[sorted]),
		       strata=as.integer(newstrat),
		       as.integer(method),
		       as.integer(error.int),
		       as.integer(nstrat),
		       as.double(ntimes.strata),
		       y=as.double(times),
		       mark=double(length(times)),
		       surv=double(length(times)),
		       varhaz=double(length(times)),
		       risksum=double(length(times)),
		       enter=double(length(times)),
		       exit.censored=double(length(times)), PACKAGE="survival")
    ntime <- length(times)
    if (error.int==1) surv$varhaz[surv$surv==0] <- NA
    ntime <- 1:ntime

    if (attr(y, 'type') == 'right') {
	if (nstrat == 1)
		temp <- list(n=n.all, # total number of obs
			     time=surv$y[ntime],
			     n.risk=surv$risksum[ntime],
			     n.event=surv$mark[ntime],
			     surv=surv$surv[ntime],
			     type=attr(y, 'type')) # type of censoring
	else {
	    temp <- surv$strata[1:nstrat]
	    tstrat <- diff(c(0, temp)) # number in each strata
	    names(tstrat) <- levels(x)[1:nstrat] # number in strata may differ
	    names(strata.all) <- levels(x)
	    temp <- list(n=n.all,
			 time=surv$y[ntime],
			 n.risk=surv$risksum[ntime],
			 n.event=surv$mark[ntime],
			 surv=surv$surv[ntime],
			 type=attr(y, 'type'),
			 ntimes.strata=ntimes.strata,
			 strata=tstrat,
			 strata.all=strata.all)
	    }
        }
    else if (attr(y, 'type') == 'counting') {
	if (nstrat == 1) 
		temp <- list(n=n.all,
			     time=surv$y[ntime],
			     n.risk=surv$risksum[ntime],
			     n.event=surv$mark[ntime],
			     surv=surv$surv[ntime],
			     type=attr(y, 'type'),
			     enter=surv$enter[ntime],
			     exit.censored=surv$exit.censored[ntime])
	else {	
	    temp <-surv$strata[1:nstrat]
	    tstrat <- diff(c(0, temp))
	    names(tstrat) <- levels(x)[1:nstrat] # number of strata may differ
	    names(strata.all) <- levels(x)
	    temp <- list(n=n.all,
			 time=surv$y[ntime],
			 ntimes.strata=ntimes.strata,
			 n.risk=surv$risksum[ntime],
			 n.event=surv$mark[ntime],
			 surv=surv$surv[ntime],
			 type=attr(y, 'type'),
			 strata=tstrat,
			 strata.all=strata.all,
			 enter=surv$enter[ntime],
			 exit.censored=surv$exit.censored[ntime])
	    }
        }

    if (!missing(new.start))
	    temp$new.start <- new.start # user defined time to start

    if (se.fit) {
	std.err <- sqrt(surv$varhaz[ntime])
	temp$std.err <- std.err
	events <- temp$n.event >0
	n.lag <- rep(c(temp$n.risk[1], temp$n.risk[events]),
		     diff(c(ntime[1], ntime[events], 1+max(ntime))))
	std.low <- switch(conf.lower,
			  'usual' = std.err,
			  'peto' = sqrt((1-temp$surv)/ temp$n.risk),
			  'modified' = std.err * sqrt(n.lag/temp$n.risk))
	zval <- qnorm(1- (1-conf.int)/2, 0,1)

	if (conf.type=='plain') {
	    temp1 <- temp$surv + zval* std.err * temp$surv
	    temp2 <- temp$surv - zval* std.low * temp$surv
	    temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
				 conf.type='plain', conf.int=conf.int))
	    }

	if (conf.type=='log') {
	    #avoid some "log(0)" messages
	    xx <- ifelse(temp$surv==0,1,temp$surv)  

	    temp1 <- ifelse(temp$surv==0, NA, exp(log(xx) + zval* std.err))
	    temp2 <- ifelse(temp$surv==0, NA, exp(log(xx) - zval* std.low))
	    temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
				 conf.type='log', conf.int=conf.int))
	    }

	if (conf.type=='log-log') {
	    who <- (temp$surv==0 | temp$surv==1) #special cases
	    temp3 <- ifelse(temp$surv==0, NA, 1)
	    xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	    temp1 <- exp(-exp(log(-log(xx)) + zval*std.err/log(xx)))
	    temp1 <- ifelse(who, temp3, temp1)
	    temp2 <- exp(-exp(log(-log(xx)) - zval*std.low/log(xx)))
	    temp2 <- ifelse(who, temp3, temp2)
	    temp <- c(temp, list(upper=temp1, lower=temp2,
				 conf.type='log-log', conf.int=conf.int))
	    }
        }
    temp
    }














