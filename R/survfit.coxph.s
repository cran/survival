# SCCS @(#)survfit.coxph.s	5.6 07/09/00

survfit.coxph <-
  function(object, newdata, se.fit=T, conf.int=.95, individual=F,
	    type, vartype,
	    conf.type=c('log', 'log-log', 'plain', 'none'),
	    call = match.call()) {

    if(!is.null((object$call)$weights))
	stop("Survfit cannot (yet) compute the result for a weighted model")
    call <- match.call()
    Terms <- terms(object)
    strat <- attr(Terms, "specials")$strata
    cluster<-attr(Terms, "specials")$cluster
    if (length(cluster)) {
	temp <- untangle.specials(Terms, 'cluster')
	Terms <- Terms[-temp$terms]
	}
    resp <-  attr(Terms, "variables")[attr(Terms, "response")]
    n <- object$n
    nvar <- length(object$coef)
    score <- exp(object$linear.predictor)
    
    temp <- c('aalen', 'kalbfleisch-prentice', 'efron',
	           'tsiatis', 'breslow', 'kaplan-meier', 'fleming-harringon',
	           'greenwood', 'exact')
    temp2 <- c(2,1,3,2,2,1,3,1,1)
    if (missing(type)) type <- object$method
    if (missing(vartype)) vartype <- type
    method <- temp2[match(match.arg(type, temp), temp)]
    if (is.na(method)) stop("Invalid survival curve type")
    vartype <- temp2[match(match.arg(vartype, temp), temp)]
    if (is.na(vartype)) stop("Invalid variance type specified")
    if (!se.fit) conf.type <- 'none'
    else conf.type <- match.arg(conf.type)

    # Recreate a copy of the data
    #  (The coxph.getdata routine never returns cluster() terms).
    data <- coxph.getdata(object, y=T, x=se.fit,
			           strata=(length(strat)))
    y <- data$y
    ny <- ncol(y)
    if (nrow(y) != n) stop ("Mismatched lengths: logic error")
    if (length(strat)) strata.all <- table(data$strata)

    # Get the sort index for the data, and add a column to y if
    #  necessary to make it of the "counting process" type  (I only
    #  wrote 1 C routine to handle both cases).
    #
    type <- attr(y, 'type')
    if (type=='counting') {
	if (method=='kaplan-meier')
	      stop ("KM method not valid for counting type data")
	if (length(strat)) ord <- order(data$strata, y[,2], -y[,3])
	else               ord <- order(y[,2], -y[,3]) 
        }
    else if (type=='right') {
	if (length(strat)) ord <- order(data$strata, y[,1], -y[,2])
	else               ord <- order(y[,1], -y[,2]) 
	miny <- min(y[,1])
	if (miny < 0) y <- cbind(2*miny -1, y)
	else          y <- cbind(-1, y)
	}
    else stop("Cannot handle \"", type, "\" type survival data")

    if (!is.null(data$x)) x <- data$x[ord,]
    else                  x <- 0
    if (is.null(object$weights))  weights <- rep(1,n)
    else                          weights <- object$weights[ord]

    # Create a 'nice' strata vector for the C code, 1 at the last obs of
    #   each strata, and 0 otherwise
    if (length(strat)) {
	newstrat <- (as.numeric(data$strata))[ord]
	newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
	}
    else newstrat <- as.integer(rep(0,n))
    newstrat[n] <- 1

    # Which type of curve do I need?
    #  1: the new data gives a covariate path for a single individual,
    #          one curve as a result
    #  2: each line of new data is "covariates over all time", and 
    #          gives rise to a separate curve
    #
    if (individual && !missing(newdata)) stype <- 1
    else {
	stype <- 2
	# Don't need (or want) strata term if it is there
	if (length(strat)) {
	    temp <- untangle.specials(Terms, 'strata')
	    Terms <- Terms[-temp$terms]
	    }
	}
    if (stype==1 && method != vartype)
	    stop("The type and vartype args must agree for individual=T")
    if (stype==1 && method==1)
	    stop("Only Aalen and F-H estimates available for individual=T")

    #
    # Get the second, "new" data set.  By default the new curve is
    #  produced around the mean of the old ones.  It the new data set
    #  is missing, use the old data set along with the mean of the old
    #  data set, but NOT the mean of the old offset variables.
    #  
    offset2 <- 0   #offset variable for the new data set
    if (!missing(newdata)) {
	m2 <- model.newframe(Terms, newdata, response=(stype==1))
	if (!inherits(m2, 'data.frame'))  {
	    x2 <- as.matrix(m2)
	    if (ncol(x2) != nvar) stop ("Wrong # of variables in new data")
	    n2 <- nrow(x2)
	    if (stype==1) stop("Program error #3")
	    }

	else  {
	    x2 <- model.matrix(delete.response(Terms), m2)[,-1,drop=F]
	    n2 <- nrow(x2)
	    offset2 <- model.extract(m2, 'offset')
	    if (is.null(offset2)) offset2 <- 0
	    if (stype==1) {
		#
		# The case of an agreg, with a multiple line newdata
		#
		strata.all <- object$n
		if (length(strat)) {
		    strata2 <- factor(x2[,strat], levels=levels(stratum))
		    x2 <- x2[, -strat, drop=F]
		    }
		else strata2 <- rep(1, nrow(x2))
		y2 <- model.extract(m2, 'response')
		if (attr(y2,'type') != type)
		    stop("Survival type of newdata does not match the fitted model")
		if (nrow(y2) != n2) stop("Wrong # of rows for Y")
		}
	    }
	}
    else x2 <- matrix(object$means, nrow=1)
    n2 <- nrow(x2)

    # Compute risk scores for the new subjects
    coef <- ifelse(is.na(object$coef), 0, object$coef)
    newrisk <- exp(c(x2 %*% coef) + offset2 - sum(coef*object$means))

    dimnames(y) <- NULL   #I only use part of Y, so names become invalid
    storage.mode(y) <- 'double'
    ndead <- sum(y[,3])
    if (stype==1) {
	surv <- .C("agsurv1", as.integer(n),
			     as.integer(nvar),
			     y[ord,],
			     as.double(score[ord]),
			     strata=as.integer(newstrat),
			     surv=double(ndead*n2),
			     varhaz=double(ndead*n2),
			     nsurv=as.integer(method==3),
			     as.double(x),
			     double(3*nvar),
			     as.double(object$var),
			     y = double(3*n*n2),
			     as.integer(n2),
			     as.double(y2),
			     as.double(x2),
			     as.double(newrisk),
			     as.integer(strata2), PACKAGE="survival" )
	ntime <- 1:surv$nsurv
	temp <- (matrix(surv$y, ncol=3))[ntime,]
	temp <- list(n=n, time = temp[,1],
		     n.risk= temp[,2],
		     n.event=temp[,3],
		     surv = surv$surv[ntime],
		     type=type)
	if (se.fit) temp$std.err <- sqrt(surv$varhaz[ntime])
	}
    else {
	surv <- .C('agsurv2', as.integer(n),
			      as.integer(nvar* se.fit),
			      y = y[ord,],
			      as.double(score[ord]),
			      strata = as.integer(newstrat),
			      surv = double(ndead*n2),
			      varhaz = double(ndead*n2),
			      as.double(x),
			      as.double(object$var),
			      nsurv = as.integer(c(method, vartype)),
			      double(3*nvar),
			      as.integer(n2),
			      as.double(x2),
			      as.double(newrisk), PACKAGE="survival")
	nsurv <- surv$nsurv[1]
	ntime <- 1:nsurv
	if (n2>1) {
	    tsurv <- matrix(surv$surv[1:(nsurv*n2)], ncol=n2)
	    tvar  <- matrix(surv$varhaz[1:(nsurv*n2)], ncol=n2)
	    dimnames(tsurv) <- list(NULL, dimnames(x2)[[1]])
	    }
	else {
	    tsurv <- surv$surv[ntime]
	    tvar  <- surv$varhaz[ntime]
	    }
	if (surv$strata[1] <=1)
	    temp _ list(n=n,time=surv$y[ntime,1],
		     n.risk=surv$y[ntime,2],
		     n.event=surv$y[ntime,3],
		     surv=tsurv,
			type=type)
	else {
	    temp <- surv$strata[1:(1+surv$strata[1])]
	    tstrat <- diff(c(0, temp[-1])) #n in each strata
	    names(tstrat) <- levels(data$strata)
	    temp _ list(n=n, time=surv$y[ntime,1],
		     n.risk=surv$y[ntime,2],
		     n.event=surv$y[ntime,3],
		     surv=tsurv,
		     strata= tstrat, ntimes.strata=tstrat,
			strata.all=strata.all,
			type=type)
	    }
	if (se.fit) temp$std.err <- sqrt(tvar)
	}

    zval _ qnorm(1- (1-conf.int)/2, 0,1)
    if (conf.type=='plain') {
	temp1 <- temp$surv + zval* temp$std * temp$surv
	temp2 <- temp$surv - zval* temp$std * temp$surv
	temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
			conf.type='plain', conf.int=conf.int))
	}
    if (conf.type=='log') {
	xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- ifelse(temp$surv==0, 0*temp$std, exp(log(xx) + zval* temp$std))
	temp2 <- ifelse(temp$surv==0, 0*temp$std, exp(log(xx) - zval* temp$std))
	temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
			conf.type='log', conf.int=conf.int))
	}
    if (conf.type=='log-log') {
	who <- (temp$surv==0 | temp$surv==1) #special cases
	xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- exp(-exp(log(-log(xx)) + zval*temp$std/log(xx)))
	temp1 <- ifelse(who, temp$surv + 0*temp$std, temp1)
	temp2 <- exp(-exp(log(-log(xx)) - zval*temp$std/log(xx)))
	temp2 <- ifelse(who, temp$surv + 0*temp$std, temp2)
	temp <- c(temp, list(upper=temp1, lower=temp2,
			conf.type='log-log', conf.int=conf.int))
	}

    temp$call <- call
    class(temp) <- c('survfit.cox', 'survfit')
    temp
    }



