# 
#  SCCS @(#)survpenal.fit.s	1.6 02/08/99
# fit a penalized parametric model
#
survpenal.fit<- function(x, y, weights, offset, init, controlvals, dist, 
		       scale=0, nstrat=1, strata, pcols, pattr, assign,
			 parms=NULL) {

    iter.max <- controlvals$iter.max
    outer.max <- controlvals$outer.max
    eps <- controlvals$rel.tol
    toler.chol <- controlvals$toler.chol
    debug <- controlvals$debug

    if (!is.matrix(x)) stop("Invalid X matrix ")
    n <- nrow(x)
    nvar <- ncol(x) 
    ny <- ncol(y)
    if (is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights)) weights<- rep(1.0,n)
    else if (any(weights<=0)) stop("Invalid weights, must be >0")

    if (scale <0) stop("Invalid scale")
    if (scale >0 && nstrat >1) 
	    stop("Cannot have both a fixed scale and strata")
    if (nstrat>1 && (missing(strata) || length(strata)!= n))
	    stop("Invalid strata variable")
    if (nstrat==1) strata <- rep(1,n)
    if (scale >0)
      nstrat2 <- 0
    else
      nstrat2 <- nstrat

    if (is.character(dist)) {
	sd <- survreg.distributions[[dist]]
	if (is.null(sd)) stop ("Unrecognized distribution")
	}
    else sd <- dist
    dnum <- match(sd$name, c("Extreme value", "Logistic", "Gaussian"))
    fdensity<-function(z) {
      stop("this is just a placeholder and should never be called")
    }
    f.expr1<-function(coef) {
      stop("this is just a placeholder and should never be called")
    }
    f.expr2<-function(coef) {
      stop("this is just a placeholder and should never be called")
    }
    if (is.na(dnum)) {
	# Not one of the "built-in distributions
	dnum <- 4
	fitter <- c('survreg3', 'survreg5')
	#Set up the callback for the sparse frailty term
	n2 <- n + sum(y[,ny]==3)
	fdensity <- function(z){
	    if (length(parms)) temp <- sd$density(z, parms)
	    else               temp <- sd$density(z)
	    
	    if (!is.matrix(temp) || any(dim(temp) != c(n2,5)))
		    stop("Density function returned an invalid matrix")
	    survlist$density <- as.vector(as.double(temp))
	    survlist
          }
	survlist <- list(z=double(n2), density=double(n2*5))
	###.C("init_survcall", as.integer(sys.nframe()), expr1)
	}
    else fitter <- c('survreg2', 'survreg4')

    # This is a subset of residuals.survreg: define the first and second
    #   derivatives at z=0 for the 4 censoring types
    #   Used below for starting estimates
    derfun <- function(y, eta, sigma, density, parms) {
	ny <<- ncol(y)
	status <- y[,ny]
	z <- (y[,1] - eta)/sigma
	dmat <- density(z,parms)
	dtemp<- dmat[,3] * dmat[,4]    #f'
	if (any(status==3)) {
	    z2 <- (y[,2] - eta)/sigma
	    dmat2 <- density(z2)
	    }
	else {
	    dmat2 <- matrix(0,1,5)   #dummy values
	    z2 <- 0
	    }
	tdenom <- ((status==0) * dmat[,2]) +
		  ((status==1) * 1 )       +
		  ((status==2) * dmat[,1]) +
		  ((status==3) * ifelse(z>0, dmat[,2]-dmat2[,2], 
		                             dmat2[,1] - dmat[,1]))
	tdenom <- 1/(tdenom* sigma)
	dg <- -tdenom   *(((status==0) * (0-dmat[,3])) +
			  ((status==1) * dmat[,4]) + 
			  ((status==2) * dmat[,3]) +
			  ((status==3) * (dmat2[,3]- dmat[,3])))

	ddg <- (tdenom/sigma)*(((status==0) * (0- dtemp)) +
			       ((status==1) * dmat[,5]) +
			       ((status==2) * dtemp) +
			       ((status==3) * (dmat2[,3]*dmat2[,4] - dtemp))) 
	list(dg = dg, ddg = ddg - dg^2)
	}
    status <- y[,ny]

    #
    # are there any sparse frailty terms?
    # 
    npenal <- length(pattr)
    if (npenal == 0 || length(pcols) != npenal)
	    stop("Invalid pcols or pattr arg")
    sparse <- sapply(pattr, function(x) !is.null(x$sparse) &&  x$sparse)
    if (sum(sparse) >1) stop("Only one sparse penalty term allowed")

    #
    # Create a marking vector for the terms, the same length as assign
    #    with pterms == 0=ordinary term, 1=penalized, 2=sparse,
    #    pindex = length of pcols = position in pterms
    # 
    # Make sure that pcols is a strict subset of assign, so that the
    #   df computation (and printing) can unambiguously decide which cols of
    #   X are penalized and which are not when doing "terms" like actions.
    # To make some downstream things easier, order pcols and pattr to be
    #   in the same relative order as the terms in 'assign' 
    #
    ## if (missing(assign)) assign <- attr(x, 'assign')
    pterms <- rep(0, length(assign))
    names(pterms) <- names(assign)
    pindex <- rep(0, npenal)
    for (i in 1:npenal) {
	temp <- unlist(lapply(assign, function(x,y) (length(x) == length(y) &&
					     all(x==y)), pcols[[i]]))
	if (sparse[i]) pterms[temp] <- 2
	else pterms[temp] <- 1
	pindex[i] <- (seq(along=temp))[temp]
	}
    if ((sum(pterms==2) != sum(sparse)) || (sum(pterms>0) != npenal))
	    stop("pcols and assign arguments disagree")
    if (any(pindex != sort(pindex))) {
	temp <- order(pindex)
	pindex <- pindex[temp]
	pcols <- pcols[temp]
	pattr <- pattr[temp]
	}
    
    # ptype= 1 or 3 if a sparse term exists, 2 or 3 if a non-sparse exists
    ptype <- any(sparse) + 2*(any(!sparse))

    if (any(sparse)) {
	sparse.attr <- (pattr[sparse])[[1]]  #can't use [[sparse]] directly
	                                     # if 'sparse' is a T/F vector
	fcol <- unlist(pcols[sparse])
	if (length(fcol) > 1) stop("Sparse term must be single column")

	# Remove the sparse term from the X matrix
	frailx <- x[, fcol]
	x <- x[, -fcol, drop=F]
	for (i in 1:length(assign)){
	    j <- assign[[i]]
	    if (j[1] > fcol) assign[[i]] <- j-1
	    }
	for (i in 1:npenal) {
	    j <- pcols[[i]]
	    if (j[1] > fcol) pcol[[i]] <- j-1
	    }

	frailx <- match(frailx, sort(unique(frailx)))
	nfrail <- max(frailx)
	nvar <- nvar - 1

	#Set up the callback for the sparse frailty term
	pfun1 <- sparse.attr$pfun
	f.expr1 <- function(coef){
            coxlist1$coef<-coef
	    if (is.null(extra1)) temp <- pfun1(coef, theta1, n.eff)
	    else  temp <- pfun1(coef, theta1, n.eff, extra1)

	    if (!is.null(temp$recenter)) 
		    coxlist1$coef <- coxlist1$coef - as.double(temp$recenter)
	    if (!temp$flag) {
		coxlist1$first <- -as.double(temp$first)
		coxlist1$second <- as.double(temp$second)
	        }
	    coxlist1$penalty <- -as.double(temp$penalty)
	    coxlist1$flag   <- as.logical(temp$flag)
	    if (any(sapply(coxlist1, length) != c(rep(nfrail,3), 1, 1)))
		    stop("Incorrect length in coxlist1")
	    coxlist1
          }
	coxlist1 <- list(coef=double(nfrail), first=double(nfrail), 
			 second=double(nfrail), penalty=0.0, flag=F)
	###.C("init_coxcall1", as.integer(sys.nframe()), expr1)
	}
    else {
	frailx <- 0
	nfrail <- 0
	}
    nvar2 <- nvar + nstrat2

    # Now the non-sparse penalties
    if (sum(!sparse) >0) {
	full.imat <- !all(unlist(lapply(pattr, function(x) x$diag)))
	ipenal <- (1:length(pattr))[!sparse]   #index for non-sparse terms
	f.expr2 <- function(coef){
            coxlist2$coef<-coef
	    pentot <- 0
	    for (i in ipenal) {
		pen.col <- pcols[[i]]
		coef<-coxlist2$coef[pen.col]
		if (is.null(extralist[[i]]))
			temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]],n.eff)
		else    temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]],
						n.eff,extralist[[i]])
		if (!is.null(temp$recenter))
		    coxlist2$coef[pen.col] <- coxlist2$coef[pen.col]- 
			                               temp$recenter
		if (temp$flag) coxlist2$flag[pen.col] <- TRUE
		else {
		    coxlist2$flag[pen.col] <- FALSE
		    coxlist2$first[pen.col] <- -temp$first
		    if (full.imat) {
			tmat <- matrix(coxlist2$second, nvar2, nvar2)
			tmat[pen.col,pen.col] <- temp$second
			coxlist2$second <- c(tmat)
		        }
		    else coxlist2$second[pen.col] <- temp$second
		    }
		pentot <- pentot - temp$penalty
	        }
	    coxlist2$penalty <- as.double(pentot)
	    if (any(sapply(coxlist2, length) != length2)) 
		    stop("Length error in coxlist2")
	    coxlist2
          }
        ##if (debug) debug(f.expr2)
	if (full.imat) {
	    coxlist2 <- list(coef=double(nvar2), first=double(nvar2), 
		    second= double(nvar2*nvar2), penalty=0.0, flag=rep(F,nvar2))
	    length2 <- c(nvar2, nvar2, nvar2*nvar2, 1, nvar2)
	    }  
	else {
	    coxlist2 <- list(coef=double(nvar2), first=double(nvar2),
		    second=double(nvar2), penalty= 0.0, flag=rep(F,nvar2))
	    length2 <- c(nvar2, nvar2, nvar2, 1, nvar2)
	    }
	###.C("init_coxcall2", as.integer(sys.nframe()), expr2)
        }
    else full.imat <- FALSE

    #
    # "Unpack" the passed in paramter list, 
    #   and make the initial call to each of the external routines
    #
    cfun <- lapply(pattr, function(x) x$cfun)
    parmlist <- lapply(pattr, function(x,eps) c(x$cparm, eps2=eps), sqrt(eps))
    extralist<- lapply(pattr, function(x) x$pparm)
    iterlist <- vector('list', length(cfun))
    thetalist <- vector('list', length(cfun))
    printfun  <- lapply(pattr, function(x) x$printfun)
    for (i in 1:length(cfun)) {
	temp <- (cfun[[i]])(parmlist[[i]], iter=0)
	if (sparse[i]) {
	    theta1 <- temp$theta
	    extra1 <- extralist[[i]]
	    }
	thetalist[[i]] <- temp$theta
	iterlist[[i]] <- temp
	}

    #
    # Manufacture the list of calls to cfun, with appropriate arguments
    #
    temp1 <- c('x', 'coef', 'plik', 'loglik', 'status', 'neff',  'df', 'trH')
    temp2 <- c('frailx', 'fcoef', 'loglik',  'fit$loglik', 'status', 'n.eff')
    temp3 <- c('x[,pen.col]', 'coef[pen.col]','loglik',
	       'fit$loglik', 'status', 'n.eff')
    calls <- vector('expression', length(cfun))
    cargs <- lapply(pattr, function(x) x$cargs)
    for (i in 1:length(cfun)) {
	tempchar <- paste("(cfun[[", i, "]])(parmlist[[", i, "]], iter,",
			  "iterlist[[", i, "]]")
	temp2b <- c(temp2, paste('pdf[', i, ']'), paste('trH[', i, ']'))
	temp3b <- c(temp3, paste('pdf[', i, ']'), paste('trH[', i, ']'))
	if (length(cargs[[i]])==0) 
	    calls[i] <- parse(text=paste(tempchar, ")"))
	else {
	    temp <- match(cargs[[i]], temp1)
	    if (any(is.na(temp))) stop(paste((cargs[[i]])[is.na(temp)],
					    "not matched"))
	    if (sparse[i]) temp4 <- paste(temp2b[temp], collapse=',')
	    else           temp4 <- paste(temp3b[temp], collapse=',')
	    
	    calls[i] <- parse(text=paste(paste(tempchar,temp4,sep=','),')'))
	    }
        }
    need.df <- any(!is.na(match(c('df', 'trH'), unlist(cargs))))#do any use df?

    #
    # Last of the setup: create the vector of variable names
    #
    varnames <- dimnames(x)[[2]]
    for (i in 1:npenal) {
	if (!is.null(pattr[[i]]$varname))
		varnames[pcols[[i]]] <- pattr[[i]]$varname
        }

    #
    # Fit the model with just a mean and scale
    #    assume initial values and penalties don't apply here
    #
    meanonly <- (nvar==1 && all(x==1) && nfrail==0)
    if (meanonly) stop("Cannot fit a penalized 'mean only' model")

    yy <- ifelse(status !=3, y[,1], (y[,1]+y[,2])/2 )
    coef <- sd$init(yy, weights,parms)
    # We sometimes get into trouble with a small estimate of sigma,
    #  (the surface isn't SPD), but never with a large one.  Double it.
    if (scale >0) vars <- log(scale)
    else vars <- log(coef[2])/2 +.7  # init returns \sigma^2, I need log(sigma)
    coef <- c(coef[1], rep(vars, nstrat))
    # get a better initial value for the mean using the "glim" trick
    deriv <- derfun(y, yy, exp(vars), sd$density, parms)
    wt <-  -1*deriv$ddg*weights
    coef[1] <- sum(weights*deriv$dg + wt*(yy -offset)) / sum(wt)

    # Now the fit proper (intercept only)
    temp <- 1 +nstrat2
    fit0 <- .C(fitter[1],
	       iter = as.integer(iter.max),
	       as.integer(n),
	       as.integer(1),
	       as.double(y),
	       as.integer(ny),
	       rep(1.0, n),
	       as.double(weights),
	       as.double(offset),
	       coef= as.double(coef),
	       as.integer(nstrat2),
	       as.integer(strata),
	       u = double(3*(temp) + temp^2),
	       var = matrix(0.0, temp, temp),
	       loglik=double(1),
	       flag=integer(1),
	       as.double(eps),
	       as.double(toler.chol), 
	       as.integer(dnum),
	       debug = as.integer(floor(debug/2)), fdensity, environment(),
                    PACKAGE="survival")

    # The "effective n" of the model
    temp <-  mean(exp(fit0$coef[-1]))   #overall sd
    n.eff <- sd$var(temp^2) * (solve(fit0$var))[1,1]

    #
    # Fit the model with all covariates
    #   Start with initial values
    #
    nvar3 <- nvar2 + nfrail
    if (is.numeric(init)) {
	if (length(init) != nvar3) {
	    if (length(init) == nvar2) init <- c(rep(0,nfrail), init)
	    else stop("Wrong length for inital values")
	    }
	if (scale >0) init <- c(init, log(scale))
	}
    else  {
	# The algebra behind the 'glim' trick just doesn't work here
	#  Use the intercept fit + zeros
	#    coef order = frailty, intercept, other covariates, sigmas
	init <- c(rep(0, nfrail), fit0$coef[1], rep(0, nvar-1), fit0$coef[-1])
	}

    #
    # Tack on the sigmas to "assign", so that the df component includes
    #   the sigmas
    if (nstrat2 >0) assign <- c(assign, list(sigma=(1+nvar):nvar2))

    iter2 <- 0
    iterfail <- NULL
    thetasave <- unlist(thetalist)
    for (iterx in 1:outer.max) {
	fit <- .C(fitter[2],
		   iter = as.integer(iter.max),
		   as.integer(n),
		   as.integer(nvar),
		   as.double(y),
		   as.integer(ny),
		   as.double(x),
	           as.double(weights),
		   as.double(offset),
		   coef= as.double(init),
	           as.integer(nstrat2),
	           as.integer(strata),
		   u = double(3*(nvar3) + nvar2*nvar3),
		   hmat = double(nvar2*nvar3),
	           hinv = double(nvar2*nvar3),
		   loglik=double(1),
		   flag=integer(1),
		   as.double(eps),
	           as.double(toler.chol), 
		   as.integer(dnum),
	           debug = as.integer(debug),
	           as.integer(ptype),
		   as.integer(full.imat),
		   as.integer(nfrail),
		   as.integer(frailx),
	           fdiag = double(nvar3),f.expr1,f.expr2,fdensity, environment(),
                    PACKAGE="survival")

	if (debug>0) browser()
	iter <- iterx
	iter2 <- iter2 + fit$iter
	if (fit$iter >=iter.max) iterfail <- c(iterfail, iter)

	if (nfrail >0) {
	    fcoef <- fit$coef[1:nfrail]
	    coef  <- fit$coef[nfrail + 1:nvar2]
	    }
	else coef <- fit$coef[1:nvar2]

	# If any penalties were infinite, the C code has made fdiag=1 out
	#  of self-preservation (0 divides).  But such coefs are guarranteed
	#  zero so the variance should be too.)
	temp <- rep(F, nvar2+nfrail)
	if (nfrail>0) temp[1:nfrail] <- coxlist1$flag
	if (ptype >1) temp[nfrail+ 1:nvar2] <- coxlist2$flag
	fdiag <- ifelse(temp, 0, fit$fdiag)

	if (need.df) {
            #get the penalty portion of the second derive matrix
	    if (nfrail>0) temp1 <- coxlist1$second
	    else 	  temp1 <- 0
	    if (ptype>1)  temp2 <- coxlist2$second
	    else          temp2 <- 0
					
	    dftemp <-coxpenal.df(matrix(fit$hmat, ncol=nvar2),  
			         matrix(fit$hinv, ncol=nvar2), fdiag, 
				 assign, ptype, nvar2,
		                 temp1, temp2, pindex[sparse])
	    df <- dftemp$df
	    var  <- dftemp$var
	    var2 <- dftemp$var2
	    pdf <- df[pterms>0]	          # df's for penalized terms
	    trH <- dftemp$trH[pterms>0]   # trace H 
	    }

	if (nfrail >0)  penalty <- -coxlist1$penalty
	else            penalty <- 0
	if (ptype >1) penalty <- penalty - coxlist2$penalty
	loglik <- fit$loglik + penalty  #C code returns PL - penalty
	if (iter==1) penalty0 <- penalty

	#
	# Call the control function(s)
	#
	done <- TRUE
	for (i in 1:length(cfun)) {
	    pen.col <- pcols[[i]]
	    temp <- eval(calls[i])
	    if (sparse[i]) theta1 <- temp$theta
	    thetalist[[i]] <- temp$theta
	    iterlist[[i]] <- temp
	    done <- done & temp$done
    	    }
	if (done) break

	# 
	# Choose starting estimates for the next iteration
	#
	if (iter==1) {
	    init <- coefsave <- fit$coef
	    thetasave <- cbind(thetasave, unlist(thetalist))
	    }
	else {
	    temp <- unlist(thetalist)
	    coefsave <- cbind(coefsave, fit$coef)
	    # temp = next guess for theta
	    # *save = prior thetas and the resultant fits
	    # choose as initial values the result for the closest old theta
	    howclose <- apply((thetasave-temp)^2,2, sum)
	    which <- min((1:iter)[howclose==min(howclose)])
	    init <- coefsave[,which]
	    thetasave <- cbind(thetasave, temp)
	    }
        }   #end of the iteration loop

    if (!need.df) {  #didn't need it iteration by iteration, but do it now
        #get the penalty portion of the second derive matrix
	if (nfrail>0) temp1 <- coxlist1$second
	else 	      temp1 <- 0
	if (ptype>1)  temp2 <- coxlist2$second
	else          temp2 <- 0
					
	dftemp <-coxpenal.df(matrix(fit$hmat,ncol=nvar2),  
			     matrix(fit$hinv,ncol=nvar2),  fdiag, 
		             assign, ptype, nvar2, 
		             temp1, temp2, pindex[sparse])
	df <- dftemp$df
	trH <- dftemp$trH
	var <- dftemp$var
	var2  <- dftemp$var2
        }

    if (iter.max >1 && length(iterfail)>0)
	    warning(paste("Inner loop failed to coverge for iterations", 
			  paste(iterfail, collapse=' ')))
    which.sing <- (fdiag[nfrail + 1:nvar] ==0)
    coef[which.sing] <- NA

    names(iterlist) <- names(pterms[pterms>0])
    cname <- varnames
    cname <- c(cname, rep("Log(scale)", nstrat2))
    dimnames(var) <- list(cname, cname)
    names(coef) <- cname

    if (nfrail >0) {
	lp <- offset + fcoef[frailx]
	lp <- lp + x %*%coef[1:nvar] 
	list(coefficients  = coef,
	     icoef = fit0$coef,
	     var    = var,
	     var2   = var2,
	     loglik = c(fit0$loglik, loglik),
	     iter   = c(iter, iter2),
	     linear.predictors = as.vector(lp),
	     frail = fcoef,
	     fvar  = dftemp$fvar,
	     df = df, 
	     penalty= c(penalty0, penalty),
	     pterms = pterms, assign2=assign,
	     history= iterlist,
	     printfun=printfun)
	}
    else {  #no sparse terms
	list(coefficients  = coef,
	     icoef = fit0$coef,
	     var    = var,
	     var2   = var2,
	     loglik = c(fit0$loglik, loglik),
	     iter   = c(iter, iter2),
	     linear.predictors = as.vector(x%*%coef[1:nvar]),
	     df = df, df2=dftemp$df2,
	     penalty= c(penalty0, penalty), 
	     pterms = pterms, assign2=assign,
	     history= iterlist,
	     printfun= printfun)
	}
    }
