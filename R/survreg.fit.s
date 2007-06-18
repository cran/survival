# 
#  SCCS @(#)survreg.fit.s	5.10 07/10/00
#
survreg.fit<- function(x, y, weights, offset, init, controlvals, dist, 
		       scale=0, nstrat=1, strata, parms=NULL) {

    controlvals<-do.call("survreg.control", controlvals)
    iter.max <- controlvals$iter.max
    eps <- controlvals$rel.tolerance
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
    if (!is.function(sd$density)) 
	stop("Missing density function in the definition of the distribution")
    dnum <- match(sd$name, c("Extreme value", "Logistic", "Gaussian"))
    if (is.na(dnum)) {
	# Not one of the "built-in distributions
	dnum <- 4
	fitter <- 'survreg3'
	#Set up the callback for the sparse frailty term
	n2 <- n + sum(y[,ny]==3)
	f.expr1 <- function(z){
            
	    if (length(parms)) temp <- sd$density(z, parms)
	    else               temp <- sd$density(z)
	    
	    if (!is.matrix(temp) || any(dim(temp) != c(n2,5)))
		    stop("Density function returned an invalid matrix")
	    list(z=z,density=as.vector(as.double(temp)))}
        
	survlist <- list(z=double(n2), density=double(n2*5))
	###.C("init_survcall", expr1, PACKAGE="survival")
	}
    else {
        fitter <- 'survreg2'
        f.expr1<-function(z) NULL
    }
    ##
    ## environment for callbacks
    rho<-environment()
    # This is a subset of residuals.survreg: define the first and second
    #   derivatives at z=0 for the 4 censoring types
    #   Used below for starting estimates
    derfun <- function(y, eta, sigma, density, parms) {
	ny <- ncol(y)
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

    #
    # Fit the model with just a mean and scale
    #    assume initial values don't apply here
    # Unless, of course, someone is fitting a mean only model!
    #
    meanonly <- (nvar==1 && all(x==1))
    if (!meanonly) {
	yy <- ifelse(y[,ny]!=3, y[,1], (y[,1]+y[,2])/2 )
	coef <- sd$init(yy, weights, parms)
	#init returns \sigma^2, I need log(sigma)
	# We sometimes get into trouble with a small estimate of sigma,
	#  (the surface isn't SPD), but never with a large one.  Double it.
	if (scale >0) vars <- log(scale)
	else vars <- log(coef[2])/2  +.7
	coef <- c(coef[1], rep(vars, nstrat))
	
	# get a better initial value for the mean using the "glim" trick
	deriv <- derfun(y, yy, exp(vars), sd$density, parms)
	wt <-  -1*deriv$ddg*weights
	coef[1] <- sum(weights*deriv$dg + wt*(yy -offset)) / sum(wt)

	# Now the fit proper (intercept only)
	nvar2 <- 1 +nstrat2
	fit0 <- .C(fitter,
		       iter = as.integer(iter.max),
		       as.integer(n),
		       as.integer(1),
		       as.double(y),
		       as.integer(ny),
		       as.double(rep(1.0, n)),
		       as.double(weights),
		       as.double(offset),
		       coef= as.double(coef),
		       as.integer(nstrat2),
		       as.integer(strata),
		       u = double(3*(nvar2) + nvar2^2),
		       var = matrix(0.0, nvar2, nvar2),
		       loglik=double(1),
		       flag=integer(1),
		       as.double(eps),
		       as.double(toler.chol), 
		       as.integer(dnum),
		       debug = as.integer(floor(debug/2)),
                       Rexpr=f.expr1, Renv=rho,
                   PACKAGE="survival")
	}

    #
    # Fit the model with all covariates
    #
    nvar2 <- nvar + nstrat2
    if (is.numeric(init)) {
	if (length(init) != nvar2) stop("Wrong length for initial parameters")
	if (scale >0) init <- c(init, log(scale))
	}
    else  {
	# Do the 'glim' method of finding an initial value of coef
	if (meanonly) {
	    yy <- ifelse(y[,ny]!=3, y[,1], (y[,1]+y[,2])/2 )
	    coef <- sd$init(yy, weights, parms)
	    if (scale >0) vars <- rep(log(scale), nstrat)
	    else vars  <- rep(log(coef[2])/2 + .7, nstrat)  
	    }
	else vars <- fit0$coef[-1]
	eta <- yy - offset     #what would be true for a 'perfect' model

	deriv <- derfun(y, yy, exp(vars[strata]), sd$density, parms)
	wt <-  -1*deriv$ddg*weights
	coef <- coxph.wtest(t(x)%*% (wt*x), 
		       c((wt*eta + weights*deriv$dg)%*% x),
			    toler.chol=toler.chol)$solve
	init <- c(coef, vars)
	}

    # Now for the fit in earnest

    fit <- .C(fitter,
		   iter = as.integer(iter.max),
		   n = as.integer(n),
		   as.integer(nvar),
		   as.double(y),
		   as.integer(ny),
		   as.double(x),
	           as.double(weights),
		   as.double(offset),
		   coef= as.double(init),
	           as.integer(nstrat2),
	           as.integer(strata),
		   u = double(3*(nvar2) + nvar2^2),
		   var = matrix(0.0, nvar2, nvar2),
		   loglik=double(1),
		   flag=integer(1),
		   as.double(eps),
	           as.double(toler.chol), 
		   as.integer(dnum),
                   debug = as.integer(debug),
              Rexpr=f.expr1, Renv=rho,
              PACKAGE="survival")

    if (debug>0) browser()
    if (iter.max >1 && fit$flag > nvar2) {
	if (controlvals$failure==1)
	       warning("Ran out of iterations and did not converge")
	else if (controlvals$failure==2)
	       return("Ran out of iterations and did not converge")
	}

    cname <- dimnames(x)[[2]]
    if (is.null(cname)) cname <- paste("x", 1:ncol(x))
    if (scale==0) cname <- c(cname, rep("Log(scale)", nstrat))
    dimnames(fit$var) <- list(cname, cname)
    if (scale>0) fit$coef <- fit$coef[1:nvar2]
    names(fit$coef) <- cname

    if (meanonly) {
	coef0 <- fit$coef
	loglik <- rep(fit$loglik,2)
	}
    else {
	coef0 <- fit0$coef
	names(coef0) <- c("Intercept", rep("Log(scale)", nstrat))
	loglik <- c(fit0$loglik, fit$loglik)
	}
    temp <- list(coefficients   = fit$coef,
		 icoef  = coef0, 
		 var    = fit$var,
		 loglik = loglik, 
		 iter   = fit$iter,
		 linear.predictors = c(x %*% fit$coef[1:nvar]),	
		 df     = length(fit$coef)
		 )
    if (debug>0) {
	temp$u <- fit$u[1:nvar2]
	JJ     <- matrix(fit$u[-seq(1, 3*nvar2)], nvar2, nvar2)
	temp$JJ <- JJ
	temp$var2 <- fit$var %*% JJ %*% fit$var
	}

    temp
    }
