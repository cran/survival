# SCCS @(#)residuals.survreg.s	4.12 02/06/99
# 
#  Residuals for survreg objects
residuals.survreg <- function(object, type=c('response', 'deviance',
		      'dfbeta', 'dfbetas', 'working', 'ldcase',
		      'ldresp', 'ldshape', 'matrix'), 
		      rsigma =TRUE, collapse=FALSE, weighted=FALSE, ...) {
    type <-match.arg(type)
    n <- length(object$linear.predictors)
    Terms <- object$terms
    if(!inherits(Terms, "terms"))
	    stop("invalid terms component of  object")

    strata <- attr(Terms, 'specials')$strata
    coef <- object$coefficients
    intercept <- attr(Terms, "intercept")
    response  <- attr(Terms, "response")
    weights <- object$weights
    if (is.null(weights)) weighted <- FALSE

    #
    # What do I need to do the computations?
    #
    if (type=='response' || type=='deviance') need.x <-FALSE
    else need.x <- TRUE

    if (type=='ldshape' || type=='ldcase') need.strata <- TRUE
    else need.strata <- FALSE
    
    need.y <- TRUE

    # grab what I need
    if (need.strata || (need.y && is.null(object$y)) || 
	               (need.x && is.null(object$x)) ) {
	# I need the model frame
	if (is.null(object$model)) m <- model.frame(object)
	else m <- object$model
	}

    
    if (need.strata && !is.null(strata)) {
	temp <- untangle.specials(Terms, 'strata', 1)
	Terms2 <- Terms[-temp$terms]
	if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
	strata <- as.numeric(strata.keep)
	nstrata <- max(strata)
	sigma <- object$scale[strata]
	}
    else {
	strata <- rep(1,n)
	nstrata <- 1
	sigma <- object$scale
	Terms2 <- Terms
	}

    if (need.x){
	x <- object$x
	if (is.null(x)) x <- model.matrix(Terms2, m)
	}
	
    if (need.y) {
	y <- object$y
	if (is.null(y)) y <- model.extract(m, 'response')
	status <- y[,ncol(y)]
	}
    #
    # Grab the distribution
    #
    if (is.character(object$dist)) 
	    dd <- survreg.distributions[[object$dist]]
    else dd <- object$dist
    if (is.null(dd$itrans)) {
	itrans <- dtrans <-function(x)x
	}
    else {
	itrans <- dd$itrans
	dtrans <- dd$dtrans
	}
    if (!is.null(dd$dist))  dd <- survreg.distributions[[dd$dist]]
    deviance <- dd$deviance
    dens <- dd$density

    nvar <- length(object$coef)
    if (rsigma) vv <- object$var
    else        vv <- object$var[1:nvar, 1:nvar]

    # Create the matrix of derivatives
    #  The "density" function returns F, 1-F, f, f'/f, and f''/f
    if (type != 'response') {
	status <- y[,ncol(y)]
	eta <- object$linear.predictor
	z <- (y[,1] - eta)/sigma
	dmat <- dens(z, object$parms)
	dtemp<- dmat[,3] * dmat[,4]    #f'
	if (any(status==3)) {
	    z2 <- (y[,2] - eta)/sigma
	    dmat2 <- dens(z2, object$parms)
	    }
	else {
	    dmat2 <- dmat   #dummy values
	    z2 <- 0
	    }
	deriv <- matrix(n,6)
	tdenom <- ((status==0) * dmat[,2]) +
		  ((status==1) * 1 )       +
		  ((status==2) * dmat[,1]) +
		  ((status==3) * ifelse(z>0, dmat[,2]-dmat2[,2], 
		                             dmat2[,1] - dmat[,1]))
	g <- log(ifelse(status==1, dmat[,3]/sigma, tdenom))
	tdenom <- 1/(tdenom* sigma)
	dg <- -tdenom   *(((status==0) * (0-dmat[,3])) +
			  ((status==1) * dmat[,4]) + 
			  ((status==2) * dmat[,3]) +
			  ((status==3) * (dmat2[,3]- dmat[,3])))

	ddg <- (tdenom/sigma)*(((status==0) * (0- dtemp)) +
			       ((status==1) * dmat[,5]) +
			       ((status==2) * dtemp) +
			       ((status==3) * (dmat2[,3]*dmat2[,4] - dtemp))) 
	td2 <- tdenom * sigma
	ds  <- ifelse(status<3, dg * sigma * z,
		                td2*(z2*dmat2[,3] - z*dmat[,3]))
	dds <- ifelse(status<3, ddg* (sigma*z)^2,
		                td2*(z2*z2*dmat2[,3]*dmat2[,4] -
				     z * z*dmat[,3] * dmat[,4]))
	dsg <- ifelse(status<3, ddg* sigma*z,
		      td2 *(z2*dmat2[,3]*dmat2[,4] - z*dtemp))
	deriv <- cbind(g, dg, ddg=ddg- dg^2, 
		       ds = ifelse(status==1, ds-1, ds), 
		       dds=dds - ds*(1+ds), 
		       dsg=dsg - dg*(1+ds))
	}

    #
    # Now, lay out the code one case at a time.
    #  There is some repetition this way, but otherwise the code gets
    #    too complicated.
    #
    if (type=='response') {
	yhat0 <- deviance(y, object$scale[strata], object$parms)
	rr <-  itrans(yhat0$center) - itrans(object$linear.predictor)
	}

    else if (type=='deviance') {
	yhat0 <- deviance(y, object$scale[strata], object$parms)
	rr <- (-1)*deriv[,2]/deriv[,3]  #working residuals
	rr <- sign(rr)* sqrt(2*(yhat0$loglik - deriv[,1]))
	}
    
    else if (type=='dfbeta' || type== 'dfbetas') {
	score <- deriv[,2] * x  # score residuals
	if (rsigma) score <- cbind(score, deriv[,4])
	rr <-    score %*% vv
	if (type=='dfbetas') rr <- rr %*% diag(1/sqrt(diag(vv)))
	}

    else if (type=='working') rr <- (-1)*deriv[,2]/deriv[,3]

    else if (type=='ldcase'){
	score <- deriv[,2] * x 
	if (rsigma) score <- cbind(score, deriv[,4])
	dfbeta<- score %*% vv
	rr <- apply(dfbeta*score,1,sum)
	}

    else if (type=='ldresp') {
	rscore <-  deriv[,3] *  (x * sigma)
	if (rsigma) rscore <- cbind(rscore, deriv[,6] * sigma)
	temp <-  rscore %*% vv
	rr <- apply(rscore * temp, 1 , sum)
	}

    else if (type=='ldshape') {
	sscore <- deriv[,6] *x
	if (rsigma) sscore <- cbind(sscore,  deriv[,5]) 
	temp <- sscore %*% vv
	rr <- apply(sscore * temp, 1, sum)
	}

   else {  #type = matrix
	rr <- deriv
	}

    #
    # Multiply up by case weights, if requested
    #
    if (weighted) rr <- rr * weights

    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	rr <- naresid(object$na.action, rr)
	if (is.matrix(rr)) n <- nrow(rr)
	else               n <- length(rr)
	}

    # Collapse if desired
    if (!missing(collapse)) {
	if (length(collapse) !=n) stop("Wrong length for 'collapse'")
	rr <- rowsum(rr, collapse)
	}

    rr
    }

	







