#
# SCCS @(#)survreg.s	5.8 07/10/00
#  The newest version of survreg, that accepts penalties and strata
#

survreg <- function(formula=formula(data), data=parent.frame(),
	weights, subset, na.action, dist='weibull', 
	init=NULL,  scale=0, control=survreg.control(), parms=NULL, 
	model=FALSE, x=FALSE, y=TRUE, robust=FALSE, ...) {

    call <- match.call()
    m <- match.call(expand=FALSE)
    temp <- c("", "formula", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]
    m[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    Terms <- if(missing(data)) terms(formula, special)
             else              terms(formula, special, data=data)
    m$formula <- Terms
    m <- eval(m, parent.frame())
    ### I commented this out last time -- don't know why
    ###Terms <- attr(m, 'terms')

    weights <- model.extract(m, 'weights')
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")

    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
        if (missing(robust)) robust <- TRUE
        tempc <- untangle.specials(Terms, 'cluster', 1:10)
        ord <- attr(Terms, 'order')[tempc$terms]
        if (any(ord>1)) stop ("Cluster can not be used in an interaction")
        cluster <- strata(m[,tempc$vars], shortlabel=TRUE)  #allow multiples
        dropx <- tempc$terms
        }
    if (length(strats)) {
        temp <- untangle.specials(Terms, 'strata', 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
        strata <- as.numeric(strata.keep)
	nstrata <- max(strata)
        }
    else {
	nstrata <- 1
	strata <- 0
	}

    if (length(dropx)) newTerms<-Terms[-dropx]
    else               newTerms<-Terms
    X<-model.matrix(newTerms,m)
    
    n <- nrow(X)
    nvar <- ncol(X)

    offset<- attr(Terms, "offset")
    if (!is.null(offset)) offset <- as.numeric(m[[offset]])
    else                  offset <- rep(0, n)

    if (is.character(dist)) {
	dlist <- survreg.distributions[[dist]]
	if (is.null(dlist)) stop(paste(dist, ": distribution not found"))
	}
    else if (is.list(dist)) dlist <- dist
    else stop("Invalid distribution object")
    if (is.null(dlist$dist)) {
	if (is.character(dlist$name) && is.function(dlist$init) &&
	    is.function(dlist$deviance)) {}
	else stop("Invalid distribution object")
	}
    else {
	if (!is.character(dlist$name) || is.null(dlist$dist) ||
	    !is.function(dlist$trans) || !is.function(dlist$dtrans))
		stop("Invalid distribution object")
	}	

    type <- attr(Y, "type")
    if (type== 'counting') stop ("Invalid survival type")
    
    logcorrect <- 0   #correction to the loglik due to transformations
    if (!is.null(dlist$trans)) {
	tranfun <- dlist$trans
	exactsurv <- Y[,ncol(Y)] ==1
	if (any(exactsurv)) logcorrect <-sum(log(dlist$dtrans(Y[exactsurv,1])))

	if (type=='interval') {
	    if (any(Y[,3]==3))
		    Y <- cbind(tranfun(Y[,1:2]), Y[,3])
	    else Y <- cbind(tranfun(Y[,1]), Y[,3])
	    }
	else if (type=='left')
	     Y <- cbind(tranfun(Y[,1]), 2-Y[,2])
	else     Y <- cbind(tranfun(Y[,1]), Y[,2])
	if (!all(is.finite(Y))) 
	    stop("Invalid survival times for this distribution")
	}
    else {
	if (type=='left') Y[,2] <- 2- Y[,2]
	else if (type=='interval' && all(Y[,3]<3)) Y <- Y[,c(1,3)]
	}

    if (is.null(dlist$itrans)) itrans <- function(x) x
    else itrans <- dlist$itrans

    if (!is.null(dlist$scale)) {
	if (!missing(scale)) warning(paste(dlist$name, 
			   "has a fixed scale, user specified value ignored"))
	scale <- dlist$scale
	}
    if (!is.null(dlist$dist)){
        if (is.atomic(dlist$dist))
            dlist <- survreg.distributions[[dlist$dist]]
        else
            dlist<-dlist$dist #<TSL>
    }
    if (missing(control)) control <- survreg.control(...)

    if (scale < 0) stop("Invalid scale value")
    if (scale >0 && nstrata >1) 
	    stop("Cannot have multiple strata with a fixed scale")

    # Check for penalized terms
    pterms <- sapply(m, inherits, 'coxph.penalty')
    if (any(pterms)) {
	pattr <- lapply(m[pterms], attributes)
	# 
	# the 'order' attribute has the same components as 'term.labels'
	#   pterms always has 1 more (response), sometimes 2 (offset)
	# drop the extra parts from pterms
	temp <- c(attr(Terms, 'response'), attr(Terms, 'offset'))
	if (length(dropx)) temp <- c(temp, dropx+1)
	pterms <- pterms[-temp]
	temp <- match((names(pterms))[pterms], attr(Terms, 'term.labels'))
	ord <- attr(Terms, 'order')[temp]
	if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')
	##pcols <- (attr(X, 'assign')[-1])[pterms]
        assign<-attrassign(X,newTerms)
        pcols<-assign[-1][pterms]
  
        fit <- survpenal.fit(X, Y, weights, offset, init=init,
				controlvals = control,
			        dist= dlist, scale=scale,
			        strata=strata, nstrat=nstrata,
				pcols, pattr,assign, parms=parms)
	}
    else fit <- survreg.fit(X, Y, weights, offset, 
			    init=init, controlvals=control,
			    dist= dlist, scale=scale, nstrat=nstrata, 
			    strata, parms=parms)

    if (is.character(fit))  fit <- list(fail=fit)  #error message
    else {
	if (scale==0) {
	    nvar <- length(fit$coef) - nstrata
	    fit$scale <- exp(fit$coef[-(1:nvar)])
	    if (nstrata==1) names(fit$scale) <- NULL
	    else names(fit$scale) <- levels(strata.keep)
	    fit$coefficients  <- fit$coefficients[1:nvar]
	    fit$idf  <- 1 + nstrata
	    }
	else {
	    fit$scale <- scale
	    fit$idf  <- 1
	    }
	fit$loglik <- fit$loglik + logcorrect
	}

    
    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    fit$df.residual <- n - sum(fit$df)
#   fit$fitted.values <- itrans(fit$linear.predictors)
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$means <- apply(X,2, mean)
    fit$call <- call
    fit$dist <- dist
    fit$df.resid<n-sum(fit$df) ##used for anova.survreg
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y)     fit$y <- Y
    if (length(parms)) fit$parms <- parms
    if (any(pterms)) class(fit)<- c('survreg.penal', 'survreg')
    else	     class(fit) <- 'survreg'

    if (robust){
        fit$naive.var<-fit$var
        if (length(cluster))
            fit$var<-crossprod(rowsum(resid(fit,"dfbeta"), cluster))
        else
            fit$var<-crossprod(rowsum(resid(fit,"dfbeta")))
    }
    fit
}
