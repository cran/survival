#SCCS  @(#)coxph.s	5.12 06/12/00
# Version with general penalized likelihoods

coxph <- function(formula=formula(data), data=parent.frame(),
	weights, subset, na.action, init,
	control, method= c("efron", "breslow", "exact"),
	singular.ok =TRUE, robust=FALSE,
	model=FALSE, x=FALSE, y=TRUE, ...) {

    method <- match.arg(method)
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    temp <- c("", "formula", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]
    special <- c("strata", "cluster")
    Terms <- if(missing(data)) terms(formula, special)
	     else              terms(formula, special, data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    if (NROW(m)==0)
        stop("No (non-missing) observations")
    
    if (missing(control)) control <- coxph.control(...)
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    weights <- model.extract(m, 'weights')
    offset<- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if(tt == 0)
		    rep(0, nrow(Y))
	      else if(tt == 1)
		      m[[offset]]
	      else {
		    ff <- m[[offset[1]]]
		    for(i in 2:tt)
			    ff <- ff + m[[offset[i]]]
		    ff
		    }

    attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
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
	strats <- as.numeric(strata.keep)
	}

    ##if (length(dropx)) X <- model.matrix(Terms[-dropx], m)[,-1,drop=F]
    ##else               X <- model.matrix(Terms, m)[,-1,drop=F]
    ### this is inefficient, but subscripting loses the assign attribute
    if (length(dropx))
        newTerms<-Terms[-dropx]
    else
        newTerms<-Terms
    X<-model.matrix(newTerms,m)
    assign<-lapply(attrassign(X,newTerms)[-1],function(x) x-1)
    X<-X[,-1,drop=FALSE]
    
    
    type <- attr(Y, "type")
    if (type!='right' && type!='counting')
	stop(paste("Cox model doesn't support \"", type,
			  "\" survival data", sep=''))
    if (missing(init)) init <- NULL

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
        pcols<-assign[pterms]
        
	#penalized are hard sometimes	
	if (control$eps.miss)   control$eps <- 1e-7
	if (control$iter.miss)  control$iter.max <- 20  

        fit <- coxpenal.fit(X, Y, strats, offset, init=init,
				control,
				weights=weights, method=method,
				row.names(m), pcols, pattr, assign)
	}
    else {
	if( method=="breslow" || method =="efron") {
	    if (type== 'right')  fitter <- get("coxph.fit")
	    else                 fitter <- get("agreg.fit")
	    }
	else if (method=='exact') fitter <- get("agexact.fit")
	else stop(paste ("Unknown method", method))

	fit <- fitter(X, Y, strats, offset, init, control, weights=weights,
			    method=method, row.names(m))
	}

    if (is.character(fit)) {
	fit <- list(fail=fit)
	class(fit) <- 'coxph'
	}
    else {
	if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
	   vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
	   msg <-paste("X matrix deemed to be singular; variable",
			   paste(vars, collapse=" "))
	   if (singular.ok) warning(msg)
	   else             stop(msg)
	   }
	fit$n <- nrow(Y)
	class(fit) <- fit$method
	fit$terms <- Terms
	##fit$assign <- attr(X, 'assign')
        fit$assign<-assign
	if (robust) {
	    fit$naive.var <- fit$var
	    fit$method    <- method
	    # a little sneaky here: by calling resid before adding the
	    #   na.action method, I avoid having missings re-inserted
	    # I also make sure that it doesn't have to reconstruct X and Y
	    fit2 <- c(fit, list(x=X, y=Y, weights=weights))
	    if (length(strats)) fit2$strata <- strata.keep
	    if (length(cluster)) {
		temp <- residuals.coxph(fit2, type='dfbeta', collapse=cluster,
					  weighted=TRUE)
		# get score for null model
		if (is.null(init))
			fit2$linear.predictors <- 0*fit$linear.predictors
		else fit2$linear.predictors <- c(X %*% init)
		temp0 <- residuals.coxph(fit2, type='score', collapse=cluster,
					 weighted=TRUE)
		}
	    else {
		temp <- residuals.coxph(fit2, type='dfbeta', weighted=TRUE)
		fit2$linear.predictors <- 0*fit$linear.predictors
		temp0 <- residuals.coxph(fit2, type='score', weighted=TRUE)
	        }
	    fit$var <- t(temp) %*% temp
	    u <- apply(as.matrix(temp0), 2, sum)
	    fit$rscore <- coxph.wtest(t(temp0)%*%temp0, u, control$toler.chol)$test
	    }
	#Wald test
	if (length(fit$coefficients) && is.null(fit$wald.test)) {
	    #not for intercept only models, or if test is already done
	    nabeta <- !is.na(fit$coefficients)
	    if (is.null(init)) temp <- fit$coefficients[nabeta]
	    else temp <- (fit$coefficients - init)[nabeta]
	    fit$wald.test <-  coxph.wtest(fit$var[nabeta,nabeta], temp,
					  control$toler.chol)$test
	    }
	na.action <- attr(m, "na.action")
	if (length(na.action)) fit$na.action <- na.action
	if (model) fit$model <- m
        ## else { ## we might want model=T and X=T
        if (x)  {
            fit$x <- X
            if (length(strats)) fit$strata <- strata.keep
        }
        if (y)     fit$y <- Y
        ##    }
    }
    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights
    
    ##fit$formula <- as.vector(attr(Terms, "formula"))
    fit$formula<-formula(Terms) ## get the environments right
    fit$call <- call
    fit$method <- method
    fit
    }
