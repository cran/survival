#SCCS @(#)survfit.s	4.19 09/08/00
survfit <- function (formula, data, weights, subset, na.action, ...) {
    call <- match.call()
    # Real tricky -- find out if the first arg is "Surv(...)" without
    #  evaluating it.  If this is so, or it is a survival object, turn it
    #  into a formula
    # (This allows people to leave off the "~1" from a formula)
    if ((mode(call[[2]]) == 'call' &&  call[[2]][[1]] == as.name('Surv'))
		|| inherits(formula, 'Surv'))  {
        formula<-eval(parse(text=paste(deparse(call[[2]]),1,sep="~")))
        ## need to add back the formula environment
        environment(formula)<-parent.frame()
	}

    # if the first object is a Cox model, call survfit.coxph
    if (!inherits(formula, 'formula')) temp <- UseMethod("survfit")
    else {
	# Ok, I have a formula
        # grab the data and process it
	m <- match.call(expand.dots=FALSE)
	m$... <- NULL

	Terms <- terms(formula, 'strata')
	ord <- attr(Terms, 'order')
	if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")
	m$formula <- Terms
	m[[1]] <- as.name("model.frame")
	m <- eval(m, parent.frame())

	n <- nrow(m)
	Y <- model.extract(m, "response")
	if (!is.Surv(Y)) stop("Response must be a survival object")

	casewt <- model.extract(m, "weights")
	# The second line below works around a bug in Splus 3.0.1, which later
	#    went away, i.e., casewt is returned as an unevaluated arg.
	if (is.null(casewt)) casewt <- rep(1,n)
	##else if (mode(casewt)=='argument') casewt <- eval(casewt[[1]])

	if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

	ll <- attr(Terms, 'term.labels')
	if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
	else X <- strata(m[ll])

	temp <- survfit.km(X, Y, casewt, ...)
	class(temp) <- "survfit"
	if (!is.null(attr(m, 'na.action'))) 
		temp$na.action <- attr(m, 'na.action')
	}
    temp$call <- call
    # added to change INFs to NAs for C program - cmb 8/25/2000
    ## I don't think we want this <TSL>
    ##if (any(is.inf(temp$std.err))) temp$std.err[is.inf(temp$std.err)] <- NA
    temp
    }

# The subscript function is bundled in here, although used most
#  often in plotting

"[.survfit" <- function(x, ..., drop=FALSE) {
    if (missing(..1)) i<- NULL  else i <- sort(..1)
    if (missing(..2)) j<- NULL  else j <- ..2
    if (is.null(x$strata)) {
	if (is.matrix(x$surv)) {
	    x$surv <- x$surv[,i,drop=drop]
	    if (!is.null(x$std.err)) x$std.err <- x$std.err[,i,drop=drop]
	    if (!is.null(x$upper)) x$upper <- x$upper[,i,drop=drop]
	    if (!is.null(x$lower)) x$lower <- x$lower[,i,drop=drop]
	    }
	else warning("Survfit object has only a single survival curve")
	}
    else {
	if (is.null(i)) keep <- seq(along.with=x$time)
	else {
	    if (is.null(x$ntimes.strata)) strata.var <- x$strata
	    else strata.var <- x$ntimes.strata
	    if (is.character(i)) strat <- rep(names(x$strata), strata.var)
	    else                 strat <- rep(1:length(x$strata), strata.var)
	    keep <- seq(along.with=strat)[match(strat, i, nomatch=0)>0]
	    if (length(i) <=1) x$strata <- NULL
	    else               x$strata  <- x$strata[i]
	    if (!is.null(x$ntimes.strata)) {
		x$strata.all <- x$strata.all[i]
		x$ntimes.strata <- x$ntimes.strata[i]
	        }
	    x$time    <- x$time[keep]
	    x$n.risk  <- x$n.risk[keep]
	    x$n.event <- x$n.event[keep]
	    }
	if (is.matrix(x$surv)) {
	    if (is.null(j)) {
		x$surv <- x$surv[keep,,drop=drop]
		if (!is.null(x$std.err)) 
			x$std.err <- x$std.err[keep,,drop=drop]
		if (!is.null(x$upper)) x$upper <-x$upper[keep,,drop=drop]
		if (!is.null(x$lower)) x$lower <-x$lower[keep,,drop=drop]
		}
	    else {
		x$surv <- x$surv[keep,j]
		if (!is.null(x$std.err)) x$std.err <- x$std.err[keep,j]
		if (!is.null(x$upper)) x$upper <- x$upper[keep,j]
		if (!is.null(x$lower)) x$lower <- x$lower[keep,j]
		}
	    }
	else {
	    x$surv <- x$surv[keep]
	    if (!is.null(x$enter)) x$enter <- x$enter[keep]
	    if (!is.null(x$exit.censored))
		    x$exit.censored <- x$exit.censored[keep]
	    if (!is.null(x$std.err)) x$std.err <- x$std.err[keep]
	    if (!is.null(x$upper)) x$upper <- x$upper[keep]
	    if (!is.null(x$lower)) x$lower <- x$lower[keep]
	    }
	}
    x
    }

basehaz<-function(fit,centered=TRUE){
    if(!inherits(fit,"coxph"))
        stop("must be a coxph object")

    sfit<-survfit(fit)

    H<- -log(sfit$surv)

    strata<-sfit$strata
    if (!is.null(strata))
        strata<-rep(names(strata),strata)
    
    if (!centered){
        z0<-fit$means
        bz0<-sum(z0*coef(fit))
        H<- H*exp(-bz0)
    }

    if (is.null(strata))
      return(data.frame(hazard=H,time=sfit$time))
    else
      return(data.frame(hazard=H,time=sfit$time,strata=strata))

}
