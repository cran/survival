# SCCS @(#)ridge.s	1.1 12/22/98
ridge <- function(..., theta, df=nvar/2, eps=.1, scale=T) {
    x <- cbind(...)
    nvar <- ncol(x)
    xname <- as.character(parse(text=substitute(cbind(...))))[-1]
    vars <- apply(x, 2, function(z) var(z[!is.na(z)]))
    class(x) <- 'coxph.penalty'

    if (!missing(theta) && !missing(df))
	    stop("Only one of df or theta can be specified")

    if (scale) 
	    pfun <- function(coef,theta, ndead, scale) {
		list(penalty= sum(coef^2 *scale)*theta/2,
		     first  = theta*coef*scale,
		     second = theta*scale,
		     flag=F)
		}
    else
	    pfun <- function(coef,theta, ndead, scale) {
		list(penalty= sum(coef^2)*theta/2,
		     first  = theta*coef,
		     second = theta,
		     flag=F)
		}


    if (!missing(theta)) {
	temp <- list(pfun=pfun,
		     diag=T,
		     cfun=function(parms, iter, history) {
				list(theta=parms$theta, done=T) }, 
		     cparm=list(theta= theta),
		     pparm= vars,
		     varname=paste('ridge(', xname, ')', sep=''))
	}
    else {
	temp <- list(pfun=pfun,
		     diag=T,
		     cfun=frailty.controldf,
		     cargs = 'df',
		     cparm=list(df=df, eps=eps, thetas=0, dfs=nvar,
		         guess=1),
		     pparm= vars,
		     varname=paste('ridge(', xname, ')', sep=''))
	}
	
    attributes(x) <- c(attributes(x), temp)
    x
    }
