# SCCS @(#)ridge.s	1.1 03/26/98
ridge <- function(..., theta=1) {
    x <- cbind(...)
    xname <- as.character(parse(text=substitute(cbind(...))))[-1]
    vars <- apply(x, 2, function(z) var(z[!is.na(z)]))
    class(x) <- 'coxph.penalty'
    temp <- list( pfun=function(coef,theta, ndead, scale) {
		       list(coef=coef, 
				    penalty= sum(coef^2 *scale)*theta/2,
				    first  = theta*coef*scale,
			       	    second = theta*scale,
			            flag=F)},
		  diag=T,
		  cfun=function(parms, iter, history) {
				list(theta=parms$theta, done=T) }, 
		  cparm=list(theta= theta),
		  pparm= vars,
		  varname=paste('ridge(', xname, ')', sep=''))

    attributes(x) <- c(attributes(x), temp)
    x
    }
