#  SCCS  @(#)survreg.old.s	1.4 02/21/99
# Map the argument list of the old survreg to the new one
#
survreg.old <- function(formula, data=sys.frame(sys.parent()), ...,
        link=c('log',"identity"),
        dist=c("extreme", "logistic", "gaussian", "exponential",
               "rayleigh", "weibull"),
	fixed=list()) {
    
    dist <- match.arg(dist)
    link <- match.arg(link)
    
    if ((dist!='weibull' && dist != 'rayleigh') && link=='log') {
	if (dist=='extreme') dist <- 'weibull'
	else dist <- paste('log', dist, sep='')
	}
    if (is.null(fixed$scale)) scale <- 0
    else scale <- fixed$scale

    survreg(formula, data, ..., dist=dist, scale=scale)
    }
