#
# SCCS @(#)is.ratetable.s	4.6 08/01/98
#
is.ratetable <- function(x, verbose=F) {
    if (!verbose) {
	if (!inherits(x, 'ratetable')) return(F)
	att <- attributes(x)
	if (any(is.na(match(c("dim", "dimnames", "dimid", "factor", "cutpoints"),
			    names(att))))) return(F)
	nd <- length(att$dim)
	if (length(x) != prod(att$dim)) return(F)
	if (!(is.list(att$dimnames) && is.list(att$cutpoints)))
		 return(F)
	if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
			 length(att$cutpoints)!=nd) return(F)
	fac <- as.numeric(att$factor)
	if (any(is.na(fac))) return(F)
	if (any(fac <0)) return(F)
	for (i in 1:nd) {
	    n <- att$dim[i]
	    if (length(att$dimnames[[i]]) !=n) return(F)
	    if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) return(F)
	    if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) return(F)
	    if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  return(F)
	    if (fac[i]>1 && i<nd) return(F)
	    }
	return(T)
	}

    #verbose return messages, useful for debugging
    msg <- NULL
    if (!inherits(x, 'ratetable')) msg <- c(msg, "wrong class")
    att <- attributes(x)
    if (any(is.na(match(c("dim", "dimnames", "dimid", "factor", "cutpoints"),
			names(att))))) msg <- c(msg, 'missing an attribute')
    nd <- length(att$dim)
    if (length(x) != prod(att$dim)) msg <- c(msg, 'dims dont match length')
    if (!(is.list(att$dimnames) && is.list(att$cutpoints)))
	     msg <- c(msg, 'dimnames or cutpoints not a list')
    if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
		     length(att$cutpoints)!=nd) msg <- c(msg, 'bad lengths')
    fac <- as.numeric(att$factor)
    if (any(is.na(fac))) msg <- c(msg, 'a missing factor')
    if (any(fac <0)) msg <- c(msg, 'factor <0')
    for (i in 1:nd) {
	n <- att$dim[i]
	if (length(att$dimnames[[i]]) !=n) 
		msg <- c(msg, 'dimnames wrong length')
	if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) 
		msg <- c(msg, 'cutpnt missing')
	if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) 
		msg <- c(msg, 'unsorted cutpoints')
	if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  
		msg <- c(msg, 'cutpnt should be null')
	if (fac[i]>1 && i<nd) 
		msg <- c(msg, 'only the last dim can be interpolated')
	}
    if (length(msg)==0) T
    else msg
    }
