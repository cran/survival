# SCCS @(#)ratetable.s	5.2 09/25/98
#
# This is a 'specials' function for pyears
#   it is a stripped down version of as.matrix(data.frame(...))
# There is no function to create a ratetable.
# This function has a class, only so that data frame subscripting will work
#
ratetable <- function(...) {
    args <- list(...)
    nargs <- length(args)
    ll <- sapply(args, length)
    n <- max(ll)
    levlist <- vector("list", nargs)
    x <- matrix(0,n,nargs)
    dimnames(x) <- list(1:n, names(args))
    for (i in 1:nargs) {
	if (ll[i] ==n) {
	    if (!is.numeric(args[[i]])) args[[i]] <- factor(args[[i]])
	    if (is.factor(args[[i]])) {
		levlist[[i]] <- levels(args[[i]])
		x[,i] <- c(args[[i]])
		}
	    else x[,i] <- args[[i]]
	    }
	else if (ll[i] ==1) levlist[i] <- args[i]
	else stop("All arguments to ratetable() must be the same length")
	}
    attr(x, "constants") <- (ll==1) & (n>1)
    attr(x, "levlist")   <- levlist
    class(x)  <- "ratetable2"
    x
    }

# The two functions below should only be called internally, when missing
#   values cause model.frame to drop some rows
is.na.ratetable2 <- function(x) {
    attributes(x) <- list(dim=dim(x))
    as.vector((1 * is.na(x)) %*% rep(1, ncol(x)) >0)
    }
"[.ratetable2" <- function(x, rows, cols, drop=FALSE) {
    if (!missing(cols)) {
       stop("This should never be called!")
       }
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- x[rows,,drop=FALSE]
    attr(y,'constants') <- aa$constants
    attr(y,'levlist')   <- aa$levlist
    class(y) <- 'ratetable2'
    y
    }

#
# Functions to manipulate rate tables
#

"[.ratetable" <- function(x, ..., drop=TRUE) {
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- NextMethod("[", drop=FALSE)
    newdim <- attr(y, 'dim')
    if (is.null(newdim)) stop("Invalid subscript")
    dropped <- (newdim==1)
    if (drop)  change <- (newdim!=aa$dim & !dropped)
    else       change <- (newdim!=aa$dim)

    if (any(change)) {  #dims that got smaller, but not dropped
	newcut <- aa$cutpoints
	for (i in (1:length(change))[change])
	    if (!is.null(newcut[[i]])) newcut[[i]] <-
		(newcut[[i]])[match(dimnames(y)[[i]], aa$dimnames[[i]])]
	aa$cutpoints <- newcut
	}
    if (drop && any(dropped)){
	if (all(dropped)) as.numeric(y)   #single element
	else {
	    #Note that we have to drop the summary function
	    attributes(y) <- list( dim = dim(y)[!dropped],
				   dimnames = dimnames(y)[!dropped],
				   dimid = aa$dimid[!dropped],
				   factor = aa$factor[!dropped],
				   cutpoints =aa$cutpoints[!dropped])
	    class(y) <- 'ratetable'
	    y
	    }
	}
    else {
	aa$dim <- aa$dimnames <- NULL
	attributes(y) <- c(attributes(y), aa)
	y
	}
    }

is.na.ratetable  <- function(x)
    structure(is.na(as.vector(x)), dim=dim(x), dimnames=dimnames(x))

Math.ratetable <- function(x, ...) {
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    NextMethod(.Generic)
    }

Ops.ratetable <- function(e1, e2) {
    #just treat it as an array
    if (nchar(.Method[1])) attributes(e1) <- attributes(e1)[c("dim","dimnames")]
    if (nchar(.Method[2])) attributes(e2) <- attributes(e2)[c("dim","dimnames")]
    NextMethod(.Generic)
    }

as.matrix.ratetable <- function(x) {
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    x
    }

print.ratetable <- function(x, ...)  {
    cat ("Rate table with dimension(s):", attr(x, 'dimid'), "\n")
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    NextMethod()
    }
