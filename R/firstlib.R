.First.lib <- function(lib, pkg) {
          library.dynam("survival", pkg, lib)
        }

## survfit.print.n=="start" is compatible with previous R and with BDR
if (is.null(getOption("survfit.print.n")))
    options(survfit.print.n="start")


"[.terms" <-function (termobj, i) {
        resp <- if (attr(termobj, "response")) 
                termobj[[2]]
        else NULL
        newformula <- attr(termobj, "term.labels")[i]
        if (length(newformula) == 0) 
                newformula <- 1
        newformula <- reformulate(newformula, resp)
        environment(newformula)<-environment(termobj)
        terms(newformula, specials = names(attr(termobj, "specials")))

}

is.category <- function(x) inherits(x,"factor") || is.factor(x)


if (R.version$minor<3 && R.version$major==1){
    ## these are from survival4, now (>1.3.0) base R
    naprint.omit <- function(x)
        paste(length(x), "deleted due to missing")
    
    ## Put the missing values back into a vector.
    ##   And be careful about the labels too.
    naresid.omit <- function(omit, x) {
        if (length(omit)==0 || !is.numeric(omit))
            stop("Invalid argument for 'omit'")
        
        if (is.matrix(x)) {
            n <- nrow(x)
            keep <- rep(NA,n+ length(omit))
            keep[-omit] <- 1:n
            x <- x[keep,,drop=FALSE]
            temp <- dimnames(x)[[1]]
            if (length(temp)) {
                temp[omit] <- names(omit)
                dimnames(x) <- list(temp, dimnames(x)[[2]])
	    }
            x
	}
        else {
            n <- length(x)
            keep <- rep(NA,n+ length(omit))
            keep[-omit] <- 1:n
            x <- x[keep]
            temp <- names(x)
            if (length(temp)) {
                temp[omit] <- names(omit)
                names(x) <- temp
	    }
            x
	}
    }
    naprint.omit <- function(x)
        paste(length(x), "observations deleted due to missing")
    naprint <- function(x, ...)
        UseMethod("naprint",x,...)
    naresid <- function(x, ...)
        UseMethod("naresid",x,...)
    
    naprint.default <- function(...)
        return("")
###
}

labels.survreg <- function(object, ...) attr(object,"term.labels")

