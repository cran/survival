.First.lib <- function(lib, pkg) {
          library.dynam("survival", pkg, lib)
        }

## survfit.print.n=="start" is compatible with previous R and with BDR
if (is.null(getOption("survfit.print.n")))
    options(survfit.print.n="start")
## survfit.print.mean==TRUE is compatible with previous R/SPLUS (but is silly)
if (is.null(getOption("survfit.print.mean")))
    options(survfit.print.mean=TRUE)


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



labels.survreg <- function(object, ...) attr(object,"term.labels")

