.onLoad <- function(lib, pkg) {
	 ## moved to NAMESPACE
          ##library.dynam("survival", pkg, lib)
          ## survfit.print.n=="start" is compatible with previous R
          ##     and with MASS
          if (is.null(getOption("survfit.print.n")))
              options(survfit.print.n="start")
          ## survfit.print.mean==TRUE is compatible with previous R/SPLUS
          ##     (but is silly)
          if (is.null(getOption("survfit.print.mean")))
              options(survfit.print.mean=TRUE)
      }



is.category <- function(x) inherits(x,"factor") || is.factor(x)



labels.survreg <- function(object, ...) attr(object,"term.labels")

