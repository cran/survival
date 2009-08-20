# $Id: labels.survreg.S 11059 2008-10-23 12:32:50Z therneau $
labels.survreg <- function(object, ...)
	attr(object$terms, "term.labels")

