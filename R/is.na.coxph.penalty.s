#  SCCS  @(#)is.na.coxph.penalty.s	1.4 02/21/99
# The subscript function for coxph.penalty objects
#  without it the "subset" arg of a model statement tosses
#  away all of the attributes
#
"[.coxph.penalty" <- function(x, ..., drop=F) {
    attlist <- attributes(x)
    attributes(x) <- attlist[match(c('dim', 'dimnames'), names(attlist), 0)] 
    x <- NextMethod('[')  #let the default method do actual subscripting

    # Tack back on all of the old attributes except dim and dimnames,
    #   which will have been properly modified by the standard [ method
    attributes(x) <- c(attributes(x),attlist[is.na(match(names(attlist),
						      c("dim", "dimnames")))])
    x
    }
			  

is.na.coxph.penalty <- function(x) {
    if (is.matrix(x)) is.na(c(unclass(x) %*% rep(1,ncol(x))))
    else is.na(unclass(x))
    }
