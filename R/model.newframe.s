#SCCS 04/14/92 @(#)model.newframe.s	4.3
# This function is called if you want to get a new data frame,
#   usually for prediction.  It's main problem is to "glue" any
#   transform specific information back onto the formula, so that
#   data dependent transforms work as they used to.
# It only works if the data dependent functions are not inside another one,
#   so  sqrt(age - min(age)) is out of luck.  It also only works for those
#   transforms that support it by adding data dependent info as an attribute
#   of their output.
# If you know this isn't so, then safe=T uses a method that is much longer,
#   but is guarranteed to work, see predict.gam

model.newframe <- function(object, newdata, safe=F, response=F, ...) {
    if (inherits(object, 'terms'))  Terms <- object
    else {
	Terms <- object$terms
	if (!inherits(Terms, 'terms'))
	    stop ("Invalid terms component of object")
	}
    offset <- attr(Terms, 'offset')

    # First, is newdata just a list of numbers?
    if (is.numeric(newdata)) {
	nvar <- length(attr(Terms,"term.labels")) + length(offset)
	if (length(newdata)>1  || newdata!=floor(newdata)  || newdata<0){ #It's not just a frame number
	    if (is.matrix(newdata) && ncol(newdata) == nvar)
		   return(newdata)
	    else if (length(newdata) == nvar)
		   return(matrix(newdata,1,nvar))
	    else stop("Argument \"newdata\" cannot be coerced to an appropriate model matrix")
	    }
	}

    # newdata is a list, data frame, or frame number
    if (!safe) {
	#augment the arguments with extra parameters
	  #someday
	if (!response) Terms <- delete.response(Terms)
	model.frame(Terms, newdata, ...)
	}
    else {
	#Do a safe call, by building up a brand new model frame
	Call <- object$call
	Call[[1]] <- as.name("model.frame")
	Call$formula <- terms.inner(formula(object))
   #might need to tack on the response here!
	if (response) stop("Not implimented yet for safe=T, response=T")
	Call$na.action <- function(x)  x
	Call <- Call[match(c("", "formula", "data", "subset", "na.action"),
	    names(Call), 0)]
	data <- eval(Call)
	attr(data, "terms") <- NULL
	Call$subset <- NULL
	Call$data <- substitute(newdata)
	newdata <- eval(Call)
	attr(newdata, "terms") <- NULL
	d2 <- dim(newdata)
	if(d2[1] < 1)
	    stop("0 rows in newdata")
	d1 <- dim(data)
	if(d1[2] != d2[2])  #newdata missing some variables
	    data <- data[, names(newdata), drop = F]
	data[seq(d1[1] + 1, d1[1] + d2[1]),  ] <- newdata  #rbind the new on
	attr(data, "row.names") <- c(rep("OLD DATA",d1[1]), row.names(newdata))
	#Now compute the combined model frame, excluding the response
	na.action <- eval(object$call$na.action)
	Terms <- object$terms
	Terms <- delete.response(Terms)
	model.frame(Terms, data, na.action = na.action)
	}
    }
