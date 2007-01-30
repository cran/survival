# SCCS @(#)Surv.s	5.5 07/09/00
# Package up surivival type data as a structure
#
Surv <- function(time, time2, event,
	      type=c('right', 'left', 'interval', 'counting', 'interval2'),
		       origin=0) {
    nn <- length(time)
    ng <- nargs()
    if (missing(type)) {
	if (ng==1 || ng==2) type <- 'right'
	else if (ng==3)     type <- 'counting'
	else stop("Invalid number of arguments")
	}
    else {
	type <- match.arg(type)
	ng <- ng-1
	if (ng!=3 && (type=='interval' || type =='counting'))
		stop("Wrong number of args for this type of survival data")
	if (ng!=2 && (type=='right' || type=='left' ||  type=='interval2'))
		stop("Wrong number of args for this type of survival data")
	}
    who <- !is.na(time)

    if (ng==1) {
	if (!is.numeric(time)) stop ("Time variable is not numeric")
	ss <- cbind(time, 1)
	dimnames(ss) <- list(NULL, c("time", "status"))
	}
    else if (type=='right' || type=='left') {
	if (!is.numeric(time)) stop ("Time variable is not numeric")
	if (length(time2) != nn) stop ("Time and status are different lengths")
	if (is.logical(time2)) status <- 1*time2
	    else  if (is.numeric(time2)) {
		who2 <- !is.na(time2)
		if (max(time2[who2]) ==2) status <- time2 -1
		else status <- time2
		if (any(status[who2] !=0  & status[who2]!=1))
				stop ("Invalid status value")
		}
	    else stop("Invalid status value")
	 ss <- cbind(time, status)
	 dimnames(ss) <- list(NULL, c("time", "status"))
	}
    else  if (type=='counting') {
	if (length(time2) !=nn) stop ("Start and stop are different lengths")
	if (length(event)!=nn) stop ("Start and event are different lengths")
	if (!is.numeric(time))stop("Start time is not numeric")
	if (!is.numeric(time2)) stop("Stop time is not numeric")
	who3 <- who & !is.na(time2)
	if (any (time[who3]>= time2[who3]))stop("Stop time must be > start time")
	if (is.logical(event)) status <- 1*event
	    else  if (is.numeric(event)) {
		who2 <- !is.na(event)
		if (max(event[who2])==2) status <- event - 1
		else status <- event
		if (any(status[who2] !=0  & status[who2]!=1))
				stop("Invalid status value")
		}
	    else stop("Invalid status value")
	ss <- cbind(time-origin, time2-origin, status)
	}

    else {  #interval censored data
	if (type=='interval2') {
	    event <- ifelse(is.na(time), 2,
		     ifelse(is.na(time2),0,
		     ifelse(time==time2, 1,3)))
	    if (any(time[event==3] > time2[event==3]))
		stop("Invalid interval: start > stop")
	    time <- ifelse(event!=2, time, time2)
	    type <- 'interval'
	    }
	else {
	    temp <- event[!is.na(event)]
	    if (!is.numeric(temp)) stop("Status indicator must be numeric")
	    if (length(temp)>0 && any(temp!= floor(temp) | temp<0 | temp>3))
		stop("Status indicator must be 0, 1, 2 or 3")
	    }
	status <- event
	ss <- cbind(time, ifelse(!is.na(event) & event==3, time2, 1),
			    status)
	}

    attr(ss, "type")  <- type
    class(ss) <- 'Surv'
    ss
    }

print.Surv <- function(x, quote=FALSE, ...)
    invisible(print(as.character.Surv(x), quote=quote, ...))

as.character.Surv <- function(x, ...) {
    class(x) <- NULL
    type <- attr(x, 'type')
    if (type=='right') {
	temp <- x[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	paste(format(x[,1]), temp, sep='')
	}
    else if (type=='counting') {
	temp <- x[,3]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	paste('(', format(x[,1]), ',', format(x[,2]), temp,
			 ']', sep='')
	}
    else if (type=='left') {
	temp <- x[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "<"," "))
	paste(temp, format(x[,1]), sep='')
	}
    else {   #interval type
	stat <- x[,3]
	temp <- c("+", "", "-", "]")[stat+1]
	temp2 <- ifelse(stat==3,
			 paste("[", format(x[,1]), ", ",format(x[,2]), sep=''),
			 format(x[,1]))
	ifelse(is.na(stat), as.character(NA), paste(temp2, temp, sep=''))
	}
    }

## the ... handling here works in R1.2 but not R1.1.1
##
##"[.Surv" <- function(x, ..., drop=F) {
##    # If only 1 subscript is given, the result will still be a Surv object
##    #  If the second is given extract the relevant columns as a matrix
##    if (missing(..2)) {
##	temp <- class(x)
##	type <- attr(x, "type")
##	class(x) <- NULL
##	x <- x[..., drop=F]
##	class(x) <- temp
##	attr(x, "type") <- type
##	x
##	}
##    else {
##	class(x) <- NULL
##	NextMethod("[")
##	}
##    }

"[.Surv" <- function(x, i,j, drop=FALSE) {
    # If only 1 subscript is given, the result will still be a Surv object
    #  If the second is given extract the relevant columns as a matrix
    if (missing(j)) {
	temp <- class(x)
	type <- attr(x, "type")
	class(x) <- NULL
	x <- x[i, , drop=FALSE]
	class(x) <- temp
	attr(x, "type") <- type
	x
	}
    else {
	class(x) <- NULL
	NextMethod("[")
	}
    }

is.na.Surv <- function(x) {
    as.vector( (1* is.na(unclass(x)))%*% rep(1, ncol(x)) >0)
    }

Math.Surv <- function(...)  stop("Invalid operation on a survival time")
Ops.Surv  <- function(...)  stop("Invalid operation on a survival time")
Summary.Surv<-function(...) stop("Invalid operation on a survival time")
is.Surv <- function(x) inherits(x, 'Surv')
