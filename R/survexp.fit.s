#  SCCS @(#)survexp.fit.s	5.2 11/03/98
#  Actually compute the expected survival for one or more cohorts
#    of subjects.  If each subject is his/her own group, it gives individual
#    survival
survexp.fit <- function(x, y, times, death, ratetable) {
    if (!is.matrix(x)) stop("x must be a matrix")
    if (ncol(x) != (1+length(dim(ratetable))))
	stop("x matrix does not match the rate table")
    atts <- attributes(ratetable)
    rfac <- atts$factor
    if (length(rfac) != ncol(x)-1) stop("Wrong length for rfac")
    ngrp <- max(x[,1])
    times <- sort(unique(times))
    if (any(times <0)) stop("Negative time point requested")
    if (missing(y))  y <- rep(max(times), nrow(x))
    ntime <- length(times)
    if (!is.logical(death)) stop("Invalid value for death indicator")

    cuts <- atts$cutpoints
    us.special <- (rfac >1)
    if (any(us.special)) {  #special handling for US pop tables
	if (sum(us.special) >1)
	    stop("Two columns marked for special handling as a US rate table")
	#slide entry date so that it appears that they were born on Jan 1
	cols <- 1+match(c("age", "year"), attr(ratetable, "dimid"))
	if (any(is.na(cols))) stop("Ratetable does not have expected shape")
	temp <- date.mdy(as.date(x[,cols[2]]) -x[,cols[1]])
	x[,cols[2]] <- x[,cols[2]] - mdy.date(temp$month, temp$day, 1960)
	# Doctor up "cutpoints"
	temp <- (1:length(rfac))[us.special]
	nyear <- length(cuts[[temp]])
	nint <- rfac[temp]       #intervals to interpolate over
	cuts[[temp]] <- round(approx(nint*(1:nyear), cuts[[temp]],
					nint:(nint*nyear))$y - .0001)
	}
    temp <- .C('pyears3',
		    as.integer(death),
		    as.integer(nrow(x)),
		    as.integer(length(atts$dim)),
		    as.integer(rfac),
		    as.integer(atts$dim),
		    as.double(unlist(cuts)),
		    ratetable,
		    as.double(x),
		    as.double(y),
		    as.integer(ntime),
		    as.integer(ngrp),
		    as.double(times),
		    surv = double(ntime * ngrp),
		    n   = integer(ntime *ngrp),
               PACKAGE="survival")
    if (ntime==1) list(surv=temp$surv, n=temp$n)
    else if (ngrp >1)
	 list(surv=apply(matrix(temp$surv, ntime, ngrp),2,cumprod),
		 n=   matrix(temp$n, ntime, ngrp))
    else list(surv=cumprod(temp$surv), n=temp$n)
    }
