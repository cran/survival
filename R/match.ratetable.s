# SCCS @(#)match.ratetable.s	4.5 11/18/97
# Do a set of error checks on whether the ratetable() vars match the
#   actual ratetable
# This is called by pyears and survexp, but not by users
#
# Returns a subscripting vector and a call
#
match.ratetable <- function(R, ratetable) {
    attR <- attributes(R)
    attributes(R) <- attR['dim']     #other attrs get in the way later
    if (!is.ratetable(ratetable)) stop("Invalid rate table")
    dimid <- attr(ratetable, 'dimid')
    ord <- match(attR$dimnames[[2]], dimid)
    if (any(is.na(ord)))
       stop(paste("Argument '", (attR$dimnames[[2]])[is.na(ord)],
	    "' in ratetable()",
	    " does not match the given table of event rates", sep=''))
    nd <- length(ord)
    if (nd != length(dimid))
	stop("The ratetable() call has the wrong number of arguments")
    ord[ord] <- 1:nd   #reverse the index, so "ord" can be on the right-hand
    R <- R[,ord,drop=FALSE]

    # Check out the dimensions of R --
    const <- attR[["constants"]][ord]
    call <- "ratetable["
    levlist <- attR[['levlist']][ord]
    dtemp <-dimnames(ratetable)
    efac  <- attr(ratetable, 'factor')
    for (i in (1:nd)) {
	if (const[i]) {   #user put in a constant
	    temp <- match(levlist[[i]], dtemp[[i]])
	    if (is.na(temp)) {
		temp <- as.numeric(levlist[[i]])
		if (is.na(temp))
		       stop(paste("Invalid value in ratetable() for variable",
				 dimid[i]))
		if (efac[i]==1) {  # this level is a factor
		    if (temp<=0 || temp!=floor(temp) || temp >length(dtemp[[i]]))
		       stop(paste("Invalid value in ratetable() for variable",
				 dimid[i]))
		    }
		else stop(paste("Invalid value in ratetable() for variable",
					dimid[i]))
		}
	    R[,i] <- temp
	    call <- paste(call, temp)
	    }
	else if (length(levlist[[i]]) >0) {  #factor or character variable
	    if (efac[i]!=1) stop(paste("In ratetable(),", dimid[i],
				     "must be a continuous variable"))
	    temp <- match(levlist[[i]], dtemp[[i]])
	    if (any(is.na(temp)))
		stop(paste("Levels do not match for ratetable() variable",
			    dimid[i]))
	    R[,i] <- temp[R[,i]]
	    }
	else {   # ratetable() thinks it is a continuous variable
	    if (efac[i]==1) {   #but it's not-- make sure it is an integer
		temp <- R[,i]
		if (any(floor(temp)!=temp) || any(temp<=0) ||
			    max(temp) > length(dtemp[[i]]))
		stop(paste("In ratetable(),",dimid[i],"is out of range"))
		}
	    }
	if (i==nd) call <- paste(call, "]")
	else       call <- paste(call, ",")
	}

    summ <- attr(ratetable, 'summary')
    if (is.null(summ))
	 list(R= R[,!const, drop=FALSE], call={if(any(const)) call else NULL})
    else list(R= R[,!const, drop=FALSE], call={if(any(const)) call else NULL},
		summ=summ(R))
    }
