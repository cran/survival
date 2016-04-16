# Automatically generated from the noweb directory
# Methods for survfitms objects
summary.survfit <- function(object, times, censored=FALSE, 
                            scale=1, extend=FALSE, 
                            rmean=getOption('survfit.rmean'),
                            ...) {
    fit <- object
    if (!inherits(fit, 'survfit'))
            stop("summary.survfit can only be used for survfit objects")

    # The print.rmean option is depreciated, it is still listened
    #   to in print.survfit, but ignored here
    if (is.null(rmean)) rmean <- "common"
    if (is.numeric(rmean)) {
        if (is.null(object$start.time)) {
            if (rmean < min(object$time)) 
                stop("Truncation point for the mean is < smallest survival")
        }
        else if (rmean < object$start.time)
            stop("Truncation point for the mean is < smallest survival")
    }
    else {
        rmean <- match.arg(rmean, c('none', 'common', 'individual'))
        if (length(rmean)==0) stop("Invalid value for rmean option")
    }

    temp <- survmean(fit, scale=scale, rmean)  
    table <- temp$matrix  #for inclusion in the output list
    rmean.endtime <- temp$end.time

    if (!missing(times)) {
        if (!is.numeric(times)) stop ("times must be numeric")
        times <- sort(times)
    }

    # The fit$surv object is sometimes a vector and sometimes a
    #  matrix.  We calculate row indices first, and then deal
    #  with the cases at the end.
    nsurv <- length(fit$time)
    if (is.null(fit$strata)) {
        nstrat <- 1
        stemp <- rep(1L, nsurv)
        strata.names <- ""
        }
    else   {
        nstrat <- length(fit$strata)
        stemp <- rep(1:nstrat, fit$strata)
        strata.names <- names(fit$strata)
    }

    if (missing(times)) {
        # just pick off the appropriate rows of the output
        # For a survfitms object n.event is a matrix, pick off all rows with an
        #  event for some endpoint.
        if (censored) indx1 <- seq(along=fit$time)
        else indx1 <- which(rowSums(as.matrix(fit$n.event)) >0)
        indx2 <- indx1
    }
    else { 
        find2 <- function(x, vec, left.open=FALSE, ...) {
            if (!left.open) findInterval(x, vec, ...)
            else length(vec) - findInterval(-x, rev(-vec), ...)
        }
        # Process the curves one at a time, adding them to the two lists
        ilist1 <- ilist2 <- ilist3 <- vector('list', nstrat)
        newtime <- ilist1
        n <- length(stemp)
        for (i in 1:nstrat) {
            who <- (1:n)[stemp==i]  # the rows of the object for this strata
            stime <- fit$time[who]

            # First, toss any printing times that are outside our range
            if (is.null(fit$start.time)) mintime <- min(stime, 0)
            else                         mintime <- fit$start.time
            ptimes <- times[times >= mintime]

            if (!extend) {
                maxtime <- max(stime)
                ptimes <- ptimes[ptimes <= maxtime]
                }
            j <- find2(ptimes, stime) 
            ilist1[[i]] <- c(0, who)[1+ j]
            ilist2[[i]] <- c(0, who)[1+ find2(ptimes, stime, left.open=TRUE)]
            ilist3[[i]] <- j  #index within a group
            newtime[[i]] <- ptimes
        }
        indx1 <- unlist(ilist1)
        indx2 <- unlist(ilist2)

        # All of the indices (ilist1, indx1, ...) contain 0 for a time point that
        #  is prior to the first observed time in the curve.  Times that
        #  are >= to the last observed time will point to that last observed
        #  time.  Variable ilist3 contains indices that are relative to the
        #  start of a curve, all other indices point to row numbers in the
        #  entire object.
        #  
        cfun <- function(x, init=0) {  #cumulative counts over a time interval
            tlist <- vector("list", nstrat)
            if (is.matrix(x)) {
                for (i in 1:nstrat) {
                    # stemp is 1,1,1,....2,2,2,,.. to mark curves
                    x2 <- x[stemp==i,]  # all those in the group
                    j  <- c(0, ilist3[[i]])
                    tlist[[i]] <- apply(rbind(0, x2), 2, function(z) {
                        diff(cumsum(z)[1+j])})
                }
                matrix(unlist(lapply(tlist, t)), byrow=T, ncol=ncol(x))
            } 
            else {
                for (i in 1:nstrat) {
                    x2 <- x[stemp==i] 
                    j  <- c(0, ilist3[[i]])
                    tlist[[i]] <- diff(cumsum(c(0,x2))[1+j])
                }
                unlist(tlist)
            }
        }
    }

    # Create an output structure
    temp <- object
    temp$table <- table
    if (length(rmean.endtime)>0  && !is.na(rmean.endtime)) 
            temp$rmean.endtime <- rmean.endtime
    if (length(indx1)==length(fit$time) && all(indx1 == seq(along=fit$time))) {
        temp$time <- temp$time/scale
        if (!is.null(temp$strata))
            temp$strata <- factor(stemp, labels=strata.names)

    }
    else if (missing(times)) {  #default censor=FALSE case
        temp$time <- temp$time[indx1]/scale
        for (j in c("n.risk", "n.event", "n.censor", "n.enter", "prev",
                    "surv", "std.err", "cumhaz", "lower", "upper")) {
            zed <- temp[[j]]
            if (!is.null(zed)) {
                if (is.matrix(zed)) temp[[j]] <- zed[indx1,,drop=FALSE]
                else temp[[j]] <- zed[indx1]
            }
        }
        if (!is.null(temp$strata))
            temp$strata <- factor(stemp[indx1], levels=1:nstrat,
                                  labels=strata.names)
    }
    else { #times argument was given
        temp$time <- unlist(newtime)/scale
        tfun <- function(x, init=0, index=indx1) {
             if (is.matrix(x)) 
                rbind(rep(init, ncol(x)), x)[1+index,,drop=FALSE]
             else c(init, x)[1 + index]
        }
        tfun2 <- function(x, end=0, index=indx2) {
             if (is.matrix(x)) 
                rbind(x, rep(end, ncol(x)))[1+index,,drop=FALSE]
             else c(x,end)[1 + index]
        }
        temp$surv  <- tfun(temp$surv, 1)
        temp$n.risk <- tfun2(temp$n.risk)
        for (j in c("std.err", "cumhaz", "lower", "upper")) {
            if (!is.null(temp[[j]])) temp[[j]] <- tfun(temp[[j]])
        }
        for (j in c("n.event", "n.censor", "n.enter")){
            zed <- temp[[j]]
            if (!is.null(zed)) temp[[j]] <- cfun(zed)
        }
        
        if (!is.null(fit$strata)) {
            scount <- unlist(lapply(ilist1, length))
            temp$strata <- factor(rep(1:nstrat, scount), levels=1:nstrat,
                                  labels=strata.names)
        }
    }

    # An ordinary survfit object contain std(cum hazard), change scales
    if (!is.null(temp$std.err)) temp$std.err <- temp$std.err * temp$surv 
    class(temp) <- 'summary.survfit'
    temp
}
summary.survfitms <- function(object, times, censored=FALSE, 
                            scale=1, extend=FALSE, 
                            rmean= getOption("survfit.rmean"),
                            ...) {
    fit <- object
    if (!inherits(fit, 'survfitms'))
            stop("summary.survfitms can only be used for survfitms objects")

    # The print.rmean option is depreciated, it is still listened
    #   to in print.survfit, but ignored here
    if (is.null(rmean)) rmean <- "common"
    if (is.numeric(rmean)) {
        if (is.null(object$start.time)) {
            if (rmean < min(object$time)) 
                stop("Truncation point for the mean is < smallest survival")
        }
        else if (rmean < object$start.time)
            stop("Truncation point for the mean is < smallest survival")
    }
    else {
        rmean <- match.arg(rmean, c('none', 'common', 'individual'))
        if (length(rmean)==0) stop("Invalid value for rmean option")
    }

    temp <- survmean2(fit, scale=scale, rmean)  
    table <- temp$matrix  #for inclusion in the output list
    rmean.endtime <- temp$end.time

    if (!missing(times)) {
        if (!is.numeric(times)) stop ("times must be numeric")
        times <- sort(times)
    }

    # The fit$prev object is sometimes a vector and sometimes a
    #  matrix.  We calculate row indices first, and then deal
    #  with the cases at the end.
    nsurv <- length(fit$time)
    if (is.null(fit$strata)) {
        nstrat <- 1
        stemp <- rep(1L, nsurv)
        strata.names <- ""
        }
    else   {
        nstrat <- length(fit$strata)
        stemp <- rep(1:nstrat, fit$strata)
        strata.names <- names(fit$strata)
    }

    if (missing(times)) {
        # just pick off the appropriate rows of the output
        # For a survfitms object n.event is a matrix, pick off all rows with an
        #  event for some endpoint.
        if (censored) indx1 <- seq(along=fit$time)
        else indx1 <- which(rowSums(as.matrix(fit$n.event)) >0)
        indx2 <- indx1
    }
    else { 
        find2 <- function(x, vec, left.open=FALSE, ...) {
            if (!left.open) findInterval(x, vec, ...)
            else length(vec) - findInterval(-x, rev(-vec), ...)
        }
        # Process the curves one at a time, adding them to the two lists
        ilist1 <- ilist2 <- ilist3 <- vector('list', nstrat)
        newtime <- ilist1
        n <- length(stemp)
        for (i in 1:nstrat) {
            who <- (1:n)[stemp==i]  # the rows of the object for this strata
            stime <- fit$time[who]

            # First, toss any printing times that are outside our range
            if (is.null(fit$start.time)) mintime <- min(stime, 0)
            else                         mintime <- fit$start.time
            ptimes <- times[times >= mintime]

            if (!extend) {
                maxtime <- max(stime)
                ptimes <- ptimes[ptimes <= maxtime]
                }
            j <- find2(ptimes, stime) 
            ilist1[[i]] <- c(0, who)[1+ j]
            ilist2[[i]] <- c(0, who)[1+ find2(ptimes, stime, left.open=TRUE)]
            ilist3[[i]] <- j  #index within a group
            newtime[[i]] <- ptimes
        }
        indx1 <- unlist(ilist1)
        indx2 <- unlist(ilist2)

        # All of the indices (ilist1, indx1, ...) contain 0 for a time point that
        #  is prior to the first observed time in the curve.  Times that
        #  are >= to the last observed time will point to that last observed
        #  time.  Variable ilist3 contains indices that are relative to the
        #  start of a curve, all other indices point to row numbers in the
        #  entire object.
        #  
        cfun <- function(x, init=0) {  #cumulative counts over a time interval
            tlist <- vector("list", nstrat)
            if (is.matrix(x)) {
                for (i in 1:nstrat) {
                    # stemp is 1,1,1,....2,2,2,,.. to mark curves
                    x2 <- x[stemp==i,]  # all those in the group
                    j  <- c(0, ilist3[[i]])
                    tlist[[i]] <- apply(rbind(0, x2), 2, function(z) {
                        diff(cumsum(z)[1+j])})
                }
                matrix(unlist(lapply(tlist, t)), byrow=T, ncol=ncol(x))
            } 
            else {
                for (i in 1:nstrat) {
                    x2 <- x[stemp==i] 
                    j  <- c(0, ilist3[[i]])
                    tlist[[i]] <- diff(cumsum(c(0,x2))[1+j])
                }
                unlist(tlist)
            }
        }
    }

    # Create an output structure
    temp <- object
    temp$table <- table
    if (length(rmean.endtime)>0  && !is.na(rmean.endtime)) 
            temp$rmean.endtime <- rmean.endtime
    if (length(indx1)==length(fit$time) && all(indx1 == seq(along=fit$time))) {
        temp$time <- temp$time/scale
        if (!is.null(temp$strata))
            temp$strata <- factor(stemp, labels=strata.names)

    }
    else if (missing(times)) {  #default censor=FALSE case
        temp$time <- temp$time[indx1]/scale
        for (j in c("n.risk", "n.event", "n.censor", "n.enter",
                    "prev", "std.err", "cumhaz", "lower", "upper")) {
            zed <- temp[[j]]
            if (!is.null(zed)) {
                if (is.matrix(zed)) temp[[j]] <- zed[indx1,,drop=FALSE]
                else temp[[j]] <- zed[indx1]
            }
        }
        if (!is.null(temp$strata))
            temp$strata <- factor(stemp[indx1], levels=1:nstrat,
                                  labels=strata.names)
    }
    else { #times argument was given
        temp$time <- unlist(newtime)/scale
        tfun <- function(x, init=0, index= indx1) {
            if (is.matrix(x)) 
                rbind(rep(init, ncol(x)), x)[1+index,,drop=FALSE]
            else c(init, x)[1 + index]
        }
        tfun2 <- function(x, end=0, index=indx2) {
             if (is.matrix(x)) 
                rbind(x, rep(end, ncol(x)))[1+index,,drop=FALSE]
             else c(x,end)[1 + index]
        }
        temp$prev <- tfun(temp$prev, 0)
        # fix up the initial states
        if (any(indx1==0)) {
            if (nstrat==1) temp$prev[indx1==0,] <- temp$p0
            else {
                ninit <- sapply(ilist1, function(x) sum(x==0))
                zz <- rep(1:nstrat, ninit)
                temp$prev[indx1==0,] <- temp$p0[zz,]
            }
        }
        temp$n.risk <- tfun2(temp$n.risk)
        for (j in c("std.err", "cumhaz", "lower", "upper")) {
            if (!is.null(temp[[j]])) temp[[j]] <- tfun(temp[[j]])
        }
        for (j in c("n.event", "n.censor", "n.enter")){
            zed <- temp[[j]]
            if (!is.null(zed)) temp[[j]] <- cfun(zed)
        }
        
        if (!is.null(fit$strata)) {
            scount <- unlist(lapply(ilist1, length))
            temp$strata <- factor(rep(1:nstrat, scount), levels=1:nstrat,
                                  labels=strata.names)
        }
    }
    class(temp) <- "summary.survfitms"
    temp
}

print.survfitms <- function(x, scale=1,
                            rmean = getOption("survfit.rmean"), ...) {
    if (!is.null(cl<- x$call)) {
        cat("Call: ")
        dput(cl)
        cat("\n")
        }        
    omit <- x$na.action
    if (length(omit)) cat("  ", naprint(omit), "\n")

    if (is.null(rmean)) rmean <- "common"
    if (is.numeric(rmean)) {
        if (is.null(x$start.time)) {
            if (rmean < min(x$time)) 
                stop("Truncation point for the mean is < smallest survival")
        }
        else if (rmean < x$start.time)
            stop("Truncation point for the mean is < smallest survival")
    }
    else {
        rmean <- match.arg(rmean, c('none', 'common', 'individual'))
        if (length(rmean)==0) stop("Invalid value for rmean option")
    }

    temp <- survmean2(x, scale=scale, rmean)
    if (is.null(temp$end.time)) print(temp$matrix, ...)
    else {
        etime <- temp$end.time
        dd <- dimnames(temp$matrix)
        cname <- dd[[2]]
        cname[length(cname)] <- paste0(cname[length(cname)], '*')
        dd[[2]] <- cname
        dimnames(temp$matrix) <- dd
        print(temp$matrix, ...)
        if (length(etime) ==1)
             cat("   *mean time in state, restricted (max time =", 
                 format(etime, ...), ")\n")
        else cat("   *mean time in state, restricted (per curve cutoff)\n")
    }
    invisible(x)
}
survmean2 <- function(x, scale, rmean) {
    nstate <- length(x$states)  #there will always be at least 1 state
    ngrp   <- max(1, length(x$strata))
    if (ngrp >1)  {
        igrp <- rep(1:ngrp, x$strata)
        rname <- names(x$strata)
        }
    else {
        igrp <- rep(1, length(x$time))
        rname <- NULL
        }

    # The n.event matrix may not have nstate columms.  Its
    #  colnames are the first elements of states, however
    if (is.matrix(x$n.event)) {
        nc <- ncol(x$n.event)
        nevent <- tapply(x$n.event, list(rep(igrp, nc), col(x$n.event)), sum)
        dimnames(nevent) <- list(rname, x$states[1:nc])
        }
    else {
        nevent <- tapply(x$n.event, igrp, sum)
        names(nevent) <- rname
        }

    outmat <- matrix(0., nrow=nstate*ngrp , ncol=2)
    outmat[,1] <- rep(x$n, nstate)
    outmat[1:length(nevent), 2] <- c(nevent)
  
    if (ngrp >1) 
        rowname <- c(outer(rname, x$states, paste, sep=", "))
    else rowname <- x$states

    # Caculate the mean time in each state
    if (rmean != "none") {
        if (is.numeric(rmean)) maxtime <- rep(rmean, ngrp)
        else if (rmean=="common") maxtime <- rep(max(x$time), ngrp)
        else maxtime <- tapply(x$time, igrp, max)
    
        meantime <- matrix(0., ngrp, nstate)
        p0 <- matrix(x$p0, nrow=ngrp)  #in case there is only one row
        for (i in 1:ngrp) {
            if (is.matrix(x$prev))
                temp <- rbind(p0[i,], x$prev[igrp==i,, drop=FALSE])
            else temp <- matrix(c(p0[i,], x$prev[igrp==i]), ncol=1)

            if (is.null(x$start.time)) tt <- c(0, x$time[igrp==i])
            else tt <- c(x$start.time, x$time[igrp==i])

            # Now cut it off at maxtime
            delta <- diff(c(tt[tt<maxtime[i]], maxtime[i]))
            if (length(delta) > nrow(temp)) delta <- delta[1:nrow(temp)]
            if (length(delta) < nrow(temp))
                delta <- c(delta, rep(0, nrow(temp) - length(delta)))
            meantime[i,] <- colSums(delta*temp)
        }

        outmat <- cbind(outmat, c(meantime)/scale)
        cname <- c("n", "nevent", "mean")
        # report back a single time, if there is only one
        if (all(maxtime == maxtime[1])) maxtime <- maxtime[1]
    }
    else cname <- c("n", "nevent")
    dimnames(outmat) <- list(rowname, cname)

    if (rmean=='none') list(matrix=outmat)
    else list(matrix=outmat, end.time=maxtime/scale)
}
"[.survfitms" <- function(x, ..., drop=TRUE) {
    nmatch <- function(indx, target) { 
        # This function lets R worry about character, negative, or logical subscripts
        #  It always returns a set of positive integer indices
        temp <- 1:length(target)
        names(temp) <- target
        temp[indx]
    }
        
    if (missing(..1)) i<- NULL  else i <- sort(..1)
    if (missing(..2)) j<- NULL  else j <- ..2
    n <- length(x$time)

    if (is.null(x$strata) && is.matrix(x$prev)) {
        # No strata, but a matrix of prevalence values
        #  In this case, allow them to use a single i subscript as well
        if (is.null(j) && !is.null(i)) {
            j <- i
            i <- NULL
        }
    }
    if (is.null(i)) {
        i2 <- 1:n
        if (is.null(strata)) i <- 1
        else i <- seq(along=strata)
    }
    else {
        if (is.null(x$strata) && (length(i) > 1 || i != 1))
            stop("subscript out of bounds")
        indx <- nmatch(i, names(x$strata)) #strata to keep
        if (any(is.na(indx))) 
            stop(paste("strata", 
                       paste(i[is.na(indx)], collapse=' '),
                       'not matched'))
        # Now, i may not be in order: a user has curve[3:2] to reorder 
        #  a plot.  Hence the "unlist(lapply(" construct which will reorder
        #  the data in the curves
        temp <- rep(1:length(x$strata), x$strata)
        keep <- unlist(lapply(i, function(x) which(temp==x)))

        if (length(i) <=1 && drop) x$strata <- NULL
        else               x$strata  <- x$strata[indx]
        i2 <- keep
    }

    if (!is.null(j)) {
        indx <- nmatch(j, x$states)
        if (any(is.na(indx)))
            stop("subscript out of bounds", j[is.na(indx)])
        else j <- as.vector(indx)
    }

    if (length(i2) ==1 && !is.null(j) && missing(drop)) drop <- FALSE
 
    # all the elements that can have "nstate" elements or columns
    #  The n.event variable can have fewer
    temp <- c("states", "n.risk", "n.event", "n.censor", "prev", 
              "cumhaz", "std.err", "lower", "upper")
    sfun <- function(z) {
        if (is.null(j)) {
            if (is.array(z)) {
                if (length(dim(z)) > 2) z[,,i2, drop=drop]  
                else z[i2,,drop=drop]
            }
            else z
        }
        else {
            if (is.array(z)) {
                if (length(dim(z)) > 2) z[j,j,i2, drop=drop]  
                else z[i2,j, drop=drop]
            }
            else z[j]
        }
    }
    for (k in temp) x[[k]] <- sfun(x[[k]])
    x$n <- x$n[i]
    x$time <- x$time[i2]
    x$transitions <- NULL  # this is incorrect after subscripting

    if (is.null(j)) x$p0<- x$p0[i,]
    else x$p0 <- x$p0[i,j]
    
    x
}
