#
# Compute the K-M for the Turnbull data
#      via a slow EM calculation
#

emsurv <- function(time, status, wt, verbose=T) {
    left.cen <- (status==2)
    if (!any(left.cen)) stop("No left censored data!")
    if (!any(status==1))stop("Must have some exact death times")

    tempy <- Surv(time[!left.cen], status[!left.cen])
    ww <- wt[!left.cen]
    tempx <- factor(rep(1, sum(!left.cen)))
    tfit <- survfit.km(tempx, tempy, casewt=ww)
    if (verbose)
       cat("Iteration 0, survival=", format(round(tfit$surv[tfit$n.event>0],3)),
		       "\n")

    stimes <- tfit$time[tfit$n.event>0]
    ltime <- time[left.cen]
    lwt   <- wt[left.cen]
    tempx <- factor(rep(1, length(stimes) + sum(!left.cen)))
    tempy <- Surv(c(time[!left.cen], stimes),
		  c(status[!left.cen], rep(1, length(stimes))))
    for (iter in 1:4) {
	wt2 <- stimes*0
	ssurv <- tfit$surv[tfit$n.event>0]
	sjump <- diff(c(1, ssurv))
	for (j in 1:(length(ltime))) {
	    k <- sum(ltime[j]>=stimes)   #index of the death time
	    if (k==0)
		stop("Left censored observation before the first death")
	    wt2[1:k] <- wt2[1:k] + lwt[j]*sjump[1:k] /(ssurv[k]-1)
	    }
	tfit <- survfit.km(tempx, tempy, casewt=c(ww, wt2))
	if (verbose) {
	   cat("Iteration", iter, "survival=",
		 format(round(tfit$surv[tfit$n.event>0],3)),  "\n")
	   cat("             weights=", format(round(wt2,3)), "\n")
	   }
	}
    survfit(tempy ~ tempx, weights=c(ww, wt2))
    }

temp <-emsurv(turnbull$time, turnbull$status, turnbull$n)
print(summary(temp))
