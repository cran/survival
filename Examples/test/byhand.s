#
# Compute some answers "by hand"
#     For test1 and test2 data
byhand1 <- function(coef, newx) {
    s <- exp(coef)
    loglik <- 2*coef - (log(3*s+3) + 2*log(s+3))
    u <- (6 + 3*s - s^2) / ((s+1)*(s+3))
    imat <- s/(s+1)^2 + 6*s/(s+3)^2

    x <- c(1,1,1,0,0,0)
    status <- c(1,0,1,1,0,1)
    xbar <- c(s/(s+1), s/(s+3), 0, 0)
    haz <- c(1/(3*s+3), 2/(s+3), 0, 1 )
    ties <- c(1,1,2,2,3,4)
    wt <- c(s,s,s,1,1,1)
    mart <- c(1,0,1,1,0,1) -  wt* (cumsum(haz))[ties]

    score <- rep(6,0)
    for (i in 1:6) {
	j <- ties[i]
	score[i] <-  -wt[i]*(cumsum((x[i]-xbar) * haz))[j]
	score[i] <- score[i] + wt[i]* ((x[i]-xbar)*status[i])[j]
	}

    scho <- c(1/(s+1), (3-s)/(3+s), 0)

    surv  <- exp(-cumsum(haz)* exp(coef*newx))
    varhaz.g <- cumsum(c(1/(3*s+3)^2, 2/(s+3)^2, 0, 1 ))

    varhaz.d <- cumsum((newx-xbar) * haz)

    varhaz <- (varhaz.g + varhaz.d0^2/ imat) * exp(2*coef*newx)

    names(xbar) <- names(haz) <- 1:4
    names(surv) <- names(varhaz) <- 1:4
    list(loglik=loglik, u=u, imat=imat, xbar=xbar, haz=haz,
	     mart=mart,  score=score,
		scho=scho, surv=surv, var=varhaz,
		varhaz.g=varhaz.g, varhaz.d=varhaz.d)
    }


byhand2 <- function(coef, newx) {
    s <- exp(coef)

    loglik <- 4*coef - log(s+1) - log(s+2) - 3*log(3*s+2) - 2*log(3*s+1)
    u <- 1/(s+1) +  1/(3*s+1) + 4/(3*s+2) -
		 ( s/(s+2) +3*s/(3*s+2) + 3*s/(3*s+1))
    imat <- s/(s+1)^2 + 2*s/(s+2)^2 + 6*s/(3*s+2)^2 +
	    3*s/(3*s+1)^2 + 3*s/(3*s+1)^2 + 12*s/(3*s+2)^2

    hazard <-c( 1/(s+1), 1/(s+2), 1/(3*s+2), 1/(3*s+1), 1/(3*s+1), 2/(3*s+2) )
    xbar <- c(s/(s+1), s/(s+2), 3*s/(3*s+2), 3*s/(3*s+1), 3*s/(3*s+1),
		3*s/(3*s+2))

    var.g <- cumsum(hazard*hazard /c(1,1,1,1,1,2))
    var.d <- cumsum( (xbar-newx)*hazard)

    surv <- exp(-cumsum(hazard) * exp(coef*newx))
    varhaz <- (var.g + var.d^2/imat)* exp(2*coef*newx)

    list(loglik=loglik, u=u, imat=imat, hazard=hazard,
	    xbar=xbar, surv=surv, varhaz=varhaz, var.g=var.g, var.d=var.d)
    }

byhand3 <- function(coef) {
    #Hard coded -- what is found in the Agreg.3 comments file
    s <- as.vector(exp(coef))  #kill the names attr

    imat <- s/(s+1)^2 + 2*s/(s+2)^2 + 6*s/(3*s+2)^2 +
	    3*s/(3*s+1)^2 + 3*s/(3*s+1)^2 + 12*s/(3*s+2)^2

    hazard <-c( 1/(s+1), 1/(s+2), 1/(3*s+2), 1/(3*s+1), 1/(3*s+1), 2/(3*s+2) )
    xbar <- c(s/(s+1), s/(s+2), 3*s/(3*s+2), 3*s/(3*s+1), 3*s/(3*s+1),
		3*s/(3*s+2))


    newx <- c(0,0,1,1,1,0,0, 2,2,2,2)
    wt <- exp(coef*newx)
    indx <- c(1,2,4,5,6,1,2,3,4,5,6)

    var.g <- hazard*hazard /c(1,1,1,1,1,2)
    surv <- exp(-cumsum(hazard[indx]*wt))
    var1 <- cumsum(var.g[indx]*wt*wt)
    d    <- cumsum( (xbar[indx] - newx)* hazard[indx] * wt)
    var2 <- d^2/imat
    names(surv) <- names(var1) <-names(var2) <- NULL
    list(time= c(2,3,7,8,9,12,13,16,17,18,19),
	 surv=surv, std= sqrt(var1 + var2), var.g=var1, var.d=var2, d=d,
		hazard=hazard, wt=wt)
    }

