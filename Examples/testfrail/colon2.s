#
#  Look at the Akaike information criteria for the colon data
#

ctheta <- 1:20/2
cdf <- cpl <- cpen <- ctheta
initials <- rep(0, 4+ length(unique(colon$id)))
for (i in seq(along=ctheta)) {
    tfit <- coxph(Surv(time, status) ~ rx + extent + node4 + 
		    + frailty(id, theta=ctheta[i]) + strata(etype), colon,
		  init=initials, iter=30)
    cat(i)
    assign(paste("tfit", i, sep='.'), tfit)
    cdf[i] <- sum(tfit$df)
    cpl[i] <- tfit$plik
    cpen[i] <- tfit$penalty
    initials <- c(tfit$coef, tfit$frail)
    }

