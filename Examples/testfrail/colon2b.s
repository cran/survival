#
# A plot of the colon data
#
cgam <- cdf
log0 <- tfit.1$loglik[1]
tempcoef <- matrix(0, 20, 4)
for (i in 1:20) {
    tfit <- get(paste('tfit',i, sep='.'))
    cgam[i] <- tfit$loglik[2] - log0
    tempcoef[i,] <- tfit$coef
    }

matplot(ctheta, cbind(cgam, 2*((cpl- log0)-cdf) - 400),
	type='l', lty=1:3,
	xlab="Theta", ylab="Corrected Partial Liklihood")
legend(4,300, c("marginal Gamma frailty", "Gamma AIC -400"), 
              lty=1:3, col=1:3, bty='n')
title(main="Colon Cancer Data")
