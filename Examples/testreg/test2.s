#
# Look at the loglik for the Peterson data, in the region of the
#  initial values
# 
x0 <- seq(2.05, 2.1, length=15)
sig<- seq(.30, .35, length=15)
x1 <- c(0.696208,0.703780 ,-7.838072 ,-0.542997,0.480374)

plog <- matrix(0,15,15) 
for (i in 1:15) {
    for (j in 1:15) {
	fit <- survreg(Surv(time, status)~factor(grp), peterson,
		       init=c(x0[i], x1, sig[j]),
		       control=list(maxiter=0))
	plog[i,j] <- fit$log[1]
	}
    }

#cline <- -c(34:38, 40, 45, 50)
contour(x0, sig, plog, levels=cline)

