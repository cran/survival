#
# Draw a contour map of the test1 data, to show the non SPD area of
#  the loglik
#

x0 <- seq(1, 3, length=20)
x1 <- -.745
sig<- seq(-2, 0, length=20)
lmat <- matrix(0, length(x0), length(sig))
for (i in 1:length(x0)) {
    for (j in 1:length(sig)) {
        fit<- survreg(Surv(time, status)~x, test1, 
		      init=c(x0[i], x1, sig[j]), control=list(maxiter=0))
	lmat[i,j] <- fit$loglik[1]
	}
    }

x1 <- -.735
lmat2 <- matrix(0, length(x0), length(sig))
for (i in 1:length(x0)) {
    for (j in 1:length(sig)) {
        fit<- survreg(Surv(time, status)~x, test1, 
		      init=c(x0[i], x1, sig[j]), control=list(maxiter=0))
	lmat2[i,j] <- fit$loglik[1]
	}
    }

clev <- -1*c(2.5,3:6, 8, 10, 15, 25)
contour(x0, sig, lmat, levels=clev)
contour(x0, sig, lmat2,levels=clev)
