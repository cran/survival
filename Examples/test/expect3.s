#
# Test out expected survival, when the parent pop is another Cox model
#
fit <- coxph(Surv(time, status) ~x, test1, method='breslow')

dummy <- data.frame(time=c(.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
		    status=c(1,0,1,0,1,0,1,1,1), x=(-4:4)/2)

efit <- survexp(time ~ ratetable(x=x), dummy, ratetable=fit, cohort=F)

#
# Now, compare to the true answer, which is known to us
#
ss <- exp(fit$coef)
haz <- c( 1/(3*ss+3), 2/(ss+3), 1) #truth at time 0,1,2,4+
chaz <- cumsum(c(0,haz))
chaz2 <- chaz[c(1,2,2,3,3,3,3,4,4)]

risk <- exp(fit$coef*dummy$x)
efit2 <- exp(-risk*chaz2)

all.equal(as.vector(efit), as.vector(efit2))  #ignore mismatched name attrib

#
# Now test the direct-adjusted curve (Ederer)
#
efit <- survexp( ~ ratetable(x=x), dummy, ratetable=fit, se=F)
direct <- survfit(fit, newdata=dummy)$surv

chaz <- chaz[-1]                  #drop time 0
d2 <- exp(outer(-chaz, risk))
all.equal(as.vector(direct), as.vector(d2))   #this tests survfit

all.equal(as.vector(efit$surv), as.vector(apply(direct,1,mean)))  #direct


#
# Now test out the Hakulinen method (Bonsel's method)
#  By construction, we have a large correlation between x and censoring
#
# In theory, hak1 and hak2 would be the same.  In practice, like a KM and
#   F-H, they differ when n is small.
#
efit <- survexp( time ~ ratetable(x=x), dummy, ratetable=fit, se=F)

surv  <- wt <- rep(1,9)
tt <- c(1,2,4)
hak1 <- hak2 <- NULL
for (i in 1:3) {
    wt[dummy$time < tt[i]]  <- 0
    hak1 <- c(hak1,  exp(-sum(haz[i]*risk*surv*wt)/sum(surv*wt)))
    hak2 <- c(hak2,  sum(exp(-haz[i]*risk)*surv*wt)/sum(surv*wt))
    surv <- surv * exp(-haz[i]*risk)
    }

all.equal(as.vector(efit$surv), as.vector(cumprod(hak2)))

#
#  Now do the conditional estimate
#
efit <- survexp( time ~ ratetable(x=x), dummy, ratetable=fit, se=F,
			conditional=T)
wt <- rep(1,9)
cond <- NULL
for (i in 1:3) {
    wt[dummy$time < tt[i]]  <- 0
    cond <- c(cond,  exp(-sum(haz[i]*risk*wt)/sum(wt)))
    }

all.equal(as.vector(efit$surv), as.vector(cumprod(cond)))

rm(wt, cond, efit, tt, surv, hak1, hak2)
rm(fit, dummy, ss, efit2, chaz, chaz2, risk)
rm(d2, direct)
