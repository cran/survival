#
#   Do the test on the simple agreg data set
#
fit <-coxph(Surv(start, stop, event)~ x, test2, method='breslow')
fit
fit0 <-coxph(Surv(start, stop, event)~ x, test2, iter=0)
fit0$coef
coxph(Surv(start, stop, event)~ x, test2, iter=1, method='breslow')$coef
coxph(Surv(start, stop, event)~ x, test2, iter=2, method='breslow')$coef
coxph(Surv(start, stop, event)~ x, test2, iter=3, method='breslow')$coef

coxph(Surv(start, stop, event) ~ x, test2, method='efron')
coxph(Surv(start, stop, event) ~ x, test2, method='exact')

resid(fit0)
resid(fit0, 'scor')
resid(fit0, 'scho')

resid(fit)
resid(fit, 'scor')
resid(fit, 'scho')

resid(coxph(Surv(start, stop, event)~ x, test2, iter=0, init=log(2)), 'score')

sfit <-survfit(fit)
sfit
summary(sfit)

# make a doubled data set
temp <- rbind(test2, test2)
temp <- data.frame(temp, x2=c(test2$x, test2$x^2), 
		         ss=c(rep(0, nrow(test2)), rep(1, nrow(test2))))
fitx <- coxph(Surv(start, stop, event) ~ x2 * strata(ss), data=temp,
	      method='breslow')

sfit <- survfit(fitx, c(fitx$means[1], 0) )
sfit
summary(sfit)
#
# Even though everyone in strata 1 has x2==0, I won't get the same survival
#  curve above if survfit is called without forcing predicted x2 to be
#  zero-- otherwise I am asking for a different baseline than the
#  simple model did.  In this particular case coef[2] is nearly zero, so
#  the curves are the same, but the variances differ.
#

# This mimics 'byhand3' and the documentation
fit <- coxph(Surv(start, stop, event) ~x, test2, method='breslow')
tdata <- data.frame( start=c(0,20, 6,0,5),
		     stop =c(5,23,10,5,15),
		     event=rep(0,5),
		     x=c(0,0,1,0,2) )

temp <- survfit(fit, tdata, individual=T)
temp2 <- byhand3(fit$coef)
all.equal(temp$surv, temp2$surv)
all.equal(temp2$std, temp$std.err)

temp2
