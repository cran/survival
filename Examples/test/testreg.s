#
# Run a test that can be verified using SAS's LIFEREG
#
fit1w <- survreg(Surv(time, status) ~x, test1, dist='weibull')
fit1w
summary(fit1w)

fit1e <- survreg(Surv(time, status) ~x, test1, dist='exp')
fit1e
summary(fit1e)

fit1l <- survreg(Surv(time, status) ~x, test1, dist='loglogistic')
fit1l
summary(fit1l)

fit1g <- survreg(Surv(time, status) ~x, test1, dist='lognormal')
summary(fit1g)
#
#  Do a test with the ovarian data
#
fitfw <- survreg.old(Surv(futime, fustat) ~ age + ecog.ps, ovarian,
	link='log', dist='extreme')
fitfw

fitfl <- survreg.old(Surv(futime, fustat) ~ age + ecog.ps, ovarian,
	link='log', dist='logistic')
fitfl

flem2 <- scan("fleming.data2", what=list(ltime=0, rtime=0))
flsurv<- Surv(flem2$ltime, flem2$rtime, type='interval2')

fitfw2 <- survreg(flsurv ~ age + ecog.ps, ovarian, dist='weibull')
summary(fitfw2)

fitfl2 <- survreg(flsurv ~ age + ecog.ps, ovarian, dist='loglogistic')
summary(fitfl2)

fitfg2 <- survreg(flsurv ~ age + ecog.ps, ovarian, dist='lognormal')
summary(fitfg2)
