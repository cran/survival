fit1 <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian)
fit2 <- censorReg(censor(futime, fustat) ~ age + ecog.ps, ovarian)
fit3 <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian,
		iter=0, init=c(fit2$coef,   log(fit2$scale)))
fit4 <- survreg(Surv(log(futime), fustat) ~age + ecog.ps, ovarian,
		dist='extreme',
		iter=0, init=c(fit2$coef,   log(fit2$scale)))
		

aeq(resid(fit2, type='working')[,1], resid(fit3, type='working'))
aeq(resid(fit2, type='response')[,1], resid(fit3, type='response'))

temp <- sign(resid(fit3, type='working'))
aeq(resid(fit2, type='deviance')[,1], temp*abs(resid(fit3, type='deviance')))
