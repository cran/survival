#
#  Check out the various residuals under an Efron approximation
#
fit0 <- coxph(Surv(time, status)~ x, test1, iter=0)
fit  <- coxph(Surv(time, status) ~x, test1)
fit0b <- coxph(Surv(start, stop, status) ~ x, test1b, iter=0)
fitb  <- coxph(Surv(start, stop, status) ~x, test1b)
fitc  <- coxph(Surv(time, status) ~ offset(fit$coef*x), test1)
fitd  <- coxph(Surv(start, stop, status) ~ offset(fit$coef*x), test1b)

resid(fit0)
resid(fit0b, collapse=test1b$id)
resid(fit)
resid(fitb, collapse=test1b$id)
resid(fitc)
resid(fitd, collapse=test1b$id)
all.equal(resid(fitc), resid(fit))

resid(fit0, type='score')
resid(fit0b, type='score', collapse=test1b$id)
resid(fit, type='score')
resid(fitb, type='score', collapse=test1b$id)

resid(fit0, type='scho')
resid(fit0b, type='scho', collapse=test1b$id)
resid(fit, type='scho')
resid(fitb, type='scho', collapse=test1b$id)

summary(survfit(fit0, list(x=0)))
summary(survfit(fit, list(x=0)))
summary(survfit(fitb,list(x=0)))
summary(survfit(fit))
