#
# Test the coxph program on the Ovarian data
#

xx _ order(ovarian$futime)   #put data in same order as SAS green book
temp <- ovarian[xx,]
attach(temp)

# List the data
temp

summary(survfit(Surv(futime, fustat)), censor=T)

# Various models
coxph(Surv(futime, fustat)~ age)
coxph(Surv(futime, fustat)~ resid.ds)
coxph(Surv(futime, fustat)~ rx)
coxph(Surv(futime, fustat)~ ecog.ps)

coxph(Surv(futime, fustat)~ resid.ds + rx + ecog.ps)
coxph(Surv(futime, fustat)~ age + rx + ecog.ps)
coxph(Surv(futime, fustat)~ age + resid.ds + ecog.ps)
coxph(Surv(futime, fustat)~ age + resid.ds + rx)

# Residuals
fit <- coxph(Surv(futime, fustat)~ age + resid.ds + rx + ecog.ps )
resid(fit)
resid(fit, 'dev')
resid(fit, 'scor')
resid(fit, 'scho')

fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps + strata(rx))
summary(fit)
summary(survfit(fit))
sfit <- survfit(fit, list(age=c(30,70), ecog.ps=c(2,3)))  #two columns
sfit
summary(sfit)
detach(w=2)


# Test the robust=T option of coxph
fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps + rx, ovarian, robust=T)
rr <- resid(fit, type='dfbeta')
all.equal(as.vector(t(rr) %*% rr), as.vector(fit$var))
