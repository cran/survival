#lfit1 <- censorReg(censor(time, status) ~ age + ph.ecog + strata(sex),lung)
lfit2 <- survreg(Surv(time, status) ~ age + ph.ecog + strata(sex), lung)
lfit3 <- survreg(Surv(time, status) ~ sex + (age+ph.ecog)*strata(sex), lung)

lfit4 <-  survreg(Surv(time, status) ~ age + ph.ecog , lung,
		  subset=(sex==1))
lfit5 <- survreg(Surv(time, status) ~ age + ph.ecog , lung,
		  subset=(sex==2))

aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
#aeq(lfit4$coef, lfit1[[1]]$coef)
#aeq(lfit4$scale, lfit1[[1]]$scale)
aeq(c(lfit4$scale, lfit5$scale), lfit3$scale )
aeq(c(lfit4$scale, lfit5$scale), sapply(lfit1, function(x) x$scale))

