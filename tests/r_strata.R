options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Test out the strata capabilities
#
tol <- survreg.control()$rel.tolerance
aeq <- function(x,y,...) all.equal(as.vector(x), as.vector(y), ...)

# intercept only models
fit1 <- survreg(Surv(time, status) ~ strata(sex), lung)
fit2 <- survreg(Surv(time, status) ~ strata(sex) + sex, lung)
fit3a<- survreg(Surv(time,status) ~1, lung, subset=(sex==1))
fit3b<- survreg(Surv(time,status) ~1, lung, subset=(sex==2))

fit1
fit2
aeq(fit2$scale, c(fit3a$scale, fit3b$scale), tolerance=tol)
aeq(fit2$loglik[2], (fit3a$loglik + fit3b$loglik)[2], tolerance=tol)
aeq(fit2$coef[1] + 1:2*fit2$coef[2], c(fit3a$coef, fit3b$coef), tolerance=tol)

#penalized models
fit1 <- survreg(Surv(time, status) ~ pspline(age, theta=.92)+strata(sex), lung)
fit2 <- survreg(Surv(time, status) ~  pspline(age, theta=.92)+ 
		strata(sex) + sex, lung)
fit1
fit2

age1 <- ifelse(lung$sex==1, lung$age, mean(lung$age))
age2 <- ifelse(lung$sex==2, lung$age, mean(lung$age))
fit3 <- survreg(Surv(time,status) ~ pspline(age1, theta=.92) +
		pspline(age2, theta=.95) + sex + strata(sex), lung,
		rel.tol=1e-6)
fit3a<- survreg(Surv(time,status) ~pspline(age, theta=.92), lung, 
		    subset=(sex==1))
fit3b<- survreg(Surv(time,status) ~pspline(age, theta=.95), lung, 
		     subset=(sex==2))

# relax the tolerance a little, since the above has lots of parameters
#  I still don't exactly match the second group, but very close
aeq(fit3$scale, c(fit3a$scale, fit3b$scale), tolerance=tol*10)
aeq(fit3$loglik[2], (fit3a$loglik + fit3b$loglik)[2], tolerance=tol*10)
pred <- predict(fit3)
aeq(pred[lung$sex==1] , predict(fit3a), tolerance=tol*10)
aeq(pred[lung$sex==2],  predict(fit3b), tolerance=tol*10)




