#
# Simplest printout to debug frailty models
#
fit1 <- survreg(Surv(time, status) ~ frailty(group, theta=.1, sparse=F), 
		data=leukemia)

fit2 <-survreg(Surv(time, status) ~ frailty(group, theta=.1, sparse=F),
	       data=leukemia, debug=2, init=c(fit1$coef, log(fit1$scale)),
	       iter=0)

fit3 <-survreg(Surv(time, status) ~ frailty(group, theta=.1, sparse=T),
	       data=leukemia, debug=2, 
	       init=c(fit1$coef[c(2,3,1)], log(fit1$scale)), iter=0)

