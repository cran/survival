# Tests of the weighted Cox model
#
# Similar data set to test1, but add weights,
#                                    a double-death/censor tied time
#                                    a censored last subject
# The latter two are cases covered only feebly elsewhere.
testw1 <- data.frame(time=  c(1,1,2,2,2,2,3,4,5),
		    status= c(1,0,1,1,1,0,0,1,0),
		    x=      c(2,0,1,1,0,1,0,1,0),
		    wt =    c(1,2,3,4,3,2,1,2,1))

fit0 <- coxph(Surv(time, status) ~x, testw1, weights=wt,
		    method='breslow', iter=0)
fit  <- coxph(Surv(time, status) ~x, testw1, weights=wt, method='breslow')
fit0
summary(fit)
resid(fit0, type='mart')
resid(fit0, type='score')
resid(fit0, type='scho')
fit0 <- coxph(Surv(time, status) ~x,testw1, weights=wt, iter=0)
resid(fit0, 'mart')
resid(coxph(Surv(time, status) ~1, testw1, weights=wt))  #Null model

resid(fit, type='mart')
resid(fit, type='score')
resid(fit, type='scho')

fit  <- coxph(Surv(time, status) ~x, testw1, weights=wt, method='efron')
fit
resid(fit, type='mart')
resid(fit, type='score')
resid(fit, type='scho')
