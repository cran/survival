# Tests of the weighted Cox model, AG form of the data
#   Same solution as doweight1.s
#
testw2 <- data.frame(id  =  c( 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 7, 8, 8, 9),
		     begin= c( 0, 5, 0, 0,10,15, 0, 0,14, 0, 0, 0,23, 0),
		     time=  c( 5,10,10,10,15,20,20,14,20,20,30,23,40,50),
		    status= c( 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0),
		    x=      c( 2, 2, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0),
		    wt =    c( 1, 1, 2, 3, 3, 3, 4, 3, 3, 2, 1, 2, 2, 1))

fit0 <- coxph(Surv(begin,time, status) ~x, testw2, weights=wt,
		    method='breslow', iter=0)
fit  <- coxph(Surv(begin,time, status) ~x, testw2, weights=wt, method='breslow')
fit0
summary(fit)
resid(fit0, type='mart', collapse=testw2$id)
resid(fit0, type='score', collapse=testw2$id)
resid(fit0, type='scho')

resid(fit, type='mart', collapse=testw2$id)
resid(fit, type='score', collapse=testw2$id)
resid(fit, type='scho')
fit0 <- coxph(Surv(begin, time, status) ~x,testw2, weights=wt, iter=0)
resid(fit0, 'mart', collapse=testw2$id)
resid(coxph(Surv(begin, time, status) ~1, testw2, weights=wt)
		      , collapse=testw2$id)  #Null model

fit  <- coxph(Surv(begin,time, status) ~x, testw2, weights=wt, method='efron')
fit
resid(fit, type='mart', collapse=testw2$id)
resid(fit, type='score', collapse=testw2$id)
resid(fit, type='scho')
