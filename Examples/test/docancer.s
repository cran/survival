#
# Test out all of the routines on a more complex data set
#
temp <- survfit(Surv(time, status) ~ ph.ecog, cancer)
summary(temp, times=c(30*1:11, 365*1:3))
print(temp[2:3])

temp <- survfit(Surv(time, status), cancer, type='fleming',
		   conf.int=.9, conf.type='log-log', error='tsiatis')
summary(temp, times=30 *1:5)

temp <- survdiff(Surv(time, status) ~ inst, cancer, rho=.5)
print(temp, digits=6)

temp <- coxph(Surv(time, status) ~ ph.ecog + ph.karno + pat.karno + wt.loss 
	      + sex + age + meal.cal + strata(inst), cancer)
summary(temp)
cox.zph(temp)
cox.zph(temp, transform='identity')

coxph(Surv(rep(0,length(time)), time, status) ~ ph.ecog + ph.karno + pat.karno
		+ wt.loss + sex + age + meal.cal + strata(inst), cancer)
