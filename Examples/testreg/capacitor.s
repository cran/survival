#
# This data set caused problems for Splus 3.4 due to a mistake in
#  my initial value code.  Data courtesy Bob Treder at Statsci
#
capacitor <- read.table('data.capacitor', row.names=1,
			col.names=c('', 'days', 'event', 'voltage'))

fitig <- survreg(Surv(days, event)~voltage, 
	dist = "gaussian", data = capacitor)
summary(fitig)

fitix <- survreg(Surv(days, event)~voltage, 
	dist = "extreme", data = capacitor)
summary(fitix)

fitil <- survreg(Surv(days, event)~voltage, 
	dist = "logistic", data = capacitor)
summary(fitil)

rm(fitil, fitig, fitix)
