#
# Fit the kidney data using AIC
#

# gamma, corrected aic
coxph(Surv(time, status) ~ age + sex + frailty(id, method='aic', caic=T), 
      kidney)

coxph(Surv(time, status) ~ age + sex + frailty(id, dist='t'), kidney)
coxph(Surv(time, status) ~ age + sex + frailty(id, dist='gauss', method='aic',
					       caic=T), kidney)


# uncorrected aic
coxph(Surv(time, status) ~ age + sex + frailty(id, method='aic', caic=F), 
      kidney)

coxph(Surv(time, status) ~ age + sex + frailty(id, dist='t', caic=F), kidney)
