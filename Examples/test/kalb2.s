jasa1 _ stan1(jasa)
attach(jasa1)
# Now fit the 6 models found in Kalbfleisch and Prentice, p139
options(contrasts=c("contr.treatment", "contr.treatment"))
sfit.1 _ coxph(Surv(start, stop, event)~ (age + surgery)*transplant,
				method='breslow')
sfit.2 _ coxph(Surv(start, stop, event)~ year*transplant,
				method='breslow')
sfit.3 _ coxph(Surv(start, stop, event)~ (age + year)*transplant,
				method='breslow')
sfit.4 _ coxph(Surv(start, stop, event)~ (year +surgery) *transplant,
				method='breslow')
sfit.5 _ coxph(Surv(start, stop, event)~ (age + surgery)*transplant + year ,
				method='breslow')
sfit.6 _ coxph(Surv(start, stop, event)~ age*transplant + surgery + year,
				method='breslow')

summary(sfit.1)
sfit.2
summary(sfit.3)
sfit.4
sfit.5
sfit.6

