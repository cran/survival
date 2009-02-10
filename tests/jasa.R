options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

expect <- survexp(futime ~ ratetable(age=(accept.dt - birth.dt), sex=1,
		year=accept.dt, race='white'), jasa, cohort=F, 
                  ratetable=survexp.usr)

survdiff(Surv(jasa$futime, jasa$fustat) ~ offset(expect))
# Now fit the 6 models found in Kalbfleisch and Prentice, p139
sfit.1 <- coxph(Surv(start, stop, event)~ (age + surgery)*transplant,
				jasa1, method='breslow')
sfit.2 <- coxph(Surv(start, stop, event)~ year*transplant,
				jasa1, method='breslow')
sfit.3 <- coxph(Surv(start, stop, event)~ (age + year)*transplant,
				jasa1, method='breslow')
sfit.4 <- coxph(Surv(start, stop, event)~ (year +surgery) *transplant,
				jasa1, method='breslow')
sfit.5 <- coxph(Surv(start, stop, event)~ (age + surgery)*transplant + year ,
				jasa1, method='breslow')
sfit.6 <- coxph(Surv(start, stop, event)~ age*transplant + surgery + year,
				jasa1, method='breslow')

summary(sfit.1)
sfit.2
summary(sfit.3)
sfit.4
sfit.5
sfit.6

# Survival curve for the "average" subject
summary(survfit(sfit.1))

# Survival curve for a subject of age 50, with prior surgery, tx at 6 months
data <- data.frame(start=c(0,183), stop=c(183,3*365), event=c(1,1),
		   age=c(50,50),  surgery=c(1,1), transplant=c(0,1))
summary(survfit(sfit.1, data, individual=T))

# These should all give the same answer
j.age <- jasa$age -48
fit1 <- coxph(Surv(futime, fustat) ~ j.age, data=jasa)
fit2 <- coxph(Surv(futime, fustat) ~ j.age, jasa, init=fit1$coef, iter=0)
fit3 <- coxph(Surv(start, stop, event) ~ age, jasa1)
fit4 <- coxph(Surv(start, stop, event) ~ offset((age-fit3$means)*fit1$coef),
              jasa1)
s1 <- survfit(fit1, fit3$means)
s2 <- survfit(fit2, fit3$means)
s3 <- survfit(fit3)
s4 <- survfit(fit4)

all.equal(s1$surv, s2$surv)
all.equal(s1$surv, s3$surv)
all.equal(s1$surv, s4$surv)

# Still the same answer, fit multiple strata at once
#  Strata 1 has independent coefs of strata 2, so putting in
#    the other data should not affect it
attach(jasa1)
ll <- length(start)
ss <- rep(0:1, c(ll,ll))
tdata <- data.frame(start=rep(start,2), stop=rep(stop,2),
		    event=rep(event,2), ss=ss, age=rep(age,2),
		    age2 = (rep(age,2))^2 * ss)
fit <- coxph(Surv(start, stop, event) ~ age*strata(ss) + age2, tdata)
#  Above replaced these 2 lines, which kill Splus5 as of 8/98
#    Something with data frames, I expect.
#fit <- coxph(Surv(rep(start,2), rep(stop,2), rep(event,2)) ~
#			rep(age,2)*strata(ss) + I(rep(age,2)^2*ss) )
all.equal(fit$coef[1], fit3$coef)
s5 <- survfit(fit, c(fit3$means, 0,0))
all.equal(s5$surv[1:(s5$strata[1])],  s3$surv)
detach("jasa1")

rm(s1, s2, s3, s4, s5, tdata, ll, ss, data)
rm(fit1, fit2, fit3, fit4, fit, j.age)
rm(sfit.1, sfit.2, sfit.3, sfit.4, sfit.5, sfit.6)
