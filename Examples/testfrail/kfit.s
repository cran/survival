# From:	McGilchrist and Aisbett, Biometrics 47, 461-66, 1991
# Data on the recurrence times to infection, at the point of insertion of
#  the catheter, for kidney patients using portable dialysis equipment.
#  Catheters may be removed for reasons other than infection, in which case
#  the observation is censored.  Each patient has exactly 2 observations.

# Variables: patient, time, status, age, 
#	   sex (1=male, 2=female),
#	   disease type (0=GN, 1=AN, 2=PKD, 3=Other)
#	   author's estimate of the frailty

# I don't match their answers, and I think that I'm right

kidney <- read.table('data.kidney', col.names=c("id", "time", "status",
				      "age", "sex", "disease", "frail"))
kidney$disease <- factor(kidney$disease, levels=c(3, 0:2),
			 labels=c('Other', 'GN', 'AN', "PKD"))

kfit <- coxph(Surv(time, status)~ age + sex + disease + frailty(id), kidney)
kfit1<- coxph(Surv(time, status) ~age + sex + disease +
	      frailty(id, theta=1), kidney, iter=20)
kfit0 <- coxph(Surv(time, status)~ age + sex + disease, kidney)
temp <-  coxph(Surv(time, status) ~age + sex + disease +
	      frailty(id, theta=1, sparse=F), kidney)


# Check out the EM based score equations
#  temp1 and kfit1 should have essentially the same coefficients
#  temp2 should equal kfit1$frail
# equality won't be exact because of the different iteration paths
temp1 <- coxph(Surv(time, status) ~ age + sex + disease +
	       offset(kfit1$frail[id]), kidney)
rr <- tapply(resid(temp1), kidney$id, sum)
temp2 <- log(rr/1 +1)
all.equal(temp1$coef, kfit1$coef) 
all.equal(temp2, kfit1$frail)



kfit
kfit1
kfit0
temp

#
# Now fit the data using REML
#
kfitm1 <- coxph(Surv(time,status) ~ age + sex + disease + 
		frailty(id, dist='gauss'), kidney)
kfitm2 <- coxph(Surv(time,status) ~ age + sex + disease + 
		      frailty(id, dist='gauss', sparse=F), kidney)
kfitm1
summary(kfitm2)
