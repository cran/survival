#
# Test out subscripting in the case of a coxph survival curve
#
fit <- coxph(Surv(time, status) ~ age + sex + meal.cal + strata(ph.ecog),
		data=cancer)

surv1 <- survfit(fit)
temp <- surv1[2:3]

which <- cumsum(surv1$strata)
zed   <- (which[1]+1):(which[3])
all.equal(surv1$surv[zed], temp$surv)
all.equal(surv1$time[zed], temp$time)

#
# Now a result with a matrix of survival curves
#
dummy <- data.frame(age=c(30,40,60), sex=c(1,2,2), meal.cal=c(500, 1000, 1500))
surv2 <- survfit(fit, newdata=dummy)

zed <- 1:which[1]
all.equal(surv2$surv[zed,1], surv2[1,1]$surv)
all.equal(surv2$surv[zed,2], surv2[1,2]$surv)
all.equal(surv2$surv[zed,3], surv2[1,3]$surv)
all.equal(surv2$surv[zed, ], surv2[1,1:3]$surv)
all.equal(surv2$surv[zed],   (surv2[1]$surv)[,1])
all.equal(surv2$surv[zed, ], surv2[1, ]$surv)
