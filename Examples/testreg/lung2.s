#
# Test out ridge regression and splines
#
lfit0 <- survreg(Surv(time, status) ~1, lung)
lfit1 <- survreg(Surv(time, status) ~ age + ridge(ph.ecog, theta=5), lung)
lfit2 <- survreg(Surv(time, status) ~ sex + ridge(age, ph.ecog, theta=1), lung)
lfit3 <- survreg(Surv(time, status) ~ sex + age + ph.ecog, lung)

lfit0
lfit1
lfit2
lfit3


xx <- pspline(lung$age, nterm=3, theta=.3)
xx <- matrix(unclass(xx), ncol=ncol(xx))   # the raw matrix
lfit4 <- survreg(Surv(time, status) ~xx, lung)
lfit5 <- survreg(Surv(time, status) ~age, lung)

lfit6 <- survreg(Surv(time, status)~pspline(age, df=2), lung)
plot(lung$age, predict(lfit6), xlab='Age', ylab="Spline prediction")
title("Lung Data")
     
lfit7 <- survreg(Surv(time, status) ~ offset(lfit6$lin), lung)

lfit4
lfit5
lfit6
lfit7$coef

rm(lfit1, lfit2, lfit3, lfit4, lfit5, lfit6, lfit7)
rm(xx, lfit0)
