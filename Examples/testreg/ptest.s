#
# Look at predicted values
#
ofit1 <- survreg(Surv(futime, fustat) ~ age + ridge(ecog.ps, rx), ovarian)

predict(ofit1)
predict(ofit1, type='response')
predict(ofit1, type='terms', se=T)

temp1 <- predict(ofit1, se=T)
temp2 <- predict(ofit1, type= 'response', se=T)
all.equal(temp2$se.fit, temp1$se.fit* sqrt(exp(temp1$fit)))
