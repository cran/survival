#
# Stanford predictions in other ways
#
fit2 <- survreg(Surv(time, status) ~ poly(age,2), stanford2,
		dist='lognormal')

p1 <- predict(fit1, type='response')
p2 <- predict(fit2, type='response')
aeq(p1, p2)

p3 <- predict(fit2, type='terms', se=T)
p4 <- predict(fit2, type='lp', se=T)
p5 <- predict(fit1, type='lp', se=T)
aeq(p3$fit + attr(p3$fit, 'constant'), p4$fit)
aeq(p4$fit, p5$fit)
aeq(p3$se.fit, p4$se.fit)  #this one should be false
aeq(p4$se.fit, p5$se.fit)  #this one true

