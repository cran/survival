#
# The special case of a single sparse frailty
#

kfit1 <- coxph(Surv(time, status) ~ frailty(id, dist='gauss'), kidney)
tempf <- predict(kfit1, type='terms')
temp  <- kfit1$frail[match(kidney$id, sort(unique(kidney$id)))]
all.equal(as.vector(tempf), as.vector(temp))

# Now fit a model with explicit offset
kfitx <- coxph(Surv(time, status) ~ offset(tempf),kidney, eps=1e-7)

all.equal(resid(kfit1), resid(kfitx))
all.equal(resid(kfit1, type='deviance'), resid(kfitx, type='deviance'))

#
# Some tests of predicted values
#
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
aeq(predict(kfit1, type='expected'), predict(kfitx, type='expected'))
aeq(predict(kfit1, type='lp'), predict(kfitx, type='lp'))

temp1 <- predict(kfit1, type='terms', se.fit=T)
all.equal(temp1$fit, kfitx$linear)
all.equal(temp1$se.fit^2, 
	  kfit1$fvar[match(kidney$id, sort(unique(kidney$id)))])

temp1
kfit1


