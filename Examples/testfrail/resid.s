#
# The residual methods treat a sparse frailty as a fixed offset with
#   no variance
#

kfit1 <- coxph(Surv(time, status) ~ age + sex + 
	           frailty(id, dist='gauss'), kidney)
tempf <- predict(kfit1, type='terms')[,3]
temp  <- kfit1$frail[match(kidney$id, sort(unique(kidney$id)))]
#all.equal(unclass(tempf), unclass(temp))
all.equal(as.vector(tempf), as.vector(temp))

# Now fit a model with explicit offset
kfitx <- coxph(Surv(time, status) ~ age + sex + offset(tempf),kidney,
	       eps=1e-7)

# These are not precisely the same, due to different iteration paths
all.equal(kfitx$coef, kfit1$coef)

# This will make them identical
kfitx <- coxph(Surv(time, status) ~ age + sex  + offset(temp),kidney,
	       iter=0, init=kfit1$coef)
all.equal(resid(kfit1), resid(kfitx))
all.equal(resid(kfit1, type='score'), resid(kfitx, type='score'))
all.equal(resid(kfit1, type='schoe'), resid(kfitx, type='schoe'))

# These are not the same, due to a different variance matrix
#  The frailty model's variance is about 2x the naive "assume an offset" var
# The score residuals are equal, however.
all.equal(resid(kfit1, type='dfbeta'), resid(kfitx, type='dfbeta'))
zed <- kfitx
zed$var <- kfit1$var
all.equal(resid(kfit1, type='dfbeta'), resid(zed, type='dfbeta'))


temp1 <- resid(kfit1, type='score')
temp2 <- resid(kfitx, type='score')
all.equal(temp1, temp2)

#
# Now for some tests of predicted values
#
all.equal(predict(kfit1, type='expected'), predict(kfitx, type='expected'))
all.equal(predict(kfit1, type='lp'), predict(kfitx, type='lp'))

temp1 <- predict(kfit1, type='terms', se.fit=T)
temp2 <- predict(kfitx, type='terms', se.fit=T)
all.equal(temp1$fit[,1:2], temp2$fit)
all.equal(temp1$se.fit[,1:2], temp2$se.fit)  #should be false
mean(temp1$se.fit[,1:2]/ temp2$se.fit)
all.equal(as.vector(temp1$se.fit[,3])^2, 
	  as.vector(kfit1$fvar[match(kidney$id, sort(unique(kidney$id)))]))

print(temp1)
kfit1
kfitx

rm(temp1, temp2, kfitx, zed, tempf)
