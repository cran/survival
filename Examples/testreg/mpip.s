temp <- matrix(scan("data.mpip", skip=23), ncol=13, byrow=T)
dimnames(temp) <- list(NULL, c('ved', 'angina', 'education', 'prior.mi',
                     'nyha', 'rales', 'ef', 'ecg', 'angina2', 'futime', 
                     'status', 'admit', 'betab'))
 
mpip <- data.frame(temp)
lved <- log(mpip$ved + .02)

fit1 <- coxph(Surv(futime, status) ~ pspline(lved) + factor(nyha) + 
	      rales + pspline(ef), mpip)

temp <- predict(fit1, type='terms', se.fit=T)
yy <- cbind(temp$fit[,4], temp$fit[,4] + 1.96*temp$se[,4],
	                  temp$fit[,4] - 1.96*temp$se[,4])
index <- order(mpip$ef)
matplot(mpip$ef[index], yy[index,], type='l', lty=c(1,2,2), col=1)
title(xlab="Ejection Fraction", ylab="Cox model risk", 
      main="Post-Infarction Survival")

fit2 <- coxph(Surv(futime, status) ~ lved + factor(nyha) + rales +
	      pspline(ef, df=0), mpip)
temp <- predict(fit2, type='terms', se.fit=T)
yy <- cbind(temp$fit[,4], temp$fit[,4] + 1.96*temp$se[,4],
	                  temp$fit[,4] - 1.96*temp$se[,4])
matplot(mpip$ef[index], yy[index,], type='l', lty=c(1,2,2), col=1)
title(xlab="Ejection Fraction", ylab="Cox model risk", 
      main="Post-Infarction Survival, AIC")


fit3 <- survreg(Surv(futime, status) ~ lved + factor(nyha) + rales +
		pspline(ef, df=2), mpip, dist='lognormal')
temp <- predict(fit3, type='terms', se.fit=T)
yy <- cbind(temp$fit[,4], temp$fit[,4] + 1.96*temp$se[,4],
	                  temp$fit[,4] - 1.96*temp$se[,4])
matplot(mpip$ef[index], yy[index,], type='l', lty=c(1,2,2), col=1)
title(xlab="Ejection Fraction", ylab="Log-normal model predictor", 
      main="Post-Infarction Survival")
