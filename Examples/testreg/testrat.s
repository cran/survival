#
# Some tests using the rat data
#
rats <- read.table('../testfrail/data.rats', 
		   col.names=c('litter', 'rx', 'time', 'status'))

rfitnull <- survreg(Surv(time, status) ~1, rats)
temp <- rfitnull$scale^2 * pi^2/6
cat("Effective n =", round(temp*(solve(rfitnull$var))[1,1],1), "\n")

rfit0 <- survreg(Surv(time, status) ~ rx , rats)
print(rfit0)

rfit1 <- survreg(Surv(time, status) ~ rx + factor(litter), rats)
temp <- rbind(c(rfit0$coef, rfit0$scale), c(rfit1$coef[1:2], rfit1$scale))
dimnames(temp) <- list(c("rfit0", "rfit1"), c("Intercept", "rx", "scale"))
temp


rfit2a <- survreg(Surv(time, status) ~ rx +
		  frailty.gaussian(litter, df=13, sparse=F), rats )
rfit2b <- survreg(Surv(time, status) ~ rx +
		  frailty.gaussian(litter, df=13, sparse=T), rats )

rfit3a <- coxph(Surv(time,status) ~ rx + 
		  frailty.gaussian(litter, df=13, sparse=F), rats )
rfit3b <- coxph(Surv(time,status) ~ rx + 
		frailty(litter, df=13, dist='gauss'), rats)

temp <- cbind(rfit2a$coef[3:52], rfit2b$frail, rfit3a$coef[2:51], rfit3b$frail)
dimnames(temp) <- list(NULL, c("surv","surv.sparse","cox","cox.sparse"))
pairs(temp)
apply(temp,2,var)/c(rfit2a$scale, rfit2b$scale, 1,1)^2
apply(temp,2,mean)

# The parametric model gives the coefficients less variance for the
#  two fits, for the same df, but the scaled results are similar.
# 13 df is near to the rmle for the rats

rm(temp, rfit2a, rfit2b, rfit3a, rfit3b, rfitnull, rfit0, rfit1)

