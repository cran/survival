#
# An example of how to manage the residuals
#
attach(jasa1)
fit <- coxph(Surv(start, stop, event) ~ (age + surgery)* transplant)
#
# The true residual for a subject is the sum of the rows that were used to
#   represent him/her to agreg.
mresid <- resid(fit, collapse=jasa1$id)
mresid

#
# The score resid has multiple columns
sresid <- resid(fit, type='score', collapse=id)
# Standardized score resids: percent change in coef if this obs were dropped
-100 * sresid %*% fit$var %*% diag(1/fit$coef)

#
# Do a real live jackknife on one of the Stanford models
#   The "acid test" of the score residuals
#
ag.diff <- matrix(double(5*103), ncol=5)
for (i in 1:103) {
    who <- !(id==i)
    tfit <- coxph(Surv(start, stop, event) ~ (age + surgery)*transplant,
			  subset=who)$coef
    ag.diff[i,] <- fit$coef - tfit
    }
# agdiff now contains the change in beta due to adding that obs to the
#   data; a reasonable measure of influence.

sresid <- sresid %*% fit$var
# Best check-  commented out so I can script this
#  for (i in 1:5) {
#     plot(ag.diff[,i], sresid[,i])
#     abline(0,1)
#     }
temp <- sqrt(apply(ag.diff, 2, var))
temp <- (ag.diff - sresid) %*% diag(1/temp)   #scaled difference
apply(abs(temp), 2, mean)   #it's about 4% error

detach("jasa1")
