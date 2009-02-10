options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Check out using a "user specified" distribution
#
mydist <- c(survreg.distributions$extreme, survreg.distributions$weibull[-1])
mydist$name <- "Weibull2"
mydist$dist <- NULL

fit1 <- survreg(Surv(time, status) ~ age + ph.ecog, lung)
fit2 <- survreg(Surv(time, status) ~ age + ph.ecog, lung, dist=mydist)

all.equal(fit1$coef, fit2$coef)
all.equal(fit1$var, fit2$var)
