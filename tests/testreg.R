options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Run a test that can be verified using other packages (we used SAS)
#
test1 <- data.frame(time=  c(9, 3,1,1,6,6,8),
                    status=c(1,NA,1,0,1,1,0),
                    x=     c(0, 2,1,1,1,0,0))
fit1w <- survreg(Surv(time, status) ~x, test1, dist='weibull')
fit1w
summary(fit1w)

fit1e <- survreg(Surv(time, status) ~x, test1, dist='exponential')
fit1e
summary(fit1e)

fit1l <- survreg(Surv(time, status) ~x, test1, dist='loglogistic')
fit1l
summary(fit1l)

fit1g <- survreg(Surv(time, status) ~x, test1, dist='lognormal')
summary(fit1g)
#
#  Do a test with the ovarian data
#
fitfw <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian,
	dist='weibull')
fitfw

fitfl <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian,
	dist='loglogistic')
fitfl

flem2 <- scan("data.fleming2", what=list(ltime=0, rtime=0))
flsurv<- Surv(flem2$ltime, flem2$rtime, type='interval2')

fitfw2 <- survreg(flsurv ~ age + ecog.ps, ovarian, dist='weibull')
summary(fitfw2)

fitfl2 <- survreg(flsurv ~ age + ecog.ps, ovarian, dist='loglogistic')
summary(fitfl2)

fitfg2 <- survreg(flsurv ~ age + ecog.ps, ovarian, dist='lognormal')
summary(fitfg2)
#
# Check out the survreg density and probability functions
#

# Gaussian
x <- -10:10
p <- seq(.1, .95, length=25)
all.equal(dsurvreg(x, 1, 5, 'gaussian'), dnorm(x, 1, 5))
all.equal(psurvreg(x, 1, 5, 'gaussian'), pnorm(x, 1, 5))
all.equal(qsurvreg(p, 1, 5, 'gaussian'), qnorm(p, 1, 5))

# Lognormal
x <- 1:10
all.equal(dsurvreg(x, 1, 5, 'lognormal'), dlnorm(x, 1, 5))
all.equal(psurvreg(x, 1, 5, 'lognormal'), plnorm(x, 1, 5))
all.equal(qsurvreg(p, 1, 5, 'lognormal'), qlnorm(p, 1, 5))

# Weibull
lambda <- exp(-2)
rho    <- 1/3
temp <- (lambda*x)^rho
all.equal(psurvreg(x, 2, 3), 1- exp(-temp))
all.equal(dsurvreg(x, 2, 3), lambda*rho*(lambda*x)^(rho-1)*exp(-temp))
