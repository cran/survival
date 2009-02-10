options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# The Stanford data from 1980 is used in Escobar and Meeker
#	t5 = T5 mismatch score
#  Their case numbers correspond to a data set sorted by age
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

stanford2$t5 <- ifelse(stanford2$t5 <0, NA, stanford2$t5)
stanford2 <- stanford2[order(stanford2$age, stanford2$time),]
stanford2$time <- ifelse(stanford2$time==0, .5, stanford2$time)

cage <- stanford2$age - mean(stanford2$age)
fit1 <- survreg(Surv(time, status) ~ cage + I(cage^2), stanford2,
		dist='lognormal')
fit1
ldcase <- resid(fit1, type='ldcase')
ldresp <- resid(fit1, type='ldresp')
print(ldresp)
# The ldcase and ldresp should be compared to table 1 in Escobar and 
#  Meeker, Biometrics 1992, p519; the colum they label as (1/2) A_{ii}

plot1 <- function() {
    # make their figure 1, 2, and 6
    plot(stanford2$age, stanford2$time, log='y', xlab="Age", ylab="Days",
	 ylim=c(.01, 10^6))
    temp <- predict(fit1, type='response', se.fit=T) 
    matlines(stanford2$age, cbind(temp$fit, temp$fit-1.96*temp$se.fit,
				            temp$fit+1.96*temp$se.fit),
	     lty=c(1,2,2))
    # these are the wrong CI lines, he plotted std dev, I plotted std err
    # here are the right ones
    #  Using uncentered age gives different coefs, but makes prediction over an
    #    extended range somewhat simpler 
    refit <- survreg(Surv(time,status)~ age + age^2, stanford2,
		     dist='lognormal')
    plot(stanford2$age, stanford2$time, log='y', xlab="Age", ylab="Days",
	 ylim=c(.01, 10^6), xlim=c(0,75))
    temp2 <- predict(refit, list(age=1:75), type='quantile', p=c(.05, .5, .95))
    matlines(1:75, temp2, lty=c(1,2,2), col=2)

    n <- length(ldcase)
    plot(1:n, ldcase, xlab="Case Number", ylab="(1/2) A")
    title (main="Case weight pertubations")
    plot(1:n, ldresp, xlab="Case Number", ylab="(1/2) A")
    title(main="Response pertubations")
    }

postscript('meekerplot.ps')
plot1()
dev.off()
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
# aeq(p3$fit + attr(p3$fit, 'constant'), p4$fit)  #R is missing the attribute
aeq(p4$fit, p5$fit)
aeq(p3$se.fit, p4$se.fit)  #this one should be false
aeq(p4$se.fit, p5$se.fit)  #this one true

