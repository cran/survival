#
# Fit the models found in Wei et. al.
#
options(contrasts='contr.treatment')
wfit <- coxph(Surv(stop, event) ~ (rx + size + number)* strata(enum) +
		 cluster(id), bladder, method='breslow')
wfit

# Check the rx coefs versus Wei, et al, JASA 1989
rx <- c(1,4,5,6)  # the treatment coefs above
cmat <- diag(4); cmat[1,] <- 1;          #contrast matrix
wfit$coef[rx] %*% cmat                   # the coefs in their paper (table 5)
t(cmat) %*% wfit$var[rx,rx] %*% cmat  # var matrix (eqn 3.2)

# Anderson-Gill fit
fita <- coxph(Surv(start, stop, event) ~ rx + size + number + cluster(id),
		  bladder2,  method='breslow')
fita

# Prentice fits.  Their model 1 a and b are the same
fit1p  <- coxph(Surv(stop, event) ~ rx + size + number, bladder2,
		subset=(enum==1), method='breslow')
fit2pa <- coxph(Surv(stop, event) ~ rx + size + number, bladder2,
		subset=(enum==2), method='breslow')
fit2pb <- coxph(Surv(stop-start,  event) ~ rx + size + number, bladder2,
		   subset=(enum==2), method='breslow')
fit3pa <- coxph(Surv(stop, event) ~ rx + size + number, bladder2,
		subset=(enum==3), method='breslow')
 #and etc.
fit1p
fit2pa
fit2pb
fit3pa
