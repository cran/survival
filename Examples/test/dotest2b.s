#
# Test out the survival curve and variance, in gory detail
#
xdata <- data.frame(rbind(test1,test1), ss = rep(1:2, rep(nrow(test1),2)))
fit <- coxph(Surv(time, status)~x + strata(ss), xdata, method='breslow')
sfit <- survfit(fit, list(x=0))	 #type='aalen' is default

# From the hand worked notes
bb <- as.vector(exp(fit$coef))
realhaz <- cumsum(c(1/(3*bb+3), 2/(bb+3), 1))
all.equal(sfit$surv, exp(-rep(realhaz,2)))

realvar <- cumsum(c( (1/(3*bb+3))^2, 2/(bb+3)^2, 1))
dd <- cumsum(c( (bb/(bb+1))*(1/(3*bb+3)),
	        (bb/(bb+3))*(2/(bb+3)),  0*1))
realvar <- realvar + dd^2 * fit$var
all.equal(sfit$std, sqrt(rep(realvar,2)))
summary(sfit)

# Get the Kalbfleisch-Prentice survival with Greenwood variance
#  
sfit <- survfit(fit, list(x=0), type='kalb')	 

# second term is the solution to an equation, in the neighborhood of .6
tfun <- function(alpha) (bb/(1-alpha^bb) + 1/(1-alpha) - (bb+3))^2
temp <- nlminb(0.6, tfun, lower=.1, upper=.9)

realkm <- cumprod(c((1- bb/(3*bb+3))^(1/bb), temp$par, 0))
all.equal(sfit$surv, rep(realkm,2))

realvar <- cumsum(c( 1/((3*bb+3)*(2*bb+3)), 2/((bb+3)*2), 1))
dd <- cumsum(c( (bb/(bb+1))*(1/(3*bb+3)),
	        (bb/(bb+3))*(2/(bb+3)),  0*1))
realvar <- realvar + dd^2 * fit$var
all.equal(sfit$std, sqrt(rep(realvar,2)))
summary(sfit)


#
# Repeat with the Efron approximation for Cox, Efron estimates
#
fit <- coxph(Surv(time, status)~x + strata(ss), xdata)
sfit <- survfit(fit, list(x=0))	 
bb <- as.vector(exp(fit$coef))

realhaz <- cumsum(c(1/(3*bb+3), 1/(bb+3) + 2/(bb+5), 1))
all.equal(sfit$surv, exp(-rep(realhaz,2)))

realvar <- cumsum(c( (1/(3*bb+3))^2, 1/(bb+3)^2 + 4/(bb+5)^2, 1))
dd <- cumsum(c( (bb/(bb+1))*(1/(3*bb+3)),
	        (bb/(bb+3))*(1/(bb+3))+ (bb/(bb+5))*(2/(bb+5)),  0*1))
realvar <- realvar + dd^2 * fit$var
all.equal(sfit$std, sqrt(rep(realvar,2)))
summary(sfit)


rm(bb, sfit, realhaz, realvar, fit, realkm, xdata)
