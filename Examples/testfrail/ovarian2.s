# Simulation for the ovarian data set
#
fit1 <- coxph(Surv(futime, fustat) ~ rx + ridge(age, ecog.ps, theta=1),
	      ovarian)

dfs <- eigen(solve(fit1$var, fit1$var2))$values

temp <- matrix(rnorm(30000), ncol=3)
temp2 <- apply((temp^2) %*% dfs, 1, sum)

round(rbind(quantile(temp2, c(.8, .9, .95, .99)), 
	     qchisq( c(.8, .9, .95, .99), sum(fit1$df))), 3)
