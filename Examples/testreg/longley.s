#
# Do a classic ridge regression on the longley data
#
longley <- read.table('data.longley',
		      col.names=c('employed', 'deflator', 'gnp', 'unemployed',
			'military', 'pop', 'year'))
lfit1 <- lm(employed ~ deflator + gnp + unemployed + military + pop + year,
	    longley)
lfit2 <- survreg(Surv(employed) ~ ., longley, dist='gaussian', 
		 toler.chol=1e-15)

all.equal(lfit1$coef, lfit2$coef, toler=1e-5)

longley2 <- longley
for (i in 2:7) longley2[[i]] <- longley2[[i]] - mean(longley2[[i]])
lfit3 <- lm(employed~., longley2)
lfit4 <- survreg(Surv(employed)~., longley2, dist='gaussian',
		 toler.chol=1e-13)

lfun <- function(theta, data) {
    assign('tempxxx', theta, frame=0)
    survreg(Surv(employed) ~ ridge(deflator, gnp, unemployed, military,
				   pop, year, theta=tempxxx), data=data,
	    dist='gaussian', toler.chol=1e-15)
    }


longley3 <- longley2
longley3$employed <- longley2$employed + sample(c(-2,2), 15, repl=T)


theta<- seq(-14, 4, length=20)
cmat1 <- matrix(0, 20, length(lfit1$coef))
cmat2 <- matrix(0, 20, length(lfit1$coef))
clog <- double(20)
for (i in 1:20) {
    temp  <- lfun(10^(theta[i]), longley2)
    clog[i] <- temp$loglik[2]
    cmat1[i,] <- temp$coef
    temp  <- lfun(10^(theta[i]), longley3)
    cmat2[i,] <- temp$coef
    }

matplot(10^theta, (cmat1[,-1]- cmat2[,-1]) %*% diag(1/lfit3$coef[-1]), 
	log='x', type='l')


