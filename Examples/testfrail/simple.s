#
# Test the logic of the new program, by fitting some no-frailty models
#  (theta=0).  It should give exactly the same answers as 'ordinary' coxph.
# By default frailty models run with eps=1e-7, ordinary with 1e-4.  I match
#   these to get the same number of iterations.
#
test1 <- data.frame(time=  c(4, 3,1,1,2,2,3),
                    status=c(1,NA,1,0,1,1,0),
                    x=     c(0, 2,1,1,1,0,0))

test2 <- data.frame(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                    stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                    event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                    x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )

zz <- rep(0, nrow(test1))
tfit1 <- coxph(Surv(time,status) ~x, test1, eps=1e-7)
tfit2 <- coxph(Surv(time,status) ~x + frailty(zz, theta=0, sparse=T), test1)
tfit3 <- coxph(Surv(zz,time,status) ~x + frailty(zz, theta=0,sparse=T), test1)

temp <- c('coefficients', 'var', 'loglik', 'linear.predictors',
	  'means', 'n')

all.equal(tfit1[temp], tfit2[temp])
all.equal(tfit1[temp], tfit3[temp])

zz <- rep(0, nrow(test2))
tfit1 <- coxph(Surv(start, stop, event) ~x, test2, eps=1e-7)
tfit2 <- coxph(Surv(start, stop, event) ~ x + frailty(zz, theta=0, sparse=T), 
	       test2)
all.equal(tfit1[temp], tfit2[temp])



