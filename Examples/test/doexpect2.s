#  Tests of expected survival

#
# Simple case 2: make sure that the averages are correct, for Ederer method
#
#  Compute the 1, 5, 10 and 12 year expected survival

temp1 <- mdy.date(1:6,6:11,36:41)
temp2 <- mdy.date(6:1,11:6,55:50)
age <- temp2 - temp1

exp1 <- survexp(~ratetable(year=temp2, age=(temp2-temp1), sex=1),
		       times=c(366, 1827, 3653, 4383))
exp2 <- survexp(~ratetable(year=temp2, age=(temp2-temp1), sex=1) + I(1:6),
			times=c(366, 1827, 3653, 4383))

print(all.equal(exp1$surv, apply(exp2$surv, 1, mean)))

#
# Check that adding more time points doesn't change things
#
exp3 <- survexp(~ratetable(year=temp2, age=(temp2-temp1), sex=1),
		times=sort(c(366, 1827, 3653, 4383, 30*(1:100))))

exp1$surv
exp3$surv[match(exp1$time, exp3$time, nomatch=0)]


#
# Now test Hakulinen's method, assuming an analysis date of 3/1/57
#
futime <- mdy.date(3,1,57) - temp2
xtime  <- sort(c(futime, 30, 60, 185, 365))

exp1 <- survexp(futime ~ ratetable(year=temp2, age=(temp2-temp1), sex=1),
		times=xtime, conditional=F)
exp2 <- survexp(~ratetable(year=temp2, age=(temp2-temp1), sex=1) + I(1:6),
			times=futime)

wt <- rep(1,6)
con <- double(6)
for (i in 1:6) {
    con[i] <- sum(exp2$surv[i,i:6])/sum(wt[i:6])
    wt <- exp2$surv[i,]
    }

exp1$surv[match(futime, xtime)]
cumprod(con)                  #should be equal


#
# Now for the conditional method
#
exp1 <- survexp(futime ~ ratetable(year=temp2, age=(temp2-temp1), sex=1),
		times=xtime, conditional=T)

cond <- exp2$surv
for (i in 6:2) cond[i,] _ (cond[i,]/cond[i-1,])  #conditional survival
for (i in 1:6) con[i] _ exp(mean(log(cond[i, i:6])))

all.equal(exp1$surv[match(futime, xtime)], cumprod(con))
cumprod(con)

rm(con, cond, exp1, exp2, wt, temp1, temp2, age)
rm(exp3, futime, xtime)
