# 
# Simple case: a single male subject, born 6/6/36 and entered on study 6/6/55.
#

temp1 <- mdy.date(6,6,36)
temp2 <- mdy.date(6,6,55)# Now compare the results from person-years
#
temp.age <- tcut(temp2-temp1, floor(c(-1, (18:31 * 365.24))),
	labels=c('0-18', paste(18:30, 19:31, sep='-')))
temp.yr  <- tcut(temp2, mdy.date(1,1,1954:1965), labels=1954:1964)
temp.time <- 3700   #total days of fu
py1 <- pyears(temp.time ~ temp.age + temp.yr, scale=1) #output in days

# The subject should appear in 20 cells
#    6/6/55 - 12/31/55, 209 days, age 19-20, 1955
#    1/1/56 -  6/ 4/56, 156 days, age 19-20, 1956
#    6/5/56 - 12/31/56, 210 days, age 20-21, 1956   (a leap year, and his
#			birthday computes one day earlier)
#    1/1/57 -  6/ 5/57, 156 days, age 20-21, 1957
#    6/6/57 - 12/31/57, 209 days, age 21-22, 1957
# and etc
#   with 203 days "off table", ie, beyond the last cell of the table
#
# It is a nuisance, but tcut follows 'cut' in that we give the ENDS of
#  the intervals, whereas the survival tables use the starts of intervals.
#  Thus this breakdown does not match that in doexpect.s
#
xx <- matrix(0, nrow=14, ncol=11)
xx[cbind(3:11, 3:11)] <- 156
xx[cbind(3:12, 2:11)] <- c(209, 210, rep(c(209, 209, 209, 210),2))
dimnames(xx) <- list(c('0-18', paste(18:30, 19:31, sep='-')), 1954:1964)
all.equal(xx, py1$pyears)
all.equal(203, py1$offtable)
all.equal(1*(xx>0), py1$n)

#
# Now with expecteds
#
py2 <- pyears(temp.time ~ temp.age + temp.yr
		+ ratetable(age=temp2-temp1, year=temp2, sex=1),
	     scale=1, ratetable=survexp.uswhite ) #output in days
all.equal(xx, py2$pyears)
all.equal(203, py2$offtable)
all.equal(1*(xx>0), py2$n)


py3 <-  pyears(temp.time ~ temp.age + temp.yr
		+ ratetable(age=temp2-temp1, year=temp2, sex=1),
	     scale=1, ratetable=survexp.uswhite , expect='pyears')
all.equal(py2$n, py3$n)
all.equal(py2$pyear, py3$pyear)
all.equal(py3$n, 1*(py3$expect>0))

# Now, compute the py3 result "by hand".  Since there is only one person
#   it can be derived from py2.
#
xx1 <- py2$expect[py2$n>0]   		# the hazard over each interval
cumhaz <- cumsum(c(0, xx1[-length(xx1)]))     # the cumulative hazard	
xx2 <- py3$expect[py3$n>0]   		# the expected number of person days
xx3 <- py3$pyears[py3$n>0]   		# the potential number of person days

# This is the integral of the curve "exp(-haz *t)" over the interval
integral <- xx3 * exp(-cumhaz)* (1- exp(-xx1))/ xx1
# They might not be exactly equal, since the C code tracks changes in the
#   rate tables that occur -within- an interval.  So try for 7 digits
all.equal(round(integral,4), round(xx2,4))

# Cut off the bottom of the table, instead of the side
temp.age <- tcut(temp2-temp1, floor(c(-1, (18:27 * 365.24))),
	labels=c('0-18', paste(18:26, 19:27, sep='-')))

py4 <- eval(py3$call)
all.equal(py4$pyear, py3$pyear[1:10,])
all.equal(py4$expect, py3$expect[1:10,])

rm(xx1, xx2, xx3, cumhaz, temp1, temp2)
rm(temp.age, temp.yr, temp.time, xx)
