#  Tests of expected survival

#
# Simple case 1: a single male subject, born 6/6/36 and entered on study 6/6/55.
#
#  Compute the 1, 5, 10 and 12 year expected survival

temp1 <- mdy.date(6,6,36)
temp2 <- mdy.date(6,6,55)
exp1 <- survexp(~ratetable(year=temp2, age=(temp2-temp1), sex=1),
		    ratetable=survexp.uswhite,times=c(366, 1827, 3653, 4383))

# Well, almost easy.  The uswhite table assumes that someone age 19 will have
#   seen 5 leap years-- but this lad has only seen 4.  Thus the routine sees
#   him as 1 day shy of his birthday.
#   (The age 19 category starts at 6940, but he is 6939 days old at entry).
# Thus his passage through the table is a bit more complicated

# First epoch:  1 day  at age 18, 1954    6/6/55
#             365 days at age 19, 1955    6/7/55 - 6/5/56

# Second      365 days at age 20, 1956    6/6/56 - 6/5/57
#             365 days at age 21, 1957    6/6/57 - 6/5/58
#             366 days at age 22, 1958    6/6/58 - 6/6/59
#             365 days at age 23, 1959    6/7/59 - 6/5/60

# Third       365 days at age 24, 1960    6/6/60 - 6/5/61
#             365 days at age 25, 1961    6/6/61 - 6/5/62
#             366 days at age 26, 1962    6/6/62 - 6/6/63
#             365 days at age 27, 1963    6/7/63 - 6/5/64
#             365 days at age 28, 1964    6/6/64 - 6/5/64

# Fourth      365 days at age 29, 1965    6/6/65 - 6/5/66
#             365 days at age 30, 1966    6/6/66 - 6/5/67

# remember, a first subscript of "1" is for age 0

xx <- survexp.uswhite[,1,]
check <- c(       (.6*xx[19,1] + .4*xx[19,2])    +
	     365* (.5*xx[20,1] + .5*xx[20,2]) ,

	     365* (.4*xx[21,1] + .6*xx[21,2]) +
	     365* (.3*xx[22,1] + .7*xx[22,2]) +
	     366* (.2*xx[23,1] + .8*xx[23,2]) +
	     365* (.1*xx[24,1] + .9*xx[24,2]) ,

	     365* (   xx[25,2]              ) +
	     365* (.9*xx[26,2] + .1*xx[26,3]) +
	     366* (.8*xx[27,2] + .2*xx[27,3]) +
	     365* (.7*xx[28,2] + .3*xx[28,3]) +
	     365* (.6*xx[29,2] + .4*xx[29,3]) ,

	     365* (.5*xx[30,2] + .5*xx[30,3]) +
	     365* (.4*xx[31,2] + .6*xx[31,3]) )

print(exp1$surv)
print(exp(-cumsum(check)))

# This does not pass the "all.equal" test.  Because of leap year (again),
#  the internal S code does not do exactly what is stated above.  US
#  rate tables are special: the entry for age 20 in 1950 is the probability
#  that someone who becomes 20 years old in 1950 will reach his 21st birthday.
#  In order to fit this into the general "cutpoints" calculations of
#  person-years, dates are adjusted so that everyone appears to have a
#  birthday on Jan 1.  But because of leap years, a birthday then moves to
#  Dec 31 sometimes --  this happens in each of the 366 day intervals above,
#  e.g., 1 day at age 22 with 1957 rates, 365 at age 22 with 1958 rates.
#  Upshot: they only agree to 6 decimals.  Close enough for me!

rm(temp1, temp2, exp1, check, xx)
