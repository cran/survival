#
# Create a "user defined" rate table, using the smoking data
#
temp <- scan("data.smoke")/100000
temp <-  matrix(temp, ncol=8, byrow=T)
smoke.rate <- c(rep(temp[,1],6), rep(temp[,2],6), temp[,3:8])
attributes(smoke.rate) <- list(
	dim=c(7,2,2,6,3),
	dimnames=list(c("45-49","50-54","55-59","60-64","65-69","70-74","75-79"),
		      c("1-20", "21+"),
		      c("Male","Female"),
		      c("<1", "1-2", "3-5", "6-10", "11-15", ">=16"),
		      c("Never", "Current", "Former")),
	dimid=c("age", "amount", "sex", "duration", "status"),
	factor=c(0,1,1,0,1),
	cutpoints=list(c(45,50,55,60,65,70,75),NULL, NULL,
				     c(0,1,3,6,11,16),NULL),
	class='ratetable'
	)
rm(temp)

is.ratetable(smoke.rate)
summary(smoke.rate)
print(smoke.rate)

summary(smoke.rate[1:3,,1,,])  #test subscripting
