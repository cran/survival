#
# The AML data, from Miller, "Survival Analysis", page 49.
#
aml <- list(time=   c( 9, 13, 13, 18, 23, 28, 31, 34, 45, 48, 161,
			   5,  5,  8,  8, 12, 16, 23, 27, 30, 33, 43, 45),
	    status= c( 1,1,0,1,1,0,1,1,0,1,0, 1,1,1,1,1,0,1,1,1,1,1,1),
	    x     = as.factor(c(rep("Maintained", 11),
				rep("Nonmaintained", 12) )))

aml <- data.frame(aml)
