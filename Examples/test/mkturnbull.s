#
# The test data set from Turnbull, JASA 1974, 169-73.
#
#  status  0=right censored
#          1=exact
#          2=left censored
#

turnbull <- data.frame( time  =c( 1,1,1, 2,2,2, 3,3,3, 4,4,4),
			status=c( 1,0,2, 1,0,2, 1,0,2, 1,0,2),
			  n   =c(12,3,2, 6,2,4, 2,0,2, 3,3,5))
