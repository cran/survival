# Create a "counting process" version of the simplest test data set
#
test1b<- list(start= c(0, 3,  0,  0, 5,  0, 6,14,  0,  0, 10,20,30, 0),
	      stop = c(3,10, 10,  5,20,  6,14,20, 30,  10,20,30,40, 10),
	      status=c(0, 1,  0,  0, 1,  0, 0, 1,  0,   0, 0, 0, 1,  0),
	      x=     c(1, 1,  1,  1, 1,  0, 0, 0,  0,   0, 0, 0, 0,  NA),
	      id =   c(3, 3,  4,  5, 5,  6, 6, 6,  7,   1, 1, 1, 1,   2))
