These are the test examples for the S-PLUS version. The file Rtest in each
directory contains all the examples that run in R (on my machine) and
Rtest.out.save contains a copy of the output they produce with

	R --vanilla <Rtest >Rtest.out 
	diff Rtest.out.save Rtest.out

should produce nothing except for date changes in the startup message.

The testreg examples take quite a long time to run, since they include a
contout plot of the loglikelihood for a particularly nasty data set,
created by calling survreg() once for each grid point.

