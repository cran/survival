#
# Do a real live jackknife on one of the Stanford models
#
xxx <- jasa1$x[, c(1,3,4,5,7)]
ag.diff <- matrix(double(5*103), ncol=5)
for (i in 1:103) {
    who <- !(jasa1$id==i)
    tfit <- agreg(start[who], stopp[who], event[who], xxx[who,])$coef
    ag.diff[i,] <- sfit.1$coef - tfit
    }
# agdiff now contains the change in beta due to adding that obs to the
#   data; a reasonable measure of influence.
