# Do a Stanford heart transplant data the way that K&P do
#
# Input - a data frame containing the raw data
# Output- a data frame with multiple obs for the transplanted subjects
#
#  There are more efficient ways to do this than the "for" loop below, but
#    this script is much more readable than any of them.
#
stan1 <- function(jasa) {
    id <- row.names(jasa)
    tx _  1*(jasa$transplant)
    covar _ cbind(jasa$age/365.25 -48,
		  (jasa$accept.dt - mdy.date(10,1,1967))/365.25,
		  jasa$surgery)
    n _ length(tx)     #number in study
    ntx _ sum(tx)      #number transplanted
    # Per paragraph 2, p138, the patient who died on the day of transplant is
    #   treated differently, for all others deaths are assumed to occur "earlier
    #   in the day" than transplants.
    special <- id ==38
    wait _ jasa$wait.time
    wait[special] _ wait[special] - .5

    age <- year <- surgery <- transplant <- id2 <- double(n+ntx)
    start <- stop <- event <- double(n+ntx)
    ii <- 1
    for (i in 1:n) {
	age[ii] <- covar[i,1]
	year[ii] <- covar[i,2]
	surgery[ii] <- covar[i,3]
	transplant[ii] <- 0
	id2[ii] <- id[i]

	if (tx[i])  { #transplanted  - 2 lines if data
	    start[ii] <- 0
	    stop[ii]  <- wait[i]
	    event[ii] <- 0

	    ii <- ii+1
	    start[ii] <- wait[i]
	    stop[ii]  <- jasa$futime[i]
	    event[ii] <- jasa$fustat[i]
	    age[ii] <- covar[i,1]
	    year[ii] <- covar[i,2]
	    surgery[ii] <- covar[i,3]
	    transplant[ii] <- 1
	    id2[ii] <- id[i]
	    }
	else {                             # one line of data
	    start[ii] <-0
	    stop[ii] <- jasa$futime[i]
	    event[ii]<- jasa$fustat[i]
	    }
	ii <- ii+1
	}

    data.frame(start, stop, event, age, year, surgery,
		    transplant=factor(transplant), id=id2)
    }

