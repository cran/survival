# SCCS @(#)summary.ratetable.s	1.1 11/03/98
#
# Print out information about a rate table: it's dimensions and keywords
#
summary.ratetable <- function(object, ...) {
    rtable<-object
    if (!inherits(rtable, 'ratetable')) stop("Argument is not a rate table")

    att <- attributes(rtable)
    ncat <- length(dim(rtable))
    cat (" Rate table with", ncat, "dimensions:\n")
    for (i in 1:ncat) {
	if (att$factor[i]==0) {
	    cat("\t", att$dimid[i], " ranges from ", 
		format(min(att$cutpoints[[i]])), " to ", 
		format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
		" categories\n", sep='')
	    }
	else if(att$factor[i]==1) {
	     cat("\t", att$dimid[i], " has levels of: ",
		 paste(att$dimnames[[i]], collapse=' '), "\n", sep='')
	     }
	else {
	    cat("\t", att$dimid[i], " ranges from " , 
		format(min(att$cutpoints[[i]])), " to ", 
		format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
		" categories,\n\t\tlinearly interpolated in ",
		att$factor[i], " steps per division\n", sep='')
	    }
	}
    invisible(att)
    }

