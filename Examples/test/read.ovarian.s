#
# Read in the ovarian data
#

ovarian <- read.table("data.ovarian", row.names=NULL,
	    col.names= c("futime", "fustat", "age", "resid.ds", "rx", "ecog.ps"))
