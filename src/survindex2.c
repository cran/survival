/*     SCCS @(#)survindex2.c	5.2 10/27/98*/
/* A subroutine for surv.fit.print
**
** Input --
**      n:      number of survival times
**      stime:  survival times, must be >0
**      strata: strata values   data is sorted by time within strata
**      ntime:  number of time values
**      time:   time values for printout, must be >=0 and in increasing order
**      nstrat: number of strata.
**
** Output --
**      indx:   a list of indices for the first strata, then for the second,
**                 and so on.  For each strata, for each time value, the
**                 index of the last survival time with stime <= time is
**                 entered.  (Remember that S uses indices starting at 1, not
**                 zero, and enter the appropriate S index).
**               If time is less than any survival time within the strata,
**                 then the index of the first obs of the next strata is
**                 returned. The proper thing to print will be
**                 a survival value of 1.
**               If there are no survival times >=time, a -1 is entered.
**                 There is  no extrapolation beyond the end of
**                 the K-M curve.
**      indx2:  when =1 indicates a time less than any survival.  When =2
**                 indicates an exact tie.
*/
#include "survS.h"
#include "survproto.h"

void survindex2(int   *n,     double *stime,   int   *strata, 
		int   *ntime, double *time,    int   *nstrat, 
		int   *indx,  int   *indx2)
    {
    register int i,j;
    int nn;
    int current_strata;
    double start_time;

    current_strata = strata[0];
    nn=0;
    start_time = -1;
    j=0;

    for (i=0; i< *nstrat * *ntime; i++) indx[i] = -1;

    for (i=0; i<*n; i++) {
	if (strata[i] != current_strata) {
	    start_time= -1;
	    current_strata = strata[i];
	    nn += *ntime- j;
	    j=0;
	    }
	for (; j< *ntime && time[j] <= stime[i]; j++) {
	    if (start_time < time[j]) {
		if (time[j] < stime[i]) {
		    if (start_time >0) indx[nn++] =i;
		    else {
			indx[nn] =i+1;
			indx2[nn++] = 1;
			}
		    }
		else {
		    indx2[nn] = 2;
		    indx[nn++] = i+1;
		    }
		}
	    }
	start_time = stime[i];
	}
    }
