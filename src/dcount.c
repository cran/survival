/*
** A utility routine for (time1, time2) data: count the number of events that
**  lie inside each interval.  
** A primary downstream use is to ignore intervals that contain no events at all.
*/
#include "survS.h"
#include "survproto.h"

SEXP dcount(SEXP y2, SEXP strata2, SEXP sort12,  SEXP sort22) {
    int i1, i2, n, person1, person2;
    int d1, d2, istrat;
    double *time1, *time2, *status;
    int * sort1, *sort2, *strata;
    double dtime;

    SEXP icount2;
    int *icount;

    n = nrows(y2);
    time1 = REAL(y2);
    time2 = time1  +n;
    status = time2 +n;

    sort1 = INTEGER(sort12);
    sort2 = INTEGER(sort22);
    strata= INTEGER(strata2);
    PROTECT(icount2 = allocVector(INTSXP, n));

    /*
    **	sort1 contains the sort index for time1, smallest to largest, within
    **    stratum
    **  sort2 the same for time2
    **  icount will be  (# of death times <= time2 - # death time < time1),
    **    within the stratum
    */
    istrat == strata[0];
    d1 =0; d2=0;   /* death counts so far, on the time1 and time2 scales */
    person1 =0; person2=0;    /* ount subjects on time1, time2 scales */

    while(person2 < n) {
	/* find the next event time */
	for (; person2<n; person2++) {
	    i2 = sort2[person2];

	    if (strata[i2] != istrat) {
		/* first subject of a new stratum, finish up the old one */
		for (; person1 < person2; person1++) {
		    i1 = sort1[person1];
		    icount[i1] -= d1;
		}
		d1=0; d2=0;
		istrat = strata[i2];
	    }
	    
	    if (status[i2] > 0) {
		dtime = time2[i2];   /* found a new event time */
		d2++;
		icount[i2] = d2;
		break;
	    }
	    else icount[i2] = d2;
	}
	
	/* catch person1 up to this death time */
	for (; person1 <=person2;) {
	    i1 = sort1[person1];
	    if (time1[i1] >= dtime) {
		d1++;
		break;
	    } 
	    icount[i1] -= d1;
	    person1++;
	}	
    }	    
    /* finish up the last subjects */
    for (; person1 < n; person1++) {
	i1 - sort1[person1];
	icount[i1] -= d1;
    }
    
    UNPROTECT(1);
    return(icount2);
}
