/*  SCCS @(#)agfit_null.c	5.2 10/27/98*/
/*
** Fit a "null" model.  We just need the loglik
**
** Input
**      n       number of subjects
**      method  ==1 for efron
**      start   start times
**      stop    stop times
**      event   =1 if there was an event at time 'stop'
**      offset  the vector of linear predictors
**      weights case weights
**      strata  is =1 for the last obs of a strata
**
** Output
**      loglik  (Breslow approx)
**
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

void agfit_null(int   *n,      int   *method,   double *start, double *stop, 
		int   *event,  double * offset,  double *weights,
		int   *strata, double loglik[2])
    {
    int i,k;
    double denom;
    double e_denom;
    double temp;
    double time;
    int deaths;
    double itemp;
    double meanwt;

    loglik[0]=0;
    for (i=0; i<*n; ) {
	if (event[i] ==1) {
	    /*
	    ** compute the sum of weights over the risk set
	    **   and count the deaths
	    */
	    denom =0;
	    e_denom =0;
	    meanwt =0;
	    deaths =0;
	    time = stop[i];
	    for (k=i; k<*n; k++) {
		if (start[k] < time) denom += exp(offset[k]);
		if (stop[k]==time && event[k]==1) {
		    deaths ++;
		    e_denom += exp(offset[k]) * weights[k];
		    loglik[0] += offset[k] * weights[k];
		    meanwt += weights[k];
		    }
		if (strata[k]==1) break;
		}

	    itemp =0;
	    meanwt /= deaths;
	    for (k=i; k<*n && stop[k]==time; k++) {
		if (event[k]==1) {
		    temp = *method * itemp / deaths;
		    loglik[0] -= meanwt *log(denom - temp*e_denom);
		    itemp++;
		    }
		i++;
		if (strata[k]==1) break;
		}
	    }
	else i++;
	}
    }
