/* SCCS @(#)survfit2.c	5.5 07/09/00
** Fit the survival curve
**  Input
**    sn = number of subjects
**    y[ny,n] =  matrix of time and status values
**    wt[n] = vector of case weights
**    strata[n] = is equal to 1 at the last obs of each strata
**    method = 1= km, 2= fleming-harrington
**    error = 1= Greenwood, 2=Tsiatis
**    mark[n], surv[n], varh[n], risksum[n] = work arrays
**
** Output
**    y[] = contains the survival times
**    strata[] = number of observations in each strata
**    mark[] = number of (weighted) events at each time point
**    surv[] = the survival for each time point
**    varh[] = the variance of the hazard function for each time point
**    risksum[] = total (weighted) number of observations still in the 
**                study at each time point
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"
void survfit2(int   *sn,      double *y,       double *wt,
	      int   *strata,  int   *method,  int   *error, 
	      double *mark,    double *surv,    double *varh, 
	      double *risksum  )
{
    int i,j;
    double hazard, varhaz;
    double sum=0, km;/*-Wall*/
    double *time, *status;
    double temp;
    int n;
    int nsurv, nstrat;

    n = *sn;
    time =y;
    status = y+n;
    strata[n-1] =1;   /*just in case the parent routine forgot */

    /*
    **  initialize a couple of arrays
    **    mark(i) contains the number of deaths at this particular time point
    **    risksum contains the running # at risk
    */
    for (i=0; i<(n-1); i++)  /*first pass, just mark the ties */
	if (time[i]==time[i+1] && strata[i]==0)  mark[i+1]= -1;
	else                                     mark[i+1]=  1;
    mark[0]=1;

    temp =0;                  /* second pass */
    for (i=n-1; i>=0; i--) {
	if (strata[i]==1) sum =0;
	sum  += wt[i];
	temp += status[i] * wt[i];

	if (mark[i] == 1) {
	    mark[i] = temp;
	    risksum[i] = sum;
	    temp =0;
	    }
	}

    /*
    ** the hazard starts out at zero;
    */
    nsurv=0;
    nstrat=0;
    km =1;
    hazard  =0;
    varhaz  =0;
    for(i=0; i<n; i++) {
	if (mark[i] >0) {
	    if (*method==1) {
		km *= (risksum[i] - mark[i]) / risksum[i];
		if (*error==1 )
		     varhaz += mark[i]/(risksum[i]*(risksum[i]-mark[i]));
		else varhaz += mark[i]/(risksum[i]*risksum[i]);
		}
	    else  if (*method==2) {
		hazard += mark[i]/risksum[i];
		km = exp(-hazard);
		if (*error==1 )
		     varhaz += mark[i]/(risksum[i]*(risksum[i]-mark[i]));
		else varhaz += mark[i]/(risksum[i]*risksum[i]);
		}

	    else  if (*method==3) {
		for (j=0; j<mark[i]; j++) {
		    temp = risksum[i] -j;
		    hazard += 1/temp;
		    if (*error==1 )
			 varhaz += 1/(temp*(temp-1));
		    else varhaz += 1/(temp*temp);
		    }
		km = exp(-hazard);
		}
	    }

	if (mark[i] >=0) {
	    time[nsurv] = time[i];
	    mark[nsurv] = mark[i];
	    risksum[nsurv] = risksum[i];
	    surv[nsurv] = km;;
	    varh[nsurv] = varhaz;
	    nsurv++;
	    }

	if (strata[i]==1) {
	    strata[nstrat]= nsurv;
	    nstrat++;
	    if (nsurv<n) {
		surv[nsurv] =1;
		varh[nsurv] = 0;
		}
	    km=1;
	    hazard  =0;
	    varhaz  =0;
	    }
	}
    }

