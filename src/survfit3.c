/* SCCS @(#)survfit3.c	1.3 07/09/00
** Fit the survival curve
**  Input
**    sn = number of subjects
**    y[ny,n] =  matrix of time and status values
**    wt[n] = vector of case weights
**    strata[n] = is equal to 1 at the last obs of each strata
**    method = 1 = km, 2 = fleming-harrington
**    error = 1 = Greenwood, 2=Tsiatis
**    nstrat = number of strata
**    ntimes_strata = number of times in each strata
**    timelist = list of times
**    weighted_event[n], surv[n], varh[n], risknum[n], enter[n],
**    exit_censored[n] = work arrays
**
** Output
**    strata[] = number of observations in each strata
**    y[] = contains the survival times
**    mark[] = number of (weighted) events at each time point 
**    surv[] = the survival for each time point
**    varh[] = the variance of the hazard function for each time point
**    risksum[] = total (weighted) number of observations still in the
**                study at each time point
**    enter[] = number of patients who entered the study at each time point
**    exit_censored[] = number of patients with no events at each time point
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"
void survfit3(int   *sn,              double *y,               double *wt,
	      int   *strata,          int   *method,          int   *error,
	      int   *nstrat,          double *ntimes_strata,  
	      double *timelist,	       double *weighted_event,  double *surv,
	      double *varh,	       double *risksum,         double *enter,
	      double *exit_censored)

{
  int i,j=0, grp, k;/*-Wall*/
  double hazard, varhaz;
  double km, c;
  double enter_count, exit_censored_count, risk_count, 
         event_count, weighted_event_count;
  double *status, *start_time, *stop_time;
  double temp, current_time;
  int n, pgm_nstrat, timelist_start, timeloop;
  int nsurv, loop_nstrat;

  n = *sn;
  pgm_nstrat = *nstrat;
  start_time = y;
  stop_time = y+n;
  status = y+(2*n);
  strata[n-1] =1;   /*just in case the parent routine forgot */

  timelist_start = 0;
  timeloop = 0;  /** keep track of total list of times **/

  loop_nstrat = 0;
  for (grp=0; grp<pgm_nstrat; grp++) {

    /** reset survival variables for each group **/
    km = 1;
    hazard = 0;
    varhaz = 0;

    for (i=0; i<ntimes_strata[grp]; i++) {
      /** reset counts for each time **/
      enter_count = 0; 
      exit_censored_count = 0; 
      event_count = 0; 
      risk_count = 0; 
      weighted_event_count = 0;

      current_time = timelist[timeloop];
      /** for each start and stop time, look at the current time from **/
      for (j=timelist_start; j<n; j++) {
	/** count number at risk at current time **/
	if (start_time[j] < current_time && stop_time[j] >= current_time)
	  risk_count += wt[j];
	  
	/** count number entered at current time **/
	if (start_time[j] == current_time)
	  enter_count++;

	/** checking where to start looking at start and stop times **/
	/** for next current_time **/
	if (stop_time[j] <= current_time)
	  timelist_start++;

	/** counting number of events or non-events at current time **/
	if (stop_time[j] == current_time) {
	  if (status[j] == 0)
	    exit_censored_count++;
	  else {
	    event_count++;
	    weighted_event_count += status[j]*wt[j];
	    }
	  }

	/** checking if at end of times for current group **/
	if (strata[j] == 1)
	  break;
        }

      if (weighted_event_count > 0) {
	if (*method==1) {
	  km *= (risk_count - weighted_event_count) / risk_count;
	  if (*error==1)
	    varhaz += weighted_event_count / 
	              (risk_count * (risk_count - weighted_event_count));
	  else 
	    varhaz += weighted_event_count / 
	              (risk_count * risk_count);
	  }
	else  if (*method==2) {
	  hazard += weighted_event_count / risk_count;
	  km = exp(-hazard);
	  if (*error==1)
	    varhaz += weighted_event_count / 
	              (risk_count * (risk_count - weighted_event_count));
	  else
	    varhaz += weighted_event_count / 
	              (risk_count * risk_count);
	  }
	else  if (*method==3) {
	  for (k=0; k<weighted_event_count; k++) {
	    c = weighted_event_count / event_count;
	    temp = risk_count - k * c;
	    hazard += 1 / temp;
	    if (*error==1)
	      varhaz += 1 / (temp * (risk_count - (k+1)*c));
	    else
	      varhaz += 1 / (temp * temp);
  	    }
	  km = exp(-hazard);
	  }
	
	weighted_event[timeloop] = weighted_event_count;
	risksum[timeloop] = risk_count;
	enter[timeloop] = enter_count;
	exit_censored[timeloop] = exit_censored_count;
	surv[timeloop] = km;
	varh[timeloop] = varhaz;
	
        }
      else {
	if (i == 0) {
	  weighted_event[timeloop] = 0;
	  risksum[timeloop] = 0;
	  enter[timeloop] = enter_count;
	  exit_censored[timeloop] = exit_censored_count;
	  surv[timeloop] = 1;
	  varh[timeloop] = 0;
	  }
	else {
	  weighted_event[timeloop] = weighted_event_count;
	  risksum[timeloop] = risk_count;
	  enter[timeloop] = enter_count;
	  exit_censored[timeloop] = exit_censored_count;
	  surv[timeloop] = surv[timeloop-1];
	  varh[timeloop] = varh[timeloop-1];
	  }
        }
      /** advance one in the timelist - points to new current time **/
      timeloop++;
      }
    timelist_start = ++j;
    }
  
  for (nsurv=0; nsurv<n; nsurv++) {
    if (strata[nsurv] == 1 ) {
      strata[loop_nstrat] = nsurv;
      loop_nstrat ++;
      }
    }  
  }











