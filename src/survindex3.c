/*    SCCS @(#)survindex3.c	1.3 07/09/00
** A subroutine for summary.survfit
**
** Input --
**      n:      number of survival times
**      stime:  survival times, must be >0
**      strata: strata values   data is sorted by time within strata
**      ntime:  number of time values specified
**      time:   time values for printout, must be >=0 and in increasing order
**      nstrat: number of strata
**      o_n_risk: original total (weighted) number of observations still in
**                the study at each time point
**      o_n_entered: original number of people that came into the study at 
**                   each time point
**      o_n_censored: original number of people that left the study without 
**                    an event at each time point
**      o_n_event: original number of (weighted) events at each time point
**      o_surv: original survival at each time point
**      o_std_err: null or original calculated standard error for each time 
**                 point
**      o_upper: null or original upper confidence interval for each time point
**      o_lower: null or original lower confidence interval for each time point
**      n_risk,n_entered, n_censored, n_event, surv, std_err, upper, 
**      lower: work arrays
**      new_start: if user specified a time to start other than first time
**      num_extend: if should calculate points beyond the curve
**      times_strata: number of times for each group, used if num_extend=0
**      temp_times: time values used for each group, used if num_extend=0
**
** Output --
**      n_risk: (weighted) number of observations still in study for specified
**              time points
**      n_entered: number of people that came in the study for the specified
**                 time points
**      n_censored: number of people that left the study without an event for
**                  the specified time points
**      n_event: (weighted) events for specified time points
**      surv: survival for specified time points
**      std_err: standard error for specified time points
**      upper: upper confidence interval for specified time points
**      lower: lower confidence interval for specified time points
*/
#include "survS.h"
#include "survproto.h"
void survindex3(int   *n,          double *stime,        int   *strata,
		int   *ntime,      double *time,         int   *nstrat, 
		int   *o_n_risk,   int   *o_n_entered,  int   *o_n_censored,
		int   *o_n_event,  double *o_surv,       double *o_std_err,
		double *o_upper,    double *o_lower,      int   *n_risk, 
		int   *n_entered,  int   *n_censored,   int   *n_event,
		double *surv,       double *std_err,      double *upper,
		double *lower,      double *new_start, 	  int   *num_extend,
		int   *times_strata,                     double *temp_times)

{
  int i,j,k;
  int nn, cc, current_strata, sum_enter, sum_censor, sum_event, /*last, last_time,-Wall unused */ 
      sum_times, strata_count;

  double start_time, starting_time;

  current_strata = strata[0];
  starting_time = stime[0] - 1;
  start_time = starting_time;
  j=0;
  nn = 0;
  cc = 0;
  sum_enter = 0;
  sum_censor = 0;
  sum_event = 0;
  strata_count = 0;
  sum_times = 0;
  
  for (i=0; i<*n; i++) { /* for i 1 */

    if (strata[i] != current_strata) { /* if 1 */
      starting_time = stime[i] - 1;
      start_time = starting_time;
      current_strata = strata[i];
      times_strata[strata_count] = sum_times;
      sum_enter = 0;
      sum_censor = 0;
      sum_event = 0;
      sum_times = 0;
      j=0;
      strata_count++;
      } /* end if 1 */

    if (stime[i] >= *new_start) {
      sum_enter += o_n_entered[i];
      sum_censor += o_n_censored[i];
      sum_event += o_n_event[i];
    
      for (; j< *ntime && time[j] <= stime[i]; j++) { /* for j 1 */
	temp_times[cc] = time[j];
	sum_times = sum_times + 1;
	if (start_time < time[j]) { /* if 2 */
	  if (time[j] < stime[i]) { /* if 3 */
	    if (start_time > starting_time) { /* if 4 */
	      n_entered[nn] = sum_enter - o_n_entered[i];
	      n_censored[nn] = sum_censor - o_n_censored[i];
	      n_event[nn] = sum_event - o_n_event[i];
	      n_risk[nn] = o_n_risk[i];
	      surv[nn] = o_surv[i-1];
	      std_err[nn] = o_std_err[i-1];
	      upper[nn] = o_upper[i-1];
	      lower[nn] = o_lower[i-1];
	      }
	    else { /* end if 4, start else 4 */
	      n_entered[nn] = 0;
	      n_censored[nn] = 0;
	      n_event[nn] = 0;
	      n_risk[nn] = o_n_risk[i];
	      surv[nn] = 1;
	      std_err[nn] = 0;
	      upper[nn] = 1;
	      lower[nn] = 1;
	      } /* end else 4 */
	    }
	  else { /* end if 3, start else 3 */
	    n_entered[nn] = sum_enter;
	    n_censored[nn] = sum_censor;
	    n_event[nn] = sum_event;
	    n_risk[nn] = o_n_risk[i];
	    surv[nn] = o_surv[i];
	    std_err[nn] = o_std_err[i];
	    upper[nn] = o_upper[i];
	    lower[nn] = o_lower[i];
	    } /* end else 3 */
	  } /* end if 2 */
	  else if (j == 0) {
	    n_entered[nn] = 0;
	    n_censored[nn] = 0;
	    n_event[nn] = 0;
	    n_risk[nn] = 0;
	    surv[nn] = 1;
	    std_err[nn] = 0;
	    upper[nn] = 1;
	    lower[nn] = 1;
	    }
	nn++; cc++;
        } /* end for j 1 */
      
      if (i+1 != *n) { /* if 5 */
	if (strata[i+1] != strata[i]) { /* if 6 */
	  for (k=0; k < *ntime; k++) { /* for k 1 */
	    if (time[k] > stime[i] && (*num_extend == 1)) { /* if 7 */
	      n_entered[nn] = sum_enter;
	      n_censored[nn] = sum_censor;
	      n_event[nn] = sum_event;
	      n_risk[nn] = 0;
	      surv[nn] = o_surv[i];
	      std_err[nn] = o_std_err[i];
	      upper[nn] = o_upper[i];
	      lower[nn] = o_lower[i];
	      nn++;
	      } /* end if 7 */
	    } /* end for k 1 */
	  } /* end if 6 */
        } 
      else { /* end if 5, start else 5 */
	for (k=0; k < *ntime; k++) { /* for k 2 */
	  if (time[k] > stime[i] && (*num_extend == 1)) { /* if 8 */
	    n_entered[nn] = sum_enter;
	    n_censored[nn] = sum_censor;
	    n_event[nn] = sum_event;
	    n_risk[nn] = 0;
	    surv[nn] = o_surv[i];
	    std_err[nn] = o_std_err[i];
	    upper[nn] = o_upper[i];
	    lower[nn] = o_lower[i];
	    nn++;
	    } /* end if 8 */
  	  } /* end for k 2 */
        } /* end else 5 */
      }
    start_time = stime[i];
    } /* end for i 1 */
  times_strata[strata_count] = sum_times;
  }





