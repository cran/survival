/*  SCCS @(#)coxdetail.c	5.2 10/27/98
** Return all of the internal peices of a cox model
**
**  the input parameters are
**
**       nused        :number of people
**       nvar         :number of covariates
**       ndead        :total number of deaths
**       y(3,n)       :start, stop, and event for each subject
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       score(n)     :the risk score for the subject
**       weights(n)   :case weights
**       means        :first element contains the method
**
**  returned parameters
**       ndead        :the number of unique death times
**       score        :the indices of the unique time points
**       y[1, ]       :the number of deaths at each time point
**       y[2, ]       :the number at risk at each time point
**       y[3, ]       :the increment in the cum -hazard at t
**       weights[]    :the increment in the variance of the cum-haz at t
**       means(nv,nd) :the matrix of weighted means, one col per unique event
**                                              time
**       u(nv,nd)     :the score vector components, one per unique event time
**       var(nd,nv,nv):components of the information matrix
**
**  work arrays
**       a(nvar)
**       a2(nvar)
**       cmat(nvar,nvar)       ragged array
**       cmat2(nvar,nvar)
**       wmeans(nvar)
**
**  the 5 arrays a, a2, cmat, cmat2 and wmeans are passed as a single
**    vector of storage, and then broken out.
**
**  the data must be sorted by ascending time within strata, deaths before
**          living within tied times.
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

void coxdetail(int   *nusedx,   int   *nvarx,    int   *ndeadx, 
	       double *y,        double *covar2,   int   *strata,  
	       double *score,    double *weights,  double *means2, 
	       double *u2,       double *var,      double *work)
{
    register int i,j,k,person;
    int     nused, nvar;
    int     nrisk, ndead;
    double **covar, **cmat;    /*ragged arrays */
    double **means;
    double **u;
    double *a;
    double *a2, **cmat2;
    double *wmeans;
    double  denom;
    double  time;
    double  temp, temp2, temp3;
    double     method;
    double  hazard;
    double  varhaz;
    int     itemp, deaths;
    int     ideath;
    double  efron_wt, d2;
    double risk;
    double  meanwt;
    double  *start,
	    *stop,
	    *event;

    nused = *nusedx;
    nvar  = *nvarx;
    method= *means2;
    ndead = *ndeadx;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(covar2, nused, nvar);
    means= dmatrix(means2, ndead, nvar);
    u    = dmatrix(u2, ndead, nvar);
    cmat = dmatrix(work,   nvar, nvar);
    cmat2= dmatrix(work + nvar*nvar, nvar, nvar);
    a = work + 2*nvar*nvar;
    a2= a+nvar;
    wmeans = a2+nvar;
    start =y;
    stop  =y + nused;
    event =y + nused +nused;

    /*
    ** Subtract the mean from each covar, as this makes the variance calc
    **  much more stable
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	wmeans[i] = temp;
	for (person=0; person<nused; person++) covar[i][person] -=temp;
	}

    /*
    ** Zero out some arrays
    */
    for (i=0; i<ndead*nvar; i++) {
	u2[i]=0;
	means2[i] =0;
	}
    for (i=0; i<ndead*nvar*nvar; i++) var[i]=0;

    /*
    ** Now walk through the data
    */
    ideath=0;
    for (person=0; person<nused;) {
	if (event[person]==0) person++;
	else {
	    /*
	    ** compute the mean and covariance over the risk set (a and c)
	    */
	    denom =0;
	    efron_wt =0;
	    meanwt =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		a2[i]=0;
		for (j=0; j<nvar; j++) {
		    cmat[i][j]=0;
		    cmat2[i][j]=0;
		    }
		}
	    time = stop[person];
	    deaths=0;
	    nrisk =0;
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    nrisk++;
		    risk = score[k] * weights[k];
		    denom += risk;
		    for (i=0; i<nvar; i++) {
			a[i] += risk*covar[i][k];
			for (j=0; j<=i; j++)
			    cmat[i][j] += risk*covar[i][k]*covar[j][k];
			}
		    if (stop[k]==time && event[k]==1) {
			deaths += 1;
			efron_wt += risk*event[k];
			meanwt += weights[k];
			for (i=0; i<nvar; i++) {
			    a2[i]+= risk*covar[i][k];
			    for (j=0; j<=i; j++)
				cmat2[i][j] += risk*covar[i][k]*covar[j][k];
			    }
			}
		     }
		if (strata[k]==1) break;
		}

	    /*
	    ** Add results into u and var for all events at this time point
	    */
	    itemp = -1;
	    hazard =0;
	    varhaz =0;
	    meanwt /= deaths;
	    for (k=person; k<nused && stop[k]==time; k++) {
		if (event[k]==1) {
		    itemp++;
		    temp = itemp*method/deaths;
		    d2 = denom - temp*efron_wt;
		    hazard += meanwt/d2;
		    varhaz += meanwt*meanwt/(d2*d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/d2;
			means[i][ideath] += (wmeans[i] +temp2)/deaths;
			u[i][ideath] += weights[k]*covar[i][k] - meanwt*temp2;
			for (j=0; j<=i; j++) {
			    temp3 =((cmat[i][j] - temp*cmat2[i][j]) -
					       temp2*(a[j]-temp*a2[j]))/d2;
			    temp3 *= meanwt;
			    var[i + j*nvar + ideath*nvar*nvar] +=temp3;
			    if (j<i)
				var[j + i*nvar + ideath*nvar*nvar] +=temp3;
			    }
			}
		    }
		person++;
		if (strata[k]==1) break;
		}
	    score[ideath] = person;
	    start[ideath] = deaths;
	    stop[ideath]  = nrisk;
	    event[ideath] = hazard;
	    weights[ideath]=varhaz;
	    ideath++;
	    }
	}   /* end  of accumulation loop */
    *ndeadx = ideath;
    }
