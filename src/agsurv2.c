/*  $Id: agsurv2.c 11080 2008-10-24 03:47:51Z therneau $*/
/*
** Fit the survival curve, the special case of an Anderson-Gill style data
**   This program differs from survfit in several key ways:
**       Only returns data at the event times, not at censoring times
**       Fewer work arrays, but it is slower.
**
**   This is similar to survfit in that a complete curve is produced for
**        each strata.  If there are multiple 'subjects' in the newdata
**        list, then a matrix of survival curves is produced, with one
**        column for each input vector.
**
**  Input
**    n=# of subjects
**    nvar - number of vars in xmat -- will be zero if se is not desired
**             (the calling routine knows that xmat is only needed for the
**              correct second term of the se)
**    y - 3 column matrix containing strart, stop, event
**    score[n] - vector of rscores
**    strata[n] - ==1 at the last obs of each strata
**    xmat   = data matrix that generated the Cox fit
**    varcov   = variance matrix of the coefs
**    nsurv[2] = the methods for estimate and variance:
**		1=Kalbfleisch/Prentice  2= Aalen/Breslow
**                              3= Tsiatis, Efron approx
**
**    ncurve  = # of curves to produce
**    newx(nvar, ncurve) =  new subject x matrix
**    newrisk(ncurve)  = rscores for the new subjects
**
** Output
**    surv  - the survival - of length ncurve*nsurv
**    varh  - the variance of the hazard function
**    nsurv - returned, number of survival time points
**    y[1,] - contains the survival times
**    y[2,] - the number of subjects at risk at that survival time
**    y[3,]  - the number of events at that time
**    strata[0]= # of strata, strata[1:n]= last obs strata 1,2, etc
**
**  Work
**    d[3*nvar]
**
**  Input must be sorted by (event before censor) within stop time within strata,
*/
#include <stdio.h>
#include <math.h>
#include "survS.h"
#include "survproto.h"

void agsurv2(Sint   *sn,      Sint   *snvar,    double *y, 
	     double *score,   Sint   *strata,   double *wt, double *surv, 
	     double *varh,    double *xmat,     double *varcov, 
	     Sint   *snsurv,  double *d,        Sint   *sncurve,
             double *newx,    double *newrisk)
{
    int i,j,k,l;
    double hazard, varhaz;
    double *start, *stop, *event;
    int n, nvar;
    int nsurv, type, vartype;
    int kk=0, psave;
    double deaths;
    double *a, *a2;
    int ncurve;
    double **covar,
	   **imat,
	   **covar2;
    int column,
	nrisk,
	nstrat,
	nsave,
	person;
    double time,
	   rscore =0,
	   e_denom,
	   denom;
    double crisk,
	   guess, inc,
	   sumt,
	   km;
    double temp,
	   downwt,
	   d2;

    n = *sn;  nvar = *snvar;
    ncurve = *sncurve;
    type = snsurv[0];  vartype=snsurv[1];
    start =y;
    stop  = y+n;
    event = y+n+n;
    a = d+nvar;
    a2 = a+nvar;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(xmat, n, nvar);
    imat = dmatrix(varcov,  nvar, nvar);
    covar2 = dmatrix(newx, ncurve, nvar);
    nsurv =0;
    nstrat =0;
    for (column=0; column<ncurve; column++) {
	crisk = newrisk[column];
	hazard  =0;
	varhaz  =0;
	km =1;
	for (i=0; i<nvar; i++) d[i] =0;
	nsave = nsurv;
	for (person=0; person<n;) {
	    if (event[person]==0) person++;
	    else {
		/*
		** compute the mean and denominator over the risk set
		*/
		denom =0;
		e_denom=0;
		for(i=0; i<nvar; i++){
		    a[i] =0;
		    a2[i]=0;
		    }
		time = stop[person];
		nrisk =0;
		deaths=0;
		for (k=person; k<n; k++) {
		    if (start[k] < time) {
			nrisk++;
			rscore = wt[k] *score[k]/crisk;
			denom += rscore;
			for (i=0; i<nvar; i++) {
			    a[i] += rscore*(covar[i][k]- covar2[i][column]);
			    }
			 }
		    if (stop[k]==time && event[k]==1) {
			deaths += wt[k];
			e_denom += rscore;
			for (i=0; i<nvar; i++) {
			    a2[i] += rscore*(covar[i][k]- covar2[i][column]);
			    }
			}
		    if (strata[k]==1) break;
		    }

		/*
		** Add results all events at this time point
		*/
		psave = person;  /* for KM case below */
		temp =0;
		for (k=person; k<n && stop[k]==time; k++) {
		    if (event[k]==1) {
			kk =k ;      /*save for km case */
			downwt = temp/deaths;
			if (type==3) {
			    d2 = (denom - downwt*e_denom);
			    hazard += wt[k]/d2;
			    }
			else  hazard += wt[k]/denom;
			
			if (vartype==3) {
			    d2 = (denom - downwt*e_denom);
			    varhaz += wt[k]/(d2*d2);
			    for (i=0; i<nvar; i++)
				d[i] += wt[k]*(a[i]- downwt*a2[i])/ (d2*d2);
			    }
			else {
			    if (vartype==2) varhaz += wt[k]/(denom*denom);
			    for (i=0; i<nvar; i++)
				d[i] += wt[k]* a[i]/(denom*denom);
			    }
			temp++;
			}
		    person++;
		    if (strata[k]==1) break;
		    }

		if (vartype==1) {
		    if (denom >e_denom) 
			varhaz += deaths/(denom*(denom-e_denom));
		    else varhaz +=deaths; /* a hack to keep from zero divide */
		    }
		if (type==1) {
		    /*
		    ** kalbfleisch estimator is harder
		    **   (But still using score/crisk as a "new" score).
		    */
		    if (deaths ==1) {
			km *= pow(1- wt[kk]*score[kk]/(crisk*denom), 
				     crisk/score[kk]);
			}
		    else {           /*find the zero of an equation */
			guess = .5;
			inc = .25;
			for (l=0; l<35; l++) { /* bisect it to death */
			    sumt =0;
			    for (k=psave; k<person; k++) {
				if (event[k] ==1) {
				    temp = score[k]/crisk;
				    sumt +=  wt[k]*temp/(1-pow(guess, temp));
				    }
				}
			    if (sumt < denom)  guess += inc;
				 else          guess -= inc;
			    inc = inc/2;
			    }
			km *= guess;
			}
		    surv[nsurv] = km;;
		    }
		else surv[nsurv] = exp(-hazard);

		temp =0;
		for (i=0; i<nvar; i++)
		    for (j=0; j< nvar; j++)
			temp += d[i]*d[j]*imat[i][j];

		varh[nsurv] = varhaz + temp;
		if (column==(ncurve-1)) {
		    /* on the last pass, I can overwrite old data with new */
		    i = nsurv - nsave;
		    start[i] = time;
		    stop[i] = nrisk;
		    event[i]= deaths;
		    }
		nsurv++;
		}

	    if (strata[person-1]==1) {
		if (column==(ncurve-1)) {
		    nstrat++;
		    strata[nstrat]= nsurv-nsave;
		    }
		km=1;
		hazard  =0;
		varhaz  =0;
		for (i=0; i<nvar; i++) d[i] =0;
		}
	    }
	}
    *snsurv = nsurv/ ncurve;
    strata[0] = nstrat;
    }



