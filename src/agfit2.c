/*  SCCS @(#)agfit2.c	5.2 10/27/98
** Anderson-Gill formulation of the cox Model
**
**  the input parameters are
**
**       maxiter      :number of iterations
**       nused        :number of people
**       nvar         :number of covariates
**       start(n)     :each row covers the time interval (start,stop]
**       stop(n)      :
**       event(n)     :was there an event at 'stop':1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       offset(n)    :linear offset
**       weights(n)   :case weights
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       tol_chol     : tolerance for the Cholesky routine
**
**  returned parameters
**       means(nv)    :column means of the X matrix
**       beta(nv)     :the vector of answers (at start contains initial est)
**       u(nv)        :score vector
**       imat(nv,nv)  :the variance matrix at beta=final, also a ragged array
**                      if flag<0, imat is undefined upon return
**       loglik(2)    :loglik at beta=initial values, at beta=final
**       sctest       :the score test at beta=initial
**       flag         :success flag  1000  did not converge
**                                   1 to nvar: rank of the solution
**       maxiter      :actual number of iterations used
**
**  work arrays
**       score(n)              the score exp(beta*z)
**       end(n)                how far to look
**       a(nvar)
**       a2(nvar)
**       cmat(nvar,nvar)       ragged array
**       cmat2(nvar,nvar)
**       newbeta(nvar)         always contains the "next iteration"
**
**  the 6 arrays score, a, cmat, and newbeta are passed as a single
**    vector of storage, and then broken out.
**
**  calls functions:  cholesky2, chsolve2, chinv2
**
**  the data must be sorted by ascending time within strata, deaths before
**          living within tied times.
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

void agfit2( int   *maxiter,  int   *nusedx,  int   *nvarx, 
	     double *start,    double *stop,    int   *event, 
	     double *covar2,   double *offset,  double *weights,
	     int   *strata,   double *means,   double *beta, 
	     double *u,        double *imat2,   double loglik[2], 
	     int   *flag,     double *work,    int   *end,
	     double *eps,      double *tol_chol, double *sctest)
{
    int i,j,k,person;
    int     iter;
    int     nused, nvar;
    int     endp;

    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *newbeta;
    double *a2, **cmat2;
    double *score;
    double  denom, zbeta, risk;
    double  time;
    double  temp, temp2;
    double  newlk=0; /*-Wall*/
    int     halving;    /*are we doing step halving at the moment? */
    double     method;
    double  meanwt;
    int itemp, deaths;
    double efron_wt, d2;

    nused = *nusedx;
    nvar  = *nvarx;
    method= *sctest;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(covar2, nused, nvar);
    imat = dmatrix(imat2,  nvar, nvar);
    cmat = dmatrix(work,   nvar, nvar);
    cmat2= dmatrix(work + nvar*nvar, nvar, nvar);
    a = work + 2*nvar*nvar;
    a2= a+nvar;
    newbeta = a2 + nvar;
    score   = newbeta + nvar;

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	means[i] = temp;
	for (person=0; person<nused; person++) covar[i][person] -=temp;
	}

    /*
    ** do the initial iteration step
    */
    loglik[1] =0;
    for (i=0; i<nvar; i++) {
	u[i] =0;
	for (j=0; j<nvar; j++)
	    imat[i][j] =0 ;
	}

    for (person=0; person<nused; person++) {
	zbeta = 0;      /* form the term beta*z   (vector mult) */
	for (i=0; i<nvar; i++)
	    zbeta += beta[i]*covar[i][person];
	score[person] = zbeta + offset[person];
        }

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
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    end[person] = k;     /*speed up -- last obs in risk set */
		    risk = exp(score[k]) * weights[k];
		    denom += risk;
		    for (i=0; i<nvar; i++) {
			a[i] += risk*covar[i][k];
			for (j=0; j<=i; j++)
			    cmat[i][j] += risk*covar[i][k]*covar[j][k];
			}
		    if (stop[k]==time && event[k]==1) {
			deaths += event[k];
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
	    ** Add results into u and imat for all events at this time point
	    */
	    meanwt /= deaths;
	    itemp = -1;
	    endp = end[person];
	    for (k=person; k<= endp && stop[k]==time; k++) {
		if (event[k]==1) {
		    itemp++;
		    temp = itemp*method/deaths;
		    d2 = denom - temp*efron_wt;
		    loglik[1] +=  weights[k]*score[k] -meanwt *log(d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/d2;
			u[i] += weights[k]*covar[i][k] - meanwt*temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] += meanwt* (
					(cmat[i][j] - temp*cmat2[i][j])/d2-
					   temp2*(a[j]-temp*a2[j])/d2);
			}
		    }
		person++;
		}
	    }
	}   /* end  of accumulation loop */

    loglik[0] = loglik[1];   /* save the loglik for iteration zero  */

    /* am I done?
    **   update the betas and test for convergence
    */

    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    *flag = cholesky2(imat, nvar, *tol_chol);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */

    *sctest=0;
    for (i=0; i<nvar; i++)
	*sctest +=  u[i]*a[i];

    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar; i++) {
	newbeta[i] = beta[i] + a[i];
	}
    if (*maxiter==0) {
	chinv2(imat,nvar);
	for (i=1; i<nvar; i++)
	    for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	return;   /* and we leave the old beta in peace */
	}

    /*
    ** here is the main loop
    */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=1; iter<=*maxiter; iter++) {
	newlk =0;
	for (i=0; i<nvar; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++)
		imat[i][j] =0;
	    }

	for (person=0; person<nused; person++) {
	    zbeta = 0;      /* form the term beta*z   (vector mult) */
	    for (i=0; i<nvar; i++)
		zbeta += newbeta[i]*covar[i][person];
	    score[person] = zbeta + offset[person];
	    }

	for (person=0; person<nused; ) {
	    if (event[person]==0) person++;
	    else {
		endp = end[person];  /* shorter loops than 1 to n */
		/*
		** compute the mean and covariance over the risk set (a and c)
		*/
		efron_wt =0;
		denom =0;
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
		for (k=person; k<=endp; k++) {
		    if (start[k] < time) {
			risk = exp(score[k]) * weights[k];
			denom += risk;
			for (i=0; i<nvar; i++) {
			    a[i] += risk*covar[i][k];
			    for (j=0; j<=i; j++)
				cmat[i][j] += risk*covar[i][k]*covar[j][k];
			    }
			if (stop[k]==time && event[k]==1) {
			    deaths += event[k];
			    efron_wt += risk*event[k];
			    meanwt += weights[k];
			    for (i=0; i<nvar; i++) {
				a2[i]+= risk*covar[i][k];
				for (j=0; j<=i; j++)
				    cmat2[i][j] += risk*covar[i][k]*covar[j][k];
				}
			    }
			}
		    }

		itemp = -1;
		meanwt /= deaths;
		for (k=person; k<=endp && stop[k]==time; k++) {
		    if (event[k]==1) {
			itemp++;
			temp = itemp*method/deaths;
			d2 = denom - temp*efron_wt;
			newlk +=  weights[k]*score[k] -meanwt *log(d2);
			for (i=0; i<nvar; i++) {
			    temp2 = (a[i] - temp*a2[i])/d2;
			    u[i] += weights[k]*covar[i][k] - meanwt*temp2;
			    for (j=0; j<=i; j++)
				imat[j][i] += meanwt*(
					 (cmat[i][j] - temp*cmat2[i][j])/d2-
					       temp2*(a[j]-temp*a2[j])/d2);
			    }
			}

		    person++;
		    }
		}
	    }   /* end  of accumulation loop */

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky2(imat, nvar, *tol_chol);

	if (fabs(1-(loglik[1]/newlk))<=*eps ) { /* all done */
	    loglik[1] = newlk;
	    chinv2(imat, nvar);     /* invert the information matrix */
	    for (i=1; i<nvar; i++)
		for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	    for (i=0; i<nvar; i++)
		beta[i] = newbeta[i];
	    if (halving==1) *flag= 1000; /*didn't converge after all */
	    *maxiter = iter;
	    return;
	    }

	if (iter==*maxiter) break;  /*skip the step halving and etc */

	if (newlk < loglik[1])   {    /*it is not converging ! */
		halving =1;
		for (i=0; i<nvar; i++)
		    newbeta[i] = (newbeta[i] + beta[i]) /2; /*half of old increment */
		}
	    else {
		halving=0;
		loglik[1] = newlk;
		chsolve2(imat,nvar,u);

		j=0;
		for (i=0; i<nvar; i++) {
		    beta[i] = newbeta[i];
		    newbeta[i] = newbeta[i] +  u[i];
		    }
		}
	}   /* return for another iteration */

    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=1; i<nvar; i++)
	for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
    for (i=0; i<nvar; i++)
	beta[i] = newbeta[i];
    *flag= 1000;
    }
