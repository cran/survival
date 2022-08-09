/* A reentrant version of the agfit program, for penalized effects modeling
**   with reasonable efficiency (I hope).  The important arrays are saved
**   from call to call so as to speed up the process.  The x-matrix itself
**   is the most important of these.
** This is the version with the "smart sorting" speedup.
**
** agfit5a: Entry and intial iteration step for beta=initial, theta=0
**              (no frailty)
**            Most of the same arguments as agfit2.
**            Allocate and save arrays in static locations.
** agfit5b: Iterate to convergence given an initial value.
** agfit5c: Compute residuals and release the saved memory.
*
**  the input parameters are
**
**       maxiter      :number of iterations
**       nused        :number of people
**       nvar         :number of covariates
**       yy[3,n]      :row 1: start time of event or censoring for person i
**		      :row 2: stop time
**                    :row 3: status for the ith person    1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(nstrat):sizes of the strata, cumulative
**       sort 	     : two column matrix:
**                      sort order for the obs, using stop time, last to first
**                      sort order for the obs, using start time, last to first
**       offset(n)    :offset for the linear predictor
**       weights(n)   :case weights
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       tolerch      :tolerance for the Cholesky routines
**       method       : Method 0=Breslow, 1=Efron
**       ptype        : 1 or 3 -- there is a sparse term
**                    : 2 or 3 -- there is a non-sparse term in the model
**       nfrail       : number of frailty groups (sparse terms), 0 if there are
**                        none
**       frail        : a vector containing the frailty groups
**       fbeta        : initial frailty estimates
**       pdiag        : if 0, then for the non-sparse terms only the diagonal
**                        of the variance matrix is penalized, otherwise the
**                        full matrix is used.
**
**  returned parameters
**       means(nv)    : vector of column means of X
**       beta(nv)     : the vector of answers (at start contains initial est)
**       u(nv)        : score vector
**       imat(nv,nv)  : the variance matrix at beta=final
**                      if flag<0, imat is undefined upon return
**       loglik       :loglik at beta=final
**       flag         :success flag  1000  did not converge
**                                   1 to nvar: rank of the solution
**       maxiter      :actual number of iterations used
**       fbeta(nfrail): fitted frailty values
**       fdiag(nfrail + nvar): diagonal of cholesky of the full inverse
**       jmat         : inverse of the cholesky
**       imat         : cholesky of the information matrix
**       expect       : contains the "expected" for each subject
** 
**  work arrays
**       score(n)     
**       a(nvar+ nfrail), a2(nvar+nfrail)
**       cmat(nvar,nvar+nfrail)       ragged array
**       cmat2(nvar,nvar+nfrail)
**       fdiag                         the diagonal of the sparse information
**       oldbeta(nvar + nfrail)         always contains the "last iteration"
**
**  the work arrays are passed as a single
**    vector of storage, and then broken out.
**
**  calls functions:  cholesky3, chsolve3, chinv2
**
**  the data must be sorted by ascending time within strata
*/
#include <math.h>
#include <stdio.h>
#include "survS.h"
#include "survproto.h"

static double **covar, **cmat, **cmat2;
static double *a, *oldbeta, *a2;
static double *offset, *weights;
static int    *event, *frail = NULL;
static double *score, *start, *stop;
static int    *sort1, *sort2;
static double *tmean;
static int    ptype, pdiag;
static double *ipen, *upen, logpen;
static int    *zflag;

static double **cmatrix(double *, int, int);

void agfit5a(int  *nusedx, int  *nvarx, double *yy, 
	       double *covar2, double *offset2,
	       double *weights2, 
	       int    *strata,  int    *sort,
	       double *means, double *beta, double *u, 
	       double *loglik, 
	       int  *methodx, int  *ptype2, int  *pdiag2,
	       int  *nfrail,  int  *frail2,
               void *fexpr1, void *fexpr2, void *rho,
	       int  *docenter) {

    int i,j,k, person;
    int     nused, nvar;
    int    nf, nvar2;
    int  deaths, itemp;
    int  istrat, indx2, p, ksave;  

    double  denom, zbeta, risk;
    double  temp;
    double  d2, efron_wt;
    double  method;
    double  meanwt, time;

    nused = *nusedx;
    nvar  = *nvarx;
    nf= *nfrail;
    method= *methodx;
    nvar2 = nvar + nf;
    ptype = *ptype2;
    pdiag = *pdiag2;

    /*
    **  Allocate storage for the arrays and vectors
    **  Since they will be used later, sizes are based on what will be
    **    needed with the frailty terms.
    */
    if (nvar >0) {
	covar= cmatrix(covar2, nused, nvar);
	cmat = cmatrix(0, nvar2, nvar+1);
	cmat2= cmatrix(0, nvar2, nvar+1);
        }

    a = CALLOC(4*nvar2 + 5*nused , double);
    oldbeta = a + nvar2;
    a2 =  oldbeta + nvar2;
    weights = a2+ nvar2;
    offset  = weights + nused;
    score   = offset + nused;
    tmean   = score + nused;
    start   = tmean + nvar2;
    stop    = start + nused;
    
    event  = CALLOC(3*nused, int);
    sort1   = event + nused;
    sort2   = sort1 + nused;

    for (i=0; i<nused; i++) {
	weights[i] = weights2[i];
	offset[i]  = offset2[i];
	event[i]  =  yy[nused + nused +i];
	sort1[i]  = sort[i];
	sort2[i]  = sort[nused+i];
	start[i]  = yy[i];
	stop[i]   = yy[nused+i];
        }

    /* scratch space for penalty 
    **    upen needs to be max(nvar, nfrail), 
    **    ipen max(nfrail, nvar(if pdiag=0) or nvar^2 )
    */
    if (nf > nvar) i=nf; else i=nvar;
    if (nf > nvar*nvar) j=nf; else j=nvar*nvar;
    if (pdiag==0)  upen = CALLOC(2*i, double);
    else           upen = CALLOC(i+j, double);
    ipen = upen + i;
    if (ptype>1)  zflag = CALLOC(nvar, int );
    else          zflag = CALLOC(2, int );

    if (nf>0) {
	frail = CALLOC(nused, int);
	for (i=0; i<nused; i++) frail[i] = frail2[i];
        }

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  more stable
    */
    for (i=0; i<nvar; i++) {
	if (docenter[i] ==0) means[i] =0;
	else {
	    temp=0;
	    for (person=0; person<nused; person++) temp += covar[i][person];
	    temp /= nused;
	    means[i] = temp;
	    for (person=0; person<nused; person++) covar[i][person] -=temp;
	}
    }	
	
    /*
    ** Find the loglik of the initial model
    **   (actually, just a no-sparse-terms model) -- loglik only
    */
    *loglik = 0;

    for (person=0; person<nused; person++) {
        zbeta = 0;      /* form the term beta*z   (vector mult) */
        for (i=0; i<nvar; i++)
            zbeta += beta[i]*covar[i][person];
        score[person] = coxsafe(zbeta + offset[person]);  /* save this away */
        }

    /*
    **  'person' walks through the stop times from largest to smallest
    **     (sort1[0] points to the largest stop time, sort1[1] the next, ...)
    **  'time' is the time of current interest, which also goes from smallest
    **      to largest.  
    **  'indx2' walks through the start times.  It will be smaller than 
    **    'person': if person=27 that means that 27 subjects have stop >=time,
    **    and are thus potential members of the risk set.  If 'indx2' =9,
    **    that means that 9 subjects have start >=time and thus are NOT part
    **    of the risk set.  (stop > start for each subject guarrantees that
    **    the 9 are a subset of the 27). 
    **  Basic algorithm: move 'person' forward, adding the new subject into
    **    the risk set.  If this is a new, unique death time, take selected
    **    old obs out of the sums, add in obs tied at this time, then
    **    add terms to the loglik, etc.
    */
    istrat=0;
    indx2 =0;
    denom =0;
    for (person=0; person<nused;) {
	p = sort1[person];
	if (event[p]==0){
	    risk = exp(score[p]) * weights[p];
	    denom += risk;
	    person++;
	    }
	else {
	    time = stop[p];
	    /*
	    ** subtract out the subjects whose start time is to the right
	    */
	    for (; indx2<strata[istrat]; indx2++) {
		p = sort2[indx2];
		if (start[p] < time) break;
		risk = exp(score[p]) * weights[p];
		denom -= risk;
		}

	    /*
	    ** compute the averages over subjects with
	    **   exactly this death time (a2 & c2)
	    */
	    efron_wt =0;
	    meanwt =0;
	    deaths=0;
	    for (k=person; k<strata[istrat]; k++) {
		p = sort1[k];
		if (stop[p] < time) break;
		risk = exp(score[p]) * weights[p];
		denom += risk;

		if (event[p]==1) {
		    deaths += event[p];
		    efron_wt += risk*event[p];
		    meanwt += weights[p];
		    }
		}
	    ksave = k;

	    /*
	    ** Now add it into the loglik
	    */
	    meanwt /= deaths;
	    itemp = -1;
	    for (; person<ksave; person++) {
		p = sort1[person];
		if (event[p]==1) {
		    itemp++;
		    temp = itemp*method/deaths;
		    d2 = denom - temp*efron_wt;
		    *loglik +=  weights[p]*score[p] -meanwt *log(d2);
		    }
		}
	    }

	if (person == strata[istrat]) {
	    istrat++;
	    denom =0;
	    indx2 = person;
	    }
	}   /* end  of accumulation loop */

    /*
    ** add in the penalty terms
    */
    if (ptype==2 || ptype==3) {
	/* there are non-sparse terms */
	cox_callback(2, beta, upen, ipen, &logpen, zflag, nvar, fexpr2, rho);
	*loglik += logpen;
        }
    }

/********************************************************************/

/*
** This call is used for iteration
*/

void agfit5b(int  *maxiter, int  *nusedx, int  *nvarx, 
	       int  *strata, double *beta, double *u,
	       double *imat2,  double *jmat2, double *loglik, 
	       int  *flag,  double *eps, double *tolerch, int  *methodx, 
	       int  *nfrail, double *fbeta, double *fdiag,
               void *fexpr1, void *fexpr2, void *rho)
{
    int i,j,k, person;
    int ii;
    int     iter;
    int     nused, nvar;
    int    nf, nvar2;
    int    fgrp;
    int    halving;
    int    itemp, deaths;
    int  istrat, indx2, p, ksave;  

    double  denom, zbeta, risk;
    double  temp, temp2;
    double  newlk =0;
    double  d2, efron_wt;
    double  meanwt, time;
    double  method;
    double  **imat, **jmat;

    nused = *nusedx;
    nvar  = *nvarx;
    nf= *nfrail;
    method= *methodx;
    nvar2 = nvar + nf;
    if (nvar >0) {
        imat  = dmatrix(imat2, nvar2, nvar);
	jmat  = dmatrix(jmat2, nvar2, nvar);
        }
    else {
	imat = 0;   /*never used, but passed as dummy to chol */
	jmat = 0;
        }

    for (i=0; i<nf; i++) oldbeta[i] = fbeta[i];
    for (i=0; i<nvar; i++) oldbeta[i+nf] = beta[i];

    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=0; iter<=*maxiter; iter++) {
	newlk = 0;
	for (i=0; i<nf; i++) fdiag[i] =0;
	for (i=0; i<nvar2; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++) jmat[j][i] =0 ;
            }

        for (person=0; person<nused; person++) {
	    if (nf>0) {
		fgrp = frail[person] -1;
		zbeta = offset[person] + fbeta[fgrp];
	        }
	    else zbeta = offset[person];
	    for (i=0; i<nvar; i++)
		zbeta += beta[i]*covar[i][person];
	    zbeta = coxsafe(zbeta);
	    score[person] = zbeta;
	    zbeta = coxsafe(zbeta);
#if(0)
	    /* I believe this is unnecessary after adding coxsafe calls */
	    if (zbeta > 20 && *maxiter >1) {
		/*
		** If the above happens, then 
		**   1. There is a real chance for catastrophic cancellation
		**       in the computation of "denom", which leads to
		**       numeric failure via log(neg number) -> inf loglik
		**   2. A risk score for one person of exp(20) > 400 million
		**       is either an infinite beta, in which case any
		**       reasonable coefficient will do, or a big overreach
		**       in the Newton-Raphson step.
		** In either case, a good solution is step halving.  However,
		**   if the user asked for exactly 1 iteration, we should
		**   just return what they asked.
		** 
		** Why 20?  Most machines have about 16 digits of precision,
		**   and this preserves approx 7 digits in the subtraction
		**   when a high risk score person leaves the risk set.
		**   (Because of centering, the average risk score is about 0).
		**   Second, if eps is small and beta is infinite, we rarely
		**   get a value above 16.  So a 20 is usually a NR overshoot.
		** A data set with zbeta=54 on iter 1 led to this fix, the
		**   true final solution had max values of 4.47.    
		*/
		halving=1;
		for (i=0; i<nvar; i++)
		    beta[i] = (oldbeta[i+nf] + beta[i]) /2; 
		for (i=0; i<nf; i++)
		    fbeta[i] = (oldbeta[i] + fbeta[i])/2;
		person = -1;  /* force the loop to start over */
		}
#endif
	    }

        istrat=0;
        indx2 =0;
	denom =0;
	for (i=0; i<nvar2; i++) {
	    a[i] =0;
	    for (j=0; j<nvar; j++) {
		cmat[j][i]=0;
		}
	    }

	for (person=0; person<nused; ) {
	    p = sort1[person];
	    if (nf>0)  fgrp = frail[p] -1;
	    else       fgrp = -1;
	    if (event[p]==0){
		risk = exp(score[p]) * weights[p];
		denom += risk;
		if (fgrp >=0) a[fgrp] += risk;
		for (i=0; i<nvar; i++) {
		    a[i+nf] += risk*covar[i][p];
		    if (fgrp >=0) cmat[i][fgrp] += risk * covar[i][p];
		    for (j=0; j<=i; j++)
			cmat[i][j+nf] += risk*covar[i][p]*covar[j][p];
		    }
		person++;
		}
	    else {
		time = stop[p];
		/*
		** subtract out the subjects whose start time is to the right
		*/
		for (; indx2<strata[istrat]; indx2++) {
		    p = sort2[indx2];
		    if (start[p] < time) break;
		    risk = exp(score[p]) * weights[p];
		    denom -= risk;
		    if (nf >0) fgrp = frail[p] - 1;
		    else       fgrp = -1;
		    if (fgrp >=0) a[fgrp] -= risk;
		    for (i=0; i<nvar; i++) {
			a[i+nf] -= risk*covar[i][p];
			if (fgrp >=0) cmat[i][fgrp] -= risk* covar[i][p];
			for (j=0; j<=i; j++)
			    cmat[i][j+nf] -= risk*covar[i][p]*covar[j][p];
			}
		    }
		/*
		** compute the averages over this death time (a2 & c2)
		*/
		efron_wt =0;
		meanwt =0;
		for (i=0; i<nvar2; i++) {
		    a2[i]=0;
		    for (j=0; j<nvar; j++) {
			cmat2[j][i]= 0;
			}
		    }	
		deaths=0;
		for (k=person; k<strata[istrat]; k++) {
		    p = sort1[k];
		    if (stop[p] < time) break;
		    risk = exp(score[p]) * weights[p];
		    denom += risk;
		    if (nf>0) {
			fgrp = frail[p] -1;
			if (fgrp>=0) a[fgrp] += risk;
			}
		    else fgrp = -1;
		    for (i=0; i<nvar; i++) {
			a[i+nf] += risk*covar[i][p];
			if (fgrp>=0) cmat[i][fgrp] += risk*covar[i][p];
			for (j=0; j<=i; j++)
			    cmat[i][j+nf] += risk*covar[i][p]*covar[j][p];
			}
		    if (event[p]==1) {
			deaths += event[p];
			efron_wt += risk* weights[p];
			meanwt += weights[p];
			if (fgrp >= 0) {
			    u[fgrp] += weights[p];
			    a2[fgrp] += risk;
			    }
			for (i=0; i<nvar; i++) {
			    a2[i+nf]+= risk*covar[i][p];
			    u[i+nf] += weights[p]* covar[i][p];
			    if (fgrp >=0) cmat2[i][fgrp] += risk*covar[i][p];
			    for (j=0; j<=i; j++)
				cmat2[i][j+nf] += risk*covar[i][p]*covar[j][p];
			    }
			}
		    }	
		ksave =k;

		/* add results into u and imat */
		itemp = -1;
		meanwt /= deaths;
		for (; person<ksave; person++) {
		    p = sort1[person];
		    if (event[p]==1) {
			itemp++;
			temp = itemp*method/deaths;
			d2 = denom - temp*efron_wt;
			newlk +=  weights[p]*score[p] -meanwt *log(d2);
			for (i=0; i<nvar2; i++) {  /* by row of full matrix */
			    temp2 = (a[i] - temp*a2[i])/d2;
			    tmean[i] = temp2;
			    u[i] -= meanwt*temp2;
			    if (i<nf) fdiag[i] += temp2 * (1-temp2);
			    else {
				ii = i-nf;  /*actual row in c/j storage space*/
				for (j=0; j<=i; j++) 
				    jmat[ii][j] +=  meanwt*(
					(cmat[ii][j] - temp*cmat2[ii][j]) /d2 -
					temp2*tmean[j]);
				}
			    }
			}
		    }
		}

	    if (person == strata[istrat]) {
		istrat++;
		denom =0;
		indx2 = person;
		for (i=0; i<nvar2; i++) {
		    a[i] =0;
		    for (j=0; j<nvar; j++) {
			cmat[j][i]=0;
			}
		    }
		}
	    }   /* end  of accumulation loop */

	/*
	** Add in the penalty terms
	*/
	if (ptype==1 || ptype==3) {
	    /* there are sparse terms */
	    cox_callback(1, fbeta, upen, ipen, &logpen, zflag, nf, fexpr1,rho); 
	    if (zflag[0] ==1) {  /* force terms to zero */
		for (i=0; i<nf; i++) {
		    u[i]=0;
		    fdiag[i] =1;
		    for (j=0; j<nvar; j++) jmat[j][i]=0;
		    }
		}
	    else {
		for (i=0; i<nf; i++) {
		    u[i] += upen[i];
		    fdiag[i] += ipen[i];
		    }
		newlk += logpen;
		}
	    }

	if (ptype==2 || ptype==3) {
	    /* there are non-sparse terms */
	    cox_callback(2, beta, upen, ipen, &logpen, zflag, nvar, fexpr2, rho);
	    newlk += logpen;
	    if (pdiag==0) {
		for (i=0; i<nvar; i++) {
		    u[i+nf] += upen[i];
		    jmat[i][i+nf] += ipen[i];
		    }
		}
	    else {
		k =0;
		for (i=0; i<nvar; i++) {
		    u[i+nf] += upen[i];
		    for (j=nf; j<nvar2; j++) jmat[i][j] += ipen[k++];
		    }
		}
	    for (i=0; i<nvar; i++) {
		if (zflag[i] ==1) {
		    u[i+nf]=0;
		    for (j=0; j<i; j++) jmat[i][j+nf]=0;
		    jmat[i+nf][i] =1;
		    }
		}
	    }

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky3(jmat, nvar2, nf, fdiag, *tolerch);
	if (fabs(newlk) < *eps ||
	    (fabs(1-(*loglik/newlk))<=*eps && halving==0)) { /* all done */
	    *loglik = newlk;
	    for (i=0; i<nvar; i++) {
		for (j=0; j<nvar2; j++)  imat[i][j] = jmat[i][j];
		}
	    chinv3(jmat, nvar2, nf, fdiag);
	    for (i=nf; i<nvar2; i++) {       /*nicer output for S user */
		fdiag[i] = jmat[i-nf][i];
		jmat[i-nf][i] =1;
		imat[i-nf][i] =1;
		for (j=i+1; j<nvar2; j++) {
		    jmat[i-nf][j] = 0;
		    imat[i-nf][j] = 0;
		    }
		}

	    *maxiter = iter;
	    return;
	    }

	if (iter==*maxiter) break;  /*skip the step halving and etc */

	if (iter>0 && newlk < *loglik)   {    /*it is not converging ! */
	    halving =1;
	    for (i=0; i<nvar; i++)
		beta[i] = (oldbeta[i+nf] + beta[i]) /2; 
	    for (i=0; i<nf; i++)
		fbeta[i] = (oldbeta[i] + fbeta[i])/2;
	    }
	else {
	    halving=0;
	    *loglik = newlk;
	    chsolve3(jmat,nvar2, nf, fdiag, u);

	    j=0;
	    for (i=0; i<nvar; i++) {
		oldbeta[i+nf] = beta[i];
		beta[i] += u[i+nf];
		}
	    for (i=0; i<nf; i++) {
		oldbeta[i] = fbeta[i];
		fbeta[i] += u[i];
		}
	    }
	}   /* return for another iteration */

    *loglik = newlk;
    for (i=0; i<nvar; i++) 
	for (j=0; j<nvar2; j++) {
            imat[i][j] = jmat[i][j];
        }
    chinv3(jmat, nvar2, nf, fdiag);
    for (i=nf; i<nvar2; i++) {       /*nicer output for S user */
	fdiag[i] = jmat[i-nf][i];
	jmat[i-nf][i] =1;
  	imat[i-nf][i] =1;
  	for (j=i+1; j<nvar2; j++) {
	    jmat[i-nf][j] = 0;
  	    imat[i-nf][j] = 0;
  	    }
        }
    *flag= 1000;
    return;
    }

static double **cmatrix(double *data, int ncol, int nrow)
    {
    int i,j;
    double **pointer;
    double *temp;
 
    pointer = CALLOC(nrow, double *);
    temp =    CALLOC(nrow*ncol, double);
    if (data==0){
	for (i=0; i<nrow; i++) {
	    pointer[i] = temp;
	    temp += ncol;
            }
        }
    else {
	for (i=0; i<nrow; i++) {
	    pointer[i] = temp;
	    for (j=0; j<ncol; j++) *temp++ = *data++;
	    }
        }
    return(pointer);
	}

static void cmatrix_free(double **data) {
    FREE(*data);
    FREE(data);
    }

void agfit5c(int  *nvar) {
    /*
    ** FREE up the extra memory
    */
    FREE(zflag);
    FREE(upen);
    FREE(event);
    FREE(a);
    if (frail != NULL) FREE(frail);
    if (*nvar > 0) {
	cmatrix_free(cmat2);
	cmatrix_free(cmat);
	cmatrix_free(covar);
        }
    }
