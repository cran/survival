/* SCCS @(#)survreg3.c	1.1 02/06/99
** This is a version of survreg2.c, but for the "user written" distributions
**   It get the derivative information for all observations at once.
**
** Input
**      maxiter - max # of iterations allowed
**      n       - number of subjects
**      nvar    - number of variables in the x matrix
**      y       - matrix of start, stop, event
**      ny      - # columns of y.  =3 if there is interval censored data
**      event   - 1=exact, 0= right censored, 2=left censored, 3=interval
**      covar   - covariates for patient i
**                   Note that S sends this in column major order
**      wt      - vector of case weights (usually 1)
**      offset  - offset vector (usually 0)
**      beta    - initial values for the parameters, of length 1+nvar
**		     last element contains the scale
**      nstrat  - an indicator: 0= scale is fixed, 1= estimate scale,
**		    >1: estimate multiple scales (strata)
**      strat   - if nstrat>0, contains the strata number for each subject
**      eps     - tolerance for convergence.  Iteration continues until the
**                  relative change in the deviance is <= eps.
**      tol_chol- tolerance for Cholesky decomposition
**      dist    -  1=extreme value, 2=logistic, 3=gaussian, 4=user
**
**  Output
**      beta    - the final coef vector
**      maxiter - the number of iterations consumed
**      imat    - the information matrix (nvar+1) by (nvar+1)
**      loglik  - the final log-liklihood
**      flag    - success flag  0 =ok
**                              -1= did not converge
**
**  Work arrays
**      newbeta(nvar)- always contains the "next iteration"
**      u(nvar)      - first deriv of the loglik
**      savediag(nvar)- for rnewton
**      JJ = the approx variance matrix J'J, guarranteed non-singular
**    these are all passed as a single vector, and then broken out.
**
**  Work array "dist" is allocated here using Salloc.  This is so that
**    all of the arguments can be the same as for the usual call.
*/
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "survS.h"
#include "survproto.h"

static int    nvar, nvar2, nstrat;
static double **covar;
static int   *strat ;
static double *time2, *time1, *status;
static double *offset;
static double **imat, **JJ;
static double *u, *wt;
static double scale;
static double **funs, *z;

static double dolik(int, double *, int, void *, void *);


static int debug;
void survreg3(int   *maxiter,   int   *nx,    int   *nvarx, 
	     double *y,          int   *ny,    double *covar2, double *wtx,
	     double *offset2,    double *beta,  int   *nstratx, 
	     int   *stratax,    double *ux,    double *imatx, 
	     double *loglik,     int   *flag,  double *eps,
	     double *tol_chol,   int   *dist,  int   *ddebug,
	      void *density, void *rho
) {
    int i,j;	
    int n;
    double *newbeta,	   *savediag;
    /*double temp; -Wall unused*/
    int halving, iter;
    double newlk;


    n = *nx;
    nvar = *nvarx;
    debug = *ddebug;
    offset = offset2;
    nstrat = *nstratx;
    strat  = stratax;
    wt = wtx;
    
    covar = dmatrix(covar2, n, nvar);
 
    /*
    ** nvar = # of "real" x variables, for iteration
    ** nvar2= # of parameters
    ** nstrat= # of strata, where 0== fixed sigma
    */
    nstrat = *nstratx;
    nvar2 = nvar + nstrat;   /* number of coefficients */
    if (nstrat==0) scale = exp(beta[nvar]);

    imat = dmatrix(imatx, nvar2, nvar2);
    u = ux;
    newbeta = u+nvar2;
    savediag= newbeta + nvar2;
    JJ  = dmatrix(savediag+nvar2, nvar2, nvar2);

    if (*ny==2) {
	time1=y;
	status = y+n;
	}
    else {
	time1=y;
	time2 = time1 + n;
	status = time2 +n;
	}

    /* count up the number of interval censored obs 
    **  and allocate memory for the callback arrarys
    */
    j =0;    for (i=0; i<n; i++)  
	if (status[i]==3) j++;
    j = j+n;
    funs  = dmatrix((double *)ALLOC(j*5, sizeof(double)), j, 5);
    z     = (double *)ALLOC(j, sizeof(double));

    /*
    ** do the initial iteration step
    */
    *loglik = dolik(n, beta, 0,density,rho); 
    if (debug >0) {
	fprintf(stderr, "nvar=%d, nvar2=%d, nstrat=%d\n", nvar, nvar2, nstrat);
        fprintf(stderr, "iter=0, loglik=%f\n", loglik[0]);
	}

    *flag= cholesky2(imat, nvar2, *tol_chol);
    if (*flag < 0) {
	i = cholesky2(JJ, nvar2, *tol_chol);
	chsolve2(JJ, nvar2, u);
	if (debug>0) fprintf(stderr, " Alternate step, flag=%d\n", i);
	}
    else chsolve2(imat,nvar2,u);        /* a replaced by  a *inverse(i) */
    if (debug>0) {
	fprintf(stderr, " flag=%d, Increment:", *flag);
	for (i=0; i<nvar2; i++) fprintf(stderr, " %f", u[i]);
	fprintf(stderr, "\n");
	}
    if (debug >2) {
	fprintf(stderr, "Imat after inverse\n");
	for (i=0; i<nvar2; i++) {
	    for (j=0; j<nvar2; j++) fprintf(stderr,"  %f", imat[i][j]);
	    fprintf(stderr, "\n");
	    }
	}
    
    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar2; i++) {
	newbeta[i] = beta[i] + u[i];
	}
    if (*maxiter==0) {
	chinv2(imat,nvar2);
	for (i=1; i<nvar2; i++)
	    for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	return;   /* and we leave the old beta in peace */
	}


    /*
    ** here is the main loop
    */
    halving =0 ;             /* >0 when in the midst of "step halving" */
    newlk = dolik(n, newbeta, 0,density,rho); 

    /* put in a call to simplex if in trouble */

    for (iter=1; iter<=*maxiter; iter++) {
	if (debug>0) fprintf(stderr, "---\niter=%d, loglik=%f\n\n", iter, 
			     newlk);

	/* 
	**   Am I done?  Check for convergence, then update betas
	*/
	if (fabs(1-(*loglik/newlk))<=*eps ) { /* all done */
	    *loglik = newlk;
	    *flag = cholesky2(imat, nvar2, *tol_chol);
	    if (debug==0) {
		chinv2(imat, nvar2);     /* invert the information matrix */
		for (i=1; i<nvar2; i++)
		    for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
		}
	    for (i=0; i<nvar2; i++)
		beta[i] = newbeta[i];
	    if (halving==1) *flag= 1000; /*didn't converge after all */
	    *maxiter = iter;
	    return;
	    }

	if (newlk < *loglik)   {    /*it is not converging ! */
	    for (j=0; j<5 && newlk < *loglik; j++) {
		halving++;
		for (i=0; i<nvar2; i++)
		    newbeta[i] = (newbeta[i] + beta[i]) /2; 
		/*
		** Special code for sigmas.  Often, they are the part
		**  that gets this routine in trouble.  The prior NR step
		**  may have decreased one of them by a factor of >10, in which
		**  case step halving isn't quite enough.  Make sure the new
		**  try differs from the last good one by no more than 1/3
		**  approx log(3) = 1.1
		**  Step halving isn't enough of a "back away" when a
		**  log(sigma) goes from 0.5 to -3, or has become singular.
		*/
		if (halving==1) {  /* only the first time */
		    for (i=0; i<nstrat; i++) {
			if ((beta[nvar+i]-newbeta[nvar+i])> 1.1)
			    newbeta[nvar+i] = beta[nvar+i] - 1.1;  
			}
		    }
		newlk = dolik(n, newbeta, 1,density,rho);
		}
	    if (debug>0) {
		fprintf(stderr,"   Step half -- %d steps, newlik=%f\n", 
			halving, newlk);
		fflush(stderr);
		}
	    }

	else {    /* take a standard NR step */
	    halving=0;
	    *loglik = newlk;
	    *flag = cholesky2(imat, nvar2, *tol_chol);
	    if (debug >2) {
		fprintf(stderr, "Imat after inverse\n");
		for (i=0; i<nvar2; i++) {
		    for (j=0; j<nvar2; j++) fprintf(stderr,"  %f", imat[i][j]);
		    fprintf(stderr, "\n");
		    }
		}
	    if (*flag < 0) {
		i = cholesky2(JJ, nvar2, *tol_chol);
		chsolve2(JJ, nvar2, u);
		if (debug>0) fprintf(stderr, " Alternate step, flag=%d\n", i);
		}
	    else chsolve2(imat,nvar2,u);
	    if (debug>1) {
		fprintf(stderr, " flag=%d, Increment:", *flag);
		for (i=0; i<nvar2; i++) fprintf(stderr, " %f", u[i]);
		fprintf(stderr, "\n");
		}
	    for (i=0; i<nvar2; i++) {
		beta[i] = newbeta[i];
		newbeta[i] = newbeta[i] +  u[i];
		}
	    }

	newlk = dolik(n, newbeta, 0,density, rho);
	}   /* return for another iteration */

    *loglik = newlk;
    if (debug==0) {
	cholesky2(imat, nvar2, *tol_chol);
	chinv2(imat, nvar2); 
	for (i=1; i<nvar2; i++) {
	    for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	    }
	}
    for (i=0; i<nvar2; i++)
	beta[i] = newbeta[i];
    *flag= 1000;
    return;
    }
    

/*
** This routine calculates the loglik, but more important
**   it sets up the arrays for the main computation
*/
static double dolik(int n, double *beta, int whichcase, void *density, void*rho) {
    int person, i,j,k;
    int strata;
    double  eta,
	    sigma;
    int     icount;  /* running count of # of interval censored */
    double  loglik,
	    temp;
    double  temp1, temp2;
    double  sz, zz, zu;
    double  sig2;
    double g=0, dg=0, ddg=0, dsig=0, ddsig=0, dsg=0;/*-Wall*/


    for (i=0; i<nvar2; i++) {
	u[i] =0;
	for (j=0; j<nvar2; j++) {
	    imat[i][j] =0 ;
	    JJ[i][j] =0;
	    }
	}

    strata =0;
    if (nstrat ==0) sigma = scale;   /* fixed scale */
    else            sigma = exp(beta[nvar]);
    sig2  = 1/(sigma*sigma);
    loglik =0;

    /*
    ** First, get the array of distribution values
    ** We get them all at once to minimize the S-callback overhead
    */
    icount =n;
    for (person=0; person<n; person++) {
	if (nstrat>1) {
	    strata= strat[person] -1; /*S likes to start counting at 1 */
	    sigma = exp(beta[strata+nvar]);
	    }
	eta =0;
	for (i=0; i<nvar; i++) eta += beta[i] * covar[i][person];
	eta += offset[person];
	z[person] = (time1[person] - eta)/sigma;  
	if (status[person]==3) {
	    z[icount] = (time2[person] - eta)/sigma;
	    icount++;
	    }
	}

    surv_callback(z, funs[0], n, density, rho);  /* treat them both as a vector */

    /*
    ** calculate the first and second derivative wrt eta,
    **   then the derivatives of the loglik (u, imat, JJ)
    */
    icount =n;
    for (person=0; person<n; person++) {
	if (nstrat>1) {
	    strata= strat[person] -1; /*S likes to start counting at 1 */
	    sigma = exp(beta[strata+nvar]);
	    sig2  = 1/(sigma*sigma);
	    }

	zz = z[person];
	sz = zz * sigma;
	j = status[person];       /*convert to integer */
	switch(j) {
	    case 1:                             /* exact */
		if (funs[2][person] <=0) {
		    /* off the probability scale -- avoid log(0), and set the
		    **  derivatives to gaussian limits (almost any deriv will
		    **  do, since the function value triggers step-halving).
		    */
		    g = -FLT_MAX;
		    dg = -zz/sigma;
		    ddg = -1/sigma;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(funs[2][person])  - log(sigma);
		    temp1 = funs[3][person]/sigma;
		    temp2 = funs[4][person]*sig2;
		    dg = -temp1;
		    dsig= -(sz*temp1 +1);
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(1- sz*temp1);
		    ddsig = sz*sz*temp2 + sz*temp1*(1- sz*temp1);
		    }
		break;
	    case 0:                             /* right censored */
		if (funs[1][person] <=0) {
		    g = -FLT_MAX;
		    dg = zz/sigma;
		    ddg =0;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(funs[1][person]);
		    temp1 = -funs[2][person]/(funs[1][person]*sigma);
		    temp2 = -funs[3][person]*funs[2][person]*sig2/
			               funs[1][person];
		    dg = -temp1;
		    dsig= -sz * temp1;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(1+dsig);
		    ddsig = sz*sz*temp2 - dsig*(1+dsig);
		    }
		break;
	    case 2:                             /* left censored */
		if (funs[2][person] <=0) {
		    /* off the probability scale -- avoid log(0) */
		    g = -FLT_MAX;
		    dg = -zz/sigma;
		    dsig =0; ddsig=0; dsg=0;
		    ddg =0;
		    }
		else {
		    g = log(funs[0][person]);
		    temp1 = funs[2][person]/(funs[0][person]*sigma);
		    temp2 = funs[3][person]*funs[2][person]*sig2/
			                      funs[0][person];
		    dg= -temp1;
		    dsig= -sz * temp1;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(1+dsig);
		    ddsig = sz*sz*temp2 - dsig*(1+dsig);
		    }
		break;
	    case 3:                             /* interval censored */
                zu = z[icount];
                /*stop roundoff in tails*/
		if (zz>0)  temp = funs[1][person] - funs[1][icount]; 
		else       temp = funs[0][icount] - funs[0][person];
		if (temp <=0) {
		    /* off the probability scale -- avoid log(0) */
		    g = -FLT_MAX;
		    dg = 1; 
		    ddg =0;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(temp);
		    dg  = -(funs[2][icount] -funs[2][person])/(temp*sigma);
		    ddg = (funs[3][icount] -funs[3][person])*sig2/temp - dg*dg;
		    dsig = (zz*funs[2][person] - zu*funs[2][icount])/temp;
		    ddsig= (zu*zu*funs[3][icount] - zz*zz*funs[3][person])
			          /temp - dsig*(1+dsig);
		    dsg = (zu*funs[3][icount] - zz*funs[3][person])/
			       (temp*sigma)  - dg *(1+dsig);
		    }
		icount++;
		break;
	    }
	loglik += g * wt[person];

	/*
	** Now the derivs wrt loglik
	*/
	if (whichcase==1) continue;     /*only needed the loglik */
	for (i=0; i<nvar; i++) {
	    temp = wt[person] * covar[i][person];
	    u[i] += temp * dg;
	    for (j=0; j<=i; j++) {
		imat[j][i] -= temp *covar[j][person] *ddg;
		JJ[j][i] += temp * covar[j][person] * dg * dg;
		}
	    }

	if (nstrat!=0) {   /* need derivative wrt log sigma */
	    k = strata+nvar;
	    u[k] += wt[person]* dsig;
	    for (i=0; i<nvar; i++) {
		imat[i][k] -= dsg * covar[i][person] * wt[person];
		JJ[i][k]   += dsig* covar[i][person] *dg * wt[person];
		}
	    imat[k][k] -=  ddsig * wt[person];
	    JJ[k][k] += dsig*dsig * wt[person];
	    }
	}

    if (debug >0) {
	fprintf(stderr, "coef" );
	if (nvar2==1) j=2; else j=nvar2;
	for (i=0; i<j; i++) fprintf(stderr," %f", beta[i]);
	if (whichcase==0) {
	    fprintf(stderr, "\nU   ");
	    for (i=0; i<nvar2; i++) fprintf(stderr," %f", u[i]);
	    }
	fprintf(stderr, "\n");
	}
    if (debug >1) {
	fprintf(stderr, "Imat\n");
	for (i=0; i<nvar2; i++) {
	    for (j=0; j<nvar2; j++) fprintf(stderr,"  %f", imat[i][j]);
	    fprintf(stderr, "\n");
	    }
	}
    return(loglik);
    }


