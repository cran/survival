/* SCCS @(#)survreg2.c	1.5 02/06/99
** Fit one of several censored data distributions
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
**      dist    -  1=extreme value, 2=logistic, 3=gaussian, 4=cauchy
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
**
**  For the calculations involving sigma and an interval censored datum, I
**    can't calculate the sigma derivatives simply from the eta derivs.
**
**  To add a new distribution:
**              add a new "static void" declaration
**              add it to the "switch(*dist)" list, (2 places)
**              add the new subroutine to the bottom of the code, see
**                      logist_d as an example
*/
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "survS.h"
#include "survproto.h"

#define  SPI    2.506628274631001     /* sqrt(2*pi) */
#define  ROOT_2 1.414213562373095

static void exvalue_d(double z, double ans[4], int j);
static void logistic_d(double z, double ans[4], int j);
static void gauss_d(double z, double ans[4], int j);
static void cauchy_d(double z, double ans[4], int j);
static void (*sreg_gg)();
static double dolik(int n, double *beta, int whichcase);

static int    nvar, nvar2, nstrat;
static double **covar;
static int   *strat ;
static double *time2, *time1, *status;
static double *offset;
static double **imat, **JJ;
static double *u, *wt;
static double scale;

static int debug;
void survreg2(int   *maxiter,   int   *nx,    int   *nvarx, 
	     double *y,          int   *ny,    double *covar2, double *wtx,
	     double *offset2,    double *beta,  int   *nstratx, 
	     int   *stratax,    double *ux,    double *imatx, 
	     double *loglik,     int   *flag,  double *eps,
	     double *tol_chol,   int   *dist,  int   *ddebug,
	      /* need same prototype as survreg3 */
	      void *placeholder1, void *placeholder2)
{
    int i,j;	
    int n;
    double *newbeta,
	   *savediag;
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

    switch(*dist) {
	case 1: sreg_gg = exvalue_d;  break;
	case 2: sreg_gg = logistic_d; break;
	case 3: sreg_gg = gauss_d;    break;
	case 4: sreg_gg = cauchy_d;   break;
	}

    /*
    ** do the initial iteration step
    */
    *loglik = dolik(n, beta, 0); 
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
    newlk = dolik(n, newbeta, 0); 

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
		newlk = dolik(n, newbeta, 1);
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

	newlk = dolik(n, newbeta, 0);
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
static double dolik(int n, double *beta, int whichcase) {
    int person, i,j,k;
    int strata;
    double  eta,
	    sigma;
    double  z, zu,
	    loglik,
	    temp, temp2;
    double  sz;
    double  sig2;
    static double  funs[4], ufun[4];
    double g=0, dg=0, ddg=0, dsig=0, ddsig=0, dsg=0; /*-Wall*/


    for (i=0; i<nvar2; i++) {
	u[i] =0;
	for (j=0; j<nvar2; j++) {
	    imat[i][j] =0 ;
	    JJ[i][j] =0;
	    }
	}

    /*
    ** calculate the first and second derivative wrt eta,
    **   then the derivatives of the loglik (u, imat, JJ)
    */
    strata =0;
    if (nstrat ==0) sigma = scale;   /* fixed scale */
    else            sigma = exp(beta[nvar]);
    sig2  = 1/(sigma*sigma);
    loglik =0;
    for (person=0; person<n; person++) {
	if (nstrat>1) {
	    strata= strat[person] -1; /*S likes to start counting at 1 */
	    sigma = exp(beta[strata+nvar]);
	    sig2  = 1/(sigma*sigma);
	    }
	eta =0;
	for (i=0; i<nvar; i++) eta += beta[i] * covar[i][person];
	eta += offset[person];
	sz = (time1[person] - eta);  /*   sigma * z  */
	z = sz /sigma;

	j = status[person];       /*convert to integer */
	switch(j) {
	    case 1:                             /* exact */
		(*sreg_gg)(z, funs,1);
		if (funs[1] <=0) {
		    /* off the probability scale -- avoid log(0), and set the
		    **  derivatives to gaussian limits (almost any deriv will
		    **  do, since the function value triggers step-halving).
		    */
		    g = -FLT_MAX;
		    dg = -z/sigma;
		    ddg = -1/sigma;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(funs[1])  - log(sigma);
		    temp = funs[2]/sigma;
		    temp2= funs[3]*sig2;
		    dg = -temp;
		    dsig= -temp*sz;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(dsig +1);
		    ddsig = sz*sz* temp2 - dsig*(1+dsig);
		    dsig -= 1;
		    }
		break;
	    case 0:                             /* right censored */
		(*sreg_gg)(z, funs,2);
		if (funs[1] <=0) {
		    g = -FLT_MAX;
		    dg = z/sigma;
		    ddg =0;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(funs[1]);
		    temp = -funs[2]/(funs[1]*sigma);
		    temp2= -funs[3]*sig2/funs[1];
		    dg = -temp;
		    dsig= -temp*sz;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(dsig +1);
		    ddsig = sz*sz* temp2 - dsig*(1+dsig);
		    }
		break;
	    case 2:                             /* left censored */
		(*sreg_gg)(z, funs,2);
		if (funs[0] <=0) {
		    /* off the probability scale -- avoid log(0) */
		    g = -FLT_MAX;
		    dg = -z/sigma;
		    dsig =0; ddsig=0; dsg=0;
		    ddg =0;
		    }
		else {
		    g = log(funs[0]);
		    temp = funs[2]/(funs[0]*sigma);
		    temp2= funs[3]*sig2/funs[0];
		    dg = -temp;
		    dsig= -temp*sz;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(dsig +1);
		    ddsig = sz*sz* temp2 - dsig*(1+dsig);
		    }
		break;
	    case 3:                             /* interval censored */
		zu = (time2[person] - eta)/sigma;  /*upper endpoint */
		(*sreg_gg)(z, funs, 2);
		(*sreg_gg)(zu,ufun ,2);
		if (z>0)  temp = funs[1] - ufun[1]; /*stop roundoff in tails*/
		else      temp = ufun[0] - funs[0];
		if (temp <=0) {
		    /* off the probability scale -- avoid log(0) */
		    g = -FLT_MAX;
		    dg = 1; 
		    ddg =0;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(temp);
		    dg  = -(ufun[2] - funs[2])/(temp*sigma);
		    ddg = (ufun[3] - funs[3])*sig2/temp - dg*dg;
		    dsig = (z*funs[2] - zu*ufun[2])/temp;
		    ddsig= ((zu*zu*ufun[3] - z*z*funs[3])/temp) -
			                dsig*(1+dsig);
		    dsg = ((zu*ufun[3] - z*funs[3])/ (temp*sigma)) -
				      dg * (dsig +1);
		    }
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


/*
**  Case      ans[0]    ans[1]       ans[2]     ans[3]
**   1                    f          f'/f        f''/ f
**   2          F        1-F         f           f'
**
**  We do both F and 1-F to avoid the error in (1-F) for F near 1
*/

static void logistic_d(double z, double ans[4], int j)
    {
    double w, temp;
    int    sign, ii;

    /*
    ** The symmetry of the logistic allows me to be careful, and never take
    **  exp(large number).  This routine should be very accurate.
    */
    if (z>0)  {
	w = exp(-z);
	sign = -1;
	ii=0;
	}
    else {
	w = exp(z);
	sign = 1;
	ii=1;
	}
    temp = 1+w;
    switch(j) {
	case 1:  ans[1] = w/(temp*temp);
		 ans[2] = sign*(1-w)/temp;
		 ans[3] = (w*w -4*w +1)/(temp*temp);
		 break;
	case 2:  ans[1-ii] = w/temp;
		 ans[ii]   = 1/temp;
		 ans[2] = w/(temp*temp);
		 ans[3] = sign*ans[2]*(1-w)/temp;
		 break;
	}
    }

static void gauss_d(double z, double ans[4], int j)
    {
    double f;

    f = exp(-z*z/2) /SPI;
    switch(j) {
	case 1: ans[1] =f;
		ans[2] = -z;
		ans[3] = z*z -1;
		break;
	case 2: if (z>0) {
	            ans[0] = pnorm5(z,0,1,1,0); /*(1 + erf(z/ROOT_2))/2;*/
		    ans[1] = pnorm5(z,0,1,0,0); /* erfc(z/ROOT_2) /2;*/
		    }
		else {
		    ans[1] = pnorm5(-z,0,1,1,0); /*(1 + erf(-z/ROOT_2))/2;*/
		    ans[0] = pnorm5(-z,0,1,0,0); /*  erfc(-z/ROOT_2) /2;*/
		    }
		ans[2] = f;
		ans[3] = -z*f;
		break;
	}
    }

/*
** In the Gaussian and logistic cases, I could avoid numeric disaster by only
**   evaluating exp(x) for x<0.  By symmetry, I could get what I need for
**   x >0.  The extreme value dist is asymmetric, and I don't yet see the
**   numeric tricks that I need.
** Probobly, a Taylor series will need to be used for large z.
*/

static void exvalue_d(double z, double ans[4], int j)
    {
    double temp;
    double w;
    w = exp(z);
    if (w > FLT_MAX) w=FLT_MAX;  /* stop infinite answers */
    temp = exp(-w);
    switch(j) {
	case 1:  ans[1] = w*temp;
		 ans[2] = 1-w;
		 ans[3] = w*(w-3) +1;
		 break;
	case 2:  ans[0] = 1-temp;
		 ans[1] = temp;
		 ans[2] = w*temp;
		 ans[3] = w*temp*(1-w);
		 break;
	}
    }

static void cauchy_d(double z, double ans[4], int j)
    {
    double temp;

    temp = 1/(1 + z*z);
    switch(j) {
	case 1:  ans[1] = temp/PI;
		 ans[2] = -2*z*temp;
		 ans[3] = (6*z*z -2) * temp * temp;
		 break;
	case 2:  ans[0] = 0.5 + atan(z)/PI;
		 ans[1] = 1 - ans[0];
		 ans[2] = temp/PI;
		 ans[3] = -2*z*temp*temp/PI;
		 break;
	}
    }

