/* SCCS @(#)survreg5.c	1.1 02/06/99
** The variant of survreg4 for user-written distributions, penalized models
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
**      wt      - case weight
**      offset  - offset vector (usually 0)
**      beta    - initial values for the parameters: frailties, then
**		    ordinary coefs, then sigmas
**		     an extra element contains the scale when nstrat=0
**      nstrat  - an indicator: 0= scale is fixed, 1= estimate scale,
**		    >1: estimate multiple scales (strata)
**      strat   - if nstrat>0, contains the strata number for each subject
**      eps     - tolerance for convergence.  Iteration continues until the
**                  relative change in the deviance is <= eps.
**      tol_chol- tolerance for Cholesky decomposition
**      dist    -  1=extreme value, 2=logistic, 3=gaussian, 4=cauchy
**      ptype        : 1 or 3 -- there is a sparse term
**                   : 2 or 3 -- there is a non-sparse term in the model
**      nfrail       : number of frailty groups (sparse terms), 0 if there are
**                       none
**      frail        : a vector containing the frailty groups
**      fbeta        : initial frailty estimates
**      pdiag        : if 0, then for the non-sparse terms only the diagonal
**                       of the variance matrix is penalized, otherwise the
**                       full matrix is used.
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

static double **funs, *z;
static double dolik(int n, double *beta, int whichcase, void *fexpr1, void *fexpr2, void *density, void *rho);

static int    nvar0, nvar, nvar2, nstrat;
static double **covar, *wt;
static int   *strat , *frail;
static double *time2, *time1, *status;
static double *offset;
static double **imat, **JJ, **jmat;
static double *u;
static int    nf, ptype, pdiag;
static double *ipen, *upen, logpen;
static int   *zflag;
static double *fdiag, *jdiag;
static double scale;

static int debug;

void survreg5(int   *maxiter,   int   *nx,       int   *nvarx, 
	      double *y,         int   *ny,       double *covar2, 
	      double *wt2,       double *offset2,  double *beta,  
	      int   *nstratx,    int   *stratax,  double *ux,    
	      double *imatx,      double *jmatx,
	      double *loglik,     int   *flag,     double *eps,
	      double *tol_chol,   int   *dist,     int   *ddebug,
	      int *ptype2,  	 int   *pdiag2,
	      int *nfrail2,      int   *frail2,   double *fdiag2,
	      /* R callback data */
	      void *fexpr1, void *fexpr2, void *density, void *rho
    )  {

    int i,j;	
    int n;
    double *newbeta;
/*    double *savediag;
      double temp; -Wall unused*/
    int halving, iter;
    double newlk;
    int lastchance=0;

    n = *nx;
    nvar = *nvarx;
    debug = *ddebug;
    offset = offset2;
    nstrat = *nstratx;
    strat  = stratax;
    ptype  = *ptype2;
    pdiag  = *pdiag2;
    nf     = *nfrail2;
    frail  = frail2;
    fdiag  = fdiag2;
    wt     = wt2;
    /*
    ** nvar0 = # of "real" x variables, for iteration
    ** nvar  = # of non-frailty vars = nvar + #sigmas
    ** nvar2= nvar + nf = other dim of u and jmat and JJ
    ** nstrat= # of strata, where 0== fixed sigma
    */
    nstrat = *nstratx;
    nvar0 = nvar;
    nvar  = nvar + nstrat;
    nvar2 = nvar  + nf;
    if (nstrat==0) scale = exp(beta[nvar]); 
    
    covar = dmatrix(covar2, n, nvar0);
    if (nvar >0) {
        imat  = dmatrix(imatx, nvar2, nvar);
	jmat  = dmatrix(jmatx, nvar2, nvar);
        }
    else {
	imat = 0;   /*never used, but passed as dummy to chol */
	jmat = 0;
        }

    u = ux;
    newbeta = u+nvar2;
    jdiag= newbeta + nvar2;
    JJ  = dmatrix(jdiag+nvar2, nvar2, nvar);

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

    /* scratch space for penalty 
    **    upen needs to be max(nvar, nfrail), 
    **    ipen max(nfrail, nvar(if pdiag=0) or nvar^2 )
    */
    if (nf > nvar) i=nf; else i=nvar;
    if (nf > nvar*nvar) j=nf; else j=nvar*nvar;
    if (pdiag==0)  upen = Calloc(2*i, double);
    else           upen = Calloc(i+j, double);
    ipen = upen + i;
    if (ptype>1)  zflag = Calloc(nvar, int);
    else          zflag = Calloc(2, int);

    if (debug>0) {
	fprintf(stderr, "\n----------Enter survreg4-----------\n");
	fprintf(stderr, "nvar=%d, nvar2=%d, nstrat=%d\n", nvar, nvar2, nstrat);
	if (nstrat==0) fprintf(stderr, "  log(scale)=%f\n", log(scale));
	}

    /*
    ** do the initial iteration step
    */
    *loglik = dolik(n, beta, 0, fexpr1,fexpr2,density,rho); 
    if (debug >0) {
        fprintf(stderr, "iter=0, loglik=%f\n", loglik[0]);
	}

    *flag = cholesky3(jmat, nvar2, nf, fdiag, *tol_chol);
    if (*flag < 0) {
	i = cholesky3(JJ, nvar2, nf, jdiag, *tol_chol);
	chsolve3(JJ, nvar2, nf, jdiag, u);
	if (debug>0) fprintf(stderr, " Alternate step, flag=%d\n", i);
	}
    else chsolve3(jmat,nvar2, nf, fdiag, u);

    if (debug>0) {
	fprintf(stderr, " flag=%d, Increment:", *flag);
	for (i=0; i<nvar2; i++) fprintf(stderr, " %f", u[i]);
	fprintf(stderr, "\n");
	}
    if (debug >2) {
	fprintf(stderr, "Imat after inverse\n");
	for (i=0; i<nvar2; i++) {
	    for (j=0; j<nvar2; j++) fprintf(stderr,"  %f", jmat[i][j]);
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
	chinv3(jmat, nvar2, nf, fdiag);
	for (i=0; i<nvar; i++) {
	    for (j=0; j<nvar2; j++)  imat[i][j] = jmat[i][j];
	    }
	for (i=nf; i<nvar2; i++) {       /*nicer output for S user */
	    fdiag[i] = jmat[i-nf][i];
	    jmat[i-nf][i] =1;
	    imat[i-nf][i] =1;
	    for (j=i+1; j<nvar2; j++) {
		jmat[i-nf][j] = 0;
		imat[i-nf][j] = 0;
		}
	    }
	return;   /* and we leave the old beta in peace */
	}


    /*
    ** here is the main loop
    */
    halving =0 ;             /* >0 when in the midst of "step halving" */
    newlk = dolik(n, newbeta, 0, fexpr1,fexpr2,density,rho); 

    /* some day put in a call to simplex if in trouble */
    for (iter=1; iter<=*maxiter; iter++) {
	if (debug>0) fprintf(stderr, "iter=%d, loglik=%f\n\n", iter, newlk);

	/* 
	**   Am I done?  Check for convergence, then update betas
	*/
	if (fabs(1-(*loglik/newlk))<=*eps ) { /* all done */
	    *loglik = newlk;
	    *flag = cholesky3(jmat, nvar2, nf, fdiag, *tol_chol);
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
	    for (i=0; i<nvar2; i++) beta[i] = newbeta[i];
	    if (halving >0) *flag += 1000; /*didn't converge after all */
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

		newlk = dolik(n, newbeta, 1, fexpr1,fexpr2,density,rho);
		}
	    if (debug>0) {
		fprintf(stderr,"   Step half -- %d steps, newlik=%f\n", 
			halving, newlk);
		fflush(stderr);
		}
	    }
	else {    /* take a standard NR step */
	    halving=0; lastchance=0;
	    *loglik = newlk;
	    *flag = cholesky3(jmat, nvar2, nf, fdiag, *tol_chol);
	    if (debug >2) {
		fprintf(stderr, "Imat after inverse\n");
		for (i=0; i<nvar2; i++) {
		    for (j=0; j<nvar2; j++) fprintf(stderr,"  %f", jmat[i][j]);
		    fprintf(stderr, "\n");
		    }
		}
	    if (*flag < 0) {
		i = cholesky3(JJ, nvar2, nf, jdiag, *tol_chol);
		chsolve3(JJ, nvar2, nf, jdiag, u);
		if (debug>0) fprintf(stderr, " Alternate step, flag=%d\n", i);
		lastchance =1;
		}
	    else chsolve3(jmat,nvar2, nf, fdiag, u);
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

	if (halving > 9 && newlk < *loglik) {
	    if (lastchance ==0){
		/* Iteration is in trouble: the last NR step is no
		**  good even after 10 step halvings.
		** Try a Fisher step instead  -- 
		*/
		for (i=0; i<nvar2; i++) newbeta[i] = beta[i]; /* go back */
		newlk = dolik(n, newbeta, 0, fexpr1,fexpr2,density,rho);   /*recreate u, fdiag, etc */
		i = cholesky3(JJ, nvar2, nf, jdiag, *tol_chol);
		chsolve3(JJ, nvar2, nf, jdiag, u);
		if (debug>0)
		    fprintf(stderr, " Alternate step b, flag=%d loglik=%f\n",
			    i, newlk);

		for (i=0; i<nvar2; i++) 
		    newbeta[i] = beta[i] +  u[i];
		halving =0;
		lastchance=1;
		}
	    else {
		/* REAL trouble: a Fisher step didn't work
		**   I hope this never happens...
		*/
		*flag=1000;
		for (i=0; i<nvar2; i++) newbeta[i] = beta[i]; /* go back */
		newlk = dolik(n, newbeta, 0, fexpr1,fexpr2,density,rho);   /*recreate u, fdiag, etc */
		break;
		}
	    }
	newlk = dolik(n, newbeta, 0, fexpr1,fexpr2,density,rho);
	}   /* return for another iteration */

    *loglik = newlk;
    *flag = cholesky3(jmat, nvar2, nf, fdiag, *tol_chol);
    chinv3(jmat, nvar2, nf, fdiag);
    for (i=0; i<nvar; i++) {
	for (j=0; j<nvar2; j++)  imat[i][j] = jmat[i][j];
	}
    for (i=nf; i<nvar2; i++) {       /*nicer output for S user */
	fdiag[i] = jmat[i-nf][i];
	jmat[i-nf][i] =1;
	imat[i-nf][i] =1;
	for (j=i+1; j<nvar2; j++) {
	    jmat[i-nf][j] = 0;
	    imat[i-nf][j] = 0;
	    }
	}

    for (i=0; i<nvar2; i++)
	beta[i] = newbeta[i];
    *flag += 1000;
    return;
    }
    

/*
** This routine calculates the loglik, but more important
**   it sets up the arrays for the main computation
*/
static double dolik(int n, double *beta, int whichcase, 
		    void *fexpr1, void *fexpr2, void *density, void *rho) {
    int person, i,j,k;
    int strata;
    double  eta,
	    sigma;
    double  zz, zu,
	    loglik,
	    temp;
    double  temp1, temp2;
    double  sz;
    double  sig2;
    double  w;
    double g=0, dg=0, ddg=0, dsig=0, ddsig=0, dsg=0;
    int fgrp=0;
    int icount;

    loglik=0;
    for (i=0; i<nf; i++) fdiag[i] =0;
    for (i=0; i<nvar2; i++) {
	u[i] =0;
	for (j=0; j<nvar; j++) {
	    jmat[j][i] =0 ;
	    JJ[j][i] =0;
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
	if (nf>0) {
	    fgrp = frail[person] -1;
	    eta  = offset[person] + beta[fgrp];
	    }
	else eta = offset[person];
	for (i=0; i<nvar; i++) eta += beta[i] * covar[i][person];
	z[person] = (time1[person] - eta)/sigma;  
	if (status[person]==3) {
	    z[icount] = (time2[person] - eta)/sigma;
	    icount++;
	    }
	}

    surv_callback(z, funs[0], n,density, rho);  /* treat them both as a vector */

    /*
    ** calculate the first and second derivative wrt eta,
    **   then the derivatives of the loglik (u, imat, JJ)
    */
    icount =n;
    for (person=0; person<n; person++) {
	if (nstrat>1) {
	    strata= strat[person] -1; /*S likes to start counting at 1 */
	    sigma = exp(beta[strata+nvar0+nf]);
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
	if (debug>3) {
	    fprintf(stderr," z=%f g=%f, dg=%f, wt=%f\n", z[icount], g, dg, wt[person]);
	    fflush(stderr);
	    }
     
	/*
	** Now the derivs wrt loglik
	**   Sparse frailties are a covariate of "1"
	*/
	if (whichcase==1) continue;     /*only needed the loglik */
	w = wt[person];
	if (nf>0) {
	    u[fgrp] += dg * w;
	    fdiag[fgrp] -= ddg * w;
	    jdiag[fgrp] += dg*dg *w;
	    }
	for (i=0; i<nvar0; i++) {
	    temp = dg * covar[i][person] *w;
	    u[i+nf] += temp;
	    for (j=0; j<=i; j++) {
		jmat[i][j+nf] -= covar[i][person] *covar[j][person] *ddg *w;
		JJ[i][j+nf]   += temp * covar[j][person] * dg;
		}
	    if (nf>0) {
		jmat[i][fgrp] -= covar[i][person] * ddg * w;
		JJ  [i][fgrp] += temp * dg;
		}
	    }

	if (nstrat!=0) {   /* need derivative wrt log sigma */
	    k = strata+nvar0;
	    u[k+nf] += w* dsig;
	    for (i=0; i<nvar0; i++) {
 		jmat[k][i+nf] -= dsg * covar[i][person] * w;
		JJ[k][i+nf]   += dsig* covar[i][person] *dg * w;
		}
	    jmat[k][k+nf] -=  ddsig * w;
	    JJ[k][k+nf] += dsig*dsig * w;

	    if (nf>0) {
		jmat[k][fgrp] -= dsg * w;
 		JJ  [k][fgrp] += dsig *dg *w;
 		}
	    }
	}
    if (debug >1) {
	fprintf(stderr," U (no penalty)");
	for (i=0; i<nvar2; i++) fprintf(stderr,"  %f", u[i]);
	fprintf(stderr, "\n");
	fflush(stderr);
	}

    /*
    ** Add in the penalty terms
    */
    if (ptype==1 || ptype==3) {
	/* there are sparse terms */
	cox_callback(1, beta, upen, ipen, &logpen, zflag,nf, fexpr1, rho); 
	if (zflag[0] ==1) {  /* force terms to zero */
	    for (i=0; i<nf; i++) {
		u[i]=0;
		fdiag[i] =1;
		jdiag[i] =1;
		for (j=0; j<nvar; j++) {
		    JJ[j][i] =0;
		    jmat[j][i]=0;
		    }
		}
	    }
	else {
	    for (i=0; i<nf; i++) {
		u[i] += upen[i];
		fdiag[i] += ipen[i];
		jdiag[i] += ipen[i];
		}
	    loglik += logpen;
	    }
	}

    if (ptype==2 || ptype==3) {
	/* there are non-sparse terms */
	cox_callback(2, beta+nf, upen, ipen, &logpen, zflag, nvar, fexpr2, rho);
	loglik += logpen;
	if (pdiag==0) {
	    for (i=0; i<nvar; i++) {
		u[i+nf] += upen[i];
		jmat[i][i+nf] += ipen[i];
		JJ[i][i+nf] += ipen[i];
		}
	    }	
	else {
	    k =0;
	    for (i=0; i<nvar; i++) {
		u[i+nf] += upen[i];
		for (j=nf; j<nvar2; j++) {
		    jmat[i][j] += ipen[k];
		    JJ  [i][j] += ipen[k];
		    k++;
		    }
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
    if (debug >0 && whichcase==0) {
	fprintf(stderr, "coef" );
	for (i=0; i<nvar2; i++) fprintf(stderr," %f", beta[i]);
	fprintf(stderr, "\nU   ");
	for (i=0; i<nvar2; i++) fprintf(stderr," %f", u[i]);
	fprintf(stderr, "\n");
	}
    if (debug >1 && whichcase==0) {
	fprintf(stderr, "Imat\n");
	for (i=0; i<nvar; i++) {
	    for (j=0; j<nvar2; j++) fprintf(stderr,"  %f", jmat[i][j]);
	    fprintf(stderr, "\n");
	    }
	fprintf(stderr,"fdiag\n");
	for (i=0; i<nvar2; i++) fprintf(stderr, "  %f", fdiag[i]);
	fprintf(stderr, "\n"); fflush(stderr);
	}
    return(loglik);
    }
