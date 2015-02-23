/* Automatically generated from all.nw using noweb */
#include <math.h>
#include "survS.h" 
#include "survproto.h"

SEXP agfit4(SEXP surv2,      SEXP covar2,    SEXP strata2,
            SEXP weights2,  SEXP offset2,    SEXP ibeta2,
            SEXP sort12,     SEXP sort22,    SEXP method2,
            SEXP maxiter2,   SEXP  eps2,     SEXP tolerance2) { 

    int i,j,k,person;
    int indx2, istrat, p;
    int ksave, nrisk, ndeath;
    int nused, nvar;

    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *oldbeta, *maxbeta;
    double *a2, **cmat2;
    double *eta;
    double  denom, zbeta, risk;
    double  time;
    double  temp, temp2;
    double  newlk =0;
    int     halving;    /*are we doing step halving at the moment? */
    double  tol_chol, eps;
    double  meanwt;
    int itemp, deaths;
    double efron_wt, d2, meaneta;

    /* inputs */
    double *start, *stop, *event;
    double *weights, *offset;
    int *sort1, *sort2, maxiter;
    int *strata;
    double method;  /* saving this as double forces some double arithmetic */

    /* returned objects */
    SEXP imat2, means2, beta2, u2, loglik2;
    double *beta, *u, *loglik, *means;
    SEXP sctest2, flag2, iter2;
    double *sctest;
    int *flag, *iter;
    SEXP rlist;
    static const char *outnames[]={"coef", "u", "imat", "loglik", "means",
                                   "sctest", "flag", "iter", ""};
    int nprotect;  /* number of protect calls I have issued */

    /* get sizes and constants */
    nused = nrows(covar2);
    nvar  = ncols(covar2);
    method= asInteger(method2);
    eps   = asReal(eps2);
    tol_chol = asReal(tolerance2);
    maxiter = asInteger(maxiter2);
  
    /* input arguments */
    start = REAL(surv2);
    stop  = start + nused;
    event = stop + nused;
    weights = REAL(weights2);
    offset = REAL(offset2);
    sort1  = INTEGER(sort12);
    sort2  = INTEGER(sort22);
    strata = INTEGER(strata2);

    /*
    ** scratch space
    **  nvar: a, a2, newbeta, maxbeta
    **  nvar*nvar: cmat, cmat2
    **  n: eta
    */
    eta = (double *) R_alloc(nused + 4*nvar + 2*nvar*nvar, sizeof(double));
    a = eta + nused;
    a2 = a +nvar;
    maxbeta = a2 + nvar;
    oldbeta = maxbeta + nvar;

    /*
    **  Set up the ragged arrays
    **  covar2 might not need to be duplicated, even though
    **  we are going to modify it, due to the way this routine was
    **  was called.  In this case NAMED(covar2) will =0
    */
    PROTECT(imat2 = allocVector(REALSXP, nvar*nvar));
    nprotect =1;
    if (NAMED(covar2)>0) {
        PROTECT(covar2 = duplicate(covar2)); 
        nprotect++;
        }
    covar= dmatrix(REAL(covar2), nused, nvar);
    imat = dmatrix(REAL(imat2),  nvar, nvar);
    cmat = dmatrix(oldbeta+ nvar,   nvar, nvar);
    cmat2= dmatrix(oldbeta+ nvar + nvar*nvar, nvar, nvar);

    /*
    ** create the output structures
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    nprotect++;
    beta2 = SET_VECTOR_ELT(rlist, 0, duplicate(ibeta2));
    beta  = REAL(beta2);

    u2 =    SET_VECTOR_ELT(rlist, 1, allocVector(REALSXP, nvar));
    u = REAL(u2);

    SET_VECTOR_ELT(rlist, 2, imat2);
    loglik2 = SET_VECTOR_ELT(rlist, 3, allocVector(REALSXP, 2)); 
    loglik  = REAL(loglik2);

    means2 = SET_VECTOR_ELT(rlist, 4, allocVector(REALSXP, nvar));
    means  = REAL(means2);

    sctest2 = SET_VECTOR_ELT(rlist, 5, allocVector(REALSXP, 1));
    sctest =  REAL(sctest2);
    flag2  =  SET_VECTOR_ELT(rlist, 6, allocVector(INTSXP, 1));
    flag   =  INTEGER(flag2);
    iter2  =  SET_VECTOR_ELT(rlist, 7, allocVector(INTSXP, 1));
    iter = INTEGER(iter2);
    
    /*
    ** Subtract the mean from each covar, as this makes the variance
    **  computation much more stable
    */
    for (i=0; i<nvar; i++) {
        temp=0;
        for (person=0; person<nused; person++) temp += covar[i][person];
        temp /= nused;
        means[i] = temp;
        for (person=0; person<nused; person++)
            covar[i][person] -=temp;
        }
    ndeath =0;
    for (i=0; i<nused; i++) ndeath += event[i];
    
    /* First iteration, which has different ending criteria */
    for (i=0; i<nvar; i++) {
        u[i] =0;
        a[i] =0;
        for (j=0; j<nvar; j++) {
            imat[i][j] =0 ;
            cmat[i][j] =0;
        }
    }

    for (person=0; person<nused; person++) {
        zbeta = 0;      /* form the term beta*z   (vector mult) */
        for (i=0; i<nvar; i++)
            zbeta += beta[i]*covar[i][person];
        eta[person] = zbeta + offset[person];
    }

    /*
    **  'person' walks through the the data from 1 to n,
    **     sort1[0] points to the largest stop time, sort1[1] the next, ...
    **  'time' is a scratch variable holding the time of current interest
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
    meaneta =0;
    nrisk =0;
    newlk =0;
    for (person=0; person<nused;) {
        p = sort1[person];
        if (event[p]==0){
            nrisk++;
            meaneta += eta[p];
            risk = exp(eta[p]) * weights[p];
            denom += risk;
            for (i=0; i<nvar; i++) {
                a[i] += risk*covar[i][p];
                for (j=0; j<=i; j++)
                    cmat[i][j] += risk*covar[i][p]*covar[j][p];
                }
            person++;
            /* nothing more needs to be done for this obs */
        }
        else {
            time = stop[p];
            /*
            ** subtract out the subjects whose start time is to the right
            */
            for (; indx2<strata[istrat]; indx2++) {
                p = sort2[indx2];
                if (start[p] < time) break;
                nrisk--;
                meaneta -= eta[p];
                risk = exp(eta[p]) * weights[p];
                denom -= risk;
                for (i=0; i<nvar; i++) {
                    a[i] -= risk*covar[i][p];
                    for (j=0; j<=i; j++)
                        cmat[i][j] -= risk*covar[i][p]*covar[j][p];
                    }
                }

            /*
            ** compute the averages over subjects with
            **   exactly this death time (a2 & c2)
            ** (and add them into a and cmat while we are at it).
            */
            efron_wt =0;
            meanwt =0;
            for (i=0; i<nvar; i++) {
                a2[i]=0;
                for (j=0; j<nvar; j++) {
                    cmat2[i][j]=0;
                    }
                }
            deaths=0;
            for (k=person; k<strata[istrat]; k++) {
                p = sort1[k];
                if (stop[p] < time) break;
                risk = exp(eta[p]) * weights[p];
                denom += risk;
                nrisk++;
                meaneta += eta[p];

                for (i=0; i<nvar; i++) {
                    a[i] += risk*covar[i][p];
                    for (j=0; j<=i; j++)
                        cmat[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                if (event[p]==1) {
                    deaths += event[p];
                    efron_wt += risk*event[p];
                    meanwt += weights[p];
                    for (i=0; i<nvar; i++) {
                        a2[i]+= risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat2[i][j] += risk*covar[i][p]*covar[j][p];
                        }
                }
                }
            ksave = k;
                
            /* 
            ** If the average eta value has gotton out of hand, fix it.
            ** We must avoid overflow in the exp function (~750 on Intel)
            ** and want to act well before that, but not take action very often.  
            ** One of the case-cohort papers suggests an offset of -100 meaning
            ** that etas of 50-100 can occur in "ok" data, so make it larger than this.
            */
            if (fabs(meaneta) > (nrisk *110)) {  
                meaneta = meaneta/nrisk;
                for (i=0; i<nused; i++) eta[i] -= meaneta;
                temp = exp(-meaneta);
                denom *= temp;
                for (i=0; i<nvar; i++) {
                    a[i] *= temp;
                    a2[i] *= temp;
                    for (j=0; j<nvar; j++) {
                        cmat[i][j]*= temp;
                        cmat2[i][j] *= temp;
                    }
                }
                meaneta =0;
            }
                
            /*
            ** Add results into u and imat for all events at this time point
            */
            meanwt /= deaths;
            itemp = -1;
            for (; person<ksave; person++) {
                p = sort1[person];
                if (event[p]==1) {
                    itemp++;
                    temp = itemp*method/(double) deaths;
                    d2 = denom - temp*efron_wt;
                    newlk +=  weights[p]*eta[p] -meanwt *log(d2);

                    for (i=0; i<nvar; i++) {
                        temp2 = (a[i] - temp*a2[i])/d2;
                        u[i] += weights[p]*covar[i][p] - meanwt*temp2;
                        for (j=0; j<=i; j++)
                            imat[j][i] += meanwt* (
                                        (cmat[i][j] - temp*cmat2[i][j])/d2-
                                           temp2*(a[j]-temp*a2[j])/d2);
                        }
                    }
                }
        }

        if (person == strata[istrat]) {
            istrat++;
            denom =0;
            meaneta=0;
            nrisk =0;
            indx2 = person;
            for (i=0; i<nvar; i++) {
                a[i] =0;
                for (j=0; j<nvar; j++) {
                    cmat[i][j]=0;
                    }
                }
        }
    }   /* end  of accumulation loop */
    loglik[0] = newlk;   /* save the loglik for iteration zero  */
    loglik[1] = newlk;

    /* Use the initial variance matrix to set a maximum coefficient */
    for (i=0; i<nvar; i++) 
        maxbeta[i] = 23/ sqrt(imat[i][i]/ndeath);

    /* Calculate the score test */
    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
        a[i] = u[i];
    *flag = cholesky2(imat, nvar, tol_chol);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */
    *sctest=0;
    for (i=0; i<nvar; i++)
        *sctest +=  u[i]*a[i];

    if (maxiter ==0) {
        *iter =0;
        loglik[1] = newlk;
        chinv2(imat, nvar);
        for (i=1; i<nvar; i++)
            for (j=0; j<i; j++)  imat[i][j] = imat[j][i];

        UNPROTECT(nprotect);
        return(rlist);
    }
    else {  
        /* Update beta for the next iteration
        **  Never complain about convergence on this first step or impose step
        **  halving.  That way someone can force one iter at a time.
        */
        for (i=0; i<nvar; i++) {
            oldbeta[i] = beta[i];
            beta[i] = beta[i] + a[i];
        }
    }
    /* main loop */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (*iter=1; *iter<= maxiter; (*iter)++) {
        for (i=0; i<nvar; i++) {
            u[i] =0;
            a[i] =0;
            for (j=0; j<nvar; j++) {
                imat[i][j] =0 ;
                cmat[i][j] =0;
            }
        }

        for (person=0; person<nused; person++) {
            zbeta = 0;      /* form the term beta*z   (vector mult) */
            for (i=0; i<nvar; i++)
                zbeta += beta[i]*covar[i][person];
            eta[person] = zbeta + offset[person];
        }

        /*
        **  'person' walks through the the data from 1 to n,
        **     sort1[0] points to the largest stop time, sort1[1] the next, ...
        **  'time' is a scratch variable holding the time of current interest
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
        meaneta =0;
        nrisk =0;
        newlk =0;
        for (person=0; person<nused;) {
            p = sort1[person];
            if (event[p]==0){
                nrisk++;
                meaneta += eta[p];
                risk = exp(eta[p]) * weights[p];
                denom += risk;
                for (i=0; i<nvar; i++) {
                    a[i] += risk*covar[i][p];
                    for (j=0; j<=i; j++)
                        cmat[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                person++;
                /* nothing more needs to be done for this obs */
            }
            else {
                time = stop[p];
                /*
                ** subtract out the subjects whose start time is to the right
                */
                for (; indx2<strata[istrat]; indx2++) {
                    p = sort2[indx2];
                    if (start[p] < time) break;
                    nrisk--;
                    meaneta -= eta[p];
                    risk = exp(eta[p]) * weights[p];
                    denom -= risk;
                    for (i=0; i<nvar; i++) {
                        a[i] -= risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat[i][j] -= risk*covar[i][p]*covar[j][p];
                        }
                    }

                /*
                ** compute the averages over subjects with
                **   exactly this death time (a2 & c2)
                ** (and add them into a and cmat while we are at it).
                */
                efron_wt =0;
                meanwt =0;
                for (i=0; i<nvar; i++) {
                    a2[i]=0;
                    for (j=0; j<nvar; j++) {
                        cmat2[i][j]=0;
                        }
                    }
                deaths=0;
                for (k=person; k<strata[istrat]; k++) {
                    p = sort1[k];
                    if (stop[p] < time) break;
                    risk = exp(eta[p]) * weights[p];
                    denom += risk;
                    nrisk++;
                    meaneta += eta[p];

                    for (i=0; i<nvar; i++) {
                        a[i] += risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat[i][j] += risk*covar[i][p]*covar[j][p];
                        }
                    if (event[p]==1) {
                        deaths += event[p];
                        efron_wt += risk*event[p];
                        meanwt += weights[p];
                        for (i=0; i<nvar; i++) {
                            a2[i]+= risk*covar[i][p];
                            for (j=0; j<=i; j++)
                                cmat2[i][j] += risk*covar[i][p]*covar[j][p];
                            }
                    }
                    }
                ksave = k;
                    
                /* 
                ** If the average eta value has gotton out of hand, fix it.
                ** We must avoid overflow in the exp function (~750 on Intel)
                ** and want to act well before that, but not take action very often.  
                ** One of the case-cohort papers suggests an offset of -100 meaning
                ** that etas of 50-100 can occur in "ok" data, so make it larger than this.
                */
                if (fabs(meaneta) > (nrisk *110)) {  
                    meaneta = meaneta/nrisk;
                    for (i=0; i<nused; i++) eta[i] -= meaneta;
                    temp = exp(-meaneta);
                    denom *= temp;
                    for (i=0; i<nvar; i++) {
                        a[i] *= temp;
                        a2[i] *= temp;
                        for (j=0; j<nvar; j++) {
                            cmat[i][j]*= temp;
                            cmat2[i][j] *= temp;
                        }
                    }
                    meaneta =0;
                }
                    
                /*
                ** Add results into u and imat for all events at this time point
                */
                meanwt /= deaths;
                itemp = -1;
                for (; person<ksave; person++) {
                    p = sort1[person];
                    if (event[p]==1) {
                        itemp++;
                        temp = itemp*method/(double) deaths;
                        d2 = denom - temp*efron_wt;
                        newlk +=  weights[p]*eta[p] -meanwt *log(d2);

                        for (i=0; i<nvar; i++) {
                            temp2 = (a[i] - temp*a2[i])/d2;
                            u[i] += weights[p]*covar[i][p] - meanwt*temp2;
                            for (j=0; j<=i; j++)
                                imat[j][i] += meanwt* (
                                            (cmat[i][j] - temp*cmat2[i][j])/d2-
                                               temp2*(a[j]-temp*a2[j])/d2);
                            }
                        }
                    }
            }

            if (person == strata[istrat]) {
                istrat++;
                denom =0;
                meaneta=0;
                nrisk =0;
                indx2 = person;
                for (i=0; i<nvar; i++) {
                    a[i] =0;
                    for (j=0; j<nvar; j++) {
                        cmat[i][j]=0;
                        }
                    }
            }
        }   /* end  of accumulation loop */

        *flag = cholesky2(imat, nvar, tol_chol);
        if (fabs(1-(loglik[1]/newlk))<= eps  && halving==0){ /* all done */
            loglik[1] = newlk;
            chinv2(imat, nvar);
            for (i=1; i<nvar; i++)
                for (j=0; j<i; j++)  imat[i][j] = imat[j][i];

            UNPROTECT(nprotect);
            return(rlist);
        }

        if (*iter < maxiter) { /*update beta */
            if (newlk < loglik[1])   {    /*it is not converging ! */
                halving =1;
                for (i=0; i<nvar; i++)
                    beta[i] = (oldbeta[i] + beta[i]) /2; /*half of old increment */
            }
            else {
                halving=0;
                loglik[1] = newlk;
                chsolve2(imat,nvar,u);

                for (i=0; i<nvar; i++) {
                    oldbeta[i] = beta[i];
                    beta[i] = beta[i] +  u[i];
                    if (beta[i]> maxbeta[i]) beta[i] = maxbeta[i];
                    else if (beta[i] < -maxbeta[i]) beta[i] = -maxbeta[i];
                }
            }
        }  
        R_CheckUserInterrupt();  /* be polite -- did the user hit cntrl-C? */
    } /*return for another iteration */
    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=1; i<nvar; i++)
        for (j=0; j<i; j++)  imat[i][j] = imat[j][i];

    UNPROTECT(nprotect);
    return(rlist);
}
