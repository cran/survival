
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "survS.h"
#include "survproto.h"
#include <math.h>
        



SEXP coxdetail_wrapper(SEXP   nusedx,   SEXP   nvarx,    SEXP   ndeadx,
               SEXP yx,        SEXP covar2x,   SEXP   stratax,
               SEXP scorex,    SEXP weightsx,  SEXP means2x,
               SEXP u2x,       SEXP varx,      SEXP   rmatx,
               SEXP nrisk2x,   SEXP workx){


    	PROTECT(nusedx = AS_INTEGER(nusedx));
    	PROTECT(nvarx = AS_INTEGER(nvarx));
    	PROTECT(ndeadx= AS_INTEGER(ndeadx));
    	PROTECT(yx = AS_NUMERIC(yx));
    	PROTECT(covar2x = AS_NUMERIC(covar2x));
    	PROTECT(stratax = AS_INTEGER(stratax));
    	PROTECT(scorex = AS_NUMERIC(scorex));
    	PROTECT(weightsx = AS_NUMERIC(weightsx));
    	PROTECT(means2x = AS_NUMERIC(means2x));
    	PROTECT(u2x = AS_NUMERIC(u2x));
    	PROTECT(varx = AS_NUMERIC(varx));
    	PROTECT(rmatx = AS_NUMERIC(rmatx));
    	PROTECT(nrisk2x = AS_NUMERIC(nrisk2x));
    	PROTECT(workx = AS_NUMERIC(workx));

        Sint *nusedxx = INTEGER_POINTER(nusedx);
        Rprintf("Size of nusedx  is %ld\n",sizeof(nusedxx));
        Sint *nvarxx = INTEGER_POINTER(nvarx);
        Sint *ndeadxx = INTEGER_POINTER(ndeadx);
        double *yxx = NUMERIC_POINTER(yx);
        double *covar2xx = NUMERIC_POINTER(covar2x);
        Sint *strataxx = INTEGER_POINTER(stratax);
        double *scorexx = NUMERIC_POINTER(scorex);
        double *weightsxx = NUMERIC_POINTER(weightsx);
        double *means2xx = NUMERIC_POINTER(means2x);
        double *u2xx = NUMERIC_POINTER(u2x);
        double *varxx = NUMERIC_POINTER(varx);
        double *rmatxx = NUMERIC_POINTER(rmatx);
        double *nrisk2xx = NUMERIC_POINTER(nrisk2x);
        double *workxx = NUMERIC_POINTER(workx);
        Rprintf("Size of rmat_double  is %ld\n with index 0 as %lf",sizeof(rmatxx),rmatxx[0]);


/*  */
        long i,j,k,person=0;
        long     nused, nvar;
        long     nrisk, ndead;
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
        long     itemp, deaths;
        long     ideath;
        double  efron_wt, d2;
        double risk;
        double  meanwt;
        double  wdeath;
        double  *start,
    	    *stop,
    	    *event;
        int     rflag;

        nused = *nusedxx;
        nvar  = *nvarxx;
        method= *means2xx;
        ndead = *ndeadxx;
        rflag = 1- rmatxx[0];

        /*
        **  Set up the ragged arrays
        */
        Rprintf("\n Trying to assign covar dmatrix \n");
        covar= dmatrix(covar2xx, nused, nvar);
        Rprintf("\n Trying to assign means dmatrix \n");
        means= dmatrix(means2xx, ndead, nvar);
        u    = dmatrix(u2xx, ndead, nvar);
        cmat = dmatrix(workxx,   nvar, nvar);
        cmat2= dmatrix(workxx + nvar*nvar, nvar, nvar);
        a = workxx + 2*nvar*nvar;
        a2= a+nvar;
        wmeans = a2+nvar;
        start =yxx;
        stop  =yxx + nused;
        event =yxx + nused +nused;

        /*
        ** Subtract the mean from each covar, as this makes the variance calc
        **  much more stable
        */
        Rprintf("\n Going thru loop of %ld iterations to set covar matrix \n",nvar); 
        for (i=0; i<nvar; i++) {
        	temp=0;
        	for (person=0; person<nused; person++) 
                    temp += covar[i][person];
        	temp /= nused;
        	wmeans[i] = temp;
        	for (person=0; person<nused; person++) 
                    covar[i][person] -=temp;
        }

        /*
        ** Zero out some arrays
        */
        for (i=0; i<ndead*nvar; i++) {
        	u2xx[i]=0;
        	means2xx[i] =0;
        }
        for (i=0; i<ndead*nvar*nvar; i++) varxx[i]=0;

        /*
        ** Now walk through the data
        */
        ideath=0;
        Rprintf("\n Walking thru the data of %ld values \n",nused); 
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
            	    deaths=0; wdeath=0;
            	    nrisk =0;
            	    for (k=person; k<nused; k++) {
                		if (start[k] < time) {
                		    nrisk++;
                		    if (rflag) rmatxx[ideath*nused +k] =1;
                		    risk = scorexx[k] * weightsxx[k];
                		    denom += risk;
                		    for (i=0; i<nvar; i++) {
                			a[i] += risk*covar[i][k];
                			for (j=0; j<=i; j++)
                			    cmat[i][j] += risk*covar[i][k]*covar[j][k];
                			}
                		    if (stop[k]==time && event[k]==1) {
                    			deaths += 1;
                    			wdeath += weightsxx[k];
                    			efron_wt += risk*event[k];
                    			meanwt += weightsxx[k];
                    			for (i=0; i<nvar; i++) {
                    			    a2[i]+= risk*covar[i][k];
                    			    for (j=0; j<=i; j++)
                    				    cmat2[i][j] += risk*covar[i][k]*covar[j][k];
                    			}
                			}
                		}
                		if (strataxx[k]==1) break;
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
                    			u[i][ideath] += weightsxx[k]*covar[i][k] - meanwt*temp2;
                    			for (j=0; j<=i; j++) {
                    			    temp3 =((cmat[i][j] - temp*cmat2[i][j]) -
                    					       temp2*(a[j]-temp*a2[j]))/d2;
                    			    temp3 *= meanwt;
                    			    varxx[i + j*nvar + ideath*nvar*nvar] +=temp3;
                    			    if (j<i)
                    				varxx[j + i*nvar + ideath*nvar*nvar] +=temp3;
                    			    }
                    			}
                		    }
                		person++;
                		if (strataxx[k]==1) break;
            		}
            	    strataxx[ideath]= person;     /* index of the death */
            	    scorexx[ideath] = wdeath;     /* weighted number of events */
            	    start[ideath] = deaths;     /* number of deaths */
            	    stop[ideath]  = nrisk;      /* number at risk   */
            	    event[ideath] = hazard;     /* increment to the hazard */
            	    weightsxx[ideath]=varhaz;     /* increment to the hazard variance */
            	    nrisk2xx[ideath]= denom ;     /* weighted number at risk */
            	    ideath++;
        	    }
        }   /* end  of accumulation loop */
        *ndeadxx = ideath;        
        Rprintf("Finished coxdetail\n");
        double *pout;
        long n = length(rmatxx);
        Rprintf("Length of rmat is %ld\n",n);
        SEXP out = PROTECT(allocVector(REALSXP, n));		         
        pout = REAL(out);
        long ii=0;
        long cntOnes = 0;
        while (ii < n){
            pout[ii] = rmatxx[ii];
            ii += 1;
            if (rmatxx[ii] == 1.0)
               cntOnes += 1;
        }
        Rprintf("size of rmatx is %ld with %ld one values",n,cntOnes);
        UNPROTECT(15);
        return out;
}
 


        
