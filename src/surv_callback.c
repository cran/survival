/*   SCCS @(#)surv_callback.c	1.1 01/31/99
** callback routines for the survreg "other" distributions
** This is modeled directly on the interface code for glm/gam
*/
/**
*** Rewritten for R 
*** pass a function closure down to survreg, evaluate it, and
*** find the answers in the list it returns
**/
#include "S.h"
#include "Rinternals.h"


/*
** This part is called by the survreg5 function, to get the distribution
*/
void surv_callback (double *z, double *dist, int n, SEXP fn, SEXP rho) {
    SEXP survlist, temp, data;
    int i;
    /* copy in the argument */
    PROTECT(data=allocVector(REALSXP,n));
    for (i=0;i<n;i++){
      REAL(data)[i]=z[i];
	}
    /* evaluate the function */
    PROTECT(survlist=eval(lang2(fn,data),rho));
    UNPROTECT(2);
    PROTECT(survlist);
    
    /* Grab the updated values from the list */
    PROTECT(temp=lang3(install("[["),survlist,install("density")));
    PROTECT(data=eval(temp,rho));
    if (!isNumeric(data))
                error("density:invalid type\n");
    for (i=0;i<length(data);i++){
      dist[i]=REAL(data)[i];
    }
    /* tidy up */
    UNPROTECT(3);	

    }
