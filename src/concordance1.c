/* Automatically generated from all.nw using noweb */
#include "survS.h"
SEXP concordance1(SEXP y, SEXP wt2,  SEXP indx2, SEXP ntree2) {
    int i, j, k, index;
    int child, parent;
    int n, ntree;
    double *time, *status;
    double *twt, *nwt, *count;
    
    SEXP count2;
    double *wt;
    int    *indx;
    
    n = nrows(y);
    ntree = asInteger(ntree2);
    wt = REAL(wt2);
    indx = INTEGER(indx2);
    
    time = REAL(y);
    status = time + n;
    PROTECT(count2 = allocVector(REALSXP, 4));
    count = REAL(count2);
    twt = (double *) R_alloc(2*ntree, sizeof(double));
    nwt = twt + ntree;
    for (i=0; i< 2*ntree; i++) twt[i] =0.0;
    for (i=0; i<4; i++) count[i]=0.0;
    for (i=n-1; i>=0; ) {
        if (status[i]==1) { /* process all tied deaths at this point */
            for (j=i; j>=0 && status[j]==1 && time[j]==time[i]; j--) {
                index = indx[j];
                for (k=i; k>j; k--) count[3] += wt[j]*wt[k]; /* tied on time */
                count[2] += wt[j] * nwt[index];              /* tied on x */
                child = (2*index) +1;  /* left child */
                if (child < ntree)
                    count[0] += wt[j] * twt[child];  /*left children */
                child++;
                if (child < ntree)
                    count[1] += wt[j] * twt[child]; /*right children */
                
                while (index >0) {  /* walk up the tree  */
                    parent = (index-1)/2;
                    if (index & 1)   /* I am the right child */
                        count[1] += wt[j] * (twt[parent] - twt[index]);
                    else count[0] += wt[j] * (twt[parent] - twt[index]);
                    index = parent;
                    }
                }
            }                    
        else j = i-1;
        
        /* Add the weights for these obs into the tree */
        for (; i>j; i--) {
            index = indx[i];
            nwt[index] += wt[i];
            twt[index] += wt[i];
            while (index >0) {
                index = (index-1)/2;
                twt[index] += wt[i];
                }
            }
        }
        
    UNPROTECT(1);
    return(count2);
}
SEXP concordance2(SEXP y,     SEXP wt2,  SEXP indx2, SEXP ntree2,
                  SEXP sortstop, SEXP sortstart) {
    int i, j, k, index;
    int child, parent;
    int n, ntree;
    int istart, iptr, jptr;
    double *time1, *time2, *status, dtime;
    double *twt, *nwt, *count;
    int *sort1, *sort2;
    
    SEXP count2;
    double *wt;
    int    *indx;
    
    n = nrows(y);
    ntree = asInteger(ntree2);
    wt = REAL(wt2);
    indx = INTEGER(indx2);
    sort2 = INTEGER(sortstop);
    sort1 = INTEGER(sortstart);
    
    time1 = REAL(y);
    time2 = time1 + n;
    status= time2 + n;
    PROTECT(count2 = allocVector(REALSXP, 4));
    count = REAL(count2);
    twt = (double *) R_alloc(2*ntree, sizeof(double));
    nwt = twt + ntree;
    for (i=0; i< 2*ntree; i++) twt[i] =0.0;
    for (i=0; i<4; i++) count[i]=0.0;
    istart = 0;  /* where we are with start times */
    for (i=0; i<n; ) {
        iptr = sort2[i];  /* In  reverse death time order */
        if (status[iptr]==1) {
            /* Toss people out of the tree */
            dtime = time2[iptr];
            for (; istart < n && time1[sort1[istart]] >= dtime; istart++) {
                jptr <- sort1[istart];
                index = indx[jptr];
                nwt[index] -= wt[jptr];
                twt[index] -= wt[jptr];
                while (index >0) {
                    index = (index-1)/2;
                    twt[index] -= wt[jptr];
                    }
                }
                
            /* Process deaths */
            for (j=i; j <n && status[sort2[j]]==1 && time2[sort2[j]]==dtime; j++) {
                jptr =  sort2[j];
                index = indx[jptr];
                for (k=i; k<j; k++) count[3] += wt[jptr]*wt[sort2[k]]; 
                count[2] += wt[jptr] * nwt[index];            /* tied on x */
                child = (2*index) +1;   /* left child */
                if (child < ntree) count[0] += wt[jptr] * twt[child];
                child++;
                if (child < ntree) count[1] += wt[jptr] * twt[child];
     
                while (index >0) {  /* walk up the tree  */
                    parent = (index-1)/2;
                    if (index &1)   /* I am the right child */
                         count[1] += wt[jptr] * (twt[parent] - twt[index]);
                    else count[0] += wt[jptr] * (twt[parent] - twt[index]);
                    index = parent;
                    }
                }                    
            }
        else j = i+1;

        /* Add the weights for these obs into the tree */
        for (; i<j; i++) {
            iptr = sort2[i];
            index = indx[iptr];
            nwt[index] += wt[iptr];
            twt[index] += wt[iptr];
            while (index >0) {
                index = (index-1)/2;
                twt[index] += wt[iptr];
                }
            }
        }
        
    UNPROTECT(1);
    return(count2);
}
