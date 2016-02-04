
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "survS.h"
#include "survproto.h"
SEXP coxdetail_wrapper(SEXP   nusedx,   SEXP   nvarx,    SEXP   ndeadx,
               SEXP y,        SEXP covar2,   SEXP   strata,
               SEXP score,    SEXP weights,  SEXP means2,
               SEXP u2,       SEXP var,      SEXP   rmat,
               SEXP nrisk2,   SEXP work){

        Sint *nusedx_int = INTEGER_POINTER(nusedx);
        Sint *nvarx_int = INTEGER_POINTER(nvarx);
        Sint *ndeadx_int = INTEGER_POINTER(ndeadx);
        double *y_double = NUMERIC_POINTER(y);
        double *covar2_double = NUMERIC_POINTER(covar2);
        Sint *strata_int = INTEGER_POINTER(strata);
        double *score_double = NUMERIC_POINTER(score);
        double *weights_double = NUMERIC_POINTER(weights);
        double *means2_double = NUMERIC_POINTER(means2);
        double *u2_double = NUMERIC_POINTER(u2);
        double *var_double = NUMERIC_POINTER(var);
        double *rmat_double = NUMERIC_POINTER(rmat);
        double *nrisk2_double = NUMERIC_POINTER(nrisk2);
        double *work_double = NUMERIC_POINTER(work);

        coxdetail(nusedx_int, nvarx_int, ndeadx_int, y_double,
                  covar2_double, strata_int, score_double, weights_double,
                  means2_double, u2_double, var_double,rmat_double, nrisk2_double,work_double);

	int n = length(rmat_double);
        double *pout;

	SEXP out = PROTECT(allocVector(REALSXP, n));
        pout = REAL(out);
        for (int i = 0; i < n; i++) {
            pout[i] = rmat_double[i];
        }
        UNPROTECT(1);
        return out;
}
