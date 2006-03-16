/** Main R changes 
*** long to int everywhere.  Penalised functions get
*** some extra arguments, init_callback functions not needed
**/

/*
**  SCCS @(#)survproto.h	5.2 10/28/98
** prototypes of all the survival functions
**  along with a few macros
*/ 

void agexact(int *maxiter,  int *nusedx,   int *nvarx,   double *start, 
	     double *stop,   int *event,    double *covar2,double *offset, 
	     int   *strata, double *means,  double *beta,  double *u, 
	     double *imat2,  double loglik[2], int *flag,  double *work, 
	     int   *work2,  double *eps,    double *tol_chol, double *sctest);

void agfit2( int   *maxiter,  int   *nusedx,  int   *nvarx, 
	     double *start,    double *stop,    int   *event, 
	     double *covar2,   double *offset,  double *weights,
	     int   *strata,   double *means,   double *beta, 
	     double *u,        double *imat2,   double loglik[2], 
	     int   *flag,     double *work,    int   *end,
	     double *eps,      double *tol_chol,double *sctest);

void agfit3( int   *maxiter,  int   *nusedx,  int   *nvarx, 
	     double *start,    double *stop,    int   *event, 
	     double *covar2,   double *offset,  double *weights,
	     int   *nstrat,   int   *strata,  int   *sort1,
	     int   *sort2,    double *means,   double *beta, 
	     double *u,        double *imat2,   double loglik[2], 
	     int   *flag,     double *work,   
	     double *eps,      double *tol_chol, double *sctest);

void agfit5_a(int *nusedx, int *nvarx, double *yy, 
	       double *covar2, double *offset2,
	       double *weights2, 
	       int   *strata,  int   *sort,
	       double *means, double *beta, double *u, 
	       double *loglik, 
	       int *methodx, int *ptype2, int *pdiag2,
	      int *nfrail,  int *frail2,
	      void *fexpr1,void *fexpr2, void *rho);

void agfit5_b(int *maxiter, int *nusedx, int *nvarx, 
	       int *strata, double *beta, double *u,
	       double *imat2,  double *jmat2, double *loglik, 
	       int *flag,  double *eps, double *tolerch, int *methodx, 
	       int *nfrail, double *fbeta, double *fdiag,
	      /* R callback information */
	      void *fexpr1, void *fexpr2, void *rho);

void agfit5_c(int *nusedx, int *nvar, int *strata,
	      int *methodx, double *expect);

void agfit_null(int   *n,      int   *method,   double *start, double *stop, 
		int   *event,  double * offset,  double *weights,
		int   *strata, double loglik[2]);
/*
void aghaz2(int   *n,     double *start,   double *stop,   int   *event, 
	    double *score, int   * strata, double *hazard, double * cumhaz);
*/
void agmart(int   *n,     int   *method,  double *start,   double *stop, 
	    int   *event, double *score,   double *wt,      int   *strata, 
	    double *resid);

void agmart2(int   *n,     int   *method,  double *start,   double *stop, 
	    int   *event,  int   *nstrat,  int *strata,    int *sort1,
	    int   *sort2,  double *score,   double *wt,      
	     double *resid,  double *haz);

void agres12(int   *nx,     int   *nvarx,   double *y,    double *covar2, 
	     int   *strata, double *score,   int *method, double *resid2, 
	     double *a);

void agscore(int   *nx,       int   *nvarx,      double *y,
	     double *covar2,   int   *strata,     double *score,
	     double *weights,  int   *method,     double *resid2, double *a);

void agsurv1(int   *sn,     int   *snvar,  double *y,      double *score, 
	     int   *strata, double *surv,   double *varh,   int   *snsurv,
	     double *xmat,   double *d,      double *varcov, double *yy,
	     int   *shisn,  double *hisy,   double *hisxmat,double *hisrisk, 
	     int   *hisstrat);

void agsurv2(int   *sn,      int   *snvar,    double *y, 
	     double *score,   int   *strata,   double *surv, 
	     double *varh,    double *xmat,     double *varcov, 
	     int   *snsurv,  double *d,        int   *sncurve,
             double *newx,    double *newrisk);

void agsurv3(int   *sn,    int   *snvar,    int   *sncurve, 
	     int   *snpt,  int   *sse,      double *score, 
	     double *sy,    double *r,        double *coef, 
	     double *var,   double *cmean,    int   *scn, 
	     double *cy,    double *cx,       double *ssurv,
	     double *varh,  double *sused,    int   *smethod);

void chinv2  (double **matrix, int n);
int cholesky2(double **matrix, int n, double toler);
void chsolve2(double **matrix, int n, double *y);
void chinv3(double **matrix , int n, int m, double *fdiag);
int cholesky3(double **matrix, int n, int m, double *diag, double toler);
void chsolve3(double **matrix, int n, int m, double *diag, double *y);

void coxdetail(int   *nusedx,   int   *nvarx,    int   *ndeadx, 
	       double *y,        double *covar2,   int   *strata,  
	       double *score,    double *weights,  double *means2, 
	       double *u2,       double *var,      double *work);

void coxfit2(int   *maxiter,   int   *nusedx,    int   *nvarx, 
	     double *Time,      int   *status,    double *covar2, 
	     double *offset,	double *weights,   int   *strata,
	     double *means,     double *beta,      double *u, 
	     double *imat2,     double loglik[2],  int   *flag, 
	     double *work,	double *eps,       double *tol_chol,
	     double *sctest);

void coxfit5_a(int *nusedx, int *nvarx, double *yy, 
	       double *covar2, double *offset2,
	       double *weights2, int *strata,  int *sorted,
	       double *means, double *beta, double *u, 
	       double *loglik, 
	       int *methodx, int *ptype2, int *pdiag2,
	       int *nfrail,  int *frail2,
	       void *fexpr1, void *fexpr2, void *rho
    ) ;
void coxfit5_b(int *maxiter, int *nusedx, int *nvarx, 
	       int *strata, double *beta, double *u,
	       double *imat2,  double *jmat2, double *loglik, 
	       int *flag,  double *eps, double *tolerch, int *methodx, 
	       int *nfrail, double *fbeta, double *fdiag,
	       void *fexpr1, void *fexpr2, void *rho);

void coxfit5_c (int *nusedx, int *nvar, int *strata, int *methodx, 
		double *expect);

void coxfit_null(int   *nusedx,    int   *method,   double *Time, 
		 int   *status,    double *score,    double *weights, 
		 int   *strata,    double *loglik, double *resid);

void coxhaz2(int   *n,      double *score,   int   *mark, 
	     int   *strata, double *hazard,  double *cumhaz);

void coxmart(int   *sn,     int   *method,    double *Time, 
	     int   *status, int   * strata,   double *score, 
	     double *wt,     double *expect);

void coxph_wtest(int *nvar2, int *ntest, double *var, double *b,
                 double *scratch, double *tolerch);

void coxscho(int   *nusedx,    int   *nvarx,    double *y, 
	     double *covar2,    double *score,    int   *strata,  
	     int   *method2,   double *work);

void coxscore(int   *nx,      int   *nvarx,    double *y, 
	      double *covar2,  int   *strata,   double *score, 
	      double *weights, int   *method,   double *resid2,
	      double *scratch);

double **dmatrix(double *array, int ncol, int nrow);

void init_doloop(int min, int max);
int doloop      (int nloops, int *index);

/* void init_coxcall1(long *ptr1, vector **ptr2);
void init_coxcall2(long *ptr1, vector **ptr2);
*/
void cox_callback (int which, double *coef, double *first, 
                   double *second, double *penalty, int *flag, int p, void *fexpr, void *rho);

void pyears1(int   *sn,      int   *sny,      int   *sdoevent, 
	     double *sy,      double *weight,       
             int   *sedim,   int   *efac, 
	     int   *edims,   double *secut,    double *expect, 
	     double *sedata,  int   *sodim,    int   *ofac, 
	     int   *odims,   double *socut,    int   *smethod, 
	     double *sodata,  double *pyears,   double *pn, 
	     double *pcount,  double *pexpect,  double *offtable);

void pyears2(int   *sn,      int   *sny,   int   *sdoevent, 
	     double *sy,      double *wt,    int   *sodim,    int   *ofac, 
	     int   *odims,   double *socut, double *sodata,
	     double *pyears,  double *pn,    double *pcount, 
	     double *offtable);

void pyears3(int   *sdeath,    int   *sn,    int   *sedim, 
	     int   *efac,      int   *edims, double *secut, 
	     double *expect,    double *sx,    double *y, 
	     int   *sntime,    int   *sngrp, double *times,
	     double *esurv,     int   *nsurv);


double pystep(int nc,        int  *index,  int  *index2,   double *wt, 
	      double *data,  int *fac,    int *dims,     double **cuts, 
	      double step,   int  edge);

/** This doesn't seem to exist, but since nothing calls it we probably don't care.
    Does this mean we could lose the space allocated to savediag?
int rnewton(int    *maxiter,   int  n,        int  nvar,        double *beta, 
	    double *u,         double **imat, double loglik[2], double eps,
	    void (*dolk)(),    void (*doimat)(),  double *tol_chol,
	    double *newbeta,   double *savediag,  int debug);
*/

void surv_callback(double *z, double *dist, int n, void *fn, void *rho);

void survdiff2(int   *nn,     int   *nngroup,    int   *nstrat, 
	       double *rho,    double *Time,       int   *status, 
	       int   *group,  int   *strata,	   double *obs, 
	       double *exp,    double *var,        double *risk, 
	       double *kaplan);

void survfit2(int   *sn,      double *y,       double *wt,
	      int   *strata,  int   *method,  int   *error, 
	      double *mark,    double *surv,    double *varh, 
	      double *risksum  );

void survfit3(int   *sn,              double *y,               double *wt,
	      int   *strata,          int   *method,          int   *error,
	      int   *nstrat,          double *ntimes_strata,  
	      double *timelist,	       double *weighted_event,  double *surv,
	      double *varh,	       double *risksum,         double *enter,
	      double *exit_censored);

/* void survindex2(int   *n,          double *stime,      int   *strata, */
/* 		int   *ntime,      double *time,       int   *nstrat, */
/* 		int   *o_n_risk,   int   *o_n_event,  double *o_surv, */
/* 		double *o_std_err,  double *o_upper,    double *o_lower,  */
/* 		int   *n_risk,     int   *n_event,    double *surv, */
/* 		double *std_err,    double *upper,      double *lower, */
/* 		double *new_start,  int   *num_extend, int   *times_strata, */
/* 		double *temp_times); */
void survindex2(int   *n,     double *stime,   int   *strata, 
		int   *ntime, double *time,    int   *nstrat, 
		int   *indx,  int   *indx2);
 
void survreg2(int   *maxiter,   int   *nx,    int   *nvarx, 
	     double *y,          int   *ny,    double *covar2, double *wtx,
	     double *offset2,    double *beta,  int   *nstratx, 
	     int   *stratax,    double *ux,    double *imatx, 
	     double *loglik,     int   *flag,  double *eps,
	     double *tol_chol,   int   *dist,  int   *ddebug,
	      void *placeholder1, void *placeholder2);

void survreg4(int   *maxiter,   int   *nx,       int   *nvarx, 
	      double *y,         int   *ny,       double *covar2, 
	      double *wt2,       double *offset2,  double *beta,  
	     int   *nstratx,    int   *stratax,  double *ux,    
	     double *imatx,      double *jmatx,
	     double *loglik,     int   *flag,     double *eps,
	     double *tol_chol,   int   *dist,     int   *ddebug,
             int *ptype2,  	 int   *pdiag2,
	     int *nfrail2,      int   *frail2,   double *fdiag2,
	     void *fexpr1, void *fexpr2, void *placeholder, void *rho);

void survreg3(int   *maxiter,   int   *nx,    int   *nvarx, 
	     double *y,          int   *ny,    double *covar2, double *wtx,
	     double *offset2,    double *beta,  int   *nstratx, 
	     int   *stratax,    double *ux,    double *imatx, 
	     double *loglik,     int   *flag,  double *eps,
	     double *tol_chol,   int   *dist,  int   *ddebug,
	      void* density, void *rho);

void survreg5(int   *maxiter,   int   *nx,       int   *nvarx, 
	      double *y,         int   *ny,       double *covar2, 
	      double *wt2,       double *offset2,  double *beta,  
	     int   *nstratx,    int   *stratax,  double *ux,    
	     double *imatx,      double *jmatx,
	     double *loglik,     int   *flag,     double *eps,
	     double *tol_chol,   int   *dist,     int   *ddebug,
             int *ptype2,  	 int   *pdiag2,
	     int *nfrail2,      int   *frail2,   double *fdiag2,
	      void *fexpr1, void *fexpr2, void *density, void *rho );

void char_date(int *n, int *order, char *cdate, int *month, int*day, int *year);
