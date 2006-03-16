#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "survproto.h"


static R_NativePrimitiveArgType agexact_t[]={INTSXP,  INTSXP,   INTSXP,   REALSXP, 
					     REALSXP,   INTSXP,    REALSXP,REALSXP, 
					     INTSXP, REALSXP,  REALSXP,  REALSXP, 
					     REALSXP, REALSXP, INTSXP,  REALSXP, 
					     INTSXP,  REALSXP,   REALSXP, REALSXP};

static R_NativePrimitiveArgType agfit2_t[]={ INTSXP,  INTSXP,  INTSXP, 
					     REALSXP,    REALSXP,    INTSXP, 
					     REALSXP,   REALSXP,  REALSXP,
					     INTSXP,   REALSXP,   REALSXP, 
					     REALSXP,        REALSXP, REALSXP, 
					     INTSXP,     REALSXP,    INTSXP,
					     REALSXP,      REALSXP, REALSXP};

static R_NativePrimitiveArgType agfit3_t[]={ INTSXP,  INTSXP,  INTSXP, 
					     REALSXP,    REALSXP,    INTSXP, 
					     REALSXP,   REALSXP,  REALSXP,
					     INTSXP,   INTSXP,  INTSXP,
					     INTSXP,    REALSXP,   REALSXP, 
					     REALSXP,        REALSXP, REALSXP,
					     INTSXP,     REALSXP,   
					     REALSXP,      REALSXP, REALSXP};

static R_NativePrimitiveArgType agfit5_a_t[]={INTSXP, INTSXP, REALSXP, 
					      REALSXP, REALSXP,
					      REALSXP,
					      INTSXP,  INTSXP,
					      REALSXP, REALSXP, REALSXP, 
					      REALSXP, 
					      INTSXP, INTSXP, INTSXP,
					      INTSXP,  INTSXP,
					      CLOSXP, CLOSXP, ENVSXP};

static R_NativePrimitiveArgType agfit5_b_t[]={INTSXP, INTSXP, INTSXP, 
					      INTSXP, REALSXP, REALSXP,
					      REALSXP,  REALSXP, REALSXP, 
					      INTSXP,  REALSXP, REALSXP, INTSXP, 
					      INTSXP, REALSXP, REALSXP,
					      CLOSXP,CLOSXP, ENVSXP};

static R_NativePrimitiveArgType agfit5_c_t[]={INTSXP, INTSXP, INTSXP,
					      INTSXP, REALSXP};

static R_NativePrimitiveArgType agfit_null_t[]={INTSXP,      INTSXP,   REALSXP, REALSXP, 
						INTSXP,  REALSXP,  REALSXP,
						INTSXP, REALSXP};


static R_NativePrimitiveArgType agmart2_t[]={INTSXP,     INTSXP,  REALSXP,   REALSXP, 
					    INTSXP, INTSXP, INTSXP, INTSXP,
					    INTSXP, REALSXP,   REALSXP,  
					    REALSXP, REALSXP};

static R_NativePrimitiveArgType agmart_t[]={INTSXP,     INTSXP,  REALSXP,   REALSXP, 
					    INTSXP, REALSXP,   REALSXP,      INTSXP, 
					    REALSXP};

static R_NativePrimitiveArgType agscore_t[]={INTSXP,       INTSXP,      REALSXP,
					     REALSXP,   INTSXP,     REALSXP,
					     REALSXP,  INTSXP,     REALSXP, REALSXP};

static R_NativePrimitiveArgType agsurv1_t[]={INTSXP,     INTSXP,  REALSXP,      REALSXP, 
					     INTSXP, REALSXP,   REALSXP,   INTSXP,
					     REALSXP,   REALSXP,      REALSXP, REALSXP,
					     INTSXP,  REALSXP,   REALSXP,REALSXP, 
					     INTSXP};

static R_NativePrimitiveArgType agsurv2_t[]={INTSXP,      INTSXP,    REALSXP, 
					     REALSXP,   INTSXP,   REALSXP, 
					     REALSXP,    REALSXP,     REALSXP, 
					     INTSXP,  REALSXP,        INTSXP,
					     REALSXP,   REALSXP};

static R_NativePrimitiveArgType agsurv3_t[]={INTSXP,    INTSXP,    INTSXP, 
					     INTSXP,  INTSXP,      REALSXP, 
					     REALSXP,    REALSXP,        REALSXP, 
					     REALSXP,   REALSXP,    INTSXP, 
					     REALSXP,    REALSXP,       REALSXP,
					     REALSXP,  REALSXP,    INTSXP};

static R_NativePrimitiveArgType coxdetail_t[]={INTSXP,   INTSXP,    INTSXP, 
					       REALSXP,        REALSXP,   INTSXP,  
					       REALSXP,    REALSXP,  REALSXP, 
					       REALSXP,       REALSXP,      REALSXP};

static R_NativePrimitiveArgType coxfit2_t[]={INTSXP,   INTSXP,    INTSXP, 
					     REALSXP,      INTSXP,    REALSXP, 
					     REALSXP,	REALSXP,   INTSXP,
					     REALSXP,     REALSXP,      REALSXP, 
					     REALSXP, REALSXP,  INTSXP, 
					     REALSXP,	REALSXP,  REALSXP,
					     REALSXP};

static R_NativePrimitiveArgType coxfit5_a_t[]={INTSXP, INTSXP, REALSXP, 
					       REALSXP, REALSXP,
					       REALSXP, INTSXP,  INTSXP,
					       REALSXP, REALSXP, REALSXP, 
					       REALSXP, 
					       INTSXP, INTSXP, INTSXP,
					       INTSXP,  INTSXP,
					       CLOSXP, CLOSXP,ENVSXP};

static R_NativePrimitiveArgType coxfit5_b_t[]={INTSXP, INTSXP, INTSXP, 
					       INTSXP, REALSXP, REALSXP,
					       REALSXP,  REALSXP, REALSXP, 
					       INTSXP,  REALSXP, REALSXP, INTSXP, 
					       INTSXP, REALSXP, REALSXP,
					       CLOSXP, CLOSXP,ENVSXP};

static R_NativePrimitiveArgType coxfit5_c_t[]={INTSXP, INTSXP, INTSXP, INTSXP, 
					       REALSXP};

static R_NativePrimitiveArgType coxmart_t[]={INTSXP,     INTSXP,    REALSXP, 
					     INTSXP, INTSXP,   REALSXP, 
					     REALSXP,     REALSXP};

static R_NativePrimitiveArgType coxph_wtest_t[]={INTSXP, INTSXP, REALSXP, REALSXP,
						 REALSXP, REALSXP};

static R_NativePrimitiveArgType coxscho_t[]={INTSXP,    INTSXP,    REALSXP, 
					     REALSXP,    REALSXP,    INTSXP,  
					     INTSXP, REALSXP};

static R_NativePrimitiveArgType coxscore_t[]={INTSXP,      INTSXP,    REALSXP, 
					      REALSXP,  INTSXP,   REALSXP, 
					      REALSXP, INTSXP,   REALSXP,
					      REALSXP};

static R_NativePrimitiveArgType pyears1_t[]={INTSXP,      INTSXP,      INTSXP, 
					     REALSXP,      REALSXP,       
					     INTSXP,   INTSXP, 
					     INTSXP,   REALSXP,    REALSXP, 
					     REALSXP,  INTSXP,    INTSXP, 
					     INTSXP,   REALSXP,    INTSXP, 
					     REALSXP,  REALSXP,   REALSXP, 
					     REALSXP,  REALSXP,  REALSXP};

static R_NativePrimitiveArgType pyears2_t[]={INTSXP,      INTSXP,   INTSXP, 
					     REALSXP,      REALSXP,    INTSXP,    INTSXP, 
					     INTSXP,   REALSXP, REALSXP,
					     REALSXP,  REALSXP,    REALSXP, 
					     REALSXP};

static R_NativePrimitiveArgType pyears3_t[]={INTSXP,    INTSXP,    INTSXP, 
					     INTSXP,      INTSXP, REALSXP, 
					     REALSXP,    REALSXP,    REALSXP, 
					     INTSXP,    INTSXP, REALSXP,
					     REALSXP,    INTSXP};



static R_NativePrimitiveArgType survdiff2_t[]={INTSXP,     INTSXP,    INTSXP, 
					       REALSXP,    REALSXP,       INTSXP, 
					       INTSXP,  INTSXP,	   REALSXP, 
					       REALSXP,    REALSXP,        REALSXP, 
					       REALSXP};

static R_NativePrimitiveArgType survfit2_t[]={INTSXP,      REALSXP,       REALSXP,
					      INTSXP,  INTSXP,  INTSXP, 
					      REALSXP,    REALSXP,    REALSXP, 
					      REALSXP};

static R_NativePrimitiveArgType survfit3_t[]={INTSXP, REALSXP, REALSXP,
					      INTSXP, INTSXP, INTSXP,
					      INTSXP, REALSXP,
					      REALSXP, REALSXP, REALSXP,
					      REALSXP, REALSXP, REALSXP,
					      REALSXP};

static R_NativePrimitiveArgType survindex2_t[]={INTSXP,     REALSXP,   INTSXP, 
						INTSXP, REALSXP,    INTSXP, 
						INTSXP,  INTSXP};
 
static R_NativePrimitiveArgType survreg2_t[]={INTSXP,   INTSXP,    INTSXP, 
					      REALSXP,          INTSXP,    REALSXP, REALSXP,
					      REALSXP,    REALSXP,  INTSXP, 
					      INTSXP,    REALSXP,    REALSXP, 
					      REALSXP,     INTSXP,  REALSXP,
					      REALSXP,   INTSXP,  INTSXP,
					      CLOSXP, ENVSXP};

static R_NativePrimitiveArgType survreg4_t[]={INTSXP,   INTSXP,       INTSXP, 
					      REALSXP,         INTSXP,       REALSXP, 
					      REALSXP,       REALSXP,  REALSXP,  
					      INTSXP,    INTSXP,  REALSXP,    
					      REALSXP,      REALSXP,
					      REALSXP,     INTSXP,     REALSXP,
					      REALSXP,   INTSXP,     INTSXP,
					      INTSXP,  	 INTSXP,
					      INTSXP,      INTSXP,   REALSXP,
					      CLOSXP, CLOSXP, CLOSXP, ENVSXP};

static R_NativePrimitiveArgType survreg3_t[]={INTSXP,   INTSXP,    INTSXP, 
					      REALSXP,          INTSXP,    REALSXP, REALSXP,
					      REALSXP,    REALSXP,  INTSXP, 
					      INTSXP,    REALSXP,    REALSXP, 
					      REALSXP,     INTSXP,  REALSXP,
					      REALSXP,   INTSXP,  INTSXP,
					      CLOSXP, ENVSXP};

static R_NativePrimitiveArgType survreg5_t[]={INTSXP,   INTSXP,       INTSXP, 
					      REALSXP,         INTSXP,       REALSXP,
					      REALSXP,       REALSXP,  REALSXP,  
					      INTSXP,    INTSXP,  REALSXP,    
					      REALSXP,      REALSXP,
					      REALSXP,     INTSXP,     REALSXP,
					      REALSXP,   INTSXP,     INTSXP,
					      INTSXP,  	 INTSXP,
					      INTSXP,      INTSXP,   REALSXP,
					      CLOSXP, CLOSXP, CLOSXP, ENVSXP };


static R_NativePrimitiveArgType char_date_t[]={INTSXP, INTSXP, STRSXP, INTSXP, INTSXP, INTSXP};

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

static const R_CMethodDef CEntries[] ={
	CDEF(agexact),
	CDEF(agfit2),
	CDEF(agfit3),
	CDEF(agfit5_a),
	CDEF(agfit5_b),
	CDEF(agfit5_c),
	CDEF(agfit_null),
	CDEF(agmart),
	CDEF(agmart2),
	CDEF(agscore),
	CDEF(agsurv1),
	CDEF(agsurv2),
	CDEF(agsurv3),
	CDEF(coxdetail),
	CDEF(coxfit2),
	CDEF(coxfit5_a),
	CDEF(coxfit5_b),
	CDEF(coxfit5_c),
	CDEF(coxmart),
	CDEF(coxph_wtest),
	CDEF(coxscho),
	CDEF(coxscore),
	CDEF(pyears1),
	CDEF(pyears2),
	CDEF(pyears3),
	CDEF(survdiff2),
	CDEF(survfit2),
	CDEF(survfit3),
	CDEF(survindex2),
	CDEF(survreg2),
	CDEF(survreg4),
	CDEF(survreg3),
	CDEF(survreg5),
	CDEF(char_date),
	{NULL,NULL,0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_survival(DllInfo *dll)
{
	R_registerRoutines(dll, CEntries, NULL,NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
