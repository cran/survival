#include "R.h"
#include "R_ext/Memory.h"
#define ALLOC(a,b) S_alloc(a,b)


#include "Rmath.h"
#ifndef erf
#define erf(x) (2*pnorm5((x)*M_SQRT2,0,1,1,0)-1)
#endif
#ifndef erfc
#define erfc(x) (2*pnorm5(-(x)*M_SQRT2,0,1,1,0))
#endif

