/*  SCCS @(#)dmatrix.c	5.2 10/27/98
** set up ragged arrays, with #of columns and #of rows
*/
#include "survS.h"
#include "survproto.h"

double **dmatrix(double *array, int ncol, int nrow)
    {
    register int i;
    register double **pointer;

    pointer = (double **) ALLOC(nrow, sizeof(double *));
    for (i=0; i<nrow; i++) {
	pointer[i] = array;
	array += ncol;
	}
    return(pointer);
    }
