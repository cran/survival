/*
 * Tear down a date string into its component parts
 */

#include <stdio.h>
#include <string.h>

char *id();

void
char_date(n, order, cdate, month, day, year)
    int  *n,
    order[],			/* 1=year, 2=month, 3=day */
    month[],
    day[],
    year[];
    char *cdate[];

{
    register int i,k, ii;
    register char *j;
    register char *cc;
    
    int what[3];
    int len;
    char tdate[10];

    for (i=0; i< *n; i++) {
	cc = cdate[i];
	for (j=cc; *j != '\0'; j++)  /* upper case to lower case */
	    if (strchr("ABCDEFGHIJKLMNOPQRSTUVWXYZ", *j) !=NULL)
		*j += 'a' - 'A';

	/*
	** If it is pure numeric, put in some delimiters based on
	**  assumptions
	*/
	len = strlen(cc);
	for (k=0; k<len; k++)
	    if ((cc[k] < '0') || (cc[k] > '9')) break;
	if (k>=len && (len>=5 && len<=8)) {
	    if (len==5)
		sprintf(tdate, "0%c/%c%c/%c%c", cc[0], cc[1], cc[2], cc[3],
						cc[4]);
	    else if (len==6)
		sprintf(tdate, "%c%c/%c%c/%c%c", cc[0], cc[1], cc[2], cc[3],
						 cc[4], cc[5]);
	    else {
		if (len==7) {
		    for (ii=7; ii>0; ii--) cc[ii]= cc[ii-1];
		    cc[0] = '0';
		    }
		if (order[0]==1)
		    sprintf(tdate, "%c%c%c%c/%c%c/%c%c", cc[0], cc[1], cc[2],
				    cc[3], cc[4], cc[5], cc[6], cc[7]);
		else if (order[1]==1)
		    sprintf(tdate, "%c%c/%c%c%c%c/%c%c", cc[0], cc[1], cc[2],
				    cc[3], cc[4], cc[5], cc[6], cc[7]);
		else
		    sprintf(tdate, "%c%c/%c%c/%c%c%c%c", cc[0], cc[1], cc[2],
				    cc[3], cc[4], cc[5], cc[6], cc[7]);
		}
	    cc = tdate;
	    }

	cc = id(cc, what, 0);
	cc = id(cc, what, 1);
	cc = id(cc, what, 2);
	if (*cc != '\0')  what[2] =0;

	if (what[0] <0) {
	    month[i] = -1*what[0];
	    day[i]   =  what[1];
	    year[i]  =  what[2];
	    }
	else if (what[1] <0) {
	    month[i] = -1*what[1];
	    day[i]   =  what[0];
	    year[i]  =  what[2];
	    }
	else for (k=0; k<3; k++) {
	    switch (order[k]) {
		case 1: year[i] = what[k]; break;
		case 2: month[i]= what[k]; break;
		case 3: day[i]  = what[k]; break;
		}
	     }
	}
    }


char *
id(str, array, k)
    char *str;
    int  array[];
    int  k;
{
    register int i;
    register char *j;

    /*skip any delimiters (leading blanks are always ok) */
    while (*str==' ') str++;
    if (k>0)
	if (strchr(" -/,", *str) !=NULL) str++;
    while (*str==' ') str++;

    if (*str=='\0') {
	array[k]=0;
	return(str);
	}

    if (strchr("0123456789", *str) ==NULL) {
	if      (strstr(str, "jan")==str) array[k] = -1;
	else if (strstr(str, "feb")==str) array[k] = -2;
	else if (strstr(str, "mar")==str) array[k] = -3;
	else if (strstr(str, "apr")==str) array[k] = -4;
	else if (strstr(str, "may")==str) array[k] = -5;
	else if (strstr(str, "jun")==str) array[k] = -6;
	else if (strstr(str, "jul")==str) array[k] = -7;
	else if (strstr(str, "aug")==str) array[k] = -8;
	else if (strstr(str, "sep")==str) array[k] = -9;
	else if (strstr(str, "oct")==str) array[k] = -10;
	else if (strstr(str, "nov")==str) array[k] = -11;
	else if (strstr(str, "dec")==str) array[k] = -12;
	else array[k] =0;

	/* pass over the rest of the string */
	while( *str!= '\0' && strchr("januaryfebmrchpilgstovd", *str)!=NULL)
		 str++;
	}
    else { /*is a number */
	i =0;
	while (*str!= '\0' &&  (j=strchr("0123456789", *str)) !=NULL) {
	    str++;
	    i = (10*i) +(*j - '0');
	    }
	array[k] = i;
	}

    return(str);
    }
