
/*
 * Code to test the downhill simplex minimiser
 *
 * Copyright 1999 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include "numlib.h"

/*       Final approximate solution:                                       */

double expect[9] = {
	-0.5706545E+00,
	-0.6816283E+00,
	-0.7017325E+00,
	-0.7042129E+00,
	-0.7013690E+00,
	-0.6918656E+00,
	-0.6657920E+00,
	-0.5960342E+00,
	-0.4164121E+00 };

double fcn(		/* Return function value */
	void *fdata,		/* Opaque data pointer */
	double tp[]);		/* Multivriate input value */

#define N 9

int main(void)
{
	double cp[N];		/* Function input values */
	double s[N];		/* Search area */
	double err;
	int j;
	int nprint = 0;		/* Itteration debugging print = off */
	int rv;

	error_program = "tdhsx";	/* Set global error reporting string */

	/*	 The following starting values provide a rough solution. */
	for (j = 0; j < N; j++) {
		cp[j] = -1.f;
		s[j] = 0.9;
	}

	nprint = 0;

	/*	 Set tol to the square root of the machine precision. */
	/*	 Unless high precision solutions are required, */
	/*	 this is the recommended setting. */

	rv = dhsx(
		&err,			/* return residual */
		N, 				/* Dimentionality */
		cp,				/* Initial starting point */
		s,				/* Size of initial search area */
		0.0000001,		/* Tollerance of error change to stop on */
		1e6,			/* Absolute error to stop on */
//		10.0,		/* Tollerance of error change to stop on */
//		0.0001,			/* Absolute error to stop on */
		1000,			/* Maximum iterations allowed */
		fcn, 			/* Error function to evaluate */
		NULL);			/* Opaque data needed by function */


	fprintf(stdout,"Return code %d, residual err = %f:\n",rv, err);
	for (j = 0; j < N; j++) {
		fprintf(stdout,"cp[%d] = %e, expect %e\n",j,cp[j],expect[j]);
	}

	return 0;
} /* main() */

/* Function being minimized */
double fcn(		/* Return function value */
void *fdata,		/* Opaque data pointer */
double tp[]			/* Multivriate input value */
) {
	double err, tt;
	double temp, temp1, temp2;
	int k;

	/* Function Body */
	err = 0.0;
	for (k = 0; k < N; ++k) {
		temp = (3.0 - 2.0 * tp[k]) * tp[k];
		temp1 = 0.0;
		if (k != 0) {
			temp1 = tp[k-1];
		}
		temp2 = 0.0;
		if (k != ((N)-1))
			temp2 = tp[k+1];
		tt = temp - temp1 - 2.0 * temp2 + 1.0;
		err += tt * tt;
	}
	err = sqrt(err);
// printf("Returning %16.14f\n",err);
	return err;
}












