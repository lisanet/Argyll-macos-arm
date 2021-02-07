#ifndef GNEWT_H
#define GNEWT_H

/* Global newton non-linear equation solver */

/*
 * Copyright 2018 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifdef __cplusplus
	extern "C" {
#endif

/*

	return values:

 	0  Success

	1  Ran out of itterations

	2  Inverting Jacobian failed

 */

int gnewt(
	void *fdata,	/* Opaque pointer to pass to fcn() and jac() */
	void (*fcn)(void *fdata, int n, double *x, double *fvec),
					/* Pointer to function we are solving */
	void (*jac)(void *fdata, int n, double *x, double **fjac),
					/* Function to compute jacobian */
	int n,			/* Number of functions and variables */
	double x[],		/* Initial solution estimate, returns final solution */
	double rfvec[],	/* Optionaly return soln. function values */
	double xtol,	/* Desired tollerance of root */
	double ftol,	/* Desired tollerance of the solution */
	int maxfcn,		/* Maximum number of function itterations */
	int maxjac		/* Maximum number of jacobian itterations */
);

#ifdef __cplusplus
	}
#endif

#endif /* GNEWT_H */










