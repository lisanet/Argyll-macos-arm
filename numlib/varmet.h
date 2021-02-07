#ifndef VARMET_H
#define VARMET_H

/* Variable Metric multivariate minimiser */

/*
 * Copyright 2000, 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifdef __cplusplus
	extern "C" {
#endif

/* Variable Metric Gradient optimiser */
/* return 0 on sucess, 1 on failure due to excessive itterations */
/* Result will be in cp */
int varmet(
double *rv,				/* If not NULL, return the residual error */
int di,					/* Dimentionality */
double cp[],			/* Initial starting point */
double s[],				/* Size of initial search area */
double ftol,			/* Tollerance of error change to stop on */
int maxit,				/* Maximum iterations allowed */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
double (*dfunc)(void *fdata, double dp[], double tp[]),		/* Gradient function to evaluate */
void *fdata				/* Opaque data needed by function */
);

/* Example user function declarations */
double varmet_funk(		/* Return function value */
	void *fdata,		/* Opaque data pointer */
	double tp[]);		/* Multivriate input value */

/* Line in multi-dimensional space minimiser */
double brentnd(			/* vector multiplier return value */
double ax,				/* Minimum of multiplier range */
double bx,				/* Starting point multiplier of search */
double cx,				/* Maximum of multiplier range */
double ftol,			/* Tollerance to stop search */
double *xmin,			/* Return value of multiplier at minimum */		
int n,					/* Dimensionality */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
void *fdata,			/* Opaque data */
double pcom[],			/* Base vector point */
double xicom[]);		/* Vector that will be multiplied and added to pcom[] */

#ifdef __cplusplus
	}
#endif

#endif /* VARMET_H */
