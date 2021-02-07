#ifndef POWELL_H
#define POWELL_H

/* Powell and Conjugate Gradient multivariate minimiser */

/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifdef __cplusplus
	extern "C" {
#endif

#define DFUNC_NRV	1e100

/* Standard interface for powell function */
/* return 0 on sucess, 1 on failure due to excessive itterations */
/* Result will be in cp */
/* Arrays start at 0 */
int powell(
double *rv,				/* If not NULL, return the residual error */
int di,					/* Dimentionality */
double cp[],			/* Initial starting point */
double s[],				/* Size of initial search area */
double ftol,			/* Tollerance of error change to stop on */
int maxit,				/* Maximum iterations allowed */
double (*funk)(void *fdata, double tp[]),		/* Error function to evaluate */
void *fdata,			/* Opaque data needed by func() */
void (*prog)(void *pdata, int perc),		/* Optional progress percentage callback */
void *pdata				/* Opaque data needed by prog() */
);

/* Conjugate Gradient optimiser using partial derivatives. */
/* return 0 on sucess, 1 on failure due to excessive itterations */
/* Result will be in cp */
int conjgrad(
double *rv,				/* If not NULL, return the residual error */
int di,					/* Dimentionality */
double cp[],			/* Initial starting point */
double s[],				/* Size of initial search area */
double ftol,			/* Tollerance of error change to stop on */
int maxit,				/* Maximum iterations allowed */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
double (*dfunc)(void *fdata, double dp[], double tp[]),		/* Gradient & function to evaluate */
						/* dfunc() should return DFUNC_NRV if it doesn't return function value */
void *fdata,			/* Opaque data needed by function */
void (*prog)(void *pdata, int perc),		/* Optional progress percentage callback */
void *pdata				/* Opaque data needed by prog() */
);

/* Example user function declarations */
double powell_funk(		/* Return function value */
	void *fdata,		/* Opaque data pointer */
	double tp[]);		/* Multivriate input value */

/* Line bracketing and minimisation routine. */
/* Return value at minimum. */
double linmin(
double cp[],		/* Start point, and returned value */
double xi[],		/* Search vector */
int di,				/* Dimensionality */
#ifdef ABSTOL
double ftol,		/* Absolute tolerance to stop on */
#else
double ftol,		/* Relative tolerance to stop on */
#endif
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
void *fdata			/* Opaque data for func() */
);

/* Line bracketing and minimisation routine using derivatives. */
/* Return value at minimum. */
double linmind(
double cp[],		/* Start point, and returned value */
double xi[],		/* Search vector */
int di,				/* Dimensionality */
#ifdef ABSTOL
double ftol,		/* Absolute tolerance to stop on */
#else
double ftol,		/* Relative tolerance to stop on */
#endif
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
double (*dfunc)(void *fdata, double dp[], double tp[]),		/* Gradient & function to evaluate */
void *fdata			/* Opaque data for func() */
);

#ifdef __cplusplus
	}
#endif

#endif /* POWELL_H */
