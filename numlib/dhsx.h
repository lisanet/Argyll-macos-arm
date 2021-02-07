#ifndef DHSX_H
#define DHSX_H

/* A general purpose downhill simplex multivariate optimser, */
/* as described in */
/* "Numerical Recipes in C", by W.H.Press, B.P.Flannery, */
/* S.A.Teukolsky & W.T.Vetterling. */

#ifdef __cplusplus
	extern "C" {
#endif

/* Down hill simplex function */
/* return 0 on sucess, 1 on failure due to excessive itterations */
/* Result will be in cp */
int dhsx(
 double *rv,			/* If not NULL, return the residual error */
 int di,				/* Dimentionality */
 double *cp,			/* Initial starting point, return minimum */
 double *s,				/* Size of initial search area */
 double ftol,			/* Finishing tollerance of error change */
 double athr,			/* Absolute return value threshold. (Set high to not use) */
 int maxit,				/* Maximum iterations allowed */
 double (*funk)(void *fdata, double *tp),		/* Error function to evaluate */
 void *fdata				/* Data needed by function */
);

double dhsx_funk(		/* Return function value */
 void *fdata,		/* Opaque data pointer */
 double tp[]);		/* Multivriate input value */

#ifdef __cplusplus
	}
#endif

#endif /* DHSX_H */
