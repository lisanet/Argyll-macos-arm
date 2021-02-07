
/* Multi-dimentional minizer using Variable Metric method */
/* This is good for smoother, well behaved functions. */

/* Code is an original expression of the algorithms decsribed in */
/* "Numerical Recipes in C", by W.H.Press, B.P.Flannery, */
/* S.A.Teukolsky & W.T.Vetterling. */

/*
 * Copyright 2000, 2006, 2007, 2017 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD:
   Fix error handling to return status (malloc, excessive itters)
   Create to "safe" library ?
   Make standalone - ie remove numsup ?
 */

/* Note that all arrays are indexed from 0 */

#include "numsup.h"
#include "powell.h"
#include "varmet.h"

#undef DEBUG					/* Debug */

#if defined(PDEBUG) || defined(CDEBUG)
# undef DBG
# define DBG(xxx) printf xxx ;
#else
# undef DBG
# define DBG(xxx) 
#endif

#define FMAX(A,B) ((A) > (B) ? (A) : (B))
#define EPS   1.0e-10	/*  Machine precision. */
#define TOLX (4 * EPS)	/*  X value stop value */
#define MAXLEN 100.0	/*  Maximum step length */

void linesearch(int di, double cpold[], double fpold, double g[], double p[], double cpnew[],
	double *pfp, double maxstep, double (*func)(void *fdata, double tp[]), void *fdata);

/* return 0 on sucess, 1 on failure due to excessive itterions */
/* Result will be in cp */
/* Note that we could use gradient in line minimiser, */
/* but haven't bothered yet. */
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
) {
	int iter;
	double fp, sumsq, maxstep;
	double *sdir, sumsdir;			/* Search direction */
	double *dp, *lastdp;
	double **hessian;				/* Hessian matrix */
	double *hlastdp;				/* Hessian times lastdp */
	double *cpnew;					/* new cp value from linemin */
	double test;

	double den, fac, fad, fae;
	double sumdg;
	int i, j;

	sdir   = dvector(0, di-1);
	dp     = dvector(0, di-1);
	lastdp = dvector(0, di-1);
	hessian = dmatrix(0, di-1, 0, di-1);
	hlastdp = dvector(0, di-1);
	cpnew   = dvector(0, di-1);

	fp = (*dfunc)(fdata, dp, cp);

	/* Initial line direction and pde squared */
	sumsq = 0.0;
	for (i = 0; i < di ;i++) {
		sdir[i] = -dp[i];				
		sumsq += cp[i] * cp[i];
	}

	DBG((" initial fp %f dp %s\n", fp, debPdv(di, dp)));

	/* Initialize inverse Hessian to unity */
	for (i = 0; i < di ;i++) {
	 	for (j = 0; j < di ; j++) {
			if (i == j)
				hessian[i][j] = 1.0;
			else
				hessian[i][j] = 0.0;
		}
	}

	/* Maximum line search step size */
	maxstep = MAXLEN * FMAX(sqrt(sumsq), (double)di);

	/* Until we give up */
	for (iter = 0; iter < maxit; iter++) {

		/* Search in direction sdir */
		linesearch(di, cp, fp, dp, sdir, cpnew, &fp, maxstep, func, fdata);

		for (i = 0; i < di; i++) {
			sdir[i] = cpnew[i] - cp[i];			/* Update the line direction, */
			cp[i] = cpnew[i];					/* and the current point. */
		}

		/* Test for convergence on x */
		test = 0.0;
		for (i = 0 ; i < di; i++) {
			double tt = fabs(sdir[i]) / FMAX(fabs(cp[i]), 1.0);
			if (tt > test)
				test = tt;
		}

		if (test < TOLX) {
			break;
		}

		for (i = 0; i < di; i++)		/* Save previous partial deriv. */
			lastdp[i] = dp[i];

		(*dfunc)(fdata, dp, cp);		/* and get the new gradient. */

		test = 0.0;					/* Test for convergence on zero gradient. */
		den = FMAX(fp, 1.0);

		for (i = 0; i < di; i++) {
			double tt = fabs(dp[i]) * FMAX(fabs(cp[i]),1.0) / den;
			if (tt > test)
				test = tt;
		}

		if (test < ftol) {
			break;
		}

	for (i = 0 ; i < di; i++)
		lastdp[i] = dp[i] - lastdp[i];	/* Compute diference of gradients, */

	for (i = 0; i < di; i++) {	/* and difernce times current matrix. */
		hlastdp[i] = 0.0;
		for (j = 0; j < di; j++)
			hlastdp[i] += hessian[i][j] * lastdp[j];
	}

	/* Calculate dot products for the denominator */
	fac = fae = sumdg = sumsdir = 0.0;
	for (i = 0; i < di; i++) { 
		fac += lastdp[i] * sdir[i];
		fae += lastdp[i] * hlastdp[i];
		sumdg += lastdp[i] * lastdp[i];
		sumsdir += sdir[i] * sdir[i];
	}
	if (fac > sqrt(EPS * sumdg * sumsdir)) { 	/* Skip update if fac not sufficiently posive */
		fac = 1.0/fac;
		fad = 1.0/fae;
		/* The vector that makes BFGS diferent from DFP: */
		for (i = 0; i < di;i++)
			lastdp[i] = fac * sdir[i] - fad * hlastdp[i];
		for (i = 0; i < di;i++) {	/* The BFGS updating formula: */
			for (j = i; j < di; j++) {
				hessian[i][j] += fac * sdir[i] * sdir[j]
				              - fad * hlastdp[i] * hlastdp[j] + fae * lastdp[i] * lastdp[j];
				hessian[j][i] = hessian[i][j];
			}
		}
	}
	for (i = 0; i < di; i++) {	/* Now calculate the next direction to go, */
		sdir[i] = 0.0;
		for (j = 0; j < di; j++)
			sdir[i] -= hessian[i][j] * dp[j];
		}
	}

	free_dvector(cpnew, 0, di-1);
	free_dvector(hlastdp, 0, di-1);
	free_dmatrix(hessian, 0, di-1, 0, di-1);
	free_dvector(lastdp, 0, di-1);
	free_dvector(dp, 0, di-1);
	free_dvector(sdir, 0, di-1);
	if (rv != NULL)
		*rv = fp;
	return 0;
}

#define ALPHA 1.0e-4 			/* Ensures sucesscient decrease in function value. */
#define LTOLX 1.0e-7		/* Convergence criterion on linesearch. */

void linesearch(
int di,
double cpold[],
double fpold,
double dp[],		/* Partial derivative */
double sdir[],		/* Current value */
double cpnew[],
double *pfp,		/* Return value */
double maxstep,
double (*func)(void *fdata, double tp[]),
void *fdata
) {
	double sum, slope;
	double slen, slen_min;
	double test, fval;
	int i;

	for (sum = 0.0, i = 0; i < di; i++)
		sum += sdir[i] * sdir[i];
	sum = sqrt(sum);

	if (sum > maxstep) {
		for (i = 0; i < di; i++)
			sdir[i] *= maxstep/sum; 		/* Scale if attempted step is too big. */
	}
	for (slope = 0.0, i = 0; i < di; i++)
		slope += dp[i] * sdir[i];

	if (slope >= 0.0)
		error("Roundoff problem in linesearch.");

	test = 0.0;
	for (i = 0;i < di; i++) {
		double tt = fabs(sdir[i])/FMAX(fabs(cpold[i]), 1.0);
		if (tt > test)
			test = tt;
	}

	slen_min = TOLX/test;
	slen = 1.0;				/* Try full step */

	/* Start of iteration loop. */
	for (;;) {
		double slen_2 = slen, slen_t;

		for (i = 0; i < di;i++)
			cpnew[i] = cpold[i] + slen * sdir[i];

		*pfp = (*func)(fdata, cpnew);

		if (slen < slen_min) {
			for (i = 0; i < di; i++)
				cpnew[i] = cpold[i];
			return;

		} else if (*pfp <= (fpold + ALPHA * slen * slope))
			return;	

		/* Backtracking */
		else {
			/* First time through */
			if (slen == 1.0)
				slen_t = -slope/(2.0 * (*pfp - fpold - slope));
			/* 2nd and subsequent times through */
			else {
				double aa, bb;
				double rhs_1, rhs_2;

				rhs_1 = *pfp - fpold - slen * slope;
				rhs_2 = fval - fpold - slen_2 * slope;
				aa = (rhs_1/(slen * slen) - rhs_2/(slen_2 * slen_2))/(slen - slen_2);
				bb = (-slen_2 * rhs_1/(slen * slen)+slen * rhs_2/(slen_2 * slen_2))/(slen - slen_2);
				if (aa == 0.0)
					slen_t = -slope/(2.0 * bb);
				else {
					double dd = bb * bb - 3.0 * aa * slope;
					if (dd < 0.0)
						slen_t = 0.5 * slen;
					else if (bb <= 0.0)
						slen_t = (-bb + sqrt(dd))/(3.0 * aa);
					else
						slen_t = -slope/(bb + sqrt(dd));
				}
				if (slen_t > 0.5 * slen)
					slen_t = 0.5 * slen;
			}
		}
		slen_2 = slen;
		fval = *pfp;
		slen = FMAX(slen_t, 0.1 * slen);
	}
}
