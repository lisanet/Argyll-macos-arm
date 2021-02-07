
/*
 * Copyright 2018 Graeme Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include "numsup.h"
#include "ludecomp.h"
#include "gnewt.h"		/* Public interface definitions */

#undef DEBUG

#ifdef DEBUG
# define DBG(xx)	printf xx;
#else
# define DBG(xx)
#endif


#define TOLX 1.0e-7 	/* Convergence criterion on delx */
#define STPMX 100.0		/* Maximum step multiplier */

static void apxjac(int n, double *x, double *fvec, double **df, void *fdata,
	void (*fcn)(void *fdata, int n, double *x, double *fvec));

static int linesearch(int n, double *xold, double fold, double *delf, double *delx,
	double *x, double *fvec, double *fp, double maxstep, void *fdata,
	void (*fcn)(void *fdata, int n, double *x, double *f), int *pfit, int maxjac, int it);

#define FMAX(A, B) ((A) > (B) ? (A) : (B))

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
) {
	int i, j, it, fit, jit, *pivx, _pivx[10];
	double f, fold;				/* half magnitide squared of fvec[] */
	double *delf, _delf[10];	/* del f where f = 0.5 F.F */
	double *fvec, _fvec[10];	/* F(x) */
	double **fjac, *_fjac[11], __fjac[10 * 10];
	double *xold, _xold[10];
	double bigfx, bigx, maxstep;
	double *delx, _delx[10];	/* Full step delta x */
	double sum;
	int rv = 0;
	
#ifdef DEBUG
	double *fvec_check;
	double **fjac_check;
#endif

	DBG(("gnewt:\n"))

	fit = jit = 0;

	/* Do local vector/array allocations */
	if (n <= 10) {
		pivx = _pivx;
		if (rfvec == NULL) {
			fvec = _fvec;
		} else
			fvec = rfvec;
		_fjac[0] = __fjac;
		fjac = _fjac+1;		/* dmatrix_reset() will setup fjac */
		xold = _xold;
		delf = _delf;
		delx = _delx;
	} else {
		pivx = ivector(0, n-1);		/* LU decomp. pivod record */
		if (rfvec == NULL) {
			fvec = dvector(0, n-1);				/* Function value */
		} else
			fvec = rfvec;
		fjac = dmatrix(0, n-1, 0, n-1);		/* Jacobian matrix */
		xold = dvector(0, n-1);				/* Previous value of x[] */
		delf = dvector(0, n-1);				/* del f */
		delx = dvector(0, n-1);				/* Full step delta x */
	}
#ifdef DEBUG
	fvec_check = dvector(0, n-1);
	fjac_check = dmatrix(0, n-1, 0, n-1);
#endif

	/* Initial function value */
	fcn(fdata, n, x, fvec);
	fit++;

	DBG((" x %s\n",debPdv(n,x)))
	DBG((" fvec %s\n",debPdv(n,fvec)))

	/* Compute half magnitide squared of function value at x */
	for (sum = 0.0, i = 0; i < n; i++)
		sum += fvec[i] * fvec[i];
	f = 0.5 * sum;
	DBG((" f %f\n",f))

	/* test for initial value being a root */
	for (bigfx = 0.0, i = 0; i < n; i++) {
		double tt = fabs(fvec[i]);
		if (tt > bigfx)
			bigfx = tt;
	}
	if (bigfx < (0.01 * ftol)) {
		goto done;
	}

	/* Compute line search x maximum step size */
	for (sum = 0.0, i = 0 ; i < n; i++)
		sum += x[i] * x[i];
	maxstep = STPMX * FMAX(sqrt(sum), (double)n);
	DBG((" maxstep %f\n",maxstep))

	/* Until we are done */
	for (it = 0; fit < maxfcn && jit < maxjac; it++) {
		double rip;

		DBG((" fit %d jit %d\n",fit,jit))

		/* Compute Jacobian matrix */
		if (jac != NULL) {
			/* lu_decomp may have swapped rows - so fix it */
			dmatrix_reset(fjac, 0, n-1, 0, n-1);
			jac(fdata, n, x, fjac);				/* User function */
		} else {
			apxjac(n, x, fvec, fjac, fdata, fcn);	/* Numerical aproximation */
		}
		jit++;

#ifdef DEBUG
		copy_dmatrix(fjac_check, fjac, 0, n-1, 0, n-1);

		DBG((" fjac = \n"))
		for (i = 0; i < n; i++)
			DBG(("  %d: %s\n",i,debPdv(n, fjac[i])))
		DBG(("\n"))
#endif

		/* Compute del f for the line search. */
		for (i = 0; i < n; i++) {
			for (sum = 0.0, j = 0; j < n; j++)
				sum += fjac[j][i] * fvec[j];		/* Hmm. df/dx . f */
			delf[i] = sum;
		}

		/* Save current values of x and f to be able to monitor progres */
		for (i = 0; i < n; i++)
			xold[i] = x[i];
		fold = f;

		/* Desired delta f to make F(x) == 0 */
		for (i = 0; i < n; i++)
			delx[i] = -fvec[i];
		DBG((" -fvec %s\n",debPdv(n,delx)))

		/* Solve for delta x using Jacobian and desired delta f */
		if (lu_decomp(fjac, n, pivx, &rip)) {
			rv = 2;
			goto done;
		}
		lu_backsub(fjac, n, pivx, delx);
		DBG((" delx %s\n",debPdv(n,delx)))

#ifdef DEBUG
		matrix_vect_mult(fvec_check, n, fjac_check, n, n, delx, n);
		DBG((" check -fvec : %s\n\n",debPdv(n,fvec_check)))
#endif

		if ((rv = linesearch(n, xold, fold, delf, delx, x, fvec, &f, maxstep, fdata, fcn,
		                                                          &fit, maxfcn, it)) != 0) {
			if (rv != 1) {	/* Not run out of itterations error */
			DBG((" linesearch failed with %d\n",rv))
				goto done;
			}
		}

		DBG((" after linesearch:\n"))
		DBG((" x %s\n",debPdv(n,x)))
		DBG((" fvec %s\n",debPdv(n,fvec)))

		/* See if f() has converged */
		for (bigfx = 0.0, i = 0; i < n; i++) {
			if (fabs(fvec[i]) > bigfx)
				bigfx = fabs(fvec[i]);
		}
		DBG((" bigfx %f ftol %f\n",bigfx,ftol))
		if (bigfx < ftol) {
			goto done;
		}

		/* Could check for zero gradient problem here... */

		/* See if x[] has converged */
		for (bigx = 0.0, i = 0; i < n; i++) {
			double tt = (fabs(x[i] - xold[i]))/FMAX(fabs(x[i]), 1.0);
			if (tt > bigx)
				bigx = tt;
		}
		DBG((" bigx %f xtol %f\n",bigx,xtol))
		if (bigx < xtol)
			goto done;
	}

	rv = 1;

  done:;

	if (n > 10) {
		if (fvec != rfvec)
			free_dvector(fvec, 0, n-1);
		free_dvector(xold, 0, n-1);
		free_dvector(delx, 0, n-1);
		free_dvector(delf, 0, n-1);
		free_dmatrix(fjac, 0, n-1, 0, n-1);
		free_ivector(pivx, 0, n-1);
	}
#ifdef DEBUG
	free_dvector(fvec_check,0, n-1);
	free_dmatrix(fjac_check,0, n-1, 0, n-1);
#endif

	return rv;
}

/* - - - - - - - - */

#define ALF 1.0e-4	/* Ensures sufficient decrease in function value. */

/* Search for a step size that makes progress */
/* Return nz on error */
static int linesearch(
	int n,
	double *xold,
	double fold,
	double *delf,	/* del f */
	double *delx,	/* full step delta x[] */
	double *x,		/* in/out current x[] */
	double *fvec,	/* return fvec at x[] */
	double *fp,		/* in/out f value */
	double maxstep,	/* maximum x step */
	void *fdata,	/* Context for fcn */
	void (*fcn)(void *fdata, int n, double *x, double *f),
	int *pfit,		/* Inc number of itts */
	int maxfcn,		/* Max function its */
	int it			/* Caller iteration count */
) {
	int i;
	double f = *fp, f2; 
	double lmda1, lmda2, min_lmda;
	double sum, slope, bigx;

	DBG(("linesearch:\n"))

	/* Comute magnitude of step */
	for (sum = 0.0, i = 0; i < n; i++)
		sum += delx[i] * delx[i];
	sum = sqrt(sum);

	/* re-scale if step is too big */
	if (sum > maxstep) {
		for (i = 0; i < n; i++)
			delx[i] *= maxstep/sum;
	}

	for (slope = 0.0, i = 0; i < n; i++)
		slope += delf[i] * delx[i];
	if (slope >= 0.0) {
		DBG((" slope %f >= 0.0\n",slope))
		return 3;
	}

	bigx = 0.0;
	for (i = 0;i < n; i++) {
		double tt = fabs(delx[i])/FMAX(fabs(xold[i]), 1.0);
		if (tt > bigx)
			bigx = tt;
	}
	min_lmda = TOLX/bigx;

	/* Try full Newton step first */
	lmda1 = 1.0;

	DBG((" lmda1 %f min_lmda %f\n",lmda1, min_lmda))

	/* Top of loop */
	for (; *pfit < maxfcn; it++) {
		double tmp_lmda;

		DBG(("  lmda1 %f\n",lmda1))

		/* Take step */
		for (i = 0;i < n;i++)
			x[i] = xold[i] + lmda1 * delx[i];

		/* Compute f = 0.5 F.F at x */
		fcn(fdata, n, x, fvec);
		(*pfit)++;
		DBG(("  x %s\n",debPdv(n,x)))
		DBG(("  fvec %s\n",debPdv(n,fvec)))
		for (sum = 0.0, i = 0; i < n; i++)
			sum += fvec[i] * fvec[i];
		f = 0.5 * sum;

//if (it == 0) printf(" linesearch: At 1st full step f %f -> %f\n", *fp, f);

		/* Convergence on delx. */
		if (lmda1 < min_lmda) {

			for (i = 0; i < n; i++)
				x[i] = xold[i];
			return 0;

		} else if (f <= fold + ALF * lmda1 * slope) {
			*fp = f;
			return 0;			/* Sufficient function decrease */

		} else {				/* Backtrack. */
			if (lmda1 == 1.0)	/* First time */
				tmp_lmda = -slope/(2.0 * (f - fold-slope));

			else {			/* Subsequent backtracks */
				double c, d, e;
				double a, b, rhs1, rhs2;

				rhs1 = f - fold - slope * lmda1;
				rhs2 = f2 - fold - slope * lmda2;
				c = rhs1/(lmda1 * lmda1);
				d = rhs2/(lmda2 * lmda2);
				e = lmda1 - lmda2;
				a = (c - d)/e;
				b = (-lmda2 * c + lmda1 * d)/e;
				if (a == 0.0)
					tmp_lmda = -slope/(2.0 * b);
				else {
					double disc = b * b - 3.0 * a * slope;

					if (disc < 0.0)
						tmp_lmda = 0.5 * lmda1;
					else if (b <= 0.0)
						tmp_lmda = (-b + sqrt(disc))/(3.0 * a);
					else
						tmp_lmda = -slope/(b + sqrt(disc));
				}
				if (tmp_lmda > 0.5 * lmda1)
					tmp_lmda = 0.5 * lmda1;
			}
		}
		lmda2 = lmda1;
		lmda1 = FMAX(tmp_lmda, lmda1 * 0.1);
		f2 = f;
	}

	*fp = f;

	return 1;
}

/* - - - - - - - - */

/* Compute forward difference as aprox. Jacobian matrix */

#define JEPS 1.0e-8	/* Aprox. sqrt of machine precision */

static void apxjac(
	int n,					/* Dimensions */
	double *x,				/* Location x to compute Jacobian */
	double *fvec,			/* Function value at x */
	double **df,			/* Return Jacobian */
	void *fdata,			/* fcn() context */
	void (*fcn)(void *fdata, int n, double *x, double *fvec)
) {
	int i, j;
	double h, temp, *f, _f[10];

	if (n <= 10)
		f = _f;
	else
		f = dvector(0, n);

	for (j = 0; j < n; j++) {
		temp = x[j];

		h = JEPS * fabs(temp);
		if (h == 0.0)
			h = JEPS;
		x[j] = temp + h;		/* Add delta */
		h = x[j] - temp;		/* Actual delta with fp precision limits */

		fcn(fdata, n, x, f);

		x[j] = temp;			/* Restore value */

		for (i = 0; i < n; i++)
			df[i][j] = (f[i] - fvec[i])/h;
	}

	if (f != _f)
		free_dvector(f, 0, n-1);
}


