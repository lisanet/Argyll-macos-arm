
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

/* Note that all arrays are indexed from 0 */

#include "numsup.h"
#include "powell.h"
#include "varmet.h"

#undef VDEBUG					/* Varmet debug */
#undef LDEBUG					/* Line min debug */
#undef PLOTL					/* Plot line search response */

#if defined(VDEBUG)
# undef VDBG
# define VDBG(xxx) printf xxx ;
#else
# undef VDBG
# define VDBG(xxx) 
#endif

#if defined(LDEBUG)
# undef LDBG
# define LDBG(xxx) printf xxx ;
#else
# undef LDBG
# define LDBG(xxx) 
#endif


#define FMAX(A,B) ((A) > (B) ? (A) : (B))
#define EPS   1.0e-10	/*  Machine precision. */
#define XTOL (4 * EPS)	/*  X value stop value */
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
double ftol,			/* relative function value change tolleranc to stop on */
//double xtol,			/* relative cp change tolleranc to stop on */
//double gtol,			/* Gradient value tollerance to stop on */
int maxit,				/* Maximum iterations allowed */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
double (*dfunc)(void *fdata, double dp[], double tp[]),		/* Gradient function to evaluate */
void *fdata				/* Opaque data needed by function */
) {
	int iter, fails;
	double fp, sumsq, maxstep;
	double *sdir, sumsdir;			/* Search direction */
	double *dp, *lastdp;
	double **hessian;				/* Hessian matrix */
	double *hlastdp;				/* Hessian times lastdp */
	double *cpnew;					/* new cp value from linemin */
	double *dels;					/* Delta's from each step */
	double test;

	double den, fac, fad, fae;
	double sumdg;
	int i, j;

	double xtol = XTOL;			/* relative cp change tolleranc to stop on */

	double pfp = 1e38, stopth, curdel;

	sdir    = dvector(0, di-1);
	dp	    = dvector(0, di-1);
	lastdp  = dvector(0, di-1);
	hessian = dmatrix(0, di-1, 0, di-1);
	hlastdp = dvector(0, di-1);
	cpnew   = dvector(0, di-1);
	dels    = dvector(0, di-1);

	VDBG((" doing partial deriv.\n"));
	fp = (*dfunc)(fdata, dp, cp);
    if (fp == DFUNC_NRV)
		fp = (*func)(fdata, cp);

	/* Initial line direction and pde squared */
	sumsq = 0.0;
	for (i = 0; i < di ;i++) {
		sdir[i] = -dp[i];				
		sumsq += cp[i] * cp[i];
	}

	VDBG((" initial fp %f dp %s\n", fp, debPdv(di, dp)));

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
	VDBG((" maxstep %g\n", maxstep));

	/* Until we give up */
	for (fails = iter = 0; fails < di && iter < maxit; iter++) {

		/* Search in direction sdir */
		VDBG((" it %d: doing line search\n",iter));
		linesearch(di, cp, fp, dp, sdir, cpnew, &fp, maxstep, func, fdata);

		for (i = 0; i < di; i++) {
			sdir[i] = cpnew[i] - cp[i];			/* Compare search direction, */
			cp[i] = cpnew[i];					/* and the current value. */
		}

		/* Check the relative change in x, and stop if it's below */
		/* the machine precision */
		for (test = 0.0, i = 0 ; i < di; i++) {
			double tt = fabs(sdir[i]) / FMAX(fabs(cp[i]), 1.0);
			if (tt > test)
				test = tt;
		}

		if (test < xtol) {
//			VDBG((" converged because test %g < xtol %g\n",test,xtol));
			VDBG((" failed to make progres (iter %d < maxit %d): test %g < xtol %g\n",iter,maxit,test,xtol));

			/* Try a move in the dp direction */
			sumsq = 0.0;
			for (i = 0; i < di ;i++) {
				sdir[i] = -dp[i];				
				sumsq += cp[i] * cp[i];
			}

			maxstep = MAXLEN * FMAX(sqrt(sumsq), (double)di);

			VDBG((" retrying\n"));

			fails++;

			continue;
//			break;
		}
		fails = 0;
		VDBG((" not converged because test %g >= xtol %g\n",test,xtol));

		/* Check relative change in function value */
		curdel = fabs(pfp - fp);
		dels[iter % di] = curdel;

		if (iter > di) {		/* Enough to compute a moving average del */
			double avgdel;

			for (avgdel = 0.0, i = 0 ; i < di; i++)
				avgdel += dels[i] * dels[i];
			avgdel = sqrt(avgdel/(double)di);	/* Moving RMS average */

			stopth = ftol * 0.5 * (fabs(pfp) + fabs(fp) + DBL_EPSILON);

			/* If average relative change is below threshold */
			if (avgdel <= stopth) {
					VDBG(("Reached stop tollerance avgdel %g <= stopth %g\n",avgdel,stopth))
					break;
			} else {
				VDBG(("Not stopping because avgdel %g > stopth %g\n",avgdel,stopth))
			}
		}
	
		pfp = fp;

		for (i = 0; i < di; i++)		/* Save previous partial deriv. */
			lastdp[i] = dp[i];

		VDBG((" it %d: doing partial deriv\n",iter));
		(*dfunc)(fdata, dp, cp);		/* and get the new gradient. */

#ifdef NEVER		/* This doesn't seem useful */
		/* Check the largest relative gradient, and stop if it's below */
		/* our stopping tollerance */
		den = FMAX(fp, 1.0);
		for (test = 0.0, i = 0; i < di; i++) {
			double tt = fabs(dp[i]) * FMAX(fabs(cp[i]),1.0) / den;
			if (tt > test)
				test = tt;
		}

		if (test < gtol) {
			VDBG((" converged because test %g < ftol %g\n",test,gtol));
			break;
		}
		VDBG((" Not converged because test %g >= ftol %g\n",test,gtol));
#endif

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
			/* The vector that makes BFGS different from DFP: */
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

	VDBG((" Returning %g\n",fp));

	free_dvector(dels, 0, di-1);
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

#ifdef PLOTL		/* Plot the step direction */
# include "aconfig.h"
# include "plot.h"
#endif

#define ALPHA 1.0e-4 		/* Ensures sufficient decrease in function value. */
#define LXTOL 1.0e-7		/* [1e-7] Min step length termination criteria on linesearch. */

void linesearch(
int di,
double cpold[],		/* Incoming (old) current value */
double fpold,		/* Incomming current function value */
double dp[],		/* Partial derivative */
double sdir[],		/* Search direction */
double cpnew[],		/* Return (new) value */
double *pfp,		/* Return objective value */
double maxstep,
double (*func)(void *fdata, double tp[]),
void *fdata
) {
	double sum, slope;
	double slen, slen_2 = 0.0, slen_min;
	double test, fp_2 = 0.0;
	int i;

	for (sum = 0.0, i = 0; i < di; i++)
		sum += sdir[i] * sdir[i];
	sum = sqrt(sum);

	LDBG(("  fpold %f, sdir length %f\n",fpold,sum));

	if (sum > maxstep) {
		LDBG(("  sdir scaled to maxstep %f\n",maxstep));
		for (i = 0; i < di; i++)
			sdir[i] *= maxstep/sum; 		/* Scale if attempted step is too big. */
	}
	for (slope = 0.0, i = 0; i < di; i++)
		slope += dp[i] * sdir[i];

	LDBG(("  dp . sdir = slope = %f\n",slope));

	if (slope >= 0.0) {		/* Hmm. */
		warning("varmet:linesearch: slope is >= 0");
#ifdef PLOTL
		goto do_return;
#else
		return;	
#endif
	}

	/* Compute largest ratio of sdir component to max(current value,1.0) */
	for (test = 0.0, i = 0;i < di; i++) {
		double tt = fabs(sdir[i])/FMAX(fabs(cpold[i]), 1.0);
		if (tt > test)
			test = tt;
	}

	slen_min = LXTOL/test;	/* min. step length termination criteria */
	slen = 1.0;				/* Try full step */

	/* Start of iteration loop. */
	for (;;) {
		double slen_t;

		LDBG(("  top of loop: slen %f slen_min %g slen_2 %f\n",slen, slen_min, slen_2));

		/* trial point */
		for (i = 0; i < di;i++)
			cpnew[i] = cpold[i] + slen * sdir[i];

		*pfp = (*func)(fdata, cpnew);
		LDBG(("   func %f, fpold %f %s\n",*pfp, fpold, *pfp < fpold ? "BETTER" : "same/worse"));

		if (slen < slen_min) {
			LDBG(("   return slen %f fp %f because slen %g < slen_min %g (cpold)\n",slen_2, fp_2, slen, slen_min));
			for (i = 0; i < di; i++)
				cpnew[i] = cpold[i];
			*pfp = fp_2;
			slen = slen_min;

#ifdef PLOTL
			goto do_return;
#else
			return;	
#endif

		} else if (*pfp <= (fpold + ALPHA * slen * slope)) {
			LDBG(("   return slen %f fp %f because func %f <= %f = fpold %f + %f * slen %f * slope %f (cpnew)\n",slen, *pfp, *pfp, fpold + ALPHA * slen * slope, fpold, ALPHA, slen, slope));
#ifdef PLOTL
			goto do_return;
#else
			return;	
#endif

		/* Backtracking */
		} else {
			LDBG(("   slen %g >= slen_min %g\n",slen, slen_min));

			LDBG(("   backtracking:\n"));

			/* First time through model as a quadratic */
			if (slen == 1.0) {
				slen_t = -slope/(2.0 * (*pfp - fpold - slope));
				LDBG(("   1: slen_t %f = -slope %f/(2.0 * (*pfp %f - fpold %f - slope %f)\n",slen_t,-slope, *pfp, fpold, slope));

			/* 2nd and subsequent times through model as a cubic */
			} else {
				double aa, bb;
				double rhs_1, rhs_2;

				LDBG(("   2:\n"));

				/* Components of cubic solution */
				rhs_1 = *pfp - fpold - slen * slope;
				rhs_2 = fp_2 - fpold - slen_2 * slope;
				aa = (rhs_1/(slen * slen) - rhs_2/(slen_2 * slen_2))/(slen - slen_2);
				bb = (-slen_2 * rhs_1/(slen * slen)+slen * rhs_2/(slen_2 * slen_2))/(slen - slen_2);
				LDBG(("   2: rhs_1 %f rhs_2 %f aa %f bb %f)\n",rhs_1,rhs_2,aa,bb));

				/* Denominator is zero */
				if (aa == 0.0) {
					slen_t = -slope/(2.0 * bb);
					LDBG(("   3: aa == 0, slen_t %f = -slope %f/(2.0 * bb %f)\n",slen_t,-slope, bb));
				} else {
					double dd = bb * bb - 3.0 * aa * slope;
					LDBG(("   2: dd == %f)\n",dd));
					if (dd < 0.0) {
						slen_t = 0.5 * slen;
						LDBG(("   4: dd %f < 0, slen_t %f = 0.5 * slen %f)\n",dd, slen_t, slen));
					} else if (bb <= 0.0) {
						slen_t = (-bb + sqrt(dd))/(3.0 * aa);
						LDBG(("   5: bb %f <= 0.0, slen_t %f = (-bb %f + sqrt(dd %f))/(3.0 * aa %f)\n",bb, slen_t, -bb, dd, aa));
					} else {
						slen_t = -slope/(bb + sqrt(dd));
						LDBG(("   6: bb %f > 0.0, slen_t %f = -slope %f /(bb %f + sqrt(dd %f)\n",bb, slen_t, -slope, bb, dd));
					}
				}
				if (slen_t > 0.5 * slen) {
					LDBG(("   7: slen_t %f > (0.5 * slen %f = %f)\n",slen_t, slen, 0.5 * slen));
					slen_t = 0.5 * slen;
				}
			}
		}
		fp_2 = *pfp;
		slen_2 = slen;
		slen = FMAX(slen_t, 0.1 * slen);
	}

#ifdef PLOTL		/* Plot the step direction */
	do_return:;
	{
# define RES 51
		double x1, x2;
		double xx[RES];
		double y1[RES];
		int j;

		/* Plot range */
		x1 = log10(0.00001);
		x2 = log10(1.0);

		printf("Computing plot points:\n");
		for (j = 0; j < RES; j++) {
			double cp[24];
			double vv = j/(RES-1.0) * (x2 - x1) + x1;
			double vl;

			vl = pow(10.0, vv);

			for (i = 0; i < di;i++)
				cp[i] = cpold[i] + vl * sdir[i];

			xx[j] = vv;
			y1[j] = (*func)(fdata, cp);
		}
		printf(" log10 step = %f:\n",log10(slen));
		do_plot(xx, y1, NULL, NULL, RES);
	}
#endif /* PLOTL */

}

