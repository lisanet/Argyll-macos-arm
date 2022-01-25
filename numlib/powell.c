
/* Multi-dimentional minizer using Powell or Conjugate Gradient methods */
/* This is good for smoother, well behaved functions. */

/* Code is an original expression of the algorithms decsribed in */
/* "Numerical Recipes in C", by W.H.Press, B.P.Flannery, */
/* S.A.Teukolsky & W.T.Vetterling. */

/*
 * Copyright 2000, 2006 Graeme W. Gill
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

/*

	Idea for improving progress accounting:

	count number of itterations already done (pitter)
	estimate number yet needed (fitter)
	progress = pitter/(pitter + fitter)

	Number yet needed estimated by progress of retval delta
	againsth threshold target. 

	ie  fitters = (lastdel - curdel)/(curdel - stopth)

*/

/* Note that all arrays are indexed from 0 */

#include "numsup.h"
#include "powell.h"

#undef SLOPE_SANITY_CHECK		/* [und] expermental */
#undef ABSTOL					/* [und] Make tollerance absolute */
#undef USE_LINMIND				/* [und] Use limind for conjgrad (typically slower) */

								/* Some debugging printfs (not comprehensive): */
#undef PDEBUG					/* Powell debug */
#undef CDEBUG					/* Conjgrad debug */
#undef LDEBUG					/* Line min debug */

#if defined(PDEBUG)
# undef PDBG
# define PDBG(xxx) printf xxx ;
#else
# undef PDBG
# define PDBG(xxx) 
#endif

#if defined(CDEBUG)
# undef CDBG
# define CDBG(xxx) printf xxx ;
#else
# undef CDBG
# define CDBG(xxx) 
#endif

#if defined(LDEBUG)
# undef LDBG
# define LDBG(xxx) printf xxx ;
#else
# undef LDBG
# define LDBG(xxx) 
#endif

/* -------------------------------------- */
/* Standard interface for powell function */
/* return 0 on sucess, 1 on failure due to excessive itterions */
/* Result will be in cp */
int powell(
double *rv,				/* If not NULL, return the residual error */
int di,					/* Dimentionality */
double cp[],			/* Initial starting point */
double s[],				/* Size of initial search area */
#ifdef ABSTOL
double ftol,			/* Absolute tollerance of error change to stop on */
#else
double ftol,			/* Relative tollerance of error change to stop on */
#endif
int maxit,				/* Maximum iterations allowed */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
void *fdata,			/* Opaque data needed by function */
void (*prog)(void *pdata, int perc),		/* Optional progress percentage callback */
void *pdata				/* Opaque data needed by prog() */
) {
	int i;
	double **dmtx, *_dmtx[10], __dmtx[10 * 10] = { 0.0 };		/* Direction vectors */
	double *spt, _spt[10];			/* Sarting point before exploring all the directions */
	double *xpt, _xpt[10];			/* Extrapolated point */
	double *svec, _svec[10];		/* Search vector */
	int    iter;
	double retv; 			/* Returned function value at p */
	double stopth;			/* Current stop threshold */
	double startdel = -1.0;	/* Initial change in function value */
	double curdel;			/* Current change in function value */
	int pc = 0;				/* Percentage complete */

	if (di <= 10) {
		int j;
		for (j = i = 0; i < di; i++, j += di)
			_dmtx[i] = __dmtx + j;
		dmtx = _dmtx;
		spt  = _spt;
		xpt  = _xpt;
		svec = _svec;
	} else {
		dmtx = dmatrixz(0, di-1, 0, di-1);	/* Zero filled */
		spt  = dvector(0, di-1);
		xpt  = dvector(0, di-1);
		svec = dvector(0, di-1);
	}

	/* Create initial direction matrix by */
	/* placing search start on diagonal */
	for (i = 0; i < di; i++)
		dmtx[i][i] = s[i];
	
	/* Save the starting point */
	for (i = 0; i < di; i++)
		spt[i] = cp[i];

	if (prog != NULL)		/* Report initial progress */
		prog(pdata, pc);

	/* Initial function evaluation */
	retv = (*func)(fdata, cp);

//printf("~1 ### initial retv = %f\n",retv);
	/* Itterate until we converge on a solution, or give up. */
	for (iter = 1; iter < maxit; iter++) {
		int j;
		double lretv;			/* Last function return value */
		int    ibig = 0;		/* Index of biggest delta */
		double del = 0.0;		/* Biggest function value decrease */
		double pretv;			/* Previous function return value */

		pretv = retv;			/* Save return value at top of itteration */

		/* Loop over all directions in the set */
		for (i = 0; i < di; i++) {

			PDBG(("Looping over direction %d\n",i))

			for (j = 0; j < di; j++)	/* Extract this direction to make search vector */
				svec[j] = dmtx[j][i];

//printf("~1 ### chosen dir = %f %f\n", svec[0],svec[1]);
			/* Minimize in that direction */
			lretv = retv;
			retv = linmin(cp, svec, di, ftol, func, fdata);

			/* Record bigest function decrease, and direction it occured on */
			if (fabs(lretv - retv) > del) {
				del = fabs(lretv - retv);
				ibig = i;
			}
		}

//printf("~1 ### biggest change was dir %d by %f\n", ibig, del);

#ifdef ABSTOL
		stopth = ftol;				/* Absolute tollerance */
#else
		stopth = ftol * 0.5 * (fabs(pretv) + fabs(retv) + DBL_EPSILON);
#endif
		curdel = fabs(pretv - retv);
		if (startdel < 0.0) {
			startdel = curdel;
		} else {
			int tt;
			tt = (int)(100.0 * pow((log(curdel) - log(startdel))/(log(stopth) - log(startdel)), 4.0) + 0.5);
			if (tt > pc && tt < 100) {
				pc = tt;
				if (prog != NULL)		/* Report initial progress */
					prog(pdata, pc);
			}

		}

		/* If we have made at least one pass through all directions and */
		/* reached a suitable tollerance, then finish */
		if (iter > 1 && curdel <= stopth) {
//printf("~1 ### stopping on itter %d because curdel %f <= stopth %f\n",iter, curdel,stopth);
			PDBG(("Reached stop tollerance because curdel %f <= stopth %f\n",curdel,stopth))
			break;
		}
		PDBG(("Not stopping because curdel %f > stopth %f\n",curdel,stopth))

//printf("~1 ### recomputing direction\n");
		/* Compute overall direction minimization moved in */
		for (i = 0; i < di; i++) {
			svec[i] = cp[i] - spt[i];	/* Average direction moved after minimization round */
			xpt[i]  = cp[i] + svec[i];	/* Extrapolated point after round of minimization */
			spt[i]  = cp[i];			/* New start point for next round */
		}
//printf("~1 ### new dir = %f %f\n", svec[0],svec[1]);

		/* Function value at extrapolated point in overall direction moved in */
		lretv = (*func)(fdata, xpt);

		if (lretv < pretv) {			/* If extrapolation is an improvement */
			double t, t1, t2;

//printf("~1 ### extrap is improvement\n");
			t1 = pretv - retv - del;
			t2 = pretv - lretv;
			t = 2.0 * (pretv -2.0 * retv + lretv) * t1 * t1 - del * t2 * t2;
			if (t < 0.0) {
//printf("~1 ### move to min in new direction\n");
				/* Move to the minimum of the new direction */
				retv = linmin(cp, svec, di, ftol, func, fdata);

				for (i = 0; i < di; i++) { 		/* Save the new direction */
					dmtx[i][ibig] = svec[i];	/* by replacing best previous */
				}
			}
		}
	}		/* Continue itterating */

//printf("~1 iters = %d\n",iter);

	/* Free up all the temporary vectors and matrix */
	if (di > 10) {
		free_dvector(svec, 0, di-1);
		free_dvector(xpt, 0, di-1);
		free_dvector(spt, 0, di-1);
		free_dmatrix(dmtx, 0, di-1, 0, di-1);
	}

	if (prog != NULL)		/* Report final progress */
		prog(pdata, 100);

	if (rv != NULL)
		*rv = retv;

	if (iter < maxit)
		return 0;

	PDBG(("powell: returning 1 due to excessive itterations\n"))
	return 1;		/* Failed due to execessive itterations */
}

/* - - - - - - - - - - - - - - - - - */

#define POWELL_GOLD 1.618034
#define POWELL_CGOLD 0.3819660
#define POWELL_MAXIT 100

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
void *fdata)		/* Opaque data for func() */
{
	int i;
	double ax, xx, bx;	/* Search vector multipliers */
	double af, xf, bf;	/* Function values at those points */
	double *xt, _xt[10];	/* Trial point */

	if (di <= 10)
		xt = _xt;
	else
		xt = dvector(0, di-1);			/* Vector for trial point */

	/* -------------------------- */
	/* First bracket the solution */

	LDBG((" linmin: Bracketing solution\n"))

	/* The line is measured as startpoint + offset * search vector. */
	/* (Search isn't symetric, but it seems to depend on cp being */
	/* best current solution ?) */
	ax = 0.0;
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + ax * xi[i];
	af = (*func)(fdata, xt);

	/* xx being vector offset 0.618 */
	xx =  1.0/POWELL_GOLD;
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + xx * xi[i];
	xf = (*func)(fdata, xt);

	LDBG((" linmin: Initial points a:%f:%f -> b:%f:%f\n",ax,af,xx,xf))

	/* Fix it so that we are decreasing from point a -> x */
	if (xf > af) {
		double tt;
		tt = ax; ax = xx; xx = tt;
		tt = af; af = xf; xf = tt;
	}
	LDBG((" linmin: Ordered Initial points a:%f:%f -> b:%f:%f\n",ax,af,xx,xf))

	bx = xx + POWELL_GOLD * (xx-ax);	/* Guess b beyond a -> x */
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + bx * xi[i];
	bf = (*func)(fdata, xt);

	LDBG((" linmin: Initial bracket a:%f:%f x:%f:%f b:%f:%f\n",ax,af,xx,xf,bx,bf))

#ifdef SLOPE_SANITY_CHECK
	/* If we're not seeing a slope indicitive of progress */
	/* of order ftol, give up straight away */
	if (2000.0 * fabs(xf - bf) <= ftol * (fabs(xf) + fabs(bf))
	 && 2000.0 * fabs(af - xf) <= ftol * (fabs(af) + fabs(xf))) {
		LDBG((" linmin: giving up because slope is too shallow\n"))
		if (xt != _xt)
			free_dvector(xt,0,di-1);

		if (bf < xf) {
			xf = bf;
			xx = bx;
		}
		goto done;
	}
#endif /* SLOPE_SANITY_CHECK */

	/* While not bracketed */
	while (xf > bf) {
		double ulim, ux, uf;
		double tt, r, q;

		LDBG((" linmin: Not bracketed because xf %f > bf %f\n",xf, bf))
		LDBG(("        ax = %f, xx = %f, bx = %f\n",ax,xx,bx))

		/* Compute ux by parabolic interpolation from a, x & b */
		q = (xx - bx) * (xf - af);
		r = (xx - ax) * (xf - bf);
		tt = q - r;
		if (tt >= 0.0 && tt < 1e-20)				/* If +ve too small */
			tt = 1e-20;
		else if (tt <= 0.0 && tt > -1e-20)		/* If -ve too small */
			tt = -1e-20;
		ux = xx - ((xx - bx) * q - (xx - ax) * r) / (2.0 * tt);
		ulim = xx + 100.0 * (bx - xx);			/* Extrapolation limit */

//printf("~1 ux = %f, ulim = %f\n",ux,ulim);
		if ((xx - ux) * (ux - bx) > 0.0) {		/* u is between x and b */

			for (i = 0; i < di; i++)			/* Evaluate u */
				xt[i] = cp[i] + ux * xi[i];
			uf = (*func)(fdata, xt);

//printf("~1 u is between x and b, uf = %f\n",uf);

			if (uf < bf) {						/* Minimum is between x and b */
//printf("~1 min is between x and b\n");
				ax = xx; af = xf;
				xx = ux; xf = uf;
				break;
			} else if (uf > xf) {				/* Minimum is between a and u */
//printf("~1 min is between a and u\n");
				bx = ux; bf = uf;
				break;
			}

			/* Parabolic fit didn't work, look further out in direction of b */
			ux = bx + POWELL_GOLD * (bx-xx);
//printf("~1 parabolic fit didn't work,look further in direction of b (%f)\n",ux);

		} else if ((bx - ux) * (ux - ulim) > 0.0) {	/* u is between b and limit */
			for (i = 0; i < di; i++)			/* Evaluate u */
				xt[i] = cp[i] + ux * xi[i];
			uf = (*func)(fdata, xt);

//printf("~1 u is between b and limit uf = %f\n",uf);
			if (uf > bf) {						/* Minimum is between x and u */
//printf("~1 min is between x and uf\n");
				ax = xx; af = xf;
				xx = bx; xf = bf;
				bx = ux; bf = uf;
				break;
			}
			xx = bx; xf = bf;					/* Continue looking */
			bx = ux; bf = uf;
			ux = bx + POWELL_GOLD * (bx - xx);	/* Test beyond b */
//printf("~1 continue looking beyond b (%f)\n",ux);

		} else if ((ux - ulim) * (ulim - bx) >= 0.0) {	/* u is beyond limit */
			ux = ulim;
//printf("~1 use limit\n");
		} else {							/* u is to left side of x ? */
			ux = bx + POWELL_GOLD * (bx-xx);
//printf("~1 look gold beyond b (%f)\n",ux);
		}
		/* Evaluate u, and move into place at b */
		for (i = 0; i < di; i++)
			xt[i] = cp[i] + ux * xi[i];
		uf = (*func)(fdata, xt);
//printf("~1 lookup ux %f value uf = %f\n",ux,uf);
		ax = xx; af = xf;
		xx = bx; xf = bf;
		bx = ux; bf = uf;
//printf("~1 move along to the right (a<-x, x<-b, b-<u)\n");
	}
	LDBG((" linmin: Got bracket a:%f:%f x:%f:%f b:%f:%f\n",ax,af,xx,xf,bx,bf))
	/* Got bracketed minimum between a -> x -> b */
//printf("~1 got bracketed minimum at %f (%f), %f (%f), %f (%f)\n",ax,af,xx,xf,bx,bf);

	/* --------------------------------------- */
	/* Now use brent minimiser bewteen a and b */
	{
		/* a and b bracket solution */
		/* x is best function value so far */
		/* w is second best function value so far */
		/* v is previous second best, or third best */
		/* u is most recently tested point */
		double wx, vx, ux;			/* Search vector multipliers */
		double wf, vf = 0.0, uf;	/* Function values at those points */
		int iter;
		double de = 0.0;	/* Distance moved on previous step */
		double e = 0.0;		/* Distance moved on 2nd previous step */

		/* Make sure a and b are in ascending order */
		if (ax > bx) {
			double tt;
			tt = ax; ax = bx; bx = tt;
			tt = af; af = bf; bf = tt;
		}

		wx = vx = xx;	/* Initial values of other center points */
		wf = xf = xf;

		for (iter = 1; iter <= POWELL_MAXIT; iter++) {
			double mx = 0.5 * (ax + bx);		/* m is center of bracket values */
#ifdef ABSTOL
			double tol1 = ftol;			/* Absolute tollerance */
#else
			double tol1 = ftol * fabs(xx) + 1e-10;
#endif
			double tol2 = 2.0 * tol1;

			LDBG((" linmin it %d: Got bracket a:%f:%f x:%f:%f b:%f:%f\n",iter,ax,af,xx,xf,bx,bf))

			/* See if we're done */
//printf("~1 linmin check %f <= %f\n",fabs(xx - mx), tol2 - 0.5 * (bx - ax));
			if (fabs(xx - mx) <= (tol2 - 0.5 * (bx - ax))) {
				LDBG((" linmin: We're done because %e <= %e\n",fabs(xx - mx), tol2 - 0.5 * (bx - ax)))
				break;
			}

			LDBG((" linmin: e %e tol2 %e\n",e,tol1))

			if (fabs(e) > tol1) {			/* Do a trial parabolic fit */
				double te, p, q, r;
				r = (xx-wx) * (xf-vf);
				q = (xx-vx) * (xf-wf);
				p = (xx-vx) * q - (xx-wx) * r;
				q = 2.0 * (q - r);
				if (q > 0.0)
					p = -p;
				else
					q = -q;
				te = e;				/* Save previous e value */
				e = de;				/* Previous steps distance moved */

				LDBG((" linmin: Trial parabolic fit\n" ))

				if (fabs(p) >= fabs(0.5 * q * te) || p <= q * (ax-xx) || p >= q * (bx-xx)) {
					/* Give up on the parabolic fit, and use the golden section search */
					e = ((xx >= mx) ? ax-xx : bx-xx);	/* Override previous distance moved */
					de = POWELL_CGOLD * e;
					LDBG((" linmin: Moving to golden section search\n" ))
				} else {	/* Use parabolic fit */
					de = p/q;			/* Change in xb */
					ux = xx + de;		/* Trial point according to parabolic fit */
					if ((ux - ax) < tol2 || (bx - ux) < tol2) {
						if ((mx - xx) > 0.0)	/* Don't use parabolic, use tol1 */
							de = tol1;			/* tol1 is +ve */
						else
							de = -tol1;
					}
					LDBG((" linmin: Using parabolic fit\n" ))
				}
			} else {	/* Keep using the golden section search */
				e = ((xx >= mx) ? ax-xx : bx-xx);	/* Override previous distance moved */
				de = POWELL_CGOLD * e;
				LDBG((" linmin: Continuing golden section search\n" ))
			}

			if (fabs(de) >= tol1) {		/* If de moves as much as tol1 would */
				ux = xx + de;			/* use it */
				LDBG((" linmin: ux = %f = xx %f + de %f\n",ux,xx,de))
			} else {					/* else move by tol1 in direction de */
				if (de > 0.0) {
					ux = xx + tol1;
					LDBG((" linmin: ux = %f = xx %f + tol1 %e\n",ux,xx,tol1))
				} else {
					ux = xx - tol1;
					LDBG((" linmin: ux = %f = xx %f - tol1 %f\n",ux,xx,tol1))
				}
			}

			/* Evaluate function */
			for (i = 0; i < di; i++)
				xt[i] = cp[i] + ux * xi[i];
			uf = (*func)(fdata, xt);

			if (uf <= xf) {					/* Found new best solution */
				LDBG((" linmin: found new best solution at %f val %f\n",ux,uf))
				if (ux >= xx) {	
					ax = xx; af = xf;		/* New lower bracket */
				} else {
					bx = xx; bf = xf;		/* New upper bracket */
				}
				vx = wx; vf = wf;			/* New previous 2nd best solution */
				wx = xx; wf = xf;			/* New 2nd best solution from previous best */
				xx = ux; xf = uf;			/* New best solution from latest */
			} else {						/* Found a worse solution */
				LDBG((" linmin: found new worse solution at %f val %f\n",ux,uf))
				LDBG((" linmin:             current best at %f val %f\n",xx,xf))
				if (ux < xx) {
					ax = ux; af = uf;		/* New lower bracket */
				} else {
					bx = ux; bf = uf;		/* New upper bracket */
				}
				if (uf <= wf || wx == xx) {	/* New 2nd best solution, or equal best */
					vx = wx; vf = wf;		/* New previous 2nd best solution */
					wx = ux; wf = uf;		/* New 2nd best from latest */
				} else if (uf <= vf || vx == xx || vx == wx) {	/* New 3rd best, or equal 1st & 2nd */
					vx = ux; vf = uf;		/* New previous 2nd best from latest */
				}
			}
		}
		/* !!! should do something if iter > POWELL_MAXIT !!!! */
		/* Solution is at xx, xf */
	}

  done:;

	/* Compute solution vector at xx */
	LDBG((" linmin: computing soln from best at %f val %f\n",xx,xf))
	for (i = 0; i < di; i++) 
		cp[i] += xx * xi[i];

	if (xt != _xt)
		free_dvector(xt,0,di-1);
//printf("~~~ line minimizer returning %e\n",xf);
	return xf;
}

#undef POWELL_GOLD
#undef POWELL_CGOLD
#undef POWELL_MAXIT

/* -------------------------------------- */
/* Conjugate Gradient optimiser using partial derivatives. */
/* return 0 on sucess, 1 on failure due to excessive itterions */
/* Result will be in cp */
/* Note that we could use gradient in line minimiser, */
/* but this seems to be slower, so we don't use it. */
int conjgrad(
double *rv,				/* If not NULL, return the residual error */
int di,					/* Dimentionality */
double cp[],			/* Initial starting point and return value */
double s[],				/* Size of initial search area */
#ifdef ABSTOL
double ftol,			/* Absolute tollerance of error change to stop on */
#else
double ftol,			/* Relative tollerance of error change to stop on */
#endif
int maxit,				/* Maximum iterations allowed */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
double (*dfunc)(void *fdata, double dp[], double tp[]),		/* Gradient & function to evaluate */
						/* dfunc() should return DFUNC_NRV if it doesn't return function value */
void *fdata,			/* Opaque data needed by function */
void (*prog)(void *pdata, int perc),		/* Optional progress percentage callback */
void *pdata				/* Opaque data needed by prog() */
) {
	int i, iter;
	double *svec, _svec[10];	/* Search vector */
	double *ssvec, _ssvec[10];	/* s[] scaled search vector */
	double *gvec, _gvec[10];	/* G direction vector */
	double *hvec, _hvec[10];	/* H direction vector */
	double retv; 			/* Returned function value at p */
	double stopth;			/* Current stop threshold */
	double startdel = -1.0;	/* Initial change in function value */
	double curdel;			/* Current change in function value */
	double brat;			/* svec to s[] ratio */
	double svec_sca;		/* svec scale factor */
	int pc = 0;				/* Percentage complete */

	if (di <= 10) {
		svec = _svec;
		ssvec = _ssvec;
		gvec  = _gvec;
		hvec  = _hvec;
	} else {
		svec = dvector(0, di-1);
		ssvec = dvector(0, di-1);
		gvec  = dvector(0, di-1);
		hvec  = dvector(0, di-1);
	}

	CDBG(("conjgrad with di %d\n", di))
	CDBG((" cp = %s\n", debPdv(di,cp)))
	CDBG(("  s = %s\n", debPdv(di,s)))

	if (prog != NULL)		/* Report initial progress */
		prog(pdata, pc);

	/* Initial function and gradient evaluation */
	CDBG((" calling dfunc\n"))
	retv = (*dfunc)(fdata, svec, cp);

	CDBG((" returned %e and d %s\n",retv,debPdv(di,svec)))

	if (retv == DFUNC_NRV) { 
		CDBG((" calling func\n"))
		retv = (*func)(fdata, cp);
	} else {
		CDBG((" dfunc returns func value %f\n",retv))
	}

	/* svec[] seems to be large after this. Compute scaled version that */
	/* has maximum of s[] so that line search is guided by the search radius. */
	for (brat = 0.0, i = 0; i < di; i++) {
		double rat = fabs(svec[i]) / fabs(s[i]);
		if (rat > brat)
			brat = rat;
	}

	svec_sca = 1.0;
	if (brat > DBL_EPSILON) 
		svec_sca /= brat;

	CDBG((" svec_sca = %f\n",svec_sca))

	/* Initial vector setup */
	for (i = 0; i < di; i++) {
		gvec[i] = -svec[i];
		svec[i] = hvec[i] = gvec[i];	/* Set G & H to -ve gradient */
		ssvec[i] = svec[i] * svec_sca;	/* Scale the search vector to s[] size */
	}

	CDBG((" initial dir = %s\n", debPdv(di, ssvec)));
	CDBG((" initial retv = %f\n",retv));

	/* Itterate untill we converge on a solution, or give up. */
	for (iter = 1; iter < maxit; iter++) {
		double gamden, gamnum, gam;
		double pretv;			/* Previous function return value */

		CDBG(("conjrad it %d: about to do linmind\n",iter))
		pretv = retv;
#ifdef USE_LINMIND
		retv = linmind(cp, ssvec, di, ftol, func, dfunc, fdata);
#else
		retv = linmin(cp, ssvec, di, ftol, func, fdata);
#endif

#ifdef ABSTOL
		stopth = ftol;				/* Absolute tollerance */
#else
		stopth = ftol * 0.5 * (fabs(pretv) + fabs(retv) + DBL_EPSILON);
#endif
		curdel = fabs(pretv - retv);
		CDBG((" this retv = %f, pretv = %f, curdel = %f\n",retv,pretv,curdel));
		if (startdel < 0.0) {
			startdel = curdel;
		} else if (prog != NULL) {	/* Update percentage */
			int tt;
			tt = (int)(100.0 * pow((log(curdel) - log(startdel))/(log(stopth) - log(startdel)), 4.0) + 0.5);
			if (tt > pc && tt < 100) {
				pc = tt;
				prog(pdata, pc); /* Report initial progress */
			}
		}

		/* If we have had at least one change of direction and */
		/* reached a suitable tollerance, then finish */
		if (iter > 1 && curdel <= stopth) {
			CDBG((" stopping on itter %d because curdel %f <= stopth %f\n",iter, curdel,stopth));
			break;
		}
		CDBG((" not stopping on itter %d because curdel %f > stopth %f\n",iter, curdel,stopth));

		CDBG(("conjrad: recomputing direction\n"))
		(*dfunc)(fdata, svec, cp);		/* (Don't use retv as it wrecks stop test) */
		CDBG((" pderiv = %s\n", debPdv(di, svec)));

		/* Compute gamma */
		for (gamnum = gamden = 0.0, i = 0; i < di; i++) {
			gamden += gvec[i] * gvec[i];
//			gamnum += svec[i] * svec[i];						/* Flecher-Reeves */
			gamnum += svec[i] * (gvec[i] + svec[i]);			/* Polak-Ribiere */
		}

		CDBG((" gamnum = %f, gamden = %f\n", gamnum,gamden));
		if (fabs(gamden) < DBL_EPSILON) {		/* Gradient is exactly zero */
			CDBG(("conjrad: gradient is exactly zero\n"))
			break;
		}

		gam = gamnum/gamden;
		CDBG(("conjrad: gamma = %f = %f/%f\n",gam,gamnum,gamden));
		CDBG((" gvec[] = %s, hvec = %s\n", debPdv(di,gvec),debPdv(di,hvec)));

		/* Adjust seach direction */
		for (i = 0; i < di; i++) {
			gvec[i] = -svec[i];
			svec[i] = hvec[i] = gvec[i] + gam * hvec[i];
		}

		/* svec[] seems to be large after this. Compute scaled version that */
		/* has maximum of s[] so that line search is guided by the search radius. */
		for (brat = 0.0, i = 0; i < di; i++) {
			double rat = fabs(svec[i]) / fabs(s[i]);
			if (rat > brat)
				brat = rat;
		}
		svec_sca = 1.0/brat;
		for (i = 0; i < di; i++)
			ssvec[i] = svec[i] * svec_sca;

		CDBG((" ssvec = %s\n", debPdv(di,ssvec)));
	}
	/* Free up all the temporary vectors and matrix */
	if (di > 10) {
		free_dvector(hvec, 0, di-1);
		free_dvector(gvec, 0, di-1);
		free_dvector(ssvec, 0, di-1);
		free_dvector(svec, 0, di-1);
	}

	if (prog != NULL)		/* Report final progress */
		prog(pdata, 100);

	if (rv != NULL)
		*rv = retv;

	CDBG((" conjgrad returning = %f\n", retv));

	if (iter < maxit)
		return 0;

	return 1;		/* Failed due to execessive itterations */
}

#define POWELL_GOLD 1.618034
#define POWELL_MAXIT 100

/* Line bracketing and minimisation routine using derivatives */
/* This is not used, because it typically makes it slower */
/* - it may take less itterations, but each itteration uses */
/* a func() and dfunc() call, at least doubling itter overhead. */
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
double (*func)(void *fdata, double tp[]),				/* Error function to evaluate */
double (*dfunc)(void *fdata, double dp[], double tp[]),	/* Gradient function to evaluate */
					/* dfunc() should return DFUNC_NRV if it doesn't return function value */
void *fdata)		/* Opaque data for func() */
{
	int i;
	double ax, xx, bx;	/* Search vector multipliers */
	double af, xf, bf;	/* Function values at those points */
	double *xt, _xt[10];	/* Trial point */
	double *df, _df[10];	/* Derivative vector */

	if (di <= 10) {
		xt = _xt;
		df = _df;
	} else {
		xt = dvector(0, di-1);			/* Vector for trial point */
		df = dvector(0, di-1);			/* Vector for trial point */
	}

	/* -------------------------- */
	/* First bracket the solution */

	LDBG((" linmind: Bracketing solution\n"))

	/* The line is measured as startpoint + offset * search vector. */
	/* (Search isn't symetric, but it seems to depend on cp being */
	/* best current solution ?) */
	ax = 0.0;
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + ax * xi[i];
	af = (*func)(fdata, xt);

	/* xx being vector offset 0.618 */
	xx =  1.0/POWELL_GOLD;
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + xx * xi[i];
	xf = (*func)(fdata, xt);

	LDBG((" linmind: Initial points a:%f:%f -> b:%f:%f\n",ax,af,xx,xf))

	/* Fix it so that we are decreasing from point a -> x */
	if (xf > af) {
		double tt;
		tt = ax; ax = xx; xx = tt;
		tt = af; af = xf; xf = tt;
	}
	LDBG((" linmind: Ordered Initial points a:%f:%f -> b:%f:%f\n",ax,af,xx,xf))

	bx = xx + POWELL_GOLD * (xx-ax);	/* Guess b beyond a -> x */
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + bx * xi[i];
	bf = (*func)(fdata, xt);

	LDBG((" linmind: Initial bracket a:%f:%f x:%f:%f b:%f:%f\n",ax,af,xx,xf,bx,bf))

#ifdef SLOPE_SANITY_CHECK
	/* If we're not seeing a slope indicitive of progress */
	/* of order ftol, give up straight away */
	if (2000.0 * fabs(xf - bf) <= ftol * (fabs(xf) + fabs(bf))
	 && 2000.0 * fabs(af - xf) <= ftol * (fabs(af) + fabs(xf))) {
		LDBG((" linmind: giving up because slope is too shallow\n"))
		if (di > 10) {
			free_dvector(df, 0, di-1);
			free_dvector(xt, 0, di-1);
		}

		if (bf < xf) {
			xf = bf;
			xx = bx;
		}
		goto done;
	}
#endif /* SLOPE_SANITY_CHECK */

	/* While not bracketed */
	while (xf > bf) {
		double ulim, ux, uf;
		double tt, r, q;

		LDBG((" linmind: Not bracketed because xf %f > bf %f\n",xf, bf))
		LDBG(("        ax = %f, xx = %f, bx = %f\n",ax,xx,bx))

		/* Compute ux by parabolic interpolation from a, x & b */
		q = (xx - bx) * (xf - af);
		r = (xx - ax) * (xf - bf);
		tt = q - r;
		if (tt >= 0.0 && tt < 1e-20)				/* If +ve too small */
			tt = 1e-20;
		else if (tt <= 0.0 && tt > -1e-20)		/* If -ve too small */
			tt = -1e-20;
		ux = xx - ((xx - bx) * q - (xx - ax) * r) / (2.0 * tt);
		ulim = xx + 100.0 * (bx - xx);			/* Extrapolation limit */

//printf("~1 ux = %f, ulim = %f\n",ux,ulim);
		if ((xx - ux) * (ux - bx) > 0.0) {		/* u is between x and b */

			for (i = 0; i < di; i++)			/* Evaluate u */
				xt[i] = cp[i] + ux * xi[i];
			uf = (*func)(fdata, xt);

//printf("~1 u is between x and b, uf = %f\n",uf);

			if (uf < bf) {						/* Minimum is between x and b */
//printf("~1 min is between x and b\n");
				ax = xx; af = xf;
				xx = ux; xf = uf;
				break;
			} else if (uf > xf) {				/* Minimum is between a and u */
//printf("~1 min is between a and u\n");
				bx = ux; bf = uf;
				break;
			}

			/* Parabolic fit didn't work, look further out in direction of b */
			ux = bx + POWELL_GOLD * (bx-xx);
//printf("~1 parabolic fit didn't work,look further in direction of b (%f)\n",ux);

		} else if ((bx - ux) * (ux - ulim) > 0.0) {	/* u is between b and limit */
			for (i = 0; i < di; i++)			/* Evaluate u */
				xt[i] = cp[i] + ux * xi[i];
			uf = (*func)(fdata, xt);

//printf("~1 u is between b and limit uf = %f\n",uf);
			if (uf > bf) {						/* Minimum is between x and u */
//printf("~1 min is between x and uf\n");
				ax = xx; af = xf;
				xx = bx; xf = bf;
				bx = ux; bf = uf;
				break;
			}
			xx = bx; xf = bf;					/* Continue looking */
			bx = ux; bf = uf;
			ux = bx + POWELL_GOLD * (bx - xx);	/* Test beyond b */
//printf("~1 continue looking beyond b (%f)\n",ux);

		} else if ((ux - ulim) * (ulim - bx) >= 0.0) {	/* u is beyond limit */
			ux = ulim;
//printf("~1 use limit\n");
		} else {							/* u is to left side of x ? */
			ux = bx + POWELL_GOLD * (bx-xx);
//printf("~1 look gold beyond b (%f)\n",ux);
		}
		/* Evaluate u, and move into place at b */
		for (i = 0; i < di; i++)
			xt[i] = cp[i] + ux * xi[i];
		uf = (*func)(fdata, xt);
//printf("~1 lookup ux %f value uf = %f\n",ux,uf);
		ax = xx; af = xf;
		xx = bx; xf = bf;
		bx = ux; bf = uf;
//printf("~1 move along to the right (a<-x, x<-b, b-<u)\n");
	}
	LDBG((" linmind: Got bracket a:%f:%f x:%f:%f b:%f:%f\n",ax,af,xx,xf,bx,bf))
	/* Got bracketed minimum between a -> x -> b */
//printf("~1 got bracketed minimum at %f (%f), %f (%f), %f (%f)\n",ax,af,xx,xf,bx,bf);

	/* --------------------------------------- */
	/* Now use brent minimiser bewteen a and b */
	{
		/* a and b bracket solution */
		/* x is best function value so far */
		/* w is second best function value so far */
		/* v is previous second best, or third best */
		/* u is most recently tested point */

		double wx, vx, ux;			/* Search vector multipliers */
		double wf, vf = 0.0, uf;	/* Function values at those points */
		double xd, wd, vd, ud;		/* Derivative values at those points */
		int iter;
		double de = 0.0;	/* Distance moved on previous step */
		double e = 0.0;		/* Distance moved on 2nd previous step */

		/* Make sure a and b are in ascending order */
		if (ax > bx) {
			double tt;
			tt = ax; ax = bx; bx = tt;
			tt = af; af = bf; bf = tt;
		}

		wx = vx = xx;	/* Initial values of other center points */
		wf = xf = xf;

		/* Lookup derivative at x (we already have xf from bracketing) */
		for (i = 0; i < di; i++)
			xt[i] = cp[i] + xx * xi[i];
		(*dfunc)(fdata, df, xt);
		for (xd = 0.0, i = 0; i < di; i++)
			xd += xi[i] * df[i];
		wd = ud = xd;

		LDBG((" linmind: xx %f, deriv. xd %f\n",xx,xd))

		for (iter = 1; iter <= POWELL_MAXIT; iter++) {
			double mx = 0.5 * (ax + bx);		/* m is center of bracket values */
#ifdef ABSTOL
			double tol1 = ftol;			/* Absolute tollerance */
#else
			double tol1 = ftol * fabs(xx) + 1e-10;
#endif
			double tol2 = 2.0 * tol1;

			LDBG((" linmind it %d: Got bracket a:%f:%f x:%f:%f b:%f:%f\n",iter, ax,af,xx,xf,bx,bf))

			/* See if we're done */
			if (fabs(xx - mx) <= (tol2 - 0.5 * (bx - ax))) {
				LDBG((" linmind: We're done because %e <= %e\n",fabs(xx - mx), tol2 - 0.5 * (bx - ax)))
				break;
			}

			LDBG((" linmind: e %f tol2 %f\n",e,tol1))

			if (fabs(e) > tol1) {			/* Do a trial secant fit */
				double te;
				double dx1, dx2;				/* Secant extrapolation points */
				double ux1, ux2;
				int ch1, ch2;

				LDBG((" linmind: Doing trial secant fit\n" ))

				dx2 = dx1 = 2.0 * (bx - ax);	/* Default to values out of the ax..bx bracket */
				
				/* Extrapolated points from last two points (secant method) */
				if (wd != xd)
					dx1 = (wx - xx) * xd/(xd - wd);
				if (vd != xd)
					dx2 = (vx - xx) * xd/(xd - vd);

				ux1 = xx + dx1;
				ux2 = xx + dx2;

				/* Check which one is reasonable */
				ch1 = (ax - ux1)  * (ux1 - bx) > 0.0 && xd * dx1 < 0.0;
				ch2 = (ax - ux2)  * (ux2 - bx) > 0.0 && xd * dx2 < 0.0;

				LDBG((" linmind: Doing dx1 %f dx2 %f ux1 %f ux2 %f ch1 %d ch2 %d\n",dx1,dx2,ux1,ux2,ch1,ch2))

				te = e;				/* Save previous e value */
				e = de;				/* Previous steps distance moved */

				if (!ch1 && !ch2)
					goto bisect;

				/* Use smallest or the one that's valid */
				if (ch1 && ch2)
					de = fabs(dx1) < fabs(dx2) ? dx1 : dx2;
				if (ch1)
					de = dx1;
				else if (ch2)
					de = dx2;

				LDBG((" linmind: set de %f\n",de))

				if (fabs(de) > fabs(0.5 * te)) {
					LDBG((" linmind: abs(de) %f > abs(te/2 = %f)\n",fabs(de),fabs(0.5 * te)))
					goto bisect;
				}

				ux = xx + de;

				if ((ux - ax) < tol2 || (bx - ux) < tol2) {
					if ((mx - xx) < 0.0)
						de = -fabs(tol1);
					else
						de = fabs(tol1);
					LDBG((" linmind: Set de to tol1 %f\n",de))
				}
#ifdef LDEBUG
				  else {
					LDBG((" linmind: Using secant fit de %f\n",de))
				}
#endif

			/* else bisect picking side using sign of derivative */
			} else {
		  bisect:
				e = (xd >= 0.0 ? ax - xx : bx  -xx);
				de = 0.5 * e;
				LDBG((" linmind: Continuing bisection search de %f\n",de))
			}

			if (fabs(de) >= tol1) {		/* If de moves as much as tol1 or more */
				ux = xx + de;			/* use it */

				/* Evaluate function */
				for (i = 0; i < di; i++)
					xt[i] = cp[i] + ux * xi[i];
				uf = (*func)(fdata, xt);

				LDBG((" linmind: ux = %f = xx %f + de %f, uf %f\n",ux,xx,de,uf))

			} else {					/* else move by tol1 in direction de */
				if (de > 0.0) {
					ux = xx + tol1;
					LDBG((" linmind: ux = %f = xx %f + tol1 %f\n",ux,xx,tol1))
				} else {
					ux = xx - tol1;
					LDBG((" linmind: ux = %f = xx %f - tol1 %f\n",ux,xx,tol1))
				}
				/* Evaluate function */
				for (i = 0; i < di; i++)
					xt[i] = cp[i] + ux * xi[i];
				uf = (*func)(fdata, xt);

				LDBG((" linmind: uf %f\n",uf))

				if (uf > xf) {		/* If tol1 step downhill takes us uphill, we're done */
					goto done;
				}
			}

			/* Evaluate derivative at trial point */
			(*dfunc)(fdata, df, xt);
			for (ud = 0.0, i = 0; i < di; i++)
				ud += xi[i] * df[i];

			LDBG((" linmind: ux %f, deriv. ud %f\n",ux,ud))

			/* Houskeeping: */
			if (uf <= xf) {					/* Found new best solution */
				LDBG((" linmind: found new best solution at %f val %f dval %f\n",ux,uf,ud))
				if (ux >= xx) {	
					ax = xx; af = xf; 		/* New lower bracket */
				} else {
					bx = xx; bf = xf;		/* New upper bracket */
				}
				vx = wx; vf = wf; vd = wd;		/* New previous 2nd best solution */
				wx = xx; wf = xf; wd = xd;		/* New 2nd best solution from previous best */
				xx = ux; xf = uf; xd = ud;		/* New best solution from latest */

			} else {			/* Found a worse solution */
				LDBG((" linmind: found new worse solution at %f val %f dval %f\n",ux,uf,ud))
				LDBG((" linmind:             current best at %f val %f dval %f\n",xx,xf,xd))
				if (ux < xx) {
					ax = ux; af = uf;		/* New lower bracket */
				} else {
					bx = ux; bf = uf;		/* New upper bracket */
				}
				if (uf <= wf || wx == xx) {	/* New 2nd best solution, or equal best */
					vx = wx; vf = wf; vd = wd;	/* New previous 2nd best solution */
					wx = ux; wf = uf; wd = ud;	/* New 2nd best from latest */
				} else if (uf <= vf || vx == xx || vx == wx) {	/* New 3rd best, or equal 1st & 2nd */
					vx = ux; vf = uf; vd = ud;	/* New previous 2nd best from latest */
				}
			}
		}	/* Next itter */

		/* !!! should do something if iter > POWELL_MAXIT ???  */
		/* Solution is at xx, xf */

	  done:;
		if (di > 10) {
			free_dvector(df, 0, di-1);
			free_dvector(xt, 0, di-1);
		}
		/* Compute solution vector */
		LDBG((" linmind: computing soln from best at %f val %f dval %f\n",xx,xf,xd))
		for (i = 0; i < di; i++) 
			cp[i] += xx * xi[i];

	}	/* Minimizer context */

	return xf;
}

#undef POWELL_GOLD
#undef POWELL_MAXIT

