
/* This is better for stocastic optimisation, where the function */
/* being evaluated may have a random component, or is not smooth. */

/*
 * Copyright 1999 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* A general purpose downhill simplex multivariate optimser, */
/* based on the Nelder and Mead algorithm. */
/* Code is an original expression of the algorithms decsribed in */
/* "Numerical Recipes in C", by W.H.Press, B.P.Flannery, */
/* S.A.Teukolsky & W.T.Vetterling. */

#include "numsup.h"

#define DEBUG

int dhsx_debug = 0;

static void simplexinit(int di, double *cp, double **p, double *r, double sv, int ii);
static double trypoint(int di,double *cp, double **p, double *y, int hix, double hpfac,
                double (*funk)(void *fdata, double *tp), void *fdata, double *tryp);

#ifdef NEVER	/* Experimental */

#define ALPHA 0.7			/* Extrapolate hight point through oposite face factor */
#define GAMMA 1.4			/* Aditional extrapolation if ALPHA is good */ 
#define BETA 0.4			/* One dimensional contraction factor (smaller is more) */
#define DELTA 0.5			/* Multi dimensional contraction factor (smaller is more) */
#define NONEXP 2			/* non expanding passes */

#else			/* Standard tuning values */

#define ALPHA 1.0			/* [1.0] Extrapolate hight point through oposite face factor */
#define GAMMA 2.0			/* [2.0] Aditional extrapolation if ALPHA is good */ 
#define BETA 0.4			/* [0.5] One dimensional contraction factor (smaller is more) */
#define DELTA 0.4			/* [0.5] Multi dimensional contraction factor (smaller is more) */
#define NONEXP 3			/* [3] non expanding passes */

#endif

/* Down hill simplex function */
/* return 0 on sucess, 1 on failure due to excessive itterations */
/* Result will be in cp */
int dhsx(
double *rv,				/* If not NULL, return the residual error */
int di,					/* Dimentionality */
double *cp,				/* Initial starting point, return minimum */
double *s,				/* Size of initial search area */
double ftol,			/* Finishing tollerance of error change */
double athr,			/* Absolute return value threshold. (Set high to not use) */
int maxit,				/* Maximum iterations allowed */
double (*funk)(void *fdata, double *tp),		/* Error function to evaluate */
void *fdata				/* Data needed by function */
) {
	int ii = 0;			/* Initial simplex orientation */
	int i, j;
	int nit;			/* Number of iterations */
	int nsp = di+1;		/* Number of simplex verticy points */
	double tryy, ysave;
	double tol;
	double **p;			/* Current simplex array */
	double *y;			/* Values of func at verticies */
	double **p2;		/* Trial simplex array */
	double *y2;			/* Trial values of func at verticies */
	int lox, hix, nhix;	/* Lowest point index, highest point, next highest point */
	double *tryp;		/* Temporary used by trypoint() */

	/* Allocate array arrays */
	tryp = dvector(0, di-1);		/* Trial value */		
	p = dmatrix(0, nsp-1, 0, di-1);	/* Vertex array of dimentions */
	y = dvector(0, nsp-1);			/* Value of function at verticies */
	p2 = dmatrix(0, nsp-1, 0, di-1);	/* Trial vertex array of dimentions */
	y2 = dvector(0, nsp-1);			/* Trial value of function at verticies */
	
	/* Init the search simplex */
	simplexinit(di, cp, p, s, 1.0, ii);

	/* Compute initial y (function) values at simplex verticies */
	for (i = 0; i < nsp; i++)			/* For all verticies */
		y[i] = (*funk)(fdata, p[i]);	/* Compute error function */

	/* Locate verticy with best value */
	lox = 0;
	for (i = 0; i < nsp; i++) {
		if (y[i] < y[lox]) 
			lox = i;
	}
	tryy = (*funk)(fdata, cp);	/* Value at initial point */
#ifdef DEBUG
	if (dhsx_debug) printf(" initial point %s = %e\n",debPdv(di,cp),tryy);
#endif /* DEBUG */

	/* If our initial point is better than any of the simplex verticies */
	if (y[lox] > tryy) {
#ifdef DEBUG
		if (dhsx_debug) printf(" initial point is better than surrounding simplex\n");
#endif /* DEBUG */
		/* Move all the verticies to match moving lox to cp */
		for (i = 0; i < nsp; i++) {
			if (i == lox)
				continue;
			for (j = 0; j < di; j++)
				p[i][j] += cp[j] - p[lox][j]; 
			y[i] = (*funk)(fdata, p[i]);	/* Compute error function */
		}
		/* Make lox be the input point */
		for (j = 0; j < di; j++)
			p[lox][j] = cp[j];
		y[lox] = tryy;
	} 

	/* Compute current center point location as sum of verticies. */
	/* (We use this to compute moves) */
	for (j = 0; j < di; j++) {					/* For all dimensions */
		double sum;
		for (i = 0, sum = 0.0; i < nsp; i++)	/* For all verticies */
			sum += p[i][j];
		cp[j] = sum;
	}

	/* Untill we find a solution or give up */
	for (nit = 0; ; nit++) {

		/* Find highest, next highest and lowest vertex */
		lox = nhix = hix = 0;
		for (i = 0; i < nsp; i++) {
			if (y[i] < y[lox]) 
				lox = i;
			if (y[i] > y[hix]) {
				nhix = hix;
				hix = i;
			} else if (y[i] > y[nhix]) {
				nhix = i;
			}
		}

		tol = y[hix] - y[lox];

#ifdef DEBUG
		if (dhsx_debug) {
			printf("Current vs =\n");
			for (i = 0; i < nsp; i++)
				printf(" %d: %s\n",i,debPdv(di, p[i]));
			printf("Current errs = %s\n",debPdv(nsp,y));
			printf("Current y[lox] = %e, y[hix] = %e\n",y[lox], y[hix]);
		}
#endif /* DEBUG */

		/* If we look like we are about to finish, */
		/* see if we should re-start with a new simplex. */
		if (tol < ftol && y[lox] < athr 	/* Found an adequate solution */
		 && nit < maxit) {
			double scale = 0.0;
			int lox2;

#ifdef DEBUG
			if (dhsx_debug) printf(" nit %d, tol %e\n",nit, tol);
#endif /* DEBUG */

			/* compute center location */
			tryy = 1.0/nsp;
			for (j = 0; j < di; j++)		/* For all dimensions */
				cp[j] *= tryy;				/* Set cp to center point of simplex */

			/* Compute scaled distance of vertexes from center */
			for (i = 0; i < nsp; i++) {
				double dist = 0.0;
				for (j = 0; j < di; j++) {
					double tt = (cp[j] - p[i][j])/s[j];
					dist += tt * tt;
				}
				scale += sqrt(dist);
			}
			scale /= (double)nsp;		/* Average scale compared to starting simplex */
#ifdef DEBUG
			if (dhsx_debug) printf(" ave scale = %f\n",scale);
#endif /* DEBUG */

			/* Enlarge search space, but not more than initial */
			scale *= 10.0;
			if (scale > 1.0)
				scale = 1.0;

			/* Compute trial simplex with different orientation */
			if (++ii >= (di+1))
				ii = 0;

			/* Init the search simplex */
			simplexinit(di, cp, p2, s, scale, ii);

			/* Compute y (function) values at simplex verticies */
			for (i = 0; i < nsp; i++)			/* For all verticies */
				y2[i] = (*funk)(fdata, p2[i]);	/* Compute error function */

			/* Locate verticy with best value */
			lox2 = 0;
			for (i = 0; i < nsp; i++) {
				if (y2[i] < y2[lox2]) 
					lox2 = i;
			}
#ifdef DEBUG
			if (dhsx_debug) printf(" y2lox %f ylox %f\n",y2[lox2], y[lox]);
#endif /* DEBUG */

			/* If any of its vertexes are better than current best, switch */
			/* to it and continue (i.e. re-start) */
			if (y2[lox2] < y[lox]) {

#ifdef DEBUG
				if (dhsx_debug) printf(" restarting\n");
#endif /* DEBUG */

				for (i = 0; i < nsp; i++) {
					for (j = 0; j < di; j++)
						p[i][j] = p2[i][j];
					y[i] = y2[i];
				}

				/* Compute current center point location as sum of verticies. */
				/* (We use this to compute moves) */
				for (j = 0; j < di; j++) {					/* For all dimensions */
					double sum;
					for (i = 0, sum = 0.0; i < nsp; i++)	/* For all verticies */
						sum += p[i][j];
					cp[j] = sum;
				}

				/* Find highest, next highest and lowest vertex */
				lox = nhix = hix = 0;
				for (i = 0; i < nsp; i++) {
					if (y[i] < y[lox]) 
						lox = i;
					if (y[i] > y[hix]) {
						nhix = hix;
						hix = i;
					} else if (y[i] > y[nhix]) {
						nhix = i;
					}
				}

				tol = y[hix] - y[lox];
			}
		}

		if ((tol < ftol && y[lox] < athr) 	/* Found an adequate solution */
		 || ((nit+1) >= maxit)) {				/* Or we are about to fail */

			/* convert cp[] to center point location, */
			/* and use best out of it and any simplex verticy. */
			tryy = 1.0/nsp;
			for (j = 0; j < di; j++)		/* For all dimensions */
				cp[j] *= tryy;				/* Set cp to center point of simplex */
#ifdef DEBUG
			if (dhsx_debug) printf("C point = %s\n",debPdv(di,cp));
#endif
			tryy = (*funk)(fdata, cp);		/* Compute error function */

			if (tryy > y[lox]) {			/* Center point is not the best */
#ifdef DEBUG
				if (dhsx_debug) printf("C point val %f is not best, using sx %d val %f instead\n",tryy,lox,y[lox]);
#endif
				tryy = y[lox];
				for (j = 0; j < di; j++)
					cp[j] = p[lox][j];
			}
#ifdef DEBUG
			else if (dhsx_debug) printf("C point val %f is best\n",tryy);
#endif
			free_dvector(y2, 0, nsp-1);
			free_dmatrix(p2, 0, nsp-1, 0, di-1);
			free_dvector(y, 0, nsp-1);
			free_dmatrix(p, 0, nsp-1, 0, di-1);
			free_dvector(tryp, 0, di-1);
#ifdef DEBUG
			if (dhsx_debug) printf("Total itterations = %d\n",nit);
#endif
			if (rv != NULL)
				*rv = tryy;
			if ((nit+1) >= maxit)
				return 1;	/* Failed */
			return 0;
		}

		/* Only try expanding after a couple of iterations */
		if (nit > NONEXP) {
			/* Try moving the high point through the oposite face by ALPHA */
#ifdef DEBUG
			if (dhsx_debug) printf("dhsx: try moving high point %d through oposite face",hix);
#endif
			tryy = trypoint(di, cp, p, y, hix, -ALPHA, funk, fdata, tryp);
		}

		/* If gave good result, continue on in that direction */
		if (nit > NONEXP && tryy <= y[lox]) {
#ifdef DEBUG
			if (dhsx_debug) printf("dhsx: moving high through oposite face worked");
#endif
			tryy = trypoint(di, cp, p, y, hix, GAMMA, funk, fdata, tryp);


		/* else if ALPHA move made things worse, do a one dimensional */
		/* contraction by a factor BETA */
		} else if (nit <= NONEXP || tryy >= y[nhix]) {

#ifdef DEBUG
			if (dhsx_debug) printf("dhsx: else try moving contracting point %d, y[ini] = %f",hix,y[hix]);
#endif
			ysave = y[hix];
			tryy = trypoint(di, cp, p, y, hix, BETA, funk, fdata, tryp);
			
			if (tryy >= ysave) {
#ifdef DEBUG
				if (dhsx_debug) printf("dhsx: contracting didn't work, try contracting other points to low");
#endif
				/* That still didn't help us, so move all the */
				/* other points towards the low point */
				for (i = 0; i < nsp; i++) {	/* For all verts except low */
					if (i != lox) {
						for (j = 0; j < di; j++) 	/* For all dimensions */
							p[i][j] = DELTA * p[i][j] + (1.0 - DELTA) * p[lox][j];
						y[i] = (*funk)(fdata, p[i]);	/* Compute function value for new point */
					}
				}
				/* Re-compute current center point location */
				for (j = 0; j < di; j++) {
					double sum;
					for (i = 0,sum = 0.0;i<nsp;i++)
						sum += p[i][j];
					cp[j] = sum;
				}
			} else {
#ifdef DEBUG
				if (dhsx_debug) printf("dhsx: contracting point %d worked, tryy = %e, ysave = %e",hix,tryy,ysave);
#endif
			}
		}
	}
}

/* Try moving the high point through the opposite face */
/* by a factor of fac, and replaces the high point if */
/* that proves to be better.  Return the failed or new */
/* function value. */
static double trypoint(
int di,					/* Dimentionality */
double *cp,				/* nsp * center coord/Returned coordinate */
double **p,				/* Starting/Current simplex (modified by dhsx) */
double *y,				/* values of func at verticies */
int hix,				/* Index of high point we are moving */
double hpfac,			/* factor to move high point */
double (*funk)(void *fdata, double tp[]),		/* Error function to evaluate */
void *fdata,			/* Data needed by function */
double *tryp			/* temporary array of size di-1 */
) {
	int j;
	double tt, tryy;

	/* Compute trial high point */
	tt = (1.0 - hpfac)/di;
	for (j = 0; j < di; j++) 
		tryp[j] = cp[j] * tt - p[hix][j] * (tt - hpfac);

	/* Evaluate trial point */
	tryy = (*funk)(fdata, tryp);	/* Compute error function */

	/* If new high point pos. is better */
	if (tryy < y[hix]) {
#ifdef DEBUG
		if (dhsx_debug) printf("Try gave improved %e from sx %d",tryy,hix);
#endif
		y[hix] = tryy;		/* Replace func val of hi with trial */

		for (j = 0; j < di; j++) {
			cp[j] += tryp[j] - p[hix][j];	/* Recompute cp */
			p[hix][j] = tryp[j];	/* Replace co-ords of hi with trial */
		}
	} else {
#ifdef DEBUG
		if (dhsx_debug) printf("Try gave worse %e from sx %d",tryy,hix);
#endif
	}
	return tryy;		/* Function value of trial point */
}


/* Make up an initial simplex for dhsx routine */
static void
simplexinit(
int di,			/* Dimentionality */
double *cp,		/* Initial solution location */
double **p,		/* Simplex to initialize */
double *s,		/* initial radius for each dimention */
double sv,		/* Radius scaling value */
int ii			/* Coordinate to start with */
) {
	double bb;
	double hh = 0.5;			/* Constant */
	double rr = sqrt(3.0)/2.0;	/* Constant */
	int i, j;
	for (i = 0; i < (di+1); i++) {	/* For each vertex */
		/* The bounding points form a equalateral simplex */
		/* whose vertexes are on a sphere about the data */
		/* point center. The coordinate sequence is: */
		/*  B = sphere radius */
		/*	H = 0.5         */
		/*	R = sqrt(3)/2   */
		/*   0  0  0 +B     */
		/*   0  0  0 -B     */

		/*   0  0   0  +B   */
		/*   0  0  +RB -HB  */
		/*   0  0  -RB -HB  */

		/*   0  0      0   +B   */
		/*   0  0     +RB  -HB  */
		/*   0  +RRb  -HRB -HB  */
		/*   0  -RRb  -HRB -HB  */

		/*   0       0      0   +B   */
		/*   0       0     +RB  -HB  */
		/*   0      +RRb   -HRB -HB  */
		/*   +RRRb  -HRRb  -HRB -HB  */
		/*   -RRRb  -HRRb  -HRB -HB  */

		/*      etc.     */
		bb = 1.0;		/* Initial unscaled radius */
		for (j = 0; j < di; j++) {	/* For each coordinate in vertex */
			if (j > ii)
				p[i][j] = cp[j] + sv * s[j] * 0.0;		/* If beyond last */
			else if (j == ii)		/* If last non-zero */
				p[i][j] = cp[j] + sv * s[j] * bb;
			else if (ii == di && j == (di-1)) /* If last of last */
				p[i][j] = cp[j] + sv * s[j] * -1.0 * bb;
			else					/* If before last */
				p[i][j] = cp[j] + sv * s[j] * -hh * bb;
			bb *= rr;
		}
		/* Increment coordinate number with wrap around */
		if (++ii >= (di+1))
			ii = 0;
	}
#ifdef DEBUG
	if (dhsx_debug) {
		for (i = 0; i < (di+1); i++)
			printf(" p[%d] = %s\n",i,debPdv(di,p[i]));
	}
#endif
}

