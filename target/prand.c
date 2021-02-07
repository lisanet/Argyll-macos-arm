
/* 
 * Argyll Color Correction System
 *
 * Perceptual space random test point class
 *
 * Author: Graeme W. Gill
 * Date:   12/9/2004
 *
 * Copyright 2004, 2009 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/* TTBD:

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include "aconfig.h"
#include "numlib.h"
#include "sort.h"
#include "icc.h"
#include "xicc.h"
#include "xcolorants.h"
#include "targen.h"
#include "prand.h"

static int prand_from_percept( prand *s, double *p, double *v);

/* ----------------------------------------------------- */

/* Default convert the nodes device coordinates into approximate perceptual coordinates */
/* (usually overriden by caller supplied function) */
static void
default_prand(void *od, double *p, double *d) {
	prand *s = (prand *)od;
	int e;

	/* Default Do nothing - copy device to perceptual. */
	for (e = 0; e < s->di; e++) {
		p[e] = d[e] * 100.0;
	}
}

/* Return the largest distance of the point outside the device gamut. */
/* This will be 0 if inside the gamut, and > 0 if outside.  */
static double
prand_in_dev_gamut(prand *s, double *d) {
	int e;
	int di = s->di;
	double tt, dd = 0.0;
	double ss = 0.0;

	for (e = 0; e < di; e++) {
		ss += d[e];

		tt = 0.0 - d[e];
		if (tt > 0.0) {
			if (tt > dd)
				dd = tt;
		}
		tt = d[e] - 1.0; 
		if (tt > 0.0) {
			if (tt > dd)
				dd = tt;
		}
	}
	tt = ss - s->ilimit;
	if (tt > 0.0) {
		if (tt > dd)
			dd = tt;
	}
	return dd;
}

/* --------------------------------------------------- */
/* Seed the object with the initial fixed points */

static void
prand_add_fixed(
prand *s,
fxpos *fxlist,			/* List of existing fixed points */
int fxno				/* Number in fixed list */
) {
	int e, di = s->di;
	int i;

	/* Add fixed points if there are any */
	if (fxno > 0) {

		for (i = 0; (i < fxno) && (i < s->tinp); i++) {
			prnode *p = &s->n[i];	/* Destination for point */

			for (e = 0; e < di; e++)
				p->p[e] = fxlist[i].p[e];

			p->fx = 1;			/* is a fixed point */
			s->percept(s->od, p->v, p->p);
			s->np = s->fnp = i+1;
		}
	}
}

/* Seed the object with the perceptual space random points. */
static void
prand_seed(prand *s) {
	int e, di = s->di;

	printf("\n");

	/* Seed the non-fixed points */
	for (; s->np < s->tinp;) {
		prnode *p = &s->n[s->np];		/* Next node */

		for (e = 0; e < di; e++) {
			if (e == 1 || e == 2)
				p->v[e] = d_rand(-128.0, 128.0);
			else
				p->v[e] = d_rand(0.0, 100.0);
		}
		if (prand_from_percept(s, p->p, p->v) == 0) {
			s->np++;
			printf("%cAdded %d/%d",cr_char,s->np,s->tinp); fflush(stdout);
		}
	}
	printf("\n");
}

/* Seed the object with the perceptual space quasi random points. */
static void
pqrand_seed(prand *s) {
	int e, di = s->di;
	sobol *sl = NULL;

	if ((sl = new_sobol(di)) == NULL)
		error("Creating sobol sequence generator failed");

	printf("\n");

	/* Seed the non-fixed points */
	for (; s->np < s->tinp;) {
		prnode *p = &s->n[s->np];		/* Next node */

		if (sl->next(sl, p->v))
			error("Run out of sobol random numbers!");

		for (e = 0; e < di; e++) {
			if (e == 1 || e == 2)
				p->v[e] = p->v[e] * 256.0 - 128.0;
			else
				p->v[e] *= 100.0;
		}
		if (prand_from_percept(s, p->p, p->v) == 0) {
			s->np++;
			printf("%cAdded %d/%d",cr_char,s->np,s->tinp); fflush(stdout);
		}
	}
	printf("\n");
	sl->del(sl);
}

/* --------------------------------------------------- */
/* Support accessing the list of generated sample points */

/* Rest the read index */
static void
prand_reset(prand *s) {
	s->rix = 0;
}

/* Read the next non-fixed point value */
/* Return nz if no more */
static int
prand_read(prand *s, double *p, double *f) {
	int e;

	/* Advance to next non-fixed point */
	while(s->rix < s->np && s->n[s->rix].fx)
		s->rix++;
	
	if (s->rix >= s->np)
		return 1;

	/* Return point info to caller */
	for (e = 0; e < s->di; e++) {
		if (p != NULL)
			p[e] = s->n[s->rix].p[e];
		if (f != NULL)
			f[e] = s->n[s->rix].v[e];
	}
	s->rix++;

	return 0;
}

/* --------------------------------------------------- */
/* Main object creation/destruction */

static void init_pmod(prand *s);

/* Destroy ourselves */
static void
prand_del(prand *s) {
	free(s->n);

	if (s->pmod != NULL)
		free(s->pmod);

	free (s);
}

/* Creator */
prand *new_prand(
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
int tinp,				/* Total number of points to generate, including fixed */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
int quasi,				/* nz to use quasi random (sobol) */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
) {
	prand *s;

	if ((s = (prand *)calloc(sizeof(prand), 1)) == NULL)
		error ("prand: malloc failed");

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	s->di = di;

	if (tinp < fxno)	/* Make sure we return at least the fixed points */
		tinp = fxno;

	s->tinp = tinp;		/* Target total number of points */
	s->ilimit = ilimit;

	/* Init method pointers */
	s->reset = prand_reset;
	s->read  = prand_read;
	s->del   = prand_del;

	/* If no perceptual function given, use default */
	if (percept == NULL) {
		s->percept = default_prand;
		s->od = s;
	} else {
		s->percept = percept;
		s->od = od;
	}

	/* Init the inverse perceptual function lookup */
	init_pmod(s);

	/* Allocate the space for the target number of points */
	if ((s->n = (prnode *)calloc(sizeof(prnode), s->tinp)) == NULL)
		error ("prand: malloc failed on sample nodes");
	s->np = s->fnp = 0;

	/* Setup the fixed points */
	prand_add_fixed(s, fxlist, fxno);

	if (tinp > fxno) { /* Create the perceptual space random points */
		if (quasi)
			pqrand_seed(s);
		else
			prand_seed(s);
	}

	prand_reset(s);		/* Reset read index */

	return s;
}

/* =================================================== */
/* Compute a simple but unbounded model of the */
/* perceptual function, used by inversion. We use the */
/* current vertex values to setup the model */
/* (Perhaps this should be moved to targen ?) */

/* A vertex point */
struct _vxpt {
	double p[MXTD];	/* Device position */
	double v[MXTD];		/* Perceptual value */
}; typedef struct _vxpt vxpt;

/* Structure to hold data for unbounded optimization function */
struct _ubfit {
	prand *s;				/* prand structure */
	vxpt *vxs;				/* List of vertex values */
	int _nvxs, nvxs;
}; typedef struct _ubfit ubfit;

/* Matrix optimisation function handed to powell() */
static double xfitfunc(void *edata, double *x) {
	ubfit *uf = (ubfit *)edata;
	prand *s = uf->s;
	int i, e, di = s->di;
	double rv = 0.0;

	/* For all the vertexes */
	for (i = 0; i < uf->nvxs; i++) {
		double v[MXTD], ev;

		/* Apply matrix cube interpolation */
		icxCubeInterp(x, di, di, v, uf->vxs[i].p);

		/* Evaluate the error */
		for (ev = 0.0, e = 0; e < di; e++) {
			double tt;
			tt = uf->vxs[i].v[e] - v[e];
			ev += tt * tt;
		}
		rv += ev;
	}
	return rv;
}

/* Fit the unbounded perceptual model to the perceptual function */
static void init_pmod(prand *s) {
	int i, ee, e, k, di = s->di;
	double *sa;
	double rerr;
	ubfit uf;

	uf.s = s;
	uf.vxs = NULL;
	uf.nvxs = uf._nvxs = 0;

	/* Allocate space for parameters */
	if ((s->pmod = malloc(di * (1 << di) * sizeof(double))) == NULL)
		error("Malloc failed for pmod");
	if ((sa = malloc(di * (1 << di) * sizeof(double))) == NULL)
		error("Malloc failed for pmod sa");

	/* Create a list of vertex values for the colorspace */
	/* Use in gamut vertexes, and compute clipped edges */
	for (ee = 0; ee < (1 << di); ee++) {
		double p[MXTD], ss;

		for (ss = 0.0, e = 0; e < di; e++) {
			if (ee & (1 << e))
				p[e] = 1.0;
			else
				p[e] = 0.0;
			ss += p[e];
		}
		if (ss < s->ilimit) {		/* Within gamut */
			if (uf.nvxs >= uf._nvxs) {
				uf._nvxs = 5 + uf._nvxs * 2;
				if ((uf.vxs = (vxpt *)realloc(uf.vxs, sizeof(vxpt) * uf._nvxs)) == NULL)
					error ("Failed to malloc uf.vxs");
			}
			for (k = 0; k < di; k++)
				uf.vxs[uf.nvxs].p[k] = p[k];
			uf.nvxs++;
		} else if ((ss - 1.0) < s->ilimit) { 		/* far end of edge out of gamut */
			double max = s->ilimit - (ss - 1.0);	/* Maximum value of one */
			for (e = 0; e < di; e++) {
				if ((ee & (1 << e)) == 0)
					continue;
				p[e] = max;
				if (uf.nvxs >= uf._nvxs) {
					uf._nvxs = 5 + uf._nvxs * 2;
					if ((uf.vxs = (vxpt *)realloc(uf.vxs, sizeof(vxpt) * uf._nvxs)) == NULL)
						error ("Failed to malloc uf.vxs");
				}
				for (k = 0; k < di; k++)
					uf.vxs[uf.nvxs].p[k] = p[k];
				uf.nvxs++;

				p[e] = 1.0;		/* Restore */
			}
		}	/* Else whole edge is out of gamut */
	}

	/* Lookup perceptual values */
	for (i = 0; i < uf.nvxs; i++) {
		s->percept(s->od, uf.vxs[i].v, uf.vxs[i].p);
//printf("~1 vtx %d: dev %f %f %f, perc %f %f %f\n",i, uf.vxs[i].p[0], uf.vxs[i].p[1], uf.vxs[i].p[2], uf.vxs[i].v[0], uf.vxs[i].v[1], uf.vxs[i].v[2]);
	}

	/* Setup matrix to be closest values initially */
	for (e = 0; e < (1 << di); e++) {	/* For each colorant combination */
		int j, f;
		double bdif = 1e6;
		double ov[MXTD];
		int bix = -1;
	
		/* Search the vertex list to find the one closest to this input combination */
		for (i = 0; i < uf.nvxs; i++) {
			double dif = 0.0;

			for (j = 0; j < di; j++) {
				double tt;
				if (e & (1 << j))
					tt = 1.0 - uf.vxs[i].p[j];
				else
					tt = 0.0 - uf.vxs[i].p[j];
				dif += tt * tt;
			}
			if (dif < bdif) {		/* best so far */
				bdif = dif;
				bix = i;
				if (dif < 0.001)
					break;			/* Don't bother looking further */
			}
		}
		for (f = 0; f < di; f++)
			 s->pmod[f * (1 << di) + e] = uf.vxs[bix].v[f];
	}

	for (e = 0; e < (di * (1 << di)); e++)
		sa[e] = 10.0;

	if (powell(&rerr, di * (1 << di), s->pmod, sa, 0.001, 1000,
	                                    xfitfunc, (void *)&uf, NULL, NULL) != 0) {
		warning("Powell failed to converge, residual error = %f",rerr);
	}

#ifdef DEBUG
	printf("Perceptual model fit residual = %f\n",sqrt(rerr));
#endif
	s->pmod_init = 1;

	free(sa);
}

/* Clip a device value to the gamut */
static int
prand_clip_point(prand *s, double *cd, double *d) {
	int e, di = s->di;
	double ss = 0.0;
	int rv = 0;

	for (e = 0; e < di; e++) {
		ss += d[e];						
		cd[e] = d[e];
		if (cd[e] < 0.0) {
			cd[e] = 0.0;
			rv |= 1;
		} else if (cd[e] > 1.0) {
			cd[e] = 1.0;
			rv |= 1;
		}										\
	}

	if (ss > s->ilimit) {
		ss = (ss - s->ilimit)/s->di;
		for (e = 0; e < di; e++)
			cd[e] -= ss;
		rv |= 1;
	}
	return rv;
}

/* Unbounded perceptual lookup. */
/* return nz if it was actually clipped and extended */
static int prand_cc_percept(prand *s, double *v, double *p) {
	double cp[MXTD];
	int clip;

	clip = prand_clip_point(s, cp, p);

	s->percept(s->od, v, cp); 

	/* Extend perceptual value using matrix model */
	if (clip) {
		int e, di = s->di;
		double mcv[MXTD], zv[MXTD];

#ifdef DEBUG
		if (s->pmod_init == 0)
			error("ofps_cc_percept() called before pmod has been inited");
#endif
		/* Lookup matrix mode of perceptual at clipped device */
		icxCubeInterp(s->pmod, di, di, mcv, cp);

		/* Compute a correction factor to add to the matrix model to */
		/* give the actual perceptual value at the clipped location */
		for (e = 0; e < di; e++)
			zv[e] = v[e] - mcv[e]; 

		/* Compute the unclipped matrix model perceptual value */
		icxCubeInterp(s->pmod, di, di, v, p);

		/* Add the correction value to it */
		for (e = 0; e < di; e++)
			v[e] += zv[e]; 
	}
	return clip;
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Reverse lookup function :- perceptual to device coordinates */
/* Using dnsq */

/* Structure to hold data for optimization function */
struct _rdatas {
	prand *s;			/* prand structure */
	double *ptv;		/* Perceptual target value */
}; typedef struct _rdatas rdatas;


/* calculate the functions at x[] */
int prand_dnsq_solver(	/* Return < 0 on abort */
	void *fdata,	/* Opaque data pointer */
	int n,			/* Dimenstionality */
	double *x,		/* Multivariate input values */
	double *fvec,	/* Multivariate output values */
	int iflag		/* Flag set to 0 to trigger debug output */
) {
	rdatas *ed = (rdatas *)fdata;
	prand *s = ed->s;
	double v[MXTD];
	int e, di = s->di;

	prand_cc_percept(s, v, x);

	for (e = 0; e < di; e++)
		fvec[e] = ed->ptv[e] - v[e];
		
//printf("~1 %f %f %f from %f %f %f\n", fvec[0], fvec[1], fvec[2], x[0], x[1], x[2]);
	return 0;
}

/* Given a point in perceptual space, an approximate point */
/* in device space, return the device value corresponding to */
/* the perceptual value, plus the clipped perceptual value. */
/* Return 1 if the point is out of gamut or dnsq failed. */
static int
prand_from_percept(
prand *s,
double *p,			/* return (clipped) device position */
double *v			/* target perceptual */
) {
	int e, di = s->di;
	rdatas ed;
	double ss;		/* Initial search area */
	double fvec[MXTD];	/* Array that will be RETURNed with thefunction values at the solution */
	double dtol;	/* Desired tollerance of the solution */
	double tol;		/* Desired tollerance of root */
	int maxfev;		/* Maximum number of function calls. set to 0 for automatic */
	int rv;

//printf("~1 percept2 called with %f %f %f\n", v[0], v[1], v[2]);
	ed.s = s;
	ed.ptv = v;			/* Set target perceptual point */

	for (e = 0; e < di; e++)
		p[e] = 0.3;				/* Start location */
	ss = 0.1;
	dtol = 1e-6;
	tol = 1e-8;
	maxfev = 1000;

	rv = dnsqe((void *)&ed, prand_dnsq_solver, NULL, di, p, ss, fvec, dtol, tol, maxfev, 0);

	if (rv != 1 && rv != 3) { /* Fail to converge */
//printf("~1 failed with rv %d\n",rv);
		return 1;
	}

//printf("~1 got soln %f %f %f\n", p[0], p[1], p[2]);
	if (prand_clip_point(s, p, p)) {
//printf("~1 clipped\n");
		return 1;
	}

	return 0;
}
















