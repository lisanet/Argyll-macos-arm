
// ppoint7c
// Approach that picks poorly supprted points with maximum interpolation
// error each time. Version that creates a candidate list when adding
// previous points to the distance grid.
// Development of version that uses interpolation error and perceptual
// distance to nearest sample point driven point placement metric, this
// one usin incremental rspl for interpolation estimation.

/* 
 * Argyll Color Correction System
 *
 * Perceptually distributed point class
 *
 * Author: Graeme W. Gill
 * Date:   5/10/96
 *
 * Copyright 1996 - 2004 Graeme W. Gill
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
#include "numlib.h"
#include "rspl.h"
#include "sort.h"
#include "icc.h"
#include "xcolorants.h"
#include "targen.h"
#include "ppoint.h"
#ifdef DUMP_PLOT
# include "plot.h"
# include "ui.h"
#endif

#undef DEBUG
#define DUMP_PLOT		/* Show on screen plot */
#define PERC_PLOT 0		/* Emit perceptive space plots */
#define DO_WAIT 1		/* Wait for user key after each plot */

#define ALWAYS
#undef NEVER

#ifdef NEVER
#ifdef	__STDC__
#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);
#else
#include <varargs.h>
void error(), warning(), verbose();
#endif
#endif	/* NEVER */

#ifdef STANDALONE_TEST
#ifdef DUMP_PLOT
static void dump_image(ppoint *s, int pcp);
#endif
#endif

static void add_dist_points(ppoint *s, co *pp, int nn);
//static double far_dist(ppoint *s, double *p);

/* Default convert the nodes device coordinates into approximate perceptual coordinates */
/* (usually overriden by caller supplied function) */
static void
default_ppoint_to_percept(void *od, double *p, double *d) {
	ppoint *s = (ppoint *)od;
	int e;

#ifndef NEVER
	/* Default Do nothing - copy device to perceptual. */
	for (e = 0; e < s->di; e++) {
		double tt = d[e];
		if (e == 0)
			tt = pow(tt, 2.0);
		else
			tt = pow(tt, 0.5);
		p[e] = tt * 100.0;
	}
#else
	for (e = 0; e < s->di; e++) {
		double tt = d[e];
		/* Two slopes with a sharp turnover in X */
		if (e == 0) {
			if (tt < 0.5)
				tt = tt * 0.3/0.5;
			else
				tt = 0.3 + ((tt-0.5) * 0.7/0.5);
		}
		p[e] = tt * 100.0;
	}
#endif
}

/* return the distance of the device value from the device gamut */
/* This will be -ve if the point is outside */
/* If bvp is non-null, the index of the closest dim times 2 */
/* will be returned for the 0.0 boundary, dim * 2 + 1 for the 1.0 */
/* boundary, and di * 2 for the ink limit boundary. */
static double
ppoint_in_dev_gamut(ppoint *s, double *d, int *bvp) {
	int e;
	int di = s->di;
	double tt, dd = 1.0;
	double ss = 0.0;
	int bv = di;
	for (e = 0; e < di; e++) {
		tt = d[e];
		if (tt < dd) {
			dd = tt;
			bv = e * 2;
		}
		tt = 1.0 - d[e];
		if (tt < dd) {
			dd = tt;
			bv = e * 2 + 1;
		}
		ss += d[e];
	}
	ss = (s->ilimit-ss)/di;	/* Axis aligned distance to ink limit */
	tt = sqrt((double)di) * ss;	/* Diagonal distance to ink limit */
	if (tt < dd) {
		dd = tt;
		bv = di * 2;
	}
	if (bvp != NULL)
		*bvp = bv;
	return dd;
}

#ifdef NEVER	/* Not currently used */
/* Given the new intended device coordinates, */
/* clip the new position to the device gamut edge */
/* return non-zero if the point was clipped */
static int
ppoint_clip_point(ppoint *s, double *d) {
	int e;
	double ss = 0.0;
	int rv = 0;
	for (e = 0; e < s->di; e++) {
		if (d[e] < 0.0) {
			d[e] = 0.0;
			rv |= 1;
		} else if (d[e] > 1.0) {
			d[e] = 1.0;
			rv |= 1;
		}
		ss += d[e];
	}
	if (ss > s->ilimit) {
		ss = (ss - s->ilimit)/s->di;
		for (e = 0; e < s->di; e++)
			d[e] -= ss;
		rv |= 1;
	}
	return rv;
}
#endif /* NEVER */

/* --------------------------------------------------- */
/* Locate the best set of points to add */

/* Definition of the optimization functions handed to powell(.) */
/* Return distance error to be minimised (maximises distance from */
/* an existing sample point) */
static double efunc1(ppoint *s, double p[]) {
	double rv = 0.0;		/* return value */

//printf("\n~1 p = %f %f\n",p[0],p[1]);
	if ((rv = (ppoint_in_dev_gamut(s, p, NULL))) < 0.0) {
		rv = rv * -500.0 + 50000.0;		/* Discourage being out of gamut */
//printf("~1 out of gamut, rv = %f\n",rv);

	} else {
		int e, di = s->di;
		double vf[MXPD];		/* Perceptual value of reference */
		co tp;					/* Lookup from interpolation grid */
		double ierr;			/* Interpolation error */
		double cdist;			/* closest distance to point */
		double errd;			/* Overall error/distance to maximise */

		for (e = 0; e < di; e++)
			tp.p[e] = p[e];

		s->pd->interp(s->pd, &tp);	/* Lookup current closest distance value */
		cdist = tp.v[di];
		if (cdist >= 10000.0)		/* Initial value */
			cdist = 0.0;

//printf("~1 min pdist = %f\n",cdist);

		/* Not quite sure which is best here. */
		/* Using percept() is slower, and has more point placement artefacts, */
		/* but seems to arrive at a better result. */
#ifdef NEVER
		for (e = 0; e < di; e++)
			vf[e] = tp.v[e];		/* Use interpolated perceptual value */
#else
		s->percept(s->od, vf, p);	/* Lookup perceptual value */
#endif
		s->g->interp(s->g, &tp);	/* Lookup current interpolation */

//printf("~1 interp %f %f, percept %f %f\n",tp.v[0],tp.v[1],vf[0],vf[1]);
		for (ierr = 0.0, e = 0; e < di; e++) {
			double tt = tp.v[e] - vf[e]; 
			ierr += tt * tt;
		}
		ierr = sqrt(ierr);
//printf("~1 interp error = %f\n",ierr);

		/* The ratio of interpolation error to support distance affects */
		/* peak vs. average error in final result. */
#ifdef NEVER
		/* Weighted squares */
		errd = ierr * ierr + DWEIGHT * cdist * cdist;
#else
		/* Linear weighted then squared */
		errd = ierr + DWEIGHT * cdist;
		errd = errd * errd;
#endif

		/* Convert max error to min return value */
		rv = 1000.0/(0.1 + errd);
//printf("~1 err val %f\n",rv);

	}

//printf("~1 efunc1 returning %f from %f %f\n",rv,p[0],p[1]);
	return rv;
}


/* return the interpolation error at the given device location */
static double
ppoint_ierr(
ppoint *s,
double *p
) {
	int e, di = s->di;
	double vf[MXPD];	/* Perceptual value of reference */
	double err;
	co tp;				/* Perceptual value of test point */

	for (e = 0; e < di; e++)
		tp.p[e] = p[e];
	s->g->interp(s->g, &tp);

	s->percept(s->od, vf, p);

	for (err = 0.0, e = 0; e < di; e++) {
		double tt = tp.v[e] - vf[e]; 
		err += tt * tt;
	}
	err = sqrt(err);

	return err;
}

/* Find the next set of points to add to our test sample set. */
/* Both device and perceptual value are returned. */
/* We try and do a batch of points because adding points to the rspl interpolation */
/* is a high cost operation. The main trap is that we may add points that are almost identical, */
/* since we don't know the effect of adding other points in this batch. */
/* To try and counter this, points are rejected that are two close together in this group. */

/* Candidate points are located that have amongst the largest distances to existing */
/* points (measured in a device/perceptual distance mix), and from those points, */
/* the ones with the highest current interpolation mis-prediction error are selected. */
/* In this way a well spread set of samples is hoped to be gemerated, but favouring */
/* those that best reduce overall interpolation error. */
static int
ppoint_find_worst(
ppoint *s,
co *p,			/* return device values */
int tnn			/* Number to return */
) {
	co *fp = s->fwfp;	/* Copy of s-> info, stored in s because of size. */
	int nfp;			/* Current number in fp[] */
	int opoints;
	int e, di = s->di;
	double sr[MXPD];	/* Search radius */
	int i, j;

	for (e = 0; e < di; e++)
		sr[e] = 0.01;			/* Device space search radius */

//printf("~1 currently %d points in fp list\n",s->nfp);

	/* The distance grid functions will have a list of the FPOINTS best */
	/* grid points to start from. Make a copy of it */
	for (nfp = 0; nfp < s->nfp; nfp++) {
		fp[nfp] = s->fp[nfp];		/* Structure copy */
		fp[nfp].v[0] = efunc1(s, fp[nfp].p);	/* Compute optimiser error value */
	}

	/* If list is not full, fill with random numbers: */
	if (nfp < FPOINTS) {
//printf("~1 not full, so adding %d random points\n",FPOINTS-nfp);
//		for (; nfp < FPOINTS; nfp++) {
		for (; nfp < tnn; nfp++) {
			double sum;

			for (;;) {		/* Find an in-gamut point */
				for (sum = 0.0, e = 0; e < di; e++)
					sum += fp[nfp].p[e] = d_rand(0.0, 1.0);
				if (sum <= s->ilimit)
				break;
			}
			fp[nfp].v[0] = efunc1(s, fp[nfp].p);	/* Compute optimiser dist error value */
		}
	}

	/* Sort them by derr, smallest to largest */
#define HEAP_COMPARE(A,B) ((A).v[0] < (B).v[0])
	HEAPSORT(co, fp, nfp);
#undef HEAP_COMPARE

	opoints = nfp < OPOINTS ? nfp : OPOINTS;

	/* Optimise best portion of the list of starting points, according to */
	/* interpolation error weighted distance. */
	for (i = 0; i < opoints; i++) {
		double mx;

		if (powell(&mx, di, fp[i].p, sr,  0.001, 1000, 
		(double (*)(void *, double *))efunc1, (void *)s, NULL, NULL) != 0 || mx >= 50000.0) {
#ifdef ALWAYS
			printf("ppoint powell failed, tt = %f\n",mx);
#endif
		}
		fp[i].v[0] = mx;
//printf("~1 optimised point %d to %f %f derr %f\n",i,fp[i].p[0],fp[i].p[1],mx);

		/* Check if this duplicates a previous point */
		for (j = 0; j < i; j++) {

			double ddif = 0.0;
			for (e = 0; e < di; e++) {
				double tt = fp[i].p[e] - fp[j].p[e];
				ddif += tt * tt;
			}
			ddif = sqrt(ddif);			/* Device value difference */
			if (ddif < CLOSED) {
//printf("~1 duplicate of %d, so marked\n",j);
				fp[i].v[0] = 50000.0;	/* Mark so it won't be used */
				break;		/* too close */
			}
		}
	}

//printf("~1 derr sorted list:\n");
//for (i = 0; i < opoints; i++)
//	printf("~1 %d: loc %f %f derr %f\n", i, fp[i].p[0],fp[i].p[1],fp[i].v[0]);

	/* Compute the interpolation error for the points of interest */
	for (i = 0; i < opoints; i++) {
		if (fp[i].v[0] >= 50000.0)		/* Duplicate or failed to optimis point */
			fp[i].v[0] = -1.0;			/* Impossibly low interpolation error */
		else
			fp[i].v[0] = ppoint_ierr(s, fp[i].p);
	}

	/* Sort them by ierr, largest to smallest */
#define HEAP_COMPARE(A,B) ((A).v[0] > (B).v[0])
	HEAPSORT(co, fp, opoints);
#undef HEAP_COMPARE

//printf("~1 ierr sorted list:\n");
//for (i = 0; i < OPOINTS; i++)
//	printf("~1 %d: loc %f %f ierr %f\n", i, fp[i].p[0],fp[i].p[1],fp[i].v[0]);

	/* Return the best tnn as next points */
	for (j = i = 0; j < tnn && i < opoints; i++) {
		if (fp[i].v[0] < 0.0)
			continue;					/* Skip marked points */
		for (e = 0; e < di; e++)
			p[j].p[e] = fp[i].p[e];
		s->percept(s->od, p[j].v, p[j].p);
		j++;
	}
//printf("~1 returning %d points\n",j);
	return j;
}


/* --------------------------------------------------- */

/* determine the errors between the rspl and 100000 random test points */
static void
ppoint_stats(
ppoint *s
) {
	int i, n;
	int e, di = s->di;
	double mx = -1e80, av = 0.0, mn = 1e80;

	for (i = n = 0; i < 100000; i++) {
		co tp;				/* Perceptual value of test point */
		double vf[MXPD];	/* Perceptual value of reference */
		double sum, err;

		for (sum = 0.0, e = 0; e < di; e++)
			sum += tp.p[e] = d_rand(0.0, 1.0);

		if (sum <= s->ilimit) {

			/* rspl estimate of expected profile interpolation */
			s->g->interp(s->g, &tp);

			/* Target values */
			s->percept(s->od, vf, tp.p);

			for (err = 0.0, e = 0; e < di; e++) {
				double tt = tp.v[e] - vf[e]; 
				err += tt * tt;
			}
			err = sqrt(err);
			if (err > mx)
				mx = err;
			if (err < mn)
				mn = err;
			av += err;
			n++;
		}
	}
	av /= (double)n;

	printf("~1 Random check errors max %f, avg %f, min %f\n",mx,av,mn);
}

/* --------------------------------------------------- */
/* Support for maintaining the device/perceptual distance grid */
/* as well as keeping the far point candidate list up to date. */

/* Structure to hold data for callback function */
struct _pdatas {
	ppoint *s;			/* ppoint structure */
	int init;			/* Initialisation flag */
	co *pp;				/* List of new points */
	int nn;				/* Number of points */
}; typedef struct _pdatas pdatas;

/* rspl set callback function for maintaining perceptual distance information */
static void
pdfunc1(
	void *ctx,			/* Context */
	double *out,		/* output value, = di percept + distance */
	double *in			/* inut value */
) {
	pdatas *pp = (pdatas *)ctx;
	ppoint *s = pp->s;
	int e, di = s->di;

	if (pp->init) {
		s->percept(s->od, out, in);		/* Lookup perceptual value */
		out[di] = 10000.0;				/* Set to very high distance */

	} else {	/* Adding some points */
		int i;
		double sd = 1e80;

		/* Find smallest distance from this grid point to any of the new points */
		for (i = 0; i < pp->nn; i++) {
			double ddist, pdist;
			double dist;			/* Combined distance */

			/* Compute device and perceptual distance */
			for (ddist = pdist = 0.0, e = 0; e < di; e++) {
				double tt = out[e] - pp->pp[i].v[e];
				pdist += tt * tt;
				tt = 100.0 * (in[e] - pp->pp[i].p[e]);
				ddist += tt * tt;
			}
			dist = DDMIX * ddist + (1.0-DDMIX) * pdist;	/* Combine both */
			if (dist < sd)
				sd = dist;
		}

		sd = sqrt(sd);
		if (sd < out[di])
			out[di] = sd;

		/* Update far point candidate list */
		if (s->nfp < FPOINTS) {				/* List isn't full yet */
			for (e = 0; e < di; e++)
				s->fp[s->nfp].p[e] = in[e];
			s->fp[s->nfp].v[0] = sd;		/* store distance here */

			if (sd > s->wfpd) {				/* If this is the worst */
				s->wfpd = sd;
				s->wfp = s->nfp;
			}
			s->nfp++;

		} else if (sd < s->wfpd) {			/* Found better, replace current worst */

			for (e = 0; e < di; e++)
				s->fp[s->wfp].p[e] = in[e];
			s->fp[s->wfp].v[0] = sd;		/* store distance here */

			/* Locate the next worst */
			s->wfpd = -1.0;
			for (i = 0; i < s->nfp; i++) {
				if (s->fp[i].v[0] > s->wfp) {
					s->wfp = i;
					s->wfpd = s->fp[i].v[0];
				}
			}
		}
	}
}

/* Add a list of new points to the perceptual distance grid */
/* (Can change this to just adding 1 point) */
static void add_dist_points(
ppoint *s,
co *pp,		/* List of points including device and perceptual values */
int nn		/* Number in the list */
) {
	pdatas pdd;			/* pd callback context */

	pdd.s = s;
	pdd.init = 0;		/* Initialise values in the grid */
	pdd.pp = pp;
	pdd.nn = nn;

	/* let callback do all the work */
	s->pd->re_set_rspl(s->pd,
	           0,					/* No special flags */
	           &pdd,				/* Callback function context */
	           pdfunc1);			/* Callback function */
}

#ifdef NEVER	/* Not currently used */
/* Return the farthest distance value for this given location */
static double far_dist(ppoint *s, double *p) {
	int e, di = s->di;
	double cdist;
	co tp;
	
	for (e = 0; e < di; e++)
		tp.p[e] = p[e];

	s->pd->interp(s->pd, &tp);	/* Lookup current closest distance value */
	cdist = tp.v[di];
	if (cdist >= 10000.0)		/* Initial value */
		cdist = 0.0;
	return cdist;
}
#endif /* NEVER */

/* --------------------------------------------------- */
/* Seed the whole thing with points */

static void
ppoint_seed(
ppoint *s,
fxpos *fxlist,			/* List of existing fixed points */
int fxno				/* Number in fixed list */
) {
	int e, di = s->di;
	int i, j;

	if (fxno > 0) {
		co *pp;

		/* Place all the fixed points at the start of the list */
		if ((pp = (co *)malloc(fxno * sizeof(co))) == NULL)
			error ("ppoint: malloc failed on %d fixed nodes",fxno);

		for (i = 0; (i < fxno) && (i < s->tinp); i++) {
			node *p = &s->list[i];	/* Destination for point */

			for (e = 0; e < di; e++)
				p->p[e] = fxlist[i].p[e];

			p->fx = 1;			/* is a fixed point */
			s->percept(s->od, p->v, p->p);

			for (e = 0; e < di; e++) {
				pp[i].p[e] = p->p[e];
				pp[i].v[e] = p->v[e];
			}
		}
		s->np = s->fnp = i;

		/* Add new points to rspl interpolation */
		s->g->add_rspl(s->g, 0, pp, i);

		free(pp);
	}

	/* Seed the remainder points randomly */
	i = 0;
	while(s->np < s->tinp) {


#ifdef NEVER
		node *p = &s->list[s->np];
		double sum;

		/* Add random points */
		for (sum = 0.0, e = 0; e < di; e++)
			sum += p->p[e] = d_rand(0.0, 1.0);

		if (sum > s->ilimit)
			continue;
		s->np++;
		i++;
		printf("%cAdded: %d",cr_char,i);
#else 

#ifdef NEVER
		int nn;
		co pp[WPOINTS];		/* Space for return values */

		/* Add points at location with the largest error */
		nn = WPOINTS;

		if ((s->np + nn) > s->tinp)		/* Limit to desired value */
			nn = s->tinp - s->np;
		nn = ppoint_find_worst(s, pp, nn);

		/* Add new points to rspl interpolation and far field */
		s->g->add_rspl(s->g, 0, pp, nn);
		add_dist_points(s, pp, nn);
#else
		/* Diagnostic version */
		int nn;
		co pp[WPOINTS];		/* Space for return values */
		double err1[WPOINTS];
		double err2[WPOINTS];

		nn = WPOINTS;

		if ((s->np + nn) > s->tinp)		/* Limit to desired value */
			nn = s->tinp - s->np;
		nn = ppoint_find_worst(s, pp, nn);

		for (j = 0; j < nn; j++)
			err1[j] = ppoint_ierr(s, pp[j].p);

		/* Add new points to rspl interpolation and far field */
		s->g->add_rspl(s->g, 0, pp, nn);
		add_dist_points(s, pp, nn);

		for (j = 0; j < nn; j++)
			err2[j] = ppoint_ierr(s, pp[j].p);

		for (j = 0; j < nn; j++)
			printf("~1 improvement after adding point is %f to %f\n",err1[j],err2[j]);
#endif
		/* Copy points into ppoint */
		for (j = 0; j < nn; j++) {
			for (e = 0; e < di; e++) {
				s->list[s->np].p[e] = pp[j].p[e]; 
				s->list[s->np].v[e] = pp[j].v[e]; 
			}
			s->np++;
		}
		i += nn;
		printf("%cAdded: %d",cr_char,i);
#endif
	}
	printf("\n");		/* Finish "Added:" */
}

/* --------------------------------------------------- */

/* Rest the read index */
static void
ppoint_reset(ppoint *s) {
	s->rix = 0;
}

/* Read the next non-fixed point value */
/* Return nz if no more */
static int
ppoint_read(ppoint *s, double *p, double *f) {
	int e;

	/* Advance to next non-fixed point */
	while(s->rix < s->np && s->list[s->rix].fx)
		s->rix++;
	
	if (s->rix >= s->np)
		return 1;

	/* Return point info to caller */
	for (e = 0; e < s->di; e++) {
		if (p != NULL)
			p[e] = s->list[s->rix].p[e];
		if (f != NULL)
			f[e] = s->list[s->rix].v[e];
	}
	s->rix++;

	return 0;
}

/* Destroy ourselves */
static void
ppoint_del(ppoint *s) {

	/* Free our nodes */
	free(s->list);

	/* Free our rspl interpolation */
	s->g->del(s->g);

	/* Free our perceptual distance grid */
	s->pd->del(s->pd);

	free (s);
}

/* Creator */
ppoint *new_ppoint(
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
int tinp,				/* Total number of points to generate, including fixed */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
) {
	ppoint *s;

	// ~~~99 Info for logging
	fprintf(stderr, "WPOINTS = %d\n",WPOINTS);
	fprintf(stderr, "FPOINTS = %d\n",FPOINTS);
	fprintf(stderr, "OPOINTS = %d\n",OPOINTS);
	fprintf(stderr, "DDMIX   = %f\n",DDMIX);
	fprintf(stderr, "DWEIGHT = %f\n",DWEIGHT);
	fprintf(stderr, "CLOSED  = %f\n",CLOSED);

	if ((s = (ppoint *)calloc(sizeof(ppoint), 1)) == NULL)
		error ("ppoint: malloc failed");

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	if (di > MXPD)
		error ("ppoint: Can't handle di %d",di);

	s->di = di;

	if (tinp < fxno)	/* Make sure we return at least the fixed points */
		tinp = fxno;

	s->tinp = tinp;		/* Target total number of points */
	s->ilimit = ilimit;

	/* Init method pointers */
	s->reset = ppoint_reset;
	s->read  = ppoint_read;
	s->stats = ppoint_stats;
	s->del   = ppoint_del;

	/* If no perceptual function given, use default */
	if (percept == NULL) {
		s->percept = default_ppoint_to_percept;
		s->od = s;
	} else {
		s->percept = percept;
		s->od = od;
	}

	/* Allocate the list of points */
	s->np = 0;
	
	if ((s->list = (node *)calloc(sizeof(node), tinp)) == NULL)
		error ("ppoint: malloc failed on nodes");

	/* Setup the interpolation and perceptual distance rspls */
	{
		int e;
		int tres, gres[MXDI];
		datai pl,ph;
		datai vl,vh;
		double avgdev[MXDO];
		pdatas pdd;			/* pd callback context */

#ifndef NEVER	/* High res. */
		if (di <= 2)
			tres = 41;		/* Make depend on no points and dim ? */
		else if (di <= 3)
			tres = 33;		/* Make depend on no points and dim ? */
		else
			tres = 15;
#else
		if (di <= 2)
			tres = 3;		/* Make depend on no points and dim ? */
		else if (di <= 3)
			tres = 17;		/* Make depend on no points and dim ? */
		else
			tres = 9;
#endif

		/* The interpolation grid mimics the operation of the profile */
		/* package creating a device to CIE mapping for the device from */
		/* the given test points. */
		s->g = new_rspl(RSPL_NOFLAGS, di, di);

		for (e = 0; e < di; e++) {
			pl[e] = 0.0;
			ph[e] = 1.0;
			if (e == 1 || e == 2) {		/* Assume Lab */
				vl[e] = -128.0;
				vh[e] = 128.0;
			} else {
				vl[e] = 0.0;
				vh[e] = 100.0;
			}
			gres[e] = tres;
			avgdev[e] = 0.005;
		}

		/* Setup other details of rspl */
		s->g->fit_rspl(s->g,
		           RSPL_INCREMENTAL |
		           /* RSPL_EXTRAFIT | */	/* Extra fit flag */
		           0,
		           NULL,				/* No test points initialy */
		           0,					/* No test points */
		           pl, ph, gres,		/* Low, high, resolution of grid */
		           vl, vh,				/* Data scale */
		           0.3,					/* Smoothing */
		           avgdev,				/* Average Deviation */
		           NULL);


		/* Track closest perceptual distance to existing test points. */
		/* To save looking up the perceptual value for every grid location */
		/* every time a point is added, cache this values in the grid too. */
		s->pd = new_rspl(RSPL_NOFLAGS, di, di+1);

		/* Initialise the pd grid ready for the first points. */
		pdd.s = s;
		pdd.init = 1;		/* Initialise values in the grid */

		s->pd->set_rspl(s->pd,
		           0,					/* No special flags */
		           &pdd,				/* Callback function context */
		           pdfunc1,				/* Callback function */
		           pl, ph, gres,		/* Low, high, resolution of grid */
		           vl, vh);				/* Data scale */

		s->wfpd = -1.0;		/* Impossibly good worst point distance */
	}

	/* Create the points */
	ppoint_seed(s, fxlist, fxno);

	/* Print some stats */
	ppoint_stats(s);

	ppoint_reset(s);		/* Reset read index */

	return s;
}

/* =================================================== */

#ifdef STANDALONE_TEST

/* Graphics Gems curve */
static double gcurve(double vv, double g) {
	if (g >= 0.0) {
		vv = vv/(g - g * vv + 1.0);
	} else {
		vv = (vv - g * vv)/(1.0 - g * vv);
	}
	return vv;
}

#ifdef NEVER
static void sa_percept(void *od, double *out, double *in) {
	double lab[3];
	
	clu->dev_to_rLab(clu, lab, in);

	out[0] = lab[0];
//	out[1] = (lab[1]+100.0)/2.0;
	out[1] = (lab[2]+100.0)/2.0;
}
#else

static void sa_percept(void *od, double *p, double *d) {

#ifndef NEVER
	/* Default Do nothing - copy device to perceptual. */
	p[0] = 100.0 * gcurve(d[0], -4.5);
	p[1] = 100.0 * gcurve(d[1], 2.8);
	p[1] = 0.8 * p[1] + 0.2 * p[0];
#else
	for (e = 0; e < di; e++) {
		double tt = d[e];
		/* Two slopes with a sharp turnover in X */
		if (e == 0) {
			if (tt < 0.5)
				tt = tt * 0.3/0.5;
			else
				tt = 0.3 + ((tt-0.5) * 0.7/0.5);
		}
		p[e] = tt * 100.0;
	}
#endif
}
#endif


int
main(argc,argv)
int argc;
char *argv[];
{
	int npoints = 21;
	ppoint *s;
	long stime,ttime;
	error_program = argv[0];

	printf("Standalone test of ppoint, argument is number of points, default %d\n",npoints);

	if (argc > 1)
		npoints = atoi(argv[1]);

	/* Create the required points */
	stime = clock();
	s = new_ppoint(2, 1.5, npoints, NULL, 0, sa_percept, (void *)NULL);

	ttime = clock() - stime;
	printf("Execution time = %f seconds\n",ttime/(double)CLOCKS_PER_SEC);

#ifdef DUMP_PLOT
	printf("Perceptual plot:\n");
	dump_image(s, 1);

	printf("Device plot:\n");
	dump_image(s, 0);
#endif /* DUMP_PLOT */

	s->del(s);

	return 0;
}

#ifdef NEVER
/* Basic printf type error() and warning() routines */
#ifdef	__STDC__
void
error(char *fmt, ...)
#else
void
error(va_alist) 
va_dcl
#endif
{
	va_list args;
#ifndef	__STDC__
	char *fmt;
#endif

	fprintf(stderr,"ppoint: Error - ");
#ifdef	__STDC__
	va_start(args, fmt);
#else
	va_start(args);
	fmt = va_arg(args, char *);
#endif
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	fflush(stdout);
	exit (-1);
}
#endif /* NEVER */
#endif /* STANDALONE_TEST */


#ifdef STANDALONE_TEST
#ifdef DUMP_PLOT

/* Dump the current point positions to a plot window file */
void
static dump_image(ppoint *s, int pcp) {
	int i;
	double minx, miny, maxx, maxy;
	static double *x1a = NULL;
	static double *y1a = NULL;

	if (pcp != 0) {	/* Perceptual range */
		minx = 0.0;	/* Assume */
		maxx = 100.0;
		miny = 0.0;
		maxy = 100.0;
	} else {
		minx = 0.0;	/* Assume */
		miny = 0.0;
		maxx = 1.0;
		maxy = 1.0;
	}
	
	if (x1a == NULL) {
		if ((x1a = (double *)malloc(s->np * sizeof(double))) == NULL)
			error ("ppoint: malloc failed");
		if ((y1a = (double *)malloc(s->np * sizeof(double))) == NULL)
			error ("ppoint: malloc failed");
	}

	for (i = 0; i < s->np; i++) {
		node *p = &s->list[i];
		
		if (pcp != 0) {
			x1a[i] = p->v[0];
			y1a[i] = p->v[1];
		} else {
			x1a[i] = p->p[0];
			y1a[i] = p->p[1];
		}
	}

	/* Plot the vectors */
	do_plot_vec(minx, maxx, miny, maxy, 
				x1a, y1a, x1a, y1a, s->np, DO_WAIT, NULL, NULL, NULL, NULL, 0);
}

#endif /* DUMP_PLOT */
#endif /* STANDALONE_TEST */















