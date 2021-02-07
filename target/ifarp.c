
/* 
 * Argyll Color Correction System
 *
 * Incremental far point class
 *
 * Author: Graeme W. Gill
 * Date:   6/11/2002
 *
 * Copyright 2002 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
   Algorithm:

   Starting with a previous test point as a seed, use a random starte point
   and minimisation algorithm to locate another point that is as far as
   possible from the nearest existing test point in perceptual space,
   while remaining in gamut at all tines. This means that ideally each
   point "fills in" the gaps in the existing distribution, while starting
   from an existing point.

   The performance is still not very good, as the inner loop involves
   locating the nearest existing point, as well as converting from
   device coordinates to perceptual space. If the powell search radius
   is reduced too much the uniformity of the distribution suffers.

 */

/* TTBD:

	It would probably help the uniformity of distribution if we could
    aproximately locate the next seed point as the one with the
    biggest adjoing "gap", and this may speed things up by allowing us
    to reduce the powel search radius. 

	Perhaps switching to a balltree indexing structure would speed up
    nearest ppoint finding as well as providing a mechanism to quickly
    locate the nearest "void".

	Subsequent experience indicates that furthest distance in perceptual
	space may not be the best strategy, but furthest distance in device
	space may be. Add #define allowing this to be tested ?? 

 */

#undef DEBUG
#define PERC_PLOT 1		/* Emit perceptive space plots (if DEBUG) */
#define DO_WAIT 1		/* Wait for user key after each plot */

#define ASSERTS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#ifdef DEBUG
# include "plot.h"
# include "ui.h"
#endif
#include "numlib.h"
#include "sort.h"
#include "icc.h"
#include "xcolorants.h"
#include "targen.h"
#include "ifarp.h"
#include "sort.h"		/* Heap sort */

#ifdef DEBUG
static void dump_image(ifarp *s, int pcp);
static void dump_image_final(ifarp *s, int pcp);
#endif

#define MAX_TRIES   30		/* Maximum itterations */


/* nn functions */
static double nearest(ifarp *s, double *q);
static void init_nn(ifarp *s);
static void add_nn(ifarp *s);
static void del_nn(ifarp *s);

/* ----------------------------------------------------- */
/* Default convert the nodes device coordinates into approximate perceptual coordinates */
static void
ifarp_to_percept(void *od, double *p, double *d) {
	ifarp *s = (ifarp *)od;
	int e;

	/* Do nothing - copy device to perceptual. */
	for (e = 0; e < s->di; e++) {
		p[e] = d[e] * 100.0;
	}
}


/* Return the largest distance of the point outside the device gamut. */
/* This will be 0 if inside the gamut, and > 0 if outside.  */
static double
ifarp_in_dev_gamut(ifarp *s, double *d) {
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

/* Snap a point to the device gamut boundary. */
/* Return nz if it has been snapped. */
static int snap_to_gamut(ifarp *s, double *d) {
	int e;
	int di = s->di;
	double dd;			/* Smallest distance */
	double ss;			/* Sum */
	int rv = 0;

	/* Snap to ink limit first */
	for (ss = 0.0, e = 0; e < di; e++)
		ss += d[e];
	dd = fabs(ss - s->ilimit);

	if (dd < 0.0) {
		int j;
		for (j = 0; j < di; j++) 
			d[j] *= s->ilimit/ss;	/* Snap to ink limit */
		rv = 1;
	}

	/* Now snap to any other dimension */
	for (e = 0; e < di; e++) {

		dd = fabs(d[e] - 0.0);
		if (dd < 0.0) {
			d[e] = 0.0;		/* Snap to orthogonal boundary */
			rv = 1;
		}
		dd = fabs(1.0 - d[e]); 
		if (dd < 0.0) {
			d[e] = 1.0;		/* Snap to orthogonal boundary */
			rv = 1;
		}
	}

	return rv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Reverse lookup function :- perceptual to device coordinates */

/* Definition of the optimization functions handed to powell() */

/* Return metric to be minimised, and */
/* an error >= 50000 on being out of device gamut */
static double efunc(void *edata, double p[]) {
	ifarp *s = (ifarp *)edata;
	double rv;
	if ((rv = (ifarp_in_dev_gamut(s, p))) > 0.0) {
		rv = rv * 500.0 + 500.0;		/* Discourage being out of gamut */
	} else {
		double v[MXTD];
		s->percept(s->od, v, p);
		rv = 500.0 - nearest(s, v);
	}
//printf("~1 rv = %f from %f %f\n",rv,p[0],p[1]);
	return rv;
}

/* Given a point in device space, optimise it to be */
/* within the device gamut, as well as being as far as */
/* possible from the nearest point in perceptual space. */
/* return nz if powell failed */
static int
optimise_point(
ifarp *s,
double *d		/* starting and returned device position */
) {
	int e, di = s->di;
	double sr[MXTD];	/* Search radius in each device dimension */
	double drad = 1.0;	/* Search Radius (affects fill evenness) */
	double ptol = 0.001;	/* Tolerance */
	double tt;

// ~~99
	for (e = 0; e < di; e++)
		sr[e] = drad;			/* Device space search radius */
	if (powell(&tt, di, d, sr,  ptol, 500, efunc, (void *)s, NULL, NULL) != 0 || tt >= 50000.0) {
#ifdef DEBUG
		warning("ifarp: powell failed, tt = %f",tt);
#endif
		return 1;
	}
	snap_to_gamut(s, d);
	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Create a new node. */
/* Return current number of nodes */
static int new_node(
ifarp *s,
int ix		/* Index of point to start from */
) {
	int di = s->di;
	int e;

// ~~99
	/* Retry if powell failes */
	for (;;) {

		/* Create the new point by cloning the existing point */
		s->nodes[s->np].fx = 0;			/* Not a fixed/pre-existing node */
		for (e = 0; e < di; e++) {
			s->nodes[s->np].p[e] = s->nodes[ix].p[e];
		}
		/* Compute new point location that is farthest from nearest existing point */
		if (optimise_point(s, s->nodes[s->np].p) == 0)
			break;
	}

	/* compute perceptual location */
	s->percept(s->od, s->nodes[s->np].v, s->nodes[s->np].p);

#ifdef DEBUG
printf("Added node %d at perc %f %f, dev %f %f\n",
s->np,
s->nodes[s->np].v[0],
s->nodes[s->np].v[1],
s->nodes[s->np].p[0],
s->nodes[s->np].p[1]);
#endif

	/* Add the node to our current list */
	s->nodes[s->np].touch = s->tbase;
	s->np++;
	add_nn(s);

	return s->np;
}

/* ============================================= */
/* Main object functions */

/* Initialise, ready to read out all the points */
static void ifarp_reset(ifarp *s) {
	s->rix = 0;
}

/* Read the next set of non-fixed points values */
/* return non-zero when no more points */
static int ifarp_read(
ifarp *s,
double *d,		/* Device position */
double *p		/* Perceptual value */
) {
	int j;

	for (; s->rix < s->np; s->rix++) {

		if (s->nodes[s->rix].fx == 0) {
			for (j = 0; j < s->di; j++) {
				if (d != NULL)
					d[j] = s->nodes[s->rix].p[j];
				if (p != NULL)
					p[j] = s->nodes[s->rix].v[j];
			}
			s->rix++;
			return 0;
		}
	}
	return 1;
}

/* Destroy ourselves */
static void
ifarp_del(ifarp *s) {

	if (s->nodes != NULL)
		free(s->nodes);

	free (s);
}

/* Constructor */
ifarp *new_ifarp(
int verb,				/* Verbosity */
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
int inp,				/* Number of points to generate */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
) {
	ifarp *s;
	int e, i;

#ifdef DEBUG
	printf("new_ifarp called with di %d, inp %d, fxno = %d\n",di,inp,fxno);
#endif

	if ((s = (ifarp *)calloc(sizeof(ifarp), 1)) == NULL)
		error ("ifarp: ifarp malloc failed");

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	s->reset = ifarp_reset;
	s->read  = ifarp_read;
	s->del   = ifarp_del;

	/* If no perceptual function given, use default */
	if (percept == NULL) {
		s->percept = ifarp_to_percept;
		s->od = s;
	} else {
		s->percept = percept;
		s->od = od;
	}

	s->ilimit = ilimit;

	s->inp = inp;	/* Intended number of points */
	s->np = 0;

	if (di > MXTD)
		error ("ifarp: Can't handle di %d",di);
	s->di = di;
	s->tbase = 0;
	
	/* Initial alloc of nodes */
	if ((s->nodes = (ifpnode *)malloc(s->inp * sizeof(ifpnode))) == NULL)
		error ("ifarp: nodes malloc failed");

	/* Copy fixed nodes */
	for (i = 0; (i < fxno) && (s->np < s->inp); i++) {
		s->nodes[s->np].fx = 1;
		for (e = 0; e < di; e++)
			s->nodes[s->np].p[e] = fxlist[i].p[e];
		s->percept(s->od, s->nodes[i].v, s->nodes[i].p);
		s->nodes[s->np].touch = s->tbase;
		s->np++;
	}

	/* Create at least one seed point */
	if (s->np == 0) {
		s->nodes[s->np].fx = 0;

		for (e = 0; e < di; e++)
			s->nodes[s->np].p[e] = 0.0;		/* This is assumed to be in gamut */
		s->percept(s->od, s->nodes[i].v, s->nodes[i].p);
		s->nodes[s->np].touch = s->tbase;
		s->np++;
	}

	/* Setup initial nearest point acceleration structure */
	init_nn(s);

	/* Create initial patches */
// ~~99

	if (verb)
		printf("Full points:\n");

	for (i = 0; s->np < s->inp; i += 17) {
		i %= s->np;
		new_node(s, i);
		if (verb) {
			int pc = (int)(100.0 * s->np/s->inp + 0.5);
			printf("  % 3d%%%c",pc,cr_char); fflush(stdout);
		}
	}

	if (verb)
		printf("\n");

	/* We're done with acceleration structure */
	del_nn(s);

	return s;
}

/* --------------------------------------------------- */
/* (This code is has been copied from gamut/gamut.c) */

#ifdef DEBUG
#define NN_INF 100000.0
#else
#define NN_INF 1e307
#endif

/* Given a point, */
/* return the nearest existint test point. */
static double
nearest(
ifarp *s,
double *q		/* Target point location */
) {
	int e, i, k;
	int di = s->di;
	int wex[MXTD * 2];		/* Current window edge indexes */
	double wed[MXTD * 2];	/* Current window edge distances */
							/* Indexes are axis * 2 +0 for lower edge, */
							/* +1 for upper edge of search box. */
							/* We are comparing lower edge of search box */
							/* with upper edge of bounding box etc. */ 

//printf("~1 nearest called\n");

	/* We have to find out which existing point the point will be nearest */

	if ((s->tbase + di) < s->tbase) {	/* Overflow of touch count */
		for (i = 0; i < s->np; i++)
			s->sax[0][i]->touch = 0;		/* reset it in all the objects */
		s->tbase = 0;
	}
	s->ttarget = s->tbase + di;		/* Target touch value */

//printf("\n");
//printf("Query point is %f %f\n",q[0], q[1]);

	/* Find starting indexes within axis arrays */
	for (e = 0; e < (2 * di); e++) {	/* For all axes min & max */
		int f = e/2;			/* Axis */
		int i0, i1, i2;			/* Search indexes */
		double v0, v1, v2;		/* Box */
		double qf, ww;

		/* Binary search this edge */
		qf = q[f]; 		/* strength reduced q[f] */

//printf("\n");
//printf("isearching axis %d %s for %f\n",f, e & 1 ? "max" : "min", qf);
		i0 = 0;
		i2 = s->np - 1;
		v0 = s->sax[f][i0]->v[f];
		v2 = s->sax[f][i2]->v[f];
//printf("start points %d - %d, bound %f - %f\n",i0, i2, v0, v2);

		if (qf <= v0) {
			i2 = i0;
			v2 = v0;
		} else if (qf >= v2) {
			i0 = i2;
			v0 = v2;
		} else {
			do {
				i1 = (i2 + i0)/2;		/* Trial point */
				v1 = s->sax[f][i1]->v[f];	/* Value at trial */
				if (v1 < qf) {
					i0 = i1;			/* Take top half */
					v0 = v1;
				} else {
					i2 = i1;			/* Take bottom half */
					v2 = v1;
				}
//printf("current point %d - %d, bound %f - %f\n",i0, i2, v0, v2);
			} while ((i2 - i0) > 1);
		}

		if (e & 1) {			/* Max side of window */
			int tc;				/* total object count */

			ww = v2 - qf;
			wed[e] = fabs(ww) * ww;
			wex[e] = i2;

			/* Check that min and max together will cover at least s->np objects */
			tc = s->np - i2 + wex[e ^ 1] + 1;
//printf("got %d, expected %d\n",tc, s->np);

			/* (I don't really understand why this works!) */
			if (tc < s->np) {		/* We haven't accounted for all the objects */
				int el = e ^ 1;		/* Low side sax */
				int ti0, ti2;
				double tv0, tv2;

				ti0 = wex[el];
				ti2 = i2;
//printf("We have straddling objects, initial indexes are %d - %d\n",ti0, ti2);

				/* While straddling objects remain undiscovered: */
				while (tc < s->np) {
					tv0 =  NN_INF;		/* Guard values */
					tv2 = -NN_INF;

					/* Increment low side until we find a straddler */
					while (ti0 < (s->np-1)) {
						ww = s->sax[f][++ti0]->v[f];	/* Position of the other end */
						if (ww < qf) {
//printf("found low object %d at index %d that straddles\n",s->sax[f][ti0]-s->nodes,ti0);
							tv0 = qf - s->sax[f][ti0]->v[f];
							break;
						}
					}

					/* Decrement high side until we find a straddler */
					while (ti2 > 0) {
						ww = s->sax[f][--ti2]->v[f];	/* Position of the other end */
						if (ww > qf) {
//printf("found high object %d at index %d that straddles\n",s->sax[f][ti2]-s->nodes,ti2);
							tv2 = s->sax[f][ti2]->v[f] - qf;
							break;
						}
					}
					/* Choose the closest */
					if (tv0 > tv2) {
						wed[el] = fabs(tv0) * tv0;
						wex[el] = ti0;
						tc++;
					} else {
						wed[e] = fabs(tv2) * tv2;
						wex[e] = ti2;
						tc++;
					}
				}
//printf("After correction we have %d - %d\n",wex[e^1], wex[e]);
			}
		} else {				/* Min side of window */
			ww = q[f] - v0;
			wed[e] = fabs(ww) * ww;
			wex[e] = i0;
		}
	}

	/* Expand a di dimenstional cube centered on the target point, */
	/* jumping to the next nearest point on any axis, discovering */
	/* any bounding boxes that are within the expanding window */
	/* by checking their touch count. */

	/* The first point found establishes the initial best distance. */
	/* When the window expands beyond the point where it can improve */
	/* the best distance, stop */

	{
		double bw = 0.0;		/* Current window distance */
		double bdist = NN_INF;	/* Best possible distance to an object outside the window */
		int bix;				/* Index of best point */

		/* Until we're done */
		for (;;) {
			int ee;				/* Axis & expanding box edge */
			int ff;				/* Axis */
			int ii;				/* Index of chosen point */
			ifpnode *ob;		/* Current object */
			unsigned int ctv;	/* Current touch value */
//printf("\n");
//printf("wwidth = %f, bdist = %f, window = %d-%d, %d-%d\n",
//bw, bdist, wex[0], wex[1], wex[2], wex[3]);
//printf("window edge distances are = %f-%f, %f-%f\n",
//wed[0], wed[1], wed[2], wed[3]);

			/* find next (smallest) window increment axis and direction */
			ee = 0;
			ii = wex[ee];
			bw = wed[ee];
			for (e = 1; e < (2 * di); e++) {
				if (wed[e] < bw) {
					ee = e;
					ii = wex[e];
					bw = wed[e];
				}
			}
//printf("Next best is axisdir %d, object %d, axis index %d, best possible dist %f\n",
//ee, s->sax[ee/2][ii] - s->nodes, ii, bw);

			if (bw == NN_INF || bw > bdist) {
				break;		/* Can't go any further, or further points will be worse */
			}

#ifdef ASSERTS
if (ii < 0 || ii >= s->np) {
printf("Assert: went out of bounds of sorted axis array\n");
exit(0);
}
#endif
			/* Chosen point on ee axis/direction, index ii */
			ff = ee / 2;			/* Axis only */

			ob = s->sax[ff][ii];

			/* Touch value of current object */
			ctv = ob->touch;

			if (ctv < s->ttarget) {		/* Not been dealt with before */

				/* Touch this new window boundary point */
				ob->touch = ctv = ((ctv < s->tbase) ? s->tbase : ctv) + 1;

//printf("New touch count on %d is %d, target %d\n",
//ob - s->nodes, s->sax[ff][ii]->touch, s->ttarget);

				/* Check the point out */
				if (ctv == (s->tbase + di)) {	/* Is within window on all axes */
					double tdist = 0.0;

					/* Compute distance from query point to this object */
					for (k = 0; k < di; k++) {
						double tt = ob->v[k] - q[k];
						tdist += tt * tt;
					}
					
//printf("Got new best point %d, dist %f\n",ob-s->nodes,sqrt(tdist));
					if (tdist < bdist) {	/* New closest distance */
						bdist = tdist;
						bix = ob - s->nodes;
					}
				}
			}

			/* Increment next window edge candidate, and figure new edge distance */
			if (ee & 1) {					/* Top */
				if (++wex[ee] >= s->np) {
					wed[ee] = NN_INF;
					wex[ee]--;
				} else {
					double ww = s->sax[ff][wex[ee]]->v[ff] - q[ff];
					wed[ee] = fabs(ww) * ww;
				}
			} else {
				if (--wex[ee] < 0) {
					wed[ee] = NN_INF;
					wex[ee]++;
				} else {
					double ww = q[ff] - s->sax[ff][wex[ee]]->v[ff];
					wed[ee] = fabs(ww) * ww;
				}
			}
		}

		s->tbase += di;			/* Next touch */

//printf("~1 returning closest to node %d distance %f\n",bix,sqrt(bdist));
		return sqrt(bdist);	/* Return nearest distance */
	}
}


/* Setup the nearest function acceleration structure */
/* with the existing points */
static void
init_nn(
ifarp *s
) {
	int di = s->di;
	int i, k;
	int np = s->np;		/* Existing number of points */

//printf("~9 init_nn called\n");

	s->tbase = 0;		/* Initialse touch flag */

	/* Allocate the arrays spaces for intended number of points */
	for (k = 0; k < di; k++) {
		if ((s->sax[k] = (ifpnode **)malloc(sizeof(ifpnode *) * s->inp)) == NULL)
			error("Failed to allocate sorted index array");
	}

	/* Add each existing test point to the axis lists. */
	for (i = 0; i < np; i++) {
		for (k = 0; k < di; k++)
			s->sax[k][i] = &s->nodes[i];
	}

	/* Sort the axis arrays */
	for (k = 0; k < di; k++) {
			/* Sort nodes edge list */
#define 	HEAP_COMPARE(A,B) (A->v[k] < B->v[k])
			HEAPSORT(ifpnode *, &s->sax[k][0], np)
#undef HEAP_COMPARE
	}
//printf("~9 init_nn done\n");
}


#ifdef NEVER		/* Slower but simpler version */

/* Add the last point to the acceleration structure */
static void
add_nn(
ifarp *s
) {
	int di = s->di;
	int i, k;
	int np = s->np;		/* Existing number of points */
	int ap = np - 1;	/* Index of point ot add */

//printf("~9 add_nn called with point ix %d, pos %f %f\n",ap, s->nodes[ap].v[0],s->nodes[ap].v[1]);

	for (k = 0; k < di; k++) {
		s->sax[k][ap] = &s->nodes[ap];
	}

	/* Sort the axis arrays */
	for (k = 0; k < di; k++) {
			/* Sort nodes edge list */
#define 	HEAP_COMPARE(A,B) (A->v[k] < B->v[k])
			HEAPSORT(ifpnode *, &s->sax[k][0], np)
#undef HEAP_COMPARE
	}
}

#else

/* Add the last point to the acceleration structure */
static void
add_nn(
ifarp *s
) {
	int di = s->di;
	int e;
	int np = s->np;		/* Existing number of points */
	int ap = np - 1;	/* Index of point to add */

//printf("~9 add_nn called with point ix %d, pos %f %f\n",ap, s->nodes[ap].v[0],s->nodes[ap].v[1]);

	for (e = 0; e < di; e++) {	/* For all axes */
		int i0, i1, i2;			/* Search indexes */
		double v0, v1, v2;		/* Box */
		double qf;

		qf = s->nodes[ap].v[e];	/* value to be insertion sorted */

//printf("isearching axis %d for %f\n",e, qf);

		/* Find index of lowest value that is greater than target */
		i0 = 0;
		i2 = ap - 1;
		v0 = s->sax[e][i0]->v[e];
		v2 = s->sax[e][i2]->v[e];
//printf("start points %d - %d, bound %f - %f\n",i0, i2, v0, v2);

		if (qf <= v0) {
			i1 = i0;
		} else if (qf >= v2) {
			i1 = ap;
		} else {
			do {
				i1 = (i2 + i0)/2;			/* Trial point */
				v1 = s->sax[e][i1]->v[e];	/* Value at trial */

				if (qf > v1) {
					i0 = i1;			/* Take top half */
					v0 = v1;
				} else { /* qf <= v1 */
					i2 = i1;			/* Take bottom half */
					v2 = v1;
				}
//printf("current point %d - %d, bound %f - %f\n",i0, i2, v0, v2);
			} while ((i2 - i0) > 1);

			i1 = i0;
			v1 = s->sax[e][i1]->v[e];

			/* Ensure we're > than target */
			while (v1 <= qf) {
				i1++;
				if (i1 < ap)
					v1 = s->sax[e][i1]->v[e];
				else 
					break;
			}
		}

		/* Make room */
		if (i1 < ap) {
			memmove((void *)&s->sax[e][i1+1], (void *)&s->sax[e][i1], (ap - i1) * sizeof(ifpnode *));
		}
		/* Insert */
		s->sax[e][i1] = &s->nodes[ap];
	}
}

#endif

/* Free everything to do with the nn */
static void del_nn(ifarp *s) {
	int di = s->di;
	int k;

	for (k = 0; k < di; k++) {
		free (s->sax[k]);
	}
}

/* =================================================== */

#ifdef STANDALONE_TEST

icxColorantLu *clu;

void sa_percept(void *od, double *out, double *in) {
	
#ifdef NEVER
	double lab[3];
	clu->dev_to_rLab(clu, lab, in);

	out[0] = lab[0];
//	out[1] = (lab[1]+100.0)/2.0;
	out[1] = (lab[2]+100.0)/2.0;
#else

	out[0] = in[0] * 100.0;
	out[1] = in[1] * 100.0;

#endif
}

int
main(argc,argv)
int argc;
char *argv[];
{
	int npoints = 500;
	ifarp *s;
	int mask = ICX_BLACK | ICX_GREEN;
	int di = 2;
	
	error_program = argv[0];
	check_if_not_interactive();

	if (argc > 1)
		npoints = atoi(argv[1]);

	if ((clu = new_icxColorantLu(mask)) == NULL)
		error ("Creation of xcolorant lu object failed");

	/* Create the required points */
	s = new_ifarp(1, di, 1.5, npoints, NULL, 0, sa_percept, (void *)NULL);

#ifdef DEBUG
	/* Dump perceptual map */
	dump_image(s, PERC_PLOT);
#endif /* DEBUG */

	s->del(s);

	return 0;
}

#endif /* STANDALONE_TEST */



#ifdef DEBUG

/* Dump the current point positions to a plot window file */
static void
dump_image(ifarp *s, int pcp) {
	double minx, miny, maxx, maxy;
	double *x1a = NULL;
	double *y1a = NULL;
	double *x2a = NULL;
	double *y2a = NULL;
	double *x3a = NULL;
	double *y3a = NULL;

	int i, nu;
	ifpnode *p;

	if (s->np == 0)
		return;

	if (pcp) {	/* Perceptual range */
		minx = 0.0;	/* Assume */
		miny = 0.0;
		maxx = 100.0;
		maxy = 100.0;
	} else {
		minx = 0.0;	/* Assume */
		miny = 0.0;
		maxx = 1.0;
		maxy = 1.0;
	}
	
	if ((x1a = (double *)malloc(s->np * sizeof(double))) == NULL)
		error ("ifarp: plot malloc failed %d",s->np);
	if ((y1a = (double *)malloc(s->np * sizeof(double))) == NULL)
		error ("ifarp: plot malloc failed %d",s->np);
	if ((x2a = (double *)malloc(s->np * sizeof(double))) == NULL)
		error ("ifarp: plot malloc failed %d",s->np);
	if ((y2a = (double *)malloc(s->np * sizeof(double))) == NULL)
		error ("ifarp: plot malloc failed %d",s->np);

	for (nu = i = 0; i < s->np; i++) {
		p = &s->nodes[i];

		if (pcp) {
			x1a[nu] = p->v[0];
			y1a[nu] = p->v[1];
			x2a[nu] = p->v[0];
			y2a[nu] = p->v[1];
		} else {
			x1a[nu] = p->p[0];
			y1a[nu] = p->p[1];
			x2a[nu] = p->p[0];
			y2a[nu] = p->p[1];
		}
		nu++;
	}

	/* Plot the vectors */
	do_plot_vec(minx, maxx, miny, maxy, 
				x1a, y1a, x2a, y2a, nu, DO_WAIT, x3a, y3a, 0);

	free(x1a);
	free(y1a);
	free(x2a);
	free(y2a);
}

#endif /* DEBUG */





