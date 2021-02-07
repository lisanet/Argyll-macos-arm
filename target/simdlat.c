
/* 
 * Argyll Color Correction System
 *
 * Device space latice test point generator class,
 * set to generate a body centered cubic lattice.
 *
 * Author: Graeme W. Gill
 * Date:   30/8/2004
 *
 * Copyright 2002 - 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Based on simplat.c
 */

/* TTBD:

	This is not the most efficient way to generate a body centered
	cubic lattice, a simpler grid counter would work faster, without
	having to track what points have been generated.

	This seems too inexact/slow to read a specified number of test
	points for use in higher dimensions.

 */

#undef DEBUG
#undef DUMP_PLOT		/* Show on screen plot */
#define PERC_PLOT 0		/* Emit perceptive space plots */
#define DO_WAIT 1		/* Wait for user key after each plot */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#if defined(DEBUG) || defined(DUMP_PLOT)
# include "plot.h"
# include "ui.h"
#endif
#include "numlib.h"
#include "sort.h"
#include "icc.h"
#include "xcolorants.h"
#include "targen.h"
#include "simdlat.h"

#if defined(DEBUG) || defined(DUMP_PLOT)
static void dump_image(simdlat *s, int pcp);
static void dump_image_final(simdlat *s, int pcp);
#endif

#define SNAP_TOL	0.02	/* Snap to gamut boundary tollerance */
#define MAX_TRIES   30		/* Maximum itterations */


/* ----------------------------------------------------- */

/* Default convert the nodes device coordinates into approximate perceptual coordinates */
/* (usually overriden by caller supplied function) */
static void
default_simdlat_to_percept(void *od, double *p, double *d) {
	simdlat *s = (simdlat *)od;
	int e;

	/* Default Do nothing - copy device to perceptual. */
	for (e = 0; e < s->di; e++) {
		p[e] = d[e] * 100.0;
	}
}


#ifdef NEVER /* Not currently used */
/* Return the largest distance of the point outside the device gamut. */
/* This will be 0 if inside the gamut, and > 0 if outside.  */
static double
simdlat_in_dev_gamut(simdlat *s, double *d) {
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
#endif /* NEVER */

/* Snap a point to the device gamut boundary. */
/* Return nz if it has been snapped. */
static int snap_to_gamut(simdlat *s, double *d) {
	int e;
	int di = s->di;
	double dd;			/* Smallest distance */
	double ss;			/* Sum */
	int rv = 0;

	/* Snap to ink limit first */
	for (ss = 0.0, e = 0; e < di; e++)
		ss += d[e];
	dd = ss - s->ilimit;

	if (dd >= -s->tol) {	/* Within tol or beyond limit */
		int j;
		for (j = 0; j < di; j++) 
			d[j] *= s->ilimit/ss;	/* Snap to ink limit */
		rv = 1;
	}

	/* Now snap to any other dimension */
	for (e = 0; e < di; e++) {

		dd = 0.0 - d[e];
		if (dd >= -s->tol) {
			d[e] = 0.0;		/* Snap to orthogonal boundary */
			rv = 1;
		}
		dd = d[e] - 1.0; 
		if (dd >= -s->tol) {
			d[e] = 1.0;		/* Snap to orthogonal boundary */
			rv = 1;
		}
	}

	return rv;
}

/* Snap a point to the gamut boundary if it is close enough. */
/* Return 1 if the point has been clipped. */
/* Return 2 if the point has been clipped by a dia. */
static int snap_to_gamut2(simdlat *s, double *d) {
	int rv = 0;
	double ud[MXTD];
	int e, di = s->di;

//printf("\n~1 snap_to_gamut2() called with %f %f\n",d[0],d[1]);

	for (e = 0; e < di; e++)
		ud[e] = d[e];			/* save unclipped location */ 

	if (snap_to_gamut(s, d)) {
		double tt;
	
		tt = 0.0;
		for (e = 0; e < di; e++) {
			double t = ud[e] - d[e];
			tt += t * t;
		}
		tt = sqrt(tt);
//printf("~1 Got snapped to %f %f by dist %f\n",d[0], d[1], tt);
		if (tt > (0.5 * s->dia))
			rv = 2;				/* invalid & !explore */
		else
			rv = 1;				/* Valid & explore */
	}
//else
//printf("~1 Didn't get snapped\n");

//printf("~1 snap_to_gamut2() on %f %f returning %d\n",d[0],d[1],rv);
	return rv;
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Compute the equalateral simplex basis vectors */
static void comp_basis(
simdlat *s,
int type,		/* 0 = body centered cubic */
				/* 1 = equilateral simplex */
				/* 2 = face centered cubic */
double dia,		/* Diameter of simplex circumspehere */
double off,		/* Starting offset in device space */
double angle	/* Rotation angle 0.0 - 1.0 */
) {
	int  i, j, di = s->di;
	double sx[MXTD+1][MXTD];	/* Simplex verticies */

	switch (type) {

		case SIMDLAT_BCC:
			/* Create the node positions for body centered */
			/* cubic latice simplex basis vectors */
			/* Body centered places points at locations where the */
			/* lattice integer coordinates are all even or all odd. */
			for (i = 0; i < (di+1); i++) {
			
				for (j = 0; j < di; j++)
					sx[i][j] = -0.5;
				if (i < di) {
					if (i > 0)
						sx[i][i-1] += dia;
				} else {
					for (j = 0; j < di; j++)
						sx[i][j] += 0.5 * dia;
				}
			}
			break;

		case SIMDLAT_EQSPLX:
			/* Create the node positions for the */
			/* equalateral simplex basis vectors */
			for (i = 0; i < (di+1); i++) {
				double rr = 1.0;		/* Current radius squared */
				double ss = dia / sqrt(3.0);	/* Scale */
			
				/* The bounding points form a equalateral simplex */
				/* whose vertexes are on a sphere about the data */
				for (j = 0; j < di; j++) {
					double ddi;
					double hh = 1.0/(di-j);	/* Weight for remaining points */
			
					if (j > i)
						sx[i][j] = 0.0;				/* If beyond last */
					else if (j == i)				/* If last non-zero */
						sx[i][j] = ss * sqrt(rr);
					else							/* If before last */
						sx[i][j] = -hh * ss * sqrt(rr);
			
					ddi = (double)(di - j);
					rr *= (ddi * ddi - 1.0)/(ddi * ddi);
				}
			}
			break;

		case SIMDLAT_FCC:
			/* Create the node positions for the */
			/* face centered arrangement */
			/* Face centered places points at locations where the */
			/* sum of the lattice integer coordinates is even. */
			for (i = 0; i < (di+1); i++) {
			
				for (j = 0; j < di; j++)
					sx[i][j] = 0.0 * -0.5;

				if (i > 0 && i < di) {
					sx[i][i-1] += 0.5 * dia;
					sx[i][di-1] += 0.5 * dia;
				} else if (i == di) {
					sx[i][di-1] += dia;
				}
			}
			break;
	}

	/* Apply a rotation to avoid possible alignment with */
	/* the device axes */
	{
		int m, k;
		int ldi = di-1;			/* Last dimension */
		double a, b;

		b = angle;
		a = sqrt(1.0 - b * b);
		
		/* Apply rotation to all except last dimension */
		for (m = 0; m < ldi; m++) {			/* Dimension being rotated */

			for (i = 0; i < (di+1); i++) {	/* Node being rotated */
				double out[MXTD];

				for (j = 0; j < di; j++) {	/* Coord being produced */
					out[j] = 0.0;

					for (k = 0; k < di; k++) {	/* Coord being used */
						if ((j == m   && k == m)
						 || (j == ldi && k == ldi))
							out[j] += a * sx[i][k];	/* Diagonal multiplier */
						else if (j == m && k == ldi)
							out[j] += b * sx[i][k];
						else if (j == ldi && k == m)
							out[j] -= b * sx[i][k];
						else if (j == k)
							out[j] += sx[i][k];
					}
				}
				for (j = 0; j < di; j++)
					sx[i][j] = out[j];			/* Transfer result */
			}
		}
	}

#ifdef DEBUG	/* Dump stats on verticies */
for(i = 0; i < (di+1); i++) {
	double val = 0.0;
	printf("vert %d = ",i);
	for(j = 0; j < di; j++) {
		val += sx[i][j] * sx[i][j];
		printf("%f ",sx[i][j]);
	}
	printf(" (%f)\n",sqrt(val));
}

for(i = 0; i < di; i++) {
	for (j = i+1; j < (di+1); j++) {
		int e;
		double val;

		/* Distance between nodes */
		for (val = 0.0, e = 0; e < di; e++) {
			double tt = sx[i][e] - sx[j][e];
			val += tt * tt;
		}
		val = sqrt(val);
		printf("dist %d %d = %f\n",i,j,val);
		}
}
#endif	/* DEBUG */

	/* Convert from di+1 verticies to di base vectors */
	for (i = 0; i < di; i++) {
		for (j = 0; j < di; j++) {
			s->bv[i][j] = sx[i+1][j] - sx[i][j];
		}
	}

	/* Establish the basis origin */
	{
		for (j = 0; j < di; j++)
			s->bo[j] = off * s->ilimit/di;
	}
}

/* Compute the hash */
static int comp_hash(
simdlat *s,
int    *x		/* Index */
) {
	int j, di = s->di;
	unsigned long hash;

	for (hash = 0, j = 0; j < di; j++)
		hash = hash * 7 + x[j];
	hash %= SPT_HASHSIZE;

	return hash;
}

/* Check if a node already exists. Return -1 if not, */
/* or node index if it does. */
static int check_exists(
simdlat *s,
int     *x,		/* Index */
int     hash	/* Hash */
) {
	int di = s->di;
	int hp;		/* node index */
	int j;

	for (hp = s->hash[hash]; hp >= 0; hp = s->nodes[hp].hp) {

		/* Check if we have a match */
		for (j = 0; j < di; j++) {
			if (s->nodes[hp].x[j] != x[j])
				break;
		}
		if (j >= di)
			break;						/* Found a match */
	}

	return hp;
}

/* Create a new node. We assume it doesn't already exist */
/* Return its index */
static int new_node(
simdlat *s,
int    *x,		/* Index */
int    hash		/* Hash */
) {
	int di = s->di;
	int b = 0; 		/* NZ if a boundary point */
	int nn;			/* New node index */
	int hp;			/* Hash chain index */
	int i, j;

	/* Make room for it */
	if ((s->np+1) >= s->np_a) {
		s->np_a *= 2;
		if ((s->nodes = (sdtnode *)realloc(s->nodes, s->np_a * sizeof(sdtnode))) == NULL)
			error ("simdlat: node realloc failed");
	}

	nn = s->np++;				/* Add the new point */

	/* Compute the intended device value */
	for (j = 0; j < di; j++)
		s->nodes[nn].p[j] = s->bo[j];

	for (i = 0; i < di; i++) {
		for (j = 0; j < di; j++) {
			s->nodes[nn].p[j] += x[i] * s->bv[i][j];	/* Sum basis vector product */
		}
	}

	/* See whether we are well outside the gamut or not */
	b = snap_to_gamut2(s, s->nodes[nn].p);

	s->percept(s->od, s->nodes[nn].v, s->nodes[nn].p);	/* Compute perceptual */

	/* Store node information */
	for (j = 0; j < di; j++)
		s->nodes[nn].x[j] = x[j];

	s->nodes[nn].b = b;
	if (b < 2) {
		s->nodes[nn].vald = 1;			/* Valid if within or on gamut */
		s->nvp++;						/* Got another valid one */
	} else
		s->nodes[nn].vald = 0;			/* Not valid if it's a boundary point */
	s->nodes[nn].expm[0] =
	s->nodes[nn].expm[1] = (1 << di)-1;	/* Assum all dimensions need exploring */
	s->nodes[nn].hp = s->nodes[nn].up = -1;			/* Linked list indexes */

	/* Add an entry in the hash table */
	if (s->hash[hash] < 0)
		s->hash[hash] = nn;			/* We are the only entry */
	else {
		hp = s->hash[hash];
		while (s->nodes[hp].hp >= 0)
			hp = s->nodes[hp].hp;		/* Follow chain */
		s->nodes[hp].hp = nn;			/* Add at the end of the chain */
	}

	return nn;
}

/* ============================================= */
/* Main object functions */

/* Initialise, ready to read out all the points */
static void simdlat_reset(simdlat *s) {
	s->rix = 0;
}

/* Read the next set of non-fixed points values */
/* return non-zero when no more points */
static int simdlat_read(
simdlat *s,
double *d,		/* Device position */
double *p		/* Perceptual value */
) {
	int j;

	for (; s->rix < s->bnp; s->rix++) {

		if (s->bnodes[s->rix].vald != 0) {
			for (j = 0; j < s->di; j++) {
				if (d != NULL)
					d[j] = s->bnodes[s->rix].p[j];
				if (p != NULL)
					p[j] = s->bnodes[s->rix].v[j];
			}
			s->rix++;
			return 0;
		}
	}
	return 1;
}

/* Do a pass of seed filling the whole gamut, given a simplex dia. */
/* Return the number of nodes produced */
static int do_pass(
simdlat *s,
double dia		/* Simplex diameter to try */
) {
	int di = s->di;
	int hash;
	int i, j, k;
	int x[MXTD];
	int nn;			/* New nodes index */
	int np;

	/* Rest the current list */
	s->np = 0;
	s->nvp = 0;
	for (i = 0; i < SPT_HASHSIZE; i++)
		s->hash[i] = -1;

	/* Initial alloc of nodes */
	if (s->nodes == NULL) {
		s->np_a = 10;
		if ((s->nodes = (sdtnode *)malloc(s->np_a * sizeof(sdtnode))) == NULL)
		error ("simdlat: nodes malloc failed");
	}

	/* Compute the simplex basis vectors */
	/* arguments: simplex diameter, device space offset, angle to skew grid */
//	comp_basis(s, s->type, dia, 0.5, s->angle);
	comp_basis(s, s->type, dia, 0.4 + dia/150.0, s->angle);

	s->dia = dia;

	/* Add an initial seed point */
	for (j = 0; j < di; j++)
		x[j] = 0;
	hash = comp_hash(s, x);
	nn = new_node(s, x, hash);

	if (s->nodes[nn].b > 1) {
		error("simdlat: initial seed point is not within gamut");
	}
	
	s->unex = nn;		/* Initial entry in unexplored list */

//printf("~1 seed node is [%d %d]\n",s->nodes[nn].x[0], s->nodes[nn].x[1]);

	/* While there is more unexplored area */
	/* and we arn't finding a ridiculous number of points */
	while(s->unex >= 0 && (s->nvp < 3 * s->inp)) {
		int pos;				/* Positive or -ve direction */
		nn = s->unex;			/* Node we're looking at */
		s->unex = s->nodes[nn].up;		/* remove from unexplored list */

//printf("\n~1 exploring beyond node [%d %d]\n",s->nodes[nn].x[0], s->nodes[nn].x[1]);

		if (s->nodes[nn].b > 1)
			continue;		/* Don't look at boundary points */

		/* For all unexplored directions */
		for (i = 0; i < di; i++) {
			for (pos = 0; pos < 2; pos++) {
				int on;			/* Other node index */

//printf("~1 checking direction dim %d, sign %d, [%d %d]\n",i,pos,x[0],x[1]);

				if (((1 << i) & s->nodes[nn].expm[pos]) == 0) {
//printf("~1 that direction has been explored\n");
					continue;		/* Try next direction */
				}

				/* Check out that direction */
				for (j = 0; j < di; j++)
					x[j] = s->nodes[nn].x[j];
				x[i] += pos ? 1 : -1;

				/* If that node already exists */
				hash = comp_hash(s, x);
				if ((on = check_exists(s, x, hash)) >= 0) {
					/* back direction doesn't need checking */
					s->nodes[on].expm[pos ^ 1] &= ~(1 << i);
//printf("~1 that node already exists\n");
					continue;		/* Try next direction */
				}

				/* Create a new node in that direction */
				on = new_node(s, x, hash);
				
				if (s->nodes[on].b > 1) {	/* If new node is boundary, don't explore beyond it */
//printf("~1 added new boundary node [%d %d]\n",x[0],x[1]);
					continue;
				}
				/* back direction on new node doesn't need checking */
				s->nodes[on].expm[pos ^ 1] &= ~(1 << i);

//printf("~1 added new internal node [%d %d] **\n",x[0],x[1]);
				s->nodes[on].up = s->unex;		/* Add this node to unexplored list */
				s->unex = on;
			} 
		}
	}


	/* Rationalise cooincident points, and count final valid */
	s->nvp =  0;
	for (i = 0; i < s->np; i++) {

//printf("~1 rationalising %d, = %f %f\n",i, s->nodes[i].p[0], s->nodes[i].p[1]);
		if (s->nodes[i].vald == 0) {
//printf("~1 point %d is not valid\n",i);
			continue;
		}

		/* First against fixed points in device space */
		for (k = 0; k < s->fxno; k++) {
			double dd;

			/* Compute distance */
			dd = 0.0;
			for (j = 0; j < di; j++) {
				double tt = s->nodes[i].p[j] - s->fxlist[k].p[j];
				dd += tt * tt;
			}
			dd = sqrt(dd);

			if (dd < s->tol) {
				s->nodes[i].vald = 0;	/* Ignore this point */
//printf("~1 point %d matches input point %d\n",i, k);
				break;
			}
		}

		if (s->nodes[i].vald == 0)
			continue;

		/* Then against all the other points */
		for (k = i+1; k < s->np; k++) {
			double dd;

			if (s->nodes[k].vald == 0)
				continue;

			/* Compute distance */
			dd = 0.0;
			for (j = 0; j < di; j++) {
				double tt = s->nodes[i].p[j] - s->nodes[k].p[j];
				dd += tt * tt;
			}
			dd = sqrt(dd);

			if (dd < s->tol) {
				s->nodes[i].vald = 0;	/* Ignore this point */
//printf("~1 point %d matches other point %d\n",i, k);
				break;
			}
		}

		if (s->nodes[i].vald != 0)
			s->nvp++;		/* Found a valid one */
	}

#ifdef DUMP_PLOT
	/* Dump plot */
	dump_image(s, PERC_PLOT);
#endif /* DUMP_PLOT */

//printf("~1 got %d valid out of %d total\n",s->nvp, s->np);
	np = s->nvp;

	/* If we have a new best */
	if (s->nvp <= s->inp && (s->inp - s->nvp) < (s->inp - s->bnvp)) {
		sdtnode *tnodes;
		int tnp_a;
		tnodes = s->bnodes;		/* Swap them */
		tnp_a = s->bnp_a;
		s->bnp    = s->np;
		s->bnvp   = s->nvp;
		s->bnodes = s->nodes;
		s->bnp_a  = s->np_a;
		s->bdia   = s->dia;
		s->nodes  = tnodes;
		s->np_a   = tnp_a;
		s->np = s->nvp = 0;		/* Zero current */
	}

	return np;
}

/* Destroy ourselves */
static void
simdlat_del(simdlat *s) {

	if (s->nodes != NULL)
		free(s->nodes);
	if (s->bnodes != NULL)
		free(s->bnodes);

	free (s);
}

/* Constructor */
simdlat *new_simdlat(
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
int inp,				/* Number of points to generate */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
int type,				/* type of geometry, 0 = body centered cubic, */
						/* 1 = equilateral simplex, 2 = face centered cubic */
double angle,			/* Angle to orient grid at (0.0 - 0.5 typical) */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
) {
	int i;
	double ctol;
	double hdia, ldia, dia;
	int    hnp,  lnp,  np;
	simdlat *s;

#ifdef DEBUG
	printf("new_simdlat called with di %d, inp %d\n",di,inp);
#endif

	if ((s = (simdlat *)calloc(sizeof(simdlat), 1)) == NULL)
		error ("simdlat: simdlat malloc failed");

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	s->reset = simdlat_reset;
	s->read  = simdlat_read;
	s->del   = simdlat_del;

	/* If no perceptual function given, use default */
	/* (Not that it's used here) */
	if (percept == NULL) {
		s->percept = default_simdlat_to_percept;
		s->od = s;
	} else {
		s->percept = percept;
		s->od = od;
	}

	s->ilimit = ilimit;

	s->inp = inp - fxno;	/* Intended number of points */
	s->angle = angle;		/* desired grid angle */

	s->tol = SNAP_TOL;

	ctol = 0.6/pow((double)s->inp, 1.0/di);
	if (ctol < s->tol) {
		s->tol = ctol;
	}
//printf("~1 tol = %f\n",s->tol);

	s->fxlist = fxlist;		/* remember fixed points */
	s->fxno   = fxno;

	/* Compute perceptual values in fixed list */
	for (i = 0; i < s->fxno; i++)
		s->percept(s->od, s->fxlist[i].v, s->fxlist[i].p);

	if (di > MXTD)
		error ("simdlat: Can't handle di %d",di);
	s->di = di;

	/* We need to do a binary search to establish the desired */
	/* latice spacing. */

	/* Do an initial stab */
#ifdef NEVER
	dia = 0.3;
#else
	{	/* For body centered cubic */
		double vol = (double)ilimit/(double)di;
		double cellvol = (2.0 * vol * di)/(double)inp;
		dia = pow(cellvol, 1.0/di);
		printf("~1 initial dia = %f\n",dia);
	}
#endif
	np = do_pass(s, dia);
	if (np == 0)
		error("simdlat: First pass gave 0 points!");

//printf("~1 first cut dia %f ang %f gave %d points, target = %d\n",dia, s->angle, np, s->inp);

	if (np < s->inp) {		/* Low count */
		ldia = dia;
		lnp = np;
		for(;;) {
			dia = pow(np/(1.5 * s->inp), 1.0/di) * dia;
//printf("~1 next try dia %f in hope of %f\n",dia, 1.5 * s->inp);

			np = do_pass(s, dia);
//printf("~1 second cut dia %f ang %f gave %d points, target = %d\n",dia, s->angle, np,s->inp);
			if (np >= s->inp)
				break;
			ldia = dia;		/* New low count */
			lnp = np;
		} 
		hdia = dia;
		hnp = np;
	} else {
		hdia = dia;			/* High count */
		hnp = np;
		for(;;) {
			dia = pow(np/(0.6 * s->inp), 1.0/di) * dia;
//printf("~1 next try dia %f in hope of %f\n",dia, 0.6 * s->inp);
			np = do_pass(s, dia);
//printf("~1 second cut dia %f ang %f gave %d points, target = %d\n",dia, s->angle, np,s->inp);
			if (np <= s->inp)
				break;
			hdia = dia;			/* new high count */
			hnp = np;
		}
		ldia = dia;
		lnp = np;
	}

	/* Now zoom into correct number, with linear interp. binary search. */
	for (i = 0; s->bnvp != s->inp && i < MAX_TRIES; i++) {
		double ratio;
	
		/* Bail out early if we're close enough */
		if ((3 * i) > MAX_TRIES) {
			if (((double)s->bnvp/(double)s->inp) > 0.99)
				break;
		}

		ratio = ((double)s->inp - lnp)/(hnp - lnp);	/* Distance between low and high */
		dia = ratio * (hdia - ldia) + ldia;
		np = do_pass(s, dia);

//printf("~1 try %d, cut dia %f ang %f gave %d points, target = %d\n",i, dia, s->angle, np,s->inp);
		if (np > s->inp) {
			hdia = dia;
			hnp = np;
		} else {
			ldia = dia;
			lnp = np;
		}
	}
	
	simdlat_reset(s);

//printf("~1 total of %d patches\n",s->bnvp);

	return s;
}

/* =================================================== */

#ifdef STANDALONE_TEST

//#define ANGLE 0.33333
#define ANGLE 0.0

icxColorantLu *clu;

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
	int e, di = 2;

#ifndef NEVER
	/* Default Do nothing - copy device to perceptual. */
	for (e = 0; e < di; e++) {
		double tt = d[e];
		if (e == 0)
			tt = pow(tt, 2.0);
		else
			tt = pow(tt, 0.5);
		p[e] = tt * 100.0;
	}
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
	int npoints = 50;
	simdlat *s;
	int mask = ICX_BLACK | ICX_GREEN;
	
	error_program = argv[0];

	if (argc > 1)
		npoints = atoi(argv[1]);

	if ((clu = new_icxColorantLu(mask)) == NULL)
		error ("Creation of xcolorant lu object failed");

	/* Create the required points */
	s = new_simdlat(2, 1.5, npoints, NULL, 0, SIMDLAT_BCC, ANGLE, sa_percept, NULL);

#ifdef DUMP_PLOT
	printf("Perceptual plot:\n");
	dump_image_final(s, 1);

	printf("Device plot:\n");
	dump_image_final(s, 0);
#endif /* DUMP_PLOT */

	s->del(s);

	return 0;
}

#endif /* STANDALONE_TEST */



#if defined(DEBUG) || defined(DUMP_PLOT)

/* Dump the current point positions to a plot window file */
static void
dump_image(simdlat *s, int pcp) {
	double minx, miny, maxx, maxy;
	double *x1a = NULL;
	double *y1a = NULL;
	double *x2a = NULL;
	double *y2a = NULL;
	double *x3a = NULL;
	double *y3a = NULL;

	int i, nu;
	sdtnode *p;

	if (s->nvp == 0)
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
	
	if ((x1a = (double *)malloc(s->nvp * sizeof(double))) == NULL)
		error ("simdlat: plot malloc failed %d",s->nvp);
	if ((y1a = (double *)malloc(s->nvp * sizeof(double))) == NULL)
		error ("simdlat: plot malloc failed %d",s->nvp);
	if ((x2a = (double *)malloc(s->nvp * sizeof(double))) == NULL)
		error ("simdlat: plot malloc failed %d",s->nvp);
	if ((y2a = (double *)malloc(s->nvp * sizeof(double))) == NULL)
		error ("simdlat: plot malloc failed %d",s->nvp);

	for (nu = i = 0; i < s->np; i++) {
		p = &s->nodes[i];

		if (p->vald == 0)
			continue;
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

/* Dump the final point positions to a plot window file */
static void
dump_image_final(simdlat *s, int pcp) {
	double minx, miny, maxx, maxy;
	double *x1a = NULL;
	double *y1a = NULL;
	double *x2a = NULL;
	double *y2a = NULL;
	double *x3a = NULL;
	double *y3a = NULL;

	int i, nu;
	sdtnode *p;

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
	
	if ((x1a = (double *)malloc(s->bnvp * sizeof(double))) == NULL)
		error ("simdlat: plot malloc failed");
	if ((y1a = (double *)malloc(s->bnvp * sizeof(double))) == NULL)
		error ("simdlat: plot malloc failed");
	if ((x2a = (double *)malloc(s->bnvp * sizeof(double))) == NULL)
		error ("simdlat: plot malloc failed");
	if ((y2a = (double *)malloc(s->bnvp * sizeof(double))) == NULL)
		error ("simdlat: plot malloc failed");

	for (nu = i = 0; i < s->bnp; i++) {
		p = &s->bnodes[i];

		if (p->vald == 0)
			continue;
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





