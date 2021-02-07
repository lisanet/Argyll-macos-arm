
#ifndef IFARP_H

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

/* A sample point node */
struct _ifpnode {
	int    fx;			/* nz if point is fixed (existing) */
	double p[MXTD];		/* Device coordinate position */
	double v[MXTD];		/* Subjective value (Labnnn..) */

	unsigned int touch;	/* nn: Per value touch count */
}; typedef struct _ifpnode ifpnode;


/* Main simplex latice object */
struct _ifarp {
/* private: */
	int di;			/* Point dimensionality */
	double ilimit;	/* Ink limit - limit on sum of p[] */
	int inp;		/* Intended number of points in list */
	int np;			/* Number of point nodes in list */
	ifpnode *nodes;	/* Current array of nodes */
	int rix;		/* Next read index */

	/* Perceptual function */
	void (*percept)(void *od, double *out, double *in);
	void *od;		/* Opaque data for perceptual point */
	
	/* nn support */
	ifpnode **sax[MXTD];	/* Sorted axis pointers, one for each direction */
	unsigned int tbase;		/* Touch base value for this pass */
	unsigned int ttarget;	/* Touch target value for this pass */

/* public: */
	/* Initialise, ready to read out all the points */
	void (*reset)(struct _ifarp *s);

	/* Read the next set of non-fixed points values */
	/* return non-zero when no more points */
	int (*read)(struct _ifarp *s, double *d, double *p);

	/* Destroy ourselves */
	void (*del)(struct _ifarp *s);

}; typedef struct _ifarp ifarp;

/* Constructor */
extern ifarp *new_ifarp(int verb, int di, double ilimit, int npoints,
	fxpos *fxlist, int fxno, 
	void (*percept)(void *od, double *out, double *in), void *od
);

#define IFARP_H
#endif /* IFARP_H */
