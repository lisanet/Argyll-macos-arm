
#ifndef PPOINT_H

/* 
 * Argyll Color Correction System
 *
 * Perceptually distributed point class
 *
 * Author: Graeme W. Gill
 * Date:   16/10/96
 *
 * Copyright 1996 - 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#define MXPD 4			/* Maximum ppoint dimentionality */
#define POW2MXPD 16		/* 2 ^ MXPD */
#define POW3MXPD 81		/* 3 ^ MXPD */
#define MXNP (MXPD + 1 + 20)	/* Maximum near points */

/* tuning parameters */
#define WPOINTS 20		/* Points returned per group */
#define FPOINTS 5000	/* Number of far points to track - more is better. */
#define OPOINTS  250	/* Number of optimsed far points to use - more is slower */
#define DDMIX 0.75		/* Device distance to perceptual ratio in distance computation */
#define DWEIGHT 0.05	/* Distance factor weight added to maximum error in opt func */
#define CLOSED 0.05		/* Too close criteria */

/* A sample point node */
struct _node {
	int    fx;			/* nz if point is fixed */
	double p[MXPD];		/* Device coordinate position */
	double v[MXPD];		/* Subjective value (Labk) */
}; typedef struct _node node;

/* Main perceptual point object */
struct _ppoint {
/* private: */
	int di;			/* Point dimensionality */
	double ilimit;	/* Ink limit - limit on sum of p[] */
	int fnp;		/* Number of existing fixed points in list */
	int tinp;		/* target number of total points in list */

	node *list;		/* tinp list of points */
	int np;			/* Number of points currently in list */

	/* Perceptual function handed in */
	void (*percept)(void *od, double *out, double *in);
	void *od;		/* Opaque data for perceptual point */
	
	/* Progressive interpolation grid */
	rspl *g;

	/* Perceptual distance map */
	rspl *pd;

	/* Candidate far point starting values */
	co fp[FPOINTS];		/* Candidate points */
	int nfp;			/* Current number in fp[] */
	int wfp;			/* Index of current worst far point */
	double wfpd;		/* worst far point distance */ 
	co fwfp[FPOINTS];	/* Working space for find_worst() */

	/* Other info */
	int rix;			/* Next read index */
//	double mn,mx,av;	/* Perceptual distance stats */
	
/* public: */
	/* return non-zero if the perceptual point is within the device gammut */
	int (*pig)(struct _ppoint *s, double *p);

	/* Initialise, ready to read out all the points */
	void (*reset)(struct _ppoint *s);

	/* Read the next set of non-fixed points values */
	/* return non-zero when no more points */
	int (*read)(struct _ppoint *s, double *p, double *f);

	/* Calculate and print stats */
	void (*stats)(struct _ppoint *s);

	/* Destroy ourselves */
	void (*del)(struct _ppoint *s);

	}; typedef struct _ppoint ppoint;


/* Constructor */
extern ppoint *new_ppoint(int di, double ilimit, int npoints,
	fxpos *fxlist, int fxno, 
	void (*percept)(void *od, double *out, double *in), void *od);

/* ------------------------------------------------------- */
/* Macros for a di dimensional counter */
/* Declare the counter name nn, dimensions di, & count */

#define DCOUNT(nn, di, start, reset, count) 				\
	int nn[MXPD];	/* counter value */						\
	int nn##_di = (di);		/* Number of dimensions */		\
	int nn##_stt = (start);	/* start count value */			\
	int nn##_rst = (reset);	/* reset on carry value */		\
	int nn##_res = (count);	/* last count +1 */				\
	int nn##_e				/* dimension index */

/* Set the counter value to 0 */
#define DC_INIT(nn) 								\
{													\
	for (nn##_e = 0; nn##_e < nn##_di; nn##_e++)	\
		nn[nn##_e] = nn##_stt;						\
	nn##_e = 0;										\
}

/* Increment the counter value */
#define DC_INC(nn)									\
{													\
	for (nn##_e = 0; nn##_e < nn##_di; nn##_e++) {	\
		nn[nn##_e]++;								\
		if (nn[nn##_e] < nn##_res)					\
			break;	/* No carry */					\
		nn[nn##_e] = nn##_rst;						\
	}												\
}

/* After increment, expression is TRUE if counter is done */
#define DC_DONE(nn)									\
	(nn##_e >= nn##_di)
	
/* ------------------------------------------------------- */

#define PPOINT_H
#endif /* PPOINT_H */
