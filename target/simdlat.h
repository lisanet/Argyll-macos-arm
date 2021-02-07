
#ifndef SIMDLAT_H
/* 
 * Argyll Color Correction System
 *
 * Simplex device space latice test point class
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
 * Based on simplat.h
 */

#define SPT_HASHSIZE 4463	/* Hash index size */

/* A sample point node */
struct _sdtnode {
	int    vald;		/* Valid flag */
	int    x[MXTD];		/* Index */
	double p[MXTD];		/* Device coordinate position */
	double v[MXTD];		/* Subjective value (Labk) */
	int    b;			/* 1 if a gamut point, 2 if a boundary point */
	unsigned long expm[2];	/* -ve/+ve dimention explored flags */

	int hp;				/* Hash linked list index */
	int up;				/* Unexplored linked list index */
}; typedef struct _sdtnode sdtnode;


/* Main simplex latice object */
struct _simdlat {
/* private: */
	int di;			/* Point dimensionality */
	double ilimit;	/* Ink limit - limit on sum of p[] */
	int inp;		/* Intended number of points in list */
	int type;		/* Type of point geometry, 0 = body centered etc. */
	double angle;	/* Grid angle */
	double bo[MXTD];		/* Basis origin */
	double bv[MXTD][MXTD];	/* Simplex basis vectors */
	int np;			/* Number of point nodes in list */
	int nvp;		/* Number of valid point nodes in list */
	int np_a;		/* Number of points allocated */
	double dia;		/* Point spacing in latice */
	sdtnode *nodes;	/* Current array of nodes */
	int bnp;		/* Number of point nodes in best list */
	int bnvp;		/* Number of best valid point nodes in list */
	int bnp_a;		/* Number of best points allocated */
	double bdia;	/* Point spacing in best latice */
	sdtnode *bnodes;	/* Current best array of nodes */
	int hash[SPT_HASHSIZE];	/* Hash index */
	int unex;		/* Head of unexplored list index */
	double tol;		/* Snap tollerance */

	/* Perceptual function */
	void (*percept)(void *od, double *out, double *in);
	void *od;		/* Opaque data for perceptual point */
	
	/* Fixed points to avoid */
	fxpos *fxlist;			/* List of existing fixed points (may be NULL) */
	int fxno;				/* Number of existing fixes points */

	int rix;				/* Next read index */

/* public: */
	/* Initialise, ready to read out all the points */
	void (*reset)(struct _simdlat *s);

	/* Read the next set of non-fixed points values */
	/* return non-zero when no more points */
	int (*read)(struct _simdlat *s, double *d, double *p);

	/* Destroy ourselves */
	void (*del)(struct _simdlat *s);

	}; typedef struct _simdlat simdlat;

	/* geometry type */
#define SIMDLAT_BCC 0			/* Body centered Cubic */
#define SIMDLAT_EQSPLX 1		/* Equilateral Simplex */
#define SIMDLAT_FCC 2			/* Face Centered Cubic */

/* Constructor */
extern simdlat *new_simdlat(int di, double ilimit, int npoints,
	fxpos *fxlist, int fxno, int type, double angle,
	void (*percept)(void *od, double *out, double *in), void *od);

#define SIMDLAT_H
#endif /* SIMDLAT_H */
