
/* 
 * Argyll Color Correction System
 *
 * Scanrd: Scan chart reader
 * This is the core chart recognition code.
 * Private H file for scanrd.c
 *
 * Author: Graeme W. Gill
 * Date:   28/9/96
 *
 * Copyright 1995, 1996, 2008, Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include "scanrd.h"		/* Include public structure */

#include "../h/llist.h" /* Linked list macros */

#ifndef M_PI
#define M_PI		3.1415926535897932384626433832795028841971693993751
#endif
#ifndef M_PI_2
#define M_PI_2      (0.5 * M_PI)
#endif
#ifndef M_PI_4
#define M_PI_4      (0.25 * M_PI)
#endif
#ifndef M_PI_3_4
#define M_PI_3_4    (0.75 * M_PI)
#endif
#ifndef M_2_PI
#define M_2_PI      (2.0/M_PI)
#endif
#ifndef M_SQRT_2
#define M_SQRT_2	1.4142135623730950488016887242096980785696718753769
#endif

#define DEG(aa) ((aa) * 180.0/M_PI)

#define iabs(a)  ((a) < 0 ? -(a) : (a))

#define MXDE 4		/* Maximum useful pixel depth */

/* A run of points at a given y, in a group */
typedef struct {
	int y,lx,hx;	/* Region covers values of lx <= x < hx */
	} run;

/* Structure to hold the points that may make up a line */
struct _points {
	int no;		/* number of coordinates (indexed from ca) */
	int mxno;	/* Maximum number that can be allocated from ca */
	run *r;		/* point runs array */
	int pn;		/* Diagnostic serial no */
	LINKSTRUCT(struct _points);	/* Linked list structure */

	/* Line stats */
	int flag;
	int nop;			/* Number of points */
	double mw,len;		/* mean width, length */
	double mx,my;		/* Mean point */
	double a;			/* Angle */
	double ca;			/* Constrained Angle (constrained to 90 degrees about 0 */
	double x1,y1,x2,y2;	/* Start/end point of fitted line */

	double pmx, pmy;			/* Raster (Perspective affected) mean point */
	double px1, py1, px2, py2;	/* Raster (Perspective affected) line start/end points */
	struct _points *opt;		/* Next in improvement list */
}; typedef struct _points points;

#define F_LINESTATS 	0x01		/* Line stats valid */
#define F_VALID			0x02		/* Line passes valid criteria */
#define F_LONGENOUGH	0x04		/* Line is long enough to be included in angle calc */
#define F_IMPROVE	    0x08		/* Line was used to improve fit */

/* Structure to hold an aggregation region description */
struct _region {
	int lx,hx;	/* Region covers values of lx <= x < hx */
	points *p;	/* Head of points linked list associated with region */
}; typedef struct _region region;

/* Structure to hold a 2D real point */
struct _point {
	double x,y;
}; typedef struct _point point;

/* Structure to hold one entry in an edge list */
struct _epoint {
	double pos;		/* Position of entry along edge */
	double len;		/* (Maximum normalized) Total length */
	double p1,p2;	/* Start and end of line in orthogonal direction */
	double ccount;	/* (Maximum normalized) Crossing count */
	struct _points *opt;	/* Linked list of feature lines to optimize to */
	int    nopt;
}; typedef struct _epoint epoint;

/* Structure of an edge list */
struct _elist {
	epoint *a;		/* Array of edge points */
	int c;			/* Count */
	double lennorm;	/* Total of max. normalized lengths of edge list */
}; typedef struct _elist elist;

#define INIT_ELIST(e) { \
	e.a = NULL;         \
	e.c = 0;            \
	e.lennorm = 0.0;    \
}

/* An edge correlation return structure */
typedef struct {
	double cc,off,scale;
} ematch;


/* An integer point */
struct _ipoint {
	int x,y;
}; typedef struct _ipoint ipoint;

/* An edge scan structure */
struct _escan {
	int e[4];			/* Indexes of points for left or right sides */
	int i;				/* Index of current pair */
	int ev,k1,k2;		/* Bresenham constants */
	int y;				/* Current y value */
	int xi,x;			/* X increment, current x value */
}; typedef struct _escan escan;

/* Structure of a sample box */
#define SBOX_NAME_SZ 20
struct _sbox {
	int diag;						/* Non-zero if diagnostic only */
	char name[SBOX_NAME_SZ];		/* Box name (usualy letter number coordinate) */
	double xpt[3];					/* Lab expected value. L < 0 if not valid */
	double x1,y1,x2,y2;				/* Reference box corner points */
	ipoint p[4];					/* Transformed sample box coordinates */
	int active;						/* Flag to indicate box is active in scan */
	int ymin,ymax;					/* Min and nmax y values */
	escan l,r;						/* left and right edge scan structures */

	unsigned long *ps[MXDE];		/* Pixel value histogram arrays (256 or 65536) */
	/* Pixel values just scanned, or from best rotation */
	double mP[MXDE];				/* Mean Pixel values (0.0 - 255.0) */
	double sdP[MXDE];				/* Standard deviations */
	double P[MXDE];					/* Output Pixel values */
	unsigned int cnt;				/* Total pixels in cell */
	LINKSTRUCT(struct _sbox);		/* Active edge linked list structure */

	/* Pixel values for alternate rotations, created if we have expected values */
	struct {
		double mP[MXDE];				/* Mean Pixel values (0.0 - 255.0) */
		double sdP[MXDE];				/* Standard deviations */
		double P[MXDE];					/* Robust mean Output Pixel values (0.0 .. 255.0) */
		unsigned int cnt;				/* Total pixels in cell */
	} rot[4];

	}; typedef struct _sbox sbox;

/* The private scanrd object */
struct _scanrd_ {
	/* Public part of structure */
	scanrd public;
	
	/* Private variables */
	int flags;				/* Operation/diagnostic flags */
	int verb;				/* verbosity level */
							/* 0 = none */
							/* 1 = warnings */
							/* 2 = minimum */
							/* 3 = per patch */
							/* 3 = patch scan details */
							/* 5 = per line */
							/* 7 = matching attempts */
							/* 8 = per pixel */

	unsigned int errv;		/* Error value */
	char errm[200];			/* Error message */

	double gammav;			/* Approximate gamma encoding of input image */
	int width,height;		/* With and height of raster in pixels */
	int depth;				/* Useful pixel plane depth */
	int tdepth;				/* Total pixel plane depth */
	int bpp;				/* Bit precision per pixel, 8 or 16 */
	int bypp;				/* Bytes per pixel, either 1 or 2 */

	unsigned char *out;		/* Diagnostic output raster array */
	
	int noslines;			/* Number of lines with valid stats */
	int novlines;			/* Number of valid lines */
	points *gdone;			/* Head of done point linked list groups */

	double ppc[4];			/* Partial perspective correction values. */
							/* persp() applies perspective distortion, */
							/* invpersp() applies perspective correction. */

	double irot;			/* Base image rotation value in radians (clockwize) */
	int norots;				/* Number of rotations to explore */
	int crot;				/* Current or best rotation being scanned */
	struct {
		double irot;		/* Image rotation value in radians (clockwize) for this rotation */
		double ixoff,iyoff,ixscale,iyscale;	/* Image offset and scale factors */
		double cc;			/* Edge match correlation */
		double xcc;			/* Expected value correlation, smaller is better */
	} rots[4];
	
	double ptrans[8];		/* Combined transform of partial perspective, */
							/* irot, i[xy]off and i[xy]scale. */
							/* Use ptrans(ptrans[]) to transform from reference to raster. */
	double iptrans[8];		/* Inverse transform of ptrans[] */
							/* Use ptrans(iptrans[]) to transform from raster to reference */
	
	elist xelist, yelist;	/* X and Y raster edge lists array */
	elist ixelist, iyelist;	/* Inverted direction X and Y raster edge lists array */
	
	elist rxelist, ryelist;	/* X and Y .cht reference edge lists array */

	double rbox_shrink;		/* Reference box shrink factor */
	int xpt;				/* NZ if got expected reference values */
	
	double fid[8];			/* Four fiducial locations, typicall clockwise from top left */
	double fidsize;			/* Fiducial diagnostic cross size */
	double havefids;		/* NZ if there are fiducials */
	int nsbox;				/* Number of sample boxes */
	sbox *sboxes;			/* List of sample boxes */
	sbox **sbstart;			/* Sorted start list */
	sbox **sbend;			/* Sorted end list */
	int csi,cei;			/* Current start/end indexes */
	sbox *alist;			/* Active list during pixel value sampling */
	
	double adivval;			/* Overall average divider value */
	int divc;				/* Average divide count */

	int next_read;			/* Next box value to read */

	char *refname;			/* Path of reference file */
	
	int inited;				/* Gamma and regions inited */
	unsigned short gamma[256 * 256];	/* Inverse gamma lookup */
	region *vrego, *vregn;	/* Old and New region for delX or vertical lines */
	int no_vo, no_vn;		/* Number of regions in array */
	region *hrego, *hregn;	/* Old and New region for delY or horizontal lines */
	int no_ho, no_hn;		/* Number of regions in array */
	double th;				/* Color change threshold */
	double divval;			/* Current orthogonal divider value */

	/* aa line stuff */
	int aa_inited;			/* Non-zero if anti-aliased line tables are inited */
	int *coverage;			/* Coverage lookup array */
	int covercells;
	int covershift;
	int Pmax;
	int adj_pixinc[4];		/* Pixel address increments for 4 directions */
	int diag_pixinc[4];
	int orth_pixinc[4];

	/*** Callbacks ***/

	int (*read_line)(void *fdata, int y, char *dst);
							/* Read line of source file, non-zero on error */
	void *fdata;			/* Opaque data for callback */

	int (*write_line)(void *ddata, int y, char *src);
							/* Write RGB line of diag file, non-zero on error */
	void *ddata;			/* Opaque data for callback */

}; typedef struct _scanrd_ scanrd_;

/*************************************************************************/
/* Heapsort macro */

/* HEAP_COMPARE(A,B)  returns true if A < B */
#define HEAPSORT(TYPE,ARRAY,NUMBER) \
		{	\
		TYPE *hs_ncb = ARRAY;	\
		int hs_l,hs_j,hs_ir,hs_i;	\
		TYPE hs_rra;	\
		\
		if (NUMBER >= 2)	\
			{	\
			hs_l = NUMBER >> 1;	\
			hs_ir = NUMBER-1;	\
			for (;;)	\
				{	\
				if (hs_l > 0)	\
					hs_rra = hs_ncb[--hs_l];	\
				else	\
					{	\
					hs_rra = hs_ncb[hs_ir];	\
					hs_ncb[hs_ir] = hs_ncb[0];	\
					if (--hs_ir == 0)	\
						{	\
						hs_ncb[0] = hs_rra;	\
						break;	\
						}	\
					}	\
				hs_i = hs_l;	\
				hs_j = hs_l+hs_l+1;	\
				while (hs_j <= hs_ir)	\
					{	\
					if (hs_j < hs_ir && HEAP_COMPARE(hs_ncb[hs_j],hs_ncb[hs_j+1]))	\
						hs_j++;	\
					if (HEAP_COMPARE(hs_rra,hs_ncb[hs_j]))	\
						{	\
						hs_ncb[hs_i] = hs_ncb[hs_j];	\
						hs_i = hs_j;	\
						hs_j = hs_j+hs_j+1;	\
						}	\
					else	\
						hs_j = hs_ir + 1;	\
					}	\
				hs_ncb[hs_i] = hs_rra;	\
				}	\
			}	\
		}

/*************************************************************************/

