
#ifndef ALPHIX_H

/* 
 * Alphabetic indexing class
 */

/*
 * Argyll Color Correction System
 *
 * Author: Graeme W. Gill
 * Date:   22/8/2005
 *
 * Copyright 2005, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* 
 *
 *	Syntax of alphix pattern:
 *
 * First comes the definition of the symbols for each digit
 * location, LS to MS. The number of definitions declares the
 * maximum number of digits. For example, for a normal 2 digit numerical
 * sequence: "0123456789, 123456789" (note the space is significant)
 * would define 0..99 with the MS digit supressed when it is 0.
 * Ranges can be used for brevity: * "0-9, 1-9".
 * As a special case, the '@' character can be used 
 * instead of '0' to indicate suppression of the leading zero.
 * Leading ' ' characters in a generated sequence are omitted.
 *
 * Optional, and delimited by a ';' character, valid segments of the
 * index sequence can be defined. For instance, to define the index
 * range to be 1..49 one could use the pattern "0-9, 1-9;1-49"
 *
 * Of course the main reason for using alphix is to allow letter index
 * sequences. For a sequence A, B, C .. AA, AB, AC etc. (the default
 * used in Argyll), the following pattern would be used: "A-Z, A-Z"
 *
 * For a some ECI2002R charts that skip columns Y and Z, the following
 * might be used: "A-Z, 2-9;A-X,2A-9Z"
 */

#define ALPHIX_MAX 4		/* Maximum digits allowed */

/* Definition for each digit sequence */
typedef struct {
	int n;		/* Number of characters in sequence */
	char *seq;	/* Sequence of characters */
	int _n;		/* Allocation size of seq */
	int b;		/* Base of this digit */
	int z;		/* NZ if leading zero is to be supressed */
} dig;

/* Definition for each range sequence */
typedef struct {
	int r0,r1;		/* Raw index start and end of range */
	int c0,c1;		/* Cooked index start and end of range */
} rngsq;

struct _alphix {
/* private: */
	int nd;		/* Number of digits */
	dig *ds;	/* Digit sequences */
	int _nd;	/* Allocation size of ds */
	int rmct;	/* Raw maximum count */
	int cmct;	/* Cooked maximum count */

	int nr;		/* Number of ranges */
	rngsq *rs;	/* Digit sequences */
	int _nr;	/* Allocation size of rs */

/* public: */
	/* Return the maximum possible length of count */
	int (*maxlen)(struct _alphix *p);

	/* Return the alpha index for the given index number (0 based) */
	/* Return NULL on error */
	char *(*aix)(struct _alphix *p, int ix);

	/* Return the index number for the given alpha index */
	/* Return -1 on error */
	int (*nix)(struct _alphix *p, char *ax);

	/* Destroy ourselves */
	void (*del)(struct _alphix *p);

}; typedef struct _alphix alphix;

/* Constructor: */
extern alphix *new_alphix(char *pattern);


/* Utility function: */
/* Given the strip and patch alhix objects, and order flag, */
/* Return a patch location. */
/* Return NULL on error */
char *patch_location(
	alphix *saix,		/* Strip alpha index object */
	alphix *paix,		/* Patch alpha index object */
	int ixord,			/* Index order, 0 = strip then patch */
	int six,			/* Strip index (from 0) */
	int pix);			/* Patch index (from 0) */

/* Utility function: */
/* Given the strip and patch alphix objects, and order flag, */
/* and a coresonding patch location string, return an index */
/* number suitable for sorting location strings. */
/* Return -1 on error */
int patch_location_order(
	alphix *saix,		/* Strip alpha index object */
	alphix *paix,		/* Patch alpha index object */
	int ixord,			/* Index order, 0 = strip then patch */
	char *_ax);			/* Patch location string */

#define ALPHIX_H
#endif /* ALPHIX_H */
