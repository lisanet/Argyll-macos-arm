
/* 
 * Alphabetic index class.
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

/* TTBD:
	It would be nice if it was possible to add a delimeter between
    the X and Y identifiers so that number + number or alpha + alpha
    identification was possible.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "numsup.h"
#include "alphix.h"

/* Return the alpha index for the given raw index number */
/* Return NULL on error */
static char *torawix(alphix *p, int ix) {
	int i, j, n;
	int z = 0;
	char *rv, *v;

	if (ix < 0 || ix >= p->rmct)
		return NULL;		/* Index range error */

	if ((rv = malloc((p->nd+1) * sizeof(char))) == NULL)
		return NULL;

	for (v = rv, j = 0, i = p->nd-1; i >= 0; j++, i--) {
		char c;
		n = ix / p->ds[i].b;
		ix -= n * p->ds[i].b;
		c = p->ds[i].seq[n];
		if (p->ds[i].z && z == 0 && c == '0')
			c = ' ';
		if (z != 0 || c != ' ') {
			*v++ = p->ds[i].seq[n];
			z = 1;
		}
	}
	*v = '\000';

	return rv;
}

/* Return the the raw index number given the alpha index */
/* return -1 on error */
static int fromanat(alphix *p, char *ax) {
	char *v, *tb, _tb[11];
	int cl;
	int i, k, rv = 0;

	cl = strlen(ax);
	if (cl > p->nd)		/* String is too long to be our index */
		return -1;

	if (p->nd > 10) {
		if ((tb = malloc((p->nd+1) * sizeof(char))) == NULL)
			return -1;	/* Malloc error */
	} else
		tb = _tb;

	/* Pack the string out to the right number of digits with spaces. */
	for (v = tb; cl < p->nd; v++, cl++)
		*v = ' ';
	strcpy(v, ax);

	/* For each digit, convert to numerical and accumulate */
	for (v = tb, i = p->nd-1; i >= 0; v++, i--) {
		for (k = 0; k < p->ds[i].n; k++) {
			if (*v == p->ds[i].seq[k]
			 || (p->ds[i].z && *v == ' ' && p->ds[i].seq[k] == '0')) {
				rv += k * p->ds[i].b;
				break;
			}
		}
		if (k >= p->ds[i].n) {
			if (tb != _tb)
				free(tb);
			return -1;		/* Unknown digit */
		}
	}

	if (tb != _tb)
		free(tb);

	return rv;
}

/* Working from the end of the string, */
/* return a pointer to the maximum number of characters */
/* that are legaly part of this alpha index. */
static char *find_start(alphix *p, char *ax) {
	int i, k;
	char *v;

	v = ax + strlen(ax) - 1;	/* start at last character */
	for (i = 0; v >= ax && i < p->nd; v--, i++) {		/* working backwards */
		for (k = 0; k < p->ds[i].n; k++) {	/* Locate */
			if (*v == p->ds[i].seq[k])
				break;				/* Found */
		}
		if (k >= p->ds[i].n)		/* Not found, so we are at the start */
			break;
	}
	return v+1;
}

/* Convert from cooked index to raw index. */
/* Return -1 if out of range */
static int cookedtoraw(alphix *p, int ix) {
	int i;
	
	if (p->nr == 0)
		return ix;

	for (i = 0; i < p->nr; i++) {
		if (ix >= p->rs[i].c0 && ix <= p->rs[i].c1) {
			ix = ix - p->rs[i].c0 + p->rs[i].r0;
			return ix;
		}
	}
	return -1;
}

/* Convert from raw index to cooked index. */
/* Return -1 if out of range */
static int rawtocooked(alphix *p, int ix) {
	int i;
	
	if (p->nr == 0)
		return ix;

	for (i = 0; i < p->nr; i++) {
		if (ix >= p->rs[i].r0 && ix <= p->rs[i].r1) {
			ix = ix - p->rs[i].r0 + p->rs[i].c0;
			return ix;
		}
	}
	return -1;
}

/* Return the maximum possible length of count */
static int alphix_maxlen(alphix *p) {
	return p->cmct;
}

/* Return the alpha index for the given index number (0 based) */
/* Return NULL on error */
/* (Free after use) */
static char *alphix_aix(alphix *p, int ix) {
	if ((ix = cookedtoraw(p, ix)) < 0)
		return NULL;
	return torawix(p, ix);
}

/* Return the index number for the given alpha index */
int alphix_nix(alphix *p, char *ax) {
	int rv;
	if ((rv = fromanat(p, ax)) < 0)
		return -1;
	return rawtocooked(p, rv);
}

/* Destroy ourselves */
static void
alphix_del(alphix *p) {
	int i;
	for (i = 0; i < p->nd; i++)
		free(p->ds[i].seq);
	free (p->ds);
	free (p->rs);
	free (p);
}

/* Constructor: */
alphix *new_alphix(char *pattern) {
	alphix *p;
	char *pp = pattern;
	int i;

	if ((p = (alphix *)calloc(1,sizeof(alphix))) == NULL)
		error("alphix: malloc failed");

	p->maxlen = alphix_maxlen;
	p->aix    = alphix_aix;
	p->nix    = alphix_nix;
	p->del    = alphix_del;

	/* We now need to parse the pattern to setup our index sequence: */
	/* For all digits */
	for (p->nd = 0; ; p->nd++) {

		if (*pp == '\000' || *pp == ';') {
			break;
		}

		/* We are starting a new digit sequence */
		if (p->nd >= p->_nd) {		/* Allocate space for current digit */
			p->_nd = 2 + p->_nd;
			if ((p->ds = (dig *)realloc(p->ds, p->_nd * sizeof(dig))) == NULL)
				error ("alphix: realloc failed");
		}

		p->ds[p->nd].n = 0;
		p->ds[p->nd].seq = NULL;
		p->ds[p->nd]._n = 0;
		p->ds[p->nd].z = 0;

		/* For all symbols in this digit */
		for (p->ds[p->nd].n = 0; ;) {
			char c, c1, c2;

			/* We are adding another character to the sequence */
			if (*pp == '\000' || *pp == ';' || *pp == ',')
				break;						/* Done this digit */
			if (pp[1] == '-' && pp[2] != '\000' && pp[2] != ';' && pp[2] != ',') {
				c1 = pp[0];
				c2 = pp[2];
				pp += 3;
			} else {
				c1 = c2 = *pp;
				pp++;
			}
			if (c1 == '@') {
				c1 = '0';
				p->ds[p->nd].z = 1;
			}
			if (c2 == '@') {
				c2 = '0';
				p->ds[p->nd].z = 1;
			}

			/* Expand digit sequence */
			for (c = c1; c <= c2; c++) {
	
				if (p->ds[p->nd].n >= p->ds[p->nd]._n) {	/* Allocate space for next character */
					p->ds[p->nd]._n = 20 + p->ds[p->nd]._n;
					if ((p->ds[p->nd].seq = (char *)realloc(p->ds[p->nd].seq, p->ds[p->nd]._n * sizeof(char))) == NULL)
						error ("alphix: realloc failed");
				}
				p->ds[p->nd].seq[p->ds[p->nd].n++] = c;
			}
		}
		if (*pp == '\000' || *pp == ';')
			continue;
		pp++;

	}

	/* Compute the native maximum index count */
	for (p->rmct = 1, i = 0; i < p->nd; i++) {
		p->ds[i].b = p->rmct;
		p->rmct *= p->ds[i].n;
	}

	/* Search for valid ranges */
	if (*pp == ';') {
		char *v, *tb, _tb[11];

		pp++;
		if (p->nd > 10) {
			if ((tb = malloc((p->nd+1) * sizeof(char))) == NULL)
				error ("alphix: malloc failed");
		} else
			tb = _tb;

		/* For each range */
		for (p->nr = 0; ; p->nr++) {

			if (*pp == '\000' || *pp == ';')
				break;

			/* We are adding a new range */
			if (p->nr >= p->_nr) {		/* Allocate space for current range */
				p->_nr = 2 + p->_nr;
				if ((p->rs = (rngsq *)realloc(p->rs, p->_nr * sizeof(rngsq))) == NULL)
					error ("alphix: realloc failed");
			}

			/* Locate the end of the range start */
			for(v = tb; *pp != '\000' && *pp != '-' && *pp != ','; v++, pp++)
				*v = *pp;
			*v = '\000';
			p->rs[p->nr].r0 = p->rs[p->nr].r1 = fromanat(p, tb);
			if (p->rs[p->nr].r0 < 0)
				error("alphix: range start definition error on '%s'",tb);


			if (*pp != '-') {		/* oops - bad definition */
				error("alphix: range definition error - missing '-'");
			}

			/* Locate the end of the range end */
			for(v = tb, pp++; *pp != '\000' && *pp != ','; v++, pp++)
				*v = *pp;
			*v = '\000';
			p->rs[p->nr].r1 = fromanat(p, tb);
			if (p->rs[p->nr].r1 < 0)
				error("alphix: range end definition error on '%s'",tb);

			if (p->rs[p->nr].r1 < p->rs[p->nr].r0)	/* Hmm */
				error("alphix: range definition error, end < start ");

			/* Compute cooked index range */
			p->rs[p->nr].c0 = 0;
			p->rs[p->nr].c1 = p->rs[p->nr].r1 - p->rs[p->nr].r0;
			if (p->nr > 0) {
				int ofs = p->rs[p->nr-1].c1+1;
				p->rs[p->nr].c0 += ofs;
				p->rs[p->nr].c1 += ofs;
			}

			if (*pp == '\000' || *pp == ';')
				continue;
			pp++;
		}

		if (tb != _tb)
			free(tb);
	}

	/* Compute cooked actual range */
	p->cmct = p->rmct;
	if (p->nr > 0) {
		p->cmct = p->rs[p->nr-1].c1+1;
	}

#ifdef NEVER
	// ###########################################################
	/* Debug stuff */
	printf("~1 There are %d digits, %d ranges, raw max = %d, cooked max = %d\n",
		p->nd,p->nr,p->rmct, p->maxlen(p));

	/* Digit sequences */
	for (i = 0; i < p->nd; i++) {
		printf(" ~1 Digit %d has %d symbols\n",i,p->ds[i].n);
		for (j = 0; j < p->ds[i].n; j++)
			printf("  ~1 symbol %d = '%c'\n",j,p->ds[i].seq[j]);
	}

	/* Ranges */
	for (i = 0; i < p->nr; i++) {
		printf(" ~1 Range %d is raw %d - %d, cooked %d - %d\n",
		       i,p->rs[i].r0, p->rs[i].r1, p->rs[i].c0, p->rs[i].c1);
	}

	/* Itterate the sequence */
	for (i = 0; i < p->cmct; i++) {
		char *v;
		v = p->aix(p, i);
		printf("ix %d -> '%s'\n",i,v);
		free(v);
	}
	// ###########################################################
#endif /* NEVER */

	return p;
}


/* ==================================================================== */
/* Utility function: */
/* Given the strip and patch alphix objects, and order flag, */
/* Return a patch location. Free the returned string after use. */
/* Return NULL on error */
char *patch_location(
	alphix *saix,		/* Strip alpha index object */
	alphix *paix,		/* Patch alpha index object */
	int ixord,			/* Index order, 0 = strip then patch */
	int six,			/* Strip index */
	int pix				/* Patch index */
) {
	char *sl, *pl, *rv;
	int ll;

	if ((sl = saix->aix(saix, six)) == NULL)
		return NULL;

	if ((pl = paix->aix(paix, pix)) == NULL) {
		free (sl);
		return NULL;
	}
	ll = strlen(sl) + strlen(pl) + 1;
	if ((rv = malloc(ll * sizeof(char))) == NULL) {
		free (pl);
		free (sl);
		return NULL;
	}

	if (ixord == 0) {
		strcpy (rv, sl);
		strcat (rv, pl);
	} else {
		strcpy (rv, pl);
		strcat (rv, sl);
	}
	return rv;
}

/* ==================================================================== */
/* Utility function: */
/* Given the strip and patch alphix objects, and order flag, */
/* and a corresonding patch location string, return an index */
/* number suitable for sorting location strings. */
/* The sort order is in order of strips, then patches within each strip */
/* Return -1 on error */
int patch_location_order(
	alphix *saix,		/* Strip alpha index object */
	alphix *paix,		/* Patch alpha index object */
	int ixord,			/* Index order, 0 = strip then patch */
	char *_ax			/* Patch location string */
) {
	char *ax;			/* Copy of input string */
	char *v;
	alphix *rh;			/* Least significant, right hand alphix */
	alphix *lh;			/* Most significant, left hand alphix */
	int ri, li;			/* Right hand and left hand index numbers */

	if ((ax = malloc(strlen(_ax)+1)) == NULL)
		return -1;
	strcpy(ax,_ax);

	if (ixord == 0) {
		lh = saix;		/* Strip is left hand */
		rh = paix;		/* Patch is right hand */
	} else {
		rh = saix;		/* Strip is right hand */
		lh = paix;		/* Patch is left hand */
	}

	/* We need to identify the boundary between */
	/* the right hand and left hand indexes. */
	/* We assume that the sequences are distinguishable ... */
	v = find_start(rh, ax);

	if (*v == '\000')	/* Nothing looks like a alphix */
		return -1;

	ri = rh->nix(rh, v);
	*v = '\000';
	li = lh->nix(lh, ax);
	free(ax);
	if (ri < 0 || li < 0)
		return -1;

	if (ixord == 0)		/* Strip is left hand */
		return li * rh->cmct + ri;
	else
		return ri * lh->cmct + li;
}










