#ifndef NAMEDC_H
#define NAMEDC_H

/* 
 * Argyll Color Correction System
 * Named color set support.
 *
 * Author: Graeme W. Gill
 * Date:   3/12/2013
 *
 * Copyright 2013 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License2.txt file for licencing details.
 *
 */

/*
 * This class provides storage for a list of named colors
 * in a file format independent fashiom.
 *
 * Currently supports CxF2, CxF3 & ICC named color profiles.
 *
 * Currently these are all assumed to be reflective colors.
 *
 * Current implementation is read only, but it could be converted
 * to support creating and writing too.
 *
 */

/* ------------------------------------------------------------------------------ */

/* An individual named color */

struct _nce {
	char *name;				/* Name in ASCII/UTF-8 */

	double Lab[3];			/* D50 L*a*b* */
	int Lab_v;				/* nz if Lab is valid */
	
	xspect *sp;				/* Non-NULL if valid spectral reflectance, norm 100% */

	icColorSpaceSignature devSig;	/* Device colorspace signature, icMaxEnumData if invalid */
	int dev_n;				/* Number of channels */
	double dev[MAX_CHAN];	/* Optional device values, ie CMYK in % */

}; typedef struct _nce nce;


struct _namedc {

	char *creator;		/* Creator name in ASCII/UTF-8 */
	char *description;	/* Description in ASCII/UTF-8 */

	int hash;			/* hash of filename and description */

	icxIllumeType ill;		/* Illuminant values was measured under */
	icxObserverType obs;	/* Observer model */

    unsigned int count;		/* Count of named colors */
    unsigned int count_a;	/* Allocated number of colors */
    nce  *data;				/* Array of [count] color values */

  /* Public: */
	void (*del)(struct _namedc *p);

#define NAMEDC_OP_NONE   0x0000
#define NAMEDC_OP_NODATA 0x0001			/* Don't load any data, just description */
#define NAMEDC_OP_NOSPEC 0x0002			/* Don't load spectral data */

	/* Read a cxf 3 format named color file */
	/* return nz on error */
	int (*read_cxf)(struct _namedc *p, const char *filename, int options);

	/* Read an ICC format named color file */
	/* return nz on error */
	int (*read_icc)(struct _namedc *p, const char *filename, int options);

	/* Read any format named color files */
	int (*read)(struct _namedc *p, const char *filename, int options);

	/* Return the index of the best mataching color, -1 on error. */
	/* Lab[] is assumed to be D50, 2 degree standard observer based CIE value, */
	/* and the spec value should only be provided if this is a reflective or */
	/* transmissive measurement, NULL if emissive. */
	/* If named color library is expects other than D50, 2 degree, then */
	/* it will use the spectral value if not NULL, or chromatically */
	/* adapt the Lab value. */
	/* deType == 0 DE76 */
	/* deType == 1 DE94 */
	/* deType == 2 DE2000 */
	/* if de != NULL, return the delta E */
	int (*match)(struct _namedc *p, double *de, double *Lab, xspect *spect, int deType);

	/* Houskeeping - should switch this to a1log ? */
#define NAMEDC_ERRL 1000
	int errc;				/* Error code */
	char err[NAMEDC_ERRL];	/* Error message */

  /* Private: */
	a1log *log;
	char *filename;			/* So we can lazy read the colors */
	int format;				/* 0 = unknown, 1 = cxf, 2 = ICC */
	int options;
	int indata;				/* State flag for sax_cb() */

	char pfx[100];			/* Prefix to apply */
#define NAMEDC_PLEN 500
	char prefix[NAMEDC_PLEN];	/* Temporary buffers to use */

	/* Color conversions */
	xsp2cie *sp2cie;		/* Reflectance or Transmittance to this namedc space */
	double chrom[3][3];		/* Chromatic transform to this namedc space */
	icmXYZNumber dXYZ;		/* Named color white point */

}; typedef struct _namedc namedc;

/* Create a new, uninitialised namedc */
namedc *new_namedc(a1log *log);

#endif /* NAMEDC_H */




































