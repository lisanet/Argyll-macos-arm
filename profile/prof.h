#ifndef PROF_H
#define PROF_H
/* 
 * ICC Profile creation library.
 *
 * Author:  Graeme W. Gill
 * Date:    11/10/00
 * Version: 1.00
 *
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This library provide high level routines to create device ICC
 * profiles from argyll cgats patch test data.
 */


/* Profile algorithm type */
typedef enum {
	prof_default          = 0,		/* Default for type of device */
	prof_clutLab          = 1,		/* Lab clut. */
	prof_clutXYZ          = 2,		/* XYZ clut. */
	prof_gammat           = 3,		/* XYZ gamut + matrix */
	prof_shamat           = 4,		/* XYZ shaper + matrix */
	prof_gam1mat          = 5,		/* XYZ shared TRC gamut + matrix */
	prof_sha1mat          = 6,		/* XYZ shared TRC shaper + matrix */
	prof_matonly          = 7		/* XYZ matrix, linear */
} prof_atype;

/* Output or Display device */
void make_output_icc(
	prof_atype ptype,		/* Profile output type */
	int mtxtoo,				/* NZ if matrix tags should be created for Display XYZ cLUT */
	icmICCVersion iccver,	/* ICC profile version to create */
	int verb,				/* Vebosity level, 0 = none */
	int iquality,			/* A2B table quality, 0..2 */
	int oquality,			/* B2A table quality, 0..2 */
	int noiluts,			/* nz to supress creation of input (Device) shaper luts */
	int noisluts,			/* nz to supress creation of input sub-grid (Device) shaper luts */
	int nooluts,			/* nz to supress creation of output (PCS) shaper luts */
	int nocied,				/* nz to supress inclusion of .ti3 data in profile */
	int noptop,				/* nz to use colorimetic source gamut to make perceptual table */ 
	int nostos,				/* nz to use colorimetic source gamut to make perceptual table */
	int gamdiag,			/* Make gamut mapping diagnostic wrl plots */
	int verify,				/* nz to print verification */
	int clipprims,			/* Clip white, black and primaries */
	double wpscale,			/* >= 0.0 for media white point scale factor */
//	double *bpo,			/* != NULL for XYZ black point override */
	icxInk *ink,			/* Ink limit/black generation setup */
	char *in_name,			/* input .ti3 file name */
	char *file_name,		/* output icc name */
	cgats *icg,				/* input cgats structure */
	int spec,				/* Use spectral data flag */
	icxIllumeType tillum,	/* Target/simulated instrument illuminant, if set. */ 
	xspect *cust_tillum,	/* Custom target/simulated illumination spectrum */
	                        /* if tillum == icxIT_custom */
	icxIllumeType illum,	/* CIE calc. illuminant spectrum, and FWA inst. */
							/* illuminant if tillum not set. */
	xspect *cust_illum,		/* Custom CIE illumination spectrum if illum == icxIT_custom */
	icxObserverType obType,	/* CIE calc. observer */
	xspect custObserver[3],	/* If obType = icxOT_custom */
	int fwacomp,			/* FWA compensation requested */
	double smooth,			/* RSPL smoothing factor, -ve if raw */
	double avgdev,			/* reading Average Deviation as a proportion of the input range */
	double demph,			/* Emphasise dark region grid resolution in cLUT */
	int gcompr,				/* Gamut compression % if > 0 rather than ipname */
	int gexpr,				/* Gamut saturation expansion % if gcompr > 0 rather */
	char *ipname,			/* input icc profile - enables gamut map, NULL if none */
	char *sgname,			/* source image gamut - NULL if none */
	char *absname[3],		/* abstract profile name for each table */
							/* may be duplicated, NULL if none */
	int sepsat,				/* Create separate Saturation B2A */
	icxViewCond *ivc_p,		/* Input Viewing Parameters for CIECAM97s */
	icxViewCond *ovc_p,		/* Output Viewing Parameters for CIECAM97s (enables CAM clip) */
	int ivc_e,				/* Input Enumerated viewing condition */
	int ovc_e,				/* Output Enumerated viewing condition */
	icxGMappingIntent *pgmi,/* Perceptual gamut mapping intent */
	icxGMappingIntent *sgmi,/* Saturation gamut mapping intent */
	profxinf *pi			/* Optional Profile creation extra data */
);

/* Input device */
void make_input_icc(
	prof_atype ptype,		/* Profile algorithm type */
	icmICCVersion iccver,	/* ICC profile version to create */
	int verb,
	int iquality,			/* A2B table quality, 0..3 */
	int oquality,			/* B2A table quality, 0..3 */
	int noisluts,			/* nz to supress creation of input (Device) shaper luts */
	int noipluts,			/* nz to supress creation of input (Device) position luts */
	int nooluts,			/* nz to supress creation of output (PCS) shaper luts */
	int nocied,				/* nz to supress inclusion of .ti3 data in profile */
	int verify,
	int autowpsc,			/* nz for Auto scale the WP to prevent clipping above WP patch */
	int clipovwp,			/* nz for Clip cLUT values above WP */
	double wpscale,			/* >= 0.0 for media white point scale factor */
	int dob2a,				/* nz to create a B2A table as well */
	int extrap,				/* nz to create extra cLUT interpolation points */
	int clipprims,			/* Clip white, black and primaries */
	char *in_name,			/* input .ti3 file name */
	char *file_name,		/* output icc name */
	cgats *icg,				/* input cgats structure */
	int emis,				/* emissive reference data */
	int spec,				/* Use spectral data flag */
	icxIllumeType illum,	/* Spectral illuminant */
	xspect *cust_illum,		/* Possible custom illumination */
	icxObserverType obType,	/* Spectral observer */
	xspect custObserver[3],	/* If obType = icxOT_custom */
	double smooth,			/* RSPL smoothing factor, -ve if raw */
	double avgdev,			/* reading Average Deviation as a proportion of the input range */
	profxinf *xpi			/* Optional Profile creation extra data */
);

#endif /* PROF_H */
