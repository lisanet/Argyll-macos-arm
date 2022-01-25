
/* 
 * Argyll Color Management System
 * Input device profile creator.
 *
 * Author: Graeme W. Gill
 * Date:   11/10/00
 *
 * Copyright 2000 - 2011 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in the scattered test chart
 * points, and interpolates them into a gridded
 * forward ICC device profile.
 * It also creates (at the moment) a limited backward
 * profile based on the forward grid.
 *
 */

/*
 * TTBD:
 *      Need to make this more of a library:
 *  ** By default limit matrix primaries to have +ve XYZ
 *     Add flag to override this.
 *  Fix error handling
 *  fix verbose output
 *  hand icc object back rather than writing file ?
 */

#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include "aconfig.h"
#include "counters.h"
#include "numlib.h"
#include "icc.h"
#include "cgats.h"
#include "xicc.h"
#include "rspl.h"
#include "prof.h"

#define DOB2A					/* [def] Create B2A table as well */
#define NO_B2A_PCS_CURVES       /* [def] PCS curves seem to make B2A less accurate. Why ? */
#define USE_CAM_CLIP_OPT        /* [def] Clip out of gamut in CAM space rather than XYZ or L*a*b* */
#undef USE_EXTRA_FITTING       	/* [und] Turn on data point error compensation */
#define USE_2ASS_SMOOTHING      /* [def] Turn on Gaussian smoothing */
#undef WARN_CLUT_CLIPPING		/* [und] Print warning if setting clut clips */
#define EXTRAP_MAXPNTS 10		/* [10] Maximum no. of extra extrapolated points per direction */
#define EXTRAP_WEIGHT 1.0		/* [1.0] Extra extrapolated point weighting */

/*
   Basic algorithm outline:

 Scanner:

   Figure out the input curves to give
   the flattest grid.

   Figure out the grid values.

   Use them to generate the A2B table.

   Do all the calculations in Lab space,
   but represent the profile in XYZ space, so that
   the white/black point normalisation doesn't cause
   the clut values to be clipped.

   This leads to a poorer accuracy as an XYZ profile,
   but can then be compensated for using the ICX_MERGE_CLUT flag
   together with a PCS override.

   Note we're hard coded as RGB device space, so we're not coping
   with grey scale or CMY.
*/

#ifdef DEBUG
#undef DBG
#define DBG(xxx) printf xxx ;
#else
#undef DBG
#define DBG(xxx) 
#endif

/* ---------------------------------------- */
#ifdef DOB2A

/* structure to support output icc B2A Lut initialisation calbacks */
/* Note that we don't cope with a LUT matrix - assume it's unity. */

typedef struct {
	int verb;
	int total, count, last;	/* Progress count information */
	int noPCScurves;		/* Flag set if we don't want PCS curves */
	icColorSpaceSignature pcsspace;	/* The PCS colorspace */
	icColorSpaceSignature devspace;	/* The device colorspace */
	icxLuLut *x;			/* A2B icxLuLut we are inverting in std PCS */

	double swxyz[3];		/* Source white point in XYZ */

	int wantLab;			/* 0 if want is XYZ PCS, 1 want is Lab PCS */
} in_b2a_callback;


/* --------------------------------------------------------- */

/* Extra non-linearity applied to BtoA XYZ PCS */
/* This distributes the LUT indexes more evenly in */
/* perceptual space, greatly improving the B2A accuracy of XYZ LUT */
static void xyzcurve(double *out, double *in) {
	int i;
	double sc = 65535.0/32768.0;

	/* Use an L* like curve, scaled to the maximum XYZ valu */
	out[0] = in[0]/sc;
	out[1] = in[1]/sc;
	out[2] = in[2]/sc;
	for (i = 0; i < 3; i++) {
		if (out[i] > 0.08)
			out[i] = pow((out[i] + 0.16)/1.16, 3.0);
		else
			out[i] = out[i]/9.032962896;
	}
	out[0] = out[0] * sc;
	out[1] = out[1] * sc;
	out[2] = out[2] * sc;
}

static void invxyzcurve(double *out, double *in) {
	int i;
	double sc = 65535.0/32768.0;

	out[0] = in[0]/sc;
	out[1] = in[1]/sc;
	out[2] = in[2]/sc;
	for (i = 0; i < 3; i++) {
		if (out[i] > 0.008856451586)
			out[i] = 1.16 * pow(out[i],1.0/3.0) - 0.16;
		else
			out[i] = 9.032962896 * out[i];
	}
	out[0] = out[0] * sc;
	out[1] = out[1] * sc;
	out[2] = out[2] * sc;
}

/* --------------------------------------------------------- */
/* NOTE :- the assumption that each stage of the BtoA is a mirror */
/* of the AtoB makes for inflexibility. */
/* Perhaps it would be better to remove this asumption from the */
/* in_b2a_clut processing ? */
/* To do this we then need inv_in_b2a_input(), and */
/* inv_in_b2a_output(), and we need to clearly distinguish between */
/* AtoB PCS' & DEV', and BtoA PCS' & DEV', since they are not */
/* necessarily the same... */


/* B2A Input table is the inverse of the AtoB output table */
/* Input PCS output PCS'' */
void in_b2a_input(void *cntx, double out[3], double in[3]) {
	in_b2a_callback *p = (in_b2a_callback *)cntx;

	DBG(("out_b2a_input got         PCS %f %f %f\n",in[0],in[1],in[2]))

	/* PCS to PCS' */
	if (p->noPCScurves) {
		out[0] = in[0];
		out[1] = in[1];
		out[2] = in[2];
	} else {
		if (p->x->inv_output(p->x, out, in) > 1)
			error("%d, %s",p->x->pp->errc,p->x->pp->err);
	}
	/* PCS' to PCS'' */
	if (p->pcsspace == icSigXYZData)	/* Apply XYZ non-linearity curve */
		invxyzcurve(out, out);

	DBG(("in_b2a_input returning PCS'' %f %f %f\n",out[0],out[1],out[2]))
}

/* clut - multitable */
/* Input PCS' output Dev' */
/* We're applying any abstract profile after gamut mapping, */
/* on the assumption is is primarily being used to "correct" the */
/* output device. Ideally the gamut mapping should take the change */
/* the abstract profile has on the output device into account, but */
/* currently we're not doing this.. */
void in_b2a_clut(void *cntx, double *out, double in[3]) {
	in_b2a_callback *p = (in_b2a_callback *)cntx;
	double in1[3];

	in1[0] = in[0];		/* in[] may be aliased with out[] */
	in1[1] = in[1];		/* so take a copy.  */
	in1[2] = in[2];

	DBG(("in_b2a_clut got       PCS' %f %f %f\n",in[0],in[1],in[2]))

	if (p->pcsspace == icSigXYZData)		/* Undo effects of extra XYZ non-linearity curve */
		xyzcurve(in1, in1);

	if (p->noPCScurves) {	/* We were given PCS or have converted to PCS */

		/* PCS to PCS' */
		if (p->x->inv_output(p->x, in1, in1) > 1)
			error("%d, %s",p->x->pp->errc,p->x->pp->err);

		DBG(("convert to PCS' got         %f %f %f\n",in1[0],in1[1],in1[2]))
	}

	/* Invert AtoB clut (PCS' to Dev') Colorimetric */
	/* to producte the colorimetric tables output. */
	if (p->x->inv_clut(p->x, out, in1) > 1)
		error("%d, %s",p->x->pp->errc,p->x->pp->err);

	DBG(("convert PCS' to DEV' got    %f %f %f %f\n",out[0],out[1],out[2],out[3]))
	DBG(("in_b2a_clut returning DEV' %f %f %f\n",out[0],out[1],out[2]))

	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = (int)(p->count * 100.0/p->total + 0.5);
		if (pc != p->last) {
			printf("%c%2d%%",cr_char,pc); fflush(stdout);
			p->last = pc;
		}
	}
}

/* Output table is the inverse of the AtoB input table */
/* Input Dev' output Dev */
void in_b2a_output(void *cntx, double out[4], double in[4]) {
	in_b2a_callback *p = (in_b2a_callback *)cntx;

	DBG(("in_b2a_output got      DEV' %f %f %f\n",in[0],in[1],in[2]))

	if (p->x->inv_input(p->x, out, in) > 1)
		error("%d, %s",p->x->pp->errc,p->x->pp->err);

	DBG(("in_b2a_output returning DEV %f %f %f\n",out[0],out[1],out[2]))
}

#else /* !DOB2A */
# pragma message("!!!!!!! profin DOB2A not defined !!!!!!!!")
#endif /* !DOB2A */
/* ---------------------------------------- */

/* Make an input device profile, where we create an A2B lut */
/* directly from the scattered input data. */
void
make_input_icc(
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
	int autowpsc,			/* 1 for Auto scale the WP to prevent clipping above WP patch */
	                        /* 2 = Force absolute colorimetric */
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
) {
	icmFile *wr_fp;
	icc *wr_icco;
	int npat;				/* Number of patches */
	int npxpat = 0;			/* Number of possible extrap extrapolation patches */
	int nxpat = 0;			/* Number of extrap extrapolation patches */
	cow *tpat;				/* Patch input values */
	int i, rv = 0;
	int isLab = 0;			/* 0 if input is XYZ, 1 if input is Lab */
	int wantLab = 0;		/* 0 if want is XYZ, 1 want is Lab. */
							/* Values will be wantLab after reading */
	int isLut = 0;			/* 0 if shaper+ matrix, 1 if lut type */
	int isShTRC = 0;		/* 0 if separate gamma/shaper TRC, 1 if shared */

	if (ptype == prof_clutLab) {		/* Lab lut */
		wantLab = 1;
		isLut = 1;
	} else if (ptype == prof_clutXYZ) {	/* XYZ lut */
		wantLab = 0;
		isLut = 1;
	} else {
		wantLab = 0;			/* gamma/shaper + matrix profile must be XYZ */
		isLut = 0;
		extrap = 0;

		if (ptype == prof_gam1mat	
		 || ptype == prof_sha1mat
		 || ptype == prof_matonly) {
			isShTRC = 1;		/* Single curve */
		}
	}

	/* Open up the file for writing */
	if ((wr_fp = new_icmFileStd_name(file_name,"w")) == NULL)
		error ("Write: Can't open file '%s'",file_name);

	if ((wr_icco = new_icc()) == NULL)
		error ("Write: Creation of ICC object failed");

	/* Add all the tags required */

	/* The header: */
	{
		icmHeader *wh = wr_icco->header;

		/* Values that must be set before writing */
		wh->deviceClass     = icSigInputClass;
    	wh->colorSpace      = icSigRgbData;				/* It's an RGB profile */
		if (wantLab)
	    	wh->pcs         = icSigLabData;
		else
	    	wh->pcs         = icSigXYZData;

		if (xpi->default_ri != icMaxEnumIntent)
	    	wh->renderingIntent = xpi->default_ri;
		else
	    	wh->renderingIntent = icRelativeColorimetric;

		/* Values that should be set before writing */
		if (xpi != NULL && xpi->manufacturer != 0L)
			wh->manufacturer = xpi->manufacturer;
		else
			wh->manufacturer = icmSigUnknownType;

		if (xpi != NULL && xpi->model != 0L)
			wh->model = xpi->model;
		else
	    	wh->model = icmSigUnknownType;

		/* Values that may be set before writing */
		if (xpi != NULL && xpi->creator != 0L)
			wh->creator = xpi->creator;
#ifdef NT
		wh->platform = icSigMicrosoft;
#endif
#ifdef __APPLE__
		wh->platform = icSigMacintosh;
#endif
#if defined(UNIX) && !defined(__APPLE__)
		wh->platform = icmSig_nix;
#endif

		if (xpi != NULL && xpi->transparency)
			wh->attributes.l |= icTransparency;
		if (xpi != NULL && xpi->matte)
			wh->attributes.l |= icMatte;
		if (xpi != NULL && xpi->negative)
			wh->attributes.l |= icNegative;
		if (xpi != NULL && xpi->blackandwhite)
			wh->attributes.l |= icBlackAndWhite;
	}
	/* Profile Description Tag: */
	{
		icmTextDescription *wo;
		char *dst;			/* description */

		if (xpi != NULL && xpi->profDesc != NULL)
			dst = xpi->profDesc;
		else {
			dst = "This is a Lut style RGB - XYZ Input Profile";
		}

		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigProfileDescriptionTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* Copyright Tag: */
	{
		icmText *wo;
		char *crt;

		if (xpi != NULL && xpi->copyright != NULL)
			crt = xpi->copyright;
		else
			crt = "Copyright, the creator of this profile";

		if ((wo = (icmText *)wr_icco->add_tag(
		           wr_icco, icSigCopyrightTag,	icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->data, crt);		/* Copy the text in */
	}
	/* Device Manufacturers Description Tag: */
	if (xpi != NULL && xpi->deviceMfgDesc != NULL) {
		icmTextDescription *wo;
		char *dst = xpi->deviceMfgDesc;

		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigDeviceMfgDescTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* Model Description Tag: */
	if (xpi != NULL && xpi->modelDesc != NULL) {
		icmTextDescription *wo;
		char *dst = xpi->modelDesc;

		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigDeviceModelDescTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* White Point Tag: */
	{
		icmXYZArray *wo;
		/* Note that tag types icSigXYZType and icSigXYZArrayType are identical */
		if ((wo = (icmXYZArray *)wr_icco->add_tag(
		           wr_icco, icSigMediaWhitePointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0].X = 0.9642;		/* Set a default value - D50 */
		wo->data[0].Y = 1.0000;
		wo->data[0].Z = 0.8249;
	}
	/* Black Point Tag: */
	{
		icmXYZArray *wo;
		if ((wo = (icmXYZArray *)wr_icco->add_tag(
		           wr_icco, icSigMediaBlackPointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0].X = 0.00;			/* Set default perfect black */
		wo->data[0].Y = 0.00;
		wo->data[0].Z = 0.00;
	}

	if (isLut == 0) {	/* shaper + matrix type */

		/* Red, Green and Blue Colorant Tags: */
		{
			icmXYZArray *wor, *wog, *wob;
			if ((wor = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigRedColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);
			if ((wog = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigGreenColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);
			if ((wob = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigBlueColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);

			wor->size = wog->size = wob->size = 1;
			wor->allocate((icmBase *)wor);	/* Allocate space */
			wog->allocate((icmBase *)wog);
			wob->allocate((icmBase *)wob);

			/* Setup some sane dummy values */
			/* icxMatrix will override these later */
			wor->data[0].X = 1.0; wor->data[0].Y = 0.0; wor->data[0].Z = 0.0;
			wog->data[0].X = 0.0; wog->data[0].Y = 1.0; wog->data[0].Z = 0.0;
			wob->data[0].X = 0.0; wob->data[0].Y = 0.0; wob->data[0].Z = 1.0;
		}

		/* Red, Green and Blue Tone Reproduction Curve Tags: */
		{
			icmCurve *wor, *wog, *wob;
			if ((wor = (icmCurve *)wr_icco->add_tag(
			           wr_icco, icSigRedTRCTag, icSigCurveType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);

			if (isShTRC) {	/* Make all TRCs shared */
				if ((wog = (icmCurve *)wr_icco->link_tag(
				           wr_icco, icSigGreenTRCTag, icSigRedTRCTag)) == NULL) 
					error("link_tag failed: %d, %s",rv,wr_icco->err);
				if ((wob = (icmCurve *)wr_icco->link_tag(
				           wr_icco, icSigBlueTRCTag, icSigRedTRCTag)) == NULL) 
					error("link_tag failed: %d, %s",rv,wr_icco->err);

			} else {		/* Else individual */
				if ((wog = (icmCurve *)wr_icco->add_tag(
				           wr_icco, icSigGreenTRCTag, icSigCurveType)) == NULL) 
					error("add_tag failed: %d, %s",rv,wr_icco->err);
				if ((wob = (icmCurve *)wr_icco->add_tag(
				           wr_icco, icSigBlueTRCTag, icSigCurveType)) == NULL) 
					error("add_tag failed: %d, %s",rv,wr_icco->err);
			}
	
			if (ptype == prof_shamat || ptype == prof_sha1mat) {	/* Shaper */
				wor->flag = wog->flag = wob->flag = icmCurveSpec; 
				wor->size = wog->size = wob->size = 256;			/* Number of entries */
			} else {						/* Gamma */
				wor->flag = wog->flag = wob->flag = icmCurveGamma;
				wor->size = wog->size = wob->size = 1;				/* Must be 1 for gamma */
			}
			wor->allocate((icmBase *)wor);	/* Allocate space */
			wog->allocate((icmBase *)wog);
			wob->allocate((icmBase *)wob);

			/* icxMatrix will set curve values */
		}

	} else {		/* Lut type profile */

		/* 16 bit dev -> pcs lut: */
		{
			icmLut *wo;

			/* Only A2B0, no intent */
			if ((wo = (icmLut *)wr_icco->add_tag(
			           wr_icco, icSigAToB0Tag,	icSigLut16Type)) == NULL) 
				error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			wo->inputChan = 3;
			wo->outputChan = 3;
			if (iquality >= 3) {
		    	wo->clutPoints = 45;
		    	wo->inputEnt = 2048;
		    	wo->outputEnt = 2048;
			} else if (iquality == 2) {
		    	wo->clutPoints = 33;
		    	wo->inputEnt = 2048;
		    	wo->outputEnt = 2048;
			} else if (iquality == 1) {
		    	wo->clutPoints = 17;
		    	wo->inputEnt = 1024;
		    	wo->outputEnt = 1024;
			} else {
		    	wo->clutPoints = 9;
		    	wo->inputEnt = 512;
		    	wo->outputEnt = 512;
			}

			wo->allocate((icmBase *)wo);/* Allocate space */

			/* icxLuLut will set tables values */
		}

#ifdef DOB2A
		/* 16 bit pcs -> dev lut: */
		if (dob2a) {
			icmLut *wo;

			/* Only B2A0, no intent */
			if ((wo = (icmLut *)wr_icco->add_tag(
			           wr_icco, icSigBToA0Tag,	icSigLut16Type)) == NULL) 
				error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			wo->inputChan = 3;
			wo->outputChan = 3;
			if (oquality >= 3) {
		    	wo->clutPoints = 45;
		    	wo->inputEnt = 2048;
		    	wo->outputEnt = 2048;
			} else if (oquality == 2) {
		    	wo->clutPoints = 33;
		    	wo->inputEnt = 2048;
		    	wo->outputEnt = 2048;
			} else if (oquality == 1) {
		    	wo->clutPoints = 17;
		    	wo->inputEnt = 1024;
		    	wo->outputEnt = 1024;
            } else if (oquality >= 0) {
		    	wo->clutPoints = 9;
		    	wo->inputEnt = 512;
		    	wo->outputEnt = 512;
            } else {                /* Special, Extremely low quality */
		    	wo->clutPoints = 3;
		    	wo->inputEnt = 64;
		    	wo->outputEnt = 64;
            }

			wo->allocate((icmBase *)wo);/* Allocate space */

			/* We set the tables below */
		}
#endif /* DOB2A */

	}

	/* Sample data use to create profile: */
	if (nocied == 0) {
		icmText *wo;
		char *crt;
		FILE *fp;

		if ((wo = (icmText *)wr_icco->add_tag(
		           wr_icco, icmMakeTag('t','a','r','g'), icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

#if defined(O_BINARY) || defined(_O_BINARY)
	    if ((fp = fopen(in_name, "rb")) == NULL)
#else
	    if ((fp = fopen(in_name, "r")) == NULL)
#endif
			error("Unable to open input file '%s' for reading",in_name);

		if (fseek(fp, 0, SEEK_END))
			error("Unable to seek to end of file '%s'",in_name);
		wo->size = ftell(fp) + 1;		/* Size needed + null */
		wo->allocate((icmBase *)wo);/* Allocate space */

		if (fseek(fp, 0, SEEK_SET))
			error("Unable to seek to end of file '%s'",in_name);

		if (fread(wo->data, 1, wo->size-1, fp) != wo->size-1)
			error("Failed to read file '%s'",in_name);
		wo->data[wo->size-1] = '\000';
		fclose(fp);

		/* Duplicate for compatibility */
		if (wr_icco->link_tag(
		         wr_icco, icmMakeTag('D','e','v','D'), icmMakeTag('t','a','r','g')) == NULL) 
			error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
		if (wr_icco->link_tag(
		         wr_icco, icmMakeTag('C','I','E','D'), icmMakeTag('t','a','r','g')) == NULL) 
			error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
	}

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	if (verb) {
		fprintf(verbo,"No of test patches = %d\n",npat);
	}

	if (extrap) {
		npxpat = 4 * EXTRAP_MAXPNTS;		/* Allow for up to 20 extra patches */
	}

	/* Allocate arrays to hold test patch input and output values */
	if ((tpat = (cow *)malloc(sizeof(cow) * (npat + npxpat))) == NULL)
		error("Malloc failed - tpat[]");

	/* Read in the CGATs fields */
	{
		int ti;
		int Xi, Yi, Zi;
		int ri, gi, bi;

		/* Check that we handle the color space */
		if ((ti = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
			error ("Input file doesn't contain keyword COLOR_REPS");
		if (strcmp(icg->t[0].kdata[ti],"LAB_RGB") == 0) {
			isLab = 1;
		} else {
			if (strcmp(icg->t[0].kdata[ti],"XYZ_RGB") == 0) {
				isLab = 0;
			} else {
				error ("Input device input file has unhandled color representation");
			}
		}

		if ((ri = icg->find_field(icg, 0, "RGB_R")) < 0)
			error ("Input file doesn't contain field RGB_R");
		if (icg->t[0].ftype[ri] != r_t)
			error ("Field RGB_R is wrong type - expect float");
		if ((gi = icg->find_field(icg, 0, "RGB_G")) < 0)
			error ("Input file doesn't contain field RGB_G");
		if (icg->t[0].ftype[gi] != r_t)
			error ("Field RGB_G is wrong type - expect float");
		if ((bi = icg->find_field(icg, 0, "RGB_B")) < 0)
			error ("Input file doesn't contain field RGB_B");
		if (icg->t[0].ftype[bi] != r_t)
			error ("Field RGB_B is wrong type - expect float");

		if (spec == 0) {        /* Using instrument tristimulous value */

			if (isLab) {
				if ((Xi = icg->find_field(icg, 0, "LAB_L")) < 0)
					error ("Input file doesn't contain field LAB_L");
				if (icg->t[0].ftype[Xi] != r_t)
					error ("Field LAB_L is wrong type - expect float");
				if ((Yi = icg->find_field(icg, 0, "LAB_A")) < 0)
					error ("Input file doesn't contain field LAB_A");
				if (icg->t[0].ftype[Yi] != r_t)
					error ("Field LAB_A is wrong type - expect float");
				if ((Zi = icg->find_field(icg, 0, "LAB_B")) < 0)
					error ("Input file doesn't contain field LAB_B");
				if (icg->t[0].ftype[Zi] != r_t)
					error ("Field LAB_B is wrong type - expect float");
			} else {
				if ((Xi = icg->find_field(icg, 0, "XYZ_X")) < 0)
					error ("Input file doesn't contain field XYZ_X");
				if (icg->t[0].ftype[Xi] != r_t)
					error ("Field XYZ_X is wrong type - expect float");
				if ((Yi = icg->find_field(icg, 0, "XYZ_Y")) < 0)
					error ("Input file doesn't contain field XYZ_Y");
				if (icg->t[0].ftype[Yi] != r_t)
					error ("Field XYZ_Y is wrong type - expect float");
				if ((Zi = icg->find_field(icg, 0, "XYZ_Z")) < 0)
					error ("Input file doesn't contain field XYZ_Z");
				if (icg->t[0].ftype[Zi] != r_t)
					error ("Field XYZ_Z is wrong type - expect float");
			}

			for (i = 0; i < npat; i++) {
				tpat[i].w = 1.0;
				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ri]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][gi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][bi]) / 100.0;
				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0) {
					error("At %d device values %f %f %f field exceeds 100.0!",i,100.0 * tpat[i].p[0],100.0 * tpat[i].p[1],100.0 * tpat[i].p[2]);
				}
				tpat[i].v[0] = *((double *)icg->t[0].fdata[i][Xi]);
				tpat[i].v[1] = *((double *)icg->t[0].fdata[i][Yi]);
				tpat[i].v[2] = *((double *)icg->t[0].fdata[i][Zi]);
				if (!isLab) {
					tpat[i].v[0] /= 100.0;		/* Normalise XYZ to range 0.0 - 1.0 */
					tpat[i].v[1] /= 100.0;
					tpat[i].v[2] /= 100.0;
				}
				if (!isLab && wantLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
					icmXYZ2Lab(&icmD50, tpat[i].v, tpat[i].v);
				} else if (isLab && !wantLab) {
					icmLab2XYZ(&icmD50, tpat[i].v, tpat[i].v);
				}
			}

		} else {		/* Using spectral data */
			int j, ii;
			xspect sp;
			char buf[100];
			int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
			xsp2cie *sp2cie;	/* Spectral conversion object */

			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_BANDS")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_BANDS");
			sp.spec_n = atoi(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_START_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_START_NM");
			sp.spec_wl_short = atof(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_END_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_END_NM");
			sp.spec_wl_long = atof(icg->t[0].kdata[ii]);
			sp.norm = 100.0;

			/* Find the fields for spectral values */
			for (j = 0; j < sp.spec_n; j++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
				            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);

				if ((spi[j] = icg->find_field(icg, 0, buf)) < 0)
					error("Input file doesn't contain field %s",buf);

				if (icg->t[0].ftype[spi[j]] != r_t)
					error("Field %s is wrong type - expect float",buf);
			}

			/* Create a spectral conversion object */
			if (emis) {
				illum = icxIT_none;
				cust_illum = NULL;
			}
			if ((sp2cie = new_xsp2cie(illum, 0.0, cust_illum, obType, custObserver,
			                          wantLab ? icSigLabData : icSigXYZData, icxClamp)) == NULL)
				error("Creation of spectral conversion object failed");

			for (i = 0; i < npat; i++) {
				tpat[i].w = 1.0;
				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ri]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][gi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][bi]) / 100.0;
				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0) {
					error("Patch %d device value field exceeds 100.0!",i+1);
				}

				/* Read the spectral values for this patch */
				for (j = 0; j < sp.spec_n; j++) {
					sp.spec[j] = *((double *)icg->t[0].fdata[i][spi[j]]);
				}

				/* Convert it to CIE space */
				sp2cie->convert(sp2cie, tpat[i].v, &sp);
			}

			sp2cie->del(sp2cie);		/* Done with this */

		}
	}	/* End of reading in CGATs file */

	if (isLut == 0) { /* Gamma/Shaper + matrix profile */
		xicc *wr_xicc;			/* extention object */
		icxLuBase *xluo;		/* Forward ixcLu */
		int flags = 0;

		/* Wrap with an expanded icc */
		if ((wr_xicc = new_xicc(wr_icco)) == NULL)
			error("Creation of xicc failed");
		
		if (verb)
			flags |= ICX_VERBOSE;

		if (ptype == prof_matonly)
			flags |= ICX_NO_IN_SHP_LUTS;	/* Make it linear */

		if (clipprims)
			flags |= ICX_CLIP_WB | ICX_CLIP_PRIMS;
				
        flags |= ICX_SET_BLACK;		/* Compute & use black */
		flags |= ICX_SET_WHITE;		/* Compute & use white */

		/* ICX_SET_WHITE_C isn't applicable to matrix profiles */
		if (autowpsc == 1)
	        flags |= ICX_SET_WHITE_US;	/* Compute & use white without scaling to L */
		else if (autowpsc == 2)
	        flags |= ICX_SET_WHITE_ABS;	/* Set dummy D50 white point to force absolute intent */

        flags |= ICX_WRITE_WBL;		/* Matrix: write white/black/luminence */

		/* Setup Device -> XYZ conversion (Fwd) object from scattered data. */
		if ((xluo = wr_xicc->set_luobj(
//		               wr_xicc, icmFwd, icRelativeColorimetric,
		               wr_xicc, icmFwd, icmDefaultIntent,
		               icmLuOrdNorm,
		               flags, 		/* Flags */
		               npat, npat, tpat, NULL, 0.0, wpscale,
//		               NULL,		/* bpo */
			           smooth, avgdev, 1.0,
		               NULL, NULL, NULL, iquality)) == NULL)
			error("%d, %s",wr_xicc->errc, wr_xicc->err);

		/* Free up xicc stuff */
		xluo->del(xluo);
		wr_xicc->del(wr_xicc);

	} else {		/* cLUT based profile */
		int flags = 0;
		icxMatrixModel *mm = NULL;

		xicc *wr_xicc;			/* extention object */
		icxLuBase *AtoB;		/* AtoB ixcLu */

		if (extrap) {
			cow *mpat;
			int nmpat;
			int j;

			if (verb) printf("Creating extrapolation black and white points:\n");

			if ((mpat = (cow *)malloc(sizeof(cow) * npat)) == NULL)
				error("Malloc failed - mpat[]");

			/* Weight points from full set to build matrix model */
			/* to extrapolate the neutral axis */
			for (nmpat = j = 0; j < npat; j++) {
				double mnp, mxp;
				int k;

				icmCpy3(mpat[nmpat].p, tpat[j].p);
				icmCpy3(mpat[nmpat].v, tpat[j].v);

				/* Locate largest/smallest RGB value */
				mxp = -1e6, mnp = 1e6;
				for (k = 0; k < 3; k++) {
					if (tpat[j].p[k] > mxp)
						mxp = tpat[j].p[k];
					if (tpat[j].p[k] < mnp)
						mnp = tpat[j].p[k];
				}
				mxp -= mnp;			/* Spread; 0 for R=G=B */

				mpat[nmpat].w = pow(1.1 - mxp, 2.0);
//printf("~1 added value %d: %f %f %f -> %f %f %f wt %f\n",j, mpat[nmpat].p[0], mpat[nmpat].p[1], mpat[nmpat].p[2], mpat[nmpat].v[0], mpat[nmpat].v[1], mpat[nmpat].v[2],mpat[nmpat].w);
				nmpat++;
			}
	
			/* Create gamma/matrix model to extrapolate with. */
			/* (Use ofset & gain, gamma curve as 0th order with 1 harmonic, */
			/* and smooth it.) */
			if ((mm = new_MatrixModel(wr_icco, verb, nmpat, mpat, wantLab,
				      /* quality */ -1, /* isLinear */ ptype == prof_matonly,
				      /* isGamma */ 0, /* isShTRC */ 0,
				      /* shape0gam */ 1, /* clipbw */ 0, /* clipprims */ 0,
//				      /* smooth */ 1.0, /* scale */ 0.7)) == NULL) {
				      /* smooth */ 1.0, /* scale */ 1.0)) == NULL) {
				error("Creating extrapolation matrix model failed - memory ?");
			}

#ifdef NEVER	/* Plot Lab of model */
{
	#define	XRES 100
	double xx[XRES];
	double y0[XRES];
	double y1[XRES];
	double y2[XRES];

	/* Display the result fit */
	for (i = 0; i < XRES; i++) {
		double rgb[3], lab[3];
		xx[i] = rgb[0] = rgb[1] = rgb[2] = i/(double)(XRES-1);
		mm->lookup(mm, lab, rgb);
		if (wantLab)
			icmLab2XYZ(&icmD50,lab,lab);
		y0[i] = lab[0];
		y1[i] = lab[1];
		y2[i] = lab[2];
	}
	do_plot(xx,y0,y1,y2,XRES);
}
#endif /* DEBUG_PLOT */
		}

		if (extrap) {
			int ii, wix = 0, j;
			int pcsy;						/* Effective PCS L or Y chanel index */
			double wpy = -1e60;
			double dwhite[MXDI];  /* Device white */
			double mxdw;
			double avgdist;		/* Average distance between points */

			/* Figure out the device white point. */
			/* Note that this is duplicating code in xicc/xmatrix.c */
			/* and xfit.c */

			if (wantLab)
				pcsy = 0;	/* L or Lab */
			else
				pcsy = 1;	/* Y of XYZ */

			for (i = 0; i < npat; i++) {
				double labv[3], yv;

				/* Create D50 Lab to allow some chromatic sensitivity */
				/* in picking the white point */
				if (wantLab)
					icmCpy3(labv, tpat[i].v);
				else
					icmXYZ2Lab(&icmD50, labv, tpat[i].v);

				/* Tilt things towards D50 neutral white patches */
				yv = labv[0] - 0.3 * sqrt(labv[1] * labv[1] + labv[2] * labv[2]);
				if (yv > wpy) {
					for (j = 0; j < 3; j++)
						dwhite[j] = tpat[i].p[j];
					wpy = yv;
					wix = i;
				}
			}

			/* Fix extrapolation matrix to be perfect at white point */
			mm->force(mm, tpat[wix].v, tpat[wix].p);

			/* Scale the white point to make one dev value 1.0 */
			mxdw = -1;
			for (j = 0; j < 3; j++) {
				if (dwhite[j] > mxdw)
					mxdw = dwhite[j];
			}
			for (j = 0; j < 3; j++) {
				dwhite[j] /= mxdw;
			}

			avgdist = pow(1.0/(double)npat, 1.0/3.0);
			if (avgdist < 0.001)
				avgdist = 0.001;
			else if (avgdist > 0.3)
				avgdist = 0.3;
//printf("~1 avgdist = %f\n",avgdist);

			/* For points with white point device ratio, */
			/* and points with R=G=B ratio, create extrapolation points. */
			for (ii = 0; ii < 2; ii++) {

				if (ii > 1)
					dwhite[0] = dwhite[1] = dwhite[2] = 1.0;

				/* Create a series of black and white patch */
				for (i = 0; i < 2; i++) {
					int cix;					/* Closest point index */
					int eix;					/* End point index */
					double cde = 1e60;			/* Closest point distance */
					double tt;
					double corr[3], cwt;		/* Correction */

					eix = npat + nxpat;

					icmScale3(tpat[eix].p, dwhite, (double)i);

					/* Locate closest point */
					for (j = 0; j < npat; j++) {
						double mnp, mxp;
						int k;
		
						/* Locate largest/smallest RGB value */
						mxp = -1e6, mnp = 1e6;
						for (k = 0; k < 3; k++) {
							if (tpat[j].p[k] > mxp)
								mxp = tpat[j].p[k];
							if (tpat[j].p[k] < mnp)
								mnp = tpat[j].p[k];
						}
						mxp -= mnp;			/* Spread; 0 for R=G=B */

						tt = icmNorm33(tpat[eix].p, tpat[j].p);
						tt += mxp;

						if (tt < cde) {
							cde = tt;
							cix = j;
						}
					}

//printf("~1 closest %d: de %f, %f %f %f -> %f %f %f\n",cix, cde, tpat[cix].p[0], tpat[cix].p[1], tpat[cix].p[2], tpat[cix].v[0], tpat[cix].v[1], tpat[cix].v[2]);

//{
//double val[3];
//mm->lookup(mm, val, tpat[cix].p);
//printf("~1 closest gam/matrix -> %f %f %f\n",val[0],val[1],val[2]);
//}

					/* Lookup matrix value for our new point */
					mm->lookup(mm, tpat[eix].v, tpat[eix].p);
//printf("~1 got value %d: %f %f %f -> %f %f %f\n",i, tpat[eix].p[0], tpat[eix].p[1], tpat[eix].p[2], tpat[eix].v[0], tpat[eix].v[1], tpat[eix].v[2]);
					/* Weight the extra point so that it doesn't overpower the */
					/* nearest real point to it too much. */
					tt = cde;
					if (tt > avgdist)		/* Distance at which sythetic point has 100% weight */
						tt = avgdist;
					tpat[eix].w = 0.5 * EXTRAP_WEIGHT * tt/avgdist;	
//printf("~1 weight %f\n",tpat[eix].w);
					if (verb)
						printf("Added synthetic point @ %f %f %f, val %f %f %f, weight %f\n",tpat[eix].p[0], tpat[eix].p[1], tpat[eix].p[2], tpat[eix].v[0], tpat[eix].v[1], tpat[eix].v[2],tpat[eix].w);
					nxpat++;
					
					/* If there is a lot of space, add a second intemediate point */
//printf("~1 cde = %f, avgdist = %f\n",cde,avgdist);
					if (cde >= (0.5 * avgdist)) {
						int nxps;				/* Number of extra points including end point */
						nxps = 1 + (int)(cde/(0.5 * avgdist));
						if (nxps > EXTRAP_MAXPNTS)
							nxps = EXTRAP_MAXPNTS;

//printf("~1 nxps = %d\n",nxps);
						for (j = 1; j < nxps; j++) {
							double bl, ipos;
		
							bl = j/(nxps + 1.0);

							ipos = (1.0 - bl) * tpat[eix].p[0]
							     +        bl * (tpat[cix].p[0] + tpat[cix].p[1] + tpat[cix].p[1])/3.0;
							icmScale3(tpat[npat + nxpat].p, dwhite, ipos);
			
							/* Lookup matrix value for our new point */
							mm->lookup(mm, tpat[npat + nxpat].v, tpat[npat + nxpat].p);
			
							/* Weight the extra point so that it doesn't overpower the */
							/* nearest real point to it too much. */
							cde = icmNorm33(tpat[cix].p, tpat[npat + nxpat].p);
			
							if (cde > avgdist)		/* Distance at which sythetic point has 100% weight */
								cde = avgdist;
							tpat[npat + nxpat].w = 0.5 * EXTRAP_WEIGHT * cde/avgdist;	
							if (verb)
								printf("Added synthetic point @ %f %f %f, val %f %f %f, weight %f\n",tpat[npat + nxpat].p[0], tpat[npat + nxpat].p[1], tpat[npat + nxpat].p[2], tpat[npat + nxpat].v[0], tpat[npat + nxpat].v[1], tpat[npat + nxpat].v[2],tpat[npat + nxpat].w);
							nxpat++;
						}
					}
				}
			}
		}

		/* Wrap with an expanded icc */
		if ((wr_xicc = new_xicc(wr_icco)) == NULL)
			error ("Creation of xicc failed");

		flags |= ICX_CLIP_NEAREST;      /* This will avoid clip caused rev setup */

		if (noisluts)
			flags |= ICX_NO_IN_SHP_LUTS;

		if (noipluts)
			flags |= ICX_NO_IN_POS_LUTS;

		if (nooluts)
			flags |= ICX_NO_OUT_LUTS;

		if (verb)
			flags |= ICX_VERBOSE;

		if (clipprims)
			flags |= ICX_CLIP_WB;
				
        flags |= ICX_SET_BLACK;		/* Compute & use black */
		flags |= ICX_SET_WHITE;		/* Compute & use white */
		if (clipovwp)
	        flags |= ICX_SET_WHITE_C;	/* Compute & use white and clip cLUT over D50 */
		else if (autowpsc == 1)
	        flags |= ICX_SET_WHITE_US;	/* Compute & use white without scaling to L */
		else if (autowpsc == 2)
	        flags |= ICX_SET_WHITE_ABS;	/* Set dummy D50 white point to force absolute intent */

		/* Setup RGB -> Lab conversion object from scattered data. */
		/* Note that we've layered it on a native XYZ icc profile. */
		/* (The skeleton model is not used - it doesn't seem to help) */
		if ((AtoB = wr_xicc->set_luobj(
		               wr_xicc, icmFwd, icmDefaultIntent,
		               icmLuOrdNorm,
#ifdef USE_EXTRA_FITTING
			               ICX_EXTRA_FIT |
#endif
#ifdef USE_2PASSSMTH
			               ICX_2PASSSMTH |
#endif
		               flags, 		/* Flags */
		               npat + nxpat, npat, tpat, NULL, 0.0, wpscale,
//			           NULL,	/* bpo */
			           smooth, avgdev, 1.0,
			           NULL, NULL, NULL, iquality)) == NULL)
			error ("%d, %s",wr_xicc->errc, wr_xicc->err);

		if (mm != NULL)
			mm->del(mm);

		/* Free up xicc stuff */
		AtoB->del(AtoB);

#ifdef DOB2A
		if (dob2a) {
			icmLut *wo;

			in_b2a_callback cx;

			if (verb)
				printf("Setting up B to A table lookup\n");

			/* Get a suitable forward conversion object to invert. */
			/* By creating a separate one to the one created using scattered data, */
			/* we ge the chance to set ICX_CAM_CLIP. It is always set to Lab 'PCS' */
			{
				int flags = 0;
	
				if (verb)
					flags |= ICX_VERBOSE;
	
				flags |= ICX_CLIP_NEAREST;		/* Not vector clip */

#ifdef USE_CAM_CLIP_OPT
				flags |= ICX_CAM_CLIP;			/* Clip in CAM Jab space rather than Lab */
#else
				warning("!!!! USE_CAM_CLIP_OPT in profout.c is off !!!!");
#endif
				if ((AtoB = wr_xicc->get_luobj(wr_xicc, flags, icmFwd,
				                  icmDefaultIntent,
				                  wantLab ? icSigLabData : icSigXYZData,
                                  icmLuOrdNorm, NULL, NULL)) == NULL)
					error ("%d, %s",wr_xicc->errc, wr_xicc->err);
			}

			/* setup context ready for B2A table setting */
			cx.verb = verb;
			cx.pcsspace = wantLab ? icSigLabData : icSigXYZData;
			cx.wantLab = wantLab;			/* Copy PCS flag over */
#ifdef NO_B2A_PCS_CURVES
			cx.noPCScurves = 1;		/* Don't use PCS curves */
#else
			cx.noPCScurves = 0;
#endif
			cx.devspace = icSigRgbData;
			cx.x = (icxLuLut *)AtoB;		/* A2B icxLuLut created from scattered data */

			if ((wo = (icmLut *)wr_icco->read_tag(
			           wr_icco, icSigBToA0Tag)) == NULL) 
				error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			/* We now setup an exact inverse, colorimetric style */
			/* Use helper function to do the hard work. */

			if (cx.verb) {
				unsigned int ui;
				int extra;
				cx.count = 0;
				cx.last = -1;
				for (cx.total = 1, ui = 0; ui < wo->inputChan; ui++, cx.total *= wo->clutPoints)
					; 
				/* Add in cell center points */
				for (extra = 1, ui = 0; ui < wo->inputChan; ui++, extra *= (wo->clutPoints-1))
					;
				cx.total += extra;
				printf("Creating B to A tables\n");
				printf(" 0%%"); fflush(stdout);
			}

			if (icmSetMultiLutTables(
			        1,
			        &wo,
					ICM_CLUT_SET_APXLS,			/* Use least squared aprox. */
					&cx,						/* Context */
					cx.pcsspace,				/* Input color space */
					icSigRgbData,				/* Output color space */
					in_b2a_input,				/* Input transform PCS->PCS' */
					NULL, NULL,					/* Use default Lab' range */
					in_b2a_clut,				/* Lab' -> Device' transfer function */
					NULL, NULL,					/* Use default Device' range */
					in_b2a_output,				/* Output transfer function, Device'->Device */
					NULL, NULL) != 0)			/* Use default APXLS range */
				error("Setting 16 bit PCS->Device Lut failed: %d, %s",wr_icco->errc,wr_icco->err);
			if (cx.verb) {
				printf("\n");
			}
#ifdef WARN_CLUT_CLIPPING
			if (wr_icco->warnc)
				warning("Values clipped in setting LUT");
#endif /* WARN_CLUT_CLIPPING */

			if (verb)
				printf("Done B to A table\n");
			AtoB->del(AtoB);
		}
#endif /* DOB2A */
		wr_xicc->del(wr_xicc);

	}

	/* Write the file (including all tags) out */
	if ((rv = wr_icco->write(wr_icco,wr_fp,0)) != 0)
		error ("Write file: %d, %s",rv,wr_icco->err);

	/* Close the file */
	wr_icco->del(wr_icco);
	wr_fp->del(wr_fp);

	/* Check the profile accuracy against the data points */
	if (verb || verify) {
		icmFile *rd_fp;
		icc *rd_icco;
		icmLuBase *luo;
		double merr = 0.0;
		double aerr = 0.0;
		double nsamps = 0.0;

		/* Open up the file for reading */
		if ((rd_fp = new_icmFileStd_name(file_name,"r")) == NULL)
			error ("Write: Can't open file '%s'",file_name);

		if ((rd_icco = new_icc()) == NULL)
			error ("Write: Creation of ICC object failed");

		/* Read the header and tag list */
		if ((rv = rd_icco->read(rd_icco,rd_fp,0)) != 0)
			error ("Read: %d, %s",rv,rd_icco->err);

		/* ~~ should use an xluobj with merge output ~~~ */
		/* Get the A2B table */
		if ((luo = rd_icco->get_luobj(rd_icco, icmFwd,
                           icAbsoluteColorimetric, icSigLabData, icmLuOrdNorm)) == NULL) {
			error ("%d, %s",rd_icco->errc, rd_icco->err);
		}

		for (i = 0; i < npat; i++) {
			double out[3], ref[3];
			double mxd;

			if (luo->lookup(luo, out, tpat[i].p) > 1)
				error ("%d, %s",rd_icco->errc,rd_icco->err);
		
			/* Our tpat data might be in XYZ, so generate an Lab ref value */
			if (!wantLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
				icmXYZ2Lab(&icmD50, ref, tpat[i].v);

			} else {
				ref[0] = tpat[i].v[0];
				ref[1] = tpat[i].v[1];
				ref[2] = tpat[i].v[2];
			}

			if (verb && verify) {
				printf("[%f] %f %f %f -> %f %f %f should be %f %f %f\n",
				       icmLabDE(ref, out),
				       tpat[i].p[0],tpat[i].p[1],tpat[i].p[2],
				       out[0],out[1],out[2],
				       ref[0],ref[1],ref[2]);
			}

			/* Check the result */
			mxd = icmLabDE(ref, out);
			if (mxd > merr)
				merr = mxd;

			aerr += mxd;
			nsamps++;
		}
		printf("Profile check complete, peak err = %f, avg err = %f\n",merr,aerr/nsamps);

		/* Done with lookup object */
		luo->del(luo);

		/* Close the file */
		rd_icco->del(rd_icco);
		rd_fp->del(rd_fp);
	}

	free(tpat);
}


