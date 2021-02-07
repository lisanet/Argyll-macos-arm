
/* 
 * Argyll Color Correction System
 * Test target chart Generator.
 *
 * Author: Graeme W. Gill
 * Date:   28/9/96
 *
 * Copyright 1996 - 2004, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This program generates a CGATS.5 compatibe file, that */
/* containing device color test patch values. */

/* TTBD:

	Include command line aruments in resulting .ti1 file.

	Should add an option to generate other PCS based pattern
	test points based on the previous profile. (i.e. near neutral
	clusters ?)
	How about an option to read in an CGATS file containing
	PCS or device values ? How is the black level chosen for PCS though ?

	Would be nice to be able to take a previous .ti3 and
	then suppliment the measured patches. Would have to add another
	set of measurement columns to .ti1 & .ti2 to carry the
	already measured values through, or do clumbsy post merge ? 
	(Latter is easiest).

	Would be nice to be able to generate secondary
	color ramps (ie. CMY for RGB space, RGB for CMYK space.)

	Using adaptive patch creation for grey colorspace is broken.
	This should be fixed.

	Would be nice to have a generator of "well behaved" device
	gamut surface points. Would need a way of creating
	a gamut surface from the .ti3 file though.

 */

/* NOTE:

	xpow is applied over the top of the normal or supplied
	device model, so it effectively becomes a modification
	of the device model.

	The general filter is applied prior to the xpow being applied.
	Many of the test patch types do take it into account
	when computing the ink limit.
	The ones that don't are the more complicated full spread patches.

 */

/* Description:

   >> THIS NEEDS REVISION <<

   Nearly all current Color correction systems generate test charts (or
   device characterisation target charts) by laying out a regular rectangular
   grid of test points in device space (Targen will do this if you feed it a non-zero
   m option). On some consideration, this approach is far from optimal. Not only
   is a regular grid inefficent in packing the multidimentional device space,
   but if the points are spaced evenly in device space, they will be poorly
   spaced in human perceptual space, and errors in perceptual space are
   the ultimate arbiter of the end profiles accuracy. Some commercial
   color systems tackle the latter problem by "pre-linearising" the device,
   which amounts to distorting the regular device space grid points with
   a perceptual inverse per device chanel lookup curve.

   The approach I have taken with Argyll, is a little different. By
   using an iterative sphere packing algorithm, I constrain the given
   number of test points to the devices physical gamut (including an
   ink limit for a printer device), and then try and pack the points
   evenly in human perceptual space, or even space them to minimise
   curvature approximation errors. Because the packing is a stocastic
   process, the resulting points are distributed without evident
   patterns.

#ifdef NEVER
   For higher dimensional spaces, where the aim is to create a
   more aproximate device profile, I've used a "perfect simplex
   latice" generator to lay out perfectly packed sample points
   in perceptual space. The latice spacing is sized by an
   iterative search to (hopefully) create the right number of
   test points.
#else
   For higher dimensional spaces, where the aim is to create a
   more aproximate device profile, I've used an "incremental far
   point" point generator, that for each added point, locates
   the device values that result in a percetual value farthest
   from any existing points in the test set.
#endif

   Another issue with laying test points out in regular grids, is
   that this means that the device response is poorly sampled
   (since the grids are usually coarse), and this can make it
   impossible to create detailed device linearisation "shaper"
   curves from the resulting data !
   Ideally, in any colorspace (input or output), when viewed from
   any possible angle, none of the test data points should appear
   to line up. The Argyll target generator seems to acheive this goal.

   A final problem withe regular grids is that they can lead to
   Runge's phenomenon, depending on the nature of the interpolation
   algorithm used.

 */

#undef DEBUG

#define VRML_DIAG		/* Enable option to dump a VRML/X3D of the resulting full spread points */
#undef ADDRECCLIPPOINTS	/* Add ink limited clipping points to regular grid */
#define EMPH_NEUTRAL	/* Emphasise neutral axis, like CIE94 does */
#define NEMPH_DEFAULT 0.5	/* Default neutral axis emphasis == 2 x CIE94 */
#define XPOW_DEFAULT 1.0	/* Default extra device power value = none */
#define DEMPH_DEFAULT 1.0	/* Default dark region emphasis == none */
#define DEFANGLE 0.3333	/* For simdlat and simplat */
#define SIMDLAT_TYPE SIMDLAT_BCC	/* Simdlat geometry type */
#define MATCH_TOLL 1e-4	/* Tollerance of device value to consider a patch a duplicate */

/* Display rise and fall time delay model. This is CRT like */
#define DISPLAY_RISE_TIME DISPTECH_WORST_RISE		/* Assumed rise time to 90% of target level */ 
#define DISPLAY_FALL_TIME DISPTECH_WORST_FALL		/* Assumed fall time to 90% of target level */
#define DISPLAY_SETTLE_AIM 0.1		/* Aim for 0.2 Delta E */

#ifdef NEVER	/* Old delay time code */
#define DISPLAY_SETTLE_AIM2 0.01		/* Aim for 1% of true level */
#define DISPLAY_ABS_AIM2 0.0001		/* Aim for .01% of true absolute level */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "vrml.h"
#include "rspl.h"
#include "cgats.h"
#include "icc.h"
#include "xicc.h"
#include "disptechs.h"
#include "targen.h"
//#include "ppoint.h"
#include "ofps.h"
#include "ifarp.h"
#include "simplat.h"
#include "simdlat.h"
#include "prand.h"
#include "ui.h"

#include <stdarg.h>

#define min2(a,b) ((a) < (b) ? (a) : (b))
#define min3(a,b,c)  (min2((a), min2((b),(c))))
#define max2(a,b) ((a) > (b) ? (a) : (b))
#define max3(a,b,c)  (max2((a), max2((b),(c))))

/* 32 bit pseudo random sequencer */
#define PSRAND32(S) (((S) & 0x80000000) ? (((S) << 1) ^ 0xa398655d) : ((S) << 1))

/* ---------------------------- */
/* The perception function data */
/* (Used for test point distribution) */
struct _pcpt {
/* public: */
	void (*del)(struct _pcpt *s);	/* We're done with it */

	int (*is_specific)(struct _pcpt *s);	/* Is a specific model, not defaulte */

	/* Conversions */
	void (*dev_to_perc)(struct _pcpt *s, double *out, double *in);	/* N-chan Perceptual */
	void (*dev_to_XYZ)(struct _pcpt *s, double *out, double *in);	/* Absolute XYZ */
	void (*dev_to_rLab)(struct _pcpt *s, double *out, double *in);	/* Relative Lab */
	void (*den_to_dev)(struct _pcpt *s, double *out, double *in);	/* Density to device */
	void (*rLab_to_dev)(struct _pcpt *s, double *out, double *in);	/* Lab to device */

	/* !!! Should add perc_to_dev using code from prand that uses dnsq !!! */

/* private: */
	inkmask xmask;		/* external xcolorants inkmask */
	inkmask nmask;		/* internal xcolorants inkmask */
	int di;				/* Number of Device dimensions */

	/* Tuning parameters */
	double nemph;		/* neutral emphasis, 0.0 - 1.0. Default 0.35 for == CIE94 */
	double idemph;		/* inv. of dark emphasis, 1.0 - 4.0. Default 1.0 == none */
	double ixpow;		/* inv. extra power Default 1.0 == none */

	/* ICC profile based */
	icmFile *fp;
	icc   *icco;
	icmLuBase *luo;		/* Device -> rLab conversion */
	icmLuBase *luo2;	/* Device -> XYZ conversion */

	/* MPP profile based */
	mpp *mlu;			/* Device -> XYZ */

	/* Xcolorants model based */
	icxColorantLu *clu;	/* Device -> CIE */

	rspl *nlin[MXTD - 3];	/* Perceptual linearisation for other chanels */
	int e;					/* Chanel being set */

	/* Reverse lookup support */
	double ilimit;			/* Ink limit (scale 1.0) */
	double den[3];			/* Target density or Lab */
	double uniform;			/* NZ if target is uniform */	
	int kchan;				/* Set to the K chanel (-1 if none) */

}; typedef struct _pcpt pcpt;


/* Absolute XYZ conversion function */
/* Internal device values 0.0 - 1.0 are converted into XYZ values */
/* (Used for downstream checking) */
static void
pcpt_to_XYZ(pcpt *s, double *out, double *in) {
	int e; 
	double inv[MXTD];

	if (s->xmask == s->nmask) {
		for (e = 0; e < s->di; e++)
			inv[e] = icx_powlike(in[e], s->ixpow);
	} else {
		for (e = 0; e < s->di; e++)
			inv[e] = 1.0 - icx_powlike(in[e], s->ixpow);
	}
	if (s->luo2 != NULL)
		s->luo2->lookup(s->luo2, out, inv);
	else if (s->mlu != NULL)
		s->mlu->lookup(s->mlu, out, inv);
	else  if (s->clu != NULL)
		s->clu->dev_to_XYZ(s->clu, out, inv);
	else {	/* Linear conversion */
		out[0] = 100.0 * inv[0];
		out[1] = 100.0 * inv[1] - 50.0;
		out[2] = 100.0 * inv[2] - 50.0;
		icmLab2XYZ(&icmD50, out, out);
	}
}


/* Relative Lab conversion function */
/* Internal device values 0.0 - 1.0 are converted into Lab values */
/* (Used for VRML/X3D visualisation checking) */
static void
pcpt_to_rLab(pcpt *s, double *out, double *in) {
	int e; 
	double inv[MXTD];

	if (s->xmask == s->nmask) {
		for (e = 0; e < s->di; e++)
			inv[e] = icx_powlike(in[e], s->ixpow);
	} else {
		for (e = 0; e < s->di; e++)
			inv[e] = 1.0 - icx_powlike(in[e], s->ixpow);
	}
	if (s->luo != NULL)
		s->luo->lookup(s->luo, out, inv);
	else if (s->mlu != NULL) {
		s->mlu->lookup(s->mlu, out, inv);
		icmXYZ2Lab(&icmD50, out, out);
	} else if (s->clu != NULL) 
		s->clu->dev_to_rLab(s->clu, out, inv);
	else {	/* Linear conversion */
		out[0] = 100.0 * inv[0];
		out[1] = 100.0 * inv[1] - 50.0;
		out[2] = 100.0 * inv[2] - 50.0;
	}
}

/* Perceptual conversion function */
/* Internal device values 0.0 - 1.0 are converted into perceptually uniform 0.0 - 100.0 */
/* This is used by optimal spread functions ? */
static void
pcpt_to_nLab(pcpt *s, double *out, double *in) {
	int e;
	double inv[MXTD];

	if (s->xmask == s->nmask) {
		for (e = 0; e < s->di; e++)
			inv[e] = icx_powlike(in[e], s->ixpow);
	} else {
		for (e = 0; e < s->di; e++)
			inv[e] = 1.0 - icx_powlike(in[e], s->ixpow);
	}

	/* If we have some sort of perceptual conversion */
	if (s->luo != NULL || s->mlu != NULL || s->clu != NULL) {
		double lab[3];

		if (s->luo != NULL)
			s->luo->lookup(s->luo, lab, inv);
		else if (s->mlu != NULL) {
			s->mlu->lookup(s->mlu, lab, inv);
			icmXYZ2Lab(&icmD50, lab, lab);
		} else { 
			s->clu->dev_to_rLab(s->clu, lab, inv);
		}

#ifdef EMPH_NEUTRAL		/* Emphasise neutral axis, like CIE94 does */
		{
			double c;		/* Chromanance */

			c = sqrt(lab[1] * lab[1] + lab[2] * lab[2]);	/* Compute chromanance */

//			c = 2.6624 / (1.0 + 0.013 * c);		/* Full strength scale factor */
			c = 3.0 / (1.0 + 0.03 * c);			/* Full strength scale factor */
			c = 1.0 + s->nemph * (c - 1.0);		/* Reduced strength scale factor */

			lab[1] *= c;			/* scale a & b */
			lab[2] *= c;
		}
#endif

		/* Dark emphasis */
		/* This doesn't actually match how demph is applied to device values... */
		if (s->idemph < 1.0) {
			double vv = lab[0];
			lab[0] = 100.0 * pow(lab[0]/100.0, s->idemph);
		}

		/* Copy Lab values to output */
		for (e = 0; e < (s->di < 3 ? s->di : 3); e++)
			out[e] = lab[e];

		/* Lookup perceptual linearised auxiliary values */
		for (e = 0; e < (s->di-3); e++) {
			co cc;
			cc.p[0] = inv[3 + e];
			s->nlin[e]->interp(s->nlin[e], &cc);
			out[3 + e] = cc.v[0];
		}

	} else {
		/* Default linear in Device space */

		for (e = 0; e < s->di; e++)
			out[e] = 100.0 * inv[e];
			if (e == 1 || e == 2)
				out[e] -= 50.0;		/* Make it Lab like */
	}
}


/* Return the largest distance of the point outside the device gamut. */
/* This will be 0 if inside the gamut, and > 0 if outside.  */
static double
pcpt_in_dev_gamut(pcpt *s, double *d) {
	int e;
	int di = s->di;
	double tt, dd = 0.0;
	double ss = 0.0;
	double id[MXTD];

	if (s->xmask == s->nmask) {
		for (e = 0; e < s->di; e++)
			id[e] = d[e];
	} else {
		for (e = 0; e < s->di; e++)
			id[e] = 1.0 - d[e];
	}

	for (e = 0; e < di; e++) {
		ss += id[e];

		tt = 0.0 - id[e];
		if (tt > 0.0) {
			if (tt > dd)
				dd = tt;
		}
		tt = id[e] - 1.0; 
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

/* optimisation function to find device values */
/* for a target density value. */
static double efunc(void *edata, double p[]) {
	pcpt *s = (pcpt *)edata;
	int e, di = s->di;
	double rv, xyz[3], den[4];

//printf("~1 efunc got dev %f %f %f %f\n",p[0],p[1],p[2],p[3]);
	pcpt_to_XYZ(s, xyz, p);		/* Convert device to XYZ */
//printf("~1 efunc got XYZ %f %f %f\n",xyz[0],xyz[1],xyz[2]);
	icx_XYZ2Tdens(den, xyz);	/* Convert XYZ to approx statusT density */
//printf("~1 efunc got density %f %f %f %f\n",den[0],den[1],den[2],den[3]);

//printf("~1 efunc got in_dev_gamut %f\n",pcpt_in_dev_gamut(s, p));

	/* Penalise for being out of gamut */
	rv = 5000.0 * pcpt_in_dev_gamut(s, p);

	/* Error to density target */
	{
		double ss = 0.0;
		for (e = 0; e < 3; e++) {
			double tt;
			tt = s->den[e] - den[e];
			ss += tt * tt;
		}
		rv += ss;
//printf("~1 efunc target den %f %f %f, err = %f, toterr %f\n",s->den[0],s->den[1],s->den[2],ss, rv);
	}

	{
		/* Minimise all channels beyond the */
		/* (assumed) first primary 3, but don't count black. */
		/* Minimise all channels except black if nchan >= 4 and uniform target */
		double ss = 0.0;

		for (e = 0; e < di; e++) {
			double tt = 0.0;

			if (di >= 4 && s->uniform && e < 3 && e != s->kchan)
				tt = p[e];			/* Minimise primary non-black if uniform */
//			else if (!s->uniform && (e < 3 || e == s->kchan))
//				tt = p[e];			/* Minimise sum of primaries & black if uniform */
			else if (e >= 3 && e != s->kchan)
				tt = 3.0 * p[e]; 	/* Suppress non-primary, non-black */

			ss += tt;
		}
		rv += 1.5 * ss * ss;
//printf("~1 efunc sum err = %f, toterr %f\n",ss, rv);
	}

//printf("~1 returning %f\n",rv);
	return rv;
}

/* Given target CMY densities, return a suitable device value */ 
static void
pcpt_den_to_dev(pcpt *s, double *out, double *in) {
	int e, di = s->di;
	double tt, sr[MXTD];	/* Search radius */

//printf("\n");
//printf("~1 targen density = %f %f %f\n",in[0],in[1],in[2]);
//printf("~1 di = %d, ilimit = %f\n",s->di,s->ilimit);

	for (e = 0; e < 3; e++)
		s->den[e] = in[e];

	for (e = 0; e < di; e++) {
		sr[e] = 0.5;			/* Device space search radius */
		out[e] = 0.5;
	}

	if (fabs(in[0] - in[1]) < 0.1
	 && fabs(in[0] - in[2]) < 0.1
	 && fabs(in[1] - in[2]) < 0.1) {
		s->uniform = 1;
//printf("~1 uniform set\n");
	} else
		s->uniform = 0;

	if (powell(&tt, di, out, sr,  0.0001, 2000, efunc, (void *)s, NULL, NULL) != 0 || tt >= 50000.0) {
		error("targen: powell failed, tt = %f\n",tt);
	}

	/* Filter out silly values */
	for (e = 0; e < di; e++) {
		if (out[e] < 0.001)
			out[e] = 0.0;
		else if (out[e] > 0.999)
			out[e] = 1.0;
	}
//printf("~1 returning device values %f %f %f\n",out[0],out[1],out[2]);
}

/* Optimisation function to find device values */
/* for a target Lab value. */
static double efunc2(void *edata, double p[]) {
	pcpt *s = (pcpt *)edata;
	int e, di = s->di;
	double rv, lab[3];

//printf("~1 efunc2 got dev %f %f %f %f\n",p[0],p[1],p[2],p[3]);
//printf("~1 efunc2 got dev %f %f %f %f %f %f\n",p[0],p[1],p[2],p[3],p[4],p[5]);

	pcpt_to_rLab(s, lab, p);		/* Convert device to rLab */
//printf("~1 efunc2 got Lab %f %f %f\n",lab[0],lab[1],lab[2]);

//printf("~1 efunc2 got in_dev_gamut %f\n",pcpt_in_dev_gamut(s, p));

	/* Penalise for being out of gamut */
	rv = 5000.0 * pcpt_in_dev_gamut(s, p);

	/* Error to Lab target */
	{
		double ss = 0.0;
		for (e = 0; e < 3; e++) {
			double tt;
			tt = s->den[e] - lab[e];
			ss += tt * tt;
		}
		rv += ss;
//printf("~1 efunc2 target Lab %f %f %f, err = %f, toterr %f\n",s->den[0],s->den[1],s->den[2],ss, rv);
	}

	{
		int f;

		/* Minimise all channels except K, and especially any */
		/* beyond the first primary 3 or 4. */
		double ss = 0.0;

		if ((s->nmask & ICX_CMYK) == ICX_CMYK)
			f = 4;
		else 
			f = 3;
		for (e = 0; e < di; e++) {
			if (e >= f)
				ss += 10.0 * p[e]; 	/* Suppress non-primary */
			else if (e < 3)
				ss += 0.05 * p[e]; 	/* Suppress first 3 primary slightly */
		}
		rv += ss * ss;
//printf("~1 efunc2 sum err = %f, toterr %f\n",ss, rv);
	}

//printf("~1 efunc2 returning %f\n\n",rv);
	return rv;
}

/* Given target Lab densities, return a suitable device value */ 
static void
pcpt_rLab_to_dev(pcpt *s, double *out, double *in) {
	int e, di = s->di;
	double tt, sr[MXTD];	/* Search radius */

//printf("\n");
//printf("#######################3\n");
//printf("~1 targen Lab = %f %f %f\n",in[0],in[1],in[2]);
//printf("~1 di = %d, ilimit = %f\n",s->di,s->ilimit);

	for (e = 0; e < 3; e++)
		s->den[e] = in[e];

	for (e = 0; e < di; e++) {
		sr[e] = 0.5;			/* Device space search radius */
		out[e] = 0.5;
	}

	if (powell(&tt, di, out, sr,  0.0001, 2000, efunc2, (void *)s, NULL, NULL) != 0 || tt >= 50000.0) {
		error("targen: powell failed, tt = %f\n",tt);
	}

	/* Filter out silly values */
	for (e = 0; e < di; e++) {
		if (out[e] <= 0.005)
			out[e] = 0.0;
		else if (out[e] >= 0.995)
			out[e] = 1.0;
	}
//printf("~1 returning device values %f %f %f\n",out[0],out[1],out[2]);
}

/* Callback to setup s->nlin[e] mapping */
static void set_nlin(void *cbntx, double *out, double *in) {
	pcpt *s = (pcpt *)cbntx;	/* Object we're setting up from */
	int e, di = s->di;
	double dev[MXTD];
	double lab[3];

	/* Just input extra channel into perceptual type lookup */
	if (s->xmask == s->nmask) {
		for (e = 0; e < di; e++) 
			dev[e] = 0.0;
		dev[3 + s->e] = in[0];
	} else {						/* Fake RGB */
		for (e = 0; e < di; e++) 
			dev[e] = 1.0;
		dev[3 + s->e] = 1.0 - in[0];
	}

	if (s->luo != NULL) {
		s->luo->lookup(s->luo, lab, dev);
	} else if (s->mlu != NULL) {
		s->mlu->lookup(s->mlu, lab, dev);
		icmXYZ2Lab(&icmD50, lab, lab);
	} else  if (s->clu != NULL) {
		s->clu->dev_to_rLab(s->clu, lab, dev);
	} else {
		lab[0] = 100.0 * in[0];
	}

	/* ~~~ should we make this delta lab along locus, rather than L value ??? */
	out[0] = lab[0];
}

/* Is a specific model, not default */
int pcpt_is_specific(pcpt *s) {
	if (s->luo2 != NULL || s->mlu != NULL)
		return 1;
	return 0;
}

/* Free the pcpt */
static void pcpt_del(pcpt *s) {

	if (s != NULL) {
		int e;

		if (s->luo != NULL) {
			s->luo->del(s->luo);
			s->luo2->del(s->luo2);
			s->icco->del(s->icco);
			s->fp->del(s->fp);
		}
		if (s->mlu != NULL) {
			s->mlu->del(s->mlu);
		}
		if (s->clu != NULL) {
			s->clu->del(s->clu);
		}
		for (e = 0; e < (s->di-3); e++) {
			if (s->nlin[e] != NULL)
				s->nlin[e]->del(s->nlin[e]);
		}
		
		free(s);
	}
}

/* Create a pcpt conversion class */
pcpt *new_pcpt(
char *profName,			/* ICC or MPP profile path, NULL for default, "none" for linear */
inkmask xmask,			/* external xcolorants mask */
inkmask nmask,			/* internal xcolorants mask */
double *ilimit,			/* ink sum limit (scale 1.0) input and return, -1 if default */
double *uilimit,		/* underlying ink sum limit (scale 1.0) input and return, -1 if default */
double nemph,			/* Neutral emphasis, 0.0 - 1.0. < 0.0 for default == CIE94 */
double demph,			/* Dark emphasis, 1.0 - 4.0. < 0.0 for default == none */
double xpow				/* Extra device power, default = none */
) {
	int e;
	pcpt *s;

	if ((s = (pcpt *)calloc(1, sizeof(pcpt))) == NULL) {
		fprintf(stderr,"targen: malloc failed allocating pcpt object\n");
		exit(-1);
	}
	
	/* Initialise methods */
	s->del          = pcpt_del;
	s->is_specific  = pcpt_is_specific;
	s->dev_to_perc  = pcpt_to_nLab;
	s->dev_to_XYZ   = pcpt_to_XYZ;
	s->dev_to_rLab  = pcpt_to_rLab;
	s->den_to_dev   = pcpt_den_to_dev;
	s->rLab_to_dev  = pcpt_rLab_to_dev;

	s->xmask = xmask;
	s->nmask = nmask;
	s->di = icx_noofinks(nmask);

	if (nemph < 0.0)
		nemph = NEMPH_DEFAULT;
	s->nemph = nemph;

	if (demph < 0.0)
		demph = DEMPH_DEFAULT;
	s->idemph = 1.0/demph;

	if (xpow < 0.0)
		xpow = XPOW_DEFAULT;
	s->ixpow = 1.0/xpow;

	/* See if we have a profile */
	if (profName != NULL
	 && profName[0] != '\000'
	 && strcmp(profName, "none") != 0
	 && strcmp(profName, "NONE") != 0) {
		int rv = 0;
		
		/* Try and open the file as an ICC profile */
		if ((s->fp = new_icmFileStd_name(profName,"r")) == NULL)
			error ("Can't open device profile '%s'",profName);
	
	
		if ((s->icco = new_icc()) == NULL)
			error ("Creation of ICC object failed");
	
		if ((rv = s->icco->read(s->icco,s->fp,0)) == 0) {
			icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
			xcal *cal = NULL;					/* Device calibration curves */

			/* Get a conversion object for relative Lab */
			if ((s->luo = s->icco->get_luobj(s->icco, icmFwd, icRelativeColorimetric,
			                                 icSigLabData, icmLuOrdNorm)) == NULL) {
				if ((s->luo = s->icco->get_luobj(s->icco, icmFwd, icmDefaultIntent,
				                                 icSigLabData, icmLuOrdNorm)) == NULL) {
					error ("%d, %s",s->icco->errc, s->icco->err);
				}
			}

			/* Get a conversion object for absolute XYZ */
			if ((s->luo2 = s->icco->get_luobj(s->icco, icmFwd, icAbsoluteColorimetric,
			                                 icSigXYZData, icmLuOrdNorm)) == NULL) {
				if ((s->luo2 = s->icco->get_luobj(s->icco, icmFwd, icmDefaultIntent,
				                                 icSigXYZData, icmLuOrdNorm)) == NULL) {
					error ("%d, %s",s->icco->errc, s->icco->err);
				}
			}
		
			/* Get details of conversion (Arguments may be NULL if info not needed) */
			s->luo->spaces(s->luo, &ins, NULL, &outs, NULL, NULL, NULL, NULL, NULL, NULL);

//printf("~1 xmask = 0x%x, ins = %s\n",xmask,icm2str(icmColorSpaceSignature, ins));
			if (icx_colorant_comb_match_icc(xmask, ins) == 0) {

				/* Should really see if ICC profile has ColorantTable tag, */
				/* and match them against targen specs. For now, */
				/* simply make sure the channel counts match and issue */
				/* a warning. */
				if (icx_noofinks(xmask) != icmCSSig2nchan(ins)) {
					s->luo->del(s->luo);
					error("ICC profile doesn't match device!");
				} else {
					warning("Profile '%s' no. channels match, but colorant types have not been checked",profName);
				}
			}

			/* Grab any device calibration curves */
			cal = xiccReadCalTag(s->icco);

			/* Set the default ink limits if not set by user */
			if (*ilimit < 0.0) {

				if (cal != NULL) {
					*ilimit = s->icco->get_tac(s->icco, NULL, xiccCalCallback, (void *)cal);
					*uilimit = s->icco->get_tac(s->icco, NULL, NULL, NULL);
				} else {
					*uilimit = *ilimit = s->icco->get_tac(s->icco, NULL, NULL, NULL);
				}
				*ilimit += 0.1;		/* + 10% */
				*uilimit += 0.1;		/* + 10% */

			/* Convert the user limit to a maximum underlying limit */
			} else if (cal != NULL && *ilimit < (double)s->di) {
				*uilimit = icxMaxUnderlyingLimit(cal, *ilimit);
			}

		} else {	/* Not a valid ICC */
			/* Close out the ICC profile */
			s->icco->del(s->icco);
			s->icco = NULL;
			s->fp->del(s->fp);
			s->fp = NULL;
		}

		/* If we don't have an ICC lookup object, look for an MPP */
		if (s->luo == NULL) {
			inkmask imask;
			double dlimit = 0.0;

			if ((s->mlu = new_mpp()) == NULL)
				error ("Creation of MPP object failed");

			if ((rv = s->mlu->read_mpp(s->mlu, profName)) != 0)
				error ("%d, %s",rv,s->mlu->err);

			s->mlu->get_info(s->mlu, &imask, NULL, &dlimit, NULL, NULL, NULL, NULL, NULL);

			if (xmask != imask) {
				s->mlu->del(s->mlu);
				error("MPP profile doesn't match device!");
			}
			if (*ilimit < 0.0 && dlimit > 0.0)	 {/* If not user specified, use MPP inklimit */
				*uilimit = *ilimit = dlimit + 0.1;
			}
		}
	}

	/* Fall back on an xcolorants model */
	if (s->luo == NULL
	 && s->mlu == NULL
	 && strcmp(profName, "none") != 0
	 && strcmp(profName, "NONE") != 0) {
		if ((s->clu = new_icxColorantLu(xmask)) == NULL)
			error ("Creation of xcolorant lu object failed");
	}
	/* else leave pointers NULL */

	if (*ilimit < 0.0)
		s->ilimit = (double)s->di;	/* Default to no limit */
	else
		s->ilimit = *ilimit/100.0;

	if (s->di > 1)
		s->kchan = icx_ink2index(xmask, ICX_BLACK);
	else
		s->kchan = -1;

	/* Create extra chanel linearisation lookups */
	for (e = 0; e < (s->di-3); e++) {
		double inmin = 0.0, inmax = 1.0;
		double outmax = 100.0;
		int gres = 256;

		if ((s->nlin[e] = new_rspl(RSPL_NOFLAGS, 1, 1)) == NULL)
			error("RSPL creation failed");

		s->e = e;	/* Chanel to set */
		s->nlin[e]->set_rspl(s->nlin[e], 0, s, set_nlin,
		                     &inmin, &inmax, &gres, &inmax, &outmax);
	}

	return s;
}

/* ------------------------------------ */

void
usage(int level, char *diag, ...) {
	int i;
	fprintf(stderr,"Generate Target deviceb test chart color values, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: targen [options] outfile\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr," -v [level]      Verbose mode [optional level 1..N]\n");
	fprintf(stderr," -d col_comb     choose colorant combination from the following:\n");
	for (i = 0; ; i++) {
		char *desc; 
		if (icx_enum_colorant_comb(i, &desc) == 0)
			break;
		fprintf(stderr,"                 %d: %s\n",i,desc);
	}
	fprintf(stderr," -D colorant     Add or delete colorant from combination:\n");
	if (level == 0)
		fprintf(stderr,"                 (Use -?? to list known colorants)\n");
	else {
		fprintf(stderr,"                 %d: %s\n",0,"Additive");
		for (i = 0; ; i++) {
			char *desc; 
			if (icx_enum_colorant(i, &desc) == 0)
				break;
			fprintf(stderr,"                 %d: %s\n",i+1,desc);
		}
	}
	fprintf(stderr," -G               Generate good optimized points rather than Fast\n");
	fprintf(stderr," -e patches       White test patches (default 4)\n");
	fprintf(stderr," -B patches       Black test patches (default 4 Grey/RGB, else 0)\n");
	fprintf(stderr," -s steps         Single channel steps (default grey 50, color 0)\n");
	fprintf(stderr," -g steps         Grey axis RGB or CMY steps (default 0)\n");
	fprintf(stderr," -n steps         Neutral axis steps (based on profile, default 0)\n");
	fprintf(stderr," -m steps         Multidimensional device space cube steps (default 0)\n");
	fprintf(stderr," -M steps         Multidimensional device space cube surface steps (default 0)\n");
	fprintf(stderr," -b steps         Multidimensional body centered cubic steps (default 0)\n");
	fprintf(stderr," -f patches       Add iterative & adaptive full spread patches to total (default grey 0, color 836)\n");
	fprintf(stderr,"                  Default is Optimised Farthest Point Sampling (OFPS)\n");
	fprintf(stderr,"  -t              Use incremental far point for full spread\n");
	fprintf(stderr,"  -r              Use device space random for full spread\n");
	fprintf(stderr,"  -R              Use perceptual space random for full spread\n");
	fprintf(stderr,"  -q              Use device space-filling quasi-random for full spread\n");
	fprintf(stderr,"  -Q              Use perceptual space-filling quasi-random for full spread\n");
	fprintf(stderr,"  -i              Use device space body centered cubic grid for full spread\n");
	fprintf(stderr,"  -I              Use perceptual space body centered cubic grid for full spread\n");
	fprintf(stderr,"  -a angle        Simplex grid angle 0.0 - 0.5 for B.C.C. grid, default %f\n",DEFANGLE);
	fprintf(stderr,"  -A adaptation   Degree of adaptation of OFPS 0.0 - 1.0 (default 0.1, -c profile used 1.0)\n");
/* Research options: */
/*	fprintf(stderr,"  -A pPERCWGHT    Device (0.0) ... Perceptual (1.0) weighting\n"); */
/*	fprintf(stderr,"  -A cCURVEWGHT   Curvature weighting  0.0 = none ... "); */
	fprintf(stderr," -l ilimit        Total ink limit in %% (default = none)\n");
	fprintf(stderr," -p power         Optional power-like value applied to all device values.\n");
	fprintf(stderr," -c profile       Optional device ICC or MPP pre-conditioning profile filename\n");
	fprintf(stderr,"                  (Use \"none\" to turn off any conditioning)\n");
	fprintf(stderr," -N nemphasis     Degree of neutral axis patch concentration 0.0-1.0 (default %.2f)\n",NEMPH_DEFAULT);
	fprintf(stderr," -V demphasis     Degree of dark region patch concentration 1.0-4.0 (default %.2f = none)\n",DEMPH_DEFAULT);
	fprintf(stderr," -F L,a,b,rad     Filter out samples outside Lab sphere.\n");
	fprintf(stderr," -O               Don't re-order display RGB patches for minimum delay\n");
	fprintf(stderr," -U               Don't filter out duplicate patches\n");
#ifdef VRML_DIAG
	fprintf(stderr," -w               Dump diagnostic outfilel%s file (Lab locations)\n",vrml_ext());
	fprintf(stderr," -W               Dump diagnostic outfiled%s file (Device locations)\n",vrml_ext());
#endif /* VRML_DIAG */
	fprintf(stderr," outfile          Base name for output(.ti1)\n");
	exit(1);
}

/* Test if outside filter sphere. Return nz if it is */
int dofilt(
	pcpt *pdata,		/* Perceptual conversion routine */
	double *filt,		/* Filter sphere definition */
	double *dev			/* Device values to check */
) {
	int i;
	double Lab[3], rr;
	pdata->dev_to_rLab(pdata, Lab, dev);
	for (rr = 0.0, i = 0; i < 3; i++) {
		double tt = Lab[i] - filt[i];
		rr += tt * tt;
	}
	if (rr > (filt[3] * filt[3])) {
//printf("~1 rejecting rad %f of %f %f %f <=> %f %f %f\n",sqrt(rr),Lab[0],Lab[1],Lab[2],filt[0],filt[1],filt[2]);
		return 1;
	}
	return 0;
}


static double disprespt(cgats *pp, int p1, int p2);

int main(int argc, char *argv[]) {
	int i, j, k;
	int fa, nfa, mfa;		/* current argument we're looking at */
	int verb = 0;			/* Verbose flag */
#ifdef VRML_DIAG
	int dumpvrml = 0;		/* Dump diagnostic VRML/X3D file */
#endif /* VRML_DIAG */
	inkmask xmask = 0;		/* External ink mask combination */
	inkmask nmask = 0;		/* Working ink mask combination (ie. CMY for printer external sRGB) */
	int di = 0;				/* Output dimensions */
	char *ident;			/* Ink combination identifier (includes possible leading 'i') */
	int good = 0;			/* 0 - fast, 1 = good */
	int esteps = -1;		/* White color patches */
	int Bsteps = -1;		/* Black color patches */
	int ssteps = -1;		/* Single channel steps */
	double xpow = 1.0;		/* Power to apply to all device values created */
	int gsteps = -1;		/* Composite grey wedge steps */
	int nsteps = -1;		/* Neutral wedge steps */
	int msteps = 0;			/* Regular grid multidimensional steps */
	int msurf = 0;			/* If nz, make msteps just on device surface */
	int bsteps = 0;			/* Regular body centered cubic grid multidimensional steps */
	int fsteps = -1;		/* Fitted Multidimensional patches */
	int uselat = 0;			/* Use incremental far point alg. for full spread points */
	int userand = 0;		/* Use random for full spread points, 2 = perceptual */
	int useqrand = 0;		/* Use sobol for full spread points, 2 = perceptual */
	int usedsim = 0;		/* Use device space simplex grid */
	int usepsim = 0;		/* Use perceptual space simplex grid */
	double simangle = DEFANGLE;	/* BCC grid angle */
	double dadapt = -2.0;	/* Degree of iterative adaptation */
	double perc_wght = 0.0;	/* Perceptual weighting */
	double curv_wght = 0.0;	/* Curvature weighting */
	double ilimit = -1.0;	/* Ink limit (scale 1.0) (default none) */
	double uilimit = -1.0;	/* Underlying (pre-calibration, scale 1.0) ink limit */
	double nemph = NEMPH_DEFAULT;
	double demph = DEMPH_DEFAULT;
	int dontdedupe = 0;		/* Don't filter duplicate samples */
	int dontreorder = 0;	/* Don't re-order RGB display patches for min delay */
	int filter = 0;			/* Filter values outside given sphere */
	double filt[4] = { 50,0,0,0 };	
	static char fname[MAXNAMEL+1] = { 0 };		/* Output file base name */
	static char pname[MAXNAMEL+1] = { 0 };		/* Device profile name */
	static char wdname[MAXNAMEL+1] = { 0 };		/* Device diagnostic .wrl/.x3d name */
	static char wlname[MAXNAMEL+1] = { 0 };		/* Lab diagnostic .wrl/.x3d name */
	char buf[500];			/* Genaral use text buffer */
	int id = 1;				/* Sample ID */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	cgats *pp;				/* cgats structure */
	long stime,ttime;
	pcpt *pdata;			/* Space linearisation callback struct */
	fxpos *fxlist = NULL;	/* Fixed point list for full spread */
	int fxlist_a = 0;		/* Fixed point list allocation */
	int fxno = 0;			/* The number of fixed points */

#ifdef NUMSUP_H
	error_program = "targen";
#endif
	check_if_not_interactive();

	if (argc <= 1)
		usage(0,"Too few arguments, got %d expect at least %d",argc-1,1);

	/* Process the arguments */
	mfa = 1;        /* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1+mfa) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?' || argv[fa][1] == '-') {
				if (argv[fa][2] == '?' || argv[fa][2] == '-')
					usage(1, "Extended usage requested");
				usage(0, "Usage requested");
			}

			else if (argv[fa][1] == 'v') {
				verb = 1;
				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					verb = atoi(na);
					fa = nfa;
				}
			}

			/* Select the ink enumeration */
			else if (argv[fa][1] == 'd') {
				fa = nfa;
				if (na == NULL) usage(0,"Expect argument after -d");
				i = atoi(na);
				if (i == 0 && na[0] != '0')
					usage(0,"Expect number argument after -d");
				if ((xmask = icx_enum_colorant_comb(i, NULL)) == 0)
					usage(0,"Argument to -d is not recognized");
			}
			/* Toggle the colorant in ink combination */
			else if (argv[fa][1] == 'D') {
				int tmask;
				fa = nfa;
				if (na == NULL) usage(0,"Expect argument after -D");
				i = atoi(na);
				if (i == 0 && na[0] != '0')
					usage(0,"Expect number argument after -D");
				if (i == 0)
					tmask = ICX_ADDITIVE;
				else
					if ((tmask = icx_enum_colorant(i-1, NULL)) == 0)
						usage(0,"Argument to -D is not recognized");
				xmask ^= tmask;
			}
			/* Good rather than fast */
			else if (argv[fa][1] == 'G') {
				good = 1;
			}
			/* White color patches */
			else if (argv[fa][1] == 'e') {
				int tt;
				if (na == NULL) usage(0,"Expect argument after -e");
				if ((tt = atoi(na)) >= 0)
					esteps = tt;
				fa = nfa;
			}
			/* Black color patches */
			else if (argv[fa][1] == 'B') {
				int tt;
				if (na == NULL) usage(0,"Expect argument after -B");
				if ((tt = atoi(na)) >= 0)
					Bsteps = tt;
				fa = nfa;
			}
			/* Individual chanel steps */
			else if (argv[fa][1] == 's') {
				int tt;
				if (na == NULL) usage(0,"Expect argument after -s");
				if ((tt = atoi(na)) >= 0)
					ssteps = tt;
				fa = nfa;
			}
			/* PCS based neutral wedge steps */
			else if (argv[fa][1] == 'n') {
				int tt;
				if (na == NULL) usage(0,"Expect argument after -n");
				if ((tt = atoi(na)) >= 0)
					nsteps = tt;
				fa = nfa;
			}
			/* RGB or CMY grey wedge steps */
			else if (argv[fa][1] == 'g') {
				int tt;
				if (na == NULL) usage(0,"Expect argument after -g");
				if ((tt = atoi(na)) >= 0)
					gsteps = tt;
				fa = nfa;
			}
			/* Multidimentional cube steps */
			else if (argv[fa][1] == 'm'
			      || argv[fa][1] == 'M') {
				int tt;
				if (na == NULL) usage(0,"Expect argument after -m");
				if ((tt = atoi(na)) >= 0) {
					msteps = tt;
					if (msteps == 1)
						msteps = 2;
				}
				if (argv[fa][1] == 'M')
					msurf = 1;
				fa = nfa;
			}
			/* Multidimentional body centered cube steps */
			else if (argv[fa][1] == 'b') {
				int tt;
				if (na == NULL) usage(0,"Expect argument after -b");
				if ((tt = atoi(na)) >= 0) {
					bsteps = tt;
					if (bsteps == 1)
						bsteps = 2;
				}
				fa = nfa;
			}
			/* Full even spread Multidimentional patches */
			else if (argv[fa][1] == 'f') {
				int tt;
				if (na == NULL) usage(0,"Expect argument after -f");
				if ((tt = atoi(na)) >= 0)
					fsteps = tt;
				fa = nfa;
			}

			/* Use incremental far point algorithm for full spread */
			else if (argv[fa][1] == 't') {
				uselat = 1;
				userand = 0;
				useqrand = 0;
				usedsim = 0;
				usepsim = 0;
			}

			/* Random requested */
			else if (argv[fa][1] == 'r'
				  || argv[fa][1] == 'R') {
				uselat = 0;
				if (argv[fa][1] == 'R')
					userand = 2;
				else
					userand = 1;
				useqrand = 0;
				usedsim = 0;
				usepsim = 0;
			}

			/* Space filling quasi-random requested */
			else if (argv[fa][1] == 'q'
				  || argv[fa][1] == 'Q') {
				uselat = 0;
				userand = 0;
				if (argv[fa][1] == 'Q')
					useqrand = 2;
				else
					useqrand = 1;
				usedsim = 0;
				usepsim = 0;
			}


			/* Device simplex grid requested */
			else if (argv[fa][1] == 'i') {
				uselat = 0;
				userand = 0;
				useqrand = 0;
				usedsim = 1;
				usepsim = 0;
			}

			/* Perceptual simplex grid requested */
			else if (argv[fa][1] == 'I') {
				uselat = 0;
				userand = 0;
				useqrand = 0;
				usedsim = 0;
				usepsim = 1;
			}

			/* Simplex grid angle */
			else if (argv[fa][1] == 'a') {
				if (na == NULL) usage(0,"Expect argument after -a");
				simangle = atof(na);
				fa = nfa;
			}

			/* Degree of iterative adaptation */
			else if (argv[fa][1] == 'A') {
				if (na == NULL) usage(0,"Expected argument to average deviation flag -A");
				if (na[0] == 'p') {			/* (relative, for verification) */
					perc_wght = atof(na+1);
					if (perc_wght < 0.0 || perc_wght > 1.0)
						usage(0,"Perceptual weighting argument %f to '-Ap' must be between 0.0 and 1.0",perc_wght);
					dadapt = -1.0;
				} else if (na[0] == 'c') {	/* (absolute, for testing) */
					curv_wght = atof(na+1);
					if (curv_wght < 0.0 || curv_wght > 100.0)
						usage(0,"Curvature weighting argument %f to '-Ac' must be between 0.0 and 100.0",curv_wght);
					dadapt = -1.0;
				} else {
					dadapt = atof(na);
					if (dadapt < 0.0 || dadapt > 1.0)
						usage(0,"Average Deviation argument %f must be between 0.0 and 1.0",dadapt);
				}
				fa = nfa;
			}

			/* Ink limit percentage */
			else if (argv[fa][1] == 'l') {
				double tt;
				if (na == NULL) usage(0,"Expect argument after -l");
				if ((tt = atof(na)) > 0.0)
					uilimit = ilimit = 0.01 * tt;
				fa = nfa;
			}

			/* Extra device power-like to use */
			else if (argv[fa][1] == 'p') {
				double tt;
				if (na == NULL) usage(0,"Expect argument after -p");
				if ((tt = atof(na)) > 0.0)
					xpow = tt;
				fa = nfa;
			}

			/* ICC profile for perceptual linearisation */
			else if (argv[fa][1] == 'c') {
				if (na == NULL) usage(0,"Expect argument after -c");
				strncpy(pname,na,MAXNAMEL-1); pname[MAXNAMEL-1] = '\000';
				fa = nfa;
			}

			/* Degree of neutral axis emphasis */
			else if (argv[fa][1] == 'N') {
				if (na == NULL) usage(0,"Expected argument to neutral emphasis flag -N");
				nemph = atof(na);
				if (nemph < 0.0 || nemph > 10.0)
					usage(0,"Neautral weighting argument %f to '-N' is out of range",nemph);
				fa = nfa;
			}

			/* Degree of dark region emphasis */
			else if (argv[fa][1] == 'V') {
				if (na == NULL) usage(0,"Expected argument to dark emphasis flag -V");
				demph = atof(na);
				if (demph < 1.0 || demph > 4.0)
					usage(0,"Dark weighting argument %f to '-V' is out of range",demph);
				fa = nfa;
			}

			/* Filter out samples outside given sphere */
			else if (argv[fa][1] == 'F') {
				if (na == NULL) usage(0,"Expect argument after -F");
				if (sscanf(na, " %lf,%lf,%lf,%lf ",&filt[0], &filt[1], &filt[2], &filt[3]) != 4)
					usage(0,"Argument to -F '%s' isn't correct",na);
				filter = 1;
				fa = nfa;
			}

			/* Don't re-order RGB display patches for best speed */
			else if (argv[fa][1] == 'O') {
				dontreorder = 1;
			}

			/* Don't filter out redundant patches */
			else if (argv[fa][1] == 'U') {
				dontdedupe = 1;
			}

#ifdef VRML_DIAG
			else if (argv[fa][1] == 'w')		/* Lab */
				dumpvrml |= 1;
			else if (argv[fa][1] == 'W')		/* Device */
				dumpvrml |= 2;
#endif /* VRML_DIAG */
			else 
				usage(0,"Unknown flag '%c'",argv[fa][1]);
		}
		else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage(0,"Expect file base name argument");
	strncpy(fname,argv[fa],MAXNAMEL-4); fname[MAXNAMEL-4] = '\000';
	strcat(fname,".ti1");

	strncpy(wdname,argv[fa],MAXNAMEL-5); wdname[MAXNAMEL-5] = '\000';
	strcat(wdname,"d");

	strncpy(wlname,argv[fa],MAXNAMEL-5); wlname[MAXNAMEL-5] = '\000';
	strcat(wlname,"l");

	/* Set default colorant combination as CMYK */
	if (xmask == 0)
		xmask = ICX_CMYK;
	
	nmask = xmask;

	/* Deal with fake printer RGB, where we use CMY internally and invert all */
	/* the resulting device values. */
	if (xmask & ICX_INVERTED) {
		if (xmask != ICX_IRGB)
			error("Don't know how to deal with inverted colorant combination 0x%x\n",xmask);
		nmask = ICX_CMY;	/* Internally treat it as CMY and invert the result */
	}

	ident = icx_inkmask2char(xmask, 1); 
	di = icx_noofinks(nmask);	/* Lookup number of dimensions */
	stime = clock();

	/* Implement some defaults */
	if (esteps < 0)
		esteps = 4;
	if (gsteps < 0)
		gsteps = 0;
	if (nsteps < 0)
		nsteps = 0;

	if (Bsteps < 0) {
		if (xmask == ICX_W || xmask == ICX_K || xmask == ICX_RGB || xmask == ICX_IRGB)
			Bsteps = 4;
		else
			Bsteps = 0; 
	}

	if (di == 1) {
		if (ssteps < 0)
			ssteps = 50;
		if (fsteps < 0)
			fsteps = 0;
	} else {
		if (ssteps < 0)		/* Defaults */
			ssteps = 0;
		if (fsteps < 0)
			fsteps = 836;
	}

	/* Do some sanity checking */
	if (di == 1) {
		if (ssteps == 0 && fsteps == 0 && msteps == 0 && bsteps == 0)
			error ("Must have some Gray steps");
		if (gsteps > 0) {
			warning ("Composite grey steps ignored for monochrome output");
			gsteps = 0;
		}
		if (nsteps > 0) {
			warning ("Neutral steps ignored for monochrome output");
			nsteps = 0;
		}
	} else if (di == 3) {
			if (ssteps == 0 && fsteps == 0 && msteps == 0 && bsteps == 0
			 && gsteps == 0 && nsteps == 0)
				error ("Must have some single or multi dimensional RGB or CMY steps");
		} else {
			if (ssteps == 0 && fsteps == 0 && msteps == 0 && bsteps == 0
			    && gsteps == 0 && nsteps == 0)
				error ("Must have some single or multi dimensional steps");
		}

		/* Deal with ICC, MPP or fallback profile */
		if ((pdata = new_pcpt(pname, xmask, nmask, &ilimit, &uilimit, nemph, demph, xpow)) == NULL) {
			error("Perceptual lookup object creation failed");
		}

		/* Set default adapation level */
		if (dadapt < -1.5) {		/* Not set by user */
			if (pname[0] != '\000')
				dadapt = 1.0;
			else
				dadapt = 0.1;
		}

		if (verb) {
			printf("%s test chart\n",ident);

			if (esteps > 0)
				printf("White patches = %d\n",esteps);
			if (Bsteps > 0)
				printf("Black patches = %d\n",Bsteps);
			if (ssteps > 0)
				printf("Single channel steps = %d\n",ssteps);
			if (gsteps > 0)
				printf("Composite Grey steps = %d\n",gsteps);
			if (nsteps > 0)
				printf("Neutral steps = %d\n",nsteps);
			if (fsteps > 0)
				printf("Full spread patches = %d\n",fsteps);
			if (msteps > 0) {
				if (msurf)
					printf("Multi-dimention cube surface steps = %d\n",msteps);
				else
					printf("Multi-dimention cube steps = %d\n",msteps);
			}
			if (bsteps > 0)
				printf("Multi-dimention body centered cube steps = %d\n",bsteps);
			if (ilimit >= 0.0)
				printf("Ink limit = %.1f%% (underlying %.1f%%)\n",ilimit * 100.0, uilimit * 100.0);
			if (filter) {
				printf("Filtering out samples outside sphere at %f %f %f radius %f\n",
						filt[0], filt[1], filt[2], filt[3]);
			}
		}
		pp = new_cgats();	/* Create a CGATS structure */
		pp->add_other(pp, "CTI1"); 	/* our special type is Calibration Target Information 1 */

		pp->add_table(pp, tt_other, 0);	/* Add the first table for target points */
		pp->add_table(pp, tt_other, 0);	/* Add the second table for density pre-defined device values */
		pp->add_table(pp, tt_other, 0);	/* Add the second table for device pre-defined device values */
		pp->add_kword(pp, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 1",NULL);
		pp->add_kword(pp, 1, "DESCRIPTOR", "Argyll Calibration Target chart information 1",NULL);
		pp->add_kword(pp, 2, "DESCRIPTOR", "Argyll Calibration Target chart information 1",NULL);
		pp->add_kword(pp, 0, "ORIGINATOR", "Argyll targen", NULL);
		pp->add_kword(pp, 1, "ORIGINATOR", "Argyll targen", NULL);
		pp->add_kword(pp, 2, "ORIGINATOR", "Argyll targen", NULL);
		atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
		pp->add_kword(pp, 0, "CREATED",atm, NULL);

		/* Make available the aproximate white point to allow relative */
		/* interpretation of the aproximate XYZ values */
		{
			int e;
			double val[MXTD], XYZ[3];

			/* Setup device white */
			if (nmask & ICX_ADDITIVE)
				for (e = 0; e < di; e++)
					val[e] = 1.0;
			else
				for (e = 0; e < di; e++)
					val[e] = 0.0;
			pdata->dev_to_XYZ(pdata, XYZ, val);		/* Lookup white XYZ */

			sprintf(buf,"%f %f %f", 100.0  * XYZ[0], 100.0 * XYZ[1], 100.0 * XYZ[2]);
			pp->add_kword(pp, 0, "APPROX_WHITE_POINT", buf, NULL);
		}

		pp->add_field(pp, 0, "SAMPLE_ID", cs_t);
		pp->add_field(pp, 1, "INDEX", i_t);			/* Index no. 0..7 in second table */
		pp->add_field(pp, 2, "INDEX", i_t);			/* Index no. 0..7 in third table */

		/* Setup CGATS fields */
		{
			int j;
			char c_ilimit[20];
			char *bident = icx_inkmask2char(xmask, 0); 

			for (j = 0; j < di; j++) {
				int imask;
				char fname[100];

				imask = icx_index2ink(xmask, j);
				sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : bident,
									  icx_ink2char(imask));

				pp->add_field(pp, 0, fname, r_t);
				pp->add_field(pp, 1, fname, r_t);
				pp->add_field(pp, 2, fname, r_t);
			}

			pp->add_kword(pp, 0, "COLOR_REP", ident, NULL);

			if (ilimit >= 0.0) {
				sprintf(c_ilimit,"%5.1f",ilimit * 100.0);
				pp->add_kword(pp, 0, "TOTAL_INK_LIMIT", c_ilimit, NULL);
			}
			free(bident);
		}

		/* ilimit is assumed to be in a valid range from here on */
		if (ilimit < 0.0) {
			uilimit = ilimit = di;	/* default is no limit */
		}

		/* Add expected XYZ values to aid previews, scan recognition & strip recognition */
		pp->add_field(pp, 0, "XYZ_X", r_t);
		pp->add_field(pp, 0, "XYZ_Y", r_t);
		pp->add_field(pp, 0, "XYZ_Z", r_t);
		pp->add_field(pp, 1, "XYZ_X", r_t);
		pp->add_field(pp, 1, "XYZ_Y", r_t);
		pp->add_field(pp, 1, "XYZ_Z", r_t);
		pp->add_field(pp, 2, "XYZ_X", r_t);
		pp->add_field(pp, 2, "XYZ_Y", r_t);
		pp->add_field(pp, 2, "XYZ_Z", r_t);

		/* Note if the expected values are expected to be accurate */
		if (pdata->is_specific(pdata))
			pp->add_kword(pp, 0, "ACCURATE_EXPECTED_VALUES", "true", NULL);

		if (xpow != 1.0) {
			sprintf(buf,"%f",xpow);
			pp->add_kword(pp, 0, "EXTRA_DEV_POW",buf, NULL);
		}

		if (demph > 1.0) {
			sprintf(buf,"%f",demph);
			pp->add_kword(pp, 0, "DARK_REGION_EMPHASIS",buf, NULL);
		}

		/* Only use optimsed full spread if <= 4 dimensions, else use ifarp */
		if (di > 4 
		 && userand == 0		/* Not other high D useful method */
		 && useqrand == 0
		 && usedsim == 0
		 && usepsim == 0)
			uselat = 1;

		/* Allocate space to record fixed steps */
		{
			fxlist_a = 4;
			if ((fxlist = (fxpos *)malloc(sizeof(fxpos) * fxlist_a)) == NULL)
			error ("Failed to malloc fxlist");
	}

	/* White color patches */
	if (esteps > 0)	{
		int j, e;

		sprintf(buf,"%d",esteps);
		pp->add_kword(pp, 0, "WHITE_COLOR_PATCHES",buf, NULL);
	
		for (j = 0; j < esteps; j++) {
			double val[MXTD], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];
	
			if (nmask & ICX_ADDITIVE) {
				for (e = 0; e < di; e++) {
					val[e] = 1.0;			/* White is full colorant */
				}
			} else {
				for (e = 0; e < di; e++) {
					val[e] = 0.0;			/* White is no colorant */
				}
			}

			/* Apply general filter */
			if (filter && dofilt(pdata, filt, val))
				continue;

			sprintf(buf,"%d",id++);
			ary[0].c = buf;

			if (xmask == nmask) {
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * val[e];
			} else {
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * (1.0 - val[e]);
			}

			pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
			ary[1 + di + 0].d = 100.0 * XYZ[0];
			ary[1 + di + 1].d = 100.0 * XYZ[1];
			ary[1 + di + 2].d = 100.0 * XYZ[2];
	
			pp->add_setarr(pp, 0, ary);
	
			if (fxlist != NULL) {		/* Note in fixed list */
				if (fxno >= fxlist_a) {
					fxlist_a *= 2;
					if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
						error ("Failed to malloc fxlist");
				}
				for (e = 0; e < di; e++)
					fxlist[fxno].p[e] = val[e];
				fxlist[fxno].eloc = pp->t[0].nsets;
				fxno++;
			}
		}
	}

	/* Black color patches */
	if (Bsteps > 0)	{
		int j, k, e;

		for (j = k = 0; j < Bsteps; j++) {
			double val[MXTD], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];
	
			if (nmask & ICX_ADDITIVE) {
				for (e = 0; e < di; e++) {
					val[e] = 0.0;			/* Black is no colorant */
				}
			} else {
				for (e = 0; e < di; e++) {
					val[e] = 1.0;			/* Black is full colorant */
				}
			}
	
			/* Apply general filter */
			if (filter && dofilt(pdata, filt, val))
				continue;

			/* Do a simple ink limit */
			if (uilimit < (double)di) {
				double tot = 0.0;
				for (e = 0; e < di; e++)
					tot += val[e];
				if (tot > uilimit) {
					for (e = 0; e < di; e++)
						val[e] *= uilimit/tot;
				}
			}

			sprintf(buf,"%d",id++);
			ary[0].c = buf;

			if (xmask == nmask) {
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * val[e];
			} else {
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * (1.0 - val[e]);
			}

			pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
			ary[1 + di + 0].d = 100.0 * XYZ[0];
			ary[1 + di + 1].d = 100.0 * XYZ[1];
			ary[1 + di + 2].d = 100.0 * XYZ[2];
	
			pp->add_setarr(pp, 0, ary);
	
			if (fxlist != NULL) {		/* Note in fixed list */
				if (fxno >= fxlist_a) {
					fxlist_a *= 2;
					if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
						error ("Failed to malloc fxlist");
				}
				for (e = 0; e < di; e++)
					fxlist[fxno].p[e] = val[e];
				fxlist[fxno].eloc = pp->t[0].nsets;
				fxno++;
			}
			k++;
		}
		sprintf(buf,"%d",k);
		pp->add_kword(pp, 0, "BLACK_COLOR_PATCHES",buf, NULL);
	
	}

	/* Primary (single channel) wedge steps */
	if (ssteps > 0)	{
		sprintf(buf,"%d",ssteps);
		pp->add_kword(pp, 0, "SINGLE_DIM_STEPS",buf, NULL);
		for (j = 0; j < di; j++) {
			for (i = 0; i < ssteps; i++) {
				int addp, e;
				double val[MXTD], XYZ[3];
				cgats_set_elem ary[1 + MXTD + 3];

				addp = 1;			/* Default add the point */

				for (e = 0; e < di; e++) {
					if (e == j)
						val[e] = (double)i/(ssteps-1);
					else
						val[e] = 0.0;
				}

				/* Extra power and dark emphasis */
				for (e = 0; e < di; e++)
					val[e] = icx_powlike(val[e], xpow * demph);

				/* See if it is already in the fixed list */
				if (!dontdedupe && fxlist != NULL) {
					int k;
					for (k = 0; k < fxno; k++) {
						for (e = 0; e < di; e++) {
							double tt;
							tt = fabs(fxlist[k].p[e] - val[e]);
							if (tt > MATCH_TOLL)
								break;			/* Not identical */
						}
						if (e >= di)
							break;				/* Was identical */
					}
					if (k < fxno)				/* Found an identical patch */
						addp = 0;				/* Don't add the point */
				}

				/* Apply general filter */
				if (filter && dofilt(pdata, filt, val))
					addp = 0;

				if (addp) {
		
					sprintf(buf,"%d",id++);
					ary[0].c = buf;

					if (xmask == nmask) {
						for (e = 0; e < di; e++)
							ary[1 + e].d = 100.0 * val[e];
					} else {
						for (e = 0; e < di; e++)
							ary[1 + e].d = 100.0 * (1.0 - val[e]);
					}

					pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
					ary[1 + di + 0].d = 100.0 * XYZ[0];
					ary[1 + di + 1].d = 100.0 * XYZ[1];
					ary[1 + di + 2].d = 100.0 * XYZ[2];
	
					pp->add_setarr(pp, 0, ary);
		
					if (fxlist != NULL) {		/* Note in fixed list */
						if (fxno >= fxlist_a) {
							fxlist_a *= 2;
							if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
								error ("Failed to malloc fxlist");
						}
						for (e = 0; e < di; e++)
							fxlist[fxno].p[e] = val[e];
						fxlist[fxno].eloc = -1;
						fxno++;
					}
				}
			}
		}
	}

	/* Composite gray wedge steps */
	if (gsteps > 0) {
		int cix[3];		/* Composite indexes */ 

		sprintf(buf,"%d",gsteps);
		pp->add_kword(pp, 0, "COMP_GREY_STEPS",buf, NULL);

		if (nmask & ICX_ADDITIVE) {	/* Look for the RGB */
			cix[0] = icx_ink2index(nmask, ICX_RED);
			cix[1] = icx_ink2index(nmask, ICX_GREEN);
			cix[2] = icx_ink2index(nmask, ICX_BLUE);

		} else {					/* Look for the CMY */
			cix[0] = icx_ink2index(nmask, ICX_CYAN);
			cix[1] = icx_ink2index(nmask, ICX_MAGENTA);
			cix[2] = icx_ink2index(nmask, ICX_YELLOW);
		}
		if (cix[0] < 0 || cix[1] < 0 || cix[2] < 0)
			error("Composite grey wedges aren't appropriate for %s device\n",ident);

		for (i = 0; i < gsteps; i++) {
			int addp, e;
			double sum, val[MXTD], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];

			addp = 1;			/* Default add the point */

			for (e = 0; e < di; e++) {
				if (e == cix[0] || e == cix[1] || e == cix[2])
					val[e] = (double)i/(gsteps-1);
				else
					val[e] = 0.0;
			}

			/* Extra power and dark emphasis */
			for (e = 0; e < di; e++)
				val[e] = icx_powlike(val[e], xpow * demph);

			/* Apply general filter */
			if (filter && dofilt(pdata, filt, val))
				addp = 0;

			/* Check if over ink limit */
			for (sum = 0.0, e = 0; e < di; e++)
				sum += val[e];
			if (sum > uilimit)
				addp = 0;

			/* See if it is already in the fixed list */
			if (!dontdedupe && fxlist != NULL) {
				int k;
				for (k = 0; k < fxno; k++) {
					for (e = 0; e < di; e++) {
						double tt;
						tt = fabs(fxlist[k].p[e] - val[e]);
						if (tt > MATCH_TOLL)
							break;			/* Not identical */
					}
					if (e >= di)
						break;				/* Was identical */
				}
				if (k < fxno)				/* Found an identical patch */
					addp = 0;				/* Don't add the point */
			}

			if (addp) {
	
				sprintf(buf,"%d",id++);
				ary[0].c = buf;

				if (xmask == nmask) {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * val[e];
				} else {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * (1.0 - val[e]);
				}

				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];

				pp->add_setarr(pp, 0, ary);
	
				if (fxlist != NULL) {		/* Note in fixed list */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxlist[fxno].eloc = -1;
					fxno++;
				}
			}
		}
	}

	/* Neutral wedge steps */
	if (nsteps > 0) {
		double lab[3];

		sprintf(buf,"%d",nsteps);
		pp->add_kword(pp, 0, "NEUTRAL_STEPS",buf, NULL);

		lab[1] = lab[2] = 0.0;
		for (i = 0; i < nsteps; i++) {
			int addp, e;
			double sum, val[MXTD], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];

			addp = 1;			/* Default add the point */

			lab[0] = 100.0 * (double)i/(nsteps-1);

			pdata->rLab_to_dev(pdata, val, lab);

			/* Apply general filter */
			if (filter && dofilt(pdata, filt, val))
				addp = 0;

			/* Check if over ink limit */
			for (sum = 0.0, e = 0; e < di; e++)
				sum += val[e];
			if (sum > uilimit)
				addp = 0;

			/* See if it is already in the fixed list */
			if (!dontdedupe && fxlist != NULL) {
				int k;
				for (k = 0; k < fxno; k++) {
					for (e = 0; e < di; e++) {
						double tt;
						tt = fabs(fxlist[k].p[e] - val[e]);
						if (tt > MATCH_TOLL)
							break;			/* Not identical */
					}
					if (e >= di)
						break;				/* Was identical */
				}
				if (k < fxno)				/* Found an identical patch */
					addp = 0;				/* Don't add the point */
			}

			if (addp) {
	
				sprintf(buf,"%d",id++);
				ary[0].c = buf;

				if (xmask == nmask) {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * val[e];
				} else {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * (1.0 - val[e]);
				}

				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];

				pp->add_setarr(pp, 0, ary);
	
				if (fxlist != NULL) {		/* Note in fixed list */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxlist[fxno].eloc = -1;
					fxno++;
				}
			}
		}
	}

	/* Regular Gridded Multi dimension steps */
	if (msteps > 0) {
		int gc[MXTD];			/* Grid coordinate */

		sprintf(buf,"%d",msteps);
		pp->add_kword(pp, 0, "MULTI_DIM_STEPS",buf, NULL);

		for (j = 0; j < di; j++)
			gc[j] = 0;			/* init coords */
			
		for (;;) {	/* For all grid points */
			double sum, val[MXTD], XYZ[3];
			int addp, e;

			/* If just surface points, reject points not on (2D) surface */
			if (msurf) {
				for (k = e = 0; e < di; e++) {
					if (gc[e] != 0 && gc[e] != (msteps-1))
						k++;
				}
				if (k > 2)
					goto next_point;
			}

			addp = 1;			/* Default add the point */

			for (e = 0; e < di; e++)
				val[e] = (double)gc[e]/(msteps-1);

			/* Extra power and dark emphasis */
			for (e = 0; e < di; e++)
				val[e] = icx_powlike(val[e], xpow * demph);

			/* Apply general filter */
			if (filter && dofilt(pdata, filt, val))
				addp = 0;

			/* Check if over ink limit */
			for (sum = 0.0, e = 0; e < di; e++)
				sum += val[e];
			if (sum > uilimit)
				addp = 0;		/* Don't add patches over ink limit */

			/* See if it is already in the fixed list */
			if (addp && !dontdedupe && fxlist != NULL) { 
				int k;
				for (k = 0; k < fxno; k++) {
					for (e = 0; e < di; e++) {
						double tt;
						tt = fabs(fxlist[k].p[e] - val[e]);
						if (tt > MATCH_TOLL)
							break;			/* Not identical */
					}
					if (e >= di)
						break;				/* Was identical */
				}
				if (k < fxno)				/* Found an identical patch */
					addp = 0;				/* Don't add the point */
			}

			/* Add patch to list if OK */
			if (addp) {
				cgats_set_elem ary[1 + MXTD + 3];

				sprintf(buf,"%d",id++);
				ary[0].c = buf;

				if (xmask == nmask) {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * val[e];
				} else {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * (1.0 - val[e]);
				}

				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];

				pp->add_setarr(pp, 0, ary);

				if (fxlist != NULL) {		/* Note in fixed list */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxlist[fxno].eloc = -1;
					fxno++;
				}
			}

			/* Increment grid index and position */
		  next_point:;
			for (j = 0; j < di; j++) {
				gc[j]++;
				if (gc[j] < msteps)
					break;	/* No carry */
				gc[j] = 0;
			}
			if (j >= di)
				break;		/* Done grid */
		}

#ifdef ADDRECCLIPPOINTS
		/* Add extra points that intersect */
		/* grid, and lie on ink limit plane */
		/* !!!!!!!!!! this doesn't cope with xpow !!!!!!!!!!! */
		if (uilimit < (di * 100.0)) {
			double val[MXTD], tv;
			double XYZ[3];
			for (k = 0; k < di; k++) {	/* dimension not on grid */
				for (j = 0; j < di; j++)
					gc[j] = 0;			/* init coords */
					
				for (;;) {	/* Until done */
					for (tv = 0.0, j = 0; j < di; j++) {
						if (j != k)
							tv += val[j] = (double)gc[j]/(msteps-1);
					}
					if (tv <= uilimit && (tv + 1.0) >= uilimit) {	/* Will intersect */
						double fr;
						val[k] = uilimit - tv; 	/* Point of intersection */
						fr = fmod((val[k] * msteps), 1.0);
						if (fr > 0.05 && fr < 0.95) {		/* Not within 5% of a grid point */ 
							int addp, e;
							cgats_set_elem ary[1 + MXTD + 3];
			
							addp = 1;			/* Default add the point */

							/* See if it is already in the fixed list */
							if (!dontdedupe && fxlist != NULL) {
								int k;
								for (k = 0; k < fxno; k++) {
									for (e = 0; e < di; e++) {
										double tt;
										tt = fabs(fxlist[k].p[e] - val[e]);
										if (tt > MATCH_TOLL)
											break;			/* Not identical */
									}
									if (e >= di)
										break;				/* Was identical */
								}
								if (k < fxno) {				/* Found an identical patch */
									addp = 0;				/* Don't add the point */
								}
							}

							/* Apply general filter */
							if (filter && dofilt(pdata, filt, val)) {
								addp = 0;
							}

							if (addp) {
								sprintf(buf,"%d",id++);
								ary[0].c = buf;
								if (xmask == nmask) {
									for (e = 0; e < di; e++)
										ary[1 + e].d = 100.0 * val[e];;
								} else {
									for (e = 0; e < di; e++)
										ary[1 + e].d = 100.0 * (1.0 - val[e]);
								}
								pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
								ary[1 + di + 0].d = 100.0 * XYZ[0];
								ary[1 + di + 1].d = 100.0 * XYZ[1];
								ary[1 + di + 2].d = 100.0 * XYZ[2];

								pp->add_setarr(pp, 0, ary);
	
								if (fxlist != NULL) {		/* Note in fixed list */
									if (fxno >= fxlist_a) {
										fxlist_a *= 2;
										if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
											error ("Failed to malloc fxlist");
									}
									for (e = 0; e < di; e++)
										fxlist[fxno].p[e] = val[e];
									fxlist[fxno].eloc = -1;
									fxno++;
								}
							}
						}
					}

					/* Increment grid index and position */
					for (j = 0; j < di; j++) {
						gc[j]++;
						if (j != k && gc[j] < msteps)
							break;	/* No carry */
						gc[j] = 0;
					}
					if (j >= di)
						break;		/* ALL done */
				}
			}
		}
#endif /* ADDRECCLIPPOINTS */
	}


	/* Regular body centered cubic gridded Multi dimension steps */
	if (bsteps > 0) {
		int gc[MXTD];			/* Grid coordinate */
		int pass = 0;			/* 0 = outer grid, 1 = inner grid */

		sprintf(buf,"%d",bsteps);
		pp->add_kword(pp, 0, "MULTI_DIM_BCC_STEPS",buf, NULL);

		for (pass = 0; pass < 2; pass++) {

			for (j = 0; j < di; j++)
				gc[j] = 0;			/* init coords */
				
			for (;;) {	/* For all grid points */
				double sum, val[MXTD], XYZ[3];
				int addp, e;

				addp = 1;			/* Default add the point */

				for (e = 0; e < di; e++)
					val[e] = (double)(pass * 0.5 + gc[e])/(bsteps-1);

				/* Extra power and dark emphasis */
				for (e = 0; e < di; e++)
					val[e] = icx_powlike(val[e], xpow * demph);

				/* Apply general filter */
				if (filter && dofilt(pdata, filt, val))
					addp = 0;

				/* Check if over ink limit */
				for (sum = 0.0, e = 0; e < di; e++)
					sum += val[e];
				if (sum > uilimit)
					addp = 0;		/* Don't add patches over ink limit */

				/* See if it is already in the fixed list */
				if (addp && !dontdedupe && fxlist != NULL) { 
					int k;
					for (k = 0; k < fxno; k++) {
						for (e = 0; e < di; e++) {
							double tt;
							tt = fabs(fxlist[k].p[e] - val[e]);
							if (tt > MATCH_TOLL)
								break;			/* Not identical */
						}
						if (e >= di)
							break;				/* Was identical */
					}
					if (k < fxno)				/* Found an identical patch */
						addp = 0;				/* Don't add the point */
				}

				/* Add patch to list if OK */
				if (addp) {
					cgats_set_elem ary[1 + MXTD + 3];

					sprintf(buf,"%d",id++);
					ary[0].c = buf;
					if (xmask == nmask) {
						for (e = 0; e < di; e++)
							ary[1 + e].d = 100.0 * val[e];
					} else {
						for (e = 0; e < di; e++)
							ary[1 + e].d = 100.0 * (1.0 - val[e]);
					}

					pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
					ary[1 + di + 0].d = 100.0 * XYZ[0];
					ary[1 + di + 1].d = 100.0 * XYZ[1];
					ary[1 + di + 2].d = 100.0 * XYZ[2];

					pp->add_setarr(pp, 0, ary);

					if (fxlist != NULL) {		/* Note in fixed list */
						if (fxno >= fxlist_a) {
							fxlist_a *= 2;
							if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
								error ("Failed to malloc fxlist");
						}
						for (e = 0; e < di; e++)
							fxlist[fxno].p[e] = val[e];
						fxlist[fxno].eloc = -1;
						fxno++;
					}
				}

				/* Increment grid index and position */
				for (j = 0; j < di; j++) {
					gc[j]++;
					if ((pass == 0 && gc[j] < bsteps)
					 || (pass == 1 && gc[j] < (bsteps-1)))
						break;	/* No carry */
					gc[j] = 0;
				}
				if (j >= di)
					break;		/* Done grid */
			}
		}
	}

	if (fsteps > fxno) { /* Top up with full spread (perceptually even) and other patch types */

		/* Generate device random numbers. Don't check for duplicates */
		if (userand == 1 || useqrand == 1) {
			int i, j;
			sobol *sl = NULL;

			if (useqrand) {
				if ((sl = new_sobol(di)) == NULL)
					error("Creating sobol sequence generator failed");
			}

			/* Create more points up to fsteps */
			if (verb)
				printf("\n");
			for (j = 0, i = fxno; i < fsteps;) {
				int e;
				double sum;
				double val[MXTD], XYZ[3];
				cgats_set_elem ary[1 + MXTD + 3];

				if (sl != NULL) {
					if (sl->next(sl, val))
						error("Run out of sobol random numbers!");
				} else {	/* else uniform random distribution */
					for (e = 0; e < di; e++)
						val[e] = d_rand(0.0, 1.0);
				}

				/* Extra power and dark emphasis */
				for (e = 0; e < di; e++)
					val[e] = icx_powlike(val[e], xpow * demph);

				/* Apply general filter */
				if (filter && dofilt(pdata, filt, val))
					continue;

				/* Check if over ink limit */
				for (sum = 0.0, e = 0; e < di; e++)
					sum += val[e];
				if (sum > uilimit)
					continue;

				sprintf(buf,"%d",id++);
				ary[0].c = buf;
				if (xmask == nmask) {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * val[e];
				} else {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * (1.0 - val[e]);
				}

				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];
	
				pp->add_setarr(pp, 0, ary);

				if (fxlist != NULL) {		/* Note in fixed list to allow stats later */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxlist[fxno].eloc = -1;
					fxno++;
				}

				if (verb) {
					printf("%cAdded %d/%d",cr_char,i+1,fxno); fflush(stdout);
				}
				i++, j++;
			}
			if (verb)
				printf("\n");

			sprintf(buf,"%d",j);
			if (sl != NULL)
				pp->add_kword(pp, 0, "SPACEFILLING_RANDOM_PATCHES", buf, NULL);
			else
				pp->add_kword(pp, 0, "RANDOM_DEVICE_PATCHES", buf, NULL);

			if (sl != NULL)
				sl->del(sl);

		} else {
//			ppoint *s = NULL;
			ofps *s = NULL;
			ifarp *t = NULL;
			simdlat *dx = NULL;
			simplat *px = NULL;
			prand *rx = NULL;

			/* (Note that the ink limit for these algorithms won't take into account the xpow), */
			/* and that we're not applying the filter until after generation, so the */
			/* number of patches won't reach the target. This could be fixed fairly easily */
			/* for some of these (new_prand).  */
			if (uselat)	 {
				/* A "greedy"/incremental far point approach */
				t = new_ifarp(verb, di, uilimit, fsteps, fxlist, fxno,
				                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "IFP_PATCHES", buf, NULL);
			} else if (usedsim) {
				/* Device space simplex latice test points */
				dx = new_simdlat(di, uilimit, fsteps, fxlist, fxno, SIMDLAT_TYPE, simangle,
			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "SIMPLEX_DEVICE_PATCHES", buf, NULL);
			} else if (usepsim) {
				/* Perceptual space simplex latice test points */
				px = new_simplat(di, uilimit, fsteps, fxlist, fxno, simangle,
			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "SIMPLEX_PERCEPTUAL_PATCHES", buf, NULL);
			} else if (userand == 2 || useqrand == 2) {
				/* Perceptual random test points */
				rx = new_prand(di, uilimit, fsteps, fxlist, fxno, useqrand == 2 ? 1 : 0,
			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "RANDOM_PERCEPTUAL_PATCHES", buf, NULL);

			} else {		/* Default full spread algorithm */
				/* Optimised Farthest Point Sampling */
				s = new_ofps(verb, di, uilimit, fsteps, good,
				            dadapt, 1.0 - perc_wght, perc_wght, curv_wght, fxlist, fxno,
			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "OFPS_PATCHES", buf, NULL);
			}

	
			for (;;) {
				int e;
				double XYZ[3], val[MXTD];
				cgats_set_elem ary[1 + MXTD + 3];
				if (( s ? s->read(s, val, NULL) :
				      t ? t->read(t, val, NULL) :
				     dx ? dx->read(dx, val, NULL) :
				     rx ? rx->read(rx, val, NULL) :
				          px->read(px, val, NULL)))
					break;
	
				/* Filter out silly values from ppoint */
				for (e = 0; e < di; e++) {
					if (val[e] < 0.001)
						val[e] = 0.0;
					else if (val[e] > 0.999)
						val[e] = 1.0;
				}
			
				/* Apply general filter */
				if (filter && dofilt(pdata, filt, val))
					continue;

				/* Do a simple ink limit */
				if (uilimit < (double)di) {
					double tot = 0.0;
					for (e = 0; e < di; e++)
						tot += val[e];
					if (tot > uilimit) {
						for (e = 0; e < di; e++)
							val[e] = val[e] * uilimit/tot;
					}
				}

				sprintf(buf,"%d",id++);
				ary[0].c = buf;
				if (xmask == nmask) {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * val[e];
				} else {
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * (1.0 - val[e]);
				}

				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];
	
				pp->add_setarr(pp, 0, ary);

				if (fxlist != NULL) {		/* Note in fixed list to allow stats later */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxlist[fxno].eloc = -1;
					fxno++;
				}
			}
			(s ? s->del(s) : t ? t->del(t) : dx ? dx->del(dx) : rx ? rx->del(rx) : px->del(px));
		}
	}

	/* Even the location of marked patches into sequence */
	{
		int ii, p1, p2, t1;

		/* For each patch to be dispersed */
		for (ii = 0; ii < fxno; ii++) {
			if (fxlist[ii].eloc >= 0) {
				p1 = fxlist[ii].eloc;

				for (k = 0; k < 10; k++) {		/* Retry 10 times */

					/* Pick a random patch to exchange it with */
					p2 = i_rand(0, pp->t[0].nsets-1);
	
					/* Check it isn't one of our patches to be dispersed */
					for (i = 0; i < fxno; i++) {
						if (fxlist[i].eloc == p2)
							break;
					}
					if (i < fxno)
						continue;		/* Try another patch to exchange with */

					/* Swap */
					for (j = 1; j < (1 + di + 3); j++) {
						double tt = *((double *)pp->t[0].fdata[p1][j]);
						*((double *)pp->t[0].fdata[p1][j]) = *((double *)pp->t[0].fdata[p2][j]);
						*((double *)pp->t[0].fdata[p2][j]) = tt;
					}
					fxlist[ii].eloc = p2;

					break;
				}
			}
		}
	}

	/* If this seems to be for a CRT, optimise the patch order to minimise the */
	/* response time delays */
	if (nmask == ICX_RGB && pp->t[0].nsets > 1 && dontreorder == 0) {
		int npat = pp->t[0].nsets;
		char *nm;						/* Don't move array */
		double udelay, *delays, adelay;
		double temp, trate;	/* Annealing temperature & rate */
		double tstart, tend;/* Annealing chedule range */

		if ((nm = (char *)malloc(sizeof(char) * npat)) == NULL)
			error ("Failed to malloc nm array");
		if ((delays = (double *)malloc(sizeof(double) * npat)) == NULL)
			error ("Failed to malloc delay array");

		/* Set nm[] to mark patches that shouldn't be moved */
		for (i = 0; i < npat; i++)
			nm[i] = 0;
		for (i = 0; i < fxno; i++) {
			if (fxlist[i].eloc >= 0)
				nm[fxlist[i].eloc] = 1;
		}

#ifdef NEVER
		/* Randomly shuffle patches */
		{
			int p1, p2;

			for (p1 = 0; p1 < npat; p1++) {
	
				if (nm[p1])
					continue;

				p2 = i_rand(0, npat-1);
				if (nm[p2])
					continue;
				for (j = 1; j < (1 + di + 3); j++) {
					double tt = *((double *)pp->t[0].fdata[p1][j]);
					*((double *)pp->t[0].fdata[p1][j]) = *((double *)pp->t[0].fdata[p2][j]);
					*((double *)pp->t[0].fdata[p2][j]) = tt;
				}
			}
		}
#endif
#ifdef NEVER
		/* Simple sort by brightness */
		{
			int p1, p2;
			double rgb1, rgb2;

			for (p1 = 0; p1 < (npat-1); p1++) {
	
				if (nm[p1])
					continue;

				rgb1 = pow(*((double *)pp->t[0].fdata[p1][1 + 0]), 2.2)
				     + pow(*((double *)pp->t[0].fdata[p1][1 + 1]), 2.2)
				     + pow(*((double *)pp->t[0].fdata[p1][1 + 2]), 2.2);
	
				for (p2 = p1 + 1; p2 < npat; p2++) {

					if (nm[p2])
						continue;

					rgb2 = pow(*((double *)pp->t[0].fdata[p2][1 + 0]), 2.2)
					     + pow(*((double *)pp->t[0].fdata[p2][1 + 1]), 2.2)
					     + pow(*((double *)pp->t[0].fdata[p2][1 + 2]), 2.2);
		
					if (rgb2 < rgb1) {
						for (j = 1; j < (1 + di + 3); j++) {
							double tt = *((double *)pp->t[0].fdata[p1][j]);
							*((double *)pp->t[0].fdata[p1][j]) = *((double *)pp->t[0].fdata[p2][j]);
							*((double *)pp->t[0].fdata[p2][j]) = tt;
						}
						rgb1 = rgb2;
					}
				}
			}
		}
#endif
		/* Compute the current overall update delay */
		udelay = 0.0;
		for (i = 1; i < npat; i++) {
			double xdelay;

			xdelay = disprespt(pp, i-1, i);

			delays[i] = xdelay;
//printf("~1 delay[%d] = %f\n",i,xdelay);
			udelay += xdelay;
		}

		if (verb)
			printf("Extra display response delay = %f sec., optimizing....\n",udelay);

		{
			int nchunks, chsize;
			int chstart, chend;

			if (verb)
				printf("%c%2d%%",cr_char,0); fflush(stdout);

			/* We'll do this in chunks of 500 to make it linear time overall, */
			/* at the cost of the best possible optimisation. */
			nchunks = (int)ceil(npat/500.0); 
			chsize = (int)ceil(npat/nchunks);
			for (chstart = 0; chstart < npat; chstart += chsize) {
				int p1, p2, bp2;
				double p1d, p2d, p1d1, p2d1;
				double p1nd, p2nd, p1nd1, p2nd1;
				double tdelay, bdelay, de;
				int noswapped;

				chend = chstart + chsize+2;
				if (chend > npat)
					chend = npat;
				noswapped = chend - chstart;
//printf("~1 chstart %d, chend %d, size %d\n",chstart,chend, chend - chstart);

				/* While we are still improving, and the improvement was significant */
				for (;noswapped > 2;) {
					noswapped = 0;

					for (p1 = chstart + 1; p1 < chend; p1++) {
						if (nm[p1])
							continue;

						p1d = delays[p1];

						/* Locate the patch ahead of us that is best to swap with */
						bp2 = -1;
						bdelay = udelay;
						for (p2 = p1 + 2; p2 < chend; p2++) {

							if (nm[p2])
								continue;

							/* Compute effect of a swap on the total delay */
							p2d = delays[p2];
							p1nd = disprespt(pp, p2-1, p1);
							p2nd = disprespt(pp, p1-1, p2);
							p1d1 = p1nd1 = 0.0;
							if ((p1+1) < chend) {
								p1d1 = delays[p1+1];
								p1nd1 = disprespt(pp, p2, p1+1);
							}
							p2d1 = p2nd1 = 0.0;
							if ((p2+1) < chend) {
								p2d1 = delays[p2+1];
								p2nd1 = disprespt(pp, p1, p2+1);
							}
			
							tdelay = udelay - p1d - p2d - p1d1 - p2d1 + p1nd + p2nd + p1nd1 + p2nd1;

							if (tdelay < bdelay)		/* Improve it */
//							if (tdelay > bdelay)		/* Make it worse */
							{
								bp2 = p2;
								bdelay = tdelay;
							}
						}
						if (bp2 < 0) {
							continue;
						} 

						/* Swap the patches */
						noswapped++;

						p2 = bp2;

						p2d = delays[p2];
						p1nd = disprespt(pp, p2-1, p1);
						p2nd = disprespt(pp, p1-1, p2);
						p1d1 = p1nd1 = 0.0;
						if ((p1+1) < chend) {
							p1d1 = delays[p1+1];
							p1nd1 = disprespt(pp, p2, p1+1);
						}
						p2d1 = p2nd1 = 0.0;
						if ((p2+1) < chend) {
							p2d1 = delays[p2+1];
							p2nd1 = disprespt(pp, p1, p2+1);
						}

						tdelay = udelay - p1d - p2d - p1d1 - p2d1 + p1nd + p2nd + p1nd1 + p2nd1;

						/* Swap the values */
						udelay = tdelay;
						delays[p2] = p1nd;
						delays[p1] = p2nd;
						if (p1 < (chend-1))
							delays[p1+1] = p1nd1;
						if (p2 < (chend-1))
							delays[p2+1] = p2nd1;

						for (j = 1; j < (1 + di + 3); j++) {
							double tt = *((double *)pp->t[0].fdata[p1][j]);
							*((double *)pp->t[0].fdata[p1][j]) = *((double *)pp->t[0].fdata[p2][j]);
							*((double *)pp->t[0].fdata[p2][j]) = tt;
						}
//printf("~1 swaping %d and %d, udelay %f\n",p1,p2,udelay);
					}
//printf("~1 udelay %f\n",udelay);
					if (verb) {
						printf("%c%2d%%",cr_char,(int)(100.0 * (chend-1 - noswapped)/(npat-1.0)));
						fflush(stdout);
					}
				}
			}
			if (verb)
				printf("%c%2d%%",cr_char,100); fflush(stdout);
		}
		if (verb)
			printf("\nOptimised display response delay = %f sec.\n",udelay);

		free(delays);
		free(nm);
	}

	/* Use ofps to measure the stats of the points */
	/* Note that if new_ofps() fails it will exit() */
	if (verb > 1
     && di <= 4
	 && (userand || useqrand || usedsim || usepsim || uselat)) {
		ofps *s;
		printf("Computing device space point stats:\n");
		if ((s = new_ofps(verb, di, uilimit, fxno, 0, 0.0, 0.0, 0.0, 0.0, fxlist, fxno,
	           (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata)) == NULL) {
			printf("Failed to compute stats\n");
		} else {
			s->stats(s);
			printf("Max distance stats: Min = %f, Average = %f, Max = %f\n",s->mn,s->av,s->mx);
			s->del(s);
		}
	}

	/* Add the eight entries in the second table. */
	/* These are legal device values that we think may */
	/* give all combinations of min/max CMY density values. */
	/* These are typically used for DTP51 and DTP41 patch separators. */
	{
		int i;

		pp->add_kword(pp, 1, "DENSITY_EXTREME_VALUES", "8", NULL);

		for (i = 0; i < 8; i++) {
			int e;
			double den[4], val[MXTD], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];

			/* Setup target density combination */
			for (e = 0; e < 3; e++) {
				if (i & (1 << e)) 
					den[e] = 2.5;
				else
					den[e] = -0.5;
			}

			/* Lookup device values for target density */
			pdata->den_to_dev(pdata, val, den);	

			/* Do a simple ink limit */
			if (uilimit < (double)di) {
				double tot = 0.0;
				for (e = 0; e < di; e++)
					tot += val[e];
				if (tot > uilimit) {
					for (e = 0; e < di; e++)
						val[e] *= uilimit/tot;
				}
			}

			ary[0].i = i;

			if (xmask == nmask) {
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * val[e];
			} else {
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * (1.0 - val[e]);
			}
			pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
			ary[1 + di + 0].d = 100.0 * XYZ[0];
			ary[1 + di + 1].d = 100.0 * XYZ[1];
			ary[1 + di + 2].d = 100.0 * XYZ[2];

			pp->add_setarr(pp, 1, ary);

		}
	}

	/* Add the nine entries in the third table. */
	/* These are legal device values that we calculate */
	/* give all combinations of typical CMY device values + 50% CMY */
	/* These are typically use for DTP20 bar coding. */
	{
		int i;
	
		icxColorantLu *ftarg = NULL;

		/* If not possible to use native space, use fake CMY */
		if ((nmask & ICX_CMYK) != ICX_CMYK
		 && (nmask & ICX_CMY) != ICX_CMY
		 && (nmask & ICX_RGB) != ICX_RGB) {
			if ((ftarg = new_icxColorantLu(ICX_CMY)) == NULL)
				error ("Creation of xcolorant lu object failed");
		}
		
		pp->add_kword(pp, 2, "DEVICE_COMBINATION_VALUES", "9", NULL);

		for (i = 0; i < 9; i++) {
			int e;
			double val[MXTD], lab[3], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];

			for (e = 0; e < di; e++)
				val[e] = 0.0;

			/* Setup target device combination */
			/* Order must be White, Cyan, Magenta, Blue Yellow Green Red Black */
			if (ftarg != NULL || (nmask & ICX_CMY) == ICX_CMY) {
				for (e = 0; e < 3; e++) {
					if (i & (1 << e)) 
						val[e] = 1.0;
					else
						val[e] = 0.0;
				}
				if (i == 7)
					val[3] = 1.0;
				if (i == 8) {			/* Special 50/50/50 grey for DTP20 */
					val[0] = val[1] = val[2] = 0.5;
					val[3] = 0.0;
				}

			} else {	/* RGB like */
				for (e = 0; e < 3; e++) {
					if (i & (1 << e)) 
						val[e] = 0.0;
					else
						val[e] = 1.0;
				}
				if (i == 8)
					val[0] = val[1] = val[2] = 0.5;
			}
			
			/* Apply extra power to device values (??) */
			for (e = 0; e < di; e++)
				val[e] = icx_powlike(val[e], xpow);

			/* If target space isn't something we recognise, convert it */
			if (ftarg != NULL) {
				ftarg->dev_to_rLab(ftarg, lab, val);
				pdata->rLab_to_dev(pdata, val, lab);	
			}

			/* Do a simple ink limit */
			if (uilimit < (double)di) {
				double tot = 0.0;
				for (e = 0; e < di; e++)
					tot += val[e];
				if (tot > uilimit) {
					for (e = 0; e < di; e++)
						val[e] *= uilimit/tot;
				}
			}

			ary[0].i = i;
			if (xmask == nmask) {
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * val[e];
			} else {
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * (1.0 - val[e]);
			}

			pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
			ary[1 + di + 0].d = 100.0 * XYZ[0];
			ary[1 + di + 1].d = 100.0 * XYZ[1];
			ary[1 + di + 2].d = 100.0 * XYZ[2];

			pp->add_setarr(pp, 2, ary);
		}

		if (ftarg != NULL)
			ftarg->del(ftarg);
	}

	ttime = clock() - stime;
	if (verb) {
		printf("Total number of patches = %d\n",pp->t[0].nsets);
		if (id < (1 + (1 << di)))
			printf("WARNING : not enough patches for %d channels, need at least %d\n",di,(1 + (1 << di)));
		printf("Execution time = %f seconds\n",ttime/(double)CLOCKS_PER_SEC);
	}

	if (pp->write_name(pp, fname))
		error("Write error : %s",pp->err);

#ifdef VRML_DIAG		/* Dump a VRML/X3D of the resulting points */
	if (dumpvrml & 1) {	/* Lab space */
		vrml *wrl;
		int nsets = pp->t[0].nsets;
		double rad;
		double dev[MXTD], Lab[3], col[3];
		int doaxes = 1;					/* Do axes */

		if ((wrl = new_vrml(wlname, doaxes, 0)) == NULL)
			error("new_vrml failed for '%s%s'",wlname,vrml_ext());

		/* Fudge sphere diameter */
		rad = 15.0/pow(nsets, 1.0/(double)(di <= 3 ? di : 3));

		for (i = 0; i < nsets; i++) {
			/* Re-do any inversion before using dev_to_rLab() */
			if (xmask == nmask) {
				for (j = 0; j < di; j++)
					dev[j] = 0.01 * *((double *)pp->t[0].fdata[i][j + 1]);
			} else {
				for (j = 0; j < di; j++)
					dev[j] = 0.01 * (100.0 - *((double *)pp->t[0].fdata[i][j + 1]));
			}
			pdata->dev_to_rLab(pdata, Lab, dev);
			wrl->Lab2RGB(wrl, col, Lab);

			wrl->add_marker(wrl, Lab, col, rad);
		}
		wrl->del(wrl);		/* Write file and delete */
	}
	if (dumpvrml & 2) {	/* Device space */
		vrml *wrl;
		int nsets = pp->t[0].nsets;
		double rad;
		double dev[MXTD], idev[MXTD], Lab[3], col[3];

		wrl = new_vrml(wdname, 0, 0);

		/* Fudge sphere diameter */
		rad = 15.0/pow(nsets, 1.0/(double)(di <= 3 ? di : 3));

		for (i = 0; i < nsets; i++) {

			/* Re-do any inversion before using dev_to_rLab() */
			if (xmask == nmask) {
				for (j = 0; j < di; j++)
					idev[j] = dev[j] = 0.01 * *((double *)pp->t[0].fdata[i][j + 1]);
			} else {
				for (j = 0; j < di; j++) {
					dev[j] = 0.01 * *((double *)pp->t[0].fdata[i][j + 1]);
					idev[j] = 1.0 - dev[j];
				}
			}

			pdata->dev_to_rLab(pdata, Lab, idev);
			wrl->Lab2RGB(wrl, col, Lab);

			/* Fudge device locations into "Lab" space */
			Lab[0] = 100.0 * dev[0];
			Lab[1] = 100.0 * dev[1] - 50.0;
			Lab[2] = 100.0 * dev[2] - 50.0;

			wrl->add_marker(wrl, Lab, col, rad);
		}
		wrl->del(wrl);		/* Write file and delete */
	}
#endif /* VRML_DIAG */

	pdata->del(pdata);	/* Cleanup perceptual conversion */

	free(ident);
	pp->del(pp);		/* Clean up */
	if (fxlist != NULL)
		free(fxlist);

	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - -  */

#ifdef NEVER

/* Compte the display response time */
static double disprespt2(cgats *pp, int p1, int p2) {
	double kr, kf;
	double orgb[3], rgb[3];
	double xdelay = 0.0;
	int j;

	kr = DISPLAY_RISE_TIME/log(1 - 0.9);	/* Exponent constant */
	kf = DISPLAY_FALL_TIME/log(1 - 0.9);	/* Exponent constant */

	orgb[0] = *((double *)pp->t[0].fdata[p1][1 + 0]) / 100.0;
	orgb[1] = *((double *)pp->t[0].fdata[p1][1 + 1]) / 100.0;
	orgb[2] = *((double *)pp->t[0].fdata[p1][1 + 2]) / 100.0;

	rgb[0] = *((double *)pp->t[0].fdata[p2][1 + 0]) / 100.0;
	rgb[1] = *((double *)pp->t[0].fdata[p2][1 + 1]) / 100.0;
	rgb[2] = *((double *)pp->t[0].fdata[p2][1 + 2]) / 100.0;

	for (j = 0; j < 3; j++) {
		double el, dl, n, t;

		el = pow(rgb[j], 2.2);
		dl = el - pow(orgb[j], 2.2);	/* Change in level */
		if (fabs(dl) > 0.01) {		/* More than 1% change in level */
			n = DISPLAY_SETTLE_AIM2 * el;
			if (n < DISPLAY_ABS_AIM2)
				n = DISPLAY_ABS_AIM2;
			if (dl > 0.0)
				t = kr * log(n/dl);
			else
				t = kf * log(n/-dl);

			if (t > xdelay)
				xdelay = t;
		}
	}
	return xdelay;
}

#endif

/* Compte the display response time */
static double disprespt(cgats *pp, int p1, int p2) {
	double orgb[3], nrgb[3];
	double xdelay = 0.0;

	orgb[0] = *((double *)pp->t[0].fdata[p1][1 + 0]) / 100.0;
	orgb[1] = *((double *)pp->t[0].fdata[p1][1 + 1]) / 100.0;
	orgb[2] = *((double *)pp->t[0].fdata[p1][1 + 2]) / 100.0;

	nrgb[0] = *((double *)pp->t[0].fdata[p2][1 + 0]) / 100.0;
	nrgb[1] = *((double *)pp->t[0].fdata[p2][1 + 1]) / 100.0;
	nrgb[2] = *((double *)pp->t[0].fdata[p2][1 + 2]) / 100.0;

	xdelay = disp_settle_time(orgb, nrgb, DISPLAY_RISE_TIME, DISPLAY_FALL_TIME, DISPLAY_SETTLE_AIM);

	return xdelay;
}

