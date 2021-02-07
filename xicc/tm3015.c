
/* 
 * Author:  Graeme W. Gill
 * Date:    4/1/19
 * Version: 1.00
 *
 * Copyright 2019 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 2 or later :-
 * see the License2.txt file for licencing details.
 *
 * Based on icx_CIE1995_CRI() and icx_EBU2012_TLCI()
 */

/*
 * These functions compute IES TM-30-15,
 * "IES Method for Evaluating Light Source Color Rendition"
 */

/*
 * TTBD:
 *
 */

#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <math.h>
# include "aconfig.h"
#ifndef SALONEINSTLIB
# include "numlib.h"
#  include "plot.h"			/* For debugging */
#  include "ui.h"
#else
# include "numsup.h"
# include "sa_conv.h"
#endif
#include "cgats.h"
#include "conv.h"
#include "xspect.h"

#undef DEBUG
#undef TRACE_CAT
#undef USE_5NM		/* To match spreadsheet */

#ifdef DEBUG
# define DBGA g_log, 0 		/* First argument to DBGF() */
# define DBGF(xx)	a1logd xx ;
#else
# define DBG(xx)
# define DBGF(xx)
#endif

#ifdef DEBUG
# define REFTRACE(xxxx) fprintf xxxx ;
# define REFOUT stderr
#else
# define REFTRACE(xxxx) 
#endif

#ifndef STANDALONE_TEST

static double MCAT02[3][3] = {
	{  0.7328, 0.4296, -0.1624 },
	{ -0.7036, 1.6975,  0.0061 },
	{  0.0030, 0.0136,  0.9834 }
};

static double MHPE[3][3] = {
	{  0.38971, 0.68898, -0.07868 },
	{ -0.22981, 1.18340,  0.04641 },
	{  0.00000, 0.00000,  1.00000 }
};

/* Create the CIECAM02 conversion from XYZ10 to R'G'B' given the XYZ white point */
/* Returns Aw */
static double comp_rgbd_mtx(double mat[3][3], double wXYZ[3]) {
	double wRGB[3];			/* White RGB values */
	double vkmat[3][3];		/* Von Kries matrix */
	double iMCAT02[3][3];
	double wRGBd[3], wRGBda[3];
	double Aw;
	int j;

	// mat = xyz -> RGB /= RwGwBw -> inv MCAT02 -> MHPE

	REFTRACE((REFOUT, "comp_rgbd_mtx with XYZ %f %f %f\n", wXYZ[0], wXYZ[1], wXYZ[2]))

	// Start with CAT matrix
	icmCpy3x3(mat, MCAT02);

	/* Setup the Von Kries white point adaptation matrix */
	icmSetUnity3x3(vkmat);
	icmMulBy3x3(wRGB, MCAT02, wXYZ); // Wxyz -> MCAT02 -> RwGwBw
	vkmat[0][0] = 1.0/wRGB[0];
 	vkmat[1][1] = 1.0/wRGB[1];
 	vkmat[2][2] = 1.0/wRGB[2];
	icmMul3x3(mat, vkmat);

	icmInverse3x3(iMCAT02, MCAT02);
	icmMul3x3(mat, iMCAT02);

	icmMul3x3(mat, MHPE);

	/* White XYZ to chromaticly adjusted RGB' */
	icmMulBy3x3(wRGBd, mat, wXYZ);

	for (j = 0; j < 3; j++) {
		double tt = pow(0.7937 * wRGBd[j], 0.42);
		wRGBda[j] = (400.0 * tt) / (27.13 + tt) + 0.1;
	}

	/* Achromatic response of white */
	Aw = (2.0 * wRGBda[0] + wRGBda[1] + (1.0/20.0) * wRGBda[2] - 0.305) * 1.0003;

	REFTRACE((REFOUT, " wRGBda %f %f %f Aw %f\n", wRGBda[0], wRGBda[1], wRGBda[2], Aw))

	return Aw;
}

#ifdef TRACE_CAT
static void apply_mtx_trace(double rgb[3], double inXYZ[3], double wXYZ[3]) {
	double wRGB[3];			/* White RGB values */
	double iMCAT02[3][3];
	double sXYZ[3];			/* Sample RGB */
	double sRGB[3];			/* Sample RGB */

	printf("\napply_mtx_trace:\n");

	printf(" wXYZ %f %f %f\n", wXYZ[0], wXYZ[1], wXYZ[2]);

	// Wxyz -> MCAT02 -> RwGwBw
	icmMulBy3x3(wRGB, MCAT02, wXYZ);

	printf(" wRGB %f %f %f\n", wRGB[0], wRGB[1], wRGB[2]);

	// ---------------------------

	//mat =  xyz -> RGB /= RwGwBw -> inv MCAT02 -> MHPE

	icmCpy3(sXYZ, inXYZ);

	printf(" sXYZ %f %f %f\n", sXYZ[0], sXYZ[1], sXYZ[2]);

	icmMulBy3x3(sRGB, MCAT02, sXYZ);

	printf(" sRGB %f %f %f\n", sRGB[0], sRGB[1], sRGB[2]);

	sRGB[0] *= 1.0/wRGB[0];
 	sRGB[1] *= 1.0/wRGB[1];
 	sRGB[2] *= 1.0/wRGB[2];

	printf(" CAT sRGB %f %f %f\n", sRGB[0], sRGB[1], sRGB[2]);

	icmInverse3x3(iMCAT02, MCAT02);
	icmMulBy3x3(sXYZ, iMCAT02, sRGB);

	printf(" CAT sXYZ %f %f %f\n", sXYZ[0], sXYZ[1], sXYZ[2]);

	icmMulBy3x3(rgb, MHPE, sXYZ);

	printf(" sRGBd %f %f %f\n", rgb[0], rgb[1], rgb[2]);
	printf("\n");
}
#endif /* TRACE_CAT */

/* Convert XYZ to CIECAM02 UCS */
/* mat is the XYZ10 to R'G'B' for the given white (see above) */
/* Aw is the acromatic response of white (see above) */
static void ciecam02_ucs(double Jabd[3], double Aw, double mat[3][3], double inXYZ[3]) {
	double iRGBd[3];
	double iRGBda[3];
	double a, b;
	double rS, ttd;
	double A, J, h, e, ss, C, M;
	double Jd, ad, bd, Md;
	int j;

	REFTRACE((REFOUT, "ciecam02_ucs with XYZ %f %f %f\n", inXYZ[0], inXYZ[1], inXYZ[2]))

	/* XYZ to chromaticly adjusted RGB' */
	icmMulBy3x3(iRGBd, mat, inXYZ);

	REFTRACE((REFOUT, " iRGBd %f %f %f\n", iRGBd[0], iRGBd[1], iRGBd[2]))

	/* Post-adapted cone response of white */
	for (j = 0; j < 3; j++) {
		double tt = pow(0.7937 * iRGBd[j], 0.42);
		iRGBda[j] = (400.0 * tt) / (27.13 + tt) + 0.1;
	}

	REFTRACE((REFOUT, " iRGBda %f %f %f\n", iRGBda[0], iRGBda[1], iRGBda[2]))

	/* Preliminary red-green & yellow-blue opponent dimensions */
	a    = iRGBda[0] - 12.0 * iRGBda[1]/11.0 + iRGBda[2]/11.0;
    b    = (1.0/9.0) * (iRGBda[0] + iRGBda[1] - 2.0 * iRGBda[2]);
	rS   = sqrt(a * a + b * b);		/* Normalised a, b */

	/* Preliminary Saturation */
	ttd = iRGBda[0] + iRGBda[1] + (21.0/20.0) * iRGBda[2];

	/* Achromatic response */
	/* Note that the minimum values of rgba[] for XYZ = 0 is 0.1, */
	/* hence magic 0.305 below comes from the following weighting of rgba[], */
	/* to base A at 0.0 */
	A = (2.0 * iRGBda[0] + iRGBda[1] + (1.0/20.0) * iRGBda[2] - 0.305) * 1.0003;

	REFTRACE((REFOUT,  " a = %f, b = %f, ttd = %f, rS = %f, A = %f\n", a, b, ttd, rS, A))

	/* Lightness */
	J = pow(A/Aw, 0.69 * 1.9272);		/* J/100  - keep Sign */

	/* Hue angle */
    h = (180.0/DBL_PI) * atan2(b, a);
	h = (h < 0.0) ? h + 360.0 : h;

	REFTRACE((REFOUT, " J = %f, h = %f\n", J, h))

	/* Eccentricity factor */
	e = 0.25 * (cos(h * DBL_PI/180.0 + 2.0) + 3.8);

	ss = (50000.0/13.0 * 1.0 * 1.0003 * e * rS) / ttd;

	/* Chroma */
	C = pow(ss, 0.9) * sqrt(J) * 0.895217884;	/* 0.895217884 = (1.64 - 0.29^0.2)^.73 */

	M = C * 0.943874156;		/* 0.943874156 = 0.7937^0.25 */
	
	REFTRACE((REFOUT, " e = %f ss = %f C = %f M = %f\n", e, ss, C, M))

	J *= 100.0;		/* Scale J to 100 */

	Md = (1.0/ 0.0228) * log(1.0 + 0.0228 * M);

	REFTRACE((REFOUT, " J = %f Md = %f\n", J, Md))

	/* Compute CIECAM02 UCS value */
	Jd = (1.0 + 100.0 * 0.007) * J/(1.0 + 0.007 * J);
	ad = Md * cos(h * DBL_PI/180.0);
	bd = Md * sin(h * DBL_PI/180.0);

	REFTRACE((REFOUT, " Jabd %f %f %f\n",Jd, ad, bd))

	Jabd[0] = Jd;
	Jabd[1] = ad;
	Jabd[2] = bd;
}

/* IES TM-30-15 evaluation samples */
static xspect TM3015_ECS[IES_TM_30_15_ESAMPLES];	/* Forward declaration */

/* Comute IES TM-30-15 */
/* Return 1 on invalid, 2 on error */
/* Invalid is when sample is not white enough. */
int icx_IES_TM_30_15(
double *pRf,			/* Return Rf */
double *pRg,			/* Return Rg */
double *pcct,			/* Return correlated color temperature */
double *pdc,			/* Return signed dU'V' to locus */
double pbins[IES_TM_30_15_BINS][2][3],		/* If not NULL, return ref & tsamp Jab */
xspect *tsamp			/* Illuminant test sample to compute TLCI of */
) {
	// Test source = input test sample
	// Ref. source = reference Plankian/D type illuminant

	// Use 380 - 780 nm, 1nm spectral calculation

	// XYZ values are normalized to have Y == 1

	// Normalize plankian and D type to have Y == 100 before blending
	// and scale result so that Y == 100

	// Compute 10 degree XYZ values for all CES for test & ref. illuminants,
	// and normalize each to Y = 100.
	
	// Convert to CIECAM02 Jab using the fom detailed in 3.7.1 thru 3.7.23

	// Compute Rf according to 3.9.1 & 3.9.2
	// Compute Rg according to 3.10.1

	int i;
	double cct;			/* CCT of test sample */
	xsp2cie *tocie2;	/* spectral illuminant to XYZ conversion 2 degree */
	xsp2cie *tocie10;	/* spectral illuminant to XYZ conversion 10 degree */
	xsp2cie *tocie10r, *tocie10s;	/* spectral reflectance to XYZ conversion 10 degree */
	xspect rspec;		/* Reference white spectrum */
	double tiXYZ[3];	/* Sample illuminant XYZ */
	double riXYZ[3];	/* Reference illuminant XYZ */
	double tUCS[3];		/* Sample CIE 1960 UCS */
	double rUCS[3];		/* Reference CIE 1960 UCS */
	double tYxy[3], rYxy[3];
	double tsampnorm;
	double dc;			/* delta of test sample to reference white in 1960 UCS */

	double rcat[3][3], tcat[3][3];	/* Ref White chromatic adapation matrix */
	double rAw, tAw;	/* acromatic response of white */
	double mdE;			/* Mean delta E */
	double Rfd;

	double esamps[IES_TM_30_15_ESAMPLES][2][3];	/* Evaluation sample Jab values */

	int bcount[IES_TM_30_15_BINS];
	double _bins[IES_TM_30_15_BINS][2][3];
	double (*bins)[2][3];

	double rArea, tArea;

	double Rf = -1.0, Rg = -1.0;
	int rv = 0;

	DBGF((DBGA,"icx_IES_TM_30_15 called\n"))

	bins = (pbins != NULL) ? pbins : _bins;

	// Use CIE1931 2 deg for CCT, xy anxd u'v'
	/* Create spectral to XYZ for CCT etc. */
	if ((tocie2 = new_xsp2cie(icxIT_none, 0.0, NULL, icxOT_CIE_1931_2, NULL, icSigXYZData, 1)) == NULL) {
		//DBGF((DBGA,"Ref new_xsp2cie 2 degreefailed\n"))
		return 2;   
	}

#ifdef USE_5NM
#   pragma message("###### tm3015.c set to 5nm !!! ######")
	tocie2->set_int_steps(tocie2, 5.0, 380.0, 780.0);
#endif

	/* Compute the 2 degree XYZ of the test sample */
	tocie2->convert(tocie2, tiXYZ, tsamp);

	/* Find the standard 2 degree observer CIE 15:2004 Plankian CCT (exp 1.4388e-2) */
	if ((cct = icx_XYZ2ill_ct(NULL, icxIT_Ptemp, icxOT_CIE_1931_2, NULL, tiXYZ, NULL, 0)) < 0.0)
	{
		//DBGF((DBGA,"Ref icx_XYZ2ill_ct failed\n"))
		return 2;
	}
	DBGF((DBGA,"CCT = %f\n", cct))

	// For CCT <= 4500 generate plankian using 3.3.1
	// For CCT >= 5500 generate Daylight according to 3.3.3
	// For 4500 < CCT < 5500 then blend according to 3.3.9

	/* Hmm. Note we're using Appendix C improved Daylight, not the IES TM-30-15 3.3.3 */

	/* Create a reference white spectrum with the same CCT. */
	if (cct <= 4500.0) {
		if (standardIlluminant(&rspec, icxIT_Ptemp, cct)) {
			//DBGF((DBGA,"planckian_il failed\n"))
			return 2;
		}
	} else if (cct >= 5500.0) {
		if (standardIlluminant(&rspec, icxIT_Dtemp, cct)) {
			//DBGF((DBGA,"daylight_il failed\n"))
			return 2;
		}

	/* It's a linear blend between plankian(4500) and daylight(5500) */
	} else {
		double dwt;			/* Blend weight */
		xspect dspec;		/* Daylight spectrum */
		double pXYZ[3], dXYZ[3];		/* Plankian and Dalight Y values */

		if (standardIlluminant(&rspec, icxIT_Ptemp, 4500.0)) {
			//DBGF((DBGA,"planckian_il failed\n"))
			return 2;
		}
		tocie2->convert(tocie2, pXYZ, &rspec);

		DBGF((DBGA," plank Y %f\n", pXYZ[1]))
		
		if (standardIlluminant(&dspec, icxIT_Dtemp, 5500.0)) {
			//DBGF((DBGA,"daylight_il failed\n"))
			return 2;
		}
		tocie2->convert(tocie2, dXYZ, &dspec);

		DBGF((DBGA," D Y %f\n", dXYZ[1]))

		dwt = (cct - 4500)/(5500.0 - 4500.0);
		DBGF((DBGA,"creating hybrid spectrum with %f Plankian & %f Daylight\n",(1.0 - dwt),dwt))

		/* Blend Y normalized values */
		for (i = 0; i < rspec.spec_n; i++) {
			double wl = XSPECT_XWL(&rspec, i);			/* Wavelength in meters */
			rspec.spec[i] =        dwt  * 1.0/dXYZ[1] * value_xspect(&dspec, wl)
			              + (1.0 - dwt) * 1.0/pXYZ[1] * rspec.spec[i];
		}
	}

	/* Compute the 2 degree XYZ of the reference white */
	tocie2->convert(tocie2, riXYZ, &rspec);

	DBGF((DBGA,"ref. XYZ = %f %f %f\n",riXYZ[0],riXYZ[1],riXYZ[2]))
	DBGF((DBGA,"test XYZ = %f %f %f\n",tiXYZ[0],tiXYZ[1],tiXYZ[2]))

	/* Normalize the spectra so as to create a Y = 1.0 normalized white */
	rspec.norm *= riXYZ[1];
	tsampnorm = tsamp->norm;		/* Save this so we can restore it before returning */
	tsamp->norm *= tiXYZ[1];
	tocie2->convert(tocie2, riXYZ, &rspec);
	tocie2->convert(tocie2, tiXYZ, tsamp);

	DBGF((DBGA,"norm XYZ ref. = %f %f %f\n",riXYZ[0],riXYZ[1],riXYZ[2]))
	DBGF((DBGA,"norm XYZ test = %f %f %f\n",tiXYZ[0],tiXYZ[1],tiXYZ[2]))

	/* Convert to perceptual CIE 1960 UCS */
	icmXYZ21960UCS(rUCS, riXYZ);		/* 1960 UCS Yuv reference white */
	icmXYZ21960UCS(tUCS, tiXYZ);		/* 1960 UCS Yuv test sample white */

	DBGF((DBGA,"UCS ref. = %f %f %f\n",rUCS[0],rUCS[1],rUCS[2]))
	DBGF((DBGA,"UCS test = %f %f %f\n",tUCS[0],tUCS[1],tUCS[2]))

	/* Distance in UCS uv from test sample to reference illuminant */
	dc = icmNorm22(&rUCS[1], &tUCS[1]);

	DBGF((DBGA,"dc = %f (fraction of thresh. %f)\n",dc,dc/0.0054))
#ifdef DEBUG
	if (dc > 0.0054)
		DBGF((DBGA,"TM-30-15 is invalid\n")) 
#endif

	/* If dc > 0.0054 we should abort computing the CRI, */
	/* but this means we fail on lots of real world lighting. */
	if (dc > 0.0054)
		rv = 1;
		
	/* Comute dc sign */
	icmXYZ2Yxy(rYxy, riXYZ);
	icmXYZ2Yxy(tYxy, tiXYZ);

	/* Dot product of vector from aprox. locus curve "center" */
	/* xy 0.5, 0.25 with vector from	*/
	if ((tYxy[1] - rYxy[1]) * (rYxy[1] - 0.5)
	  + (tYxy[2] - rYxy[2]) * (rYxy[2] - 0.25) < 0.0)
		dc = -dc;
			
	// Use CIE1931 10 deg for comparison tristimulus values
	if ((tocie10 = new_xsp2cie(icxIT_none, 0.0, &rspec, icxOT_CIE_1964_10, NULL, icSigXYZData, 0)) == NULL) {
		//DBGF((DBGA,"Ref new_xsp2cie 10 degree failed\n"))
		return 2;   
	}

	/* Create spectral to XYZ with illuminants for delta E's */
	if ((tocie10r = new_xsp2cie(icxIT_custom, 0.0, &rspec, icxOT_CIE_1964_10, NULL, icSigXYZData, 0)) == NULL) {
		//DBGF((DBGA,"Ref new_xsp2cie 10 degree ref. failed\n"))
		return 2;   
	}
	if ((tocie10s = new_xsp2cie(icxIT_custom, 0.0, tsamp, icxOT_CIE_1964_10, NULL, icSigXYZData, 0)) == NULL) {
		//DBGF((DBGA,"Ref new_xsp2cie 10 degree ref. failed\n"))
		return 2;   
	}

#ifdef USE_5NM
# pragma message("###### tm3015.c set to 5nm !!! ######")
	tocie10->set_int_steps(tocie10, 5.0, 380.0, 780.0);
	tocie10r->set_int_steps(tocie10r, 5.0, 380.0, 780.0);
	tocie10s->set_int_steps(tocie10s, 5.0, 380.0, 780.0);
#endif

	/* Compute the 10 degree XYZ of the test sample and reference */
	rspec.norm = 1.0;
	tsamp->norm = 1.0;

	tocie10->convert(tocie10, tiXYZ, tsamp);
	tocie10->convert(tocie10, riXYZ, &rspec);

	/* Normalize the spectra so as to create a Y = 1.0 normalized white */
	rspec.norm *= riXYZ[1];
	tsamp->norm *= tiXYZ[1];

	tocie10->convert(tocie10, riXYZ, &rspec);
	tocie10->convert(tocie10, tiXYZ, tsamp);

	/* Compute CIECAM02 constants */
	rAw = comp_rgbd_mtx(rcat, riXYZ);
	tAw = comp_rgbd_mtx(tcat, tiXYZ);

	DBGF((DBGA,"rAw %f tAw\n",rAw, tAw))

	/* ------------------------------------------------------ */
	/* Compute delta Jab for each evaluation sample from reference white */
	mdE = 0.0;
	for (i = 0; i < IES_TM_30_15_ESAMPLES; i++) {
		double rXYZ[3];
		double sXYZ[3];
		double dE;

		tocie10r->convert(tocie10r, rXYZ, &TM3015_ECS[i]);

#ifdef TRACE_CAT
		apply_mtx_trace(esamps[i][0], rXYZ, riXYZ);
#endif /* TRACE_CAT */

		ciecam02_ucs(esamps[i][0], rAw, rcat, rXYZ);

		tocie10s->convert(tocie10s, sXYZ, &TM3015_ECS[i]);

#ifdef TRACE_CAT
		apply_mtx_trace(esamps[i][1], sXYZ, tiXYZ);
#endif /* TRACE_CAT */

		ciecam02_ucs(esamps[i][1], tAw, tcat, sXYZ);

		dE = icmNorm33(esamps[i][0], esamps[i][1]);

		DBGF((DBGA," Sample %d XYZ ref. %f %f %f test %f %f %f\n", i, rXYZ[0], rXYZ[1], rXYZ[2], sXYZ[0], sXYZ[1], sXYZ[2]))
		DBGF((DBGA,"          Jab ref. %f %f %f test %f %f %f dE %f\n", esamps[i][0][0], esamps[i][0][1], esamps[i][0][2], esamps[i][1][0], esamps[i][1][1], esamps[i][1][2], dE))

		mdE += dE;
	}

	/* Note that CIE244 uses a scaling factor of 6.73 rather than 7.54 */
	mdE /= (double)i;
	Rfd = 100.0 - 7.54 * mdE;
	if (Rfd < 0.0)
		Rfd = 0.0;
	Rf = 10.0 * log(exp(Rfd/10.0) + 1.0);

	/* ------------------------------------------------------ */
	/* Average evaluation sample Jab values into the hue bins */
	for (i = 0; i < IES_TM_30_15_BINS; i++) {
		icmSet3(bins[i][0], 0.0);
		icmSet3(bins[i][1], 0.0);
		bcount[i] = 0;
	}

	for (i = 0; i < IES_TM_30_15_ESAMPLES; i++) {
		double h;
		int bix;

		h = atan2(esamps[i][0][2], esamps[i][0][1]) / (2.0 * DBL_PI);	/* Reference ab */
		h = (h < 0.0) ? h + 1.0 : h;

		bix = (int)floor(IES_TM_30_15_BINS * h);
		if (bix > (IES_TM_30_15_BINS-1))		/* Just in case */
			bix = 0;
		
		DBGF((DBGA," esamp %d Jab ref. %f %f %f test %f %f %f Bin %d\n",i, esamps[i][0][0], esamps[i][0][1], esamps[i][0][2], esamps[i][1][0], esamps[i][1][1], esamps[i][1][2], bix))

		icmAdd3(bins[bix][0], bins[bix][0], esamps[i][0]);
		icmAdd3(bins[bix][1], bins[bix][1], esamps[i][1]);
		bcount[bix]++;

//DBGF((DBGA," Bin %d: ab ref. %f %f test %f %f cnt %d\n", bix, bins[bix][0][1], bins[bix][0][2], bins[bix][1][1], bins[bix][1][2], bcount[bix]))
	}

	for (i = 0; i < IES_TM_30_15_BINS; i++) {
		if (bcount[i] > 0) {
			icmScale3(bins[i][0], bins[i][0], 1.0/bcount[i]);
			icmScale3(bins[i][1], bins[i][1], 1.0/bcount[i]);
		}
		DBGF((DBGA," Bin %d: ab ref. %f %f test %f %f\n", i, bins[i][0][1], bins[i][0][2], bins[i][1][1], bins[i][1][2]))
	}

	/* Compute areas */
	rArea = tArea = 0.0;
	for (i = 0; i < IES_TM_30_15_BINS; i++) {
		int j = i < (IES_TM_30_15_BINS-1) ? i + 1 : 0;

		rArea += bins[i][0][1] * bins[j][0][2] - bins[j][0][1] * bins[i][0][2];
		tArea += bins[i][1][1] * bins[j][1][2] - bins[j][1][1] * bins[i][1][2];
	}
	rArea *= 0.5;
	tArea *= 0.5;

	DBGF((DBGA," Area: ref. %f test %f\n", rArea, tArea))

	Rg = 100.0 * tArea/rArea;

	tsamp->norm = tsampnorm;		/* Restore this */

	tocie2->del(tocie2);
	tocie2 = NULL;
	tocie10->del(tocie10);
	tocie10r->del(tocie10r);
	tocie10s->del(tocie10s);

	//DBGF((DBGA,"returning Rf %f Rg %f\n",Rf,Rg))

	if (pRf != NULL)
		*pRf = Rf;

	if (pRg != NULL)
		*pRg = Rg;

	if (pcct != NULL)
		*pcct = cct;

	if (pdc != NULL)
		*pdc = dc;

	return rv;
}

#else	/* STANDALONE_TEST */

static void tm3015_plot(double bins[IES_TM_30_15_BINS][2][3]);

int
main() {
	int rv;
	double Rf = -1.0, Rg = -1.0;
	double CCT = -1.0, duv = -1.0;
	double bins[IES_TM_30_15_BINS][2][3];

	/* Spreadsheet example */
	xspect f4 = { 81, 380.0, 780.0, 1.0, { 0.0188243065, 0.0231175694, 0.0287318362, 0.0323645971, 0.0663804491, 0.4540951123, 0.0643989432, 0.0525099075, 0.0581241744, 0.0637384412, 0.0693527081, 1, 0.2651915456, 0.0842140026, 0.0891677675, 0.0931307794, 0.0961030383, 0.0987450462, 0.1003963012, 0.1017173052, 0.1020475561, 0.1020475561, 0.1036988111, 0.1010568032, 0.0990752972, 0.0984147952, 0.0994055482, 0.1036988111, 0.1126155878, 0.1287978864, 0.1548877147, 0.1918758256, 0.2417437252, 0.7460369881, 0.499009247, 0.4583883752, 0.5392998679, 0.6169088507, 0.6816380449, 0.8018494055, 0.8672391017, 0.7688243065, 0.7575957728, 0.7311756935, 0.6905548217, 0.641677675, 0.5858652576, 0.5284015852, 0.4762219287, 0.4147952444, 0.3609643329, 0.3143989432, 0.2701453104, 0.2315059445, 0.1981505945, 0.1687582563, 0.143989432, 0.1218626156, 0.1033685601, 0.0871862616, 0.0739762219, 0.0630779392, 0.0561426684, 0.0459048877, 0.0389696169, 0.034015852, 0.0290620872, 0.0244385733, 0.0211360634, 0.0178335535, 0.0161822985, 0.0151915456, 0.0138705416, 0.0122192867, 0.0122192867, 0.0108982827, 0.0115587847, 0.0118890357, 0.0102377807, 0.0085865258, 0.0062747688 } };

	xspect f5 = {
	81, 380.0, 780.0,	/* 109 bands from 300 to 830 nm in 5nm steps */
	1.0,		/* Arbitrary scale factor */
	{
/* 380 */	1.87, 2.35, 2.92, 3.45, 5.10, 18.91, 6.00, 6.11, 6.85, 7.58,
/* 430 */	8.31, 40.76, 16.06, 10.32, 10.91, 11.40, 11.83, 12.17, 12.40, 12.54,
/* 480 */	12.58, 12.52, 12.47, 12.20, 11.89, 11.61, 11.33, 11.10, 10.96, 10.97,
/* 530 */	11.16, 11.54, 12.12, 27.78, 17.73, 14.47, 15.20, 15.77, 16.10, 18.54,
/* 580 */	19.50, 15.39, 14.64, 13.72, 12.69, 11.57, 10.45, 9.35, 8.29, 7.32,
/* 630 */	6.41, 5.63, 4.90, 4.26, 3.72, 3.25, 2.83, 2.49, 2.19, 1.93,
/* 680 */	1.71, 1.52, 1.48, 1.26, 1.13, 1.05, 0.96, 0.85, 0.78, 0.72,
/* 730 */	0.68, 0.67, 0.65, 0.61, 0.62, 0.59, 0.62, 0.64, 0.55, 0.47,
/* 780 */	0.40,
	}
};

/* CIE F8 */
/* Fluorescent, Wide band 5000K, CRI 95 */
	xspect f8 = {
	81, 380.0, 780.0,	/* 109 bands from 300 to 830 nm in 5nm steps */
	1.0,		/* Arbitrary scale factor */
	{
/* 380 */	1.21, 1.5, 1.81, 2.13, 3.17, 13.08, 3.83, 3.45, 3.86, 4.42,
/* 430 */	5.09, 34.10, 12.42, 7.68, 8.6, 9.46, 10.24, 10.84, 11.33, 11.71,
/* 480 */	11.98, 12.17, 12.28, 12.32, 12.35, 12.44, 12.55, 12.68, 12.77, 12.72,
/* 530 */	12.60, 12.43, 12.22, 28.96, 16.51, 11.79, 11.76, 11.77, 11.84, 14.61,
/* 580 */	16.11, 12.34, 12.53, 12.72, 12.92, 13.12, 13.34, 13.61, 13.87, 14.07,
/* 630 */	14.20, 14.16, 14.13, 14.34, 14.50, 14.46, 14.00, 12.58, 10.99, 9.98,
/* 680 */	9.22, 8.62, 8.07, 7.39, 6.71, 6.16, 5.63, 5.03, 4.46, 4.02,
/* 730 */	3.66, 3.36, 3.09, 2.85, 2.65, 2.51, 2.37, 2.15, 1.89, 1.61,
/* 780 */	1.32,
	}
};

	printf("Checking IES TM-30-15 CRI:\n");

	rv = icx_IES_TM_30_15(&Rf, &Rg, &CCT, &duv, bins, &f4);
	printf(" f4 returned %d Rf %f Rg %f CCT %f duv %f\n",rv,Rf,Rg,CCT,duv);
	printf("        expect Rf 51.5      Rg 83.5\n\n",rv,Rf,Rg,CCT,duv);
	tm3015_plot(bins);
	

	rv = icx_IES_TM_30_15(&Rf, &Rg, &CCT, &duv, bins, &f5);
	printf(" f5 returned %d Rf %f Rg %f CCT %f duv %f\n",rv,Rf,Rg,CCT,duv);
	printf("        expect Rf 75.2      Rg 87.4\n\n",rv,Rf,Rg,CCT,duv);
	tm3015_plot(bins);

	rv = icx_IES_TM_30_15(&Rf, &Rg, &CCT, &duv, bins, &f8);
	printf(" f8 returned %d Rf %f Rg %f CCT %f duv %f\n",rv,Rf,Rg,CCT,duv);
	printf("        expect Rf 95.9      Rg 101.0\n\n",rv,Rf,Rg,CCT,duv);
	tm3015_plot(bins);

	return 0;
}


#define SSAMP 4

// Note that bins is modified
static void tm3015_plot(double bins[IES_TM_30_15_BINS][2][3]) {
	double rvecs[SSAMP * 16][2];	// DEBUG ref light circle vectors
	double tvecs[SSAMP * 16][2];	// test light circle vectors
	double shvec[2 * 16][2];		// Shift vectors
	int i, j, k;
	double maxr = 0.0;
	plot_g gg = { 0 };
	float lblack[3] = { 0.5, 0.5, 0.5 };
	float lred[3]   = { 1.0, 0.5, 0.5 };
	float black[3] = { 0.0, 0.0, 0.0 };
	float red[3]   = { 1.0, 0.0, 0.0 };
	float blue[3] = { 0.2, 0.2, 1.0 };

#ifdef NEVER
	clear_g(&gg);

	for (i = 0; i < 16; i++) {
		int ip1 = i < 15 ? i+1 : 0; 
		
		add_vec_g(&gg, bins[i][0][1], bins[i][0][2], bins[ip1][0][1], bins[ip1][0][2], lblack);
		add_vec_g(&gg, bins[i][1][1], bins[i][1][2], bins[ip1][1][1], bins[ip1][1][2], lred);
	}
	do_plot_g(&gg, 0.0, 0.0, 0.0, 0.0, 1.0, 0, 1);
#endif

	clear_g(&gg);

	/* Copy raw shift vectors */
	for (i = 0; i < 16; i++) {
		shvec[2 * i + 0][0] = bins[i][0][1];
		shvec[2 * i + 0][1] = bins[i][0][2];
		shvec[2 * i + 1][0] = bins[i][1][1];
		shvec[2 * i + 1][1] = bins[i][1][2];
	}

	// Convert to angle & radius 
	for (i = 0; i < 16; i++) {
		for (j = 0; j < 2; j++) {
			double x = bins[i][j][1];
			double y = bins[i][j][2];

			// atan2 returns -pi < val <= +pi
			bins[i][j][1] = atan2(y, x);
			bins[i][j][2] = sqrt(x * x + y * y);
		}
	}

	// Interpolate test points, and normalize them by interp. reference points */
	for (i = 0; i < 16; i++) {
		int ip1 = i < 15 ? i+1 : 0; 
		int k = i, kp1 = ip1;
		for (j = 0; j < SSAMP; j++) {	// 4 x interpolation
			double ra, rr, ta, tr;
			double ra0, ra1, ta0, ta1;
			double bl = (double)j/SSAMP;
			double nf = 1.0;

//printf("i %d j %d\n",i,j);

			ra0 = bins[i][0][1];
			ra1 = bins[ip1][0][1];
			ta0 = bins[i][1][1];
			ta1 = bins[ip1][1][1];

			// Make sure pairs are adjacent
			if ((ra1 - ra0) > (1.0 * DBL_PI))
				ra0 += 2.0 * DBL_PI;
			else if ((ra1 - ra0) < -(1.0 * DBL_PI))
				ra0 -= 2.0 * DBL_PI;

			if ((ta1 - ta0) > (1.0 * DBL_PI))
				ta0 += 2.0 * DBL_PI;
			else if ((ta1 - ta0) < -(1.0 * DBL_PI))
				ta0 -= 2.0 * DBL_PI;

			// Interpolated values
			ra = (1.0 - bl) * ra0 + bl * ra1;
			rr = (1.0 - bl) * bins[i][0][2] + bl * bins[ip1][0][2];
			ta = (1.0 - bl) * ta0 + bl * ta1;
			tr = (1.0 - bl) * bins[i][1][2] + bl * bins[ip1][1][2];

//printf(" bl %f ar %f rr %f tr %f rr %f\n",bl, ra, rr, tr, rr);

			nf = rr;		// Normalization factor

			tr /= nf;
			rr /= nf;

			if (tr > maxr)
				maxr = tr;

//printf(" norm ta %f tr %f\n",ta,tr);
			rvecs[SSAMP * i + j][0] = ra;
			rvecs[SSAMP * i + j][1] = rr; 

			tvecs[SSAMP * i + j][0] = ta;
			tvecs[SSAMP * i + j][1] = tr; 

			// Set shift vectors
			if (j == 0) {
				shvec[2 * i + 0][0] = ra;
				shvec[2 * i + 0][1] = rr;
				shvec[2 * i + 1][0] = ta; 
				shvec[2 * i + 1][1] = tr; 
//printf(" Shift ra %f rr %f ta %f tr %f\n",bins[i][0][1],1.0,tr,ta);
			}
//printf("\n");
		}
	}

//printf("Done all, convert back\n");

	// Convert back to cartesian
	for (i = 0; i < (SSAMP * 16); i++) {
		double a = tvecs[i][0];
		double r = tvecs[i][1];

		tvecs[i][0] = r * cos(a); 
		tvecs[i][1] = r * sin(a);
//printf(" tvecs[%d] a %f r %f x %f y %f\n",i,a,r,tvecs[i][0],tvecs[i][1]);
	}

	for (i = 0; i < (SSAMP * 16); i++) {
		double a = rvecs[i][0];
		double r = rvecs[i][1];

		rvecs[i][0] = r * cos(a); 
		rvecs[i][1] = r * sin(a);
//printf(" rvecs[%d] a %f r %f x %f y %f\n",i,a,r,tvecs[i][0],tvecs[i][1]);
	}

	for (i = 0; i < (2 * 16); i++) {
		double a = shvec[i][0];
		double r = shvec[i][1];

		shvec[i][0] = r * cos(a); 
		shvec[i][1] = r * sin(a);
//printf(" shvec[%d] a %f r %f x %f y %f\n",i,a,r,shvec[i][0],shvec[i][1]);
	}

//printf("About to plot\n");

	for (i = 0; i < (SSAMP * 16); i++) {
		int ip1 = i < (SSAMP * 16 -1) ? i+1 : 0; 
		double x0, y0, x1, y1, r0, r1;

		x0 = tvecs[i][0];
		y0 = tvecs[i][1];
		x1 = tvecs[ip1][0];
		y1 = tvecs[ip1][1];

		r0 = sqrt(x0 * x0 + y0 * y0);
		r1 = sqrt(x1 * x1 + y1 * y1);

		add_vec_g(&gg, x0/r0, y0/r0, x1/r1, y1/r1, black);
//		add_vec_g(&gg, rvecs[i][0], rvecs[i][1], rvecs[ip1][0], rvecs[ip1][1], black);

		add_vec_g(&gg, x0, y0, x1, y1, red);

	}

	for (i = 0; i < 16; i++) {
		add_vec_g(&gg, shvec[2 * i + 0][0], shvec[2 * i + 0][1],
		               shvec[2 * i + 1][0], shvec[2 * i + 1][1], blue);
	}

	do_plot_g(&gg, 0.0, 0.0, 0.0, 0.0, 1.0, 0, 1);
	
}

#undef SSAMP

#endif /* STANDALONE_TEST */

/* ------------------------------------------------------ */

#ifndef STANDALONE_TEST

/* IES TM-30-15 evaluation samples */
/* 99 samples, 360 .. 830 at 5nm intervals */
static xspect TM3015_ECS[IES_TM_30_15_ESAMPLES] = {
	{ 95, 360.0, 830.0, 1.0, { 0.61947, 0.61947, 0.61947, 0.61947, 0.61947, 0.62357, 0.62766, 0.63172, 0.63590, 0.64015, 0.64540, 0.64888, 0.65070, 0.65147, 0.65080, 0.64823, 0.64400, 0.63825, 0.63010, 0.61895, 0.60630, 0.59361, 0.58020, 0.56548, 0.55150, 0.54016, 0.53000, 0.51942, 0.50950, 0.50163, 0.49580, 0.49171, 0.48940, 0.48884, 0.48940, 0.49068, 0.49380, 0.50011, 0.51030, 0.52488, 0.54420, 0.56832, 0.59630, 0.62679, 0.65780, 0.68743, 0.71490, 0.73974, 0.76170, 0.78072, 0.79720, 0.81156, 0.82390, 0.83427, 0.84280, 0.84971, 0.85560, 0.86098, 0.86570, 0.86955, 0.87280, 0.87578, 0.87850, 0.88094, 0.88330, 0.88579, 0.88840, 0.89091, 0.89260, 0.89420, 0.89573, 0.89723, 0.89872, 0.90019, 0.90164, 0.90307, 0.90449, 0.90588, 0.90726, 0.90861, 0.90995, 0.91128, 0.91258, 0.91387, 0.91514, 0.91514, 0.91514, 0.91514, 0.91514, 0.91514, 0.91514, 0.91514, 0.91514, 0.91514, 0.91514 } },
	{ 95, 360.0, 830.0, 1.0, { 0.25173, 0.25173, 0.25173, 0.25173, 0.25173, 0.25413, 0.25655, 0.25898, 0.26150, 0.26409, 0.26700, 0.26834, 0.26840, 0.26779, 0.26630, 0.26349, 0.25900, 0.25241, 0.24310, 0.23086, 0.21750, 0.20471, 0.19180, 0.17813, 0.16560, 0.15602, 0.14810, 0.14026, 0.13270, 0.12604, 0.12090, 0.11769, 0.11610, 0.11569, 0.11610, 0.11710, 0.11870, 0.12114, 0.12520, 0.13198, 0.14330, 0.16092, 0.18550, 0.21707, 0.25420, 0.29505, 0.33760, 0.37991, 0.42050, 0.45828, 0.49320, 0.52533, 0.55400, 0.57872, 0.60030, 0.61963, 0.63670, 0.65151, 0.66510, 0.67850, 0.69160, 0.70422, 0.71710, 0.73050, 0.74200, 0.74976, 0.75710, 0.76684, 0.77460, 0.78201, 0.78903, 0.79587, 0.80256, 0.80907, 0.81542, 0.82161, 0.82763, 0.83349, 0.83919, 0.84473, 0.85011, 0.85534, 0.86042, 0.86534, 0.87012, 0.87012, 0.87012, 0.87012, 0.87012, 0.87012, 0.87012, 0.87012, 0.87012, 0.87012, 0.87012 } },
	{ 95, 360.0, 830.0, 1.0, { 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00663, 0.01694, 0.01773, 0.01167, 0.00511, 0.00379, 0.00874, 0.01888, 0.02972, 0.04199, 0.04896, 0.05457, 0.05674, 0.05410, 0.05487, 0.05908, 0.06075, 0.06259, 0.06304, 0.06069, 0.05967, 0.05944, 0.05601, 0.05392, 0.05202, 0.04895, 0.04762, 0.04647, 0.04398, 0.04102, 0.04046, 0.04187, 0.04197, 0.04123, 0.04185, 0.04305, 0.04484, 0.04429, 0.04632, 0.04914, 0.05151, 0.05305, 0.05463, 0.05863, 0.06604, 0.07427, 0.08516, 0.09790, 0.10928, 0.12064, 0.12749, 0.13322, 0.14094, 0.15003, 0.15656, 0.16207, 0.17095, 0.18004, 0.18664, 0.18611, 0.18650, 0.19529, 0.20210, 0.20393, 0.21328, 0.21966, 0.22617, 0.23281, 0.23959, 0.24651, 0.25356, 0.26073, 0.26804, 0.27548, 0.28305, 0.29074, 0.29855, 0.30648, 0.31452, 0.32268, 0.32268, 0.32268, 0.32268, 0.32268, 0.32268, 0.32268, 0.32268, 0.32268, 0.32268, 0.32268 } },
	{ 95, 360.0, 830.0, 1.0, { 0.47376, 0.47376, 0.47376, 0.47376, 0.47376, 0.47793, 0.46406, 0.44102, 0.42157, 0.40559, 0.38950, 0.37567, 0.36501, 0.35611, 0.34936, 0.34200, 0.33382, 0.32411, 0.31282, 0.29816, 0.27982, 0.26042, 0.24078, 0.22206, 0.20462, 0.18884, 0.17488, 0.16337, 0.15560, 0.15190, 0.15143, 0.15412, 0.16036, 0.17075, 0.18586, 0.20558, 0.22892, 0.25477, 0.28272, 0.31202, 0.34200, 0.37116, 0.39782, 0.42123, 0.44136, 0.45795, 0.47105, 0.48107, 0.48842, 0.49386, 0.49824, 0.50229, 0.50567, 0.50819, 0.51009, 0.51142, 0.51262, 0.51402, 0.51522, 0.51566, 0.51527, 0.51410, 0.51322, 0.51343, 0.51337, 0.51311, 0.51261, 0.51208, 0.51191, 0.51291, 0.51485, 0.51667, 0.51849, 0.52041, 0.52228, 0.52415, 0.52602, 0.52789, 0.52976, 0.53163, 0.53350, 0.53537, 0.53723, 0.53910, 0.54096, 0.54096, 0.54096, 0.54096, 0.54096, 0.54096, 0.54096, 0.54096, 0.54096, 0.54096, 0.54096 } },
	{ 95, 360.0, 830.0, 1.0, { 0.13351, 0.13351, 0.13351, 0.13351, 0.13351, 0.13540, 0.13565, 0.13540, 0.13579, 0.13757, 0.14001, 0.14195, 0.14194, 0.13887, 0.13327, 0.12595, 0.11738, 0.10794, 0.09819, 0.08868, 0.07978, 0.07182, 0.06508, 0.05976, 0.05586, 0.05323, 0.05140, 0.04991, 0.04869, 0.04783, 0.04761, 0.04831, 0.05011, 0.05302, 0.05654, 0.06008, 0.06335, 0.06628, 0.06946, 0.07376, 0.08059, 0.09115, 0.10539, 0.12247, 0.13976, 0.15471, 0.16698, 0.17763, 0.19118, 0.21263, 0.24557, 0.29221, 0.35063, 0.41726, 0.48613, 0.55140, 0.61007, 0.66051, 0.70370, 0.74105, 0.77303, 0.79995, 0.82248, 0.84133, 0.85712, 0.87050, 0.88222, 0.89292, 0.90256, 0.91091, 0.91776, 0.92304, 0.92737, 0.93156, 0.93639, 0.93985, 0.94353, 0.94700, 0.95027, 0.95335, 0.95625, 0.95897, 0.96153, 0.96394, 0.96620, 0.96620, 0.96620, 0.96620, 0.96620, 0.96620, 0.96620, 0.96620, 0.96620, 0.96620, 0.96620 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04525, 0.04525, 0.04525, 0.04525, 0.04525, 0.05337, 0.06285, 0.07388, 0.08249, 0.10382, 0.11965, 0.13053, 0.13699, 0.13975, 0.14015, 0.13946, 0.13788, 0.13541, 0.13217, 0.12835, 0.12428, 0.12021, 0.11600, 0.11157, 0.10743, 0.10413, 0.10171, 0.09999, 0.09856, 0.09712, 0.09600, 0.09563, 0.09619, 0.09777, 0.10033, 0.10388, 0.10861, 0.11472, 0.12221, 0.13104, 0.14123, 0.15275, 0.16518, 0.17796, 0.19031, 0.20157, 0.21160, 0.22030, 0.22708, 0.23139, 0.23328, 0.23318, 0.23230, 0.23215, 0.23457, 0.24132, 0.25369, 0.27287, 0.30030, 0.33695, 0.38151, 0.43149, 0.48165, 0.52720, 0.56789, 0.60430, 0.63589, 0.66181, 0.68123, 0.70786, 0.73017, 0.75138, 0.77145, 0.79034, 0.80806, 0.82462, 0.84003, 0.85433, 0.86754, 0.87973, 0.89094, 0.90122, 0.91063, 0.91923, 0.92706, 0.92706, 0.92706, 0.92706, 0.92706, 0.92706, 0.92706, 0.92706, 0.92706, 0.92706, 0.92706 } },
	{ 95, 360.0, 830.0, 1.0, { 0.06001, 0.06001, 0.06001, 0.06001, 0.06001, 0.05958, 0.05916, 0.05873, 0.05830, 0.05786, 0.05730, 0.05692, 0.05680, 0.05694, 0.05720, 0.05744, 0.05760, 0.05763, 0.05730, 0.05646, 0.05530, 0.05406, 0.05290, 0.05192, 0.05110, 0.05039, 0.04980, 0.04936, 0.04900, 0.04866, 0.04840, 0.04828, 0.04820, 0.04807, 0.04790, 0.04775, 0.04770, 0.04780, 0.04800, 0.04829, 0.04880, 0.04973, 0.05130, 0.05390, 0.05850, 0.06636, 0.07920, 0.09894, 0.12770, 0.16683, 0.21430, 0.26677, 0.31900, 0.36580, 0.40400, 0.43172, 0.45030, 0.46188, 0.46860, 0.47230, 0.47360, 0.47296, 0.47150, 0.47030, 0.46950, 0.46891, 0.46790, 0.46607, 0.46450, 0.46297, 0.46148, 0.45999, 0.45850, 0.45702, 0.45553, 0.45404, 0.45256, 0.45108, 0.44959, 0.44811, 0.44663, 0.44515, 0.44367, 0.44220, 0.44072, 0.44072, 0.44072, 0.44072, 0.44072, 0.44072, 0.44072, 0.44072, 0.44072, 0.44072, 0.44072 } },
	{ 95, 360.0, 830.0, 1.0, { 0.05818, 0.05818, 0.05818, 0.05818, 0.05818, 0.05961, 0.05984, 0.05941, 0.05886, 0.05858, 0.05836, 0.05797, 0.05769, 0.05777, 0.05795, 0.05790, 0.05764, 0.05729, 0.05707, 0.05715, 0.05725, 0.05711, 0.05685, 0.05663, 0.05647, 0.05633, 0.05609, 0.05569, 0.05529, 0.05505, 0.05496, 0.05495, 0.05491, 0.05477, 0.05459, 0.05446, 0.05441, 0.05442, 0.05444, 0.05444, 0.05458, 0.05501, 0.05571, 0.05682, 0.05915, 0.06379, 0.07216, 0.08577, 0.10608, 0.13417, 0.16962, 0.21130, 0.25684, 0.30366, 0.34959, 0.39266, 0.43125, 0.46417, 0.49160, 0.51406, 0.53220, 0.54676, 0.55870, 0.56897, 0.57812, 0.58662, 0.59483, 0.60303, 0.61109, 0.61878, 0.62574, 0.63173, 0.63709, 0.64228, 0.64779, 0.65285, 0.65804, 0.66319, 0.66830, 0.67337, 0.67840, 0.68339, 0.68834, 0.69325, 0.69811, 0.69811, 0.69811, 0.69811, 0.69811, 0.69811, 0.69811, 0.69811, 0.69811, 0.69811, 0.69811 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04157, 0.04157, 0.04157, 0.04157, 0.04157, 0.04172, 0.04187, 0.04202, 0.04225, 0.04229, 0.04246, 0.04275, 0.04318, 0.04374, 0.04441, 0.04516, 0.04596, 0.04681, 0.04771, 0.04862, 0.04936, 0.04981, 0.05008, 0.05044, 0.05142, 0.05340, 0.05585, 0.05805, 0.05925, 0.05899, 0.05771, 0.05602, 0.05420, 0.05240, 0.05060, 0.04872, 0.04668, 0.04446, 0.04225, 0.04038, 0.03936, 0.03960, 0.04081, 0.04258, 0.04483, 0.04750, 0.05039, 0.05323, 0.05554, 0.05696, 0.05771, 0.05813, 0.05843, 0.05870, 0.05894, 0.05917, 0.05966, 0.06065, 0.06193, 0.06320, 0.06420, 0.06474, 0.06492, 0.06489, 0.06471, 0.06441, 0.06399, 0.06344, 0.06276, 0.06228, 0.06173, 0.06118, 0.06064, 0.06010, 0.05957, 0.05904, 0.05852, 0.05800, 0.05748, 0.05697, 0.05646, 0.05596, 0.05546, 0.05496, 0.05447, 0.05447, 0.05447, 0.05447, 0.05447, 0.05447, 0.05447, 0.05447, 0.05447, 0.05447, 0.05447 } },
	{ 95, 360.0, 830.0, 1.0, { 0.17602, 0.17602, 0.17602, 0.17602, 0.17602, 0.18172, 0.18756, 0.19354, 0.19983, 0.20628, 0.21152, 0.20944, 0.20352, 0.19850, 0.19473, 0.19149, 0.18842, 0.18529, 0.18231, 0.17978, 0.17804, 0.17718, 0.17637, 0.17467, 0.17192, 0.16834, 0.16488, 0.16250, 0.16135, 0.16120, 0.16111, 0.16047, 0.16075, 0.16357, 0.16904, 0.17689, 0.18673, 0.19850, 0.21355, 0.23345, 0.25938, 0.29341, 0.34157, 0.40830, 0.48749, 0.57053, 0.64908, 0.71596, 0.76830, 0.80508, 0.82849, 0.84185, 0.84979, 0.85637, 0.86217, 0.86695, 0.87064, 0.87337, 0.87575, 0.87822, 0.88011, 0.88074, 0.88075, 0.88086, 0.88083, 0.88044, 0.88050, 0.88167, 0.88291, 0.88412, 0.88529, 0.88645, 0.88760, 0.88874, 0.88987, 0.89099, 0.89210, 0.89320, 0.89429, 0.89537, 0.89644, 0.89750, 0.89855, 0.89959, 0.90062, 0.90062, 0.90062, 0.90062, 0.90062, 0.90062, 0.90062, 0.90062, 0.90062, 0.90062, 0.90062 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03951, 0.03951, 0.03951, 0.03951, 0.03951, 0.04449, 0.05006, 0.05629, 0.06081, 0.07233, 0.08032, 0.08522, 0.08742, 0.08748, 0.08644, 0.08520, 0.08358, 0.08121, 0.07806, 0.07423, 0.06998, 0.06554, 0.06091, 0.05612, 0.05174, 0.04833, 0.04593, 0.04441, 0.04337, 0.04243, 0.04169, 0.04131, 0.04139, 0.04207, 0.04366, 0.04679, 0.05302, 0.06419, 0.08220, 0.10861, 0.14340, 0.18552, 0.23131, 0.27651, 0.31696, 0.34933, 0.37353, 0.39051, 0.40201, 0.40996, 0.41601, 0.42172, 0.42833, 0.43684, 0.44755, 0.46045, 0.47495, 0.49066, 0.50855, 0.52962, 0.55369, 0.57982, 0.60514, 0.62726, 0.64752, 0.66760, 0.68675, 0.70360, 0.71681, 0.73386, 0.74863, 0.76285, 0.77650, 0.78958, 0.80209, 0.81403, 0.82541, 0.83623, 0.84651, 0.85625, 0.86547, 0.87419, 0.88242, 0.89018, 0.89748, 0.89748, 0.89748, 0.89748, 0.89748, 0.89748, 0.89748, 0.89748, 0.89748, 0.89748, 0.89748 } },
	{ 95, 360.0, 830.0, 1.0, { 0.24167, 0.24167, 0.24167, 0.24167, 0.24167, 0.21622, 0.18687, 0.16204, 0.14514, 0.12981, 0.11536, 0.10584, 0.10430, 0.09871, 0.09483, 0.09767, 0.09829, 0.09474, 0.09395, 0.09191, 0.08603, 0.08133, 0.07541, 0.06828, 0.06159, 0.05719, 0.05416, 0.05087, 0.04867, 0.04995, 0.05235, 0.05565, 0.05980, 0.06338, 0.07119, 0.08203, 0.09658, 0.11374, 0.13250, 0.15516, 0.18331, 0.21658, 0.25334, 0.29135, 0.33392, 0.37966, 0.42691, 0.47137, 0.51260, 0.54874, 0.57960, 0.60302, 0.62230, 0.63769, 0.64759, 0.65705, 0.66606, 0.67239, 0.67883, 0.67562, 0.67555, 0.67885, 0.67821, 0.67843, 0.68161, 0.68412, 0.68273, 0.68325, 0.68676, 0.68639, 0.68745, 0.68851, 0.68957, 0.69063, 0.69169, 0.69275, 0.69380, 0.69485, 0.69590, 0.69695, 0.69799, 0.69903, 0.70007, 0.70111, 0.70215, 0.70215, 0.70215, 0.70215, 0.70215, 0.70215, 0.70215, 0.70215, 0.70215, 0.70215, 0.70215 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04337, 0.04337, 0.04337, 0.04337, 0.04337, 0.04333, 0.04328, 0.04324, 0.04327, 0.04311, 0.04307, 0.04320, 0.04357, 0.04419, 0.04496, 0.04572, 0.04625, 0.04641, 0.04644, 0.04655, 0.04635, 0.04557, 0.04496, 0.04556, 0.04863, 0.05497, 0.06332, 0.07191, 0.07909, 0.08373, 0.08654, 0.08841, 0.08912, 0.08836, 0.08654, 0.08414, 0.08118, 0.07767, 0.07423, 0.07149, 0.06947, 0.06819, 0.06818, 0.07069, 0.07949, 0.09754, 0.12216, 0.14933, 0.17516, 0.19672, 0.21476, 0.23037, 0.24254, 0.25030, 0.25485, 0.25781, 0.26021, 0.26266, 0.26487, 0.26641, 0.26735, 0.26785, 0.26795, 0.26770, 0.26735, 0.26713, 0.26696, 0.26668, 0.26616, 0.26598, 0.26566, 0.26535, 0.26503, 0.26471, 0.26440, 0.26408, 0.26376, 0.26345, 0.26313, 0.26282, 0.26250, 0.26219, 0.26187, 0.26156, 0.26124, 0.26124, 0.26124, 0.26124, 0.26124, 0.26124, 0.26124, 0.26124, 0.26124, 0.26124, 0.26124 } },
	{ 95, 360.0, 830.0, 1.0, { 0.33532, 0.33532, 0.33532, 0.33532, 0.33532, 0.33683, 0.33833, 0.33984, 0.34140, 0.34301, 0.34510, 0.34676, 0.34820, 0.34976, 0.35130, 0.35261, 0.35380, 0.35499, 0.35620, 0.35738, 0.35850, 0.35955, 0.36070, 0.36205, 0.36340, 0.36454, 0.36560, 0.36676, 0.36820, 0.37009, 0.37270, 0.37625, 0.38080, 0.38614, 0.39110, 0.39457, 0.39660, 0.39753, 0.39760, 0.39728, 0.39790, 0.40103, 0.40830, 0.42058, 0.43550, 0.45009, 0.46220, 0.47043, 0.47560, 0.47895, 0.48120, 0.48286, 0.48410, 0.48498, 0.48560, 0.48607, 0.48670, 0.48768, 0.48870, 0.48941, 0.48980, 0.48999, 0.49030, 0.49098, 0.49190, 0.49283, 0.49370, 0.49444, 0.49490, 0.49534, 0.49577, 0.49619, 0.49662, 0.49704, 0.49747, 0.49789, 0.49832, 0.49874, 0.49917, 0.49959, 0.50002, 0.50044, 0.50087, 0.50129, 0.50172, 0.50172, 0.50172, 0.50172, 0.50172, 0.50172, 0.50172, 0.50172, 0.50172, 0.50172, 0.50172 } },
	{ 95, 360.0, 830.0, 1.0, { 0.16992, 0.16992, 0.16992, 0.16992, 0.16992, 0.17295, 0.17602, 0.17913, 0.18141, 0.18593, 0.18883, 0.19137, 0.19479, 0.20029, 0.20890, 0.22103, 0.23492, 0.24851, 0.26083, 0.27135, 0.28027, 0.28799, 0.29480, 0.30093, 0.30660, 0.31208, 0.31800, 0.32476, 0.33168, 0.33767, 0.34119, 0.34098, 0.33712, 0.33041, 0.32312, 0.31756, 0.31476, 0.31513, 0.31810, 0.32241, 0.32510, 0.32407, 0.32228, 0.32449, 0.33775, 0.36682, 0.40525, 0.44457, 0.47934, 0.50614, 0.52637, 0.54227, 0.55469, 0.56419, 0.57172, 0.57820, 0.58395, 0.58911, 0.59367, 0.59765, 0.60150, 0.60563, 0.60997, 0.61431, 0.61843, 0.62214, 0.62533, 0.62795, 0.62993, 0.63284, 0.63536, 0.63788, 0.64039, 0.64289, 0.64538, 0.64787, 0.65035, 0.65282, 0.65528, 0.65773, 0.66018, 0.66262, 0.66504, 0.66746, 0.66987, 0.66987, 0.66987, 0.66987, 0.66987, 0.66987, 0.66987, 0.66987, 0.66987, 0.66987, 0.66987 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03466, 0.03466, 0.03466, 0.03466, 0.03466, 0.03537, 0.03609, 0.03683, 0.03767, 0.03829, 0.03910, 0.04003, 0.04104, 0.04212, 0.04349, 0.04531, 0.04737, 0.04944, 0.05166, 0.05419, 0.05696, 0.05989, 0.06299, 0.06628, 0.06973, 0.07323, 0.07667, 0.07996, 0.08320, 0.08647, 0.08973, 0.09289, 0.09576, 0.09824, 0.10056, 0.10293, 0.10535, 0.10777, 0.11015, 0.11250, 0.11485, 0.11730, 0.12026, 0.12415, 0.12934, 0.13622, 0.14527, 0.15696, 0.17140, 0.18852, 0.20764, 0.22810, 0.24970, 0.27223, 0.29483, 0.31669, 0.33780, 0.35841, 0.37905, 0.40003, 0.42039, 0.43917, 0.45673, 0.47371, 0.49042, 0.50703, 0.52319, 0.53846, 0.55239, 0.56804, 0.58283, 0.59748, 0.61196, 0.62624, 0.64030, 0.65413, 0.66770, 0.68100, 0.69401, 0.70672, 0.71912, 0.73119, 0.74292, 0.75432, 0.76537, 0.76537, 0.76537, 0.76537, 0.76537, 0.76537, 0.76537, 0.76537, 0.76537, 0.76537, 0.76537 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03080, 0.03080, 0.03080, 0.03080, 0.03080, 0.03310, 0.03556, 0.03820, 0.04031, 0.04443, 0.04751, 0.04967, 0.05105, 0.05186, 0.05253, 0.05346, 0.05470, 0.05623, 0.05805, 0.06016, 0.06249, 0.06492, 0.06732, 0.06952, 0.07145, 0.07306, 0.07441, 0.07559, 0.07648, 0.07701, 0.07747, 0.07817, 0.07924, 0.08086, 0.08358, 0.08811, 0.09531, 0.10580, 0.11916, 0.13450, 0.15040, 0.16541, 0.17859, 0.18928, 0.19761, 0.20396, 0.20914, 0.21381, 0.21771, 0.22041, 0.22166, 0.22146, 0.22077, 0.22081, 0.22303, 0.22881, 0.23900, 0.25451, 0.27704, 0.30800, 0.34672, 0.39138, 0.43730, 0.48016, 0.51979, 0.55666, 0.58977, 0.61775, 0.63924, 0.66823, 0.69292, 0.71656, 0.73905, 0.76037, 0.78045, 0.79930, 0.81691, 0.83330, 0.84849, 0.86252, 0.87545, 0.88732, 0.89819, 0.90812, 0.91717, 0.91717, 0.91717, 0.91717, 0.91717, 0.91717, 0.91717, 0.91717, 0.91717, 0.91717, 0.91717 } },
	{ 95, 360.0, 830.0, 1.0, { 0.06498, 0.06498, 0.06498, 0.06498, 0.06498, 0.06692, 0.06890, 0.07094, 0.07310, 0.07538, 0.07850, 0.08130, 0.08410, 0.08752, 0.09180, 0.09696, 0.10250, 0.10789, 0.11290, 0.11747, 0.12180, 0.12611, 0.13040, 0.13466, 0.13900, 0.14355, 0.14840, 0.15359, 0.15900, 0.16433, 0.16880, 0.17166, 0.17280, 0.17246, 0.17160, 0.17132, 0.17250, 0.17559, 0.17960, 0.18329, 0.18580, 0.18685, 0.18810, 0.19170, 0.19990, 0.21409, 0.23220, 0.25140, 0.26940, 0.28450, 0.29690, 0.30720, 0.31570, 0.32272, 0.32890, 0.33486, 0.34080, 0.34674, 0.35240, 0.35758, 0.36260, 0.36783, 0.37320, 0.37853, 0.38360, 0.38829, 0.39300, 0.39783, 0.40120, 0.40448, 0.40765, 0.41084, 0.41403, 0.41723, 0.42043, 0.42364, 0.42686, 0.43008, 0.43332, 0.43655, 0.43979, 0.44304, 0.44629, 0.44955, 0.45281, 0.45281, 0.45281, 0.45281, 0.45281, 0.45281, 0.45281, 0.45281, 0.45281, 0.45281, 0.45281 } },
	{ 95, 360.0, 830.0, 1.0, { 0.10927, 0.10927, 0.10927, 0.10927, 0.10927, 0.11190, 0.11458, 0.11733, 0.12020, 0.12316, 0.12590, 0.12589, 0.12450, 0.12359, 0.12320, 0.12298, 0.12300, 0.12334, 0.12380, 0.12417, 0.12450, 0.12484, 0.12510, 0.12518, 0.12510, 0.12501, 0.12560, 0.12801, 0.13480, 0.14805, 0.16640, 0.18711, 0.20520, 0.21662, 0.22330, 0.22905, 0.23930, 0.25936, 0.29250, 0.33928, 0.39140, 0.43925, 0.47690, 0.50097, 0.51460, 0.52210, 0.52600, 0.52810, 0.52900, 0.52902, 0.52860, 0.52810, 0.52750, 0.52676, 0.52600, 0.52534, 0.52460, 0.52361, 0.52260, 0.52181, 0.52110, 0.52029, 0.51950, 0.51887, 0.51820, 0.51729, 0.51620, 0.51512, 0.51440, 0.51370, 0.51303, 0.51235, 0.51168, 0.51100, 0.51033, 0.50965, 0.50897, 0.50830, 0.50762, 0.50695, 0.50627, 0.50560, 0.50492, 0.50425, 0.50357, 0.50357, 0.50357, 0.50357, 0.50357, 0.50357, 0.50357, 0.50357, 0.50357, 0.50357, 0.50357 } },
	{ 95, 360.0, 830.0, 1.0, { 0.02509, 0.02509, 0.02509, 0.02509, 0.02509, 0.02701, 0.02908, 0.03129, 0.03197, 0.03714, 0.03952, 0.03916, 0.03614, 0.03094, 0.02572, 0.02242, 0.02035, 0.01859, 0.01787, 0.01943, 0.02472, 0.03577, 0.05650, 0.08996, 0.13394, 0.18288, 0.22320, 0.24320, 0.24664, 0.24070, 0.23085, 0.22122, 0.21248, 0.20455, 0.19769, 0.19234, 0.18925, 0.18913, 0.19252, 0.19974, 0.21059, 0.22470, 0.24147, 0.26062, 0.28317, 0.31075, 0.34583, 0.39045, 0.44412, 0.50458, 0.56506, 0.61852, 0.66127, 0.69166, 0.71270, 0.72818, 0.74021, 0.75035, 0.75977, 0.76939, 0.77972, 0.79098, 0.80256, 0.81379, 0.82440, 0.83423, 0.84297, 0.85030, 0.85588, 0.86331, 0.86964, 0.87572, 0.88156, 0.88715, 0.89252, 0.89765, 0.90257, 0.90728, 0.91178, 0.91609, 0.92020, 0.92413, 0.92788, 0.93146, 0.93487, 0.93487, 0.93487, 0.93487, 0.93487, 0.93487, 0.93487, 0.93487, 0.93487, 0.93487, 0.93487 } },
	{ 95, 360.0, 830.0, 1.0, { 0.12356, 0.12356, 0.12356, 0.12356, 0.12356, 0.13071, 0.13822, 0.14608, 0.14874, 0.16578, 0.17316, 0.17444, 0.17316, 0.17220, 0.17167, 0.17114, 0.17068, 0.17047, 0.17068, 0.17150, 0.17316, 0.17586, 0.17962, 0.18430, 0.18935, 0.19425, 0.19927, 0.20499, 0.21238, 0.22278, 0.23849, 0.26174, 0.29340, 0.33353, 0.38028, 0.43126, 0.48394, 0.53577, 0.58432, 0.62791, 0.66763, 0.70476, 0.73852, 0.76792, 0.79343, 0.81553, 0.83334, 0.84628, 0.85638, 0.86565, 0.87355, 0.87932, 0.88378, 0.88789, 0.89152, 0.89442, 0.89689, 0.89928, 0.90145, 0.90319, 0.90443, 0.90521, 0.90572, 0.90617, 0.90672, 0.90745, 0.90830, 0.90916, 0.90989, 0.91076, 0.91157, 0.91237, 0.91316, 0.91394, 0.91472, 0.91550, 0.91626, 0.91703, 0.91778, 0.91853, 0.91927, 0.92001, 0.92074, 0.92146, 0.92218, 0.92218, 0.92218, 0.92218, 0.92218, 0.92218, 0.92218, 0.92218, 0.92218, 0.92218, 0.92218 } },
	{ 95, 360.0, 830.0, 1.0, { 0.00012, 0.00012, 0.00012, 0.00012, 0.00012, 0.01062, 0.02084, 0.03107, 0.04157, 0.05232, 0.06229, 0.07034, 0.07607, 0.07961, 0.08233, 0.08564, 0.08976, 0.09456, 0.09986, 0.10554, 0.11177, 0.11877, 0.12662, 0.13529, 0.14440, 0.15365, 0.16337, 0.17397, 0.18569, 0.19876, 0.21356, 0.23039, 0.24906, 0.26935, 0.29154, 0.31594, 0.34258, 0.37131, 0.40142, 0.43235, 0.46450, 0.49808, 0.53141, 0.56281, 0.59236, 0.62034, 0.64620, 0.66928, 0.68945, 0.70682, 0.72194, 0.73531, 0.74695, 0.75689, 0.76575, 0.77407, 0.78137, 0.78720, 0.79221, 0.79723, 0.80262, 0.80842, 0.81375, 0.81785, 0.82133, 0.82495, 0.82862, 0.83218, 0.83587, 0.83979, 0.84305, 0.84486, 0.84588, 0.84712, 0.84960, 0.85043, 0.85190, 0.85335, 0.85480, 0.85624, 0.85766, 0.85907, 0.86047, 0.86186, 0.86323, 0.86323, 0.86323, 0.86323, 0.86323, 0.86323, 0.86323, 0.86323, 0.86323, 0.86323, 0.86323 } },
	{ 95, 360.0, 830.0, 1.0, { 0.37853, 0.37853, 0.37853, 0.37853, 0.37853, 0.40061, 0.41716, 0.43082, 0.44274, 0.45364, 0.46388, 0.47366, 0.48294, 0.49171, 0.50041, 0.50937, 0.51814, 0.52630, 0.53427, 0.54260, 0.55151, 0.56099, 0.57041, 0.57915, 0.58732, 0.59526, 0.60357, 0.61278, 0.62303, 0.63425, 0.64599, 0.65789, 0.67036, 0.68406, 0.69997, 0.71856, 0.73784, 0.75548, 0.77019, 0.78135, 0.78986, 0.79689, 0.80310, 0.80879, 0.81341, 0.81640, 0.81801, 0.81876, 0.81934, 0.82031, 0.82163, 0.82311, 0.82467, 0.82636, 0.82855, 0.83152, 0.83481, 0.83794, 0.84114, 0.84471, 0.84844, 0.85209, 0.85572, 0.85943, 0.86329, 0.86723, 0.87082, 0.87373, 0.87657, 0.87984, 0.88271, 0.88441, 0.88579, 0.88770, 0.88946, 0.89030, 0.89074, 0.89141, 0.89200, 0.89212, 0.89200, 0.89196, 0.89200, 0.89204, 0.89200, 0.89200, 0.89200, 0.89200, 0.89200, 0.89200, 0.89200, 0.89200, 0.89200, 0.89200, 0.89200 } },
	{ 95, 360.0, 830.0, 1.0, { 0.20026, 0.20026, 0.20026, 0.20026, 0.20026, 0.21220, 0.22466, 0.23763, 0.25150, 0.26605, 0.28120, 0.28563, 0.28380, 0.28249, 0.28200, 0.28135, 0.28100, 0.28161, 0.28320, 0.28567, 0.28900, 0.29330, 0.29910, 0.30695, 0.31700, 0.32939, 0.34470, 0.36343, 0.38540, 0.41039, 0.43860, 0.47037, 0.50620, 0.54626, 0.58940, 0.63410, 0.67880, 0.72176, 0.76060, 0.79314, 0.81870, 0.83727, 0.85000, 0.85830, 0.86340, 0.86637, 0.86780, 0.86812, 0.86770, 0.86688, 0.86590, 0.86494, 0.86390, 0.86263, 0.86100, 0.85898, 0.85700, 0.85543, 0.85410, 0.85264, 0.85040, 0.84708, 0.84390, 0.84212, 0.84150, 0.84134, 0.84070, 0.83899, 0.83740, 0.83582, 0.83428, 0.83273, 0.83116, 0.82958, 0.82800, 0.82640, 0.82478, 0.82316, 0.82152, 0.81988, 0.81822, 0.81655, 0.81486, 0.81317, 0.81146, 0.81146, 0.81146, 0.81146, 0.81146, 0.81146, 0.81146, 0.81146, 0.81146, 0.81146, 0.81146 } },
	{ 95, 360.0, 830.0, 1.0, { 0.10281, 0.10281, 0.10281, 0.10281, 0.10281, 0.09292, 0.08404, 0.07621, 0.06814, 0.06240, 0.06158, 0.05854, 0.05464, 0.05502, 0.05832, 0.05840, 0.05462, 0.05290, 0.05304, 0.05210, 0.05075, 0.05128, 0.05300, 0.05316, 0.05188, 0.05240, 0.05510, 0.05851, 0.06343, 0.07062, 0.08218, 0.10589, 0.14388, 0.19506, 0.25332, 0.30750, 0.35271, 0.38845, 0.41659, 0.43753, 0.45206, 0.46315, 0.47193, 0.47900, 0.48644, 0.49282, 0.49629, 0.50020, 0.50456, 0.50626, 0.50779, 0.51042, 0.51055, 0.50984, 0.50953, 0.50935, 0.51270, 0.51496, 0.51194, 0.51094, 0.51107, 0.50985, 0.50758, 0.50201, 0.49937, 0.50130, 0.50444, 0.51390, 0.51821, 0.52735, 0.53455, 0.54175, 0.54893, 0.55608, 0.56322, 0.57032, 0.57740, 0.58445, 0.59146, 0.59844, 0.60537, 0.61227, 0.61911, 0.62592, 0.63267, 0.63267, 0.63267, 0.63267, 0.63267, 0.63267, 0.63267, 0.63267, 0.63267, 0.63267, 0.63267 } },
	{ 95, 360.0, 830.0, 1.0, { 0.18151, 0.18151, 0.18151, 0.18151, 0.18151, 0.17829, 0.17511, 0.17198, 0.16880, 0.16558, 0.16160, 0.15900, 0.15810, 0.15881, 0.16130, 0.16571, 0.17170, 0.17895, 0.18760, 0.19777, 0.20900, 0.22112, 0.23560, 0.25362, 0.27340, 0.29329, 0.31500, 0.34048, 0.36920, 0.40017, 0.43290, 0.46717, 0.50320, 0.54088, 0.57840, 0.61400, 0.64770, 0.67979, 0.70970, 0.73677, 0.76070, 0.78141, 0.79940, 0.81523, 0.82920, 0.84155, 0.85240, 0.86185, 0.86990, 0.87659, 0.88230, 0.88735, 0.89160, 0.89487, 0.89740, 0.89946, 0.90110, 0.90232, 0.90330, 0.90419, 0.90500, 0.90570, 0.90630, 0.90680, 0.90720, 0.90752, 0.90790, 0.90840, 0.90880, 0.90918, 0.90956, 0.90993, 0.91030, 0.91066, 0.91103, 0.91140, 0.91176, 0.91212, 0.91249, 0.91285, 0.91320, 0.91356, 0.91392, 0.91427, 0.91463, 0.91463, 0.91463, 0.91463, 0.91463, 0.91463, 0.91463, 0.91463, 0.91463, 0.91463, 0.91463 } },
	{ 95, 360.0, 830.0, 1.0, { 0.06462, 0.06462, 0.06462, 0.06462, 0.06462, 0.06500, 0.06620, 0.06758, 0.06839, 0.07185, 0.08117, 0.09389, 0.11434, 0.14548, 0.17980, 0.20758, 0.22587, 0.24072, 0.25486, 0.26272, 0.26805, 0.27812, 0.28681, 0.29048, 0.29200, 0.29516, 0.30288, 0.31468, 0.32907, 0.34463, 0.35691, 0.36974, 0.38422, 0.39686, 0.40592, 0.41045, 0.41237, 0.41646, 0.42207, 0.42657, 0.42941, 0.43090, 0.43186, 0.43417, 0.43799, 0.44134, 0.44276, 0.44410, 0.44574, 0.44656, 0.44617, 0.44794, 0.44909, 0.44919, 0.45155, 0.45456, 0.45681, 0.45528, 0.45003, 0.44666, 0.44579, 0.44445, 0.44430, 0.44272, 0.44242, 0.44704, 0.45358, 0.46121, 0.46357, 0.47205, 0.47816, 0.48429, 0.49042, 0.49655, 0.50268, 0.50881, 0.51494, 0.52106, 0.52718, 0.53329, 0.53939, 0.54548, 0.55155, 0.55761, 0.56365, 0.56365, 0.56365, 0.56365, 0.56365, 0.56365, 0.56365, 0.56365, 0.56365, 0.56365, 0.56365 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04580, 0.04580, 0.04580, 0.04580, 0.04580, 0.04489, 0.04394, 0.04317, 0.04467, 0.04650, 0.04793, 0.04699, 0.04621, 0.04400, 0.04456, 0.04653, 0.04711, 0.04794, 0.04801, 0.04537, 0.04556, 0.04574, 0.04704, 0.04889, 0.04999, 0.05113, 0.05184, 0.05266, 0.05780, 0.06583, 0.07923, 0.09316, 0.10443, 0.11237, 0.11717, 0.11799, 0.11843, 0.11869, 0.11877, 0.11740, 0.11503, 0.12031, 0.12027, 0.12087, 0.12009, 0.12214, 0.12237, 0.12116, 0.11951, 0.11896, 0.11864, 0.11826, 0.11813, 0.11937, 0.11881, 0.11634, 0.11621, 0.11637, 0.11656, 0.11571, 0.11537, 0.11267, 0.11420, 0.11624, 0.11636, 0.11320, 0.11456, 0.11781, 0.11627, 0.11443, 0.11559, 0.11697, 0.11711, 0.11871, 0.11979, 0.11884, 0.11780, 0.11833, 0.11880, 0.11906, 0.12030, 0.12106, 0.12099, 0.12120, 0.12123, 0.12123, 0.12123, 0.12123, 0.12123, 0.12123, 0.12123, 0.12123, 0.12123, 0.12123, 0.12123 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03981, 0.03981, 0.03981, 0.03981, 0.03981, 0.04155, 0.04336, 0.04525, 0.04798, 0.04887, 0.05115, 0.05478, 0.05972, 0.06611, 0.07483, 0.08683, 0.10270, 0.12281, 0.14711, 0.17523, 0.20611, 0.23842, 0.27043, 0.30066, 0.32892, 0.35533, 0.37997, 0.40328, 0.42723, 0.45369, 0.48277, 0.51378, 0.54473, 0.57367, 0.60006, 0.62365, 0.64365, 0.65956, 0.67224, 0.68280, 0.69184, 0.69975, 0.70674, 0.71305, 0.71910, 0.72521, 0.73114, 0.73663, 0.74196, 0.74751, 0.75350, 0.75999, 0.76667, 0.77308, 0.77882, 0.78351, 0.78709, 0.78962, 0.79148, 0.79321, 0.79597, 0.80060, 0.80638, 0.81236, 0.81832, 0.82419, 0.82976, 0.83482, 0.83915, 0.84423, 0.84885, 0.85335, 0.85775, 0.86203, 0.86621, 0.87028, 0.87424, 0.87810, 0.88185, 0.88551, 0.88906, 0.89252, 0.89589, 0.89916, 0.90233, 0.90233, 0.90233, 0.90233, 0.90233, 0.90233, 0.90233, 0.90233, 0.90233, 0.90233, 0.90233 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03043, 0.03043, 0.03043, 0.03043, 0.03043, 0.03198, 0.03135, 0.02956, 0.02799, 0.02943, 0.03159, 0.03054, 0.03068, 0.03510, 0.04618, 0.05858, 0.07378, 0.09711, 0.12843, 0.16071, 0.18945, 0.21301, 0.23025, 0.23860, 0.24313, 0.24970, 0.25715, 0.26432, 0.27011, 0.27155, 0.27200, 0.27353, 0.27630, 0.27819, 0.27862, 0.27889, 0.27955, 0.28173, 0.28524, 0.28833, 0.29177, 0.29664, 0.30154, 0.30507, 0.30871, 0.31340, 0.31903, 0.32509, 0.32974, 0.33388, 0.33894, 0.34220, 0.34451, 0.34661, 0.34779, 0.35040, 0.35438, 0.35485, 0.35338, 0.35403, 0.35561, 0.35852, 0.35971, 0.35739, 0.35449, 0.35373, 0.35540, 0.36075, 0.36514, 0.36959, 0.37419, 0.37882, 0.38348, 0.38815, 0.39284, 0.39756, 0.40229, 0.40704, 0.41181, 0.41660, 0.42140, 0.42622, 0.43105, 0.43589, 0.44074, 0.44074, 0.44074, 0.44074, 0.44074, 0.44074, 0.44074, 0.44074, 0.44074, 0.44074, 0.44074 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04401, 0.04401, 0.04401, 0.04401, 0.04401, 0.04373, 0.04499, 0.04761, 0.05142, 0.05624, 0.06182, 0.06809, 0.07569, 0.08532, 0.09712, 0.11108, 0.12710, 0.14513, 0.16551, 0.18850, 0.21364, 0.24081, 0.27189, 0.30813, 0.34621, 0.38299, 0.42055, 0.46097, 0.50101, 0.53702, 0.56908, 0.59786, 0.62260, 0.64264, 0.65921, 0.67373, 0.68657, 0.69790, 0.70804, 0.71734, 0.72615, 0.73467, 0.74241, 0.74895, 0.75484, 0.76070, 0.76648, 0.77190, 0.77655, 0.78017, 0.78326, 0.78641, 0.78978, 0.79338, 0.79694, 0.80012, 0.80264, 0.80445, 0.80640, 0.80916, 0.81170, 0.81302, 0.81388, 0.81524, 0.81722, 0.81966, 0.82215, 0.82429, 0.82611, 0.82769, 0.82899, 0.82993, 0.83048, 0.83060, 0.83026, 0.83060, 0.83067, 0.83074, 0.83081, 0.83088, 0.83095, 0.83102, 0.83109, 0.83116, 0.83123, 0.83123, 0.83123, 0.83123, 0.83123, 0.83123, 0.83123, 0.83123, 0.83123, 0.83123, 0.83123 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04628, 0.04628, 0.04628, 0.04628, 0.04628, 0.04542, 0.04457, 0.04374, 0.04290, 0.04206, 0.04110, 0.04069, 0.04100, 0.04209, 0.04400, 0.04679, 0.05050, 0.05522, 0.06130, 0.06913, 0.07890, 0.09077, 0.10500, 0.12180, 0.14120, 0.16307, 0.18690, 0.21218, 0.23870, 0.26636, 0.29520, 0.32517, 0.35580, 0.38644, 0.41620, 0.44402, 0.46830, 0.48756, 0.50130, 0.50941, 0.51240, 0.51113, 0.50720, 0.50235, 0.49800, 0.49512, 0.49310, 0.49093, 0.48760, 0.48248, 0.47630, 0.47030, 0.46630, 0.46623, 0.47200, 0.48518, 0.50610, 0.53447, 0.56880, 0.60699, 0.64570, 0.68173, 0.71380, 0.74151, 0.76620, 0.78956, 0.81310, 0.83652, 0.85260, 0.86699, 0.87979, 0.89151, 0.90221, 0.91196, 0.92083, 0.92887, 0.93615, 0.94273, 0.94867, 0.95403, 0.95885, 0.96318, 0.96707, 0.97057, 0.97370, 0.97370, 0.97370, 0.97370, 0.97370, 0.97370, 0.97370, 0.97370, 0.97370, 0.97370, 0.97370 } },
	{ 95, 360.0, 830.0, 1.0, { 0.17315, 0.17315, 0.17315, 0.17315, 0.17315, 0.17955, 0.18551, 0.19147, 0.19787, 0.20470, 0.21023, 0.21321, 0.21595, 0.22112, 0.22921, 0.24063, 0.25764, 0.28283, 0.31841, 0.36487, 0.41631, 0.46632, 0.51290, 0.55466, 0.58809, 0.61048, 0.62424, 0.63310, 0.64079, 0.65018, 0.66068, 0.67101, 0.68048, 0.68881, 0.69672, 0.70489, 0.71287, 0.72001, 0.72588, 0.73048, 0.73522, 0.74123, 0.74699, 0.75091, 0.75383, 0.75692, 0.76021, 0.76361, 0.76763, 0.77265, 0.77797, 0.78271, 0.78627, 0.78844, 0.79028, 0.79278, 0.79528, 0.79706, 0.79877, 0.80118, 0.80412, 0.80718, 0.80995, 0.81215, 0.81403, 0.81593, 0.81807, 0.82055, 0.82326, 0.82593, 0.82804, 0.82924, 0.83027, 0.83213, 0.83581, 0.83693, 0.83904, 0.84114, 0.84321, 0.84525, 0.84728, 0.84929, 0.85127, 0.85323, 0.85517, 0.85517, 0.85517, 0.85517, 0.85517, 0.85517, 0.85517, 0.85517, 0.85517, 0.85517, 0.85517 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04354, 0.04354, 0.04354, 0.04354, 0.04354, 0.04363, 0.04371, 0.04380, 0.04389, 0.04397, 0.04404, 0.04401, 0.04400, 0.04411, 0.04433, 0.04464, 0.04511, 0.04589, 0.04721, 0.04923, 0.05164, 0.05427, 0.05796, 0.06361, 0.07121, 0.08063, 0.09207, 0.10573, 0.12138, 0.13906, 0.16048, 0.18727, 0.21932, 0.25474, 0.28642, 0.30780, 0.31975, 0.32497, 0.32621, 0.32579, 0.32446, 0.32253, 0.32020, 0.31762, 0.31486, 0.31189, 0.30849, 0.30443, 0.29971, 0.29443, 0.28917, 0.28446, 0.28043, 0.27707, 0.27414, 0.27141, 0.26885, 0.26649, 0.26449, 0.26303, 0.26213, 0.26186, 0.26248, 0.26410, 0.26609, 0.26786, 0.26981, 0.27228, 0.27417, 0.27602, 0.27782, 0.27963, 0.28144, 0.28326, 0.28509, 0.28692, 0.28876, 0.29061, 0.29247, 0.29433, 0.29620, 0.29808, 0.29996, 0.30185, 0.30375, 0.30375, 0.30375, 0.30375, 0.30375, 0.30375, 0.30375, 0.30375, 0.30375, 0.30375, 0.30375 } },
	{ 95, 360.0, 830.0, 1.0, { 0.05265, 0.05265, 0.05265, 0.05265, 0.05265, 0.05260, 0.05255, 0.05250, 0.05245, 0.05239, 0.05217, 0.05170, 0.05108, 0.05042, 0.04980, 0.04932, 0.04899, 0.04883, 0.04885, 0.04903, 0.04929, 0.04960, 0.05029, 0.05178, 0.05440, 0.05841, 0.06371, 0.06999, 0.07648, 0.08244, 0.08763, 0.09226, 0.09785, 0.10570, 0.11485, 0.12372, 0.13045, 0.13368, 0.13428, 0.13355, 0.13216, 0.13054, 0.12886, 0.12720, 0.12554, 0.12381, 0.12196, 0.11992, 0.11767, 0.11522, 0.11283, 0.11074, 0.10894, 0.10738, 0.10598, 0.10468, 0.10349, 0.10242, 0.10154, 0.10088, 0.10041, 0.10012, 0.10019, 0.10071, 0.10131, 0.10168, 0.10213, 0.10293, 0.10364, 0.10434, 0.10502, 0.10571, 0.10640, 0.10709, 0.10779, 0.10850, 0.10920, 0.10991, 0.11063, 0.11135, 0.11207, 0.11279, 0.11353, 0.11426, 0.11500, 0.11500, 0.11500, 0.11500, 0.11500, 0.11500, 0.11500, 0.11500, 0.11500, 0.11500, 0.11500 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03431, 0.03431, 0.03431, 0.03431, 0.03431, 0.03708, 0.04283, 0.05358, 0.07506, 0.09776, 0.12729, 0.16371, 0.20727, 0.25618, 0.30712, 0.35739, 0.40413, 0.44397, 0.47549, 0.49917, 0.51519, 0.52400, 0.52781, 0.52870, 0.52832, 0.52853, 0.52959, 0.53055, 0.53092, 0.53193, 0.53347, 0.53402, 0.53283, 0.53039, 0.52860, 0.52848, 0.52890, 0.52837, 0.52700, 0.52649, 0.52843, 0.53210, 0.53516, 0.53670, 0.53798, 0.54007, 0.54310, 0.54624, 0.54817, 0.54899, 0.55035, 0.55266, 0.55518, 0.55746, 0.55801, 0.55702, 0.55717, 0.55869, 0.55967, 0.55989, 0.56082, 0.56192, 0.56165, 0.56041, 0.55912, 0.55486, 0.54941, 0.54330, 0.53563, 0.52630, 0.51786, 0.51116, 0.50353, 0.49575, 0.48819, 0.48065, 0.47311, 0.46558, 0.45806, 0.45057, 0.44310, 0.43562, 0.42817, 0.42076, 0.41339, 0.41339, 0.41339, 0.41339, 0.41339, 0.41339, 0.41339, 0.41339, 0.41339, 0.41339, 0.41339 } },
	{ 95, 360.0, 830.0, 1.0, { 0.05028, 0.05028, 0.05028, 0.05028, 0.05028, 0.05023, 0.05019, 0.05014, 0.05000, 0.05009, 0.05000, 0.04991, 0.05000, 0.05039, 0.05100, 0.05177, 0.05300, 0.05490, 0.05700, 0.05875, 0.06000, 0.06073, 0.06100, 0.06109, 0.06200, 0.06454, 0.06800, 0.07115, 0.07200, 0.07038, 0.07400, 0.08994, 0.11500, 0.14309, 0.16700, 0.18168, 0.19200, 0.20293, 0.21000, 0.20811, 0.19900, 0.18615, 0.17300, 0.16242, 0.15500, 0.15106, 0.15200, 0.15771, 0.16100, 0.15559, 0.14600, 0.13792, 0.13100, 0.12386, 0.11700, 0.11101, 0.10500, 0.09848, 0.09400, 0.09471, 0.10310, 0.12080, 0.14660, 0.17877, 0.21630, 0.25833, 0.30380, 0.35162, 0.40070, 0.45484, 0.50886, 0.56267, 0.61505, 0.66489, 0.71130, 0.75367, 0.79164, 0.82512, 0.85421, 0.87917, 0.90035, 0.91817, 0.93303, 0.94536, 0.95553, 0.95553, 0.95553, 0.95553, 0.95553, 0.95553, 0.95553, 0.95553, 0.95553, 0.95553, 0.95553 } },
	{ 95, 360.0, 830.0, 1.0, { 0.11323, 0.11323, 0.11323, 0.11323, 0.11323, 0.11525, 0.10979, 0.09541, 0.08689, 0.08224, 0.08142, 0.08681, 0.10151, 0.12741, 0.16431, 0.21014, 0.26265, 0.31889, 0.37622, 0.43060, 0.47844, 0.51829, 0.55029, 0.57599, 0.59562, 0.61049, 0.62077, 0.62698, 0.63270, 0.64011, 0.64719, 0.65171, 0.65373, 0.65447, 0.65538, 0.65739, 0.65954, 0.66027, 0.65923, 0.65779, 0.65799, 0.65940, 0.66049, 0.66094, 0.66122, 0.66190, 0.66315, 0.66458, 0.66544, 0.66585, 0.66684, 0.66881, 0.67060, 0.67138, 0.67119, 0.67049, 0.67041, 0.67108, 0.67106, 0.67024, 0.66944, 0.66846, 0.66851, 0.67028, 0.67149, 0.67258, 0.67428, 0.67638, 0.67857, 0.68190, 0.68592, 0.68952, 0.69328, 0.69711, 0.70083, 0.70452, 0.70819, 0.71183, 0.71544, 0.71903, 0.72258, 0.72615, 0.72967, 0.73316, 0.73661, 0.73661, 0.73661, 0.73661, 0.73661, 0.73661, 0.73661, 0.73661, 0.73661, 0.73661, 0.73661 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04566, 0.04566, 0.04566, 0.04566, 0.04566, 0.04573, 0.04580, 0.04588, 0.04606, 0.04597, 0.04606, 0.04634, 0.04678, 0.04741, 0.04823, 0.04925, 0.05039, 0.05160, 0.05297, 0.05453, 0.05595, 0.05698, 0.05791, 0.05932, 0.06224, 0.06742, 0.07389, 0.08031, 0.08563, 0.08916, 0.09130, 0.09256, 0.09254, 0.09086, 0.08800, 0.08461, 0.08110, 0.07777, 0.07461, 0.07165, 0.06935, 0.06811, 0.06770, 0.06784, 0.06863, 0.07018, 0.07224, 0.07444, 0.07646, 0.07805, 0.07924, 0.08017, 0.08089, 0.08147, 0.08192, 0.08233, 0.08295, 0.08401, 0.08532, 0.08661, 0.08759, 0.08807, 0.08821, 0.08823, 0.08821, 0.08818, 0.08811, 0.08796, 0.08769, 0.08760, 0.08743, 0.08727, 0.08711, 0.08694, 0.08678, 0.08662, 0.08645, 0.08629, 0.08613, 0.08597, 0.08581, 0.08565, 0.08549, 0.08533, 0.08517, 0.08517, 0.08517, 0.08517, 0.08517, 0.08517, 0.08517, 0.08517, 0.08517, 0.08517, 0.08517 } },
	{ 95, 360.0, 830.0, 1.0, { 0.07500, 0.07500, 0.07500, 0.07500, 0.07500, 0.07502, 0.07503, 0.07504, 0.07522, 0.07498, 0.07502, 0.07534, 0.07595, 0.07683, 0.07801, 0.07950, 0.08130, 0.08347, 0.08625, 0.08964, 0.09243, 0.09371, 0.09480, 0.09794, 0.10655, 0.12302, 0.14427, 0.16617, 0.18600, 0.20198, 0.21465, 0.22435, 0.22825, 0.22393, 0.21393, 0.20158, 0.18858, 0.17596, 0.16384, 0.15244, 0.14313, 0.13716, 0.13396, 0.13286, 0.13458, 0.13978, 0.14736, 0.15579, 0.16364, 0.16980, 0.17456, 0.17846, 0.18157, 0.18391, 0.18569, 0.18725, 0.18930, 0.19238, 0.19610, 0.19976, 0.20259, 0.20400, 0.20434, 0.20413, 0.20352, 0.20254, 0.20115, 0.19927, 0.19682, 0.19524, 0.19333, 0.19144, 0.18956, 0.18769, 0.18584, 0.18400, 0.18218, 0.18037, 0.17857, 0.17679, 0.17502, 0.17327, 0.17153, 0.16981, 0.16809, 0.16809, 0.16809, 0.16809, 0.16809, 0.16809, 0.16809, 0.16809, 0.16809, 0.16809, 0.16809 } },
	{ 95, 360.0, 830.0, 1.0, { 0.37973, 0.37973, 0.37973, 0.37973, 0.37973, 0.38712, 0.39455, 0.40204, 0.40980, 0.41775, 0.42690, 0.43209, 0.43680, 0.44546, 0.45770, 0.47231, 0.48950, 0.50977, 0.53330, 0.55996, 0.58880, 0.61854, 0.64750, 0.67396, 0.69650, 0.71398, 0.72620, 0.73336, 0.73640, 0.73638, 0.73420, 0.73061, 0.72610, 0.72102, 0.71550, 0.70967, 0.70390, 0.69860, 0.69410, 0.69062, 0.68810, 0.68630, 0.68460, 0.68243, 0.67990, 0.67729, 0.67490, 0.67300, 0.67180, 0.67137, 0.67140, 0.67164, 0.67250, 0.67464, 0.67890, 0.68616, 0.69710, 0.71218, 0.73100, 0.75284, 0.77640, 0.80025, 0.82300, 0.84343, 0.86090, 0.87518, 0.88720, 0.89754, 0.90410, 0.91008, 0.91555, 0.92071, 0.92558, 0.93018, 0.93451, 0.93859, 0.94243, 0.94604, 0.94945, 0.95264, 0.95565, 0.95847, 0.96112, 0.96360, 0.96594, 0.96594, 0.96594, 0.96594, 0.96594, 0.96594, 0.96594, 0.96594, 0.96594, 0.96594, 0.96594 } },
	{ 95, 360.0, 830.0, 1.0, { 0.06598, 0.06598, 0.06598, 0.06598, 0.06598, 0.06537, 0.06477, 0.06418, 0.06374, 0.06292, 0.06235, 0.06213, 0.06235, 0.06307, 0.06414, 0.06540, 0.06682, 0.06845, 0.07050, 0.07328, 0.07725, 0.08275, 0.08966, 0.09786, 0.10793, 0.12099, 0.13980, 0.16761, 0.20821, 0.26399, 0.33133, 0.40352, 0.46746, 0.51124, 0.53408, 0.53861, 0.52991, 0.51304, 0.49059, 0.46466, 0.43777, 0.41211, 0.38793, 0.36527, 0.34523, 0.32877, 0.31534, 0.30404, 0.29410, 0.28480, 0.27563, 0.26625, 0.25696, 0.24839, 0.24177, 0.23840, 0.23909, 0.24436, 0.25398, 0.26709, 0.28089, 0.29267, 0.30194, 0.30939, 0.31803, 0.33103, 0.34980, 0.37529, 0.40848, 0.43191, 0.45996, 0.48826, 0.51664, 0.54491, 0.57290, 0.60042, 0.62733, 0.65347, 0.67872, 0.70296, 0.72611, 0.74810, 0.76889, 0.78845, 0.80677, 0.80677, 0.80677, 0.80677, 0.80677, 0.80677, 0.80677, 0.80677, 0.80677, 0.80677, 0.80677 } },
	{ 95, 360.0, 830.0, 1.0, { 0.08862, 0.08862, 0.08862, 0.08862, 0.08862, 0.08277, 0.07727, 0.07211, 0.06710, 0.06228, 0.05620, 0.05196, 0.04970, 0.04885, 0.04890, 0.04950, 0.05090, 0.05356, 0.05820, 0.06544, 0.07510, 0.08761, 0.10660, 0.13515, 0.17100, 0.21203, 0.26190, 0.32330, 0.38910, 0.45024, 0.49970, 0.53263, 0.55070, 0.55722, 0.55550, 0.54854, 0.53800, 0.52506, 0.51020, 0.49382, 0.47680, 0.46002, 0.44390, 0.42861, 0.41390, 0.39905, 0.38200, 0.36163, 0.34190, 0.32717, 0.31840, 0.31512, 0.31440, 0.31336, 0.31170, 0.30968, 0.30730, 0.30349, 0.29310, 0.27208, 0.24480, 0.21777, 0.19760, 0.18894, 0.18850, 0.19437, 0.21810, 0.26567, 0.30730, 0.35107, 0.39594, 0.44264, 0.49037, 0.53828, 0.58549, 0.63118, 0.67464, 0.71528, 0.75271, 0.78669, 0.81713, 0.84409, 0.86771, 0.88823, 0.90592, 0.90592, 0.90592, 0.90592, 0.90592, 0.90592, 0.90592, 0.90592, 0.90592, 0.90592, 0.90592 } },
	{ 95, 360.0, 830.0, 1.0, { 0.04163, 0.04163, 0.04163, 0.04163, 0.04163, 0.04219, 0.04274, 0.04331, 0.04390, 0.04451, 0.04520, 0.04555, 0.04570, 0.04589, 0.04610, 0.04629, 0.04650, 0.04676, 0.04700, 0.04714, 0.04720, 0.04724, 0.04740, 0.04779, 0.04830, 0.04879, 0.04930, 0.04986, 0.05040, 0.05082, 0.05110, 0.05122, 0.05120, 0.05107, 0.05090, 0.05076, 0.05060, 0.05036, 0.05010, 0.04987, 0.04960, 0.04922, 0.04890, 0.04879, 0.04880, 0.04880, 0.04880, 0.04882, 0.04880, 0.04868, 0.04860, 0.04868, 0.04890, 0.04918, 0.04950, 0.04984, 0.05020, 0.05057, 0.05100, 0.05150, 0.05200, 0.05242, 0.05270, 0.05278, 0.05270, 0.05252, 0.05230, 0.05211, 0.05200, 0.05189, 0.05179, 0.05168, 0.05158, 0.05147, 0.05137, 0.05127, 0.05116, 0.05106, 0.05096, 0.05085, 0.05075, 0.05065, 0.05055, 0.05045, 0.05034, 0.05034, 0.05034, 0.05034, 0.05034, 0.05034, 0.05034, 0.05034, 0.05034, 0.05034, 0.05034 } },
	{ 95, 360.0, 830.0, 1.0, { 0.05540, 0.05540, 0.05540, 0.05540, 0.05540, 0.05519, 0.05609, 0.05723, 0.05760, 0.05773, 0.05841, 0.05900, 0.06076, 0.06181, 0.06364, 0.06536, 0.06841, 0.07063, 0.07333, 0.07831, 0.08470, 0.09311, 0.10313, 0.11523, 0.13064, 0.14961, 0.17359, 0.20369, 0.23996, 0.27414, 0.30443, 0.32537, 0.33756, 0.34253, 0.34239, 0.33797, 0.33044, 0.32147, 0.31229, 0.29914, 0.28479, 0.27639, 0.26307, 0.24861, 0.23426, 0.21943, 0.20576, 0.19114, 0.17990, 0.17159, 0.16571, 0.16047, 0.15739, 0.15454, 0.15364, 0.15204, 0.15024, 0.15039, 0.15037, 0.15041, 0.15234, 0.15367, 0.15924, 0.16394, 0.16737, 0.16914, 0.17330, 0.17773, 0.17987, 0.18020, 0.18000, 0.17773, 0.17586, 0.17644, 0.17847, 0.18116, 0.18747, 0.19526, 0.20319, 0.20930, 0.21407, 0.21864, 0.21997, 0.21784, 0.21926, 0.21926, 0.21926, 0.21926, 0.21926, 0.21926, 0.21926, 0.21926, 0.21926, 0.21926, 0.21926 } },
	{ 95, 360.0, 830.0, 1.0, { 0.09982, 0.09982, 0.09982, 0.09982, 0.09982, 0.10287, 0.10599, 0.10920, 0.11260, 0.11615, 0.12020, 0.12231, 0.12350, 0.12532, 0.12800, 0.13152, 0.13620, 0.14245, 0.15070, 0.16141, 0.17520, 0.19261, 0.21370, 0.23822, 0.26510, 0.29289, 0.31950, 0.34297, 0.36240, 0.37735, 0.38800, 0.39477, 0.39830, 0.39931, 0.39860, 0.39682, 0.39400, 0.39003, 0.38490, 0.37870, 0.37190, 0.36490, 0.35750, 0.34938, 0.34030, 0.33000, 0.31810, 0.30447, 0.29020, 0.27656, 0.26440, 0.25433, 0.24630, 0.24004, 0.23490, 0.23025, 0.22600, 0.22218, 0.21890, 0.21633, 0.21470, 0.21422, 0.21480, 0.21628, 0.21840, 0.22098, 0.22420, 0.22799, 0.23080, 0.23356, 0.23624, 0.23895, 0.24168, 0.24443, 0.24720, 0.24999, 0.25280, 0.25564, 0.25849, 0.26137, 0.26426, 0.26718, 0.27011, 0.27307, 0.27605, 0.27605, 0.27605, 0.27605, 0.27605, 0.27605, 0.27605, 0.27605, 0.27605, 0.27605, 0.27605 } },
	{ 95, 360.0, 830.0, 1.0, { 0.14254, 0.14254, 0.14254, 0.14254, 0.14254, 0.14305, 0.14356, 0.14407, 0.14460, 0.14514, 0.14580, 0.14646, 0.14840, 0.15294, 0.16010, 0.16970, 0.18200, 0.19737, 0.21620, 0.23886, 0.26560, 0.29650, 0.33100, 0.36800, 0.40480, 0.43811, 0.46400, 0.47933, 0.48480, 0.48229, 0.47450, 0.46402, 0.45220, 0.43995, 0.42770, 0.41577, 0.40450, 0.39412, 0.38450, 0.37543, 0.36690, 0.35895, 0.35180, 0.34575, 0.34130, 0.33870, 0.33700, 0.33506, 0.33220, 0.32808, 0.32320, 0.31848, 0.31560, 0.31636, 0.32210, 0.33395, 0.35250, 0.37809, 0.41040, 0.44884, 0.49220, 0.53897, 0.58690, 0.63375, 0.67790, 0.71834, 0.75580, 0.78958, 0.81150, 0.83107, 0.84844, 0.86431, 0.87876, 0.89187, 0.90371, 0.91438, 0.92396, 0.93256, 0.94024, 0.94710, 0.95321, 0.95864, 0.96347, 0.96776, 0.97155, 0.97155, 0.97155, 0.97155, 0.97155, 0.97155, 0.97155, 0.97155, 0.97155, 0.97155, 0.97155 } },
	{ 95, 360.0, 830.0, 1.0, { 0.08776, 0.08776, 0.08776, 0.08776, 0.08776, 0.08843, 0.08933, 0.09022, 0.09089, 0.09134, 0.09245, 0.09578, 0.10461, 0.12154, 0.14474, 0.17179, 0.20227, 0.23616, 0.27292, 0.31153, 0.34940, 0.38408, 0.41516, 0.44282, 0.46753, 0.48971, 0.50918, 0.52583, 0.54042, 0.55383, 0.56646, 0.57838, 0.58866, 0.59648, 0.60236, 0.60683, 0.60921, 0.60857, 0.60427, 0.59607, 0.58490, 0.57159, 0.55521, 0.53495, 0.51211, 0.48821, 0.46354, 0.43824, 0.41298, 0.38846, 0.36486, 0.34231, 0.32117, 0.30177, 0.28420, 0.26844, 0.25429, 0.24158, 0.23040, 0.22080, 0.21248, 0.20511, 0.19862, 0.19300, 0.18827, 0.18444, 0.18137, 0.17893, 0.17713, 0.17599, 0.17550, 0.17564, 0.17652, 0.17830, 0.18111, 0.18243, 0.18437, 0.18633, 0.18831, 0.19030, 0.19230, 0.19433, 0.19637, 0.19842, 0.20049, 0.20049, 0.20049, 0.20049, 0.20049, 0.20049, 0.20049, 0.20049, 0.20049, 0.20049, 0.20049 } },
	{ 95, 360.0, 830.0, 1.0, { 0.00057, 0.00057, 0.00057, 0.00057, 0.00057, 0.00739, 0.01389, 0.02038, 0.02720, 0.03434, 0.04051, 0.04446, 0.04633, 0.04675, 0.04681, 0.04748, 0.04889, 0.05103, 0.05418, 0.05879, 0.06559, 0.07524, 0.08795, 0.10363, 0.12163, 0.14087, 0.15921, 0.17469, 0.18701, 0.19627, 0.20241, 0.20542, 0.20553, 0.20313, 0.19895, 0.19369, 0.18754, 0.18056, 0.17279, 0.16436, 0.15574, 0.14741, 0.13947, 0.13197, 0.12523, 0.11951, 0.11449, 0.10977, 0.10532, 0.10120, 0.09750, 0.09436, 0.09190, 0.09028, 0.08950, 0.08951, 0.09004, 0.09084, 0.09188, 0.09319, 0.09465, 0.09614, 0.09758, 0.09892, 0.10023, 0.10164, 0.10336, 0.10553, 0.10789, 0.11019, 0.11246, 0.11484, 0.11757, 0.12090, 0.12509, 0.12825, 0.13195, 0.13574, 0.13962, 0.14359, 0.14766, 0.15182, 0.15608, 0.16044, 0.16489, 0.16489, 0.16489, 0.16489, 0.16489, 0.16489, 0.16489, 0.16489, 0.16489, 0.16489, 0.16489 } },
	{ 95, 360.0, 830.0, 1.0, { 0.07883, 0.07883, 0.07883, 0.07883, 0.07883, 0.08225, 0.08580, 0.08949, 0.09339, 0.09727, 0.10142, 0.10575, 0.11016, 0.11458, 0.11919, 0.12413, 0.12931, 0.13467, 0.14052, 0.14714, 0.15412, 0.16118, 0.16930, 0.18000, 0.19580, 0.21832, 0.24453, 0.27019, 0.29107, 0.30412, 0.31122, 0.31477, 0.31449, 0.30979, 0.30149, 0.29047, 0.27638, 0.25896, 0.23937, 0.21933, 0.20136, 0.18742, 0.17665, 0.16780, 0.16107, 0.15676, 0.15442, 0.15328, 0.15223, 0.15040, 0.14826, 0.14662, 0.14618, 0.14745, 0.15005, 0.15358, 0.15829, 0.16436, 0.17109, 0.17743, 0.18201, 0.18381, 0.18359, 0.18240, 0.18052, 0.17803, 0.17486, 0.17093, 0.16613, 0.16279, 0.15896, 0.15521, 0.15154, 0.14793, 0.14440, 0.14093, 0.13754, 0.13421, 0.13096, 0.12777, 0.12464, 0.12159, 0.11859, 0.11566, 0.11280, 0.11280, 0.11280, 0.11280, 0.11280, 0.11280, 0.11280, 0.11280, 0.11280, 0.11280, 0.11280 } },
	{ 95, 360.0, 830.0, 1.0, { 0.09950, 0.09950, 0.09950, 0.09950, 0.09950, 0.10406, 0.10880, 0.11373, 0.11889, 0.12416, 0.12971, 0.13543, 0.14122, 0.14703, 0.15303, 0.15937, 0.16583, 0.17226, 0.17913, 0.18680, 0.19471, 0.20244, 0.21118, 0.22283, 0.24056, 0.26640, 0.29663, 0.32597, 0.34883, 0.36119, 0.36570, 0.36586, 0.36203, 0.35401, 0.34248, 0.32800, 0.30983, 0.28742, 0.26249, 0.23757, 0.21634, 0.20167, 0.19213, 0.18575, 0.18280, 0.18373, 0.18766, 0.19304, 0.19699, 0.19725, 0.19521, 0.19302, 0.19233, 0.19418, 0.19788, 0.20264, 0.20890, 0.21709, 0.22617, 0.23465, 0.24066, 0.24281, 0.24224, 0.24044, 0.23778, 0.23432, 0.22994, 0.22450, 0.21783, 0.21320, 0.20789, 0.20268, 0.19757, 0.19256, 0.18764, 0.18282, 0.17810, 0.17347, 0.16894, 0.16450, 0.16016, 0.15591, 0.15175, 0.14769, 0.14371, 0.14371, 0.14371, 0.14371, 0.14371, 0.14371, 0.14371, 0.14371, 0.14371, 0.14371, 0.14371 } },
	{ 95, 360.0, 830.0, 1.0, { 0.05860, 0.05860, 0.05860, 0.05860, 0.05860, 0.05842, 0.05824, 0.05806, 0.05802, 0.05765, 0.05750, 0.05753, 0.05771, 0.05801, 0.05853, 0.05937, 0.06059, 0.06229, 0.06471, 0.06784, 0.07048, 0.07175, 0.07306, 0.07679, 0.08687, 0.10600, 0.13056, 0.15571, 0.17817, 0.19575, 0.20908, 0.21853, 0.22052, 0.21204, 0.19630, 0.17770, 0.15900, 0.14200, 0.12623, 0.11115, 0.09810, 0.08841, 0.08151, 0.07650, 0.07296, 0.07063, 0.06935, 0.06891, 0.06884, 0.06872, 0.06863, 0.06872, 0.06894, 0.06923, 0.06966, 0.07040, 0.07182, 0.07418, 0.07708, 0.07997, 0.08223, 0.08342, 0.08378, 0.08364, 0.08306, 0.08200, 0.08048, 0.07849, 0.07605, 0.07441, 0.07251, 0.07065, 0.06884, 0.06707, 0.06535, 0.06366, 0.06202, 0.06041, 0.05885, 0.05732, 0.05583, 0.05438, 0.05296, 0.05158, 0.05023, 0.05023, 0.05023, 0.05023, 0.05023, 0.05023, 0.05023, 0.05023, 0.05023, 0.05023, 0.05023 } },
	{ 95, 360.0, 830.0, 1.0, { 0.07278, 0.07278, 0.07278, 0.07278, 0.07278, 0.07402, 0.07528, 0.07656, 0.07790, 0.07930, 0.08120, 0.08291, 0.08480, 0.08736, 0.09050, 0.09403, 0.09800, 0.10263, 0.10860, 0.11658, 0.12670, 0.13900, 0.15380, 0.17152, 0.19280, 0.21822, 0.24790, 0.28105, 0.31370, 0.34111, 0.35850, 0.36251, 0.35540, 0.34058, 0.32030, 0.29646, 0.27060, 0.24412, 0.21820, 0.19387, 0.17170, 0.15204, 0.13470, 0.11940, 0.10590, 0.09403, 0.08390, 0.07567, 0.06940, 0.06504, 0.06220, 0.06044, 0.05940, 0.05878, 0.05840, 0.05812, 0.05790, 0.05775, 0.05790, 0.05855, 0.05960, 0.06093, 0.06250, 0.06428, 0.06610, 0.06780, 0.06950, 0.07121, 0.07240, 0.07356, 0.07471, 0.07586, 0.07704, 0.07823, 0.07943, 0.08066, 0.08190, 0.08316, 0.08443, 0.08572, 0.08704, 0.08837, 0.08971, 0.09108, 0.09247, 0.09247, 0.09247, 0.09247, 0.09247, 0.09247, 0.09247, 0.09247, 0.09247, 0.09247, 0.09247 } },
	{ 95, 360.0, 830.0, 1.0, { 0.40857, 0.40857, 0.40857, 0.40857, 0.40857, 0.42408, 0.43973, 0.45551, 0.47190, 0.48883, 0.51160, 0.53114, 0.54880, 0.56754, 0.58620, 0.60333, 0.62040, 0.63907, 0.65890, 0.67891, 0.69800, 0.71523, 0.73050, 0.74379, 0.75470, 0.76285, 0.76840, 0.77166, 0.77300, 0.77289, 0.77220, 0.77164, 0.77090, 0.76950, 0.76730, 0.76421, 0.76000, 0.75437, 0.74680, 0.73691, 0.72520, 0.71226, 0.69820, 0.68280, 0.66490, 0.64336, 0.61800, 0.58959, 0.56170, 0.53802, 0.51990, 0.50780, 0.50090, 0.49798, 0.49760, 0.49839, 0.49950, 0.50065, 0.50320, 0.50874, 0.51800, 0.53104, 0.54610, 0.56122, 0.57550, 0.58862, 0.60150, 0.61435, 0.62320, 0.63169, 0.63982, 0.64787, 0.65583, 0.66371, 0.67150, 0.67920, 0.68680, 0.69430, 0.70170, 0.70899, 0.71618, 0.72326, 0.73023, 0.73709, 0.74383, 0.74383, 0.74383, 0.74383, 0.74383, 0.74383, 0.74383, 0.74383, 0.74383, 0.74383, 0.74383 } },
	{ 95, 360.0, 830.0, 1.0, { 0.30073, 0.30073, 0.30073, 0.30073, 0.30073, 0.32661, 0.34992, 0.37230, 0.39416, 0.41553, 0.43632, 0.45652, 0.47668, 0.49725, 0.51775, 0.53774, 0.55789, 0.57898, 0.60125, 0.62458, 0.64804, 0.67049, 0.69066, 0.70754, 0.72132, 0.73246, 0.74115, 0.74752, 0.75151, 0.75312, 0.75259, 0.75010, 0.74537, 0.73810, 0.72853, 0.71689, 0.70284, 0.68608, 0.66694, 0.64607, 0.62480, 0.60444, 0.58562, 0.56856, 0.55242, 0.53611, 0.51865, 0.49975, 0.48196, 0.46791, 0.45788, 0.45146, 0.44792, 0.44644, 0.44631, 0.44699, 0.44862, 0.45160, 0.45656, 0.46396, 0.47319, 0.48341, 0.49376, 0.50347, 0.51212, 0.51931, 0.52450, 0.52735, 0.52854, 0.52899, 0.52954, 0.53109, 0.53478, 0.54158, 0.55168, 0.56548, 0.58514, 0.61045, 0.62991, 0.63369, 0.62991, 0.62865, 0.62991, 0.63117, 0.62991, 0.62991, 0.62991, 0.62991, 0.62991, 0.62991, 0.62991, 0.62991, 0.62991, 0.62991, 0.62991 } },
	{ 95, 360.0, 830.0, 1.0, { 0.15521, 0.15521, 0.15521, 0.15521, 0.15521, 0.17199, 0.19018, 0.20980, 0.23150, 0.25503, 0.28201, 0.29616, 0.30297, 0.31130, 0.32158, 0.33253, 0.34437, 0.35808, 0.37612, 0.40032, 0.42836, 0.45815, 0.49262, 0.53399, 0.57658, 0.61368, 0.64244, 0.66137, 0.67063, 0.67106, 0.66461, 0.65332, 0.63836, 0.62052, 0.59994, 0.57664, 0.55071, 0.52241, 0.49259, 0.46218, 0.43191, 0.40231, 0.37328, 0.34465, 0.31651, 0.28904, 0.26263, 0.23787, 0.21609, 0.19856, 0.18547, 0.17653, 0.17054, 0.16621, 0.16281, 0.15981, 0.15704, 0.15455, 0.15295, 0.15286, 0.15443, 0.15772, 0.16295, 0.16994, 0.17687, 0.18216, 0.18699, 0.19243, 0.19647, 0.20045, 0.20436, 0.20833, 0.21235, 0.21642, 0.22056, 0.22475, 0.22899, 0.23330, 0.23765, 0.24207, 0.24654, 0.25106, 0.25564, 0.26027, 0.26496, 0.26496, 0.26496, 0.26496, 0.26496, 0.26496, 0.26496, 0.26496, 0.26496, 0.26496, 0.26496 } },
	{ 95, 360.0, 830.0, 1.0, { 0.13470, 0.13470, 0.13470, 0.13470, 0.13470, 0.13497, 0.13470, 0.13370, 0.13470, 0.13984, 0.14600, 0.14999, 0.15360, 0.15934, 0.16770, 0.17857, 0.19160, 0.20654, 0.22380, 0.24394, 0.26740, 0.29427, 0.32340, 0.35336, 0.38270, 0.40983, 0.43250, 0.44870, 0.45800, 0.46055, 0.45710, 0.44861, 0.43620, 0.42088, 0.40300, 0.38264, 0.35940, 0.33313, 0.30510, 0.27694, 0.25030, 0.22651, 0.20550, 0.18674, 0.16940, 0.15247, 0.13470, 0.11572, 0.09900, 0.08790, 0.08160, 0.07827, 0.07630, 0.07448, 0.07290, 0.07173, 0.07010, 0.06692, 0.06130, 0.05304, 0.04440, 0.03781, 0.03390, 0.03263, 0.03320, 0.03593, 0.04650, 0.06786, 0.08660, 0.09028, 0.08660, 0.08561, 0.08660, 0.08686, 0.08660, 0.08653, 0.08660, 0.08662, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660, 0.08660 } },
	{ 95, 360.0, 830.0, 1.0, { 0.10606, 0.10606, 0.10606, 0.10606, 0.10606, 0.11035, 0.11514, 0.12020, 0.12531, 0.13033, 0.13550, 0.14114, 0.14747, 0.15464, 0.16259, 0.17128, 0.18103, 0.19242, 0.20673, 0.22518, 0.24804, 0.27486, 0.30339, 0.33116, 0.35681, 0.37901, 0.39557, 0.40468, 0.40690, 0.40341, 0.39562, 0.38483, 0.37179, 0.35707, 0.34094, 0.32357, 0.30490, 0.28488, 0.26380, 0.24215, 0.22077, 0.20033, 0.18041, 0.16051, 0.14091, 0.12209, 0.10463, 0.08914, 0.07639, 0.06696, 0.06048, 0.05632, 0.05385, 0.05246, 0.05170, 0.05122, 0.05094, 0.05089, 0.05139, 0.05276, 0.05499, 0.05802, 0.06171, 0.06589, 0.07020, 0.07429, 0.07806, 0.08141, 0.08396, 0.08534, 0.08554, 0.08476, 0.08362, 0.08285, 0.08317, 0.08202, 0.08147, 0.08092, 0.08037, 0.07983, 0.07929, 0.07876, 0.07822, 0.07769, 0.07717, 0.07717, 0.07717, 0.07717, 0.07717, 0.07717, 0.07717, 0.07717, 0.07717, 0.07717, 0.07717 } },
	{ 95, 360.0, 830.0, 1.0, { 0.18471, 0.18471, 0.18471, 0.18471, 0.18471, 0.23474, 0.29345, 0.35993, 0.43450, 0.51357, 0.60420, 0.65038, 0.66450, 0.67117, 0.67560, 0.67830, 0.68050, 0.68332, 0.68690, 0.69119, 0.69650, 0.70332, 0.71260, 0.72513, 0.74050, 0.75796, 0.77650, 0.79473, 0.81000, 0.81990, 0.82440, 0.82407, 0.81960, 0.81160, 0.80030, 0.78591, 0.76890, 0.74963, 0.72790, 0.70380, 0.67930, 0.65657, 0.63660, 0.61990, 0.60610, 0.59451, 0.58410, 0.57394, 0.56400, 0.55461, 0.54680, 0.54140, 0.53770, 0.53473, 0.53190, 0.52908, 0.52750, 0.52840, 0.53160, 0.53599, 0.53820, 0.53583, 0.53260, 0.53260, 0.53530, 0.53839, 0.53720, 0.52892, 0.52050, 0.51216, 0.50406, 0.49595, 0.48784, 0.47974, 0.47166, 0.46358, 0.45553, 0.44749, 0.43949, 0.43152, 0.42358, 0.41568, 0.40783, 0.40002, 0.39226, 0.39226, 0.39226, 0.39226, 0.39226, 0.39226, 0.39226, 0.39226, 0.39226, 0.39226, 0.39226 } },
	{ 95, 360.0, 830.0, 1.0, { 0.43710, 0.43710, 0.43710, 0.43710, 0.43710, 0.47272, 0.50312, 0.53129, 0.55735, 0.58083, 0.60183, 0.62063, 0.63777, 0.65381, 0.66908, 0.68370, 0.69721, 0.70918, 0.71973, 0.72912, 0.73751, 0.74500, 0.75159, 0.75729, 0.76220, 0.76636, 0.76953, 0.77145, 0.77193, 0.77082, 0.76796, 0.76324, 0.75680, 0.74895, 0.74036, 0.73159, 0.72237, 0.71241, 0.70206, 0.69189, 0.68268, 0.67515, 0.66956, 0.66589, 0.66336, 0.66107, 0.65854, 0.65553, 0.65249, 0.65002, 0.64874, 0.64903, 0.65032, 0.65199, 0.65430, 0.65764, 0.66202, 0.66738, 0.67372, 0.68090, 0.68803, 0.69419, 0.69911, 0.70280, 0.70579, 0.70842, 0.70991, 0.70960, 0.70853, 0.70807, 0.70918, 0.71261, 0.71859, 0.72710, 0.73766, 0.75043, 0.76850, 0.79265, 0.81153, 0.81521, 0.81153, 0.81030, 0.81153, 0.81276, 0.81153, 0.81153, 0.81153, 0.81153, 0.81153, 0.81153, 0.81153, 0.81153, 0.81153, 0.81153, 0.81153 } },
	{ 95, 360.0, 830.0, 1.0, { 0.28873, 0.28873, 0.28873, 0.28873, 0.28873, 0.30831, 0.32860, 0.34955, 0.36977, 0.39383, 0.41611, 0.43654, 0.45502, 0.47151, 0.48618, 0.49915, 0.51039, 0.51992, 0.52806, 0.53497, 0.53987, 0.54220, 0.54344, 0.54563, 0.55128, 0.56200, 0.57549, 0.58884, 0.60070, 0.61031, 0.61757, 0.62184, 0.61965, 0.60824, 0.59038, 0.56965, 0.54731, 0.52399, 0.50027, 0.47704, 0.45660, 0.44088, 0.42892, 0.41930, 0.41175, 0.40621, 0.40242, 0.39998, 0.39815, 0.39634, 0.39478, 0.39387, 0.39388, 0.39495, 0.39686, 0.39943, 0.40282, 0.40710, 0.41175, 0.41604, 0.41909, 0.42028, 0.42008, 0.41917, 0.41780, 0.41607, 0.41393, 0.41127, 0.40798, 0.40562, 0.40289, 0.40017, 0.39746, 0.39475, 0.39205, 0.38936, 0.38667, 0.38399, 0.38132, 0.37865, 0.37599, 0.37334, 0.37070, 0.36806, 0.36544, 0.36544, 0.36544, 0.36544, 0.36544, 0.36544, 0.36544, 0.36544, 0.36544, 0.36544, 0.36544 } },
	{ 95, 360.0, 830.0, 1.0, { 0.17009, 0.17009, 0.17009, 0.17009, 0.17009, 0.17930, 0.18890, 0.19888, 0.20849, 0.22043, 0.23133, 0.24199, 0.25317, 0.26551, 0.27899, 0.29357, 0.30976, 0.32796, 0.34749, 0.36733, 0.38621, 0.40281, 0.41600, 0.42493, 0.42990, 0.43128, 0.42890, 0.42262, 0.41302, 0.40101, 0.38820, 0.37623, 0.36636, 0.35908, 0.35246, 0.34451, 0.33558, 0.32660, 0.31870, 0.31303, 0.31076, 0.31238, 0.31572, 0.31825, 0.31870, 0.31651, 0.31274, 0.30881, 0.30579, 0.30449, 0.30480, 0.30627, 0.30778, 0.30843, 0.30877, 0.30964, 0.31175, 0.31579, 0.32267, 0.33330, 0.34848, 0.36893, 0.39515, 0.42758, 0.46663, 0.51216, 0.56194, 0.61321, 0.66321, 0.70794, 0.74967, 0.78723, 0.82050, 0.84956, 0.87463, 0.89604, 0.91415, 0.92936, 0.94204, 0.95256, 0.96125, 0.96840, 0.97427, 0.97907, 0.98299, 0.98299, 0.98299, 0.98299, 0.98299, 0.98299, 0.98299, 0.98299, 0.98299, 0.98299, 0.98299 } },
	{ 95, 360.0, 830.0, 1.0, { 0.02353, 0.02353, 0.02353, 0.02353, 0.02353, 0.02714, 0.03128, 0.03604, 0.03956, 0.04880, 0.05538, 0.05973, 0.06229, 0.06356, 0.06429, 0.06519, 0.06650, 0.06835, 0.07080, 0.07393, 0.07801, 0.08318, 0.08883, 0.09427, 0.09914, 0.10306, 0.10525, 0.10508, 0.10285, 0.09908, 0.09414, 0.08838, 0.08232, 0.07644, 0.07090, 0.06586, 0.06169, 0.05878, 0.05728, 0.05732, 0.05899, 0.06227, 0.06670, 0.07168, 0.07651, 0.08049, 0.08302, 0.08366, 0.08262, 0.08028, 0.07721, 0.07397, 0.07080, 0.06788, 0.06539, 0.06346, 0.06189, 0.06055, 0.05979, 0.05995, 0.06109, 0.06312, 0.06589, 0.06936, 0.07381, 0.07946, 0.08592, 0.09264, 0.09904, 0.10694, 0.11492, 0.12341, 0.13244, 0.14202, 0.15216, 0.16290, 0.17424, 0.18619, 0.19876, 0.21197, 0.22580, 0.24026, 0.25534, 0.27102, 0.28730, 0.28730, 0.28730, 0.28730, 0.28730, 0.28730, 0.28730, 0.28730, 0.28730, 0.28730, 0.28730 } },
	{ 95, 360.0, 830.0, 1.0, { 0.16989, 0.16989, 0.16989, 0.16989, 0.16989, 0.18211, 0.19500, 0.20857, 0.22326, 0.23893, 0.25743, 0.26842, 0.27548, 0.28427, 0.29488, 0.30635, 0.31906, 0.33377, 0.35150, 0.37310, 0.39842, 0.42668, 0.45565, 0.48271, 0.50514, 0.52064, 0.52858, 0.52895, 0.52250, 0.51016, 0.49295, 0.47185, 0.44761, 0.42095, 0.39281, 0.36403, 0.33487, 0.30556, 0.27667, 0.24893, 0.22344, 0.20115, 0.18221, 0.16643, 0.15310, 0.14145, 0.13099, 0.12147, 0.11320, 0.10657, 0.10171, 0.09859, 0.09675, 0.09567, 0.09500, 0.09454, 0.09440, 0.09475, 0.09577, 0.09754, 0.09981, 0.10231, 0.10484, 0.10725, 0.10932, 0.11084, 0.11159, 0.11151, 0.11115, 0.11079, 0.11044, 0.11008, 0.10973, 0.10938, 0.10903, 0.10869, 0.10834, 0.10799, 0.10765, 0.10730, 0.10696, 0.10662, 0.10628, 0.10594, 0.10560, 0.10560, 0.10560, 0.10560, 0.10560, 0.10560, 0.10560, 0.10560, 0.10560, 0.10560, 0.10560 } },
	{ 95, 360.0, 830.0, 1.0, { 0.00602, 0.00602, 0.00602, 0.00602, 0.00602, 0.00905, 0.01360, 0.02039, 0.02949, 0.04609, 0.06811, 0.09024, 0.10713, 0.11533, 0.11875, 0.12207, 0.12550, 0.12850, 0.13206, 0.13733, 0.14467, 0.15417, 0.16581, 0.17902, 0.19103, 0.19880, 0.20027, 0.19443, 0.18339, 0.16987, 0.15589, 0.14313, 0.13275, 0.12537, 0.12004, 0.11552, 0.11101, 0.10609, 0.10157, 0.09857, 0.09810, 0.10068, 0.10475, 0.10827, 0.10922, 0.10612, 0.09949, 0.09047, 0.08072, 0.07169, 0.06355, 0.05616, 0.04935, 0.04295, 0.03694, 0.03133, 0.02631, 0.02207, 0.01867, 0.01613, 0.01430, 0.01301, 0.01211, 0.01150, 0.01112, 0.01094, 0.01092, 0.01106, 0.01132, 0.01139, 0.01154, 0.01169, 0.01184, 0.01200, 0.01216, 0.01232, 0.01248, 0.01265, 0.01281, 0.01298, 0.01315, 0.01332, 0.01350, 0.01368, 0.01386, 0.01386, 0.01386, 0.01386, 0.01386, 0.01386, 0.01386, 0.01386, 0.01386, 0.01386, 0.01386 } },
	{ 95, 360.0, 830.0, 1.0, { 0.09963, 0.09963, 0.09963, 0.09963, 0.09963, 0.10544, 0.11155, 0.11796, 0.12490, 0.13231, 0.14130, 0.14721, 0.15180, 0.15780, 0.16550, 0.17467, 0.18540, 0.19798, 0.21300, 0.23107, 0.25240, 0.27660, 0.30110, 0.32295, 0.33970, 0.34943, 0.35190, 0.34737, 0.33650, 0.32016, 0.29980, 0.27692, 0.25270, 0.22820, 0.20430, 0.18174, 0.16090, 0.14202, 0.12520, 0.11054, 0.09830, 0.08866, 0.08130, 0.07575, 0.07150, 0.06811, 0.06530, 0.06290, 0.06090, 0.05933, 0.05820, 0.05751, 0.05710, 0.05680, 0.05660, 0.05650, 0.05650, 0.05658, 0.05680, 0.05720, 0.05770, 0.05822, 0.05870, 0.05910, 0.05940, 0.05961, 0.05970, 0.05967, 0.05960, 0.05953, 0.05947, 0.05940, 0.05933, 0.05927, 0.05920, 0.05914, 0.05907, 0.05900, 0.05894, 0.05887, 0.05881, 0.05874, 0.05868, 0.05861, 0.05854, 0.05854, 0.05854, 0.05854, 0.05854, 0.05854, 0.05854, 0.05854, 0.05854, 0.05854, 0.05854 } },
	{ 95, 360.0, 830.0, 1.0, { 0.10927, 0.10927, 0.10927, 0.10927, 0.10927, 0.11775, 0.12680, 0.13644, 0.14700, 0.15842, 0.17220, 0.18102, 0.18720, 0.19467, 0.20420, 0.21573, 0.22910, 0.24431, 0.26200, 0.28281, 0.30680, 0.33330, 0.35940, 0.38192, 0.39890, 0.40878, 0.41030, 0.40282, 0.38790, 0.36754, 0.34330, 0.31661, 0.28890, 0.26145, 0.23510, 0.21043, 0.18760, 0.16666, 0.14770, 0.13090, 0.11670, 0.10545, 0.09680, 0.09021, 0.08510, 0.08095, 0.07740, 0.07421, 0.07140, 0.06907, 0.06740, 0.06650, 0.06610, 0.06589, 0.06570, 0.06547, 0.06540, 0.06569, 0.06630, 0.06714, 0.06800, 0.06869, 0.06910, 0.06920, 0.06910, 0.06891, 0.06860, 0.06816, 0.06780, 0.06745, 0.06712, 0.06679, 0.06646, 0.06613, 0.06580, 0.06548, 0.06516, 0.06483, 0.06451, 0.06419, 0.06387, 0.06356, 0.06324, 0.06293, 0.06262, 0.06262, 0.06262, 0.06262, 0.06262, 0.06262, 0.06262, 0.06262, 0.06262, 0.06262, 0.06262 } },
	{ 95, 360.0, 830.0, 1.0, { 0.22149, 0.22149, 0.22149, 0.22149, 0.22149, 0.23594, 0.25104, 0.26676, 0.28245, 0.30036, 0.31769, 0.33467, 0.35149, 0.36826, 0.38457, 0.39997, 0.41415, 0.42690, 0.43836, 0.44837, 0.45536, 0.45810, 0.45804, 0.45745, 0.45897, 0.46435, 0.47144, 0.47732, 0.47999, 0.47816, 0.47257, 0.46378, 0.44970, 0.42856, 0.40271, 0.37520, 0.34809, 0.32279, 0.29914, 0.27704, 0.25793, 0.24313, 0.23196, 0.22349, 0.21784, 0.21519, 0.21485, 0.21593, 0.21763, 0.21934, 0.22104, 0.22281, 0.22444, 0.22571, 0.22681, 0.22807, 0.23021, 0.23373, 0.23804, 0.24224, 0.24535, 0.24669, 0.24669, 0.24599, 0.24484, 0.24338, 0.24154, 0.23921, 0.23629, 0.23428, 0.23192, 0.22958, 0.22726, 0.22495, 0.22267, 0.22039, 0.21814, 0.21590, 0.21368, 0.21147, 0.20928, 0.20711, 0.20495, 0.20282, 0.20069, 0.20069, 0.20069, 0.20069, 0.20069, 0.20069, 0.20069, 0.20069, 0.20069, 0.20069, 0.20069 } },
	{ 95, 360.0, 830.0, 1.0, { 0.13225, 0.13225, 0.13225, 0.13225, 0.13225, 0.13933, 0.14672, 0.15443, 0.16220, 0.17098, 0.17961, 0.18831, 0.19733, 0.20679, 0.21630, 0.22540, 0.23381, 0.24136, 0.24814, 0.25399, 0.25741, 0.25724, 0.25494, 0.25267, 0.25277, 0.25679, 0.26298, 0.26902, 0.27359, 0.27587, 0.27606, 0.27413, 0.26802, 0.25609, 0.24031, 0.22327, 0.20640, 0.19064, 0.17600, 0.16260, 0.15168, 0.14438, 0.14004, 0.13802, 0.13932, 0.14481, 0.15313, 0.16242, 0.17116, 0.17820, 0.18373, 0.18816, 0.19156, 0.19396, 0.19558, 0.19682, 0.19847, 0.20117, 0.20444, 0.20756, 0.20980, 0.21071, 0.21063, 0.21003, 0.20908, 0.20784, 0.20630, 0.20442, 0.20218, 0.20051, 0.19863, 0.19676, 0.19491, 0.19307, 0.19124, 0.18942, 0.18762, 0.18584, 0.18406, 0.18230, 0.18055, 0.17882, 0.17709, 0.17538, 0.17369, 0.17369, 0.17369, 0.17369, 0.17369, 0.17369, 0.17369, 0.17369, 0.17369, 0.17369, 0.17369 } },
	{ 95, 360.0, 830.0, 1.0, { 0.31794, 0.31794, 0.31794, 0.31794, 0.31794, 0.34608, 0.37534, 0.40554, 0.43537, 0.46845, 0.49988, 0.52990, 0.55872, 0.58644, 0.61261, 0.63672, 0.65857, 0.67805, 0.69495, 0.70890, 0.71896, 0.72438, 0.72576, 0.72392, 0.71916, 0.71161, 0.70123, 0.68800, 0.67187, 0.65288, 0.63147, 0.60796, 0.58201, 0.55343, 0.52327, 0.49283, 0.46320, 0.43523, 0.40899, 0.38458, 0.36303, 0.34528, 0.33078, 0.31880, 0.30924, 0.30208, 0.29698, 0.29355, 0.29142, 0.29027, 0.28997, 0.29041, 0.29121, 0.29206, 0.29307, 0.29450, 0.29688, 0.30057, 0.30502, 0.30941, 0.31264, 0.31392, 0.31378, 0.31299, 0.31182, 0.31034, 0.30852, 0.30629, 0.30358, 0.30160, 0.29933, 0.29708, 0.29484, 0.29260, 0.29038, 0.28817, 0.28597, 0.28377, 0.28159, 0.27942, 0.27726, 0.27511, 0.27296, 0.27083, 0.26871, 0.26871, 0.26871, 0.26871, 0.26871, 0.26871, 0.26871, 0.26871, 0.26871, 0.26871, 0.26871 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03646, 0.03646, 0.03646, 0.03646, 0.03646, 0.03952, 0.04283, 0.04641, 0.05040, 0.05485, 0.06170, 0.06943, 0.07860, 0.09002, 0.10340, 0.11814, 0.13360, 0.14925, 0.16520, 0.18141, 0.19670, 0.20953, 0.21820, 0.22140, 0.21950, 0.21333, 0.20390, 0.19206, 0.17800, 0.16197, 0.14520, 0.12916, 0.11530, 0.10467, 0.09690, 0.09100, 0.08500, 0.07718, 0.06770, 0.05734, 0.04730, 0.03901, 0.03430, 0.03478, 0.04080, 0.05134, 0.06120, 0.06540, 0.06410, 0.05901, 0.05300, 0.04859, 0.04590, 0.04443, 0.04370, 0.04331, 0.04310, 0.04314, 0.04410, 0.04681, 0.05200, 0.06016, 0.07080, 0.08316, 0.09630, 0.10963, 0.12420, 0.14004, 0.15140, 0.16305, 0.17496, 0.18755, 0.20083, 0.21480, 0.22946, 0.24480, 0.26083, 0.27752, 0.29486, 0.31280, 0.33133, 0.35040, 0.36995, 0.38994, 0.41031, 0.41031, 0.41031, 0.41031, 0.41031, 0.41031, 0.41031, 0.41031, 0.41031, 0.41031, 0.41031 } },
	{ 95, 360.0, 830.0, 1.0, { 0.37756, 0.37756, 0.37756, 0.37756, 0.37756, 0.39319, 0.40905, 0.42509, 0.44022, 0.45815, 0.47443, 0.48877, 0.50091, 0.51073, 0.51874, 0.52556, 0.53162, 0.53716, 0.54203, 0.54598, 0.54903, 0.55122, 0.55254, 0.55300, 0.55264, 0.55149, 0.54955, 0.54681, 0.54326, 0.53894, 0.53389, 0.52809, 0.52121, 0.51298, 0.50390, 0.49461, 0.48566, 0.47742, 0.46979, 0.46263, 0.45609, 0.45034, 0.44537, 0.44115, 0.43795, 0.43598, 0.43496, 0.43460, 0.43496, 0.43610, 0.43774, 0.43955, 0.44114, 0.44227, 0.44310, 0.44397, 0.44537, 0.44767, 0.45052, 0.45332, 0.45526, 0.45577, 0.45536, 0.45475, 0.45433, 0.45436, 0.45475, 0.45534, 0.45598, 0.45652, 0.45709, 0.45767, 0.45825, 0.45883, 0.45941, 0.45998, 0.46056, 0.46114, 0.46172, 0.46230, 0.46288, 0.46346, 0.46404, 0.46462, 0.46519, 0.46519, 0.46519, 0.46519, 0.46519, 0.46519, 0.46519, 0.46519, 0.46519, 0.46519, 0.46519 } },
	{ 95, 360.0, 830.0, 1.0, { 0.16160, 0.16160, 0.16160, 0.16160, 0.16160, 0.18016, 0.20034, 0.22217, 0.23998, 0.27370, 0.29807, 0.32072, 0.34930, 0.38927, 0.43737, 0.48744, 0.53051, 0.55904, 0.57419, 0.57916, 0.57668, 0.56907, 0.55741, 0.54267, 0.52663, 0.51095, 0.49575, 0.48059, 0.46418, 0.44544, 0.42486, 0.40323, 0.38078, 0.35769, 0.33451, 0.31157, 0.28794, 0.26281, 0.23701, 0.21187, 0.18905, 0.16994, 0.15459, 0.14279, 0.13464, 0.13012, 0.12858, 0.12914, 0.13086, 0.13275, 0.13364, 0.13271, 0.13047, 0.12777, 0.12550, 0.12447, 0.12520, 0.12807, 0.13295, 0.13949, 0.14695, 0.15462, 0.16234, 0.17089, 0.18428, 0.20617, 0.23561, 0.27055, 0.30889, 0.34890, 0.39166, 0.43615, 0.48169, 0.52754, 0.57293, 0.61712, 0.65946, 0.69940, 0.73652, 0.77057, 0.80140, 0.82900, 0.85348, 0.87497, 0.89371, 0.89371, 0.89371, 0.89371, 0.89371, 0.89371, 0.89371, 0.89371, 0.89371, 0.89371, 0.89371 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03388, 0.03388, 0.03388, 0.03388, 0.03388, 0.03483, 0.03581, 0.03682, 0.03783, 0.03893, 0.04001, 0.04112, 0.04230, 0.04360, 0.04508, 0.04678, 0.04885, 0.05133, 0.05391, 0.05623, 0.05798, 0.05898, 0.05937, 0.05937, 0.05927, 0.05923, 0.05888, 0.05785, 0.05629, 0.05448, 0.05282, 0.05164, 0.05083, 0.05013, 0.04915, 0.04764, 0.04597, 0.04463, 0.04398, 0.04429, 0.04547, 0.04724, 0.04885, 0.04958, 0.04934, 0.04829, 0.04686, 0.04551, 0.04448, 0.04393, 0.04378, 0.04390, 0.04408, 0.04420, 0.04428, 0.04442, 0.04478, 0.04555, 0.04706, 0.04965, 0.05351, 0.05890, 0.06642, 0.07676, 0.09035, 0.10762, 0.12907, 0.15524, 0.18665, 0.22083, 0.26013, 0.30371, 0.35111, 0.40165, 0.45437, 0.50813, 0.56171, 0.61388, 0.66357, 0.70988, 0.75220, 0.79016, 0.82368, 0.85284, 0.87789, 0.87789, 0.87789, 0.87789, 0.87789, 0.87789, 0.87789, 0.87789, 0.87789, 0.87789, 0.87789 } },
	{ 95, 360.0, 830.0, 1.0, { 0.02919, 0.02919, 0.02919, 0.02919, 0.02919, 0.03355, 0.03853, 0.04422, 0.05090, 0.05866, 0.06910, 0.07815, 0.08670, 0.09650, 0.10760, 0.11962, 0.13220, 0.14512, 0.15870, 0.17307, 0.18690, 0.19843, 0.20550, 0.20636, 0.20120, 0.19094, 0.17740, 0.16234, 0.14650, 0.13040, 0.11480, 0.10054, 0.08860, 0.07973, 0.07360, 0.06947, 0.06600, 0.06197, 0.05720, 0.05188, 0.04670, 0.04255, 0.04050, 0.04145, 0.04540, 0.05156, 0.05680, 0.05823, 0.05620, 0.05198, 0.04730, 0.04369, 0.04130, 0.03997, 0.03950, 0.03970, 0.04030, 0.04114, 0.04260, 0.04518, 0.04930, 0.05525, 0.06270, 0.07116, 0.08000, 0.08885, 0.09840, 0.10867, 0.11600, 0.12346, 0.13106, 0.13904, 0.14743, 0.15623, 0.16546, 0.17512, 0.18522, 0.19576, 0.20675, 0.21819, 0.23008, 0.24241, 0.25519, 0.26840, 0.28204, 0.28204, 0.28204, 0.28204, 0.28204, 0.28204, 0.28204, 0.28204, 0.28204, 0.28204, 0.28204 } },
	{ 95, 360.0, 830.0, 1.0, { 0.21199, 0.21199, 0.21199, 0.21199, 0.21199, 0.22951, 0.24803, 0.26752, 0.28860, 0.31111, 0.33960, 0.36066, 0.37770, 0.39675, 0.41750, 0.43808, 0.45680, 0.47192, 0.48120, 0.48293, 0.47790, 0.46766, 0.45430, 0.43959, 0.42360, 0.40603, 0.38690, 0.36621, 0.34350, 0.31848, 0.29190, 0.26503, 0.24020, 0.21929, 0.20150, 0.18521, 0.16840, 0.14935, 0.12800, 0.10516, 0.08340, 0.06563, 0.05420, 0.05104, 0.05680, 0.07068, 0.08720, 0.10057, 0.10850, 0.11045, 0.10940, 0.10857, 0.10870, 0.10986, 0.11190, 0.11481, 0.11940, 0.12646, 0.13580, 0.14672, 0.15740, 0.16579, 0.17020, 0.16964, 0.16570, 0.15994, 0.15120, 0.13917, 0.12970, 0.12095, 0.11294, 0.10540, 0.09831, 0.09165, 0.08540, 0.07953, 0.07404, 0.06889, 0.06408, 0.05959, 0.05539, 0.05147, 0.04781, 0.04440, 0.04123, 0.04123, 0.04123, 0.04123, 0.04123, 0.04123, 0.04123, 0.04123, 0.04123, 0.04123, 0.04123 } },
	{ 95, 360.0, 830.0, 1.0, { 0.02805, 0.02805, 0.02805, 0.02805, 0.02805, 0.03537, 0.05422, 0.07830, 0.10593, 0.13046, 0.15591, 0.18140, 0.20586, 0.22735, 0.24649, 0.26234, 0.27425, 0.28313, 0.29062, 0.29480, 0.29517, 0.29327, 0.28875, 0.28176, 0.27334, 0.26346, 0.25141, 0.23809, 0.22535, 0.21297, 0.20098, 0.18934, 0.17750, 0.16598, 0.15478, 0.14464, 0.13534, 0.12632, 0.11734, 0.10870, 0.10113, 0.09491, 0.08967, 0.08526, 0.08169, 0.07914, 0.07764, 0.07707, 0.07730, 0.07825, 0.08023, 0.08357, 0.08886, 0.09681, 0.10784, 0.12216, 0.13992, 0.16133, 0.18574, 0.21193, 0.23877, 0.26496, 0.28852, 0.30737, 0.32249, 0.33247, 0.33782, 0.33897, 0.33667, 0.33133, 0.32425, 0.31788, 0.31165, 0.30494, 0.29854, 0.29222, 0.28597, 0.27981, 0.27373, 0.26774, 0.26182, 0.25589, 0.25007, 0.24437, 0.23878, 0.23878, 0.23878, 0.23878, 0.23878, 0.23878, 0.23878, 0.23878, 0.23878, 0.23878, 0.23878 } },
	{ 95, 360.0, 830.0, 1.0, { 0.12947, 0.12947, 0.12947, 0.12947, 0.12947, 0.14167, 0.15481, 0.16893, 0.18120, 0.20172, 0.21754, 0.23408, 0.25676, 0.28940, 0.32944, 0.37206, 0.40977, 0.43567, 0.44790, 0.44654, 0.43439, 0.41472, 0.39001, 0.36265, 0.33550, 0.31099, 0.28943, 0.27027, 0.25170, 0.23219, 0.21258, 0.19421, 0.17803, 0.16458, 0.15330, 0.14323, 0.13285, 0.12098, 0.10823, 0.09570, 0.08469, 0.07633, 0.07069, 0.06766, 0.06722, 0.06932, 0.07348, 0.07900, 0.08479, 0.08966, 0.09254, 0.09274, 0.09095, 0.08820, 0.08559, 0.08415, 0.08450, 0.08701, 0.09155, 0.09774, 0.10495, 0.11256, 0.12024, 0.12850, 0.14079, 0.16039, 0.18676, 0.21842, 0.25388, 0.29220, 0.33405, 0.37868, 0.42547, 0.47364, 0.52230, 0.57054, 0.61747, 0.66232, 0.70442, 0.74331, 0.77869, 0.81043, 0.83857, 0.86323, 0.88465, 0.88465, 0.88465, 0.88465, 0.88465, 0.88465, 0.88465, 0.88465, 0.88465, 0.88465, 0.88465 } },
	{ 95, 360.0, 830.0, 1.0, { 0.56749, 0.56749, 0.56749, 0.56749, 0.56749, 0.57857, 0.58956, 0.60046, 0.60662, 0.62422, 0.63353, 0.63805, 0.64127, 0.64576, 0.65031, 0.65295, 0.65239, 0.64837, 0.64415, 0.64297, 0.64445, 0.64684, 0.64673, 0.64119, 0.63104, 0.61778, 0.60176, 0.58330, 0.56373, 0.54443, 0.52600, 0.50863, 0.49145, 0.47370, 0.45601, 0.43928, 0.42404, 0.41066, 0.39922, 0.38961, 0.38145, 0.37420, 0.36715, 0.35995, 0.35394, 0.35060, 0.35017, 0.35247, 0.35672, 0.36207, 0.36794, 0.37408, 0.38115, 0.39035, 0.40398, 0.42415, 0.45114, 0.48460, 0.52342, 0.56583, 0.60811, 0.64650, 0.67910, 0.70496, 0.72516, 0.74134, 0.75525, 0.76868, 0.78344, 0.79525, 0.80745, 0.81908, 0.83016, 0.84069, 0.85069, 0.86016, 0.86913, 0.87760, 0.88559, 0.89313, 0.90023, 0.90690, 0.91317, 0.91906, 0.92458, 0.92458, 0.92458, 0.92458, 0.92458, 0.92458, 0.92458, 0.92458, 0.92458, 0.92458, 0.92458 } },
	{ 95, 360.0, 830.0, 1.0, { 0.31684, 0.31684, 0.31684, 0.31684, 0.31684, 0.34433, 0.36972, 0.39599, 0.42314, 0.45025, 0.47575, 0.49823, 0.51751, 0.53377, 0.54733, 0.55847, 0.56720, 0.57340, 0.57683, 0.57723, 0.57441, 0.56815, 0.55801, 0.54360, 0.52497, 0.50256, 0.47780, 0.45229, 0.42720, 0.40344, 0.38140, 0.36127, 0.34289, 0.32612, 0.31138, 0.29903, 0.28877, 0.28021, 0.27333, 0.26820, 0.26479, 0.26307, 0.26304, 0.26462, 0.26728, 0.27052, 0.27437, 0.27890, 0.28390, 0.28921, 0.29508, 0.30178, 0.30938, 0.31796, 0.32802, 0.34004, 0.35416, 0.37048, 0.38943, 0.41133, 0.43586, 0.46260, 0.49134, 0.52196, 0.55435, 0.58815, 0.62191, 0.65419, 0.68463, 0.71299, 0.73843, 0.76029, 0.77927, 0.79615, 0.81061, 0.82241, 0.83267, 0.84211, 0.84839, 0.84957, 0.84839, 0.84800, 0.84839, 0.84878, 0.84839, 0.84839, 0.84839, 0.84839, 0.84839, 0.84839, 0.84839, 0.84839, 0.84839, 0.84839, 0.84839 } },
	{ 95, 360.0, 830.0, 1.0, { 0.12067, 0.12067, 0.12067, 0.12067, 0.12067, 0.17602, 0.23298, 0.27447, 0.30442, 0.32357, 0.33682, 0.34501, 0.34763, 0.35069, 0.35500, 0.35271, 0.34519, 0.33814, 0.33335, 0.32585, 0.31343, 0.30296, 0.29343, 0.27907, 0.26179, 0.24644, 0.23007, 0.21652, 0.20356, 0.18994, 0.17481, 0.15838, 0.14046, 0.12315, 0.11069, 0.10618, 0.10310, 0.09997, 0.09749, 0.09121, 0.08211, 0.07440, 0.06799, 0.06738, 0.07171, 0.08124, 0.09024, 0.09630, 0.10233, 0.10517, 0.10010, 0.09536, 0.09978, 0.11368, 0.13297, 0.16167, 0.19616, 0.23898, 0.29173, 0.34345, 0.37835, 0.40403, 0.42472, 0.43600, 0.44674, 0.45591, 0.46076, 0.46165, 0.46531, 0.47145, 0.47513, 0.47881, 0.48250, 0.48618, 0.48987, 0.49356, 0.49725, 0.50094, 0.50463, 0.50832, 0.51201, 0.51570, 0.51938, 0.52307, 0.52675, 0.52675, 0.52675, 0.52675, 0.52675, 0.52675, 0.52675, 0.52675, 0.52675, 0.52675, 0.52675 } },
	{ 95, 360.0, 830.0, 1.0, { 0.70772, 0.70772, 0.70772, 0.70772, 0.70772, 0.71503, 0.72222, 0.72930, 0.73650, 0.74377, 0.75320, 0.76059, 0.76660, 0.77259, 0.77790, 0.78171, 0.78450, 0.78682, 0.78820, 0.78804, 0.78620, 0.78261, 0.77710, 0.76964, 0.76100, 0.75186, 0.74180, 0.73044, 0.71880, 0.70796, 0.69790, 0.68834, 0.67900, 0.66980, 0.66150, 0.65484, 0.64970, 0.64587, 0.64380, 0.64381, 0.64510, 0.64674, 0.64830, 0.64956, 0.65070, 0.65200, 0.65380, 0.65629, 0.65910, 0.66177, 0.66390, 0.66528, 0.66630, 0.66737, 0.66840, 0.66966, 0.67350, 0.68201, 0.69420, 0.70871, 0.72580, 0.74588, 0.76840, 0.79213, 0.81420, 0.83253, 0.84990, 0.86826, 0.88140, 0.89315, 0.90355, 0.91304, 0.92167, 0.92952, 0.93663, 0.94307, 0.94889, 0.95414, 0.95888, 0.96314, 0.96698, 0.97044, 0.97354, 0.97632, 0.97882, 0.97882, 0.97882, 0.97882, 0.97882, 0.97882, 0.97882, 0.97882, 0.97882, 0.97882, 0.97882 } },
	{ 95, 360.0, 830.0, 1.0, { 0.56349, 0.56349, 0.56349, 0.56349, 0.56349, 0.57154, 0.57956, 0.58753, 0.59570, 0.60396, 0.61320, 0.61739, 0.61800, 0.61762, 0.61590, 0.61213, 0.60710, 0.60152, 0.59440, 0.58478, 0.57360, 0.56173, 0.54790, 0.53111, 0.51360, 0.49766, 0.48270, 0.46771, 0.45300, 0.43916, 0.42670, 0.41592, 0.40640, 0.39769, 0.39000, 0.38378, 0.37970, 0.37817, 0.37830, 0.37928, 0.38180, 0.38643, 0.39180, 0.39639, 0.40010, 0.40309, 0.40520, 0.40630, 0.40670, 0.40672, 0.40620, 0.40499, 0.40330, 0.40177, 0.40240, 0.40703, 0.41560, 0.42779, 0.44420, 0.46526, 0.48980, 0.51664, 0.54610, 0.57781, 0.60700, 0.63023, 0.65360, 0.68162, 0.70290, 0.72287, 0.74141, 0.75913, 0.77600, 0.79202, 0.80717, 0.82147, 0.83492, 0.84755, 0.85937, 0.87042, 0.88072, 0.89031, 0.89921, 0.90746, 0.91511, 0.91511, 0.91511, 0.91511, 0.91511, 0.91511, 0.91511, 0.91511, 0.91511, 0.91511, 0.91511 } },
	{ 95, 360.0, 830.0, 1.0, { 0.03375, 0.03375, 0.03375, 0.03375, 0.03375, 0.03555, 0.03745, 0.03943, 0.04138, 0.04379, 0.04605, 0.04830, 0.05071, 0.05334, 0.05597, 0.05829, 0.06004, 0.06105, 0.06153, 0.06171, 0.06143, 0.06054, 0.05935, 0.05833, 0.05816, 0.05926, 0.06103, 0.06267, 0.06351, 0.06315, 0.06193, 0.06027, 0.05806, 0.05512, 0.05170, 0.04814, 0.04456, 0.04107, 0.03781, 0.03496, 0.03285, 0.03177, 0.03156, 0.03212, 0.03394, 0.03742, 0.04208, 0.04716, 0.05170, 0.05496, 0.05726, 0.05907, 0.06044, 0.06132, 0.06183, 0.06214, 0.06252, 0.06316, 0.06391, 0.06457, 0.06500, 0.06512, 0.06500, 0.06478, 0.06451, 0.06423, 0.06391, 0.06349, 0.06292, 0.06257, 0.06214, 0.06171, 0.06128, 0.06085, 0.06043, 0.06001, 0.05959, 0.05918, 0.05877, 0.05836, 0.05795, 0.05755, 0.05715, 0.05675, 0.05635, 0.05635, 0.05635, 0.05635, 0.05635, 0.05635, 0.05635, 0.05635, 0.05635, 0.05635, 0.05635 } },
	{ 95, 360.0, 830.0, 1.0, { 0.01425, 0.01425, 0.01425, 0.01425, 0.01425, 0.01292, 0.01731, 0.02444, 0.03831, 0.04940, 0.06324, 0.07955, 0.09523, 0.11233, 0.13166, 0.14950, 0.15820, 0.16275, 0.17107, 0.17496, 0.17395, 0.17370, 0.17216, 0.16727, 0.16007, 0.15576, 0.14927, 0.14321, 0.13436, 0.12644, 0.12123, 0.11599, 0.11004, 0.10400, 0.09716, 0.09359, 0.09213, 0.09261, 0.09066, 0.08902, 0.08907, 0.08792, 0.08794, 0.08588, 0.08463, 0.08725, 0.09088, 0.09758, 0.10274, 0.10874, 0.10959, 0.11171, 0.11603, 0.12276, 0.12916, 0.13826, 0.14464, 0.15049, 0.16785, 0.19265, 0.20784, 0.22794, 0.25255, 0.26974, 0.28737, 0.30481, 0.31630, 0.32417, 0.33878, 0.34908, 0.36054, 0.37215, 0.38392, 0.39582, 0.40785, 0.41999, 0.43223, 0.44455, 0.45694, 0.46938, 0.48186, 0.49436, 0.50687, 0.51938, 0.53185, 0.53185, 0.53185, 0.53185, 0.53185, 0.53185, 0.53185, 0.53185, 0.53185, 0.53185, 0.53185 } },
	{ 95, 360.0, 830.0, 1.0, { 0.50252, 0.50252, 0.50252, 0.50252, 0.50252, 0.51400, 0.52546, 0.53689, 0.54636, 0.56056, 0.57138, 0.57981, 0.58686, 0.59284, 0.59520, 0.59176, 0.58448, 0.57575, 0.56542, 0.55287, 0.53821, 0.52174, 0.50366, 0.48437, 0.46484, 0.44597, 0.42761, 0.40934, 0.39068, 0.37145, 0.35265, 0.33549, 0.32078, 0.30890, 0.29884, 0.28927, 0.27879, 0.26664, 0.25456, 0.24483, 0.23947, 0.23969, 0.24364, 0.24875, 0.25258, 0.25331, 0.25139, 0.24833, 0.24771, 0.25339, 0.26836, 0.29491, 0.33339, 0.38325, 0.44221, 0.50670, 0.56969, 0.62438, 0.66837, 0.70087, 0.72318, 0.73733, 0.74601, 0.75185, 0.75644, 0.76102, 0.76637, 0.77313, 0.78195, 0.78747, 0.79409, 0.80056, 0.80687, 0.81303, 0.81904, 0.82490, 0.83060, 0.83616, 0.84157, 0.84683, 0.85195, 0.85693, 0.86177, 0.86647, 0.87103, 0.87103, 0.87103, 0.87103, 0.87103, 0.87103, 0.87103, 0.87103, 0.87103, 0.87103, 0.87103 } },
	{ 95, 360.0, 830.0, 1.0, { 0.10613, 0.10613, 0.10613, 0.10613, 0.10613, 0.11586, 0.12635, 0.13764, 0.14596, 0.16476, 0.17753, 0.18653, 0.19401, 0.20193, 0.21099, 0.22114, 0.23055, 0.23700, 0.23849, 0.23384, 0.22469, 0.21322, 0.20076, 0.18829, 0.17614, 0.16448, 0.15340, 0.14290, 0.13255, 0.12197, 0.11120, 0.10054, 0.09085, 0.08296, 0.07705, 0.07288, 0.06921, 0.06482, 0.05957, 0.05375, 0.04826, 0.04403, 0.04140, 0.04111, 0.04617, 0.05927, 0.07933, 0.10352, 0.12550, 0.13984, 0.14804, 0.15301, 0.15638, 0.15927, 0.16204, 0.16492, 0.16820, 0.17215, 0.17684, 0.18200, 0.18627, 0.18830, 0.18796, 0.18571, 0.18319, 0.18185, 0.18120, 0.18026, 0.17803, 0.17751, 0.17632, 0.17513, 0.17395, 0.17278, 0.17161, 0.17045, 0.16929, 0.16814, 0.16700, 0.16586, 0.16473, 0.16361, 0.16249, 0.16138, 0.16027, 0.16027, 0.16027, 0.16027, 0.16027, 0.16027, 0.16027, 0.16027, 0.16027, 0.16027, 0.16027 } },
	{ 95, 360.0, 830.0, 1.0, { 0.42640, 0.42640, 0.42640, 0.42640, 0.42640, 0.44033, 0.45435, 0.46845, 0.47490, 0.50054, 0.51283, 0.51584, 0.51363, 0.50938, 0.50280, 0.49313, 0.48116, 0.46804, 0.45465, 0.44164, 0.42903, 0.41662, 0.40401, 0.39068, 0.37591, 0.35912, 0.34056, 0.32103, 0.30254, 0.28662, 0.27176, 0.25598, 0.23869, 0.22042, 0.20493, 0.19580, 0.19252, 0.19315, 0.19401, 0.19190, 0.18706, 0.18151, 0.18071, 0.18983, 0.20940, 0.23795, 0.27076, 0.30302, 0.33272, 0.35883, 0.38147, 0.40076, 0.41563, 0.42508, 0.42953, 0.42990, 0.42774, 0.42497, 0.42427, 0.42833, 0.43926, 0.45877, 0.48761, 0.52562, 0.56992, 0.61695, 0.66286, 0.70375, 0.73574, 0.77275, 0.80334, 0.83071, 0.85496, 0.87625, 0.89481, 0.91086, 0.92467, 0.93649, 0.94656, 0.95511, 0.96235, 0.96846, 0.97360, 0.97793, 0.98156, 0.98156, 0.98156, 0.98156, 0.98156, 0.98156, 0.98156, 0.98156, 0.98156, 0.98156, 0.98156 } },
	{ 95, 360.0, 830.0, 1.0, { 0.07275, 0.07275, 0.07275, 0.07275, 0.07275, 0.08208, 0.09249, 0.10408, 0.11546, 0.13191, 0.14757, 0.16011, 0.16724, 0.16828, 0.16927, 0.17537, 0.18160, 0.18230, 0.17906, 0.17458, 0.16856, 0.16017, 0.14940, 0.13679, 0.12423, 0.11336, 0.10344, 0.09358, 0.08459, 0.07741, 0.07195, 0.06776, 0.06420, 0.06078, 0.05768, 0.05526, 0.05381, 0.05358, 0.05483, 0.05763, 0.06145, 0.06564, 0.06961, 0.07304, 0.07664, 0.08120, 0.08693, 0.09386, 0.10181, 0.11045, 0.11913, 0.12714, 0.13391, 0.13922, 0.14390, 0.14909, 0.15603, 0.16583, 0.17906, 0.19607, 0.21697, 0.24170, 0.26986, 0.30107, 0.33549, 0.37309, 0.41253, 0.45219, 0.49039, 0.53155, 0.57117, 0.60989, 0.64728, 0.68294, 0.71658, 0.74797, 0.77696, 0.80350, 0.82757, 0.84925, 0.86864, 0.88587, 0.90110, 0.91449, 0.92621, 0.92621, 0.92621, 0.92621, 0.92621, 0.92621, 0.92621, 0.92621, 0.92621, 0.92621, 0.92621 } },
	{ 95, 360.0, 830.0, 1.0, { 0.22178, 0.22178, 0.22178, 0.22178, 0.22178, 0.23127, 0.24103, 0.25108, 0.26170, 0.27271, 0.28380, 0.28598, 0.28220, 0.27711, 0.27030, 0.26063, 0.24900, 0.23663, 0.22390, 0.21107, 0.19870, 0.18718, 0.17590, 0.16417, 0.15200, 0.13979, 0.12870, 0.11977, 0.11260, 0.10642, 0.10020, 0.09321, 0.08600, 0.07945, 0.07430, 0.07114, 0.06990, 0.07019, 0.07090, 0.07109, 0.07110, 0.07175, 0.07450, 0.08084, 0.09190, 0.10811, 0.12760, 0.14791, 0.16650, 0.18130, 0.19230, 0.20011, 0.20570, 0.20999, 0.21340, 0.21626, 0.21890, 0.22164, 0.22460, 0.22768, 0.23010, 0.23137, 0.23290, 0.23606, 0.24010, 0.24393, 0.24720, 0.24973, 0.25120, 0.25261, 0.25397, 0.25533, 0.25670, 0.25807, 0.25945, 0.26083, 0.26222, 0.26362, 0.26501, 0.26642, 0.26782, 0.26923, 0.27065, 0.27207, 0.27350, 0.27350, 0.27350, 0.27350, 0.27350, 0.27350, 0.27350, 0.27350, 0.27350, 0.27350, 0.27350 } },
	{ 95, 360.0, 830.0, 1.0, { 0.02966, 0.02966, 0.02966, 0.02966, 0.02966, 0.03379, 0.04777, 0.07917, 0.13075, 0.17691, 0.22905, 0.28446, 0.34279, 0.40110, 0.45479, 0.50150, 0.53963, 0.56762, 0.58708, 0.60061, 0.60850, 0.61064, 0.60879, 0.60406, 0.59734, 0.59009, 0.58275, 0.57511, 0.56692, 0.55949, 0.55297, 0.54597, 0.53793, 0.52943, 0.52207, 0.51669, 0.51244, 0.50807, 0.50371, 0.50097, 0.50109, 0.50342, 0.50601, 0.50809, 0.51072, 0.51513, 0.52154, 0.52925, 0.53705, 0.54511, 0.55458, 0.56537, 0.57661, 0.58749, 0.59586, 0.60146, 0.60694, 0.61324, 0.61886, 0.62331, 0.62790, 0.63150, 0.63258, 0.63234, 0.63235, 0.62884, 0.62321, 0.61633, 0.60749, 0.59715, 0.58830, 0.58146, 0.57354, 0.56538, 0.55745, 0.54949, 0.54151, 0.53350, 0.52548, 0.51744, 0.50940, 0.50135, 0.49329, 0.48523, 0.47718, 0.47718, 0.47718, 0.47718, 0.47718, 0.47718, 0.47718, 0.47718, 0.47718, 0.47718, 0.47718 } },
	{ 95, 360.0, 830.0, 1.0, { 0.00526, 0.00526, 0.00526, 0.00526, 0.00526, 0.02448, 0.05165, 0.07373, 0.09670, 0.12048, 0.14705, 0.17355, 0.19903, 0.22337, 0.24497, 0.26321, 0.27932, 0.29203, 0.29990, 0.30270, 0.30154, 0.29826, 0.29386, 0.28892, 0.28247, 0.27445, 0.26514, 0.25488, 0.24532, 0.23679, 0.22756, 0.21696, 0.20634, 0.19676, 0.18816, 0.18045, 0.17340, 0.16725, 0.16263, 0.16000, 0.15893, 0.15868, 0.15918, 0.16048, 0.16283, 0.16661, 0.17222, 0.18041, 0.19233, 0.20953, 0.23350, 0.26384, 0.29850, 0.33518, 0.37373, 0.41399, 0.45456, 0.49236, 0.52439, 0.55001, 0.57049, 0.58799, 0.60570, 0.62287, 0.63523, 0.64717, 0.65931, 0.67212, 0.68575, 0.70060, 0.71303, 0.72275, 0.73370, 0.74439, 0.75455, 0.76442, 0.77403, 0.78335, 0.79239, 0.80115, 0.80963, 0.81816, 0.82632, 0.83411, 0.84154, 0.84154, 0.84154, 0.84154, 0.84154, 0.84154, 0.84154, 0.84154, 0.84154, 0.84154, 0.84154 } },
	{ 95, 360.0, 830.0, 1.0, { 0.23077, 0.23077, 0.23077, 0.23077, 0.23077, 0.24352, 0.25636, 0.26920, 0.28195, 0.29460, 0.30754, 0.32130, 0.33662, 0.35357, 0.36925, 0.38105, 0.39054, 0.39975, 0.40839, 0.41539, 0.41877, 0.41719, 0.41283, 0.40771, 0.39971, 0.38682, 0.37157, 0.35700, 0.34356, 0.33108, 0.31952, 0.30897, 0.30003, 0.29316, 0.28788, 0.28365, 0.28066, 0.27914, 0.27880, 0.27935, 0.28107, 0.28419, 0.28827, 0.29295, 0.29889, 0.30683, 0.31670, 0.32832, 0.34189, 0.35764, 0.37558, 0.39567, 0.41783, 0.44204, 0.46840, 0.49679, 0.52606, 0.55502, 0.58337, 0.61088, 0.63663, 0.65968, 0.67984, 0.69709, 0.71148, 0.72321, 0.73309, 0.74189, 0.74957, 0.75581, 0.76013, 0.76226, 0.76306, 0.76369, 0.76530, 0.76570, 0.76658, 0.76746, 0.76833, 0.76920, 0.77007, 0.77094, 0.77181, 0.77267, 0.77353, 0.77353, 0.77353, 0.77353, 0.77353, 0.77353, 0.77353, 0.77353, 0.77353, 0.77353, 0.77353 } },
	{ 95, 360.0, 830.0, 1.0, { 0.06102, 0.06102, 0.06102, 0.06102, 0.06102, 0.06701, 0.07354, 0.08065, 0.08810, 0.09694, 0.10607, 0.11496, 0.12312, 0.13016, 0.13629, 0.14162, 0.14547, 0.14718, 0.14700, 0.14528, 0.14180, 0.13637, 0.12955, 0.12194, 0.11362, 0.10476, 0.09637, 0.08937, 0.08341, 0.07795, 0.07299, 0.06865, 0.06493, 0.06185, 0.05952, 0.05800, 0.05727, 0.05727, 0.05788, 0.05896, 0.06033, 0.06187, 0.06370, 0.06599, 0.06891, 0.07264, 0.07748, 0.08381, 0.09218, 0.10311, 0.11679, 0.13341, 0.15364, 0.17791, 0.20530, 0.23477, 0.26624, 0.29975, 0.33484, 0.37096, 0.40773, 0.44475, 0.48134, 0.51686, 0.55127, 0.58455, 0.61630, 0.64600, 0.67316, 0.70121, 0.72707, 0.75148, 0.77439, 0.79576, 0.81558, 0.83389, 0.85071, 0.86610, 0.88012, 0.89286, 0.90440, 0.91481, 0.92418, 0.93259, 0.94014, 0.94014, 0.94014, 0.94014, 0.94014, 0.94014, 0.94014, 0.94014, 0.94014, 0.94014, 0.94014 } },
	{ 95, 360.0, 830.0, 1.0, { 0.07372, 0.07372, 0.07372, 0.07372, 0.07372, 0.09065, 0.11100, 0.13525, 0.16428, 0.19678, 0.23541, 0.27695, 0.31816, 0.35643, 0.39144, 0.42267, 0.44627, 0.45905, 0.46390, 0.46464, 0.46288, 0.45926, 0.45279, 0.44269, 0.43067, 0.41862, 0.40662, 0.39437, 0.38166, 0.36832, 0.35445, 0.34032, 0.32703, 0.31562, 0.30624, 0.29880, 0.29299, 0.28848, 0.28494, 0.28227, 0.28117, 0.28238, 0.28596, 0.29206, 0.30186, 0.31662, 0.33692, 0.36293, 0.39399, 0.42915, 0.46746, 0.50750, 0.54614, 0.58068, 0.61187, 0.64119, 0.66965, 0.69796, 0.72611, 0.75367, 0.77900, 0.80048, 0.81742, 0.82988, 0.83974, 0.84881, 0.85666, 0.86235, 0.86492, 0.87151, 0.87627, 0.88089, 0.88535, 0.88966, 0.89384, 0.89787, 0.90177, 0.90553, 0.90917, 0.91267, 0.91606, 0.91932, 0.92247, 0.92551, 0.92844, 0.92844, 0.92844, 0.92844, 0.92844, 0.92844, 0.92844, 0.92844, 0.92844, 0.92844, 0.92844 } },
	{ 95, 360.0, 830.0, 1.0, { 0.33851, 0.33851, 0.33851, 0.33851, 0.33851, 0.36145, 0.38504, 0.40918, 0.43394, 0.45859, 0.48408, 0.50808, 0.52830, 0.54340, 0.55592, 0.56810, 0.57702, 0.57960, 0.57732, 0.57245, 0.56571, 0.55731, 0.54695, 0.53443, 0.52046, 0.50593, 0.49141, 0.47729, 0.46359, 0.45025, 0.43720, 0.42459, 0.41355, 0.40517, 0.39949, 0.39625, 0.39511, 0.39575, 0.39816, 0.40244, 0.40897, 0.41800, 0.42915, 0.44206, 0.45737, 0.47574, 0.49692, 0.52053, 0.54655, 0.57500, 0.60596, 0.63909, 0.67241, 0.70387, 0.73294, 0.75943, 0.78318, 0.80405, 0.82201, 0.83704, 0.84892, 0.85751, 0.86339, 0.86749, 0.87174, 0.87761, 0.88387, 0.88860, 0.88988, 0.89554, 0.89928, 0.90290, 0.90640, 0.90979, 0.91307, 0.91624, 0.91930, 0.92226, 0.92513, 0.92789, 0.93056, 0.93314, 0.93563, 0.93803, 0.94035, 0.94035, 0.94035, 0.94035, 0.94035, 0.94035, 0.94035, 0.94035, 0.94035, 0.94035, 0.94035 } },
	{ 95, 360.0, 830.0, 1.0, { 0.11931, 0.11931, 0.11931, 0.11931, 0.11931, 0.13183, 0.14544, 0.16019, 0.17660, 0.19451, 0.21470, 0.22438, 0.22750, 0.23063, 0.23420, 0.23707, 0.23810, 0.23625, 0.23090, 0.22182, 0.20980, 0.19595, 0.18170, 0.16814, 0.15480, 0.14093, 0.12640, 0.11163, 0.09880, 0.08993, 0.08480, 0.08234, 0.08030, 0.07659, 0.07100, 0.06414, 0.05820, 0.05623, 0.06330, 0.08396, 0.11870, 0.16498, 0.21210, 0.24927, 0.27350, 0.28466, 0.28610, 0.28181, 0.27470, 0.26721, 0.26090, 0.25685, 0.25500, 0.25496, 0.25620, 0.25835, 0.26200, 0.26778, 0.27570, 0.28525, 0.29450, 0.30142, 0.30490, 0.30450, 0.30140, 0.29673, 0.28970, 0.28021, 0.27280, 0.26568, 0.25889, 0.25222, 0.24566, 0.23922, 0.23290, 0.22669, 0.22060, 0.21463, 0.20878, 0.20305, 0.19743, 0.19193, 0.18655, 0.18129, 0.17614, 0.17614, 0.17614, 0.17614, 0.17614, 0.17614, 0.17614, 0.17614, 0.17614, 0.17614, 0.17614 } },
	{ 95, 360.0, 830.0, 1.0, { 0.27425, 0.27425, 0.27425, 0.27425, 0.27425, 0.28297, 0.29185, 0.30090, 0.31040, 0.32028, 0.33277, 0.34184, 0.34836, 0.35424, 0.35855, 0.35975, 0.35650, 0.34789, 0.33440, 0.31699, 0.29728, 0.27687, 0.25673, 0.23768, 0.22049, 0.20561, 0.19230, 0.17970, 0.16765, 0.15642, 0.14727, 0.14136, 0.13849, 0.13788, 0.13781, 0.13679, 0.13518, 0.13405, 0.13543, 0.14163, 0.15519, 0.17835, 0.21210, 0.25629, 0.30765, 0.36205, 0.41520, 0.46276, 0.50047, 0.52541, 0.54006, 0.54796, 0.55154, 0.55303, 0.55484, 0.55912, 0.56677, 0.57838, 0.59454, 0.61548, 0.63996, 0.66655, 0.69439, 0.72255, 0.74925, 0.77306, 0.79489, 0.81494, 0.82814, 0.84022, 0.85125, 0.86163, 0.87141, 0.88058, 0.88919, 0.89725, 0.90478, 0.91182, 0.91838, 0.92450, 0.93019, 0.93548, 0.94040, 0.94496, 0.94920, 0.94920, 0.94920, 0.94920, 0.94920, 0.94920, 0.94920, 0.94920, 0.94920, 0.94920, 0.94920 } },
	{ 95, 360.0, 830.0, 1.0, { 0.16958, 0.16958, 0.16958, 0.16958, 0.16958, 0.17882, 0.18846, 0.19848, 0.20920, 0.22044, 0.23190, 0.23457, 0.23160, 0.22794, 0.22340, 0.21695, 0.20920, 0.20100, 0.19260, 0.18414, 0.17580, 0.16773, 0.15980, 0.15183, 0.14380, 0.13581, 0.12840, 0.12204, 0.11650, 0.11137, 0.10610, 0.10032, 0.09440, 0.08894, 0.08470, 0.08227, 0.08150, 0.08193, 0.08270, 0.08309, 0.08340, 0.08441, 0.08770, 0.09508, 0.10830, 0.12914, 0.15950, 0.20020, 0.24750, 0.29719, 0.34770, 0.39782, 0.44510, 0.48712, 0.52280, 0.55159, 0.57370, 0.58975, 0.60130, 0.60992, 0.61620, 0.62057, 0.62380, 0.62656, 0.62880, 0.63041, 0.63170, 0.63294, 0.63380, 0.63463, 0.63542, 0.63622, 0.63702, 0.63782, 0.63861, 0.63941, 0.64020, 0.64099, 0.64179, 0.64258, 0.64337, 0.64416, 0.64495, 0.64574, 0.64652, 0.64652, 0.64652, 0.64652, 0.64652, 0.64652, 0.64652, 0.64652, 0.64652, 0.64652, 0.64652 } }
};

#endif /* !STANDALONE_TEST */








































