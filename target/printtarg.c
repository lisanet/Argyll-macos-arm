
/* 
 * Argyll Color Correction System
 * PostScript print chart generator module.
 *
 * Author: Graeme W. Gill
 * Date:   28/9/96
 *
 * Copyright 1996 - 2009 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* 

	TTBD:

	Add option to return the number of patches that will
	exactly fit the given number of pages.

	Add optional device link processing support (same slot
	as -K file.cal) to permit a smoother proofing verification workflow. 

	Add independent w & h patch scaling option.

	Allow scaling minimum leading/trailing white space.

	Add option to omit labelling.

	Add "single pixel patch" mode, for pure digital processing for
	abstract profile creation.

	Add -h2 flag for Munki for super high-res chart ?
	Note:   i1Pro:  Illum spot: 3.5mm Aperture: 4.5mm, Physical aperture: 4.55mm

            Munki:  Illum spot: 8.0mm Aperture: 6.0mm, Physical aperture: 7.63mm
			Min patch size is 6mm x 6mm, below that increases delta E.

	Add an option that allows including a scale gauge, to detect
	accidental re-scaling.

	Make it aportion extra space evenly around the chart
	rather than at the trailing edges.

	Add direct PDF support, including NChannel output.

	Add option to apply a scale to counteract the Adobe utility problem.
*/

/* This program generates a PostScript or TIFF print target file, */
/* containing color test patches, given the .ti1 file specifying */
/* what the colors are. */

/* The output is designed to suite a general XY spectrometer (such as */
/* the Gretag SpectrScan), a handheld, manual instrument, or */
/* an Xrite DTP51, DTP41 or Eye-One strip spectrometer. */

/* Description:

	This program simply generates a PostScripto or TIFF file containing
	the patches layed out for an Xrite DTP20/DTP22/DTP51/DTP41/SpectroScan/i1pro/Munki.
	It allows them to be layed out on a choice of paper sizes,
	with the appropriate contrasting color spacers between
	each patch for the strip reading instruments. Unlike other
	charts, Argyll charts are generated as required, rather
	that being fixed. Also unlike most other strip reading charts,
	the spacers may colored, so that the density contrast ratio is
	guaranteed, even when two patches are about 50% density.

	Another feature is the pseudo random patch layout. This has
	three purposes. One is to try and average out any variation
	in the device response in relationship to the location of
	the patch on the paper. Color copiers and printing presses
	(for instance), are notorious in having side to side density
	variations.

	Another purpose of the random patch layout, is that it gives
	the reading program a good mechanism for detecting user error.
	It can guess the expected values, compare them to the readings,
	and complain if it seems that the strip is probably the wrong
	one. It can also be used to identify and rectify a strip
    that has been read in backwards.

	The final purpose of the random patch layout is to optimse the
	contrast between patches in a strip, to improve the robustness
	of the strip reading, and to be able to distinguish the directin
	a strip has been read in. Using this, small charts may be even be
	generated without any gaps between the test patches.

 */

/*
 * Nomencalture:
 *
 *	Largely due to how the strip readers name things, the following terms
 *  are used for how patches are grouped:
 *
 *  Pass, Row: One row of patches in a strip. A pass is usually labeled
 *             with a unique alphabetic label. 
 *  Strip:     A group of passes that can be read by a strip reader.
 *             For an XY instrument, the strip is a complete sheet, and
 *             each pass is one column. The rows of an XY chart are
 *             the step numbers within a pass.
 *  Step:      One test patch in a pass.
 *  Sheet:     One sheet of paper, containing full and partial strips.
 *             For an XY instrument, there will be only one strip per sheet.
 *           
 */

/* TTBD:
 *
 *    Improve EPS support to add a preview to each eps file.
 */

#undef DEBUG			/* Print edge details to stderr */
#undef FORCEN			/* For testing, force DeviceN */
#define DEN_COMPRESS	/* Compress density estimates > 1.0 */
						/* - this biases it towards white spacers */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <ctype.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "cgats.h"
#include "icc.h"
#include "xicc.h"
#include "insttypes.h"
#include "render.h"
#include "randix.h"
#include "alphix.h"
#include "rspl.h"
#include "sort.h"
#include "ui.h"

#include <stdarg.h>

/* Convert inches into mm */
#define inch2mm(xx) ((xx) * 25.4)

/* Convert mm into points */
#define mm2pnt(xx) ((xx) * 72.0/25.4)

/* A color structure */
struct _col {
	int nmask;		/* colorant mask */
	int altrep;		/* alternate grey or CMY representation type 0..8 */
	int i;			/* cols list index */
	int ix;			/* random list index */
	char *id;		/* Id string */
	char loc[10];	/* Location ID string */
	int t;			/* Tag */
#define T_XYZ    0x0001
#define T_LAB    0x0002
#define T_DEN    0x0004
#define T_RGB    0x0008
#define T_N      0x0010
#define T_NFB    0x2000			/* DeviceN fallback enabled */
#define T_PRESET 0x4000			/* A preset color rather than a test patch */
#define T_PAD    0x8000			/* A padding color patch */
	double XYZ[3];				/* Aproximate XYZ */
	double Lab[3];				/* Aproximate Lab */
	double den[4];				/* Approx statusT density + visual density */
	int    dtp20_octval;		/* DTP20 octal value */
	double dtp20_psize;			/* DTP20 patch width */
	double rgb[3];				/* Aproximate sRGB */
	int n;						/* Number of colorants */
	double dev[ICX_MXINKS];		/* Value of colorants */
	
	struct _col *nc[2];			/* Neigborhood colors */
	struct _col *oc;			/* Opposite direction color */
	double       wnd;			/* Worst neigborhood contrast density */
}; typedef struct _col col;

#define min2(a,b) ((a) < (b) ? (a) : (b))
#define min3(a,b,c)  (min2((a), min2((b),(c))))
#define max2(a,b) ((a) > (b) ? (a) : (b))
#define max3(a,b,c)  (max2((a), max2((b),(c))))

/* Declare edge tracking functions */
void et_init(void);
void et_height(double height);
void et_media(double *rgb);
void et_color(double *rgb);
void et_edge(int isx, int negh, double mj, double mi0, double mi1);
void et_patch(char *id, double xo, double yo, double w, double h);
void et_fiducial(double x, double y);
void et_write(char *fname, col *cols, int *rix, int si, int ei);
void et_clear(void);

/* ====================================================== */
/* Calibration Target rendering class */
/* Outputs either PostScript or TIFF raster */
/* test charts. */
/* We just do an error() if something goes wrong */

/* Common class structure */
#define TREND_STRUCT																	\
	/* Start a page */ 																	\
	void (*startpage)(struct _trend *s, int pn);										\
	/* End a page */ 																	\
	void (*endpage)(struct _trend *s);													\
	/* set the color */ 																\
	void (*setcolor)(struct _trend *s, xcal *cal, col *c);								\
	/* A rectangle, with optional edge tracking */										\
	void (*rectangle)(struct _trend *s,					/* Render a rectangle */		\
		double x, double y,		/* Top left corner of rectangle in mm from origin */	\
		double w, double h,		/* Width and height */									\
		char *id,				/* Patch id, NULL if not a diagnostic mark */			\
		int et					/* nz if use edge tracking on this */					\
	);																					\
	/* A testpad hexagon. */															\
	void (*hexagon)(struct _trend *s,													\
		double x, double y,		/* Top left vertex of hex mm from origin */				\
		double w, double h,		/* Width and height */									\
		int step,               /* Step number from 0 to figure odd/even */				\
		char *id				/* Patch id, NULL if not a diagnostic mark */			\
	);																					\
	/* A centered string */																\
	void (*string)(struct _trend *s,													\
		double x, double y,		/* Bot Left Corner of rectangle in mm from origin */	\
		double w, double h,		/* Width and height */									\
		char *str				/* String */											\
	);																					\
	/* A vertically centered string, rendered from bottom to top */						\
	void (*vstring)(struct _trend *s,													\
		double x, double y,		/* Bot Right Corner of rectangle in mm from origin */	\
		double w, double h,		/* Width and height */									\
		char *str				/* String */											\
	);																					\
	/* A dotted line */																	\
	void (*dline)(struct _trend *s,														\
		double x0, double y0,		/* Start of line */									\
		double x1, double y1,		/* End of line */									\
		double w					/* Width */											\
	);																					\
	/* Delete the object */ 															\
	void (*del)(struct _trend *s);														\

struct _trend {
	TREND_STRUCT
}; typedef struct _trend trend;

/* ==================================== */
/* PostScript output class */

struct _ps_trend {
	TREND_STRUCT
	FILE *of;			/* Postscript output file */
	int eps;			/* EPS flag */
	char *fname;
}; typedef struct _ps_trend ps_trend;

/* Start a page */
static void ps_startpage(trend *ss, int pagen) {
	ps_trend *s = (ps_trend *)ss; 

	fprintf(s->of,"%%%%Page: (Page %d) %d\n",pagen,pagen);
}

/* End a page */
static void ps_endpage(trend *ss) {
	ps_trend *s = (ps_trend *)ss; 

	fprintf(s->of,"showpage\n");
	fprintf(s->of,"\n");
}

/* Set a device N color with fallback */
static void
gen_ncolor(ps_trend *s, col *c) {
	int i;

	/* define the colorspace */
	fprintf(s->of,"[ /DeviceN [ ");
	for (i = 0; i < c->n; i++) {
		int imask = icx_index2ink(c->nmask, i);
		fprintf(s->of,"/%s ", icx_ink2psstring(imask));
	}

	if (c->t & T_NFB) {		/* Use color fallback */
		fprintf(s->of,"] /DeviceRGB ");	/* Fallback to RGB */
		fprintf(s->of,"{ ");
		for (i = 0; i < c->n; i++)		/* Remove N values */
			fprintf(s->of,"pop ");
		for (i = 0; i < 3; i++)			/* Set RGB values */
			fprintf(s->of,"%f ",c->rgb[i]);
	} else {
		fprintf(s->of,"] /DeviceGray ");	/* Fallback to Gray */
		fprintf(s->of,"{ ");
		for (i = 0; i < c->n; i++)		/* Remove N values */
			fprintf(s->of,"pop ");
		fprintf(s->of,"%f ",(c->rgb[0] + c->rgb[1] + c->rgb[2])/3.0); /* Set Gray value */
	}

	fprintf(s->of," } ] setcolorspace\n");

	/* Set the color */
	for (i = 0; i < c->n; i++)
		fprintf(s->of,"%f ",c->dev[i]);
	fprintf(s->of,"setcolor\n");
}


/* Set a device color */
/* Set it by the rep with most components */
static	void ps_setcolor(trend *ss, xcal *cal, col *c) {
	ps_trend *s = (ps_trend *)ss; 
	double cdev[ICX_MXINKS];     /* Calibrated device color */

	if (cal != NULL)
		cal->interp(cal, cdev, c->dev);	
	else {
		int j;
		for (j = 0; j < c->n; j++)
			cdev[j] = c->dev[j];
	}

	if ((c->t & T_N) == 0)
		error("ps_setcolor with no device values set");

#ifndef FORCEN
	if (c->nmask == ICX_W) {
		if ((c->t & T_PRESET) == 0)
			fprintf(s->of,"%% Ref %s %s %f\n",c->id, c->loc, 100.0 * cdev[0]);

		if (c->altrep == 0) {	/* DeviceGray */
			fprintf(s->of,"%f setgray\n",cdev[0]);
		} else if (c->altrep == 4) {	/* DeviceRGB */
			fprintf(s->of,"%f %f %f setrgbcolor\n",cdev[0],cdev[0],cdev[0]);
		} else if (c->altrep == 5) {	/* Separation */
			fprintf(s->of,"[ /Separation (White) /DeviceGray { pop %f } ] setcolorspace\n",
			           cdev[0]);
			fprintf(s->of,"%f setcolor\n",cdev[0]);
		} else if (c->altrep == 6) {	/* DeviceN */
			gen_ncolor(s, c);
		} else {
			error("Device white encoding not approproate!");
		}

	} else if (c->nmask == ICX_K) {
		if ((c->t & T_PRESET) == 0)
			fprintf(s->of,"%% Ref %s %s %f\n",c->id, c->loc, 100.0 * cdev[0]);
		if (c->altrep == 0) {	/* DeviceGray */
			fprintf(s->of,"%f setgray\n",1.0 - cdev[0]);
		} else if (c->altrep == 1) {	/* DeviceCMYK */
			fprintf(s->of,"0.0 0.0 0.0 %f setcmykcolor\n",cdev[0]);
		} else if (c->altrep == 2) {	/* Separation */
			fprintf(s->of,"[ /Separation (Black) /DeviceGray { pop %f } ] setcolorspace\n",
			           1.0 - cdev[0]);
			fprintf(s->of,"%f setcolor\n",cdev[0]);
		} else if (c->altrep == 3) {	/* DeviceN */
			gen_ncolor(s, c);
		} else {
			error("Device black encoding not approproate!");
		}

	} else if (c->nmask == ICX_CMY) {
		if ((c->t & T_PRESET) == 0)
			fprintf(s->of,"%% Ref %s %s %f %f %f\n", c->id, c->loc,
			        100.0 * cdev[0], 100.0 * cdev[1], 100.0 * cdev[2]);

		if (c->altrep == 0) {			/* DeviceCMYK */
			fprintf(s->of,"%f %f %f 0.0 setcmykcolor\n",cdev[0],cdev[1],cdev[2]);
		} else if (c->altrep == 7) {	/* Inverted DeviceRGB */
			fprintf(s->of,"%f %f %f setrgbcolor\n",1.0-cdev[0],1.0-cdev[1],1.0-cdev[2]);
		} else if (c->altrep == 8) {	/* DeviceN */
			gen_ncolor(s, c);
		} else {
			error("Device CMY encoding not approproate!");
		}

	} else if (c->nmask == ICX_RGB || c->nmask == ICX_IRGB) {
		if ((c->t & T_PRESET) == 0)
			fprintf(s->of,"%% Ref %s %s %f %f %f\n",c->id, c->loc,
			       100.0 * cdev[0], 100.0 *cdev[1], 100.0 *cdev[2]);
		fprintf(s->of,"%f %f %f setrgbcolor\n",cdev[0],cdev[1],cdev[2]);

	} else if (c->nmask == ICX_CMYK) {
		if ((c->t & T_PRESET) == 0)
			fprintf(s->of,"%% Ref %s %s %f %f %f %f\n", c->id, c->loc,
			        100.0 * cdev[0], 100.0 * cdev[1], 100.0 * cdev[2], 100.0 * cdev[3]);
		fprintf(s->of,"%f %f %f %f setcmykcolor\n",cdev[0],cdev[1],cdev[2],cdev[3]);

	} else
#endif /* !FORCEN */
	       {	/* Device N */
		int i;
		if ((c->t & T_PRESET) == 0) {
			fprintf(s->of,"%% Ref %s %s",c->id, c->loc);
			for (i = 0; i < c->n; i++)
				fprintf(s->of,"%f ", 100.0 * cdev[i]);
			fprintf(s->of,"\n");
		}
		gen_ncolor(s, c);
	}

	/* Remember edge tracking color */
	et_color(c->rgb);
}

/* Generate a rectangle, with optional edge tracking */
/* Note the page coordinate origin is bottom left. */
static void ps_rectangle(trend *ss,
	double x, double y,		/* Top left corner of rectangle in mm from origin */
	double w, double h,		/* Width and height */
	char *id,				/* Patch id, NULL if not a diagnostic mark */
	int et					/* nz if use edge tracking on this */
) {
	ps_trend *s = (ps_trend *)ss; 

	if (w < 1e-6 || h < 1e-6)
		return;			/* Skip zero sized rectangle */
	y -= h;				/* Convert to bottom left corner */
	x = mm2pnt(x);
	y = mm2pnt(y);
	w = mm2pnt(w);
	h = mm2pnt(h);
	fprintf(s->of,"%f %f %f %f rect\n",w,h,x,y);

	if (et) {
		et_patch(id, x, y, w, h);
		et_edge(1, 0, x, y, y + h);
		et_edge(1, 1, x + w, y, y + h);
		et_edge(0, 0, y, x, x + w); 
		et_edge(0, 1, y + h, x, x + w);
	}
}

/* Generate one testpad hexagon. */
/* Note the page coordinate origin is bottom left. */
/* The hex always has left/right sides */
/* and peaks at the top and the bottom. */
static void ps_hexagon(trend *ss,
	double x, double y,		/* Top left vertex of hex mm from origin */
	double w, double h,		/* Width and height */
	int step,               /* Step number from 0 to figure odd/even */
	char *id				/* Patch id, NULL if not a diagnostic mark */
) {
	ps_trend *s = (ps_trend *)ss; 

	if (w < 1e-6 || h < 1e-6)
		return;			/* Skip zero sized rectangle */
	if ((step & 1) == 0) /* Even so left side of stagger */
		x -= 0.25 * w;
	else 				 /* Odd so right side of stagger */
		x += 0.25 * w;
	y = y - 5.0/6.0 * h;
	h *= 2.0/3.0;		/* Convert to hex side length */
	x = mm2pnt(x);
	y = mm2pnt(y);
	w = mm2pnt(w);
	h = mm2pnt(h);
	fprintf(s->of,"%f %f %f %f hex\n",w,h,x,y);
}

/* A centered string */
static void ps_string(trend *ss,
	double x, double y,		/* Bot Left Corner of rectangle in mm from origin */
	double w, double h,		/* Width and height */
	char *str				/* String */
) {
	ps_trend *s = (ps_trend *)ss; 

	if (fabs(w) < 1e-6 || fabs(h) < 1e-6)
		return;			/* Skip zero sized string */
	x = mm2pnt(x);
	y = mm2pnt(y);
	w = mm2pnt(w);
	h = mm2pnt(h);
	fprintf(s->of,"%f scaleTimes\n",h * 0.75);
	fprintf(s->of,"(%s) %f %f centerShow\n",str,x+w/2.0,y+h/2.0);
}

/* A vertically centered string */
static void ps_vstring(trend *ss,
	double x, double y,		/* Bot Right Corner of rectangle in mm from origin */
	double w, double h,		/* Width and height */
	char *str				/* String */
) {
	ps_trend *s = (ps_trend *)ss; 

	if (fabs(w) < 1e-6 || fabs(h) < 1e-6)
		return;			/* Skip zero sized string */
	x = mm2pnt(x);
	y = mm2pnt(y);
	w = mm2pnt(w);
	h = mm2pnt(h);
	fprintf(s->of,"%f scaleTimes\n",w * 0.75);
	fprintf(s->of,"(%s) %f %f vcenterShow\n",str,x-w/2.0,y+h/2.0);
}

/* A dotted line */
static void ps_dline(trend *ss,
	double x1, double y1,		/* Start of line */
	double x2, double y2,		/* End of line */
	double w					/* Width */
) {
	ps_trend *s = (ps_trend *)ss; 

	if (fabs(w) < 1e-6 || fabs(y2 - y1) < 1e-6)
		return;			/* Skip zero sized line */
	x1 = mm2pnt(x1);
	x2 = mm2pnt(x2);
	y1 = mm2pnt(y1);
	y2 = mm2pnt(y2);
	w = mm2pnt(w);
	fprintf(s->of,"[%f %f] %f setdash\n",mm2pnt(1.0),mm2pnt(2.0),mm2pnt(0.0));
	fprintf(s->of,"%f setlinewidth\n",w);
	fprintf(s->of,"newpath %f %f moveto %f %f lineto stroke\n",x1,y1,x2,y2);
}

/* Complete operations, then delete the object */ 
static void ps_del(trend *ss) {
	ps_trend *s = (ps_trend *)ss; 

	if (s->of) {
		fprintf(s->of,"\n");
		fprintf(s->of,"%%%%EOF\n");

		if (fclose(s->of))
			error ("Unable to close output file '%s'",s->fname);
	}
	if (s->fname)
		free(s->fname);

	free(s);
}

trend *new_ps_trend(
	char *fname,			/* File name */
	int npages,				/* Number of pages needed */
	int nmask,				/* Non zero if we are doing a DeviceN chart */
	double pw, double ph,	/* Page width and height in mm */
	int eps,				/* EPS flag */
	int nocups,				/* NZ to supress cups job ticket */
	int rand,				/* randomize */
	int rstart				/* Random start number/chart ID */
) {
	ps_trend *s;

	if ((s = (ps_trend *)calloc(1, sizeof(ps_trend))) == NULL) {
		return NULL;
	}

	s->startpage = ps_startpage;
	s->endpage = ps_endpage;
	s->setcolor = ps_setcolor;
	s->rectangle = ps_rectangle;
	s->hexagon = ps_hexagon;
	s->string = ps_string;
	s->vstring = ps_vstring;
	s->dline = ps_dline;
	s->del = ps_del;

	s->eps = eps;

	if ((s->of = fopen(fname,"w")) == NULL) {
		error ("Unable to open output file '%s'",fname);
	}

	if ((s->fname = strdup(fname)) == NULL)
		error("stdup of fname failed");
	
	/* Generate PS file prolog */
	{
		time_t clk = time(0);
		struct tm *tsp = localtime(&clk);
		int ipw, iph;

		ipw = (int)ceil(mm2pnt(pw));
		iph = (int)ceil(mm2pnt(ph));
		if (eps)
			fprintf(s->of,"%%!PS-Adobe-3.0 EPSF-3.0\n");
		else
			fprintf(s->of,"%%!PS-Adobe-3.0\n");
		fprintf(s->of,"%%%%Title: Argyll Color Calibration Target\n");
#ifdef FORCEN
		if (1) {
#else
		if (nmask != ICX_W		/* If not a Gray, RGB or CMYK device space */
		 && nmask != ICX_K
		 && nmask != ICX_RGB
		 && nmask != ICX_IRGB
		 && nmask != ICX_CMYK) {
#endif
			fprintf(s->of,"%%%%LanguageLevel: 3\n");
		} else {
			fprintf(s->of,"%%%%LanguageLevel: 1\n");
			if (nmask == ICX_CMYK)
				fprintf(s->of,"%%%%Extensions: CMYK\n");
		}

		fprintf(s->of,"%%%%Creator: Argyll target chart generator\n");
		fprintf(s->of,"%%%%For: The user who wants accurate color\n");
//	fprintf(s->of,"%%%%Version: %s\n",VERSION);
		fprintf(s->of,"%%%%CreationDate: %s",asctime(tsp));
		fprintf(s->of,"%%%%DocumentData: Clean7Bit\n");
		if (eps)
			fprintf(s->of,"%%%%Pages: %d\n",1);
		else
			fprintf(s->of,"%%%%Pages: %d\n",npages);
		fprintf(s->of,"%%%%PageOrder: Ascend\n");
		fprintf(s->of,"%%%%BoundingBox: %d %d %d %d\n",0,0,ipw-1,iph-1);
		fprintf(s->of,"%%%%Orientation: Portrait\n");		/* Rows are always virtical */
		if (!nocups)
			fprintf(s->of,"%%cupsJobTicket: cups-disable-cmm\n");
		fprintf(s->of,"%%%%EndComments\n");
		fprintf(s->of,"\n");
		if (!eps) {
			fprintf(s->of,"<< /PageSize [%d %d] >> setpagedevice\n",ipw, iph);
			fprintf(s->of,"\n");
		}
		fprintf(s->of,"%%%%BeginProlog\n\n");
#ifdef NEVER
		fprintf(s->of,"%% Duplicate nth element of stack\n");
		fprintf(s->of,"%% arguments: n, the offset from the element bellow the n\n");
		fprintf(s->of,"/dupn { 2 add dup -1 roll dup 3 -1 roll 1 roll } bind def\n");
		fprintf(s->of,"\n");
#endif
		fprintf(s->of,"%% arbitrary rectangle\n");
		fprintf(s->of,"%% arguments: w h x y\n");
		fprintf(s->of,"/rect { gsave \n");
		fprintf(s->of,"newpath\n");
		fprintf(s->of,"moveto\n");
		fprintf(s->of,"dup 0.0 exch rlineto\n");
		fprintf(s->of,"exch 0.0  rlineto\n");
		fprintf(s->of,"0.0 exch neg rlineto\n");
		fprintf(s->of,"closepath\n");
		fprintf(s->of,"fill\n");
		fprintf(s->of,"grestore } bind def\n");
		fprintf(s->of,"\n");
		fprintf(s->of,"%% hexagon with bottom left origin\n");
		fprintf(s->of,"%% arguments: w h x y\n");
		fprintf(s->of,"/hex { gsave \n");
		fprintf(s->of,"newpath\n");
		fprintf(s->of,"moveto\n");
		fprintf(s->of,"0 1 index rlineto\n");
		fprintf(s->of,"1 index 2 div 1 index 2 div rlineto\n");
		fprintf(s->of,"1 index 2 div 1 index 2 div neg rlineto\n");
		fprintf(s->of,"0 1 index neg rlineto\n");
		fprintf(s->of,"1 index 2 div neg 1 index 2 div neg rlineto\n");
		fprintf(s->of,"pop pop\n");
		fprintf(s->of,"closepath\n");
		fprintf(s->of,"fill\n");
		fprintf(s->of,"grestore } bind def\n");
		fprintf(s->of,"\n");
		fprintf(s->of,"%% set times-roman font\n");
		fprintf(s->of,"%% argument: scale\n");
		fprintf(s->of,"/scaleTimes {\n");
		fprintf(s->of,"/Times-Roman findfont\n");
		fprintf(s->of,"exch scalefont\n");
		fprintf(s->of,"setfont } bind def\n");
		fprintf(s->of,"\n");
		fprintf(s->of,"%% Print a centered string\n");
		fprintf(s->of,"%% argument: string, x, y\n");
		fprintf(s->of,"/centerShow {\n");
		fprintf(s->of,"gsave translate\n");
		fprintf(s->of,"newpath 0.0 0.0 moveto dup true charpath pathbbox\n");
		fprintf(s->of,"3 -1 roll sub exch 3 -1 roll sub\n");
		fprintf(s->of,"-0.5 mul exch -0.5 mul\n");
		fprintf(s->of,"moveto show grestore} bind def\n");
		fprintf(s->of,"\n");
		fprintf(s->of,"%% Print a vertically centered string\n");
		fprintf(s->of,"%% argument: string, x, y\n");
		fprintf(s->of,"/vcenterShow {\n");
		fprintf(s->of,"gsave translate 90.0 rotate\n");
		fprintf(s->of,"newpath 0.0 0.0 moveto dup true charpath pathbbox\n");
		fprintf(s->of,"3 -1 roll sub exch 3 -1 roll sub\n");
		fprintf(s->of,"-0.5 mul exch -0.5 mul\n");
		fprintf(s->of,"moveto show grestore} bind def\n");

		fprintf(s->of,"%%%%EndProlog\n");
		fprintf(s->of,"\n");

		if (rand != 0)
			fprintf(s->of,"%% RandomStart %d\n",rstart);
		else
			fprintf(s->of,"%% ChartID %d\n",rstart);
		fprintf(s->of,"\n");
	}

	return (trend *)s;
}

/* ==================================== */
/* TIFF raster output class */
/* We use the render library to do all the hard work. */

struct _tiff_trend {
	TREND_STRUCT
	
	render2d *r;		/* Raster renderer object */
	char *fname;
	color2d c;			/* Last set color */
	int comp;			/* Flag, use compression */

}; typedef struct _tiff_trend tiff_trend;

/* Start a page */
static void tiff_startpage(trend *ss, int pn) {
	tiff_trend *s = (tiff_trend *)ss; 

	/* Nothing to do */
}

/* End a page */
static void tiff_endpage(trend *ss) {
	tiff_trend *s = (tiff_trend *)ss; 

	/* Nothing to do */
}

/* set the color */
static	void tiff_setcolor(trend *ss, xcal *cal, col *c) {
	tiff_trend *s = (tiff_trend *)ss; 
	double cdev[ICX_MXINKS];     /* Calibrated device color */

	if (cal != NULL)
		cal->interp(cal, cdev, c->dev);	
	else {
		int j;
		for (j = 0; j < c->n; j++)
			cdev[j] = c->dev[j];
	}

	if ((c->t & T_N) == 0)
		error("tiff_setcolor with no device values set");

	if (c->nmask == ICX_W) {
		if (c->altrep == 0) {	/* DeviceGray */
			s->c[0] = cdev[0];
		} else if (c->altrep == 4) {	/* DeviceRGB */
			s->c[0] = cdev[0];
			s->c[1] = cdev[0];
			s->c[2] = cdev[0];
		} else if (c->altrep == 5) {	/* Separation */
			s->c[0] = cdev[0];
		} else if (c->altrep == 6) {	/* DeviceN single channel */
			s->c[0] = cdev[0];
		} else {
			error("Device white encoding not approproate!");
		}

	} else if (c->nmask == ICX_K) {
		if (c->altrep == 0) {	/* DeviceGray */
			s->c[0] = cdev[0];
		} else if (c->altrep == 1) {	/* DeviceCMYK */
			s->c[0] = 0.0;
			s->c[1] = 0.0;
			s->c[2] = 0.0;
			s->c[3] = cdev[0];
		} else if (c->altrep == 2) {	/* Separation */
			s->c[0] = cdev[0];
		} else if (c->altrep == 3) {	/* DeviceN single channel */
			s->c[0] = cdev[0];
		} else {
			error("Device black encoding not approproate!");
		}

	} else if (c->nmask == ICX_CMY) {
		if (c->altrep == 0) {			/* DeviceCMYK */
			s->c[0] = cdev[0];
			s->c[1] = cdev[1];
			s->c[2] = cdev[2];
			s->c[3] = 0.0;
		} else if (c->altrep == 7) {	/* Inverted DeviceRGB */
			s->c[0] = 1.0-cdev[0];
			s->c[1] = 1.0-cdev[1];
			s->c[2] = 1.0-cdev[2];
		} else if (c->altrep == 8) {	/* DeviceN three channel */
			s->c[0] = cdev[0];
			s->c[1] = cdev[1];
			s->c[2] = cdev[2];
		} else {
			error("Device CMY encoding not approproate!");
		}

	} else {
		int j;
		for (j = 0; j < s->r->ncc; j++)
			s->c[j] = cdev[j];
	}

	/* Remember edge tracking color */
	et_color(c->rgb);
}

/* A rectangle, with optional edge tracking */
static void tiff_rectangle(trend *ss,
	double x, double y,		/* Top left corner of rectangle in mm from origin */
	double w, double h,		/* Width and height */
	char *id,				/* Patch id, NULL if not a diagnostic mark */
	int et					/* nz if use edge tracking on this */
) {
	tiff_trend *s = (tiff_trend *)ss; 

	y -= h;				/* Convert to bottom left corner */
	s->r->add(s->r, new_rect2d(s->r, x, y, w, h, s->c));

	if (et) {
		et_patch(id, x, y, w, h);
		et_edge(1, 0, x, y, y + h);
		et_edge(1, 1, x + w, y, y + h);
		et_edge(0, 0, y, x, x + w); 
		et_edge(0, 1, y + h, x, x + w);
	}
}

/* A testpad hexagon. */
static void tiff_hexagon(trend *ss,
	double x, double y,		/* Top left vertex of hex mm from origin */
	double w, double h,		/* Width and height */
	int step,               /* Step number from 0 to figure odd/even */
	char *id				/* Patch id, NULL if not a diagnostic mark */
) {
	tiff_trend *s = (tiff_trend *)ss; 
	double vv[3][2];
	color2d cc[3];
	int i, j;

	if ((step & 1) == 0) /* Even so left side of stagger */
		x -= 0.25 * w;
	else 				 /* Odd so right side of stagger */
		x += 0.25 * w;
	y = y - 5.0/6.0 * h;
	h *= 2.0/3.0;		/* Convert to hex side length */

	/* Triangle color */
	for (i = 0; i < 3; i++)
		for (j = 0; j < s->r->ncc; j++)
			cc[i][j] = s->c[j];

	/* Top triangle */
	vv[0][0] = x;
	vv[0][1] = y + h;
	vv[1][0] = x + w;
	vv[1][1] = y + h;
	vv[2][0] = x + 0.5 * w;
	vv[2][1] = y + 1.5 * h;
	s->r->add(s->r, new_trivs2d(s->r, vv, cc));
	
	/* Center rectangle */
	s->r->add(s->r, new_rect2d(s->r, x, y, w, h, s->c));

	/* Bottom triangle */
	vv[0][0] = x;
	vv[0][1] = y;
	vv[1][0] = x + w;
	vv[1][1] = y;
	vv[2][0] = x + 0.5 * w;
	vv[2][1] = y - 0.5 * h;
	s->r->add(s->r, new_trivs2d(s->r, vv, cc));
}

/* A centered string */
static void tiff_string(trend *ss,
	double x, double y,		/* Bot Left Corner of rectangle in mm from origin */
	double w, double h,		/* Width and height */
	char *str				/* String */
) {
	tiff_trend *s = (tiff_trend *)ss; 
	double sw = 0.0, sh = 0.0;

	sh = h * 0.57;
	meas_string2d(s->r,&sw,NULL,timesr_b,str,sh,0);
	/* Center the string within the recangle */
	x += 0.5 * (w - sw);
	y += 0.5 * (h - sh);
	add_string2d(s->r,NULL,NULL,timesr_b,str,x,y,sh,0,s->c);
}

/* A vertically centered string */
static void tiff_vstring(trend *ss,
	double x, double y,		/* Bot Right Corner of rectangle in mm from origin */
	double w, double h,		/* Width and height */
	char *str				/* String */
) {
	tiff_trend *s = (tiff_trend *)ss; 
	double sw = 0.0, sh = 0.0;

	x -= w;		/* Make it bot left corner */
	sw = w * 0.57;
	meas_string2d(s->r,NULL, &sh,timesr_b,str,sw,3);
	/* Center the string within the recangle */
	x += 0.5 * (w + sw);
	y += 0.5 * (h - sh);
	add_string2d(s->r,NULL,NULL,timesr_b,str,x,y,sw,3,s->c);
}

/* A dotted line */
static void tiff_dline(trend *ss,
	double x0, double y0,		/* Start of line */
	double x1, double y1,		/* End of line */
	double w					/* Width */
) {
	tiff_trend *s = (tiff_trend *)ss; 

	add_dashed_line2d(s->r,x0,y0,x1,y1,w,1.0,2.0,0,s->c);
}

/* Complete operations, then delete the object */ 
static void tiff_del(trend *ss) {
	tiff_trend *s = (tiff_trend *)ss; 

	if (s->r != NULL) {
		s->r->write(s->r, s->fname, s->comp, NULL, NULL, tiff_file);
		s->r->del(s->r);
	}
	if (s->fname != NULL)
		free(s->fname);
	free(s);
}

static trend *new_tiff_trend(
	char *fname,				/* File name */
	int nmask,					/* Non zero if we are doing a DeviceN chart */
	depth2d dpth,				/* 8 or 16 bit */
	double pw, double ph,		/* Page width and height in mm */
	double marg,				/* Page margine in mm */
	double hres, double vres,	/* Resolution */
	int altrep,					/* printer grey/CMY representation type 0..8 */
	int ncha,					/* flag, use nchannel alpha */
	int comp,					/* flag, use compression */
	int dith					/* flag, 1 = use 8 bit stocastic dithering */
) {
	tiff_trend *s;
	color2d c;					/* Background color */
	colort2d csp;
	double ma[4];				/* Page margins */
	int nc = 0;
	int j;

	if ((s = (tiff_trend *)calloc(1, sizeof(tiff_trend))) == NULL) {
		error("Failed to create a tiff target rendering object");
	}

	s->startpage = tiff_startpage;
	s->endpage = tiff_endpage;
	s->setcolor = tiff_setcolor;
	s->rectangle = tiff_rectangle;
	s->hexagon = tiff_hexagon;
	s->string = tiff_string;
	s->vstring = tiff_vstring;
	s->dline = tiff_dline;
	s->del = tiff_del;

	if (nmask == ICX_W) {
		if (altrep == 0				/* DeviceGray */
		 || altrep == 5) {			/* Separation single channel */
			csp = w_2d;
			nc = 1;
		} else if (altrep == 4) {	/* DeviceRGB */
			csp = rgb_2d;
			nc = 3;
		} else if (altrep == 6) {	/* DeviceN single channel */
			csp = ncol_2d;
			nc = icx_noofinks(nmask);
			nc = 1;
		} else {
			error("Device white encoding not approproate");
		}

	} else if (nmask == ICX_K) {
		if (altrep == 0				/* DeviceGray */
		 || altrep == 2) {			/* Separation single channel */
			csp = k_2d;
			nc = 1;
		} else if (altrep == 1) {	/* DeviceCMYK */
			csp = cmyk_2d;
			nc = 4;
		} else if (altrep == 3) {	/* DeviceN single channel */
			csp = ncol_2d;
			nc = icx_noofinks(nmask);
			nc = 1;
		} else {
			error("Device black encoding not approproate");
		}

	} else if (nmask == ICX_RGB || nmask == ICX_IRGB) {
		csp = rgb_2d;
		nc = 3;

	} else if (nmask == ICX_CMY) {
		if (altrep == 0) {			/* DeviceCMYK */
			csp = cmyk_2d;
			nc = 4;
		} else if (altrep == 7) {	/* Inverted DeviceRGB */
			csp = rgb_2d;
			nc = 3;
		} else if (altrep == 8) {	/* DeviceN three channel */
			csp = ncol_2d;
			nc = icx_noofinks(nmask);
		} else {
			error("Device CMY encoding not approproate");
		}

	} else if (nmask == ICX_CMYK) {
		csp = cmyk_2d;
		nc = 4;
	} else {	/* Device N */
		if (ncha)
			csp = ncol_a_2d;
		else
			csp = ncol_2d;
		nc = icx_noofinks(nmask);
	}

	if (marg > 0)
		ma[0] = ma[1] = ma[2] = ma[3] = marg;
	else
		ma[0] = ma[1] = ma[2] = ma[3] = 0;

	if ((s->r = new_render2d(pw, ph, ma, hres, vres,  csp, nc, dpth, dith, NULL, NULL, 0.0)) == NULL) {
		error("Failed to create a render2d object for tiff output");
	} 

	/* We're going to assume this is all printed output, so */
	/* the background should be white. */
	if ((nmask & ICX_ADDITIVE)
	  || (nmask == ICX_CMY && altrep == 7)) {	/* CMY as inverted RGB */
		for (j = 0; j < nc; j++)
			c[j] = 1.0;
	} else {
		for (j = 0; j < nc; j++)
			c[j] = 0.0;
	}
	s->r->set_defc(s->r, c);

	s->comp = comp;

	if ((s->fname = strdup(fname)) == NULL)
		error("stdup of fname failed");

	return (trend *)s;
}


/* ====================================================== */

/* Convert XYZ represention into Lab, XYZ density and RGB */
void
col_convert(col *cp, double *wp) {

	if ((cp->t & T_XYZ) == 0 || (cp->t & T_N) == 0)
		error("gen_color needs XYZ and device colors set !");

	if ((cp->t & T_LAB) == 0) {
		icmXYZNumber w;
		w.X = wp[0];
		w.Y = wp[1];
		w.Z = wp[2];
		icmXYZ2Lab(&w, cp->Lab, cp->XYZ);
		cp->t |= T_LAB;
	}
	if ((cp->t & T_DEN) == 0) {
		icx_XYZ2Tdens(cp->den, cp->XYZ);

#ifdef DEN_COMPRESS
		/* Compress densities > 1.0, 2:1 */
		{
			int e;
			for (e = 0; e < 3; e++) {
				if (cp->den[e] > 1.0)
					cp->den[e] = 1.0 + (cp->den[e] - 1.0) * 0.5;
			}
		}
#endif /* DEN_COMPRESS */
		cp->t |= T_DEN;
	}

	if ((cp->t & T_RGB) == 0) {
		icx_XYZ2sRGB(cp->rgb, wp, cp->XYZ);
		cp->t |= T_RGB;
	}
}

/* return the middle of 3 values */
double mid3(double a, double b, double c) {

	if ((a < b && a > c) || (a < c && a > b))
		return a;
	if ((b < a && b > c) || (b < c && b > a))
		return b;
	return c;
}

/* return the vector difference of 3 values */
double vec3(double a, double b, double c) {

	return sqrt(a * a + b * b + c * c);
}

/* Type of density rule to use with color patches/spacers. */
/* This really depends on the algorithm used by the strip */
/* reading instrument. For an Xrite DTP41, it appears that */
/* using max3 is the best. This agrees with their documentation, */
/* that it is looking for the largest change in one of the three */
/* channels. Another make of instrument might use a different */
/* algorithm. */
/* Choices are: */
/* max3 - aim for largest change to be greatest */
/* mid3 - aim for middle change to be greatest */
/* min3 - aim for minimum change to be greatest */
/* vec3 - aim for the vector change to be greatest */

#define RULE3 max3

/* Setup a suitable spacer color, and */
/* Return the worst case density difference. */
double setup_spacer(
col **psc,		/* Return pointer to spacer color */
col *pp,		/* Previous patch color */
col *cp,		/* Current patch color */
col *pcol,		/* 8 pre-defined spacer colors */
int sptype,		/* Spacer type code, -1 = none, */
				/* 0 = No spacer, 1 = b&W spacer, */
				/* 2 = colored */
int usede		/* Aim for maximum delta E rather than density */
) {
	col *sc = NULL;	/* Spacer chosen */
	double dd, pdd;

//printf("~1\n");
//printf("~1 setting spacer between %s (%s) and %s (%s)\n",pp->loc, pp->id, cp->loc, cp->id);

	/* Compute contrast between the patches */
	if (usede) {
		pdd = icmLabDE(pp->Lab, cp->Lab); 
	} else {
		/* return the density contrast between the patches */
		if (pp->nmask == ICX_W
		 || pp->nmask == ICX_K) {	/* If only capable of single density */

			pdd = fabs(pp->den[3] - cp->den[3]);
//printf("~1 computed B&W diff of %f\n",dd);

		} else {
			pdd = RULE3(fabs(pp->den[0] - cp->den[0]),
			            fabs(pp->den[1] - cp->den[1]),
			            fabs(pp->den[2] - cp->den[2]));

//printf("~1 computed color diff of %f\n",dd);
		}
	}

	if (sptype <= 0)	/* No spacers */
		return pdd;

	if (pp->nmask == ICX_W
	 || pp->nmask == ICX_K) {	/* If only capable of single density */
		sptype = 1;				/* Force to B&W spacer */
	}

	if (sptype == 1) {			/* B&W spacer */
		double d1, d2;

		if (usede) {
			/* Choose whether space should be white or black */
			/* Shoose color that has greatest worst contrast */
			d1 = min2(icmLabDE(pcol[0].Lab, pp->Lab), icmLabDE(pcol[0].Lab, cp->Lab));
			d2 = min2(icmLabDE(pcol[7].Lab, pp->Lab), icmLabDE(pcol[7].Lab, cp->Lab));

//printf("~1 worst difference to white = %f\n", d1);
//printf("~1 worst difference to black = %f\n", d2);

			if (d1 > d2) {
//printf("~1 chosen white\n");
				sc = &pcol[0];
				dd = d1;
			} else {
//printf("~1 chosen black\n");
				sc = &pcol[7];
				dd = d2;
			}

		} else {
			/* Choose whether space should be white or black */
			/* Shoose color that has greatest worst contrast */
			d1 = min2(fabs(pcol[0].den[3] - pp->den[3]), fabs(pcol[0].den[3] - cp->den[3]));
			d2 = min2(fabs(pcol[7].den[3] - pp->den[3]), fabs(pcol[7].den[3] - cp->den[3]));

//printf("~1 worst difference to white = %f\n", d1);
//printf("~1 worst difference to black = %f\n", d2);

			if (d1 > d2) {
//printf("~1 chosen white\n");
				sc = &pcol[0];
				dd = d1;
			} else {
//printf("~1 chosen black\n");
				sc = &pcol[7];
				dd = d2;
			}
		}

	} else {				/* else colored spacer */
		int ii, i;

		/* Check out all the possible space values for the one that gives the best */
		/* and second best contrast to each edge */

		if (usede) {

			/* for all possible spacer colors */
			dd = -1.0;
			for (i = 0; i < 8; i++) {
				double bb;

				bb = min2(icmLabDE(pcol[i].Lab, pp->Lab), icmLabDE(pcol[i].Lab, cp->Lab));

				/* Worst of two edges best is better than any previous */
				if (bb > dd) {
					dd = bb;		/* Worst color of worst edge */
					ii = i;
					sc = &pcol[i];
				}
			}

		} else {
			/* for all possible spacer colors */
			dd = -1.0;
			for (i = 0; i < 8; i++) {
				double b1, b2, bb;

				b1 = RULE3(fabs(pcol[i].den[0] - pp->den[0]),
				          fabs(pcol[i].den[1] - pp->den[1]),
				          fabs(pcol[i].den[2] - pp->den[2]));

				b2 = RULE3(fabs(pcol[i].den[0] - cp->den[0]),
				          fabs(pcol[i].den[1] - cp->den[1]),
				          fabs(pcol[i].den[2] - cp->den[2]));

				bb = min2(b1, b2);	/* Worst of two edges */

				/* Worst of two edges best is better than any previous */
				if (bb > dd) {
					dd = bb;		/* Worst color of worst edge */
					ii = i;
					sc = &pcol[i];
				}
			}
		}
	}

//printf("~1 returning spacer contrast %f + patch contrast %f\n",dd,pdd);
	if (psc != NULL)
		*psc = sc;			/* Return pointer to chosen spacer */

	return 0.6 * dd + 0.4 * pdd;	/* Return propotion of spacer and patch contrast */
}

/* Given two patches, compute the density difference between them */
double density_difference(
col *p1,		/* Previous patch color */
col *p2,		/* Current patch color */
int usede		/* Aim for maximum delta E rather than density */
) {
	double dd, pdd;

	/* Compute contrast between the patches */
	if (usede) {
		pdd = icmLabDE(p1->Lab, p2->Lab); 
	} else {
		/* return the density contrast between the patches */
		if (p1->nmask == ICX_W
		 || p1->nmask == ICX_K) {	/* If only capable of single density */

			pdd = fabs(p1->den[3] - p2->den[3]);

		} else {
			pdd = RULE3(fabs(p1->den[0] - p2->den[0]),
			            fabs(p1->den[1] - p2->den[1]),
			            fabs(p1->den[2] - p2->den[2]));

		}
	}

	return pdd;
}

/* Given a number, return the DTP20 SID patch encoding. */
/* return nz on error */
static int dtp20_enc(
	col *pcol,			/* pcol[8] of spacer patch colors */
	int ndig,			/* Number of octal digits/patches */
	int lend,			/* NZ if little endian order for SID row index */
	col **ppcol,		/* return the pointers to pcol for the encoding */
	unsigned int rix	/* Number to encode */
) {
	int si, ei, ii;
	int i, wv = rix;

	if (lend) 			/* Little endian order */
		si = 0, ei = ndig, ii = 1;
	else
		si = ndig-1, ei = -1, ii = -1;

	/* LSB to MSB octal */
	for (i = si; i != ei; i += ii) {
		int j, k;
		j = wv % 8;		/* Digit value */
		wv /= 8;
		for (k = 0; k < 8; k++) {
			if (pcol[k].dtp20_octval == j)
				break;
		}
		if (k >= 8)
			return 1;		/* Something weird happened */
		ppcol[i] = &pcol[k];
	}
	if (wv != 0)
		return 1;			/* Number is too big to be encoded */
//printf("~1dtp20_enc %d ->",rix);
//for (i = 0; i < ndig; i++)
//	printf(" %d",ppcol[i]->dtp20_octval);
//printf("\n");
	return 0;
}


/* Sort function for stree */
static int cmp_eperr(const void *p1, const void *p2) {
	return ((col *)p1)->wnd == ((col *)p2)->wnd ? 0 :
	                           (((col *)p1)->wnd < ((col *)p2)->wnd ? -1 : 1);
}

#define SYMWT 1.3	/* Amount to discount direction delta E compared to spacers */

/* Setup the randomised index. */
/* The index only covers test sample patches, not TID or max/min/SID patches */
void setup_randix(
int *rix,			/* Index lookup array to fill in */
int npat,			/* Number of test targets needed */
int rand,			/* Randomise flag */
int rstart,			/* Starting index for random */
int verb,			/* Verbose flag */
col *cols,			/* Array of colors to be put on target chart */
col *pcol,			/* 8 spacer colors */
int tpprow,			/* Test sample patches per row */
int spacer,			/* Spacer code, 0 = None, 1 = b&w, 2 = colored */
int needpc,			/* Need patch to patch contrast in a row */
int domaxmin,		/* Top and tail strip with max/min or SID.  0 = DTP51, 2 = DTP20 SID */
col *media,			/* Alias for media color */
int usede			/* NZ to use delta E rather than density */
) {
	int i;
	randix *r = NULL;	/* Random index order object */

	if (rand)
		r = new_randix(npat, rstart);

	/* Setup initial linear or randomised index */
	for (i = 0; i < npat; i++) {
		if (rand) {
			rix[i] = r->next(r);
		} else {
			rix[i] = i;
		}
		cols[rix[i]].ix = i;	 /* This colors random index */
	}
	rix[i] = 0;		/* npat+1 may be read, so fill extra entry at end. */

	if (rand)
		r->del(r);

	/* Setup initial contrast check */
	{
		col *pp, *cp, *np, *op;	/* Previous, current, next and opposite patch */
		col *maxd;			/* Alias for maximum density  */
		col *mind;			/* Alias for minimum density */
		aat_atree_t *stree;	/* Tree holding colors sorted by worst case contrast */
		aat_atrav_t *aat_tr;	/* Tree accessor */
		double temp, trate;	/* Annealing temperature & rate */
		double tstart, tend;/* Annealing chedule range */

		mind = &pcol[0];		/* White */
		maxd = &pcol[7];		/* Black */

		/* Create sorted tree so as to be able to locate the largest wnd */
		if ((stree = aat_anew(cmp_eperr)) == NULL)
			error("Allocating aat tree for colors failed");
		if ((aat_tr = aat_atnew()) == NULL)
			error("aat_atnew returned NULL");

		for (i = 0; i < npat; i++) {
			int j, k;		/* Patch index, Row index */
			int tpitr;		/* Test patches in this row */
			double tt;
	
			j = (i % tpprow);
			k = i / tpprow;

			/* Figure previous patch */
			if (j == 0) {			/* First in row */
				if (domaxmin == 1)
					pp = maxd;		/* Maxd will be before first patch */
				else if (domaxmin == 2) {
					col *ppcol[3];
					/* Compute the SID colors the DTP20 will use before row */
					if (dtp20_enc(pcol, 3, 1, ppcol, k+1) != 0)
						error("Internal, dtp20 SID row id failed, val %d, digits %d",k+1,3);
					pp = ppcol[2];	/* Barcode patch next to first test patch */
				} else
					pp = media;		/* Media will be before first patch */
			} else {
				pp = &cols[rix[i-1]];
			}

			/* Current patch */
			cp = &cols[rix[i]];

			/* Next patch */
			if (j == (tpprow-1) || i == (npat-1)) { /* Last in row or last patch */
				if ((j == (tpprow-1) || i == (npat-1)) && domaxmin == 1)
					np = mind;
				else if ((j == (tpprow-1) || i == (npat-1)) && domaxmin == 2)
					np = mind;		/* DTP20 has white trailing patch */
				else
					np = media;
			} else {
				np = &cols[rix[i+1]];
			}

			/* Opposite patch */
			
			if (k == npat/tpprow)	/* Last row */
				tpitr = npat - tpprow * (npat/tpprow);
			else
				tpitr = tpprow;
			op = &cols[rix[ k * tpprow + (tpitr -1 - j)]];

			/* Setup pointers and worst case contrast */
			cp->nc[0] = pp;
			cp->nc[1] = np;
			cp->oc = op;
			cp->wnd = setup_spacer(NULL, pp, cp, pcol, spacer, usede);
			tt      = setup_spacer(NULL, cp, np, pcol, spacer, usede);
			if (tt < cp->wnd)
				cp->wnd = tt;
			if (cp != op) {
				tt = SYMWT * density_difference(cp, op, usede);
				if (tt < cp->wnd)
					cp->wnd = tt;
			}

			/* Insert it into the sorted list */
			if ((aat_ainsert(stree, (void *)cp)) == 0)
				error("aat_ainsert color %d failed",i);
		}

		if (verb) {
			double wrdc = 1e300;

			if ((cp = aat_atfirst(aat_tr, stree)) == NULL)
				error("There seem to be no colors in the tree");
			if (usede)
				printf("Worst case delta E = %f\n", cp->wnd);
			else
				printf("Worst case density contrast = %f\n", cp->wnd);

			/* Evaluate each strips direction confusion */
			for (i = 0; ; i++) {
				int j, tpitr;		/* Test patches in this row */
				double tot = 0.0;
				if ((i * tpprow) >= npat)
					break;
				if (i == npat/tpprow)	/* Last row */
					tpitr = npat - tpprow * (npat/tpprow);
				else
					tpitr = tpprow;
	
				for (tot = 0.0, j = 0; j < tpitr; j++) {
					tot += density_difference(&cols[rix[ i * tpprow + j]],
					                          &cols[rix[ i * tpprow + (tpitr -1 - j)]], usede);
 				}
				tot /= tpitr;
				if (tot < wrdc)
					wrdc = tot;
			}
			if (usede)
				printf("Worst case direction distinction delta E = %f\n", wrdc);
			else
				printf("Worst case direction distinction density contrast = %f\n", wrdc);
		}

		if (needpc == 0 || !rand || npat < 3)
			return;		/* Current order is sufficient */

		if (verb) {
			printf("Optimising layout for strip reader:\n");
			printf(" 0%%"); fflush(stdout);
		}

		if (spacer == 2) {	/* Colored spacer, don't optimise too hard */
			tstart = 0.4;
			tend   = 0.00001;
			trate  = 0.87;
		} else {			/* No spacer or B&W spacer, do more optimisation */
			tstart = 0.4;
			tend   = 0.000005;
			trate  = 0.95;
		}
		/* Simulated anealing */
		for (temp = tstart; temp > tend; temp *= trate) {

			int ii, itlim;	/* Maximum passes at a temperature */
			int nsuc = 0;	/* Number that succeed */
			int suclim;		/* Number of successful changes before continuing */

			if (spacer == 2) {	/* Colored spacer, don't optimise too hard */
				itlim = npat * 10;
				suclim = npat;
			} else {
				itlim = npat * 14;
				suclim = npat;
			}

			if (verb) {		/* Output percent intervals */
				double pc;
	
				pc = (log(temp) - log(tstart))/(log(tend) - log(tstart));
				printf("%c%2d%%",cr_char,(int)(100.0 * pc+0.5)); fflush(stdout);
			}

			/* Improve the ordering */
			for (ii = 0; ii < itlim ; ii++) { 
				col *p1, *p2;
				double tt, de;

				/* Chose another patch to try swapping worst with */
				p1 = aat_atfirst(aat_tr, stree);	/* worst */
				for (;;) {
					tt = d_rand(0.0, 1.0);
					p2 = &cols[(int)(tt * (npat-1.0))];	/* Swap candidate */
					if (p1 != p2 && p2 != p1->oc)
						break;		/* Swap is not the worst or opposite */
				}

				/* Check p1 in p2's place */
				de = setup_spacer(NULL, p2->nc[0], p1, pcol, spacer, usede);
				tt = setup_spacer(NULL, p1, p2->nc[1], pcol, spacer, usede);
				if (tt < de)
					de = tt;
				tt = SYMWT * density_difference(p2->oc, p1, usede);
				if (tt < de)
					de = tt;

				/* Check p2 in p1's place */
				tt = setup_spacer(NULL, p1->nc[0], p2, pcol, spacer, usede);
				if (tt < de)
					de = tt;
				tt = setup_spacer(NULL, p2, p1->nc[1], pcol, spacer, usede);
				if (tt < de)
					de = tt;
				tt = SYMWT * density_difference(p1->oc, p2, usede);
				if (tt < de)
					de = tt;

				de = de - p1->wnd;		/* Increase in worst difference */

				/* If this swap will improve things, or temp is high enough, */
				/* then actually do the swap. */
				if (de > 0.0
				   || d_rand(0.0, 1.0) < exp(de/temp)) {
					int t;
					col *tp0, *tp1;

					nsuc++;

//printf("~1 temp = %f, ii = %d, swapping %d and %d\n",temp,ii,p1->i, p2->i);
 
					/* Remove them from the tree */
					if ((aat_aerase(stree, (void *)p1)) == 0)
						error("aat_aerase failed to find color  no %d", p1->i);
					if ((aat_aerase(stree, (void *)p2)) == 0)
						error("aat_aerase failed to find color  no %d", p2->i);

					/* Swap them in random index list */
					rix[p1->ix] = p2->i; 
					rix[p2->ix] = p1->i; 
					t = p1->ix; 
					p1->ix = p2->ix; 
					p2->ix = t;

					/* Swap their neighbors, taking care */
					/* of the situation if they are neigbors */
					tp0 = p1->nc[0];
					tp1 = p2->nc[0];
					if (tp0 == p1)
						tp0 = p2;
					else if (tp0 == p2)
						tp0 = p1;
					if (tp1 == p1)
						tp1 = p2;
					else if (tp1 == p2)
						tp1 = p1;
					p2->nc[0] = tp0;
					p1->nc[0] = tp1;

					tp0 = p1->nc[1];
					tp1 = p2->nc[1];
					if (tp0 == p1)
						tp0 = p2;
					else if (tp0 == p2)
						tp0 = p1;
					if (tp1 == p1)
						tp1 = p2;
					else if (tp1 == p2)
						tp1 = p1;
					p2->nc[1] = tp0;
					p1->nc[1] = tp1;

					/* Swap their opposites (they cannot be opposites of each other) */
					p1->oc->oc = p2;
					p2->oc->oc = p1;
					tp0 = p1->oc;
					p1->oc = p2->oc;
					p2->oc = tp0;

					/* Reset backwards references */
					p1->nc[0]->nc[1] = p1;
					p1->nc[1]->nc[0] = p1;
					p2->nc[0]->nc[1] = p2;
					p2->nc[1]->nc[0] = p2;

					/* re-compute contrast to neighbors */
					p1->wnd = setup_spacer(NULL, p1->nc[0], p1, pcol, spacer, usede);
					tt      = setup_spacer(NULL, p1, p1->nc[1], pcol, spacer, usede);
					if (tt < p1->wnd)
						p1->wnd = tt;
					if (p1 != p1->oc) {
						tt = SYMWT * density_difference(p1, p1->oc, usede);
						if (tt < p1->wnd)
							p1->wnd = tt;
					}

					p2->wnd = setup_spacer(NULL, p2->nc[0], p2, pcol, spacer, usede);
					tt      = setup_spacer(NULL, p2, p2->nc[1], pcol, spacer, usede);
					if (tt < p2->wnd)
						p2->wnd = tt;
					if (p2 != p2->oc) {
						tt = SYMWT * density_difference(p2, p2->oc, usede);
						if (tt < p2->wnd)
							p2->wnd = tt;
					}

					/* (!!! We haven't recomputed the possible change in the ->oc's */
					/* ->wnd due to it's opposite haveing changed. !!!) */

					/* Add them back to the tree */
					if ((aat_ainsert(stree, (void *)p1)) == 0)
						error("aat_ainsert color no %d failed",p1->i);
					if ((aat_ainsert(stree, (void *)p2)) == 0)
						error("aat_ainsert color no %d failed",p2->i);

#ifdef NEVER
printf("~1 current list = \n");
for (cp = aat_atfirst(aat_tr, stree); cp != NULL; cp = aat_atnext(aat_tr))
	printf("%d: %f\n",cp->i,cp->wnd);
cp = aat_atfirst(aat_tr, stree);
printf("~1 worst case contrast = %f, list index %d, id '%s', loc '%s'\n",
cp->wnd, cp->i,cp->id,cp->loc);
printf("~1 neighbor list index %d, id '%s', loc '%s'\n",
cp->nc[0]->i,cp->nc[0]->id,cp->nc[0]->loc);
printf("~1 neighbor list index %d, id '%s', loc '%s'\n",
cp->nc[1]->i,cp->nc[1]->id,cp->nc[1]->loc);
#endif

					if (nsuc > suclim)
						break;
				}
			}
		}

		if (verb) {
			double wrdc = 1e300;
			if ((cp = aat_atfirst(aat_tr, stree)) == NULL)
				error("There seem to be no colors in the tree");
			if (usede)
				printf("%c100%%\nAfter optimisation, worst delta E = %f\n",cr_char,cp->wnd);
			else
				printf("%c100%%\nAfter optimisation, density contrast = %f\n",cr_char,cp->wnd);

			/* Evaluate each strips direction confusion */
			for (i = 0; ; i++) {
				int j, tpitr;		/* Test patches in this row */
				double tot = 0.0;
				if ((i * tpprow) >= npat)
					break;
				if (i == npat/tpprow)	/* Last row */
					tpitr = npat - tpprow * (npat/tpprow);
				else
					tpitr = tpprow;
	
				for (tot = 0.0, j = 0; j < tpitr; j++) {
					tot += density_difference(&cols[rix[ i * tpprow + j]],
					                          &cols[rix[ i * tpprow + (tpitr -1 - j)]], usede);
 				}
				tot /= tpitr;
				if (tot < wrdc)
					wrdc = tot;
			}
			if (usede)
				printf("Worst case direction distinction delta E = %f\n", wrdc);
			else
				printf("Worst case direction distinction density contrast = %f\n", wrdc);
		}


		aat_atdelete(aat_tr);
		aat_adelete(stree);

	}

}

#define MAXPPROW 500		/* Absolute maximum patches per pass/row permitted */
#define MAXROWLEN 2000.0	/* Absolute maximum row length */

void
generate_file(
instType itype,		/* Target instrument type */
int itype_mod,		/* Target instrument type modifier */
char *bname,		/* Output file basename */
col *cols,			/* Array of colors to be put on target chart */
int npat,			/* Number of test targets needed */
xcal *cal,			/* Optional printer calibration, NULL if none */
char *label,		/* Per strip label */
double pw,			/* Page width */
double ph,			/* Page height */
double bord, 		/* Border margin in mm */
int nosubmarg,		/* NZ if bord is not to be subtracted from raster */
int nollimit,		/* NZ to not limit the strip length */
int nolpcbord,		/* NZ to suppress left paper clip border */
int rand,			/* Randomise flag */
int rstart,			/* Starting index for random */
alphix *saix,		/* Strip alpha index object */
alphix *paix,		/* Patch alpha index object */
int ixord,			/* Index order, 0 = strip then patch */
double pscale,		/* Test patch & spacers scale factor */
double sscale,		/* Spacers scale factor */
int hflag,			/* Spectroscan/Munki high density modified */
int verb,			/* Verbose flag */
int scanc,			/* Scan compatible bits, 1 = .cht gen, 2 = wide first row */
int oft,			/* PS/EPS/TIFF select (0,1,2) */
int nocups,			/* NZ to supress cups job ticket in PS/EPS */
depth2d tiffdpth,	/* TIFF pixel depth */
double tiffres,		/* TIFF resolution in DPI */
int ncha,			/* flag, use nchannel alpha */
int tiffdith,		/* flag, nz to use TIFF 8 bit dithering */
int tiffcomp,		/* flag, nz to use TIFF compression */
int spacer,			/* Spacer code, -1 = default, 0 = None, 1 = b&w, 2 = colored */
int nmask,			/* DeviceN mask */
int altrep,			/* printer grey/CMY representation type 0..8 */
col *pcol,			/* 8 spacer colors or 8 barcode colors for DTP20 */
double *wp,			/* Approximate white XYZ point */
int *ptpprow,		/* Return Test sample patches per row */
unsigned char **pprps,	/* Return malloced array holding Passes per strip */
double *p_patchlen,	/* Return patch length in mm */
double *p_gaplen,	/* Return gap length in mm */
double *p_taplen,	/* Return trailer length in mm */
int *p_npat			/* Return number of patches including padding */
) {
	char psname[MAXNAMEL+20];		/* Name of output file */
	trend *tro = NULL;		/* Target rendering object */
	double x1, y1, x2, y2;	/* Bounding box in mm */
	double iw, ih;			/* Imagable areas width and height in mm */
	double arowl;			/* Available row length */

	int hex = 0;			/* Hexagon patches flag (Spectrolino) */
	double stagger = 0.0;	/* Stagger alternate rows by half patch length */ 

	/* Chart definition variables. Set to default and override */
	/* for particular instruments. */
	double lbord = 0.0;	/* Additional left border */
	int domaxmin = 0;	/* if == 1, Top and tail strip with max and min values (DTP51) */
						/* if == 2, Top and tail with  DTP20 strip identification (SID) */
	int nmaxp = 0;		/* Number of max (header) patches for max/min/sid */
	int nminp = 0;		/* Number of min (trailer) patches for max/min/sid */
	int nextrap = 0;	/* Number of extra patches for max and min = nmaxp + nminp */
	int needpc = 0;		/* Need patch to patch contrast in a row */
	int dorspace = 0;	/* Do a rrsp from center of last patch to cut line & print label. */
	int dopglabel = 0;	/* Write a per page label */	
	double pglth = 0.0;	/* Page Label text height */
	int padlrow = 0;	/* flag - Pad the last row with white */
	double lspa = 0.0;	/* Leader space before first patch containint border, label, SID etc. */
	double lcar = 0.0;	/* Leading clear area before first patch. Will be white */
	double plen = 0.0;	/* Patch min length */
	double tidplen = 0.0;	/* TID Patch min length */
	double pspa = 0.0;	/* Patch to patch spacer */
	double tspa = 0.0;	/* Clear space after last patch */
	double txhi = 0.0;	/* Step/cut, row label text height */
	double txhisl = 0.0;/* Strip/column label text height */
	int docutmarks = 0;	/* Generate strip cut marks */
	double clwi = 0.0;	/* Cut line width */

	int dorowlabel = 0;	/* Generate a set of row labels */
	double rlwi = 0.0;	/* Row label test width */

	double hxew = 0.0;	/* Hexagon chart extra width padding to the right of patches */
	double hxeh = 0.0;	/* Hexagon chart extra height padding around patches */

	double pwid = 0.0;	/* Patch min width */
	double rrsp = 0.0; 	/* Row to row spacing */
	double pwex = 0.0;	/* Patch width expansion between rows of a strip */

	int mxpprow = 0;	/* Maximum patches per row permitted (including min/max patches) */
	double mxrowl = 0;	/* Maximum row length for patchs (including min/max patches) */
						/* Number of patches is strip by whichever is shorter */
	int tidrows = 0;	/* Rows on first page for target ID (ie. DTP20) */
	int tidtype = 0;	/* Target ID type. 0 = DTP20 */
	int tidminp = 0;	/* Target ID minumum number of patches */
	int tidpad = 0;		/* Initial padding to place TID near middle */
	int tidnpat = npat;	/* Number of test targets needed, including TID row */
	int rpstrip = 0;	/* Rows per strip */
	int usede = 0;		/* Use delta E to maximize patch/spacer conrast rather than density */

	double mints, minbs;	/* Minimum top & bottom space from paper edges */
	double amints, aminbs;	/* Actual mints & minbs, allowing for unused space */
	double swid;	/* Whole strip width */
	int pprow;		/* patches per row (inc. max/min/sid patches) */
	int tidpprow;	/* TID patches per row (inc. max/min/sid patches) */
	int tpprow;		/* test patches per row (excludes max/min/sid patches) */
	int sppage;		/* whole & partial strips per page */
	int rppstrip;	/* rows per partial strip on whole page */
	int npages;		/* Number whole & partial pages */
	int lsppage;	/* Last page whole & partial strips per page */
	int lrpstrip;	/* Last strips whole & partial rows per strip */
	int lpprow;		/* last row patches per row (inc. max/min/sid patches) */
	int ppstrip;	/* Real patches per whole strip */
	int pppage;		/* Real patches per whole page */
	int rem;		/* temporary */
	double sxwi = 0.0;	/* Scan compatible first row extra width */

	int *rix;	/* Random index lookup (Logical -> patch index) */
	int i;		/* Logical patch in target list */
	int l_si;	/* Last start of page value of i */
	int ix;		/* Patch index in target list */
	int pir;	/* Patch in row */
	int ris;	/* Row in strip */
	int sip;	/* Set in page (including partial strip) */
	int pif;	/* Page in file */
	double x = 0.0, y = 0.0;	/* Current position */
	col *pp = NULL, *cp = NULL;	/* Previous and current printed patch colors */
	int cpf = 0;				/* Current patch special flag. 1 = start bit, 2 = stop bit */
	int slix;	/* Strip label index, -1 for TID  */
	char *slab = NULL;	/* Strip label string */
	unsigned char *rpsp;	/* Rows per strip, pointer */
	col *mark;	/* Alias for mark color */
	col *media;	/* Alias for media color */
	col *maxd;	/* Alias for minimum density  */
	col *mind;	/* Alias for maximum density */
	col *sc;	/* Alias for current spacer color */

	/* Note pcol[] is setup by targen to be:
		0 = white
		1 = Cyan
		2 = Magenta
		3 = Blue
		4 = Yellow
		5 = Green
		6 = Red
		7 = Black
		8 = 50/50/50 CMY Grey
		(Should switch to symbols for these ??)
	 */

	/* We assume that since this is intended for a printer, */
	/* the media is always white. This may not be the case */
	/* on other output media. */
	if (itype == instDTP20)
		mark  = &pcol[8];		/* Grey */
	else
		mark  = &pcol[7];		/* Black */
	media = &pcol[0];		/* White */
	mind = &pcol[0];		/* White */
	maxd = &pcol[7];		/* Black */

	/* Setup DTP20 bar code encoding for spacers. */
	pcol[0].dtp20_octval = 0;	/* White */
	pcol[0].dtp20_psize  = 6.5;
	pcol[1].dtp20_octval = 4;	/* Cyan */
	pcol[1].dtp20_psize  = 10.0;
	pcol[2].dtp20_octval = 2;	/* Magenta */
	pcol[2].dtp20_psize  = 0.0;
	pcol[3].dtp20_octval = 6;	/* Blue */
	pcol[3].dtp20_psize  = 12.5;
	pcol[4].dtp20_octval = 1;	/* Yellow */
	pcol[4].dtp20_psize  = 7.0;
	pcol[5].dtp20_octval = 5;	/* Green */
	pcol[5].dtp20_psize  = 0.0;
	pcol[6].dtp20_octval = 3;	/* Red */
	pcol[6].dtp20_psize  = 0.0;
	pcol[7].dtp20_octval = 7;	/* Black */
	pcol[7].dtp20_psize  = 13.0;

	/* Setup .cht edge tracking information */
	et_init();
	et_height(oft != 2 ? mm2pnt(ph) : ph);
	et_media(media->rgb);

	/* Set Instrument specific parameters */
	if (itype == instDTP20) { 	/* Xrite DTP20 */
		hex = 0;				/* No hex for strip instruments */
		hxew = hxeh = 0.0;		/* No extra padding because no hex */
		domaxmin = 2;			/* Print SID patches */
		nmaxp = 4;				/* Extra header patches */
		nminp = 1;				/* Extra trailer patches */
		nextrap = nmaxp + nminp;/* Number of extra patches for max and min */
		dorspace = 0;			/* Maximise number of rows */
		dopglabel = 1;			/* Write a per page label */
		padlrow = 1;			/* Pad the last row with white */
		spacer = 0;				/* No Spacer between patches */
		pspa  = 0.0;			/* No spacer width */
		usede = 1;				/* Use delta E to maximize patch/spacer conrast */
		needpc = 1;				/* Helps to have patch to patch contrast in a row ? */
		lspa  = bord + 5.0 + 5.0;	/* Leader space before first patch = bord + pcar + yxhi */
		lcar  = 5.0;			/* Leading clear area before first patch */
		plen  = pscale * (6.5);	/* Patch min length */
		if(plen <= 6.75)		/* Patch length must be one of 5 lengths */
			plen = 6.5;
		else if(plen <= 8.0)
			plen = 7.0;
		else if(plen <= 11.25)
			plen = 10.0;
		else if(plen <= 12.75)
			plen = 12.5;
		else 
			plen = 13.0;
		tidplen = 6.0;			/* TID Patch length. Can't vary. */
		tspa  = 5.0;			/* Clear space after last patch */
		pwid  = 10.0;			/* Patch min width. (The guide slot is 12mm ?) */
		if (plen > pwid)
			pwid = plen;		/* Make patch at least as wide as long */
		rrsp  = pwid;			/* Row center to row center spacing */
		pwex  = 0.0;			/* Patch width expansion between rows of a strip */
		if (nollimit == 0) {
			mxpprow = MAXPPROW;		/* Maximum patches per row permitted (set by length) */
			mxrowl = (240.0 - lcar - tspa);	/* Maximum row length */
		} else {
			mxpprow = MAXPPROW;		/* Maximum patches per row permitted (set by length) */
			mxrowl = MAXROWLEN;				/* No row length */
		}
		tidrows = 1;			/* Rows on first page for target ID */
		tidtype = 0;			/* Target ID type. 0 = DTP20 */
		tidminp = 21;			/* Target ID minumum number of patches (+ white) */
		rpstrip = 999;			/* Rows per strip */
		txhi   = 5.0;			/* Label Text Height */
		txhisl = 2.0;			/* Strip Label Text Height */
		docutmarks = 0;			/* Don't generate strip cut marks */
		clwi  = 0.0;			/* Cut line width */
		dorowlabel = 0;			/* Don't generate row labels */
		rlwi = 0.0;				/* Row label width */
		pglth = 5.0;			/* Page Label text height */

		if (npat > 4095)
			error ("Number of patchs %d exceeds maximum of 4095 for DTP20");

	} else if (itype == instDTP22 ) {	/* X-Rite DTP22 Digital Swatchbook */
		hex = hflag ? 1 : 0;	/* Hex if requestested */
		domaxmin = 0;			/* Don't print max and min patches */
		nextrap = 0;			/* Number of extra patches for max and min */
		nmaxp = nminp = 0;		/* Extra max/min patches */
		nextrap = nmaxp + nminp;/* Number of extra patches for max and min */
		dorspace = 0;			/* Do a rrsp from center of last patch to cut line */
		dopglabel = 1;			/* Write a per page label */
		padlrow = 0;			/* Pad the last row with white */
		spacer = 0;				/* No spacer */
		usede = 1;				/* Use delta E to maximize patch/spacer conrast */
		needpc = 1;				/* Need patch to patch contrast in a row */
		lspa  = bord + 8.0;		/* Leader space before first patch = border + text */
		lcar  = 0.0;			/* Leading clear area before first patch */
		if (hex) {
			plen = pscale * sqrt(0.75) * 8.0;	/* Patch min length */
			hxeh = 1.0/6.0 * plen;				/* Extra border for hex tops & bottoms */
			hxew = pscale * 0.25 * 8.0;			/* Extra border for hex sides */
		} else {
			plen = pscale * 8.0;	/* Patch min length */
			hxew = hxeh = 0.0;		/* No extra padding because no hex */
		}
		pspa  = 0.0;			/* Inbetween Patch spacer */
		tspa  = 0.0;			/* Clear space after last patch */
		pwid  = pscale * 8.0;	/* Patch min width */
		rrsp  = pscale * 8.0;	/* Row center to row center spacing */
		pwex  = 0.0;			/* Patch width expansion between rows of a strip */
		mxpprow = MAXPPROW;		/* Maximum patches per row permitted */
		mxrowl = MAXROWLEN;		/* Maximum row length */
		tidrows = 0;			/* No rows on first page for target ID */
		rpstrip = 999;			/* Rows per strip */
		txhi = txhisl = 5.0;	/* Text Height */
		docutmarks = 0;			/* Don't generate strip cut marks */
		clwi  = 0.0;			/* Cut line width */
		dorowlabel = 1;			/* Generate row labels */
		rlwi = 8.0;				/* Row label width */
		pglth = 5.0;			/* Page Label text height */

	} else if (itype == instDTP41) {	/* Xrite DTP41 */
		hex = 0;				/* No hex for strip instruments */
		hxew = hxeh = 0.0;		/* No extra padding because no hex */
		domaxmin = 0;			/* Don't print max and min patches */
		nmaxp = nminp = 0;		/* Extra max/min patches */
		nextrap = nmaxp + nminp;/* Number of extra patches for max and min */
		dorspace = 0;			/* Maximise number of rows */
		dopglabel = 1;			/* Write a per page label */
		padlrow = 1;			/* Pad the last row with white */
		if (spacer < 0)
			spacer = 2;			/* Colored Spacer */
		needpc = 1;				/* Need patch to patch contrast in a row */
		lspa  = inch2mm(1.5);	/* Leader space before first patch */
		lcar  = inch2mm(0.5);	/* Leading clear area before first patch */
		plen  = pscale * inch2mm(0.29);	/* Patch min length (should be 7.0 mm min.) */
		if (spacer > 0)
			pspa  = pscale * sscale * inch2mm(0.08);	/* Inbetween Patch spacer (should be 2.0 mm min.)*/
		else
			pspa  = 0.0;		/* No spacer */
		tspa  = 2 * (plen + pspa);	/* Clear space after last patch */
		pwid  = inch2mm(0.5);	/* Patch min width */
		rrsp  = inch2mm(0.5);	/* Row center to row center spacing */
		pwex  = (rrsp - pwid)/2.0;	/* Patch width expansion between rows of a strip */
		if (nollimit == 0) {
			mxpprow = 100;			/* Maximum patches per row permitted */
			mxrowl = inch2mm(55.0);	/* Maximum row length */
		} else {
			mxpprow = MAXPPROW;		/* Maximum patches per row permitted */
			mxrowl = MAXROWLEN;		/* Maximum row length */
		}
		tidrows = 0;			/* No rows on first page for target ID */
		rpstrip = 8;			/* Rows per strip */
		txhi = txhisl = 5.0;	/* Text Height */
		docutmarks = 1;			/* Generate strip cut marks */
		clwi  = 0.3;			/* Cut line width */
		dorowlabel = 0;			/* Don't generate row labels */
		rlwi = 0.0;				/* Row label width */
		pglth = 5.0;			/* Page Label text height */

	} else if (itype == instDTP51) { 	/* Xrite DTP51 */
		hex = 0;				/* No hex for strip instruments */
		hxew = hxeh = 0.0;		/* No extra padding because no hex */
		domaxmin = 1;			/* Print max and min patches */
		nmaxp = nminp = 1;		/* Extra max/min patches */
		nextrap = nmaxp + nminp;/* Number of extra patches for max and min */
		dorspace = 1;			/* Do a rrsp from center of last patch to cut line */
		dopglabel = 0;			/* No need for a per page label */
		padlrow = 1;			/* Pad the last row with white */
		if (spacer < 0)
			spacer = 2;			/* Colored Spacer */
		needpc = 1;				/* Need patch to patch contrast in a row */
		lspa  = inch2mm(1.2);	/* Leader space before first patch */
		lcar  = inch2mm(0.25);	/* Leading clear area before first patch */
		plen  = pscale * inch2mm(0.4);	/* Patch min length */
		if (spacer > 0)
			pspa  = pscale * sscale * inch2mm(0.07);	/* Inbetween Patch spacer */
		else
			pspa  = 0.0;		/* No spacer */
		tspa  = inch2mm(0.0);	/* Clear space after last patch */
		pwid  = inch2mm(0.4);	/* Patch min width */
		rrsp  = inch2mm(0.5);	/* Row center to row center spacing */
		pwex  = (rrsp - pwid)/2.0;	/* Patch width expansion between rows of a strip */
		if (nollimit == 0) {
			mxpprow = 72;			/* Maximum patches per row permitted */
			mxrowl = inch2mm(40.0);	/* Maximum row length */
		} else {
			mxpprow = MAXPPROW;		/* Maximum patches per row permitted */
			mxrowl = MAXROWLEN;		/* Maximum row length */
		}
		tidrows = 0;			/* No rows on first page for target ID */
		rpstrip = 6;			/* Rows per strip */
		txhi = txhisl = 5.0;	/* Text Height */
		docutmarks = 1;			/* Generate strip cut marks */
		clwi  = 0.3;			/* Cut line width */
		dorowlabel = 0;			/* Don't generate row labels */
		rlwi = 0.0;				/* Row label width */
		pglth = 5.0;			/* Page Label text height */

	} else if (itype == instSpectroScan ) {	/* GretagMacbeth SpectroScan */
		hex = hflag ? 1 : 0;	/* Hex if requestested */
		domaxmin = 0;			/* Don't print max and min patches */
		nextrap = 0;			/* Number of extra patches for max and min */
		nmaxp = nminp = 0;		/* Extra max/min patches */
		nextrap = nmaxp + nminp;/* Number of extra patches for max and min */
		dorspace = 0;			/* Maximise number of rows */
		dopglabel = 1;			/* Write a per page label */
		padlrow = 0;			/* Pad the last row with white */
		spacer = 0;				/* No spacer */
		needpc = 0;				/* Don't need patch to patch contrast in a row */
		lspa  = bord + 7.0;		/* Leader space before first patch = border + text */
		lcar  = 0.0;			/* Leading clear area before first patch */
		if (hex) {
			plen = pscale * sqrt(0.75) * 7.0;	/* Patch min length */
			hxeh = 1.0/6.0 * plen;				/* Extra border for hex tops & bottoms */
			hxew = pscale * 0.25 * 7.0;			/* Extra border for hex sides */
		} else {
			plen = pscale * 7.0;	/* Patch min length */
			hxew = hxeh = 0.0;		/* No extra padding because no hex */
		}
		pspa  = 0.0;			/* Inbetween Patch spacer */
		tspa  = 0.0;			/* Clear space after last patch */
		pwid  = pscale * 7.0;	/* Patch min width */
		rrsp  = pscale * 7.0;	/* Row center to row center spacing */
		pwex  = 0.0;			/* Patch width expansion between rows of a strip */
		mxpprow = MAXPPROW;		/* Maximum patches per row permitted */
		mxrowl = MAXROWLEN;		/* Maximum row length */
		tidrows = 0;			/* No rows on first page for target ID */
		rpstrip = 999;			/* Rows per strip */
		txhi = txhisl = 5.0;	/* Row/Column Text Height */
		docutmarks = 0;			/* Don't generate strip cut marks */
		clwi  = 0.0;			/* Cut line width */
		dorowlabel = 1;			/* Generate row labels */
		rlwi = 7.5;				/* Row label width */
		pglth = 5.0;			/* Page Label text height */

	} else if (itype == instI1Pro ) {	/* GretagMacbeth/X-Rite Eye-One Pro */
		if (nolpcbord == 0 && bord < 26.0)
			lbord = 26.0 - bord;	/* need this for holder to grip paper and plastic spacer */
		hex = 0;				/* No hex for strip instruments */
		hxew = hxeh = 0.0;		/* No extra padding because no hex */
		domaxmin = 0;			/* Don't print max and min patches */
		nextrap = 0;			/* Number of extra patches for max and min */
		nmaxp = nminp = 0;		/* Extra max/min patches */
		nextrap = nmaxp + nminp;/* Number of extra patches for max and min */
		dorspace = 0;			/* Maximise number of rows by having no space between them */
		dopglabel = 1;			/* Write a per page label */
		padlrow = 1;			/* Don't need to pad the last row for the i1, */
								/* but the strip read logic can't handle it. */
		if (spacer < 0)
			spacer = 2;			/* Colored Spacer */
		needpc = 1;				/* Need patch to patch contrast in a row */
		usede = 1;				/* Use delta E to maximize patch/spacer conrast */

		/* 5 mm aperture */
		if (itype_mod == 0) {
			lcar  = 10.0;			/* Leading clear area before first patch */
			plen  = pscale * 10.00;	/* Patch min length - total 10 mm (absolute limit is 10) */
			if (spacer > 0)
				pspa  = pscale * sscale * 1.00;	/* Inbetween Patch spacer 1mm */
			else
				pspa  = 0.0;		/* No spacer */
			tspa  = 10.0; 			/* Clear space after last patch - run off */
			pwid  = pscale * 8.0;	/* Patch min width (absollute limit is 7) */
			rrsp  = pscale * 8.0;	/* Row center to row center spacing */
			pwex  = 0.0;			/* Patch width expansion between rows of a strip */

		/* 8 mm aperture */
		} else {
			lcar  = 20.0;			/* Leading clear area before first patch */
			plen  = pscale * 20.00;	/* Patch min length - total 20 mm (absolute limit 16 with zr) */
			if (spacer > 0)
				pspa  = pscale * sscale * 2.0;	/* Inbetween Patch spacer 2.0 mm */
			else
				pspa  = 0.0;		/* No spacer */
			tspa  = 20.0; 			/* Clear space after last patch - run off */
			pwid  = pscale * 16.0;	/* Patch min width (absollute limit is 14) */
			rrsp  = pscale * 16.0;	/* Row center to row center spacing */
			pwex  = 0.0;			/* Patch width expansion between rows of a strip */
		}
		txhi = txhisl = 7.0;		/* Text Height */
		rlwi = 0.0;					/* Row label width */
		pglth = 5.0;				/* Page Label text height */
		lspa = bord + txhisl + lcar;	/* Leader space before first patch */

		if (nollimit == 0) {
			mxpprow = MAXPPROW;		/* Maximum patches per row permitted (set by length) */
			mxrowl = (260.0 - lcar - tspa);	/* Maximum holder row length */
		} else {
			mxpprow = MAXPPROW;		/* Maximum */
			mxrowl = MAXROWLEN;		/* Maximum */
		}
		tidrows = 0;			/* No rows on first page for target ID */
		rpstrip = 999;			/* Rows per strip */
		docutmarks = 0;			/* Don't generate strip cut marks */
		clwi  = 0.0;			/* Cut line width */
		dorowlabel = 0;			/* Don't generate row labels */

	} else if (itype == instColorMunki ) {	/* X-Rite ColorMunki */
		hex = 0;				/* No hex for strip instruments */
		hxew = hxeh = 0.0;		/* No extra padding because no hex */
		domaxmin = 0;			/* Don't print max and min patches */
		nextrap = 0;			/* Number of extra patches for max and min */
		nmaxp = nminp = 0;		/* Extra max/min patches */
		nextrap = nmaxp + nminp;/* Number of extra patches for max and min */
		dorspace = 0;			/* Put spaces between rows for guidance */
		dopglabel = 1;			/* Write a per page label */
		padlrow = 1;			/* Don't need to pad the last row for the Munki, */
								/* but the strip read logic can't handle it. */
		if (spacer < 0)
			spacer = 2;			/* Colored Spacer */
		needpc = 1;				/* Need patch to patch contrast in a row */
		usede = 1;				/* Use delta E to maximize patch/spacer conrast */
		lspa  = bord + 7.0 + 20.0;	/* Leader space before first patch = bord + txhisl + lcar */
		lcar  = 20.0;			/* Leading clear area before first patch */
		plen  = pscale * 14.00;	/* Patch min length - total 15 mm */
		if (spacer > 0)
			pspa  = pscale * sscale * 1.0;	/* Inbetween Patch spacer 1mm */
		else
			pspa  = 0.0;		/* No spacer */
		tspa  = 25.0; 			/* Clear space after last patch - run off */
		if (hflag) {			/* High density */
			pwid  = pscale * 13.7;	/* Patch min width */
			rrsp  = pscale * 13.7;	/* Row center to row center spacing */
			hxeh = 0.25 * plen;		/* Extra space for row stagger */
			stagger = 0.5 * (plen + 0.5 * pspa);	/* Do stagger */
		} else {
			pwid  = pscale * 28.0;	/* Patch min width */
			rrsp  = pscale * 28.0;	/* Row center to row center spacing */
		}
		pwex  = 0.0;			/* Patch width expansion between rows of a strip */
		mxpprow = MAXPPROW;		/* Maximum patches per row permitted (set by length) */
		mxrowl = MAXROWLEN;		/* Maximum row length */
		tidrows = 0;			/* No rows on first page for target ID */
		rpstrip = 999;			/* Rows per strip */
		txhi = txhisl = 7.0;	/* Text Height */
		docutmarks = 0;			/* Don't generate strip cut marks */
		clwi  = 0.0;			/* Cut line width */
		dorowlabel = 0;			/* Don't generate row labels */
		rlwi = 0.0;				/* Row label width */
		pglth = 5.0;			/* Page Label text height */


	} else {
		error("Unsupported intrument type");
	}
	
	/* Compute page limits */
	x1 = bord + lbord;	/* Bounding box in mm */
	y1 = bord;
	x2 = pw - bord;
	y2 = ph - bord;
	iw = x2 - x1;		/* Imagable areas width and height in mm */
	ih = y2 - y1;

	*p_patchlen = plen;		/* Return patch lenth in mm */
	*p_gaplen = pspa;		/* Return gap lenth in mm */
	*p_taplen = tspa;		/* Return trailer lenth in mm */

	if (scanc & 2) 			/* Scan compatiblity */
		sxwi = pwid/2.0;	/* First row patches extra width */

	/* Compute limits for this page size */
	/* Figure the available space for patches */
	mints = bord + txhisl + lcar;		/* Minimum top space due to border, text and clear area */
	if (mints < lspa)
		mints = lspa;				/* Minimum top space due to leader */
	minbs = bord;					/* Minimum botom space due to border */
	if (minbs < tspa)
		minbs = tspa;				/* Minimum botom space due to trailer */
	arowl = ph - mints - minbs - 2.0 * hxeh;	/* Available space for printing test patches */
	if (arowl > mxrowl)
		arowl = mxrowl;				/* Limit maximum row length */

	/* We are assuming that every patch may be surounded by a spacer */
	/* (ie. there are always pprow+1 gaps spacers are used) */
	pprow = (int)((arowl - pspa)/(plen + pspa));	/* Raw Patches per row */
	if (pprow > mxpprow)			/* Limit to maximum */
		pprow = mxpprow;
	
	tidpprow = 0;
	if (tidminp > 0 && tidplen > 0.0) {
		tidpprow = (int)((arowl - pspa)/(tidplen + pspa));	/* Raw TID Patches per row */
		if (tidpprow < tidminp)			/* TID doesn't use nextrap */
			error("Paper size not long enough for target identification row (need %.1f mm, got %.1f mm)!",tidminp * (tidplen + pspa) - pspa, arowl);
	}

	tidpad = (pprow - tidminp)/2;	/* Center TID */

	if (pprow < (1+nextrap))
		error("Paper size not long enought for a single patch per row!");

	*ptpprow = tpprow = pprow - nextrap;	/* Test sample patches per row */

	tidnpat = npat + (tidrows * tpprow);	/* Total patches including TID row, but not max/min/sid */

	/* (Sample patches per row including TID) */
	if ((*pprps = (unsigned char *)malloc(sizeof(unsigned char) * (2 + (tidnpat/tpprow)))) == NULL)
		error("Malloc failed!");
	rpsp = *pprps;
	
	/* Compute actual lowest coordinate used */
	aminbs = ph - mints - pspa - pprow * (plen + pspa);
	amints = mints + 0.5 * (aminbs - minbs); 	/* Distribute extra space */
	aminbs = minbs + 0.5 * (aminbs - minbs);

	/* Compute whole strip width */
	if (dorspace)
		swid = rpstrip * rrsp + pwid/2.0;			/* set gutter is rrsp - pwid/2 wide */
	else
		swid = (rpstrip-1) * rrsp + pwid + clwi;	/* set gutter is 0, but allow for cut line */
		
	/* Compute strips per page.  Number of whole strips + partial strips */
	sppage = (int)((iw - rlwi - sxwi - 2.0 * hxew - (dopglabel ? pglth : 0.0))/swid) + 1;

	/* Compute rows per partial strip on whole page */
	if (dorspace)
		rppstrip = (int)((iw - rlwi - sxwi - 2.0 * hxew - (dopglabel ? pglth : 0.0) - swid * (sppage-1) - pwid/2.0)/rrsp);
	else
		rppstrip = (int)((iw - rlwi - sxwi - 2.0 * hxew - (dopglabel ? pglth : 0.0) - swid * (sppage-1) - pwid + rrsp)/rrsp);
	if (rppstrip < 0)
		rppstrip = 0;
	if (rppstrip == 0) {	/* Make last partial strip a full strip */
		sppage--;
		rppstrip = rpstrip;
	}
	
	if (sppage <= 0)
		error("Not enough width for even one row!");

	/* The number of pages needed */
	pppage = tpprow * ((sppage-1) * rpstrip + rppstrip);/* Real patches per page */
	npages = (tidnpat + pppage -1)/pppage;				/* whole & partial pages */
	ppstrip = tpprow * rpstrip;							/* Real patches per full strip */

	rem = tidnpat;						/* Total test patches to print */
	rem -= (npages-1) * pppage;			/* Remaining patches to be printed on last page */

	lsppage = (rem + ppstrip -1)/ppstrip;	/* Last pages whole & partial strips per page */
	rem -= (lsppage - 1) * ppstrip;			/* remaining patches to be printed in last strip */

	lrpstrip = (rem + tpprow - 1)/tpprow;
										/* Last strips whole & partial rows per strip */

	rem -= (lrpstrip - 1) * tpprow;		/* remaining patches to be printed in last row */

	lpprow = rem + nextrap;				/* Patches in last row of last strip of last page */

	if (verb) {
		fprintf(stderr,"Patches = %d\n",npat);
		fprintf(stderr,"Test patches per row = %d\n",tpprow);
		if (sppage == 1)
			fprintf(stderr,"Rows per page = %d, patches per page = %d\n",rppstrip, pppage);
		else
			fprintf(stderr,"Strips per page = %d, rows per partial strip = %d, patches per page = %d\n",sppage, rppstrip, pppage);
		if (tidrows > 0)
			fprintf(stderr,"Target ID rows in first page = %d\n", tidrows);
		fprintf(stderr,"Rows in last strip = %d, patches in last row = %d\n", lrpstrip, lpprow-nextrap);
		fprintf(stderr,"Total pages needed = %d\n",npages);
	}

	if (padlrow) {		/* Add in extra padding patches */
		int i;
		for (i = 0; lpprow < pprow; lpprow++, npat++, tidnpat++, i = (i + 1) & 7) {
#ifdef NEVER
			if (needpc && rand)
				cols[npat] = pcol[i];		/* structure copy */
			else
				cols[npat] = *media;		/* structure copy */
#else
			cols[npat] = *media;		/* structure copy */
#endif
			cols[npat].i = npat;			/* Now has an index */
			cols[npat].t &= ~T_PRESET;
			cols[npat].t |= T_PAD; 
			cols[npat].id = "0";			/* Padding identification */
		}
	}
	*p_npat = npat;			/* Return number of patches including padding */
//printf("~1 padded no of patches = %d\n", npat);

	/* Setup logical to random patch mapping */
	if ((rix = (int *)malloc(sizeof(int) * (npat + 1))) == NULL)
		error("Malloc failed!");

	setup_randix(rix, npat, rand, rstart, verb, cols, pcol,
	             tpprow, spacer, needpc, domaxmin, media, usede);
	rix[npat] = -1;		/* Shouldn't use this */

	/* Init everything */
	l_si = i = 0;	/* Physical test patch index. */
	ix = rix[i];	/* First index */

	pir = 1;		/* Starting patch in row (includes max/min patches) */
	ris = 1;		/* Starting row in strip. */
	sip = 1;		/* Starting strip in page */
	pif = 1;		/* Starting page in file */

						/* slix is 0..n but is pre-incremented, so start at -1 */
	slix = -1 -tidrows;	/* Start at -2 if there is a TID  */

	/* Until there are no more patches to do */
	for (;;) {
		char *sp = NULL;	/* String pointer - label */
		double w;			/* Width of next patch */
		int flags;			/* flags for current patch */
#define IS_FPIR   0x00001	/* Is first patch in row (possibly max density patch) */
#define IS_XPAT   0x00002	/* Is max density/SID patch */
#define IS_NPAT   0x00004	/* Is min density/SID patch */
#define IS_LPIR   0x00008	/* Is last patch in row (possibly min density patch) */
#define IS_FRIS   0x00010	/* Is first row in strip */
#define IS_LRIS   0x00020	/* Is last row in strip */
#define IS_FSIP   0x00040	/* Is first strip in page */
#define IS_LSIP   0x00080	/* Is last strip in page */
#define IS_FPIF   0x00100	/* Is first page in file */
#define IS_LPIF   0x00200	/* Is last page in file */
#define IS_PAD    0x04000	/* Is fake padding patch, to round out very last row */

		/* Init flags for this patch */
		flags = 0;
		if (pir == 1)
			flags |= IS_FPIR;
		if (pir <= nmaxp)
			flags |= IS_XPAT;

		if (slix < -1) {		/* If TID row */

			if (pir == tidpprow)
				flags |= IS_LPIR;

		} else {

			if (pir > (pprow - nminp))
				flags |= IS_NPAT;
	
			if (pir == pprow)
				flags |= IS_LPIR;
		}

		if (ris == 1)
			flags |= IS_FRIS;
		if (ris == rpstrip)
			flags |= IS_LRIS;

		if (sip == 1)
			flags |= IS_FSIP;
		if (sip == sppage) {
			flags |= IS_LSIP;
			if (ris == rppstrip)
				flags |= IS_LRIS;
		}
		if (pif == 1)
			flags |= IS_FPIF;
		if (pif == npages) {
			flags |= IS_LPIF;
			if (sip == lsppage) {
				flags |= IS_LSIP;
				if (ris == lrpstrip) {
					flags |= IS_LRIS;

					if (padlrow) {				/* Last row in chart may be a runt */
						if (pir > lpprow)
							flags |= IS_PAD;
					} else {
						if (pir > (lpprow - nminp))
							flags |= IS_NPAT;
						if (pir == lpprow)
							flags |= IS_LPIR;
					}
				}
			}
		}

//printf("~1 pir %d, ris %d, sip %d, pif %d\n", pir, ris, sip, pif);

		/* Set initial patch width */
		w = pwid;
		if (!(flags & IS_FRIS))
			w += pwex;		/* Extend patch into previous row to row gap */
		if (!(flags & IS_LRIS))
			w += pwex;		/* Extend patch into next row to row gap */

		if ((flags & IS_FRIS) && (flags & IS_FSIP))	/* First row on page */
			w += sxwi;		/* Make first row fatter for scan compatiblity */

		if (flags & IS_FPIR) {			/* Start of row */
			y = ph - amints;			/* Start ready for spacer or patch */ 

			if (flags & IS_FRIS) {							/* Start of strip */
				if (flags & IS_FSIP) {						/* Start of page */
					x = x1;									/* Start at leftmost position */
					if (oft == 0) {	/* PS */
						if (flags & IS_FPIF) {					/* First page */
							sprintf(psname,"%s.ps",bname);
							if ((tro = new_ps_trend(psname,npages,nmask,pw,ph,oft,nocups,rand,rstart)) == NULL)
								error ("Unable to create output rendering object file '%s'",psname);
							if (verb)
								printf("Creating file '%s'\n",psname);
						}
					} else if (oft == 1) {	/* EPS */
						if (npages > 1)
							sprintf(psname,"%s_%02d.eps",bname,pif);
						else
							sprintf(psname,"%s.eps",bname);
						if ((tro = new_ps_trend(psname,npages,nmask,pw,ph,oft,rand,rstart,nocups)) == NULL)
							error ("Unable to create output rendering object file '%s'",psname);
						if (verb)
							printf("Creating file '%s'\n",psname);
					} else {	/* TIFF */
						double res;		/* pix/mm */
						if (npages > 1)
							sprintf(psname,"%s_%02d.tif",bname,pif);
						else
							sprintf(psname,"%s.tif",bname);

						res = tiffres/25.4; 
						if ((tro = new_tiff_trend(psname,nmask,tiffdpth,pw,ph,
						    nosubmarg ? 0 : bord, res,res,altrep,ncha,tiffcomp, tiffdith)) == NULL)
							error ("Unable to create output rendering object file '%s'",psname);
						if (verb)
							printf("Creating file '%s'\n",psname);
					}
					tro->startpage(tro,pif);

					/* Print all the row labels */
					if (dorowlabel) {
						double ty = y;		/* Temp y coord */
						int tpir;			/* Temp patch in row */

						for (tpir = 0; tpir < pprow; tpir++) {
							/* If we're within test sample patch range */
							if (tpir >= nmaxp && tpir < (pprow - nminp)) {
								char *rlabl;
								if ((rlabl = paix->aix(paix, tpir - nmaxp)) == NULL)
									error ("Patch in row label %d out of range",tpir);
								tro->setcolor(tro, cal, mark);
								tro->string(tro, x, ty-plen, rlwi, plen, rlabl);
								free(rlabl);
							}
						
							ty -= plen + pspa;
						}

						x += rlwi;
					}

					x += hxew; /* Extra space on left for bits of hex */

					/* Clear edge list tracking */
					et_clear();
				}
			}

			/* Increment strip label 0..n */
			slix++;

			/* Print strip label */
			if ((lspa - lcar - bord) >= txhisl) {	/* There is room for label */
				if (slix < 0) {	/* TID */
					tro->setcolor(tro, cal, mark);
					tro->string(tro,x,y2-txhisl,w,txhisl,"TID");
				} else {		/* Not TID */
					if (slab != NULL)
						free(slab);
					if ((slab = saix->aix(saix, slix)) == NULL)
						error("strip index %d out of range",slix);
		
					tro->setcolor(tro, cal, mark);
					tro->string(tro,x,y2-txhisl,w,txhisl,slab);
				}
			}

			/* Start with background = media */
			pp = media;

			/* Stagger rows */
			if (stagger > 0.0) {
				if (slix & 1)
					y -= stagger;
			}
		}

		/* Figure the current patch color */
		cpf = 0;
		if (slix < 0) {		/* TID */
			cp = mind;		/* Default padding color is white */
			if (pir > tidpad && pir <= (tidpad + tidminp) ) {
				int opir = pir - tidpad -1;	/* TID index 0 .. tidminp-1 */
				
				if (opir == 0) {
					cp = &pcol[1];	/* Cyan */

				} else if (opir >= 1 && opir <= 2) {	/* Patches in each strip */
					col *ppcol[2];
					if (dtp20_enc(pcol, 2, 0, ppcol, tpprow) != 0)
						error("Internal, dtp20 TID ppstrip failed, val %d, digits %d",tpprow,2);
					cp = ppcol[opir-1];

				} else if (opir >= 3 && opir <= 6) {	/* Total patches, including padding */
					col *ppcol[4];
					if (dtp20_enc(pcol, 4, 0, ppcol, npat) != 0)
						error("Internal, dtp20 tot patches failed, val %d, digits %d",npat,4);
					cp = ppcol[opir-3];

				} else if (opir == 7) {	/* Patch size */
					int j;
					for (j = 0; j < 8; j++) {
						if (fabs(pcol[j].dtp20_psize - plen) < 0.001)
							break;
					}
					if (j >= 8)
						error("Can't encode patch length for DTP20");
					cp = &pcol[j];

				} else if (opir == 8) {	/* Patch spacer width */
					int j, k;
					k = (int)(pspa / 0.5 + 0.5);
					for (j = 0; j < 8; j++) {
						if (pcol[j].dtp20_octval == k)
							break;
					}
					if (j >= 8)
						error("Can't encode spacer length for DTP20");
					cp = &pcol[j];

				} else if (opir >= 9 && opir <= 18) {	/* User defined */

					if (opir == 9) {			/* Indicate what user defined is used for */
						cp = &pcol[0];			/* 0 = random seed format */
					} else if (opir >= 10 && opir <= 13) {	/* Random seend value */
						col *ppcol[4];
						if (dtp20_enc(pcol, 4, 0, ppcol, rstart) != 0)
							error("Internal, dtp20 chart id failed, val %d, digits %d",rstart,4);
						cp = ppcol[opir-10];

					} else {	/* 14 .. 18 */
						cp = &pcol[0];		/* Currently unused */
					}

				} else if (opir == 19) {
					cp = &pcol[4];	/* Yellow */
				} else if (opir == 20) {
					cp = &pcol[2];	/* Magenta */
				}
			}
		} else {
			if (flags & IS_XPAT) {						/* Max or SID at start of row */
				sp = NULL;								/* Not a test patch (no label) */
				if (domaxmin == 1) {
					cp = maxd;							/* Maximum density patch at start */
				} else if (domaxmin == 2) {
					if (pir == 1) {						/* At very start */
						cpf = 1;						/* Create start bit */
						cp = &pcol[8];					/* Starts with 50/50/50 DTP20 Grey */
					} else {
						col *ppcol[3];
						/* Compute the patch colors the DTP20 will use before row */
						if (dtp20_enc(pcol, 3, 1, ppcol, slix+1) != 0)
							error("Internal, dtp20 SID row id failed, val %d, digits %d",slix+1,3);
						cp = ppcol[pir-2];	
					}
				}

			} else if (flags & IS_NPAT) { 				/* Min at end of rows or stop bit */
				sp = NULL;								/* Not a test patch (no label) */
				if (domaxmin == 1) {
					cp = mind;
				} else if (domaxmin == 2) {				/* DTP20 end patch */
					cpf = 2;							/* Create stop bit */
					cp = mind;							/* Starts with mind */
				}

			} else if (flags & IS_PAD) {				/* Fake blank patch */
				sp = NULL;								/* Not a test patch (no label) */
				cp = media;

			} else {	/* set test sample patch location and color */
				int apir = pir - nmaxp;			/* Adjusted pir for max/min patches */
				if (sp != NULL)
					free(sp);
				if ((sp = patch_location(saix, paix, ixord, slix, apir-1)) == NULL)
					error ("Patch location out of range, strip %d patch %d",slix,apir-1);
				if (ix < 0) {
					error("Internal, got -ve patch index for generating patch %d",i);
				}
				strcpy(cols[ix].loc, sp);		/* Record location */
				cp = &cols[ix];					/* Get color for this patch */

				i++;							/* Consumed a test patch */
				if (i > npat)
					error("Internal - ran out of test patches !");
					 
				ix = rix[i];		/* Next patch index */
			}
		}

		/* Print a spacer in front of patch if requested */
		if (spacer > 0) {
			setup_spacer(&sc, pp, cp, pcol, spacer, usede);
			tro->setcolor(tro, cal, sc);
			tro->rectangle(tro, x, y, w, pspa, NULL,1);
			y -= pspa;
		}

		/* Print patch */
		{
			double wplen = plen;

			if (slix < 0) {
				wplen = tidplen;	/* TID can have a different length patch */
			}

			tro->setcolor(tro, cal, cp);	/* Patch color set above */
			if (hex) {
				int apir = pir - nmaxp;		/* Adjusted pir for max/min patches */
				tro->hexagon(tro, x, y, w, wplen, apir-1, sp);
			} else {
				/* We hack in the twin patches for the DTP20 start and stop */
				/* Initial color is set above (as for regular patches) */
				if (cpf == 1) {				/* DTP20 start bit */
					tro->rectangle(tro, x, y, w, 1.0, sp,1);	/* 50/50/50 set above */
					tro->setcolor(tro, cal, mind);		/* White */
					tro->rectangle(tro, x, y - 1.0, w, wplen - 1.0, sp,1);
				} else if (cpf == 2) {		/* DTP20 stop bit */
					tro->rectangle(tro, x, y, w, wplen - 3.0, sp,1);	/* mind set above */
					tro->setcolor(tro, cal, &pcol[8]);		/* 50/50/50 Grey for DTP20 */
					tro->rectangle(tro, x, y - wplen + 3.0, w, 3.0, sp,1);
				} else {					/* Normal patch */
					tro->rectangle(tro, x, y, w, wplen, sp,1);
				}
			}
			y -= wplen;
		}
		if (sp != NULL) {		/* Done with sp for the moment */
			free(sp);
			sp = NULL;
		}

		/* Advance the patch count */
		pir++;
		pp = cp;	/* Current color becomes last color */

		/* If this is the last patch in the row, */
		/* print a possible last spacer. */
		if (flags & IS_LPIR) {								/* End of row */
			cp = media;
			if (spacer > 0) {
				setup_spacer(&sc, pp, cp, pcol, spacer, usede);
				tro->setcolor(tro, cal, sc);
				tro->rectangle(tro, x, y, w, pspa, NULL,1);
				y -= pspa;
			}
		}

		if (flags & IS_LPIR) {								/* Last patch in row */
			pir = 1;
			ris++;

			/* First step to the middle of the patch */
			if (flags & IS_FRIS)
				x += pwid/2.0;
			else
				x += pwid/2.0 + pwex;

			/* Then step to the start of the next patch */
			if (flags & IS_LRIS) {
				if (dorspace)
					x += rrsp;				/* row to row space, making gutter */
				else
					x += pwid/2.0 + clwi;	/* no gutter, but room for cut line */
			} else
				x += pwid/2.0 + pwex;

			if ((flags & IS_FRIS) && (flags & IS_FSIP))	/* First row on page */
				x += sxwi;					/* Allow for scan compatible fatter first row */

			if (flags & IS_LRIS) {							/* End of strip */
				/* Ignore TID rows */
				if ((flags & IS_FPIF) && (flags & IS_FSIP))
					*rpsp++ = ris-tidrows-1;	/* Record how many rows in this strip */
				else
					*rpsp++ = ris-1;	/* Record how many rows in this strip */
				ris = 1;
				sip++;
				
				/* Print end of strip crop line */
				tro->setcolor(tro, cal, mark);
				if (docutmarks)			/* Generate strip cut marks */
					tro->dline(tro,x-0.3/2.0,y1,x-0.3/2.0,y2,0.3); /* 0.3 wide dotted line */

				/* Print end of strip identification if we've allowed space */
				if (dorspace)
					tro->vstring(tro,x,y1,rrsp-pwid/2.0-pwex,y2-y1,label);

				if (flags & IS_LSIP) {						/* End of page */
					sip = 1;

					x += hxew;		/* Allow space for extra bits of hexagons */

					/* Print per page label if we've allowed for it */
					if (dopglabel) {
						//tro->vstring(tro,x2,y1,pglth,y2-y1,label); /* At end of page */
						tro->vstring(tro,x+pglth,y1,pglth,y2-y1,label);	/* After last strip */
					}

					/* If we expect to scan this chart in, add some fiducial marks at the corners */
					if (scanc & 1) {
						double lw = 0.5;			/* Line width */
						double ll = 5.0;			/* Line length */

//						fprintf(of,"%% Fiducial marks\n");

						tro->setcolor(tro, cal, mark);
				
						tro->rectangle(tro, x1, y2, ll, lw, NULL, 0);			/* Top left */
						tro->rectangle(tro, x1, y2, lw, ll, NULL, 0);
						if (oft != 2)
							et_fiducial(mm2pnt(x1 + 0.5 * lw), mm2pnt(y2 - 0.5 * lw));
						else
							et_fiducial(x1 + 0.5 * lw, y2 - 0.5 * lw);
				
						tro->rectangle(tro, x2 - ll, y2, ll, lw, NULL, 0);		/* Top right */
						tro->rectangle(tro, x2 - lw, y2, lw, ll, NULL, 0);
						if (oft != 2)
							et_fiducial(mm2pnt(x2 - 0.5 * lw), mm2pnt(y2 - 0.5 * lw));
						else
							et_fiducial(x2 - 0.5 * lw, y2 - 0.5 * lw);
				
						tro->rectangle(tro, x2 - ll, y1 + lw, ll, lw, NULL, 0);	/* Bottom right */
						tro->rectangle(tro, x2 - lw, y1 + ll, lw, ll, NULL, 0);
						if (oft != 2)
							et_fiducial(mm2pnt(x2 - 0.5 * lw), mm2pnt(y1 + 0.5 * lw));
						else
							et_fiducial(x2 - 0.5 * lw, y1 + 0.5 * lw);
				
						tro->rectangle(tro, x1, y1 + lw, ll, lw, NULL, 0);		/* Bottom left */
						tro->rectangle(tro, x1, y1 + ll, lw, ll, NULL, 0);
						if (oft != 2)
							et_fiducial(mm2pnt(x1 + 0.5 * lw), mm2pnt(y1 + 0.5 * lw));
						else
							et_fiducial(x1 + 0.5 * lw, y1 + 0.5 * lw);
					}

					tro->endpage(tro);
					if (oft != 0) {		/* EPS or TIFF */
						tro->del(tro);
					}
					if (flags & IS_LPIF) {					/* Last page in file */
						if (oft == 0) {	/* PS */
							tro->del(tro);
						}
					}

					/* If we are anticipating that scanin may be used with the */
					/* chart, create the scanin recognition file for this page. */
					if (scanc & 1) {
						char chtname[MAXNAMEL+20];	/* Name of .cht file */

						if (npages > 1)
							sprintf(chtname,"%s_%02d.cht",bname,pif);
						else
							sprintf(chtname,"%s.cht",bname);
						et_write(chtname, cols, rix, l_si, i);
					}

					l_si = i;		/* New last start i */

					if (flags & IS_LPIF) {					/* Last page in file */
						break;								/* Done */
					}
					pif++;
				}
			}
		}
	}
	if (slab != NULL)
		free(slab);
	free(rix);

	*rpsp++ = 0;	/* End of rows per strip stuff */

	et_clear();		/* Cleanup edge list structures */
}

/* A paper size structure */
typedef struct {
	char *name;			/* User name (lower case) */
	double w,h;			/* Width and height in mm */
	int def;			/* Non-zero if default */
} paper;

static paper psizes[] = {
	{ "A4", 	 210.0,	297.0, 0 },
	{ "A4R", 	 297.0,	210.0, 0 },
	{ "A3", 	 297.0,	420.0, 1 },
	{ "A2", 	 420.0,	594.0, 0 },
	{ "Letter",	 215.9,	279.4, 0 },
	{ "LetterR", 279.4,	215.9, 0 },
	{ "Legal",	 215.9,	355.6, 0 },
	{ "4x6",	 101.6,	152.4, 0 },
	{ "11x17",	 279.4,	431.8, 0 },
	{ NULL, 0.0, 0.0, 0 }
};

#define DEF_MARGINE 6.0

/* Case independent string compare */
int
cistrcmp(char *s1, char *s2) {
	for (;;s1++, s2++) {
		if (tolower(*s1) != tolower(*s2))
			return 1;
		if (*s1 == '\000')
			return 0;
	}
}

#define DEF_SIXPAT	"A-Z, A-Z"				/* Default strip index pattern */		
#define DEF_PIXPAT	"0-9,@-9,@-9;1-999"		/* Default patch index pattern */		

void usage(char *diag, ...) {
	paper *pp;
	fprintf(stderr,"Generate Target PostScrip file, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: printtarg [-v] [-i instr] [-r] [-s] [-p size] basename\n");
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -i 20 | 22 | 41 | 51 | SS | i1 | p3 | CM Select target instrument (default DTP41)\n");
	fprintf(stderr,"                 i1 = i1Pro, 3p = i1Pro3+, CM = ColorMunki\n");
	fprintf(stderr,"                 20 = DTP20, 22 = DTP22, 41 = DTP41, 51 = DTP51,\n");
	fprintf(stderr,"                 SS = SpectroScan\n");
	fprintf(stderr," -h              Use hexagon patches for SS, double density for CM\n");
	fprintf(stderr," -a scale        Scale patch size and spacers by factor (e.g. 0.857 or 1.5 etc.)\n");
	fprintf(stderr," -A scale        Scale spacers by additional factor (e.g. 0.857 or 1.5 etc.)\n");
	fprintf(stderr," -r              Don't randomize patch location\n");
	fprintf(stderr," -s              Create a scan image recognition (.cht) file\n");
	fprintf(stderr," -S              Same as -s, but don't generate wide orientation strip.\n");
	fprintf(stderr," -c              Force colored spacers\n");
	fprintf(stderr," -b              Force B&W spacers\n");
	fprintf(stderr," -n              Force no spacers\n");
	fprintf(stderr," -f              Create PostScript DeviceN Color fallback\n");
	fprintf(stderr," -w g|r|s|n      White colorspace encoding: DeviceGray (def), DeviceRGB, Separation or DeviceN\n");            
	fprintf(stderr," -k g|c|s|n      Black colorspace encoding: DeviceGray (def), DeviceCMYK, Separation or DeviceN\n");            
	fprintf(stderr," -o k|r|n        CMY colorspace encoding DeviceCMYK (def), inverted DeviceRGB or DeviceN\n");            
	fprintf(stderr," -e              Output EPS compatible file\n");
	fprintf(stderr," -t [res]        Output 8 bit TIFF raster file, optional res DPI (default 100)\n");
	fprintf(stderr," -T [res]        Output 16 bit TIFF raster file, optional res DPI (default 100)\n");
	fprintf(stderr," -C              Don't use TIFF compression\n");
	fprintf(stderr," -N              Use TIFF alpha N channels more than 4\n");
	fprintf(stderr," -D              Dither 8 bit TIFF values down from 16 bit\n");
	fprintf(stderr," -Q nbits        Quantize test values to fit in nbits\n");
	fprintf(stderr," -R rsnum        Use given random start number\n");
	fprintf(stderr," -K file.cal     Apply printer calibration to patch values and include in .ti2\n");
	fprintf(stderr," -I file.cal     Include calibration in .ti2 (but don't apply it)\n");
	fprintf(stderr," -x pattern      Use given strip indexing pattern (Default = \"%s\")\n",DEF_SIXPAT);
	fprintf(stderr," -y pattern      Use given patch indexing pattern (Default = \"%s\")\n",DEF_PIXPAT);
	fprintf(stderr," -m margin       Set a page margin in mm (default %3.1f mm)\n",DEF_MARGINE);       
	fprintf(stderr," -M margin       Set a page margin in mm and include it in TIFF\n");       
	fprintf(stderr," -P              Don't limit strip length\n");       
	fprintf(stderr," -L              Suppress any left paper clip border\n");       
	fprintf(stderr," -U              Suppress CUPS cupsJobTicket: cups-disable-cmm in PS & EPS files\n");
	fprintf(stderr," -p size         Select page size from:\n");
	for (pp = psizes; pp->name != NULL; pp++)
		fprintf(stderr,"                 %-8s [%.1f x %.1f mm]%s\n", pp->name, pp->w, pp->h,
		pp->def ? " (default)" : "");
	fprintf(stderr," -p WWWxHHH      Custom size, WWW mm wide by HHH mm high\n");
	fprintf(stderr," basname         Base name for input(.ti1), output(.ti2) and output(.ps/.eps/.tif)\n");
	exit(1);
	}

int
main(argc,argv)
int argc;
char *argv[];
{
	int fa, nfa, mfa;		/* argument we're looking at */
	int verb = 0;
	int hflag = 0;			/* Hexagon patches for SS, high density for CM */
	double pscale = 1.0;	/* Patch size scale */
	double sscale = 1.0;	/* Spacer size scale */
	int rand = 1;
	int qbits = 0;			/* Quantization bits */
	int oft = 0;			/* Ouput File type, 0 = PS, 1 = EPS , 2 = TIFF */
	int nocups = 0;			/* Supress CUPS PS/EPS job ticket */
	depth2d tiffdpth = bpc8_2d;	/* TIFF pixel depth */
	double tiffres = 100.0;	/* TIFF resolution in DPI */
	int ncha = 0;			/* flag, use nchannel alpha */
	int tiffdith = 0;		/* flag, use TIFF 8 bit dithering */
	int tiffcomp = 1;		/* flag, use TIFF compression */
	int spacer = -1;		/* -1 = default for instrument */
							/* 0 = forse no spacer, 1 = Force B&W spacers */
							/* 2 = Force colored spacer */
	int rstart = -1;		/* Random sequence start value */
	char *sixpat = DEF_SIXPAT;	/* Strip index pattern */		
	char *pixpat = DEF_PIXPAT;	/* Patch index pattern */		
	alphix *saix, *paix;	/* Strip and Patch index generators */
	int ixord = 0;			/* Index order, 0 = strip then patch */
	int scanc = 0;			/* Scan compatible bits, 1 = .cht, 2 = wide first row */
	int devnfb = 0;			/* Add device N fallback colors */
	int altrep = 0;			/* Device K/W/CMY color type 0..8 */
	int applycal = 0;		/* NZ to apply calibration */
	static char inname[MAXNAMEL+20] = { 0 };	/* Input cgats file name */
	static char calname[MAXNAMEL+1] = { 0 };	/* Input printer calibration */
	static char psname[MAXNAMEL+1] = { 0 };		/* Output postscrip file base name */
	static char outname[MAXNAMEL+20] = { 0 };	/* Output cgats file name */
	cgats *icg;				/* input cgats structure */
	cgats *ocg;				/* output cgats structure */
	xcal *cal = NULL;		/* printer calibration */ 
	instType itype = instDTP41;		/* Default target instrument */
	int itype_mod = 0;		/* Instrument type modifier. 1 = i1Pro3+ */
	int nmask = 0;			/* Device colorant mask */
	int nchan = 0;			/* Number of device chanels */
	int i;
	int si, fi, wi;			/* sample id index, field index, keyWord index */
	char label[400];		/* Space for chart label */
	double marg = DEF_MARGINE;	/* Margin from paper edge in mm */
	int nosubmarg = 0;		/* Don't subtract it from raster */
	int nolpcbord = 0;		/* NZ to suppress left paper clip border */
	int nollimit = 0;		/* NZ to release any strip length limits */
	paper *pap = NULL;		/* Paper size pointer, NULL if custom */
	double cwidth, cheight;	/* Custom paper width and height in mm */
	col *cols;				/* test patch colors */
	int npat;				/* Number of patches */
	int nppat;				/* Number of patches including padding */
	col pcold[8];			/* pre-defined density extreme colors */
	int pcolvv = 0;			/* pcolv valid if nz */
	col pcolv[9];			/* pre-defined device combination ecolors */
	col *pcol;				/* Chosen color patches for device */
	double wp[3];			/* Approximate XYZ white point */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	int sip;				/* Steps in Pass */
	unsigned char *pis;		/* Passes in strip array */
	double plen, glen, tlen;/* Patch, gap and trailer length in mm */
	char *bp, buf[500];		/* general sprintf buffer */

	error_program = "printtarg";
	check_if_not_interactive();

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	if (argc <= 1)
		usage("Not enough arguments");

#ifdef DEBUG
# pragma message("######### printtarg DEBUG is #define ########")
	fprintf(stderr,"target: DEBUG is #defined\n");
#endif

	/* Find the default paper size */
	for (pap = psizes; pap->name != NULL; pap++) {
		if (pap->def != 0)
			break;
	}
	if (pap->name == NULL)
		error ("Internal - can't find default paper size");

	/* Process the arguments */
	mfa = 1;						/* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')	{	/* Look for any flags */
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

			if (argv[fa][1] == '?')
				usage("Requested usage");

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;

			/* hflag patches */
			else if (argv[fa][1] == 'h' || argv[fa][1] == 'H')
				hflag = 1;

			/* Patch scaling */
			else if (argv[fa][1] == 'a') {
				fa = nfa;
				if (na == NULL) usage("Expected scale factor to -a");
				pscale = atof(na);
				if (pscale < 0.1 || pscale > 4.0)
					usage("Scale factor %f is outside expected range 0.1 - 4.0",pscale);
			}

			/* Spacer scaling */
			else if (argv[fa][1] == 'A') {
				fa = nfa;
				if (na == NULL) usage("Expected scale factor to -a");
				sscale = atof(na);
				if (sscale < 0.1 || sscale > 8.0)
					usage("Scale factor %f is outside expected range 0.1 - 8.0",sscale);
			}

			/* Scan compatible */
			else if (argv[fa][1] == 's')
				scanc = 3;

			else if (argv[fa][1] == 'S')
				scanc = 1;

			/* Force colored spacer */
			else if (argv[fa][1] == 'c')
				spacer = 2;

			/* Force B&W spacer */
			else if (argv[fa][1] == 'b' || argv[fa][1] == 'B')
				spacer = 1;

			/* No spacer */
			else if (argv[fa][1] == 'n')
				spacer = 0;

			/* Randomisation off */
			else if (argv[fa][1] == 'r')
				rand = 0;

			/* Specify random seed */
			else if (argv[fa][1] == 'R') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -R");
				rstart = atoi(na);
				if (rstart < 0)
					usage("Argument to -R must be positive");
			}

			/* Enable DeviceN color fallback */
			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F')
				devnfb = 1;

			/* Select the printer W color representation */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -w");
				switch(na[0]) {
					case 'g':
					case 'G':
						altrep = 0;
						break;
					case 'r':
					case 'R':
						altrep = 4;
						break;
					case 's':
					case 'S':
						altrep = 5;
						break;
					case 'n':
					case 'N':
						altrep = 6;
						break;
					default:
						usage("Unexpected argument to -w");
				}
			}

			/* Select the printer K color representation */
			else if (argv[fa][1] == 'k') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -k");
				switch(na[0]) {
					case 'g':
					case 'G':
						altrep = 0;
						break;
					case 'c':
					case 'C':
						altrep = 1;
						break;
					case 's':
					case 'S':
						altrep = 2;
						break;
					case 'n':
					case 'N':
						altrep = 3;
						break;
					default:
						usage("Unexpected argument to -k");
				}
			}

			/* Select the printer CMY color representation */
			else if (argv[fa][1] == 'o') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -o");
				switch(na[0]) {
					case 'k':
					case 'K':
						altrep = 0;
						break;
					case 'r':
					case 'R':
						altrep = 7;
						break;
					case 'n':
					case 'N':
						altrep = 8;
						break;
					default:
						usage("Unexpected argument to -o");
				}
			}

			/* EPS */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E')
				oft = 1;

			/* TIFF */
			else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				oft = 2;
				if (argv[fa][1] == 'T')
					tiffdpth = bpc16_2d;
				else
					tiffdpth = bpc8_2d;

				if (na != NULL) {	/* Found an optional resolution */
					fa = nfa;
					tiffres = atof(na);
					if (tiffres <= 1.0 || tiffres > 1e6)
						usage("TIFF resolution is out of range");
				}
			}
			/* use Nchannel alpha for TIFF  */
			else if (argv[fa][1] == 'N') {
				ncha = 1;
			}
			/* use 16->8 bit dithering for 8 bit TIFF  */
			else if (argv[fa][1] == 'D') {
				tiffdith = 1;
			}
			/* Don't use TIFF compression */
			else if (argv[fa][1] == 'C') {
				tiffcomp = 0;
			}

			/* Specify quantization bits */
			else if (argv[fa][1] == 'Q') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -Q");
				qbits = atoi(na);
				if (qbits < 1 || qbits > 32)
					usage("Argument to -Q must be between 1 and 32");
			}

			/* Specify strip index pattern */
			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -x");
				sixpat = na;
			}

			/* Specify patch index pattern */
			else if (argv[fa][1] == 'y' || argv[fa][1] == 'Y') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -y");
				pixpat = na;
			}

			/* Border margin */
			else if (argv[fa][1] == 'm' || argv[fa][1] == 'M') {
				if (argv[fa][1] == 'M')
					nosubmarg = 1;
				fa = nfa;
				if (na == NULL) usage("Expected border margine argument to -m");
				marg = atof(na);
				if (marg < 0.0 || marg > 50.0)
					usage("Border margin %f is outside expected range",marg);
			}

			/* Don't limit the strip length */
			else if (argv[fa][1] == 'P') {
				nollimit = 1;
			}

			/* Suppress left paper clip border */
			else if (argv[fa][1] == 'L') {
				nolpcbord = 1;
			}

			/* Suppress CUPS job ticket */
			else if (argv[fa][1] == 'U') {
				nocups = 1;
			}

			/* Page size */
			else if (argv[fa][1] == 'p') {
				fa = nfa;
				if (na == NULL) usage("Expected an argument to -p");
				for (pap = psizes; pap->name != NULL; pap++) {
					if (cistrcmp(na, pap->name) == 0)
						break;
				}
				
				if (pap->name == NULL) {	/* See if it matches a custom size */
					if (sscanf(na,"%lfx%lf",&cwidth, &cheight) == 2) {
						pap = NULL;			/* Indicate custom */
						if (cwidth < 1.0 || cwidth > 4000.0
						 || cheight < 1.0 || cheight > 4000.0)
							usage("Argument to -p was of unexpected size");		/* Sanity check */
					} else {
						usage("Failed to recognise argument to -p");
					}
				}
			}
			/* Target Instrument type */
			else if (argv[fa][1] == 'i') {
				fa = nfa;
				if (na == NULL) usage("Expected an argument to -i");

				if (strcmp("20", na) == 0)
					itype = instDTP20;
				else if (strcmp("22", na) == 0)
					itype = instDTP22;
				else if (strcmp("41", na) == 0)
					itype = instDTP41;
				else if (strcmp("51", na) == 0)
					itype = instDTP51;
				else if (strcmp("SS", na) == 0 || strcmp("ss", na) == 0)
					itype = instSpectroScan;
				else if (strcmp("i1", na) == 0)
					itype = instI1Pro;
				else if (strcmp("3p", na) == 0) {
					itype = instI1Pro;
					itype_mod = 1;
				} else if (strcmp("cm", na) == 0 || strcmp("CM", na) == 0)
					itype = instColorMunki;
				else
					usage("Argument to -i wasn't recognised");

			/* Printer calibration */
			} else if (argv[fa][1] == 'K' || argv[fa][1] == 'I') {
				if (argv[fa][1] == 'K')
					applycal = 1;
				else
					applycal = 0;
				fa = nfa;
				if (na == NULL) usage("Expected an argument to -%c",argv[fa][1]);
				strncpy(calname,na,MAXNAMEL); calname[MAXNAMEL] = '\000';
			} else 
				usage("Unknown flag");
		}
		else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage("Expecting basename argument");
	strncpy(inname,argv[fa],MAXNAMEL); inname[MAXNAMEL] = '\000';
	strcpy(outname,inname);
	strcpy(psname,inname);
	strcat(inname,".ti1");
	strcat(outname,".ti2");

	if (calname[0] != '\000') {
		if ((cal = new_xcal()) == NULL)
			error("new_xcal failed");
		if ((cal->read(cal, calname)) != 0)
			error("%s",cal->err);
	}

	/* Set default qantization for known output */
	if (qbits == 0 && oft == 2) { 
		if (tiffdpth == bpc16_2d || tiffdith != 0)
			qbits = 16;
		else if (tiffdpth == bpc8_2d)
			qbits = 8;
	}

	if (itype == instSpectroScan) {
		if (scanc) {
			if (verb)
				printf("Can only select hexagonal patches if no scan recognition is needed - ignored!\n");
			hflag = 0;
		}
	} else if (itype == instColorMunki) {
		/* OK */
	} else if (hflag) {
		if (verb)
			printf("Can only select h flag for SpectrScan or ColorMunki - ignored!\n");
		hflag = 0;
	}

	if ((saix = new_alphix(sixpat)) == NULL)
		error("Strip indexing pattern '%s' doesn't parse",sixpat);

	if ((paix = new_alphix(pixpat)) == NULL)
		error("Patch in strip indexing pattern '%s' doesn't parse",pixpat);

	icg = new_cgats();	/* Create a CGATS structure */
	icg->add_other(icg, "CTI1"); 	/* our special input type is Calibration Target Information 1 */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI1 format file");
	if (icg->ntables < 2 || icg->ntables > 3)
		error ("Input file doesn't contain two or three tables");

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	/* Allocate room for test patches and maximum padding patches */
	if ((cols = (col *)malloc(sizeof(col) * (npat + MAXPPROW))) == NULL)
		error("Malloc failed!");

	/* Setup output cgats file */
	ocg = new_cgats();	/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI2"); 	/* our special type is Calibration Target Information 2 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 2",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll printtarg", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);

	/* Note what instrument the chart is setup for */
	ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", inst_name(itype) , NULL);

	/* Copy various parameters through */
	if ((wi = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[wi], NULL);
	
	if ((wi = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[wi], NULL);
	
	if ((wi = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[wi], NULL);

	if ((wi = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
		ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[wi], NULL);

	if ((wi = icg->find_kword(icg, 0, "ACCURATE_EXPECTED_VALUES")) >= 0)
		ocg->add_kword(ocg, 0, "ACCURATE_EXPECTED_VALUES",icg->t[0].kdata[wi], NULL);

	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
	ocg->add_field(ocg, 0, "SAMPLE_LOC", cs_t);

	if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
		error ("Input file '%s' doesn't contain field SAMPLE_ID in first table",inname);
	if (icg->t[0].ftype[si] != nqcs_t)
		error ("Field SAMPLE_ID is wrong type");

	/* Read the approximate white point */
	if ((fi = icg->find_kword(icg, 0, "APPROX_WHITE_POINT")) < 0)
		error ("Input file doesn't contain keyword APPROX_WHITE_POINT");
	if (sscanf(icg->t[0].kdata[fi], "%lf %lf %lf", &wp[0], &wp[1], &wp[2]) != 3)
		error ("Couldn't parse the white point data correctly");
	wp[0] /= 100.0; wp[1] /= 100.0; wp[2] /= 100.0;
	ocg->add_kword(ocg, 0, "APPROX_WHITE_POINT",icg->t[0].kdata[fi], NULL);

//printf("~1 got approx white point of %f %f %f\n",wp[0],wp[1],wp[2]);

	/* Figure out the color space */
	if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file '%s' doesn't contain keyword COLOR_REPS",inname);

	if ((nmask = icx_char2inkmask(icg->t[0].kdata[fi])) != 0) {
		int i, j, ii;
		int chix[ICX_MXINKS];	/* Device chanel indexes */
		int xyzix[3];			/* XYZ chanel indexes */
		char *ident;			/* Full ident */
		char *bident;			/* Base ident */
		char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };
		double qscale = (1 << qbits) - 1.0;

		if (cal != NULL && nmask != cal->devmask)
			error ("Calibration colorspace %s doesn't match .ti1 %s",icx_inkmask2char(cal->devmask, 1),icx_inkmask2char(nmask, 1));

		if ((ii = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0)
			ocg->add_kword(ocg, 0, "TOTAL_INK_LIMIT",icg->t[0].kdata[ii], NULL);

		nchan = icx_noofinks(nmask);
		ident = icx_inkmask2char(nmask, 1); 
		bident = icx_inkmask2char(nmask, 0); 

		for (j = 0; j < nchan; j++) {
			int imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : bident,
			                      icx_ink2char(imask));

			if ((ii = icg->find_field(icg, 0, fname)) < 0)
				error ("Input file '%s' doesn't contain field %s in first table",inname,fname);
			if (icg->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type",fname);
	
			ocg->add_field(ocg, 0, fname, r_t);
			chix[j] = ii;
		}

		for (j = 0; j < 3; j++) {
			if ((ii = icg->find_field(icg, 0, xyzfname[j])) < 0)
				error ("Input '%s' file doesn't contain field %s in first table",inname,xyzfname[j]);
			if (icg->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type",xyzfname[j]);
	
			ocg->add_field(ocg, 0, xyzfname[j], r_t);
			xyzix[j] = ii;
		}

		ocg->add_kword(ocg, 0, "COLOR_REP", ident, NULL);

		/* Read all the test patches in, and quantize them */
		for (i = 0; i < npat; i++) {
			cols[i].i = i;
			cols[i].t = T_N | T_XYZ;
			if (devnfb)
				cols[i].t |= T_NFB;
			cols[i].nmask = nmask;
			cols[i].altrep = altrep;
			cols[i].n  = nchan;
			cols[i].id = ((char *)icg->t[0].fdata[i][si]);
			sprintf(cols[i].loc, "???");
			for (j = 0; j < nchan; j++) {
				double vr, vv = *((double *)icg->t[0].fdata[i][chix[j]]) / 100.0;
				if (qbits > 0) {
					vv *= qscale;
					vr = floor(vv + 0.5);
					if ((vr - vv) == 0.5 && (((int)vr) & 1) != 0) /* Round to even */
						vr -= 1.0;
					vv = vr/qscale;
				}
				cols[i].dev[j] = vv;
			}
			for (j = 0; j < 3; j++)
				cols[i].XYZ[j] = *((double *)icg->t[0].fdata[i][xyzix[j]]) / 100.0;
			col_convert(&cols[i], wp);	/* Ensure other representations */
		}

		free(ident);
		free(bident);
	} else
		error ("Input file keyword COLOR_REPS has unknown value");

	/* Load up the pre-defined density extreme spacer colors */
	{
		int i, j, ii;
		int nsp;
		int chix[ICX_MXINKS];	/* Device chanel indexes */
		int xyzix[3];			/* XYZ chanel indexes */
		char *bident;
		char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };

		nchan = icx_noofinks(nmask);
		bident = icx_inkmask2char(nmask, 0); 

		if ((nsp = icg->t[1].nsets) <= 0)
			error ("No sets of data in second table");

		for (j = 0; j < nchan; j++) {
			int imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : bident,
			                      icx_ink2char(imask));

			if ((ii = icg->find_field(icg, 1, fname)) < 0)
				error ("Input file '%s' doesn't contain field %s in second table",inname,fname);
			if (icg->t[1].ftype[ii] != r_t)
				error ("Field %s is wrong type",fname);
			chix[j] = ii;
		}

		for (j = 0; j < 3; j++) {
			if ((ii = icg->find_field(icg, 1, xyzfname[j])) < 0)
				error ("Input file '%s' doesn't contain field %s in second table",inname,xyzfname[j]);
			if (icg->t[1].ftype[ii] != r_t)
				error ("Field %s is wrong type",xyzfname[j]);
			xyzix[j] = ii;
		}

		if (nsp != 8)
			error ("Expect second set of data to have 8 sets, found %d",nsp);

		/* Read all the density spacer patches in */
		for (i = 0; i < nsp; i++) {
			pcold[i].i = -1;
			pcold[i].t = T_N | T_XYZ | T_PRESET;
			if (devnfb)
				pcold[i].t |= T_NFB;
			pcold[i].nmask = nmask;
			pcold[i].altrep = altrep;
			pcold[i].n  = nchan;
			pcold[i].id = "";
			sprintf(cols[i].loc, "???");
			for (j = 0; j < nchan; j++)
				pcold[i].dev[j] = *((double *)icg->t[1].fdata[i][chix[j]]) / 100.0;
			for (j = 0; j < 3; j++)
				pcold[i].XYZ[j] = *((double *)icg->t[1].fdata[i][xyzix[j]]) / 100.0;
			col_convert(&pcold[i], wp);	/* Ensure other representations */
		}

		free(bident);
	}

	/* Load up the pre-defined device combination barcode colors */
	if (icg->ntables >= 3) {
		int i, j, ii;
		int nsp;
		int chix[ICX_MXINKS];	/* Device chanel indexes */
		int xyzix[3];			/* XYZ chanel indexes */
		char *bident;			/* Base ident */
		char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };

		nchan = icx_noofinks(nmask);
		bident = icx_inkmask2char(nmask, 0); 

		if ((nsp = icg->t[2].nsets) > 0) {

			for (j = 0; j < nchan; j++) {
				int imask;
				char fname[100];

				imask = icx_index2ink(nmask, j);
				sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : bident,
				                      icx_ink2char(imask));

				if ((ii = icg->find_field(icg, 2, fname)) < 0)
					error ("Input file '%s' doesn't contain field %s in third table",inname,fname);
				if (icg->t[2].ftype[ii] != r_t)
					error ("Field %s is wrong type",fname);
				chix[j] = ii;
			}

			for (j = 0; j < 3; j++) {
				if ((ii = icg->find_field(icg, 2, xyzfname[j])) < 0)
					error ("Input file '%s' doesn't contain field %s in third table",inname,xyzfname[j]);
				if (icg->t[2].ftype[ii] != r_t)
					error ("Field %s is wrong type",xyzfname[j]);
				xyzix[j] = ii;
			}

			if (nsp != 9)
				error ("Expect third set of data to have 9 sets, found %d",nsp);

			/* Read all the barcode CMY color patches in */
			for (i = 0; i < nsp; i++) {
				pcolv[i].i = -1;
				pcolv[i].t = T_N | T_XYZ | T_PRESET;
				if (devnfb)
					pcolv[i].t |= T_NFB;
				pcolv[i].nmask = nmask;
				pcolv[i].altrep = altrep;
				pcolv[i].n  = nchan;
				pcolv[i].id = "";
				sprintf(cols[i].loc, "???");
				for (j = 0; j < nchan; j++)
					pcolv[i].dev[j] = *((double *)icg->t[2].fdata[i][chix[j]]) / 100.0;
				for (j = 0; j < 3; j++)
					pcolv[i].XYZ[j] = *((double *)icg->t[2].fdata[i][xyzix[j]]) / 100.0;
				col_convert(&pcolv[i], wp);	/* Ensure other representations */
			}

			free(bident);
			pcolvv = 1;
		}
	}

	if (pap != NULL) {
		sprintf(buf, "%.1fx%.1f",pap->w, pap->h);
		if (verb)
			printf("Paper chosen is %s	[%.1f x %.1f mm]\n", pap->name, pap->w, pap->h);
	} else {
		sprintf(buf, "%.1fx%.1f",cwidth, cheight);
		if (verb)
			printf("Paper chosen is custom %.1f x %.1f mm\n", cwidth, cheight);
	}
	ocg->add_kword(ocg, 0, "PAPER_SIZE", buf, NULL);

	if (rstart == -1) {
		rstart = clk % npat;
	} else {
		rstart = rstart % npat;
	}
	sprintf(buf,"%d",rstart);
	if (rand)
		ocg->add_kword(ocg, 0, "RANDOM_START", buf, NULL);
	else
		ocg->add_kword(ocg, 0, "CHART_ID", buf, NULL);

	if (itype == instSpectroScan && hflag) {
		ocg->add_kword(ocg, 0, "HEXAGON_PATCHES", "True", NULL);
	}

	if (itype == instDTP20) {
		if (pcolvv == 0)
			error("Input file '%s' doesn't contain device combination table needed for DTP20",inname);
		pcol = pcolv;		/* Barcode color values */
	} else
		pcol = pcold;		/* Density spacer alues */

	
	sprintf(label, "ArgyllCMS - Chart \"%s\" (%s %d) %s",
	               psname, rand ? "Random Start" : "Chart ID", rstart, atm);
	generate_file(itype, itype_mod, psname, cols, npat, applycal ? cal : NULL, label,
	            pap != NULL ? pap->w : cwidth, pap != NULL ? pap->h : cheight,
	            marg, nosubmarg, nollimit, nolpcbord, rand, rstart, saix, paix,	ixord,
	            pscale, sscale, hflag, verb, scanc, oft, nocups, tiffdpth, tiffres, ncha, tiffdith,
	            tiffcomp, spacer, nmask, altrep, pcol, wp,
	            &sip, &pis, &plen, &glen, &tlen, &nppat);

	if (itype == instDTP20
	 || itype == instDTP41) {	/* DTP20/41 needs this */
		sprintf(buf,"%f",plen);
		ocg->add_kword(ocg, 0, "PATCH_LENGTH", buf, NULL);
		sprintf(buf,"%f",glen);
		ocg->add_kword(ocg, 0, "GAP_LENGTH", buf, NULL);
		if (itype == instDTP41) {	/* DTP41 needs this */
			sprintf(buf,"%f",tlen);
			ocg->add_kword(ocg, 0, "TRAILER_LENGTH", buf, NULL);
		}
	}

	sprintf(buf,"%d",sip);
	ocg->add_kword(ocg, 0, "STEPS_IN_PASS", buf, NULL);

	/* Convert pass in strips count to base 62 */
	buf[0] = '\000';
	bp = buf;
	for (i = 0; ;i++) {
		if (pis[i] == 0)
			break;
		sprintf(bp, "%s%d", i > 0 ? "," : "", pis[i]);
		bp += strlen(bp);
	}
	ocg->add_kword(ocg, 0, "PASSES_IN_STRIPS2", buf, NULL);

	/* Output the default Argyll style strip and patch numbering */
	ocg->add_kword(ocg, 0, "STRIP_INDEX_PATTERN", sixpat, NULL);
	ocg->add_kword(ocg, 0, "PATCH_INDEX_PATTERN", pixpat, NULL);
	ocg->add_kword(ocg, 0, "INDEX_ORDER", ixord ? "PATCH_THEN_STRIP" : "STRIP_THEN_PATCH", NULL);

	/* Write out the patch info to the output CGATS file */
	for (i = 0; i < nppat; i++) {
		cgats_set_elem ary[2 + ICX_MXINKS + 3];
		int j;

		if (strcmp(cols[i].loc, "???") == 0)
			warning ("Internal, patch %s (%d) wasn't given a valid location string",cols[i].id,i+1);
		ary[0].c = cols[i].id;
		ary[1].c = cols[i].loc;
		for (j = 0; j < nchan; j++)
			ary[2 + j].d = 100.0 * cols[i].dev[j];
		for (j = 0; j < 3; j++)
			ary[2 + nchan + j].d = 100.0 * cols[i].XYZ[j];
		ocg->add_setarr(ocg, 0, ary);
	}

	/* If there is a calibration, append it to the .ti2 file */
	if (cal != NULL) {
		if (cal->write_cgats(cal, ocg) != 0)
			error("Writing cal error : %s",cal->err);
	}

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	if (cal != NULL)
		cal->del(cal);
	paix->del(paix);
	saix->del(saix);
	free(pis);
	free(cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */

	return 0;
}

/******************************************************************/
/* Edge tracking support, for generating the scanner image        */
/* recognition reference chart file. */

/* Establish width and height to convert between topleft and */
/* bottom left origin ??? ~~9 */

/*
	Basic algorithm strategy:

	First we simply accumulate the raw recognition and patch
	identification information. Once the chart is generated, we:
		sort into horizontal and vertical half edges
		sort into +ve and -ve edges
		match +ve and -ve edges
		for each match, generate a delta edge segment
		Assume any non-matched edges are against the media.
		Coalesce delta edges into X & Y edge lists
		Compute normalised strength.
		Compute crossing count.
		Figure average box size, and compute shrink.
		
*/

/* A half edge structure */
/* coordinate origin is top left */
struct _hedge {
	int ix;			/* Index for debug id */
	double rgb[3];	/* Color this half edge transitions to */
	int negh;		/* 1 if this is a -ve major coordinate side half edge, 0 otherwise */
	double mj;		/* Major coordinate offset (ie. X coord for vertical edge) */
	double mi0;		/* Minor coordinate smaller value (ie. Y for vertical edge) */
	double mi1;		/* Minor coordinate larger value (ie. Y for vertical edge) */
	struct _hedge *next;	/* Next in linked list */
}; typedef struct _hedge hedge;

/* A patch identifier */
/* coordinate origin is top left */
struct _patch {
	char id[20];	/* ID string, Zero length if a diagnostic rectangle */
	double xo;		/* Location of the rectangle origin (bottom left ???) */
	double yo;
	double w;		/* Size of the patch */
	double h;
	struct _patch *next;	/* Next in linked list */
}; typedef struct _patch patch;


/* Structure to one edge */
struct _edge {
	double p1,p2;	/* Start and end of line in orthogonal direction */
	struct _edge *next;	/* next in the linked list */
}; typedef struct _edge edge;

/* Structure of an edge list */
struct _elist {
	double pos;		/* Position of edges along major axis */
	double len;		/* Total length of edges atthis position */
	double cc;		/* Crossing count */
	int ne;			/* Count of edges */
	edge *e;		/* Head of linked list of edges at this position */
	struct _elist *next;	/* Next in linked list */
}; typedef struct _elist elist;

/* - - - - - - - - - - - - - - - - - - - */
/* Structure to track recognition edges */
struct {
	double height;	/* Height of the page */
	double mrgb[3];	/* Media RGB */	
	double rgb[3];	/* Currently set RGB */	

	int nfid;		/* Number of fiducial marks. Must be 4 to cause output */
	double fx[4];	/* Fiducial mark X coordinates */
	double fy[4];	/* Fiducial mark Y coordinates */

	/* Raw half edge lists, [vertical, horizontal] */
	int nhe[2];
	hedge *he[2];	/* Pointer to start of hedhe linked list */

	/* Patch identity information */
	int npatches;
	patch *patches;

	/* Processed information */
	hedge **she[2];	/* Sorted half edges */

	int nel[2];		/* Number of edge positions */
	elist *el[2];	/* Head of edge linked list */
	elist **nelp;	/* Next edge to append to */
} et;

/* Initialise the structure */
void et_init(void) {
	memset(&et, 0, sizeof(et));
}

/* Tell et of the height, so the Y coordinate can be flipped */
void et_height(double height) {
	et.height = height;
//printf("~1 media height set to %f\n",height);
}

/* Tell et of the media color */
void et_media(double *rgb) {
	int e;
	for (e = 0; e < 3; e++)
		et.mrgb[e] = rgb[e];
}

/* Track the current GC color */
void et_color(
double *rgb		/* New RGB values */
) {
	int e;
	for (e = 0; e < 3; e++)
		et.rgb[e] = rgb[e];
}

/* Track a drawn object half edge */
/* We assume that no object is written over any other object, */
/* and that each half edge has a perfect opposite edge (ie. same */
/* mi0 and mi1), or no matching half edge (it is over the media) - */
/* ie. no partialy overlapping half edges. */
/* The arguments origin is assumed to be bottom left */
void et_edge(
int isx,		/* NZ if this is a vertical edge */
int negh,		/* NZ if this is a -ve major coordinate side half edge */
double mj,		/* Major coordinate offset (ie. X coord for vertical edge) */
double mi0,		/* Minor coordinate smaller value (ie. Y for vertical edge) */
double mi1		/* Minor coordinate larger value (ie. Y for vertical edge) */
) {
	int e, h;
	hedge *he;

	if (mi1 < mi0)
		error ("et_edge, minor coords wern't ordered");

	if ((he = (hedge *)calloc(sizeof(hedge), 1)) == NULL)
		error("Malloc of half edge structure failed");

	for (e = 0; e < 3; e++)
		he->rgb[e] = et.rgb[e];

	/* Flip the Y coordinate */
	if (isx) {
		double tmi0, tmi1;
		tmi0 = et.height - mi1;		/* swap to keep smallest small */
		tmi1 = et.height - mi0;
		mi0 = tmi0;
		mi1 = tmi1;
	} else {
		mj = et.height - mj;
	}

	he->negh = negh ? 1 : 0;
	he->mj   = mj;
	he->mi0  = mi0;
	he->mi1  = mi1;

	/* Add half edges to the list */
	h = isx ? 0 : 1;
	et.nhe[h]++;
	he->next = et.he[h];
	et.he[h] = he;
}

/* Track a patch identity */
/* The arguments origin is assumed to be bottom left */
void et_patch(
char  *id,		/* ID string, NULL if a diagnostic mark */
double xo,		/* Bottom left of the rectangle */
double yo,
double w,		/* Size of the patch */
double h
) {
	patch *p;

//printf("~1 got patch at %f %f, w %f, h %f\n", xo, yo, w, h);

	if ((p = (patch *)calloc(sizeof(patch), 1)) == NULL)
		error("Malloc of patch structure failed");

	/* Flip Y */
	yo = et.height - (yo + h);

	if (id != NULL) {
		strncpy(p->id, id, 19);
		p->id[19] = '\000';
	} else {
		p->id[0] = '\000';
	}
	p->xo = xo;
	p->yo = yo;
	p->h  = h;
	p->w  = w;

	/* Add patch to list */
	et.npatches++;
	p->next = et.patches;
	et.patches = p;
}

/* Add a fiducial mark location */
/* It is an error to add more than 4, */
/* and exactly 4 have to be added to cause fiducials to be output. */
void et_fiducial(
double x,		/* Bottom left of the rectangle */
double y
) {
	if (et.nfid >= 4)
		error("et_fiducial: too many fiducial marks");
	et.fx[et.nfid] = x;
	et.fy[et.nfid] = et.height - y;		/* Flip Y */
	et.nfid++;
}

/* Compute the image recognition information, and write the */
/* .cht file. */
void et_write(char *fname, col *cols, int *rix, int si, int ei) {
	FILE *of;
	hedge *ep0, *ep1, *epe;
	int i, h;

//printf("~1 et has %d vertical and %d horizontal half edges\n", et.nhe[0], et.nhe[1]);
//printf("~1 et has %d patches\n", et.npatches);

	/* Do X then Y */
	for (h = 0; h < 2; h++) {

		/* Create sorted list of vertical half edges */
		if ((et.she[h] = (hedge **)malloc(sizeof(patch*) * et.nhe[h])) == NULL)
			error("Malloc of array of vertical halfedge pointers failed");
	
		for (ep0 = et.he[h], i = 0; i < et.nhe[h]; i++, ep0 = ep0->next)
			et.she[h][i] = ep0;
	
		/* Sort helf edges by their X location, then their Y0 location */
#define HEAP_COMPARE(A,B) (fabs(A->mj - B->mj) < 1e-6 ? A->mi0 < B->mi0 : A->mj < B->mj)
		HEAPSORT(hedge *, et.she[h], et.nhe[h]);
#undef HEAP_COMPARE

		/* Re-create the linked list in sorted order */
		et.he[h] = NULL;
		for (i = et.nhe[h]-1; i >= 0; i--) {
			et.she[h][i]->next = et.he[h];
			et.he[h] = et.she[h][i];
			et.he[h]->ix = i;
		}
		
		free(et.she[h]);
		et.she[h] = NULL;

#ifdef DEBUG
		fprintf(stderr,"Sorted %s half edges:\n",h ? "Vertical" : "Horizontal");
		for (ep0 = et.he[h]; ep0 != NULL; ep0 = ep0->next) {
			fprintf(stderr,"%s %d at %c = %f from %c = %f to %f\n",
			h == 0 ? "Vert" : "Horiz", ep0->ix,
			h == 0 ? 'X' : 'Y', ep0->mj,
			h == 0 ? 'Y' : 'X', ep0->mi0, ep0->mi1);
		}
#endif /* DEBUG */

		/* Do a first pass to locate and split any part overlapping pairs of edges */
		for (ep0 = et.he[h]; ep0 != NULL; ep0 = epe) {

			/* Locate the end of the half edges at the same position */
			for (epe = ep0; epe != NULL; epe = epe->next) {
				if (epe == NULL || fabs(ep0->mj - epe->mj) > 1e-6)
					break;
			}

#ifdef DEBUG
			fprintf(stderr,"Doing group from %d to %d\n",ep0->ix, epe ? epe->ix : -1);
#endif

			/* Look for overlapping half edges, and split them up so */
			/* there are no overlaps */
			for (; ep0 != epe; ep0 = ep0->next) {

				for (ep1 = ep0->next; ep1 != epe; ep1 = ep1->next) {

					if (ep1->mi0 > ep0->mi1)	/* Out of range for overlap */
						break;

					/* If partial overlap */
					if (ep0->mi0 < (ep1->mi0-1e-6)
					 && ep0->mi1 > (ep1->mi0+1e-6)
					 && ep0->mi1 < (ep1->mi1-1e-6)) {
						hedge *ep0b, *ep1b;
						
#ifdef DEBUG
						fprintf(stderr,"Half edges partial overlap:\n");
						fprintf(stderr,"i = %d, j = %d\n",ep0->ix,ep1->ix);
						fprintf(stderr,"%s %d at %c = %f from %c = %f to %f, half %s\n",
						h == 0 ? "Vert" : "Horiz", ep0->ix,
						h == 0 ? 'X' : 'Y', ep0->mj,
						h == 0 ? 'Y' : 'X', ep0->mi0, ep0->mi1,
						ep0->negh ? "Neg" : "Pos");
						fprintf(stderr,"%s %d at %c = %f from %c = %f to %f, half %s\n",
						h == 0 ? "Vert" : "Horiz", ep1->ix,
						h == 0 ? 'X' : 'Y', ep1->mj,
						h == 0 ? 'Y' : 'X', ep1->mi0, ep1->mi1,
						ep1->negh ? "Neg" : "Pos");
#endif
						/* Split up the two edges so that we have four edges */
						
						if ((ep0b = (hedge *)calloc(sizeof(hedge), 1)) == NULL)
							error("Malloc of half edge structure failed");
						memcpy(ep0b, ep0, sizeof(hedge));

						if ((ep1b = (hedge *)calloc(sizeof(hedge), 1)) == NULL)
							error("Malloc of half edge structure failed");
						memcpy(ep1b, ep1, sizeof(hedge));

						ep0b->mi0 = ep1->mi0;
						ep1b->mi1 = ep0->mi1;
						ep0->mi1 = ep1->mi0;
						ep1->mi0 = ep0->mi1;

						/* Insert them in order into linked list */
						ep1b->next = ep1;
						ep0b->next = ep1b;
						ep0->next = ep0b;

						et.nhe[h] += 2;

					/* If full overlap */
					} else if (ep0->mi0 < (ep1->mi0-1e-6)
					 && ep0->mi1 > (ep1->mi1+1e-6)) {
						hedge *ep0b, *ep0c;
						
#ifdef DEBUG
						fprintf(stderr,"Half edges full overlap:\n");
						fprintf(stderr,"i = %d, j = %d\n",ep0->ix,ep1->ix);
						fprintf(stderr,"%s %d at %c = %f from %c = %f to %f, half %s\n",
						h == 0 ? "Vert" : "Horiz", ep0->ix,
						h == 0 ? 'X' : 'Y', ep0->mj,
						h == 0 ? 'Y' : 'X', ep0->mi0, ep0->mi1,
						ep0->negh ? "Neg" : "Pos");
						fprintf(stderr,"%s %d at %c = %f from %c = %f to %f, half %s\n",
						h == 0 ? "Vert" : "Horiz", ep1->ix,
						h == 0 ? 'X' : 'Y', ep1->mj,
						h == 0 ? 'Y' : 'X', ep1->mi0, ep1->mi1,
						ep1->negh ? "Neg" : "Pos");
#endif
						/* Split up the first edge so that we have four edges */
						
						if ((ep0b = (hedge *)calloc(sizeof(hedge), 1)) == NULL)
							error("Malloc of half edge structure failed");
						memcpy(ep0b, ep0, sizeof(hedge));

						if ((ep0c = (hedge *)calloc(sizeof(hedge), 1)) == NULL)
							error("Malloc of half edge structure failed");
						memcpy(ep0c, ep0, sizeof(hedge));

						ep0b->mi0 = ep1->mi0;
						ep0b->mi1 = ep1->mi1;
						ep0c->mi0 = ep1->mi1;
						ep0->mi1 = ep1->mi0;

						/* Insert them in order into linked list */
						ep0c->next = ep1->next;
						ep1->next = ep0c;
						ep0b->next = ep1;
						ep0->next = ep0b;
						et.nhe[h] += 2;
					}
				}
			}
		}

		/* Now we can assume that edges match completely or not at all */

		/* Re-create sorted list of vertical half edges */
		/* (Instead of this we could convert the half edge matching loops */
		/*  below to use the linked list, like the above code.) */

		if ((et.she[h] = (hedge **)malloc(sizeof(patch*) * et.nhe[h])) == NULL)
			error("Malloc of array of vertical halfedge pointers failed");
	
		for (ep0 = et.he[h], i = 0; i < et.nhe[h]; i++, ep0 = ep0->next)
			et.she[h][i] = ep0;
	
		/* Sort helf edges by their X location, then their Y0 location */
#define HEAP_COMPARE(A,B) (fabs(A->mj - B->mj) < 1e-6 ? A->mi0 < B->mi0 : A->mj < B->mj)
		HEAPSORT(hedge *, et.she[h], et.nhe[h]);
#undef HEAP_COMPARE

		et.nel[h] = 0;
		et.nelp = &et.el[h];		/* Append next edge list here */
		*et.nelp = NULL;			/* No edge list at this position yet */

		/* Create the edge list information */
		for (i = 0; i < et.nhe[h];) {
			int j, ii, nj;
			double *rgb = NULL;	/* Contrast RGB */
			elist *el;		/* Current elist */

			el = *et.nelp;	/* current end of lits */

			/* Locate the end of the half edges at the same position */
			for (ii = i; ii < et.nhe[h]; ii++) {
				if (fabs(et.she[h][i]->mj - et.she[h][ii]->mj) > 1e-6)
					break;
			}

#ifdef DEBUG
			fprintf(stderr,"Doing group from %d to %d\n",i, ii);
#endif

			/* Find half edge pairs */
			/* Note that we assume that the half edges match perfectly, */
			/* or not at all. */
			for (j = i; j < ii; j = nj, j++) {
				int e, k = j+1;
				double vv;

				if (k < ii
				 && et.she[h][j]->negh != et.she[h][k]->negh
				 && fabs(et.she[h][j]->mi0 - et.she[h][k]->mi0) < 1e-5
				 && fabs(et.she[h][j]->mi1 - et.she[h][k]->mi1) < 1e-5) {
					/* Found a matching half edge */

					nj = k;
					rgb = et.she[h][k]->rgb;

				} else if (k < ii	/* Assert */
				 && (   (et.she[h][j]->mi0+1e-6) < et.she[h][k]->mi1
				     && et.she[h][j]->mi1 > (et.she[h][k]->mi0+1e-6))) {

					/* Found an overlapping non-matching edge */
					nj = k;

#ifdef DEBUG
					fprintf(stderr,"Half edges overlap but don't match:\n");
					fprintf(stderr,"i = %d, j = %d\n",i,j);
					fprintf(stderr,"%s %d at %c = %f from %c = %f to %f, half %s\n",
					h == 0 ? "Vert" : "Horiz", i,
					h == 0 ? 'X' : 'Y', et.she[h][j]->mj,
					h == 0 ? 'Y' : 'X', et.she[h][j]->mi0, et.she[h][j]->mi1,
					et.she[h][j]->negh ? "Neg" : "Pos");
					fprintf(stderr,"%s %d at %c = %f from %c = %f to %f, half %s\n",
					h == 0 ? "Vert" : "Horiz", j,
					h == 0 ? 'X' : 'Y', et.she[h][k]->mj,
					h == 0 ? 'Y' : 'X', et.she[h][k]->mi0, et.she[h][k]->mi1,
					et.she[h][k]->negh ? "Neg" : "Pos");
#endif /* DEBUG */
						error("Internal - half edges don't match");

				} else {
					/* Must be a non-matching edge against the media */
					nj = j;
					rgb = et.mrgb;
				}

				/* Compute vector delta in rgb */
				/* Add entry to edge list */
				for (e = 0, vv = 0.0; e < 3; e++) {
					double tt = rgb[e] - et.she[h][j]->rgb[e];
					vv += tt * tt;
				}

//printf("~1 h %d, mj %f, mi %f, vv = %f\n",h, et.she[h][j]->mj, et.she[h][j]->mi0, vv);
				/* If edge is of sufficient magnitude */
				// ~~99
				if (vv > 0.2) {
					edge *ep;

					/* See if we need to add a first elist for this position */
					if (el == NULL) {		/* We do */
						if ((el = (elist *)calloc(sizeof(elist), 1)) == NULL)
							error("Malloc of elist structure failed");
						*et.nelp = el;
						el->pos = et.she[h][j]->mj;
						et.nel[h]++;
					}

					/* Add another edge entry */
					if ((ep = (edge *)calloc(sizeof(edge), 1)) == NULL)
						error("Malloc of edge structure failed");

					ep->next = el->e;
					ep->p1 = et.she[h][j]->mi0;
					ep->p2 = et.she[h][j]->mi1;

					el->e = ep;		/* Add to edge list */
					el->ne++;
				}
			}

			if (el != NULL) {
				/* We've done that position, so get ready for next */
				et.nelp = &el->next;		/* Append any more positions here */
				*et.nelp = NULL;
			}
			i = ii;		/* Start search for next group here */
		}
	}

	/* Figure the crossing count */
	for (h = 0; h < 2; h++) {
		int oh = 1 - h;	/* The other list */
		elist *el, *pl;		/* Current, previous elist */

		for (pl = NULL, el = et.el[h]; el != NULL; pl = el, el = el->next) {
			edge *ep;
			double pp, np;	/* Window in pos direction for crossing */

			if (pl != NULL)
				pp = (pl->pos + el->pos)/2.0;	/* Half distance to next line */
			else
				pp = -1e6;

			if (el->next != NULL)
				np = (el->next->pos + el->pos)/2.0;	/* Half distance to next line */
			else
				np = 1e6;

			/* For each edge on this edge position */
			for (ep = el->e; ep != NULL; ep = ep->next) {
				elist *oel;		/* Other edge list pointer */

				/* For each edge in other list */
				for (oel = et.el[oh]; oel != NULL; oel = oel->next) {
					edge *oep;

					if (oel->pos < ep->p1 || oel->pos > ep->p2)
						continue;	/* Other edge doesn't intersect this segment */

					for (oep = oel->e; oep != NULL; oep = oep->next) {

						/* If crosses on this line within +-0.5 of line each side */
						if (oep->p1 <= np && oep->p2 >= pp) {
							el->cc++;	/* Count crossing */
						}

					}
				}
			}
		}
	}

	/* Compute and normalise the length (strength) & crossing count of each edge */
	for (h = 0; h < 2; h++) {
		elist *el;		/* Current elist */
		double maxlen = 0.0;  
		double maxcc = 0.0;  

		for (el = et.el[h]; el != NULL; el = el->next) {
			edge *ep;
			double tlen;

			for (tlen = 0.0, ep = el->e; ep != NULL; ep = ep->next) {
				tlen += ep->p2 - ep->p1;
			}
			el->len = tlen;
			if (maxlen < tlen)
				maxlen = tlen;
			if (maxcc < el->cc)
				maxcc = el->cc;
		}

		/* Normalise */
		for (el = et.el[h]; el != NULL; el = el->next) {
			el->len /= maxlen;
			el->cc  /= maxcc;
		}
	}

	/* Output the .cht file */
	if ((of = fopen(fname,"w")) == NULL)
			error ("Unable to open output file '%s' for writing",fname);
	
	fprintf(of,"\n\n");
	fprintf(of, "BOXES %d\n",et.npatches);


	/* Locate fiducials if they've been added to chart */
	if (et.nfid == 4) {
		fprintf(of, "  F _ _ %f %f %f %f %f %f %f %f\n",
		        et.fx[0], et.fy[0], et.fx[1], et.fy[1], et.fx[2], et.fy[2], et.fx[3], et.fy[3]);
	}

	{
		int fidc = 0;
		patch *pp;
		double mins;		/* Minimum sample box size */

		for (pp = et.patches; pp != NULL; pp = pp->next) {

			if (pp->id[0] == '\000') {
				fprintf(of, "  D MARK%d MARK%d _ _ %f %f %f %f 0 0\n",
				        fidc, fidc, pp->w, pp->h, pp->xo, pp->yo);
				fidc++;
			}
		}
		
		mins = 1e6;
		for (pp = et.patches; pp != NULL; pp = pp->next) {

			if (pp->id[0] != '\000') {
				fprintf(of, "  X %s %s _ _ %f %f %f %f 0 0\n",
				        pp->id, pp->id, pp->w, pp->h, pp->xo, pp->yo);
				if (mins > pp->w)
					mins = pp->w;
				if (mins > pp->h)
					mins = pp->h;
			}
		}
		fprintf(of,"\n");
		
		/* Use a 15% box shrink */
		fprintf(of, "BOX_SHRINK %f\n", mins * 0.15);
		fprintf(of,"\n");

	}

	fprintf(of,"REF_ROTATION %f\n", 0.0);
	fprintf(of,"\n");

	{
		elist *el;		/* Current elist */

		fprintf(of,"XLIST %d\n",et.nel[0]);
		for (el = et.el[0]; el != NULL; el = el->next)
			fprintf(of,"  %f %f %f\n",el->pos, el->len, el->cc);
		fprintf(of,"\n");
	
		fprintf(of,"YLIST %d\n",et.nel[1]);
		for (el = et.el[1]; el != NULL; el = el->next)
			fprintf(of,"  %f %f %f\n",el->pos, el->len, el->cc);
		fprintf(of,"\n");
	
		fprintf(of,"\n");
	}


	fprintf(of, "EXPECTED XYZ %d\n",ei - si);

	for (i = si; i < ei; i++) {
		int ix = rix[i];
		fprintf(of, "  %s %f %f %f\n", cols[ix].loc, 100.0 * cols[ix].XYZ[0], 100.0 * cols[ix].XYZ[1], 100.0 * cols[ix].XYZ[2]);
	}
	fprintf(of,"\n");

	if (fclose(of))
		error ("Unable to close output file '%s'",fname);
}


/* Cleanup any allocation */
void et_clear(void) {
	int h;
	patch *p;

	et.nfid = 0;

	for (h = 0; h < 2; h++) {
		hedge *he;
		elist *el;

		/* Free up half edges */
		he = et.he[h];
		while (he != NULL) {
			hedge *ne = he->next;
			free(he);
			he = ne;
		}
		et.nhe[h] = 0;
		et.he[h] = NULL;
	
		/* Free up sorted half edge lists */
		if (et.she[h] != NULL)
			free(et.she[h]);
		et.she[h] = NULL;

		/* Free up edge lists and edges */
		el = et.el[h];
		while (el != NULL) {
			elist *nel;
			edge  *ep;

			ep = el->e;
			while (ep != NULL) {
				edge *nep;

				nep = ep->next;
				free(ep);
				ep = nep;
			}
			el->ne = 0;
			el->e = NULL;

			nel = el->next;
			free(el);
			el = nel;
		}
		et.nel[h] = 0;
		et.el[h] = NULL;
	}
	et.nelp = NULL;

	p = et.patches;
	while (p != NULL) {
		patch *np = p->next;
		free(p);
		p = np;
	}
	et.patches = NULL;
	et.npatches = 0;
}








