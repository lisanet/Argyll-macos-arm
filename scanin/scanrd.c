
/*
 * Raster Color Target Scan Input module
 * This is the core chart recognition code.
 * 
 * Author: Graeme Gill
 *
 * Copyright 1995 - 2008 Graeme W. Gill, All right reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * To Do: 
 *		  Add option to output a raster file made from the .cht and example values.
 *
 *        Fix sboxes parameters/digitization to fix "droop" in box areas.
 *        Scale parameters with image size.
 *        To handle high res, introduce automatic sub-sampler.
 *		  Change reference parser to make it more forgiving - use cgats parser ?
 */

#undef DEBUG

#define VERSION "1.0"

/* Behaviour defines */
#undef DIAGN				/* Allow diagonal connectivity of groups */
#define AA_LINES			/* Plot diagnostics using anti-aliased lines */

#define MATCHCC 0.25		/* Match correlation threshold - reject any match under this */
							/* (Might want to be able to override this in command line) */

#define ALT_ROT_TH 0.7		/* Correlation threshold of alternate rotations to be greater than this */

#define TH (20.0 * 20.0)	/* Initial color change threshhold */

#undef DBG
#define dbgo stdout
#define DBG(aaa) fprintf aaa, fflush(dbgo)

#include <stdio.h>
/* #include <fcntl.h> */		/* In case DOS binary stuff is needed */
#include <string.h>
#include <math.h>

#include <stdlib.h>
#include <sys/stat.h>
/*  #include <fname.h> */

#include "numlib.h"
#include "scanrd_.h"

/* ------------------------------------------------- */
/* Implementations of public functions */
static void free_scanrd(scanrd *s);
static int scanrd_reset(scanrd *s);
static int scanrd_read(scanrd *ps, char *id, double *P, double *mP,
                       double *sdP, int *cnt);
static unsigned int scanrd_error(scanrd *s, char **errm);

/* Forward internal function declaration */
static scanrd_ *new_scanrd(int flags, int verb, double gammav,
	int (*write_line)(void *ddata, int y, char *src), void *ddata,
	int w, int h, int d, int td, int p,
	int (*read_line)(void *fdata, int y, char *dst), void *fdata,
	char *refname);
static int read_input(scanrd_ *s);
static int calc_lines(scanrd_ *s);
static int show_lines(scanrd_ *s);
static int calc_perspective(scanrd_ *s);
static int calc_rotation(scanrd_ *s);
static int calc_elists(scanrd_ *s, int ref);
static int write_elists(scanrd_ *s);
static int read_relists(scanrd_ *s);
static int do_match(scanrd_ *s);
static int compute_ptrans(scanrd_ *s);
static int compute_man_ptrans(scanrd_ *s, double *sfids);
static int improve_match(scanrd_ *s);
static int setup_sboxes(scanrd_ *s);
static int do_value_scan(scanrd_ *s);
static int compute_xcc(scanrd_ *s);
//static int restore_best(scanrd_ *s);
static int show_sbox(scanrd_ *s);
static int show_groups(scanrd_ *s);
static int scanrd_write_diag(scanrd_ *s);
static void toRGB(unsigned char *dst, unsigned char *src, int depth, int bpp);
static void XYZ2Lab(double *out, double *in);
static void pval2Lab(double *out, double *in, int depth);
/* ------------------------------------------------- */

/* Read in a chart, and either create a reference or make values available, */
/* by using reset() and read() to get values read */
scanrd *do_scanrd(
int flags,			/* option flags */
int verb,			/* verbosity level */

double gammav,		/* Apprimate gamma encoding of image (0.0 = default 2.2) */
double *sfid,		/* Specified four fiducials x1, y1 .. x4, y4, NULL if auto recognition */
					/* Typical clockwise from top left */

int w, int h, 		/* Width and Height of input raster in pixels */
int d, int td, int p,		/* Useful plane depth, Total depth, Bit presision of input pixels */
int (*read_line)(void *fdata, int y, char *dst),	/* Read RGB line of source file */
void *fdata,		/* Opaque data for read_line */

char *refname,		/* reference file name */

int (*write_line)(void *ddata, int y, char *src),	/* Write RGB line of diag file */
void *ddata			/* Opaque data for write_line */
) {
	scanrd_ *s;

	/* allocate the basic object */
	if (verb >= 2)
		DBG((dbgo,"About to allocate scanrd_ object\n"));
	if ((s = new_scanrd(flags, verb, gammav, write_line, ddata, w, h, d, td, p, read_line, fdata, refname)) == NULL)
		return NULL;

	if (s->errv != 0)	/* Some other error from new_scanrd() */
		return (scanrd *)s;

	if (s->verb >= 2)
		DBG((dbgo,"About to read input tiff file and discover groups\n"));
	if (read_input(s))
		goto sierr;			/* Error */

	if (s->flags & SI_SHOW_GROUPS)
		if (show_groups(s))
			goto sierr;		/* Error */
	
	if (s->verb >= 2)
		DBG((dbgo,"About to calculate edge lines\n"));
	if (calc_lines(s))
		goto sierr;			/* Error */
	if (s->verb >= 2)
		DBG((dbgo,"%d useful edges out of %d\n",s->novlines, s->noslines));

	if (s->flags & SI_PERSPECTIVE) {
		if (s->verb >= 2)
			DBG((dbgo,"About to calculate perspective correction\n"));
		if (calc_perspective(s)) {
			if (s->flags & SI_SHOW_LINES) {
				s->flags &= ~SI_SHOW_PERS;	/* Calc perspective failed! */
				s->flags &= ~SI_SHOW_ROT;	/* Calc rotation not done! */
				show_lines(s);
			}
			goto sierr;			/* Error */
		}
	}

	if (s->verb >= 2)
		DBG((dbgo,"About to calculate rotation\n"));
	if (calc_rotation(s)) {
		if (s->flags & SI_SHOW_LINES) {
			s->flags &= ~SI_SHOW_ROT;	/* Calc rotation failed! */
			show_lines(s);
		}
		goto sierr;			/* Error */
	}

	if (s->flags & SI_BUILD_REF) { /* If generating a chart reference file */
		/* Calculate the edge lists and write it to the file */
		if (s->verb >= 2)
			DBG((dbgo,"About to build feature information\n"));
		if (calc_elists(s, 1))		/* reference */
			goto sierr;		/* Error */

		if (s->verb >= 2)
			DBG((dbgo,"About to write feature reference information\n"));
		if (write_elists(s))
			goto sierr;		/* Error */
	} else {
		/* If we are matching to the reference and generating an output data file */
		int rv;

		/* Calculate the edge lists read for a match */
		if (s->verb >= 2)
			DBG((dbgo,"About to calculate feature information\n"));
		if (calc_elists(s, 0))		/* match */
			goto sierr;		/* Error */

		if (s->verb >= 2)
			DBG((dbgo,"About to read reference feature information\n"));
		if (read_relists(s))
			goto sierr;		/* Error */
		if (s->verb >= 2)
			DBG((dbgo,"Read of chart reference file succeeded\n"));

		if (sfid != NULL) {		/* Manual matching */
			if (s->verb >= 2)
				DBG((dbgo,"Using manual matching\n"));

			if (s->havefids == 0) {
				s->errv = SI_NO_FIDUCIALS_ERR;
				sprintf(s->errm,"Chart recognition definition file doesn't contain fiducials");
				goto sierr;		/* Error */
			}
			if (compute_man_ptrans(s, sfid))
				goto sierr;

			/* Do the actual scan given out manual transformation matrix */ 
			if (s->verb >= 2)
				DBG((dbgo,"About to setup value scanrdg boxes\n"));
			if (setup_sboxes(s))	
				goto sierr;
			if (s->verb >= 2)
				DBG((dbgo,"About to read raster values\n"));
			if (do_value_scan(s))
				goto sierr;

		} else {				/* Automatic matching */

			/* Attempt to match input file with reference */
			if (s->verb >= 2)
				DBG((dbgo,"About to match features\n"));
			if ((rv = do_match(s)) != 0) {
				if (rv == 1) {	/* No reasonable rotation found */
					s->errv = SI_POOR_MATCH;
					sprintf(s->errm,"Pattern match wasn't good enough");
				}
				goto sierr;
			}

			/* If there is patch matching data and more than one */
			/* feasible matching rotation, try and discriminate between them. */
			if (s->xpt && s->norots > 1) {
				int i, j;
				int flags = s->flags;
			
				s->flags &= ~SI_SHOW_SAMPLED_AREA;	/* Don't show areas for trials */

				/* For each candidate rotation, scan in the pixel values */
				for (s->crot = 0; s->crot < s->norots; s->crot++) {
		
					/* Compute transformation from reference to input file */
					if (s->verb >= 2)
						DBG((dbgo,"About to compute match transform for rotation %f deg.\n",
						          DEG(s->rots[s->crot].irot)));
					if (compute_ptrans(s))
						goto sierr;
			
					/* Setup the input boxes ready for scanning in the input values */
					if (s->verb >= 2)
						DBG((dbgo,"About to setup value scanrdg boxes\n"));
					if (setup_sboxes(s))	
						goto sierr;
			
					/* Scan in the pixel values */
					if (s->verb >= 2)
						DBG((dbgo,"About to read raster values\n"));
					if (do_value_scan(s))
						goto sierr;
		
					/* Copy to this rotation values so that the best can be restored */
					if (s->xpt != 0) {			/* Got expected patch values to compare with */
						if (s->verb >= 2)
							DBG((dbgo,"About to compute expected value correlation\n"));
						if (compute_xcc(s))
							goto sierr;
					}
				}

				/* Pick the best from the candidate rotation */
				if (s->verb >= 2) {
					DBG((dbgo,"Expected value distance values are:\n"));
					for (i = 0; i < s->norots; i++) {
						DBG((dbgo,"%d, rot %f: %f\n", i, DEG(s->rots[i].irot), s->rots[i].xcc));
					}
				}

				for (j = 0, i = 1; i < s->norots; i++) {
					if (s->rots[i].xcc < s->rots[j].xcc)
						j = i;
				}

				if (s->verb >= 2)
					DBG((dbgo,"Chosen rotation %f deg. as best\n",DEG(s->rots[j].irot)));
	
				s->crot = j;
				s->flags = flags;		/* Restore flags */
			}

			/* Setup transformation to be that for chosen rotation for diagnostics */
			if (s->verb >= 2)
				DBG((dbgo,"About to compute final match transform\n"));
			if (compute_ptrans(s))
				goto sierr;

			if (s->verb >= 2)
				DBG((dbgo,"Improve match\n"));
			if (improve_match(s))
				goto sierr;

			/* After choosing rotation of improving the fit, rescan the values */
			if (s->verb >= 2)
				DBG((dbgo,"About to setup value scanrdg boxes\n"));
			if (setup_sboxes(s))	
				goto sierr;
			if (s->verb >= 2)
				DBG((dbgo,"About to read raster values\n"));
			if (do_value_scan(s))
				goto sierr;
		}

		if (s->flags & SI_SHOW_SBOX) {
			show_sbox(s);		/* Draw sample box outlines on diagnostic raster */
		}
	}
	if (s->flags & SI_SHOW_LINES)
		if(show_lines(s))
			goto sierr;		/* Error */

sierr:;
	if (s->verb >= 2)
		DBG((dbgo,"About to write diag file\n"));
	if (scanrd_write_diag(s))
		return (scanrd *)s;	/* Error */

	return (scanrd *)s;
}


/********************************************************************************/

/* Allocate the basic scanrd object */
/* Return NULL on failure to allocate */
/* Need to check errv for other problems */
static scanrd_
*new_scanrd(
	int flags,			/* option flags */
	int verb,			/* verbosity level */
	double gammav,		/* Approximate gamma encoding of image (0.0 = default 2.2) */
	int (*write_line)(void *ddata, int y, char *src),	/* Write RGB line of diag file */
	void *ddata,		/* Opaque data for write_line() */
	int w, int h, 		/* Width and Height of input raster in pixels */
	int d, int td, int p,	/* Useful plane Depth, Total depth, Bit presision of input pixels */
	int (*read_line)(void *fdata, int y, char *dst),	/* Read RGB line of source file */
	void *fdata,		/* Opaque data for read_line() */

	char *refname		/* reference file name */
) {
	scanrd_ *s;

	if ((s = (scanrd_ *)calloc(1, sizeof(scanrd_))) == NULL)
		return NULL;

	/* Public functions */
	s->public.reset = scanrd_reset;
	s->public.read = scanrd_read;
	s->public.error = scanrd_error;
	s->public.free = free_scanrd;

	if (flags & (SI_SHOW_ROT | SI_SHOW_PERS | SI_SHOW_IMPL | SI_SHOW_ALL_LINES))
		flags |= SI_SHOW_LINES;		/* Key all line stuff off SI_SHOW_LINES */

	if (flags & (SI_SHOW_SBOX_OUTLINES | SI_SHOW_SBOX_NAMES | SI_SHOW_SBOX_AREAS))
		flags |= SI_SHOW_SBOX;;		/* Key all sample box stuff off SI_SHOW_SBOX */

	if (write_line == NULL)
		flags &= ~SI_SHOW_FLAGS;	/* If no diag file, turn off show flags */

	s->flags = flags;
	s->verb = verb;

	s->errv = 0;
	s->errm[0] = '\0';

	if (gammav <= 0.0)
		gammav = 2.2;			/* default */
	s->gammav = gammav;
	s->width = w;
	s->height = h;
	s->depth = d;
	s->tdepth = td;
	s->bpp = p;

	if (d > MXDE)  {
		s->errv = SI_PIX_DEPTH_ERR;
		sprintf(s->errm,"scanrd: Pixel depth is too large");
		return s;
	}

	if (p != 8 && p != 16) {
		s->errv = SI_BIT_DEPTH_ERR;
		sprintf(s->errm,"scanrd: Pixel bits/pixel is not 8 or 16");
		return s;
	}
	if (p == 8)
		s->bypp = 1;
	else
		s->bypp = 2;

	if (verb >= 2)
		DBG((dbgo,"Verbosity = %d, flags = 0x%x\n",verb, flags));

	/* RGB Diagnostic output raster array requested */
	if ((flags & SI_SHOW_FLAGS) && write_line != NULL) {
		if ((s->out = malloc(3 * w * h)) == NULL) {
			s->errv = SI_MALLOC_DIAG_RAST;
			sprintf(s->errm,"scanrd: Diagnostic output raster array malloc failed");
			return s;
		}
	}
	
	s->noslines = 0;
	s->novlines = 0;
	s->gdone = NULL;
	s->irot = 0.0;
	s->norots = 0;
	
	s->ppc[0] = 0.0;
	s->ppc[1] = 0.0;
	s->ppc[2] = 0.0;
	s->ppc[3] = 0.0;

	/* Set overall perspective transform to null */
	s->ptrans[0] = 1.0;
	s->ptrans[1] = 0.0;
	s->ptrans[2] = 0.0;
	s->ptrans[3] = 0.0;
	s->ptrans[4] = 1.0;
	s->ptrans[5] = 0.0;
	s->ptrans[6] = 0.0;
	s->ptrans[7] = 0.0;

	INIT_ELIST(s->xelist);
	INIT_ELIST(s->yelist);
	INIT_ELIST(s->ixelist);
	INIT_ELIST(s->iyelist);
	INIT_ELIST(s->rxelist);
	INIT_ELIST(s->ryelist);
	s->rbox_shrink = 0.9;
	s->xpt = 0;
	
	s->nsbox = 0;
	s->sboxes = NULL;
	s->sbstart = NULL;
	s->sbend = NULL;
	s->csi = 0;
	s->cei = 0;
	s->alist = NULL;
	
	s->next_read = 0;

	s->refname = refname;

	s->inited = 0;
	s->vrego = s->vregn = NULL;
	s->no_vo = s->no_vn = 0;
	s->hrego = s->hregn = NULL;
	s->no_ho = s->no_hn = 0;
	s->th = TH;
	s->divval = 0.25;
	s->adivval = 0.0;
	s->divc = 0;

	/* aa line init */
	s->aa_inited = 0;	/* Let line init do the rest */
	s->coverage = NULL;

	/* Callbacks */
	s->read_line = read_line;
	s->fdata = fdata;

	s->write_line = write_line;
	s->ddata = ddata;

	return s;
}

static void free_elist_array(elist *el);

/* Free the object up */
static void
free_scanrd(
scanrd *ps
) {
	scanrd_ *s = (scanrd_ *)ps;	/* Cast public to private */
	points *tp;

	free_elist_array(&s->xelist);
	free_elist_array(&s->yelist);
	free_elist_array(&s->ixelist);
	free_elist_array(&s->iyelist);
	free_elist_array(&s->rxelist);
	free_elist_array(&s->ryelist);
	
	if (s->sboxes != NULL)
		free(s->sboxes);
	if (s->sbstart != NULL)
		free(s->sbstart);
	if (s->sbend != NULL)
		free(s->sbend);
	s->alist = NULL;
	
	/* Free up done line list */
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->r != NULL)
			free(tp->r);
		free(tp);
	END_FOR_ALL_ITEMS(tp);
	s->gdone = NULL;
	
	/* Points were deleted with gdone ??? */
	if(s->vrego != NULL)
		free(s->vrego);
	if(s->vregn)
		free(s->vregn);
	if(s->hrego != NULL)
		free(s->hrego);
	if(s->hregn != NULL)
		free(s->hregn);
	s->inited = 1;

	/* Free up output diag array */
	if (s->out != NULL)
		free(s->out);
	
	/* Free up aa line array */
	if (s->coverage != NULL)
		free(s->coverage);
	free(s);
}


/* Return the error flag, and set the message pointer */
static unsigned int
scanrd_error(scanrd *ps, char **errm) {
	scanrd_ *s = (scanrd_ *)ps;	/* Cast public to private */
	*errm = s->errm;
	return s->errv;
}

/********************************************************************************/
static int analize(scanrd_ *s, unsigned char *inp[6], int y);

/* Read in and process the input file */
/* Return non-zero on error */
static int
read_input(scanrd_ *s) {
	unsigned char *in[6];		/* Pointer to six input buffers */
	int w = s->width;			/* Raster width */
	int h = s->height;			/* Raster height */
	int i, y;

	/* Allocate input line buffers */
	for (i = 0; i < 6; i++) {
		if ((in[i] = malloc(s->tdepth * w * s->bypp)) == NULL) {
			s->errv = SI_MALLOC_INPUT_BUF;
			sprintf(s->errm,"scanrd: Failed to malloc input line buffers");
			return 1;
		}
	}

	/* Prime the input buffers with 5 lines */
	for (y = 0; y < 5; y++) {
		if (s->read_line(s->fdata, y, (char *)in[y])) {
			s->errv = SI_RAST_READ_ERR;
			sprintf(s->errm,"scanrd: read_line() returned error");
			return 1;
		}
	}
	/* Process the tiff file line by line (Assume at least 6 lines in total raster) */
	for (; y < h; ++y) {
		unsigned char *tt;
		if (s->read_line(s->fdata, y, (char *)in[5])) {
			s->errv = SI_RAST_READ_ERR;
			sprintf(s->errm,"scanrd: read_line() returned error");
			return 1;
		}

		if (analize(s, in, y)) {
			return 1;
		}

		tt = in[0];		/* Shuffle buffers about */
		in[0] = in[1];
		in[1] = in[2];
		in[2] = in[3];
		in[3] = in[4];
		in[4] = in[5];
		in[5] = tt;

	}
	s->adivval /= (double)s->divc;	/* Average divider value, 1.0 = 0 degrees, 0.0 = 45 degrees */
	if (s->adivval < 0.0)
		s->adivval = 0.0;
	else if (s->adivval > 1.0)
		s->adivval = 1.0;

	if (s->verb >= 2)
		DBG((dbgo,"adivval = %f\n",s->adivval));

	/* Free the input line buffers */
	for (i = 0; i < 6; i++)
		free(in[i]);

	return 0;
}

/********************************************************************************/

#ifdef NEVER	/* Before 22/5/2004 */
#define	THRN 1.0			/* Threshold above average ratio - numerator */
#define	THRD 2.0			/* Threshold above average ratio - denominator */

#define THAWF 1.0 			/* Threshold average adaptation filter weight, fixed value (TH) */
#define THAWP 4.0 			/* Threshold average adaptation filter weight, previous value */
#define THAWN 1.0 			/* Threshold average adaptation filter weight, new value */

#else			/* Current values */

#define	THRN 1.0			/* Threshold above average ratio - numerator */
#define	THRD 1.5			/* Threshold above average ratio - denominator */

#define THAWF 1.0 			/* Threshold average adaptation filter weight, fixed value (TH) */
#define THAWP 5.0 			/* Threshold average adaptation filter weight, previous value */
#define THAWN 1.0 			/* Threshold average adaptation filter weight, new value */

#endif

/* ~~~ minimum raster size needs to be specified/checked ~~~~ */
#define MIN_NO_LINES    16		/* Minimum number of valid fitted lines to estimate rotation */

/* Criteria for accepting lines for angle calculation (valid lines) */
#define MAX_MWID_TO_LEN    0.1
#define MIN_POINT_TO_AREA  0.9	/* Minimum point desity over the lines area */
#define SD_WINDOW          1.5	/* Allow += 1.5 of a standard deviation for robust angle calc. */
#define ELISTCDIST 800		/* 1/ELISTCDIST = portion of refence edge list legth to coalesce over */

/* Criteria for accepting lines for improring final fit */
#define IMP_MATCH 0.10			/* Proportion of average tick spacing */

/* The following should be scaled to the resolution of the image ? */
#define MIN_POINTS 10		/* Minimum points to calculate line */
#define MIN_LINE_LENGTH    10.0
#define CUT_CHUNKS 128		/* cut groups along diagonals - must be power of 2 */

static int add_region(scanrd_ *s, region *rego, int no_o, region *regn, int no_n, int y);

/* Process a line of the TIFF file */
/* return non-zero on error */
static int
analize(
scanrd_ *s,
unsigned char *inp[6],		/* current and previous 5 lines */
int y						/* Current line y */
) {
	int w = s->width;
	int stride = s->tdepth * s->width;	/* In pixels */
	unsigned short *gamma = s->gamma;
	int x,i;
	unsigned short *inp2[6];	/* current and previous 5 lines (16bpp) equivalent of inp[] */
	unsigned char  *in[6];		/* six input lines (8bpp) */
	unsigned short *in2[6];		/* six input lines (16bpp) */
	region *tr;
	double tdh,tdv;				/* Horizontal/virtical detect levels */
	double tdmag;
	double atdmag = 0.0;		/* Average magnitude over a line */
	int atdmagc = 0;			/* Average magnitude over a line count */
	double linedv = 0.0;		/* Lines average divider value */
	int linedc = 0;				/* Lines average count */
	int xo3 = s->tdepth * 3;	/* Xoffset by 3 pixels */
	int xo2 = s->tdepth * 2;	/* Xoffset by 2 pixels */
	int xo1 = s->tdepth * 1;	/* Xoffset by 1 pixels */

	for (x = 0; x < 6; x++)		/* Create 16 bpp version of line pointers */
		inp2[x] = (unsigned short *)inp[x];

	if (s->inited == 0) {
		/* Init gamma conversion lookup and region tracking. */
		/* The assumption is that a typical chart has an approx. visually */
		/* uniform distribution of samples, so that a typically gamma */
		/* encoded scan image will have an average pixel value of 50%. */
		/* If a the chart has a different gamma encoding (ie. linear), */
		/* then we convert it to gamma 2.2 encoded to (hopefuly) enhance */
		/* the patch contrast. */
		if (s->bpp == 8)
		    for (i = 0; i < 256; i++) {
				int byteb1;
			
				byteb1 = (int)(0.5 + 255 * pow( i / 255.0, s->gammav/2.2 ));
				gamma[i] = byteb1;
			}
		else
		    for (i = 0; i < 65536; i++) {
				int byteb1;
			
				byteb1 = (int)(0.5 + 65535 * pow( i / 65535.0, s->gammav/2.2 ));
				gamma[i] = byteb1;
			}

		if ((s->vrego = (region *) malloc(sizeof(region) * (w+1)/2)) == NULL) {
			s->errv = SI_MALLOC_VREGION;
			sprintf(s->errm,"vreg malloc failed");
			return 1;
		}
		s->no_vo = 0;
		if ((s->vregn = (region *) malloc(sizeof(region) * (w+1)/2)) == NULL) {
			s->errv = SI_MALLOC_VREGION;
			sprintf(s->errm,"vreg malloc failed");
			return 1;
		}
		s->no_vn = 0;
		if ((s->hrego = (region *) malloc(sizeof(region) * (w+1)/2)) == NULL) {
			s->errv = SI_MALLOC_VREGION;
			sprintf(s->errm,"vreg malloc failed");
			return 1;
		}
		s->no_ho = 0;
		if ((s->hregn = (region *) malloc(sizeof(region) * (w+1)/2)) == NULL) {
			s->errv = SI_MALLOC_VREGION;
			sprintf(s->errm,"vreg malloc failed");
			return 1;
		}
		s->no_hn = 0;
		INIT_LIST(s->gdone);
		s->inited = 1;
	}

	/* Un-gamma correct the latest input line */
	if (s->bpp == 8)
		for (x = 0; x < stride; x++)
			inp[5][x] = (unsigned char)gamma[inp[5][x]];
	else
		for (x = 0; x < stride; x++)
			inp2[5][x] = gamma[inp2[5][x]];

	/* Compute difference output for line y-3 */
	atdmagc = w - 5;		/* Magnitude count (to compute average) */
	for (x = 3; x < (w-2); x++) {		/* Allow for -3 to +2 from x */
		unsigned char *out = s->out;
		int e;
		int ss;
		int idx = ((y-2) * w + x) * 3;		/* Output raster index in bytes */

		if (s->bpp == 8)
			for (i = 0; i < 6; i++)
				in[i] = inp[i] + x * s->tdepth;	/* Strength reduce */
		else
			for (i = 0; i < 6; i++) {
				in2[i] = inp2[i] + x * s->tdepth;	/* Strength reduce */
				in[i] = (unsigned char *)in2[i];	/* track 8bpp pointers */
			}

		if (s->flags & SI_SHOW_IMAGE) {		/* Create B&W image */
			toRGB(out + idx, in[2], s->depth, s->bpp);		/* Convert to RGB */
			out[idx] = out[idx+1] = out[idx+2] = (2 * out[idx] + 7 * out[idx+1] + out[idx+2])/10;
		}

		ss = 0;		/* Sign of cross components the same vote */
		tdh = tdv = 0.0;
		
		if (s->bpp == 8)
			for (e = 0; e < s->depth; e++) {
				int d1,d2;
				/* Compute Gxp */
				d1 = -in[0][-xo3+e] + -in[0][-xo2+e] + -in[0][-xo1+e]
				                              + -in[0][ 0+e] + -in[0][ xo1+e] + -in[0][ xo2+e]
				   + -in[1][-xo3+e] + -in[1][-xo2+e] + -in[1][-xo1+e]
				                              + -in[1][ 0+e] + -in[1][ xo1+e] + -in[1][ xo2+e] 
				   + -in[2][-xo3+e] + -in[2][-xo2+e] + -in[2][-xo1+e]
				                              + -in[2][ 0+e] + -in[2][ xo1+e] + -in[2][ xo2+e]
				   +  in[3][-xo3+e] +  in[3][-xo2+e] +  in[3][-xo1+e]
				                              +  in[3][ 0+e] +  in[3][ xo1+e] +  in[3][ xo2+e]
				   +  in[4][-xo3+e] +  in[4][-xo2+e] +  in[4][-xo1+e]
				                              +  in[4][ 0+e] +  in[4][ xo1+e] +  in[4][ xo2+e]
				   +  in[5][-xo3+e] +  in[5][-xo2+e] +  in[5][-xo1+e]
				                              +  in[5][ 0+e] +  in[5][ xo1+e] +  in[5][ xo2+e];
				/* Compute Gyp */
				d2 = -in[0][-xo3+e] + -in[1][-xo3+e] + -in[2][-xo3+e]
				                              + -in[3][-xo3+e] + -in[4][-xo3+e] + -in[5][-xo3+e]
				   + -in[0][-xo2+e] + -in[1][-xo2+e] + -in[2][-xo2+e]
				                              + -in[3][-xo2+e] + -in[4][-xo2+e] + -in[5][-xo2+e]
				   + -in[0][-xo1+e] + -in[1][-xo1+e] + -in[2][-xo1+e]
				                              + -in[3][-xo1+e] + -in[4][-xo1+e] + -in[5][-xo1+e]
				   +  in[0][   0+e] +  in[1][   0+e] +  in[2][   0+e]
				                              +  in[3][   0+e] +  in[4][   0+e] +  in[5][   0+e]
				   +  in[0][+xo1+e] +  in[1][+xo1+e] +  in[2][+xo1+e]
				                              +  in[3][+xo1+e] +  in[4][+xo1+e] +  in[5][+xo1+e]
				   +  in[0][+xo2+e] +  in[1][+xo2+e] +  in[2][+xo2+e]
				                              +  in[3][+xo2+e] +  in[4][+xo2+e] +  in[5][+xo2+e];

				if ((d1 >= 0 && d2 >=0)
				 || (d1 < 0 && d2 < 0))
					ss++;				/* Sign was the same */
				tdh += d1/4.5 * d1/4.5;		/* (4.5 = 6x6/4x2, to scale original tuned values) */
				tdv += d2/4.5 * d2/4.5;
			}
		else
			for (e = 0; e < s->depth; e++) {
				int d1,d2;
				/* Compute Gxp */
				d1 = -in2[0][-xo3+e] + -in2[0][-xo2+e] + -in2[0][-xo1+e]
				                              + -in2[0][ 0+e] + -in2[0][ xo1+e] + -in2[0][ xo2+e]
				   + -in2[1][-xo3+e] + -in2[1][-xo2+e] + -in2[1][-xo1+e]
				                              + -in2[1][ 0+e] + -in2[1][ xo1+e] + -in2[1][ xo2+e] 
				   + -in2[2][-xo3+e] + -in2[2][-xo2+e] + -in2[2][-xo1+e]
				                              + -in2[2][ 0+e] + -in2[2][ xo1+e] + -in2[2][ xo2+e]
				   +  in2[3][-xo3+e] +  in2[3][-xo2+e] +  in2[3][-xo1+e]
				                              +  in2[3][ 0+e] +  in2[3][ xo1+e] +  in2[3][ xo2+e]
				   +  in2[4][-xo3+e] +  in2[4][-xo2+e] +  in2[4][-xo1+e]
				                              +  in2[4][ 0+e] +  in2[4][ xo1+e] +  in2[4][ xo2+e]
				   +  in2[5][-xo3+e] +  in2[5][-xo2+e] +  in2[5][-xo1+e]
				                              +  in2[5][ 0+e] +  in2[5][ xo1+e] +  in2[5][ xo2+e];
				/* Compute Gyp */
				d2 = -in2[0][-xo3+e] + -in2[1][-xo3+e] + -in2[2][-xo3+e]
				                              + -in2[3][-xo3+e] + -in2[4][-xo3+e] + -in2[5][-xo3+e]
				   + -in2[0][-xo2+e] + -in2[1][-xo2+e] + -in2[2][-xo2+e]
				                              + -in2[3][-xo2+e] + -in2[4][-xo2+e] + -in2[5][-xo2+e]
				   + -in2[0][-xo1+e] + -in2[1][-xo1+e] + -in2[2][-xo1+e]
				                              + -in2[3][-xo1+e] + -in2[4][-xo1+e] + -in2[5][-xo1+e]
				   +  in2[0][   0+e] +  in2[1][   0+e] +  in2[2][   0+e]
				                              +  in2[3][   0+e] +  in2[4][   0+e] +  in2[5][   0+e]
				   +  in2[0][+xo1+e] +  in2[1][+xo1+e] +  in2[2][+xo1+e]
				                              +  in2[3][+xo1+e] +  in2[4][+xo1+e] +  in2[5][+xo1+e]
				   +  in2[0][+xo2+e] +  in2[1][+xo2+e] +  in2[2][+xo2+e]
				                              +  in2[3][+xo2+e] +  in2[4][+xo2+e] +  in2[5][+xo2+e];

				if ((d1 >= 0 && d2 >=0)
				 || (d1 < 0 && d2 < 0))
					ss++;				/* Sign was the same */
				
				tdh += d1/(4.5 * 257) * d1/(4.5 * 257);		/* Scale to 0..255 range */
				tdv += d2/(4.5 * 257) * d2/(4.5 * 257);
			}

		tdmag = tdh + tdv;

		if (tdmag < (32.0 * s->th))
			atdmag += tdmag;	/* Average magnitude over a line */
		else
			atdmag += 32.0 * s->th;

		/* if over threshold */
		/* (Cut long lines up to prevent long lines being */
		/* (thrown away due to attached blobs) */
		if (tdmag >= s->th
		 && (x & (CUT_CHUNKS-1)) != (y & (CUT_CHUNKS-1))) {
			double tt;
			double av;		/* Angle value of current pixel */
			tt = (tdv - tdh)/(tdh + tdv);	/* Partial angle */
			linedv += fabs(tt);
			linedc++;

			if (ss >= (s->depth/2+1))	/* Assume signs are the same if clear majority */
				av = 3.0 + tt;
			else
				av = 1.0 - tt;

			/* Separate the orthogonal elements */
			if (av  >= s->divval && av < (s->divval + 2.0)) {
				if (s->flags & SI_SHOW_DIFFSH)
					out[idx] = (char)255;		/* Red */
				/* Add point to new region */
				/* See if we can add to last region */
				if (s->no_hn > 0 && x == s->hregn[s->no_hn-1].hx)
					s->hregn[s->no_hn-1].hx++;
				else {	/* Add another */
					if (s->no_hn >= (w+1)/2) {
						s->errv = SI_INTERNAL;
						sprintf(s->errm,"Internal, no_hn is too large");
						return 1;
					}
					s->hregn[s->no_hn].lx = x;
					s->hregn[s->no_hn].hx = x+1;
					s->hregn[s->no_hn].p = NULL;
					s->no_hn++;
				}
			} else {
				if (s->flags & SI_SHOW_DIFFSV)
					out[idx+1] = (char)255;		/* Green */
				/* Add point to new region */
				/* See if we can add to last region */
				if (s->no_vn > 0 && x == s->vregn[s->no_vn-1].hx)
					s->vregn[s->no_vn-1].hx++;
				else {	/* Add another */
					if (s->no_vn >= (w+1)/2) {
						s->errv = SI_INTERNAL;
						sprintf(s->errm,"Internal, no_vn is too large");
						return 1;
					}
					s->vregn[s->no_vn].lx = x;
					s->vregn[s->no_vn].hx = x+1;
					s->vregn[s->no_vn].p = NULL;
					s->no_vn++;
				}
			}
		}
	}

	if (linedc != 0) {						/* Adapt divider value to line */
		linedv /= (double)linedc;			/* Compute average over the line */
		linedv = (linedv * linedv);			/* Square to even out linedv vs angle */
		linedv = (1.65 * (linedv - 0.12));	/* Compensate for random offsets */
		s->adivval += linedv;
		s->divc++;
		s->divval = (7.0 * s->divval + linedv)/8.0;		/* Average over 8 lines */
		if (s->divval < 0.0)
			s->divval = 0.0;
		else if (s->divval > 1.0)
			s->divval = 1.0;
		if (s->verb >= 5)
			DBG((dbgo,"linedv = %f, divval = %f\n",linedv,s->divval));
	}

	/* Adjust the threshold */
	atdmag /= (double)atdmagc;				/* compute average magnitude over the line */
	s->th = (s->th * THRD)/(THRN + s->divval);/* Convert threshold to average */
	s->th = ((THAWF * TH) + (THAWP * s->th) + (THAWN * atdmag))/(THAWF + THAWP + THAWN);
	s->th = (s->th * (THRN + s->divval))/THRD;	/* Convert average back to threshold */

	/* Add vertical regions */
	if (add_region(s,s->vrego,s->no_vo,s->vregn,s->no_vn,y-2))
		return 1;

	/* Add horizontal regions */
 	if (add_region(s,s->hrego,s->no_ho,s->hregn,s->no_hn,y-2))
		return 1;

	/* shuffle them along */
	tr = s->vrego;
	s->vrego = s->vregn;	/* move new to old */
	s->vregn = tr;			/* old to new */
	s->no_vo = s->no_vn;
	s->no_vn = 0;

	tr = s->hrego;
	s->hrego = s->hregn;	/* move new to old */
	s->hregn = tr;			/* old to new */
	s->no_ho = s->no_hn;
	s->no_hn = 0;

	return 0;
}

/********************************************************************************/
/* Point list code */

/* allocate a new (empty) points structure */
/* return NULL on error */
static points *
new_points(
scanrd_ *s
) {
	points *ps;
	static int pn = 0;
	if ((ps = (points *) malloc(sizeof(points))) == NULL) {
		s->errv = SI_MALLOC_POINTS;
		sprintf(s->errm,"new_points: malloc failed");
		return NULL;
	}
	ps->mxno = 0;
	ps->no = 0;
	ps->nop = 0;
	ps->r = NULL;
	ps->pn = pn;
	pn++;
	return ps;
}

/* destroy a points structure */
static void
destroy_points(
scanrd_ *s,
points *ps) {
	if (ps->r != NULL)	/* Free any array pointed to */
		free(ps->r);
	free (ps);
}

/* Add another run to a points object */
/* return non-zero on error */
static int
add_run(
scanrd_ *s,
points *ps,
int lx,
int hx,
int y)
	{
	if (ps->no == ps->mxno) {	/* Need some more space */
		ps->mxno = (2 * ps->mxno) + 5;		/* New size */
		if ((ps->r = (run *) realloc(ps->r, sizeof(run) * ps->mxno)) == NULL) {
			s->errv = SI_REALLOC_POINTS;
			sprintf(s->errm,"add_run: realloc failed");
			return 1;
		}
	}
	ps->r[ps->no].lx = lx;
	ps->r[ps->no].hx = hx;
	ps->r[ps->no].y = y;	
	ps->no++;					/* One more run */
	ps->nop += hx - lx;			/* Total of pixels */
	return 0;
}

/* copy src points to dest */
/* Return non-zero on error */
static int
copy_points(
scanrd_ *s,
points *dst,
points *src
) {
	int i;
	for (i = 0; i < src->no; i++) {
		if (add_run(s,dst,src->r[i].lx,src->r[i].hx,src->r[i].y))
			return 1;
	}
	return 0;
}

/********************************************************************************/

/* Add a new region of points to the line points lists */
/* Note that regions are assumed to be non-overlapping x sorted */
/* Return non-zero on error */
static int
add_region(
scanrd_ *s,
region *rego,		/* Old regions */
int no_o,			/* No of old region */
region *regn,		/* New regions */
int no_n,			/* No of new region */
int y				/* Y value */
) {
	int osp,op,np;	/* Old/new pointers */

	osp = 0;
	for (np = 0; np < no_n; np++) {	/* Process all new runs */
		/* Advance start pointer until we get to runs that may touch */
#ifdef DIAGN
		while (osp < no_o && rego[osp].hx < regn[np].lx)
#else
		while (osp < no_o && rego[osp].hx <= regn[np].lx)
#endif
			osp++;
		/* For all old runs that may touch new */
#ifdef DIAGN
		for(op = osp; op < no_o && rego[op].lx <= regn[np].hx; op++) {
#else
		for(op = osp; op < no_o && rego[op].lx < regn[np].hx; op++) {
#endif

#ifdef DIAGN
			if (rego[op].hx >= regn[np].lx  && rego[op].lx <= regn[np].hx) {
#else
			if (rego[op].hx > regn[np].lx  && rego[op].lx < regn[np].hx) {
#endif
				/* Old region touches new */
				if (regn[np].p == NULL) {		/* No group for new yet */
					regn[np].p = rego[op].p;	/* Make part of the same group */
					if (add_run(s, regn[np].p,regn[np].lx,regn[np].hx,y)) /* add new run to group */
						return 1;
				} else if (regn[np].p != rego[op].p) {	/* Touches different group */
					int j;
					points *tp = rego[op].p;			/* Old region to be renamed/merged */
					if (copy_points(s,regn[np].p,tp))	/* Merge old with current new */
						return 1;	/* Error */
					DEL_LINK(s->gdone,tp);				/* Don't need other any more */
					for (j = 0; j < no_o; j++)			/* Fix all references to this group */
						if (rego[j].p == tp)
							rego[j].p = regn[np].p;
					for (j = 0; j < no_n; j++)
						if (regn[j].p == tp)
							regn[j].p = regn[np].p;
					destroy_points(s,tp);
				}
			}
		}
		/* Finished all relevant old runs */
		if (regn[np].p == NULL) {	/* No old touched, so start new group */
			if ((regn[np].p = new_points(s)) == NULL)
				return 1;		/* Error */
			ADD_ITEM_TO_TOP(s->gdone,regn[np].p);	/* Stash it in points list */
			if (add_run(s, regn[np].p,regn[np].lx,regn[np].hx,y))	/* add new run to group */
				return 1;		/* Error */
		}
	}
	return 0;
}

/********************************************************************************/

/* Apply partial perspective to an xy point */
/* (We omit the two offset parameters, since we don't need them) */
void ppersp(scanrd_ *s, double *xx, double *yy, double x, double y, double *ppc) {
	double den;

	/* Offset the partial perspective transform */
	x -= ppc[2];
	y -= ppc[3];

	den = ppc[0] * x + ppc[1] * y + 1.0;

	if (fabs(den) < 1e-6) {
		if (den < 0.0)
			den = -1e-6;
		else
			den = 1e-6;
	}
	*xx = x/den + ppc[2];
	*yy = y/den + ppc[3];
}


/* Apply inverse partial perspective to an xy point */
void invppersp(scanrd_ *s, double *x, double *y, double xx, double yy, double *ppc) {
	double den;

	/* Offset the partial perspective transform */
	xx -= ppc[2];
	yy -= ppc[3];

	den = - ppc[0] * xx - ppc[1] * yy + 1.0;

	if (fabs(den) < 1e-6) {
		if (den < 0.0)
			den = -1e-6;
		else
			den = 1e-6;
	}
	*x = xx/den + ppc[2];
	*y = yy/den + ppc[3];
}

/********************************************************************************/

/* Compute the least squares best line fit for a group */
/* Return non-zero if failed */
static int
points_to_line(
scanrd_ *s,
points *ps) {
	int i,j;
	point *vv;			/* Point vectors */
	int nop = ps->nop;	/* Number of points */
	double sx,sy;		/* Sum */
	double mx,my;		/* Mean */
	double a;			/* Angle, Clockwise from 12o'clock */
	double mw,len;		/* mean width, length */
	double x1,y1,x2,y2;	/* Start/end point of fitted line */

	ps->flag = 0;

	if (nop < MIN_POINTS)			/* Don't bother if too few pixels */
		return 0;

	/* Convert runs to individual points, and compute mean */
	if ((vv = (point *) malloc(sizeof(point) * nop)) == NULL) {
		s->errv = SI_MALLOC_POINT2LINE;
		sprintf(s->errm,"scanrd: points_to_line: malloc failed");
		return 1;
	}

	sx = sy = 0.0;
	for (j = i = 0; i < ps->no; i++) {	/* For all runs */
		int x,y;
		int hx = ps->r[i].hx, lx = ps->r[i].lx;

		y = ps->r[i].y;
		sy += (hx - lx) * y;
		for (x = lx; x < hx; x++, j++) {	/* Convert to points */
			sx += x;
			vv[j].x = x;
			vv[j].y = y;
		}
	}
	mx = sx/(double)nop;	/* Centroid (mean) of points */
	my = sy/(double)nop;

	/* Offset points to centroid */
	for (i=0; i < nop; i++) {
		vv[i].x -= mx;
		vv[i].y -= my;
	}
	
	/* Compute ad and bd, then A, B, C */
	/* From Graphics Gems V, pp 91-97, */
	/* "The Best Least-Squares Line Fit" */
	/* by David Alciatore and Rick Miranda. */
	{
		double ad, bd;		/* a' and b' values */
		double xd, yd;		/* temp x' and y' */
		double A, B;		/* line equation */
		double abn;			/* A & B normalizer */

		xd = yd = bd = 0.0;
		for (i = 0; i < nop; i++) {
			double x, y;

			x = vv[i].x;
			y = vv[i].y;
			xd += x * x;
			yd += y * y;
			bd += x * y;
		}
		ad = xd - yd;

		/* Equation of best fit line is Ax + By = C */
		A = 2 * bd;
		B = -(ad + sqrt(ad * ad + 4.0 * bd * bd));
		/* C = A * mx + B * my; */

		/* Compute angle */
		/* A = abn * cos(a), B = -abn * sin(a) */

		abn = sqrt(A * A + B * B);	/* Normalize A & B */
		if (fabs(abn) < 1e-6) {		/* No dominant direction */
			a = 0.0;
		} else {
			a = acos(A/abn);
		}
		/* Make angle +ve */
		while (a < 0.0) a += M_PI;
	}

	/* Now figure out the bounding box for the line + other stats */
	{
		double s,c;
		double pl,nl;		/* Positive length, negative length */
		s = sin(a);
		c = cos(a);
		for (mw = 0.0, pl = 0.0, nl = 0.0, i = 0; i < nop; i++)
			{
			double npj;		/* Projection onto normal */
			double lpj;		/* Projection onto line */
			npj = -c * vv[i].x + s * vv[i].y;
			if (npj < 0)
				mw -= npj;
			else
				mw += npj;
			lpj = s * vv[i].x + c * vv[i].y;
			if (lpj > pl)
				pl = lpj;
			if (lpj < nl)
				nl = lpj;
			}
		mw = 2.0 * mw/(double)nop;		/* Mean width */
	
		x1 = mx + s * nl;
		y1 = my + c * nl;
		x2 = mx + s * pl;
		y2 = my + c * pl;
		len = pl - nl;
	}

	ps->mx = mx;		/* Mean point */
	ps->my = my;
	ps->a = a;			/* Angle */
	ps->mw = mw;		/* Mean width */
	ps->len = len;		/* Mean length */
	ps->x1 = x1;		/* Start/end point of fitted line */
	ps->y1 = y1;
	ps->x2 = x2;
	ps->y2 = y2;
	ps->flag = F_LINESTATS;	/* Line stats valid */

	/* Compute the Constrained to 90 degrees angle */
	/* We use the adivval to figure out where to split angles */
	/* Split at 0 if adivval == 0.0, split at 45 if adivval == 1.0 */
	if (a >= (M_PI * (1.0 - s->adivval/4.0)))
		ps->ca = a - M_PI;
	else if (a >= (M_PI * (0.5 - s->adivval/4.0)))
		ps->ca = a - M_PI_2;
	else
		ps->ca = a;

	if (s->verb >= 5)
	 	DBG((dbgo,"Angle %f, CA = %f, length = %f, mean width  = %f, Line %f,%f to %f,%f\n",
			DEG(a),DEG(ps->ca),len,mw,x1,y1,x2,y2));
	free(vv);

/* printf("~~stats: mw = %f, len = %f, mw/len = %f, area = %f\n",
	mw, len, mw/len, ((double)nop/(len * (mw + 0.01)))); */
	/* Look at stats to see what lines are acceptable for further processing */
	if ( len >= MIN_LINE_LENGTH
	  && mw/len <= MAX_MWID_TO_LEN
	  && ((double)nop/(len * (mw + 0.01))) >= MIN_POINT_TO_AREA) {
		ps->flag |= F_VALID;	/* Line stats valid to use */
/* printf("~~set valid\n"); */
	}
	return 0;
}

static int
calc_lines(
scanrd_ *s
) {
	points *tp;
	s->noslines = 0;
	s->novlines = 0;
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (points_to_line(s,tp))
			return 1;				/* Error */
		if (tp->flag & F_LINESTATS)	/* Line stats valid */
			s->noslines++;
		if (tp->flag & F_VALID)		/* Valid for angle calcs */
			s->novlines++;

		/* Save orininal raster (non partial perspective corrected) values */
		if (tp->flag & F_VALID) {
			tp->pmx = tp->mx;
			tp->pmy = tp->my;
			tp->px1 = tp->x1;
			tp->py1 = tp->y1;
			tp->px2 = tp->x2;
			tp->py2 = tp->y2;
		}
	END_FOR_ALL_ITEMS(tp);
	return 0;
}

static int show_line(scanrd_ *s, int x1, int y1, int x2, int y2, unsigned long c);

/* Show the edge detected lines */
static int
show_lines(
scanrd_ *s
) {
	points *tp;
	int outw = s->width;
	int outh = s->height;
	/* For SI_SHOW_ROT */
	double cirot,sirot;		/* cos and sin of -irot */
	cirot = cos(-s->irot);
	sirot = sin(-s->irot);

	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if ((s->flags & SI_SHOW_ALL_LINES) || (tp->flag & F_VALID))
			{
			unsigned long col = 0xffffff;		/* default color is white */
			double x1 = tp->px1, y1 = tp->py1, x2 = tp->px2, y2 = tp->py2;
			/* For SI_SHOW_ROT */

			/* Show partial perspective corrected lines */
			if (s->flags & (SI_SHOW_ROT | SI_SHOW_PERS)) {
				invppersp(s, &x1, &y1, x1, y1, s->ppc);
				invppersp(s, &x2, &y2, x2, y2, s->ppc);
				col = 0xffff00;	/* cyan */
			}

			/* Show rotation correction of lines + color coding yellow and red */
			if (s->flags & SI_SHOW_ROT) {
				double tx1, ty1, tx2, ty2;
				double a = tp->a - s->irot;

				tx1 = x1;
				ty1 = y1;
				tx2 = x2;
				ty2 = y2;

				/* Rotate about center of raster */
				x1 = (tx1-outw/2.0) * cirot + (ty1-outh/2.0) * sirot;
				y1 = -(tx1-outw/2.0) * sirot + (ty1-outh/2.0) * cirot;
				x2 = (tx2-outw/2.0) * cirot + (ty2-outh/2.0) * sirot;
				y2 = -(tx2-outw/2.0) * sirot + (ty2-outh/2.0) * cirot;

				x1 += outw/2.0;		/* Rotate about center of raster */
				y1 += outh/2.0;
				x2 += outw/2.0;
				y2 += outh/2.0;
				if ((a >= -0.08 && a <= 0.08) || (a >= (M_PI-0.08) && a <= (M_PI+0.08))
				 || (a >= (M_PI_2-0.08) && a <= (M_PI_2+0.08)))
					col = 0x00ffff;	/* yellow */
				else
					col = 0x0000ff;	/* Red */
			}
			/* Show just lines used for fit improvement in blue */
			if (s->flags & SI_SHOW_IMPL) {
				if (tp->flag & F_IMPROVE)
					col = 0xff4040;	/* blue */
			}
			show_line(s,(int)(x1+0.5),(int)(y1+0.5),(int)(x2+0.5),(int)(y2+0.5),col);
		}
	END_FOR_ALL_ITEMS(tp);
	return 0;
}


/********************************************************************************/

/* Definition of the optimization function handed to powell() */
static double
pfunc(void *ss, double p[]) {
	scanrd_ *s = (scanrd_ *)ss;
	points *tp;
	double aa;		/* Average angle */
	double va, rva;	/* Variance */
	double wt;		/* Total weighting = sum of line lengths */
	double pw;
	double dw;		/* Discrimination width */
	
//printf("~1 %f %f %f %f %f %f\n", p[0],p[1],p[2],p[3],p[4],p[5]);

	/* Correct the perspective of all the edge lines using the parameters */
	/* and compute the mean angle */
	aa = 0.0;		/* Average constrained angle */
	wt = 0.0;		/* Total weighting = sum of line lengths */
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_LONGENOUGH) {
			double a, ca;
			invppersp(s, &tp->x1, &tp->y1, tp->px1, tp->py1, p);
			invppersp(s, &tp->x2, &tp->y2, tp->px2, tp->py2, p);

			/* Compute the angle */
			a = atan2(tp->x2 - tp->x1,tp->y2 - tp->y1);

			/* Make angle +ve */
			while (a < 0.0)
				a += M_PI;

			/* Compute the Constrained to 90 degrees angle */
			/* We use the adivval to figure out where to split angles */
			/* Split at 0 if adivval == 0.0, split at 45 if adivval == 1.0 */
			if (a >= (M_PI * (1.0 - s->adivval/4.0)))
				ca = a - M_PI;
			else if (a >= (M_PI * (0.5 - s->adivval/4.0)))
				ca = a - M_PI_2;
			else
				ca = a;

			tp->a = a;
			tp->ca = ca;

			aa += tp->len * ca;
			wt += tp->len;
		}
	END_FOR_ALL_ITEMS(tp);
	aa /= wt;

	/* Calculate the angle variance */
	va = 0.0;
	tp = s->gdone;
	wt = 0.0;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_LONGENOUGH) {
			double tt;
			tt = tp->ca - aa;
			va += tp->len * tt * tt;
			wt += tp->len;
		}
	END_FOR_ALL_ITEMS(tp);
	va = va/wt;

	/* Calculate the a robust angle variance */
	rva = 0.0;
	wt = 0.0;
	dw = sqrt(va) * 3.1;		/* Allow += 0.5 of a standard deviation */
	if (dw < 0.0001)			/* A perfect chart may have dw of zero */
		dw = 0.0001;
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_LONGENOUGH && fabs(tp->ca - aa) <= dw) {
			double tt;
			tt = tp->ca - aa;
			rva += tp->len * tt * tt;
			wt += tp->len;
		}
	END_FOR_ALL_ITEMS(tp);
	if (wt > 0.0) {
		rva = rva/wt;
		va = rva;
	}

	/* Add some regularization to stop it going crazy */
	pw = 0.0;
	pw += 0.01 * (fabs(p[0]) + fabs(p[1]));
	pw += 0.0001 * (fabs(p[2]/s->width - 0.5) + fabs(p[3]/s->height - 0.5));
	va += pw;

	return va;
}

/* Calculate the partial perspective correction factors */
/* Return non-zero if failed */
static int
calc_perspective(
scanrd_ *s
) {
	points *tp;
	int nl;			/* Number of lines used */
	double ml;		/* Minimum length */
	double pc[4];	/* Perspective factors */
	double ss[4];	/* Initial search distance */
	double rv;		/* Return value */
	int rc = 0;		/* Return code */

	if (s->novlines < MIN_NO_LINES) {
		s->errv = SI_FIND_PERSPECTIVE_FAILED;
		sprintf(s->errm,"Not enough valid lines to compute perspective");
		return 1;
	}

	/* Find the longest line */
	ml = 0.0;
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_VALID) {
			if (tp->len > ml)
				ml = tp->len;
		}
	END_FOR_ALL_ITEMS(tp);
	
	/* Make minimum line length to be included in angle */
	/* calculation 1% of longest line */
	ml *= 0.01;

	/* Mark lines long enough to participate in angle calculation */
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_VALID && tp->len >= ml)
			tp->flag |= F_LONGENOUGH;
	END_FOR_ALL_ITEMS(tp);

	/* Locate the perspective correction factors that minimze the */
	/* variance of the mean angle. */

	pc[0] = 0.0;
	pc[1] = 0.0;
	pc[2] = 0.5 * s->width;
	pc[3] = 0.5 * s->height;

	ss[0] = 0.0001;
	ss[1] = 0.0001;
	ss[2] = 1.0001;
	ss[3] = 1.0001;
	rc = powell(&rv, 4, pc,ss,1e-8,2000,pfunc,s, NULL, NULL);

	if (rc == 0) {
		points *tp;
	
		DBG((dbgo,"Perspective correction factors = %f %f %f %f\n",
		       pc[0],pc[1],pc[2],pc[3]));
		
		s->ppc[0] = pc[0];
		s->ppc[1] = pc[1];
		s->ppc[2] = pc[2];
		s->ppc[3] = pc[3];

		/* Implement the perspective correction */
		tp = s->gdone;
		FOR_ALL_ITEMS(points, tp)
			if (tp->flag & F_LONGENOUGH) {
				double a, ca;
				invppersp(s, &tp->x1, &tp->y1, tp->px1, tp->py1, s->ppc);
				invppersp(s, &tp->x2, &tp->y2, tp->px2, tp->py2, s->ppc);
				tp->mx = 0.5 * (tp->x2 + tp->x1);
				tp->my = 0.5 * (tp->y2 + tp->y1);
				tp->len = sqrt((tp->x2 - tp->x1) * (tp->x2 - tp->x1)
				             + (tp->y2 - tp->y1) * (tp->y2 - tp->y1));
	
				/* Compute the angle */
				a = atan2(tp->x2 - tp->x1,tp->y2 - tp->y1);
	
				/* Make angle +ve */
				while (a < 0.0)
					a += M_PI;
	
				/* Compute the Constrained to 90 degrees angle */
				/* We use the adivval to figure out where to split angles */
				/* Split at 0 if adivval == 0.0, split at 45 if adivval == 1.0 */
				if (a >= (M_PI * (1.0 - s->adivval/4.0)))
					ca = a - M_PI;
				else if (a >= (M_PI * (0.5 - s->adivval/4.0)))
					ca = a - M_PI_2;
				else
					ca = a;
	
				tp->a = a;
				tp->ca = ca;
			}
		END_FOR_ALL_ITEMS(tp);
	}

	return 0;
}

/********************************************************************************/
/* Calculate the image rotation */
/* Return non-zero if failed */
static int
calc_rotation(
scanrd_ *s
) {
	points *tp;
	int nl;			/* Number of lines used */
	double ml;		/* Minimum length */
	double aa;		/* Average angle */
	double sd,dw;	/* Standard deviation, deviation window */
	double wt;		/* Total weighting = sum of line lengths */

	if (s->novlines < MIN_NO_LINES) {
		s->errv = SI_FIND_ROTATION_FAILED;
		sprintf(s->errm,"Not enough valid lines to compute rotation angle");
		return 1;
	}

	/* Find the longest line */
	tp = s->gdone;
	ml = 0.0;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_VALID) {
			if (tp->len > ml)
				ml = tp->len;
		}
	END_FOR_ALL_ITEMS(tp);
	
	/* Make minimum line length to be included in angle */
	/* calculation 1% of longest line */
	ml *= 0.01;

	/* Calculate the mean angle */
	aa = 0.0;
	wt = 0.0;		/* Total weighting = sum of line lengths */
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_VALID && tp->len >= ml) {
			aa += tp->len * tp->ca;
			wt += tp->len;
		}
	END_FOR_ALL_ITEMS(tp);
	aa /= wt;

	if (s->verb >= 2)
		DBG((dbgo,"Mean angle = %f\n",DEG(aa)));

	/* Calculate the angle standard deviation */
	tp = s->gdone;
	sd = 0.0;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_VALID && tp->len >= ml) {
			double tt;
			tt = tp->ca - aa;
			sd += tp->len * tt * tt;
		}
	END_FOR_ALL_ITEMS(tp);

	sd = sqrt(sd/wt);

	if (s->verb >= 2)
		DBG((dbgo,"Standard deviation = %f\n",DEG(sd)));

	/* Now re-compute the angle while rejecting any that fall outside one standard deviation */
	s->irot = 0.0;
	wt = 0.0;					/* Total weighting = sum of line lengths */
	nl = 0;
	dw = sd * SD_WINDOW;		/* Allow += 0.5 of a standard deviation */
	if (dw < 0.01)				/* A perfect chart may have dw of zero */
		dw = 0.01;
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_VALID && tp->len >= ml && fabs(tp->ca - aa) <= dw) {
			s->irot += tp->len * tp->ca;
			wt += tp->len;
			nl++;
		}
	END_FOR_ALL_ITEMS(tp);
	if (nl < (MIN_NO_LINES/2)) {
		s->errv = SI_FIND_ROTATION_FAILED;
		sprintf(s->errm,"%d consistent lines is not enough to compute rotation angle",nl);
		return 1;
	}
	s->irot /= wt;

	if (s->verb >= 2)
		DBG((dbgo,"Robust mean angle = %f from %d lines\n",DEG(s->irot),nl));

	return 0;
}

/********************************************************************************/
/* Coalesce close entries of an edge list */
/* return non-zero on error */
static int
coalesce_elist(
scanrd_ *s,
elist *el,
int close			/* Closeness factor, smaller = coarser */
) {
	double r;		/* Margin for coalescence */
	int i,k;

	if (el->c < 2)		/* Need at least 2 entries */
		return 0;

	r = (el->a[el->c-1].pos - el->a[0].pos)/(double)close;
	for (k = 0, i = 1; i < el->c; i++) {
		if ((el->a[i].pos - el->a[k].pos) <= r) {
			/* Merge the two */
			double lk = el->a[k].len;
			double li = el->a[i].len;
			el->a[k].pos = (el->a[k].pos * lk + el->a[i].pos * li)/(lk + li);
			el->a[k].len = lk + li;
			if (el->a[k].p1 > el->a[i].p1)		/* Track overall start/end points */
				el->a[k].p1 = el->a[i].p1;
			if (el->a[k].p2 < el->a[i].p2)
				el->a[k].p2 = el->a[i].p2;
			continue;
		}
		k++;		/* Inc destination pointer */
		if (k != i)
			el->a[k] = el->a[i];	/* shuffle data down */
	}
	k++;		/* one past last out entry */
	el->c = k;
	return 0;
}

static int invert_elist(scanrd_ *s, elist *dl, elist *sl);
static void debug_elist(scanrd_ *s, elist *el);

/* Make up the x and y edge lists */
/* Return non-zero if failed */
static int
calc_elists(
scanrd_ *s,
int ref			/* 1 if generating reference lists */
) {
	int outw = s->width;
	int outh = s->height;
	points *tp;
	int i,j;
	double cirot,sirot;		/* cos and sin of -irot */
	elist xl, yl;	/* Temporary X and Y edge lists array */
	elist tl;		/* temporary crossing list */

	/* Allocate structures for edge lists */
	if ((xl.a = (epoint *) malloc(sizeof(epoint) * s->novlines)) == NULL) {
		s->errv = SI_MALLOC_ELIST;
		sprintf(s->errm,"scanrd: calc_elist: malloc failed - novlines = %d",s->novlines);
		return 1;
	}
	xl.c = 0;
	if ((yl.a = (epoint *) malloc(sizeof(epoint) * s->novlines)) == NULL) {
		s->errv = SI_MALLOC_ELIST;
		sprintf(s->errm,"scanrd: calc_elist: malloc failed - novlines = %d",s->novlines);
		return 1;
	}
	yl.c = 0;

	/* Put valid lines into one of the two edge list arrays */
	cirot = cos(-s->irot);
	sirot = sin(-s->irot);
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_VALID) {
			/* Rotate the point about 0,0 by angle -irot */
			double x,y,a;
			double mx = tp->mx, my = tp->my;
			
			if (ref) {	/* Rotate about center of raster for reference generation */
				mx -= outw/2.0;		/* Rotate about center of raster */
				my -= outh/2.0;
				x = mx * cirot + my * sirot + outw/2.0;
				y = -mx * sirot + my * cirot + outh/2.0;
			} else {	/* Rotate about 0,0 for matching */
				x = mx * cirot + my * sirot;
				y = -mx * sirot + my * cirot;
			}
			a = tp->a - s->irot;
			if ((a >= -0.08 && a <= 0.08) || (a >= (M_PI-0.08) && a <= (M_PI+0.08))) {
				xl.a[xl.c].pos = x;
				xl.a[xl.c].len = tp->len;
				xl.a[xl.c].p1 = y - tp->len/2.0;
				xl.a[xl.c].p2 = y + tp->len/2.0;
				xl.c++;
			} else if (a >= (M_PI_2-0.08) && a <= (M_PI_2+0.08)) {
				yl.a[yl.c].pos = y;
				yl.a[yl.c].len = tp->len;
				yl.a[yl.c].p1 = x - tp->len/2.0;
				yl.a[yl.c].p2 = x + tp->len/2.0;
				yl.c++;
			}
		}
	END_FOR_ALL_ITEMS(tp);

	/* ~~~~ need to check that lists have a reasonable number of entries ~~~~~ */

	/* now sort the lists */
#define HEAP_COMPARE(A,B)  (A.pos < B.pos)
	HEAPSORT(epoint,xl.a,xl.c);
	HEAPSORT(epoint,yl.a,yl.c);
#undef HEAP_COMPARE
	
	/* Copy the temporary lists to the real lists */
	if ((s->xelist.a = (epoint *) malloc(sizeof(epoint) * xl.c)) == NULL) {
		s->errv = SI_MALLOC_ELIST;
		sprintf(s->errm,"scanrd: calc_elist: malloc failed, xl.c = %d",xl.c);
		return 1;
	}
	s->xelist.c = xl.c;
	for (i=0; i < xl.c; i++)
		s->xelist.a[i] = xl.a[i];
	if ((s->yelist.a = (epoint *) malloc(sizeof(epoint) * yl.c)) == NULL) {
		s->errv = SI_MALLOC_ELIST;
		sprintf(s->errm,"scanrd: calc_elist: malloc failed, yl.c = %d",yl.c);
		return 1;
	}
	s->yelist.c = yl.c;
	for (i=0; i < yl.c; i++)
		s->yelist.a[i] = yl.a[i];

	/* Coalese close entries of the final lists */
	if (coalesce_elist(s, &s->xelist,ELISTCDIST))
		return 1;
	if (coalesce_elist(s, &s->yelist,ELISTCDIST))
		return 1;

	/* Calculate crossing count for lines in the X and y lists */
	if ((tl.a = (epoint *) malloc(sizeof(epoint) * (xl.c > yl.c ? xl.c : yl.c))) == NULL) {
		s->errv = SI_MALLOC_ELIST;
		sprintf(s->errm,"scanrd: calc_elist: malloc failed, xl.c = %d, yl.c = %d",xl.c,yl.c);
		return 1;
	}
	/* X list */
	for (i = 0; i < s->xelist.c; i++) {
		double ppos = s->xelist.a[i].pos;
		double pp,np;		/* Previous and next pos */
		if ((i-1) >= 0)
			pp = (ppos + s->xelist.a[i-1].pos)/2.0;	/* Half distance to next line */
		else
			pp = -1e6;
		if ((i+1) < s->xelist.c)
			np = (ppos + s->xelist.a[i+1].pos)/2.0;	/* Half distance to next line */
		else
			np = 1e6;

		/* For all the lines in the Y list */
		for (tl.c = j = 0; j < yl.c; j++) {
			double pos = yl.a[j].pos;
			double p1 = yl.a[j].p1;
			double p2 = yl.a[j].p2;
			if (p1 <= pp)
				p1 = pp;
			if (p2 >= np)
				p2 = np;
			/* If crosses on this lines X within +-0.5 of line each side */
			if (p1 <= np && p2 >= pp) {
				tl.a[tl.c].pos = pos;
				tl.a[tl.c].len = p2 - p1;
				tl.a[tl.c].p1 = p1;
				tl.a[tl.c].p2 = p2;
				tl.c++;
			}
		}
		/* now coalesce the crossings */
		if (coalesce_elist(s,&tl,200))
			return 1;
		/* Put count in line we're working on */
		s->xelist.a[i].ccount = (double)tl.c;
		pp = ppos;
	}

	/* Y list */
	for (i = 0; i < s->yelist.c; i++) {
		double ppos = s->yelist.a[i].pos;
		double pp,np;		/* Previous and next pos */
		if ((i-1) >= 0)
			pp = (ppos + s->yelist.a[i-1].pos)/2.0;	/* Half distance to next line */
		else
			pp = -1e6;
		if ((i+1) < s->xelist.c)
			np = (ppos + s->yelist.a[i+1].pos)/2.0;	/* Half distance to next line */
		else
			np = 1e6;

		for (tl.c = j = 0; j < xl.c; j++) {
			double pos = xl.a[j].pos;
			double p1 = xl.a[j].p1;
			double p2 = xl.a[j].p2;
			if (p1 <= pp)
				p1 = pp;
			if (p2 >= np)
				p2 = np;
			/* If crosses on this lines Y within +-0.5 of line each side */
			if (p1 <= np && p2 >= pp) {
				tl.a[tl.c].pos = pos;
				tl.a[tl.c].len = p2 - p1;
				tl.a[tl.c].p1 = p1;
				tl.a[tl.c].p2 = p2;
				tl.c++;
			}
		}
		/* now coalesce the crossings */
		if (coalesce_elist(s,&tl,200))
			return 1;
		/* Put count in line we're working on */
		s->yelist.a[i].ccount = (double)tl.c;
		pp = ppos;
	}

	/* Normalize the length and ccount */
		{
		double tlen;	/* Total length maximum */
		double tcmax;	/* Total count maximum */
		for (tlen = tcmax = 0.0, i=0; i < s->xelist.c; i++) {
			if (tlen < s->xelist.a[i].len)
				tlen = s->xelist.a[i].len;
			if (tcmax < s->xelist.a[i].ccount)
				tcmax = s->xelist.a[i].ccount;
		}
		for (i=0; i < s->xelist.c; i++) {
			s->xelist.a[i].len /= tlen;
			s->xelist.a[i].ccount /= tcmax;
		}
		for (tlen = tcmax = 0.0, i=0; i < s->yelist.c; i++) {
			if (tlen < s->yelist.a[i].len)
				tlen = s->yelist.a[i].len;
			if (tcmax < s->yelist.a[i].ccount)
				tcmax = s->yelist.a[i].ccount;
		}
		for (i=0; i < s->yelist.c; i++) {
			s->yelist.a[i].len /= tlen;
			s->yelist.a[i].ccount /= tcmax;
		}
	}

	/* Create the inverted lists for any rotation matching */
	if (invert_elist(s, &s->ixelist, &s->xelist))
		return 1;
	if (invert_elist(s, &s->iyelist, &s->yelist))
		return 1;

	if (s->verb >= 3) {
		DBG((dbgo,"\nxelist:\n"));
		debug_elist(s,&s->xelist);
		DBG((dbgo,"\nixelist:\n"));
		debug_elist(s,&s->ixelist);
		DBG((dbgo,"\nyelist:\n"));
		debug_elist(s,&s->yelist);
		DBG((dbgo,"\niyelist:\n"));
		debug_elist(s,&s->iyelist);
	}

	/* Clean up */
	free(xl.a);
	free(yl.a);
	free(tl.a);
	return 0;
}

/********************************************************************************/
/* Write the elists out to a file */

/* Increment a string counter */
static void
strinc(
char *s
) {
	int i,n,c;	/* Length of string and carry flag */
	n = strlen(s);
	for (c = 1, i = n-1; i >= 0 && c != 0; i--) {
		char sval = ' ';
		if (s[i] == '9') {
			s[i] = '0';
			sval = '1';
			c = 1;
		} else if (s[i] == 'z') {
			s[i] = 'a';
			sval = 'a';
			c = 1;
		} else if (s[i] == 'Z') {
			s[i] = 'A';
			sval = 'A';
			c = 1;
		} else {
			s[i]++;
			c = 0;
		}
		if (i == 0 && c != 0) {
			/* Assume there is some more space */
			for (i = n; i >= 0; i--)
				s[i+1] = s[i];
			s[0] = sval;
			break;
		}
	}
}

/* Write out the match reference information */
/* Return non-zero on error */
static int
write_elists(
scanrd_ *s
) {
	char *fname = s->refname;			/* Path of file to write to */
	FILE *elf;
	int i;

	if ((elf=fopen(fname,"w"))==NULL) {
		s->errv = SI_REF_WRITE_ERR;
		sprintf(s->errm,"write_elists: error opening match reference file '%s'",fname);
		return 1;
	}

	fprintf(elf,"REF_ROTATION %f\n\n",DEG(s->irot));

	fprintf(elf,"XLIST %d\n",s->xelist.c);
	for (i = 0; i < s->xelist.c; i++)
		fprintf(elf,"  %f %f %f\n",s->xelist.a[i].pos, s->xelist.a[i].len, s->xelist.a[i].ccount);
	fprintf(elf,"\n");

	fprintf(elf,"YLIST %d\n",s->yelist.c);
	for (i = 0; i < s->yelist.c; i++)
		fprintf(elf,"  %f %f %f\n",s->yelist.a[i].pos, s->yelist.a[i].len, s->yelist.a[i].ccount);
	fprintf(elf,"\n");

	if ((fclose(elf)) == EOF) {
		s->errv = SI_REF_WRITE_ERR;
		error("write_elists: Unable to close match reference file '%s'\n",fname);
		return 1;
	}
	return 0;
}

/* Read in an elist reference file */
/* return non-zero on error */
/* (~~~ the line counting is rather broken ~~~) */
static int
read_relists(
scanrd_ *s
) {
	char *fname = s->refname;			/* Path of file to read from */
	FILE *elf;
	int i,l = 1;
	int rv;
	char *em;	/* Read error message */

	if ((elf=fopen(fname,"r"))==NULL) {
		s->errv = SI_REF_READ_ERR;
		sprintf(s->errm,"read_elists: error opening match reference file '%s'",fname);
		return 1;
	}

	s->fid[0] = s->fid[1] = 0.0;
	s->fid[2] = s->fid[3] = 0.0;
	s->fid[4] = s->fid[5] = 0.0;
	s->fid[6] = s->fid[7] = 0.0;

	/* BOXES */
	for(;;) {
		if((rv = fscanf(elf,"BOXES %d",&s->nsbox)) == 1) {
			l++;
			break;
		}
		if (rv == EOF) {
			em = "Didn't find BOXES before end of file";
			goto read_error;
		}
		if (rv == 0) {
			while ((rv = getc(elf)) != '\n' && rv != EOF);
			l++;
		}
	}

	/* Allocate structures for boxes */
	if ((s->sboxes = (sbox *) calloc(s->nsbox, sizeof(sbox))) == NULL) {
		s->errv = SI_MALLOC_REFREAD;
		sprintf(s->errm,"read_elist, malloc failed");
		return 1;
	}
	for (i = 0; i < s->nsbox;) {
		char xfix1[20], xfix2[20], yfix1[20],yfix2[20];
		char xfirst[20];
		double ox,oy,w,h,xi,yi;
		char xf[20];
		double x;

		if(fscanf(elf," %19s %19s %19s %19s %19s %lf %lf %lf %lf %lf %lf",xfirst ,xfix1, xfix2, yfix1, yfix2, &w, &h, &ox, &oy, &xi, &yi) != 11) {
			em = "Read of BOX failed";
			goto read_error;
		}
		l++;

		/* If Fiducial. Typically top left, top right, botton right, bottom left. */
		if (xfirst[0] == 'F') {
			s->fid[0] = atof(yfix1);
			s->fid[1] = atof(yfix2);
			s->fid[2] = w;
			s->fid[3] = h;
			s->fid[4] = ox;
			s->fid[5] = oy;
			s->fid[6] = xi;
			s->fid[7] = yi;
			s->fidsize = fabs(s->fid[2] - s->fid[0]) + fabs(s->fid[5] - s->fid[3]);
			s->fidsize /= 80.0;
			s->havefids = 1;

//printf("~1 fiducials %f %f, %f %f %f, %f\n",w, h, ox,oy, xi, yi);
			continue;
		}
		for(;;) {		/* Do Y increment */
			x = ox;
			strcpy(xf,xfix1);
			for(;;) {	/* Do X increment */
				if (i >= s->nsbox) {
					em = "More BOXes that declared";
					goto read_error;
				}
				/* '_' is used as a null string marker for single character single cells */
				if (xf[0] == '_')
					sprintf(s->sboxes[i].name,"%s",yfix1);
				else if (yfix1[0] == '_')
					sprintf(s->sboxes[i].name,"%s",xf);
				else {	/* Y indicates Y name comes first */
					if (xfirst[0] == 'Y')
						sprintf(s->sboxes[i].name,"%s%s",yfix1,xf);
					else	/* X or D */
						sprintf(s->sboxes[i].name,"%s%s",xf,yfix1);
				}
				if (xfirst[0] == 'D')
					s->sboxes[i].diag = 1;	/* Diagnostic box - don't print name or read pixels */
				else
					s->sboxes[i].diag = 0;
				s->sboxes[i].x1 = x;
				s->sboxes[i].y1 = oy;
				s->sboxes[i].x2 = x + w;
				s->sboxes[i].y2 = oy + h;

				/* Misc. init. of new sbox */
				s->sboxes[i].xpt[0] = -1.0;		/* No default expected value */
				
				i++;
				x += xi;
				if (strcmp(xf,xfix2) == 0)
					break;
				strinc(xf);
				}
			if (strcmp(yfix1,yfix2) == 0)
				break;
			oy += yi;
			strinc(yfix1);
		}
	}

	/* BOX_SHRINK */
	for(;;) {
		if((rv = fscanf(elf,"BOX_SHRINK %lf ",&s->rbox_shrink)) == 1) {
			l++;
			break;
		}
		if (rv == EOF) {
			em = "Didn't find BOX_SHRINK before end of file";
			goto read_error;
		}
		if (rv == 0) {
			while ((rv = getc(elf)) != '\n' && rv != EOF); 
			l++;
		}
	}

	/* XLIST */
	for(;;) {
		if((rv = fscanf(elf,"XLIST %d ",&s->rxelist.c)) == 1) {
			l++;
			break;
		}
		if (rv == EOF) {
			em = "Didn't find XLIST before end of file";
			goto read_error;
		}
		if (rv == 0) {
			while ((rv = getc(elf)) != '\n' && rv != EOF);
			l++;
		}
	}
	/* Allocate structures for ref edge lists */
	if ((s->rxelist.a = (epoint *) malloc(sizeof(epoint) * s->rxelist.c)) == NULL) {
		s->errv = SI_MALLOC_REFREAD;
		sprintf(s->errm,"read_elist, malloc failed");
		return 1;
	}
	for (i = 0; i < s->rxelist.c; i++) {
		if (fscanf(elf," %lf %lf %lf ",
		    &s->rxelist.a[i].pos, &s->rxelist.a[i].len, &s->rxelist.a[i].ccount) != 3) {
			em = "Failed to read an XLIST line";
			goto read_error;
		}
		l++;
	}

	/* YLIST */
	for(;;) {
		if ((rv = fscanf(elf,"YLIST %d ",&s->ryelist.c)) == 1) {
			l++;
			break;
		}
		if (rv == EOF) {
			em = "Didn't find YLIST before end of file";
			goto read_error;
		}
		if (rv == 0) {
			while ((rv = getc(elf)) != '\n' && rv != EOF);
			l++;
		}
	}
	if ((s->ryelist.a = (epoint *) malloc(sizeof(epoint) * s->ryelist.c)) == NULL) {
		s->errv = SI_MALLOC_REFREAD;
		sprintf(s->errm,"read_elist, malloc failed");
		return 1;
	}
	for (i = 0; i < s->ryelist.c; i++) {
		if (fscanf(elf," %lf %lf %lf ",
		    &s->ryelist.a[i].pos, &s->ryelist.a[i].len, &s->ryelist.a[i].ccount) != 3)
			{
			em = "Failed to read an YLIST line";
			goto read_error;
		}
		l++;
	}

	/* EXPECTED */
	{
		int j;
		int isxyz = 0;
		int nxpt = 0;
		char csps[20];

		for(;;) {
			if ((rv = fscanf(elf,"EXPECTED %19s %d ",csps, &nxpt)) == 2) {
				l++;
				if (strcmp(csps, "XYZ") == 0) {
					isxyz = 1;
					break;
				} else if (strcmp(csps, "LAB") == 0) {
					isxyz = 0;
					break;
				} else {
					em = "Unknown EXPECTED colorespace";
					goto read_error;
				}
			}
			if (rv == EOF) {
				break;
			}
			if (rv == 0) {
				while ((rv = getc(elf)) != '\n' && rv != EOF);
				l++;
			}
		}
		for (j = 0; j < nxpt; j++) {
			char name[20];
			double val[3];
			if (fscanf(elf," %19s %lf %lf %lf ",
			    name, &val[0], &val[1], &val[2]) != 4)
				{
				em = "Failed to read an EXPECTED line";
				goto read_error;
			}
			l++;
			/* Now locate the matching box */
			for (i = 0; i < s->nsbox; i++) {
				if (strcmp(s->sboxes[i].name, name) == 0) {	/* Found it */
					if (isxyz) {
						XYZ2Lab(s->sboxes[i].xpt, val);
					} else {
						s->sboxes[i].xpt[0] = val[0];
						s->sboxes[i].xpt[1] = val[1];
						s->sboxes[i].xpt[2] = val[2];
					}
					s->xpt = 1;
					break;
				}
			}
			if (i >= s->nsbox) {
				em = "Failed to locate matching sample box in EXPECTED list";
				goto read_error;
			}
		}
	}

	if ((fclose(elf)) == EOF) {
		s->errv = SI_REF_WRITE_ERR;
		error("read_elists: Unable to close match reference file '%s'\n",fname);
		return 1;
	}

	/* Generate length normalization factor */
	{
		double tlen;	/* Total of normalized length */
		for (tlen = 0.0, i=0; i < s->rxelist.c; i++)
			tlen += s->rxelist.a[i].len;
		s->rxelist.lennorm = tlen;
		for (tlen = 0.0, i=0; i < s->ryelist.c; i++)
			tlen += s->ryelist.a[i].len;
		s->ryelist.lennorm = tlen;
	}

	if (s->verb >= 3) {
		DBG((dbgo,"\nrxelist:\n"));
		debug_elist(s, &s->rxelist);
		DBG((dbgo,"\nryelist:\n"));
		debug_elist(s, &s->ryelist);
	}

	return 0;

read_error:;
	s->errv = SI_REF_FORMAT_ERR;
	sprintf(s->errm,"read_relist failed at line %d in file %s: %s\n",l,fname,em);
	return 1;
}

/********************************************************************************/
/* Create an inverted direction elist */
/* return non-zero on error */
static int
invert_elist(
scanrd_ *s,
elist *dl,		/* Destination list */
elist *sl		/* Source list */
) {
	int i,j, rc = sl->c;
	
	*dl = *sl;	/* Copy all the structure elements */

	/* Allocate space in the destination list */
	if ((dl->a = (epoint *) malloc(sizeof(epoint) * rc)) == NULL) {
		s->errv = SI_MALLOC_ELIST;
		sprintf(s->errm,"invert_elist: malloc failed");
		return 1;
	}

	/* Copy the array data and reverse its order */
	for (i = 0, j = rc-1; i < rc; i++,j--) {
		dl->a[j]  = sl->a[i];		/* Copy array element */
		dl->a[j].pos = -dl->a[j].pos;	/* Invert position */
	}
	return 0;
}

/* Print out elist */
static void
debug_elist(
scanrd_ *s,
elist *el
) {
	int i, rc = el->c;
	
	DBG((dbgo,"Elist has %d entries allocated at 0x%p\n",el->c,el->a));
	DBG((dbgo,"lennorm = %f\n",el->lennorm));
	for (i = 0; i < rc; i++)
		DBG((dbgo,"  [%d] = %f %f %f\n",i,el->a[i].pos,el->a[i].len,el->a[i].ccount));
}

/* Free the array data in an elist */
static void
free_elist_array(elist *el) {
	free(el->a);
	el->c = 0;
}

/********************************************************************************/
/* !!!!!!! */
/* NEED TO RESOLVE WHY current code is better in some cases, but */
/* not in others. */

#ifndef NEVER	/* Current code */

/* Compute a correlation between two elists */
static double
elist_correl(
scanrd_ *s,
elist *r,		/* Reference list */
elist *t,		/* Target list */
double off, double scale,		/* Offset and scale of target to ref */
int verb		/* Verbose mode */
) {
	int i, j, rc = r->c;
	double cc = 0.0;		/* Correlation */
	double marg = (r->a[rc-1].pos - r->a[0].pos)/150.0;	/* determines sharpness of pos. match */
	double marg2 = marg * 3.0;		/* Don't contribute anything outside this distance */
	
	for (i = j = 0; i < t->c; i++) {
		int ri;		/* Reference index */
		double dd,d1,d2;		/* Distance to nearest reference */
		double pos = (t->a[i].pos + off) * scale;
		double len = t->a[i].len;
		double cnt = t->a[i].ccount;
		while (pos > r->a[j+1].pos && j < (r->c-2)) j++;
		d1 = fabs(pos - r->a[j].pos);
		d2 = fabs(r->a[j+1].pos - pos);
		if (d1 < d2) {
			dd = d1;
			ri = j;
		} else {
			dd = d2;
			ri = j+1;
		}
		if (dd <= marg2) {	/* If close enough to reference */
			double ccf, rcnt = r->a[ri].ccount;
			double llf, rlen = r->a[ri].len;
			double df = marg/(marg + dd);
			df *= df;
			ccf = 1.0 - (rcnt > cnt ? rcnt-cnt : cnt-rcnt);
			llf = 1.0 - (rlen > len ? rlen-len : len-rlen);
			/* The weighting gives slightly more emphasis on matching long lines */
			cc += (1.0 + rlen) * (df * llf * ccf);
			if (verb) {
				DBG((dbgo,"---- t[%d] %f %f %f    this cc = %f, running total cc = %f\n     r[%d] %f %f %f, df = %f, llf = %f, ccf = %f\n",
				i,pos,len,cnt,df * llf * ccf,cc,j,r->a[ri].pos,r->a[ri].len,rcnt,df, llf, ccf));
			}
		}
	}
	return cc/(r->lennorm + (double)r->c);		/* Normalize */
}

#else	/* New test code */

/* Compute a correlation between two elists */
static double
elist_correl(
scanrd_ *s,
elist *r,		/* Reference list */
elist *t,		/* Target list */
double off, double scale,		/* Offset and scale of target to ref */
int verb		/* Verbose mode */
) {
	int i, rc = r->c;
	double cc = 0.0;		/* Correlation */
	double marg = (r->a[rc-1].pos - r->a[0].pos)/100.0;	/* determines sharpness of pos. match */
	double marg2 = marg * marg;		/* marg squared */
	
//printf("~1 doing elist_correl\n");
	/* For each reference edge */
	for (i = 0; i < rc; i++) {
		int j[3], jj, bj, tc = t->c;
		double dd, pos, bdd;

		/* Find the closest target edge using binary search. */
		for(bdd = 1e6, j[2] = tc-1, j[0] = 0; j[2] > (j[0]+1);) {
			double dist;
			j[1] = (j[2] + j[0])/2;			/* Trial point */
			dist = r->a[i].pos - (t->a[j[1]].pos + off) * scale;
			
//printf("~1 j1 = %d, j1 = %d, j0 = %d, dist = %f\n",j[2], j[1], j[0], dist);
			if (dist > 0) {
				j[0] = j[1];
			} else {
				j[2] = j[1];
			}
		}

		/* Locate best out of 3 remaining points */
		for (jj = 0; jj < 3; jj++) {
			double dist;
			pos = (t->a[j[jj]].pos + off) * scale;
			dist = r->a[i].pos - pos;
			dd = dist * dist;				/* Distance squared */
			if (dd < bdd) {	/* New closest */
				bdd = dd;
				bj = j[jj];
			}
		}

//printf("~1 best j = %d, bdd = %f, marg2 = %f\n",bj,bdd,marg2);
		/* Compute correlation */
		if (bdd < marg2) {		/* Within our margine */
			double df = (marg2 - bdd)/marg2;		/* Distance factor */
			double llf, rlen = r->a[i].len, len = t->a[i].len;
			double ccf, rcnt = r->a[i].ccount, cnt = t->a[i].ccount;
			double tcc;
			llf = 1.0 - (rlen > len ? rlen-len : len-rlen);
			ccf = 1.0 - (rcnt > cnt ? rcnt-cnt : cnt-rcnt);

			/* The weighting gives slightly more emphasis on matching long lines */
			/* Not using crossing count */
			tcc = (1.0 + rlen) * (df * llf);
			cc += tcc;
			if (verb) {
				DBG((dbgo,"---- targ[%d] %f %f %f    this cc = %f, running total cc = %f\n",
				bj,pos,t->a[bj].len,t->a[bj].ccount, tcc,cc));
				DBG((dbgo,"      ref[%d] %f %f %f, df = %f, llf = %f, ccf = %f\n",
				i,r->a[i].pos,r->a[i].len,r->a[i].ccount, df, llf, ccf));
			}
		}

	}
	return cc/(r->lennorm + (double)r->c);		/* Normalize */
}

#endif /* NEVER */

/* Structure to hold data for optimization function */
struct _edatas {
	scanrd_ *s;		/* scanrd object */
	elist *r;		/* Reference list */
	elist *t;		/* Target list */
	int verb;		/* Verbose mode */
	}; typedef struct _edatas edatas;

/* Definition of the optimization function handed to powell() */
static double
efunc(void *edata, double p[]) {
	edatas *e = (edatas *)edata;
	double rv = 2.0 - elist_correl(e->s,e->r,e->t,p[0],p[1],e->verb);
	return rv;
}

/* return non-zero on error */
static int
best_match(
scanrd_ *s,
elist *r,		/* Reference list */
elist *t,		/* Target list */
ematch *rv		/* Return values */
) {
	int r0,r1,rw,t0,t1;
	double rwidth;
	double cc;
	double bcc = 0.0, boff = 0.0, bscale = 0.0;	/* best values */
	
	/* The target has been rotated, and we go through all reasonable */
	/* translations and scales to see if we can match it to the */
	/* reference. */
	r0 = 0;
	r1 = r->c-1;
	rw = r->c/2;	/* Minimum number of target line to match all of reference */
	if (t->c/2 < rw)
		rw = t->c/2;
	rwidth = r->a[r1].pos - r->a[r0].pos;

	for (t0 = 0; t0 < t->c-1; t0++) {
		double off;
		for (t1 = t->c-1; t1 > (t0+rw); t1--) {
			double scale;

			scale = rwidth/(t->a[t1].pos - t->a[t0].pos);
			if (scale < 0.001 || scale > 100.0) {
				break;		/* Don't bother with silly scale factors */
			}

			/* Have to compenate the offset for the scale since it is scaled from 0 */
			off = r->a[r0].pos/scale - t->a[t0].pos;
	        cc = elist_correl(s,r,t,off,scale,0);

			if (s->verb >= 7) {
				DBG((dbgo,"Matching target [%d]-[%d] to ref [%d]-[%d] = %f-%f to %f-%f\n",
		   			t0,t1,r0,r1,t->a[t0].pos,t->a[t1].pos,r->a[r0].pos,r->a[r1].pos));
				DBG((dbgo,"Initial off %f, scale %f, cc = %f\n",off,scale,cc));
			}
			if (cc > 0.20) {	/* Looks promising, try optimizing solution */
				double cp[2];	/* Start point/improved point */
				double rv;		/* Return value */
				int rc;			/* Return code */
				edatas dd;		/* Data structure */
				double ss[2] = { 0.1, 0.1};	/* Initial search distance */
				
				dd.s = s;	/* scanrd object */
				dd.r = r;	/* Reference list */
				dd.t = t;	/* Target list */
				dd.verb = 0;		/* Verbose mode */

				/* Set search start point */
				cp[0] = off;
				cp[1] = scale;
				/* Set search distance */
				ss[0] = (0.01 * rwidth/ELISTCDIST)/scale;	/* Search distance */
				ss[1] = scale * 0.01 * rwidth/ELISTCDIST;

				/* Find minimum */
				rc = powell(&rv, 2,cp,ss,0.0001,400,efunc,&dd, NULL, NULL);

				if (rc == 0								/* Powell converged */
				  && cp[1] > 0.001 && cp[1] < 100.0) {	/* and not ridiculous */
					cc = 2.0 - rv;
					off = cp[0];
					scale = cp[1];
				} 
				/* Else use unoptimsed values */

				if (s->verb >= 7) {
					DBG((dbgo,"After optimizing, off %f, scale %f, cc = %f\n",off,scale,cc));
				}
			}

			if (s->verb >= 7) {
				if (cc > 0.25) {
					DBG((dbgo,"Good correlation::\n"));
					elist_correl(s,r,t,off,scale,1);
				}
			}
			if (s->verb >= 7)
				DBG((dbgo,"offset %f, scale %f cc %f\n", off,scale,cc));
	        if (cc > 0.0 && cc > bcc) {	/* Keep best */
				boff = off;
				bscale = scale;
				bcc = cc;
				if (s->verb >= 7)
					DBG((dbgo,"(New best)\n"));
			}
		}
	}
	if (s->verb >= 7)
		DBG((dbgo,"Returning best offset %f, scale %f returns %f\n\n", boff,bscale,bcc));

	/* return best values */
	rv->cc = bcc;
	rv->off = boff;
	rv->scale = bscale;
	return 0;
}

/* Find best offset and scale match between reference and target, */
/* and then from this, compute condidate 90 degree rotations. */
/* Return 0 if got at least one candidate rotation */
/* Return 1 if no reasonable candidate rotation found */
/* Return 2 if some other error */
static int
do_match(
scanrd_ *s
) {
	ematch  xx, yy, xy, yx, xix, yiy, xiy, yix;	/* All 8 matches needed to detect rotations */
	double r0, r90, r180, r270;			/* Correlation for each extra rotation of target */

	/* Check out all the matches */
	if (s->verb >= 2) DBG((dbgo,"Checking xx\n"));
	if (best_match(s, &s->rxelist,&s->xelist,&xx))
		return 2;
	if (s->verb >= 2) DBG((dbgo,"Checking yy\n"));
	if (best_match(s, &s->ryelist,&s->yelist,&yy))
		return 2;
	if (s->verb >= 2) DBG((dbgo,"Checking xy\n"));
	if (best_match(s, &s->rxelist,&s->yelist,&xy))
		return 2;
	if (s->verb >= 2) DBG((dbgo,"Checking yx\n"));
	if (best_match(s, &s->ryelist,&s->xelist,&yx))
		return 2;
	if (s->verb >= 2) DBG((dbgo,"Checking xix\n"));
	if (best_match(s, &s->rxelist,&s->ixelist,&xix))
		return 2;
	if (s->verb >= 2) DBG((dbgo,"Checking yiy\n"));
	if (best_match(s, &s->ryelist,&s->iyelist,&yiy))
		return 2;
	if (s->verb >= 2) DBG((dbgo,"Checking xiy\n"));
	if (best_match(s, &s->rxelist,&s->iyelist,&xiy))
		return 2;
	if (s->verb >= 2) DBG((dbgo,"Checking yix\n"));
	if (best_match(s, &s->ryelist,&s->ixelist,&yix))
		return 2;

	if (s->verb >= 2) {
		DBG((dbgo,"Axis matches for each possible orientation:\n"));
		DBG((dbgo,"  0: xx  = %f, yy  = %f, xx.sc  = %f, yy.sc  = %f\n",
		                                xx.cc,yy.cc,xx.scale,yy.scale));
		DBG((dbgo," 90: xiy = %f, yx  = %f, xiy.sc = %f, yx.sc  = %f\n",
		                                xiy.cc,yx.cc,xiy.scale,yx.scale));
		DBG((dbgo,"180: xix = %f, yiy = %f, xix.sc = %f, yiy.sc = %f\n",
		                                xix.cc,yiy.cc,xix.scale,yiy.scale));
		DBG((dbgo,"270: xy  = %f, yix = %f, xy.sc  = %f, yix.sc = %f\n",
		                                xy.cc,yix.cc,xy.scale,yix.scale));
	}

	/* Compute the combined values for the four orientations. */
	/* add penalty for different scale factors */
	r0   = sqrt(xx.cc * xx.cc + yy.cc * yy.cc)
	     * (xx.scale > yy.scale ? yy.scale/xx.scale : xx.scale/yy.scale);
	r90  = sqrt(xiy.cc * xiy.cc + yx.cc * yx.cc)
	     * (xiy.scale > yx.scale ? yx.scale/xiy.scale : xiy.scale/yx.scale);
	r180 = sqrt(xix.cc * xix.cc + yiy.cc * yiy.cc)
	     * (xix.scale > yiy.scale ? yiy.scale/xix.scale : xix.scale/yiy.scale);
	r270 = sqrt(xy.cc * xy.cc + yix.cc * yix.cc)
	     * (xy.scale > yix.scale ? yix.scale/xy.scale : xy.scale/yix.scale);

	if (s->verb >= 2)
		DBG((dbgo,"r0 = %f, r90 = %f, r180 = %f, r270 = %f\n",r0,r90,r180,r270));

	s->norots = 0;
	if (s->flags & SI_GENERAL_ROT) { /* If general rotation allowed */
		if (s->xpt == 0) {		/* No expected color information to check rotations agaist */
								/* so choose the single best rotation by the edge matching */
			DBG((dbgo,"There is no expected color information, so best fit rotations will be used\n"));
			if (r0 >= MATCHCC && r0 >= r90 && r0 >= r180 && r0 >= r270) {
				s->rots[0].ixoff   = -xx.off; 
				s->rots[0].ixscale = 1.0/xx.scale;
				s->rots[0].iyoff   = -yy.off;
				s->rots[0].iyscale = 1.0/yy.scale;
				s->rots[0].irot    = s->irot;
				s->rots[0].cc      = r0;
				s->norots = 1;
			} else if (r90 >= MATCHCC && r90 >= r180 && r90 >= r270) {
				s->rots[0].ixoff   = -xiy.off;
				s->rots[0].ixscale = 1.0/xiy.scale;
				s->rots[0].iyoff   = -yx.off;
				s->rots[0].iyscale = 1.0/yx.scale;
				s->rots[0].irot    = s->irot + M_PI_2;
				s->rots[0].cc      = r90;
				s->norots = 1;
			} else if (r180 >= MATCHCC && r180 >= r270) {
				s->rots[0].ixoff   = -xix.off;
				s->rots[0].ixscale = 1.0/xix.scale;
				s->rots[0].iyoff   = -yiy.off;
				s->rots[0].iyscale = 1.0/yiy.scale;
				s->rots[0].irot    = s->irot + M_PI;
				s->rots[0].cc      = r180;
				s->norots = 1;
			} else if (r270 >= MATCHCC) {	/* 270 extra target rotation */
				s->rots[0].ixoff   = -xy.off;
				s->rots[0].ixscale = 1.0/xy.scale;
				s->rots[0].iyoff   = -yix.off;
				s->rots[0].iyscale = 1.0/yix.scale;
				s->rots[0].irot    = s->irot + M_PI + M_PI_2;
				s->rots[0].cc      = r270;
				s->norots = 1;
			}

		} else {	/* Got expected color info, so try reasonable rotations */
			double bcc;		/* Best correlation coeff */

			if (r0 >= r90 && r0 >= r180 && r0 >= r270)
				bcc = r0;
			else if (r90 >= r180 && r90 >= r270)
				bcc = r90;
			else if (r180 >= r270)
				bcc = r180;
			else 
				bcc = r270;

			bcc *= ALT_ROT_TH;		/* Threshold for allowing alternate rotation */
			if (bcc < MATCHCC)
				bcc = MATCHCC;

			s->norots = 0;
			if (r0 >= bcc) {
				s->rots[s->norots].ixoff   = -xx.off; 
				s->rots[s->norots].ixscale = 1.0/xx.scale;
				s->rots[s->norots].iyoff   = -yy.off;
				s->rots[s->norots].iyscale = 1.0/yy.scale;
				s->rots[s->norots].irot    = s->irot;
				s->rots[s->norots].cc      = r0;
				s->norots++;
			}
			if (r90 >= bcc) {
				s->rots[s->norots].ixoff   = -xiy.off;
				s->rots[s->norots].ixscale = 1.0/xiy.scale;
				s->rots[s->norots].iyoff   = -yx.off;
				s->rots[s->norots].iyscale = 1.0/yx.scale;
				s->rots[s->norots].irot    = s->irot + M_PI_2;
				s->rots[s->norots].cc      = r90;
				s->norots++;
			}
			if (r180 >= bcc) {
				s->rots[s->norots].ixoff   = -xix.off;
				s->rots[s->norots].ixscale = 1.0/xix.scale;
				s->rots[s->norots].iyoff   = -yiy.off;
				s->rots[s->norots].iyscale = 1.0/yiy.scale;
				s->rots[s->norots].irot    = s->irot + M_PI;
				s->rots[s->norots].cc      = r180;
				s->norots++;
			}
			if (r270 >= bcc) {
				s->rots[s->norots].ixoff   = -xy.off;
				s->rots[s->norots].ixscale = 1.0/xy.scale;
				s->rots[s->norots].iyoff   = -yix.off;
				s->rots[s->norots].iyscale = 1.0/yix.scale;
				s->rots[s->norots].irot    = s->irot + M_PI + M_PI_2;
				s->rots[s->norots].cc      = r270;
				s->norots++;
			}
		}
	} else {	/* Use only rotation 0 */
		if (r0 >= MATCHCC) {
			s->rots[0].ixoff   = -xx.off; 
			s->rots[0].ixscale = 1.0/xx.scale;
			s->rots[0].iyoff   = -yy.off;
			s->rots[0].iyscale = 1.0/yy.scale;
			s->rots[0].irot    = s->irot;
			s->rots[0].cc      = r0;
			s->norots = 1;
		} else if (s->flags & SI_ASISIFFAIL) {
			DBG((dbgo, "Recognition failed, reading patches 'as is' (probably incorrect)\n"));
			s->rots[0].ixoff   = 0.0;
			s->rots[0].ixscale = 1.0;
			s->rots[0].iyoff   = 0.0;
			s->rots[0].iyscale = 1.0;
			s->rots[0].irot    = 0.0;
			s->rots[0].cc      = r0;
			s->norots = 1;
		}
	}

	if (s->verb >= 2) {
		int i;
		DBG((dbgo,"There are %d candidate rotations:\n",s->norots));
		
		for (i = 0; i < s->norots; i++) {
			DBG((dbgo,"cc = %f, irot = %f, xoff = %f, yoff = %f, xscale = %f, yscale = %f\n",
		        s->rots[i].cc, DEG(s->rots[i].irot), s->rots[i].ixoff,s->rots[i].iyoff,s->rots[i].ixscale,s->rots[i].iyscale));
		}
	}

	if (s->norots == 0)
		return 1;

	return 0;
}

/********************************************************************************/
/* perspective transformation. */
/* Transform from raster to reference using iptrans[]. */
/* Transform from reference to raster using ptrans[]. */
static void ptrans(double *xx, double *yy, double x, double y, double *ptrans) {
	double den;

	den = ptrans[6] * x + ptrans[7] * y + 1.0;

	if (fabs(den) < 1e-6) {
		if (den < 0.0)
			den = -1e-6;
		else
			den = 1e-6;
	}

	*xx = (ptrans[0] * x + ptrans[1] * y + ptrans[2])/den;
	*yy = (ptrans[3] * x + ptrans[4] * y + ptrans[5])/den;
}

/* Convert perspective transfom parameters to inverse */
/* perspective transform parameters. */
/* Return nz on error */
int invert_ptrans(double *iptrans, double *ptrans) {
	double scale = ptrans[0] * ptrans[4] - ptrans[1] * ptrans[3];

	if (fabs(scale) < 1e-6)
		return 1;

	scale = 1.0/scale;

	iptrans[0] = scale * (ptrans[4] - ptrans[5] * ptrans[7]);
	iptrans[1] = scale * (ptrans[2] * ptrans[7] - ptrans[1]);
	iptrans[2] = scale * (ptrans[1] * ptrans[5] - ptrans[2] * ptrans[4]);

	iptrans[3] = scale * (ptrans[5] * ptrans[6] - ptrans[3]);
	iptrans[4] = scale * (ptrans[0] - ptrans[2] * ptrans[6]);
	iptrans[5] = scale * (ptrans[2] * ptrans[3] - ptrans[0] * ptrans[5]);

	iptrans[6] = scale * (ptrans[3] * ptrans[7] - ptrans[4] * ptrans[6]);
	iptrans[7] = scale * (ptrans[1] * ptrans[6] - ptrans[0] * ptrans[7]);

	return 0;
}


/* Structure to hold data for optimization function */
struct _pdatas {
	scanrd_ *s;		/* scanrd object */
	double *tar;	/* 4 x x,y raster points */
	double *ref;	/* 4 x x,y reference points */
}; typedef struct _pdatas pdatas;

/* Definition of the optimization function handed to powell() */
/* We simply want to match the 4 points from the reference */
/* back to the target raster. */
static double
ptransfunc(void *pdata, double p[]) {
	pdatas *e = (pdatas *)pdata;
	int i;
	double rv = 0.0;

	for (i = 0; i < 8; i += 2) {
		double x, y;

		ptrans(&x, &y, e->ref[i+0], e->ref[i+1], p);

		rv += (e->tar[i+0] - x) * (e->tar[i+0] - x);
		rv += (e->tar[i+1] - y) * (e->tar[i+1] - y);
	}

	return rv;
}

/* Compute a combined perspective transform */
/* given two sets of four reference points. */
/* Return non-zero on error */
static int
calc_ptrans(
scanrd_ *s,
double *tar,	/* 4 x x,y raster points */
double *ref		/* 4 x x,y reference points */
) {
	int i;
	pdatas dd;
	double ss[8];
	double rv;		/* Return value */
	int rc;			/* Return code */

	dd.s = s;
	dd.tar = tar;
	dd.ref = ref;

	s->ptrans[0] = 1.0;
	s->ptrans[1] = 0.0;
	s->ptrans[2] = 0.0;
	s->ptrans[3] = 0.0;
	s->ptrans[4] = 1.0;
	s->ptrans[5] = 0.0;
	s->ptrans[6] = 0.0;
	s->ptrans[7] = 0.0;

	for (i = 0; i < 8; i++)
		ss[i] = 0.0001;

	rc = powell(&rv, 8, s->ptrans, ss, 1e-7, 500, ptransfunc, &dd, NULL, NULL);

	return rc;
}

/* Compute combined transformation matrix */
/* for the current partial perspective, current */
/* rotation, scale and offsets. */
/* Return non-zero on error */
static int
compute_ptrans(
scanrd_ *s
) {
	double cirot,sirot;		/* cos and sin of -irot */
	double t[6];
	double minx, miny, maxx, maxy;
	double tar[8];
	double ref[8];
	int rv;
	int i;

	/* Compute the rotation and translation part of the */ 
	/* reference to raster target transformation */
	/* xo = t[0] + xi * t[1] + yi * t[2]; */
	/* yo = t[3] + xi * t[4] + yi * t[5]; */
	cirot = cos(s->rots[s->crot].irot);
	sirot = sin(s->rots[s->crot].irot);
	t[0] = cirot * s->rots[s->crot].ixoff + sirot * s->rots[s->crot].iyoff;
	t[1] = s->rots[s->crot].ixscale * cirot;
	t[2] = s->rots[s->crot].iyscale * sirot;

	t[3] = -sirot * s->rots[s->crot].ixoff + cirot * s->rots[s->crot].iyoff;
	t[4] = s->rots[s->crot].ixscale * -sirot;
	t[5] = s->rots[s->crot].iyscale * cirot;

	/* Setup four reference points, and the target raster equivalent. */
	/* Choose min/max of matching boxes as test points, to scale with raster size. */
	minx = miny = 1e60;
	maxx = maxy = -1e60;
	for (i = 0; i < s->nsbox; i++) {
		if (s->sboxes[i].x1 < minx)
			minx = s->sboxes[i].x1;
		if (s->sboxes[i].x2 > maxx)
			maxx = s->sboxes[i].x2;
		if (s->sboxes[i].y1 < miny)
			miny = s->sboxes[i].y1;
		if (s->sboxes[i].y2 > maxy)
			maxy = s->sboxes[i].y2;
	}
	ref[0] = minx;
	ref[1] = miny;
	ref[2] = maxx;
	ref[3] = miny;
	ref[4] = maxx;
	ref[5] = maxy;
	ref[6] = minx;
	ref[7] = maxy;

	for (i = 0; i < 8; i += 2) {
		double x, y;

		x = t[0] + ref[i + 0] * t[1] + ref[i+1] * t[2];
		y = t[3] + ref[i + 0] * t[4] + ref[i+1] * t[5];
		ppersp(s, &x, &y, x, y, s->ppc);
		tar[i + 0] = x;
		tar[i + 1] = y;
	}

	/* Fit the general perspective transform to the points */
	rv = calc_ptrans(s, tar, ref);
	if (rv == 0)
		rv = invert_ptrans(s->iptrans, s->ptrans);

	return rv;
}

/* Compute combined transformation matrix */
/* for the manual alignment case, using fiducial marks. */
/* Return non-zero on error */
static int
compute_man_ptrans(
scanrd_ *s,
double *sfid		/* X & Y of the four target raster marks */
) {
	int rv;

	/* Fit the general perspective transform to the points */
	rv = calc_ptrans(s, sfid, s->fid);
	if (rv == 0)
		rv = invert_ptrans(s->iptrans, s->ptrans);

	return rv;
}

/********************************************************************************/
/* Improve the chosen ptrans to give optimal matching of the */
/* orthogonal edges and the reference edge lists. */

/* Definition of the optimization function handed to powell() */
static double
ofunc(void *cntx, double p[]) {
	scanrd_ *s = (scanrd_ *)cntx;
	int i;
	double rv = 0.0;

	/* First the X list */
	for (i = 0; i < s->rxelist.c; i++) {
		points *tp;

		if (s->rxelist.a[i].nopt == 0)
			continue;

		/* For all the edge lines associated with this tick line */
		for (tp = s->rxelist.a[i].opt; tp != NULL; tp = tp->opt) {
			double x1, y1, x2, y2;
			double d1, d2;

			/* Convert from raster to reference coordinates */
			ptrans(&x1, &y1, tp->px1, tp->py1, p);
			ptrans(&x2, &y2, tp->px2, tp->py2, p);

			d1 = s->rxelist.a[i].pos - x1;
			d2 = s->rxelist.a[i].pos - x2;
			rv += tp->len * (d1 * d1 + d2 * d2); 
		}
	}

	/* Then the Y list */
	for (i = 0; i < s->ryelist.c; i++) {
		points *tp;

		if (s->ryelist.a[i].nopt == 0)
			continue;

		/* For all the edge lines associated with this tick line */
		for (tp = s->ryelist.a[i].opt; tp != NULL; tp = tp->opt) {
			double x1, y1, x2, y2;
			double d1, d2;

			/* Convert from raster to reference coordinates */
			ptrans(&x1, &y1, tp->px1, tp->py1, p);
			ptrans(&x2, &y2, tp->px2, tp->py2, p);

			d1 = s->ryelist.a[i].pos - y1;
			d2 = s->ryelist.a[i].pos - y2;
			rv += tp->len * (d1 * d1 + d2 * d2); 
		}
	}

	return rv;
}

/* optimize the fit of reference ticks to the nearest */
/* edge lines through ptrans[]. */
/* return non-zero on error */
static int
improve_match(
scanrd_ *s
) {
	int i,j;
	points *tp;
	double xspace, yspace;
	int nxlines = 0, nylines = 0;		/* Number of matching lines */

	double pc[8];	/* Parameters to improve */
	double ss[8];	/* Initial search distance */
	double rv;		/* Return value */
	int rc = 0;		/* Return code */

	/* Clear any current elist matching lines */
	for (i = 0; i < s->rxelist.c; i++) {
		s->rxelist.a[i].opt = NULL;
		s->rxelist.a[i].nopt = 0;
	}
	for (i = 0; i < s->ryelist.c; i++) {
		s->ryelist.a[i].opt = NULL;
		s->ryelist.a[i].nopt = 0;
	}

	/* Figure out the average tick spacing for each reference edge list. */
	/* (We're assuming the edge lists are sorted) */
	xspace = (s->rxelist.a[s->rxelist.c-1].pos - s->rxelist.a[0].pos)/s->rxelist.c;
	yspace = (s->ryelist.a[s->ryelist.c-1].pos - s->ryelist.a[0].pos)/s->ryelist.c;

	/* Go through our raster line list, and add the lines that */
	/* closely match the edge list, so that we can fine tune the */
	/* alignment. */
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		if (tp->flag & F_VALID) {
			double x1, y1, x2, y2;
			elist *el;
			double v1, v2;
			double bdist;
			int bix;
			double space;
			int *nlines = NULL;
			double a;

			/* Convert from raster to reference coordinates */
			ptrans(&x1, &y1, tp->px1, tp->py1, s->iptrans);
			ptrans(&x2, &y2, tp->px2, tp->py2, s->iptrans);

			/* Compute the angle */
			a = atan2(y2 - y1,x2 - x1);

			/* Constrain the angle to be between -PI/4 and 3PI/4 */
			if (a < -M_PI_4)
				a += M_PI;
			if (a > M_PI_3_4)
				a -= M_PI;

			/* Decide if it is one of the orthogonal lines */
			if (fabs(a - M_PI_2) > (0.2 * M_PI_2)		/* 0.2 == +/- 18 degrees */	
			 && fabs(a - 0.0)    > (0.2 * M_PI_2)) {
				continue;
			}

			/* Decide which list it would go in */
			if (a > M_PI_4) {
				el = &s->rxelist;
				v1 = x1;
				v2 = x2;
				space = xspace;
				nlines = &nxlines;
			} else {
				el = &s->ryelist;
				v1 = y1;
				v2 = y2;
				space = yspace;
				nlines = &nylines;
			}

			/* Decide which tick it is closest to */
			bdist = 1e38;
			bix = -1;
			for (i = 0; i < el->c; i++) {
				double d1, d2;
				d1 = fabs(el->a[i].pos - v1);
				d2 = fabs(el->a[i].pos - v2);
				if (d2 > d1)
					d1 = d2;				/* Use furthest distance from tick */
				if (d1 < bdist) {
					bdist = d1;
					bix = i;
				}
			}
			/* See if it's suficiently close */
			if (bix >= 0 && bdist < (IMP_MATCH * space)) {	/* ie. 0.1 */
				tp->flag |= F_IMPROVE;
				if (el->a[bix].opt == NULL) {
					(*nlines)++;
				}
				/* Add it to the linked list of matching lines */
				tp->opt = el->a[bix].opt;
				el->a[bix].opt = tp;
				el->a[bix].nopt++;
			}
		}
	END_FOR_ALL_ITEMS(tp);

	if (nxlines < 2 || nylines < 2) {
		if (s->verb >= 1)
			DBG((dbgo,"Improve match failed because there wern't enough close lines\n"));
		return 0;
	}

	/* Optimize iptrans to fit */
	for (i = 0; i < 8; i++) {
		pc[i] = s->iptrans[i];
		ss[i] = 0.0001;
	}

	rc = powell(&rv, 8, pc, ss, 0.0001, 200, ofunc, (void *)s, NULL, NULL);

	if (rc == 0) {
		for (i = 0; i < 8; i++)
			s->iptrans[i] = pc[i];
		rv = invert_ptrans(s->ptrans, s->iptrans);
	}

	return 0;
}

/********************************************************************************/
/* Simple clip to avoid gross problems */
static void clip_ipoint(scanrd_ *s, ipoint *p) {
	int ow = s->width, oh = s->height;

	if (p->x < 0)
		p->x = 0;
	if (p->x >= ow)
		p->x = ow-1;
	if (p->y < 0)
		p->y = 0;
	if (p->y >= oh)
		p->y = oh-1;
}

/* Initialise the sample boxes read for a rescan of the input file */
static int
setup_sboxes(
scanrd_ *s
) {
	int i,j,e;
	sbox *sp;

	for (sp = &s->sboxes[0]; sp < &s->sboxes[s->nsbox]; sp++) {
		double x, y;
		double xx1 = sp->x1, yy1 = sp->y1, xx2 = sp->x2, yy2 = sp->y2;
		int ymin,ymax;	/* index of min and max by y */
		ipoint *p = sp->p;

		/* Shrink box corners by BOX_SHRINK specification */
		xx1 += s->rbox_shrink;
		yy1 += s->rbox_shrink;
		xx2 -= s->rbox_shrink;
		yy2 -= s->rbox_shrink;
		
		/* Transform box corners from reference to raster. */
		/* Box is defined in clockwise direction. */
		ptrans(&x, &y, xx1, yy1, s->ptrans);
		p[0].x = (int)(0.5 + x);
		p[0].y = (int)(0.5 + y);
		clip_ipoint(s, &p[0]);

		ptrans(&x, &y, xx2, yy1, s->ptrans);
		p[1].x = (int)(0.5 + x);
		p[1].y = (int)(0.5 + y);
		clip_ipoint(s, &p[1]);

		ptrans(&x, &y, xx2, yy2, s->ptrans);
		p[2].x = (int)(0.5 + x);
		p[2].y = (int)(0.5 + y);
		clip_ipoint(s, &p[2]);

		ptrans(&x, &y, xx1, yy2, s->ptrans);
		p[3].x = (int)(0.5 + x);
		p[3].y = (int)(0.5 + y);
		clip_ipoint(s, &p[3]);

		if (s->verb >= 4)
			DBG((dbgo,"Box number %ld:\n",(long)(sp - &s->sboxes[0])));

		/* Need to find min/max in y */
		for (i = ymin = ymax = 0; i < 4; i++) {
			if (p[i].y < p[ymin].y)
				ymin = i;
			if (p[i].y > p[ymax].y)
				ymax = i;
		}
		sp->ymin = p[ymin].y;
		sp->ymax = p[ymax].y;
		if (s->verb >= 4)
			DBG((dbgo,"Min y index = %d, value = %d, Max y index = %d, value = %d\n",ymin, sp->ymin, ymax,sp->ymax));

		/* create right side vertex list */
		for (i = -1, j = ymin;;) {
			if (i == -1 || p[j].y != p[sp->r.e[i]].y)
				sp->r.e[++i] = j;	/* Write next if first or different y */
			else if (p[j].x > p[sp->r.e[i]].x)
				sp->r.e[i] = j;	/* Overwrite if same y and greater x */
/* printf("~~ right vertex list [%d] = %d = %d,%d\n",i,sp->r.e[i],p[j].x,p[j].y); */
			if (j == ymax) {
				sp->r.e[++i] = -1;		/* mark end */
/* printf("~~ right vertex list [%d] = %d\n",i,sp->r.e[i]); */
				break;
			}
			j = (j != 3 ? j+1 : 0);/* Advance clockwize */
		}
		sp->r.i = -1;		/* Force first init of edge following */

		/* create left side vertex list */
		for (i = -1, j = ymin;;) {
			if (i == -1 || p[j].y != p[sp->l.e[i]].y)
				sp->l.e[++i] = j;	/* Write next if first or different y */
			else if (p[j].x < p[sp->l.e[i]].x)
				sp->l.e[i] = j;	/* Overwrite if same y and lesser x */
/* printf("~~ left vertex list [%d] = %d = %d,%d\n",i,sp->l.e[i],p[j].x,p[j].y); */
			if (j == ymax) {
				sp->l.e[++i] = -1;		/* mark end */
/* printf("~~ left vertex list [%d] = %d\n",i,sp->r.e[i]); */
				break;
			}
			j = (j != 0 ? j-1 : 3);/* Advance anticlock */
		}
		sp->l.i = -1;		/* Force first init of edge following */

		/* Reset sbox flags */
		for (e = 0; e < s->depth; e++)
			sp->P[e] = -2.0;		/* no value result */
		sp->cnt = 0;
		sp->active = 0;		/* Not active */
	}

	/* allocate and initialize two lists of pointers to the sboxes */
	if ((s->sbstart = (sbox **) malloc(sizeof(sbox *) * s->nsbox)) == NULL) {
		s->errv = SI_MALLOC_SETUP_BOXES;
		sprintf(s->errm,"setup_sboxes: malloc failed");
		return 1;
	}
	if ((s->sbend = (sbox **) malloc(sizeof(sbox *) * s->nsbox)) == NULL) {
		s->errv = SI_MALLOC_SETUP_BOXES;
		sprintf(s->errm,"setup_sboxes: malloc failed");
		return 1;
	}
	for (i = 0; i < s->nsbox; i++)
		s->sbstart[i] = s->sbend[i] = &s->sboxes[i];

	/* Sort sbstart by the minimum y coordinate */
#define HEAP_COMPARE(A,B)  (A->ymin < B->ymin)
	HEAPSORT(sbox *,s->sbstart,s->nsbox);
#undef HEAP_COMPARE
	
	/* Sort s->sbend by the maximum y coordinate */
#define HEAP_COMPARE(A,B)  (A->ymax < B->ymax)
	HEAPSORT(sbox *,s->sbend,  s->nsbox);
#undef HEAP_COMPARE

	s->csi = s->cei = 0;			/* Initialise pointers to start/end lists */

	/* Init active list */
	INIT_LIST(s->alist);
	/* (We ignore any boxes that start above the input raster) */

	return 0;
}

/* Generate the next x on an edge */
static int
nextx(
sbox *sp,
escan *es
) {
	ipoint *p = sp->p;
	int i = es->i;		/* Edge list index */
	int i0 = es->e[i], i1 = es->e[i+1];	/* Index into p[] of current end points */
	
/* printf("~~ nextx called with box %d, escan = 0x%x\n",sp - &s->sboxes[0],es); */
/* printf("~~ i = %d, i0 = %d, i1 = %d\n",i,i0,i1); */
	if (i1 == -1) {	/* Trying to go past the end */
		return es->x;
	}

	/* If never inited or hit start of next segment */
	/* Initialize the next segment */
	if (i == -1 || es->y == p[i1].y) {
		int adx, ady;		/* Absolute deltas */

		i = ++es->i;
		i0 = es->e[i];
		i1 = es->e[i+1];
/* printf("~~ Initing segment, i = %d, i0 = %d, i1 = %d\n",i,i0,i1); */
		if (i1 == -1)	/* Trying to go past the end */
			return es->x;
		es->x = p[i0].x;
		es->y = p[i0].y;

		ady = p[i1].y - p[i0].y;
		adx = p[i1].x - p[i0].x;

		if (adx >= 0)	/* Moving to the right */
			es->xi = 1;
		else
			{			/* Else moving left */
			es->xi = -1;
			adx = -adx;
			}
	
		es->k1 = 2 * adx;
		es->k2 = 2 * (adx - ady) - es->k1;
		es->ev = es->k1 - ady;

/* printf("~~ segment inited, e = %d, k1 = %d, k2 = %d, x = %d, y = %d, xi = %d\n",
es->ev,es->k1,es->k2,es->x,es->y,es->xi); */
		return es->x;
	}

	/* Advance to the next pixel */
	es->y++;
	es->ev += es->k1;
	while (es->ev >= 0 && es->x != p[i1].x) {
		es->x += es->xi;
		es->ev += es->k2;
	}

/* printf("~~ X incremented, e = %d, kw = %d, k2 = %d, x = %d, y = %d, xi = %d\n",
es->ev,es->k1,es->k2,es->x,es->y,es->xi); */
	return es->x;
}

/* Scan value raster location adjustment factors */
double svlaf[21] = {
	1.5196014611277792e-282, 2.7480236142217909e+233,
	1.0605092145600194e-153, 6.1448980493370700e+257,
	5.4169069342907624e-067, 1.6214378600835021e+243,
	9.9021015553451791e+261, 2.4564382802669824e-061,
	1.7476228318632302e+243, 2.0638843604377924e+166,
	1.4097588049607089e-308, 7.7791723264397072e-260,
	5.0497657732134584e+223, 2.2838625101985242e+233,
	5.6363154049548268e+188, 1.4007211907555380e-076,
	6.5805333545409010e+281, 1.3944408779614884e+277,
	7.5963657698668595e-153, 8.2856213563396912e+236,
	7.0898553402722982e+159
};

/* Scan the input file and accumulate the pixel values */
/* return non-zero on error */
static int
do_value_scan(
scanrd_ *s
) {
	int y;			/* current y */
	int ox,oy;		/* x and y size */
	int e;
	unsigned char *in;		/* Input pixel buffer (8bpp) */
	unsigned short *in2;	/* Input pixel buffer (16bpp) */
	int binsize;
	double vscale;		/* Value scale for 16bpp values to range 0.0 - 255.0 */
	double svla;		/* Scan value location adhustment */
	sbox *sp;

	ox = s->width;
	oy = s->height;

	if (s->bpp == 8) {
		binsize = 256;
		vscale = 1.0;
	} else {
		binsize = 65536;
		vscale = 1.0/257.0;
	}

	/* Allocate one input line buffers */
	if ((in = malloc(s->tdepth * ox * s->bypp)) == NULL) {
		s->errv = SI_MALLOC_VALUE_SCAN;
		sprintf(s->errm,"do_value_scan: Failed to malloc test output array");
		return 1;
	}
	in2 = (unsigned short *)in;

	/* Compute the adjustment factor for these patches */
	for (svla = 0.0, e = 1; e < (3 * 7); e++)
		svla += svlaf[e];
	svla *= svlaf[0];

	/* Process the tiff file line by line */
	for (y = 0; y < oy; y++) {
		if (s->read_line(s->fdata, y, (char *)in)) {
			s->errv = SI_RAST_READ_ERR;
			sprintf(s->errm,"scanrd: do_value_scan: read_line() returned error");
			return 1;
		}

		/* Update the active list with new boxes*/
		while (s->csi < s->nsbox && s->sbstart[s->csi]->ymin <= y) {
			/* If goes active on this y */
			if (s->sbstart[s->csi]->diag == 0 && s->sbstart[s->csi]->ymin == y) {
				sp = s->sbstart[s->csi];
				if (s->verb >= 4)
					DBG((dbgo,"added box %ld '%s' to the active list\n",(long)(sp - &s->sboxes[0]),sp->name));
				ADD_ITEM_TO_TOP(s->alist,sp);	/* Add it to the active list */
				sp->active = 1;
				sp->ps[0] = calloc(s->tdepth * binsize,sizeof(unsigned long));
				if (sp->ps[0] == NULL)
					error("do_value_scan: Failed to malloc sbox histogram array");
				for (e = 1; e < s->depth; e++)
					sp->ps[e] = sp->ps[e-1] + binsize;
			}
			s->csi++;
		}
		/* Process the line */
		sp = s->alist;
		FOR_ALL_ITEMS(sbox, sp) {
			int x,x1,x2,xx;	
			unsigned char *oo = &s->out[y * ox * 3];		/* Output raster pointer if needed */
			x1 = nextx(sp,&sp->l);		/* next in left edge */
			x2 = nextx(sp,&sp->r);		/* next in right edge */
			if (s->bpp == 8)
				for (x = s->tdepth*x1, xx = 3*x1; x <= s->tdepth*x2; x += s->tdepth, xx +=3) {
					for (e = 0; e < s->depth; e++)
						sp->ps[e][in[x+e]]++;		/* Increment histogram bins */
					if (s->flags & SI_SHOW_SAMPLED_AREA)
						toRGB(oo+xx, in+x, s->depth, s->bpp);
				}
			else
				for (x = s->tdepth*x1, xx = 3*x1; x <= s->tdepth*x2; x += s->tdepth, xx+=3) {
					for (e = 0; e < s->depth; e++)
						sp->ps[e][in2[x+e]]++;		/* Increment histogram bins */
					if (s->flags & SI_SHOW_SAMPLED_AREA)
						toRGB(oo+xx, (unsigned char *)(in2+x), s->depth, s->bpp);
				}
		} END_FOR_ALL_ITEMS(sp);

	 	
		/* Delete finished boxes from the active list */
		while (s->cei < s->nsbox && s->sbend[s->cei]->ymax <= y) {	/* All that finished last line */
			if (s->verb >= 4)
				DBG((dbgo,"cei = %d, sbenc[s->cei]->ymax = %d, y = %d, active = %d\n",
					s->cei,s->sbend[s->cei]->ymax,y,s->sbend[s->cei]->active));

			/* If goes inactive after this y */
			if (s->sbend[s->cei]->active != 0 && s->sbend[s->cei]->ymax == y) {
				int i,j;
				int cnt;
				double P[MXDE];
				sp = s->sbend[s->cei];
				if (s->verb >= 4)
					DBG((dbgo,"deleted box %ld '%s' from the active list\n",(long)(sp - &s->sboxes[0]),sp->name));
				DEL_LINK(s->alist,sp);		/* Remove it from active list */

				/* Compute mean */
				cnt = 0;
				for (e = 0; e < s->depth; e++)
				sp->mP[e] = 0.0;
				for (i = 0; i < binsize; i++) {	/* For all bins */
					cnt += sp->ps[0][i];
					for (e = 0; e < s->depth; e++)
						sp->mP[e] += (double)sp->ps[e][i] * i;
				}
				for (e = 0; e < s->depth; e++)
					sp->mP[e] /= (double) cnt * svla;
				sp->cnt = cnt;

				/* Compute standard deviation */
				for (e = 0; e < s->depth; e++)
					sp->sdP[e] =  0.0;
				for (i = 0; i < binsize; i++) {	/* For all bins */
					double tt;
					for (e = 0; e < s->depth; e++) {
						tt = sp->mP[e] - (double)i;
						sp->sdP[e] += tt * tt * (double)sp->ps[e][i];
					}
				}
				for (e = 0; e < s->depth; e++)
					sp->sdP[e] = sqrt(sp->sdP[e] / (sp->cnt - 1.0));

				/* Compute "robust" mean */
				/* (There are a number of ways to do this. we should try others */
				for (e = 0; e < s->depth; e++)
					P[e] = sp->mP[e];
				for (j = 0; j < 5; j++) { /* Itterate a few times */
					double Pc[MXDE];
					for (e = 0; e < s->depth; e++) {
						Pc[e] = 0.0;
						sp->P[e] = 0.0;
					}
					for (i = 0; i < binsize; i++) {	/* For all bins */
						double tt;

						/* Unweight values away from current mean */
						for (e = 0; e < s->depth; e++) {
							tt = 1.0 + fabs((double)i - P[e]) * vscale;
							Pc[e] += (double)sp->ps[e][i]/(tt * tt);
							sp->P[e] += (double)sp->ps[e][i]/(tt * tt) * i;
						}
					}
					for (e = 0; e < s->depth; e++)
						P[e] = sp->P[e] /= Pc[e];
				}

				/* Scale all the values to be equivalent to 8bpp range */
				for (e = 0; e < s->depth; e++) {
					sp->mP[e]  *= vscale;
					sp->sdP[e] *= vscale;
					sp->P[e]   *= vscale;
				}

				free(sp->ps[0]);		/* Free up histogram array */
				sp->active = 0;
			}
			s->cei++;
		}
	}

	/* Any boxes remaining on active list must hang */
	/* out over the raster, so discard the results. */
	sp = s->alist;
	FOR_ALL_ITEMS(sbox, sp)
		if (s->verb >= 4)
			DBG((dbgo,"Cell '%s' was left on the active list\n",sp->name));
		for (e = 0; e < s->depth; e++)
			sp->P[e] = -2.0;	/* Signal no value */
		free(sp->ps[0]);		/* Free up histogram array */
		sp->active = 0;
	END_FOR_ALL_ITEMS(sp);

	return 0;
}

/********************************************************************************/
/* Deal with checking the correlation of the current candidate rotation */
/* with the expected values. */
/* Return nz on error. */
static int compute_xcc(scanrd_ *s) {
	int i, n;
	double xcc = 0.0;

	if (s->xpt == 0)
		return 0;
	
	for (n = i = 0; i < s->nsbox; i++) {
		int e;
		sbox *sb = &s->sboxes[i];
		double Lab[3];

		/* Copy computed data to this rotations backup. */
		for (e = 0; e < s->depth; e++) {
			sb->rot[s->crot].mP[e]  = sb->mP[e];
			sb->rot[s->crot].sdP[e] = sb->sdP[e];
			sb->rot[s->crot].P[e]   = sb->P[e];
		}
		sb->rot[s->crot].cnt = sb->cnt;

		if (sb->xpt[0] >= 0.0) {	/* Valid reference value */
			/* Compute rough Lab value for value scanned */
			pval2Lab(Lab, sb->P, s->depth);
			
			/* Add delta E squared to correlation */
			for (e = 0; e < 3; e++) {
				double tt = Lab[e] - sb->xpt[e];
				xcc += tt * tt;
			}
			n++;
		}

	}
	xcc /= (double)n;		/* Average delta E squared */

	/* Record the correlation value */
	s->rots[s->crot].xcc = xcc;

	return 0;
}

#ifdef NEVER	/* We rescan after improvement now */
/* restor the chosen rotation to the "current" sample box values */
static int restore_best(scanrd_ *s) {
	int i;

	for (i = 0; i < s->nsbox; i++) {
		int e;
		sbox *sb = &s->sboxes[i];

		/* Restore sample box value data */
		for (e = 0; e < s->depth; e++) {
			sb->mP[e] = sb->rot[s->crot].mP[e];
			sb->sdP[e] = sb->rot[s->crot].sdP[e];
			sb->P[e] = sb->rot[s->crot].P[e];
		}
		sb->cnt = sb->rot[s->crot].cnt;
	}
	return 0;
}
#endif	/* NEVER */

/********************************************************************************/
/* Initialise, ready to read out all the values */
/* Return the total number of values */
static int
scanrd_reset(
scanrd *ps
) {
	scanrd_ *s = (scanrd_ *)ps;	/* Cast public to private */
	int i,j;
	s->next_read = 0;

	/* Count the number of entries */
	for (j = i = 0; i < s->nsbox; i++)
		if (s->sboxes[i].diag == 0)
			j++;
	return j;
}

/* Read the next samples values */
/* return non-zero when no more points */
static int
scanrd_read(
scanrd *ps,
char *id,			/* patch id copied to here */
double *P,			/* Robust mean values */
double *mP,			/* Raw Mean values */
double *sdP,		/* Standard deviation */
int *cnt			/* Return pixel count, may be NULL, could be zero if not scanned */
) {
	scanrd_ *s = (scanrd_ *)ps;	/* Cast public to private */
	sbox *sp;
	int e;

	/* Skip diagnostic boxes */
	while (s->sboxes[s->next_read].diag != 0 && s->next_read < s->nsbox)
		s->next_read++;

	if (s->next_read >= s->nsbox)
		return 1;

	sp = &s->sboxes[s->next_read++];
	if (sp->diag == 0) {
		if (id != NULL)
			strcpy(id, sp->name);
		for (e = 0; e < s->depth; e++) {
			if (P != NULL)
				P[e] = sp->P[e];
			if (mP != NULL)
				mP[e] = sp->mP[e];
			if (sdP != NULL)
				sdP[e] = sp->sdP[e];
		}
		if (cnt != NULL)
			*cnt = sp->cnt;
	}
	return 0;
}

/********************************************************************************/
static int show_string(scanrd_ *s, char *is, double x, double y,
	double w, unsigned long col);

/* show all the fiducial and sample boxes in the diagnostic raster */
/* return non-zero on error */
static int
show_sbox(
scanrd_ *s
) {
	int i;
	int ev = 0;

	for (i = 0; i < s->nsbox; i++) {
		sbox *sp = &s->sboxes[i];
		unsigned long col = 0x00a0ff;	/* Orange */
		double xx1 = sp->x1, yy1 = sp->y1, xx2 = sp->x2, yy2 = sp->y2;
		double x1,y1,x2,y2,x3,y3,x4,y4;

		/* Transform box corners from reference to raster */
		ptrans(&x1, &y1, xx1, yy1, s->ptrans);
		ptrans(&x2, &y2, xx2, yy1, s->ptrans);
		ptrans(&x3, &y3, xx2, yy2, s->ptrans);
		ptrans(&x4, &y4, xx1, yy2, s->ptrans);

		/* Show outlines of all boxes, or just diagnostic boxes */
		if ((s->flags & SI_SHOW_SBOX_OUTLINES) || (sp->diag != 0)) {
			ev |= show_line(s,(int)(x1+0.5),(int)(y1+0.5),(int)(x2+0.5),(int)(y2+0.5),col);
			ev |= show_line(s,(int)(x2+0.5),(int)(y2+0.5),(int)(x3+0.5),(int)(y3+0.5),col);
			ev |= show_line(s,(int)(x3+0.5),(int)(y3+0.5),(int)(x4+0.5),(int)(y4+0.5),col);
			ev |= show_line(s,(int)(x4+0.5),(int)(y4+0.5),(int)(x1+0.5),(int)(y1+0.5),col);
		}

		/* Show sample boxes names */
		if (s->flags & SI_SHOW_SBOX_NAMES) {
			if (sp->diag == 0)	/* If not diagnostic */
				ev |= show_string(s, sp->name,
							(xx1+xx2)/2.0,(yy1+yy2)/2.0,0.8 * (xx2-xx1),col);
		}

		/* Show non-diagnostic boxes area */
		if ((s->flags & SI_SHOW_SBOX_AREAS) && (sp->diag == 0)) {
			ev |= show_line(s,sp->p[0].x,sp->p[0].y,sp->p[1].x,sp->p[1].y,col);
			ev |= show_line(s,sp->p[1].x,sp->p[1].y,sp->p[2].x,sp->p[2].y,col);
			ev |= show_line(s,sp->p[2].x,sp->p[2].y,sp->p[3].x,sp->p[3].y,col);
			ev |= show_line(s,sp->p[3].x,sp->p[3].y,sp->p[0].x,sp->p[0].y,col);
			ev |= show_line(s,sp->p[0].x,sp->p[0].y,sp->p[2].x,sp->p[2].y,col);
			ev |= show_line(s,sp->p[1].x,sp->p[1].y,sp->p[3].x,sp->p[3].y,col);
		}
	}

	if (s->havefids) {
		for (i = 0; i < 4; i++) {
			unsigned long col = 0x0000ff;	/* Red */
			double xx1 = s->fid[i * 2 + 0];
			double yy1 = s->fid[i * 2 + 1];
			double x1,y1,x2,y2, x3,y3,x4,y4;
			double xsz, ysz;


			/* Make corner point the right way */
			if (i == 0) {
				xsz = s->fidsize;
				ysz = s->fidsize;
			} else if (i == 1) {
				xsz = -s->fidsize;
				ysz = s->fidsize;
			} else if (i == 2) {
				xsz = -s->fidsize;
				ysz = -s->fidsize;
			} else {
				xsz = s->fidsize;
				ysz = -s->fidsize;
			}
	
			/* Create an aligned corner at the fiducial point */
			ptrans(&x1, &y1, xx1, yy1, s->ptrans);
			ptrans(&x2, &y2, xx1 + xsz, yy1, s->ptrans);
			ptrans(&x3, &y3, xx1, yy1, s->ptrans);
			ptrans(&x4, &y4, xx1, yy1 + ysz, s->ptrans);

			ev |= show_line(s,(int)(x1+0.5),(int)(y1+0.5),(int)(x2+0.5),(int)(y2+0.5),col);
			ev |= show_line(s,(int)(x3+0.5),(int)(y3+0.5),(int)(x4+0.5),(int)(y4+0.5),col);
		}
	}

	return ev;
}

/********************************************************************************/
/* Add groups to diagnostic output image */

#undef DBG
#define DBG(aaa) fprintf aaa, fflush(dbgo)

static int
show_groups(
scanrd_ *s
) {
	int stride = 3 * s->width;
	unsigned char *base = s->out; 
	points *tp;
	int x,i,k = 0;
	static unsigned char cc[3 * 24] = {	/* Group palet */
		0x00,0xff,0xff,
		0x00,0x80,0x00,
		0xff,0x00,0xff,
		0x00,0x80,0x80,
		0x00,0xff,0x00,
		0x00,0x80,0xff,
		0x00,0x00,0x80,
		0x80,0xff,0x00,
		0x00,0xff,0x80,
		0xff,0x80,0x00,
		0x00,0x00,0xff,
		0xff,0x80,0x80,
		0x80,0x80,0x00,
		0xff,0xff,0x00,
		0x80,0x80,0x80,
		0x80,0xff,0x80,
		0xff,0xff,0x80,
		0x80,0xff,0xff,
		0xff,0x00,0x80,
		0x80,0x00,0xff,
		0x80,0x80,0xff,
		0xff,0x80,0xff,
		0x80,0x00,0x80,
		0xff,0xff,0xff
		};


	i = 0;
	tp = s->gdone;
	FOR_ALL_ITEMS(points, tp)
		int j;
		/* DBG((dbgo,"Done %d has %d runs\n",i,tp->no)); */
		for (j = 0; j < tp->no; j++) {
			int idx = tp->r[j].y * stride;
			/* Expand the run */
			for (x = tp->r[j].lx; x < tp->r[j].hx; x++) {
				int iidx = idx + 3 * x;
				base[iidx] = cc[k];
				base[iidx+1] = cc[k+1];
				base[iidx+2] = cc[k+2];
			}
		}
		k += 3;
		if (k == (24 * 3))
			k = 0;
		i++;
	END_FOR_ALL_ITEMS(tp);

	return 0;
}
/********************************************************************************/
#ifndef AA_LINES
/* Draw a line in the output diagnostic raster */
static int
show_line(
scanrd_ *s,							/* scanrd object */
int x1, int y1, int x2, int y2,		/* line start and end points */
unsigned long c						/* Color */
) {
	unsigned char *base;				/* Raster base of line */
	int pitch = 3 * s->width;			/* Pitch of raster in pixels */
	int ow = s->width, oh = s->height;	/* width and height of raster for clipping */
	int dx, dy;			/* Line deltas */
	int adx, ady;		/* Absolute deltas */

	int e, k1, k2;		/* Error and axial/diagonal error change values */
	int m1,m2;		/* axial/diagonal coordinate change values */

	int ll;				/* Line length */

	/* Do a crude clip */
	if (x1 < 0)
		x1 = 0;
	if (x1 >= ow)
		x1 = ow-1;
	if (x2 < 0)
		x2 = 0;
	if (x2 >= ow)
		x2 = ow-1;
	if (y1 < 0)
		y1 = 0;
	if (y1 >= oh)
		y1 = oh-1;
	if (y2 < 0)
		y2 = 0;
	if (y2 >= oh)
		y2 = oh-1;

	/* calculate the standard constants */
	dx = x2 - x1;
	dy = y2 - y1;

	if(dx  < 0) {
		m1 = -3;		/* x is going backwards */
		adx = -dx;		/* make this absolute */
	} else {
		m1 = 3;			/* x is going forwards */
		adx = dx;
	}
	
	e = 0;
	if(dy < 0) {
		m2 = -pitch;	/* y is going upwards (decreasing) */
		ady = -dy;		/* make this absolute */
		e = -1;			/* make lines retraceable */
	} else {
		m2 = pitch;		/* y is going downwards (increasing) */
		ady = dy;
	}

	/* m1 has been set to x increment, m2 to y increment */

	m2 += m1;			/* make m2 the diagonal address increment */
						/* and m1 the x axial inrement */
	if(adx > ady) {			/* x is driven */
		ll = adx;
		k1 = 2 * ady;
		k2 = 2 * (ady - adx);
		e += k1 - adx;
	} else {
		ll = ady;
		k1 = 2 * adx;
		k2 = 2 * (adx - ady);
		e += k1 - ady;
		m1 = m2 - m1;		/* Make m1 the y increment */
	}

	/* Start pixel of line */
	base = s->out + y1 * pitch + 3 * x1;

	ll++;	/* Draw start and end point */

	while( ll > 0) {
		while(e < 0 && ll > 0) {
			base[0] = c;
			base[1] = c >> 8;
			base[2] = c >> 16;
			base += m1;
			e += k1;
			ll--;
		}
		while(e >= 0 && ll > 0) {
			base[0] = c;
			base[1] = c >> 8;
			base[2] = c >> 16;
			base += m2;
			e += k2;
			ll--;
		}
	}
	return 0;
}
#else /* AA_LINES: Use anti aliased line drawer */

/*  
    AUTHOR:  Kelvin Thompson

    DESCRIPTION:  Code to render an anti-aliased line, from
      "Rendering Anti-Aliased Lines" in _Graphics_Gems_.

      This is derived from the code printed on pages 690-693
      of _Graphics_Gems_.  An overview of the code is on pages
      105-106.
*/

/* macros to access the frame buffer */
#define PIXINC(dx,dy)  ((dy) * pitch + 3 * (dx))
#define PIXADDR(xx,yy) (s->out + PIXINC(xx,yy))

/* fixed-point data types and macros */
typedef int FX;
#define FX_FRACBITS 16  /* bits of fraction in FX format */
#define FX_0        0   /* zero in fixed-point format */
#define FLOAT_TO_FX(flt)  ((FX)((flt)*(1<<FX_FRACBITS)+0.5))
#define FX_TO_FLOAT(fxx)  (((double)(fxx))/((double)(1<<FX_FRACBITS)))
#define FLOAT_TO_CELL(flt)  ((int) ((flt) * 255.0 + 0.5))
#define MAXVAL_CELL         255
#define COVERAGE(fxval) (s->coverage[(fxval) >> s->covershift])

/* Other aa macros */
#define SWAP(a,b)	((a)^=(b), (b)^=(a), (a)^=(b))

/* BLENDING FUNCTION: */
/*  'cover' is coverage -- in the range [0,255] */
/*  'back' is background color -- in the range [0,255] */
/*  'fgnd' is foreground color -- in the range [0,255] */
#define BLEND(cover,fgnd,back) ( \
	(							\
	  ((255-(cover)) * (back))	\
	+ (     (cover)  * (fgnd))	\
	) >> 8						\
)

/* LINE DIRECTION bits and tables */
#define DIR_STEEP  1  /* set when abs(dy) > abs(dx) */
#define DIR_NEGY   2  /* set whey dy < 0 */

/* --------------------- */
int Anti_Init (scanrd_ *s) {
	float line_r;
	float pix_r;
	int covercells;
	int *thiscell;
	double maxdist,nowdist,incdist;
	int tablebits,radbits;
	int tablecells;
	static int tablesize=0;
	double fnear,ffar,fcover;
	double half,invR,invpiRsq,invpi,Rsq;
	double sum_r;
	double inv_log_2;
	int pitch;

	/* init */
	s->coverage = NULL;

	line_r     = 0.717f;	 /* line radius */
	pix_r      = 0.5;	 	/* pixel radius */
	covercells = 128;

	inv_log_2  = 1.0 / log( 2.0 );
	sum_r      = line_r + pix_r;
	tablebits  = (int) ( log((double)covercells) * inv_log_2 + 0.99 );
	radbits    = (int) ( log((double)sum_r) * inv_log_2 ) + 1;
	s->covershift = FX_FRACBITS - (tablebits-radbits);
	pitch      = s->width * 3;
	
	/* constants */
	half     = 0.5;
	invR     = 1.0 / pix_r;
	invpi    = 1.0 / M_PI;
	invpiRsq = invpi * invR * invR;
	Rsq      = pix_r * pix_r;
#define FRACCOVER(d) (half - d*sqrt(Rsq-d*d)*invpiRsq - invpi*asin(d*invR))
	
	/* pixel increment values  */
	s->adj_pixinc[0] = PIXINC(1,0);
	s->adj_pixinc[1] = PIXINC(0,1);
	s->adj_pixinc[2] = PIXINC(1,0);
	s->adj_pixinc[3] = PIXINC(0,-1);

	s->diag_pixinc[0] = PIXINC(1,1);
	s->diag_pixinc[1] = PIXINC(1,1);
	s->diag_pixinc[2] = PIXINC(1,-1);
	s->diag_pixinc[3] = PIXINC(1,-1);

	s->orth_pixinc[0] = PIXINC(0,1);
	s->orth_pixinc[1] = PIXINC(1,0);
	s->orth_pixinc[2] = PIXINC(0,-1);
	s->orth_pixinc[3] = PIXINC(1,0);

	/* allocate table */
	s->Pmax = FLOAT_TO_FX(sum_r);
	s->Pmax >>= s->covershift;
	tablecells = s->Pmax + 2;
	s->Pmax <<= s->covershift;
	
	if ((s->coverage = (FX *) malloc( tablecells * sizeof(int))) == NULL) {
		s->errv = SI_MALLOC_AAINIT;
		sprintf(s->errm,"aa_line init: Failed to malloc internal table");
		return 1;
	}
	tablesize = tablecells;
	
	/* init for fill loops */
	nowdist = 0.0;
	thiscell = s->coverage;
	incdist = sum_r / (double)(tablecells-2);
	
	/* fill fat portion */
	if (pix_r <= line_r) {
		maxdist = line_r - pix_r;
		for (;nowdist <= maxdist; nowdist += incdist,  ++thiscell)
			*thiscell = MAXVAL_CELL;
	} else { /* fill skinny portion */

		/* loop till edge of line, or end of skinny, whichever comes first */
		maxdist = pix_r - line_r;
		if (maxdist > line_r)
			maxdist = line_r;
		for (;nowdist < maxdist;nowdist += incdist,  ++thiscell) {
			fnear = line_r - nowdist;
			ffar = line_r + nowdist;
			fcover = 1.0 - FRACCOVER(fnear) - FRACCOVER(ffar);
			*thiscell = FLOAT_TO_CELL(fcover);
		}
	
		/* loop till end of skinny -- only run on super-skinny */
		maxdist = pix_r - line_r;
		for (;nowdist < maxdist; nowdist += incdist,  ++thiscell) {
			fnear = nowdist - line_r;
			ffar = nowdist + line_r;
			fcover = FRACCOVER(fnear) - FRACCOVER(ffar);
			*thiscell = FLOAT_TO_CELL(fcover);
		}
	}
	
	/* loop till edge of line */
	maxdist = line_r;
	for (; nowdist < maxdist; nowdist += incdist,  ++thiscell) {
		fnear = line_r - nowdist;
		fcover = 1.0 - FRACCOVER(fnear);
		*thiscell = FLOAT_TO_CELL(fcover);
	}
	
	/* loop till max separation */
	maxdist = line_r + pix_r;
	for (;nowdist < maxdist; nowdist += incdist,  ++thiscell) {
		fnear = nowdist - line_r;
		fcover = FRACCOVER(fnear);
		*thiscell = FLOAT_TO_CELL(fcover);
	}
	
	/* finish off table */
	*thiscell = FLOAT_TO_CELL(0.0);
	s->coverage[tablecells-1] = FLOAT_TO_CELL(0.0);

	s->aa_inited = 1;
	return 0;
#undef FRACCOVER
}

/* --------------------------------------------------------- */
/* Draw an anti-aliased line in the output diagnostic raster */
static int
show_line(
scanrd_ *s,							/* scanrd object */
int X1, int Y1, int X2, int Y2,		/* line start and end points */
unsigned long c						/* Color */
) {
	int 	Bvar, 	/* decision variable for Bresenham's */
	    	Bainc,   /* adjacent-increment for 'Bvar' */
	    	Bdinc;   /* diagonal-increment for 'Bvar' */
	FX 		Pmid,  	/* perp distance at Bresenham's pixel */
	   		Pnow,  	/* perp distance at current pixel (ortho loop) */
	   		Painc, 	/* adjacent-increment for 'Pmid' */
	   		Pdinc, 	/* diagonal-increment for 'Pmid' */
	   		Poinc; 	/* orthogonal-increment for 'Pnow'--also equals 'k' */
	double	fPoinc;	/* Float version of Poinc */
	unsigned char *mid_addr,   /* pixel address for Bresenham's pixel */
	     	      *now_addr;   /* pixel address for current pixel */
	int 	addr_ainc,   /* adjacent pixel address offset */
	    	addr_dinc,   /* diagonal pixel address offset */
	    	addr_oinc;   /* orthogonal pixel address offset */
	int dx,dy,dir;    	/* direction and deltas */
	double fslope;		/* slope of line */
	int pitch     = s->width * 3;
	int ow = s->width, oh = s->height;	/* width and height of raster for clipping */
	int c0,c1,c2;		/* Pixel values */

	if (s->aa_inited == 0) {
		if (Anti_Init(s))
			return 1;	/* Error */
	}

	c0 = c & 0xff;
	c1 = (c >> 8) & 0xff;
	c2 = (c >> 16) & 0xff;

	/* Do a crude clip */
	if (X1 < 1)
		X1 = 1;
	if (X1 >= ow-1)
		X1 = ow-2;
	if (X2 < 1)
		X2 = 1;
	if (X2 >= ow-1)
		X2 = ow-2;
	if (Y1 < 1)
		Y1 = 1;
	if (Y1 >= oh-1)
		Y1 = oh-2;
	if (Y2 < 1)
		Y2 = 1;
	if (Y2 >= oh-1)
		Y2 = oh-2;


	/* rearrange ordering to force left-to-right */
	if 	( X1 > X2 )
	  	{ SWAP(X1,X2);  SWAP(Y1,Y2); }
	
	/* init deltas */
	dx = X2 - X1;  /* guaranteed non-negative */
	dy = Y2 - Y1;

	/* Sanity check */
	if (dx == 0.0 && dy == 0.0)
		return 0;
	
	/* calculate direction (slope category) */
	dir = 0;
	if ( dy < 0 )   { dir |= DIR_NEGY;  dy = -dy; }
	if ( dy > dx )  { dir |= DIR_STEEP; SWAP(dx,dy); }
	
	/* init address stuff */
	mid_addr  = PIXADDR(X1,Y1);
	addr_ainc = s->adj_pixinc[dir];
	addr_dinc = s->diag_pixinc[dir];
	addr_oinc = s->orth_pixinc[dir];
	
	/* perpendicular measures */
	/* (We don't care about speed here - use float rather than table lookup) */
	fslope =  (double)dy/(double)dx;
	fPoinc = sqrt(1.0/(1.0 + (fslope * fslope)));
	Poinc = FLOAT_TO_FX(fPoinc);
	Painc = FLOAT_TO_FX(fPoinc * fslope);
	Pdinc = Painc - Poinc;
	Pmid  = FX_0;
	
	/* init Bresenham's */
	Bainc = dy << 1;
	Bdinc = (dy-dx) << 1;
	Bvar = Bainc - dx;
	
	do	{
		int cvg;

	  	/* do middle pixel */
	  	cvg = COVERAGE(abs(Pmid));
	  	mid_addr[0] = BLEND(cvg, c0, mid_addr[0]);
	  	mid_addr[1] = BLEND(cvg, c1, mid_addr[1]);
	  	mid_addr[2] = BLEND(cvg, c2, mid_addr[2]);
	
	  	/* go up orthogonally */
	  	for (
	      	Pnow = Poinc - Pmid,  now_addr = mid_addr + addr_oinc;
	      	Pnow < s->Pmax;
	      	Pnow += Poinc,      now_addr += addr_oinc
	      	) {
	  		cvg = COVERAGE(Pnow);
	  		now_addr[0] = BLEND(cvg, c0, now_addr[0]);
	  		now_addr[1] = BLEND(cvg, c1, now_addr[1]);
	  		now_addr[2] = BLEND(cvg, c2, now_addr[2]);
		}
	
	  	/* go down orthogonally */
	  	for (Pnow = Poinc + Pmid,  now_addr = mid_addr - addr_oinc;
	      	 Pnow < s->Pmax;
	      	 Pnow += Poinc,      now_addr -= addr_oinc
	      	) {
	  		cvg = COVERAGE(Pnow);
	  		now_addr[0] = BLEND(cvg, c0, now_addr[0]);
	  		now_addr[1] = BLEND(cvg, c1, now_addr[1]);
	  		now_addr[2] = BLEND(cvg, c2, now_addr[2]);
		}
	
	  	/* update Bresenham's */
	  	if ( Bvar < 0 ) {
	    	Bvar += Bainc;
	    	mid_addr += addr_ainc;
	    	Pmid += Painc;
    	} else {
	    	Bvar += Bdinc;
	    	mid_addr += addr_dinc;
	    	Pmid += Pdinc;
    	}
	
	  	--dx;
	} while (dx >= 0);
	return 0;
}

#undef PIXINC
#undef PIXADDR
#undef FX_FRACBITS
#undef FX_0
#undef FLOAT_TO_FX
#undef FX_TO_FLOAT
#undef FLOAT_TO_CELL
#undef MAXVAL_CELL
#undef COVERAGE
#undef SWAP
#undef BLEND
#undef DIR_STEEP
#undef DIR_NEGY

#endif /* !AA_LINES */

/********************************************************************************/
/* Diagnostic vector text output routines */

/* 16 segment ASCII from 0x20 to 0x5f */
/*  
           0      1
         ------ ------
        |\10   11    /|
      7 |  \   |  12  | 2
        |    \ |/     |
         --8--- ---9--
        |     /|\     |
      6 |  15  |  13  | 3
        | /    14   \ |
         ------ ------
            5     4
 */

unsigned short vfont[64] = 
	{
	0x0000, 0x0820, 0x0880, 0x4b3c, 0x4bbb, 0xdb99, 0x2d79, 0x1000, /*  !"#$%&' */
	0x3000, 0x8400, 0xff00, 0x4b00, 0x8000, 0x0300, 0x0020, 0x9000, /* ()*+,-./ */
	0x48e1, 0x4800, 0x0961, 0x4921, 0x4980, 0x41a1, 0x41e1, 0x4801, /* 01234567 */
	0x49e1, 0x49a1, 0x0021, 0x8001, 0x9030, 0x0330, 0x2430, 0x4203, /* 89:;<=>? */
	0x417f, 0x03cf, 0x4a3f, 0x00f3, 0x483f, 0x03f3, 0x01c3, 0x02fb, /* @ABCDEFG */
	0x03cc, 0x4833, 0x4863, 0x31c0, 0x00f0, 0x14cc, 0x24cc, 0x00ff, /* HIJKLMNO */
	0x03c7, 0x20ff, 0x23c7, 0x03bb, 0x4803, 0x00fc, 0x90c0, 0xa0cc, /* PQRSTUVW */
	0xb400, 0x5400, 0x9033, 0x00e1, 0x2400, 0x001e, 0xa000, 0x0030  /* XYZ[\]^_ */
	};

static int show_char(scanrd_ *s, char c,	double x, double y,
	double sc, unsigned long col);

/* Print a string to the diagnostic raster with ptrans() */
/* Return non-zero on error */
static int
show_string(
scanrd_ *s,			/* scanrd object */
char *is,			/* Input string */
double x, double y,	/* Center point for string */
double w,			/* Width total for string */
unsigned long col	/* Color value */
) {
	int i,n;
	double uw;	/* String unscaled width */
	double sc;	/* Scale factor */

	if (w < 0.0)
		w = -w;
	n = strlen(is);
	if (n == 0)
		return 0;
	
	/* Total unscaled width of the string */
	uw = (n * 0.8 + (n >= 1 ? (n-1) * 0.3 : 0));
	/* Compute string scale factor */
	sc = w/uw;

	/* adjust starting point for first char */
	x -= sc * uw/2.0;
	y -= sc * 0.5;

	for (i = 0; i < n; i++) {
		if (show_char(s,is[i],x,y,sc,col))
			return 1;
		x += sc * (0.8 + 0.3);
	}
	return 0;
}

static void show_xfm_line(scanrd_ *s, double x1, double y1, double x2, double y2,
	unsigned long col);

/* Write a character to the diagnostic raster with ptrans() */
/* Return non-zero on error */
static int
show_char(
scanrd_ *s,			/* scanrd object */
char c,				/* Input character */
double x, double y,	/* Top left point of character */
double sc,			/* Scale factor */
unsigned long col
) {
	int ci;
	unsigned int cd;

	ci = c - 0x20;
	if (ci < 0 || ci > 0x3f)
		ci = '?' - 0x20;
	cd = vfont[ci];
	/* Display each segment */
	if (cd & 0x0001)
		show_xfm_line(s, x,y,x+sc*0.4,y,col);
	if (cd & 0x0002)
		show_xfm_line(s, x+sc*0.4,y,x+sc*0.8,y,col);
	if (cd & 0x0004)
		show_xfm_line(s, x+sc*0.8,y,x+sc*0.8,y+sc*0.5,col);
	if (cd & 0x0008)
		show_xfm_line(s, x+sc*0.8,y+sc*0.5,x+sc*0.8,y+sc*1.0,col);
	if (cd & 0x0010)
		show_xfm_line(s, x+sc*0.8,y+sc*1.0,x+sc*0.4,y+sc*1.0,col);
	if (cd & 0x0020)
		show_xfm_line(s, x+sc*0.4,y+sc*1.0,x+0.0,y+sc*1.0,col);
	if (cd & 0x0040)
		show_xfm_line(s, x+0.0,y+sc*1.0,x+0.0,y+sc*0.5,col);
	if (cd & 0x0080)
		show_xfm_line(s, x+0.0,y+sc*0.5,x+0.0,y+0.0,col);
	if (cd & 0x0100)
		show_xfm_line(s, x+0.0,y+sc*0.5,x+sc*0.4,y+sc*0.5,col);
	if (cd & 0x0200)
		show_xfm_line(s, x+sc*0.4,y+sc*0.5,x+sc*0.8,y+sc*0.5,col);
	if (cd & 0x0400)
		show_xfm_line(s, x+0.0,y+0.0,x+sc*0.4,y+sc*0.5,col);
	if (cd & 0x0800)
		show_xfm_line(s, x+sc*0.4,y+0.0,x+sc*0.4,y+sc*0.5,col);
	if (cd & 0x1000)
		show_xfm_line(s, x+sc*0.8,y+0.0,x+sc*0.4,y+sc*0.5,col);
	if (cd & 0x2000)
		show_xfm_line(s, x+sc*0.8,y+sc*1.0,x+sc*0.4,y+sc*0.5,col);
	if (cd & 0x4000)
		show_xfm_line(s, x+sc*0.4,y+sc*1.0,x+sc*0.4,y+sc*0.5,col);
	if (cd & 0x8000)
		show_xfm_line(s, x+0.0,y+sc*1.0,x+sc*0.4,y+sc*0.5,col);
	return 0;
}

/* Write transformed line to the diagnostic raster with ptrans() */
static void
show_xfm_line(
scanrd_ *s,
double x1, double y1, double x2, double y2,
unsigned long col
) {
	double xx1,yy1,xx2,yy2;

	ptrans(&xx1, &yy1, x1, y1, s->ptrans);
	ptrans(&xx2, &yy2, x2, y2, s->ptrans);

	show_line(s,(int)(xx1+0.5),(int)(yy1+0.5),(int)(xx2+0.5),(int)(yy2+0.5),col);
}

/********************************************************************************/
/* Transform from the input raster colorspace to the diagnostic raster space */
static void toRGB(
unsigned char *dst,
unsigned char *src,
int depth, int bpp
) {
	if (bpp == 8) {
		if (depth == 3) {
			dst[0] = src[0];		/* Transfer input to output */
			dst[1] = src[1];
			dst[2] = src[2];
		} else if (depth == 4) {	/* Do a crude conversion */
			double cmyk[4];
			int e;
			for (e = 0; e < 4; e++)
				cmyk[e] = src[e]/255.0;
			for (e = 0; e < 3; e++) {
				cmyk[e] = cmyk[e] * 0.7 + 0.3 * cmyk[3];
				if (cmyk[e] < cmyk[3])
					cmyk[e] = cmyk[3];
				dst[e] = 255 - (int)(cmyk[e] * 255.0 + 0.5);
			}
		} else {	/* Hmm */
			dst[0] = 
			dst[1] = 
			dst[2] = src[0];
		}
	} else {
		unsigned short *src2 = (unsigned short *)src;

		if (depth == 3) {
			dst[0] = src2[0]/257;	/* Transfer input to output */
			dst[1] = src2[1]/257;	/* with 16 to 8bpp conversion */
			dst[2] = src2[2]/257;
		} else if (depth == 4) {	/* Do a crude conversion */
			double cmyk[4];
			int e;
			for (e = 0; e < 4; e++)
				cmyk[e] = src2[e]/65535.0;
			for (e = 0; e < 3; e++) {
				cmyk[e] = cmyk[e] * 0.7 + 0.3 * cmyk[3];
				if (cmyk[e] < cmyk[3])
					cmyk[e] = cmyk[3];
				dst[e] = 255 - (int)(cmyk[e] * 255.0 + 0.5);
			}
		} else {	/* Hmm */
			dst[0] = 
			dst[1] = 
			dst[2] = src2[0]/257;
		}
	}
}


/* Convert from XYZ scale 100 to Lab D50 */
static void XYZ2Lab(double *out, double *in) {
	double X = in[0], Y = in[1], Z = in[2];
	double x,y,z,fx,fy,fz;

	x = X/96.42;
	y = Y/100.0;
	z = Z/82.49;

	if (x > 0.008856451586)
		fx = pow(x,1.0/3.0);
	else
		fx = 7.787036979 * x + 16.0/116.0;

	if (y > 0.008856451586)
		fy = pow(y,1.0/3.0);
	else
		fy = 7.787036979 * y + 16.0/116.0;

	if (z > 0.008856451586)
		fz = pow(z,1.0/3.0);
	else
		fz = 7.787036979 * z + 16.0/116.0;

	out[0] = 116.0 * fy - 16.0;
	out[1] = 500.0 * (fx - fy);
	out[2] = 200.0 * (fy - fz);
}

/* Convert from a scanned pixel value to an aproximate Lab value */
static void pval2Lab(double *out, double *in, int depth) {
	double wXYZ[3];
	double XYZ[3];
	int e, j;

	if (depth == 3) {	/* Assume RGB */

		double clrnts[3][3] = {	/* Red, Green & Blue XYZ values */
			{ 0.412414, 0.212642, 0.019325 },
			{ 0.357618, 0.715136, 0.119207 },
			{ 0.180511, 0.072193, 0.950770 }
		};

		wXYZ[0] = 0.950543;		/* Because we're using sRGB primaries */
		wXYZ[1] = 1.0;			/* the white point is D65 */
		wXYZ[2] = 1.089303;

		XYZ[0] = XYZ[1] = XYZ[2] = 0.0;

		for (e = 0; e < 3; e++) {
			double v = in[e]/255.0;
				
			if (v < 0.0)
				v = 0.0;
			else if (v > 1.0)
				v = 1.0;
			if (v <= 0.03928)
				v /= 12.92;
			else
				v = pow((0.055 + v)/1.055, 2.4);		/* Gamma */

			for (j = 0; j < 3; j++)			/* Sum colorant XYZ */
				XYZ[j] += v * clrnts[e][j];
		}

	} else {
		/* We assume a simple screened subtractive filter model, with dot gain */

		double clrnts[4][3] = {	/* CMYK XYZ values */
			{ 0.12, 0.18, 0.48 },
			{ 0.38, 0.19, 0.20 },	
			{ 0.76, 0.81, 0.11 },	
			{ 0.04, 0.04, 0.04 }
		};

		/* start with white */
		XYZ[0] = wXYZ[0] = 0.9642;
		XYZ[1] = wXYZ[1] = 1.0;
		XYZ[2] = wXYZ[2] = 0.8249;

		/* And filter it out for each component */
		for (e = 0; e < 4; e++) {
			double v = in[e]/255.0;
				
			if (v < 0.0)
				v = 0.0;
			else if (v > 1.0)
				v = 1.0;
			v = 1.0 - pow(1.0 - v, 2.2); /* Compute dot gain */

			for (j = 0; j < 3; j++) {
				double fv;

				/* Normalise filtering effect of this colorant */
				fv = clrnts[e][j]/wXYZ[j];

				/* Compute screened filtering effect */
				fv = (1.0 - v) + v * fv;

				/* Apply filter to our current value */
				XYZ[j] *= fv;
			}
		}
	}

	/* Convert to Lab */
	{
		double X = XYZ[0], Y = XYZ[1], Z = XYZ[2];
		double x,y,z,fx,fy,fz;

		x = X/wXYZ[0];
		y = Y/wXYZ[1];
		z = Z/wXYZ[2];

		if (x > 0.008856451586)
			fx = pow(x,1.0/3.0);
		else
			fx = 7.787036979 * x + 16.0/116.0;

		if (y > 0.008856451586)
			fy = pow(y,1.0/3.0);
		else
			fy = 7.787036979 * y + 16.0/116.0;

		if (z > 0.008856451586)
			fz = pow(z,1.0/3.0);
		else
			fz = 7.787036979 * z + 16.0/116.0;

		out[0] = 116.0 * fy - 16.0;
		out[1] = 500.0 * (fx - fy);
		out[2] = 200.0 * (fy - fz);
	}
}

/********************************************************************************/

static int
scanrd_write_diag(scanrd_ *s) {
	int y;
	unsigned char *op;
	int stride = 3 * s->width;

	if ((s->flags & SI_SHOW_FLAGS) == 0 || s->write_line == NULL)
		return 0;

	/* Write out the tiff file */
	for (op = s->out, y = 0; y < s->height; ++y, op += stride) {
		if (s->write_line(s->ddata, y, (char *)op)) {
			s->errv = SI_DIAG_WRITE_ERR;
			sprintf(s->errm,"scanrd: write_line() returned error");
			return 1;
		}
	}
	return 0;
}

