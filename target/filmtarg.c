/* 
 * Argyll Color Correction System
 * Film recorder chart generator module.
 *
 * Author: Neil Okamoto
 * Date:   1/3/2001
 *
 * Copyright 2001, DreamWorks LLC
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 2 :-
 * see the License2.txt file for licencing details.
 */

/* This program generates a set of TIFF images containing color test patches,
 * given the .ti1 file specifying what the colors are.
 *
 * The output is designed to facilitate output on a 35mm motion picture film
 * recorder such as those manufactured by Cineon, Celco, MGI, ARRI, etc. */

/* Description:
 *
 */

/* TTBD:
 *
 *  - eliminate gamma hack -- do it correctly in targen.
 *  - write a proper description
 *  - multi-patch mosaics on a single tiff image
 *  - render text labels into images
 */

#undef DEBUG

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
#include "cgats.h"
#include "randix.h"
#include "tiffio.h"

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

/* A color structure */
/* This can hold all representations simultaniously */
typedef struct {
	int t;			/* Tag */
#define T_W 0x1
#define T_RGB  0x2
#define T_CMYK 0x4
	double gy;
	double r,g,b;
	double c,m,y,k;
	char *id;		/* Id string */
	char loc[5];	/* Location ID string */
} col;


void
generate_tiffs(
char *basename,		/* Base file name for output tiffs */
col *cols,			/* Array of colors to be put on target chart */
int npat,			/* Number of test targets needed */
char *label,		/* Per set label */
int pw, int ph,		/* Patch width and height */
int rand,			/* Randomise flag */
int rstart,			/* Starting index for random */
int verb)			/* Verbose flag */
{
	randix *r = NULL;	/* Random index order object */
	int ix;		    	/* Patch index in target list */
	int i;		    	/* Logical patch in target list */
	char slab[6];		/* Strip label */
	char fname[200];	/* File name */
	TIFF *tiff;					/* TIFF file handle */
	unsigned short* scanline;	/* scanline buffer */
	int x,y;					/* position in TIFF image */

	if (rand)
		r = new_randix(npat, rstart);

	i = 0;
	while (i < npat) {

			if (rand)
				ix = r->next(r);				/* Next pseudo random index */
			else
				ix = i;

			sprintf(slab, "P%04d",i+1);
			strcpy(cols[ix].loc, slab);
			sprintf(fname, "%s.%04d.tiff", basename, i+1);

			if ((tiff = TIFFOpen(fname, "w")) == NULL) {
				fprintf(stderr,"Failed to open output %s\n", fname);
				exit(-1);
			}

			if (!(scanline = (unsigned short*)malloc(3*sizeof(short)*pw))) {
				fprintf(stderr,"Failed to malloc a scanline buffer\n");
				exit(-1);
			}
			
			TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, pw);
			TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, ph);
			TIFFSetField(tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
			TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
			TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 16);
			TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
			TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
			TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW);

			if (verb) {
				printf("writing patch %4d (%4.2f %4.2f %4.2f) to %s\n",
					   ix, cols[ix].r, cols[ix].g, cols[ix].b, fname);
			}

			/* fill scanline */
			for (x=0; x<pw; x++) {
				int xx = 3*x;
				scanline[xx++] = (unsigned short)((cols[ix].r * 65535.0) + 0.5);
				scanline[xx++] = (unsigned short)((cols[ix].g * 65535.0) + 0.5);
				scanline[xx]   = (unsigned short)((cols[ix].b * 65535.0) + 0.5);
			}
			for (y=0; y<ph; y++) {
				if (TIFFWriteScanline(tiff, (tdata_t)scanline, y, 0) < 0) {
					fprintf(stderr,"WriteScanline Failed at line %d\n",y);
					exit (-1);
				}
			}
			TIFFClose(tiff);
			free(scanline);

			i++;
	}

	if (rand)
		r->del(r);

}

/* film size structure */
typedef struct {
	char *name;			/* User name (lower case) */
	int w,h;			/* Width and height in pixels */
} film;

static film filmsizes[] = {
	{ "Acad_Full",   3656, 2664 },
	{ "Acad_Half",   1828, 1332 },
	{ "Cscope_Full", 3656, 3112 },
	{ "Cscope_Half", 1828, 1556 },
	{ "FullAp_Full", 4096, 3112 },
	{ "FullAp_Half", 2048, 1556 },
	{ "Vista_Full",  4096, 6144 },
	{ "Vista_Half",  2048, 3072 },
	{ NULL, 0, 0 }
};

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

void
usage(void) {
	film *ff;
	fprintf(stderr,"Generate Film Target image sequence, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Licensed under the GPL\n");
	fprintf(stderr,"usage: filmtarg [-v] [-r] [-g gamma] [-f format] basename\n");
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -r              Don't randomise patch location\n");
	fprintf(stderr," -g gamma        Gamma correction / linearization\n");
	fprintf(stderr," -f format       Select format from:\n");
	for (ff = filmsizes; ff->name != NULL; ff++)
	fprintf(stderr,"                 \t%s\t[%d x %d]\n", ff->name, ff->w, ff->h);
	fprintf(stderr," basename        Base name for input(.ti1)/output(.ti2)\n");
	exit(1);
}


int
main(argc,argv)
int argc;
char *argv[];
{
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int rand = 1;
	double gamma = 1.0;
	static char inname[200] = { 0 };		/* Input cgats file base name */
	static char outname[200] = { 0 };		/* Output cgats file base name */
	static char tiffname[200] = { 0 };		/* Output postscrip file base name */
	cgats *icg;			/* input cgats structure */
	cgats *ocg;			/* output cgats structure */
	int i;
	int si, fi, wi;		/* sample id index, field index, keyWord index */
	char *label = "Argyll Color Correction System - Development chart";
	film *flm = &filmsizes[1];	/* Default film format = Academy Half */
	col *cols;
	int npat;		/* Number of patches */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	char buf[100];	/* general sprintf buffer */
	int rstart = 0;	/* Random sequence start value */

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	if (argc <= 1)
		usage();

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++)
		{
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')		/* Look for any flags */
			{
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else
				{
				if ((fa+1) < argc)
					{
					if (argv[fa+1][0] != '-')
						{
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
						}
					}
				}

			if (argv[fa][1] == '?')
				usage();

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;
			/* Randomisation */
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R')
				rand = 0;
			/* Gamma */
			else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				fa = nfa;
				if (na == NULL) usage();
				gamma = atof(na);
			}
			/* Image format */
			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F') {
				fa = nfa;
				if (na == NULL) usage();
				for (flm = filmsizes; flm->name != NULL; flm++) {
					if (cistrcmp(na, flm->name) == 0)
						break;
				}
				if (flm->name == NULL)
					usage();
			} else 
				usage();
			}
		else
			break;
		}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(inname,argv[fa]);
	strcat(inname,".ti1");
	strcpy(outname,argv[fa]);
	strcat(outname,".ti2");
	strcpy(tiffname,argv[fa]);

	icg = new_cgats();	/* Create a CGATS structure */
	icg->add_other(icg, "CTI1"); 	/* our special input type is Calibration Target Information 1 */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI1 format file");
	if (icg->ntables < 1)
		error ("Input file doesn't contain at least one table");

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	if ((cols = (col *)malloc(sizeof(col) * npat)) == NULL)
		error("Malloc failed!");

	/* Setup output cgats file */
	ocg = new_cgats();	/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI2"); 	/* our special type is Calibration Target Information 2 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 2",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll filmtarg", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);

	/* Note what instrument this chart is setup for */
	ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", "GretagMacbeth SpectroScanT", NULL);

	/* Copy various parameters through */
	if ((wi = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[wi], NULL);
	
	if ((wi = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[wi], NULL);
	
	if ((wi = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[wi], NULL);

	if ((wi = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
		ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[wi], NULL);

	
	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
	ocg->add_field(ocg, 0, "SAMPLE_LOC", cs_t);

	if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
		error ("Input file doesn't contain field SAMPLE_ID");
	if (icg->t[0].ftype[si] != nqcs_t)
		error ("Field SAMPLE_ID is wrong type");

	/* Figure out the color space */
	if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file doesn't contain keyword COLOR_REPS");
	if (strcmp(icg->t[0].kdata[fi],"CMYK") == 0) {
		error ("We don't support CMYK yet");
	} else if (strcmp(icg->t[0].kdata[fi],"RGB") == 0) {
		int ri, gi, bi;
		if ((ri = icg->find_field(icg, 0, "RGB_R")) < 0)
			error ("Input file doesn't contain field RGB_R");
		if (icg->t[0].ftype[ri] != r_t)
			error ("Field RGB_R is wrong type");
		if ((gi = icg->find_field(icg, 0, "RGB_G")) < 0)
			error ("Input file doesn't contain field RGB_G");
		if (icg->t[0].ftype[gi] != r_t)
			error ("Field RGB_G is wrong type");
		if ((bi = icg->find_field(icg, 0, "RGB_B")) < 0)
			error ("Input file doesn't contain field RGB_B");
		if (icg->t[0].ftype[bi] != r_t)
			error ("Field RGB_B is wrong type");
		ocg->add_field(ocg, 0, "RGB_R", r_t);
		ocg->add_field(ocg, 0, "RGB_G", r_t);
		ocg->add_field(ocg, 0, "RGB_B", r_t);
		ocg->add_kword(ocg, 0, "COLOR_REP","RGB", NULL);
		if (gamma != 1.0) {
			sprintf(buf, "%6.4f", gamma);
			ocg->add_kword(ocg, 0, "GAMMA", buf, NULL);
		}
		for (i = 0; i < npat; i++) {
			cols[i].t = T_RGB;
			cols[i].id = ((char *)icg->t[0].fdata[i][si]);
			/* normalized and possibly gamma corrected */
			cols[i].r = pow(*((double *)icg->t[0].fdata[i][ri]) / 100.0, gamma);
			cols[i].g = pow(*((double *)icg->t[0].fdata[i][gi]) / 100.0, gamma);
			cols[i].b = pow(*((double *)icg->t[0].fdata[i][bi]) / 100.0, gamma);
		}
	} else if (strcmp(icg->t[0].kdata[fi],"W") == 0) {
		error ("We don't support GRAY yet");
	} else
		error ("Input file keyword COLOR_REPS has unknown value");


	if (verb)
		printf("Film format chosen is %s   [%d x %d]\n", flm->name, flm->w, flm->h);

	if (rand) {
		rstart = clk % npat;
		sprintf(buf,"%d",rstart);
		ocg->add_kword(ocg, 0, "RANDOM_START", buf, NULL);
	}

	generate_tiffs(tiffname, cols, npat, label, flm->w, flm->h, rand, rstart, verb);

	/* Write out the patch info to the output CGATS file */
	for (i = 0; i < npat; i++) {
		if (cols[i].t & T_CMYK)
			ocg->add_set(ocg, 0, cols[i].id, cols[i].loc, 100.0 * cols[i].c, 100.0 * cols[i].m,
			                         100.0 * cols[i].y, 100.0 * cols[i].k);
		else if (cols[i].t & T_RGB)
			ocg->add_set(ocg, 0, cols[i].id, cols[i].loc, 100.0 * cols[i].r,
			                         100.0 * cols[i].g, 100.0 * cols[i].b);
		else if (cols[i].t & T_W)
			ocg->add_set(ocg, 0, cols[i].id, cols[i].loc, 100.0 * cols[i].gy);
	}

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	free(cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */

	return 0;
}

/******************************************************************/
/* Error/debug output routines */
/******************************************************************/


/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"filmtarg: Error - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit (-1);
}

void
warning(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"filmtarg: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}







