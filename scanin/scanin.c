
/* 
 * Argyll Color Correction System
 *
 * Scanin: Input the scan of a test chart, and output cgats data
 *         Uses scanrd to do the hard work.
 *
 * Author: Graeme W. Gill
 * Date:   29/1/97
 *
 * Copyright 1995 - 2002 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
	TTBD

	Add "single pixel patch" mode, for pure digital processing for
	abstract profile creation.


 */

#include <stdio.h>
#include <fcntl.h>		/* In case DOS binary stuff is needed */
#include <ctype.h>
#include <string.h>
#include <math.h>

#include <time.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "cgats.h"
#include "icc.h"
#include "xicc.h"
#include "scanrd.h"
#include "tiffio.h"
#include "ui.h"

/* Utilities */

void fix_it8(char *o, char *i);

#ifdef NT		/* You'd think there might be some standards.... */
# ifndef __BORLANDC__
#  define stricmp _stricmp
# endif
#else
# define stricmp strcasecmp
#endif

#define TXBUF (256*1024L)

/* NOTE: We aren't handling the libtiff error/warning messages !! */

/* Read a line of data from the input Grey, RGB or CMYK tiff file */
/* return non-zero on error */
int read_line(
void *fdata,
int y,
char *dst
) {
	if (TIFFReadScanline((TIFF *)fdata, (tdata_t)dst, y, 0) < 0)
		return 1;
	return 0;
}

/* Write a line of data to the diagnostic RGB tiff file */
/* return non-zero on error */
static int
write_line(
void *ddata,
int y,
char *src
) {
	if (TIFFWriteScanline((TIFF *)ddata, (tdata_t)src, y, 0) < 0)
		return 1;
	return 0;
}

/* Compute a simple string hash value */
unsigned int shash(char *c) {
	unsigned int hash = 0;

	for (; *c != '\000'; c++)
		hash = hash * 67 + (unsigned char)(*c);

	return hash;
}

void
usage(void) {
	fprintf(stderr,"Scanin, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin [options] input.tif recogin.cht valin.cie [diag.tif]\n");
	fprintf(stderr,"   :- inputs 'input.tif' and outputs scanner 'input.ti3', or\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin -g [options] input.tif recogout.cht [diag.tif]\n");
	fprintf(stderr,"   :- outputs file 'recogout.cht', or\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin -o [options] input.tif recogin.cht [diag.tif]\n");
	fprintf(stderr,"   :- outputs file 'input.val', or\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin -c [options] input.tif recogin.cht scanprofile.[%s|mpp] pbase [diag.tif]\n",ICC_FILE_EXT_ND);
	fprintf(stderr,"   :- inputs pbase.ti2 and outputs printer pbase.ti3, or\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin -r [options] input.tif recogin.cht pbase [diag.tif]\n");
	fprintf(stderr,"   :- inputs pbase.ti2+.ti3 and outputs pbase.ti3\n");
	fprintf(stderr,"\n");
	fprintf(stderr," -g                   Generate a chart reference (.cht) file\n");
	fprintf(stderr," -o                   Output patch values in .val file\n");
	fprintf(stderr," -c                   Use image to measure color to convert printer pbase .ti2 to .ti3\n");
	fprintf(stderr," -ca                  Same as -c, but accumulates more values to pbase .ti3\n");
	fprintf(stderr,"                       from subsequent pages\n");
	fprintf(stderr," -r                   Replace device values in pbase .ti2/.ti3\n");
	fprintf(stderr,"                      Default is to create a scanner .ti3 file\n");
	fprintf(stderr," -F x1,y1,x2,y2,x3,y3,x4,y4\n");
	fprintf(stderr,"                      Don't auto recognize, locate using four fiducual marks\n");
	fprintf(stderr," -p                   Compensate for perspective distortion\n");
	fprintf(stderr," -a                   Recognise chart in normal orientation only (-A fallback as is)\n");
	fprintf(stderr,"                      Default is to recognise all possible chart angles\n");
	fprintf(stderr," -m                   Return true mean (default is robust mean)\n");
	fprintf(stderr," -G gamma             Approximate gamma encoding of image\n");
	fprintf(stderr," -v [n]               Verbosity level 0-9\n");
	fprintf(stderr," -d [ihvglLIcrsonap]    Generate diagnostic output (try -dipn)\n");
	fprintf(stderr,"     i                diag - B&W of input image\n");
	fprintf(stderr,"     h                diag - Horizontal edge/tick detection\n");
	fprintf(stderr,"     v                diag - Vertical edge/tick detection\n");
	fprintf(stderr,"     g                diag - Groups detected\n");
	fprintf(stderr,"     l                diag - Lines detected\n");
	fprintf(stderr,"     L                diag - All lines detected\n");
	fprintf(stderr,"     I                diag - lines used to improve fit\n");
	fprintf(stderr,"     c                diag - lines perspective corrected\n");
	fprintf(stderr,"     r                diag - lines rotated\n");
	fprintf(stderr,"     s                diag - diagnostic sample boxes rotated\n");
	fprintf(stderr,"     o                diag - sample box outlines\n");
	fprintf(stderr,"     n                diag - sample box names\n");
	fprintf(stderr,"     a                diag - sample box areas\n");
	fprintf(stderr,"     p                diag - pixel areas sampled\n");
	fprintf(stderr," -O outputfile Override the default output filename & extension.\n");
	exit(1);
	}


int main(int argc, char *argv[])
{
	int fa,nfa;					/* current argument we're looking at */
	static char tiffin_name[MAXNAMEL+1] = { 0 };	/* TIFF Input file name (.tif) */
	static char datin_name[MAXNAMEL+4+1] = { 0 };	/* Data input name (.cie/.q60) */
	static char datout_name[MAXNAMEL+4+1] = { 0 };	/* Data output name (.ti3/.val) */
	static char recog_name[MAXNAMEL+1] = { 0 };		/* Reference chart name (.cht) */
	static char prof_name[MAXNAMEL+1] = { 0 };		/* scanner profile name (.cht) */
	static char diag_name[MAXNAMEL+1] = { 0 };		/* Diagnostic Output (.tif) name, if used */
	int verb = 1;
	int tmean = 0;		/* Return true mean, rather than robust mean */
	int repl = 0;		/* Replace .ti3 device values from raster file */
	int outo = 0;		/* Output the values read, rather than creating scanner .ti3 */
	int colm = 0;		/* Use inage values to measure color for print profile. > 1 == append */
	int flags = SI_GENERAL_ROT;	/* Default allow all rotations */

	TIFF *rh = NULL, *wh = NULL;
	uint16 depth, bps;			/* Useful depth, bits per sample */
	uint16 tdepth;				/* Total depth including alpha */
	uint16 pconfig, photometric;
	uint16 rextrasamples;		/* Extra "alpha" samples */
	uint16 *rextrainfo;			/* Info about extra samples */
	int gotres = 0;
	uint16 resunits;
	float resx, resy;

	icColorSpaceSignature tiffs = 0;	/* Type of tiff color space */

	int i, j;
	double gamma = 0.0;		/* default */
	double _sfid[8], *sfid = NULL;		/* Specified fiducials */
	int width, height;		/* x and y size */

	scanrd *sr;				/* Scanrd object */
	int err;	
	char *errm;
	int pnotscan = 0;		/* Number of patches that wern't scanned */

	if (argc <= 1)
		usage();

	error_program = argv[0];

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
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

			if (argv[fa][1] == '?') {
				usage();
			} else if (argv[fa][1] == 'v') {
				verb = 2;
				if (na != NULL && isdigit(na[0])) {
					verb = atoi(na);
				}
			} else if (argv[fa][1] == 'm') {
				tmean = 1;

			} else if (argv[fa][1] == 'g') {
				flags |= SI_BUILD_REF;
				repl = 0;
				outo = 0;
				colm = 0;

			} else if (argv[fa][1] == 'r') {
				repl = 1;
				outo = 0;
				colm = 0;

			} else if (argv[fa][1] == 'o') {
				repl = 0;
				outo = 1;
				colm = 0;

			} else if (argv[fa][1] == 'c') {
				repl = 0;
				outo = 0;
				colm = 1;
				if (argv[fa][2] != '\000' && argv[fa][2] == 'a')
					colm = 2;

			/* Approximate gamma encoding of image */
			} else if (argv[fa][1] == 'G') {
				fa = nfa;
				if (na == NULL) usage();
				gamma = atof(na);
				if (gamma < 0.0 || gamma > 5.0)
					usage();

			/* Use specified fiducials instead of auto recognition */
			} else if (argv[fa][1] == 'F') {
				fa = nfa;
				if (na == NULL) usage();
				if (sscanf(na, " %lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf ",
				    &_sfid[0], &_sfid[1], &_sfid[2], &_sfid[3],
				    &_sfid[4], &_sfid[5], &_sfid[6], &_sfid[7]) != 8) {
					usage();
				}

				sfid = _sfid;

			/* Compensate for perspective */
			} else if (argv[fa][1] == 'p') {
				flags |= SI_PERSPECTIVE;

			/* Don't recognise rotations */
			} else if (argv[fa][1] == 'a') {
				flags &= ~SI_GENERAL_ROT;

			/* Don't recognise rotations, and read patches */
			/* anyway "as is", if everything else failes */
			} else if (argv[fa][1] == 'A') {
				flags &= ~SI_GENERAL_ROT;
				flags |= SI_ASISIFFAIL;

			} else if (argv[fa][1] == 'd' || argv[fa][1] == 'D') {
				while (na != NULL && *na != '\000') {
					switch (*na) {
						case 'i':
							flags |= SI_SHOW_IMAGE;
							break;
						case 'h':
							flags |= SI_SHOW_DIFFSH;
							break;
						case 'v':
							flags |= SI_SHOW_DIFFSV;
							break;
						case 'g':
							flags |= SI_SHOW_GROUPS;
							break;
						case 'l':
							flags |= SI_SHOW_LINES;
							break;
						case 'L':
							flags |= SI_SHOW_ALL_LINES;
							break;
						case 'I':
							flags |= SI_SHOW_IMPL;
							break;
						case 'c':
							flags |= SI_SHOW_PERS;
							break;
						case 'r':
							flags |= SI_SHOW_ROT;
							break;
						case 's':
							flags |= SI_SHOW_SBOX;
							break;
						case 'o':
							flags |= SI_SHOW_SBOX_OUTLINES;
							break;
						case 'n':
							flags |= SI_SHOW_SBOX_NAMES;
							break;
						case 'a':
							flags |= SI_SHOW_SBOX_AREAS;
							break;
						case 'p':
							flags |= SI_SHOW_SAMPLED_AREA;
							break;
						default:
							usage();
					}
					na++;
				}

			/* Output file name */
			} else if (argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage();
				strncpy(datout_name,na,MAXNAMEL); datout_name[MAXNAMEL] = '\000';

			} else 
				usage();
		} else
			break;
	}

	/* TIFF Raster input file name */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(tiffin_name,argv[fa],MAXNAMEL); tiffin_name[MAXNAMEL] = '\000';

	/* Create a desination file path and name */
	if (datout_name[0] == '\000'		/* Not been overridden */
	 && (flags & SI_BUILD_REF) == 0
	 && repl == 0 && colm == 0) {		/* Not generate ref or replacing .ti3 dev */
										// ~~~99 Hmm. Should we honour -O ??
		char *xl;
		strncpy(datout_name,argv[fa],MAXNAMEL); datout_name[MAXNAMEL] = '\000';
		if ((xl = strrchr(datout_name, '.')) == NULL)	/* Figure where extention is */
			xl = datout_name + strlen(datout_name);
		if (outo == 0)	/* Creating scan calib data */
			strcpy(xl,".ti3");
		else			/* Just outputing values for some other purpose */
			strcpy(xl,".val");
	}

	/* .cht Reference file in or out */
	if (++fa >= argc || argv[fa][0] == '-') usage();
	strncpy(recog_name,argv[fa],MAXNAMEL); recog_name[MAXNAMEL] = '\000';

	if (colm > 0) {
		if (++fa >= argc || argv[fa][0] == '-') usage();
		strncpy(prof_name,argv[fa],MAXNAMEL); prof_name[MAXNAMEL] = '\000';
	}

	/* CGATS Data file input/output */
	if ((flags & SI_BUILD_REF) == 0 && outo == 0) {	/* Not generate ref or just outputing */
		if (++fa >= argc || argv[fa][0] == '-') usage();
		if (outo == 0) {	/* Creating scan calib data */
			/* Data file */
			strncpy(datin_name,argv[fa],MAXNAMEL); datin_name[MAXNAMEL] = '\000';
		}
		if (repl != 0 || colm > 0) {	/* Color from image or replacing .ti3 device data */
			strcpy(datin_name,argv[fa]);
			strcat(datin_name,".ti2");
			strcpy(datout_name,argv[fa]);		// ~~~99 Hmm. Should we honour -O ??
			strcat(datout_name,".ti3");
		}
	}

	/* optional diagnostic file */
	if (++fa < argc) {
		if (argv[fa][0] == '-')
			usage();
		strncpy(diag_name,argv[fa],MAXNAMEL); diag_name[MAXNAMEL] = '\000';
	} else {	/* Provide a default name */
		strcpy(diag_name,"diag.tif");
	}

	if (stricmp(diag_name, tiffin_name) == 0) {
		error("Diagnostic output '%s' might overwrite the input '%s'!",diag_name,tiffin_name);
	}

	/* ----------------------------------------- */
	/* Open up input tiff file ready for reading */
	/* Got arguments, so setup to process the file */
	if ((rh = TIFFOpen(tiffin_name, "r")) == NULL)
		error("error opening read file '%s'",tiffin_name);

	TIFFGetField(rh, TIFFTAG_IMAGEWIDTH,  &width);
	TIFFGetField(rh, TIFFTAG_IMAGELENGTH, &height);

	TIFFGetField(rh, TIFFTAG_BITSPERSAMPLE, &bps);
	if (bps != 8 && bps != 16)
		error("TIFF Input file '%s' must be 8 or 16 bits/channel",tiffin_name);

	/* See if there are alpha planes */
	TIFFGetFieldDefaulted(rh, TIFFTAG_EXTRASAMPLES, &rextrasamples, &rextrainfo);

	TIFFGetField(rh, TIFFTAG_SAMPLESPERPIXEL, &depth);

	if (rextrasamples > 0 && verb)
		printf("%d extra (alpha ?) samples will be ignored\n",rextrasamples);

	tdepth = depth;
	depth = tdepth - rextrasamples;

	if (depth != 1 && depth != 3 && depth != 4)
		error("Input '%s' must be a Grey, RGB or CMYK tiff file",tiffin_name);

	TIFFGetField(rh, TIFFTAG_PHOTOMETRIC, &photometric);
	if (depth == 1 && photometric != PHOTOMETRIC_MINISBLACK
	               && photometric != PHOTOMETRIC_MINISWHITE)
		error("1 chanel input '%s' must be a Grey tiff file",tiffin_name);
	else if (depth == 3 && photometric != PHOTOMETRIC_RGB)
		error("3 chanel input '%s' must be an RGB tiff file",tiffin_name);
	else if (depth == 4 && photometric != PHOTOMETRIC_SEPARATED)
		error("4 chanel input '%s' must be a CMYK tiff file",tiffin_name);

	if (depth == 1)
		tiffs = icSigGrayData;
	else if (depth == 3) 
		tiffs = icSigRgbData;
	else if (depth == 4)
		tiffs = icSigCmykData;

	TIFFGetField(rh, TIFFTAG_PLANARCONFIG, &pconfig);
	if (pconfig != PLANARCONFIG_CONTIG)
		error("TIFF Input file '%s' must be planar",tiffin_name);


	if (TIFFGetField(rh, TIFFTAG_RESOLUTIONUNIT, &resunits) != 0) {
		TIFFGetField(rh, TIFFTAG_XRESOLUTION, &resx);
		TIFFGetField(rh, TIFFTAG_YRESOLUTION, &resy);

		if (resunits == RESUNIT_NONE		/* If it looks valid */
		 || resunits == RESUNIT_INCH
		 || resunits == RESUNIT_CENTIMETER)
			gotres = 1;
	}

	/* -------------------------- */
	/* setup the diag output file */
	if (flags & SI_SHOW_FLAGS) {
		if ((wh = TIFFOpen(diag_name, "w")) == NULL)
			error("Can\'t create TIFF file '%s'!",diag_name);
	
		TIFFSetField(wh, TIFFTAG_IMAGEWIDTH,  width);
		TIFFSetField(wh, TIFFTAG_IMAGELENGTH, height);
		TIFFSetField(wh, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
		TIFFSetField(wh, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(wh, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(wh, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(wh, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(wh, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
		if (gotres) {
			TIFFSetField(wh, TIFFTAG_RESOLUTIONUNIT, resunits);
			TIFFSetField(wh, TIFFTAG_XRESOLUTION, resx);
			TIFFSetField(wh, TIFFTAG_YRESOLUTION, resy);
		}
		TIFFSetField(wh, TIFFTAG_IMAGEDESCRIPTION, "Scanin diagnosis output");
	}

	/* -------------------------- */
	if (verb >= 2) {
		printf("Input file '%s': w=%d, h=%d, d = %d, bpp = %d\n",
		        tiffin_name, width, height, depth, bps);
		if (flags & SI_BUILD_REF)
			printf("Build Scan Chart reference file '%s'\n",recog_name);
		else {
			printf("Data input file '%s'\n",datin_name);
			printf("Data output file '%s'\n",datout_name);
			printf("Chart reference file '%s'\n",recog_name);
		}
		if (flags & SI_SHOW_FLAGS)
			printf("Creating diagnostic tiff file '%s'\n",diag_name);
	}

	/* -------------------------- */
	/* Do the operation */

	if ((sr = do_scanrd(
		flags,			/* option flags */
		verb,			/* verbosity level */

		gamma,
		sfid,			/* Specified fiducuals, if any */
		width, height, depth, tdepth, bps,	/* Width, Height and Depth of input in pixels */
		read_line,		/* Read line function */
		(void *)rh,		/* Opaque data for read_line */

		recog_name,		/* reference file name */

		write_line,		/* Write line function */
		(void *)wh		/* Opaque data for write_line */
	)) == NULL) {
		if (flags & SI_SHOW_FLAGS)
			TIFFClose(wh);
		error("Unable to allocate scanrd object");
	}

	if ((err = sr->error(sr, &errm)) != 0) {
		if ((flags & SI_SHOW_FLAGS) && err != SI_DIAG_WRITE_ERR)
			TIFFClose(wh);		/* Close diagnostic file */
		error("Scanin failed with code 0x%x, %s",err,errm);
	}

	/* Read an output the values */
	if ((flags & SI_BUILD_REF) == 0) {	/* Not generate ref */

		/* -------------------------------------------------- */
		if (outo != 0) {		/* Just output the values */
								/* Note value range is raw 0..255, */
								/* while all others output formats are out of 100 */
			cgats *ocg;			/* output cgats structure */
			time_t clk = time(0);
			struct tm *tsp = localtime(&clk);
			char *atm = asctime(tsp); /* Ascii time */
	
			/* Setup output cgats file */
			ocg = new_cgats();	/* Create a CGATS structure */
			ocg->add_other(ocg, "VALS"); 	/* Dummy type */
			ocg->add_table(ocg, tt_other, 0);	/* Start the first table */
	
			ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration raster values",NULL);
			ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll scanin", NULL);
			atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
			ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	
			ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
			if (depth == 1) {
				ocg->add_field(ocg, 0, "GREY", r_t);
			} else if (depth == 3) {
				ocg->add_field(ocg, 0, "RGB_R", r_t);
				ocg->add_field(ocg, 0, "RGB_G", r_t);
				ocg->add_field(ocg, 0, "RGB_B", r_t);
			} else if (depth == 4) {
				ocg->add_field(ocg, 0, "CMYK_C", r_t);
				ocg->add_field(ocg, 0, "CMYK_M", r_t);
				ocg->add_field(ocg, 0, "CMYK_Y", r_t);
				ocg->add_field(ocg, 0, "CMYK_K", r_t);
			}
	
			/* Initialise, ready to read out all the values */
			for (j = 0; ; j++) {
				char id[100];		/* Input patch id */
				double P[4];		/* Robust/true mean values */
				int pixcnt;			/* PIxel count */

				if (tmean) {
					if (sr->read(sr, id, NULL, P, NULL, &pixcnt) != 0)
						break;
				} else {
					if (sr->read(sr, id, P, NULL, NULL, &pixcnt) != 0)
						break;
				}
		
				if (pixcnt == 0)
					pnotscan++;

				if (depth == 1) {
					ocg->add_set( ocg, 0, id, P[0]);
				} else if (depth == 3) {
					ocg->add_set( ocg, 0, id, P[0], P[1], P[2]);
				} else if (depth == 4) {
					ocg->add_set( ocg, 0, id, P[0], P[1], P[2], P[3]);
				}
			}
	
			if (verb)
				printf("Writing output values to file '%s'\n",datout_name);
				
			if (ocg->write_name(ocg, datout_name))
				error("Write error to '%s' : %s",datout_name,ocg->err);
	
			ocg->del(ocg);		/* Clean up */

		/* -------------------------------------------------- */
		} else if (repl != 0) {	/* Replace .ti3 device values */
			cgats *icg;			/* input .ti2 cgats structure */
			cgats *ocg;			/* input/output .ti3 cgats structure */
			int npat;			/* Number of test patches */
			int dim = 0;		/* Dimenstionality of device space */
			int fi;				/* Field index */
			int isi, ili;		/* Input file sample and location indexes */
			char *dfnames[5][4] = {	/* Device colorspace names */
				{ "" },
				{ "GRAY_W" },
				{ "" },
				{ "RGB_R", "RGB_G", "RGB_B" },
				{ "CMYK_C", "CMYK_M", "CMYK_Y", "CMYK_K" }
			};
			int odim = 0;		/* Output file device dimensionality */
			int dfi[5][4];		/* Output file device colorspace indexes */
			int osi;			/* Output file sample id index */
	
			/* Setup input .ti2 file */
			icg = new_cgats();			/* Create a CGATS structure */
			icg->add_other(icg, "CTI2"); 	/* Calibration Target Information 2 */
			if (icg->read_name(icg, datin_name))
				error("CGATS file '%s' read error : %s",datin_name,icg->err);
	
			if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
				error("Input file '%s' isn't a CTI2 format file",datin_name);

			if (icg->ntables < 1)
				error("Input file '%s' doesn't contain at least one table",datin_name);
	
			if ((npat = icg->t[0].nsets) <= 0)
				error("Input file '%s' doesn't contain any data sets",datin_name);
	
			/* Figure out the color space */
			if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
				error("Input file '%s' doesn't contain keyword COLOR_REP",datin_name);

			if (strcmp(icg->t[0].kdata[fi],"CMYK") == 0) {
				dim = 4;
			} else if (strcmp(icg->t[0].kdata[fi],"RGB") == 0) {
				dim = 3;
			} else if (strcmp(icg->t[0].kdata[fi],"W") == 0) {
				dim = 1;
			} else
				error("Input file '%s' keyword COLOR_REP has unknown value",datin_name);

			/* Find fields we want in the input file */
			if ((isi = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
				error("Input file '%s' doesn't contain field SAMPLE_ID",datin_name);
			if (icg->t[0].ftype[isi] != nqcs_t)
				error("Input file '%s' Field SAMPLE_ID is wrong type",datin_name);
		
			if ((ili = icg->find_field(icg, 0, "SAMPLE_LOC")) < 0)
				error("Input file '%s' doesn't contain field SAMPLE_LOC",datin_name);
			if (icg->t[0].ftype[ili] != cs_t
			 && icg->t[0].ftype[ili] != nqcs_t)
				error("Input file '%s' Field SAMPLE_LOC is wrong type",datin_name);

			/* Setup input/output .ti3 file */
			ocg = new_cgats();			/* Create a CGATS structure */
			ocg->add_other(ocg, "CTI3"); 	/* Calibration Target Information 3 */
			if (ocg->read_name(ocg, datout_name))
				error("CGATS file '%s' read error : %s",datout_name,ocg->err);
	
			if (ocg->t[0].tt != tt_other || ocg->t[0].oi != 0)
				error("Input file '%s' isn't a CTI3 format file",datout_name);

			if (ocg->ntables < 1)
				error("Input file '%s' doesn't contain at least one table",datout_name);
	
			if (npat != ocg->t[0].nsets)
				error("Input file '%s' doesn't contain same number of data sets",datout_name);
	

			/* Find the fields we want in the output file */

			/* Figure out the color space */
			if ((fi = ocg->find_kword(ocg, 0, "COLOR_REP")) < 0)
				error("Input file '%s' doesn't contain keyword COLOR_REP",datout_name);

			if (strncmp(ocg->t[0].kdata[fi],"CMYK",4) == 0) {
				odim = 4;
			} else if (strncmp(ocg->t[0].kdata[fi],"RGB",3) == 0) {
				odim = 3;
			} else if (strncmp(ocg->t[0].kdata[fi],"W",1) == 0) {
				odim = 1;
			} else
				error("Input file '%s' keyword COLOR_REP has unknown value",datout_name);

			if (odim != dim)
				error("File '%s' has different device space to '%s'",datin_name, datout_name);

			if ((osi = ocg->find_field(ocg, 0, "SAMPLE_ID")) < 0)
				error("Input file '%s' doesn't contain field SAMPLE_ID",datout_name);
			if (ocg->t[0].ftype[osi] != nqcs_t)
				error("Input file '%s' Field SAMPLE_ID is wrong type",datout_name);
	
			for (i = 0; i < dim; i++) {
				if ((dfi[dim][i] = ocg->find_field(ocg, 0, dfnames[dim][i])) < 0)
					error("Input '%s' file doesn't contain field %s", datout_name, dfnames[dim][i]);
				if (ocg->t[0].ftype[dfi[dim][i]] != r_t)
					error("Input '%s' Field %s is wrong type",datout_name,dfnames[dim][i]);
			}

			/* Initialise, ready to read out all the values */
			for (i = sr->reset(sr); i > 0; i--) {	/* For all samples in .tiff file */
				char loc[100];		/* Target patch location */
				double P[4];		/* Robust/raw mean values */
				int pixcnt;			/* Pixel count */
				int k, e;

				if (tmean)
					sr->read(sr, loc, NULL, P, NULL, &pixcnt);
				else
					sr->read(sr, loc, P, NULL, NULL, &pixcnt);
		
				if (pixcnt == 0)
					pnotscan++;

				/* Search for this location in the .ti2 file */
				for (j = 0; j < npat; j++) {
					if (strcmp(loc, (char *)icg->t[0].fdata[j][ili]) == 0) {
						char *sidp = (char *)icg->t[0].fdata[j][isi];

						/* Search for this sample id in .ti3 file */
						for (k = 0; k < npat; k++) {
							if (strcmp(sidp, (char *)ocg->t[0].fdata[k][osi]) == 0) {
								/* Update the device values */
								for (e = 0; e < dim; e++) {
									double vv = 100.0 * P[e]/255.0;
									*((double *)ocg->t[0].fdata[k][dfi[dim][e]]) = vv;
								}
								break;
							}	
						}
						if (k >= npat && verb >= 1)
							printf("Warning: Couldn't find sample '%s' in '%s'\n",sidp,datout_name);
						break;
					}
				}
				if (j >= npat && verb >= 1)
					printf("Warning: Couldn't find location '%s' in '%s'\n",loc,datin_name);
			}
	
			/* Flush our changes */
			if (verb)
				printf("Writing output values to file '%s'\n",datout_name);
				
			if (ocg->write_name(ocg, datout_name))
				error("Write error to file '%s' : %s",datout_name,ocg->err);
	
			ocg->del(ocg);		/* Clean up */
			icg->del(icg);		/* Clean up */

		/* ---------------------------------------------------------- */
		} else if (colm > 0) {	/* Using the image to measure color */
			/* All this needs to track the code in spectro/printread.c */
			cgats *icg;			/* input .ti2 cgats structure */
			cgats *ocg;			/* input/output .ti3 cgats structure */
			icmFile *rd_fp = NULL;	/* Image to CIE lookup */
			icc *rd_icco = NULL;
			icmLuBase *luo;
			mpp *mlu = NULL;
			time_t clk = time(0);
			struct tm *tsp = localtime(&clk);
			char *atm = asctime(tsp); /* Ascii time */
			int nmask = 0;		/* Device colorant mask */
			int nchan = 0;		/* Number of device chanels */
			int npat;			/* Number of input patches (inc. padding) */
			int nopat = 0;		/* Number of output patches */
			int si;				/* Sample id index */
			int li;				/* Location id index */
			int ti;				/* Temp index */
			int fi;				/* Colorspace index */
	
			icg = new_cgats();			/* Create a CGATS structure */
			icg->add_other(icg, "CTI2"); 	/* special type Calibration Target Information 2 */
			icg->add_other(icg, "CAL"); 	/* There may be a calibration too */
		
			if (icg->read_name(icg, datin_name))
				error("CGATS file '%s' read error : %s",datin_name,icg->err);
		
			if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
				error("Input file '%s' isn't a CTI2 format file",datin_name);
			if (icg->ntables < 1)
				error("Input file '%s' doesn't contain at least one table",datin_name);
		
			if ((npat = icg->t[0].nsets) <= 0)
				error("Input file '%s' has no sets of data",datin_name);
		
			/* Setup output cgats file */
			ocg = new_cgats();			/* Create a CGATS structure */
			ocg->add_other(ocg, "CTI3");	/* special type Calibration Target Information 3 */

			if (colm > 1) {		/* Appending information to .ti3 */

				if (ocg->read_name(ocg, datout_name))
					error("CGATS file read error on '%s': %s",datout_name, ocg->err);
		
				if (ocg->t[0].tt != tt_other || ocg->t[0].oi != 0)
					error("Input file '%s' isn't a CTI3 format file",datout_name);
				if (ocg->ntables < 1)
					error("Input file '%s' doesn't at least exactly one table",datout_name);
				if ((nopat = ocg->t[0].nsets) <= 0)
					error("Input file '%s' has no existing sets of data",datout_name);

			} else {			/* Creating .ti3 */

				ocg->add_table(ocg, tt_other, 0);	/* Start the first table */
		
				ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
				ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll printread", NULL);
				atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
				ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
				ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */
				if ((ti = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
					ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[ti], NULL);
			
				if ((ti = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
					ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[ti], NULL);
		
				if ((ti = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
					ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[ti], NULL);
		
				if ((ti = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
					ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[ti], NULL);

				if ((ti = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0)
					ocg->add_kword(ocg, 0, "TOTAL_INK_LIMIT",icg->t[0].kdata[ti], NULL);

				/* See if there is a calibration in the .ti2, and copy it if there is */
				{
					int oi, tab;
			
					oi = icg->get_oi(icg, "CAL");
			
					for (tab = 0; tab < icg->ntables; tab++) {
						if (icg->t[tab].tt == tt_other && icg->t[tab].oi == oi) {
							break;
						}
					}
					if (tab < icg->ntables) {
						xcal *cal = NULL;
			
						if (verb)
							printf("Copying .cal from '%s' to '%s'\n",datin_name,datout_name);
						
						if ((cal = new_xcal()) == NULL) {
							error("new_xcal failed");
						}
						if (cal->read_cgats(cal, icg, tab, datin_name) != 0)  {
							error("%s",cal->err);
						}
			
						if (cal->write_cgats(cal, ocg)) {
							error("%s",cal->err);
						}
			
						cal->del(cal);
					}
				}
			}

			if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
				error("Input file '%s' doesn't contain field SAMPLE_ID",datin_name);
			if (icg->t[0].ftype[si] != nqcs_t)
				error("Input file '%s' Field SAMPLE_ID is wrong type",datin_name);
		
			/* Fields we want */
			if (colm > 1) {		/* Appending information to .ti3 */
				if ((ti = ocg->find_field(ocg, 0, "SAMPLE_ID")) != 0)
					error("Input file '%s' field SAMPLE_ID (%d) not in expected location (%d)",
					       datout_name, ti, 0);
			} else {
				ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
			}
		
			if ((li = icg->find_field(icg, 0, "SAMPLE_LOC")) < 0)
				error("Input file '%s' doesn't contain field SAMPLE_LOC",datin_name);
			if (icg->t[0].ftype[li] != cs_t)
				error("Input file '%s' field SAMPLE_LOC is wrong type",datin_name);
		
			/* Figure out the color space */
			if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
				error("Input file '%s' doesn't contain keyword COLOR_REPS",datin_name);
		
			if ((nmask = icx_char2inkmask(icg->t[0].kdata[fi])) != 0) {
				int i, j, ii;
				int chix[ICX_MXINKS];	/* Device chanel indexes */
				int xyzix[3];			/* XYZ chanel indexes */
				char *ident, *bident;
				char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };
		
				nchan = icx_noofinks(nmask);
				ident = icx_inkmask2char(nmask, 1); 
				bident = icx_inkmask2char(nmask, 0); 
		
				/* Device channels */
				for (j = 0; j < nchan; j++) {
					int imask;
					char fname[100];
		
					imask = icx_index2ink(nmask, j);
					sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : bident,
					                      icx_ink2char(imask));
		
					if ((ii = icg->find_field(icg, 0, fname)) < 0)
						error("Input file '%s' doesn't contain field %s",datin_name,fname);
					if (icg->t[0].ftype[ii] != r_t)
						error("Field %s is wrong type",fname);
			
					if (colm > 1) {		/* Appending information to .ti3 */
						if (ocg->find_field(ocg, 0, fname) != 1 + j)
							error("Input file '%s' field %s not in expected location",datout_name,fname);
					} else {
						ocg->add_field(ocg, 0, fname, r_t);
					}
					chix[j] = ii;
				}
		
				/* Approximate XYZ and real XYZ */
				for (j = 0; j < 3; j++) {
					if ((ii = icg->find_field(icg, 0, xyzfname[j])) >= 0) {

						if (icg->t[0].ftype[ii] != r_t)
							error("Input file '%s' field %s is wrong type",datin_name,xyzfname[j]);
					}
			
					if (colm > 1) {		/* Appending information to .ti3 */
						if (ocg->find_field(ocg, 0, xyzfname[j]) != 1 + nchan + j)
							error("Input file '%s' field %s not in expected location",
							                            datout_name,xyzfname[j]);
					} else {
						ocg->add_field(ocg, 0, xyzfname[j], r_t);
					}
					xyzix[j] = ii;
				}
		
				if (colm <= 1) {		/* Creating .ti3 */
					char fname[100];
					sprintf(fname, "%s_XYZ", ident);
					ocg->add_kword(ocg, 0, "COLOR_REP", fname, NULL);
				}
		
				if (colm > 1) {		/* Appending .ti3 data */

					/* Check that all the patches match */
					for (ii = i = 0; i < npat; i++) {

						if (strcmp(((char *)icg->t[0].fdata[i][si]), "0") == 0)
							continue;			/* Padding, so skip it */
		
						/* Id's */
						if (strcmp (((char *)icg->t[0].fdata[i][si]),
								    ((char *)ocg->t[0].fdata[ii][si])) != 0)
							error("'%s' and '%s' field id's don't match at patch %d\n",datin_name,datout_name,i+1);

						/* device values */
						for (j = 0; j < nchan; j++) {
							double ival, oval;
							ival = *((double *)icg->t[0].fdata[i][chix[j]]);
							oval = *((double *)ocg->t[0].fdata[ii][1 + j]);
							if (fabs(ival - oval) > 0.001)
								error("'%s' and '%s' device values (%f %f) don't match at patch %d %d\n",datin_name,datout_name,ival, oval, i+1, ii+1);
						}
						ii++;
					}
					if (ii != nopat)
						error("Different number of patches in '%s' (%d) to expected(%d)",datout_name,nopat,ii);

				} else { /* Read all the test patches in, and create output slots */
					cgats_set_elem *setel;	/* Array of set value elements */

					if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * (1 + nchan + 3))) == NULL)
						error("Malloc failed!");


					for (ii = i = 0; i < npat; i++) {
						int k = 0;

						if (strcmp(((char *)icg->t[0].fdata[i][si]), "0") == 0)
							continue;			/* Padding, so skip it */
		
						/* Id */
						setel[k++].c = ((char *)icg->t[0].fdata[i][si]);
					
						/* device values */
						for (j = 0; j < nchan; j++) {
							setel[k++].d = *((double *)icg->t[0].fdata[i][chix[j]]);
						}

						/* Unset XYZ values */
						setel[k++].d = -1.0;
						setel[k++].d = -1.0;
						setel[k++].d = -1.0;

						ocg->add_setarr(ocg, 0, setel);

						ii++;
					}
					nopat = ii;
					free(setel);
				}
				free(ident);
				free(bident);

			} else
				error("Input file '%s' keyword COLOR_REPS has unknown value",datin_name);
		
			/* Setup RGB to XYZ conversion */
			{
				int inn, outn;			/* Chanels for input and output spaces */
				icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
				int rv;

				/* Open up the file for reading */
				if ((rd_fp = new_icmFileStd_name(prof_name,"r")) == NULL)
					error("Write: Can't open file '%s'",prof_name);
		
				if ((rd_icco = new_icc()) == NULL)
					error("Read: Creation of ICC object failed");
		
				/* Read the header and tag list */
				if ((rv = rd_icco->read(rd_icco,rd_fp,0)) == 0) {
		
					/* Get the Fwd table, absolute with XYZ override */
					if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, icAbsoluteColorimetric,
					                              icSigXYZData, icmLuOrdNorm)) == NULL) {
						error("%d, %s",rd_icco->errc, rd_icco->err);
					}

					/* Get details of conversion */
					luo->spaces(luo, &ins, &inn, &outs, &outn, NULL, NULL, NULL, NULL, NULL);

					/* Check that it matches what we expect */

				} else {	/* Not a valid ICC */
					inkmask cnv_nmask = 0;  /* Conversion input nmask */

					/* Close out the ICC profile */
					rd_icco->del(rd_icco);
					rd_icco = NULL;
					rd_fp->del(rd_fp);
					rd_fp = NULL;

					/* If we don't have an ICC lookup object, look for an MPP */

					if ((mlu = new_mpp()) == NULL)
						error ("Creation of MPP object failed");

					if ((rv = mlu->read_mpp(mlu, prof_name)) == 0) {

						/* mlu defaults to absolute XYZ lookup */
						mlu->get_info(mlu, &cnv_nmask, &inn, NULL, NULL, NULL, NULL, NULL, NULL);
				
						outn = 3;
						outs = icSigXYZData;
				
						if ((ins = icx_colorant_comb_to_icc(cnv_nmask)) == 0)
							error ("Couldn't match MPP mask to valid ICC colorspace");

					} else {
						mlu->del(mlu);
						mlu = NULL;
						error("File '%s' failed to read as ICC or MPP profile",prof_name);
					}
				}
				if (inn != depth || tiffs != ins)
					error("%s profile '%s' doesn't match TIFF file type",luo != NULL ? "ICC" : "MPP", prof_name);
			}

			/* Initialise, ready to read out all the values */
			for (i = sr->reset(sr); i > 0; i--) {	/* For all samples in .tiff file */
				char loc[100];		/* Target patch location */
				double P[ICX_MXINKS];	/* Robust/true mean values */
				double xyz[3];			/* profile XYZ value */
				int pixcnt;				/* Pixel count */
				int k, e;

				if (tmean)
					sr->read(sr, loc, NULL, P, NULL, &pixcnt);
				else
					sr->read(sr, loc, P, NULL, NULL, &pixcnt);
		
				if (pixcnt == 0)
					pnotscan++;

				/* Search for this location in the .ti2 file */
				for (j = 0; j < npat; j++) {
					if (strcmp(loc, (char *)icg->t[0].fdata[j][li]) == 0) {	/* Got location */
						char *sidp = (char *)icg->t[0].fdata[j][si];	/* Get id */

						if (strcmp(sidp, "0") == 0)
							break;			/* Padding, so ignore it */

						/* Search for this sample id in .ti3 file */
						for (k = 0; k < nopat; k++) {
							if (strcmp(sidp, (char *)ocg->t[0].fdata[k][si]) == 0) {
								
//printf("Loc %s, ID %s got RGB value %f %f %f\n",sidp, loc, P[0], P[1], P[2]);

								/* Convert RGB to XYZ */
								for (e = 0; e < depth; e++)
									P[e] /= 255.0;			/* Convert to 0.0 .. 1.0 range */

								/* Convert to XYZ */
								if (luo != NULL)
									luo->lookup(luo, xyz, P);
								else
									mlu->lookup(mlu, xyz, P);

								/* Sanity check XYZ ? */
								// ~~~99

								/* Update the XYZ values */
								for (e = 0; e < 3; e++) {
									double ev = *((double *)ocg->t[0].fdata[k][1 + nchan + e]);

									if (ev != -1.0)
										error("Found an existing value in '%s' file (%f)",datout_name,ev);

									*((double *)ocg->t[0].fdata[k][1 + nchan + e]) = 100.0 * xyz[e];
								}
								break;
							}	
						}
						if (k >= nopat)
							error("Couldn't find sample '%s' in '%s'\n",sidp,datout_name);
						break;
					}
				}
				if (j >= npat && verb >= 1)
					error("Couldn't find location '%s' in '%s'\n",loc, datin_name);
			}
	
			/* Warn if not all patch values have been filled */
			if (verb) {
				int e, k;
				for (k = 0; k < nopat; k++) {
					for (e = 0; e < 3; e++) {
						double ev = *((double *)ocg->t[0].fdata[k][1 + nchan + e]);

						if (ev == -1.0)
							break;
					}
					if (e < 3)
						break;
				}
				if (k < nopat)
					printf("Not all sample values have been filled\n");
				else
					printf("All sample values have been filled\n");
			}

			if (verb)
				printf("Writing output values to file '%s'\n",datout_name);
				
			if (ocg->write_name(ocg, datout_name))
				error("File '%s' write error : %s",datout_name,ocg->err);
		
			if (luo != NULL)
				luo->del(luo);
			if (rd_icco != NULL)
				rd_icco->del(rd_icco);
			if (rd_fp != NULL)
				rd_fp->del(rd_fp);
			if (mlu != NULL)
				mlu->del(mlu);

			ocg->del(ocg);		/* Clean up */
			icg->del(icg);		/* Clean up */
		
		/* ----------------------------------- */
		} else {	/* Normal scan calibration */
			cgats *icg;			/* input cgats structure */
			cgats *ocg;			/* output cgats structure */
			time_t clk = time(0);
			struct tm *tsp = localtime(&clk);
			char *atm = asctime(tsp); /* Ascii time */
			int ti;				/* Temp index */
			int sx;				/* Sample id index */
			int isLab = 0;		/* D50 Lab reference */
			int Xx, Yx, Zx;		/* XYZ_X, XYZ_Y, XYZ_Z index */
			int spec_n = 0;		/* Number of spectral bands */
			double spec_wl_short;/* First reading wavelength in nm (shortest) */
			double spec_wl_long; /* Last reading wavelength in nm (longest) */
			int spi[XSPECT_MAX_BANDS];  /* CGATS indexes for each wavelength */
			int npat;			/* Number of test patches in it8 chart */
			int nsetel = 0;		/* Number of output set elements */
			cgats_set_elem *setel;  /* Array of set value elements */
			unsigned int *idhash; 	/* Array of reference id hashes */
	
			icg = new_cgats();			/* Create a CGATS structure */
			icg->add_other(icg, ""); 	/* Accept any type */
			if (icg->read_name(icg, datin_name))
				error("CGATS file '%s' read error : %s",datin_name,icg->err);
	
			/* ~~ should accept ti2 file and convert RGB to XYZ using    */
			/*    device cal., to make W/RGB/CMYK ->XYZ reading chart ~~ */
			if (icg->ntables < 1)
				error("Input file '%s' doesn't contain at least one table",datin_name);
	
			if ((npat = icg->t[0].nsets) <= 0)
				error("File '%s' no sets of data in first table",datin_name);
	
			/* Setup output cgats file */
			ocg = new_cgats();	/* Create a CGATS structure */
			ocg->add_other(ocg, "CTI3"); 	/* our special type is Calibration Target Information 3 */
			ocg->add_table(ocg, tt_other, 0);	/* Start the first table */
	
			ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
			ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll target", NULL);
			atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
			ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	
			ocg->add_kword(ocg, 0, "DEVICE_CLASS","INPUT", NULL);	/* What sort of device this is */
			ocg->add_kword(ocg, 0, "COLOR_REP","XYZ_RGB", NULL);
	
			/* Fields we want from input chart reference file */
			if ((sx = icg->find_field(icg, 0, "Sample_Name")) < 0) {
				if ((sx = icg->find_field(icg, 0, "SAMPLE_NAME")) < 0) {
					if ((sx = icg->find_field(icg, 0, "SAMPLE_LOC")) < 0) {
						if ((sx = icg->find_field(icg, 0, "SAMPLE_ID")) < 0) {
							error("Input file '%s' doesn't contain field SAMPLE_ID, Sample_Name or SAMPLE_NAME",datin_name);
						}
					}
				}
			}
			if (icg->t[0].ftype[sx] != nqcs_t && icg->t[0].ftype[sx] != cs_t)
				error("Input file '%s' field %s is wrong type", datin_name, icg->t[0].fsym[sx]);

			if ((Xx = icg->find_field(icg, 0, "XYZ_X")) < 0) {
				if ((Xx = icg->find_field(icg, 0, "LAB_L")) < 0)
					error("Input file '%s' doesn't contain field XYZ_X or LAB_L",datin_name);
				
				isLab = 1;
				if (icg->t[0].ftype[Xx] != r_t)
					error("Input file '%s' field LAB_L is wrong type",datin_name);
				if ((Yx = icg->find_field(icg, 0, "LAB_A")) < 0)
					error("Input file doesn't contain field LAB_A",datin_name);
				if (icg->t[0].ftype[Yx] != r_t)
					error("Input file '%s' field LAB_A is wrong type",datin_name);
				if ((Zx = icg->find_field(icg, 0, "LAB_B")) < 0)
					error("Input file '%s' doesn't contain field LAB_B",datin_name);
				if (icg->t[0].ftype[Zx] != r_t)
					error("Input file '%s' field LAB_B is wrong type",datin_name);
			} else {
				if (icg->t[0].ftype[Xx] != r_t)
					error("Input file '%s' field XYZ_X is wrong type",datin_name);
				if ((Yx = icg->find_field(icg, 0, "XYZ_Y")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Y",datin_name);
				if (icg->t[0].ftype[Yx] != r_t)
					error("Input file '%s' field XYZ_Y is wrong type",datin_name);
				if ((Zx = icg->find_field(icg, 0, "XYZ_Z")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Z",datin_name);
				if (icg->t[0].ftype[Zx] != r_t)
					error("Input file '%s' field XYZ_Z is wrong type",datin_name);
			}

			/* Find possible spectral fields in reference */
			if ((ti = icg->find_kword(icg, 0, "SPECTRAL_BANDS")) >= 0) {
				spec_n = atoi(icg->t[0].kdata[ti]);
				if ((ti = icg->find_kword(icg, 0, "SPECTRAL_START_NM")) < 0)
					error ("Input file '%s' doesn't contain keyword SPECTRAL_START_NM",datin_name);
				spec_wl_short = atof(icg->t[0].kdata[ti]);
				if ((ti = icg->find_kword(icg, 0, "SPECTRAL_END_NM")) < 0)
					error ("Input file '%s' doesn't contain keyword SPECTRAL_END_NM",datin_name);
				spec_wl_long = atof(icg->t[0].kdata[ti]);
		
				/* Find the fields for spectral values */
				for (i = 0; i < spec_n; i++) {
					char buf[100];
					int nm;
			
					/* Compute nearest integer wavelength */
					nm = (int)(spec_wl_short + ((double)i/(spec_n-1.0))
					            * (spec_wl_long - spec_wl_short) + 0.5);
					
					sprintf(buf,"SPEC_%03d",nm);
		
					if ((spi[i] = icg->find_field(icg, 0, buf)) < 0)
						error("Input file doesn't contain field %s",datin_name);
				}
			}

			ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
			nsetel += 1;
			ocg->add_field(ocg, 0, "XYZ_X", r_t);
			ocg->add_field(ocg, 0, "XYZ_Y", r_t);
			ocg->add_field(ocg, 0, "XYZ_Z", r_t);
			nsetel += 3;

			/* If we have spectral information, output it too */
			if (spec_n > 0) {
				char buf[100];
	
				nsetel += spec_n;       /* Spectral values */
				sprintf(buf,"%d", spec_n);
				ocg->add_kword(ocg, 0, "SPECTRAL_BANDS",buf, NULL);
				sprintf(buf,"%f", spec_wl_short);
				ocg->add_kword(ocg, 0, "SPECTRAL_START_NM",buf, NULL);
				sprintf(buf,"%f", spec_wl_long);
				ocg->add_kword(ocg, 0, "SPECTRAL_END_NM",buf, NULL);
	
				/* Generate fields for spectral values */
				for (i = 0; i < spec_n; i++) {
					int nm;
			
					/* Compute nearest integer wavelength */
					nm = (int)(spec_wl_short + ((double)i/(spec_n-1.0))
					            * (spec_wl_long - spec_wl_short) + 0.5);
					
					sprintf(buf,"SPEC_%03d",nm);
					ocg->add_field(ocg, 0, buf, r_t);
				}
			}

			if (depth == 1) {
				ocg->add_field(ocg, 0, "GREY", r_t);
				ocg->add_field(ocg, 0, "STDEV_GREY", r_t);
			} else if (depth == 3) {
				ocg->add_field(ocg, 0, "RGB_R", r_t);
				ocg->add_field(ocg, 0, "RGB_G", r_t);
				ocg->add_field(ocg, 0, "RGB_B", r_t);
				ocg->add_field(ocg, 0, "STDEV_R", r_t);
				ocg->add_field(ocg, 0, "STDEV_G", r_t);
				ocg->add_field(ocg, 0, "STDEV_B", r_t);
			} else if (depth == 4) {
				ocg->add_field(ocg, 0, "CMYK_C", r_t);
				ocg->add_field(ocg, 0, "CMYK_M", r_t);
				ocg->add_field(ocg, 0, "CMYK_Y", r_t);
				ocg->add_field(ocg, 0, "CMYK_K", r_t);
				ocg->add_field(ocg, 0, "STDEV_C", r_t);
				ocg->add_field(ocg, 0, "STDEV_M", r_t);
				ocg->add_field(ocg, 0, "STDEV_Y", r_t);
				ocg->add_field(ocg, 0, "STDEV_K", r_t);
			}
			nsetel += 2 * depth;
	
			if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * nsetel)) == NULL)
				error("Malloc failed!");

			if ((idhash = (unsigned int *)malloc(sizeof(unsigned int) * npat)) == NULL)
				error("Malloc failed!");

			/* Setup hash list of reference labels to speed comparisons */
			for (j = 0; j < npat; j++) {
				char id[100];		/* Reference patch id */

				/* Normalise reference labels */
				fix_it8(id, ((char *)icg->t[0].fdata[j][sx]));	/* Copy and fix */

				idhash[j] = shash(id);
			}

			/* Initialise, ready to read out all the values */
			for (i = sr->reset(sr); i > 0; i--) {
				char tod[100];			/* Temp output patch id */
				char od[100];			/* Output patch id */
				unsigned int odhash;	/* Chart id hashes */
				double P[4];			/* Robust/true mean values */
				double sdP[4];			/* Standard deviation */
				int pixcnt;				/* Pixel count */

				if (tmean)
					sr->read(sr, tod, NULL, P, sdP, &pixcnt);
				else
					sr->read(sr, tod, P, NULL, sdP, &pixcnt);

				fix_it8(od, tod);

				odhash = shash(od);
		
				if (pixcnt == 0)
					pnotscan++;

				/* Search for matching id in reference */
				for (j = 0; j < npat; j++) {
					char id[100];		/* Reference patch id */

					if (odhash != idhash[j]) {	/* Fast reject */
						continue;
					}

					/* Normalise reference labels */
					fix_it8(id, ((char *)icg->t[0].fdata[j][sx]));	/* Copy and fix */
	
					if (strcmp(id, od) == 0) {
						int k = 0, m;
						double XYZ[3];

						setel[k++].c = id;

				        XYZ[0] = *((double *)icg->t[0].fdata[j][Xx]);
				        XYZ[1] = *((double *)icg->t[0].fdata[j][Yx]);
				        XYZ[2] = *((double *)icg->t[0].fdata[j][Zx]);
						if (isLab) {
							icmLab2XYZ(&icmD50, XYZ, XYZ);
							XYZ[0] *= 100.0;
							XYZ[1] *= 100.0;
							XYZ[2] *= 100.0;
						}

						setel[k++].d = XYZ[0];
						setel[k++].d = XYZ[1];
						setel[k++].d = XYZ[2];

						if (spec_n > 0) {
							for (m = 0; m < spec_n; m++) {
								setel[k++].d = *((double *)icg->t[0].fdata[j][spi[m]]);
							}
						}

						for (m = 0; m < depth; m++) 
							setel[k++].d = P[m] * 100.0/255.0;
						for (m = 0; m < depth; m++) 
							setel[k++].d = sdP[m] * 100.0/255.0;

						ocg->add_setarr(ocg, 0, setel);

						break;
					}

					if (j >= npat && verb >= 1)
						printf("Warning: Couldn't match field '%s'\n",od);
				}
			}
	
			if (verb)
				printf("Writing output values to file '%s'\n",datout_name);
				
			if (ocg->write_name(ocg, datout_name))
				error("Output file '%s' write error : %s",datout_name, ocg->err);
	
			free(idhash);
			free(setel);

			ocg->del(ocg);		/* Clean up */
			icg->del(icg);		/* Clean up */

		}
	}

	if (pnotscan > 0)
		warning("A total of %d patches had no value set!",pnotscan);

	/* Clean up */
	sr->free(sr);

	TIFFClose(rh);

	if (flags & SI_SHOW_FLAGS)
		TIFFClose(wh);

	return 0;
}

/* Fix IT8 chart labels */
void fix_it8(char *o, char *i) {
	/* Some hacks for it8 */
	if (strcmp(i,"Dmin")==0) {
		strcpy(o,"GS00");
		return;
	}
	if (strcmp(i,"Dmax")==0) {
		strcpy(o,"GS23");
		return;
	}

	while (!isdigit(*i) && *i != '\000') 	/* Skip non-numbers */
		*o++ = *i++; 
	if (i[0] != '\000' && i[1] == '\000')	/* Single last digit */
		*o++ = '0';				/* Add leading zero */
	strcpy(o, i);				/* Copy remainder */
}

/********************************************************************************/
