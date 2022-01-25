
/* 
 * Argyll Color Management System
 *
 * Read in the RGB/CMYK CGATS device data from Gretag/Logo/X-Rite etc.
 * and convert it into a .ti3 CGATs format suitable for the Argyll CMS.
 *
 * Derived from  kodak2cgats.c 
 * Author: Graeme W. Gill
 * Date:   16/11/00
 *
 * Copyright 2000 - 2010, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD
	
	Need to add support for SPECTRAL_NM SPECTRAL_PCT type spectral values.
	i.e. format with wavlength then corrsponding value.
	See /src/argyll/test/JosvanRiswick/R080505.cgt

	Do we need to worry about normalising display values to Y = 100, or marking
	them not normalised ?

	Should have an option to output .ti1 files instead of .ti3 ?

	Should parse and support XRGA & Polarized tags.
	We could then add an XRGA selection for the output.

	~~

	MEASUREMENT_SOURCE
	iSis_INFO
		Illumination=D50

		ObserverAngle=2

		Filter=Unknown
		Filter=None
		Filter=No
		Filter=UVcut
		Filter=D50

		(???
			Filter=Polarizer
			Filter=D65
		???)

		MeasurementCondition=M1

	ILLUMINATION_NAME   "D50"
	OBSERVER_ANGLE  "2"
	FILTER "UVcut"

	ICCOLOR_CSV
		MeasurementOptics=Remission

	DEVCALSTD UNKNOWN
	DEVCALSTD XRDI
	DEVCALSTD GMDI
	DEVCALSTD XRGA

	Cxf:DevicePolarization
	Cxf:DeviceFilter
	SpectralPolFilter

 */

#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include "copyright.h"
#include "aconfig.h"
#include "cgats.h"
#include "xspect.h"
#include "insttypes.h"
#include "numlib.h"
#include "ui.h"

static int trans(char *dst, char *src);

void
usage(char *mes) {
	fprintf(stderr,"Convert Gretag/Logo or X-Rite ColorPport raw RGB or CMYK device profile data to Argyll CGATS data, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (mes != NULL)
		fprintf(stderr,"error: %s\n",mes);
	fprintf(stderr,"usage: txt2ti3 [-v] [-l limit] [devfile] infile [specfile] outbase\n");
/*	fprintf(stderr," -v            Verbose mode\n"); */
	fprintf(stderr," -2            Create dummy .ti2 file as well\n");
	fprintf(stderr," -l limit      set ink limit, 0 - 400%% (default max in file)\n");
	fprintf(stderr," -d            Set type of device as Display, rather than Output\n");
	fprintf(stderr," -D            Set type of device as Display and not normalized\n");
	fprintf(stderr," -i            Set type of device as Input, not Output\n");
	fprintf(stderr," -T            Transpose sample name Letters and Numbers\n");
	fprintf(stderr," [devfile]     Input Device CMYK target file (typically file.txt)\n");
	fprintf(stderr," infile        Input CIE, Spectral or Device & Spectral file (typically file.txt)\n");
	fprintf(stderr," [specfile]    Input Spectral file (typically file.txt)\n");
	fprintf(stderr," outbasename   Output file basename for .ti3 and .ti2\n");
	exit(1);
	}

int main(int argc, char *argv[])
{
	int i, j;
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int out2 = 0;			/* Create dumy .ti2 file output */
	int disp = 0;			/* 1 if this is a display device, 2 if nor normalized */
	int inp = 0;			/* nz if this is an input device */
	int transpose = 0;		/* nz to transpose letters and numbers */
	static char devname[MAXNAMEL+1] = { 0 };		/* Input CMYK/Device .txt file (may be null) */
	static char ciename[MAXNAMEL+1] = { 0 };		/* Input CIE .txt file (may be null) */
	static char specname[MAXNAMEL+1] = { 0 };		/* Input Device / Spectral .txt file */
	static char outname[MAXNAMEL+9] = { 0 };		/* Output cgats .ti3 file base name */
	static char outname2[MAXNAMEL+9] = { 0 };		/* Output cgats .ti2 file base name */
	int ti;					/* Temporary field index */
	cgats *cmy = NULL;		/* Input RGB/CMYK reference file */
	int f_id1 = -1, f_c, f_m, f_y, f_k = 0;	/* Field indexes */
	double dev_scale = 1.0;	/* Device value scaling */
	cgats *ncie = NULL;		/* Input CIE readings file (may be Dev & spectral too) */
	int f_id2, f_cie[3];	/* Field indexes */
	cgats *spec = NULL;		/* Input spectral readings (NULL if none) */
	double spec_scale = 1.0;	/* Spectral value scaling */
	int f_id3 = 0;			/* Field indexes */
	int spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
	cgats *ocg;				/* output cgats structure for .ti3 */
	cgats *ocg2;			/* output cgats structure for .ti2 */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	xcalstd calstd = xcalstd_none;	/* X-Rite Calibration standard */
	int isPolarize = 0;		/* nz if polarized measurement */
	int islab = 0;			/* CIE is Lab rather than XYZ */
	int specmin = 0, specmax = 0, specnum = 0;	/* Min and max spectral in nm, inclusive */
	int npat = 0;			/* Number of patches */
	int ndchan = 0;			/* Number of device channels, 0 = no device, RGB = 3, CMYK = 4 */
	int tlimit = -1;		/* Not set */
	double mxsum = -1.0;	/* Maximim sum of inks found in file */
	int mxsumix = 0;

	error_program = "txt2ti3";

	if (argc <= 1)
		usage("Too few arguments");

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage(NULL);

			else if (argv[fa][1] == '2')
				out2 = 1;

			else if (argv[fa][1] == 'l') {
				if (na == NULL) usage("No ink limit parameter");
				fa = nfa;
				tlimit = atoi(na);
				if (tlimit < 1)
					tlimit = -1;
			}

			else if (argv[fa][1] == 'd') {
				disp = 1;
				inp = 0;

			} else if (argv[fa][1] == 'D') {
				disp = 2;
				inp = 0;

			} else if (argv[fa][1] == 'i') {
				disp = 0;
				inp = 1;

			} else if (argv[fa][1] == 'T') {
				transpose = 1;

			} else if (argv[fa][1] == 'v')
				verb = 1;
			else 
				usage("Unknown flag");
		} else
			break;
	}

	/* See how many arguments remain */
	switch (argc - fa) {
		case 2:		/* Must be a new combined device + cie/spectral */
			if (fa >= argc || argv[fa][0] == '-') usage("Bad dev filename");
			strcpy(devname,argv[fa]);
			strcpy(ciename,argv[fa++]);
			if (fa >= argc || argv[fa][0] == '-') usage("Bad output filename");
			strcpy(outname, argv[fa++]);
			if (verb) printf("Single source file, assumed dev/cie/spectral\n");
			break;
		case 3:		/* Device + cie or spectral */
					/* or combined cie + spectral */
			if (fa >= argc || argv[fa][0] == '-') usage("Bad dev filename");
			strcpy(devname,argv[fa++]);
			if (fa >= argc || argv[fa][0] == '-') usage("Bad cie/spec filename");
			strcpy(ciename,argv[fa++]);	
			if (fa >= argc || argv[fa][0] == '-') usage("Bad output filename");
			strcpy(outname, argv[fa++]);
			if (verb) printf("Two source files, assumed dev + cie/spectral or dev/cie + spectral\n");
			break;
		case 4:		/* Device, ci and spectral */
			if (fa >= argc || argv[fa][0] == '-') usage("Bad dev filename");
			strcpy(devname,argv[fa++]);
			if (fa >= argc || argv[fa][0] == '-') usage("Bad cie filename");
			strcpy(ciename,argv[fa++]);	
			if (fa >= argc || argv[fa][0] == '-') usage("Bad spec filename");
			strcpy(specname,argv[fa++]);	
			if (fa >= argc || argv[fa][0] == '-') usage("Bad output filename");
			strcpy(outname, argv[fa++]);
			if (verb) printf("Three source files, assumed dev + cie + spectral\n");
			break;
		default:
			usage("Wrong number of filenames");
	}

	if (out2) {
		strcpy (outname2, outname);
		strcat(outname2,".ti2");
	}

	/* Convert basename into .ti3 */
	strcat(outname,".ti3");

	/* Open up the Input CMYK/RGB reference file (might be same as ncie/spec) */
	cmy = new_cgats();	/* Create a CGATS structure */
	cmy->add_other(cmy, "LGOROWLENGTH"); 	/* Gretag/Logo Target file */
	cmy->add_other(cmy, "Date:");			/* Gretag/Logo Target file */
	cmy->add_other(cmy, "ECI2002");			/* Gretag/Logo Target file */
	cmy->add_other(cmy, ""); 				/* Wildcard */
	if (cmy->read_name(cmy, devname))
		error ("Read: Can't read dev file '%s'. Unknown format, missing or corrupted file ? (%s)",devname,cmy->err);
	if (cmy->ntables != 1)
		warning("Input file '%s' doesn't contain exactly one table",devname);

	if ((npat = cmy->t[0].nsets) <= 0)
		error("No patches");

	/* See if we can figure out polarization and xrga */
	if ((ti = cmy->find_field(cmy, 0, "MEASUREMENT_SOURCE")) >= 0) {
		if (strstr(cmy->t[0].kdata[ti], "Filter=Polarizer") != NULL)
			isPolarize = 1;
	}
	if ((ti = cmy->find_field(cmy, 0, "FILTER")) >= 0) {
		if (strstr(cmy->t[0].kdata[ti], "Polariz") != NULL)
			isPolarize = 1;
	}
	if ((ti = cmy->find_field(cmy, 0, "DEVCALSTD")) >= 0) {
		calstd = str2xcalstd(cmy->t[0].kdata[ti]);
	}
	
	if ((f_id1 = cmy->find_field(cmy, 0, "SampleName")) < 0
	 && (f_id1 = cmy->find_field(cmy, 0, "Sample_Name")) < 0
	 && (f_id1 = cmy->find_field(cmy, 0, "SAMPLE_NAME")) < 0
	 && (f_id1 = cmy->find_field(cmy, 0, "SAMPLE_ID")) < 0
	 && (f_id1 = cmy->find_field(cmy, 0, "SampleID")) < 0) {
		warning("Input file '%s' doesn't contain field SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SampleID",devname);
	} else {
		if (cmy->t[0].ftype[f_id1] != nqcs_t
		 && cmy->t[0].ftype[f_id1] != cs_t
		 && cmy->t[0].ftype[f_id1] != i_t)
			error("Field SampleName (%s) from CMYK/RGB file '%s' is wrong type",cmy->t[0].fsym[f_id1],devname);
	}

	if (cmy->find_field(cmy, 0, "RGB_R") >= 0) {
		ndchan = 3;
		if (verb)  {
			if (inp || disp)
				printf("Seems to be an RGB device\n");
			else
				printf("Seems to be a psuedo-RGB device\n");
		}
	} else if (cmy->find_field(cmy, 0, "CMYK_C") >= 0) {
		ndchan = 4;
		if (verb)
			printf("Seems to be a CMYK device\n");
	} else {
		printf("No device values found - hope that's OK!\n");
	}

	if (ndchan == 3) {
		if ((f_c = cmy->find_field(cmy, 0, "RGB_R")) < 0) {
			error("Input file '%s' doesn't contain field RGB_R",devname);
		}
		if (cmy->t[0].ftype[f_c] != r_t)
			error("Field RGB_R from file '%s' is wrong type - expect float",devname);

		if ((f_m = cmy->find_field(cmy, 0, "RGB_G")) < 0)
			error("Input file '%s' doesn't contain field RGB_G",devname);
		if (cmy->t[0].ftype[f_m] != r_t)
			error("Field RGB_G from file '%s' is wrong type - expect float",devname);

		if ((f_y = cmy->find_field(cmy, 0, "RGB_B")) < 0)
			error("Input file '%s' doesn't contain field RGB_B",devname);
		if (cmy->t[0].ftype[f_y] != r_t)
			error("Field RGB_B from file '%s' is wrong type - expect float",devname);

	} else if (ndchan == 4) {
		if ((f_c = cmy->find_field(cmy, 0, "CMYK_C")) < 0) {
			error("Input file '%s' doesn't contain field CMYK_C",devname);
		}
		if (cmy->t[0].ftype[f_c] != r_t)
			error("Field CMYK_C from file '%s' is wrong type - expect float",devname);

		if ((f_m = cmy->find_field(cmy, 0, "CMYK_M")) < 0)
			error("Input file '%s' doesn't contain field CMYK_M",devname);
		if (cmy->t[0].ftype[f_m] != r_t)
			error("Field CMYK_M from file '%s' is wrong type - expect float",devname);

		if ((f_y = cmy->find_field(cmy, 0, "CMYK_Y")) < 0)
			error("Input file '%s' doesn't contain field CMYK_Y",devname);
		if (cmy->t[0].ftype[f_y] != r_t)
			error("Field CMYK_Y from file '%s' is wrong type - expect float",devname);

		if ((f_k = cmy->find_field(cmy, 0, "CMYK_K")) < 0)
			error("Input file '%s' doesn't contain field CMYK_Y",devname);
		if (cmy->t[0].ftype[f_k] != r_t)
			error("Field CMYK_K from file '%s' is wrong type - expect float",devname);
	}
	if (verb && ndchan > 0) printf("Read device values\n");

	if (cmy->find_field(cmy, 0, "XYZ_X") >= 0
	 || cmy->find_field(cmy, 0, "LAB_L") >= 0
	 || cmy->find_field(cmy, 0, "Lab_L") >= 0) {
		/* We've got a new combined device+cie file as the first one. */
		/* Shuffle it into ciename , and ciename into specname */

		strcpy(specname, ciename);
		strcpy(ciename, devname);

		if (verb) printf("We've got a combined device + instrument readings file\n");
	}

	/* Open up the input nCIE or Spectral device data file */
	ncie = new_cgats();	/* Create a CGATS structure */
	ncie->add_other(ncie, "LGOROWLENGTH"); 	/* Gretag/Logo Target file */
	ncie->add_other(ncie, "ECI2002"); 		/* Gretag/Logo Target file */
	ncie->add_other(ncie, ""); 				/* Wildcard */
	if (ncie->read_name(ncie, ciename))
		error ("Read: Can't read cie file '%s'. Unknown format, missing or corrupted file ?",ciename);
	if (ncie->ntables != 1)
		warning("Input file '%s' doesn't contain exactly one table",ciename);

	if (npat != ncie->t[0].nsets)
		error("Number of patches between '%s' and '%s' doesn't match",devname,ciename);

	/* See if we can figure out polarization and xrga */
	if ((ti = ncie->find_field(ncie, 0, "MEASUREMENT_SOURCE")) >= 0) {
		if (strstr(ncie->t[0].kdata[ti], "Filter=Polarizer") != NULL)
			isPolarize = 1;
	}
	if ((ti = ncie->find_field(ncie, 0, "FILTER")) >= 0) {
		if (strstr(ncie->t[0].kdata[ti], "Polariz") != NULL)
			isPolarize = 1;
	}
	if ((ti = ncie->find_field(ncie, 0, "DEVCALSTD")) >= 0) {
		calstd = str2xcalstd(ncie->t[0].kdata[ti]);
	}

	if ((f_id2 = ncie->find_field(ncie, 0, "SampleName")) < 0
	 && (f_id2 = ncie->find_field(ncie, 0, "Sample_Name")) < 0
	 && (f_id2 = ncie->find_field(ncie, 0, "SAMPLE_NAME")) < 0
	 && (f_id2 = ncie->find_field(ncie, 0, "SAMPLE_ID")) < 0
	 && (f_id2 = ncie->find_field(ncie, 0, "SampleID")) < 0) {
		if (strcmp(devname, ciename))
			warning("Input file '%s' doesn't contain field SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SampleID",ciename);
	} else {
		if (ncie->t[0].ftype[f_id2] != nqcs_t
		 && ncie->t[0].ftype[f_id2] != cs_t
		 && ncie->t[0].ftype[f_id2] != i_t)
			error("Field SampleName (%s) from cie file '%s' is wrong type",ncie->t[0].fsym[f_id2],ciename);
	}

	if (ncie->find_field(ncie, 0, "XYZ_X") < 0
	 && ncie->find_field(ncie, 0, "LAB_L") < 0
	 && ncie->find_field(ncie, 0, "Lab_L") < 0) {

		/* Not a cie file. See if it's a spectral file */
		if (ncie->find_field(ncie, 0, "nm500") < 0
		 && ncie->find_field(ncie, 0, "NM_500") < 0
		 && ncie->find_field(ncie, 0, "SPECTRAL_NM_500") < 0
		 && ncie->find_field(ncie, 0, "R_500") < 0
		 && ncie->find_field(ncie, 0, "SPECTRAL_500") < 0
		 && ncie->find_field(ncie, 0, "SPECTRAL_NM500") < 0)
			error("Input file '%s' doesn't contain field XYZ_X, LAB_L or spectral",ciename);	/* Nope */

		/* We have a spectral file only. Fix things and drop through */
		ncie->del(ncie);
		ncie = NULL;
		strcpy(specname, ciename);
		ciename[0] = '\000';

	} else {	/* Continue dealing with cie value file */
		int fix = 0;
		char *fields[3][3] = {
			{ "XYZ_X", "XYZ_Y", "XYZ_Z" },
			{ "LAB_L", "LAB_A", "LAB_B" },
			{ "Lab_L", "Lab_a", "Lab_b" }
		};

		if (ncie->find_field(ncie, 0, "nm500") >= 0
		 || ncie->find_field(ncie, 0, "NM_500") < 0
		 || ncie->find_field(ncie, 0, "SPECTRAL_NM_500") >= 0
		 || ncie->find_field(ncie, 0, "R_500") >= 0
		 || ncie->find_field(ncie, 0, "SPECTRAL_500") >= 0
		 || ncie->find_field(ncie, 0, "SPECTRAL_NM500") >= 0) {
			if (verb) printf("Found spectral values\n");
			/* It's got spectral data too. Make sure we read it */
			strcpy(specname, ciename);
		}

		if (ncie->find_field(ncie, 0, "LAB_L") >= 0) {
			fix = 1;
			islab = 1;
		} else if (ncie->find_field(ncie, 0, "Lab_L") >= 0) {
			fix = 2;
			islab = 1;
		}

		for (i = 0; i < 3; i++) {

			if ((f_cie[i] = ncie->find_field(ncie, 0, fields[fix][i])) < 0
			 && (f_cie[i] = ncie->find_field(ncie, 0, fields[fix+1][i])) < 0)
				error("Input file '%s' doesn't contain field %s",ciename,fields[fix][i]);

			if (ncie->t[0].ftype[f_cie[i]] != r_t)
				error("Field %s from file '%s' is wrong type - expect float",fields[fix][i], ciename);
		}
	
		if (verb) printf("Found CIE values\n");
	}

	/* Open up the input Spectral device data file */
	if (specname[0] != '\000') {
		char bufs[6][50];

		spec = new_cgats();	/* Create a CGATS structure */
		spec->add_other(spec, "LGOROWLENGTH"); 	/* Gretag/Logo Target file */
		spec->add_other(spec, "ECI2002"); 		/* Gretag/Logo Target file */
		spec->add_other(spec, ""); 				/* Wildcard */
		if (spec->read_name(spec, specname))
			error ("Read: Can't read spec file '%s'. Unknown format, missing or corrupted file ?",specname);
		if (spec->ntables != 1)
			warning("Input file '%s' doesn't contain exactly one table",specname);

		if (npat != spec->t[0].nsets)
			error("Number of patches between '%s' and '%s' doesn't match",specname);

		/* See if we can figure out polarization and xrga */
		if ((ti = spec->find_field(spec, 0, "MEASUREMENT_SOURCE")) >= 0) {
			if (strstr(spec->t[0].kdata[ti], "Filter=Polarizer") != NULL)
				isPolarize = 1;
		}
		if ((ti = spec->find_field(spec, 0, "FILTER")) >= 0) {
			if (strstr(spec->t[0].kdata[ti], "Polariz") != NULL)
				isPolarize = 1;
		}
		if ((ti = spec->find_field(spec, 0, "DEVCALSTD")) >= 0) {
			calstd = str2xcalstd(spec->t[0].kdata[ti]);
		}
	
		if ((f_id3 = spec->find_field(spec, 0, "SampleName")) < 0
		 && (f_id3 = spec->find_field(spec, 0, "Sample_Name")) < 0
		 && (f_id3 = spec->find_field(spec, 0, "SAMPLE_NAME")) < 0
		 && (f_id3 = spec->find_field(spec, 0, "SAMPLE_ID")) < 0
		 && (f_id3 = spec->find_field(spec, 0, "SampleID")) < 0) {
			if (strcmp(devname, specname))
				warning("Input file '%s' doesn't contain field SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SampleID",specname);
		} else {
			if (spec->t[0].ftype[f_id3] != nqcs_t
			 && spec->t[0].ftype[f_id3] != cs_t
			 && spec->t[0].ftype[f_id3] != i_t)
				error("Field SampleName (%s) from spec file '%s' is wrong type",spec->t[0].fsym[f_id3],specname);
		}

		/* Find the spectral readings nm range */
		for (specmin = 500; specmin >= 300; specmin -= 10) {
			sprintf(bufs[0],"nm%03d", specmin);
			sprintf(bufs[1],"NM_%03d", specmin);
			sprintf(bufs[2],"SPECTRAL_NM_%03d", specmin);
			sprintf(bufs[3],"R_%03d", specmin);
			sprintf(bufs[4],"SPECTRAL_%03d", specmin);
			sprintf(bufs[5],"SPECTRAL_NM%03d", specmin);

			if (spec->find_field(spec, 0, bufs[0]) < 0
			 && spec->find_field(spec, 0, bufs[1]) < 0
			 && spec->find_field(spec, 0, bufs[2]) < 0
			 && spec->find_field(spec, 0, bufs[3]) < 0
			 && spec->find_field(spec, 0, bufs[4]) < 0
			 && spec->find_field(spec, 0, bufs[5]) < 0)	/* Not found */
			break;
		}
		specmin += 10;
		for (specmax = 500; specmax <= 900; specmax += 10) {
			sprintf(bufs[0],"nm%03d", specmax);
			sprintf(bufs[1],"NM_%03d", specmax);
			sprintf(bufs[2],"SPECTRAL_NM_%03d", specmax);
			sprintf(bufs[3],"R_%03d", specmax);
			sprintf(bufs[4],"SPECTRAL_%03d", specmax);
			sprintf(bufs[5],"SPECTRAL_NM%03d", specmax);

			if (spec->find_field(spec, 0, bufs[0]) < 0
			 && spec->find_field(spec, 0, bufs[1]) < 0
			 && spec->find_field(spec, 0, bufs[2]) < 0
			 && spec->find_field(spec, 0, bufs[3]) < 0
			 && spec->find_field(spec, 0, bufs[4]) < 0
			 && spec->find_field(spec, 0, bufs[5]) < 0) /* Not found */
			break;
		}
		specmax -= 10;
				
		if (specmin > 420 || specmax < 680) {	/* Not enough range to be useful */
			spec->del(spec);
			spec = NULL;
			specname[0] = '\000';
			if (verb)
				printf("Not enough (%d) spectral values - discarding spectral\n",specnum);
		} else {

			specnum = (specmax - specmin)/10 + 1;

			if (verb)
				printf("Found there are %d spectral values, from %d to %d nm\n",specnum,specmin,specmax);
			

			/* Locate the fields for spectral values */
			for (j = 0; j < specnum; j++) {
				sprintf(bufs[0],"nm%03d", specmin + 10 * j);
				sprintf(bufs[1],"NM_%03d", specmin + 10 * j);
				sprintf(bufs[2],"SPECTRAL_NM_%03d", specmin + 10 * j);
				sprintf(bufs[3],"R_%03d", specmin + 10 * j);
				sprintf(bufs[4],"SPECTRAL_%03d", specmin + 10 * j);
				sprintf(bufs[5],"SPECTRAL_NM%03d", specmin + 10 * j);
	
				if ((spi[j] = spec->find_field(spec, 0, bufs[0])) < 0
				 && (spi[j] = spec->find_field(spec, 0, bufs[1])) < 0
				 && (spi[j] = spec->find_field(spec, 0, bufs[2])) < 0
				 && (spi[j] = spec->find_field(spec, 0, bufs[3])) < 0
				 && (spi[j] = spec->find_field(spec, 0, bufs[4])) < 0
				 && (spi[j] = spec->find_field(spec, 0, bufs[5])) < 0) {	/* Not found */
					
					spec->del(spec);
					spec = NULL;
					specname[0] = '\000';
					error("Failed to find spectral band %d nm in file '%s'\n",specmin + 10 * j,specname);
				} else {
					if (spec->t[0].ftype[spi[j]] != r_t)
						error("Field '%s' from file '%s' is wrong type - expect float",spec->t[0].fsym[spi[j]], specname);
				}

			}
		}
	}

	if (ciename[0] == '\000' && specname[0] == '\000')
		error("Input file doesn't contain either CIE or spectral data");

	/* Setup output cgats file */
	ocg = new_cgats();	/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI3"); 	/* our special type is Calibration Target Information 3 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll target", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	if (disp) {
		ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);	/* What sort of device this is */

		if (disp == 2)
			ocg->add_kword(ocg, 0, "NORMALIZED_TO_Y_100","NO", NULL);	/* What sort of device this is */
	}else if (inp)
		ocg->add_kword(ocg, 0, "DEVICE_CLASS","INPUT", NULL);	/* What sort of device this is */
	else
		ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */

	/* Note what instrument the chart was read with */
	/* Assume this - could try reading from file INSTRUMENTATION "SpectroScan" ?? */
	ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", inst_name(instSpectrolino) , NULL);

	if (calstd != xcalstd_none)
		ocg->add_kword(ocg, 0, "DEVCALSTD", xcalstd2str(calstd), NULL);

	if (isPolarize)
		ocg->add_kword(ocg, 0, "INSTRUMENT_FILTER", "POLARIZED", NULL);

	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
	if (f_id1 >= 0) 
		ocg->add_field(ocg, 0, "SAMPLE_LOC", cs_t);

	if (ndchan == 3) {
		ocg->add_field(ocg, 0, "RGB_R", r_t);
		ocg->add_field(ocg, 0, "RGB_G", r_t);
		ocg->add_field(ocg, 0, "RGB_B", r_t);
		if (inp) {
			if (islab)
				ocg->add_kword(ocg, 0, "COLOR_REP","LAB_RGB", NULL);
			else
				ocg->add_kword(ocg, 0, "COLOR_REP","XYZ_RGB", NULL);
		} else if (disp) {
			if (islab)
				ocg->add_kword(ocg, 0, "COLOR_REP","RGB_LAB", NULL);
			else
				ocg->add_kword(ocg, 0, "COLOR_REP","RGB_XYZ", NULL);
		} else {
			if (islab)
				ocg->add_kword(ocg, 0, "COLOR_REP","iRGB_LAB", NULL);
			else
				ocg->add_kword(ocg, 0, "COLOR_REP","iRGB_XYZ", NULL);
		}
	} else if (ndchan == 4) {
		ocg->add_field(ocg, 0, "CMYK_C", r_t);
		ocg->add_field(ocg, 0, "CMYK_M", r_t);
		ocg->add_field(ocg, 0, "CMYK_Y", r_t);
		ocg->add_field(ocg, 0, "CMYK_K", r_t);
		if (inp) {	/* Does this make any sense ? */
			if (islab)
				ocg->add_kword(ocg, 0, "COLOR_REP","LAB_CMYK", NULL);
			else
				ocg->add_kword(ocg, 0, "COLOR_REP","XYZ_CMYK", NULL);
		} else {
			if (islab)
				ocg->add_kword(ocg, 0, "COLOR_REP","CMYK_LAB", NULL);
			else
				ocg->add_kword(ocg, 0, "COLOR_REP","CMYK_XYZ", NULL);
		}
	}

	if (ncie != NULL) {
		if (islab) {
			ocg->add_field(ocg, 0, "LAB_L", r_t);
			ocg->add_field(ocg, 0, "LAB_A", r_t);
			ocg->add_field(ocg, 0, "LAB_B", r_t);
		} else {
			ocg->add_field(ocg, 0, "XYZ_X", r_t);
			ocg->add_field(ocg, 0, "XYZ_Y", r_t);
			ocg->add_field(ocg, 0, "XYZ_Z", r_t);
		}
	}

	/* Guess the device data scaling */
	if (ndchan > 0) {
		double maxv = 0.0;
		int f_dev[4] = { f_c, f_m, f_y, f_k };

		/* Guess what scale the spectral data is set to */
		for (i = 0; i < npat; i++) {
			for (j = 0; j < ndchan; j++) {
				double vv;
				vv = *((double *)cmy->t[0].fdata[i][f_dev[j]]);
				if (vv > maxv)
					maxv = vv;
			}
		}
		if (maxv < 10.0) {
			dev_scale = 100.0/1.0;
			if (verb) printf("Device max found = %f, scale by 100.0\n",maxv);
		} else if (maxv > 160.0) {
			dev_scale = 100.0/255.0;
			if (verb) printf("Device max found = %f, scale by 100/255\n",maxv);
		} else {
			dev_scale = 100.0/100.0;
			if (verb) printf("Device max found = %f, scale by 1.0\n",maxv);
		}
	}

	if (spec != NULL) {
		char buf[100];
		double maxv = 0.0;
		sprintf(buf,"%d", specnum);
		ocg->add_kword(ocg, 0, "SPECTRAL_BANDS",buf, NULL);
		sprintf(buf,"%d", specmin);
		ocg->add_kword(ocg, 0, "SPECTRAL_START_NM",buf, NULL);
		sprintf(buf,"%d", specmax);
		ocg->add_kword(ocg, 0, "SPECTRAL_END_NM",buf, NULL);

		/* Generate fields for spectral values */
		for (j = 0; j < specnum; j++) {
			sprintf(buf,"SPEC_%03d", specmin + 10 * j);
			ocg->add_field(ocg, 0, buf, r_t);
		}

		/* Guess what scale the spectral data is set to */
		for (i = 0; i < npat; i++) {
			for (j = 0; j < specnum; j++) {
				double vv;
				vv = *((double *)spec->t[0].fdata[i][spi[j]]);
				if (vv > maxv)
					maxv = vv;
			}
		}
		if (maxv < 10.0) {
			spec_scale = 100.0/1.0;
			if (verb) printf("Spectral max found = %f, scale by 100.0\n",maxv);
		} else if (maxv > 160.0) {
			spec_scale = 100.0/255.0;
			if (verb) printf("Spectral max found = %f, scale by 100/255\n",maxv);
		} else {
			spec_scale = 100.0/100.0;
			if (verb) printf("Spectral max found = %f, scale by 1.0\n",maxv);
		}
	}

	/* Write out the patch info to the output CGATS file */
	{
		cgats_set_elem *setel;	/* Array of set value elements */

		if ((setel = (cgats_set_elem *)malloc(
		     sizeof(cgats_set_elem) * ocg->t[0].nfields)) == NULL)
			error("Malloc failed!");

		if (ncie != NULL && strcmp(ciename, devname)) {
			if (f_id1 < 0)
				error("Input file '%s' needs SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SampleID to synchronize samples",devname);
			if (f_id2 < 0)
				error("Input file '%s' needs SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SampleID to synchronize samples",ciename);
		}

		if (spec != NULL && strcmp(specname, devname)) {
			if (f_id1 < 0)
				error("Input file '%s' needs SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SampleID to synchronize samples",devname);
			if (f_id3 < 0)
				error("Input file '%s' needs SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SampleID to synchronize samples",specname);
		}


		/* Write out the patch info to the output CGATS file */
		for (i = 0; i < npat; i++) {
			char id[100];
			char loc[100];
			int k = 0;

			if (ncie != NULL && f_id1 >= 0 && f_id1 >= 0) {
				if (strcmp(((char *)cmy->t[0].rfdata[i][f_id1]), 
				           ((char *)ncie->t[0].rfdata[i][f_id2])) != 0) {
					error("Patch label mismatch to CIE values, patch %d, '%s' != '%s'\n",
					       i, ((char *)cmy->t[0].rfdata[i][f_id1]), 
				              ((char *)ncie->t[0].rfdata[i][f_id2]));
				}
			}

			if (spec != NULL && f_id1 >= 0 && f_id3 >= 0) {
				if (strcmp(((char *)cmy->t[0].rfdata[i][f_id1]), 
				           ((char *)spec->t[0].rfdata[i][f_id3])) != 0) {
					error("Patch label mismatch to spectral values, patch %d, '%s' != '%s'\n",
					       i, ((char *)cmy->t[0].rfdata[i][f_id1]), 
				              ((char *)spec->t[0].rfdata[i][f_id3]));
				}
			}

			/* SAMPLE ID */
			sprintf(id, "%d", i+1);
			setel[k++].c = id;

			/* SAMPLE NAME */
			if (f_id1 >= 0) { 
				strcpy(loc, (char *)cmy->t[0].rfdata[i][f_id1]);

				/* Transpose Letters and Numbers */
				if (transpose)
					trans(loc, (char *)cmy->t[0].rfdata[i][f_id1]);

				setel[k++].c = loc;
			}
			
			if (ndchan == 3) {
				setel[k++].d = dev_scale * *((double *)cmy->t[0].fdata[i][f_c]);
				setel[k++].d = dev_scale * *((double *)cmy->t[0].fdata[i][f_m]);
				setel[k++].d = dev_scale * *((double *)cmy->t[0].fdata[i][f_y]);
			} else if (ndchan == 4){
				double sum = 0.0;
				sum += setel[k++].d = dev_scale * *((double *)cmy->t[0].fdata[i][f_c]);
				sum += setel[k++].d = dev_scale * *((double *)cmy->t[0].fdata[i][f_m]);
				sum += setel[k++].d = dev_scale * *((double *)cmy->t[0].fdata[i][f_y]);
				sum += setel[k++].d = dev_scale * *((double *)cmy->t[0].fdata[i][f_k]);
				if (sum > mxsum) {
					mxsum = sum;
					mxsumix = i;
				}
			}
	
			if (ncie != NULL) {
				setel[k++].d = *((double *)ncie->t[0].fdata[i][f_cie[0]]);	
				setel[k++].d = *((double *)ncie->t[0].fdata[i][f_cie[1]]);
				setel[k++].d = *((double *)ncie->t[0].fdata[i][f_cie[2]]);
			}
	
			if (spec) {
				for (j = 0; j < specnum; j++) {
					setel[k++].d = spec_scale * *((double *)spec->t[0].fdata[i][spi[j]]);
				}
			}

			ocg->add_setarr(ocg, 0, setel);
		}

		free(setel);
	}

	if (tlimit < 0 && mxsum > 0.0) {
		if (verb)
			printf("No ink limit given, using maximum %f found in file at %d\n",mxsum,mxsumix+1);

		tlimit = (int)(mxsum + 0.5);
	}

	if (tlimit > 0) {
		char buf[100];
		sprintf(buf, "%d", tlimit);
		ocg->add_kword(ocg, 0, "TOTAL_INK_LIMIT", buf, NULL);
	}

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	/* Create a dummy .ti2 file (used with scanin -r) */
	if (out2) {

		/* Setup output cgats file */
		ocg2 = new_cgats();	/* Create a CGATS structure */
		ocg2->add_other(ocg2, "CTI2"); 	/* our special type is Calibration Target Information 2 */
		ocg2->add_table(ocg2, tt_other, 0);	/* Start the first table */

		ocg2->add_kword(ocg2, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 2",NULL);
		ocg2->add_kword(ocg2, 0, "ORIGINATOR", "Argyll txt2ti3", NULL);
		atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
		ocg2->add_kword(ocg2, 0, "CREATED",atm, NULL);
		if (disp)
			ocg2->add_kword(ocg2, 0, "DEVICE_CLASS","DISPLAY", NULL);	/* What sort of device this is */
		else
			ocg2->add_kword(ocg2, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */
		/* Note what instrument the chart was read with */
		/* Assume this - could try reading from file INSTRUMENTATION "SpectroScan" ?? */
		ocg2->add_kword(ocg2, 0, "TARGET_INSTRUMENT", inst_name(instSpectrolino) , NULL);

		/* Fields we want */
		ocg2->add_field(ocg2, 0, "SAMPLE_ID", nqcs_t);
		ocg2->add_field(ocg2, 0, "SAMPLE_LOC", nqcs_t);

		/* We're missing lots of .ti2 stuff like: */
		/* ocg->add_kword(ocg, 0, "APPROX_WHITE_POINT",icg->t[0].kdata[fi], NULL); */
		/* ocg->add_kword(ocg, 0, "PATCH_LENGTH", buf, NULL); */
		/* ocg->add_kword(ocg, 0, "GAP_LENGTH", buf, NULL); */
		/* ocg->add_kword(ocg, 0, "TRAILER_LENGTH", buf, NULL); */
		/* ocg->add_kword(ocg, 0, "STEPS_IN_PASS", buf, NULL); */
		/* ocg->add_kword(ocg, 0, "PASSES_IN_STRIPS", pis, NULL); */
		/* ocg->add_kword(ocg, 0, "STRIP_INDEX_PATTERN", sixpat, NULL); */
		/* ocg->add_kword(ocg, 0, "PATCH_INDEX_PATTERN", pixpat, NULL); */
		/* ocg->add_kword(ocg, 0, "INDEX_ORDER", ixord ? "PATCH_THEN_STRIP" : "STRIP_THEN_PATCH", NULL); */

		if (tlimit > 0) {
			char buf[100];
			sprintf(buf, "%d", tlimit);
			ocg2->add_kword(ocg2, 0, "TOTAL_INK_LIMIT", buf, NULL);
		}

		if (ndchan == 3) {
			ocg2->add_field(ocg2, 0, "RGB_R", r_t);
			ocg2->add_field(ocg2, 0, "RGB_G", r_t);
			ocg2->add_field(ocg2, 0, "RGB_B", r_t);
			if (inp || disp) {
				ocg2->add_kword(ocg2, 0, "COLOR_REP","RGB", NULL);
			} else {
				ocg2->add_kword(ocg2, 0, "COLOR_REP","iRGB", NULL);
			}
		} else if (ndchan == 4) {
			ocg2->add_field(ocg2, 0, "CMYK_C", r_t);
			ocg2->add_field(ocg2, 0, "CMYK_M", r_t);
			ocg2->add_field(ocg2, 0, "CMYK_Y", r_t);
			ocg2->add_field(ocg2, 0, "CMYK_K", r_t);
			ocg2->add_kword(ocg2, 0, "COLOR_REP","CMYK", NULL);
		}

		if (ncie != NULL) {
			if (islab) {
				ocg2->add_field(ocg2, 0, "LAB_L", r_t);
				ocg2->add_field(ocg2, 0, "LAB_A", r_t);
				ocg2->add_field(ocg2, 0, "LAB_B", r_t);
			} else {
				ocg2->add_field(ocg2, 0, "XYZ_X", r_t);
				ocg2->add_field(ocg2, 0, "XYZ_Y", r_t);
				ocg2->add_field(ocg2, 0, "XYZ_Z", r_t);
			}
		}

		if (f_id1 < 0)
			error("Input file '%s' needs SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SampleID to be able to generate .ti2 file",devname);

		/* Write out the patch info to the output CGATS file */
		{
			cgats_set_elem *setel;	/* Array of set value elements */

			if ((setel = (cgats_set_elem *)malloc(
			     sizeof(cgats_set_elem) * (2 + (ndchan) + (ncie != NULL ? 3 : 0)))) == NULL)
				error("Malloc failed!");

			/* Write out the patch info to the output CGATS file */
			for (i = 0; i < npat; i++) {
				char id[100];
				int k = 0;

				if (ncie != NULL) {
					if (strcmp(((char *)cmy->t[0].rfdata[i][f_id1]), 
					           ((char *)ncie->t[0].rfdata[i][f_id2])) != 0) {
						error("Patch label mismatch to CIE values, patch %d, '%s' != '%s'\n",
						       i, ((char *)cmy->t[0].rfdata[i][f_id1]), 
					              ((char *)ncie->t[0].rfdata[i][f_id2]));
					}
				}

				if (spec != NULL) {
					if (strcmp(((char *)cmy->t[0].rfdata[i][f_id1]), 
					           ((char *)spec->t[0].rfdata[i][f_id3])) != 0) {
						error("Patch label mismatch to spectral values, patch %d, '%s' != '%s'\n",
						       i, ((char *)cmy->t[0].rfdata[i][f_id1]), 
					              ((char *)spec->t[0].rfdata[i][f_id3]));
					}
				}

				sprintf(id, "%d", i+1);
				setel[k++].c = id;						/* ID */
				setel[k++].c = ((char *)cmy->t[0].rfdata[i][f_id1]); 	/* Location */

				if (ndchan == 3) {
					setel[k++].d = 100.0/255.0 * *((double *)cmy->t[0].fdata[i][f_c]);
					setel[k++].d = 100.0/255.0 * *((double *)cmy->t[0].fdata[i][f_m]);
					setel[k++].d = 100.0/255.0 * *((double *)cmy->t[0].fdata[i][f_y]);
				} else if (ndchan == 4) {
					setel[k++].d = *((double *)cmy->t[0].fdata[i][f_c]);
					setel[k++].d = *((double *)cmy->t[0].fdata[i][f_m]);
					setel[k++].d = *((double *)cmy->t[0].fdata[i][f_y]);
					setel[k++].d = *((double *)cmy->t[0].fdata[i][f_k]);
				}
		
				if (ncie != NULL) {
					setel[k++].d = *((double *)ncie->t[0].fdata[i][f_cie[0]]);	
					setel[k++].d = *((double *)ncie->t[0].fdata[i][f_cie[1]]);
					setel[k++].d = *((double *)ncie->t[0].fdata[i][f_cie[2]]);
				}
				ocg2->add_setarr(ocg2, 0, setel);
			}

			free(setel);
		}

		if (ocg2->write_name(ocg2, outname2))
			error("Write error : %s",ocg2->err);

	}

	/* Clean up */
	cmy->del(cmy);
	if (ncie != NULL)
		ncie->del(ncie);
	if (spec != NULL)
		spec->del(spec);
	ocg->del(ocg);

	return 0;
}


/* Transpose Letters and Numbers of location */
/* return nz on error */
static int trans(char *dst, char *src) {
	int tt = 0, first = 0, second = 0;		/* 1..n */
	char *c, *d;		/* Dividing point */

//printf("\nGot '%s'\n",src);
	d = src;

	if (*d == '\000')
		return 1; 

	if ((*d >= 'A' && *d <= 'Z')
	 || (*d >= 'a' && *d <= 'z'))
		tt = 1;						// Initial is letters
	else
		tt = 0;
//printf("Initial is %s\n",tt ? "letters" : "numbers");

	for (; ; d++) {
		if (*d == '\000') {
//printf("Failed to find division\n");
			return 1;
		}
		if (tt && *d >= '0' && *d <= '9') {
//printf("Found number at offset %d\n",d - src);
			break;

		} else if (!tt && ((*d >= 'A' && *d <= 'Z')
			            || (*d >= 'a' && *d <= 'z'))) {
//printf("Found letter at offset %d\n",d - src);
			break;
		}
	}


//printf("2nd seg = '%s'\n",d);

	for (c = src; c < d; c++) {
		if (tt) {
			if (*c >= 'A' && *c <= 'Z')
				first = first * 26 + (*c - 'A' + 1);
			else
				first = first * 26 + (*c - 'a' + 1);
		} else {
				first = first * 10 + (*c - '0');
		}
	}
	for (; *c != '\000'; c++) {
		if (!tt) {
			if (*c >= 'A' && *c <= 'Z')
				second = second * 26 + (*c - 'A' + 1);
			else
				second = second * 26 + (*c - 'a' + 1);
		} else {
				second = second * 10 + (*c - '0');
		}
	}
//printf("First = %d, second = %d\n",first,second);

	c = dst;
	if (tt) {		/* Output letters */
		if (second > 26) {
			*c++ = 'A' + ((second-1) / 26) -1;
			second = ((second-1) % 26) + 1;
		}
		*c++ = 'A' + second -1;
	} else {		/* Output letters */
		c += sprintf(c, "%d",second);
	}
	if (tt) {		/* Output letters */
		c += sprintf(c, "%d",first);
	} else {		/* Output numbers */
		if (first > 26) {
			*c++ = 'A' + ((first-1) / 26) -1;
			first = ((first-1) % 26) + 1;
		}
		*c++ = 'A' + first -1;
	}
	*c = '\000';
//printf("Output '%s'\n",dst);

	return 0;
}

