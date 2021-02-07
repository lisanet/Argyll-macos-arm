
/* 
 * Argyll Color Correction System
 * Color Device profile checker.
 *
 * Author: Graeme W. Gill
 * Date:   15/7/2001
 *
 * Copyright 2001 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in the .ti3 scattered test chart
 * points, and checks them against an ICC profile.
 * forward ICC device profile.
 */

/*
 * TTBD:
 *		Switch to generic colorant read code rather than Grey/RGB/CMYK,
 *		and allow checking ICC profiles > 4 colors
 */

#undef DEBUG

#undef HACK					/* Print per patch XYZ differences */

#define IMP_MONO			/* Turn on development code */

#define verbo stdout

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "cgats.h"
#include "xicc.h"
#include "insttypes.h"
#include "sort.h"
#include "plot.h"
#include "vrml.h"
#include "ui.h"

void
usage(void) {
	fprintf(stderr,"Check accuracy of ICC profile, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: profcheck [-options] data.ti3 iccprofile.icm\n");
	fprintf(stderr," -v [level]      Verbosity level (default 1), 2 to print each DE\n");
	fprintf(stderr," -c              Show CIE94 delta E values\n");
	fprintf(stderr," -k              Show CIEDE2000 delta E values\n");
	fprintf(stderr," -w              create %s visualisation (iccprofile%s)\n",vrml_format(),vrml_ext());
	fprintf(stderr," -x              Use %s axes\n",vrml_format());
	fprintf(stderr," -m              Make %s lines a minimum of 0.5\n",vrml_format());
	fprintf(stderr," -e              Color vectors acording to delta E\n");
	fprintf(stderr," -h              Plot a histogram of delta E's\n");
	fprintf(stderr," -s              Sort output by delta E\n");
	fprintf(stderr," -P N.NN         Create a pruned .ti3 with points less or equal to N.NN delta E\n");
	fprintf(stderr," -d devval1,deval2,devvalN\n");
	fprintf(stderr,"                 Specify a device value to sort against\n");
	fprintf(stderr," -p              Sort device value by PCS (Lab) target\n");
	fprintf(stderr," -f [illum]      Use Fluorescent Whitening Agent compensation [opt. simulated inst. illum.:\n");
	fprintf(stderr,"                  M0, M1, M2, A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp]\n");
	fprintf(stderr," -i illum        Choose illuminant for computation of CIE XYZ from spectral data & FWA:\n");
	fprintf(stderr,"                  A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                  1931_2 (def), 1964_10, 2012_2, 2012_10, S&B 1955_2, shaw, J&V 1978_2 or file.cmf\n");
	fprintf(stderr," -I intent       r = relative colorimetric, a = absolute (default)\n");
	fprintf(stderr," data.ti3        Test data file\n");
	fprintf(stderr," iccprofile.icm  Profile to check against\n");
	exit(1);
	}

static void DE2RGB(double *out, double in);

/* Patch value type */
typedef struct {
	char sid[50];		/* sample id */
	char slo[50];		/* sample location, "" if not known */
	double p[MAX_CHAN];	/* Device value */
	double v[3];		/* CIE value */

	double pv[3];		/* Profile CIE value */
	double de;			/* Delta E to profile CIE value */

	double dp;			/* Delta p[] from target device value */
	double dv;			/* Delta E from CIE value */
} pval;

/* Histogram bin type */
typedef struct {
	int count;			/* Raw count */
	double val;			/* Normalized value */
	double min, max;	/* Bin range */
} hbin;

int main(int argc, char *argv[])
{
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int cie94 = 0;
	int cie2k = 0;
	int dovrml = 0;
	int dominl = 0;
	int doaxes = 0;
	int dodecol = 0;
	char ti3name[MAXNAMEL+1] = { 0 };	/* Input cgats file base name */
	cgats *icg;				/* input cgats structure */
	char iccname[MAXNAMEL+1] = { 0 };	/* Input icc file base name */
	icmFile *rd_fp;
	icRenderingIntent intent = icAbsoluteColorimetric;
	icc *rd_icco;
	icmLuBase *luo;
	char out_name[MAXNAMEL+1], *xl;		/* VRML/X3D name */
	vrml *wrl = NULL;

	int fwacomp = 0;			/* FWA compensation */
	int isemis = 0;				/* nz if this has emissive reference data */
	int isdisp = 0;				/* nz if this is a display device, 0 if output */
	int isdnormed = 0;      	/* Has display data been normalised to 100 ? */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType tillum = icxIT_none;	/* Target/simulated instrument illuminant */ 
	xspect cust_tillum, *tillump = NULL; /* Custom target/simulated illumination spectrum */
	icxIllumeType illum = icxIT_none;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType obType = icxOT_none;
	xspect custObserver[3];		/* If obType = icxOT_custom */

	int sortbyde = 0;			/* Sort by delta E */

	int ddevv = 0;				/* Do device value sort */
	double devval[MAX_CHAN];	/* device value to sort on */
	int sortbypcs = 0;			/* Sort by PCS error to device target */

	int dohisto = 0;			/* Plot histogram of delta E's */
	double prune = 0.0;			/* If > 0.0, created a pruned .ti3 file */

	int npat;					/* Number of patches */
	pval *tpat;					/* Patch input values */
	pval **stpat;				/* Pointers to internal sorted tpat[] */
	int i, j, rv = 0;
	icColorSpaceSignature devspace = 0;	/* The device colorspace */
	int isAdditive = 0;			/* 0 if subtractive, 1 if additive colorspace */
	int isLab = 0;				/* 0 if input is XYZ, 1 if input is Lab */
	int devchan = 0;			/* Number of device chanels */

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	error_program = "profcheck";

	if (argc <= 1)
		usage();

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
				usage();

			/* Verbosity */
			else if (argv[fa][1] == 'v') {
				verb = 1;
				if (na != NULL && isdigit(na[0])) {
					verb = atoi(na);
				}
			}

			/* VRML/X3D */
			else if (argv[fa][1] == 'w')
				dovrml = 1;

			/* Minimum line length */
			else if (argv[fa][1] == 'm')
				dominl = 1;

			/* Axes */
			else if (argv[fa][1] == 'x')
				doaxes = 1;

			/* Delta E coloring */
			else if (argv[fa][1] == 'e')
				dodecol = 1;

			else if (argv[fa][1] == 'c') {
				cie94 = 1;
				cie2k = 0;
			}

			else if (argv[fa][1] == 'k') {
				cie94 = 0;
				cie2k = 1;
			}

			/* Plot histogram */
			else if (argv[fa][1] == 'h')
				dohisto = 1;

			/* Sort by delta E */
			else if (argv[fa][1] == 's')
				sortbyde = 1;

			/* Create a pruned .ti3 */
			else if (argv[fa][1] == 'P') {
				if (na == NULL) usage();
				fa = nfa;
				prune = atof(na);
				if (prune <= 0.0)
					usage();
			}

			/* Device sort value */
			else if (argv[fa][1] == 'd') {
				char *tp, buf[200];
				int ndv;
				if (na == NULL) usage();
				fa = nfa;

				ddevv = 1;
				strcpy(buf, na);
				
				/* Replace ',' with '\000' */
				for (ndv = 1,tp = buf; *tp != '\000'; tp++) {
					if (*tp == ',') {
						*tp = '\000';
						ndv++;
					}
				}
				if (ndv >= MAX_CHAN)
					ndv = MAX_CHAN;

				for (tp = buf, i = 0; i < ndv; i++, tp += strlen(tp) + 1) {
					devval[i] = atof(tp);
				}
			}

			else if (argv[fa][1] == 'p')
				sortbypcs = 1;

			/* FWA compensation */
			else if (argv[fa][1] == 'f') {
				fwacomp = 1;

				if (na != NULL) {	/* Argument is present - target/simulated instr. illum. */
					fa = nfa;
					if (strcmp(na, "A") == 0
					 || strcmp(na, "M0") == 0) {
						spec = 1;
						tillum = icxIT_A;
					} else if (strcmp(na, "C") == 0) {
						spec = 1;
						tillum = icxIT_C;
					} else if (strcmp(na, "D50") == 0
					        || strcmp(na, "M1") == 0) {
						spec = 1;
						tillum = icxIT_D50;
					} else if (strcmp(na, "D50M2") == 0
					        || strcmp(na, "M2") == 0) {
						spec = 1;
						tillum = icxIT_D50M2;
					} else if (strcmp(na, "D65") == 0) {
						spec = 1;
						tillum = icxIT_D65;
					} else if (strcmp(na, "F5") == 0) {
						spec = 1;
						tillum = icxIT_F5;
					} else if (strcmp(na, "F8") == 0) {
						spec = 1;
						tillum = icxIT_F8;
					} else if (strcmp(na, "F10") == 0) {
						spec = 1;
						tillum = icxIT_F10;
					} else {	/* Assume it's a filename */
						inst_meas_type mt;

						spec = 1;
						tillum = icxIT_custom;
						if (read_xspect(&cust_tillum, &mt, na) != 0)
							usage();

						if (mt != inst_mrt_none
						 && mt != inst_mrt_emission
						 && mt != inst_mrt_ambient
						 && mt != inst_mrt_emission_flash
						 && mt != inst_mrt_ambient_flash) {
							error("Target illuminant '%s' is wrong measurement type",na);
						}
					}
				}
			}
			/* Spectral Illuminant type */
			else if (argv[fa][1] == 'i') {
				if (na == NULL) usage();
				fa = nfa;
				if (strcmp(na, "A") == 0) {
					spec = 1;
					illum = icxIT_A;
				} else if (strcmp(na, "C") == 0) {
					spec = 1;
					illum = icxIT_C;
				} else if (strcmp(na, "D50") == 0) {
					spec = 1;
					illum = icxIT_D50;
				} else if (strcmp(na, "D50M2") == 0) {
					spec = 1;
					illum = icxIT_D50M2;
				} else if (strcmp(na, "D65") == 0) {
					spec = 1;
					illum = icxIT_D65;
				} else if (strcmp(na, "F5") == 0) {
					spec = 1;
					illum = icxIT_F5;
				} else if (strcmp(na, "F8") == 0) {
					spec = 1;
					illum = icxIT_F8;
				} else if (strcmp(na, "F10") == 0) {
					spec = 1;
					illum = icxIT_F10;
				} else {	/* Assume it's a filename */
					inst_meas_type mt;

					spec = 1;
					illum = icxIT_custom;
					if (read_xspect(&cust_illum, &mt, na) != 0)
						usage();

					if (mt != inst_mrt_none
					 && mt != inst_mrt_emission
					 && mt != inst_mrt_ambient
					 && mt != inst_mrt_emission_flash
					 && mt != inst_mrt_ambient_flash) {
						error("CIE illuminant '%s' is wrong measurement type",na);
					}
				}
			}

			/* Spectral Observer type */
			else if (argv[fa][1] == 'o') {
				if (na == NULL) usage();
				fa = nfa;
				if (strcmp(na, "1931_2") == 0) {			/* Classic 2 degree */
					spec = 1;
					obType = icxOT_CIE_1931_2;
				} else if (strcmp(na, "1964_10") == 0) {	/* Classic 10 degree */
					spec = 1;
					obType = icxOT_CIE_1964_10;
				} else if (strcmp(na, "2012_2") == 0) {		/* Latest 2 degree */
					spec = 1;
					obType = icxOT_CIE_2012_2;
				} else if (strcmp(na, "2012_10") == 0) {	/* Latest 10 degree */
					spec = 1;
					obType = icxOT_CIE_2012_10;
				} else if (strcmp(na, "1955_2") == 0) {		/* Stiles and Burch 1955 2 degree */
					spec = 1;
					obType = icxOT_Stiles_Burch_2;
				} else if (strcmp(na, "1978_2") == 0) {		/* Judd and Voss 1978 2 degree */
					spec = 1;
					obType = icxOT_Judd_Voss_2;
				} else if (strcmp(na, "shaw") == 0) {		/* Shaw and Fairchilds 1997 2 degree */
					spec = 1;
					obType = icxOT_Shaw_Fairchild_2;
				} else {	/* Assume it's a filename */
					obType = icxOT_custom;
					if (read_cmf(custObserver, na) != 0)
						usage();
				}
			}

			/* Intent (only applies to ICC profile) */
			else if (argv[fa][1] == 'I') {
				if (na == NULL) usage();
				fa = nfa;
    			switch (na[0]) {
					case 'r':
						intent = icRelativeColorimetric;
						break;
					case 'a':
						intent = icAbsoluteColorimetric;
						break;
					default:
						usage();
				}
			}

			else 
				usage();
		} else
			break;
	}

	/* Get the file name arguments */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(ti3name,argv[fa++],MAXNAMEL); ti3name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(iccname,argv[fa++],MAXNAMEL); iccname[MAXNAMEL] = '\000';

	strncpy(out_name,iccname,MAXNAMEL-4); out_name[MAXNAMEL-4] = '\000';
	if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
		xl = out_name + strlen(out_name);
	xl[0] = '\000';			/* Remove extension */

	if (fwacomp && spec == 0)
		error ("FWA compensation only works when viewer and/or illuminant selected");

	/* Open and look at the .ti3 profile patches file */
	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI3"); 	/* our special input type is Calibration Target Information 3 */

	if (icg->read_name(icg, ti3name))
		error("CGATS file read error on '%s': %s",ti3name,icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file '%s' isn't a CTI3 format file",ti3name);
	if (icg->ntables < 1)
		error ("Input file '%s' doesn't contain at least one table",ti3name);

	/* Figure out what sort of device it is */
	{
		int ti;

		if ((ti = icg->find_kword(icg, 0, "DEVICE_CLASS")) < 0)
			error ("Input file '%s' doesn't contain keyword DEVICE_CLASS",ti3name);

		if (strcmp(icg->t[0].kdata[ti],"EMISINPUT") == 0) {
			isemis = 1;

		} else if (strcmp(icg->t[0].kdata[ti],"DISPLAY") == 0) {
			isemis = 1;
			isdisp = 1;
		}

		/* See if the CIE data has been normalised to Y = 100 */
		if ((ti = icg->find_kword(icg, 0, "NORMALIZED_TO_Y_100")) < 0
		 || strcmp(icg->t[0].kdata[ti],"NO") == 0) {
			isdnormed = 0;
		} else {
			isdnormed = 1;
		}
	}

	if (isemis && illum != icxIT_none)
		warning("-i illuminant ignored for emissive reference type");

	if (isemis & fwacomp) {
			warning("-f FWA compensation ignored for emissive reference type");
		fwacomp = 0;
		tillum = icxIT_none;
	}

	/* Set defaults */
	if (illum == icxIT_none)
		illum = icxIT_D50;
	
	if (obType == icxOT_none)
		obType = icxOT_CIE_1931_2;

	/* See if CIE is actually available - some sources of .TI3 don't provide it */
	if (!spec
	 && icg->find_field(icg, 0, "LAB_L") < 0
	 && icg->find_field(icg, 0, "XYZ_X") < 0) {

		if (icg->find_kword(icg, 0, "SPECTRAL_BANDS") < 0)
			error ("Neither CIE nor spectral data found in file '%s'",ti3name);

		/* Switch to using spectral information */
		if (verb)
			printf("No CIE data found, switching to spectral with standard observer & D50\n");
		spec = 1;
		illum = icxIT_D50;
		obType = icxOT_CIE_1931_2;
	}
	
	/* Figure out what sort of device colorspace it is */
	{
		int ti;

		if ((ti = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
			error("Input file '%s' doesn't contain keyword COLOR_REPS",ti3name);

		if (strcmp(icg->t[0].kdata[ti],"CMYK_XYZ") == 0) {
			devspace = icSigCmykData;
			devchan = 4;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMYK_LAB") == 0) {
			devspace = icSigCmykData;
			devchan = 4;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMY_XYZ") == 0) {
			devspace = icSigCmyData;
			devchan = 3;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMY_LAB") == 0) {
			devspace = icSigCmyData;
			devchan = 3;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"RGB_XYZ") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"RGB_LAB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"iRGB_XYZ") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"iRGB_LAB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
		/* Scanner .ti3 files: */
		} else if (strcmp(icg->t[0].kdata[ti],"XYZ_RGB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"LAB_RGB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
#ifdef IMP_MONO
		} else if (strcmp(icg->t[0].kdata[ti],"K_XYZ") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"K_LAB") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"W_XYZ") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"W_LAB") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 1;
			isAdditive = 1;
#endif /* IMP_MONO */

		} else 
			error("Device input file '%s' has unhandled color representation '%s'",
			                                                     ti3name, icg->t[0].kdata[ti]);
	}

	if ((npat = icg->t[0].nsets) <= 0)
		error("Input file '%s' has no sets of data",ti3name);

	if (verb) {
		fprintf(verbo,"No of test patches = %d\n",npat);
	}

	/* Allocate arrays to hold test patch input and output values */
	if ((tpat = (pval *)malloc(sizeof(pval) * npat)) == NULL)
		error("Malloc failed - tpat[]");

	if ((stpat = (pval **)malloc(sizeof(pval *) * npat)) == NULL)
		error("Malloc failed - stpat[]");

	for (i = 0; i < npat; i++)
		stpat[i] = &tpat[i];

	/* Read in the CGATs fields */
	{
		int sidx;					/* Sample ID index */
		int sloc;					/* Sample location index (if any) */
		int ti, ci, mi, yi, ki;
		int Xi, Yi, Zi;

		if ((sidx = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
			error("Input file '%s' doesn't contain field SAMPLE_ID",ti3name);
		if (icg->t[0].ftype[sidx] != nqcs_t)
			error("Input file '%s' field SAMPLE_ID is wrong type",ti3name);

		if ((sloc = icg->find_field(icg, 0, "SAMPLE_LOC")) >= 0) {
			if (icg->t[0].ftype[sloc] != cs_t)
				error("Input file '%s' field SAMPLE_LOC is wrong type",ti3name);
		}

		if (devspace == icSigGrayData) {
			if (isAdditive) {
				if ((ci = icg->find_field(icg, 0, "GRAY_W")) < 0)
					error("Input file doesn't contain field GRAY_W");
				if (icg->t[0].ftype[ci] != r_t)
					error("Field GRAY_W is wrong type - expect float");
			} else {
				if ((ci = icg->find_field(icg, 0, "GRAY_K")) < 0)
					error("Input file doesn't contain field GRAY_K");
				if (icg->t[0].ftype[ci] != r_t)
					error("Field GRAY_K is wrong type - expect float");
			}
			mi = yi = ki = ci;

		} else if (devspace == icSigRgbData) {
			if ((ci = icg->find_field(icg, 0, "RGB_R")) < 0)
				error("Input file '%s' doesn't contain field RGB_R",ti3name);
			if (icg->t[0].ftype[ci] != r_t)
				error("Input file '%s' field RGB_R is wrong type - expect float",ti3name);
			if ((mi = icg->find_field(icg, 0, "RGB_G")) < 0)
				error("Input file '%s' doesn't contain field RGB_G",ti3name);
			if (icg->t[0].ftype[mi] != r_t)
				error("Input file '%s' field RGB_G is wrong type - expect float",ti3name);
			if ((yi = icg->find_field(icg, 0, "RGB_B")) < 0)
				error("Input file '%s' doesn't contain field RGB_B",ti3name);
			if (icg->t[0].ftype[yi] != r_t)
				error("Input file '%s' field RGB_B is wrong type - expect float",ti3name);
			ki = yi;

		} else if (devspace == icSigCmyData) {

			if ((ci = icg->find_field(icg, 0, "CMY_C")) < 0)
				error("Input file '%s' doesn't contain field CMY_C",ti3name);
			if (icg->t[0].ftype[ci] != r_t)
				error("Input file '%s' field CMY_C is wrong type - expect float",ti3name);
			if ((mi = icg->find_field(icg, 0, "CMY_M")) < 0)
				error("Input file '%s' doesn't contain field CMY_M",ti3name);
			if (icg->t[0].ftype[mi] != r_t)
				error("Input file '%s' field CMY_M is wrong type - expect float",ti3name);
			if ((yi = icg->find_field(icg, 0, "CMY_Y")) < 0)
				error("Input file '%s' doesn't contain field CMY_Y",ti3name);
			ki = yi;
		} else {	/* Assume CMYK */

			if ((ci = icg->find_field(icg, 0, "CMYK_C")) < 0)
				error("Input file '%s' doesn't contain field CMYK_C",ti3name);
			if (icg->t[0].ftype[ci] != r_t)
				error("Input file '%s' field CMYK_C is wrong type - expect float",ti3name);
			if ((mi = icg->find_field(icg, 0, "CMYK_M")) < 0)
				error("Input file '%s' doesn't contain field CMYK_M",ti3name);
			if (icg->t[0].ftype[mi] != r_t)
				error("Input file '%s' field CMYK_M is wrong type - expect float",ti3name);
			if ((yi = icg->find_field(icg, 0, "CMYK_Y")) < 0)
				error("Input file '%s' doesn't contain field CMYK_Y",ti3name);
			if (icg->t[0].ftype[yi] != r_t)
				error("Input file '%s' field CMYK_Y is wrong type - expect float",ti3name);
			if ((ki = icg->find_field(icg, 0, "CMYK_K")) < 0)
				error("Input file '%s' doesn't contain field CMYK_K",ti3name);
			if (icg->t[0].ftype[ki] != r_t)
				error("Input file '%s' field CMYK_K is wrong type - expect float",ti3name);
		}

		if (spec == 0) { 		/* Using instrument tristimulous value */

			if (isLab) {		/* Expect Lab */
				if ((Xi = icg->find_field(icg, 0, "LAB_L")) < 0)
					error("Input file '%s' doesn't contain field LAB_L",ti3name);
				if (icg->t[0].ftype[Xi] != r_t)
					error("Input file '%s' field LAB_L is wrong type - expect float",ti3name);
				if ((Yi = icg->find_field(icg, 0, "LAB_A")) < 0)
					error("Input '%s' file doesn't contain field LAB_A",ti3name);
				if (icg->t[0].ftype[Yi] != r_t)
					error("Input file '%s' field LAB_A is wrong type - expect float",ti3name);
				if ((Zi = icg->find_field(icg, 0, "LAB_B")) < 0)
					error("Input file '%s' doesn't contain field LAB_B",ti3name);
				if (icg->t[0].ftype[Zi] != r_t)
					error("Input file '%s' field LAB_B is wrong type - expect float",ti3name);

			} else { 		/* Expect XYZ */
				if ((Xi = icg->find_field(icg, 0, "XYZ_X")) < 0)
					error("Input file '%s' doesn't contain field XYZ_X",ti3name);
				if (icg->t[0].ftype[Xi] != r_t)
					error("Input file '%s' field XYZ_X is wrong type - expect float",ti3name);
				if ((Yi = icg->find_field(icg, 0, "XYZ_Y")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Y",ti3name);
				if (icg->t[0].ftype[Yi] != r_t)
					error("Input file '%s' field XYZ_Y is wrong type - expect float",ti3name);
				if ((Zi = icg->find_field(icg, 0, "XYZ_Z")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Z",ti3name);
				if (icg->t[0].ftype[Zi] != r_t)
					error("Input file '%s' field XYZ_Z is wrong type - expect float",ti3name);
			}

			for (i = 0; i < npat; i++) {
				strcpy(tpat[i].sid, (char *)icg->t[0].fdata[i][sidx]);
				if (sloc >= 0) 
					strcpy(tpat[i].slo, (char *)icg->t[0].fdata[i][sloc]);
				else
					strcpy(tpat[i].slo, "");
				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ci]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][mi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][yi]) / 100.0;
				tpat[i].p[3] = *((double *)icg->t[0].fdata[i][ki]) / 100.0;
				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0
				 || tpat[i].p[3] > 1.0) {
					error("Input file '%s' device value field value exceeds 100.0 !",ti3name);
				}
				tpat[i].v[0] = *((double *)icg->t[0].fdata[i][Xi]);
				tpat[i].v[1] = *((double *)icg->t[0].fdata[i][Yi]);
				tpat[i].v[2] = *((double *)icg->t[0].fdata[i][Zi]);
				if (!isLab && (!isdisp || isdnormed != 0)) {
					tpat[i].v[0] /= 100.0;		/* Normalise XYZ to range 0.0 - 1.0 */
					tpat[i].v[1] /= 100.0;
					tpat[i].v[2] /= 100.0;
				}
				if (!isLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
					icmXYZ2Lab(&icmD50, tpat[i].v, tpat[i].v);
				}
			}

		} else { 		/* Using spectral data */
			int ii;
			xspect sp;
			char buf[100];
			int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
			xsp2cie *sp2cie;	/* Spectral conversion object */

			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_BANDS")) < 0)
				error ("Input file '%s' doesn't contain keyword SPECTRAL_BANDS",ti3name);
			sp.spec_n = atoi(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_START_NM")) < 0)
				error ("Input file '%s' doesn't contain keyword SPECTRAL_START_NM",ti3name);
			sp.spec_wl_short = atof(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_END_NM")) < 0)
				error ("Input file '%s; doesn't contain keyword SPECTRAL_END_NM",ti3name);
			sp.spec_wl_long = atof(icg->t[0].kdata[ii]);
			if (!isdisp || isdnormed != 0)
				sp.norm = 100.0;
			else
				sp.norm = 1.0;

			/* Find the fields for spectral values */
			for (j = 0; j < sp.spec_n; j++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
				            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);

				if ((spi[j] = icg->find_field(icg, 0, buf)) < 0)
					error("Input file '%s' doesn't contain field %s",ti3name,buf);

				if (icg->t[0].ftype[spi[j]] != r_t)
					error("Field %s is wrong type - expect float",buf);
			}

			if (isemis) {
				illum = icxIT_none;		/* Emissive data has no illuminant */
				fwacomp = 0;
			}

			/* Create a spectral conversion object */
			if ((sp2cie = new_xsp2cie(illum, 0.0, illum == icxIT_none ? NULL : &cust_illum,
			                          obType, custObserver, icSigLabData, icxClamp)) == NULL)
				error("Creation of spectral conversion object failed");

			if (fwacomp) {
				double nw = 0.0;		/* Number of media white patches */
				xspect mwsp;			/* Medium spectrum */
				instType itype;			/* Spectral instrument type */
				xspect insp;			/* Instrument illuminant */

				mwsp = sp;		/* Struct copy */

	 			if ((ti = icg->find_kword(icg, 0, "TARGET_INSTRUMENT")) < 0)
					error ("Input file '%s' can't find target instrument needed for FWA compensation",ti3name);

				if ((itype = inst_enum(icg->t[0].kdata[ti])) == instUnknown)
					error ("Input file '%s' unrecognised target instrument '%s'",ti3name, icg->t[0].kdata[ti]);

				if (inst_illuminant(&insp, itype) != 0)
					error ("Instrument doesn't have an FWA illuminent");

				/* Find the media white spectral reflectance */
				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] = 0.0;

				/* Compute the mean of all the media white patches */
				for (i = 0; i < npat; i++) {
					int use = 0;
	
					if (devspace == icSigGrayData) {
						if (isAdditive) {
							if (*((double *)icg->t[0].fdata[i][ci]) > (100.0 - 0.1))
								use = 1;
						} else {
							if (*((double *)icg->t[0].fdata[i][ci]) < 0.1)
								use = 1;
						}
					} else if (devspace == icSigRgbData) {
						if (*((double *)icg->t[0].fdata[i][ci]) > (100.0 - 0.1)
						 && *((double *)icg->t[0].fdata[i][mi]) > (100.0 - 0.1)
						 && *((double *)icg->t[0].fdata[i][yi]) > (100.0 - 0.1))
							use = 1;
					} else if (devspace == icSigCmyData) {
						if (*((double *)icg->t[0].fdata[i][ci]) < 0.1
						 && *((double *)icg->t[0].fdata[i][mi]) < 0.1
						 && *((double *)icg->t[0].fdata[i][yi]) < 0.1)
							use = 1;
					} else {	/* Assume CMYK */

						if (*((double *)icg->t[0].fdata[i][ci]) < 0.1
						 && *((double *)icg->t[0].fdata[i][mi]) < 0.1
						 && *((double *)icg->t[0].fdata[i][yi]) < 0.1
						 && *((double *)icg->t[0].fdata[i][ki]) < 0.1) {
							use = 1;
						}
					}

					if (use) {
						/* Read the spectral values for this patch */
						for (j = 0; j < mwsp.spec_n; j++) {
							mwsp.spec[j] += *((double *)icg->t[0].fdata[i][spi[j]]);
						}
						nw++;
					}
				}
				if (nw == 0.0) {
					warning("Input file '%s' can't find a media white patch to init FWA",ti3name);

					/* Track the maximum reflectance for any band to determine white. */
					/* This might give bogus results if there is no white patch... */
					for (i = 0; i < npat; i++) {
						for (j = 0; j < mwsp.spec_n; j++) {
							double rv = *((double *)icg->t[0].fdata[i][spi[j]]);
							if (rv > mwsp.spec[j])
								mwsp.spec[j] = rv;
						}
					}
					nw++;
				}

				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] /= nw;	/* Compute average */

				/* If we are setting a specific simulated instrument illuminant */
				if (tillum != icxIT_none) {
					tillump = &cust_tillum;
					if (tillum != icxIT_custom) {
						if (standardIlluminant(tillump, tillum, 0.0)) {
							error("simulated inst. illum. not recognised");
						}
					}
				}

				if (sp2cie->set_fwa(sp2cie, &insp, tillump, &mwsp)) 
					error ("Set FWA on sp2cie failed");

				if (verb) {
					double FWAc;
					sp2cie->get_fwa_info(sp2cie, &FWAc);
					fprintf(verbo,"FWA content = %f\n",FWAc);
				}
				
			}

			for (i = 0; i < npat; i++) {
				strcpy(tpat[i].sid, (char *)icg->t[0].fdata[i][sidx]);
				if (sloc >= 0) 
					strcpy(tpat[i].slo, (char *)icg->t[0].fdata[i][sloc]);
				else
					strcpy(tpat[i].slo, "");
				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ci]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][mi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][yi]) / 100.0;
				tpat[i].p[3] = *((double *)icg->t[0].fdata[i][ki]) / 100.0;
				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0
				 || tpat[i].p[3] > 1.0) {
					error("Input file '%s' device value field value exceeds 100.0 !",ti3name);
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
		/* Normalize display values to Y = 1.0 if needed */
		/* (re-norm spec derived, since observer may be different) */
		if (isdisp && (isdnormed == 0 || spec != 0)) {
			double scale = -1e6;
			double bxyz[3];

			/* Locate max Y */
			for (i = 0; i < npat; i++) {
				icmLab2XYZ(&icmD50, bxyz,  tpat[i].v);
				if (bxyz[1] > scale)
					scale = bxyz[1];
			}
			
			scale = 1.0/scale;

			/* Scale max Y to 1.0 */
			for (i = 0; i < npat; i++) {
				icmLab2XYZ(&icmD50, tpat[i].v, tpat[i].v);
				tpat[i].v[0] *= scale;
				tpat[i].v[1] *= scale;
				tpat[i].v[2] *= scale;
				icmXYZ2Lab(&icmD50, tpat[i].v, tpat[i].v);
			}
		}

	}	/* End of reading in CGATs file */

	/* - - - - - - - - - - */
	/* Compute the delta E of each point against the profile value */

	/* Open up the file for reading */
	if ((rd_fp = new_icmFileStd_name(iccname,"r")) == NULL)
		error("Read: Can't open file '%s'",iccname);

	if ((rd_icco = new_icc()) == NULL)
		error("Read: Creation of ICC object failed");

	/* Read the header and tag list */
	if ((rv = rd_icco->read(rd_icco,rd_fp,0)) != 0)
		error("Read: %d, %s",rv,rd_icco->err);

	/* Get the Fwd table, absolute with Lab override */
	if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, intent,
	                              icSigLabData, icmLuOrdNorm)) == NULL) {
		error("%d, %s",rd_icco->errc, rd_icco->err);
	}

	for (i = 0; i < npat; i++) {

		/* Lookup the patch value in the profile */
		if (luo->lookup(luo, tpat[i].pv, tpat[i].p) > 1)
			error("%d, %s",rd_icco->errc,rd_icco->err);

		if (cie2k)
			tpat[i].de = icmCIE2K(tpat[i].v, tpat[i].pv);
		else if (cie94)
			tpat[i].de = icmCIE94(tpat[i].v, tpat[i].pv);
		else
			tpat[i].de = icmLabDE(tpat[i].v, tpat[i].pv);
	}

	/* - - - - - - - - - - */
	/* Sort by delta E */
	if (sortbyde) {
		/* Sort the dE's */
#define HEAP_COMPARE(A,B) (A.de > B.de)
		HEAPSORT(pval, tpat, npat);
#undef HEAP_COMPARE
	}

	/* - - - - - - - - - - */
	/* Plot a dE histogram */
	if (dohisto) {
		double demax = -1e6, demin = 1e6;
		int maxbins = 50;		/* Maximum bins */
		int minbins = 20;		/* Target minimum bins (depends on distribution) */ 
		int mincount = 10;		/* Minimum number of points in a bin */
		double mbwidth;
		int nbins = 0;
		hbin *bins = NULL;
		double tval;
		double *x, *y;
		
		/* Figure out the range of dE's */
		for (i = 0; i < npat; i++) {
			double de = tpat[i].de;

			if (de > demax)
				demax = de;
			if (de < demin)
				demin = de;
		}

		if (demax < 1e-6)
			error("histogram: dE range is too small to plot");

		/* Bin width that gives maxbins */
		mbwidth = demax / maxbins;
		
		/* Reduce mincount if needed to get minbins */
		if (npat/minbins < mincount)
			mincount = npat/minbins;

		if ((bins = (hbin *)calloc(maxbins, sizeof(hbin))) == NULL)
			error("malloc of histogram bins failed");

		/* Sort the dE's */
#define HEAP_COMPARE(A,B) (A->de < B->de)
		HEAPSORT(pval *, stpat, npat);
#undef HEAP_COMPARE

		/* Create bins and add points */
		bins[0].min = 0.0;
		for (nbins = i = 0; i < npat && nbins < maxbins; i++) {
			double de = stpat[i]->de;

			/* Move on to next bin ? */
			if (bins[nbins].count >= mincount
			 && (de - bins[nbins].min) >= mbwidth) {
				if (i > 0)
					bins[nbins].max = 0.5 * (de + stpat[i-1]->de);
				else
					bins[nbins].max = de;
				nbins++;
				bins[nbins].min = bins[nbins-1].max; 
			} 
			bins[nbins].count++;
			bins[nbins].max = de;
		}
		if (bins[nbins].count != 0)
			nbins++;

		/* Compute value */
		tval = 0.0;
		for (i = 0; i < nbins; i++) {
			bins[i].val = bins[i].count/(bins[i].max - bins[i].min);
			tval += bins[i].val;
		}

		tval /= 100.0;		/* Make it % */
		for (i = 0; i < nbins; i++) {
			bins[i].val /= tval;
			if (verb) fprintf(verbo,"Bin %d, %f - %f, % 2.4f%%, count %d\n",
			             i,bins[i].min,bins[i].max,bins[i].val,bins[i].count);
		}

		/* Plot it */
		if ((x = (double *)calloc(nbins+1, sizeof(double))) == NULL)
			error("malloc of histogram x array");
		if ((y = (double *)calloc(nbins+1, sizeof(double))) == NULL)
			error("malloc of histogram y array");

		for (i = 0; i < nbins; i++) {
			x[i] = 0.5 * (bins[i].min + bins[i].max);
			y[i] = bins[i].val;
		}
		x[i] = demax;
		y[i] = 0.0;
		do_plot(x, y, NULL, NULL, nbins+1);

		free(y);
		free(x);
		free(bins);
	}

	/* - - - - - - - - - - */
	/* Create a pruned .ti3 file */
	if (prune > 0.0) {
		char *cp, outname[MAXNAMEL+31];
		cgats *ocg;
		cgats_set_elem *setel;		/* Array of set value elements */

		strcpy(outname, ti3name);
		if ((cp = strrchr(outname, '.')) == NULL)
			cp = outname + strlen(outname);
		sprintf(cp, "_p%.2f.ti3",prune);
		
		if (verb) fprintf(verbo,"Created pruned file '%s'\n",outname);

		/* Create the output files */
		if ((ocg = new_cgats()) == NULL)
			error("Failed to create cgats object");
	
		/* Duplicate the type of the file */
		if (icg->t[0].tt == cgats_X) {
			ocg->add_other(ocg, icg->cgats_type);
			ocg->add_table(ocg, tt_other, 0);
		} else if (icg->t[0].tt == tt_other) {
			ocg->add_other(ocg, icg->others[icg->t[0].oi]);
			ocg->add_table(ocg, tt_other, 0);
		} else {
			ocg->add_table(ocg, icg->t[0].tt, 0);
		}
	
		/* Duplicate all the keywords */
		for (i = 0; i < icg->t[0].nkwords; i++) {
			ocg->add_kword(ocg, 0, icg->t[0].ksym[i], icg->t[0].kdata[i], NULL);
		}
	
		/* Duplicate all of the fields */
		for (i = 0; i < icg->t[0].nfields; i++) {
			ocg->add_field(ocg, 0, icg->t[0].fsym[i], icg->t[0].ftype[i]);
		}
	
		if ((setel = (cgats_set_elem *)malloc(
		     sizeof(cgats_set_elem) * icg->t[0].nfields)) == NULL)
			error("Malloc failed!");
	
		/* Copy them approproately */
		for (i = 0; i < icg->t[0].nsets; i++) {

			if (tpat[i].de <= prune) {
				icg->get_setarr(icg, 0, i, setel);
				ocg->add_setarr(ocg, 0, setel);
			}
		}
	
		if (verb) {
			double acc;

			acc = (double)ocg->t[0].nsets/(double)icg->t[0].nsets * 100.0;
			fprintf(verbo,"%.2f%% accepted, %.3f%% rejected\n",acc, 100.0-acc);
		}

		/* Write out the file */
		if (ocg->write_name(ocg, outname))
			error("CGATS file '%s' write error : %s",outname,ocg->err);
	}

	/* - - - - - - - - - - */
	/* Display various results */
	{
		double merr = 0.0;		/* Max */
		double aerr = 0.0;		/* Avg */
		double rerr = 0.0;		/* RMS */
		double nsamps = 0.0;
		int inn, outn;			/* Chanells for input and output spaces */

		if (dovrml) {
			wrl = new_vrml(out_name, doaxes, vrml_lab);
			wrl->start_line_set(wrl, 0);
		}

		/* Get details of conversion (Arguments may be NULL if info not needed) */
		luo->spaces(luo, NULL, &inn, NULL, &outn, NULL, NULL, NULL, NULL, NULL);

		for (i = 0; i < npat; i++) {
			double de, *out;

			de = tpat[i].de;
			out = tpat[i].pv;

			if (verb > 1) {
#ifdef HACK		// Print XYZ
#pragma message("!!!!!!!!!!!!!!!!!! profile/profcheck.c HACK enabled !!!!!!!!!!!!!!!")
				double outxyz[3], vxyz[3];
				icmLab2XYZ(&icmD50, outxyz, out);
				icmLab2XYZ(&icmD50, vxyz, tpat[i].v);

				printf("[%f] %s%s%s: %s -> %f %f %f should be %f %f %f\n",
				       de,
				       tpat[i].sid,
				       tpat[i].slo[0] != '\000' ? " @ " : "",
				       tpat[i].slo,
				       icmPdv(devchan, tpat[i].p),
				       outxyz[0],outxyz[1],outxyz[2],
				       vxyz[0],vxyz[1],vxyz[2]);
#else
				printf("[%f] %s%s%s: %s -> %f %f %f should be %f %f %f\n",
				       de,
				       tpat[i].sid,
				       tpat[i].slo[0] != '\000' ? " @ " : "",
				       tpat[i].slo,
				       icmPdv(devchan, tpat[i].p),
				       out[0],out[1],out[2],
				       tpat[i].v[0],tpat[i].v[1],tpat[i].v[2]);
#endif
			}
			if (dovrml) {
				int ix[2];
				double p1[3], p2[3];

				/* Add the vertexes */
				if (dominl && de < 0.5) {		/* Make a minimum length */
					double cent[3], vec[3], vlen;

					/* Compute center */
					icmAdd3(cent, tpat[i].v, out);
					icmScale3(cent, cent, 0.5);
					if ((vlen = de) < 1e-6) {
						vec[0] = 0.25; vec[1] = 0.0; vec[2] = 0.0;
					} else {
						icmSub3(vec, tpat[i].v, out);
						icmScale3(vec, vec, 0.25/vlen);
					}
					icmSub3(p1, cent, vec);
					icmAdd3(p2, cent, vec);
				} else {
					icmCpy3(p1, tpat[i].v);
					icmCpy3(p2, out);
				}
				ix[0] = wrl->add_vertex(wrl, 0, p1);
				ix[1] = wrl->add_vertex(wrl, 0, p2);

				/* Add the line */
				if (dodecol) {		/* Lines with color determined by length */
					double rgb[3];
					DE2RGB(rgb, icmNorm33(p1, p2));
					wrl->add_col_line(wrl, 0, ix, rgb);

				} else {	/* Natural color */
					wrl->add_line(wrl, 0, ix);
				}
			}

			/* Stats */
			aerr += de;
			rerr += de * de;

			nsamps++;
			if (de > merr)
				merr = de;

		}
		if (dovrml) {
			wrl->make_lines_vc(wrl, 0, 0.0);
			wrl->del(wrl);
		}
		printf("Profile check complete, errors%s: max. = %f, avg. = %f, RMS = %f\n",
            cie2k ? "(CIEDE2000)" : cie94 ? " (CIE94)" : "", merr, aerr/nsamps, sqrt(rerr/nsamps));

		/* ------------------------------- */
		/* If we want sort by target value */
		if (ddevv) {
			double cieval[3];

			/* Lookup the CIE value of the target */
			if (luo->lookup(luo, cieval, devval) > 1)
				error("%d, %s",rd_icco->errc,rd_icco->err);

			/* Compute deltas to target value. */
			for (i = 0; i < npat; i++) {
				if (cie2k)
					tpat[i].dv = icmCIE2K(tpat[i].v, cieval);
				else if (cie94)
					tpat[i].dv = icmCIE94(tpat[i].v, cieval);
				else
					tpat[i].dv = icmLabDE(tpat[i].v, cieval);

				tpat[i].dp = 0.0;
				for (j = 0; j < inn; j++) {
					double tt;
					tt = tpat[i].p[j] - devval[j];
					tpat[i].dp += tt * tt;
				}
				tpat[i].dp = sqrt(tpat[i].dp);
			}

			if (sortbypcs) {
				/* Sort by pcs delta */
#define HEAP_COMPARE(A,B) (A->dv < B->dv)
				HEAPSORT(pval *, stpat, npat);
#undef HEAP_COMPARE
			} else {
				/* Sort by device delta */
#define HEAP_COMPARE(A,B) (A->dp < B->dp)
				HEAPSORT(pval *, stpat, npat);
#undef HEAP_COMPARE
			}

			printf("Target point:\n");
			if (devspace == icSigCmykData) {
				printf("%f %f %f %f -> %f %f %f\n",
				       devval[0],devval[1],devval[2],devval[3],
				       cieval[0],cieval[1],cieval[2]);
			} else {	/* Assume RGB/CMY */
				printf("%f %f %f -> %f %f %f\n",
				       devval[0],devval[1],devval[2],
				       cieval[0],cieval[1],cieval[2]);
			}
			printf("\n");
	
			for (i = 0; i < npat; i++) {
				if (devspace == icSigCmykData) {
					printf("%s: %f %f %f %f [%f] -> %f %f %f [%f]\n",
					       stpat[i]->sid,
					       stpat[i]->p[0],stpat[i]->p[1],stpat[i]->p[2],stpat[i]->p[3],
					       stpat[i]->dp,
					       stpat[i]->v[0],stpat[i]->v[1],stpat[i]->v[2],
					       stpat[i]->dv);
				} else {	/* Assume RGB/CMY */
					printf("%s: %f %f %f [%f] -> %f %f %f [%f]\n",
					       stpat[i]->sid,
					       stpat[i]->p[0],stpat[i]->p[1],stpat[i]->p[2],
					       stpat[i]->dp,
					       stpat[i]->v[0],stpat[i]->v[1],stpat[i]->v[2],
					       stpat[i]->dv);
				}
			}
		}

		/* Done with lookup object */
		luo->del(luo);

		/* Close the file */
		rd_icco->del(rd_icco);
		rd_fp->del(rd_fp);
	}

	icg->del(icg);		/* Clean up */

	return 0;
}


/* ------------------------------------------------ */

/* Convert a delta E value into a signal color: */
static void
DE2RGB(double *out, double in) {
	struct {
		double de;
		double r, g, b;
	} range[6] = {
		{ 10.0, 1, 1, 0 },		/* yellow */
		{ 4.0,  1, 0, 0 },		/* red */
		{ 2.0, 1, 0, 1 },		/* magenta */
		{ 1.0, 0, 0, 1 },		/* blue */
		{ 0.5, 0, 1, 1 },		/* cyan */
		{ 0.0, 0, 1, 0 }		/* green */
	};
	int i;
	double bl;

	/* Locate the range we're in */
	if (in > range[0].de) {
		out[0] = range[0].r;
		out[1] = range[0].g;
		out[2] = range[0].b;
	} else {
		for (i = 0; i < 5; i++) {
			if (in <= range[i].de && in >= range[i+1].de)
				break;
		}
		bl = (in - range[i+1].de)/(range[i].de - range[i+1].de);
		out[0] = bl * range[i].r + (1.0 - bl) * range[i+1].r;
		out[1] = bl * range[i].g + (1.0 - bl) * range[i+1].g;
		out[2] = bl * range[i].b + (1.0 - bl) * range[i+1].b;
	}
}



