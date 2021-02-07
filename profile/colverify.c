/* 
 * Argyll Color Correction System
 * Verify two sets of PCS values.
 *
 * Author: Graeme W. Gill
 * Date:   7/6/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in two CGATS files (probably but not necesserily .ti3 files) of PCS
 * values (either XYZ, L*a*b* or spectral), matches the values, and computes
 * overall errors. This is useful for verifying proofing systems.
 */

/*
 * TTBD:
 *
 *	We're no warning about setting illuminant or FWA for emissive spectral data,
 *  they are just ignored.
 */

#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include <string.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "vrml.h"
#include "cgats.h"
#include "xicc.h"
#include "insttypes.h"
#include "disptechs.h"
#include "ccmx.h"
#include "sort.h"
#include "plot.h"
#include "ui.h"

#ifdef DEBUG
#undef DBG
#define DBG(xxx) printf xxx ;
#else
#undef DBG
#define DBG(xxx) 
#endif

void
usage(void) {
	fprintf(stderr,"Verify CIE values, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: colverify [-options] target.ti3 measured.ti3\n");
	fprintf(stderr," -v [n]          Verbose mode, n >= 2 print each value\n");
	fprintf(stderr," -l              Match patches by sample location rather than id\n");
	fprintf(stderr," -n              Normalise each files reading to its white Y\n");
	fprintf(stderr," -N              Normalise each files reading to its white XYZ\n");
	fprintf(stderr," -m              Normalise each files reading to its white X+Y+Z\n");
	fprintf(stderr," -M              Normalise both files reading to mean white XYZ\n");
	fprintf(stderr," -D              Use D50 100.0 as L*a*b* white reference\n");
	fprintf(stderr," -c              Show CIE94 delta E values\n");
	fprintf(stderr," -k              Show CIEDE2000 delta E values\n");
	fprintf(stderr," -h [hist.txt]   Plot a histogram of delta E's [Optionaly save points to .txt]\n");
	fprintf(stderr," -s              Sort patch values by error\n");
	fprintf(stderr," -w              create PCS %s vector visualisation (measured%s)\n",vrml_format(),vrml_ext());
	fprintf(stderr," -W              create PCS %s marker visualisation (measured%s)\n",vrml_format(),vrml_ext());
	fprintf(stderr," -d              create Device RGB %s marker visualisation (measured%s)\n",vrml_format(),vrml_ext());
//	fprintf(stderr," -d y            create Device YCbCr %s marker visualisation (measured%s)\n",vrml_format(),vrml_ext());
	fprintf(stderr," -x              Use %s axes\n",vrml_format());
	fprintf(stderr," -f [illum]      Use Fluorescent Whitening Agent compensation [opt. simulated inst. illum.:\n");
	fprintf(stderr,"                  M0, M1, M2, A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp]\n");
	fprintf(stderr," -i illum        Choose illuminant for computation of CIE XYZ from spectral data & FWA:\n");
	fprintf(stderr,"                  A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                  1931_2 (def), 1964_10, 2012_2, 2012_10, S&B 1955_2, shaw, J&V 1978_2 or file.cmf\n");
	fprintf(stderr," -L profile.%s  Skip any first file, out of profile gamut patches\n",ICC_FILE_EXT_ND);
	fprintf(stderr," -X file.ccmx    Apply Colorimeter Correction Matrix to second file\n");
//	fprintf(stderr," -Z A|X          Just print Average|Max +tab\n");
	fprintf(stderr," target.ti3      Target (reference) PCS or spectral values.\n");
	fprintf(stderr," measured.ti3    Measured (actual) PCS or spectral values\n");
	exit(1);
	}

/* Patch value type */
typedef struct {
	char sid[50];		/* sample id */
	char loc[100];		/* sample location (empty if none) */
	double rgb[3];		/* RGB value if RGB device space present, or YCbCr if dovrml==4 */
	double ycc[3];		/* YCbCr if RGB and dovrml==4 */
	int og;				/* Out of gamut flag */
	double xyz[3];		/* XYZ value */
	double v[3];		/* Lab value */
	double de;			/* Delta E */
	double ixde[3];		/* XYZ Component DE */
	double ide[3];		/* Lab Component DE */
} pval;

/* Histogram bin type */
typedef struct {
	int count;			/* Raw count */
	double val;			/* Normalized value */
	double min, max;	/* Bin range */
} hbin;

int main(int argc, char *argv[])
{
	int fa,nfa,mfa;			/* current argument we're looking at */
	int verb = 0;       	/* Verbose level */
	int useloc = 0;			/* Match patches by sample location */
	int norm = 0;			/* 1 = norm to White Y, 2 = norm to White XYZ */
							/* 3 = norm to White X+Y+Z, 4 = norm to average XYZ */
	int usestdd50 = 0;		/* Use standard D50 instead of avg white as reference */
	int cie94 = 0;
	int cie2k = 0;
	int dovrml = 0;			/* 1 = PCS vector, 2 = PCS marker, 3 = RGB, 4 - YCbCr */
	int doaxes = 0;
	int dohisto = 0;			/* Plot histogram of delta E's */
	char histoname[MAXNAMEL+1] = "\000";  /* Optional file to save histogram points to */
	int dosort = 0;
	int dozrep = 0;				/* 1 = print average, 2 = print max */
	char ccmxname[MAXNAMEL+1] = "\000";  /* Colorimeter Correction Matrix name */
	ccmx *cmx = NULL;					/* Colorimeter Correction Matrix */
	char gprofname[MAXNAMEL+1] = "\000";  /* Gamut limit profile name */
	icmFile *fp = NULL;
	icc *icco = NULL;
	xicc *xicco = NULL;
	icxLuBase *luo = NULL;

	struct {
		char name[MAXNAMEL+1];	/* Patch filename  */
		int isemis;				/* nz if emsissive spectral reference data */
		int isdisp;				/* nz if display */
		int isdnormed;      	/* Has display data been normalised to 100 ? */
		int isrgb;				/* Is RGB device space ? */
		int npat;				/* Number of patches */
		int nig;				/* Number of patches in gamut */
		double w[3];			/* XYZ of "white" */
		double nw[3];			/* Normalised XYZ of "white" */
		pval *pat;				/* patch values */
	} cg[2];					/* Target and current patch file information */

	int *match;					/* Array mapping first list indexes to corresponding second */
	int *sort;					/* Array of first list indexes in sorted order */
	int fwacomp = 0;			/* FWA compensation */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType tillum = icxIT_none;	/* Target/simulated instrument illuminant */ 
	xspect cust_tillum, *tillump = NULL; /* Custom target/simulated illumination spectrum */
	icxIllumeType illum = icxIT_none;	/* Spectral defaults to D50 */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType obType = icxOT_none;	/* Defaults to 1931 2 degree */
	xspect custObserver[3];		/* If obType = icxOT_custom */

	icmXYZNumber labw = icmD50;	/* The Lab white reference */

	char out_name[MAXNAMEL+4+1]; /* VRML/X3D name */
	vrml *wrl = NULL;

	int i, j, n;

	if (argc <= 1)
		usage();

	/* Process the arguments */
	mfa = 2;        				/* Minimum final arguments */
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

			if (argv[fa][1] == '?')
				usage();

			/* Verbose */
			else if (argv[fa][1] == 'v') {
				verb = 1;

				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					fa = nfa;
					verb = atoi(na);
				}
			}

			/* Use location to match patches */
			else if (argv[fa][1] == 'l') {
				useloc = 1;
			}

			/* normalize */
			else if (argv[fa][1] == 'n'
			      || argv[fa][1] == 'N') {
				norm = 1;
				if (argv[fa][1] == 'N')
					norm = 2;
			}

			else if (argv[fa][1] == 'm') {
				norm = 3;
			}

			else if (argv[fa][1] == 'M') {
				norm = 4;
			}

			else if (argv[fa][1] == 'D')
				usestdd50 = 1;

			/* VRML/X3D */
			else if (argv[fa][1] == 'w')
				dovrml = 1;

			else if (argv[fa][1] == 'W')
				dovrml = 2;

			else if (argv[fa][1] == 'd') {
				dovrml = 3;
				if (na != NULL) {	/* Argument is present - RGB or YCbCr. */
					fa = nfa;
					if (strcmp(na, "y") == 0)
						dovrml = 4;
					else
						usage();
				}					
			}

			/* Axes */
			else if (argv[fa][1] == 'x')
				doaxes = 1;

			/* CIE94 delta E */
			else if (argv[fa][1] == 'c') {
				cie94 = 1;
				cie2k = 0;
			}

			else if (argv[fa][1] == 'k') {
				cie94 = 0;
				cie2k = 1;
			}

			/* Plot histogram */
			else if (argv[fa][1] == 'h') {
				dohisto = 1;

				if (na != NULL) {	/* Argument is present - file to save points to */
					fa = nfa;
					strncpy(histoname,na,MAXNAMEL); histoname[MAXNAMEL] = '\000';
				}
			}

			/* Sort */
			else if (argv[fa][1] == 's')
				dosort = 1;

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

			/* Spectral to CIE Illuminant type */
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

			/* Gamut limit profile for first file */
			else if (argv[fa][1] == 'L') {
				if (na == NULL) usage();
				fa = nfa;
				strncpy(gprofname,na,MAXNAMEL-1); gprofname[MAXNAMEL-1] = '\000';
			}

			/* Colorimeter Correction Matrix for second file */
			else if (argv[fa][1] == 'X') {
				if (na == NULL) usage();
				fa = nfa;
				strncpy(ccmxname,na,MAXNAMEL-1); ccmxname[MAXNAMEL-1] = '\000';
			}

			else if (argv[fa][1] == 'Z') {
				if (na == NULL) usage();
				fa = nfa;
				if (strcmp(na, "A") == 0) {
					dozrep = 1;
				} else if (strcmp(na, "X") == 0) {
					dozrep = 2;
				} else
					usage();

			} else 
				usage();
		} else
			break;
	}

	/* Get the file name arguments */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(cg[0].name,argv[fa++],MAXNAMEL); cg[0].name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(cg[1].name,argv[fa],MAXNAMEL); cg[1].name[MAXNAMEL] = '\000';

	/* Create VRML/X3D base name */
	{
		char *xl;
		strcpy(out_name, cg[1].name);
		if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
			xl = out_name + strlen(out_name);
		xl[0] = '\000';		/* Remove extension */
	}

	if (fwacomp && spec == 0)
		error ("FWA compensation only works when viewer and/or illuminant selected");

	/* Gamut limit profile */
	if (gprofname[0] != '\000') {
		int rv;

		if ((fp = new_icmFileStd_name(gprofname,"r")) == NULL)
			error ("Can't open file '%s'",gprofname);
	
		if ((icco = new_icc()) == NULL)
			error ("Creation of ICC object failed");
	
		if ((rv = icco->read(icco,fp,0)) != 0)
			error("Reading profile '%s' failed failed with error %d:'%s'\n",
		     	       gprofname, icco->errc,  icco->err);

		if (icco->header->deviceClass != icSigInputClass
		 && icco->header->deviceClass != icSigDisplayClass
		 && icco->header->deviceClass != icSigOutputClass)
			error("Profile '%s' must be a device profile to filter by gamut",gprofname);

		/* Wrap with an expanded icc */
		if ((xicco = new_xicc(icco)) == NULL)
			error ("Creation of xicc failed");

		/* Get a expanded color conversion object */
		if ((luo = xicco->get_luobj(xicco, ICX_CLIP_NEAREST | ICX_FAST_SETUP,
		    icmFwd, icRelativeColorimetric, icSigXYZData, icmLuOrdNorm, NULL, NULL)) == NULL)
			error ("%d, %s",xicco->errc, xicco->err);
	}

	/* Colorimeter Correction Matrix */
	if (ccmxname[0] != '\000') {
		if ((cmx = new_ccmx()) == NULL)
			error("new_ccmx failed\n");
		if (cmx->read_ccmx(cmx,ccmxname))
			error("Reading Colorimeter Correction Matrix file '%s' failed with error %d:'%s'\n",
		     	       ccmxname, cmx->errc,  cmx->err);
	}


	/* Open up each file in turn, target then measured, */
	/* and read in the CIE values. */
	for (n = 0; n < 2; n++) {
		cgats *cgf = NULL;			/* cgats file data */
		int isLab = 0;				/* 0 if file CIE is XYZ, 1 if is Lab */
		int sidx;					/* Sample ID index */
		int sldx = -1;				/* Sample location index, < 0 if invalid */
		int xix, yix, zix;
		int rgbix[3];				/* RGB field indexes (if rgb ) */
		int dti;					/* Device type index */

		/* Open CIE target values */
		cgf = new_cgats();			/* Create a CGATS structure */
		cgf->add_other(cgf, ""); 	/* Allow any signature file */
	
		DBG(("Opening file '%s'\n",cg[n].name))

		if (cgf->read_name(cgf, cg[n].name))
			error("CGATS file '%s' read error : %s",cg[n].name,cgf->err);
	
		if (cgf->ntables < 1)
			error ("Input file '%s' doesn't contain at least one table",cg[n].name);
	
		if ((dti = cgf->find_kword(cgf, 0, "DEVICE_CLASS")) < 0)
			warning("Input file '%s' doesn't contain keyword DEVICE_CLASS",cg[n].name);

		/* Figure out what sort of device it is */
		{
			int ti;
	
			cg[n].isemis = 0;
			cg[n].isdisp = 0;
			cg[n].isdnormed = 0;
			cg[n].w[0] = cg[n].w[1] = cg[n].w[2] = 0.0;

			if (dti >= 0) {
				if (strcmp(cgf->t[0].kdata[dti],"DISPLAY") == 0) {
					cg[n].isemis = 1;
					cg[n].isdisp = 1;
					cg[n].isdnormed = 1;	/* Assume display type is normalised to 100 */

				} else if (strcmp(cgf->t[0].kdata[dti],"EMISINPUT") == 0) {
					cg[n].isemis = 1;
				}
	
				if (cg[n].isdisp) {
	
					if ((ti = cgf->find_kword(cgf, 0, "LUMINANCE_XYZ_CDM2")) >= 0) {
						if (sscanf(cgf->t[0].kdata[ti], " %lf %lf %lf ",&cg[n].w[0], &cg[n].w[1], &cg[n].w[2]) != 3)
							cg[n].w[0] = cg[n].w[1] = cg[n].w[2] = 0.0;
					}
	
					/* See if there is an explicit tag indicating data has been normalised to Y = 100 */
					if ((ti = cgf->find_kword(cgf, 0, "NORMALIZED_TO_Y_100")) >= 0) {
						if (strcmp(cgf->t[0].kdata[ti],"NO") == 0) {
							cg[n].isdnormed = 0;
						} else {
							cg[n].isdnormed = 1;
						}
					}
				}
			}
		}

		/* Check if the file is suitable */
		if (!spec
		 && cgf->find_field(cgf, 0, "LAB_L") < 0
		 && cgf->find_field(cgf, 0, "XYZ_X") < 0) {
	
			if (cgf->find_kword(cgf, 0, "SPECTRAL_BANDS") < 0)
				error ("Neither CIE nor spectral data found in file '%s'",cg[n].name);
	
			/* Switch to using spectral information */
			if (verb)
				printf("No CIE data found, switching to spectral with standard observer & D50 for file '%s'\n",cg[n].name);
			spec = 1;
		}
		if (spec && cgf->find_kword(cgf, 0, "SPECTRAL_BANDS") < 0)
			error ("No spectral data data found in file '%s' when spectral expected",cg[n].name);
	
		if (!spec && cgf->find_field(cgf, 0, "LAB_L") >= 0)
			isLab = 1;
		
		cg[n].nig = cg[n].npat = cgf->t[0].nsets;		/* Number of patches */

		/* See if it has RGB device space (for -d option) */
		if ((rgbix[0] = cgf->find_field(cgf, 0, "RGB_R")) >= 0
		 && cgf->t[0].ftype[rgbix[0]] == r_t

		 &&	(rgbix[1] = cgf->find_field(cgf, 0, "RGB_G")) >= 0
		 && cgf->t[0].ftype[rgbix[1]] == r_t

		 &&	(rgbix[2] = cgf->find_field(cgf, 0, "RGB_B")) >= 0
		 && cgf->t[0].ftype[rgbix[2]] == r_t) {
			cg[n].isrgb = 1;
		} else {
			cg[n].isrgb = 0;
		}

		/* Read all the target patches */
		if (cg[n].npat <= 0)
			error("No sets of data in file '%s'",cg[n].name);
	
		if (verb && n == 0) {
			fprintf(verbo,"No of test patches = %d\n",cg[n].npat);
		}
	
		/* Allocate arrays to hold test patch input and output values */
		if ((cg[n].pat = (pval *)malloc(sizeof(pval) * cg[n].npat)) == NULL)
			error("Malloc failed - pat[]");
	
		/* Read in the CGATs fields */
		if ((sidx = cgf->find_field(cgf, 0, "SAMPLE_ID")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SampleName")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "Sample_Name")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SAMPLE_NAME")) < 0)
			error("Input file '%s' doesn't contain field SAMPLE_ID, SampleName, Sample_Name, SAMPLE_NAME",cg[n].name);
		if (cgf->t[0].ftype[sidx] != nqcs_t
		 && cgf->t[0].ftype[sidx] != cs_t)
			error("Sample ID/Name field isn't a quoted or non quoted character string");

		if ((sldx = cgf->find_field(cgf, 0, "SAMPLE_LOC")) < 0
		 || cgf->t[0].ftype[sldx] != cs_t)
			sldx = -1;

		if (spec == 0) { 		/* Using instrument tristimulous value */

			if (isLab) {		/* Expect Lab */
				if ((xix = cgf->find_field(cgf, 0, "LAB_L")) < 0)
					error("Input file '%s' doesn't contain field LAB_L",cg[n].name);
				if (cgf->t[0].ftype[xix] != r_t)
					error("Field LAB_L is wrong type - expect float");
				if ((yix = cgf->find_field(cgf, 0, "LAB_A")) < 0)
					error("Input file '%s' doesn't contain field LAB_A",cg[n].name);
				if (cgf->t[0].ftype[yix] != r_t)
					error("Field LAB_A is wrong type - expect float");
				if ((zix = cgf->find_field(cgf, 0, "LAB_B")) < 0)
					error("Input file '%s' doesn't contain field LAB_B",cg[n].name);
				if (cgf->t[0].ftype[zix] != r_t)
					error("Field LAB_B is wrong type - expect float");

			} else { 		/* Expect XYZ */
				if ((xix = cgf->find_field(cgf, 0, "XYZ_X")) < 0)
					error("Input file '%s' doesn't contain field XYZ_X",cg[n].name);
				if (cgf->t[0].ftype[xix] != r_t)
					error("Field XYZ_X is wrong type - expect float");
				if ((yix = cgf->find_field(cgf, 0, "XYZ_Y")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Y",cg[n].name);
				if (cgf->t[0].ftype[yix] != r_t)
					error("Field XYZ_Y is wrong type - expect float");
				if ((zix = cgf->find_field(cgf, 0, "XYZ_Z")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Z",cg[n].name);
				if (cgf->t[0].ftype[zix] != r_t)
					error("Field XYZ_Z is wrong type - expect float");
			}

			for (i = 0; i < cg[n].npat; i++) {
				strcpy(cg[n].pat[i].sid, (char *)cgf->t[0].fdata[i][sidx]);
				if (sldx >= 0)
					strcpy(cg[n].pat[i].loc, (char *)cgf->t[0].fdata[i][sldx]);
				else
					cg[n].pat[i].loc[0] = '\000';
				cg[n].pat[i].og = 0;
				cg[n].pat[i].xyz[0] = *((double *)cgf->t[0].fdata[i][xix]);
				cg[n].pat[i].xyz[1] = *((double *)cgf->t[0].fdata[i][yix]);
				cg[n].pat[i].xyz[2] = *((double *)cgf->t[0].fdata[i][zix]);

				if (isLab) {	/* Convert Lab to XYZ 0..100% */
					icmLab2XYZ(&icmD50_100, cg[n].pat[i].xyz, cg[n].pat[i].xyz);
				}
//printf("~1 file %d patch %d = XYZ %f %f %f\n", n,i,cg[n].pat[i].xyz[0],cg[n].pat[i].xyz[1],cg[n].pat[i].xyz[2]);

				/* restore normalised display values to absolute */
				if (cg[n].isdnormed) {
					if (cg[n].w[1] > 0.0) {		// Found absoluute display white tag
						cg[n].pat[i].xyz[0] *= cg[n].w[1]/100.0;
						cg[n].pat[i].xyz[1] *= cg[n].w[1]/100.0;
						cg[n].pat[i].xyz[2] *= cg[n].w[1]/100.0;
					}

				} else if (!cg[n].isdisp) {
					/* If reflective or transmissive that are 0..100%, */
					/* scale back to 0.. 1 */
					cg[n].pat[i].xyz[0] /= 100.0;		/* scale back to XYZ 1.0 */
					cg[n].pat[i].xyz[1] /= 100.0;
					cg[n].pat[i].xyz[2] /= 100.0;
				}

				/* Apply ccmx */
				if (n == 1 && cmx != NULL) {
					cmx->xform(cmx, cg[n].pat[i].xyz, cg[n].pat[i].xyz);
				}

				if ((dovrml == 3 || dovrml == 4) && cg[n].isrgb) {
					cg[n].pat[i].rgb[0] = 0.01 * *((double *)cgf->t[0].fdata[i][rgbix[0]]);
					cg[n].pat[i].rgb[1] = 0.01 * *((double *)cgf->t[0].fdata[i][rgbix[1]]);
					cg[n].pat[i].rgb[2] = 0.01 * *((double *)cgf->t[0].fdata[i][rgbix[2]]);

					if (dovrml == 4) {
						icmRec709_RGBd_2_YPbPr(cg[n].pat[i].ycc, cg[n].pat[i].rgb);
					}
				}
			}

		} else { 		/* Using spectral data */
			int ii;
			xspect sp;
			char buf[100];
			int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
			xsp2cie *sp2cie;	/* Spectral conversion object */

			/* Copies of global values: */
			int l_fwacomp = fwacomp;
			int l_spec = spec;
			icxIllumeType l_tillum = tillum;
			xspect *l_tillump = tillump;
			icxIllumeType l_illum = illum;
			icxObserverType l_observ = obType;

			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_BANDS")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_BANDS");
			sp.spec_n = atoi(cgf->t[0].kdata[ii]);
			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_START_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_START_NM");
			sp.spec_wl_short = atof(cgf->t[0].kdata[ii]);
			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_END_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_END_NM");
			sp.spec_wl_long = atof(cgf->t[0].kdata[ii]);
			if (!cg[n].isdisp || cg[n].isdnormed != 0)
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

				if ((spi[j] = cgf->find_field(cgf, 0, buf)) < 0)
					error("Input file doesn't contain field %s",buf);

				if (cgf->t[0].ftype[spi[j]] != r_t)
					error("Field %s is wrong type - expect float",buf);
			}

			/* Create a spectral conversion object */
			if (cg[n].isemis) {
				if (l_illum != icxIT_none)
					warning("-i illuminant ignored for emissive reference type");
				if (l_fwacomp)
					warning("-f FWA ignored for emissive reference type");
				l_illum = icxIT_none;		/* Make emissive conversion */
				l_tillum = icxIT_none;
				l_fwacomp = 0;
			} else {
				/* Set default */
				if (l_illum == icxIT_none)
					l_illum = icxIT_D50;
			}
	
			/* Set default */
			if (l_observ == icxOT_none)
				l_observ = icxOT_CIE_1931_2;

			if ((sp2cie = new_xsp2cie(l_illum, 0.0, l_illum == icxIT_none ? NULL : &cust_illum,
			                          l_observ, custObserver, icSigXYZData, icxClamp)) == NULL)
				error("Creation of spectral conversion object failed");

			if (l_fwacomp) {
				int ti;
				xspect mwsp;			/* Medium spectrum */
				instType itype;			/* Spectral instrument type */
				xspect insp;			/* Instrument illuminant */

				mwsp = sp;		/* Struct copy */

	 			if ((ti = cgf->find_kword(cgf, 0, "TARGET_INSTRUMENT")) < 0)
					error ("Can't find target instrument in '%s' needed for FWA compensation",cg[n].name);

				if ((itype = inst_enum(cgf->t[0].kdata[ti])) == instUnknown)
					error ("Unrecognised target instrument '%s'", cgf->t[0].kdata[ti]);

				if (inst_illuminant(&insp, itype) != 0)
					error ("Instrument doesn't have an FWA illuminent");

				/* Determine a media white spectral reflectance */
				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] = 0.0;

				/* Since we don't want to assume that there are any associated device */
				/* values present in each file, we can't use this as means of */
				/* determining the media color. Use an alternative approach here, */
				/* which may give slightly different results to profile. */

				/* Track the maximum reflectance for any band to determine white. */
				/* This might silently fail, if there isn't white in the sample set. */
				for (i = 0; i < cg[0].npat; i++) {
					for (j = 0; j < mwsp.spec_n; j++) {
						double rv = *((double *)cgf->t[0].fdata[i][spi[j]]);
						if (rv > mwsp.spec[j])
							mwsp.spec[j] = rv;
					}
				}

				/* If we are setting a specific simulated instrument illuminant */
				if (l_tillum != icxIT_none) {
					l_tillump = &cust_tillum;
					if (l_tillum != icxIT_custom) {
						if (standardIlluminant(l_tillump, l_tillum, 0.0)) {
							error("simulated inst. illum. not recognised");
						}
					}
				}

				if (sp2cie->set_fwa(sp2cie, &insp, l_tillump, &mwsp)) 
					error ("Set FWA on sp2cie failed");

				if (verb) {
					double FWAc;
					sp2cie->get_fwa_info(sp2cie, &FWAc);
					fprintf(verbo,"FWA content = %f\n",FWAc);
				}
			}

			for (i = 0; i < cg[0].npat; i++) {

				strcpy(cg[n].pat[i].sid, (char *)cgf->t[0].fdata[i][sidx]);
				if (sldx >= 0)
					strcpy(cg[n].pat[i].loc, (char *)cgf->t[0].fdata[i][sldx]);
				else
					cg[n].pat[i].loc[0] = '\000';
				cg[n].pat[i].og = 0;

				/* Read the spectral values for this patch */
				for (j = 0; j < sp.spec_n; j++) {
					sp.spec[j] = *((double *)cgf->t[0].fdata[i][spi[j]]);
				}

				/* Convert it to XYZ space */
				sp2cie->convert(sp2cie, cg[n].pat[i].xyz, &sp);

				/* restore normalised display values to absolute */
				if (cg[n].isdnormed) {
					if (cg[n].w[1] > 0.0) {		// Found absoluute display white tag
						cg[n].pat[i].xyz[0] *= cg[n].w[1];
						cg[n].pat[i].xyz[1] *= cg[n].w[1];
						cg[n].pat[i].xyz[2] *= cg[n].w[1];
					}

				}

				/* Apply ccmx */
				if (n == 1 && cmx != NULL) {
					cmx->xform(cmx, cg[n].pat[i].xyz, cg[n].pat[i].xyz);
				}
			}

			sp2cie->del(sp2cie);		/* Done with this */
		}	/* End of reading in CGATs file spectral data */


		/* Locate the patch with maximum Y, a possible white patch */
		/* in case we need it latter. */
		{
			int ii;

			if (cg[n].w[1] == 0.0) {	/* No display white patch tag */

				/* Locate patch with biggest Y, assume it is white... */
				for (i = 0; i < cg[n].npat; i++) {
					if (cg[n].pat[i].xyz[1] > cg[n].w[1]) {
						icmCpy3(cg[n].w, cg[n].pat[i].xyz);
						ii = i;
					}
				}
				if (verb) printf("File %d Chose patch %d as white, XYZ %f %f %f\n",
				                       n, ii+1,cg[n].w[0],cg[n].w[1],cg[n].w[2]);
			} else {
				if (verb) printf("File %d White is from display luminance ref. XYZ %f %f %f\n",
				                       n, cg[n].w[0],cg[n].w[1],cg[n].w[2]);
			}
			icmCpy3(cg[n].nw, cg[n].w);
		}
		cgf->del(cgf);		/* Clean up */
	}	/* Next file */

	if (norm == 4) {		/* Normalise to average of white XYZ of the two files */
		icmBlend3(cg[0].w, cg[0].w, cg[1].w, 0.5);
		icmCpy3(cg[1].w, cg[0].w);
//		if (verb) printf("Average White XYZ %f %f %f\n",cg[0].w[0],cg[0].w[1],cg[0].w[2]);
	}

	/* For both files */
	for (n = 0; n < 2; n++) {

		/* Normalise this file to white = 1.0 or D50 */
		if (norm) {
			int ii;

			double chmat[3][3];				/* Chromatic adapation matrix */

			DBG(("Normalizng '%s' to white\n",cg[n].name))

			if (norm == 2 || norm == 4) {		/* Norm to white XYZ */ 
				icmXYZNumber s_wp;
				icmAry2XYZ(s_wp, cg[n].w);
				icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, s_wp, chmat);
			}

			for (i = 0; i < cg[n].npat; i++) {
				if (norm == 1) {
					cg[n].pat[i].xyz[0] *= 100.0 / cg[n].w[1];
					cg[n].pat[i].xyz[1] *= 100.0 / cg[n].w[1];
					cg[n].pat[i].xyz[2] *= 100.0 / cg[n].w[1];
				} else if (norm == 2 || norm == 4) { 
					icmMulBy3x3(cg[n].pat[i].xyz, chmat, cg[n].pat[i].xyz);
				} else {
					cg[n].pat[i].xyz[0] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
					cg[n].pat[i].xyz[1] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
					cg[n].pat[i].xyz[2] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
				}
//printf("~1 file %d patch %d = norm XYZ %f %f %f\n", n,i,cg[n].pat[i].xyz[0],cg[n].pat[i].xyz[1],cg[n].pat[i].xyz[2]);
			}
			/* Compute normalised white too */
			if (norm == 1) {
				cg[n].nw[0] *= 100.0 / cg[n].w[1];
				cg[n].nw[1] *= 100.0 / cg[n].w[1];
				cg[n].nw[2] *= 100.0 / cg[n].w[1];
			} else if (norm == 2 || norm == 4) { 
				icmMulBy3x3(cg[n].nw, chmat, cg[n].w);
			} else {
				cg[n].nw[0] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
				cg[n].nw[1] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
				cg[n].nw[2] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
			}
//printf("~1 file %d norm white XYZ %f %f %f\n", n,cg[n].nw[0], cg[n].nw[1], cg[n].nw[2]);
		}
	}	/* Next file */

	if (cmx != NULL) 
		cmx->del(cmx);
	cmx = NULL;

	/* Check that the number of test patches matches */
	if (cg[0].npat != cg[1].npat)
		error("Number of patches between '%s' and '%s' doesn't match",cg[0].name,cg[1].name);
	
	/* Create a list to map the second list of patches to the first */
	if ((match = (int *)malloc(sizeof(int) * cg[0].npat)) == NULL)
		error("Malloc failed - match[]");

	/* Use location to match */
	if (useloc) {
		for (i = 0; i < cg[0].npat; i++) {
			if (cg[0].pat[0].loc[0] == '\000'
			 || cg[1].pat[0].loc[0] == '\000')
				error("Need both files to have SAMPLE_LOC field to match on location");

			for (j = 0; j < cg[1].npat; j++) {
				if (strcmp(cg[0].pat[i].loc, cg[1].pat[j].loc) == 0)
					break;			/* Found it */
			}
			if (j < cg[1].npat) {
				match[i] = j;
			} else {
				error("Failed to find matching patch to '%s'",cg[0].pat[i].loc);
			}
		}

	/* Use id */
	} else {
		for (i = 0; i < cg[0].npat; i++) {
			for (j = 0; j < cg[1].npat; j++) {
				if (strcmp(cg[0].pat[i].sid, cg[1].pat[j].sid) == 0)
					break;			/* Found it */
			}
			if (j < cg[1].npat) {
				match[i] = j;
			} else {
				error("Failed to find matching patch to '%s'",cg[0].pat[i].sid);
			}
		}
	}

	/* Figure out which patches to skip because they are out of gamut */
	if (luo != NULL) {
		double chmat[3][3];				/* Chromatic adapation matrix */
		double out[MAX_CHAN], in[3], check[3];
		icmXYZNumber s_wp;
		int rv;

		DBG(("Figuring out of gamut patches\n"))

		/* Convert sample PCS to relative */
//printf("   cg[0].w %f %f %f\n", cg[0].w[0], cg[0].w[1], cg[0].w[2]);
		icmAry2XYZ(s_wp, cg[0].w);
		s_wp.X /= s_wp.Y;		// Normalise the white to 1.0
		s_wp.Y /= s_wp.Y;		// so that matrix doesn't change magnitude
		s_wp.Z /= s_wp.Y;
//printf("   s_wp %f %f %f\n", s_wp.X, s_wp.Y, s_wp.Z);
		icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, s_wp, chmat);
//printf("~1 matrix = \n");
//printf("   %f %f %f\n", chmat[0][0], chmat[0][1], chmat[0][2]);
//printf("   %f %f %f\n", chmat[1][0], chmat[1][1], chmat[1][2]);
//printf("   %f %f %f\n", chmat[2][0], chmat[2][1], chmat[2][2]);

		for (i = 0; i < cg[0].npat; i++) {
			icmMulBy3x3(in, chmat, cg[0].pat[i].xyz);

//printf("~1 %d: xyz %f %f %f, rel %f %f %f\n", i+1, cg[0].pat[i].xyz[0], cg[0].pat[i].xyz[1], cg[0].pat[i].xyz[2], in[0], in[1], in[2]);

			if ((rv = luo->inv_lookup(luo, out, in)) > 0 || 1) {
				double de;

				luo->lookup(luo, check, out);
				de = icmXYZLabDE(&icmD50,check, in);
//printf("~1 %d: rv %d, de %f, check XYZ %f %f %f\n",i+1,rv, de, check[0],check[1],check[2]);

				if (de >= 0.01) {
					cg[0].pat[i].og = 1;
//printf("~1 Patch %d is out of gamut by DE %f RGB %f %f %f\n",i+1,de,out[0],out[1],out[2]);
					if (verb >= 3)
						printf("Patch %d is out of gamut by DE %f\n",i+1,de);
					cg[0].nig--;
				}
			}
		}
		if (verb)
			fprintf(verbo,"No of test patches in gamut = %d/%d\n",cg[0].nig,cg[0].npat);
	}

	if (cg[0].nig <= 0) {
		if (verb)
			fprintf(verbo,"No test patches in gamut - givig up\n");
		return 0;
	}

	/* Adjust the Lab reference white to be the mean of the white of the two files */
	if (norm != 0 && !usestdd50) {
		labw.X = 0.5 * (cg[0].nw[0] + cg[1].nw[0]);
		labw.Y = 0.5 * (cg[0].nw[1] + cg[1].nw[1]);
		labw.Z = 0.5 * (cg[0].nw[2] + cg[1].nw[2]);

		if (verb)
			printf("L*a*b* white reference = XYZ %f %f %f\n",labw.X,labw.Y,labw.Z);
	}

	/* labw defaults to D50 */

	/* Convert XYZ to Lab */
	for (n = 0; n < 2; n++) {
		for (i = 0; i < cg[n].npat; i++) {
			icmXYZ2Lab(&labw, cg[n].pat[i].v, cg[n].pat[i].xyz);
		}
	}

	/* Compute the delta E's */
	DBG(("Computing the delta E's\n"))
	for (i = 0; i < cg[0].npat; i++) {

		cg[0].pat[i].ixde[0] = fabs(cg[0].pat[i].xyz[0] - cg[1].pat[match[i]].xyz[0]);
		cg[0].pat[i].ixde[1] = fabs(cg[0].pat[i].xyz[1] - cg[1].pat[match[i]].xyz[1]);
		cg[0].pat[i].ixde[2] = fabs(cg[0].pat[i].xyz[2] - cg[1].pat[match[i]].xyz[2]);

		if (cie2k)
			cg[0].pat[i].de = icmCIE2K(cg[0].pat[i].v, cg[1].pat[match[i]].v);
		else if (cie94)
			cg[0].pat[i].de = icmCIE94(cg[0].pat[i].v, cg[1].pat[match[i]].v);
		else
			cg[0].pat[i].de = icmLabDE(cg[0].pat[i].v, cg[1].pat[match[i]].v);

		cg[0].pat[i].ide[0] = fabs(cg[0].pat[i].v[0] - cg[1].pat[match[i]].v[0]);
		cg[0].pat[i].ide[1] = fabs(cg[0].pat[i].v[1] - cg[1].pat[match[i]].v[1]);
		cg[0].pat[i].ide[2] = fabs(cg[0].pat[i].v[2] - cg[1].pat[match[i]].v[2]);
	}

	/* Create sorted list, from worst to best. */
	if ((sort = (int *)malloc(sizeof(int) * cg[0].npat)) == NULL)
		error("Malloc failed - sort[]");
	for (i = 0; i < cg[0].npat; i++) 
		sort[i] = i;

#define HEAP_COMPARE(A,B) (cg[0].pat[A].de > cg[0].pat[B].de)
	HEAPSORT(int, sort, cg[0].npat);
#undef HEAP_COMPARE

	/* - - - - - - - - - - */
	/* Plot a dE histogram */
	if (dohisto) {
		double demax = -1e6, demin = 1e6;
		int maxbins = 50;		/* Maximum bins */
//		int minbins = 20;		/* Target minimum bins (depends on distribution) */ 
//		int mincount = 10;		/* Minimum number of points in a bin */
		int minbins = 10;		/* Target minimum bins (depends on distribution) */ 
		int mincount = 5;		/* Minimum number of points in a bin */
		double mbwidth;
		int nbins = 0;
		hbin *bins = NULL;
		pval **stpat;          /* Pointers to sorted cg[0].pat[] */
		double tval;
		double *x, *y;
		
		DBG(("Plotting histogram\n"))

		/* Figure out the range of dE's */
		for (i = 0; i < cg[0].npat; i++) {
			double de = cg[0].pat[i].de;

			if (de > demax)
				demax = de;
			if (de < demin)
				demin = de;
		}

		if (demax < 1e-6)
			error("histogram: dE range is too small to plot");

		/* Bin width that gives maxbins */
		mbwidth = demax / maxbins;
		
#ifdef NEVER
		/* Reduce mincount if needed to get minbins */
		if (cg[0].npat/minbins < mincount)
			mincount = cg[0].npat/minbins;
#endif

		if ((bins = (hbin *)calloc(maxbins, sizeof(hbin))) == NULL)
			error("malloc of histogram bins failed");

		if ((stpat = (pval **)malloc(sizeof(pval *) * cg[0].npat)) == NULL)
			error("Malloc failed - stpat[]");

		for (i = 0; i < cg[0].npat; i++)
			stpat[i] = &cg[0].pat[i];

	  	/* Sort the dE's */
#define HEAP_COMPARE(A,B) (A->de < B->de)
		HEAPSORT(pval *, stpat, cg[0].npat);
#undef HEAP_COMPARE

		/* Create bins and add points */
		bins[0].min = 0.0;
		for (nbins = i = 0; i < cg[0].npat && nbins < maxbins; i++) {
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

		/* Save points to a file ? */
		if (histoname[0] != '\000') {
			FILE *fp;

			if ((fp = fopen(histoname, "w")) == NULL)
				error("Opening '%s' for writing failed",histoname);

			fprintf(fp, "%s\t%s\n\n",cg[0].name, cg[1].name);

 			for (i = 0; i < nbins; i++) {
				fprintf(fp, "%f\t%f\n",0.5 * (bins[i].min + bins[i].max), bins[i].val);
			}

			if (fclose(fp) != 0)
				error("Closing '%s' failed",histoname);
		}

		free(y);
		free(x);
		free(bins);
		free(stpat);
	}

	/* - - - - - - - - - - */
	/* Figure out the report */
	{
		double merr = 0.0, aerr = 0.0;
		int n90;
		double merr90 = 0.0, aerr90 = 0.0;
		int n10;
		double merr10 = 0.0, aerr10 = 0.0;
		double rad;
		double aierr[3] = { 0.0, 0.0, 0.0 };
		double aixerr[3] = { 0.0, 0.0, 0.0 };
		double red[3] = { 1.0, 0.2, 0.2 };
		double green[3] = { 0.2, 1.0, 0.2 };
		double min[3], max[3];
		double col[3];

		if (dovrml) {
			double vol;
			int k;

			wrl = new_vrml(out_name, doaxes, (dovrml == 3 || dovrml == 4) ? vrml_rgb : vrml_lab);
			wrl->start_line_set(wrl, 0);

			for (j = 0; j < 3; j++) {
				min[j] = 1e6;
				max[j] = -1e6;
			}

			/* Get bounding box */
			for (i = 0; i < cg[0].npat; i++) {
				for (k = 0; k < 2; k++) {
					for (j = 0; j < 3; j++) {
						if (dovrml == 3 || dovrml == 4) {		/* RGB or YCC device plot */
							if (cg[k].pat[i].rgb[j] > max[j])
								max[j] = cg[k].pat[i].rgb[j];
							if (cg[k].pat[i].rgb[j] < min[j])
								min[j] = cg[k].pat[i].rgb[j];
						} else {
							if (cg[k].pat[i].v[j] > max[j])
								max[j] = cg[k].pat[i].v[j];
							if (cg[k].pat[i].v[j] < min[j])
								min[j] = cg[k].pat[i].v[j];
						}
					}
				}
			}

			for (vol = 1.0, j = 0; j < 3; j++) {
//printf("~1 size[%d] = %f\n",j, max[j] - min[j]);
				vol *= (max[j] - min[j]);
			}
			vol = sqrt(vol);
//printf("~1 vol = %f\n",vol);
			rad = 0.02 * vol/pow(cg[0].npat, 1.0/3.0);
//printf("~1 rad = %f\n",rad);

			if (dovrml == 3)	// Hack
				rad = 0.02;
			else if (dovrml == 4)	// Hack
				rad = 0.015;
		}

		if (dovrml && (dovrml == 3 || dovrml == 4)) {		/* RGB/YCC device plot */
			if (!cg[0].isrgb || !cg[1].isrgb)
				error("Both files must have RGB devices space for -d option");
		}

		/* Do overall results */
		for (i = 0; i < cg[0].npat; i++) {
			double de;

			if (cg[0].pat[i].og)		/* Skip out of gamut patches */
				continue;

			if (dosort)
				j = sort[i];
			else
				j = i;

			de = cg[0].pat[j].de;
			aerr += de;

			aierr[0] += cg[0].pat[j].ide[0];
			aierr[1] += cg[0].pat[j].ide[1];
			aierr[2] += cg[0].pat[j].ide[2];

			aixerr[0] += cg[0].pat[j].ixde[0];
			aixerr[1] += cg[0].pat[j].ixde[1];
			aixerr[2] += cg[0].pat[j].ixde[2];

			if (verb >= 2) {

				printf("%s%s%s: %f %f %f <=> %f %f %f  de %f\n",
					cg[0].pat[j].sid,
					cg[0].pat[j].loc[0] != '\000' ? " " : "",
					cg[0].pat[j].loc,
					cg[0].pat[j].v[0], cg[0].pat[j].v[1], cg[0].pat[j].v[2],
					cg[1].pat[match[j]].v[0], cg[1].pat[match[j]].v[1], cg[1].pat[match[j]].v[2],
					de);

#ifdef NEVER	/* Print XYZ as well */
				printf("  %f %f %f <=> %f %f %f\n",
					cg[0].pat[j].xyz[0], cg[0].pat[j].xyz[1], cg[0].pat[j].xyz[2],
					cg[1].pat[match[j]].xyz[0], cg[1].pat[match[j]].xyz[1], cg[1].pat[match[j]].xyz[2]);
#endif
			}

			if (de > merr)
				merr = de;

			if (dovrml) {
				if ((dovrml == 3 || dovrml == 4)) {		/* RGB/YCC device plot */
					double *val1, *val2;
					int k;

					if (dovrml == 3) {
						val1 = cg[0].pat[i].rgb; 
						val2 = cg[1].pat[match[i]].rgb;
					} else {
						val1 = cg[0].pat[i].ycc; 
						val2 = cg[1].pat[match[i]].ycc;
					}

					de = icmNorm33(val1, val2);

					if (de > 1e-6) {
						wrl->add_vertex(wrl, 0, val1);
						wrl->add_vertex(wrl, 0, val2);
					}

#ifdef NEVER	// Green target
					wrl->add_marker(wrl, val1, green, rad);

#else		// Natural color
					for (k = 0; k < 3; k++)
						col[k] = 0.3 + 0.7 * (cg[0].pat[i].rgb[k] - min[k])/(max[k] - min[k]);
					wrl->add_marker(wrl, val1, col, rad);
#endif

					wrl->add_marker_trans(wrl, val2, red, 0.3, rad * 0.99);

				} else {		/* PCS */
					if (de > 1e-6) {
						wrl->add_vertex(wrl, 0, cg[0].pat[i].v);
						wrl->add_vertex(wrl, 0, cg[1].pat[match[i]].v);
					}
					if (dovrml == 2) {
						wrl->add_marker(wrl, cg[0].pat[i].v, green, rad);
						wrl->add_marker(wrl, cg[1].pat[match[i]].v, red, rad);
					}
				}
			}

		}
		if (cg[0].nig > 0) {
			aerr /= (double)cg[0].nig;
			aierr[0] /= (double)cg[0].nig;
			aierr[1] /= (double)cg[0].nig;
			aierr[2] /= (double)cg[0].nig;

			aixerr[0] /= (double)cg[0].nig;
			aixerr[1] /= (double)cg[0].nig;
			aixerr[2] /= (double)cg[0].nig;
		}

		if (dovrml) {
			wrl->make_lines(wrl, 0, 2);
			wrl->del(wrl);
			wrl = NULL;
		}

		/* Do best 90% */
		n90 = (int)(cg[0].nig * 9.0/10.0 + 0.5);
		for (i = j = 0; i < cg[0].npat; i++) {
			double de = cg[0].pat[sort[i]].de;
			if (cg[0].pat[i].og)		/* Skip out of gamut */
				continue;
			if (j >= (cg[0].nig-n90)) {	/* If in top 90% of in gamut patches */
				aerr90 += de;
				if (de > merr90)
					merr90 = de;
			}
			j++;						/* Index of within gamut patches */
		}
		if (n90 > 0)
			aerr90 /= (double)n90;

		/* Do worst 10% */
		n10 = (int)(cg[0].nig * 1.0/10.0 + 0.5);
		for (i = j = 0; i < cg[0].npat; i++) {
			double de = cg[0].pat[sort[i]].de;
			if (cg[0].pat[i].og)		/* Skip out of gamut */
				continue;
			if (j < n10) {				/* If in worst 10% of in gamut patches */
				aerr10 += de;
				if (de > merr10)
					merr10 = de;
			}
			j++;
		}
		if (n10 > 0)
			aerr10 /= (double)n10;

		if (verb) {
			fprintf(verbo,"No of test patches in worst 10%% are = %d\n",n10);
			fprintf(verbo,"No of test patches in best 90%% are = %d\n",n90);
		}

		/* Single number report */
		if (dozrep != 0) {
			if (dozrep == 1) {
				printf("  %f\t", aerr);
			} else if (dozrep == 2) {
				printf("  %f\t", merr);
			}
			fflush(stdout);

		} else {
			printf("Verify results:\n");
			if (norm == 4)
				printf("  L*a*b* ref. = average XYZ %f %f %f\n",cg[0].w[0],cg[0].w[1],cg[0].w[2]);
			else if (norm == 1) {
				printf("  File 1 L* ref. Y %f\n", cg[0].w[1]);
				printf("  File 2 L* ref. Y %f\n", cg[1].w[1]);
			} else if (norm == 2) {
				printf("  File 1 L*a*b* ref. XYZ %f %f %f\n", cg[0].w[0],cg[0].w[1],cg[0].w[2]);
				printf("  File 2 L*a*b* ref. XYZ %f %f %f\n", cg[1].w[0],cg[1].w[1],cg[1].w[2]);
			} else if (norm == 3) {
				printf("  File 1 L* ref. X+Y+Z %f %f %f\n", cg[0].w[0],cg[0].w[1],cg[0].w[2]);
				printf("  File 2 L* ref. X+Y+Z %f %f %f\n", cg[1].w[0],cg[1].w[1],cg[1].w[2]);
			}
			printf("  Total errors%s:     peak = %f, avg = %f\n", cie2k ? " (CIEDE2000)" : cie94 ? " (CIE94)" : "", merr, aerr);
			printf("  Worst 10%% errors%s: peak = %f, avg = %f\n", cie2k ? " (CIEDE2000)" : cie94 ? " (CIE94)" : "", merr10, aerr10);
			printf("  Best  90%% errors%s: peak = %f, avg = %f\n", cie2k ? " (CIEDE2000)" : cie94 ? " (CIE94)" : "", merr90, aerr90);
			printf("  avg err X  %f, Y  %f, Z  %f\n", aixerr[0], aixerr[1], aixerr[2]);
			printf("  avg err L* %f, a* %f, b* %f\n", aierr[0], aierr[1], aierr[2]);
		}
	
		free(sort);
		free(match);
		free(cg[0].pat);
		free(cg[1].pat);
	}

	if (luo != NULL)
		luo->del(luo);
	if (xicco != NULL)
		xicco->del(xicco);		/* Expansion wrapper */
	if (icco != NULL)
		icco->del(icco);		/* Icc */
	if (fp != NULL)
		fp->del(fp);

	return 0;
}


