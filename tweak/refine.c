
/* 
 * Create/modify an abstract ICC transformation that will
 * correct measured color inacuracies of a proof reproduction.
 * Input is a list of target and actual CIE values.
 *
 * Author:  Graeme W. Gill
 * Date:    12/5/05
 * Version: 1.00
 *
 * Copyright 2005, 2008 Graeme W. Gill
 * Please refer to License.txt file for details.
 */

/* TTBD:
 *
 */

/* Basic idea:

	Given two .ti3 files (or equivalent CGATS files), one containing a
	spread of target patch values (XYZ, Lab or spectral), and the
	other containing the corresponding measured values, a PCS->PCS
	correction RSPL mapping is created to adjust for any innacuracy
	in the B2A table, which is then used to refine an existing abstract
	correction profile. (Any device values in the .ti3 tables are ignored.)

	The refined abstract profile can then be used to create an adjusted
	device profile, or an adjusted device link.

	Complications are:

		Measurements are absolute, while an abstract profile
		is relative.
		- use flag to mark this when abstract is used, or mark in profile header ?

		Corrections shouldn't go outside the target devices gamut, or
		they will lead to out of control regions on the gamut surface.
		- need output profile to clip changes to gamut surface.

		Refinement feedback could be unstable.
		- use a damping factor to improve stability.

	Currently the way out of gamut value are handled is to
	allow a full attempt at correction for the first round,
	and then constrain any subsequent corrections to be
	of no greater magnitude than that first correction.
	Subsequent corrections can change the direction of
	correction, but cannot increase its magnitude. This
	seems to work OK with light out of gamut colors,
	but dark CMYK out of gamut colors sometimes regress at
	the first step, possibly because of the gamut boundary
	topolgy in those regions.

	Interestingly, for CMYK the results are most stable (in simulation)
	when applied to the simple linked device link, and tends to
	be slightly unstable when applied to inverse A2B lookups
	that are used in profile and icclink -G. This could be a symptom
	of the black generation non-uniformity problem causing instability
	in the black inversion. Moving to optimised separation CMYK
	profile generation might overcome this problem.

	NOTE:- the current value for the rspl weak default weight seems OK
	for a reasonable number of points, but if refine was to be used for
	arbitrary tweaking, it should probably be made a tunable parameter
	that affects the radius of influence of each adjustment point.
*/


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "rspl.h"
#include "xicc.h"
#include "ui.h"

#define COMPLOOKUP	/* Compound with previous in ICM lookup rather than rspl */
#undef WARN_CLUT_CLIPPING   /* [Undef] Print warning if setting clut clips */
#undef DEBUG1		/* Print each correction value */
#undef DEBUG2		/* Print each value changed */
#undef DEBUG3		/* Trace history of particular points */

#define verbo stdout

#define RSPLFLAGS (0 /* | RSPL_2PASSSMTH | RSPL_EXTRAFIT2 */)

#define DEF_DAMP1 0.95		/* Initial */
#define DEF_DAMP2 0.70		/* Subsequent */
#define DEF_CLUTRES 33
#define GAMRES 10.0
#define SMOOTHF 0.3			/* RSPL smoothing factor */
#define AVGDEV 0.003		/* Average deviation of input values */
#define WWEIGHT 1.0			/* weak default function weight */
#define WHITEWEIGHT 5.0		/* Weight to put on discovered white correction */

#ifdef DEBUG3
/* Debug points of interest */
int poi[] = {
	8,
	15,
	10,
	14,
	380,
	274,
	562,
	172,
	510,
	331,
	297,
	494,
	102,
	18,
	13,
	6,
	305,
	61,
	455,
	63,
	461,
	68,
	7,
	369,
	211,
	42,
	427,
	113,
	204,
	224,
	334,
	28,
	175,
	330,
	273,
	376,
	318,
	44,
	57,
	469,
	9,
	85,
	278,
	414,
	124,
	-1
};

#endif /* DEBUG3 */

void usage(char *diag, ...) {
	fprintf(stderr,"Create abstract correction profile given table of absolute CIE correction values, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: refine [-options] cietarget ciecurrent [outdevicc] [inabs] outabs\n");
	fprintf(stderr," -v              Verbose\n");
	fprintf(stderr," -c              Create initial abstract correction profile\n");
	fprintf(stderr," -g              Don't impose output device gamut limit\n");
	fprintf(stderr," -r res          Set abstract profile clut resolution (default %d)\n",DEF_CLUTRES);
	fprintf(stderr," -d factor       Override default damping factor (default %f, then %f)\n",DEF_DAMP1,DEF_DAMP2);
	fprintf(stderr," -R              Aim for white point relative match rather than absolute\n");
	fprintf(stderr," -f [illum]      Use Fluorescent Whitening Agent compensation [opt. simulated inst. illum.:\n");
	fprintf(stderr,"                  M0, M1, M2, A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp]\n");
	fprintf(stderr," -i illum        Choose illuminant for computation of CIE XYZ from spectral data & FWA:\n");
	fprintf(stderr,"                  A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                  1931_2, 1964_10, 2012_2 (def), 2012_10, S&B 1955_2, shaw, J&V 1978_2 or file.cmf\n");
	fprintf(stderr," cietarget       Target CIE or spectral values, CGATS file (e.g. .ti3)\n");
	fprintf(stderr," ciecurrent      Actual CIE or spectral values, CGATS file (e.g. .ti3)\n");
	fprintf(stderr," [outdevicc]     Output device ICC profile to set gamut limit (not used if -g)\n");
	fprintf(stderr," [inabs]         Previous abstract correction ICC profile (not used if -c)\n");
	fprintf(stderr," outabs          Created/refined abstract correction ICC profile\n");
	exit(1);
}

/* ------------------------------------------- */
/* structure to support icc Lut initialisation/modification calbacks */

struct _callback {
	int verb;				/* Verbosity */
	int total, count, last;	/* Progress count information */
	rspl *r;				/* correction transform */
	icmLuBase *rd_luo;		/* Existing abstract profile (NULL if none) */
	gamut *dev_gam;			/* Gamut of output device (NULL if none) */
}; typedef struct _callback callback;


/* - - - - */
/*  clut  */

/* New CLUT table */
/* Correct for PCS errors */
void PCSp_PCSp(void *cntx, double *out, double *in) {
	callback *p = (callback *)cntx;
	co pp;

#ifdef DEBUG2
	printf("Got Lab in %f %f %f\n",in[0],in[1],in[2]);
#endif

	pp.p[0] = in[0];
	pp.p[1] = in[1];
	pp.p[2] = in[2];
	p->r->interp(p->r, &pp);				/* This correction */
	out[0] = pp.v[0];
	out[1] = pp.v[1];
	out[2] = pp.v[2];

#ifdef COMPLOOKUP
	/* Compound with previous correction */
	if (p->rd_luo != NULL) {
		p->rd_luo->lookup(p->rd_luo, out, out);			/* Previous correction */
	}
#endif

#ifdef DEBUG2
	printf("Got Lab out %f %f %f\n",out[0],out[1],out[2]);
	printf("\n");
#endif

	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = p->count * 100.0/p->total + 0.5;
		if (pc != p->last) {
			printf("%c%2d%%",cr_char,pc), fflush(stdout);
			p->last = pc;
		}
	}
}

/* ------------------------------------------- */

/* Patch value type */
typedef struct {
	char sid[50];		/* sample id */
	double v[3];		/* CIE value */
	double de;			/* Delta E */
} pval;

/* Weak default function */
static void wfunc(void *cbntx, double *out, double *in) {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
}

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	struct {
		char name[MAXNAMEL+1];	/* Patch filename  */
		int npat;				/* Number of patches */
		pval *pat;				/* patch values */
	} cg[2];					/* Target and current patch file information */
	char dev_name[MAXNAMEL+1];	/* Output device ICC filename for gamut */
	char rd_name[MAXNAMEL+1];	/* Abstract profile ICC to modify */
	char wr_name[MAXNAMEL+1];	/* Modified/created abstract profile ICC */

	int dorel = 0;				/* Do white point relative match */
	int *match;					/* Array mapping first list indexes to corresponding second */
	int fwacomp = 0;			/* FWA compensation on spectral ? */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType tillum = icxIT_none;	/* Target/simulated instrument illuminant */ 
	xspect cust_tillum, *tillump = NULL; /* Custom target/simulated illumination spectrum */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;					/* Custom illumination spectrum */
	icxObserverType obType = icxOT_CIE_2012_2;
	xspect custObserver[3];		/* If obType = icxOT_custom */
	callback cb;				/* Callback support stucture for setting abstract profile */

	icmFile *rd_fp = NULL;		/* Existing abstract profile to modify */
	icc *rd_icc = NULL;

	icmFile *wr_fp;				/* Modified/created abstract profile to write */
	icc *wr_icc;

	int verb = 0;
	int nogamut = 0;					/* Don't impose a gamut limit */
	int docreate = 0;					/* Create an initial abstract correction profile */
	int clutres = DEF_CLUTRES;			/* Output abstract profile clut resolution */
	double damp1 = DEF_DAMP1;			/* Initial damping factor */
	double damp2 = DEF_DAMP2;			/* Subsequent damping factor */
	double smoothf = SMOOTHF;			/* RSPL Smoothing factor */
	double avgdev[MXDO];				/* RSPL Average Deviation */
	double wweight = WWEIGHT;			/* weak default function weight */
	int whitepatch = -1;				/* Index of white patch */
	double merr = 0.0, aerr = 0.0;		/* Stats on color change */
	int i, j, e, n, rv = 0;

	error_program = argv[0];
	check_if_not_interactive();

	if (argc < 6)
		usage("Too few arguments");

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')	{	/* Look for any flags */
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
				usage("Usage requested");

			/* Verbosity */
			else if (argv[fa][1] == 'v') {
				verb = 1;
			}
			/* Create initial abstract correction profile */
			else if (argv[fa][1] == 'c') {
				docreate = 1;
			}
			/* Don't impose a gamut limit */
			else if (argv[fa][1] == 'g') {
				nogamut = 1;
			}
			/* Override the correction clut resolution */
			else if (argv[fa][1] == 'r') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to -r");
				clutres = atoi(na);
			}
			/* Override the damping factor */
			else if (argv[fa][1] == 'd') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to -d");
				damp2 = atof(na);
			}
			/* Aim for white point relative match */
			else if (argv[fa][1] == 'R') {
				dorel = 1;
			}

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
							usage("Unable to read target spectrum '%s'",na);

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
				fa = nfa;
				if (na == NULL) usage("Expect argument to -i");
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
						usage("Unable to read custom spectrum '%s'",na);

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
				fa = nfa;
				if (na == NULL) usage("Expected argument to -o");
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
						usage(0,"Failed to read custom observer CMF from -o file '%s'",na);
				}
			}

			else 
				usage("Unrecognised flag -%c",argv[fa][1]);
		} else
			break;
	}

	/* Grab all the filenames: */

	/* The two CIE value files */
	if (fa >= argc || argv[fa][0] == '-') usage("Expected cietarget file argument");
	strncpy(cg[0].name,argv[fa++],MAXNAMEL); cg[0].name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage("Expected ciecurrent file argument");
	strncpy(cg[1].name,argv[fa++],MAXNAMEL); cg[1].name[MAXNAMEL] = '\000';

	/* Optional output device name */
	if (nogamut == 0) {
		if (fa >= argc || argv[fa][0] == '-') usage("Expected outdevicc file argument");
		strncpy(dev_name,argv[fa++],MAXNAMEL); dev_name[MAXNAMEL] = '\000';
	}

	/* Optional input abstract profile name */
	if (docreate == 0) {
		if (fa >= argc || argv[fa][0] == '-') usage("Expected inabs file argument");
		strncpy(rd_name,argv[fa++],MAXNAMEL); rd_name[MAXNAMEL] = '\000';
	}

	/* Output abstract profile name */
	if (fa >= argc || argv[fa][0] == '-') usage("Expected outabs file argument");
	strncpy(wr_name,argv[fa++],MAXNAMEL); wr_name[MAXNAMEL] = '\000';

	/* ======================= */
	/* Open up each CIE file in turn, target then measured, */
	/* and read in the CIE values. */
	for (n = 0; n < 2; n++) {
		cgats *cgf = NULL;			/* cgats file data */
		int isLab = 0;				/* 0 if file CIE is XYZ, 1 if is Lab */
		int sidx;					/* Sample ID index */
		int xix, yix, zix;

		/* Open CIE target values */
		cgf = new_cgats();			/* Create a CGATS structure */
		cgf->add_other(cgf, ""); 	/* Allow any signature file */
	
		if (cgf->read_name(cgf, cg[n].name))
			error("CGATS file '%s' read error : %s",cg[n].name,cgf->err);
	
		if (cgf->ntables < 1)
			error ("Input file '%s' doesn't contain at least one table",cg[n].name);
	
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
			illum = icxIT_D50;
			obType = icxOT_CIE_1931_2;
		}
		if (spec && cgf->find_kword(cgf, 0, "SPECTRAL_BANDS") < 0)
			error ("No spectral data data found in file '%s' when spectral expected",cg[n].name);
	
		if (!spec && cgf->find_field(cgf, 0, "LAB_L") >= 0)
			isLab = 1;
		
		cg[n].npat = cgf->t[0].nsets;		/* Number of patches */
	
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
		 && (sidx = cgf->find_field(cgf, 0, "SAMPLE_NAME")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SAMPLE_LOC")) < 0)
			error("Input file '%s' doesn't contain field SAMPLE_ID, SampleName, Sample_Name, SAMPLE_NAME or SAMPLE_LOC",cg[n].name);
		if (cgf->t[0].ftype[sidx] != nqcs_t
		 && cgf->t[0].ftype[sidx] != cs_t)
			error("Sample ID/Name field isn't a quoted or non quoted character string");

		if (spec == 0) { 		/* Using instrument tristimulous value */

			if (isLab) {		/* Expect Lab */
				if ((xix = cgf->find_field(cgf, 0, "LAB_L")) < 0)
					error("Input file '%s' doesn't contain field LAB_L",cg[n].name);
				if (cgf->t[0].ftype[xix] != r_t)
					error("Field LAB_L is wrong type");
				if ((yix = cgf->find_field(cgf, 0, "LAB_A")) < 0)
					error("Input file '%s' doesn't contain field LAB_A",cg[n].name);
				if (cgf->t[0].ftype[yix] != r_t)
					error("Field LAB_A is wrong type");
				if ((zix = cgf->find_field(cgf, 0, "LAB_B")) < 0)
					error("Input file '%s' doesn't contain field LAB_B",cg[n].name);
				if (cgf->t[0].ftype[zix] != r_t)
					error("Field LAB_B is wrong type");

			} else { 		/* Expect XYZ */
				if ((xix = cgf->find_field(cgf, 0, "XYZ_X")) < 0)
					error("Input file '%s' doesn't contain field XYZ_X",cg[n].name);
				if (cgf->t[0].ftype[xix] != r_t)
					error("Field XYZ_X is wrong type");
				if ((yix = cgf->find_field(cgf, 0, "XYZ_Y")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Y",cg[n].name);
				if (cgf->t[0].ftype[yix] != r_t)
					error("Field XYZ_Y is wrong type");
				if ((zix = cgf->find_field(cgf, 0, "XYZ_Z")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Z",cg[n].name);
				if (cgf->t[0].ftype[zix] != r_t)
					error("Field XYZ_Z is wrong type");
			}

			for (i = 0; i < cg[n].npat; i++) {
				strcpy(cg[n].pat[i].sid, (char *)cgf->t[0].fdata[i][sidx]);
				cg[n].pat[i].v[0] = *((double *)cgf->t[0].fdata[i][xix]);
				cg[n].pat[i].v[1] = *((double *)cgf->t[0].fdata[i][yix]);
				cg[n].pat[i].v[2] = *((double *)cgf->t[0].fdata[i][zix]);
				if (!isLab) {
					cg[n].pat[i].v[0] /= 100.0;		/* Normalise XYZ to range 0.0 - 1.0 */
					cg[n].pat[i].v[1] /= 100.0;
					cg[n].pat[i].v[2] /= 100.0;
				}
				if (!isLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
					icmXYZ2Lab(&icmD50, cg[n].pat[i].v, cg[n].pat[i].v);
				}
			}

		} else { 		/* Using spectral data */
			int ii;
			xspect sp;
			char buf[100];
			int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
			xsp2cie *sp2cie;	/* Spectral conversion object */

			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_BANDS")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_BANDS");
			sp.spec_n = atoi(cgf->t[0].kdata[ii]);
			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_START_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_START_NM");
			sp.spec_wl_short = atof(cgf->t[0].kdata[ii]);
			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_END_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_END_NM");
			sp.spec_wl_long = atof(cgf->t[0].kdata[ii]);
			sp.norm = 100.0;

			/* Find the fields for spectral values */
			for (j = 0; j < sp.spec_n; j++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
				            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);

				if ((spi[j] = cgf->find_field(cgf, 0, buf)) < 0)
					error("Input file doesn't contain field %s",buf);
			}

			/* Figure out what sort of device it is */
			{
				int ti;
		
				if ((ti = cgf->find_kword(cgf, 0, "DEVICE_CLASS")) < 0)
					error ("Input file '%s' doesn't contain keyword DEVICE_CLASS",cg[n].name);
		
				if (strcmp(cgf->t[0].kdata[ti],"DISPLAY") == 0) {
					illum = icxIT_none;		/* Displays are assumed to be self luminous */
				}
			}

			/* Create a spectral conversion object */
			if ((sp2cie = new_xsp2cie(illum, 0.0, illum == icxIT_none ? NULL : &cust_illum,
			                          obType, custObserver, icSigLabData, icxClamp)) == NULL)
				error("Creation of spectral conversion object failed");

			if (fwacomp) {
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

				/* Track the maximum reflectance for any band to determine white. */
				/* This might silently fail, if there isn't white in the sampe set. */
				for (i = 0; i < cg[0].npat; i++) {
					for (j = 0; j < mwsp.spec_n; j++) {
						double rv = *((double *)cgf->t[0].fdata[i][spi[j]]);
						if (rv > mwsp.spec[j])
							mwsp.spec[j] = rv;
					}
				}

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
			}

			for (i = 0; i < cg[0].npat; i++) {

				strcpy(cg[n].pat[i].sid, (char *)cgf->t[0].fdata[i][sidx]);

				/* Read the spectral values for this patch */
				for (j = 0; j < sp.spec_n; j++) {
					sp.spec[j] = *((double *)cgf->t[0].fdata[i][spi[j]]);
				}

				/* Convert it to CIE space */
				sp2cie->convert(sp2cie, cg[n].pat[i].v, &sp);
			}

			sp2cie->del(sp2cie);		/* Done with this */

		}	/* End of reading in CGATs file */
		cgf->del(cgf);		/* Clean up */
	}

	/* Check that the number of test patches matches */
	if (cg[0].npat != cg[1].npat)
		error("Number of patches between '%s' and '%s' doesn't match",cg[0].name,cg[1].name);
	
	/* Create a list to map the second list (measured) of patches to the first (target) */
	if ((match = (int *)malloc(sizeof(int) * cg[0].npat)) == NULL)
		error("Malloc failed - match[]");
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

	/* Try and figure out which is the white patch */
	{
		double hL = -1.0;
		for (i = 0; i < cg[0].npat; i++) {
			if (cg[0].pat[i].v[0] > hL) {
				hL = cg[0].pat[i].v[0];
				whitepatch = i;
			}
		}
	}

	/* If we are aiming for a white point relative match, adjust the */
	/* measured and target values to have a D50 white point */
	if (dorel) {
		for (n = 0; n < 2; n++) {
			int wpix;			/* White patch index */
			double wp_xyz[3];
			icmXYZNumber wp;	/* White value */
			double mat[3][3];	/* Chromatic transform */
			
			if (n == 0)
				wpix = whitepatch;
			else
				wpix = match[whitepatch];
		

			/* Compute a chromatic correction matrix */
			icmLab2XYZ(&icmD50, wp_xyz, cg[n].pat[wpix].v);
			icmAry2XYZ(wp, wp_xyz);

			icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, wp, mat);

			for (i = 0; i < cg[n].npat; i++) {
				icmLab2XYZ(&icmD50, cg[n].pat[i].v, cg[n].pat[i].v);
				icmMulBy3x3(cg[n].pat[i].v, mat, cg[n].pat[i].v);
				icmXYZ2Lab(&icmD50, cg[n].pat[i].v, cg[n].pat[i].v);
//printf("Table %d, patch %d, Lab %f %f %f\n",n,i,cg[n].pat[i].v[0],cg[n].pat[i].v[1],cg[n].pat[i].v[2]);
			}
		}
	}

	/* Compute the delta E's just for information */
	for (i = 0; i < cg[0].npat; i++) {
		double de = icmLabDE(cg[0].pat[i].v, cg[1].pat[match[i]].v);
		cg[0].pat[i].de = de;
		if (de > merr)
			merr = de;
		aerr += de;
	}
	if (cg[0].npat > 0)
		aerr /= (double)cg[0].npat;

	if (verb) {
		fprintf(verbo,"No of correction patches = %d\n",cg[0].npat);
		fprintf(verbo,"Average dE = %f, Maximum dE = %f\n",aerr,merr);
		fprintf(verbo,"White patch assumed to be patch %s\n",cg[0].pat[whitepatch].sid);
	}

	/* ======================= */
	/* Possible limiting gamut */
	if (nogamut == 0) {
		icmFile *dev_fp;
		icc *dev_icc;
		xicc *dev_xicc;
		icxLuBase *dev_luo;
		icxInk ink;							/* Ink parameters */

		/* Open up the device ICC profile, so that we can create a gamut */
		/* and get an absolute PCS->device conversion */
		if ((dev_fp = new_icmFileStd_name(dev_name,"r")) == NULL)
			error ("Can't open file '%s'",dev_name);
	
		if ((dev_icc = new_icc()) == NULL)
			error("Creation of ICC object failed");
	
		/* Read header etc. */
		if ((rv = dev_icc->read(dev_icc,dev_fp,0)) != 0)
			error("Reading profile '%s' failed: %d, %s",dev_name,rv,dev_icc->err);
	
		/* Check that the profile is appropriate */
		if (dev_icc->header->deviceClass != icSigInputClass
		 && dev_icc->header->deviceClass != icSigDisplayClass
		 && dev_icc->header->deviceClass != icSigOutputClass)
			error("Device Profile '%s' isn't a device profile",dev_name);

		ink.tlimit = -1.0;		/* No ink limit by default */
		ink.klimit = -1.0;

		/* Wrap with an expanded icc */
		if ((dev_xicc = new_xicc(dev_icc)) == NULL)
			error ("Creation of xicc failed");

		/* Use a heuristic to guess the ink limit */
		icxGetLimits(dev_xicc, &ink.tlimit, &ink.klimit);
		ink.tlimit += 0.05;		/* allow a slight margine */

		if (verb)
			printf("Estimated Total inklimit is %f%%, Black %f%% \n",100.0 * ink.tlimit,ink.klimit < 0.0 ? 100.0 : 100.0 * ink.klimit);

		/* Get a expanded color conversion object suitable for gamut */
		if ((dev_luo = dev_xicc->get_luobj(dev_xicc, ICX_CLIP_NEAREST, icmFwd,
		     dorel ? icRelativeColorimetric : icAbsoluteColorimetric,
		     icSigLabData, icmLuOrdNorm, NULL, &ink)) == NULL)
			error ("%d, %s",dev_xicc->errc, dev_xicc->err);
	
		/* Creat a gamut surface */
		if ((cb.dev_gam = dev_luo->get_gamut(dev_luo, GAMRES)) == NULL)
			error ("%d, %s",dev_xicc->errc, dev_xicc->err);

		dev_luo->del(dev_luo);
		dev_xicc->del(dev_xicc);
		dev_icc->del(dev_icc);
		dev_fp->del(dev_fp);
	} else {
		cb.dev_gam = NULL;
	}

	/* ======================= */
	/* Open up the existing abstract profile that is to be refined. */
	if (docreate == 0) {
		if ((rd_fp = new_icmFileStd_name(rd_name,"r")) == NULL)
			error ("Can't open file '%s'",rd_name);
	
		if ((rd_icc = new_icc()) == NULL)
			error ("Creation of ICC object failed");
	
		/* Read header etc. */
		if ((rv = rd_icc->read(rd_icc,rd_fp,0)) != 0)
			error ("%d, %s",rv,rd_icc->err);
	
		if (rd_icc->header->deviceClass != icSigAbstractClass)
			error("Input Profile '%s' isn't abstract type",rd_name);

		if ((cb.rd_luo = rd_icc->get_luobj(rd_icc, icmFwd,
		        dorel ? icRelativeColorimetric : icAbsoluteColorimetric,
		        icSigLabData, icmLuOrdNorm)) == NULL)
				error ("%d, %s",rd_icc->errc, rd_icc->err);
	} else {
		cb.rd_luo = NULL;
	}

	/* ======================= */
	/* Create refining rspl */
	{
		cow *rp;		/* rspl setup points */
		int npnts = 0;	/* Total number of test points */
		int gres[MXDI];	/* rspl grid resolution */
		double damp;
		datai mn, mx;

		if ((rp = (cow *)malloc(sizeof(cow) * cg[0].npat)) == NULL)
			error("Malloc failed - rp[]");
		
		/* Create mapping points */
		for (i = 0; i < cg[0].npat; i++) {
			double temp[3];
			double ccor[3], cmag;				/* Current correction vector */
			double ncor[3], nmag;				/* New correction vector */

			/* Input is target [0] */
			for (j = 0; j < 3; j++)
				rp[i].p[j] = cg[0].pat[i].v[j];

			/* Cull out of range points */
			if (rp[i].p[0] < 0.0 || rp[i].p[0] > 100.0
			 || rp[i].p[1] < -127.0 || rp[i].p[1] > 127.0
			 || rp[i].p[2] < -127.0 || rp[i].p[2] > 127.0) {
#ifdef DEBUG1
			printf("Ignoring %f %f %f\n",rp[i].p[0],rp[i].p[1],rp[i].p[2]);
#endif
				continue;
			}
			
#ifdef DEBUG1
			printf("%d: Target        %f %f %f\n",i,rp[i].p[0],rp[i].p[1],rp[i].p[2]);
#endif

			damp = cb.rd_luo != NULL ? damp2 : damp1;
			ccor[0] = ccor[1] = ccor[2] = 0.0;
			cmag = 0.0;

			/* Lookup the current correction applied to the target */
			if (cb.rd_luo != NULL) {		/* Subsequent pass */
				double corval[3];
				cb.rd_luo->lookup(cb.rd_luo, corval, cg[0].pat[i].v);
				icmSub3(ccor, corval, cg[0].pat[i].v);
				cmag = icmNorm3(ccor);
#ifdef DEBUG1
				printf("%d: ccor          %f %f %f, mag %f\n",i, ccor[0],ccor[1],ccor[2],cmag);
#endif
			}

			/* Create a trial 100% full correction point */
			for (j = 0; j < 3; j++)
				rp[i].v[j] = ccor[j] + 2.0 * cg[0].pat[i].v[j] - cg[1].pat[match[i]].v[j];

			/* If a first pass and the target or the correction are out of gamut, */
			/* use a damping factor of 1.0 */
			if (cb.rd_luo == NULL
			 && cb.dev_gam != NULL
			 && cb.dev_gam->nradial(cb.dev_gam, temp, rp[i].p) > 1.0
			 && cb.dev_gam->nradial(cb.dev_gam, temp, rp[i].v) > 1.0) {
				damp = 1.0;
			}

#ifdef DEBUG1
			printf("%d: damp =         %f\n",i, damp);
#endif

			/* Create a new correction that does a damped correction to the current error */
			/* [0] = target, [1] = measured */
			for (j = 0; j < 3; j++)
				ncor[j] = ccor[j] + damp * (cg[0].pat[i].v[j] - cg[1].pat[match[i]].v[j]);
			nmag = icmNorm3(ncor);

#ifdef DEBUG1
			printf("%d: ncor          %f %f %f, mag %f\n",i, ncor[0],ncor[1],ncor[2],nmag);
#endif

			/* If this is not the first pass, limit the new correction */
			/* to be 1 + damp as big as the previous correction */
			if (cb.rd_luo != NULL) {
				if ((nmag/cmag) > (1.0 + damp2)) {
#ifdef DEBUG1
					printf("%d: Limited cor mag from %f to %f\n",i, nmag, (1.0 + damp2) * cmag);
#endif
					icmScale3(ncor, ncor, (1.0 + damp2) * cmag/nmag);
				}
			}

			/* Create correction point */
			for (j = 0; j < 3; j++)
				rp[i].v[j] = cg[0].pat[i].v[j] + ncor[j];

			/* If the target point or corrected point is likely to be outside */
			/* the gamut, limit the magnitude of the correction to be the same */
			/* as the previous correction. */ 
			if (cb.rd_luo != NULL && cb.dev_gam != NULL) {
				if (cb.dev_gam->nradial(cb.dev_gam, temp, rp[i].p) > 1.0
				 || cb.dev_gam->nradial(cb.dev_gam, temp, rp[i].v) > 1.0) {
#ifdef DEBUG1
					printf("%d: Limited cor mag from %f to %f\n",i, nmag, cmag);
#endif
					icmScale3(ncor, ncor, cmag/nmag);
				}
				/* Create correction point again */
				for (j = 0; j < 3; j++)
					rp[i].v[j] = cg[0].pat[i].v[j] + ncor[j];
			}

#ifdef DEBUG1
			printf("%d: Was           %f %f %f\n",i, cg[1].pat[match[i]].v[0], cg[1].pat[match[i]].v[1], cg[1].pat[match[i]].v[2]);
			printf("%d: Correction to %f %f %f\n",i, rp[i].v[0], rp[i].v[1], rp[i].v[2]);
#endif

#ifdef COMPLOOKUP
			/* Remove current correction from new change */
			for (j = 0; j < 3; j++)
				rp[i].v[j] -= ccor[j];
#endif
			/* Set weighting */
			if (i == whitepatch)
				rp[i].w = WHITEWEIGHT;
			else
				rp[i].w = 1.0;
			npnts++;

#ifdef DEBUG3
			{
				char fname[50], tmp[50];
				FILE *lf;
				int mi = match[i];
				double tig, cig, rig;  
				double vv[3], temp[3];
				double del[3], delt;
				double corrdel[3], corrdelt;
				double pcval[3], pcorrdel[3], pcorrdelt;

				for (j = 0;; j++) {
					if (poi[j] == (i+1) || poi[j] < 0)
						break;
				}
				if (poi[j] < 0) {
					continue;
				}

#ifdef COMPLOOKUP
				/* Compute total correction point */
				for (j = 0; j < 3; j++)
					vv[j] = rp[i].v[j] + ccor[j];
#else
				for (j = 0; j < 3; j++)
					vv[j] = rp[i].v[j];
#endif
				sprintf(fname,"patch%04d.log",i+1);
				if ((lf = fopen(fname, "a")) == NULL)
					error("Unable to open debug3 log file '%s'\n",fname);

	 			cig = cb.dev_gam->nradial(cb.dev_gam, temp, cg[1].pat[mi].v) - 1.0;
				if (cig > 0.0)
					sprintf(tmp, " OUT %f",cig);
				else
					sprintf(tmp, "");
				fprintf(lf,"Currently  %f %f %f%s\n", cg[1].pat[mi].v[0], cg[1].pat[mi].v[1], cg[1].pat[mi].v[2], tmp);

	 			tig = cb.dev_gam->nradial(cb.dev_gam, temp, cg[0].pat[i].v) - 1.0;
				if (tig > 0.0)
					sprintf(tmp, " OUT %f",tig);
				else
					sprintf(tmp, "");
				fprintf(lf,"Target     %f %f %f%s\n", cg[0].pat[i].v[0], cg[0].pat[i].v[1], cg[0].pat[i].v[2], tmp);

				icmSub3(del, cg[1].pat[mi].v, cg[0].pat[i].v);
				delt = icmNorm3(del);
				fprintf(lf,"DE         %f %f %f (%f)\n", del[0], del[1], del[2], delt);
				
				rig = cb.dev_gam->nradial(cb.dev_gam, temp, vv) - 1.0;
				if (rig > 0.0)
					sprintf(tmp, " OUT %f",rig);
				else
					sprintf(tmp, "");
				fprintf(lf,"Correction %f %f %f%s\n", vv[0], vv[1], vv[2], tmp);
				icmSub3(corrdel, vv, cg[0].pat[i].v);
				corrdelt = icmNorm3(corrdel);
				fprintf(lf,"CorrDelta  %f %f %f (%f)\n", corrdel[0], corrdel[1], corrdel[2], corrdelt);
				/* Note the previous correction we're compunded with */
				if (cb.rd_luo != NULL) {
					cb.rd_luo->lookup(cb.rd_luo, pcval, cg[0].pat[i].v);
					icmSub3(pcorrdel, pcval, cg[0].pat[i].v);
					pcorrdelt = icmNorm3(pcorrdel);
					fprintf(lf,"PrevCorrDelta %f %f %f (%f)\n", pcorrdel[0], pcorrdel[1], pcorrdel[2], pcorrdelt);
				}
				fprintf(lf,"\n");

				fclose(lf);
			}
#endif /* DEBUG3 */
		}

		/* Create refining rspl */
		mn[0] =   0.0, mn[1] = mn[2] = -128.0;			/* Allow for 16 bit grid range */
		mx[0] = 100.0, mx[1] = mx[2] =  (65535.0 * 255.0)/65280.0 - 128.0;
		cb.verb = verb;
		if ((cb.r = new_rspl(RSPL_NOFLAGS, 3, 3)) == NULL)
			error("new_rspl failed");

		for (e = 0; e < 3; e++)
			gres[e] = clutres;
		for (e = 0; e < 3; e++)
			avgdev[e] = AVGDEV;

		cb.r->fit_rspl_w_df(cb.r,
		           RSPLFLAGS			/* Extra flags */
		           | verb ? RSPL_VERBOSE : 0,
		           rp,					/* Test points */
		           npnts,				/* Number of test points */
		           mn, mx, gres,		/* Low, high, resolution of grid */
		           NULL, NULL,			/* Default data scale */
		           smoothf,				/* Smoothing */
		           avgdev,				/* Average Deviation */
		           NULL,				/* Grid width */
                   wweight,				/* weak default function weight */
				   NULL,				/* No context */
		           wfunc				/* Weak function */
		);
		if (verb) printf("\n");

		/* Report how good the fit is */
		if (verb) {
			co tco;	/* Test point */
			double maxe = -1e6, avge = 0.0;

			for (i = 0; i < npnts; i++) {
				double de;

				icmAry2Ary(tco.p, rp[i].p);
				cb.r->interp(cb.r, &tco);

				de = icmLabDE(tco.v, rp[i].v);
				if (de > maxe)
					maxe = de;
				avge += de;
			}
			avge /= (double)npnts;
			printf("Refining transform has error to defining points avg: %f, max %f\n",avge,maxe);
		}
		free(rp);
	}

	/* ======================= */
	/* Create new abstract ICC profile */
	if ((wr_fp = new_icmFileStd_name(wr_name,"w")) == NULL)
		error ("Can't open file '%s' for writing",wr_name);

	if ((wr_icc = new_icc()) == NULL)
		error ("Creation of write ICC object failed");

	/* Add all the tags required */

	/* The header: */
	{
		icmHeader *wh = wr_icc->header;

		/* Values that must be set before writing */
		wh->deviceClass     = icSigAbstractClass;
    	wh->colorSpace      = icSigLabData;
    	wh->pcs             = icSigLabData;
		if (dorel)
	    	wh->renderingIntent = icRelativeColorimetric;	/* White point relative */
		else
	    	wh->renderingIntent = icAbsoluteColorimetric;	/* Instrument reading based */
	}
	/* Profile Description Tag: */
	{
		icmTextDescription *wo;
		char *dst = "Argyll refine output";
		if ((wo = (icmTextDescription *)wr_icc->add_tag(
		           wr_icc, icSigProfileDescriptionTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* Copyright Tag: */
	{
		icmText *wo;
		char *crt = "Copyright the user who created it";
		if ((wo = (icmText *)wr_icc->add_tag(
		           wr_icc, icSigCopyrightTag,	icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

		wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->data, crt);		/* Copy the text in */
	}
	/* White Point Tag: */
	{
		icmXYZArray *wo;
		/* Note that tag types icSigXYZType and icSigXYZArrayType are identical */
		if ((wo = (icmXYZArray *)wr_icc->add_tag(
		           wr_icc, icSigMediaWhitePointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0] = icmD50;			/* So absolute/relative rendering is the same */
	}
	/* 16 bit pcs -> pcs lut: */
	{
		icmLut *wo;
		int flags = ICM_CLUT_SET_EXACT;	/* Assume we're setting from RSPL's */

		/* Intent 0 = default/perceptual */
		if ((wo = (icmLut *)wr_icc->add_tag(
		           wr_icc, icSigAToB0Tag,	icSigLut16Type)) == NULL) 
			error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

		wo->inputChan = 3;
		wo->outputChan = 3;
    	wo->clutPoints = clutres;
    	wo->inputEnt = 256;				/* Not actually used */
    	wo->outputEnt = 256;
		wo->allocate((icmBase *)wo);/* Allocate space */

		/* The matrix is only applicable to XYZ input space, */
		/* so it is not used here. */

		/* Use helper function to do the hard work. */
		if (cb.verb) {
			int extra;
			for (cb.total = 1, i = 0; i < 3; i++, cb.total *= wo->clutPoints)
				; 
			/* Add in cell center points */
			for (extra = 1, i = 0; i < wo->inputChan; i++, extra *= (wo->clutPoints-1))
				;
			cb.total += extra;
			cb.count = 0;
			cb.last = -1;
			printf(" 0%%"), fflush(stdout);
		}

#ifdef COMPLOOKUP
		/* Compound with previous correction */
		if (cb.rd_luo != NULL)
			flags = ICM_CLUT_SET_APXLS;	/* Won't be least squares, so do extra sampling */
#endif

		if (wo->set_tables(wo,
				flags,
				&cb,
				icSigLabData, 			/* Input color space */
				icSigLabData, 			/* Output color space */
				NULL,					/* Linear input transform Lab->Lab' (NULL = default) */
				NULL, NULL,				/* Use default Maximum range of Lab' values */
				PCSp_PCSp,				/* Lab' -> Lab' transfer function */
				NULL, NULL,				/* Use default Maximum range of Lab' values */
				NULL,					/* Linear output transform Lab'->Lab */
				NULL, NULL				/* Default SET_APXLS range */ 
		) != 0)
			error("Setting 16 bit Lab->Lab Lut failed: %d, %s",wr_icc->errc,wr_icc->err);

		if (verb)
			printf("\n");
#ifdef WARN_CLUT_CLIPPING
		if (wr_icc->warnc)
			warning("Values clipped in setting abstract LUT");
#endif /* WARN_CLUT_CLIPPING */
		if (verb)
			printf("Done filling abstract table\n");
	}
	/* Write the file out */
	if ((rv = wr_icc->write(wr_icc,wr_fp,0)) != 0)
		error ("Write file: %d, %s",rv,wr_icc->err);
	
	/* ======================================= */
	
	/* Clean everything up */
	wr_icc->del(wr_icc);
	wr_fp->del(wr_fp);

	if (docreate == 0) {
		cb.rd_luo->del(cb.rd_luo);
		rd_icc->del(rd_icc);
		rd_fp->del(rd_fp);
	}

	if (nogamut == 0) {
		cb.dev_gam->del(cb.dev_gam);
	}

	cb.r->del(cb.r);

	free(match);
	free(cg[0].pat);
	free(cg[1].pat);

	return 0;
}

