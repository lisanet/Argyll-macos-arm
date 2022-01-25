
/* Colorimeter Correction Matrix and */
/* Colorimeter Calibration Spectral Sample creation utility */

/* 
 * Argyll Color Management System
 * Author: Graeme W. Gill
 * Date:   19/8/2010
 *
 * Copyright 2010, 2011, 2013 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This program uses display measurements from a colorimeter and */
/* a spectrometer to create a correction matrix for a particular */
/* colorimeter/display combination,. */
/* or */
/* It uses display measurements from a spectrometer to create */
/* calibration samples that can be used with a Colorimeter that */
/* knowns its own spectral sensitivity curves (ie. X-Rite i1d3, Spyder 4). */

/* Based on spotread.c, illumread.c, dispcal.c */

/* 
	TTBD:

		Would be nice to have a way of not changing the target
		instruments absolute calibration.

		Would be nice to have a veryify option that produces
		a fit report of a matrix vs. the input files.

		Would be nice to have the option of procssing a Spyder 3 correction.txt file.
		(See post from umberto.guidali@tiscali.it)

		Would be nice to be able to use an i1D3 to correct other instruments,
		or an i1D3 created .ti3 as the reference. (can do this with .ti3 files).

		Would be nice to have an option of providing two ICC profiles,
		instead of using .ti3 files (?? How well would it work though ?)
 */


#undef DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "cgats.h"
#include "xspect.h"
#include "insttypes.h"
#include "conv.h"
#include "icoms.h"
#include "inst.h"
#include "ccast.h"
#include "ui.h"
#include "dispwin.h"
#include "webwin.h"
#include "dummywin.h"
#ifdef NT
# include "madvrwin.h"
#endif
#include "dispsup.h"
#include "ccss.h"
#include "ccmx.h"
#include "instappsup.h"
#ifdef ENABLE_USB
# include "spyd2.h"
#endif

#if defined (NT)
#include <conio.h>
#endif

#define DEFAULT_MSTEPS 1
#undef SHOW_WINDOW_ONFAKE	/* Display a test window up for a fake device */
#define COMPORT 1			/* Default com port 1..4 */


#if defined(__APPLE__) && defined(__POWERPC__)
/* Workaround for a ppc gcc 3.3 optimiser bug... */
static int gcc_bug_fix(int i) {
	static int nn;
	nn += i;
	return nn;
}
#endif	/* APPLE */

/* Invoke with -dfake for testing with a fake device. */
/* Invoke with -dFAKE for automatic creation of test matrix. */
/* Will use a fake.icm/.icc profile if present, or a built in fake */
/* device behaviour if not. */

void
/* Flag = 0x0000 = default */
/* Flag & 0x0001 = list ChromCast's */
/* Flag & 0x0002 = list Technology choice */
/* Flag & 0x1xxx = -S flag */
usage(int flag, char *diag, ...) {
	disppath **dp;
	icompaths *icmps = new_icompaths(0);
	inst2_capability cap = 0;

	fprintf(stderr,"Create CCMX or CCSS, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: ccxxmake -t dtech [-options] output.%s\n", flag & 0x1000 ? "ccss" : "ccmx");
	fprintf(stderr," -v                Verbose mode\n");
	fprintf(stderr," -S                Create CCSS rather than CCMX\n");
	fprintf(stderr," -f ref.ti3[,targ.ti3]  Create from one or two .ti3 files rather than measure.\n");
#if defined(UNIX_X11)
	fprintf(stderr," -display displayname   Choose X11 display name\n");
	fprintf(stderr," -d n[,m]          Choose the display n from the following list (default 1)\n");
	fprintf(stderr,"                   Optionally choose different display m for VideoLUT access\n"); 
#else
	fprintf(stderr," -d n              Choose the display from the following list (default 1)\n");
#endif
	dp = get_displays();
	if (dp == NULL || dp[0] == NULL)
		fprintf(stderr,"    ** No displays found **\n");
	else {
		int i;
		for (i = 0; ; i++) {
			if (dp[i] == NULL)
				break;
			fprintf(stderr,"    %d name = '%s'\n",i+1,dp[i]->name);
			fprintf(stderr,"    %d = '%s'\n",i+1,dp[i]->description);
		}
	}
	free_disppaths(dp);
	fprintf(stderr," -dweb[:port]      Display via a web server at port (default 8080)\n");
	fprintf(stderr," -dcc[:n]          Display via n'th ChromeCast (default 1, ? for list)\n");
	if (flag & 0x001) {
		ccast_id **ids;
		if ((ids = get_ccids()) == NULL) {
			fprintf(stderr,"    ** Error discovering ChromCasts **\n");
		} else {
			if (ids[0] == NULL)
				fprintf(stderr,"    ** No ChromCasts found **\n");
			else {
				int i;
				for (i = 0; ids[i] != NULL; i++)
					fprintf(stderr,"    %d = '%s'\n",i+1,ids[i]->name);
				free_ccids(ids);
			}
		}
	}
#ifdef NT
	fprintf(stderr," -d madvr          Display via MadVR Video Renderer\n");
#endif
	fprintf(stderr," -d dummy          Dummy (non-existant, invisible) display\n");
//	fprintf(stderr," -d fake           Use a fake (ICC profile) display device for testing, fake%s if present\n",ICC_FILE_EXT);
	fprintf(stderr," -p                Use telephoto mode (ie. for a projector, if available)\n");
	fprintf(stderr," -a                Use ambient measurement mode (ie. for a projector, if available)\n");
	cap = inst_show_disptype_options(stderr, " -y c|l                 ", icmps, 1, 0);
	fprintf(stderr," -z disptype       Different display type for spectrometer (see -y)\n");
	fprintf(stderr," -P ho,vo,ss[,vs]  Position test window and scale it\n");
	fprintf(stderr,"                   ho,vi: 0.0 = left/top, 0.5 = center, 1.0 = right/bottom etc.\n");
	fprintf(stderr,"                   ss: 0.5 = half, 1.0 = normal, 2.0 = double etc.\n");
	fprintf(stderr," -F                Fill whole screen with black background\n");
#if defined(UNIX_X11)
	fprintf(stderr," -n                Don't set override redirect on test window\n");
#endif
	fprintf(stderr," -N                Disable initial calibration of instrument if possible\n");
	fprintf(stderr," -H                Use high resolution spectrum mode (if available)\n");
//	fprintf(stderr," -V                Use adaptive measurement mode (if available)\n");
	fprintf(stderr," -C \"command\"      Invoke shell \"command\" each time a color is set\n");
	fprintf(stderr," -M \"command\"      Invoke shell \"command\" each time a color is measured\n");
	fprintf(stderr," -o observ         Choose CIE Observer for CCMX spectrometer data:\n");
	fprintf(stderr,"                    1931_2 (def), 1964_10, 2012_2, 2012_10, S&B 1955_2, shaw, J&V 1978_2 or file.cmf\n");
	fprintf(stderr," -s steps          Override default patch sequence combination steps  (default %d)\n",DEFAULT_MSTEPS);
	fprintf(stderr," -W n|h|x          Override serial port flow control: n = none, h = HW, x = Xon/Xoff\n");
	fprintf(stderr," -D [level]        Print debug diagnostics to stderr\n");
	fprintf(stderr," -E desciption     Override the default overall description\n");
	fprintf(stderr," -I displayname    Set display make and model description (optional)\n");
	if (flag & 0x0002) {
		int i;
		disptech_info *list = disptech_get_list();
		for (i = 0; list[i].dtech != disptech_end; i++)
			fprintf(stderr," %s %s               %s\n",i == 0 ? "-t" : "  ", list[i].lsel,list[i].desc);
	} else {
		fprintf(stderr," -t dtech          Set display technology type\n");
		fprintf(stderr,"                    (Use -?? to list technology choices)\n");
	}
	fprintf(stderr," -U c              Set UI selection character(s)\n");
	fprintf(stderr," -Y r|n            Set or override refresh/non-refresh display type\n");
	fprintf(stderr," -Y R:rate         Override measured refresh rate with rate Hz\n");
	fprintf(stderr," -Y A              Use non-adaptive integration time mode (if available).\n");
	fprintf(stderr," correction.ccmx | calibration.ccss\n");
	fprintf(stderr,"                   File to save result to\n");
	if (icmps != NULL)
		icmps->del(icmps);
	exit(1);
}

typedef double ary3[3];

int main(int argc, char *argv[]) {
	int i,j;
	int fa, nfa, mfa;					/* current argument we're looking at */
	disppath *disp = NULL;				/* Display being used */
	double hpatscale = 1.0, vpatscale = 1.0;	/* scale factor for test patch size */
	double ho = 0.0, vo = 0.0;			/* Test window offsets, -1.0 to 1.0 */
	int fullscreen = 0;            		/* NZ if whole screen should be filled with black */
	int verb = 0;
	int debug = 0;
	int doccss = 0;						/* Create CCSS rather than CCMX */
	int fake = 0;						/* Use the fake device for testing, 2 for auto */
	int faketoggle = 0;					/* Toggle fake between "colorimeter" and "spectro" */
	int fakeseq = 0;					/* Fake auto CCMX sequence */
	int spec = 0;						/* Need spectral data to implement option */
	icxObserverType obType = icxOT_CIE_1931_2;
	xspect custObserver[3];				/* If obType = icxOT_custom */
	int override = 1;					/* Override redirect on X11 */
	icompaths *icmps = NULL;			/* Ports to choose from */
	int comno = COMPORT;				/* COM port used */
	flow_control fc = fc_nc;			/* Default flow control */
	int highres = 0;					/* High res mode if available */
	int ditype = 0;						/* Display kind selector, 0 = default */
	int sditype = -1;					/* Spectro display kind, -1 = use ditype */
	int refrmode = -1;					/* Refresh mode */
	double refrate = 0.0;				/* 0.0 = default, > 0.0 = override refresh rate */ 
	int cbid = 0;						/* Calibration base display mode ID */
	int nadaptive = 0;					/* Use non-adaptive mode if available */
	int tele = 0;						/* NZ if telephoto mode */
	int ambient = 0;					/* NZ if ambient mode */
	int noinitcal = 0;					/* Disable initial calibration */
	int webdisp = 0;					/* NZ for web display, == port number */
	int ccdisp = 0;			 			/* NZ for ChromeCast, == list index */
	ccast_id **ccids = NULL;
	ccast_id *ccid = NULL;
#ifdef NT
	int madvrdisp = 0;					/* NZ for MadVR display */
#endif
	int dummydisp = 0;					/* NZ for dummy display */
	char *ccallout = NULL;				/* Change color Shell callout */
	char *mcallout = NULL;				/* Measure color Shell callout */
	int msteps = DEFAULT_MSTEPS;		/* Patch surface size */
	int npat = 0;						/* Number of patches/colors */
	ary3 *refs = NULL;					/* Reference XYZ values */
	int gotref = 0;
	char *refname = NULL;				/* Name of reference instrument */
	char *reffile = NULL;				/* Name of reference file */
	ary3 *cols = NULL;					/* Colorimeter XYZ values */
	int gotcol = 0;
	char *colname = NULL;				/* Name of colorimeter instrument */
	char *colfile = NULL;				/* Name of colorimeter file */
	col *rdcols = NULL;					/* Internal storage of all the patch colors */
	int saved = 0;						/* Saved result */
	char innames[2][MAXNAMEL+1] = { "\000", "\000" };  /* .ti3 input names */
	char outname[MAXNAMEL+5+1] = "\000";  /* ccmx output file name */
	char *description = NULL;			/* Given overall description */
	char *displayname = NULL;			/* Given display name */
	disptech_info *dtinfo = NULL;		/* Display technology */
	char *uisel = NULL;					/* UI selection letters */
	int uflag = 0;						/* usage flag */
	int rv;

	set_exe_path(argv[0]);				/* Set global exe_path and error_program */
	check_if_not_interactive();

	/* Process the arguments */
	mfa = 0;        /* Minimum final arguments */
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
					usage(uflag | 2, "Extended usage requested");
				usage(uflag | 0, "Usage requested");

			} else if (argv[fa][1] == 'v') {
				verb = 1;
				g_log->verb = verb;

			} else if (argv[fa][1] == 'S') {
				doccss = 1;
				uflag |= 0x1000;

			} else if (argv[fa][1] == 'f') {
				char *cna, *f1 = NULL;
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Expect argument to input file flag -f");

				if ((cna = strdup(na)) == NULL)
					error("Malloc failed");

				/* If got just one file - enough for CCSS */
				if ((f1 = strchr(cna, ',')) == NULL) {
					strncpy(innames[0],cna,MAXNAMEL-1); innames[0][MAXNAMEL-1] = '\000';
					free(cna);

				/* Got two files - needed for CCMX */
				} else {
					*f1++ = '\000';
					strncpy(innames[0],cna,MAXNAMEL-1); innames[0][MAXNAMEL-1] = '\000';
					strncpy(innames[1],f1,MAXNAMEL-1); innames[1][MAXNAMEL-1] = '\000';
					free(cna);
				}

			/* Display number */
			} else if (argv[fa][1] == 'd') {
				if (strncmp(na,"web",3) == 0
				 || strncmp(na,"WEB",3) == 0) {
					webdisp = 8080;
					if (na[3] == ':') {
						webdisp = atoi(na+4);
						if (webdisp == 0 || webdisp > 65535)
							usage(uflag | 0,"Web port number must be in range 1..65535");
					}
					fa = nfa;
				} else if (strncmp(na,"cc",2) == 0
				 || strncmp(na,"CC",2) == 0) {
					ccdisp = 1;
					if (na[2] == ':') {
						if (na[3] < '0' || na[3] > '9')
							usage(uflag | 0x0001,"Available ChromeCasts");

						ccdisp = atoi(na+3);
						if (ccdisp <= 0)
							usage(uflag | 0,"ChromCast number must be in range 1..N");
					}
					fa = nfa;
#ifdef NT
				} else if (strncmp(na,"madvr",5) == 0
				 || strncmp(na,"MADVR",5) == 0) {
					madvrdisp = 1;
					fa = nfa;
#endif
				} else if (strncmp(na,"dummy",5) == 0
				 || strncmp(na,"DUMMY",5) == 0) {
					dummydisp = 1;
					fa = nfa;
				} else {
#if defined(UNIX_X11)
					int ix, iv;

					if (strcmp(&argv[fa][2], "isplay") == 0 || strcmp(&argv[fa][2], "ISPLAY") == 0) {
						if (++fa >= argc || argv[fa][0] == '-') usage(uflag | 0,"Parameter expected following -display");
						setenv("DISPLAY", argv[fa], 1);
					} else {
						if (na == NULL) usage(uflag | 0,"Parameter expected following -d");
						fa = nfa;
						if (strcmp(na,"fake") == 0 || strcmp(na,"FAKE") == 0) {
							fake = 1;
							if (strcmp(na,"FAKE") == 0)
								fakeseq = 1;
						} else {
							if (sscanf(na, "%d,%d",&ix,&iv) != 2) {
								ix = atoi(na);
								iv = 0;
							}
							if (disp != NULL)
								free_a_disppath(disp);
							if ((disp = get_a_display(ix-1)) == NULL)
								usage(uflag | 0,"-d parameter %d out of range",ix);
							if (iv > 0)
								disp->rscreen = iv-1;
						}
					}
#else
					int ix;
					if (na == NULL) usage(uflag | 0,"Parameter expected following -d");
					fa = nfa;
					if (strcmp(na,"fake") == 0 || strcmp(na,"FAKE") == 0) {
						fake = 1;
						if (strcmp(na,"FAKE") == 0)
							fakeseq = 1;
					} else {
						ix = atoi(na);
						if (disp != NULL)
							free_a_disppath(disp);
						if ((disp = get_a_display(ix-1)) == NULL)
							usage(uflag | 0,"-d parameter %d out of range",ix);
					}
#endif
				}
#if defined(UNIX_X11)
			} else if (argv[fa][1] == 'n') {
				override = 0;
#endif /* UNIX */

			/* COM port  */
			} else if (argv[fa][1] == 'c') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Paramater expected following -c");
				comno = atoi(na);
				if (comno < 1 || comno > 40) usage(uflag | 0,"-c parameter %d out of range",comno);

			/* Telephoto */
			} else if (argv[fa][1] == 'p') {
				tele = 1;
				ambient = 0;

			/* Ambient */
			} else if (argv[fa][1] == 'a') {
				ambient = 1;
				tele = 0;

			/* Display type */
			} else if (argv[fa][1] == 'y') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Parameter expected after -y");
				ditype = na[0];
				if (ditype == '_' && na[1] != '\000')
					ditype = ditype << 8 | na[1];

				/* For ccss, set a default */
				if (na[0] == 'r') {
					refrmode = 1;
				} else if (na[0] == 'n') {
					refrmode = 0;
				}

			/* Spectro display type */
			} else if (argv[fa][1] == 'z') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Parameter expected after -z");
				sditype = na[0];

			/* Test patch offset and size */
			} else if (argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Parameter expected after -P");
				if (sscanf(na, " %lf,%lf,%lf,%lf ", &ho, &vo, &hpatscale, &vpatscale) == 4) {
					;
				} else if (sscanf(na, " %lf,%lf,%lf ", &ho, &vo, &hpatscale) == 3) {
					vpatscale = hpatscale;
				} else {
					usage(uflag | 0,"-P parameter '%s' not recognised",na);
				}
				if (ho < 0.0 || ho > 1.0
				 || vo < 0.0 || vo > 1.0
				 || hpatscale <= 0.0 || hpatscale > 50.0
				 || vpatscale <= 0.0 || vpatscale > 50.0)
					usage(uflag | 0,"-P parameters %f %f %f %f out of range",ho,vo,hpatscale,vpatscale);
				ho = 2.0 * ho - 1.0;
				vo = 2.0 * vo - 1.0;

			/* Full screen black background */
			} else if (argv[fa][1] == 'F') {
				fullscreen = 1;

			/* No initial calibration */
			} else if (argv[fa][1] == 'N') {
				noinitcal = 1;

			/* High res spectral mode */
			} else if (argv[fa][1] == 'H') {
				highres = 1;

			/* Adaptive mode - now default, so flag is deprecated */
			} else if (argv[fa][1] == 'V') {
				warning("dispcal -V flag is deprecated");

			/* Spectral Observer type (only relevant for CCMX) */
			} else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Parameter expecte after -o");
				if (strcmp(na, "1931_2") == 0) {			/* Classic 2 degree */
					spec = 2;
					obType = icxOT_CIE_1931_2;
				} else if (strcmp(na, "1964_10") == 0) {	/* Classic 10 degree */
					spec = 2;
					obType = icxOT_CIE_1964_10;
				} else if (strcmp(na, "2012_2") == 0) {		/* Latest 2 degree */
					spec = 2;
					obType = icxOT_CIE_2012_2;
				} else if (strcmp(na, "2012_10") == 0) {	/* Latest 10 degree */
					spec = 2;
					obType = icxOT_CIE_2012_10;
				} else if (strcmp(na, "1955_2") == 0) {		/* Stiles and Burch 1955 2 degree */
					spec = 2;
					obType = icxOT_Stiles_Burch_2;
				} else if (strcmp(na, "1978_2") == 0) {		/* Judd and Voss 1978 2 degree */
					spec = 2;
					obType = icxOT_Judd_Voss_2;
				} else if (strcmp(na, "shaw") == 0) {		/* Shaw and Fairchilds 1997 2 degree */
					spec = 2;
					obType = icxOT_Shaw_Fairchild_2;
				} else {	/* Assume it's a filename */
					obType = icxOT_custom;
					if (read_cmf(custObserver, na) != 0)
						usage(uflag | 0,"Failed to read custom observer CMF from -o file '%s'",na);
				}

			} else if (argv[fa][1] == 's') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Parameter expecte after -s");
				msteps = atoi(na);
				if (msteps < 1 || msteps > 16)
					usage(uflag | 0,"-s parameter value %d is outside the range 1 to 16",msteps);

			/* Change color callout */
			} else if (argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Parameter expected after -C");
				ccallout = na;

			/* Measure color callout */
			} else if (argv[fa][1] == 'M') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Parameter expected after -M");
				mcallout = na;

			/* Serial port flow control */
			} else if (argv[fa][1] == 'W') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Paramater expected following -W");
				if (na[0] == 'n' || na[0] == 'N')
					fc = fc_None;
				else if (na[0] == 'h' || na[0] == 'H')
					fc = fc_Hardware;
				else if (na[0] == 'x' || na[0] == 'X')
					fc = fc_XonXOff;
				else
					usage(uflag | 0,"-W parameter '%c' not recognised",na[0]);

			} else if (argv[fa][1] == 'D') {
				debug = 1;
				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					debug = atoi(na);
					fa = nfa;
				}
				g_log->debug = debug;

			} else if (argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Expect argument to display description flag -I");
				displayname = strdup(na);

			} else if (argv[fa][1] == 't') {
				fa = nfa;
				if (na == NULL) usage(uflag | 2,"Expect argument to display technology flag -t");
				dtinfo = disptech_get_list();
				if (na[1] != '\000')
					usage(uflag | 2,"Expect single character argument to display technology flag -t");
				if ((dtinfo = disptech_select(dtinfo, na[0])) == NULL)
					usage(uflag | 2,"-t parameter '%c' not recognized",na[0]);

			/* Copyright string */
			} else if (argv[fa][1] == 'E') {
				fa = nfa;
				if (na == NULL) usage(uflag | 0,"Expect argument to overall description flag -E");
				description = strdup(na);

			/* Extra flags */
			} else if (argv[fa][1] == 'Y') {
				if (na == NULL)
					usage(uflag | 0,"Flag '-Y' expects extra flag");
			
				if (na[0] == 'r') {
					refrmode = 1;
				} else if (na[0] == 'n') {
					refrmode = 0;
				} else if (na[0] == 'R') {
					if (na[1] != ':')
						usage(uflag | 0,"-Y R:rate syntax incorrect");
					refrate = atof(na+2);
					if (refrate < 5.0 || refrate > 150.0)
						usage(uflag | 0,"-Y R:rate %f Hz not in valid range",refrate);
				} else if (na[0] == 'A') {
					nadaptive = 1;
				} else {
					usage(uflag | 0,"Flag '-Z %c' not recognised",na[0]);
				}
				fa = nfa;

			/* UI selection character */
			} else if (argv[fa][1] == 'U') {
				fa = nfa;
				if (na == NULL || na[0] == '\000') usage(uflag | 0,"Expect argument to flag -U");
				uisel = na;
				for (i = 0; uisel[i] != '\000'; i++) {
					if (!( (uisel[i] >= '0' && uisel[i] <= '9')
					    || (uisel[i] >= 'A' && uisel[i] <= 'Z')
					    || (uisel[i] >= 'a' && uisel[i] <= 'z'))) {
						usage(uflag | 0,"-U character(s) must be 0-9,A-Z,a-z");
					}
				}

			} else 
				usage(uflag | 0,"Flag '-%c' not recognised",argv[fa][1]);
		}
		else
			break;
	}

	/* Get the output ccmx file name argument */
	if (fa >= argc)
		usage(uflag | 0,"Output filname expected");

	strncpy(outname,argv[fa++],MAXNAMEL-1); outname[MAXNAMEL-1] = '\000';
	if (strrchr(outname, '.') == NULL)	/* no extension */
		strcat(outname, doccss ? ".ccss" : ".ccmx");

	if (fakeseq && doccss)
		error("Fake CCSS test not implemeted");

	printf("\n");

	if (dtinfo == NULL)
		error("Display technology (-t) must be set");

	/* CCSS: See if we're working from a .ti3 file */
	if (doccss && innames[0][0] != '\000') {
		cgats *cgf = NULL;			/* cgats file data */
		int sidx;					/* Sample ID index */
		int ii, ti;
		char buf[100];
		int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
		xspect sp, *samples = NULL;
		ccss *cc;
		double bigv = -1e60;

		if (innames[1][0] != '\000')
			warning("second -f .ti3 '%s' ignored! (not needed for CCSS)", innames[1]);

		/* Open spectral values file */
		cgf = new_cgats();			/* Create a CGATS structure */
		cgf->add_other(cgf, ""); 	/* Allow any signature file */
	
		if (cgf->read_name(cgf, innames[0]))
			error("CGATS file '%s' read error : %s",innames[0],cgf->err);
	
		if (cgf->ntables < 1)
			error ("Input file '%s' doesn't contain at least one table",innames[0]);
	
		if ((npat = cgf->t[0].nsets) <= 0)
			error("No sets of data in file '%s'",innames[0]);

		if ((samples = (xspect *)malloc(npat * sizeof(xspect))) == NULL)
				error("malloc failed");

		if ((ii = cgf->find_kword(cgf, 0, "TARGET_INSTRUMENT")) < 0)
			error ("Can't find keyword TARGET_INSTRUMENT in '%s'",innames[0]);

		if ((ti = cgf->find_kword(cgf, 0, "DISPLAY_TYPE_REFRESH")) >= 0) {
			if (stricmp(cgf->t[0].kdata[ti], "YES") == 0)
				refrmode = 1;
			else if (stricmp(cgf->t[0].kdata[ti], "NO") == 0)
				refrmode = 0;
		}

		if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_BANDS")) < 0)
			error ("Input file '%s' doesn't contain keyword SPECTRAL_BANDS",innames[0]);
		sp.spec_n = atoi(cgf->t[0].kdata[ii]);
		if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_START_NM")) < 0)
			error ("Input file '%s' doesn't contain keyword SPECTRAL_START_NM",innames[0]);
		sp.spec_wl_short = atof(cgf->t[0].kdata[ii]);
		if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_END_NM")) < 0)
			error ("Input file '%s' doesn't contain keyword SPECTRAL_END_NM",innames[0]);
		sp.spec_wl_long = atof(cgf->t[0].kdata[ii]);
		sp.norm = 1.0;		/* We assume emssive */

		/* Find the fields for spectral values */
		for (j = 0; j < sp.spec_n; j++) {
			int nm;
	
			/* Compute nearest integer wavelength */
			nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
			            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
			
			sprintf(buf,"SPEC_%03d",nm);

			if ((spi[j] = cgf->find_field(cgf, 0, buf)) < 0)
				error("Input file '%s' doesn't contain field %s",innames[0],buf);

			if (cgf->t[0].ftype[spi[j]] != r_t)
				error("Field %s is wrong type - expect float",buf);
		}

		/* Transfer all the spectral values */
		for (i = 0; i < npat; i++) {

			XSPECT_COPY_INFO(&samples[i], &sp);
	
			for (j = 0; j < sp.spec_n; j++) {
				samples[i].spec[j] = *((double *)cgf->t[0].fdata[i][spi[j]]);
			}
		}
		cgf->del(cgf);		/* Clean up */
		cgf = NULL;

		if (description == NULL) {
			char *disp = displayname != NULL ? displayname : dtinfo->desc;
			char *tt = "CCSS for ";
			if ((description = malloc(strlen(disp) + strlen(tt) + 1)) == NULL)
				error("Malloc failed");
			strcpy(description, tt);
			strcat(description, disp);
		}

		/* See what the highest value is */
		for (i = 0; i < npat; i++) {	/* For all grid points */

			for (j = 0; j < samples[i].spec_n; j++) {
				if (samples[i].spec[j] > bigv)
					bigv = samples[i].spec[j];
			}
		}
			
		/* Normalize the values */
		for (i = 0; i < npat; i++) {	/* For all grid points */
			double scale = 100.0;

			for (j = 0; j < samples[i].spec_n; j++)
				samples[i].spec[j] *= scale / bigv;
		}
			
		if (refrmode < 0)
			error("The display refresh mode is not known - use the -Y flag");

		if ((cc = new_ccss()) == NULL)
			error("new_ccss() failed");

		if (cc->set_ccss(cc, "Argyll ccxxmake", NULL, description, displayname,
		                 dtinfo->dtech, refrmode, uisel, refname, 0, samples, npat)) {
			error("set_ccss failed with '%s'\n",cc->err);
		}
		if(cc->write_ccss(cc, outname))
			printf("\nWriting CCSS file '%s' failed with '%s\n",outname, cc->err);
		else
			printf("\nWriting CCSS file '%s' succeeded\n",outname);
		cc->del(cc);
		free(samples);

#ifdef DEBUG
		printf("About to exit\n");
#endif
		return 0;
	}

	/* CCMX: See if we're working from two files */
	if (!doccss && innames[0][0] != '\000') {
		int n;
		char *oname = NULL;			/* Observer name */
		ccmx *cc;

		if (innames[1][0] == '\000') {
			error("Need two -f .ti3 files to create CCMX");
		}

		/* Open up each CIE file in turn, target then measured, */
		/* and read in the CIE values. */
		for (n = 0; n < 2; n++) {
			cgats *cgf = NULL;			/* cgats file data */
			int isLab = 0;				/* 0 if file CIE is XYZ, 1 if is Lab */
			double wxyz[3], scale = 1.0;/* Scale factor back to absolute */
			int sidx;					/* Sample ID index */
			int xix, yix, zix;
			ary3 *current = NULL;		/* Current value array */
			int ii, ti;
			int instspec = 0;			/* File is spectrale */

			/* Open CIE target values */
			cgf = new_cgats();			/* Create a CGATS structure */
			cgf->add_other(cgf, ""); 	/* Allow any signature file */
		
			if (cgf->read_name(cgf, innames[n]))
				error("CGATS file '%s' read error : %s",innames[n],cgf->err);
		
			if (cgf->ntables < 1)
				error ("Input file '%s' doesn't contain at least one table",innames[n]);
		
			/* Check if the file is suitable */
			if (cgf->find_field(cgf, 0, "LAB_L") < 0
			 && cgf->find_field(cgf, 0, "XYZ_X") < 0) {
		
				error ("No CIE data found in file '%s'",innames[n]);
			}
		
			if (cgf->find_field(cgf, 0, "LAB_L") >= 0)
				isLab = 1;
			
			if (cols == NULL) {
				if ((npat = cgf->t[0].nsets) <= 0)
					error("No sets of data in file '%s'",innames[n]);

				if ((refs = (ary3 *)malloc(npat * sizeof(ary3))) == NULL)
					error("malloc failed");
				if ((cols = (ary3 *)malloc(npat * sizeof(ary3))) == NULL)
					error("malloc failed");

			} else {
				if (npat != cgf->t[0].nsets)
					error ("Number of sets %d in file '%s' doesn't match other file %d",cgf->t[0].nsets,innames[n],npat);
			}

			if ((ii = cgf->find_kword(cgf, 0, "TARGET_INSTRUMENT")) < 0)
				error ("Can't find keyword TARGET_INSTRUMENT in '%s'",innames[n]);
	
			if ((ti = cgf->find_kword(cgf, 0, "INSTRUMENT_TYPE_SPECTRAL")) < 0)
				error ("Can't find keyword INSTRUMENT_TYPE_SPECTRAL in '%s'",innames[n]);
	
			if (strcmp(cgf->t[0].kdata[ti],"YES") == 0) {
				instspec = 1;		/* Currently is a spectral file */
				if (gotref)
					error("Found two spectral files - expect one colorimtric file");
				current = refs;
				refname = strdup(cgf->t[0].kdata[ii]);
				reffile = innames[n];
				gotref = 1;
				
			} else if (strcmp(cgf->t[0].kdata[ti],"NO") == 0) {
				instspec = 0;			/* Currently is not spectral file */
				if (gotcol) {
					/* Copy what we though was cols to refs */
					refname = colname;
					reffile = colfile;
					for (i = 0; i < npat; i++) {
						refs[i][0] = cols[i][0];
						refs[i][1] = cols[i][1];
						refs[i][2] = cols[i][2];
					}
					gotref = 1;
					warning("Got two colorimetric files - assuming '%s' is the refrence",innames[0]);
					refrmode = -1;
					cbid = 0;

					if (spec) {
						error("Spectral reference is required to use non-standard observer");
					}
				}
				if ((ti = cgf->find_kword(cgf, 0, "DISPLAY_TYPE_REFRESH")) >= 0) {
					if (stricmp(cgf->t[0].kdata[ti], "YES") == 0)
						refrmode = 1;
					else if (stricmp(cgf->t[0].kdata[ti], "NO") == 0)
						refrmode = 0;
				}
				if ((ti = cgf->find_kword(cgf, 0, "DISPLAY_TYPE_BASE_ID")) >= 0) {
					cbid = atoi(cgf->t[0].kdata[ti]);
				} else {
					cbid = 0;
				}
				current = cols;
				colname = strdup(cgf->t[0].kdata[ii]);
				colfile = innames[n];
				gotcol = 1;
			} else {
				error ("Unknown INSTRUMENT_TYPE_SPECTRAL value '%s'",cgf->t[0].kdata[ti]);
			}

			if ((ti = cgf->find_kword(cgf, 0, "NORMALIZED_TO_Y_100")) < 0
			 || strcmp(cgf->t[0].kdata[ti],"NO") == 0) {
				scale = 1.0;		/* Leave absolute */
			} else {
				if ((ti = cgf->find_kword(cgf, 0, "LUMINANCE_XYZ_CDM2")) < 0)
					error ("Can't find keyword LUMINANCE_XYZ_CDM2 in '%s'",innames[n]);
				if (sscanf(cgf->t[0].kdata[ti],"%lf %lf %lf", &wxyz[0], &wxyz[1], &wxyz[2]) != 3)
					error ("Unable to parse LUMINANCE_XYZ_CDM2 in '%s'",innames[n]);
				scale = wxyz[1]/100.0;	/* Convert from Y = 100 normalise back to absolute */
			}

			if (instspec && spec) {
				int ii;
				xspect sp;
				char buf[100];
				int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
				xsp2cie *sp2cie;			/* Spectral conversion object */

				if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_BANDS")) < 0)
					error ("Input file '%s' doesn't contain keyword SPECTRAL_BANDS",innames[n]);
				sp.spec_n = atoi(cgf->t[0].kdata[ii]);
				if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_START_NM")) < 0)
					error ("Input file '%s' doesn't contain keyword SPECTRAL_START_NM",innames[n]);
				sp.spec_wl_short = atof(cgf->t[0].kdata[ii]);
				if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_END_NM")) < 0)
					error ("Input file '%s' doesn't contain keyword SPECTRAL_END_NM",innames[n]);
				sp.spec_wl_long = atof(cgf->t[0].kdata[ii]);
				sp.norm = 1.0;		/* We assume emssive */

				/* Find the fields for spectral values */
				for (j = 0; j < sp.spec_n; j++) {
					int nm;
			
					/* Compute nearest integer wavelength */
					nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
					            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
					
					sprintf(buf,"SPEC_%03d",nm);

					if ((spi[j] = cgf->find_field(cgf, 0, buf)) < 0)
						error("Input file '%s' doesn't contain field %s",innames[n],buf);
				}

				/* Create a spectral conversion object */
				if ((sp2cie = new_xsp2cie(icxIT_none, 0.0, NULL, obType, custObserver, icSigXYZData, icxClamp)) == NULL)
					error("Creation of spectral conversion object failed");

				for (i = 0; i < npat; i++) {
	
					/* Read the spectral values for this patch */
					for (j = 0; j < sp.spec_n; j++) {
						sp.spec[j] = *((double *)cgf->t[0].fdata[i][spi[j]]) * scale;
					}
	
					/* Convert it to CIE space */
					sp2cie->convert(sp2cie, current[i], &sp);
				}
				sp2cie->del(sp2cie);        /* Done with this */

			/* Colorimetric file - assume it's the target */
			} else {

				if (isLab) {		/* Expect Lab */
					if ((xix = cgf->find_field(cgf, 0, "LAB_L")) < 0)
						error("Input file '%s' doesn't contain field LAB_L",innames[n]);
					if (cgf->t[0].ftype[xix] != r_t)
						error("Field LAB_L is wrong type - expect float");
					if ((yix = cgf->find_field(cgf, 0, "LAB_A")) < 0)
						error("Input file '%s' doesn't contain field LAB_A",innames[n]);
					if (cgf->t[0].ftype[yix] != r_t)
						error("Field LAB_A is wrong type - expect float");
					if ((zix = cgf->find_field(cgf, 0, "LAB_B")) < 0)
						error("Input file '%s' doesn't contain field LAB_B",innames[n]);
					if (cgf->t[0].ftype[zix] != r_t)
						error("Field LAB_B is wrong type - expect float");

				} else { 		/* Expect XYZ */
					if ((xix = cgf->find_field(cgf, 0, "XYZ_X")) < 0)
						error("Input file '%s' doesn't contain field XYZ_X",innames[n]);
					if (cgf->t[0].ftype[xix] != r_t)
						error("Field XYZ_X is wrong type - expect float");
					if ((yix = cgf->find_field(cgf, 0, "XYZ_Y")) < 0)
						error("Input file '%s' doesn't contain field XYZ_Y",innames[n]);
					if (cgf->t[0].ftype[yix] != r_t)
						error("Field XYZ_Y is wrong type - expect float");
					if ((zix = cgf->find_field(cgf, 0, "XYZ_Z")) < 0)
						error("Input file '%s' doesn't contain field XYZ_Z",innames[n]);
					if (cgf->t[0].ftype[zix] != r_t)
						error("Field XYZ_Z is wrong type - expect float");
				}

				for (i = 0; i < npat; i++) {
					current[i][0] = *((double *)cgf->t[0].fdata[i][xix]);
					current[i][1] = *((double *)cgf->t[0].fdata[i][yix]);
					current[i][2] = *((double *)cgf->t[0].fdata[i][zix]);
					if (isLab) {	/* Convert test patch Lab to XYZ scale 100 */
						icmLab2XYZ(&icmD50_100, current[i], current[i]);
					}
					/* Rescale to absolute if needed */
					current[i][0] *= scale;
					current[i][1] *= scale;
					current[i][2] *= scale;
				}
			}
			cgf->del(cgf);		/* Clean up */
			cgf = NULL;
		}

		if (spec != 0 && obType != icxOT_CIE_1931_2)
			oname = standardObserverDescription(obType);

		if (oname != NULL) {
			char *tt = colname;

			if ((colname = malloc(strlen(tt) + strlen(oname) + 3)) == NULL)
				error("Malloc failed");
			strcpy(colname, tt);
			strcat(colname, " (");
			strcat(colname, oname);
			strcat(colname, ")");
		}
		if (description == NULL) {
			char *disp = displayname != NULL ? displayname : dtinfo->desc;
			if ((description = malloc(strlen(colname) + strlen(disp) + 4)) == NULL)
				error("Malloc failed");
			strcpy(description, colname);
			strcat(description, " & ");
			strcat(description, disp);
		}

		if (refrmode < 0)
			error("The display refresh mode is not in '%s' - use the -Y flag",colfile);

		if (cbid == 0)
			error("The calibration base display mode not specified in the '%s' file",colfile);

		if ((cc = new_ccmx()) == NULL)
			error("new_ccmx() failed");

		if (cc->create_ccmx(cc, description, colname, displayname, dtinfo->dtech,
			                     refrmode, cbid, uisel, refname, 0, npat, refs, cols)) {
			error("create_ccmx failed with '%s'\n",cc->err);
		}
		if (verb) {
			printf("Fit error is max %f, avg %f DE94\n",cc->mx_err,cc->av_err);
			printf("Correction matrix is:\n");
			printf("  %f %f %f\n", cc->matrix[0][0], cc->matrix[0][1], cc->matrix[0][2]);
			printf("  %f %f %f\n", cc->matrix[1][0], cc->matrix[1][1], cc->matrix[1][2]);
			printf("  %f %f %f\n", cc->matrix[2][0], cc->matrix[2][1], cc->matrix[2][2]);
		}

		if(cc->write_ccmx(cc, outname))
			printf("\nWriting CCMX file '%s' failed with '%s'\n",outname,cc->err);
		else
			printf("\nWriting CCMX file '%s' succeeded\n",outname);
		cc->del(cc);

	/* Do interactive measurements */
	} else {

		/* No explicit display has been set */
		if (
#ifndef SHOW_WINDOW_ONFAKE
		!fake 
#endif
#ifdef NT
		 && madvrdisp == 0
#endif
		 && dummydisp == 0
		 && webdisp == 0
		 && ccdisp == 0
		 && disp == NULL) {
			int ix = 0;
#if defined(UNIX_X11)
			char *dn, *pp;

			if ((dn = getenv("DISPLAY")) != NULL) {
				if ((pp = strrchr(dn, ':')) != NULL) {
					if ((pp = strchr(pp, '.')) != NULL) {
						if (pp[1] != '\000')
							ix = atoi(pp+1);
					}
				}
			}
#endif
			if ((disp = get_a_display(ix)) == NULL)
				error("Unable to open the default display");

			if (displayname == NULL && (displayname = strdup(disp->description)) == NULL)
				error("Malloc failed");

			printf("Display description is '%s'\n",displayname);
		}
		if (fake) {
			displayname = strdup("fake display");
		}

		/* If we've requested ChromeCast, look it up */
		if (ccdisp) {
			if ((ccids = get_ccids()) == NULL)
				error("discovering ChromCasts failed");
			if (ccids[0] == NULL)
				error("There are no ChromCasts to use\n");
			for (i = 0; ccids[i] != NULL; i++)
				;
			if (ccdisp < 1 || ccdisp > i)
				error("Chosen ChromCasts (%d) is outside list (1..%d)\n",ccdisp,i);
			ccid = ccids[ccdisp-1];
		}

		/* Create grid of device test values */
		{
			int j;
			int gc[3];			/* Grid coordinate */

			if (msteps == 1)
				npat = 4;
			else
				npat = msteps * msteps * msteps;

			if ((rdcols = (col *)malloc(npat * sizeof(col))) == NULL) {
				error("malloc failed");
			}
			if ((refs = (ary3 *)malloc(npat * sizeof(ary3))) == NULL) {
				free(rdcols);
				error("malloc failed");
			}
			if ((cols = (ary3 *)malloc(npat * sizeof(ary3))) == NULL) {
				free(rdcols);
				free(refs);
				error("malloc failed");
			}

			/* RGBW */
			if (msteps == 1) {
				npat = 0;
				rdcols[npat].r = 1.0;
				rdcols[npat].g = 0.0;
				rdcols[npat].b = 0.0;
				npat++;
				rdcols[npat].r = 0.0;
				rdcols[npat].g = 1.0;
				rdcols[npat].b = 0.0;
				npat++;
				rdcols[npat].r = 0.0;
				rdcols[npat].g = 0.0;
				rdcols[npat].b = 1.0;
				npat++;
				rdcols[npat].r = 1.0;
				rdcols[npat].g = 1.0;
				rdcols[npat].b = 1.0;
				npat++;
#ifdef DEBUG
				for (j = 0; j < 4; j++) 
					printf("Dev val %f %f %f\n",rdcols[j].r,rdcols[j].g,rdcols[j].b);
#endif
			} else {
				for (j = 0; j < 3; j++)
					gc[j] = 0;			/* init coords */
					
				for (npat = 0; ;) {	/* For all grid points */

					/* Just colors with at least one channel at 100% */
					if (gc[0] == (msteps-1)
					 || gc[1] == (msteps-1)
					 || gc[2] == (msteps-1)) 
					{
						
						rdcols[npat].r = (double)gc[0]/(msteps-1);
						rdcols[npat].g = (double)gc[1]/(msteps-1);
						rdcols[npat].b = (double)gc[2]/(msteps-1);
#ifdef DEBUG
						printf("Dev val %f %f %f\n",rdcols[npat].r,rdcols[npat].g,rdcols[npat].b);
#endif
						npat++;
					}
					
					/* Increment grid index and position */
					for (j = 0; j < 3; j++) {
						gc[j]++;
						if (gc[j] < msteps)
							break;	/* No carry */
						gc[j] = 0;
					}
					if (j >= 3)
						break;		/* Done grid */
				}
			}
			if (verb)
				printf("Total test patches = %d\n",npat);
		}

		/* Until the measurements are done, or we give up */
		for (;;) {
			int c;

			/* Print the menue of adjustments */
			printf("\n");
			if (gotref)
				printf("[Got spectrometer readings]\n");
			if (gotcol)
				printf("[Got colorimeter readings]\n");
			printf("Press 1 .. 4:\n");
			{
				printf("1) Select an instrument, Currently %d (", comno);
				if (icmps == NULL)
					icmps = new_icompaths(g_log);
				else
					icmps->refresh(icmps);
				if (icmps != NULL) {
					icompath **paths;
					if ((paths = icmps->paths) != NULL) {
						int i;
						for (i = 0; ; i++) {
							if (paths[i] == NULL)
								break;
							if ((i+1) == comno) {
								printf(" '%s'",paths[i]->name);
								break;
							}
						}
					}
				}
				printf(")\n");
			}
			if (doccss)
				printf("2) Measure test patches with current (spectrometer) instrument\n");
			else
				printf("2) Measure test patches with current instrument\n");

			if (doccss) {
				if (gotref)
					printf("3) Save Colorimeter Calibration Spectral Set\n");
				else
					printf("3) [ Save Colorimeter Calibration Spectral Set ]\n");

			} else {
				if (gotref && gotcol)
					printf("3) Compute Colorimeter Correction Matrix & save it\n");
				else
					printf("3) [ Compute Colorimeter Correction Matrix & save it ]\n");
			}
			printf("4) Exit\n");

			if (fakeseq == 0) {
				empty_con_chars();
				c = next_con_char();
			} else {
				switch (fakeseq) {
					case 1:
						c = '2';
						fakeseq = 2;
						break;
					case 2:
						c = '2';
						fakeseq = 3;
						break;
					case 3:
						c = '3';
						fakeseq = 4;
						break;
					default:
						c = '4';
						break;
				}
			}
			printf("'%c'\n",c);


			/* Deal with selecting the instrument */
			if (c == '1') {
				if (icmps == NULL)
					icmps = new_icompaths(g_log);
				else
					icmps->refresh(icmps);
				if (icmps != NULL) {
					icompath **paths;
					if ((paths = icmps->paths) != NULL) {
						int i;
						for (i = 0; ; i++) {
							if (paths[i] == NULL)
								break;
							if ((paths[i]->dtype == instSpyder1 && setup_spyd2(0) == 0)
							 || (paths[i]->dtype == instSpyder2 && setup_spyd2(1) == 0))
								fprintf(stderr,"    %d = '%s' !! Disabled - no firmware !!\n",i+1,paths[i]->name);
							else
								fprintf(stderr,"    %d = '%s'\n",i+1,paths[i]->name);
						}
						printf("Select device 1 - %d: \n",i);
						empty_con_chars();
						c = next_con_char();

						if (c < '1' || c > ('0' + i)) {
							printf("'%c' is out of range - ignored !\n",c);
						} else {
							comno = c - '0'; 
						}
						
					} else {
						fprintf(stderr,"No ports to select from!\n");
					}
				}
				continue;
			}

			/* Deal with doing a measurement */
			if (c == '2') {
				int errc;				/* Return value from new_disprd() */
				disprd *dr;				/* Display patch read object */
				inst *it;				/* Instrument */
				inst_mode  cap = inst_mode_none;	/* Instrument mode capabilities */
				inst2_capability cap2 = inst2_none;	/* Instrument capabilities 2 */
				inst3_capability cap3 = inst3_none;	/* Instrument capabilities 3 */

				if (fake)
					comno = FAKE_DEVICE_PORT;
				if (icmps == NULL)
					icmps = new_icompaths(g_log);

				/* Should we use current cal rather than native ??? */
				if ((dr = new_disprd(&errc, icmps->get_path(icmps, comno),
				                     fc, ditype, sditype, 1, tele, ambient, nadaptive,
				                     noinitcal, 0, highres, refrate, 3, NULL, NULL,
					                 NULL, 0, disp, 0, fullscreen,
				                     override, webdisp, ccid,
#ifdef NT
					                 madvrdisp,
#endif
									 dummydisp, ccallout, mcallout, 0,
					                 100.0 * hpatscale, 100.0 * vpatscale, ho, vo,
					                 disptech_unknown, 0, NULL, NULL, 0, 2, icxOT_default, NULL, 
				                     0, 0, "fake" ICC_FILE_EXT, g_log)) == NULL) {
					printf("Opening the instrument failed with '%s'. Try selecting it again ?\n",
					        disprd_err(errc));
					continue;
				}
			
				it = dr->it;

				if (fake) {
					if (faketoggle)
						cap = inst_mode_spectral;
					else
						cap = inst_mode_colorimeter;
					cap2 = inst2_none;
					cap3 = inst3_none;
					refrmode = 0; 
					cbid = 1;
				} else {
					it->capabilities(it, &cap, &cap2, &cap3);
					if (!IMODETST(cap, inst_mode_spectral)) {
						dr->get_disptype(dr, &refrmode, &cbid);	/* Get the display type info */
					}
				}

				if (doccss && !IMODETST(cap, inst_mode_spectral)) {
					printf("You have to use a spectrometer to create a CCSS!\n");
					continue;
				}

				if (faketoggle)
					dr->fake2 = 1;
				else
					dr->fake2 = -1;

				/* Test the CRT with all of the test points */
				if ((rv = dr->read(dr, rdcols, npat, 1, npat, 1, 0, instClamp)) != 0) {
					dr->del(dr);
					error("disprd returned error code %d\n",rv);
				}

				if (doccss) {		/* We'll use the rdcols values */
					gotref = 1;
				} else {
					if (IMODETST(cap, inst_mode_spectral)) {
						xsp2cie *sp2cie = NULL;

						if (spec) {
							/* Create a spectral conversion object */
							if ((sp2cie = new_xsp2cie(icxIT_none, 0.0, NULL, obType, custObserver, icSigXYZData, icxClamp)) == NULL)
								error("Creation of spectral conversion object failed");
						}
						for (i = 0; i < npat; i++) {	/* For all grid points */
							if (spec) {
								if (rdcols[i].sp.spec_n <= 0)
									error("Didn't get spectral value");
								sp2cie->convert(sp2cie, refs[i], &rdcols[i].sp);
							} else {
								if (rdcols[i].XYZ_v == 0)
									error("Didn't get XYZ value");
								refs[i][0] = rdcols[i].XYZ[0];
								refs[i][1] = rdcols[i].XYZ[1];
								refs[i][2] = rdcols[i].XYZ[2];
							}
						}
						if (fake)
							refname = "fake spectrometer";
						else
							refname = inst_name(it->dtype);
						gotref = 1;
						if (sp2cie != NULL)
							sp2cie->del(sp2cie);
					} else if (IMODETST(cap, inst_mode_colorimeter)) {
						for (i = 0; i < npat; i++) {	/* For all grid points */
							if (rdcols[i].XYZ_v == 0)
								error("Didn't get XYZ value");
							cols[i][0] = rdcols[i].XYZ[0];
							cols[i][1] = rdcols[i].XYZ[1];
							cols[i][2] = rdcols[i].XYZ[2];
						}
						if (fake)
							colname = "fake colorimeter";
						else
							colname = inst_name(it->dtype);
						gotcol = 1;
					}
				}
				dr->del(dr);

				faketoggle ^= 1;

			}	/* End of take a measurement */

			if (c == '3') {		/* Compute result and save */
				/* Save the CCSS */
				if (doccss) {
					ccss *cc;
					xspect *samples = NULL;
					double bigv = -1e60;

					if (!gotref) {
						printf("You have to read the spectrometer values first!\n");
						continue;
					}

					if (description == NULL) {
						char *disp = displayname != NULL ? displayname : dtinfo->desc;
						char *tt = "CCSS for ";
						if ((description = malloc(strlen(disp) + strlen(tt) + 1)) == NULL)
							error("Malloc failed");
						strcpy(description, tt);
						strcat(description, disp);
					}

					if ((samples = (xspect *)malloc(sizeof(xspect) * npat)) == NULL)
						error("Malloc failed");
					
					/* See what the highest value is */
					for (i = 0; i < npat; i++) {	/* For all grid points */
						if (rdcols[i].sp.spec_n <= 0)
							error("Didn't get spectral values");
						for (j = 0; j < rdcols[i].sp.spec_n; j++) {
							if (rdcols[i].sp.spec[j] > bigv)
								bigv = rdcols[i].sp.spec[j];
						}
					}
						
					/* Copy all the values and normalize them */
					for (i = 0; i < npat; i++) {	/* For all grid points */
						double scale = 100.0;

						samples[i] = rdcols[i].sp;	/* Structure copy */
						for (j = 0; j < rdcols[i].sp.spec_n; j++)
							samples[i].spec[j] *= scale / bigv;
					}
						
					if (refrmode < 0)
						warning("No refresh mode specified! Assuming non-refresh !");

					if ((cc = new_ccss()) == NULL)
						error("new_ccss() failed");
	
					if (cc->set_ccss(cc, "Argyll ccxxmake", NULL, description, displayname,
					                 dtinfo->dtech, refrmode, NULL, refname, 0, samples, npat)) {
						error("set_ccss failed with '%s'\n",cc->err);
					}
					if(cc->write_ccss(cc, outname))
						printf("\nWriting CCSS file '%s' failed with '%s'\n",outname,cc->err);
					else
						printf("\nWriting CCSS file '%s' succeeded\n",outname);
					cc->del(cc);
					free(samples);
					saved = 1;

				/* Compute and save CCMX */
				} else {
					char *oname = NULL;		/* Observer desciption */
					ccmx *cc;
	
					if (!gotref) {
						printf("You have to read the spectrometer values first!\n");
						continue;
					}
					if (!gotcol) {
						printf("You have to read the colorimeter values first!\n");
						continue;
					}
	
					if (spec != 0 && obType != icxOT_CIE_1931_2)
						oname = standardObserverDescription(obType);
	
					if (oname != NULL) {		/* Incorporate observer name in colname */
						char *tt = colname;
	
						if ((colname = malloc(strlen(tt) + strlen(oname) + 3)) == NULL)
							error("Malloc failed");
						strcpy(colname, tt);
						strcat(colname, " (");
						strcat(colname, oname);
						strcat(colname, ")");
					}
					if (description == NULL) {
						char *disp = displayname != NULL ? displayname : dtinfo->desc;
						if ((description = malloc(strlen(colname) + strlen(disp) + 4)) == NULL)
							error("Malloc failed");
						strcpy(description, colname);
						strcat(description, " & ");
						strcat(description, disp);
					}

					if (refrmode < 0)
						error("Internal error - the instrument did not return a refmode");
					if (cbid == 0)
						error("Internal error - the instrument did not return a cbid");

					if ((cc = new_ccmx()) == NULL)
						error("new_ccmx() failed");
	
					if (cc->create_ccmx(cc, description, colname, displayname, dtinfo->dtech,
					                    refrmode, cbid, uisel, refname, 0, npat, refs, cols)) {
						error("create_ccmx failed with '%s'\n",cc->err);
					}
					if (verb) {
						printf("Fit error is avg %f, max %f DE94\n",cc->av_err,cc->mx_err);
						printf("Correction matrix is:\n");
						printf("  %f %f %f\n", cc->matrix[0][0], cc->matrix[0][1], cc->matrix[0][2]);
						printf("  %f %f %f\n", cc->matrix[1][0], cc->matrix[1][1], cc->matrix[1][2]);
						printf("  %f %f %f\n", cc->matrix[2][0], cc->matrix[2][1], cc->matrix[2][2]);
					}
	
					if(cc->write_ccmx(cc, outname))
						printf("\nWriting CCMX file '%s' failed with '%s'\n",outname,cc->err);
					else
						printf("\nWriting CCMX file '%s' succeeded\n",outname);
					cc->del(cc);
					saved = 1;
				}
			}

			if (c == '4' || c == 0x3) {		/* Exit */
				if (!saved) {
					printf("Not saved yet, are you sure ? (y/n): "); fflush(stdout);
					empty_con_chars();
					c = next_con_char();
					printf("\n");
					if (c != 'y' && c != 'Y')
						continue;
				}
				break;
			}

		}	/* Next command */

		free(displayname);
		if (icmps != NULL)
			icmps->del(icmps);
		free_ccids(ccids);
	}

#ifdef DEBUG
	/* Do a CCMX verification */
	if (!doccss) {
		ccmx *cc;
		double av_err, mx_err;
		int wix;
		double maxy = -1e6;
		icmXYZNumber wh;

		for (i = 0; i < npat; i++) {
			if (refs[i][1] > maxy) {
				maxy = refs[i][1];
				wix = i;
			}
		}
		wh.X = refs[wix][0];
		wh.Y = refs[wix][1];
		wh.Z = refs[wix][2];

		if ((cc = new_ccmx()) == NULL)
			error("new_ccmx() failed");
		if(cc->read_ccmx(cc, outname))
			printf("Reading file '%s' failed\n",outname);

		av_err = mx_err = 0.0;
		for (i = 0; i < npat; i++) {
			double txyz[3], tlab[3], xyz[3], _xyz[3], lab[3], de;
			icmCpy3(txyz, refs[i]);
			icmXYZ2Lab(&wh, tlab, txyz);
			icmCpy3(xyz, cols[i]);
			cc->xform(cc,_xyz, xyz);
			icmXYZ2Lab(&wh, lab, _xyz);
			de = icmCIE94(tlab, lab);
			av_err += de;
			if (de > mx_err)
				mx_err = de;
			printf("%d: txyz %f %f %f\n",i,txyz[0], txyz[1], txyz[2]);
			printf("%d:  xyz %f %f %f, _xyz %f %f %f\n",i,xyz[0], xyz[1], xyz[2], _xyz[0], _xyz[1], _xyz[2]);
			printf("%d: tlab %f %f %f, lab %f %f %f\n",i,tlab[0], tlab[1], tlab[2], lab[0], lab[1], lab[2]);
			printf("%d: de %f\n",i,de);
		}
		av_err /= npat;
		printf("Avg = %f, max = %f\n",av_err,mx_err);
		cc->del(cc);
	}
#endif

#ifdef DEBUG
	printf("About to exit\n");
#endif

	return 0;
}




