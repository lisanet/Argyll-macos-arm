/* 
 * Argyll Color Management System
 * DTP92/Spectrolino display target reader
 *
 * Author: Graeme W. Gill
 * Date:   4/10/96
 *
 * Copyright 1996 - 2013 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This program displays test patches, and takes readings from a display device */

/* TTBD

	Add support for black recalibration using i1pro or munki.
	Setable timeout ? Need to allow placing instrument back
	on screen. Need this to properly handle ss anyway ?

    Add bell at end of readings ?

	Ideally this should be changed to always create non-normalized (absolute)
	readings, since normalised readings are not natural, but everything
	that deals with .ti3 data then needs fixing to deal with non-normalized
	readings !
 */

#undef DEBUG
#undef DEBUG_OFFSET	/* Keep test window out of the way */

/* Invoke with -dfake for testing with a fake device */
/* Will use fake.icm/.icc if present */

#define COMPORT 1		/* Default com port 1..4 */

#ifdef __MINGW32__
# define WINVER 0x0500
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#if defined (NT)
#include <conio.h>
#endif
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "cgats.h"
#include "xspect.h"
#include "insttypes.h"
#include "conv.h"
#include "icoms.h"
#include "inst.h"
#include "ccmx.h"
#include "ccss.h"
#include "ccast.h"
#include "dispwin.h"
#include "dispsup.h"
#include "sort.h"
#include "instappsup.h"
#ifdef ENABLE_USB
# include "spyd2.h"
#endif
#include "ui.h"

/* ------------------------------------------------------------------- */
#if defined(__APPLE__) && defined(__POWERPC__)

/* Workaround for a PPC gcc 3.3 optimiser bug... */
/* It seems to cause a segmentation fault instead of */
/* converting an integer loop index into a float, */
/* when there are sufficient variables in play. */
static int gcc_bug_fix(int i) {
	static int nn;
	nn += i;
	return nn;
}
#endif	/* APPLE */

/* ------------------------------------------------------------------- */

/*

  Flags used:

         ABCDEFGHIJKLMNOPQRSTUVWXYZ
  upper    .... .... .. ..    .....
  lower    ..      .  . .  .  .. . 

*/

/* Flag = 0x0000 = default */
/* Flag & 0x0001 = list ChromeCast's */
void usage(int flag, char *diag, ...) {
	int i;
	disppath **dp;
	icompaths *icmps;
	inst2_capability cap2 = inst2_none;

	fprintf(stderr,"Read a Display, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: dispread [options] outfile\n");
	fprintf(stderr," -v              Verbose mode\n");
#if defined(UNIX_X11)
	fprintf(stderr," -display displayname Choose X11 display name\n");
	fprintf(stderr," -d n[,m]             Choose the display n from the following list (default 1)\n");
	fprintf(stderr,"                      Optionally choose different display m for VideoLUT access\n"); 
#else
	fprintf(stderr," -d n                 Choose the display from the following list (default 1)\n");
#endif
	dp = get_displays();
	if (dp == NULL || dp[0] == NULL)
		fprintf(stderr,"    ** No displays found **\n");
	else {
		int i;
		for (i = 0; ; i++) {
			if (dp[i] == NULL)
				break;
			fprintf(stderr,"    %d = '%s'\n",i+1,dp[i]->description);
		}
	}
	free_disppaths(dp);
	fprintf(stderr," -dweb[:port]         Display via a web server at port (default 8080)\n");
	fprintf(stderr," -dcc[:n]             Display via n'th ChromeCast (default 1, ? for list)\n");
	if (flag & 0x001) {
		ccast_id **ids;
		if ((ids = get_ccids()) == NULL) {
			fprintf(stderr,"    ** Error discovering ChromeCasts **\n");
		} else {
			if (ids[0] == NULL)
				fprintf(stderr,"    ** No ChromeCasts found **\n");
			else {
				int i;
				for (i = 0; ids[i] != NULL; i++)
					fprintf(stderr,"    %d = '%s'\n",i+1,ids[i]->name);
				free_ccids(ids);
			}
		}
	}
#ifdef NT
	fprintf(stderr," -d madvr             Display via MadVR Video Renderer\n");
#endif
	fprintf(stderr," -d dummy             Dummy (non-existant, invisible) display\n");
//	fprintf(stderr," -d fake              Use a fake (ICC Profile) display device for testing, fake%s if present\n",ICC_FILE_EXT);
	fprintf(stderr," -c listno            Set communication port from the following list (default %d)\n",COMPORT);
	if ((icmps = new_icompaths(g_log)) != NULL) {
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
		} else
			fprintf(stderr,"    ** No ports found **\n");
	}
	fprintf(stderr," -p                   Use telephoto mode (ie. for a projector, if available)\n");
	fprintf(stderr," -a                   Use ambient measurement mode (ie. for a projector, if available)\n");
	cap2 = inst_show_disptype_options(stderr, " -y                   ", icmps, 0, 0);
	fprintf(stderr," -k file.cal          Load calibration file into display while reading\n");
	fprintf(stderr," -K file.cal          Apply calibration file to test values while reading\n");
#ifdef NT
	fprintf(stderr," -V                   Enable MadVR color management (3dLut)\n");
#endif
	fprintf(stderr," -s                   Save spectral information (default don't save)\n");
	fprintf(stderr," -P ho,vo,ss[,vs]     Position test window and scale it\n");
	fprintf(stderr,"                      ho,vi: 0.0 = left/top, 0.5 = center, 1.0 = right/bottom etc.\n");
	fprintf(stderr,"                      ss: 0.5 = half, 1.0 = normal, 2.0 = double etc.\n");
	fprintf(stderr," -F                   Fill whole screen with black background\n");
#if defined(UNIX_X11)
	fprintf(stderr," -n                   Don't set override redirect on test window\n");
#endif
	fprintf(stderr," -E                   Encode the test values for video range 16..235/255\n");
	fprintf(stderr," -Z nbits             Quantize test values to fit in nbits\n");
	fprintf(stderr," -J                   Run instrument calibration first (used rarely)\n");
	fprintf(stderr," -N                   Disable initial calibration of instrument if possible\n");
	fprintf(stderr," -H                   Use high resolution spectrum mode (if available)\n");
//	fprintf(stderr," -V                   Use adaptive measurement mode (if available)\n");
	fprintf(stderr," -w                   Disable normalisation of white to Y = 100\n");
	if (cap2 & inst2_ccmx)
		fprintf(stderr," -X file.ccmx         Apply Colorimeter Correction Matrix\n");
	if (cap2 & inst2_ccss) {
		fprintf(stderr," -X file.ccss         Use Colorimeter Calibration Spectral Samples for calibration\n");
		fprintf(stderr," -Q observ            Choose CIE Observer for spectrometer or CCSS colorimeter data:\n");
		fprintf(stderr,"                      1931_2 (def), 1964_10, 2012_2, 2012_10, S&B 1955_2, shaw, J&V 1978_2, 1964_10c or file.cmf\n");
	}
	fprintf(stderr," -I b|w               Drift compensation, Black: -Ib, White: -Iw, Both: -Ibw\n");
	fprintf(stderr," -Y R:rate            Override measured refresh rate with rate Hz\n");
	fprintf(stderr," -Y A                 Use non-adaptive integration time mode (if available).\n");
	fprintf(stderr," -Y p                 Don't wait for the instrument to be placed on the display\n");
	fprintf(stderr," -Y k                 Omit the file.cal information from the .ti3 file\n"); 
	fprintf(stderr," -C \"command\"         Invoke shell \"command\" each time a color is set\n");
	fprintf(stderr," -M \"command\"         Invoke shell \"command\" each time a color is measured\n");
//	fprintf(stderr," -x [lx]              Take manually entered values, either L*a*b* (-xl) or XYZ (-xx).\n");
	fprintf(stderr," -x x                 Take manually entered XYZ values\n");
	fprintf(stderr," -W n|h|x             Override serial port flow control: n = none, h = HW, x = Xon/Xoff\n");
	fprintf(stderr," -D [level]           Print debug diagnostics to stderr\n");
	fprintf(stderr," outfile              Base name for input[ti1]/output[ti3] file\n");
	if (icmps != NULL)
		icmps->del(icmps);
	exit(1);
}

int main(int argc, char *argv[]) {
	int i, j;
	int fa, nfa, mfa;					/* current argument we're looking at */
	disppath *disp = NULL;				/* Display being used */
	double hpatscale = 1.0, vpatscale = 1.0;	/* scale factor for test patch size */
	double ho = 0.0, vo = 0.0;			/* Test window offsets, -1.0 to 1.0 */
	int out_tvenc = 0;					/* 1 to use RGB Video Level encoding */
	int qbits = 0;						/* Quantization bits, 0 = not set */
	int fullscreen = 0;            		/* NZ if whole screen should be filled with black */
	int verb = 0;
	int debug = 0;
	int fake = 0;						/* Use the fake device for testing */
	int override = 1;					/* Override redirect on X11 */
	int comport = COMPORT;				/* COM port used */
	icompaths *icmps = NULL;
	icompath *ipath = NULL;
	flow_control fc = fc_nc;			/* Default flow control */
	int docalib = 0;					/* Force a calibration */
	int highres = 0;					/* Use high res mode if available */
	double refrate = 0.0;			    /* 0.0 = default, > 0.0 = override refresh rate */ 
	int nadaptive = 0;					/* Use non-adaptive mode if available */
	int bdrift = 0;						/* Flag, nz for black drift compensation */
	int wdrift = 0;						/* Flag, nz for white drift compensation */
	int ditype = 0;						/* Display type selection charater(s) */
	int tele = 0;						/* NZ if telephoto mode */
	int ambient = 0;					/* NZ if ambient mode */
	int noautocal = 0;					/* Disable auto calibration */
	int noplace = 0;					/* Disable user instrument placement */
	int donorm = 1;						/* Enable Y = 100 normalisation */
	char ccxxname[MAXNAMEL+1] = "\000";  /* Colorimeter Correction Matrix name */
	ccmx *cmx = NULL;					/* Colorimeter Correction Matrix */
	ccss *ccs = NULL;					/* Colorimeter Calibration Spectral Samples */
	int spec = 0;						/* Don't save spectral information */
	icxObserverType obType = icxOT_default;
	xspect custObserver[3];				/* If obType = icxOT_custom */
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
	int xtern = 0;						/* Use external (user supplied) values rather than */
										/* instument read 1 = Lab, 2 = XYZ */ 
	char inname[MAXNAMEL+1] = "\000";	/* Input cgats file base name */
	char outname[MAXNAMEL+1] = "\000";	/* Output cgats file base name */
	char calname[MAXNAMEL+1] = "\000";	/* Calibration file name (if any) */
	int native = 2;						/* X0 = use current per channel calibration curve */
										/* X1 = set native linear output and use ramdac high prec */
										/* 0X = use current color management cLut (MadVR) */
										/* 1X = disable color management cLUT (MadVR) */
										/* Note that if calname, X0 loads calfile into hardware */
										/* and X1 applied soft calibration. */
	double cal[3][MAX_CAL_ENT];			/* Display calibration */
	int ncal = 256;						/* Default number of cal entries used */
	int nocaloutput = 0;				/* NZ to omit calibration info from .ti3 */
	cgats *icg;							/* input cgats structure */
	cgats *ocg;							/* output cgats structure */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp);			/* Ascii time */
	col *cols;							/* Internal storage of all the patch colors */
	int dim = 0;						/* Dimensionality - 1, 3, or 4 */
	int npat;							/* Number of patches/colors */
	int xpat = 0;						/* Set to number of extra patches */
	int wpat;							/* Set to index of white patch */
	int si;								/* Sample id index */
	int ti;								/* Temp index */
	int fi;								/* Colorspace index */
	int nsetel = 0;
	cgats_set_elem *setel;				/* Array of set value elements */
	disprd *dr;							/* Display patch read object */
	int noramdac = 0;					/* Will be set to nz if can't set ramdac */
	int nocm = 0;						/* Will be set to nz if can't set color management */
	int errc;							/* Return value from new_disprd() */
	int rv;

	set_exe_path(argv[0]);				/* Set global exe_path and error_program */
	check_if_not_interactive();

#ifdef DEBUG_OFFSET
	ho = 0.8;
	vo = -0.8;
#endif

#if defined(DEBUG) || defined(DEBUG_OFFSET)
	printf("!!!!!! Debug turned on !!!!!!\n");
#endif

	if (argc <= 1)
		usage(0,"Too few arguments");

	if (ncal > MAX_CAL_ENT)
		error("Internal, ncal = %d > MAX_CAL_ENT %d\n",ncal,MAX_CAL_ENT);

	/* Process the arguments */
	mfa = 1;        /* Minimum final arguments */
	for (fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {		/* Look for any flags */
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
				usage(0,"Usage requested");

			} else if (argv[fa][1] == 'v') {
				verb = 1;
				g_log->verb = verb;

			/* Display number */
			} else if (argv[fa][1] == 'd') {
				if (strncmp(na,"web",3) == 0
				 || strncmp(na,"WEB",3) == 0) {
					webdisp = 8080;
					if (na[3] == ':') {
						webdisp = atoi(na+4);
						if (webdisp == 0 || webdisp > 65535)
							usage(0,"Web port number must be in range 1..65535");
					}
					fa = nfa;
				} else if (strncmp(na,"cc",2) == 0
				 || strncmp(na,"CC",2) == 0) {
					ccdisp = 1;
					if (na[2] == ':') {
						if (na[3] < '0' || na[3] > '9')
							usage(0x0001,"Available ChromeCasts");

						ccdisp = atoi(na+3);
						if (ccdisp <= 0)
							usage(0,"ChromeCast number must be in range 1..N");
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
						if (++fa >= argc || argv[fa][0] == '-') usage(0,"Parameter expected following -display");
						setenv("DISPLAY", argv[fa], 1);
					} else {
						if (na == NULL) usage(0,"Parameter expected following -d");
						fa = nfa;
						if (strcmp(na,"fake") == 0) {
							fake = 1;
						} else {
							if (sscanf(na, "%d,%d",&ix,&iv) != 2) {
								ix = atoi(na);
								iv = 0;
							}
							if (disp != NULL)
								free_a_disppath(disp);
							if ((disp = get_a_display(ix-1)) == NULL)
								usage(0,"-d parameter %d out of range",ix);
							if (iv > 0)
								disp->rscreen = iv-1;
						}
					}
#else
					int ix;
					if (na == NULL) usage(0,"Parameter expected following -d");
					fa = nfa;
					if (strcmp(na,"fake") == 0) {
						fake = 1;
					} else {
						ix = atoi(na);
						if (disp != NULL)
							free_a_disppath(disp);
						if ((disp = get_a_display(ix-1)) == NULL)
							usage(0,"-d parameter %d out of range",ix);
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
				if (na == NULL) usage(0,"Paramater expected following -c");
				comport = atoi(na);
				if (comport < 1 || comport > 50) usage(0,"-c parameter %d out of range",comport);

			/* Telephoto */
			} else if (argv[fa][1] == 'p') {
				tele = 1;
				ambient = 0;

			/* Request ambient measurement */
			} else if (argv[fa][1] == 'a') {
				ambient = 1;
				tele = 0;

			/* Display type */
			} else if (argv[fa][1] == 'y') {
				fa = nfa;
				if (na == NULL) usage(0,"Parameter expected after -y");
				ditype = na[0];
				if (ditype == '_' && na[1] != '\000')
					ditype = ditype << 8 | na[1];

			/* Calibration file */
			} else if (argv[fa][1] == 'k'
			        || argv[fa][1] == 'K') {
				if (na == NULL) usage(0,"Parameter expected after -%c",argv[fa][1]);
				strncpy(calname,na,MAXNAMEL); calname[MAXNAMEL] = '\000';
				if (argv[fa][1] == 'K')
					native |= 1;			/* Use native linear & soft cal */
				else
					native &= ~1;			/* Use HW cal */
				fa = nfa;

#ifdef NT
			/* MadVR verify mode */
			} else if (argv[fa][1] == 'V') {
				native &= ~2;
#endif

			/* Save spectral data */
			} else if (argv[fa][1] == 's') {
				spec = 1;

			/* Test patch offset and size */
			} else if (argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage(0,"Parameter expected after -P");
				if (sscanf(na, " %lf,%lf,%lf,%lf ", &ho, &vo, &hpatscale, &vpatscale) == 4) {
					;
				} else if (sscanf(na, " %lf,%lf,%lf ", &ho, &vo, &hpatscale) == 3) {
					vpatscale = hpatscale;
				} else {
					usage(0,"-P parameter '%s' not recognised",na);
				}
				if (ho < 0.0 || ho > 1.0
				 || vo < 0.0 || vo > 1.0
				 || hpatscale <= 0.0 || hpatscale > 50.0
				 || vpatscale <= 0.0 || vpatscale > 50.0)
					usage(0,"-P parameters %f %f %f %f out of range",ho,vo,hpatscale,vpatscale);
				ho = 2.0 * ho - 1.0;
				vo = 2.0 * vo - 1.0;

			/* Full screen black background */
			} else if (argv[fa][1] == 'F') {
				fullscreen = 1;

			/* Video encoded output */
			} else if (argv[fa][1] == 'E') {
				out_tvenc = 1;
				if (qbits == 0)
					qbits = 8;

			/* Specify quantization bits */
			} else if (argv[fa][1] == 'Z') {
				fa = nfa;
				if (na == NULL) usage(0,"Expected argument to -Z");
				qbits = atoi(na);
				if (qbits < 1 || qbits > 32)
					usage(0,"Argument to -Q must be between 1 and 32");

			/* Force calibration */
			} else if (argv[fa][1] == 'J') {
				docalib = 1;

			/* No auto-cal */
			} else if (argv[fa][1] == 'N') {
				noautocal = 1;

			/* High res mode */
			} else if (argv[fa][1] == 'H') {
				highres = 1;

			/* Disable normalisation of values to white Y = 100 */ 
			} else if (argv[fa][1] == 'w') {
				donorm = 0;

			/* Colorimeter Correction Matrix */
			/* or Colorimeter Calibration Spectral Samples */
			} else if (argv[fa][1] == 'X') {
				int ix;
				fa = nfa;
				if (na == NULL) usage(0,"Parameter expected following -X");
				strncpy(ccxxname,na,MAXNAMEL-1); ccxxname[MAXNAMEL-1] = '\000';

			} else if (argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL || na[0] == '\000') usage(0,"Parameter expected after -I");
				for (i=0; ; i++) {
					if (na[i] == '\000')
						break;
					if (na[i] == 'b' || na[i] == 'B')
						bdrift = 1;
					else if (na[i] == 'w' || na[i] == 'W')
						wdrift = 1;
					else
						usage(0,"-I parameter '%c' not recognised",na[i]);
				}

			/* Spectral Observer type */
			} else if (argv[fa][1] == 'Q') {
				fa = nfa;
				if (na == NULL) usage(0,"Parameter expecte after -Q");
				if (strcmp(na, "1931_2") == 0) {			/* Classic 2 degree */
					obType = icxOT_CIE_1931_2;
				} else if (strcmp(na, "1964_10") == 0) {	/* Classic 10 degree */
					obType = icxOT_CIE_1964_10;
				} else if (strcmp(na, "2012_2") == 0) {		/* Latest 2 degree */
					obType = icxOT_CIE_2012_2;
				} else if (strcmp(na, "2012_10") == 0) {	/* Latest 10 degree */
					obType = icxOT_CIE_2012_10;
				} else if (strcmp(na, "1964_10c") == 0) {	/* 10 degree corrected */
					obType = icxOT_CIE_1964_10c;
				} else if (strcmp(na, "1955_2") == 0) {		/* Stiles and Burch 1955 2 degree */
					obType = icxOT_Stiles_Burch_2;
				} else if (strcmp(na, "1978_2") == 0) {		/* Judd and Voss 1978 2 degree */
					obType = icxOT_Judd_Voss_2;
				} else if (strcmp(na, "shaw") == 0) {		/* Shaw and Fairchilds 1997 2 degree */
					obType = icxOT_Shaw_Fairchild_2;
				} else {	/* Assume it's a filename */
					obType = icxOT_custom;
					if (read_cmf(custObserver, na) != 0)
						usage(0,"Failed to read custom observer CMF from -Q file '%s'",na);
				}

			/* Change color callout */
			} else if (argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage(0,"Parameter expected after -C");
				ccallout = na;

			/* Measure color callout */
			} else if (argv[fa][1] == 'M') {
				fa = nfa;
				if (na == NULL) usage(0,"Parameter expected after -M");
				mcallout = na;

			/* Request external values */
			} else if (argv[fa][1] == 'x') {
				fa = nfa;
				if (na == NULL) usage(0, "Parameter expected after -x");

				if (na[0] == 'x' || na[0] == 'X')
					xtern = 2;
//				else if (na[0] == 'l' || na[0] == 'L')
//					xtern = 1;
				else
					usage(0, "Unexpected parameter -x %c",na[0]);

			/* Serial port flow control */
			} else if (argv[fa][1] == 'W') {
				fa = nfa;
				if (na == NULL) usage(0,"Parameter expected after -W");
				if (na[0] == 'n' || na[0] == 'N')
					fc = fc_None;
				else if (na[0] == 'h' || na[0] == 'H')
					fc = fc_Hardware;
				else if (na[0] == 'x' || na[0] == 'X')
					fc = fc_XonXOff;
				else
					usage(0,"-W parameter '%s' not recognised",na);

			} else if (argv[fa][1] == 'D') {
				debug = 1;
				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					debug = atoi(na);
					fa = nfa;
				}
				g_log->debug = debug;
				callback_ddebug = 1;		/* dispwin global */

			/* Extra misc. flags */
			} else if (argv[fa][1] == 'Y') {
				if (na == NULL)
					usage(0,"Flag '-Y' expects extra flag");
			
				if (na[0] == 'R') {
					if (na[1] != ':')
						usage(0,"-Y R:rate syntax incorrect");
					refrate = atof(na+2);
					if (refrate < 5.0 || refrate > 150.0)
						usage(0,"-Y R:rate %f Hz not in valid range",refrate);
				} else if (na[0] == 'p') {
					noplace = 1;
				} else if (na[0] == 'A') {
					nadaptive = 1;
				} else if (na[0] == 'k') {
					nocaloutput = 1;
				} else {
					usage(0,"Flag '-Y %c' not recognised",na[0]);
				}
				fa = nfa;

			} else 
				usage(0,"Flag '-%c' not recognised",argv[fa][1]);
			}
		else
			break;
	}

	/* No explicit display has been set */
	if (!fake 
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
	}

	/* See if there is an environment variable ccxx */
	if (ccxxname[0] == '\000') {
		char *na;
		if ((na = getenv("ARGYLL_COLMTER_CAL_SPEC_SET")) != NULL) {
			strncpy(ccxxname,na,MAXNAMEL-1); ccxxname[MAXNAMEL-1] = '\000';

		} else if ((na = getenv("ARGYLL_COLMTER_COR_MATRIX")) != NULL) {
			strncpy(ccxxname,na,MAXNAMEL-1); ccxxname[MAXNAMEL-1] = '\000';
		}
	}

	/* Load up CCMX or CCSS */
	if (ccxxname[0] != '\000') {
		if ((cmx = new_ccmx()) == NULL
		  || cmx->read_ccmx(cmx, ccxxname)) {
			if (cmx != NULL) {
				cmx->del(cmx);
				cmx = NULL;
			}
			
			/* CCMX failed, try CCSS */
			if ((ccs = new_ccss()) == NULL
			  || ccs->read_ccss(ccs, ccxxname)) {
				if (ccs != NULL) {
					ccs->del(ccs);
					ccs = NULL;
					error("Reading CCMX/CCSS File '%s' failed\n", ccxxname);
				}
			}
		}
	}

	if (fake)
		comport = FAKE_DEVICE_PORT;
	if ((icmps = new_icompaths(g_log)) == NULL)
		error("Finding instrument paths failed");
	if ((ipath = icmps->get_path(icmps, comport)) == NULL)
		error("No instrument at port %d",comport);

	/* If we've requested ChromeCast, look it up */
	if (ccdisp) {
		if ((ccids = get_ccids()) == NULL)
			error("discovering ChromeCasts failed");
		if (ccids[0] == NULL)
			error("There are no ChromeCasts to use\n");
		for (i = 0; ccids[i] != NULL; i++)
			;
		if (ccdisp < 1 || ccdisp > i)
			error("Chosen ChromeCasts (%d) is outside list (1..%d)\n",ccdisp,i);
		ccid = ccids[ccdisp-1];
	}

	if (docalib) {
		
		if (verb) printf("Forcing a calibration\n");

		if ((rv = disprd_calibration(ipath, fc, ditype, -1, 0, tele, ambient, nadaptive, noautocal, 
			                         disp, webdisp, ccid,
#ifdef NT
			                         madvrdisp,
#endif
			                         dummydisp, out_tvenc, fullscreen, override,
			                         100.0 * hpatscale, 100.0 * vpatscale, ho, vo,
		                             g_log)) != 0) {
			error("docalibration failed with return value %d\n",rv);
		}
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage(0,"Filname parameter not found");
	strncpy(inname,argv[fa++],MAXNAMEL-4); inname[MAXNAMEL-4] = '\000';
	strcpy(outname,inname);
	strcat(inname,".ti1");
	strcat(outname,".ti3");

	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI1"); 	/* our special input type is Calibration Target Information 1 */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file '%s' isn't a CTI1 format file",inname);
	if (icg->ntables < 1)		/* We don't use second table at the moment */
		error ("Input file '%s' doesn't contain at least one table",inname);

	if ((npat = icg->t[0].nsets) <= 0)
		error ("Input file '%s' has no sets of data",inname);

	/* Setup output cgats file */
	ocg = new_cgats();					/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI3"); 		/* our special type is Calibration Target Information 3 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll dispread", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);	/* What sort of device this is */

	if ((ti = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
		ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[ti], NULL);
	
	if ((ti = icg->find_kword(icg, 0, "DARK_REGION_EMPHASIS")) >= 0)
		ocg->add_kword(ocg, 0, "DARK_REGION_EMPHASIS",icg->t[0].kdata[ti], NULL);

	if (verb) printf("Number of patches = %d\n",npat);

	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);

	if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
		error ("Input file '%s' doesn't contain field SAMPLE_ID",inname);
	if (icg->t[0].ftype[si] != nqcs_t)
		error ("Input file %s' field SAMPLE_ID is wrong type",inname);

	if ((cols = (col *)malloc(sizeof(col) * (npat+1))) == NULL)
		error("Malloc failed!");

	/* Figure out the color space */
	/* Read all the test patches in, and quantize them */
	if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file '%s' doesn't contain keyword COLOR_REP",inname);
	if (strcmp(icg->t[0].kdata[fi],"RGB") == 0) {
		int ri, gi, bi;
		double rgb[3];
		double qscale = (1 << qbits) - 1.0;
		dim = 3;
		if ((ri = icg->find_field(icg, 0, "RGB_R")) < 0)
			error ("Input file '%s' doesn't contain field RGB_R",inname);
		if (icg->t[0].ftype[ri] != r_t)
			error ("Input file '%s' field RGB_R is wrong type - expect float",inname);
		if ((gi = icg->find_field(icg, 0, "RGB_G")) < 0)
			error ("Input file '%s' doesn't contain field RGB_G",inname);
		if (icg->t[0].ftype[gi] != r_t)
			error ("Input file '%s' field RGB_G is wrong type - expect float",inname);
		if ((bi = icg->find_field(icg, 0, "RGB_B")) < 0)
			error ("Input file '%s' doesn't contain field RGB_B",inname);
		if (icg->t[0].ftype[bi] != r_t)
			error ("Input file '%s' field RGB_B is wrong type - expect float",inname);
		ocg->add_field(ocg, 0, "RGB_R", r_t);
		ocg->add_field(ocg, 0, "RGB_G", r_t);
		ocg->add_field(ocg, 0, "RGB_B", r_t);
		ocg->add_kword(ocg, 0, "COLOR_REP","RGB_XYZ", NULL);
		ocg->add_field(ocg, 0, "XYZ_X", r_t);
		ocg->add_field(ocg, 0, "XYZ_Y", r_t);
		ocg->add_field(ocg, 0, "XYZ_Z", r_t);
		for (i = 0; i < npat; i++) {
			cols[i].id = ((char *)icg->t[0].fdata[i][si]);
			rgb[0] = *((double *)icg->t[0].fdata[i][ri]) / 100.0;
			rgb[1] = *((double *)icg->t[0].fdata[i][gi]) / 100.0;
			rgb[2] = *((double *)icg->t[0].fdata[i][bi]) / 100.0;
			if (qbits > 0) {
				double vr;
				for (j = 0; j < 3; j++) {
					rgb[j] *= qscale;
					vr = floor(rgb[j] + 0.5);
					if ((vr - rgb[j]) == 0.5 && (((int)vr) & 1) != 0) /* Round to even */
						vr -= 1.0;
					rgb[j] = vr/qscale;
				}
			}
			cols[i].r = rgb[0];
			cols[i].g = rgb[1];
			cols[i].b = rgb[2];
			cols[i].XYZ[0] = cols[i].XYZ[1] = cols[i].XYZ[2] = -1.0;
		}
	} else
		error ("Input file '%s' keyword COLOR_REP has illegal value (RGB colorspace expected)",inname);

	/* Check that there is a white patch, and if not, add one, */
	/* so that we can normalize the values to white. */
	for (wpat = 0; wpat < npat; wpat++) {
		if (cols[wpat].r > 0.9999999 && 
		    cols[wpat].g > 0.9999999 && 
		    cols[wpat].b > 0.9999999) {
			break;
		}
	}
	if (wpat >= npat) {	/* Create a white patch */
		if (verb)
			printf("Adding one white patch\n");
		xpat = 1;
		cols[wpat].r = cols[wpat].g = cols[wpat].b = 1.0;
	}

	/* Setup a display calibration set if we are given one */
	/* (Should switch to xcal ?) */
	if (calname[0] != '\000') {
		cgats *ccg;			/* calibration cgats structure */
		int ii, ri, gi, bi;
		
		ccg = new_cgats();			/* Create a CGATS structure */
		ccg->add_other(ccg, "CAL"); /* our special calibration type */
	
		if (ccg->read_name(ccg, calname))
			error("CGATS calibration file read error %s on file '%s'",ccg->err,calname);
	
		if (ccg->ntables == 0 || ccg->t[0].tt != tt_other || ccg->t[0].oi != 0)
			error ("Calibration file isn't a CAL format file");
		if (ccg->ntables < 1)
			error ("Calibration file '%s' doesn't contain at least one table",calname);
	
		if ((ncal = ccg->t[0].nsets) <= 0)
			error ("No data in set of file '%s'",calname);
	
		if (ncal < 2 || ncal > MAX_CAL_ENT)
			error("Data set size %d is out of range for '%s'",ncal,calname);

		if (ncal > MAX_CAL_ENT)
			error ("Cant handle %d data sets in file '%s', max is %d",ncal,calname,MAX_CAL_ENT);
	
		if ((fi = ccg->find_kword(ccg, 0, "DEVICE_CLASS")) < 0)
			error ("Calibration file '%s' doesn't contain keyword DEVICE_CLASS",calname);
		if (strcmp(ccg->t[0].kdata[fi],"DISPLAY") != 0)
			error ("Calibration file '%s' doesn't have DEVICE_CLASS of DISPLAY",calname);

		if ((fi = ccg->find_kword(ccg, 0, "VIDEO_LUT_CALIBRATION_POSSIBLE")) >= 0) {
			if (stricmp(ccg->t[0].kdata[fi],"NO") == 0) {
				native = 1;
				if (verb) printf("Switching to soft cal because there is no access to VideoLUTs\n");
			}
		}

		if ((fi = ccg->find_kword(ccg, 0, "TV_OUTPUT_ENCODING")) >= 0) {
			if (!out_tvenc && (strcmp(ccg->t[0].kdata[fi], "YES") == 0
			 || strcmp(ccg->t[0].kdata[fi], "yes")) == 0) {
				if (verb) printf("Using Video range (16-235)/255 range encoding because cal. used it\n");
				out_tvenc = 1;
			}
		}

		if ((fi = ccg->find_kword(ccg, 0, "COLOR_REP")) < 0)
			error ("Calibration file '%s' doesn't contain keyword COLOR_REP",calname);
		if (strcmp(ccg->t[0].kdata[fi],"RGB") != 0)
			error ("Calibration file '%s' doesn't have COLOR_REP of RGB",calname);

		if ((ii = ccg->find_field(ccg, 0, "RGB_I")) < 0)
			error ("Calibration file '%s' doesn't contain field RGB_I",calname);
		if (ccg->t[0].ftype[ii] != r_t)
			error ("Field RGB_R in file '%s' is wrong type - expect float",calname);
		if ((ri = ccg->find_field(ccg, 0, "RGB_R")) < 0)
			error ("Calibration file '%s' doesn't contain field RGB_R",calname);
		if (ccg->t[0].ftype[ri] != r_t)
			error ("Field RGB_R in file '%s' is wrong type - expect float",calname);
		if ((gi = ccg->find_field(ccg, 0, "RGB_G")) < 0)
			error ("Calibration file '%s' doesn't contain field RGB_G",calname);
		if (ccg->t[0].ftype[gi] != r_t)
			error ("Field RGB_G in file '%s' is wrong type - expect float",calname);
		if ((bi = ccg->find_field(ccg, 0, "RGB_B")) < 0)
			error ("Calibration file '%s' doesn't contain field RGB_B",calname);
		if (ccg->t[0].ftype[bi] != r_t)
			error ("Field RGB_B in file '%s' is wrong type - expect float",calname);
		for (i = 0; i < ncal; i++) {
			cal[0][i] = *((double *)ccg->t[0].fdata[i][ri]);
			cal[1][i] = *((double *)ccg->t[0].fdata[i][gi]);
			cal[2][i] = *((double *)ccg->t[0].fdata[i][bi]);
		}
		ccg->del(ccg);
	} else {
		cal[0][0] = -1.0;	/* Not used */
	}

	if ((dr = new_disprd(&errc, ipath, fc, ditype, -1, 0, tele, ambient, nadaptive, noautocal, noplace,
	                     highres, refrate, native, &noramdac, &nocm, cal, ncal, disp,
		                 out_tvenc, fullscreen, override, webdisp, ccid,
#ifdef NT
						 madvrdisp,
#endif
		                 dummydisp, ccallout, mcallout, xtern,
		                 100.0 * hpatscale, 100.0 * vpatscale, ho, vo,
	                     ccs != NULL ? ccs->dtech : cmx != NULL ? cmx->dtech : disptech_unknown,
	                     cmx != NULL ? cmx->cc_cbid : 0,
	                     cmx != NULL ? cmx->matrix : NULL,
	                     ccs != NULL ? ccs->samples : NULL, ccs != NULL ? ccs->no_samp : 0,
		                 spec, obType, custObserver, bdrift, wdrift,
		                 "fake" ICC_FILE_EXT, g_log)) == NULL)
		error("new_disprd failed with '%s'\n",disprd_err(errc));

	if (cmx != NULL)
		cmx->del(cmx);
	if (ccs != NULL)
		ccs->del(ccs);

	/* Test the CRT with all of the test points */
	if ((rv = dr->read(dr, cols, npat + xpat, 1, npat + xpat, 1, 0, instNoClamp)) != 0) {
		dr->del(dr);
		error("dispd->read returned error code %d\n",rv);
	}
	/* Note what instrument the chart was read with */
	if (dr->it != NULL) {
		int refrmode, cbid;
		ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", inst_name(dr->it->get_itype(dr->it)) , NULL);
		dr->get_disptype(dr, &refrmode, &cbid);
		if (refrmode >= 0)
			ocg->add_kword(ocg, 0, "DISPLAY_TYPE_REFRESH", refrmode ? "YES" : "NO", NULL);
		if (cbid != 0) {
			char buf[100];
			sprintf(buf, "%d", cbid);
			ocg->add_kword(ocg, 0, "DISPLAY_TYPE_BASE_ID", buf, NULL);
		}
	} else {
		ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", "Fake" , NULL);
	}

	/* Note if the instrument is natively spectral or not */
	if (dr->it != NULL) {
		inst_mode cap;
		dr->it->capabilities(dr->it, &cap, NULL, NULL);

		if (dr->it != NULL && cap & inst_mode_spectral)
			ocg->add_kword(ocg, 0, "INSTRUMENT_TYPE_SPECTRAL", "YES" , NULL);
		else
			ocg->add_kword(ocg, 0, "INSTRUMENT_TYPE_SPECTRAL", "NO" , NULL);
	} else {
		ocg->add_kword(ocg, 0, "INSTRUMENT_TYPE_SPECTRAL", "NO" , NULL);
	}

	dr->del(dr);

	/* And save the result: */

	/* Convert from absolute XYZ to relative XYZ */
	if (npat > 0) {
		double nn;

		/* Make sure there is a copy of the white patch beyond npat, */
		/* so that it can be left absolute. */
		if (wpat != npat) {
			cols[npat].r = cols[npat].g = cols[npat].b = 1.0;
			cols[npat].XYZ[0] = cols[wpat].XYZ[0];
			cols[npat].XYZ[1] = cols[wpat].XYZ[1];
			cols[npat].XYZ[2] = cols[wpat].XYZ[2];
			cols[npat].XYZ_v = cols[wpat].XYZ_v; 
			wpat = npat;
		}

		if (donorm) {
			if (cols[wpat].XYZ_v == 0)
				error("XYZ of white patch is not valid!\nCan't normalise value to white",i);

			nn = 100.0 / cols[wpat].XYZ[1];		/* Normalise Y of white to 100 */
		} else
			nn = 1.0;

		for (i = 0; i < npat; i++) {

			if (cols[i].XYZ_v == 0)
				warning("XYZ patch %d is not valid!",i+1);

			for (j = 0; j < 3; j++)
				cols[i].XYZ[j] = nn * cols[i].XYZ[j];
		
			/* Keep spectral aligned with normalised XYZ */
			if (cols[i].sp.spec_n > 0) {
				for (j = 0; j < cols[i].sp.spec_n; j++)
					cols[i].sp.spec[j] *= nn;
			}
		}
	}

	nsetel += 1;		/* For id */
	nsetel += dim;		/* For device values */
	nsetel += 3;		/* For XYZ */

	/* If we have spectral information, output it too */
	if (npat > 0 && cols[0].sp.spec_n > 0) {
		char buf[100];

		nsetel += cols[0].sp.spec_n;		/* Spectral values */
		sprintf(buf,"%d", cols[0].sp.spec_n);
		ocg->add_kword(ocg, 0, "SPECTRAL_BANDS",buf, NULL);
		sprintf(buf,"%f", cols[0].sp.spec_wl_short);
		ocg->add_kword(ocg, 0, "SPECTRAL_START_NM",buf, NULL);
		sprintf(buf,"%f", cols[0].sp.spec_wl_long);
		ocg->add_kword(ocg, 0, "SPECTRAL_END_NM",buf, NULL);

		/* Generate fields for spectral values */
		for (i = 0; i < cols[0].sp.spec_n; i++) {
			int nm;
	
			/* Compute nearest integer wavelength */
			nm = (int)(cols[0].sp.spec_wl_short + ((double)i/(cols[0].sp.spec_n-1.0))
			            * (cols[0].sp.spec_wl_long - cols[0].sp.spec_wl_short) + 0.5);
			
			sprintf(buf,"SPEC_%03d",nm);
			ocg->add_field(ocg, 0, buf, r_t);
		}
	}

	if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * nsetel)) == NULL)
		error("Malloc failed!");

	/* Write out the patch info to the output CGATS file */
	for (i = 0; i < npat; i++) {
		int k = 0;

		if (cols[i].XYZ_v == 0) {
			warning("Omitting patch %d from .ti3 file!",i+1);
			continue;
		}

		setel[k++].c = cols[i].id;
		setel[k++].d = 100.0 * cols[i].r;
		setel[k++].d = 100.0 * cols[i].g;
		setel[k++].d = 100.0 * cols[i].b;

		setel[k++].d = cols[i].XYZ[0];
		setel[k++].d = cols[i].XYZ[1];
		setel[k++].d = cols[i].XYZ[2];

		for (j = 0; j < cols[i].sp.spec_n; j++) {
			setel[k++].d = cols[i].sp.spec[j];
		}

		ocg->add_setarr(ocg, 0, setel);
	}

	free(setel);

	/* If we have the absolute brightness of the display, record it */
	if (cols[wpat].XYZ_v != 0) {
		char buf[100];

		sprintf(buf,"%f %f %f", cols[wpat].XYZ[0], cols[wpat].XYZ[1], cols[wpat].XYZ[2]);
		ocg->add_kword(ocg, 0, "LUMINANCE_XYZ_CDM2",buf, NULL);
	}

	if (donorm) 
		ocg->add_kword(ocg, 0, "NORMALIZED_TO_Y_100","YES", NULL);
	else
		ocg->add_kword(ocg, 0, "NORMALIZED_TO_Y_100","NO", NULL);

	/* Write out the calibration if we have it and we want to save it */
	if (cal[0][0] >= 0.0 && !nocaloutput) {
		ocg->add_other(ocg, "CAL"); 		/* our special type is Calibration file */
		ocg->add_table(ocg, tt_other, 1);	/* Add another table for RAMDAC values */
		ocg->add_kword(ocg, 1, "DESCRIPTOR", "Argyll Device Calibration State",NULL);
		ocg->add_kword(ocg, 1, "ORIGINATOR", "Argyll dispread", NULL);
		ocg->add_kword(ocg, 1, "CREATED",atm, NULL);

		ocg->add_kword(ocg, 1, "DEVICE_CLASS","DISPLAY", NULL);
		ocg->add_kword(ocg, 1, "COLOR_REP","RGB", NULL);
		ocg->add_kword(ocg, 0, "VIDEO_LUT_CALIBRATION_POSSIBLE",noramdac ? "NO" : "YES", NULL);

		ocg->add_field(ocg, 1, "RGB_I", r_t);
		ocg->add_field(ocg, 1, "RGB_R", r_t);
		ocg->add_field(ocg, 1, "RGB_G", r_t);
		ocg->add_field(ocg, 1, "RGB_B", r_t);

		if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * 4)) == NULL)
			error("Malloc failed!");

		for (i = 0; i < ncal; i++) {
			double vv;

#if defined(__APPLE__) && defined(__POWERPC__)
			gcc_bug_fix(i);
#endif
			vv = i/(ncal-1.0);

			setel[0].d = vv;
			setel[1].d = cal[0][i];
			setel[2].d = cal[1][i];
			setel[3].d = cal[2][i];

			ocg->add_setarr(ocg, 1, setel);
		}

		free(setel);
	}

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	if (verb)
		printf("Written '%s'\n",outname);

	icmps->del(icmps);
	free(cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */
	free_a_disppath(disp);
	free_ccids(ccids);

	return 0;
}


