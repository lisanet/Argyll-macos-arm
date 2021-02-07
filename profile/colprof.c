/* 
 * Argyll Color Correction System
 * Color Device profile generator.
 *
 * Author: Graeme W. Gill
 * Date:   15/2/97
 *
 * Copyright 1996-2011 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in the scattered test chart
 * points, and interpolates them into a gridded
 * forward ICC device profile, as well as creating
 * backward conversions based on the forward grid.
 * 
 * Preview profiles are not currently generated.
 * 
 * The gamut cLUT should be implemented with xicc/rspl
 */

/*
 * TTBD:
 *
 *		Should add -ua option, to restore option to force input profile to be
 *      absolute intent for all ICC intents.
 *
 *		Should allow ICC Device attributes to be set.
 *
 *      Add Argyll private tag to record ink limit etc. to automate link parameters.
 *      Estimate ink limit from B2A tables if no private tag ?
 *      Add used option for black relative
 *      Add used option for separate high res reverse tables
 *      Fix 400% assumptions for > 4 color devices ?
 *
 *      Should allow creating profiles from .MPP directly for <= 4 dev channels.
 *      Should allow creating profiles from existing ICC profiles (deprecate revfix ?)
 *
 *      Should allow creating profiles >4 channels by providing .MPP for input,
 *      dev link .icm for psudo-dev to device & .ti3 for Pseudo-dev to PCS.
 *		Note gamut should come from psudo-dev to PCS.
 */

#undef DEBUG
#undef DO_TIME			/* Time the operation */

#define verbo stdout

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "cgats.h"
#include "xicc.h"
#include "prof.h"
#include "ui.h"

#define DEFAVGDEV 0.5		/* Default average deviation percentage */
							/* This equates to a uniform added error of +/- 1% */

#define DEMPH_DEFAULT 1.0	/* Default dark region emphasis == none */

/*

  Flags used:

         ABCDEFGHIJKLMNOPQRSTUVWXYZ
  upper  ....    . ... .. ....    .
  lower  .... .. . .. .........    

*/

void usage(char *diag, ...) {
	int i;
	fprintf(stderr,"Create ICC profile, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: %s [-options] inoutfile\n",error_program);
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -A manufacturer Manufacturer description string\n");
	fprintf(stderr," -M model        Model description string\n");
	fprintf(stderr," -D description  Profile Description string (Default \"inoutfile\")\n");
	fprintf(stderr," -C copyright    Copyright string\n");
	fprintf(stderr," -Z tmnb         Attributes: Transparency, Matte, Negative, BlackAndWhite\n");
	fprintf(stderr," -Z prsa         Default intent: Perceptual, Rel. Colorimetric, Saturation, Abs. Colorimetric\n");

	fprintf(stderr," -q lmhu         Quality - Low, Medium (def), High, Ultra\n");
//	fprintf(stderr," -q fmsu         Speed - Fast, Medium (def), Slow, Ultra Slow\n");
	fprintf(stderr," -b [lmhun]      Low quality B2A table - or specific B2A quality or none for input device\n");
//	fprintf(stderr," -b [fmsun]      B2A Speed - Fast, Medium, Slow, Ultra Slow, None, same as -q (def)\n");
	fprintf(stderr," -ni             Don't create input (Device) shaper curves\n");
	fprintf(stderr," -np             Don't create input (Device) grid position curves\n");
	fprintf(stderr," -no             Don't create output (PCS) shaper curves\n");
	fprintf(stderr," -nc             Don't put the input .ti3 data in the profile\n");
	fprintf(stderr," -k zhxr         Black Ink generation target: z = zero K,\n");
	fprintf(stderr,"                 h = 0.5 K, x = max K, r = ramp K (def.)\n");
	fprintf(stderr," -k p stle stpo enpo enle shape\n");
	fprintf(stderr,"                 stle: K level at White 0.0 - 1.0\n");
	fprintf(stderr,"                 stpo: start point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"                 enpo: End point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"                 enle: K level at Black 0.0 - 1.0\n");
	fprintf(stderr,"                 shape: 1.0 = straight, 0.0-1.0 concave, 1.0-2.0 convex\n");
	fprintf(stderr," -K parameters   Same as -k, but target is K locus rather than K value itself\n");
	fprintf(stderr," -l tlimit       override total ink limit, 0 - 400%% (default from .ti3)\n");
	fprintf(stderr," -L klimit       override black ink limit, 0 - 100%% (default from .ti3)\n");
	fprintf(stderr," -a lxXYgsmGS    Algorithm type override\n");
	fprintf(stderr,"                 l = Lab cLUT (def.), x = XYZ cLUT,\n");
	fprintf(stderr,"                 X = display XYZ cLUT + matrix, Y = display XYZ cLUT + debug matrix\n");
	fprintf(stderr,"                 g = gamma+matrix, s = shaper+matrix, m = matrix only,\n");
	fprintf(stderr,"                 G = single gamma+matrix, S = single shaper+matrix\n");
//  Development - not supported
//	fprintf(stderr," -I ver          Set ICC profile version > 2.2.0\n");
//	fprintf(stderr,"                 ver = 4, Enable ICC V4 creation\n");
	fprintf(stderr," -u              If input profile, auto scale WP to allow extrapolation\n");
	fprintf(stderr," -ua             If input profile, force absolute intent\n");
	fprintf(stderr," -uc             If input profile, clip cLUT values above WP\n");
	fprintf(stderr," -U scale        If input profile, scale media white point by scale\n");
	fprintf(stderr," -R              Restrict white <= 1.0, black and primaries to be +ve\n");
//	fprintf(stderr," -B X,Y,Z        Display Black Point override hack\n");
	fprintf(stderr," -V demphasis    Degree of dark region cLUT grid emphasis 1.0-4.0 (default %.2f = none)\n",DEMPH_DEFAULT);
	fprintf(stderr," -f [illum]      Use Fluorescent Whitening Agent compensation [opt. simulated inst. illum.:\n");
	fprintf(stderr,"                  M0, M1, M2, A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp]\n");
	fprintf(stderr," -i illum        Choose illuminant for computation of CIE XYZ from spectral data & FWA:\n");
	fprintf(stderr,"                  A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                  1931_2 (def), 1964_10, 2012_2, 2012_10, S&B 1955_2, shaw, J&V 1978_2, or file.cmf\n");
	fprintf(stderr," -r avgdev       Average deviation of device+instrument readings as a percentage (default %4.2f%%)\n",DEFAVGDEV);
/* Research options: */
/*	fprintf(stderr," -r sSMOOTH      RSPL or shaper suplimental optimised smoothing factor\n"); */
/*	fprintf(stderr," -r rSMOOTH      RSPL or shaper raw underlying smoothing factor\n"); */
	fprintf(stderr," -s src%s|cperc  Apply gamut mapping to output profile perceptual B2A table\n",ICC_FILE_EXT);
	fprintf(stderr,"                   for given source space, or compression percentage\n");
	fprintf(stderr," -S src%s|experc Apply gamut mapping to output profile perceptual and\n",ICC_FILE_EXT);
	fprintf(stderr,"                   and saturation B2A table, or expansion percentage\n");
	fprintf(stderr," -nP             Use colormetric source gamut to make output profile perceptual table\n");
	fprintf(stderr," -nS             Use colormetric source gamut to make output profile saturation table\n");
	fprintf(stderr," -g src.gam      Use source image gamut as well for output profile gamut mapping\n");
	fprintf(stderr," -p absprof,...  Incorporate abstract profile(s) into output tables\n");
	fprintf(stderr," -t intent       Override gamut mapping intent for output profile perceptual table:\n");
	fprintf(stderr," -T intent       Override gamut mapping intent for output profile saturation table:\n");
	for (i = 0; ; i++) {
		icxGMappingIntent gmi;
		if (xicc_enum_gmapintent(&gmi, i, NULL) == icxIllegalGMIntent)
			break;
		fprintf(stderr,"              %s\n",gmi.desc);
	}
	fprintf(stderr," -c viewcond     set input viewing conditions for output profile %s gamut mapping,\n",icxcam_description(cam_default));
	fprintf(stderr,"                  either an enumerated choice, or a parameter\n");
	fprintf(stderr," -d viewcond     set output viewing conditions for output profile %s gamut mapping\n",icxcam_description(cam_default));
	fprintf(stderr,"                  either an enumerated choice, or a parameter\n");
	fprintf(stderr,"                  Also sets out of gamut clipping CAM space.\n");
	fprintf(stderr,"                  either an enumerated choice, or a series of parameters:value changes\n");
	for (i = 0; ; i++) {
		icxViewCond vc;
		if (xicc_enum_viewcond(NULL, &vc, i, NULL, 1, NULL) == -999)
			break;

		fprintf(stderr,"             %s\n",vc.desc);
	}
	fprintf(stderr," -P              Create gamut gammap_p.wrl and gammap_s.wrl diagostics\n");
	fprintf(stderr," -O outputfile   Override the default output filename.\n");
	fprintf(stderr," inoutfile       Base name for input.ti3/output%s file\n",ICC_FILE_EXT);
	exit(1);
}

int main(int argc, char *argv[]) {
	int fa,nfa,mfa;				/* current argument we're looking at */
#ifdef DO_TIME			/* Time the operation */
	clock_t stime, ttime;		/* Start and total times */
#endif
	int verb = 0;
	int iquality = 1;			/* A2B quality */
	int oquality = -1;			/* B2A quality same as A2B */
	int verify = 0;				/* Not used anymore */
	int noisluts = 0;			/* No input shaper luts */
	int noipluts = 0;			/* No input position luts */
	int nooluts = 0;			/* No output shaper luts */
	int nocied = 0;				/* No .ti3 CIE data in profile */
	int noptop = 0;				/* Use colormetric source gamut to make perceptual table */
	int nostos = 0;				/* Use colormetric source gamut to make saturation table */
	int gamdiag = 0;			/* Make gamut mapping diagnostic wrl plots */
	int autowpsc = 0;			/* 1 = Auto scale the WP to prevent clipping above WP patch */
								/* 2 = Force absolute colorimetric */
	int clipovwp = 0;			/* Clip cLUT values above WP */
	int clipprims = 0;			/* Clip white, black and primaries */
//	double bpo[3] = { -1,-1,-1 };	/* Black point override hack XYZ value */
	double demph = 0.0;			/* Emphasise dark region grid resolution in cLUT */
	double iwpscale = -1.0;		/* Input white point scale factor */
	int doinextrap = 1;			/* Sythesize extra sample points for input device cLUT */
	int doinb2a = 1;			/* Create an input device B2A table */
	int inking = 3;				/* Default K target ramp K */
	int locus = 0;				/* Default K value target */
	double Kstle = 0.0, Kstpo = 0.0, Kenle = 0.0, Kenpo = 0.0, Kshap = 0.0;
	int tlimit = -1;			/* Total ink limit as a % */
	int klimit = -1;			/* Black ink limit as a % */
	int fwacomp = 0;			/* FWA compensation */
	double avgdev = DEFAVGDEV/100.0;	/* Average measurement deviation */
	double smooth = 1.0;		/* RSPL Smoothness factor (relative, for verification) */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType tillum = icxIT_none;	/* Target/simulated instrument illuminant */ 
	xspect cust_tillum;			/* Custom target/simulated illumination spectrum */
								/* xspect will use illum/cust_illum if tillum == none */
	icxIllumeType illum = icxIT_none;	/* Spectral illuminant (defaults to D50) */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType obType = icxOT_none;	/* Observer (defaults to 1931 2 degree) */
	xspect custObserver[3];		/* If obType = icxOT_custom */
	int gcompr = 0, gexpr = 0;		/* Gamut compression/expansion % instead of Input icc profile */
	char ipname[MAXNAMEL+1] = "";	/* Input icc profile - enables gamut map */
	char sgname[MAXNAMEL+1] = "";	/* Image source gamut name */
	char absstring[3 * MAXNAMEL +1];	/* Storage for absnames */
	char *absnames[3] = { NULL, NULL, NULL };	/* Abstract profile name */
	int sepsat = 0;				/* Create separate saturation B2A table */
	icxViewCond ivc_p;			/* Input Viewing Parameters for CAM */
	icxViewCond ovc_p;			/* Output Viewing Parameters for CAM (enables CAM clip) */
	int ivc_e = -1, ovc_e = -1;	/* Enumerated viewing condition */
	icxGMappingIntent pgmi;		/* default Perceptual gamut mapping intent */
	int pgmi_set = 0;			/* Set by user option */
	icxGMappingIntent sgmi;		/* default Saturation gamut mapping intent */
	int sgmi_set = 0;			/* Set by user option */
	char baname[MAXNAMEL+1] = "";	/* Input & Output base name */
	char inname[MAXNAMEL+1] = "";	/* Input cgats file base name */
	char outname[MAXNAMEL+1] = "";	/* Output cgats file base name */
	cgats *icg;					/* input cgats structure */
	int dti, ti;				/* Device type index, Temporary CGATs index */
	prof_atype ptype = prof_default;	/* Default for each type of device */
	int mtxtoo = 0;				/* 1 if matrix tags should be created for Display XYZ cLUT */
								/* 2 if debug matrix tags should be created for Display XYZ cLUT */
	icmICCVersion iccver = icmVersionDefault;	/* ICC profile version to create */
	profxinf xpi;		/* Extra profile information */

	
#ifdef DO_TIME			/* Time the operation */
	stime = clock();
#endif /* DO_TIME */
	error_program = argv[0];
	check_if_not_interactive();
	memset((void *)&xpi, 0, sizeof(profxinf));	/* Init extra profile info to defaults */
	xpi.default_ri = icMaxEnumIntent;			/* Default default */

	/* Init VC overrides so that we know when the've been set */
	ivc_p.Ev = -1;
	ivc_p.Wxyz[0] = -1.0; ivc_p.Wxyz[1] = -1.0; ivc_p.Wxyz[2] = -1.0;
	ivc_p.La = -1.0;
	ivc_p.Yb = -1.0;
	ivc_p.Lv = -1.0;
	ivc_p.Yf = -1.0;
	ivc_p.Yg = -1.0;
	ivc_p.Gxyz[0] = -1.0; ivc_p.Gxyz[1] = -1.0; ivc_p.Gxyz[2] = -1.0;
	ivc_p.hkscale = -1.0;
	ivc_p.mtaf = -1.0;
	ivc_p.Wxyz2[0] = -1.0; ivc_p.Wxyz2[1] = -1.0; ivc_p.Wxyz2[2] = -1.0;

	ovc_p.Ev = -1;
	ovc_p.Wxyz[0] = -1.0; ovc_p.Wxyz[1] = -1.0; ovc_p.Wxyz[2] = -1.0;
	ovc_p.La = -1.0;
	ovc_p.Yb = -1.0;
	ovc_p.Lv = -1.0;
	ovc_p.Yf = -1.0;
	ovc_p.Yg = -1.0;
	ovc_p.Gxyz[0] = -1.0; ovc_p.Gxyz[1] = -1.0; ovc_p.Gxyz[2] = -1.0;
	ovc_p.hkscale = -1.0;
	ovc_p.mtaf = -1.0;
	ovc_p.Wxyz2[0] = -1.0; ovc_p.Wxyz2[1] = -1.0; ovc_p.Wxyz2[2] = -1.0;

	xicc_enum_gmapintent(&pgmi, icxPerceptualGMIntent, NULL);
	xicc_enum_gmapintent(&sgmi, icxSaturationGMIntent, NULL);

	if (argc <= 1)
		usage("Too few arguments, got %d expect at least %d",argc-1,1);

	/* Process the arguments */
	mfa = 1;		/* Minimum final arguments */
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
				usage("Usage requested");

			else if (argv[fa][1] == 'v')
				verb = 1;

			/* Manufacturer description string */
			else if (argv[fa][1] == 'A') {
				if (na == NULL) usage("Expect argument to manufacturer description flag -A");
				fa = nfa;
				xpi.deviceMfgDesc = na;
			}

			/* Model description string */
			else if (argv[fa][1] == 'M') {
				if (na == NULL) usage("Expect argument to model description flag -M");
				fa = nfa;
				xpi.modelDesc = na;
			}

			/* Profile Description */
			else if (argv[fa][1] == 'D') {
				if (na == NULL) usage("Expect argument to profile description flag -D");
				fa = nfa;
				xpi.profDesc = na;
			}

			/* Copyright string */
			else if (argv[fa][1] == 'C') {
				if (na == NULL) usage("Expect argument to copyright flag -C");
				fa = nfa;
				xpi.copyright = na;
			}

			/* Attribute bits */
			else if (argv[fa][1] == 'Z') {
				int i, j;

				if (na == NULL) usage("Expect argument to attribute flag -Z");
				fa = nfa;
				for (j = i = 0; j == 0; i++) {

					switch(na[i]) {

						/* Attribute */
						case 't':
						case 'T':
							xpi.transparency = 1;
							break;
						case 'm':
						case 'M':
							xpi.matte = 1;
							break;
						case 'n':
						case 'N':
							xpi.negative = 1;
							break;
						case 'b':
						case 'B':
							xpi.blackandwhite = 1;
							break;

						/* Default intent */
						case 'p':
						case 'P':
							xpi.default_ri = icPerceptual;
							break;
						case 'r':
						case 'R':
							xpi.default_ri = icRelativeColorimetric;
							break;
						case 's':
						case 'S':
							xpi.default_ri = icSaturation;
							break;
						case 'a':
						case 'A':
							xpi.default_ri = icAbsoluteColorimetric;
							break;
						default:
							j = 1;
							break;
					}
				}
				if (na[i-1] != '\000') usage("Unknown argument '%c' to attribute flag -Z",na[i-1]);
			}

			/* Quality */
			else if (argv[fa][1] == 'q') {
//				if (na == NULL) usage("Expect argument to quality flag -q");
				if (na == NULL) usage("Expect argument to speed flag -q");
				fa = nfa;
   	 			switch (na[0]) {
					case 'f':				/* Fast */
					case 'l':
					case 'L':
						iquality = 0;
						break;
					case 'm':				/* Medium */
					case 'M':
						iquality = 1;
						break;
					case 's':				/* Slow */
					case 'h':
					case 'H':
						iquality = 2;
						break;
					case 'u':				/* Ultra Slow */
					case 'U':
						iquality = 3;
						break;
					default:
						usage("Unknown argument '%c' to quality flag -q",na[0]);
//						usage("Unknown argument '%c' to speed flag -q",na[0]);
				}
			}
			else if (argv[fa][1] == 'b') {
				if (na != NULL) {	/* Got a B2A quaiity */
					fa = nfa;
	    			switch (na[0]) {
						case 'f':				/* Fast */
						case 'l':
						case 'L':
							oquality = 0;
							break;
						case 'm':				/* Medium */
						case 'M':
							oquality = 1;
							break;
						case 's':				/* Slow */
						case 'h':
						case 'H':
							oquality = 2;
							break;
						case 'u':				/* Ultra Slow */
						case 'U':
							oquality = 3;
							break;
						case 'n':				/* No B2A for input device */
						case 'N':
							oquality = -2;
							doinb2a = 0;
							break;
						default:
							usage("Unknown argument '%c' to quality flag -q",na[0]);
					}
				} else
					oquality = 0;
			}

			/* Disable input or output luts */
			else if (argv[fa][1] == 'n') {
				if (na == NULL) {	/* Backwards compatible */
					nooluts = 1;
				} else {
					fa = nfa;
					if (na[0] == 'i')
						noisluts = 1;
					else if (na[0] == 'p')
						noipluts = 1;
					else if (na[0] == 'o')
						nooluts = 1;
					else if (na[0] == 'c')
						nocied = 1;
					else if (na[0] == 'P')
						noptop = 1;
					else if (na[0] == 'S')
						nostos = 1;
					else
						usage("Unknown argument '%c' to flag -n",na[0]);
				}
			}

			else if (argv[fa][1] == 'u') {
				autowpsc = 1;
				clipovwp = 0;
				if (argv[fa][2] == 'a') {
					autowpsc = 2;
					clipovwp = 0;
				} else if (argv[fa][2] == 'c') {
					autowpsc = 0;
					clipovwp = 1;
				} else if (argv[fa][2] != '\000') {
					usage("Unknown flag '%c' after -u",argv[fa][2]);
				}
			}
			else if (argv[fa][1] == 'U') {
				if (na == NULL) usage("Expected argument to input white point scale flag -U");
				fa = nfa;
				iwpscale = atof(na);
				if (iwpscale < 0.0 || iwpscale > 200.0)
					usage("Argument '%s' to flag -U out of range",na);
			}
			/* Clip primaries */
			else if (argv[fa][1] == 'R') {
				clipprims = 1;
			}

#ifdef NEVER	/* Prototype - not used */
			/* Black Point override hack */
			else if (argv[fa][1] == 'B') {
				if (na == NULL) usage("Expect X,Y,Z value after -B");
				fa = nfa;
				if (sscanf(na, " %lf , %lf , %lf ",&bpo[0],&bpo[1],&bpo[2]) != 3)
					usage("Couldn't parse hack black point (-B) value '%s'",na);
				if (bpo[0] < 0.0 || bpo[1] < 0.0 || bpo[1] < 0.0)
					usage("Bad hack black point (-B) value '%s'",na);
			}
#endif

			/* Degree of dark region emphasis */
			else if (argv[fa][1] == 'V') {
				if (na == NULL) usage(0,"Expected argument to dark emphasis flag -V");
				fa = nfa;
				demph = atof(na);
				if (demph < 1.0 || demph > 3.0)
					usage("Dark weighting argument %f to '-V' is out of range",demph);
			}

			/* Inking rule */
			else if (argv[fa][1] == 'k'
			      || argv[fa][1] == 'K') {
				if (argv[fa][1] == 'k')
					locus = 0;		/* Use K value target */
				else
					locus = 1;		/* Use K locus target */
				if (na == NULL) usage("Expect argument to inking flag -k");
				fa = nfa;
    			switch (na[0]) {
					case 'z':
					case 'Z':
						inking = 0;		/* Use minimum k */
						break;
					case 'h':
					case 'H':
						inking = 1;		/* Use 0.5 k */
						break;
					case 'x':
					case 'X':
						inking = 2;		/* Use maximum K */
						break;
					case 'r':
					case 'R':
						inking = 3;		/* Use ramping K */
						break;
					case 'p':
					case 'P':
						inking = 4;		/* Use parameter curve */
						++fa;
						if (fa >= argc) usage("Too few arguments to inking flag -kp");
						Kstle = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage("Too few arguments to inking flag -kp");
						Kstpo = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage("Too few arguments to inking flag -kp");
						Kenpo = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage("Too few arguments to inking flag -kp");
						Kenle = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage("Too few arguments to inking flag -kp");
						Kshap = atof(argv[fa]);
						break;
					default:
						usage("Unknown inking rule (-k) argument '%c'",na[0]);
				}
			}

			/* Total Ink Limit */
			else if (argv[fa][1] == 'l') {
				if (na == NULL) usage("Expected argument to total ink limit flag -l");
				fa = nfa;
				tlimit = atoi(na);
			}

			/* Black Ink Limit */
			else if (argv[fa][1] == 'L') {
				if (na == NULL) usage("Expected argument to black ink limit flag -L");
				fa = nfa;
				klimit = atoi(na);
			}

			/* Algorithm type */
			else if (argv[fa][1] == 'a') {
				if (na == NULL) usage("Expect argument to algorithm flag -a");
				fa = nfa;
    			switch (na[0]) {
					case 'l':
					case 'L':
						ptype = prof_clutLab;
						break;
					case 'x':
						ptype = prof_clutXYZ;
						break;
					case 'X':
						mtxtoo = 1;
						ptype = prof_clutXYZ;
						break;
					case 'Y':
						mtxtoo = 2;
						ptype = prof_clutXYZ;
						break;
					case 'g':
						ptype = prof_gammat;
						break;
					case 'G':
						ptype = prof_gam1mat;
						break;
					case 's':
						ptype = prof_shamat;
						break;
					case 'S':
						ptype = prof_sha1mat;
						break;
					case 'm':
						ptype = prof_matonly;
						break;
					default:
						usage("Unknown argument '%c' to algorithm flag -a",na[0] );
				}
			}
			/* Profile version */
			else if (argv[fa][1] == 'I') {
				if (na == NULL) usage("Expect argument to version flag -I");
				fa = nfa;
    			switch (na[0]) {
					case '4':
						iccver = icmVersion4_1;
						break;
					default:
						usage("Unknown argument '%c' to version flag -I",na[0] );
				}
			}

			/* FWA compensation */
			else if (argv[fa][1] == 'f') {
				fwacomp = 1;
				spec = 1;		/* Have to use spectral data */

				if (na != NULL) {	/* Argument is present - target/simulated instr. illum. */
					fa = nfa;
					if (strcmp(na, "A") == 0
					 || strcmp(na, "M0") == 0) {
						tillum = icxIT_A;
					} else if (strcmp(na, "C") == 0) {
						tillum = icxIT_C;
					} else if (strcmp(na, "D50") == 0
					        || strcmp(na, "M1") == 0) {
						tillum = icxIT_D50;
					} else if (strcmp(na, "D50M2") == 0
					        || strcmp(na, "M2") == 0) {
						tillum = icxIT_D50M2;
					} else if (strcmp(na, "D65") == 0) {
						tillum = icxIT_D65;
					} else if (strcmp(na, "F5") == 0) {
						tillum = icxIT_F5;
					} else if (strcmp(na, "F8") == 0) {
						tillum = icxIT_F8;
					} else if (strcmp(na, "F10") == 0) {
						tillum = icxIT_F10;
					} else {	/* Assume it's a filename */
						inst_meas_type mt;

						tillum = icxIT_custom;
						if (read_xspect(&cust_tillum, &mt, na) != 0)
							usage("Failed to read custom target illuminant spectrum in file '%s'",na);

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

			/* CIE Illuminant type */
			else if (argv[fa][1] == 'i') {
				if (na == NULL) usage("Expect argument to illuminant flag -i");
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
						usage("Failed to read custom illuminant spectrum in file '%s'",na);

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
				if (na == NULL) usage("Expect argument to observer flag -o");
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
						usage(0,"Failed to read custom observer CMF from -o file '%s'",na);
				}
			}

			/* Average Deviation percentage */
			else if (argv[fa][1] == 'r') {
				if (na == NULL) usage("Expected argument to average deviation flag -r");
				fa = nfa;
				if (na[0] == 's') {			/* (relative, for verification) */
					smooth = atof(na+1);
					if (smooth < 0.0)
						usage("Optimised smoothing factor argument to '-rs' must be over 0.0");
				} else if (na[0] == 'r') {	/* (absolute, for testing) */
					smooth = atof(na+1);
					if (smooth < 0.0)
						usage("Raw smoothing factor argument to '-rr' must be over 0.0");
					smooth = -smooth;		/* Signal raw factor */
				} else {
					avgdev = 0.01 * atof(na);
					if (avgdev < 0.0 || avgdev > 1.0)
						usage("Average Deviation argument must be between 0.0 and 100.0");
				}
			}

			/* Percetual Source Compression and Perceptual/Saturation Gamut Maping mode enable */
			/* Percetual Source Gamut and Perceptual/Saturation Gamut Maping mode enable */
			else if (argv[fa][1] == 's'
			      || argv[fa][1] == 'S') {
				int sat = 0;
				if (argv[fa][1] == 'S') {
					sepsat = 1;
					sat = 1;
				} else {
					sepsat = 0;
				}
				if (na == NULL)
					usage("Unrecognised argument to source gamut flag -%c",argv[fa][1]);
				fa = nfa;
				{
					char *cp;
					int comp = 0;
					for (cp = na; *cp != '\000'; cp++) {
						if (!isdigit(*cp))
							break;
					}
					/* If general gamut compression/expansion mode */
					if (*cp == '\000') {
						int ratio = atoi(na);
						if (ratio <= 0 || ratio > 100) {
							usage("Gamut -%c %s %d%% must be > 0 and < 100",
							       argv[fa][1], sat ? "expansion" : "compression", ratio);
						}
						if (sat) {
							gexpr = ratio;
							if (gcompr == 0)
								gcompr = ratio;
						} else {
							gcompr = ratio;
							if (gexpr == 0)
								gexpr = ratio;
						}
					}
				}
				/* Not compression %, so assume it's a filename. */
				/* We allow both filename and general compression to allow */
				/* for a non-default L mapping. */
				if (gcompr == 0) {
					strncpy(ipname,na,MAXNAMEL); ipname[MAXNAMEL] = '\000';
				}
			}

			/* Source image gamut */
			else if (argv[fa][1] == 'g') {
				if (na == NULL)
					usage("Unrecognised argument to source image gamut flag -g",argv[fa][1]);
				fa = nfa;
				strncpy(sgname,na,MAXNAMEL); sgname[MAXNAMEL] = '\000';
			}

			/* Abstract profile */
			else if (argv[fa][1] == 'p') {
				char *f1 = NULL, *f2 = NULL;
				if (na == NULL) usage("Expected abstract profile filename after -p");
				fa = nfa;
				strncpy(absstring,na,MAXNAMEL*3); absstring[MAXNAMEL*3] = '\000';
				if ((f1 = strchr(absstring, ',')) == NULL) {		/* Only one profile */
					absnames[2] = absnames[1] = absnames[0] = absstring;	/* Duplicate */
				} else {	/* At least one comma */
					*f1++ = '\000';
					if ((f2 = strchr(f1, ',')) != NULL)		/* Two commas */
						*f2++ = '\000';
					if (*absstring != '\000')
						absnames[0] = absstring;
					if (*f1 != '\000')
						absnames[1] = f1;
					if (f2 != NULL && *f2 != '\000')
						absnames[2] = f2;
				}
			}

			/* Perceptual Mapping intent override */
			else if (argv[fa][1] == 't') {
				if (na == NULL) usage("Expect argument to perceptul intent override flag -t");
				fa = nfa;
				if (xicc_enum_gmapintent(&pgmi, icxNoGMIntent, na) == icxIllegalGMIntent)
					usage("Unrecognised intent '%s' to perceptual override flag -t",na);
				pgmi_set = 1;
			}

			/* Saturation Mapping intent override */
			else if (argv[fa][1] == 'T') {
				if (na == NULL) usage("Expect argument to saturation intent override flag -T");
				fa = nfa;
				if (xicc_enum_gmapintent(&sgmi, icxNoGMIntent, na) == icxIllegalGMIntent)
					usage("Unrecognised intent '%s' to saturation override flag -T",na);
				sgmi_set = 1;
			}

			/* Viewing conditions */
			else if (argv[fa][1] == 'c'
			      || argv[fa][1] == 'd') {
				icxViewCond *vc;

				if (argv[fa][1] == 'c') {
					vc = &ivc_p;
				} else {
					vc = &ovc_p;
				}

				if (na == NULL) usage("Viewing conditions flag (-%c) needs an argument",argv[fa][1]);
				fa = nfa;
				if (na[1] != ':') {
					if (vc == &ivc_p) {
						if ((ivc_e = xicc_enum_viewcond(NULL, NULL, -2, na, 1, NULL)) == -999)
							usage("Urecognised Enumerated Viewing conditions '%s'",na);
					} else {
						if ((ovc_e = xicc_enum_viewcond(NULL, NULL, -2, na, 1, NULL)) == -999)
							usage("Urecognised Enumerated Viewing conditions '%s'",na);
					}
				} else if (na[0] == 's' || na[0] == 'S') {
					if (na[1] != ':')
						usage("Viewing conditions (-%cs) missing ':'",argv[fa][1]);
					if (na[2] == 'n' || na[2] == 'N') {
						vc->Ev = vc_none;		/* Automatic */
					} else if (na[2] == 'a' || na[2] == 'A') {
						vc->Ev = vc_average;
					} else if (na[2] == 'm' || na[2] == 'M') {
						vc->Ev = vc_dim;
					} else if (na[2] == 'd' || na[2] == 'D') {
						vc->Ev = vc_dark;
					} else if (na[2] == 'c' || na[2] == 'C') {
						vc->Ev = vc_cut_sheet;
					} else
						usage("Viewing condition (-%c) unrecognised surround '%c'",argv[fa][1],na[2]);
				} else if (na[0] == 'w' || na[0] == 'W') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc->Wxyz[0] = x; vc->Wxyz[1] = y; vc->Wxyz[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc->Wxyz[0] = x; vc->Wxyz[1] = y; vc->Wxyz[2] = -1;
					} else
						usage("Viewing condition (-%cw) unrecognised white point '%s'",argv[fa][1],na+1);
				} else if (na[0] == 'a' || na[0] == 'A') {
					if (na[1] != ':')
						usage("Viewing conditions (-%ca) missing ':'",argv[fa][1]);
					vc->La = atof(na+2);
				} else if (na[0] == 'b' || na[0] == 'B') {
					if (na[1] != ':')
						usage("Viewing conditions (-%cb) missing ':'",argv[fa][1]);
					vc->Yb = atof(na+2)/100.0;
				} else if (na[0] == 'l' || na[0] == 'L') {
					if (na[1] != ':')
						usage("Viewing conditions (-%cl) missing ':'",argv[fa][1]);
					vc->Lv = atof(na+2);
				} else if (na[0] == 'f' || na[0] == 'F') {
					if (na[1] != ':')
						usage("Viewing conditions (-%cf) missing ':'",argv[fa][1]);
					vc->Yf = atof(na+2)/100.0;
				} else if (na[0] == 'g' || na[0] == 'G') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc->Gxyz[0] = x; vc->Gxyz[1] = y; vc->Gxyz[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc->Gxyz[0] = x; vc->Gxyz[1] = y; vc->Gxyz[2] = -1;
					} else if (sscanf(na+1,":%lf",&x) == 1) {
						vc->Yg = x/100.0;
					} else
						usage("Viewing condition (-%cg) unrecognised flare '%s'",argv[fa][1],na+1);

				} else if (na[0] == 'h' || na[0] == 'H') {
					if (na[1] != ':')
						usage("Viewing conditions (-%ch) missing ':'",argv[fa][1]);
					vc->hkscale = atof(na+2);
				} else if (na[0] == 'm' || na[0] == 'M') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc->Wxyz2[0] = x; vc->Wxyz2[1] = y; vc->Wxyz2[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc->Wxyz2[0] = x; vc->Wxyz2[1] = y; vc->Wxyz2[2] = -1;
					} else if (sscanf(na+1,":%lf",&x) == 1) {
						vc->mtaf = x;
					} else
						usage("Viewing condition (-%cm) unrecognised flare '%s'",argv[fa][1],na+1);

				} else
					usage("Viewing condition (-%c) unrecognised sub flag '%c'",argv[fa][1],na[0]);
			}

			/* Gammut mapping diagnostic plots */
			else if (argv[fa][1] == 'P')
				gamdiag = 1;

			/* Output file name */
			else if (argv[fa][1] == 'O') {
				if (na == NULL) usage("Output filename override (-O) needs an argument");
				fa = nfa;
				strncpy(outname,na,MAXNAMEL); outname[MAXNAMEL] = '\000';
			}

			else 
				usage("Unknown flag '%c'",argv[fa][1]);
		} else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage("Missing input .ti3 and output ICC basename");
	strncpy(baname,argv[fa++],MAXNAMEL-4); baname[MAXNAMEL-4] = '\000';
	if (xpi.profDesc == NULL) {
		char *cp;
		if ((cp = strrchr(baname, '/')) == NULL
		 && (cp = strrchr(baname, '\\')) == NULL)
			cp = baname;
		else
			cp++;
		xpi.profDesc = cp;	/* Default description = file name of base path */
	}
	strcpy(inname, baname);
	strcat(inname,".ti3");
	if (outname[0] == '\000') {		/* If not overridden */
		strcpy(outname, baname);
		strcat(outname, ICC_FILE_EXT);
	}

	/* Issue some errors & warnings for strange combinations */
	if (fwacomp && spec == 0)
		error("FWA compensation only works when viewer and/or illuminant selected");

	if (pgmi_set && (ipname[0] == '\000' && gcompr == 0))
		warning("-t perceptual intent override only works if -s srcprof or -S srcprof is used");

	if (sgmi_set && (ipname[0] == '\000' && gcompr == 0))
		warning("-T saturation intent override only works if -S srcprof is used");

	if (sgmi_set && sepsat == 0) {	/* Won't do much otherwise */
		if (verb)
			printf("Saturation intent override was set, so adding saturation intent table\n");
		sepsat = 1;
	}

	if (gamdiag && (ipname[0] == '\000' && gcompr == 0))
		warning("No gamut mapping called for, so -P will produce nothing");
		
	if (sgname[0] != '\000' && ipname[0] == '\000')
		warning("-g srcgam will do nothing without -s srcprof or -S srcprof");

	/* Should warning input viewing condition set and ipname[0] == '\000' */

	if (oquality == -1) {		/* B2A tables will be used */
		oquality = iquality;
	}

	/* Open and look at the .ti3 profile patches file */
	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI3"); 	/* our special input type is Calibration Target Information 3 */
	icg->add_other(icg, "CAL"); 	/* our special device Calibration state */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI3 format file");
	if (icg->ntables < 1)
		error ("Input file doesn't contain at least one table");

	/* Read the device/profile type tag */
	if ((dti = icg->find_kword(icg, 0, "DEVICE_CLASS")) < 0)
		error ("Input file doesn't contain keyword DEVICE_CLASS");

	/* Issue some errors & warnings for strange combinations for profile type */

	/* Reflective options when not a reflective profile type */
	if (strcmp(icg->t[0].kdata[dti],"DISPLAY") == 0
	 || strcmp(icg->t[0].kdata[dti],"EMISINPUT") == 0) {
		if (illum != icxIT_none)
			warning("-i illuminant ignored for emissive reference type");
		if (fwacomp != 0)
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
			error("Neither CIE nor spectral data found in file '%s'",inname);

		/* Switch to using spectral information */
		if (verb)
			printf("No CIE data found, switching to spectral with standard observer & D50\n");
		spec = 1;
		illum = icxIT_D50;
		obType = icxOT_CIE_1931_2;
	}
	
	/* If we requested spectral, check that it is available */
	if (spec) {
		if (icg->find_kword(icg, 0, "SPECTRAL_BANDS") < 0)
			error ("Requested spectral interpretation when data not available");
	}

	/* If not set in arguments, set default demph from .ti1 or default none */
	if (demph == 0.0) {
		if ((ti = icg->find_kword(icg, 0, "DARK_REGION_EMPHASIS")) >= 0) {
			demph = atof(icg->t[0].kdata[ti]);
			demph = pow(demph, 0.7);		/* Reduce intensity */
			if (verb)
				printf("Dark emphasis factor %f from targen\n",demph);
		} else {
			demph = DEMPH_DEFAULT;
		}
	}

	/* check the device class, and call function to create profile. */
	if (strcmp(icg->t[0].kdata[dti],"OUTPUT") == 0) {		/* i.e. printer */
		icxInk ink;							/* Ink parameters */

//		if (bpo[1] >= 0.0)
//			error("-B option not valid for output profile");

		if ((ti = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0) {
			int imax;
			imax = atoi(icg->t[0].kdata[ti]);
			if (imax > 0 && imax <= 400.0) {
				if (tlimit > 0 && tlimit <= 400.0) {	/* User has specified limit as option */
					if (imax < tlimit) {
						warning("Ink limit greater than original chart! (%d%% > %d%%)",tlimit,imax);
					}
				} else {
					if (imax > 80.0)
						tlimit = (int)imax - 10;	/* Rule of thumb - 10% below chart maximum */
					else
						tlimit = (int)imax;
				}
			}
		}

		/* (Note that this isn't set by any of the Argyll tools currently, */
		/*  but can be set manually.) */
		if ((ti = icg->find_kword(icg, 0, "BLACK_INK_LIMIT")) >= 0) {
			int kmax;
			kmax = atoi(icg->t[0].kdata[ti]);
			if (kmax > 0 && kmax <= 100.0) {
				if (klimit >= 0 && klimit <= 100.0) {	/* User has specified limit as option */
					if (kmax < klimit) {
						warning("Black ink limit greater than original chart! (%d%% > %d%%)",klimit,kmax);
					}
				} else {
					klimit = (int)kmax;
				}
			}
		}

		if (tlimit >= 0 && tlimit < 400.0) {
			if (verb)
				printf("Total ink limit being used is %d%%\n",tlimit);
			ink.tlimit = tlimit/100.0;	/* Set a total ink limit */
		} else {
			if (verb)
				printf("No total ink limit being used\n");
			ink.tlimit = -1.0;			/* Don't use a limit */
		}

		if (klimit >= 0 && klimit < 100.0) {
			if (verb)
				printf("Black ink limit being used is %d%%\n",klimit);
			ink.klimit = klimit/100.0;	/* Set a black ink limit */
		} else {
			if (verb)
				printf("No black ink limit being used\n");
			ink.klimit = -1.0;			/* Don't use a limit */
		}

		ink.KonlyLmin = 0;				/* Use normal black Lmin for locus */
		ink.c.Ksmth = ICXINKDEFSMTH;	/* default black curve smoothing */
		ink.c.Kskew = ICXINKDEFSKEW;	/* default black curve skew */
		ink.x.Ksmth = ICXINKDEFSMTH;
		ink.x.Kskew = ICXINKDEFSKEW;

		if (inking == 0) {			/* Use minimum */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = 0.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 0.0;
			ink.c.Kshap = 1.0;
		} else if (inking == 1) {	/* Use 0.5  */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = 0.5;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 0.5;
			ink.c.Kshap = 1.0;
		} else if (inking == 2) {	/* Use maximum  */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = 1.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 1.0;
			ink.c.Kshap = 1.0;
		} else if (inking == 3) {	/* Use ramp  */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = 0.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 1.0;
			ink.c.Kshap = 1.0;
		} else {				/* Use specified curve */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = Kstle;
			ink.c.Kstpo = Kstpo;
			ink.c.Kenpo = Kenpo;
			ink.c.Kenle = Kenle;
			ink.c.Kshap = Kshap;
		}

		if (ptype == prof_default)
			ptype = prof_clutLab;
		else if (ptype != prof_clutLab && ptype != prof_clutXYZ) {
			error ("Output profile can only be a cLUT algorithm");
		}

		if (autowpsc == 1)
			error ("Input auto WP scale mode isn't applicable to an output device");
		if (autowpsc == 2)
			error ("Force absolute colorimetric isn't applicable to an output device");
		if (clipovwp)
			error ("Input cLUT clipping above WP mode isn't applicable to an output device");

		make_output_icc(ptype, 0, iccver, verb, iquality, oquality,
		                noisluts, noipluts, nooluts, nocied, noptop, nostos,
		                gamdiag, verify, clipprims, iwpscale,
//		                NULL,		/* bpo */
		                &ink, inname, outname, icg, 
		                spec, tillum, &cust_tillum, illum, &cust_illum, obType, custObserver,
						fwacomp, smooth, avgdev, 1.0,
						gcompr, gexpr,
		                ipname[0] != '\000' ? ipname : NULL,
		                sgname[0] != '\000' ? sgname : NULL,
		                absnames,
						sepsat, &ivc_p, &ovc_p, ivc_e, ovc_e,
						&pgmi, &sgmi, &xpi);

	/* Scanner, Camera with reflective target, or Camera with emissive target */
	} else if (strcmp(icg->t[0].kdata[dti],"INPUT") == 0
	        || strcmp(icg->t[0].kdata[dti],"EMISINPUT") == 0) {

		int emis = strcmp(icg->t[0].kdata[dti],"EMISINPUT") == 0;

//		if (bpo[1] >= 0.0)
//			error("-B option not valid for input profile");

		if (ptype == prof_default)
			ptype = prof_clutLab;		/* For best possible quality */

		if (clipovwp && ptype != prof_clutLab && ptype != prof_clutXYZ)
			error ("Input cLUT clipping above WP mode isn't applicable to a matrix profile");

		make_input_icc(ptype, iccver, verb, iquality, oquality, noisluts, noipluts, nooluts, nocied,
		               verify, autowpsc, clipovwp, iwpscale, doinb2a, doinextrap, clipprims,
		               inname, outname, icg, emis, spec, illum, &cust_illum, obType, custObserver,
		               smooth, avgdev, &xpi);

	/* Display or Projector */
	} else if (strcmp(icg->t[0].kdata[dti],"DISPLAY") == 0) {

		if (autowpsc)
			error ("Input auto WP scale mode isn't applicable to an output device");
		if (clipovwp)
			error ("Input cLUT clipping above WP mode isn't applicable to an output device");

		if (fwacomp)
			error ("FWA compensation isn't applicable to a display device");

		if (ptype == prof_default)
			ptype = prof_clutLab;		/* ?? or should it default to prof_shamat ?? */

		/* If a source gamut is provided for a Display, then a V2.4.0 profile will be created */
		make_output_icc(ptype, mtxtoo, iccver, verb, iquality, oquality,
		                noisluts, noipluts, nooluts, nocied, noptop, nostos,
		                gamdiag, verify, clipprims, iwpscale,
//		                bpo[1] >= 0.0 ? bpo : NULL,
		                NULL, inname, outname, icg,
		                spec, icxIT_none, NULL, illum, &cust_illum, obType, custObserver,
						0, smooth, avgdev, demph,
						gcompr, gexpr,
		                ipname[0] != '\000' ? ipname : NULL,
		                sgname[0] != '\000' ? sgname : NULL,
		                absnames,
						sepsat, &ivc_p, &ovc_p, ivc_e, ovc_e,
						&pgmi, &sgmi, &xpi);

	} else
		error ("Input file keyword DEVICE_CLASS has unknown value");

	icg->del(icg);		/* Clean up */

#ifdef DO_TIME			/* Time the operation */
	ttime = clock() - stime;
	printf("Exectution time = %f seconds\n",(double)ttime/(double)CLOCKS_PER_SEC);
#endif /* DO_TIME */

	return 0;
}










