
/* 
 * Argyll Color Correction System
 * Inverse profile checker.
 *
 * Author: Graeme W. Gill
 * Date:   1999/11/29
 *
 * Copyright 1999 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes checks the round trip errors of 
 * the colorimetric forward and inverse profile direction
 * of an ICC profile.
 * (Was called icc/fbtest.c)
 */


/* TTBD:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "icc.h"
#include "xicc.h"
#include "vrml.h"

/* Resolution of the sampling modes */
#define TRES 11
#define HTRES 27
#define UHTRES 61

/* ------------------------------------------------------- */
/* Macros for an di or fdi dimensional counter */
/* Declare the counter name nn, dimensions di, & count */

#define DCOUNT(nn, di, start, reset, count) 				\
	int nn[MAX_CHAN];	/* counter value */					\
	int nn##_di = (di);		/* Number of dimensions */		\
	int nn##_stt = (start);	/* start count value */			\
	int nn##_rst = (reset);	/* reset on carry value */		\
	int nn##_res = (count); /* last count +1 */				\
	int nn##_e				/* dimension index */

/* Set the counter value to 0 */
#define DC_INIT(nn) 								\
{													\
	for (nn##_e = 0; nn##_e < nn##_di; nn##_e++)	\
		nn[nn##_e] = nn##_stt;						\
	nn##_e = 0;										\
}

/* Increment the counter value */
#define DC_INC(nn)									\
{													\
	for (nn##_e = 0; nn##_e < nn##_di; nn##_e++) {	\
		nn[nn##_e]++;								\
		if (nn[nn##_e] < nn##_res)					\
			break;	/* No carry */					\
		nn[nn##_e] = nn##_rst;						\
	}												\
}

/* After increment, expression is TRUE if counter is done */
#define DC_DONE(nn)									\
	(nn##_e >= nn##_di)
	
/* ---------------------------------------- */

void usage(void) {
	fprintf(stderr,"Check fwd to bwd relative transfer of an ICC file, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill\n");
	fprintf(stderr,"usage: invprofcheck [-] profile.icm\n");
	fprintf(stderr," -v [level]   verbosity level (default 1), 2 to print each DE\n");
	fprintf(stderr," -l limit     set total ink limit (estimate by default)\n");
	fprintf(stderr," -L klimit    set black channel ink limit (estimate by default)\n");
	fprintf(stderr," -i intent      a = absolute, r = relative colorimetric (def.)\n");
	fprintf(stderr,"                p = perceptual, s = saturation\n");
	fprintf(stderr," -h           high res test (%d)\n",HTRES);
	fprintf(stderr," -u           Ultra high res test (%d)\n",UHTRES);
	fprintf(stderr," -R res       Specific grid resolution\n");
	fprintf(stderr," -I           Do bwd to fwd check\n");
	fprintf(stderr," -c           Show CIE94 delta E values\n");
	fprintf(stderr," -k           Show CIEDE2000 delta E values\n");
	fprintf(stderr," -w           create %s visualisation (profile%s)\n",vrml_format(),vrml_ext());
	fprintf(stderr," -x           Use %s axes\n",vrml_format());
	fprintf(stderr," -e           Color vectors acording to delta E\n");
	fprintf(stderr," profile.icm  Profile to check\n");
	exit(1);
}

static void DE2RGB(double *out, double in);

#if defined(__IBMC__) && defined(_M_IX86)
void bug_workaround(int *co) { };			/* Workaround optimiser bug */
#endif

int
main(
	int argc,
	char *argv[]
) {
	int fa,nfa;						/* argument we're looking at */
	int verb = 0;
	int cie94 = 0;
	int cie2k = 0;
	int dovrml = 0;
	int doaxes = 0;
	int dodecol = 0;
	char in_name[MAXNAMEL+1];
	char out_name[MAXNAMEL+1], *xl;		/* VRML/X3D name */
	icmFile *rd_fp;
	icc *icco;
	int rv = 0;
	int inv = 0;
	int tres = TRES;
	double tlimit = -1.0;
	double klimit = -1.0;
	icRenderingIntent intent = icRelativeColorimetric;	/* Default */
	vrml *wrl = NULL;

	error_program = "invprofcheck";

	if (argc < 2)
		usage();

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
				usage();

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
				if (na != NULL && isdigit(na[0])) {
					verb = atoi(na);
				}
			}

			/* Resolution */
			else if (argv[fa][1] == 'h' || argv[fa][1] == 'H') {
				tres = HTRES;
			}

			/* Resolution */
			else if (argv[fa][1] == 'u' || argv[fa][1] == 'U') {
				tres = UHTRES;
			}

			/* Resolution */
			else if (argv[fa][1] == 'R') {
				int res;
				if (na == NULL) usage();
				fa = nfa;
				res = atoi(na);
				if (res < 2 || res > 500)
					usage();
				tres = res;
			}

			/* Inverse */
			else if (argv[fa][1] == 'I') {
				inv = 1;
			}

			else if (argv[fa][1] == 'l') {
				int limit;
				if (na == NULL) usage();
				fa = nfa;
				limit = atoi(na);
				if (limit < 1)
					limit = 1;
				tlimit = limit/100.0;
			}

			else if (argv[fa][1] == 'L') {
				int limit;
				if (na == NULL) usage();
				fa = nfa;
				limit = atoi(na);
				if (limit < 1)
					limit = 1;
				klimit = limit/100.0;
			}

			/* Intent */
			else if (argv[fa][1] == 'i') {
				if (na == NULL) usage();
				fa = nfa;
    			switch (na[0]) {
					case 'p':
						intent = icPerceptual;
						break;
					case 'r':
						intent = icRelativeColorimetric;
						break;
					case 's':
						intent = icSaturation;
						break;
					case 'a':
						intent = icAbsoluteColorimetric;
						break;
					default:
						usage();
				}
			}

			/* VRML/X3D */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W')
				dovrml = 1;

			/* Axes */
			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X')
				doaxes = 1;

			/* Delta E coloring */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E')
				dodecol = 1;

			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				cie94 = 1;
				cie2k = 0;
			}

			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				cie94 = 0;
				cie2k = 1;
			}

			else 
				usage();
		}
		else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(in_name,argv[fa++],MAXNAMEL); in_name[MAXNAMEL] = '\000';

	strncpy(out_name,in_name,MAXNAMEL-4); out_name[MAXNAMEL-4] = '\000';
	if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
		xl = out_name + strlen(out_name);
	xl[0] = '\000';								/* Remove extension */

	/* Open up the file for reading */
	if ((rd_fp = new_icmFileStd_name(in_name,"r")) == NULL)
		error ("Read: Can't open file '%s'",in_name);

	if ((icco = new_icc()) == NULL)
		error ("Read: Creation of ICC object failed");

	/* Read the header and tag list */
	if ((rv = icco->read(icco,rd_fp,0)) != 0)
		error ("Read: %d, %s",rv,icco->err);

	/* Check the forward lookup against the bwd function */
	{
		xcal *cal = NULL;                   /* Device calibration curves */
		icColorSpaceSignature ins, outs;	/* Type of input and output spaces of fwd */
		int inn, outn;						/* Channels of fwd conversion */
		int kch;							/* Black channel, -1 if not known/applicable */
		icmLuBase *luo1, *luo2;
		double merr = 0.0;		/* Max */
		double aerr = 0.0;		/* Avg */
		double rerr = 0.0;		/* RMS */
		double nsamps = 0.0;

		/* Get a Device to PCS conversion object */
		if ((luo1 = icco->get_luobj(icco, icmFwd, intent, icSigLabData, icmLuOrdNorm)) == NULL) {
			if ((luo1 = icco->get_luobj(icco, icmFwd, icmDefaultIntent, icSigLabData, icmLuOrdNorm)) == NULL)
				error ("%d, %s",icco->errc, icco->err);
		}
		/* Get details of conversion */
		luo1->spaces(luo1, &ins, &inn, &outs, &outn, NULL, NULL, NULL, NULL, NULL);

		/* Get a PCS to Device conversion object */
		if ((luo2 = icco->get_luobj(icco, icmBwd, intent, icSigLabData, icmLuOrdNorm)) == NULL) {
			if ((luo2 = icco->get_luobj(icco, icmBwd, icmDefaultIntent, icSigLabData, icmLuOrdNorm)) == NULL)
				error ("%d, %s",icco->errc, icco->err);
		}

		if (dovrml) {
			wrl = new_vrml(out_name, doaxes, vrml_lab);
			wrl->start_line_set(wrl, 0);
		}
		
		/* Grab any device calibration curves */
		cal = xiccReadCalTag(icco);

		kch = icxGuessBlackChan(icco);

		/* Set the default ink limits if not set by user */
		if (tlimit < 0.0 || klimit < 0.0) {
			double max[MAX_CHAN], total;

			total = icco->get_tac(icco, max, cal != NULL ? xiccCalCallback : NULL, (void *)cal);

			if (tlimit < 0.0)
				tlimit = total;

			if (klimit < 0.0 && kch >= 0)
				klimit = max[kch];
		}

		if (verb) {
			printf("Grid resolution is %d\n",tres);
			if (tlimit >= 0.0)
				printf("Input total ink limit assumed is %3.1f%%\n",100.0 * tlimit);
			if (klimit >= 0.0)
				printf("Input black ink limit assumed is %3.1f%%\n",100.0 * klimit);
		}

		/* Device -> PCS -> Device */
		if (!inv) {
			double dev[MAX_CHAN], cdev[MAX_CHAN], pcsin[3], devout[MAX_CHAN], pcsout[3];
			DCOUNT(co, inn, 0, 0, tres);		/* Multi-D counter */
	
			/* Go through the chosen device grid */
			DC_INIT(co)
			for (; !DC_DONE(co);) {
				int n, rv1, rv2;
				double sum;
				double de;

				/* Check the (possibly calibrated) device values */
				/* end reject any over the limits. */
				for (sum = 0, n = 0; n < inn; n++) {
					cdev[n] = dev[n] = co[n]/(tres-1.0);
					sum += cdev[n];
				}
				if (cal != NULL) {
					cal->interp(cal, cdev, dev);
					for (sum = 0, n = 0; n < inn; n++)
						sum += cdev[n];
				}

				if ((tlimit > 0.0 && sum > tlimit)
				 || (klimit > 0.0 && kch >= 0 && cdev[kch] > klimit)) {
					DC_INC(co);
					continue;
				}

				/* Generate the in-gamut PCS test point */
				/* by converting device to pcsin */
				if ((rv1 = luo1->lookup(luo1, pcsin, dev)) > 1)
					error ("%d, %s",icco->errc,icco->err);

				/* Now do the check */
				/* PCS -> Device */
				if ((rv2 = luo2->lookup(luo2, devout, pcsin)) > 1)
					error ("%d, %s",icco->errc,icco->err);

				/* Device to PCS */
				if ((rv2 = luo1->lookup(luo1, pcsout, devout)) > 1)
					error ("%d, %s",icco->errc,icco->err);

				/* Delta E */
				if (dovrml) {
					int ix[2];

					/* Add the verticies */
					ix[0] = wrl->add_vertex(wrl, 0, pcsin);
					ix[1] = wrl->add_vertex(wrl, 0, pcsout);

					/* Add the line */
					if (dodecol) {		/* Lines with color determined by length */
						double rgb[3];
						DE2RGB(rgb, icmNorm33(pcsin, pcsout));
						wrl->add_col_line(wrl, 0, ix, rgb);
	
					} else {	/* Natural color */
						wrl->add_line(wrl, 0, ix);
					}
				}
	
				/* Check the result */
				if (cie2k)
					de = icmCIE2K(pcsout, pcsin);
				else if (cie94)
					de = icmCIE94(pcsout, pcsin);
				else
					de = icmLabDE(pcsout, pcsin);
	
				aerr += de;
				rerr += de * de;
				if (de > merr)
					merr = de;
				nsamps++;

				if (verb > 1) {
					printf("[%f] %f %f %f -> ",de, pcsin[0], pcsin[1], pcsin[2]);
					for (n = 0; n < inn; n++)
						printf("%f ",devout[n]);
					printf("-> %f %f %f\n",pcsout[0], pcsout[1], pcsout[2]);
				}

				DC_INC(co);
			}

		/* PCS -> Device -> PCS */
		} else {
			double dev[MAX_CHAN], cdev[MAX_CHAN], pcsin[3], devout[MAX_CHAN], pcsout[3];
			DCOUNT(co, 3, 0, 0, tres);		/* Multi-D counter */
	
			/* Go through the chosen Lab grid */
			DC_INIT(co)
			for (; !DC_DONE(co);) {
				int n, rv1, rv2;
				double sum;
				double de;

				pcsin[0] = 100.0 * co[0]/(tres-1.0);
				pcsin[1] = (127.0 * 2.0 * co[1]/(tres-1.0)) - 127.0;
				pcsin[2] = (127.0 * 2.0 * co[2]/(tres-1.0)) - 127.0;

				/* PCS -> Device */
				if ((rv2 = luo2->lookup(luo2, devout, pcsin)) > 1)
					error ("%d, %s",icco->errc,icco->err);

				/* Device to PCS */
				if ((rv2 = luo1->lookup(luo1, pcsout, devout)) > 1)
					error ("%d, %s",icco->errc,icco->err);

				/* Delta E */
				if (dovrml) {
//				if (fabs(pcsin[0] - 5.0) < 0.1 && dovrml) {
					int ix[2];

					/* Add the verticies */
					ix[0] = wrl->add_vertex(wrl, 0, pcsin);
					ix[1] = wrl->add_vertex(wrl, 0, pcsout);

					/* Add the line */
					if (dodecol) {		/* Lines with color determined by length */
						double rgb[3];
						DE2RGB(rgb, icmNorm33(pcsin, pcsout));
						wrl->add_col_line(wrl, 0, ix, rgb);
	
					} else {	/* Natural color */
						wrl->add_line(wrl, 0, ix);
					}
				}
	
				/* Check the result */
				if (cie2k)
					de = icmCIE2K(pcsout, pcsin);
				else if (cie94)
					de = icmCIE94(pcsout, pcsin);
				else
					de = icmLabDE(pcsout, pcsin);
	
				aerr += de;
				rerr += de * de;
				if (de > merr)
					merr = de;
				nsamps++;

				if (verb > 1) {
					printf("[%f] %f %f %f -> ",de, pcsin[0], pcsin[1], pcsin[2]);
					for (n = 0; n < inn; n++)
						printf("%f ",devout[n]);
					printf("-> %f %f %f\n",pcsout[0], pcsout[1], pcsout[2]);
				}

				DC_INC(co);
			}
		}
		if (dovrml) {
			wrl->make_lines_vc(wrl, 0, 0.0);
			wrl->del(wrl);
		}

		printf("Profile check complete, errors%s: max. = %f, avg. = %f, RMS = %f\n",
            cie2k ? "(CIEDE2000)" : cie94 ? " (CIE94)" : "", merr, aerr/nsamps, sqrt(rerr/nsamps));

		/* Done with lookup object */
		luo1->del(luo1);
		luo2->del(luo2);
	}

	icco->del(icco);
	rd_fp->del(rd_fp);

	return 0;
}


/* ------------------------------------------------ */
/* Convert a delta E value into a signal color: */
static void DE2RGB(double *out, double in) {
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

//printf("~1 input de = %f\n",in);

	/* Locate the range we're in */
	if (in > range[0].de) {
		out[0] = range[0].r;
		out[1] = range[0].g;
		out[2] = range[0].b;
//printf("~1 too big\n");
	} else {
		for (i = 0; i < 5; i++) {
			if (in <= range[i].de && in >= range[i+1].de)
				break;
		}
		bl = (in - range[i+1].de)/(range[i].de - range[i+1].de);
//printf("~1 located at ix %d, bl = %f\n",i,bl);
		out[0] = bl * range[i].r + (1.0 - bl) * range[i+1].r;
		out[1] = bl * range[i].g + (1.0 - bl) * range[i+1].g;
		out[2] = bl * range[i].b + (1.0 - bl) * range[i+1].b;
	}
//printf("~1 returning rgb %f %f %f\n",out[0],out[1],out[2]);
}


