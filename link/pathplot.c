
/* 
 * Argyll Color Management System
 *
 * Plot the L in vs. L out for a path through the
 * input device space, for a given link.
 *
 * We are assuming RGB device input space !!!!
 *
 * Author:  Graeme W. Gill
 * Date:    2002/1/28
 * Version: 1.00
 *
 * Copyright 2002 Graeme W. Gill
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
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
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "icc.h"
#include "xicc.h"
#include "plot.h"
#include "ui.h"

#define PRES 100

static double start[13][4]   = {
	{ 0.0, 0.0, 0.0, 0.0},
	{ 0.0, 0.0, 0.0, 0.0},
	{ 0.0, 0.0, 0.0, 0.0},
	{ 0.0, 0.0, 0.0, 0.0},
	{ 0.0, 0.0, 0.0, 0.0},
	{ 0.0, 0.0, 0.0, 0.0},

	{ 0.0, 0.0, 0.0, 0.0},

	{ 1.0, 0.0, 0.0, 0.0},
	{ 0.0, 1.0, 0.0, 0.0},
	{ 0.0, 0.0, 1.0, 0.0},
	{ 0.0, 1.0, 1.0, 0.0},
	{ 1.0, 0.0, 1.0, 0.0},
	{ 1.0, 1.0, 0.0, 0.0}
};

static double end[13][4]  = {
	{ 1.0, 0.0, 0.0, 0.0},
	{ 0.0, 1.0, 0.0, 0.0},
	{ 0.0, 0.0, 1.0, 0.0},
	{ 0.0, 1.0, 1.0, 0.0},
	{ 1.0, 0.0, 1.0, 0.0},
	{ 1.0, 1.0, 0.0, 0.0},

	{ 1.0, 1.0, 1.0, 0.0},

	{ 1.0, 1.0, 1.0, 0.0},
	{ 1.0, 1.0, 1.0, 0.0},
	{ 1.0, 1.0, 1.0, 0.0},
	{ 1.0, 1.0, 1.0, 0.0},
	{ 1.0, 1.0, 1.0, 0.0},
	{ 1.0, 1.0, 1.0, 0.0}
};

static char *name[13] = {
		"Black - Red", "Black - Green", "Black - Blue", 
		"Black - Cyan", "Black - Magenta", "Black - Yellow", 
		"Grey",
		"Red - White", "Green - White", "Blue - White", 
		"Cyan - White", "Magenta - White", "Yellow - White", 
	};

#define USE_JAB		/* Assume CRT in (2), print out (0) */

/* ---------------------------------------- */

void usage(void) {
	fprintf(stderr,"Plot device space path L in/out curve from an ICC link file, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill\n");
	fprintf(stderr,"usage: pathplot inprof linkprof outprof\n");
	fprintf(stderr," -v        verbose\n");
	exit(1);
}

int
main(
	int argc,
	char *argv[]
) {
	int fa, nfa;				/* argument we're looking at */
	int verb = 0;
	char in_name[100];
	char link_name[100];
	char out_name[100];
	icmFile *in_fp, *link_fp, *out_fp;
	icc *in_icco, *link_icco, *out_icco;
	xicc *in_xicco, *out_xicco;
	icmLuBase *link_lu;
	icxLuBase *in_lu, *out_lu;
	icColorSpaceSignature pcsor;	/* PCS to use */
	icxViewCond ivc[1], ovc[1];
	int rv = 0;

	error_program = argv[0];

	if (argc < 4)
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

			/* Verbosity */
			if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			}
			else if (argv[fa][1] == '?')
				usage();
			else 
				usage();
		}
		else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(in_name,argv[fa++]);

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(link_name,argv[fa++]);

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(out_name,argv[fa++]);

	/* Open up the files for reading */
	if ((in_fp = new_icmFileStd_name(in_name,"r")) == NULL)
		error ("Read: Can't open file '%s'",in_name);

	if ((in_icco = new_icc()) == NULL)
		error ("Read: Creation of ICC object failed");

	if ((rv = in_icco->read(in_icco,in_fp,0)) != 0)
		error ("Read: %d, %s",rv,in_icco->err);

	if ((in_xicco = new_xicc(in_icco)) == NULL)
		error ("Creation of input profile xicc failed");


	if ((link_fp = new_icmFileStd_name(link_name,"r")) == NULL)
		error ("Read: Can't open file '%s'",link_name);

	if ((link_icco = new_icc()) == NULL)
		error ("Read: Creation of ICC object failed");

	if ((rv = link_icco->read(link_icco,link_fp,0)) != 0)
		error ("Read: %d, %s",rv,link_icco->err);


	if ((out_fp = new_icmFileStd_name(out_name,"r")) == NULL)
		error ("Read: Can't open file '%s'",out_name);

	if ((out_icco = new_icc()) == NULL)
		error ("Read: Creation of ICC object failed");

	if ((rv = out_icco->read(out_icco,out_fp,0)) != 0)
		error ("Read: %d, %s",rv,out_icco->err);

	if ((out_xicco = new_xicc(out_icco)) == NULL)
		error ("Creation of output profile xicc failed");


#ifdef USE_JAB
	pcsor = icxSigJabData;		/* Use CIECAM as PCS */

	if (xicc_enum_viewcond(in_xicco, ivc, -2, "mt", 0, NULL) == -999)	/* Set input at monitor in typical */
		error ("%d, %s",in_xicco->errc, in_xicco->err);

	if (xicc_enum_viewcond(out_xicco, ovc, -2, "pp", 0, NULL) == -999)	/* Set output at practical reflection print */
		error ("%d, %s",out_xicco->errc, out_xicco->err);

#else
	pcsor = icSigLabData;		/* Default use Lab as PCS */
#endif	/* !USE_JAB */

	/* Device to PCS conversion object */
	if ((in_lu = in_xicco->get_luobj(in_xicco, ICX_CLIP_NEAREST, icmFwd, icAbsoluteColorimetric, pcsor, icmLuOrdNorm, ivc, NULL)) == NULL) {
		if ((in_lu = in_xicco->get_luobj(in_xicco, ICX_CLIP_NEAREST, icmBwd, icmDefaultIntent, pcsor, icmLuOrdNorm, ivc, NULL)) == NULL)
			error ("%d, %s",in_xicco->errc, in_xicco->err);
	}

	/* Get a Device to Device conversion object */
	if ((link_lu = link_icco->get_luobj(link_icco, icmFwd, icmDefaultIntent, pcsor, icmLuOrdNorm)) == NULL)
		error ("%d, %s",link_icco->errc, link_icco->err);

	/* Get a Device to PCS conversion object */
	if ((out_lu = out_xicco->get_luobj(out_xicco, ICX_CLIP_NEAREST, icmFwd, icAbsoluteColorimetric, pcsor, icmLuOrdNorm, ovc, NULL)) == NULL) {
		if ((out_lu = out_xicco->get_luobj(out_xicco, ICX_CLIP_NEAREST, icmFwd, icmDefaultIntent, pcsor, icmLuOrdNorm, ovc, NULL)) == NULL)
			error ("%d, %s",out_xicco->errc, out_xicco->err);
	}

	{
		double xx[PRES], yy[PRES];
		int k, i, j;
		double tt[10], tt2[10];

		for (k = 0; k < 13; k++) {
			printf("Doing %s\n",name[k]);

			for (i = 0; i < PRES; i++) {
				double frac = i/(PRES-1.0);

				for (j = 0; j < 4; j++)
					tt[j] = (end[k][j] - start[k][j]) * frac + start[k][j];

				if (verb)
					printf(" %f %f %f -> ",tt[0],tt[1],tt[2]);

				/* input device space to PCS */
				if ((rv = in_lu->lookup(in_lu, tt2, tt)) > 1)
					error ("%d, %s",in_icco->errc,in_icco->err);

				xx[i] = tt2[0];		/* L value */

				/* input device space to output device space */
				if ((rv = link_lu->lookup(link_lu, tt, tt)) > 1)
					error ("%d, %s",link_icco->errc,link_icco->err);

				/* output device space to PCS */
				if ((rv = out_lu->lookup(out_lu, tt, tt)) > 1)
					error ("%d, %s",out_icco->errc,out_icco->err);

				yy[i] = tt[0];		/* L value */

				if (verb)
					printf("%f %f %f\n",tt[0],tt[1],tt[2]);
			}

			if (do_plot(xx,yy,NULL,NULL,PRES) < 0)
				error("do_plot returned -1!\n");
		}

	}

	/* Done with lookup objects */
	in_lu->del(in_lu);
	link_lu->del(link_lu);
	out_lu->del(out_lu);

	in_icco->del(in_icco);
	in_fp->del(in_fp);
	link_icco->del(link_icco);
	link_fp->del(link_fp);
	out_icco->del(out_icco);
	out_fp->del(out_fp);

	return 0;
}

