
/* 
 * International Color Consortium Format Library (icclib)
 * Profile color lookup utility.
 *
 * Author:  Graeme W. Gill
 * Date:    1999/11/29
 * Version: 2.15
 *
 * Copyright 1998 - 2012 Graeme W. Gill
 *
 * This material is licensed with an "MIT" free use license:-
 * see the License4.txt file in this directory for licensing details.
 */

/* TTBD:
 *
 */

/*

	This file is intended to exercise the ability
	of the icc library to translate colors through a profile.
	It also serves as a concise example of how to do this.

 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "icc.h"

void error(char *fmt, ...), warning(char *fmt, ...);

void usage(void) {
	fprintf(stderr,"Translate colors through an ICC profile, V%s\n",ICCLIB_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill\n");
	fprintf(stderr,"usage: icclu [-v level] [-f func] [-i intent] [-o order] profile\n");
	fprintf(stderr," -v level      Verbosity level 0 - 2 (default = 1)\n");
	fprintf(stderr," -f function   f = forward, b = backwards, g = gamut, p = preview\n");
	fprintf(stderr," -i intent     p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"               s = saturation, a = absolute\n");
//	fprintf(stderr,"               P = absolute perceptual, S = absolute saturation\n");
	fprintf(stderr," -p oride      x = XYZ_PCS, l = Lab_PCS, y = Yxy,\n");
	fprintf(stderr," -o order      n = normal (priority: lut > matrix > monochrome)\n");
	fprintf(stderr,"               r = reverse (priority: monochrome > matrix > lut)\n");
	fprintf(stderr," -s scale      Scale device range 0.0 - scale rather than 0.0 - 1.0\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"    The colors to be translated should be fed into standard input,\n");
	fprintf(stderr,"    one input color per line, white space separated.\n");
	fprintf(stderr,"    A line starting with a # will be ignored.\n");
	fprintf(stderr,"    A line not starting with a number will terminate the program.\n");
	exit(1);
}

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char prof_name[500];
	icmFile *fp;
	icc *icco;
	int verb = 1;
	double scale = 0.0;		/* Device value scale factor */
	int rv = 0;
	int repYxy = 0;			/* Report Yxy */
	char buf[200];
	double oin[MAX_CHAN], in[MAX_CHAN], out[MAX_CHAN];

	icmLuBase *luo;
	icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
	int inn, outn;						/* Number of components */
	icmLuAlgType alg;					/* Type of lookup algorithm */

	/* Lookup parameters */
	icmLookupFunc     func   = icmFwd;				/* Default */
	icRenderingIntent intent = icmDefaultIntent;	/* Default */
	icColorSpaceSignature pcsor = icmSigDefaultData;	/* Default */
	icmLookupOrder    order  = icmLuOrdNorm;		/* Default */
	
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
				fa = nfa;
				if (na == NULL)
					verb = 2;
				else {
					if (na[0] == '0')
						verb = 0;
					else if (na[0] == '1')
						verb = 1;
					else if (na[0] == '2')
						verb = 2;
					else
						usage();
				}
			}

			/* function */
			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F') {
				fa = nfa;
				if (na == NULL) usage();
    			switch (na[0]) {
					case 'f':
					case 'F':
						func = icmFwd;
						break;
					case 'b':
					case 'B':
						func = icmBwd;
						break;
					case 'g':
					case 'G':
						func = icmGamut;
						break;
					case 'p':
					case 'P':
						func = icmPreview;
						break;
					default:
						usage();
				}
			}

			/* Intent */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage();
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
					/* Special function icclib intents */
					case 'P':
						intent = icmAbsolutePerceptual;
						break;
					case 'S':
						intent = icmAbsoluteSaturation;
						break;
					default:
						usage();
				}
			}

			/* Search order */
			else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage();
    			switch (na[0]) {
					case 'n':
					case 'N':
						order = icmLuOrdNorm;
						break;
					case 'r':
					case 'R':
						order = icmLuOrdRev;
						break;
					default:
						usage();
				}
			}

			/* PCS override */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage();
    			switch (na[0]) {
					case 'x':
					case 'X':
						pcsor = icSigXYZData;
						repYxy = 0;
						break;
					case 'l':
					case 'L':
						pcsor = icSigLabData;
						repYxy = 0;
						break;
					case 'y':
					case 'Y':
						pcsor = icSigXYZData;
						repYxy = 1;
						break;
					default:
						usage();
				}
			}

			/* Device scale */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				fa = nfa;
				if (na == NULL) usage();
				scale = atof(na);
				if (scale <= 0.0) usage();
			}

			else 
				usage();
		} else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(prof_name,argv[fa]);

	/* Open up the profile for reading */
	if ((fp = new_icmFileStd_name(prof_name,"r")) == NULL)
		error ("Can't open file '%s'",prof_name);

	if ((icco = new_icc()) == NULL)
		error ("Creation of ICC object failed");

	if ((rv = icco->read(icco,fp,0)) != 0)
		error ("%d, %s",rv,icco->err);

	if (icco->header->cmmId == str2tag("argl"))
		icco->allowclutPoints256 = 1;

	if (verb > 1) {
		icmFile *op;
		if ((op = new_icmFileStd_fp(stdout)) == NULL)
			error ("Can't open stdout");
		icco->header->dump(icco->header, op, 1);
		op->del(op);
	}

	/* Get a conversion object */
	if ((luo = icco->get_luobj(icco, func, intent, pcsor, order)) == NULL)
		error ("%d, %s",icco->errc, icco->err);

	/* Get details of conversion (Arguments may be NULL if info not needed) */
	luo->spaces(luo, &ins, &inn, &outs, &outn, &alg, NULL, NULL, NULL, NULL);

	if (repYxy) {	/* report Yxy rather than XYZ */
		if (ins == icSigXYZData)
			ins = icSigYxyData; 
		if (outs == icSigXYZData)
			outs = icSigYxyData; 
	}
		
	/* Process colors to translate */
	for (;;) {
		int i,j;
		char *bp, *nbp;

		/* Init default input values */
		for (i = 0; i < MAX_CHAN; i++) {
			in[i] = oin[i] = 0.0;
		}

		/* Read in the next line */
		if (fgets(buf, 200, stdin) == NULL)
			break;
		if (buf[0] == '#') {
			if (verb > 0)
				fprintf(stdout,"%s\n",buf);
			continue;
		}
		/* For each input number */
		for (nbp = buf, i = 0; i < MAX_CHAN; i++) {
			bp = nbp;
			in[i] = oin[i] = strtod(bp, &nbp);
			if (nbp == bp)
				break;			/* Failed */
		}
		if (i == 0)
			break;

		/* If device data and scale */
		if(scale > 0.0
    	 && ins != icSigXYZData
		 && ins != icSigLabData
		 && ins != icSigLuvData
		 && ins != icSigYCbCrData
		 && ins != icSigYxyData
		 && ins != icSigHsvData
		 && ins != icSigHlsData) {
			for (i = 0; i < MAX_CHAN; i++) {
				in[i] /= scale;
			}
		}

		if (repYxy && ins == icSigYxyData) {
			icmYxy2XYZ(in, in);
		}

		/* Do conversion */
		if ((rv = luo->lookup(luo, out, in)) > 1)
			error ("%d, %s",icco->errc,icco->err);

		/* Output the results */
		if (verb > 0) {
			for (j = 0; j < inn; j++) {
				if (j > 0)
					fprintf(stdout," %f",oin[j]);
				else
					fprintf(stdout,"%f",oin[j]);
			}
			printf(" [%s] -> %s -> ", icm2str(icmColorSpaceSignature, ins),
		                          icm2str(icmLuAlg, alg));
		}
		
		if (repYxy && outs == icSigYxyData) {
			icmXYZ2Yxy(out, out);
		}

		/* If device data and scale */
		if(scale > 0.0
    	 && outs != icSigXYZData
		 && outs != icSigLabData
		 && outs != icSigLuvData
		 && outs != icSigYCbCrData
		 && outs != icSigYxyData
		 && outs != icSigHsvData
		 && outs != icSigHlsData) {
			for (i = 0; i < MAX_CHAN; i++) {
				out[i] *= scale;
			}
		}

		for (j = 0; j < outn; j++) {
			if (j > 0)
				fprintf(stdout," %f",out[j]);
			else
				fprintf(stdout,"%f",out[j]);
		}
		if (verb > 0)
			printf(" [%s]", icm2str(icmColorSpaceSignature, outs));

		if (verb > 0 && rv != 0)
			fprintf(stdout," (clip)");

		fprintf(stdout,"\n");

	}

	/* Done with lookup object */
	luo->del(luo);

	icco->del(icco);

	fp->del(fp);

	return 0;
}


/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"icclu: Error - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit (-1);
}

void
warning(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"icclu: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}
