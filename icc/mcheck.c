
/* 
 * International Color Consortium Format Library (icclib)
 * Check the device chanel to PCS monotonicity.
 *
 * Author:  Graeme W. Gill
 * Date:    2000/12/11
 * Version: 2.15
 *
 * Copyright 2000 - 2012 Graeme W. Gill
 *
 * This material is licensed with an "MIT" free use license:-
 * see the License4.txt file in this directory for licensing details.
 */

/* TTBD:
 *
 * Make general device input, not just CMYK
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "vrml.h"
#include "icc.h"

void error(char *fmt, ...), warning(char *fmt, ...);

void usage(void) {
	fprintf(stderr,"Check device to PCS monotonicity of a CMYK ICC file, V%s\n",ICCLIB_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill\n");
	fprintf(stderr,"usage: mcheck [-v] [-w] infile\n");
	fprintf(stderr," -v        verbose\n");
	fprintf(stderr," -c        Check just Cyan monotonicity\n");
	fprintf(stderr," -m        Check just Magenta monotonicity\n");
	fprintf(stderr," -y        Check just Yellow monotonicity\n");
	fprintf(stderr," -k        Check just Black monotonicity\n");
	fprintf(stderr," -w        create %s visualisation\n",vrml_format());
	fprintf(stderr," -x        use %s axes\n",vrml_format());
	exit(1);
}

#define MGR 50			/* Maximum grid resolution handled */

int
main(
	int argc,
	char *argv[]
) {
	int fa,nfa;				/* argument we're looking at */
	int verb = 0;
	int cchan = -1;			/* default all */
	int dovrml = 0;
	int doaxes = 0;
	char in_name[500];
	char out_name[500], *xl;
	icmFile *rd_fp;
	icc *wr_icco, *rd_icco;		/* Keep object separate */
	int rv = 0;

	/* Check variables */
	icmLuBase *luo;
	icmLuLut *luluto;	/* Lookup xLut type object */
	int gres;			/* Grid resolution */
	icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
	int inn;							/* Number of input chanels */
	icmLuAlgType alg;
	vrml *wrl;
	int dx[4];			/* Device index mapping */
	int chan, cs, ce;
	
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

			/* Verbosity */
			if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			}
			/* VRML/X3D */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				dovrml = 1;
			}
			/* Cyan */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				cchan = 0;
			}
			/* Magenta */
			else if (argv[fa][1] == 'm' || argv[fa][1] == 'M') {
				cchan = 1;
			}
			/* Yellow */
			else if (argv[fa][1] == 'y' || argv[fa][1] == 'Y') {
				cchan = 2;
			}
			/* Black */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				cchan = 3;
			}
			else if (argv[fa][1] == 'x') {
				doaxes = 1;
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
	strcpy(in_name,argv[fa]);

	strcpy(out_name, in_name);
	if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
		xl = out_name + strlen(out_name);
	xl[0] = '\000';			/* Remove extension */

	/* Open up the file for reading */
	if ((rd_fp = new_icmFileStd_name(in_name,"r")) == NULL)
		error ("Read: Can't open file '%s'",in_name);

	if ((rd_icco = new_icc()) == NULL)
		error ("Read: Creation of ICC object failed");

	/* Read the header and tag list */
	if ((rv = rd_icco->read(rd_icco,rd_fp,0)) != 0)
		error ("Read: %d, %s",rv,rd_icco->err);

	/* Get a Device to PCS conversion object */
	if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, icRelativeColorimetric, icSigLabData, icmLuOrdNorm)) == NULL) {
		if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, icmDefaultIntent, icSigLabData, icmLuOrdNorm)) == NULL)
			error ("%d, %s",rd_icco->errc, rd_icco->err);
	}
	/* Get details of conversion */
	luo->spaces(luo, &ins, &inn, &outs, NULL, &alg, NULL, NULL, NULL, NULL);

	if (alg != icmLutType) {
		error("Expecting Lut based profile");
	}

	if (ins != icSigCmykData) {
		error("Expecting CMYK device");
	}

	if (outs != icSigLabData) {
		error("Expecting Lab PCS");
	}

	luluto = (icmLuLut *)luo;	/* Lookup xLut type object */

	gres = luluto->lut->clutPoints;
	if (gres > MGR) {
		error("Can't handle grid resolution greater than %d\n",MGR);
	}

	if (dovrml) {
		wrl = new_vrml(out_name, doaxes, vrml_lab);
		if (wrl == NULL)
			error("new_vrml for '%s%s' failed",out_name,vrml_ext());
		wrl->start_line_set(wrl, 0);
	}

	/* For all the device chanels chosen */
	if (cchan < 0) {
		cs = 0;
		ce = inn;
	} else {
		cs = cchan;
		ce = cs + 1;
	}
	for (chan = cs; chan < ce; chan++) {

		/* Check the monotonicity of the output for a given device input */
		int co[4];
		if (chan == 0) {
			dx[0] = 1;
			dx[1] = 2;
			dx[2] = 3;
			dx[3] = 0;		/* Cyan is variable */
		} else if (chan == 1) {
			dx[0] = 0;
			dx[1] = 2;
			dx[2] = 3;
			dx[3] = 1;		/* Magenta is variable */
		} else if (chan == 2) {
			dx[0] = 0;
			dx[1] = 1;
			dx[2] = 3;
			dx[3] = 2;		/* Yellow is variable */
		} else if (chan == 3) {
			dx[0] = 0;
			dx[1] = 1;
			dx[2] = 2;
			dx[3] = 3;		/* Black is variable */
		}

		/* Itterate throught the CMY clut grid points */
		for (co[0] = 0; co[0] < gres; co[0]++) {
			for (co[1] = 0; co[1] < gres; co[1]++) {
				for (co[2] = 0; co[2] < gres; co[2]++) {
					int j, k, ck, nm;
					double dev[MGR][4];
					double pcs[MGR][3];
					double apcs[3], ss;

					/* Run up the variable axis */
					for (ck = 0; ck < gres; ck++) {

						dev[ck][dx[0]] = co[0]/(gres-1.0);
						dev[ck][dx[1]] = co[1]/(gres-1.0);
						dev[ck][dx[2]] = co[2]/(gres-1.0);
						dev[ck][dx[3]] = ck/(gres-1.0);

						/* Device to PCS */
						if ((rv = luluto->clut(luluto, pcs[ck], dev[ck])) > 1)
							error ("%d, %s",rd_icco->errc,rd_icco->err);

//						if (dovrml)
//							wrl->add_vertex(wrl, 0, pcs[ck]);
					}

					/* Compute average vector direction */
					for (ss = 0.0, k = 0; k < 3; k++) {
						double tt;
						tt = pcs[gres-1][k] - pcs[0][k];
						ss += tt * tt;
						apcs[k] = tt;
					}
					for (k = 0; k < 3; k++)
						apcs[k] /= ss;

					/* Now compute the dot product for each vector, */
					/* and check for reversals. */
					j = 0;
//printf("Checking          CMYK %f %f %f %f Lab %f %f %f\n",
//       dev[j][0], dev[j][1], dev[j][2], dev[j][3],
//       pcs[j][0], pcs[j][1], pcs[j][2]);
					for (nm = 0, j = 1; j < gres; j++) {
						for (ss = 0.0, k = 0; k < 3; k++)	/* Dot product */
							ss += (pcs[j][k] - pcs[j-1][k]) * apcs[k];

//printf("Checking %f CMYK %f %f %f %f Lab %f %f %f\n",
//       ss, dev[j][0], dev[j][1], dev[j][2], dev[j][3],
//       pcs[j][0], pcs[j][1], pcs[j][2]);

						if (ss <= 0.0) {
							nm = 1;
							printf("NonMon %f at CMYK %f %f %f %f Lab %f %f %f\n",
							       ss, dev[j][0], dev[j][1], dev[j][2], dev[j][3],
							       pcs[j][0], pcs[j][1], pcs[j][2]);
						}
					}
//printf("\n");

					/* Display just the non mono threads */
					if (nm && dovrml) {
						for (j = 0; j < gres; j++)
							wrl->add_vertex(wrl, 0, pcs[j]);
					}
					if (verb) {
						printf("."); fflush(stdout);
					}
				}
			}
		}
	}

	if (dovrml) {
		wrl->make_lines(wrl, 0, gres);
		wrl->del(wrl);
	}

	/* Done with lookup object */
	luo->del(luo);

	rd_icco->del(rd_icco);
	rd_fp->del(rd_fp);

	return 0;
}

/* ------------------------------------------------ */
/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"icctest: Error - ");
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

	fprintf(stderr,"icctest: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}
