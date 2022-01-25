
/* 
 * Argyll Color Management System
 * Simple CMYK profile generator.
 *
 * Author: Graeme W. Gill
 * Date:   9/11/96
 *
 * Copyright 1996, 2002 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This program generates a simple mathematical profile for a CMYK device. */
/* It is intended for use in bootstrapping the test chart generation. */

#define VERSION "1.1"

#undef DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include "../cgats/cgats.h"
#include "../numlib/numlib.h"
#include "icc.h"

/* A color structure */
/* This holds the test patch results */
typedef struct {
	double c,m,y,k;
	double bc[8];			/* Gamma corrected blend coefficients of cmyk */
	double Lab[3];
	double err;				/* Delta E squared */
} col;

/* Structure to hold data for optimization function */
struct _edatas {
	col *cols;			/* Pointer to table of patch results */
	int npat;			/* Number of patches */
	int xyzi;			/* current xyz index */
	double gam[4];		/* Gamma values */
	double k[3][8];	/* Primary combination values */
	}; typedef struct _edatas edatas;

/* Definition of the power optimization function handed to powell() */
/* This function is for optimising the gamma function */
double efunc1(void *edata, double p[]) {
	edatas *ed = (edatas *)edata;
	double rv;
	col *cp;
	for (rv = 0.0, cp = &ed->cols[ed->npat-1]; cp >= &ed->cols[0]; cp--) {
		double cc,mm,yy,kk;
		double nc,nm,ny,nk;
		double XYZ[3], Lab[3];
		int j;

		/* Apply gamma correction to each input */
		for (j = 0; j < 4; j++) {
			if (p[j] < 0.2)
				p[j] = 0.2;
			else if (p[j] > 5.0)
				p[j] = 5.0;
		}
		cc = pow(cp->c, p[0]);
		nc = 1.0 - cc;
		mm = pow(cp->m, p[1]);
		nm = 1.0 - mm;
		yy = pow(cp->y, p[2]);
		ny = 1.0 - yy;
		kk = pow(cp->k, p[3]);
		nk = 1.0 - kk;

		/* Then interpolate between all combinations of primaries. */
		/* plus one that stands for all that are close to black */
		for (j = 0; j < 3; j++) {
			XYZ[j] = nc * nm * ny * nk * ed->k[j][0]
			       + nc * nm * yy * nk * ed->k[j][1]
			       + nc * mm * ny * nk * ed->k[j][2]
			       + nc * mm * yy * nk * ed->k[j][3]
			       + cc * nm * ny * nk * ed->k[j][4]
			       + cc * nm * yy * nk * ed->k[j][5]
			       + cc * mm * ny * nk * ed->k[j][6]
			       + (cc * mm * yy * nk + kk) * ed->k[j][7];
		}
		icmXYZ2Lab(&icmD50, Lab, XYZ);
		rv += cp->err = icmLabDEsq(Lab, cp->Lab);
	}
printf("Efunc1 returning %f\n",rv);
	return rv;
}

/* Definition of the primary coefficient optimization function handed to powell() */
/* This function is for optimising the primary values */
double efunc2(void *edata, double p[]) {
	edatas *ed = (edatas *)edata;
	int j, os = ed->xyzi;
	double tt, rv;
	col *cp;

	rv = 0.0;

	for (j = 0; j < 8; j++)  {
		if (p[j] < 0.0) {			/* Protect against silly values */
			p[j] = 0.0;
			rv += 1000.0;
		}
		else if (p[j] > 1.5) {
			p[j] = 1.5;
			rv += 1000.0;
		}
		ed->k[os][j] = p[j];	/* Load into current */
	}

	/* Compute error */
	for (cp = &ed->cols[ed->npat-1]; cp >= &ed->cols[0]; cp--) {
		double XYZ[3], Lab[3];

		for (os = 0; os < 3; os++) { 

			/* Interpolate between all combinations of primaries. */
			for (tt = 0.0, j = 0; j < 8; j++) 
				tt += cp->bc[j] * ed->k[os][j];
			XYZ[os] = tt;
		}
		icmXYZ2Lab(&icmD50, Lab, XYZ);
		rv += cp->err = icmLabDEsq(Lab, cp->Lab);
	}
printf("Efunc2 returning %f\n",rv);
	return rv;
}

/* Calculate blend coefficients */
void calc_bc(edatas *ed) {
	int j;
	col *cp;
	for (cp = &ed->cols[ed->npat-1]; cp >= &ed->cols[0]; cp--) {
		double cc,mm,yy,kk;
		double nc,nm,ny,nk;

		/* Apply gamma correction to each input, and calculate complement */
		for (j = 0; j < 4; j++) {
			if (ed->gam[j] < 0.2)
				ed->gam[j] = 0.2;
			else if (ed->gam[j] > 5.0)
				ed->gam[j] = 2.0;
		}
		cc = pow(cp->c, ed->gam[0]);
		nc = 1.0 - cc;
		mm = pow(cp->m, ed->gam[1]);
		nm = 1.0 - mm;
		yy = pow(cp->y, ed->gam[2]);
		ny = 1.0 - yy;
		kk = pow(cp->k, ed->gam[3]);
		nk = 1.0 - kk;

		/* Go through all 8 combinations */
		cp->bc[0] = nc * nm * ny * nk;
		cp->bc[1] = nc * nm * yy * nk;
		cp->bc[2] = nc * mm * ny * nk;
		cp->bc[3] = nc * mm * yy * nk;
		cp->bc[4] = cc * nm * ny * nk;
		cp->bc[5] = cc * nm * yy * nk;
		cp->bc[6] = cc * mm * ny * nk;
		cp->bc[7] = cc * mm * yy * nk + kk;
	}
}

void usage(void);

int main(int argc, char *argv[])
{
	int i,j,k;
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	static char inname[200] = { 0 };		/* Input cgats file base name */
	static char outname[200] = { 0 };		/* Output cgats file base name */
	cgats *icg;			/* input cgats structure */
	cgats *ocg;			/* output cgats structure */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	int ti;				/* Temporary index */
	edatas ed;			/* Optimising function data structure */
	double resid[4];
	double presid,dresid;
	double sarea;

	error_program = argv[0];
	if (argc <= 1)
		usage();

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++)
		{
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')		/* Look for any flags */
			{
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else
				{
				if ((fa+1) < argc)
					{
					if (argv[fa+1][0] != '-')
						{
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
						}
					}
				}

			if (argv[fa][1] == '?')
				usage();

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;
			else 
				usage();
			}
		else
			break;
		}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(inname,argv[fa]);
	strcat(inname,".ti3");
	strcpy(outname,argv[fa]);
	strcat(outname,".pr1");

	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI3"); 	/* our special input type is Calibration Target Information 3 */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI3 format file");
	if (icg->ntables != 1)
		error ("Input file doesn't contain exactly one table");

	if ((ed.npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	if (verb) {
		printf("No of test patches = %d\n",ed.npat);
	}

	if ((ed.cols = (col *)malloc(sizeof(col) * ed.npat)) == NULL)
		error("Malloc failed!");

	/* Setup output cgats file */
	/* This is a simple interpolation CMYK -> XYZ device profile */
	ocg = new_cgats();	/* Create a CGATS structure */
	ocg->add_other(ocg, "PROF1"); 		/* our special type is Profile type 1 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Device Profile Type 1",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll sprof", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);

	/* Figure out the color space */
	if ((ti = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file doesn't contain keyword COLOR_REPS");
	if (strcmp(icg->t[0].kdata[ti],"CMYK_XYZ") == 0) {
		int ci, mi, yi, ki;
		int Xi, Yi, Zi;
		if ((ci = icg->find_field(icg, 0, "CMYK_C")) < 0)
			error ("Input file doesn't contain field CMYK_C");
		if (icg->t[0].ftype[ci] != r_t)
			error ("Field CMYK_C is wrong type");
		if ((mi = icg->find_field(icg, 0, "CMYK_M")) < 0)
			error ("Input file doesn't contain field CMYK_M");
		if (icg->t[0].ftype[mi] != r_t)
			error ("Field CMYK_M is wrong type");
		if ((yi = icg->find_field(icg, 0, "CMYK_Y")) < 0)
			error ("Input file doesn't contain field CMYK_Y");
		if (icg->t[0].ftype[yi] != r_t)
			error ("Field CMYK_Y is wrong type");
		if ((ki = icg->find_field(icg, 0, "CMYK_K")) < 0)
			error ("Input file doesn't contain field CMYK_K");
		if (icg->t[0].ftype[ki] != r_t)
			error ("Field CMYK_K is wrong type");
		if ((Xi = icg->find_field(icg, 0, "XYZ_X")) < 0)
			error ("Input file doesn't contain field XYZ_X");
		if (icg->t[0].ftype[Xi] != r_t)
			error ("Field XYZ_X is wrong type");
		if ((Yi = icg->find_field(icg, 0, "XYZ_Y")) < 0)
			error ("Input file doesn't contain field XYZ_Y");
		if (icg->t[0].ftype[Yi] != r_t)
			error ("Field XYZ_Y is wrong type");
		if ((Zi = icg->find_field(icg, 0, "XYZ_Z")) < 0)
			error ("Input file doesn't contain field XYZ_Z");
		if (icg->t[0].ftype[Zi] != r_t)
			error ("Field XYZ_Z is wrong type");
		for (i = 0; i < ed.npat; i++) {
			double XYZ[3];
			ed.cols[i].c = *((double *)icg->t[0].fdata[i][ci]) / 100.0;
			ed.cols[i].m = *((double *)icg->t[0].fdata[i][mi]) / 100.0;
			ed.cols[i].y = *((double *)icg->t[0].fdata[i][yi]) / 100.0;
			ed.cols[i].k = *((double *)icg->t[0].fdata[i][ki]) / 100.0;
			XYZ[0] = *((double *)icg->t[0].fdata[i][Xi]) / 100.0;
			XYZ[1] = *((double *)icg->t[0].fdata[i][Yi]) / 100.0;
			XYZ[2] = *((double *)icg->t[0].fdata[i][Zi]) / 100.0;
			icmXYZ2Lab(&icmD50, ed.cols[i].Lab, XYZ);
		}

		/* Initialise the model */
		ed.gam[0] = 1.0;		/* First four are CMYK gamma values */
		ed.gam[1] = 1.0;
		ed.gam[2] = 1.0;
		ed.gam[3] = 1.0;

		/* Initialise interpolation end points for each combination of primary, */
		/* with all combinations close to black being represented by param[7]. */
		ed.k[0][0] = .82; ed.k[1][0] = .83; ed.k[2][0] = .75;		/* White */
		ed.k[0][1] = .66; ed.k[1][1] = .72; ed.k[2][1] = .05;		/*   Y  */
		ed.k[0][2] = .27; ed.k[1][2] = .12; ed.k[2][2] = .06;		/*  M   */
		ed.k[0][3] = .27; ed.k[1][3] = .12; ed.k[2][3] = .00;		/*  MY  */
		ed.k[0][4] = .09; ed.k[1][4] = .13; ed.k[2][4] = .44;		/* C    */
		ed.k[0][5] = .03; ed.k[1][5] = .10; ed.k[2][5] = .04;		/* C Y  */
		ed.k[0][6] = .02; ed.k[1][6] = .01; ed.k[2][6] = .05;		/* CM   */
		ed.k[0][7] = .01; ed.k[1][7] = .01; ed.k[2][7] = .01;		/* Black */

		sarea = 0.3;
		presid = dresid = 100.0;
		for (k=0; /* dresid > 0.0001 && */ k < 40; k++) {	/* Untill we're done */
			double sresid;
			double sr[8];
			double p[8];

			/* Adjust the gamma */
			for (i = 0; i < 4; i++)
				sr[i] = 0.1;			/* Device space search radius */
			if (powell(&resid[3], 4, &ed.gam[0], sr,  0.1, 1000, efunc1, (void *)&ed, NULL, NULL) != 0)
				error ("Powell failed");

			/* Adjust the primaries */
			calc_bc(&ed);		/* Calculate blend coefficients */
			for (i = 0; i < 8; i++)
				sr[i] = 0.2;			/* Device space search radius */
			sresid = 99.0;
			for (j = 0; j < 3; j++) {	/* For each of X, Y and Z */
				ed.xyzi = j;

				for (i = 0; i < 8; i++)
					p[i] = ed.k[j][i];
printf("##############\n");
printf("XYZ = %d\n",j);
				if (powell(&resid[j], 8, p, sr,  0.1, 1000, efunc2, (void *)&ed, NULL, NULL) != 0)
					error ("Powell failed");

				for (i = 0; i < 8; i++)
					ed.k[j][i] = p[i];

				if (sresid > resid[j])
					sresid = resid[j];
			}
			dresid = presid - sresid;
			if (dresid < 0.0)
				dresid = 100.0;
			presid = sresid;
printf("~1 presid = %f, sresid = %f, dresid = %f\n",presid, sresid, dresid);
		}

		/* Fields we want */
		ocg->add_kword(ocg, 0, "DSPACE","CMYK", NULL);
		ocg->add_kword(ocg, 0, "DTYPE","PRINTER", NULL);
		ocg->add_field(ocg, 0, "PARAM_ID", i_t);
		ocg->add_field(ocg, 0, "PARAM", r_t);

		/* Output model parameters */
		for (j = 0; j < 4; j++)
			ocg->add_set(ocg, 0, j, ed.gam[j]);

		for (j = 0; j < 3; j++) {
			for (i = 0; i < 8; i++)
				ocg->add_set(ocg, 0, 10 * (j + 1) + i, 100.0 * ed.k[j][i]);
		}

		if (verb) {
			double aver = 0.0;
			double maxer = 0.0;
			for (i = 0; i < ed.npat; i++) {
				double err = sqrt(ed.cols[i].err);
				if (err > maxer)
					maxer = err;
				aver += err;
			}
			aver = aver/((double)i);
			printf("Average fit error = %f, maximum = %f\n",aver,maxer);
		}
	} else if (strcmp(icg->t[0].kdata[ti],"RGB") == 0) {
		error ("We can't handle RGB !");
	} else if (strcmp(icg->t[0].kdata[ti],"W") == 0) {
		error ("We can't handle Grey !");
	} else
		error ("Input file keyword COLOR_REPS has unknown value");

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	free(ed.cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */

	return 0;
}

/******************************************************************/
/* Error/debug output routines */
/******************************************************************/

void
usage(void) {
	fprintf(stderr,"Create Simple CMYK Device Profile, Version %s\n",VERSION);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: %s [-v] outfile\n",error_program);
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," outfile         Base name for input.tr3/output.pr1 file\n");
	exit(1);
	}


