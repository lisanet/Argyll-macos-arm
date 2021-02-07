
/* 
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    2006/5/9
 * Version: 1.00
 *
 * Copyright 2006 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 */

/*
 * This is some test code to test the Daylight and Plankian spectra, 
 * Correlated and Visual Color Temperatures, and CRI.
 * and plot a spectrum, CMF or CCSS.
 */

#include <stdio.h>
#include <math.h>
#include "aconfig.h"
#include "cgats.h"
#include "xspect.h"
#include "numlib.h"
#include "plot.h"
#include "ui.h"

#define PLANKIAN
#define XRES 500

#define MAXGRAPHS 10

/* Display a spectrum etc. */
/* We are guaranteed that the x range/increments are identical, */
/* and that there is only one spectrum if douv */
static int do_spec(
	char name[MAXGRAPHS][200],
	xspect *sp,
	int nsp,				/* Number of sp */
	inst_meas_type mt,		/* Measurement type */
	int dozero,				/* Include zero in the range */
	int noplot,				/* Don't actually plot */
	icxDensityType dens,	/* Density type to print */
	int douv,				/* Do variation of added UV test */
	double uvmin,
	double uvmax
) {
	int n, i, j, k, m;
	double wl_short, wl_long;		/* Common range */
	double xyz[3];		/* Color temperature */
	double Yxy[3];
	double Lab[3];		/* D50 Lab value */
	double xx[XRES];
	double yy[10][XRES];
	double *yp[10];
	double cct, vct;
	double cct_xyz[3], vct_xyz[3];
	double cct_lab[3], vct_lab[3];
	double cdt, vdt;
	double cdt_xyz[3], vdt_xyz[3];
	double cdt_lab[3], vdt_lab[3];
	icmXYZNumber wp;
	double de;
	double uv = uvmin;
	double step = 0.1;
	xspect tsp;			/* Spectrum with possible UV added */
	char *color[] = {
		"Black", "Red", "Green", "Blue", "Yellow", "Purple", "Brown", "Orange", "Grey", "Magenta"
	};

	printf("\n");

	for (j = 0; j < 10; j++)
		yp[j] = NULL;

	if (nsp > 10)
		nsp = 10;

	m = 0;				/* offset in output array */
	n = 1;

	wl_short = 1e6;
	wl_long = -1e6;
	for (k = 0; k < nsp; k++) {
		if (sp[k].spec_wl_long > wl_long)
			wl_long = sp[k].spec_wl_long;
		if (sp[k].spec_wl_short < wl_short)
			wl_short = sp[k].spec_wl_short;
	}

	if (douv) {
		n = 1 + (int)(0.5 + (uvmax-uvmin)/0.1);
		if (n > 9)
			n = 9;		/* Don't use white */
		if (n > 1)
			step = (uvmax-uvmin)/(n-1.0); 
	}

	for (k = 0; k < nsp; k++) {
		tsp = sp[k];
		for (uv = uvmax, j = 0; j < n; j++, uv -= step) {

			if (douv) {
				printf("UV level = %f\n",uv);
				xsp_setUV(&tsp, &sp[k], uv);
			}
				
			printf("Type = %s [%s]\n",name[k], color[k]);

			if (mt == inst_mrt_none
			 || mt == inst_mrt_emission
			 || mt == inst_mrt_ambient
			 || mt == inst_mrt_emission_flash
			 || mt == inst_mrt_ambient_flash) {
	
				/* Compute absolute XYZ of illuminant */
				if (icx_ill_sp2XYZ(xyz, icxOT_CIE_1931_2, NULL, icxIT_custom, 0, &tsp, 1) != 0) 
					warning("icx_ill_sp2XYZ returned error");
	
				printf("Abs. Y = %f\n", xyz[1]);

				/* Normalise to Y = 1 */
				xyz[0] /= xyz[1];
				xyz[2] /= xyz[1];
				xyz[1] /= xyz[1];

				icmXYZ2Yxy(Yxy, xyz);
				icmXYZ2Lab(&icmD50, Lab, xyz);

				printf("XYZ = %f %f %f, x,y = %f %f\n", xyz[0], xyz[1], xyz[2], Yxy[1], Yxy[2]);
				printf("D50 L*a*b* = %f %f %f\n", Lab[0], Lab[1], Lab[2]);
				
				/* Compute CCT */
				if ((cct = icx_XYZ2ill_ct(cct_xyz, icxIT_Ptemp, icxOT_CIE_1931_2, NULL, xyz, NULL, 0)) < 0)
					warning("Got bad cct\n");
	
				/* Compute VCT */
				if ((vct = icx_XYZ2ill_ct(vct_xyz, icxIT_Ptemp, icxOT_CIE_1931_2, NULL, xyz, NULL, 1)) < 0)
					warning("Got bad vct\n");
	
				printf("CCT = %f, VCT = %f\n",cct, vct);
	
				/* Compute CDT */
				if ((cct = icx_XYZ2ill_ct(cct_xyz, icxIT_Dtemp, icxOT_CIE_1931_2, NULL, xyz, NULL, 0)) < 0)
					warning("Got bad cct\n");
	
				/* Compute VDT */
				if ((vct = icx_XYZ2ill_ct(vct_xyz, icxIT_Dtemp, icxOT_CIE_1931_2, NULL, xyz, NULL, 1)) < 0)
					warning("Got bad vct\n");
	
				printf("CDT = %f, VDT = %f\n",cct, vct);
	
				{
					int invalid = 0;
					double RR[14];
					double cri;
					cri = icx_CIE1995_CRI(&invalid, RR, &tsp);
					printf("CRI = %.1f [ R9 = %.1f ]%s\n",cri,RR[9-1],invalid ? " (Invalid)" : "");
				}
				{
					int invalid = 0;
					double tlci;
					tlci = icx_EBU2012_TLCI(&invalid, &tsp);
					printf("TLCI = %.1f%s\n",tlci,invalid ? " (Invalid)" : "");
				}
	
				/* Use modern color difference - gives a better visual match */
				icmAry2XYZ(wp, vct_xyz);
				icmXYZ2Lab(&wp, cct_lab, cct_xyz);
				icmXYZ2Lab(&wp, vct_lab, vct_xyz);
				de = icmCIE2K(cct_lab, vct_lab);
				printf("CIEDE2000 Delta E = %f\n",de);

			} else if (mt == inst_mrt_none
			 || mt == inst_mrt_reflective
			 || mt == inst_mrt_transmissive) {
				double cmyv[4];
	
				printf("CIE values under D50 illuminant:\n");

				if (icx_sp2XYZ(xyz, icxOT_CIE_1931_2, NULL, icxIT_D50, 0, NULL, &tsp) != 0) 
					warning("icx_sp2XYZ returned error");
	
				icmXYZ2Yxy(Yxy, xyz);
				icmXYZ2Lab(&icmD50, Lab, xyz);

				printf("XYZ = %f %f %f, x,y = %f %f\n", xyz[0], xyz[1], xyz[2], Yxy[1], Yxy[2]);
				printf("D50 L*a*b* = %f %f %f\n", Lab[0], Lab[1], Lab[2]);
				
				if (dens != icxDT_none) {
					/* Density */
					xsp_density(cmyv, &tsp, dens);
					printf("CMYV %s = %f %f %f %f\n", xsp_density_desc(dens), cmyv[0], cmyv[1], cmyv[2], cmyv[3]);
				}
	
			} else {
				printf("Unhandled measurement type '%s'\n",meas_type2str(mt));
			}

			/* Plot spectrum out */
			for (i = 0; i < XRES; i++) {
				double ww;

				ww = (wl_long - wl_short)
				   * ((double)i/(XRES-1.0)) + wl_short;
			
				xx[i] = ww;
				yy[(m + k + j) % 10][i] = value_xspect(&tsp, ww);
			}
			yp[(m + k + j) % 10] = &yy[(m + k + j) % 10][0];
			printf("\n");
		}
	}
	if (!noplot) {
		do_plot10(xx, yp[0], yp[1], yp[2], yp[3], yp[4], yp[5], yp[6], yp[7], yp[8], yp[9], XRES, dozero);
	}


	return 0;
}


void usage(void) {
	fprintf(stderr,"Plot spectrum and calculate CCT and VCT\n");
	fprintf(stderr,"Author: Graeme W. Gill\n");
	fprintf(stderr,"usage: specplot [infile.sp]\n");
	fprintf(stderr," -v               verbose\n");
	fprintf(stderr," -c               combine multiple files into one plot\n");
	fprintf(stderr," -z               don't make range cover zero\n");
	fprintf(stderr," -s               don't plot spectra\n");
	fprintf(stderr," -d I|A|M|T|E     print density values\n");
	fprintf(stderr," -u level         plot effect of adding estimated UV level\n");
	fprintf(stderr," -U               plot effect of adding range of estimated UV level\n");
	fprintf(stderr," [infile.sp ...]  spectrum files to plot\n");
	fprintf(stderr,"                  default is all built in illuminants\n");
	exit(1);
}

int
main(
	int argc,
	char *argv[]
) {
	int fa, nfa;			/* argument we're looking at */
	int k;
	int verb = 0;
	int comb = 0;
	int zero = 1;
	int noplot = 0;
	icxDensityType dens = icxDT_none;
	double temp;
	xspect sp[MAXGRAPHS];
	icxIllumeType ilType;
	int douv = 0;
	double uvmin = -1.0, uvmax = 1.0;
	char buf[MAXGRAPHS][200];

	error_program = argv[0];

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

			/* Show added UV  */
			if (argv[fa][1] == 'u') {
				douv = 1;

				fa = nfa;
				if (na == NULL)
					usage();

				uvmin = uvmax = atof(na);
				if (uvmin < -10.0 || uvmax > 10.0)
					usage();
			}

			else if (argv[fa][1] == 'U') {
				douv = 1;
			}

			/* Verbosity */
			else if (argv[fa][1] == 'v') {
				verb = 1;

			} else if (argv[fa][1] == 'c') {
				comb = 1;

			} else if (argv[fa][1] == 'z') {
				zero = 0;

			} else if (argv[fa][1] == 's') {
				noplot = 1;

			} else if (argv[fa][1] == 'd') {
				if (na == NULL)
					usage();

				switch(na[0]) {
					case 'N':
						dens = icxDT_none;
						break;
					case 'I':
						dens = icxDT_ISO;
						break;
					case 'A':
						dens = icxDT_A;
						break;
					case 'M':
						dens = icxDT_M;
						break;
					case 'T':
						dens = icxDT_T;
						break;
					case 'E':
						dens = icxDT_E;
						break;
					default:
						usage();
				}
				fa = nfa;

			} else {
				usage();
			}
		}
		else
			break;
	}

	if (fa < argc && argv[fa][0] != '-') {	/* Got file arguments */
		int nsp = 0;		/* Current number in sp[] */
		int soff = 0;		/* Offset within file */
		int maxgraphs = MAXGRAPHS;
		int eof;

		if (douv)
			maxgraphs = 1;
		
		nsp = 0;
	
		/* Until we run out */
		for (;;) {
			int i, nret, nreq;
			inst_meas_type mt;

			/* If we've got to the limit of each plot, */
			/* or at least one and we're not combining files and at start of a new file, */
			/* or at least one and there are no more files */
			if (nsp >= MAXGRAPHS || (nsp > 0 && ((!comb && soff == 0) || fa >= argc))) {
				/* Plot what we've got */
				do_spec(buf, sp, nsp, mt, zero, noplot, dens, douv, uvmin, uvmax);
				nsp = 0;
			}

			if (fa >= argc)		/* No more files */
				break;

			/* Read as many spectra from the file as possible */
			nreq = MAXGRAPHS - nsp;

			if (read_nxspect(&sp[nsp], &mt, argv[fa], &nret, soff, nreq, 0) != 0) {
				error ("Unable to read custom spectrum, CMF or CCSS '%s'",argv[fa]);
			}

			for (i = 0; i < nret; i++) {
				xspect_denorm(&sp[nsp + i]);
				sprintf(buf[nsp + i],"File '%s' spect %d",argv[fa], soff + i);
			}
			nsp += nret;
			soff += nret;
			if (nret < nreq) {		/* We're done with this file */
				fa++;
				soff = 0;
			}
		}

	} else {

		/* For each standard illuminant */
		for (ilType = icxIT_A; ilType <= icxIT_F10; ilType++) {
			char *inm = NULL;
	
			switch (ilType) {
			    case icxIT_A:
					inm = "A"; break;
			    case icxIT_C:
					inm = "C"; break;
			    case icxIT_D50:
					inm = "D50"; break;
			    case icxIT_D50M2:
					inm = "D50M2"; break;
			    case icxIT_D55:
					inm = "D55"; break;
			    case icxIT_D65:
					inm = "D65"; break;
			    case icxIT_D75:
					inm = "D75"; break;
			    case icxIT_E:
					inm = "E"; break;
			    case icxIT_F5:
					inm = "F5"; break;
			    case icxIT_F8:
					inm = "F8"; break;
			    case icxIT_F10:
					inm = "F10"; break;
				default:
					inm = "Unknown"; break;
					break;
			}
	
			if (standardIlluminant(&sp[0], ilType, 0) != 0)
				error ("standardIlluminant returned error for %d (%s)",ilType,inm);
		
			strcpy(buf[0],inm);
			do_spec(buf, sp, 1, inst_mrt_ambient, zero, noplot, dens, douv, uvmin, uvmax);
		}

		/* For each material and illuminant */
		for (k = 0; k < 2; k++) {
			double start, end;

			if (k == 0)
				start = 3400, end = 9000;
			else
				start = 500, end = 9000;

			for (temp = start; temp <= end; temp += 500) {
		
				ilType = k == 0 ? icxIT_Dtemp : icxIT_Ptemp;
	
				if (standardIlluminant(&sp[0], ilType, temp) != 0)
					error ("standardIlluminant returned error for %s temp %f",
					k == 0 ? "Daylight" : "Plankian", temp);
		
				sprintf(buf[0], "%s at %f", k == 0 ? "Daylight" : "Black body", temp);
	
				do_spec(buf, sp, 1, inst_mrt_ambient, zero, noplot, dens, douv, uvmin, uvmax);
			}
		}

	}
	return 0;
}








