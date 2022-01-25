
/* 
 * Argyll Color Management System
 *
 * Read in the device data from Colorblind device files,
 * and convert it into a .ti3 CGATs format suitable for
 * the Argyll CMM.
 *
 * Derived from  kodak2cgats.c 
 * Author: Graeme W. Gill
 * Date:   16/11/00
 *
 * Copyright 2000, 2010, Graeme W. Gill
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#define VERSION "1.0"

/* TTBD
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

void
usage(void) {
	fprintf(stderr,"Convert Colorblind raw device profile data to Argyll data, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: cb2ti3 [-v] [-l limit] infile outfile\n");
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -l limit        Set inklimit in .ti3 file\n");
	fprintf(stderr," infile	         Base name for input.CMY and input.nCIE file\n");
	fprintf(stderr," outfile         Base name for output.ti3 file\n");
	exit(1);
	}

int main(int argc, char *argv[])
{
	int i;
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	static char tarname[200] = { 0 };		/* Input .CMY file */
	static char inname[200] = { 0 };		/* Input .nCIE file */
	static char outname[200] = { 0 };		/* Output cgats .ti3 file base name */
	cgats *cmy;			/* Input RGB reference file */
	int f_id1, f_c, f_m, f_y;	/* Field indexes */
	cgats *ncie;		/* Inpit CIE readings file */
	int f_id2, f_xx, f_yy, f_zz;	/* Field indexes */
	cgats *ocg;			/* output cgats structure */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	int npat = 0;		/* Number of patches */

	error_program = "cb2ti3";

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
	strcpy(tarname,argv[fa++]);
	strcat(inname,".CMY");
	strcat(tarname,".nCIE");

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(outname, argv[fa++]);
	strcat(outname,".ti3");

	/* Open up the Input CMY reference file */
	cmy = new_cgats();	/* Create a CGATS structure */
	cmy->add_other(cmy, "CBTA"); 	/* Colorblind Target file */
	if (cmy->read_name(cmy, inname))
		error ("Read: Can't open file '%s'",inname);
	if (cmy->ntables == 0 || cmy->t[0].tt != tt_other || cmy->t[0].oi != 0)
		error ("Input file isn't a 'CBTA' format file");
	if (cmy->ntables != 1)
		fprintf(stderr,"Input file '%s' doesn't contain exactly one table",inname);

	if ((npat = cmy->t[0].nsets) <= 0)
		error("No patches");

	if ((f_id1 = cmy->find_field(cmy, 0, "SAMPLE_ID")) < 0)
		error("Input file doesn't contain field SAMPLE_ID");
	if (cmy->t[0].ftype[f_id1] != nqcs_t)
		error("Field SAMPLE_ID is wrong type");

	if ((f_c = cmy->find_field(cmy, 0, "C")) < 0)
		error("Input file doesn't contain field C");
	if (cmy->t[0].ftype[f_c] != r_t)
		error("Field C is wrong type");

	if ((f_m = cmy->find_field(cmy, 0, "M")) < 0)
		error("Input file doesn't contain field M");
	if (cmy->t[0].ftype[f_m] != r_t)
		error("Field M is wrong type");

	if ((f_y = cmy->find_field(cmy, 0, "Y")) < 0)
		error("Input file doesn't contain field Y");
	if (cmy->t[0].ftype[f_y] != r_t)
		error("Field Y is wrong type");

	/* Open up the input nCIE device data file */
	ncie = new_cgats();	/* Create a CGATS structure */
	ncie->add_other(ncie, "CBPR"); 	/* Colorblind Printer Response file */
	if (ncie->read_name(ncie, tarname))
		error ("Read: Can't open file '%s'",tarname);
	if (ncie->ntables == 0 || ncie->t[0].tt != tt_other || ncie->t[0].oi != 0)
		error ("Input file isn't a 'CBTA' format file");
	if (ncie->ntables != 1)
		fprintf(stderr,"Input file '%s' doesn't contain exactly one table",tarname);

	if (npat != ncie->t[0].nsets)
		error("Number of patches doesn't match");

	if ((f_id2 = ncie->find_field(ncie, 0, "SAMPLE_ID")) < 0)
		error("Input file doesn't contain field SAMPLE_ID");
	if (ncie->t[0].ftype[f_id2] != nqcs_t)
		error("Field SAMPLE_ID is wrong type");

	if ((f_xx = ncie->find_field(ncie, 0, "XYZ_X")) < 0)
		error("Input file doesn't contain field XYZ_X");
	if (ncie->t[0].ftype[f_xx] != r_t)
		error("Field XYZ_X is wrong type");

	if ((f_yy = ncie->find_field(ncie, 0, "XYZ_Y")) < 0)
		error("Input file doesn't contain field XYZ_Y");
	if (ncie->t[0].ftype[f_yy] != r_t)
		error("Field XYZ_Y is wrong type");

	if ((f_zz = ncie->find_field(ncie, 0, "XYZ_Z")) < 0)
		error("Input file doesn't contain field XYZ_Z");
	if (ncie->t[0].ftype[f_zz] != r_t)
		error("Field XYZ_Z is wrong type");

	/* Setup output cgats file */
	ocg = new_cgats();	/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI3"); 	/* our special type is Calibration Target Information 3 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll target", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */

	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);

	ocg->add_field(ocg, 0, "RGB_R", r_t);
	ocg->add_field(ocg, 0, "RGB_G", r_t);
	ocg->add_field(ocg, 0, "RGB_B", r_t);
	ocg->add_kword(ocg, 0, "COLOR_REP","RGB_XYZ", NULL);
	ocg->add_field(ocg, 0, "XYZ_X", r_t);
	ocg->add_field(ocg, 0, "XYZ_Y", r_t);
	ocg->add_field(ocg, 0, "XYZ_Z", r_t);

	/* Write out the patch info to the output CGATS file */
	for (i = 0; i < npat; i++) {
		char id[100];
		double rgb[3]; 
		double xyz[3]; 

		if (strcmp(((char *)cmy->t[0].fdata[i][f_id1]), 
		           ((char *)ncie->t[0].fdata[i][f_id2])) != 0) {
			error("Patch label mismatch, patch %d, '%s' != '%s'\n",
			       i, ((char *)cmy->t[0].fdata[i][f_id1]), 
		              ((char *)ncie->t[0].fdata[i][f_id2]));
		}

		rgb[0] = 100.0 - *((double *)cmy->t[0].fdata[i][f_c]);	/* Convert to RGB */
		rgb[1] = 100.0 - *((double *)cmy->t[0].fdata[i][f_m]);
		rgb[2] = 100.0 - *((double *)cmy->t[0].fdata[i][f_y]);

		xyz[0] = *((double *)ncie->t[0].fdata[i][f_xx]);	
		xyz[1] = *((double *)ncie->t[0].fdata[i][f_yy]);
		xyz[2] = *((double *)ncie->t[0].fdata[i][f_zz]);

		sprintf(id, "%d", i+1);
		ocg->add_set(ocg, 0, id, rgb[0], rgb[1], rgb[2], 
		                         xyz[0], xyz[1], xyz[2]);
	}

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	ncie->del(ncie);		/* Clean up */
	cmy->del(cmy);
	ocg->del(ocg);

	return 0;
}



