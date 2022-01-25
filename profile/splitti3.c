/* 
 * Argyll Color Management System
 * Split a .ti3 (or other CGATS like) file into two parts.
 *
 * Author: Graeme W. Gill
 * Date:   14/12/2005
 *
 * Copyright 2005, 2010 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in a CGATS .ti3 file, and splits it into
 * two .ti3 files, spreading the readings between them. 
 * This is intended for use in verifying the profiler.
 */

/*
 * TTBD:
	
	This doesn't pass calibration table information through.
	(ie. should copy all tables after the first.)

	Write a companion "combineti3" to merge .ti3's together.
 */

#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include <string.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include <sys/types.h>
#include <time.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "cgats.h"
#include "xicc.h"
#include "insttypes.h"
#include "sort.h"

void
usage(void) {
	fprintf(stderr,"Split a .ti3 into two, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: splitti3 [-options] input.ti3 output1.ti3 output2.ti3\n");
	fprintf(stderr," -v              Verbose - print each patch value\n");
	fprintf(stderr," -n no           Put no sets in first file, and balance in second file.\n");
	fprintf(stderr," -p percent      Put percent%% sets in first file, and balance in second file. (def. 50%%)\n");
	fprintf(stderr," -w              Put white patches in both files.\n");
	fprintf(stderr," -r seed         Use given random seed.\n");
	fprintf(stderr," input.ti3       File to be split up.\n");
	fprintf(stderr," output1.ti3     First output file\n");
	fprintf(stderr," output2.ti3     Second output file\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int numb = -1;			/* Number to put in first */
	double prop = 0.5;		/* Proportion to put in first */
	int dow = 0;			/* Put white patches in both files */
	int seed = 0x12345678;
	int doseed = 0;

	cgats *cgf = NULL;			/* cgats file data */
	char in_name[MAXNAMEL+1];	/* Patch filename  */

	cgats *cg1 = NULL;			/* cgats file data */
	char out_name1[MAXNAMEL+4+1]; /* VRML name */
	cgats *cg2 = NULL;			/* cgats file data */
	char out_name2[MAXNAMEL+4+1]; /* VRML name */

	cgats_set_elem *setel;		/* Array of set value elements */
	int *flags;					/* Point to destination of set */

	int i, j, n;

	error_program = "splitti3";

	if (argc <= 1)
		usage();

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
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

			if (argv[fa][1] == '?') {
				usage();

			} else if (argv[fa][1] == 'v') {
				verb = 1;

			} else if (argv[fa][1] == 'n') {
				if (na == NULL) usage();
				fa = nfa;
				numb = atoi(na);
				if (numb < 0) usage();

			} else if (argv[fa][1] == 'p') {
				if (na == NULL) usage();
				fa = nfa;
				prop = atoi(na);
				if (prop < 0) usage();
				prop = prop / 100.0;

			} else if (argv[fa][1] == 'w') {
				dow = 1;

			} else if (argv[fa][1] == 'r') {
				if (na == NULL) usage();
				fa = nfa;
				seed = atoi(na);
				doseed = 1;
			}

			else 
				usage();
		} else
			break;
	}

	if (doseed)
		rand32(seed);			/* Init seed deterministicaly */
	else
		rand32(time(NULL));		/* Init seed randomly */

	/* Get the file name arguments */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(in_name,argv[fa++],MAXNAMEL); in_name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(out_name1,argv[fa++],MAXNAMEL); out_name1[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(out_name2,argv[fa++],MAXNAMEL); out_name2[MAXNAMEL] = '\000';

	if ((cgf = new_cgats()) == NULL)
		error("Failed to create cgats object");
	cgf->add_other(cgf, ""); 	/* Allow any signature file */
	
	if (cgf->read_name(cgf, in_name))
		error("CGATS file '%s' read error : %s",in_name,cgf->err);
	
	if (cgf->ntables < 1)
		error ("Input file '%s' doesn't contain at least one table",in_name);
	
	/* Create the two output files */
	if ((cg1 = new_cgats()) == NULL)
		error("Failed to create cgats object");
	if ((cg2 = new_cgats()) == NULL)
		error("Failed to create cgats object");

	/* Duplicate the type of the file */
	if (cgf->t[0].tt == cgats_X) {
		cg1->add_other(cg1, cgf->cgats_type);
		cg1->add_table(cg1, tt_other, 0);
		cg2->add_other(cg2, cgf->cgats_type);
		cg2->add_table(cg2, tt_other, 0);
	} else if (cgf->t[0].tt == tt_other) {
		cg1->add_other(cg1, cgf->others[cgf->t[0].oi]);
		cg1->add_table(cg1, tt_other, 0);
		cg2->add_other(cg2, cgf->others[cgf->t[0].oi]);
		cg2->add_table(cg2, tt_other, 0);
	} else {
		cg1->add_table(cg1, cgf->t[0].tt, 0);
		cg2->add_table(cg1, cgf->t[0].tt, 0);
	}

	/* Duplicate all the keywords */
	for (i = 0; i < cgf->t[0].nkwords; i++) {
		cg1->add_kword(cg1, 0, cgf->t[0].ksym[i], cgf->t[0].kdata[i], NULL);
		cg2->add_kword(cg2, 0, cgf->t[0].ksym[i], cgf->t[0].kdata[i], NULL);
	}

	/* Duplicate all of the fields */
	for (i = 0; i < cgf->t[0].nfields; i++) {
		cg1->add_field(cg1, 0, cgf->t[0].fsym[i], cgf->t[0].ftype[i]);
		cg2->add_field(cg2, 0, cgf->t[0].fsym[i], cgf->t[0].ftype[i]);
	}

	if ((setel = (cgats_set_elem *)malloc(
	     sizeof(cgats_set_elem) * cgf->t[0].nfields)) == NULL)
		error("Malloc failed!");

	if ((flags = (int *)calloc(cgf->t[0].nsets, sizeof(int))) == NULL)
		error("Malloc failed!");
	
	if (numb < 0) {	/* Use percentage */
		numb = (int)(cgf->t[0].nsets * prop + 0.5);
	}
	if (numb > cgf->t[0].nsets)
		numb = cgf->t[0].nsets;

	if (verb)
		printf("Putting %d sets in '%s' and %d in '%s'\n",numb,out_name1,cgf->t[0].nsets-numb,out_name2);

	n = 0;

	/* If dow, add white patches to both sets */
	if (dow) {
		int ti;
		char *buf;
		char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };
		char *labfname[3] = { "LAB_L", "LAB_A", "LAB_B" };
		char *outc;
		int nmask;
		int nchan;
		char *bident;
		int chix[ICX_MXINKS];	/* Device chanel indexes */
		int pcsix[3];			/* PCS chanel indexes */
		int isin = 0;
		int isadd = 0;
		int isLab = 0;

		if ((ti = cgf->find_kword(cgf, 0, "DEVICE_CLASS")) < 0)
			error ("Input file doesn't contain keyword DEVICE_CLASS");

		if (strcmp(cgf->t[0].kdata[ti],"INPUT") == 0)
			isin = 1;

		if ((ti = cgf->find_kword(cgf, 0, "COLOR_REP")) < 0)
			error("Input file doesn't contain keyword COLOR_REPS");
		
		if ((buf = strdup(cgf->t[0].kdata[ti])) == NULL)
			error("Malloc failed");
		
		/* Split COLOR_REP into device and PCS space */
		if ((outc = strchr(buf, '_')) == NULL)
			error("COLOR_REP '%s' invalid", cgf->t[0].kdata[ti]);
		*outc++ = '\000';
		
		if (strcmp(outc, "XYZ") == 0) {
			isLab = 0;
		} else if (strcmp(outc, "LAB") == 0) {
			isLab = 1;
		} else
			error("COLOR_REP '%s' invalid (Neither XYZ nor LAB)", cgf->t[0].kdata[ti]);
	
		if ((nmask = icx_char2inkmask(buf)) == 0) {
			error ("File '%s' keyword COLOR_REPS has unknown device value '%s'",in_name,buf);
		}

		if (nmask & ICX_ADDITIVE)
			isadd = 1;

		nchan = icx_noofinks(nmask);
		bident = icx_inkmask2char(nmask, 0); 

		/* Find device fields */
		for (j = 0; j < nchan; j++) {
			int ii, imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : bident,
			                      icx_ink2char(imask));

			if ((ii = cgf->find_field(cgf, 0, fname)) < 0)
				error ("Input file doesn't contain field %s",fname);
			if (cgf->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type",fname);
			chix[j] = ii;
		}

		/* Find PCS fields */
		for (j = 0; j < 3; j++) {
			int ii;

			if ((ii = cgf->find_field(cgf, 0, isLab ? labfname[j] : xyzfname[j])) < 0)
				error ("Input file doesn't contain field %s",isLab ? labfname[j] : xyzfname[j]);
			if (cgf->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type",isLab ? labfname[j] : xyzfname[j]);
			pcsix[j] = ii;
		}

		if (isin) {
			int wix = -1;
			double wv = -1e60;
			int pcsy = 1;

			if (isLab)
				pcsy = 0;

			/* We assume that the white point is the patch with the */
			/* highest L* or Y value. */
			for (i = 0; i < cgf->t[0].nsets; i++) {
				double val;

				val = *((double *)cgf->t[0].fdata[i][pcsix[pcsy]]);
				if (val > wv) {
					wv = val;
					wix = i;
				}
			}
			if (wix > 0) {
				n++;
				flags[wix] = 3;
				if (verb)
					printf("Found input white patch index %d\n",wix);
			}

		} else {

			if (isadd) {
				for (i = 0; i < cgf->t[0].nsets; i++) {
					for (j = 0; j < nchan; j++) {
						if (*((double *)cgf->t[0].fdata[i][chix[j]]) < 99.99)
							break;
					}
					if (j >= nchan) {
						n++;
						flags[i] = 3;
						if (verb)
							printf("Found additive white patch index %d\n",i);
					}
				}
			} else {
				for (i = 0; i < cgf->t[0].nsets; i++) {
					for (j = 0; j < nchan; j++) {
						if (*((double *)cgf->t[0].fdata[i][chix[j]]) > 0.01)
							break;
					}
					if (j >= nchan) {
						n++;
						flags[i] = 3;
						if (verb)
							printf("Found subtractive white patch index %d\n",i);
					}
				}
			}
		}
		free(bident);
	}

	/* Chose which of the sets go into file 1 and 2*/
	for (;n < numb;) {
		i = i_rand(0, cgf->t[0].nsets-1);
		if (flags[i] == 0) {
			flags[i] = 1;
			n++;
		}
	}

	/* Assume any patch not flagged goes to 2 */
	for (i = 0; i < cgf->t[0].nsets; i++) {
		if (flags[i] == 0)
			flags[i] = 2;
	}

	/* Copy them approproately */
	for (i = 0; i < cgf->t[0].nsets; i++) {
		cgf->get_setarr(cgf, 0, i, setel);
		if (flags[i] & 1) {
			cg1->add_setarr(cg1, 0, setel);
		}
		if (flags[i] & 2) {
			cg2->add_setarr(cg2, 0, setel);
		}
	}

	/* Write out the files */
	if (cg1->write_name(cg1, out_name1))
		error("CGATS file '%s' write error : %s",out_name1,cg1->err);
	if (cg2->write_name(cg2, out_name2))
		error("CGATS file '%s' write error : %s",out_name2,cg2->err);
	

	free(flags);
	free(setel);
	
	return 0;
}





