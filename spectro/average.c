/* 
 * Argyll Color Management System
 * Average one or more .ti3 (or other CGATS like) file values together.
 * If just one file is supplied, all patches within it with
 * the same device value are averaged together. 
 *
 * Author: Graeme W. Gill
 * Date:   18/1/2011
 *
 * Copyright 2005, 2010, 2011 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * (based on splitti3.c)
 */

/*
 * TTBD:
	
	Should probably re-index SAMPLE_ID field

 */


#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "sort.h"
#include "cgats.h"
#include "xicc.h"
#include "insttypes.h"

static double average(double *vals, int nvals);
static double median(double *vals, int nvals);
static void geommed(double res[3], double vals[][3], int nvals);

void usage(char *diag, ...) {
	int i;
	fprintf(stderr,"Average or merge values in .ti3 like files, Version %s\n",ARGYLL_VERSION_STR);
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: average [-options] input1.ti3 input2.ti3 ... output.ti3\n");
	fprintf(stderr," -v              Verbose\n");
	fprintf(stderr," -e              Median rather than average\n");
	fprintf(stderr," -g              Geometric Median of PCS in encoded space\n");
	fprintf(stderr," -L              Geometric Median of PCS in L*a*b* space\n");
	fprintf(stderr," -X              Geometric Median of PCS in XYZ space\n");
	fprintf(stderr," -m              Merge rather than average\n");
	fprintf(stderr," input1.ti3      First input file\n");
	fprintf(stderr," input2.ti3      Second input file\n");
	fprintf(stderr," ...             etc.\n");
	fprintf(stderr," output.ti3      Resulting averaged or merged output file\n");
	exit(1);
}

/* Information about each file */
struct _inpinfo {
	char name[MAXNAMEL+1];
	cgats *c;
}; typedef struct _inpinfo inpinfo;

int main(int argc, char *argv[]) {
	int fa,nfa;					/* current argument we're looking at */
	int verb = 0;
	int domedian = 0;			/* Median rather than average */
	int dogeom = 0;				/* Do geometric median of PCS, 2 = Lab, 3 = PCS */
	int domerge = 0;			/* Merge rather than average */

	int ninps = 0;				/* Number of input files */
	inpinfo *inps;				/* Input file info. inp[ninp] == output file */
	cgats *ocg;					/* Copy of output cgats * */

	cgats_set_elem *setel;		/* Array of set value elements */
	int *flags;					/* Point to destination of set */

	int nchan = 0;				/* Number of device channels */
	int chix[ICX_MXINKS];		/* Device chanel indexes */
	int npcs = 0;

	int haspcs[2] =  { 0 };		/* Has Lab, XYZ */
	int pcsix[3][3];			/* Lab, XYZ chanel indexes */

	int i, j, n;

	error_program = "average";

	if (argc < 2)
		usage("Too few arguments (%d, minimum is 2)",argc-1);

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
				usage("Usage requested");

			/* Median */
			else if (argv[fa][1] == 'e') {
				domedian = 1;
			}

			/* Geometric Median of PCS */
			else if (argv[fa][1] == 'g') {
				dogeom = 1;
			}

			/* Geometric Median of PCS in L*a*b* */
			else if (argv[fa][1] == 'L') {
				dogeom = 2;
			}

			/* Geometric Median of PCS in XYZ */
			else if (argv[fa][1] == 'X') {
				dogeom = 3;
			}

			/* Merge */
			else if (argv[fa][1] == 'm') {
				domerge = 1;
			}

			/* Verbosity */
			else if (argv[fa][1] == 'v') {
				verb = 1;
			}

			else  {
				usage("Unknown flag '%c'",argv[fa][1]);
			}

		} else if (argv[fa][0] != '\000') {
			/* Get the next filename */
		
			if (ninps == 0)
				inps = (inpinfo *)malloc(sizeof(inpinfo));
			else
				inps = (inpinfo *)realloc(inps, (ninps+1) * sizeof(inpinfo));
			if (inps == NULL)
				error("Malloc failed in allocating space for file info.");

			memset((void *)&inps[ninps], 0, sizeof(inpinfo));
			strncpy(inps[ninps].name,argv[fa],MAXNAMEL);
			inps[ninps].name[MAXNAMEL] = '\000';

			ninps++;
		} else {
			break;
		}
	}

	if (ninps < 2)
		error("Must be at least one input and one output file specified");

	ninps--;	/* Number of inputs */

	/* Open and read each input file, and create output file */
	for (n = 0; n <= ninps; n++) {

		if ((inps[n].c = new_cgats()) == NULL)
			error("Failed to create cgats object for file '%s'",inps[n].name);
	
		if (n < ninps) {		/* If input file, read it */
			inps[n].c->add_other(inps[n].c, ""); 	/* Allow any signature file */
		
			if (inps[n].c->read_name(inps[n].c, inps[n].name))
				error("CGATS file '%s' read error : %s",inps[n].name,inps[n].c->err);
		
			if (inps[n].c->ntables < 1)
				error ("Input file '%s' doesn't contain at least one table",inps[n].name);
		}
	}
	ocg = inps[ninps].c;		/* Alias for output file */

	/* Duplicate everything from the first input file into the output file. */
	for (n = 0; n < inps[0].c->ntables; n++) {

		if (inps[0].c->t[n].tt == cgats_X) {
			ocg->add_other(ocg, inps[0].c->cgats_type);
			ocg->add_table(ocg, tt_other, 0);
		} else if (inps[0].c->t[n].tt == tt_other) {
			int oi;
			oi = ocg->add_other(ocg, inps[0].c->others[inps[0].c->t[n].oi]);
			ocg->add_table(ocg, tt_other, oi);
		} else {
			ocg->add_table(ocg, inps[0].c->t[n].tt, 0);
		}

		/* Duplicate all the keywords */
		for (i = 0; i < inps[0].c->t[n].nkwords; i++) {
//printf("~1 table %d, adding keyword '%s'\n",n,inps[0].c->t[n].ksym[i]);
			ocg->add_kword(ocg, n, inps[0].c->t[n].ksym[i], inps[0].c->t[n].kdata[i], NULL);
		}

		/* Duplicate all of the fields */
		for (i = 0; i < inps[0].c->t[n].nfields; i++) {
			ocg->add_field(ocg, n, inps[0].c->t[n].fsym[i], inps[0].c->t[n].ftype[i]);
		}

		/* If more than one file, must be merging or averaging between files  */
		if (ninps > 1) {
			/* Duplicate all of the data or first file to output file */
			if ((setel = (cgats_set_elem *)malloc(
			     sizeof(cgats_set_elem) * inps[0].c->t[n].nfields)) == NULL)
				error("Malloc failed!");
	
			for (i = 0; i < inps[0].c->t[n].nsets; i++) {
				inps[0].c->get_setarr(inps[0].c, n, i, setel);
				ocg->add_setarr(ocg, n, setel);
			}
			free(setel);
		}
	}

	/* Figure out the indexes of the device channels */
	if (inps[0].c->find_kword(inps[0].c, 0, "COLOR_REP") < 0) {
		warning("Input file '%s' doesn't contain keyword COLOR_REP", inps[0].name);
	} else {
		int ti;
		char *buf;
		char *outc;
		int nmask;
		char *bident;

		if ((ti = inps[0].c->find_kword(inps[0].c, 0, "COLOR_REP")) < 0)
			error("Input file '%s' doesn't contain keyword COLOR_REP", inps[0].name);
		
		if ((buf = strdup(inps[0].c->t[0].kdata[ti])) == NULL)
			error("Malloc failed");
		
		/* Split COLOR_REP into device and PCS space */
		if ((outc = strchr(buf, '_')) == NULL)
			error("COLOR_REP '%s' invalid", inps[0].c->t[0].kdata[ti]);
		*outc++ = '\000';
		
		if ((nmask = icx_char2inkmask(buf)) == 0) {

			/* Hmm. This might be an input reference. */
			if (strcmp(buf, "XYZ") == 0 || strcmp(buf, "LAB") == 0) {
				if ((nmask = icx_char2inkmask(outc)) == 0) {
					error ("File '%s' keyword COLOR_REP has unknown device value '%s'",inps[0].name,outc);
				}
			} else {
				error ("File '%s' keyword COLOR_REP has unknown device value '%s'",inps[0].name,buf);
			}
		}

		nchan = icx_noofinks(nmask);
		bident = icx_inkmask2char(nmask, 0); 

		/* Find device fields */
		for (j = 0; j < nchan; j++) {
			int ii, imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : bident,
			                      icx_ink2char(imask));

			if ((ii = inps[0].c->find_field(inps[0].c, 0, fname)) < 0)
				error ("Input file doesn't contain field %s",fname);
			if (inps[0].c->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type",fname);
			chix[j] = ii;
		}
	}
	
	/* Figure out the indexes of the PCS channels, if any */
	{
		int npcs;
		char *fname[2][3] =  { { "LAB_L", "LAB_A", "LAB_B" },
		                       { "XYZ_X", "XYZ_Y", "XYZ_Z" } };

		/* For Lab and XYZ */
		for (j = 0; j < 2; j++) {
			for (npcs = 0; npcs < 3; npcs++) {
				int ii;
	
				if ((ii = inps[0].c->find_field(inps[0].c, 0, fname[j][npcs])) < 0)
					break;		/* Try next or give up */
	
				if (inps[0].c->t[0].ftype[ii] != r_t)
					error ("Field %s is wrong type",fname[j][npcs]);
				pcsix[j][npcs] = ii;
			}
			if (npcs == 3)
				haspcs[j] = 1;
		}
	}
	if (!haspcs[0] && !haspcs[1])
		warning("No PCS fields found - hope that's OK!");
	
	if (!domerge && verb) {
		printf("Averaging the following fields:");
		for (j = 0; j < inps[0].c->t[0].nfields; j++) {
			int jj;

			/* Only real types */
			if (inps[0].c->t[0].ftype[j] != r_t) {
				continue;
			}

			/* Not device channels */
			for (jj = 0; jj < nchan; jj++) {
				if (chix[jj] == j)
					break;
			}
			if (jj < nchan) {
				continue;
			}

			printf(" %s",inps[0].c->t[0].fsym[j]);
		}
		printf("\n");
	}

	/* Get ready to add more values to output, */
	/* for merging or averaging within one file. */
	if ((setel = (cgats_set_elem *)malloc(
	     sizeof(cgats_set_elem) * inps[0].c->t[0].nfields)) == NULL)
		error("Malloc failed!");

	/* If averaging values within the one file */
	if (ninps == 1) {
		int *valdone;
		int npat;
		double *vlist;
		double (*v3list)[3] = NULL;
		int k, e;
		n = 0;		/* Output set index */

		if ((valdone = (int *)calloc(inps[0].c->t[0].nsets, sizeof(int))) == NULL) 
			error("Malloc failed!");

		if ((vlist = (double *)calloc(inps[0].c->t[0].nsets, sizeof(double))) == NULL) 
			error("Malloc failed!");

		if (dogeom && (haspcs[0] || haspcs[1])) {
			if ((v3list = (double (*)[3])calloc(inps[0].c->t[0].nsets, 3 * sizeof(double))) == NULL) 
				error("Malloc failed!");
		}

		/* For each patch */
		for (i = 0; i < inps[0].c->t[0].nsets; i++) {

			if (valdone[i])
				continue;

			inps[0].c->get_setarr(inps[0].c, 0, i, setel);
			ocg->add_setarr(ocg, 0, setel);

			/* For each non-device real field values */
			for (j = 0; j < inps[0].c->t[0].nfields; j++) {
				int jj;

				/* Only real types */
				if (inps[0].c->t[0].ftype[j] != r_t)
					continue;

				/* Not device channels */
				for (jj = 0; jj < nchan; jj++) {
					if (chix[jj] == j)
						break;
				}
				if (jj < nchan)
					continue;

				/* Locate any patches (including starting patch) with matching device values */
				npat = 0;
				for (k = i; k < inps[0].c->t[0].nsets; k++) {
	
					/* Check if the device values match */
					for (e = 0; e < nchan; e++) {
						double diff;
	
						diff = *((double *)inps[0].c->t[0].fdata[i][chix[e]])
						     - *((double *)inps[0].c->t[0].fdata[k][chix[e]]);
	
						if (fabs(diff) > 0.001) {
							break;
						}
					}
					if (e < nchan) {
						continue;
					}

					vlist[npat++] = *((double *)inps[0].c->t[0].fdata[k][j]);
					valdone[k] = 1;
				}
				if (domedian)
					*((double *)ocg->t[0].fdata[n][j]) = median(vlist, npat);
				else
					*((double *)ocg->t[0].fdata[n][j]) = average(vlist, npat);
			}

			/* Override per-component average/median if PCS Geometric Median */
			if (dogeom && (haspcs[0] || haspcs[1])) {
				double res[3];

				/* For Lab and XYZ */
				for (j = 0; j < 2; j++) {

					if (haspcs[j] == 0)
						continue;

					/* Locate any patches (including starting patch) with matching device values */
					npat = 0;
					for (k = i; k < inps[0].c->t[0].nsets; k++) {
		
						/* Check if the device values match */
						for (e = 0; e < nchan; e++) {
							double diff;
		
							diff = *((double *)inps[0].c->t[0].fdata[i][chix[e]])
							     - *((double *)inps[0].c->t[0].fdata[k][chix[e]]);
		
							if (fabs(diff) > 0.001) {
								break;
							}
						}
						if (e < nchan) {
							continue;
						}
	
						v3list[npat][0] = *((double *)inps[0].c->t[0].fdata[k][pcsix[j][0]]);
						v3list[npat][1] = *((double *)inps[0].c->t[0].fdata[k][pcsix[j][1]]);
						v3list[npat][2] = *((double *)inps[0].c->t[0].fdata[k][pcsix[j][2]]);

						if (j == 0 && dogeom == 3)		/* Lab and want XYZ */
							icmLab2XYZ(&icmD50_100, v3list[npat], v3list[npat]);
						else if (j == 1 && dogeom == 2)	/* XYZ and want Lab */
							icmXYZ2Lab(&icmD50_100, v3list[npat], v3list[npat]);

						npat++;
					}
					geommed(res, v3list, npat);
	
					if (j == 0 && dogeom == 3)
						icmXYZ2Lab(&icmD50_100, res, res);
					else if (j == 1 && dogeom == 2)
						icmLab2XYZ(&icmD50_100, res, res);

					*((double *)ocg->t[0].fdata[n][pcsix[j][0]]) = res[0];
					*((double *)ocg->t[0].fdata[n][pcsix[j][1]]) = res[1];
					*((double *)ocg->t[0].fdata[n][pcsix[j][2]]) = res[2];
				}
			}
			n++;		/* One more output set */
		}

		if (v3list != NULL)
			free(v3list);
		free(vlist);
		free(valdone);

	/* Averaging patches between identical files, */
	/* or concatenating (merging) several files */
	} else {

		/* Check/process all the other input files */
		for (n = 1; n < ninps; n++) {
	
			/* Check all the fields match the first file */
			if (inps[0].c->t[0].nfields != inps[n].c->t[0].nfields)
				error ("File '%s' has %d fields, file '%s has %d",
				       inps[n].name, inps[n].c->t[0].nfields, inps[0].name, inps[0].c->t[0].nfields);
			for (j = 0; j < inps[0].c->t[0].nfields; j++) {
				if (inps[0].c->t[0].ftype[j] != inps[n].c->t[0].ftype[j])
					error ("File '%s' field no. %d named '%s' doesn't match file '%s' field '%s'",
					       inps[n].name, j, inps[n].c->t[0].fsym[j], inps[0].name, inps[0].c->t[0].fsym[j]);
			}
	
			/* If merging, append all the values */
			if (domerge) {
				for (i = 0; i < inps[n].c->t[0].nsets; i++) {
					inps[n].c->get_setarr(inps[n].c, 0, i, setel);
					ocg->add_setarr(ocg, 0, setel);
				}
	
			} else {	/* Averaging */

				/* Check the number of patches matches the first file */
				if (inps[0].c->t[0].nsets != inps[n].c->t[0].nsets)
					error ("File '%s' has %d sets, file '%s has %d",
					       inps[n].name, inps[n].c->t[0].nsets, inps[0].name, inps[0].c->t[0].nsets);
				/* For all the patches: */
				for (i = 0; i < inps[n].c->t[0].nsets; i++) {
	
					/* Check that the device values match the first file */
					for (j = 0; j < nchan; j++) {
						double diff;
						diff = *((double *)inps[0].c->t[0].fdata[i][chix[j]])
						     - *((double *)inps[n].c->t[0].fdata[i][chix[j]]);
	
						if (fabs(diff) > 0.001)
							error ("File '%s' set %d has field '%s' value that differs from '%s'",
					       inps[n].name, i+1, inps[n].c->t[0].fsym[j], inps[0].name);
					}
				}
			}
		}
		
		/* If averaging */
		if (!domerge) {
			int npat;
			double *vlist;
			double (*v3list)[3] = NULL;

			if ((vlist = (double *)calloc(ninps, sizeof(double))) == NULL) 
				error("Malloc failed!");

			if (dogeom && (haspcs[0] || haspcs[1])) {
				if ((v3list = (double (*)[3])calloc(inps[0].c->t[0].nsets, 3 * sizeof(double))) == NULL) 
					error("Malloc failed!");
			}

			/* For all the non-device real field values */
			for (j = 0; j < inps[0].c->t[0].nfields; j++) {
				int jj;

				/* Only real types */
				if (inps[0].c->t[0].ftype[j] != r_t)
					continue;

				/* Not device channels */
				for (jj = 0; jj < nchan; jj++) {
					if (chix[jj] == j)
						break;
				}
				if (jj < nchan)
					continue;
	
				/* For each patch */
				for (i = 0; i < inps[n].c->t[0].nsets; i++) {

					/* For all input files */
					npat = 0;
					for (n = 0; n < ninps; n++) {
						vlist[npat++] = *((double *)inps[n].c->t[0].fdata[i][j]);
					}

					if (domedian)
						*((double *)ocg->t[0].fdata[i][j]) = median(vlist, npat);
					else
						*((double *)ocg->t[0].fdata[i][j]) = average(vlist, npat);
				}
			}

			/* Override per-component average/median if PCS Geometric Median */
			if (dogeom && (haspcs[0] || haspcs[1])) {
				double res[3];

				/* For each patch */
				for (i = 0; i < inps[n].c->t[0].nsets; i++) {

					/* For Lab and XYZ */
					for (j = 0; j < 2; j++) {

						if (haspcs[j] == 0)
							continue;
	
						/* For all input files */
						npat = 0;
						for (n = 0; n < ninps; n++) {
							v3list[npat][0] = *((double *)inps[n].c->t[0].fdata[i][pcsix[j][0]]);
							v3list[npat][1] = *((double *)inps[n].c->t[0].fdata[i][pcsix[j][1]]);
							v3list[npat][2] = *((double *)inps[n].c->t[0].fdata[i][pcsix[j][2]]);

							if (j == 0 && dogeom == 3)		/* Lab and want XYZ */
								icmLab2XYZ(&icmD50_100, v3list[npat], v3list[npat]);
							else if (j == 1 && dogeom == 2)	/* XYZ and want Lab */
								icmXYZ2Lab(&icmD50_100, v3list[npat], v3list[npat]);

							npat++;
						}
						geommed(res, v3list, npat);

						if (j == 0 && dogeom == 3)
							icmXYZ2Lab(&icmD50_100, res, res);
						else if (j == 1 && dogeom == 2)
							icmLab2XYZ(&icmD50_100, res, res);

						*((double *)ocg->t[0].fdata[i][pcsix[j][0]]) = res[0];
						*((double *)ocg->t[0].fdata[i][pcsix[j][1]]) = res[1];
						*((double *)ocg->t[0].fdata[i][pcsix[j][2]]) = res[2];
					}
				}
			}

			if (v3list != NULL)
				free(v3list);
			free(vlist);
		}
	}

	/* Write out the output and free the cgats * */
	for (n = 0; n <= ninps; n++) {

		if (n >= ninps) {		/* If ouput file, write it */
			if (inps[n].c->write_name(inps[n].c, inps[ninps].name))
				error("CGATS file '%s' write error : %s",inps[n].name,inps[n].c->err);
		}
		inps[n].c->del(inps[n].c);
	}

	free(setel);
	free(inps);
	
	return 0;
}


static double average(double *vals, int nvals) {
	double rv;
	int i;

	for (rv = 0.0, i = 0; i < nvals; i++)
		rv += vals[i];

	if (nvals > 0)
		rv /= (double)nvals;

	return rv;
}

static double median(double *vals, int nvals) {
	if (nvals < 3)
		return average(vals, nvals);

#define HEAP_COMPARE(A,B) (A < B)
	HEAPSORT(double,vals,nvals);

	if ((nvals & 1) != 0)
		return vals[nvals/2]; 
	else
		return 0.5 * (vals[nvals/2] + vals[nvals/2-1]); 
}

/* Compute Geometric Median of PCS values */
/* using Weiszfeld's algorithm. */
static void geommed(double res[3], double vals[][3], int nvals) {
	int i, j;

	/* Start with mean value */
	icmSet3(res, 0.0);
	for (i = 0; i < nvals; i++)
		icmAdd3(res, res, vals[i]);
	icmScale3(res, res, 1.0/(double)nvals);

//printf("\n~1 average = %f %f %f\n", res[0], res[1], res[2]);

	/* Itterate to approach Geometric Mean */
	for (j = 0; j < 20; j++) {
        double tv[3], tt;
		int k;

        icmSet3(tv, 0.0);
		tt = 0.0;
		for (k = i = 0; i < nvals; i++) {
			double norm = icmNorm33(vals[i], res);
			if (norm < 1e-6)
				continue;
            tv[0] += vals[i][0]/norm;
            tv[1] += vals[i][1]/norm;
            tv[2] += vals[i][2]/norm;
            tt += 1.0/norm;
			k++;
//printf("Norm = %f, tv = %f %f %f, tt = %f\n",norm, tv[0], tv[1], tv[2], tt);
        }
		if (k > 0)
			icmScale3(res, tv, 1.0/tt);
//printf("~1 res = %f %f %f\n", res[0], res[1], res[2]);
	}

//printf("~1 geomm = %f %f %f\n", res[0], res[1], res[2]);
}


