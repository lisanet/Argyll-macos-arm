
/* Integer Multi-Dimensional Interpolation */

/*
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licensed under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * Top level kernel code generator
 *
 * This module is invoked from the make system,
 * and generates all the versions and configurations of 
 * the IMDI kernel code. It includes all the generated
 * files in imdi_k.h, which also contains a table
 * so that the run time code knows what kernels
 * are available.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include "copyright.h"
#include "aconfig.h"

#include "imdi.h"
#include "imdi_tab.h"

#ifndef MAXNAMEL
# define MAXNAMEL 512	/* Maximum command line filename lengths */
#endif

#undef VERBOSE
#undef TEST1		/* Generate one test case */

/*
	Ideal grid resolutions for 8 bit precision calculations.
	See imdi_gen.c for a more detailed list.

	Grid	Jumps	
	4		0
	6		0
	16		0
	18		0
	52		0
	86		0
	256		0

	3		1
	5		1
	9		1
	17		1
	33		1
	65		1
	128		1
	129		1
	255		1

 */

/* The following structure initialisations define what kernel routines should be built */
static
gendesc descs[] = {
#ifdef TEST1
	{
		{ 3, 0 },					/* * Input dimension combinations */
		{ 33, 0 },	 				/* + Interpolation table resolutions */
		{ 8, 0 },	 				/* + Min Simplex table resolutions */
		{ 4, 0 },					/* * Output dimension combinations */
		{OOPT(oopts_check,3), oopts_none},	/* + Output channel options */
		{pixint8, 0 },				/* * Input pixel representation */
		{prec_p16, 0},				/* + Internal precision */
		{pixint8, 0},				/* + Output pixel representation */
		{opts_sort_splx, opts_end}		/* * Direction & stride combinations */
//		{opts_splx_sort, opts_end}		/* * Direction & stride combinations */
	}
#else
	/* A reasonably full set of combinations */
								/* * means multiplies combination */
								/* + means lockstep with previous line */
	{
		{ 1,    3,  4,  5,  6,  7,  8, 9, 10, 0 },	/* * Input dimension combinations */
		{ 256, 33, 18, 16, 12,  8,  7, 6, 5,  0 }, 	/* + Min Interpolation table resolutions */
		{ 1,    8, 17,  1,  1,  1,  1, 1, 1,  0 }, 	/* + Min Simplex table resolutions */

		{ 1, 3, 4, 5, 6, 7, 8, 9, 10, 0 },				/* * Output dimension combinations */
		{oopts_none, oopts_none, oopts_none, oopts_none, oopts_none, oopts_none,
		 oopts_none, oopts_none, oopts_none, oopts_none, oopts_none, oopts_none},
														/* + Output channel options */

		{pixint8, pixint16, pixint8,  pixint16, pixint16, 0 },	/* * Input pixel representation */
		{prec_p8, prec_p8,  prec_p8,  prec_p16, prec_p16, 0 },	/* + Internal precision */
		{pixint8, pixint8,  pixint16, pixint16, pixint16, 0 },	/* + Output pixel representation */

		{
		 opts_splx_sort,						/* (both, but default to simple alg, no stride) */
		 opts_istride | opts_ostride,			/* + (Sort only with stride) */
		 opts_end }								/* * Direction & stride combinations */
	}
#endif	/* !TEST1 */
};

void set_architecture(mach_arch *ar, int use64);

struct _knamestr {
	char name[100];
	char desc[100];
	struct _knamestr *next;
}; typedef struct _knamestr knamestr;

knamestr *
new_knamestr(char *name, char *desc) {
	knamestr *kn;
	
	if ((kn = (knamestr *)malloc(sizeof(knamestr))) == NULL) {
		fprintf(stderr,"new_knamestr malloc failed\n");
		exit(-1);
	}
	strcpy(kn->name, name);
	strcpy(kn->desc, desc);
	kn->next = NULL;
	return kn;
}

void usage(void) {
	fprintf(stderr,"Make imdi kernel code Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"usage: imdi_make [-i]\n");
	fprintf(stderr," -d dir        Directory to create them in (default .)\n");
	fprintf(stderr," -i            Individial Files\n");
	fprintf(stderr," -f            Force 64 bit\n");
	exit(1);
}

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	int indiv = 0;			/* Individual files */
	int rv;
	int dn, tnd;
	genspec gs, ogs;
	tabspec ts, ots;
	mach_arch ar;
	int ix = 1;				/* kernel index */
	knamestr *list = NULL, *lp = NULL;
#if defined(ALLOW64) && defined(USE64)
	int use64 = 1;
#else
	int use64 = 0;
#endif
	char dirname[MAXNAMEL+1+1] = "";   /* Output directory name */
	char temp[MAXNAMEL+100+1];		/* Buffer to compose filenames in */
	FILE *kcode = NULL;	/* Kernel routine code file */
	FILE *kheader;		/* Kernel routine header file */

	/* Zero out the gen and tabspecs, to give diff a place to start */
	memset((void *)&ogs, 0, sizeof(genspec));
	memset((void *)&gs, 0, sizeof(genspec));
	memset((void *)&ots, 0, sizeof(tabspec));
	memset((void *)&ts, 0, sizeof(tabspec));

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

			if (argv[fa][1] == '?') {
				usage();

			}
			/* Destination directory */
			else if (argv[fa][1] == 'd' || argv[fa][1] == 'D') {
				int len;
				fa = nfa;
				if (na == NULL) usage();
				strncpy(dirname,na,MAXNAMEL); dirname[MAXNAMEL] = '\000';
				len = strlen(dirname);
				if (len > 0) {
					if (dirname[len-1] != '/' && dirname[len-1] != '\\')
						strcat(dirname, "/");
				}
			}
			/* Individual files */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				indiv = 1;
			}
			/* Force 64 bit */
			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F') {
#ifdef ALLOW64
				use64 = 1;
#else
				fprintf(stderr,"ALLOW64 bits is undefined\n");
				usage();
#endif
			}
			else { 
				usage();
			}
		} else
			break;
	}

	set_architecture(&ar, use64);

	/* Open the file for kernel routine declaration header */
	sprintf(temp, "%simdi_k.h",dirname);
	if ((kheader = fopen(temp, "w")) == NULL) {
		fprintf(stderr,"imdi_make: unable to open file '%s'\n",temp);
		exit(-1);
	}

	if (!indiv) {
		sprintf(temp, "%simdi_k.c",dirname);
		if ((kcode = fopen(temp, "w")) == NULL) {
			fprintf(stderr,"imdi_make: unable to open file '%s'\n",temp);
			exit(-1);
		}
	}

	tnd = sizeof(descs)/sizeof(gendesc);	/* Total number of descriptions */
#ifdef VERBOSE
	printf("Number of descriptions = %d\n",tnd);
#endif /* VERBOSE */

	fprintf(kheader,"/* Integer Multi-Dimensional Interpolation */\n");
	fprintf(kheader,"/* Declarations for all the generated kernel functions */\n");
	fprintf(kheader,"/* This file is generated by imdi_make */\n\n");
	fprintf(kheader,"/* Copyright 2000 - 2007 Graeme W. Gill */\n");
	fprintf(kheader,"/* All rights reserved. */\n");
	fprintf(kheader,"/* This material is licensed under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :- */\n");
	fprintf(kheader,"/* see the License.txt file for licensing details.*/\n");

	fprintf(kheader,"\n");

	/* For all the descriptions */
	for (dn = 0; dn < tnd; dn++) {
		int cb, ncb;

		/* Do all combinations for that description */
		for (cb = 0, ncb = 1; cb < ncb; cb++) {
			int nalg, alg;
			char ofname[100];

			/* Compute generate spec. and number of combinations */
			ncb = set_genspec(&gs, &descs[dn], cb, &ar);

			if (indiv) {
				sprintf(temp, "%s%s.c",dirname,ofname);
				if ((kcode = fopen(temp, "w")) == NULL) {
					fprintf(stderr,"imdi_make: unable to open file '%s'\n",temp);
					exit(-1);
				}
			}

			nalg = 2;		/* By default generate just sort algorithm */
			alg = 1;

			if ((gs.opt & opts_splx_sort)
			 || (gs.opt & opts_sort_splx)) {
				alg = 0;	/* Generate both simplex and sort algorithms */			
			}
			if (gs.opt & opts_splx) {
				nalg = 1;	/* Generate just simplex algorithm */
				alg = 0;
			}

			for (; alg < nalg; alg++) {

				if (alg == 0)
					gs.opt |=  opts_splx;
				else
					gs.opt &=  ~opts_splx;

				/* Generate it */
				rv = gen_c_kernel(&gs, &ts, &ar, kcode, ix, &ogs, &ots);
				if (rv != 0 && rv != 1) {
					fprintf(stderr,"imdi_make: gen_c_kernel returned a err %d\n",rv);
					exit(-1);
				}

				/* Add the name to the list */
				if (list == NULL)
					lp = list = new_knamestr(gs.kname, gs.kdesc);
				else {
					lp->next = new_knamestr(gs.kname, gs.kdesc);
					lp = lp->next;
				}
				if (indiv) {
					if (fclose(kcode) != 0) {
						fprintf(stderr,"imdi_make: unable to close file '%s'\n",ofname);
						exit(-1);
					}
				}
				ogs = gs;	/* Structure copy */
				ots = ts;
				ix++;

				if (rv == 0)
					break;		/* there was only sort available */
			}
		}
	}

	/* Include the kernel functions in the header file */
	if (indiv) {
		for(lp = list; lp != NULL; lp = lp->next) {
			fprintf(kheader,"#include \"%s_%s.c\"\n",lp->name,lp->desc);
		}
	} else {
		fprintf(kheader,"#include \"imdi_k.c\"	/* All the kernel code */\n");
	}
	fprintf(kheader,"\n");

	/* Output function table */
	
	fprintf(kheader,
		"struct {\n"
		"	void (*interp)(imdi *s, void **outp, int ostride, void **inp, int  istride, unsigned int npix);\n"
		"	void (*gentab)(genspec *g, tabspec *t);\n"
		"} ktable[%d] = {\n",ix-1);

	for(lp = list; lp != NULL; lp = lp->next) {
		fprintf(kheader,"\t{ %s, %s_gentab }%s\n", lp->name, lp->name, 
		lp->next != NULL ? "," : "");
	}
	fprintf(kheader,"};\n");
	fprintf(kheader,"\n");
	fprintf(kheader,"int no_kfuncs = %d;\n",ix-1);
	fprintf(kheader,"\n");

	if (!indiv) {
		if (fclose(kcode) != 0) {
			fprintf(stderr,"imdi_make: unable to close file 'imdi_k.c'\n");
			exit(-1);
		}
	}

	if (fclose(kheader) != 0) {
		fprintf(stderr,"imdi_make: unable to close file 'imdi_k.h'\n");
		exit(-1);
	}

	/* Free the kname list */
	for(lp = list; lp != NULL;) {
		char *p = (char *)lp;
		lp = lp->next;
		free(p);
	}

	return 0;
}


/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Initialse the architecture structure properly. */
/* We're doing this purely at run time, on the assumption */
/* that the target machine is the one we're running on. */
/* We would have to do this differently in a cross development */
/* environment. */
void
set_architecture(
mach_arch *ar,
int use64
) {
	unsigned long etest = 0xff;
	char *machtype;		/* Environment value */

	if (*((unsigned char *)&etest) == 0xff) {
		ar->bigend = 0;		/* Little endian */
	} else {
		ar->bigend = 1;		/* Big endian endian */
	}

	machtype = getenv("MACHTYPE");

	/* Unfortunetaly many environments don't export MACHTYPE :-( */
	/* so we implement a fall back */
	if (machtype == NULL) {
#ifdef __ppc__
		machtype = "powerpc";
#endif
	}

	if (machtype != NULL && strcmp(machtype, "powerpc") == 0) {
	
		/* Section tunable for PowerPC */

		ar->uwa    = 0;		/* Use wide memory access */
		ar->shfm   = 0;		/* Use shifts to mask values */
		ar->oscale = 8;		/* Has scaled indexing up to * 8 */
		ar->smmul  = 0;		/* Doesn't have fast small multiply for index scaling */
		if (use64) {
			ar->nords  = 4;		/* Number of ord types */
			ar->nints  = 4;		/* Number of int types */
		} else {
			ar->nords  = 3;		/* Number of ord types */
			ar->nints  = 3;		/* Number of int types */
		}
		ar->natord = 2;		/* Most natural type (assume unsigned int) */
		ar->natint = 2;		/* Most natural type (assume int) */

		ar->pbits = sizeof(void *) * 8;		/* Number of bits in a pointer */

		ar->ords[0].bits = 8 * sizeof(unsigned char);
		ar->ords[0].name = "unsigned char";
		ar->ords[0].align = 1;

		ar->ords[1].bits = 8 * sizeof(unsigned short);
		ar->ords[1].name = "unsigned short";
		ar->ords[1].align = 1;

		ar->ords[2].bits = 8 * sizeof(unsigned int);
		ar->ords[2].name = "unsigned int";
		ar->ords[2].align = 1;

#ifdef ALLOW64
		ar->ords[3].bits = 8 * sizeof(unsigned longlong);
		ar->ords[3].name = "unsigned " str_longlong ;
		ar->ords[3].align = 0;
#endif /* ALLOW64 */

		ar->ints[0].bits = 8 * sizeof(signed char);
		ar->ints[0].name = "signed char";
		ar->ints[0].align = 1;

		ar->ints[1].bits = 8 * sizeof(short);
		ar->ints[1].name = "short";
		ar->ints[1].align = 1;

		ar->ints[2].bits = 8 * sizeof(int);
		ar->ints[2].name = "int";
		ar->ints[2].align = 1;

#ifdef ALLOW64
		ar->ints[3].bits = 8 * sizeof(longlong);
		ar->ints[3].name = str_longlong ;
		ar->ints[3].align = 0;
#endif /* ALLOW64 */

	} else {

		/* Currently assume x86 type */

		ar->uwa    = 0;		/* Use wide memory access */
		ar->shfm   = 0;		/* Use shifts to mask values */
		ar->oscale = 8;		/* Has scaled indexing up to * 8 */
		ar->smmul  = 0;		/* Doesn't have fast small multiply for index scaling */
		if (use64) {
			ar->nords  = 4;		/* Number of ord types */
			ar->nints  = 4;		/* Number of int types */
		} else {
			ar->nords  = 3;		/* Number of ord types */
			ar->nints  = 3;		/* Number of int types */
		}
		ar->natord = 2;		/* Most natural type (assume unsigned int) */
		ar->natint = 2;		/* Most natural type (assume int) */

		ar->pbits = sizeof(void *) * 8;		/* Number of bits in a pointer */

		ar->ords[0].bits = 8 * sizeof(unsigned char);
		ar->ords[0].name = "unsigned char";
		ar->ords[0].align = 1;

		ar->ords[1].bits = 8 * sizeof(unsigned short);
		ar->ords[1].name = "unsigned short";
		ar->ords[1].align = 1;

		ar->ords[2].bits = 8 * sizeof(unsigned int);
		ar->ords[2].name = "unsigned int";
		ar->ords[2].align = 1;

#ifdef ALLOW64
		ar->ords[3].bits = 8 * sizeof(unsigned longlong);
		ar->ords[3].name = "unsigned " str_longlong ;
		ar->ords[3].align = 0;
#endif /* ALLOW64 */

		ar->ints[0].bits = 8 * sizeof(signed char);
		ar->ints[0].name = "signed char";
		ar->ints[0].align = 1;

		ar->ints[1].bits = 8 * sizeof(short);
		ar->ints[1].name = "short";
		ar->ints[1].align = 1;

		ar->ints[2].bits = 8 * sizeof(int);
		ar->ints[2].name = "int";
		ar->ints[2].align = 1;

#ifdef ALLOW64
		ar->ints[3].bits = 8 * sizeof(longlong);
		ar->ints[3].name = str_longlong ;
		ar->ints[3].align = 0;
#endif /* ALLOW64 */
	}
}


