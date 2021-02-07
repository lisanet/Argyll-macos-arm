
/* 
 * Argyll Color Correction System
 *
 * Read in the RGB & CIE data from a LightSpace format .bcs XML file
 * and convert it into a .ti3 CGATs format suitable for the Argyll CMS.
 *
 * Derived from  txt2cgats.c 
 * Author: Graeme W. Gill
 * Date:   16/11/00
 *
 * Copyright 2000 - 2014, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD
	
 */

#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include "copyright.h"
#include "aconfig.h"
#include "cgats.h"
#include "xspect.h"
#include "insttypes.h"
#include "mxml.h"
#include "numlib.h"
#include "ui.h"

#define DEB 6

void
usage(char *mes) {
	fprintf(stderr,"Convert LightSpace raw RGB device profile data to Argyll CGATS data, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (mes != NULL)
		fprintf(stderr,"error: %s\n",mes);
	fprintf(stderr,"usage: ls2ti3 [-v] infile outbase\n");
/*	fprintf(stderr," -v            Verbose mode\n"); */
	fprintf(stderr," infile        Input LightSpace .bcs file\n");
	fprintf(stderr," outbasename   Output file basename for .ti3\n");
	exit(1);
	}

/* XML data type callback for mxmlLoadFile() */
mxml_type_t
type_cb(mxml_node_t *node) {
	mxml_node_t *parent = mxmlGetParent(node);
	const char *pname;
	const char *name = node->value.element.name;

	if (parent == NULL)
		return MXML_TEXT;

	pname = parent->value.element.name;

//	printf("~1 type_cb got pnode '%s' node '%s'\n",pname, name);

	if (strcmp(pname, "stimuli") == 0
	 && (strcmp(name, "red") == 0
	  || strcmp(name, "green") == 0
	  || strcmp(name, "blue") == 0))
		return MXML_REAL;

	if (strcmp(pname, "XYZ") == 0
	 && (strcmp(name, "X") == 0
	  || strcmp(name, "Y") == 0
	  || strcmp(name, "Z") == 0))
		return MXML_REAL;

#ifdef NEVER
	if (strcmp(pname, "FileInformation") == 0
	 && (strcmp(name, "Creator") == 0
	  || strcmp(name, "CreationDate") == 0
	  || strcmp(name, "Description") == 0))
		return MXML_OPAQUE;		/* Don't split strings up */
#endif

	return MXML_TEXT;
}

struct _patch {
	int pno;				/* Patch number */

	double dev[MAX_CHAN];	/* Device value */
	double XYZ[3];			/* XYZ value */

}; typedef struct _patch patch;


int main(int argc, char *argv[]) {
	int i, j;
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	static char inname[MAXNAMEL+1] = { 0 };		/* Input LightSpace .xml file */
	static char outname[MAXNAMEL+1] = { 0 };		/* Output cgats .ti3 file base name */
	double dev_scale = 1.0;	/* Device value scaling */

    FILE *ifp;
    mxml_node_t *tree, *pnode, *node;
    mxml_node_t *top, *data;
	const char *attr, *name;
	patch *patches;

	cgats *ocg;				/* output cgats structure for .ti3 */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */

	int islab = 0;			/* CIE is Lab rather than XYZ */

	int npat = 0;			/* Number of patches */
	int ndchan = 3;			/* Number of device channels, 0 = no device, RGB = 3, CMYK = 4 */
	cgats_set_elem *setel;	/* Array of set value elements */

	error_program = "ls2ti3";

	if (argc < 3)
		usage("Too few arguments");

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

			if (argv[fa][1] == '?')
				usage(NULL);

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;
			else 
				usage("Unknown flag");
		} else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage("input filename not found");
	strncpy(inname,argv[fa++],MAXNAMEL); inname[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage("output basename not found");
	strncpy(outname,argv[fa++],MAXNAMEL-4); outname[MAXNAMEL-4] = '\000';
	strcat(outname,".ti3");

    if ((ifp = fopen(inname, "r")) == NULL)
		error("Unable to read '%s'",inname);

#ifndef NEVER
    tree = mxmlLoadFile(NULL, ifp, type_cb);
#else
	tree = mxmlSAXLoadFd(NULL, fileno(ifp), type_cb, sax_cb, (void *)p);
#endif
    fclose(ifp);

	if (tree == NULL)
		error("Parsing '%s' failed",inname);

	if ((top = mxmlFindElement(tree, tree, NULL, NULL, NULL, MXML_DESCEND_FIRST)) == NULL
	 || mxmlGetType(top) != MXML_ELEMENT)
		error("Failed to find top element in '%s'",inname);

	if (strcmp(top->value.element.name, "builder_color_space") != 0)
		error("'%s' doesn't seem to be a LightSpace .bcs file ?",inname);

#ifdef NEVER
	/* Check that its CxF3 */
	if ((attr = mxmlElementGetAttr(top, "xmlns")) == NULL)
		error("Failed to find CxF attribute %s in '%s'","xmlns", inname);

	if (strcmp(attr, "http://colorexchangeformat.com/CxF3-core") != 0)
		error("File '%s' is not CxF format",inname);

	/* Look for header information */
	if ((pnode = mxmlFindElement(top, top, "head", NULL, NULL, MXML_DESCEND_FIRST)) == NULL)
		error("Failed to find head in '%s'",inname);

	/* Grab the creation date */
	if ((node = mxmlFindElement(pnode, pnode, "created", NULL, NULL, MXML_DESCEND_FIRST)) == NULL)
		error();

	name = mxmlGetOpaque(node);
#endif

	/* Locate the data node */
	if ((data = mxmlFindElement(top, top, "data", NULL, NULL, MXML_DESCEND_FIRST)) == NULL)
		error("Failed to find data in '%s'",inname);

	/* Get the number of patches */
	if ((attr = mxmlElementGetAttr(data, "frames")) == NULL)
		error("Failed to find data attribute 'frames' in '%s'", inname);

	npat = atoi(attr);

	if (npat <= 0)
		error("Illegal number of patches %d in '%s",inname);

	a1logd(g_log, DEB, "npat = %d\n",npat);

	if ((patches = calloc(npat, sizeof(patch))) == NULL)
		error("malloc failed");

	/* Read all the patches */
	node = mxmlFindElement(data, data, "patch", NULL, NULL, MXML_DESCEND_FIRST);
	for (i = 0; node != NULL && i < npat;) {
		int j;
	    mxml_node_t *stim, *res, *xyz, *val;
		char *rgb_key[3] = { "red", "green", "blue" };
		char *xyz_key[3] = { "X", "Y", "Z" };
		
		if (mxmlGetType(node) != MXML_ELEMENT) {
			a1logd(g_log, DEB, "skipping non element node type %d\n",mxmlGetType(node));
			goto next;
		}
		a1logd(g_log, DEB, "read node '%s'\n",node->value.element.name);

		if ((attr = mxmlElementGetAttr(node, "frame")) == NULL) {
			a1logd(g_log, DEB, "read_cxf: skipping node without frame\n");
			goto next;							/* Skip this one */
		}
		patches[i].pno = atoi(attr);

		a1logd(g_log, DEB, "got patch %d no %d\n",i,patches[i].pno);

		/* Read the RGB stimulus */
		if ((stim = mxmlFindElement(node, node, "stimuli", NULL, NULL, MXML_DESCEND_FIRST))
			                                                                       == NULL) {
			a1logd(g_log, DEB, "Can't find 'stimuli' in patch %d - skipping\n",patches[i].pno);
			goto next;
		}

		for (j = 0; j < 3; j++) {
			if ((val = mxmlFindElement(stim, stim, rgb_key[j], NULL, NULL, MXML_DESCEND_FIRST))
				                                                                       == NULL) {
				a1logd(g_log, DEB, "failed to find RGB component %s\n",rgb_key[j]);
				goto next;
			}
			patches[i].dev[j] = mxmlGetReal(val);
			a1logd(g_log, DEB, "got RGB component %s value %f\n",rgb_key[j],patches[i].dev[j]);
		}

		/* Read the XYZ results */
		if ((res = mxmlFindElement(node, node, "results", NULL, NULL, MXML_DESCEND_FIRST))
			                                                                       == NULL) {
			a1logd(g_log, DEB, "Can't find 'results' in patch %d - skipping\n",patches[i].pno);
			goto next;
		}

		if ((xyz = mxmlFindElement(res, res, "XYZ", NULL, NULL, MXML_DESCEND_FIRST))
			                                                                       == NULL) {
			a1logd(g_log, DEB, "Can't find 'XYZ' in patch %d - skipping\n",patches[i].pno);
			goto next;
		}

		for (j = 0; j < 3; j++) {
			if ((val = mxmlFindElement(xyz, xyz, xyz_key[j], NULL, NULL, MXML_DESCEND_FIRST))
				                                                                       == NULL) {
				a1logd(g_log, DEB, "failed to find XYZ component %s\n",xyz_key[j]);
				goto next;
			}
			patches[i].XYZ[j] = mxmlGetReal(val);
			a1logd(g_log, DEB, "got XYZ component %s value %f\n",xyz_key[j],patches[i].XYZ[j]);
		}
		i++;

		/* Add patch to output CGATS */

	next:;	/* Next color */
		node = mxmlGetNextSibling(node);
	}

	mxmlDelete(tree);


	/* Setup output cgats file */
	ocg = new_cgats();	/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI3"); 	/* our special type is Calibration Target Information 3 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll target", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);	/* What sort of device this is */

	/* Note what instrument the chart was read with */
	/* Assume this - could try reading from file INSTRUMENTATION "SpectroScan" ?? */
	ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", inst_name(instSpectrolino) , NULL);

	/* Don't make any assumptions about normalisation (this is default anyway) */
	ocg->add_kword(ocg, 0, "NORMALIZED_TO_Y_100","NO", NULL);

	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);

	if (ndchan == 3) {
		ocg->add_field(ocg, 0, "RGB_R", r_t);
		ocg->add_field(ocg, 0, "RGB_G", r_t);
		ocg->add_field(ocg, 0, "RGB_B", r_t);
		if (islab)
			ocg->add_kword(ocg, 0, "COLOR_REP","RGB_LAB", NULL);
		else
			ocg->add_kword(ocg, 0, "COLOR_REP","RGB_XYZ", NULL);
	}

	if (islab) {
		ocg->add_field(ocg, 0, "LAB_L", r_t);
		ocg->add_field(ocg, 0, "LAB_A", r_t);
		ocg->add_field(ocg, 0, "LAB_B", r_t);
	} else {
		ocg->add_field(ocg, 0, "XYZ_X", r_t);
		ocg->add_field(ocg, 0, "XYZ_Y", r_t);
		ocg->add_field(ocg, 0, "XYZ_Z", r_t);
	}

	/* LS uses 0.0 .. 1.0 */
	dev_scale = 100.0/1.0;

	/* Write out the patch info to the output CGATS file */
	if ((setel = (cgats_set_elem *)malloc(
	     sizeof(cgats_set_elem) * ocg->t[0].nfields)) == NULL)
		error("Malloc failed!");

	/* Write out the patch info to the output CGATS file */
	for (i = 0; i < npat; i++) {
		char id[100];
		int k = 0;

		/* SAMPLE ID */
		sprintf(id, "%d", patches[i].pno);
		setel[k++].c = id;

		if (ndchan == 3) {
			setel[k++].d = dev_scale * patches[i].dev[0];
			setel[k++].d = dev_scale * patches[i].dev[1];
			setel[k++].d = dev_scale * patches[i].dev[2];
		} 

		setel[k++].d = patches[i].XYZ[0];
		setel[k++].d = patches[i].XYZ[1];
		setel[k++].d = patches[i].XYZ[2];

		ocg->add_setarr(ocg, 0, setel);
	}

	free(setel);

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	/* Clean up */
	ocg->del(ocg);

	return 0;
}



