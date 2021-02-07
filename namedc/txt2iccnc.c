
/* 
 * Argyll Color Correction System
 * x11 style .txt to ICC Named Color converter
 *
 * Author: Graeme W. Gill
 * Date:   24/20/2014
 *
 * Copyright 2014 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * TTBD
 */

#undef DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include "copyright.h"
#include "aconfig.h"
#ifndef SALONEINSTLIB
#include "numlib.h"
#include "icc.h"
#else
#include "numsup.h"
#endif
#include "cgats.h"
#include "xspect.h"
#include "ui.h"

void usage(char *diag, ...) {
	fprintf(stderr,"x11 .txt to ICC Named Colors\n");
	fprintf(stderr,"Author: Graeme W. Gill\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: txt2iccNC [-v level] description infile outfile\n");
	exit(1);
}

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	int verb = 0;
	char desc[MAXNAMEL+1];
	char inname[MAXNAMEL+1];
	char outname[MAXNAMEL+1];
	FILE *ifp;
#define BUFSZ 512
	char buf[BUFSZ];
	icmFile *wr_fp;
	icc *wr_icco;
	int i, rv;

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
				usage(NULL);

			else if (argv[fa][1] == 'v')
				verb = 1;

			else 
				usage("Unknown option '%c'",argv[fa][1]);
		}
		else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage("Missing description");
	strncpy(desc,argv[fa++],MAXNAMEL); desc[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage("Missing input filename");
	strncpy(inname,argv[fa++],MAXNAMEL); inname[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage("Missing output filename");
	strncpy(outname,argv[fa++],MAXNAMEL); outname[MAXNAMEL] = '\000';

	if ((ifp = fopen(inname,"r"))==NULL) {
		error("Read: Can't open file '%s'",inname);
	}

	/* Open up the file for writing */
	if ((wr_fp = new_icmFileStd_name(outname,"w")) == NULL)
		error ("Write: Can't open file '%s'",outname);

	if ((wr_icco = new_icc()) == NULL)
		error ("Write: Creation of ICC object failed");

	/* Values that must be set before writing */
	wr_icco->header->deviceClass     = icSigNamedColorClass;
   	wr_icco->header->colorSpace      = icSigRgbData;
   	wr_icco->header->pcs             = icSigLabData;
   	wr_icco->header->renderingIntent = icRelativeColorimetric;

	/* Add the description tag */
	{
		icmTextDescription *wo;
		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigProfileDescriptionTag,	icSigTextDescriptionType)) == NULL) 
			error("Failed to add icmTextDescription");
	
		wo->size = strlen(desc)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, desc);		/* Copy the string in */
	}

	/* Copyright Tag: */
	{
		icmText *wo;
		char *crt;

		crt = "";

		if ((wo = (icmText *)wr_icco->add_tag(
		           wr_icco, icSigCopyrightTag,	icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->data, crt);		/* Copy the text in */
	}

	/* White Point Tag: */
	/* Use the orgininal, non PCS adapted white point */
	{
		icmXYZArray *wo;
		/* Note that tag types icSigXYZType and icSigXYZArrayType are identical */
		if ((wo = (icmXYZArray *)wr_icco->add_tag(
		           wr_icco, icSigMediaWhitePointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0] = icmD65;			/* sRGB is D65 */
	}

	/* Add a named color tag */
	{
		icmNamedColor *wo;

		if ((wo = (icmNamedColor *)wr_icco->add_tag(
	           wr_icco, icSigNamedColor2Tag, icSigNamedColor2Type)) == NULL) 
			error("Adding icmNamedColor tag failed");
	
	   	wo->count = 0;	
	
		/* Read lines and count the colors */
		for (;;) {
			double rgb[3];
			char s1[BUFSZ];
			char s2[BUFSZ];
			char s3[BUFSZ];
			
			if (fgets(buf, BUFSZ, ifp) == NULL)
				break;
	
			if ((rv = sscanf(buf, " %lf %lf %lf %s %s %s\n",&rgb[0], &rgb[1], &rgb[2], s1, s2, s3)) >= 4) {
				wo->count++;
			}
		}
	
	   	wo->nDeviceCoords =	3;	/* Num of device coordinates */
		strcpy(wo->prefix,""); /* Prefix for each color name, max 32, null terminated */
		strcpy(wo->suffix,""); /* Suffix for each color name, max 32, null terminated */
	
		wo->allocate((icmBase *)wo);	/* Allocate named color structures */
	
		if (verb)
			printf("Counted %d colors\n",wo->count);

		/* Read lines and colors */
		rewind(ifp);
		for (i = 0; i < wo->count; i++) {
			double rgb[3], lab[3];
			char s1[BUFSZ];
			char s2[BUFSZ];
			char s3[BUFSZ];
			unsigned int j;
			
			if (fgets(buf, BUFSZ, ifp) == NULL)
				break;
	
			if ((rv = sscanf(buf, " %lf %lf %lf %s %s %s\n",&rgb[0], &rgb[1], &rgb[2], s1, s2, s3)) >= 4) {
				rgb[0] /= 255.0;
				rgb[1] /= 255.0;
				rgb[2] /= 255.0;
				if (rv >= 5) {
					strcat(s1, " ");
					strcat(s1, s2);
				}
				if (rv >= 6) {
					strcat(s1, " ");
					strcat(s1, s3);
				}
				/* Convert from sRGB to Bradford adapted D50 */
				icx_sRGB2XYZ(lab, icmD50_ary3, rgb);
				icmXYZ2Lab(&icmD50, lab, lab);
	
				if (verb)
					printf("Got %f %f %f '%s'\n",rgb[0], rgb[1], rgb[2], s1);
	
				strncpy(wo->data[i].root,s1,31);
				wo->data[i].root[31] = '\000';
	
				for (j = 0; j < wo->nDeviceCoords; j++)
					wo->data[i].deviceCoords[j] = rgb[j];
				for (j = 0; j < 3; j++)
					wo->data[i].pcsCoords[j] = lab[j];
			}
		}
		fclose(ifp);
	}

	if ((rv = wr_icco->write(wr_icco, wr_fp, 0)) != 0)
		error ("Write file: %d, %s",rv,wr_icco->err);

	wr_icco->del(wr_icco);
	wr_fp->del(wr_fp);

	return 0;
}
























