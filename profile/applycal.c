
/* 
 * ArgyllCMS.
 * Apply a device calibration to an ICC profile.
 *
 * Author:  Graeme W. Gill
 * Date:    2009/8/31
 * Version: 1.00
 *
 * Copyright 2009 Graeme W. Gill
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD:
 *
 *		Would be good if remove restores shared curves, rather than
 *		leaving duplicates.
 *
 * 		Could stash the whole cal file in a text tag.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <time.h>
#include "copyright.h"
#include "aconfig.h"
#include "cgats.h"
#include "numlib.h"
#include "rspl.h"
#include "xicc.h"

#undef DEBUG

#ifdef DEBUG
#undef DBG
#define DBG(xxx) printf xxx ;
#else
#undef DBG
#define DBG(xxx) 
#endif

void usage(char *diag, ...) {
	int i;
	fprintf(stderr,"Apply device calibration to an ICC profile, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: %s [-options] [calfile.cal] inprof%s [outprof%s]\n","applycal",ICC_FILE_EXT,ICC_FILE_EXT);
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -a              Apply or re-apply calibration (default)\n");
	fprintf(stderr," -u              Remove calibration\n");
	fprintf(stderr," -c              Check calibration\n");
	fprintf(stderr," calfile.cal     Calibration file to apply\n");
	fprintf(stderr," inprof%s      ICC profile to read\n",ICC_FILE_EXT);
	fprintf(stderr," outprof%s     modified ICC profile to write\n",ICC_FILE_EXT);
	exit(2);
}

/* A primary signature and its backup */
struct _tagsigpair {
	icTagSignature prim;
	icTagSignature back;
	int chan;
	int dir;			/* 0 = none or out, 1 = in */
}; typedef struct _tagsigpair tagsigpair;

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char cal_name[MAXNAMEL+1];
	char in_name[MAXNAMEL+1];
	char out_name[MAXNAMEL+1];
	xcal *cal = NULL;			/* Calibration to apply */
	icmFile *rd_fp = NULL, *wr_fp = NULL;
	icc *icco;
	int apply = 1;
	int remove = 0;
	int check = 0;
	int verb = 0;
	int found = -1;
	int rv;

	if (argc < 3)
		usage("Too few arguments");

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

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;

			else if (argv[fa][1] == 'a' || argv[fa][1] == 'A') {
				apply = 1;
				remove = 0;
				check = 0;
			}

			else if (argv[fa][1] == 'u' || argv[fa][1] == 'U') {
				apply = 0;
				remove = 1;
				check = 0;
			}

			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				apply = 0;
				remove = 0;
				check = 1;
			}

			else
				usage("Unknown option '%s'",argv[fa][1]);

		} else
			break;
	}

	if (apply) {
		if (fa >= argc || argv[fa][0] == '-') usage("Missing calibration filename");
		strncpy(cal_name,argv[fa++],MAXNAMEL); cal_name[MAXNAMEL] = '\000';
	}

	if (fa >= argc || argv[fa][0] == '-') usage("Missing input profile filename");
	strncpy(in_name,argv[fa++],MAXNAMEL); in_name[MAXNAMEL] = '\000';

	if (apply || remove) {
		if (fa >= argc || argv[fa][0] == '-') usage("Missing output profile name");
		strncpy(out_name,argv[fa++],MAXNAMEL); out_name[MAXNAMEL] = '\000';
	}
	if (fa < argc) usage("Extra argument '%s'",argv[fa]);

	DBG(("apply %d, remove %d, check %d, calfile '%s', infile '%s', outfile '%s'\n",apply,remove,check,cal_name,in_name,out_name));

	if (apply) {
		/* Open up the calibration file */
		if ((cal = new_xcal()) == NULL)
			error("new_xcal failed");
		if ((cal->read(cal, cal_name)) != 0)
			error("%s",cal->err);
	}

	/* Open up the profile for reading */
	if ((rd_fp = new_icmFileStd_name(in_name,"r")) == NULL)
		error ("Can't open file '%s'",in_name);

	if ((icco = new_icc()) == NULL)
		error ("Creation of ICC object failed");

	/* Read header etc. */
	if ((rv = icco->read(icco,rd_fp,0)) != 0)
		error ("%d, %s",rv,icco->err);

	/* Read every tag */
	if (icco->read_all_tags(icco) != 0) {
		error("Unable to read all tags: %d, %s",icco->errc,icco->err);
	}

	rd_fp->del(rd_fp);

	/* ======================================= */
	{
		tagsigpair sigs[] = {			/* Signatures to modify and their backups */
			{ icSigProfileDescriptionTag, 	icmMakeTag('A','R','0','T'), 0, 0 },
			{ icSigAToB0Tag,				icmMakeTag('A','R','1','0'), 0, 1 },
			{ icSigAToB1Tag,				icmMakeTag('A','R','2','0'), 0, 1 },
			{ icSigAToB2Tag,				icmMakeTag('A','R','3','0'), 0, 1 },
			{ icSigBToA0Tag,				icmMakeTag('A','R','4','0'), 0, 0 },
			{ icSigBToA1Tag,				icmMakeTag('A','R','5','0'), 0, 0 },
			{ icSigBToA2Tag,				icmMakeTag('A','R','6','0'), 0, 0 },
			{ icSigRedTRCTag,				icmMakeTag('A','R','7','0'), 0, 0 },
			{ icSigGreenTRCTag,				icmMakeTag('A','R','7','1'), 1, 0 },
			{ icSigBlueTRCTag,				icmMakeTag('A','R','7','2'), 2, 0 },
			{ icSigGrayTRCTag, 				icmMakeTag('A','R','8','0'), 0, 0 },
			{ 0, 0, 0 }
		};
		tagsigpair linksigs[] = {		/* Signatures to modify */
			{ icSigProfileDescriptionTag, 	icmMakeTag('A','R','0','T'), 0 },
			{ icSigAToB0Tag,				icmMakeTag('A','R','9','1'), 1 },
			{ icSigAToB0Tag,				icmMakeTag('A','R','9','0'), 0 },
			{ 0, 0, 0 }
		};
		tagsigpair *ssigp, *sigp;
		int ntags = 0;						/* Number of tags done */
		icmBase *dtags[50];					/* Address of tags done */
		int inp = 0;						/* NZ for input calibration */
		unsigned int i, j;

		if (icco->header->deviceClass == icSigInputClass && cal->devclass != icSigInputClass) {
			warning("Non-input calibration being applied to an input profile");
		}

		if (check) {
			DBG(("Checking...\n"));

			for (sigp = linksigs; sigp->prim != 0; sigp++) {
				icmBase *primt;

				if (sigp->prim != icSigProfileDescriptionTag)
					continue;

				if ((primt = icco->read_tag(icco, sigp->prim)) == NULL)
					error("Can't find icSigProfileDescriptionTag in profile"); 

					/* See if we have a backup */
				if (icco->read_tag(icco, sigp->back) != NULL)
					found = 1;
				else
					found = 0;
			}
			DBG(("found = %d\n",found));
			if (found < 0)
				error("Internal, didn't look for ProfileDescriptionTag");

			if (verb) {
				if (found)
					printf("Profile has had calibration applied\n");
				else
					printf("Profile has NOT had calibration applied\n");
			}
		} else if (apply) {
			DBG(("Applying...\n"));

			if (icco->header->deviceClass == icSigInputClass
			 || icco->header->deviceClass == icSigDisplayClass
			 || icco->header->deviceClass == icSigOutputClass) {

				DBG(("Input, display or output profile\n"));

				/* Check colorspace is compatible */
				if (cal->colspace != icco->header->colorSpace)
					error("Calibration space %s doesn't match profile %s",
					             icm2str(icmColorSpaceSignature, cal->colspace),
					             icm2str(icmColorSpaceSignature, icco->header->colorSpace));

				ssigp = sigs;
				/* Note the cal direction */
				if (cal->devclass == icSigInputClass)
					inp = 1;
				else
					inp = 0;

			} else if (icco->header->deviceClass == icSigLinkClass) {
				DBG(("Device link profile\n"));

				/* Check colorspace is compatible */
				if (cal->colspace != icco->header->pcs)
					error("Calibration space %s doesn't match profile %s",
					             icm2str(icmColorSpaceSignature, cal->colspace),
					             icm2str(icmColorSpaceSignature, icco->header->pcs));
				ssigp = linksigs;
				/* Note the cal direction */
				if (cal->devclass == icSigInputClass)
					inp = 1;
				else
					inp = 0;
			} else {
				error("Can't apply calibration to profile of class %s",
				      icm2str(icmProfileClassSignature, icco->header->deviceClass));
			}
			DBG(("input calibration = %d\n",inp));

			/* First pass is to duplicate any linked TRCTags */
			for (sigp = ssigp; sigp->prim != 0; sigp++) { /* Process each tag */
				icmBase *primt;

				DBG(("looking for tag '%s'\n",icm2str(icmTagSignature, sigp->prim)));
				if ((primt = icco->read_tag(icco, sigp->prim)) == NULL) {
					if (sigp->prim == icSigProfileDescriptionTag)
						error("Can't find icSigProfileDescriptionTag in profile"); 
					continue;		/* Don't have this tag */
				}

				/* XXXXTRCTag */
			    if (primt->ttype == icSigCurveType) {
					icmCurve *wo, *ro = (icmCurve *)primt;

					/* See if we have a backup */
					if ((wo = (icmCurve *)icco->read_tag(icco, sigp->back)) == NULL) {

						/* If tag is shared, we need to separated it */
						if (ro->refcount > 1) {
							DBG(("tag is shared, so separate it\n"));
							if ((wo = (icmCurve *)icco->add_tag(
							           icco, sigp->back, icSigCurveType)) == NULL) 
								error("Failed to create tag '%s'\n", icm2str(icmTagSignature, sigp->back));
							wo->flag = ro->flag;
							wo->size = ro->size;
							wo->allocate((icmBase *)wo);	/* Allocate space */
							for (i = 0; i < wo->size; i++)	/* Copy the curve */
								wo->data[i] = ro->data[i];

							if (icco->delete_tag(icco, sigp->prim))
								error("Failed to delete tag '%s'",icm2str(icmTagSignature, sigp->prim));
							if (icco->rename_tag(icco, sigp->back, sigp->prim))
								error("Failed to rename tag '%s' to '%s'",
								          icm2str(icmTagSignature, sigp->back),
								          icm2str(icmTagSignature, sigp->prim));
						}
					} else {
						if (ro->refcount > 1)
							error("Found tag %s has backup, but is shared",icm2str(icmTagSignature,sigp->prim));
					}
				}
			}
			/* Second pass we create backups and calibrated curves */
			for (sigp = ssigp; sigp->prim != 0; sigp++) { /* Process each tag */
				icmBase *primt;

				DBG(("looking for tag '%s'\n",icm2str(icmTagSignature, sigp->prim)));
				if ((primt = icco->read_tag(icco, sigp->prim)) == NULL) {
					if (sigp->prim == icSigProfileDescriptionTag)
						error("Can't find icSigProfileDescriptionTag in profile"); 
					continue;		/* Don't have this tag */
				}

				/* icSigProfileDescriptionTag type */
				if (primt->ttype == icSigTextDescriptionType) {
					icmTextDescription *wo, *ro = (icmTextDescription *)primt;
					char *extra = NULL;

					/* See if we've done this tag before due to links */
					for (i = 0; i < ntags; i++) {
						if (dtags[i] == primt)
							break;				/* Yes */
					}
					if (i < ntags) {
						DBG(("Found this tag before (link)\n"));
						continue;				/* Skip tag */
					}
					if (ntags >= 50)			/* Impossible */
						error("Internal, run out of previuos tags space");
					dtags[ntags++] = primt;		/* Remember this one */

					DBG(("ProfileDescriptionTag\n"));
					DBG(("Looking for backup tag '%s'\n",icm2str(icmTagSignature, sigp->back)));

					/* See if we have a backup */
					if ((wo = (icmTextDescription *)icco->read_tag(icco, sigp->back)) == NULL) {
						DBG(("No backup, creating one\n"));
						/* No, so create one */
						if ((wo = (icmTextDescription *)icco->add_tag(
						           icco, sigp->back, icSigTextDescriptionType)) == NULL) 
							error("Failed to create tag '%s'\n", icm2str(icmTagSignature, sigp->back));

						wo->size = ro->size;
						wo->allocate((icmBase *)wo);	/* Allocate space */
						strcpy(wo->desc, ro->desc);		/* Copy the string in */
						/* Hmm. what should we do with Unicode and script ? */
					}

					if (cal->xpi.profDesc != NULL)
						extra = cal->xpi.profDesc;
					else
						extra = cal_name;

					ro->size = strlen(ro->desc) + 3 + strlen(extra) + 3;
					ro->allocate((icmBase *)ro);	/* Allocate space */
					strcpy(ro->desc, wo->desc);
					strcat(ro->desc, " [ ");
					strcat(ro->desc, extra);
					strcat(ro->desc, " ]");

					DBG(("Set tag contents to '%s'\n",ro->desc));

				/* icSigAToBXTag or icSigBToAXTag */
				} else if (primt->ttype == icSigLut8Type
				        || primt->ttype == icSigLut16Type) {
					icmLut *ro = (icmLut *)primt;	/* Modified Lut */

					/* See if we've done this tag before due to links */
					for (i = 0; i < ntags; i++) {
						if (dtags[i] == primt)
							break;				/* Yes */
					}
					if (i < ntags) {
						DBG(("Found this tag before (link)\n"));
						continue;				/* Skip tag */
					}
					if (ntags >= 50)			/* Impossible */
						error("Internal, run out of previous tags space");
					dtags[ntags++] = primt;		/* Remember this one */

					DBG(("Lut8 or Lut16\n"));

					if (sigp->dir) {
						/* Apply calibration to the input table */
						for (j = 0; j < ro->inputChan; j++) {
							icTagSignature bsig = sigp->back;
							icmCurve *wo;					/* Backup of original */

							/* Create a tag per channel */
							bsig += j; 

							DBG(("Looking for backup tag '%s'\n",icm2str(icmTagSignature, bsig)));

							/* See if we have a backup */
							if ((wo = (icmCurve *)icco->read_tag(icco, bsig)) == NULL) {
								DBG(("No backup, creating one\n"));
								/* No, so create one */
								if ((wo = (icmCurve *)icco->add_tag(
								           icco, bsig, icSigCurveType)) == NULL) 
									error("Failed to create tag '%s'\n", icm2str(icmTagSignature, bsig));
								wo->flag = icmCurveSpec; 		/* Specified version */
								wo->size = ro->inputEnt;
								wo->allocate((icmBase *)wo);	/* Allocate space */
								for (i = 0; i < wo->size; i++)	/* Copy the curve */
									wo->data[i] = ro->inputTable[j * ro->inputEnt + i];
							}

							/* Create new input curve from inv cal + orginal curve */
							for (i = 0; i < ro->inputEnt; i++) {
								double val;
								val = i/(ro->inputEnt-1.0);
								if (inp)
									val = cal->interp_ch(cal, j, val);		/* Do calibration */
								else
									val = cal->inv_interp_ch(cal, j, val);	/* Undo calibration */
								wo->lookup_fwd(wo, &val, &val);		/* Original curve */
								ro->inputTable[j * ro->inputEnt + i] = val;
							}
							DBG(("Created calibrated input curve\n"));
						}
					} else {
						/* Apply calibration to the output table */
						for (j = 0; j < ro->outputChan; j++) {
							icTagSignature bsig = sigp->back;
							icmCurve *wo;					/* Backup of original */

							/* Create a tag per channel */
							bsig += j; 

							DBG(("Looking for backup tag '%s'\n",icm2str(icmTagSignature, bsig)));
							/* See if we have a backup */
							if ((wo = (icmCurve *)icco->read_tag(icco, bsig)) == NULL) {
								DBG(("No backup, creating one\n"));
								/* No, so create one */
								if ((wo = (icmCurve *)icco->add_tag(
								           icco, bsig, icSigCurveType)) == NULL) 
									error("Failed to create tag '%s'\n", icm2str(icmTagSignature, bsig));
						
								wo->flag = icmCurveSpec; 		/* Specified version */
								wo->size = ro->outputEnt;
								wo->allocate((icmBase *)wo);	/* Allocate space */
								for (i = 0; i < wo->size; i++)	/* Copy the curve */
									wo->data[i] = ro->outputTable[j * ro->outputEnt + i];
							}

							/* Create new output curve from original + cal */
							for (i = 0; i < ro->outputEnt; i++) {
								double val;
								val = i/(ro->outputEnt-1.0);
								wo->lookup_fwd(wo, &val, &val);			/* Original curve */
								if (inp)
									val = cal->interp_ch(cal, j, val);		/* Undo calibration */
								else
									val = cal->interp_ch(cal, j, val);		/* Do calibration */
								ro->outputTable[j * ro->outputEnt + i] = val;
							}
							DBG(("Created calibrated output curve\n"));
						}
					}

				/* XXXXTRCTag */
				} else if (primt->ttype == icSigCurveType) {
					icmCurve *wo, *ro = (icmCurve *)primt;

					DBG(("CurveType\n"));

					DBG(("Looking for backup tag '%s'\n",icm2str(icmTagSignature, sigp->back)));

					/* See if we have a backup */
					if ((wo = (icmCurve *)icco->read_tag(icco, sigp->back)) == NULL) {

						DBG(("No backup, creating one\n"));
						/* No, so create one */
						if ((wo = (icmCurve *)icco->add_tag(
						           icco, sigp->back, icSigCurveType)) == NULL) 
							error("Failed to create tag '%s'\n", icm2str(icmTagSignature, sigp->back));
				
						wo->flag = ro->flag;
						wo->size = ro->size;
						wo->allocate((icmBase *)wo);	/* Allocate space */
						for (i = 0; i < wo->size; i++)	/* Copy the curve */
							wo->data[i] = ro->data[i];

						/* Change type & size of ro if necessary */
						if (ro->flag != icmCurveSpec || wo->size < 256) {
							ro->flag = icmCurveSpec;
							ro->size = 256;
							ro->allocate((icmBase *)wo);	/* Allocate space */
						} 
					}

					/* Create new forward direction curve from cal + orginal curve */
					j = sigp->chan;
					for (i = 0; i < ro->size; i++) {
						double val;
						val = i/(ro->size-1.0);
//printf("~1 Input val %f", val);
						val = cal->inv_interp_ch(cal, j, val);	/* Inverse output calibration */
//printf(", after inv curve %f", val);
						wo->lookup_fwd(wo, &val, &val);			/* Original curve */
//printf(", after orig %f\n", val);
						ro->data[i] = val;
					}
					DBG(("Created calibrated %s curve for chan %d\n",inp ? "input" : "output",j));
				} else {
					error("Tag %s is type %s we don't know how to handle",
					      icm2str(icmTagSignature, sigp->prim),
					      icm2str(icmTypeSignature, primt->ttype));
				}
			}

		} else if (remove) {
			int k;
			DBG(("Removing...\n"));
			for (k = 0; k < 2; k++) {
				if (k == 0)
					ssigp = sigs;
				else if (k == 1)
					ssigp = linksigs;

				for (sigp = ssigp; sigp->prim != 0; sigp++) { /* Process each tag */
					icmBase *backt, *primt;

					DBG(("Looking for baclup tag '%s'\n",icm2str(icmTagSignature, sigp->back)));
					if ((backt = icco->read_tag(icco, sigp->back)) == NULL)
						continue;		/* Don't have this backup tag */
					
					DBG(("Looking for primary tag '%s'\n",icm2str(icmTagSignature, sigp->prim)));
					if ((primt = icco->read_tag(icco, sigp->prim)) == NULL) {
						error("Can't find primary tag %s for backup %s",
						                  icm2str(icmTagSignature, sigp->prim),
						                  icm2str(icmTagSignature, sigp->back));
					}

					/* icSigProfileDescriptionTag type */
					if (primt->ttype == icSigTextDescriptionType) {
						icmTextDescription *wo, *ro = (icmTextDescription *)primt;

						DBG(("ProfileDescriptionTag\n"));

						wo = (icmTextDescription *)backt;

						/* Restore primary table */
						ro->size = wo->size;
						ro->allocate((icmBase *)ro);	/* Reallocate space */
						strcpy(ro->desc, wo->desc);		/* Restore description */

						/* delete backup */
						if (icco->delete_tag(icco, sigp->back))
							error("Failed to delete tag '%s'",icm2str(icmTagSignature, sigp->prim));
						DBG(("Restored primary and deleted backup\n"));

					/* icSigAToBXTag or icSigBToAXTag */
					} else if (primt->ttype == icSigLut8Type
					        || primt->ttype == icSigLut16Type) {
						icmLut *ro = (icmLut *)primt;	/* Modified Lut */

						if (sigp->dir) {
							/* Restore the input table */
							for (j = 0; j < ro->inputChan; j++) {
								icTagSignature bsig = sigp->back;
								icmCurve *wo;					/* Backup of original */

								/* Create a tag per channel */
								bsig += j; 

								DBG(("Looking for backup tag '%s'\n",icm2str(icmTagSignature, bsig)));

								/* See if we have a backup */
								if ((wo = (icmCurve *)icco->read_tag(icco, bsig)) == NULL)
									error("Can't find original table data in tag %s",
									                 icm2str(icmTagSignature, sigp->back));

								/* Restore primary table */
								for (i = 0; i < wo->size; i++)	/* Copy the curve */
									ro->inputTable[j * ro->inputEnt + i] = wo->data[i];
			
								/* delete backup */
								if (icco->delete_tag(icco, bsig))
									error("Failed to delete tag '%s'",icm2str(icmTagSignature, bsig));
								DBG(("Restored primary and deleted backup\n"));
							}
						} else {
							/* Restore the output table */
							for (j = 0; j < ro->outputChan; j++) {
								icTagSignature bsig = sigp->back;
								icmCurve *wo;					/* Backup of original */

								/* Create a tag per channel */
								bsig += j; 

								DBG(("Looking for backup tag '%s'\n",icm2str(icmTagSignature, bsig)));

								/* See if we have a backup */
								if ((wo = (icmCurve *)icco->read_tag(icco, bsig)) == NULL)
									error("Can't find original table data in tag %s",
									                 icm2str(icmTagSignature, sigp->back));

								/* Restore primary table */
								for (i = 0; i < wo->size; i++)	/* Copy the curve */
									ro->outputTable[j * ro->outputEnt + i] = wo->data[i];
			
								/* delete backup */
								if (icco->delete_tag(icco, bsig))
									error("Failed to delete tag '%s'",icm2str(icmTagSignature, bsig));
								DBG(("Restored primary and deleted backup\n"));
							}
						}

					/* XXXXTRCTag */
					} else if (primt->ttype == icSigCurveType) {
						icmCurve *wo, *ro = (icmCurve *)primt;

						DBG(("CurveType\n"));
						DBG(("Looking for backup tag '%s'\n",icm2str(icmTagSignature, sigp->back)));

						/* See if we have a backup */
						wo = (icmCurve *)backt;

						/* Restore primary table */
						ro->flag = wo->flag;
						ro->size = wo->size;
						ro->allocate((icmBase *)ro);	/* Allocate space */
						for (i = 0; i < wo->size; i++)	/* Copy the curve */
							ro->data[i] = wo->data[i];

						/* delete backup */
						if (icco->delete_tag(icco, sigp->back))
							error("Failed to delete tag '%s'",icm2str(icmTagSignature, sigp->back));
						DBG(("Restored primary and deleted backup\n"));

					} else {
						error("Tag %s is type %s we don't know how to handle",
						      icm2str(icmTagSignature, sigp->prim),
						      icm2str(icmTypeSignature, primt->ttype));
					}
				}
			}
		}
	}
	/* ======================================= */
	
	if (apply || remove) {
		/* Open up the other profile for writing */
		if ((wr_fp = new_icmFileStd_name(out_name,"w")) == NULL)
			error ("Can't open file '%s'",out_name);
	
		if ((rv = icco->write(icco,wr_fp,0)) != 0)
			error ("Write file: %d, %s",rv,icco->err);
		wr_fp->del(wr_fp);
	}
	
	if (cal != NULL)
		cal->del(cal);
	icco->del(icco);

	if (found == 1)
		return 1;
	return 0;
}

