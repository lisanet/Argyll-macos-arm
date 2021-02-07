
/* 
 * Argyll Color Correction System
 * CxF to CIE/TI1/TI2/TI3 conversion utility.
 *
 * Author: Graeme W. Gill
 * Date:   22/1/2019
 *
 * Copyright 2019 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Based on namedc.c
 */

/*
 * The idea is that this will handle all sorts of CxF3 to ArgyllCMS 
 * conversions, but support will have to be added on a case-by-case
 * basis.
 * 
 * Currently supports CxF3 Input Chart references (Outputs .cie)
 * 
 * May go through the motions of generating an RGB or CMYK .ti3, but
 * some tags will be missing.
 *
 */

/*
 * TTBD:
 *
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
#include "conv.h"
#include "cgats.h"
#include "xspect.h"
#include "mxml.h"
#ifdef STANDALONE_TEST
#include "ui.h"
#endif /* STANDALONE_TEST */

#ifdef NT       /* You'd think there might be some standards.... */
# ifndef __BORLANDC__
#  define stricmp _stricmp
# endif
#else
# define stricmp strcasecmp
#endif

/* Context for parsing file */
typedef struct {
	int indata;				/* State flag for sax_cb() (Not used) */
    char pfx[100];          /* Prefix to apply */
#define NAMEDC_PLEN 500
	char prefix[NAMEDC_PLEN];   /* Temporary buffers to use */

} cxfctx;

/* Compare while ignoring possible prefix on s1 */
static int pfxcmp(const char *s1, const char *s2) {
	const char *cp;

	if ((cp = strchr(s1, ':')) == NULL)
		cp = s1;
	else
		cp++;

	return strcmp(cp, s2);
}

/* XML data type callback for mxmlLoadFile() */
static mxml_type_t
type_cb(mxml_node_t *node) {
	const char *name = node->value.element.name;
	mxml_node_t *parent = mxmlGetParent(node);
	const char *pname = NULL;

	if (parent != NULL)
		pname = parent->value.element.name;

//	printf("~1 type_cb got node named '%s', parent '%s'\n",name,pname);

	if (pname == NULL)
		return MXML_TEXT;

	// ReflectanceSpectrum
	// Example doesn't have NumPoints, Increment
	// <ReflectanceSpectrum MeasureDate="2003-09-28T12:15:33-05:00"_ColorSpecification="CSD65-2" Name="45/0 Spectral" StartWL="400" NumPoints="31" Increment="10">
	// 0.0580 0.0594 0.0594 0.0584 0.0581 0.0591 0.0599 0.0601 0.0603 0.0610 0.0634 0.0695 0.0760 0.0786 0.0798 0.0826 0.0897 0.1024 0.1197 0.1350 0.1434 0.1455 0.1499 0.1594 0.1721 0.1842 0.1913 0.1928 0.1878 0.1734 0.1704

	if (pfxcmp(pname, "ColorValues") == 0
	 && pfxcmp(name, "ReflectanceSpectrum") == 0)
		return MXML_REAL;
//		return MXML_OPAQUE;

	if ((pfxcmp(pname, "ColorCIELab") == 0
	  || pfxcmp(pname, "ColorSpaceCIELab") == 0)
	 && (pfxcmp(name, "L") == 0
	  || pfxcmp(name, "A") == 0
	  || pfxcmp(name, "B") == 0))
		return MXML_REAL;

	if ((pfxcmp(pname, "ColorCIEXYZ") == 0
	  || pfxcmp(pname, "ColorIEXYZ") == 0
	  || pfxcmp(pname, "ColorSpaceCIEXYZ") == 0)
	 && (pfxcmp(name, "X") == 0
	  || pfxcmp(name, "Y") == 0
	  || pfxcmp(name, "Z") == 0))
		return MXML_REAL;

	if ((pfxcmp(pname, "ColorSRGB") == 0
	  || pfxcmp(pname, "ColorSpaceSRGB") == 0)
	 && (pfxcmp(name, "R") == 0
	  || pfxcmp(name, "G") == 0
	  || pfxcmp(name, "B") == 0))
		return MXML_REAL;

	if ((pfxcmp(pname, "ColorCMYK") == 0
	  || pfxcmp(pname, "ColorSpaceCMYK") == 0)
	 && (pfxcmp(name, "Cyan") == 0
	  || pfxcmp(name, "Magenta") == 0
	  || pfxcmp(name, "Yellow") == 0
	  || pfxcmp(name, "Black") == 0))
		return MXML_REAL;

	if ((pfxcmp(pname, "FileInformation") == 0
	  || pfxcmp(pname, "Header") == 0)
	 && (pfxcmp(name, "Creator") == 0
	  || pfxcmp(name, "CreationDate") == 0
	  || pfxcmp(name, "Description") == 0))
		return MXML_OPAQUE;		/* Don't split strings up */

	if ((pfxcmp(pname, "IlluminationOptions") == 0
	  || pfxcmp(pname, "TristimulusSpec") == 0
	  || pfxcmp(pname, "ColorSpaceSpecificationSpectrumTristimulus") == 0)
	 && (pfxcmp(name, "Illuminant") == 0
	  || pfxcmp(name, "Observer") == 0
	  || pfxcmp(name, "FieldOfView") == 0))
		return MXML_OPAQUE;		/* Don't split strings up */

	return MXML_TEXT;
}

static void sax_cb(mxml_node_t *node, mxml_sax_event_t event, void *data) {
	cxfctx *p = (cxfctx *)data;

	if (event == MXML_SAX_ELEMENT_OPEN) {
		const char *pname = mxmlGetElement(node);

		if (pfxcmp(pname, "Resources") == 0) {
//printf("~1 open of Resources\n");
			p->indata = 1;
		}

		mxmlRetain(node);

	} else if (event == MXML_SAX_DIRECTIVE) {
		mxmlRetain(node);

	} else if (event == MXML_SAX_DATA) {
		/* If the parent was retained, then retain this data node as well. */
		if (mxmlGetRefCount(mxmlGetParent(node)) > 1) {
			mxmlRetain(node);
		}
	} else if (event == MXML_SAX_ELEMENT_CLOSE) {
		const char *pname = mxmlGetElement(node);
		if(pfxcmp(pname, "Resources") == 0) {
//printf("~1 close of Resources\n");
			p->indata = 0;
		}
	}
}


/* Return a temporary key with prefix */
static char *pfx(cxfctx *p, char *key) {
	if (p->pfx[0] == '\000')
		return key;

	snprintf(p->prefix, NAMEDC_PLEN, "%s:%s", p->pfx, key);

	return p->prefix;
}

/* Return a temporary key with prefix at front and after each '/' */
static char *pfxp(cxfctx *p, char *path) {
	char *cp = p->prefix;
	char *sp, *ep;
	int plen, slen;

	if (p->pfx[0] == '\000')
		return path;

	plen = strlen(p->pfx);

	/* Figure out the first source key */
	sp = path;
	if ((ep = strchr(sp, '/')) != NULL) {
		slen = ep - sp;
		ep++;				/* Point beyond '/' */
	} else {
		slen = strlen(sp);
		ep = sp + slen;		/* point at '\000' */
	}

	for (;;) {
		if ((NAMEDC_PLEN -1 -(cp - p->prefix)) < (plen + 1 + slen))
			break;			/* No room for prefix + path element */

		strcpy(cp, p->pfx);
		cp += plen;
		*cp++ = ':';
		strncpy(cp, sp, slen);
		cp += slen;
		
		/* Figure out the next source key */
		sp = ep;
		
		if ((ep = strchr(sp, '/')) != NULL) {
			slen = ep - sp;
			ep++;				/* Point beyond '/' */
		} else {
			slen = strlen(sp);
			ep = sp + slen;		/* point at '\000' */
		}

		if (slen <= 0)			/* No more path elements */
			break;

		if ((NAMEDC_PLEN -1 -(cp - p->prefix)) < 1)
			break;			/* No room for '/' */
		*cp++ = '/';
	}
	*cp++ = '\000';

	return p->prefix;
}

/* Return a temporary key with suffix */
static char *sfx(cxfctx *p, char *key) {
	if (p->pfx[0] == '\000')
		return key;

	snprintf(p->prefix, 200, "%s:%s", key, p->pfx);

	return p->prefix;
}

void usage(char *diag, ...) {
	fprintf(stderr,"cxf2ti3\n");
	fprintf(stderr,"Author: Graeme W. Gill\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: cxf2ti3 [-v level] infile.cxf outbase\n");
	fprintf(stderr," -v level               Verbosity level 1-9\n");
	fprintf(stderr," -D level               Debugging level 1-9\n");
	fprintf(stderr," infile.cxf             CxF3 file to convert\n");
	fprintf(stderr," outbase                Basename of output file\n");
	exit(1);
}

int
main(int argc, char *argv[]) {
	int fa,nfa,mfa;				/* argument we're looking at */
	cxfctx ctx;
	static char inname[MAXNAMEL+1];
	static char outname[MAXNAMEL+9] = { 0 };		/* Output cgats .ti3 file base name */
	int verb = 0;
	int debug = 0;

    FILE *fp;				/* Input file */
    mxml_node_t *tree, *cxf, *pnode, *node;
	const char *attr, *name;
	int cxf2 = 0;
	int i, j;
	int phase;

	int ti;					/* Temporary field index */
	cgats *ocg = NULL;			/* output cgats structure for .ti2 */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	char buf[200];

	int fix_id = 0;			/* Field id's */
	int fix_loc = 0;
	int fix_lab = 0;
	int fix_xyz = 0;
	int fix_spect = 0;
	int fix_rgb = 0;
	int fix_cmyk = 0;

	int id;
	int has_loc = 0;	
	int has_lab = 0;	
	int has_xyz = 0;	
	int has_spect = 0;	
	int has_rgb = 0;	
	int has_cmyk = 0;	
	int specmin = 0, specinc = 0, specmax = 0, specnum = 0;	/* Min and max spectral in nm, inclusive */
	int npat = 0;			/* Number of patches */
	int got_descriptor = 0;
	int got_orginator = 0;
	int got_created = 0;
	int nsetel = 0;

	error_program = "CxF2ti3";

	/* Process the arguments */
	mfa = 2;						/* Minium final arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')	{	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1+mfa) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage(NULL);

			/* Verbose level */
			else if (argv[fa][1] == 'v') {
				if (na != NULL) { 
					verb = atoi(na);
					fa = nfa;
				} else
					verb = 1;
			}

			/* Debug level */
			else if (argv[fa][1] == 'D') {
				if (na != NULL) { 
					debug = atoi(na);
					fa = nfa;
				} else
					debug = 1;
			}

			else 
				usage("Unknown option '%c'",argv[fa][1]);
		}
		else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage("Missing input filename");
	strncpy(inname,argv[fa++],MAXNAMEL); inname[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage("Missing output basename");
	strncpy(outname,argv[fa++],MAXNAMEL); outname[MAXNAMEL] = '\000';

	g_log->verb = verb;
	g_log->debug = debug;

	/* - - - - - - - - - - - - */

	a1logd(g_log, 1, "cxf2ti3: file '%s\n",inname);
	
    if ((fp = fopen(inname, "r")) == NULL) {
		error("Opening CxF file '%s' failed",inname);
	}

	memset((void *)&ctx, 0, sizeof(cxfctx));

#ifdef NEVER
    tree = mxmlLoadFile(NULL, fp, type_cb);
#else
	tree = mxmlSAXLoadFd(NULL, fileno(fp), type_cb, sax_cb, (void *)&ctx);
#endif
    fclose(fp);

	if (tree == NULL) {
		error("Parsing XML file '%s' failed",inname);
	}

	if ((cxf = mxmlFindElement(tree, tree, NULL, NULL, NULL, MXML_DESCEND_FIRST)) == NULL
	 || mxmlGetType(cxf) != MXML_ELEMENT) {
		if (cxf == NULL)
			error("Failed to find top element in '%s'",inname);
		else
			error("Top element is not type MXML_ELEMENT in '%s'",inname);
	}
	name = cxf->value.element.name;

	if ((attr = strchr(cxf->value.element.name, ':')) != NULL) {
		int len = attr - name;
		if (len > 99)
			len = 99;
		strncpy(ctx.pfx, name, len);
		ctx.pfx[len] = '\000';
	} else {
		ctx.pfx[0] = '\000';
	}

	a1logd(g_log, 4, "cxf2ti3: prefix '%s'\n",ctx.pfx);

	if (strcmp(name, pfx(&ctx, "CxF")) != 0) {
		error("Top element not called CxF in '%s'",inname);
	}

	/* Check that its CxF3 or CxF2 */
	if ((attr = mxmlElementGetAttr(cxf, sfx(&ctx, "xmlns"))) == NULL) {
		error("Failed to find CxF attribute %s in '%s'",sfx(&ctx, "xmlns"), inname);
	}

	if (strcmp(attr, "http://colorexchangeformat.com/CxF3-core") != 0) {
		if (strcmp(attr, "http://colorexchangeformat.com/v2") != 0) {
			error("cxf2ti3: %s\n","File '%s' is not CxF format",inname);
		} else {
			cxf2 = 1;
			a1logd(g_log, 1, "cxf2ti3: This is a CxF2 file\n");
		}
	} else {
		a1logd(g_log, 1, "cxf2ti3: This is a CxF3 file\n");
	}

	ocg = new_cgats();	/* Create a CGATS structure */

	ocg->add_table(ocg, tt_none, 0);	/* Start the first table, type to be determined */

	/* Field we always want */
	fix_id = ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
	nsetel++;

	/* - - - - - - - - - - - - */

	/* Copy the description */
	if ((node = mxmlFindPathNode(cxf, pfxp(&ctx,"FileInformation/Description"))) != NULL) {
		name = mxmlGetOpaque(node);
		a1logd(g_log, 2, "cxf2ti3: got description '%s'\n",name);
		ocg->add_kword(ocg, 0, "DESCRIPTOR", name, NULL);
		got_descriptor = 1;
	}
		
	/* Grab the creator */
	if ((node = mxmlFindPathNode(cxf, pfxp(&ctx,"FileInformation/Creator"))) != NULL) {
		name = mxmlGetOpaque(node);
		a1logd(g_log, 2, "cxf2ti3: got creator '%s'\n",name);
		ocg->add_kword(ocg, 0, "ORIGINATOR", name, NULL);
		got_orginator = 1;
	}

	if ((node = mxmlFindPathNode(cxf, pfxp(&ctx,"FileInformation/CreationDate"))) != NULL) {
		name = mxmlGetOpaque(node);
		a1logd(g_log, 2, "cxf2ti3: got creation date '%s'\n",name);
		ocg->add_kword(ocg, 0, "CREATED", name, NULL);
		got_created = 1;
	}

	/* Look through the color specifications and see if there are spectral details */
	pnode = mxmlFindPathNode(cxf, pfxp(&ctx,"Resources/ColorSpecificationCollection/ColorSpecification/WavelengthRange"));
	if (pnode == NULL)
		pnode = mxmlFindPathNode(cxf, pfxp(&ctx,"Resources/ColorSpecificationCollection/ColorSpecification/MeasurementSpec/WavelengthRange"));
	while (pnode != NULL) {
		name = mxmlElementGetAttr(pnode, "StartWL");
		if (name != NULL) {
			specmin = atoi(name);
			a1logd(g_log, 2, "cxf2ti3: got StartWL %d\n",specmin);
		}
		name = mxmlElementGetAttr(pnode, "Increment");
		if (name != NULL) {
			specinc = atoi(name);
			a1logd(g_log, 2, "cxf2ti3: got Increment %d\n",specinc);
		}
		pnode = mxmlGetNextSibling(pnode);
	}

	/* Setup, and then read all the color data */
	/* phase == 0 setup */
	/* phase == 1 read */
	for (phase = 0; phase < 2; phase++) {
		const char *LabSpecification = NULL;
		const char *XYZSpecification = NULL;
		const char *SpectSpecification = NULL;
		const char *RGBSpecification = NULL;
		const char *CMYKSpecification = NULL;
		cgats_set_elem *setel = NULL;			/* Array of set value elements */

		if (phase == 1) {
			if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * nsetel)) == NULL)
				error("Failed to malloc setel length %d",nsetel);
		}

		/* Locate the Resources/ObjectCollection node */
		if ((pnode = mxmlFindPathNode(cxf, pfxp(&ctx,"Resources/ObjectCollection"))) == NULL) {
			error("Failed to find Resources/ObjectCollection in '%s'",inname);
		}

		/* - - - - - - - - - - - - */
		/* Setup or Read all the colors. */
		/* Start with first node */
		node = mxmlFindElement(pnode, pnode, pfx(&ctx,"Object"), NULL, NULL, MXML_DESCEND_FIRST);
		id = 1;
		for (i = 0; node != NULL; i++) {
			int j;
		    mxml_node_t *ppvals, *pvals, *val;
			const char *name;
			double Lab[3];
			double XYZ[3];
			xspect spect;
			double RGB[3];
			double CMYK[4];
			int six = 0;			/* setel index */
			
			if (mxmlGetType(node) != MXML_ELEMENT) {
				a1logd(g_log, 6, "cxf2ti3: skipping non element node type %d\n",mxmlGetType(node));
				goto next;
			}
			a1logd(g_log, 6, "cxf2ti3: read node '%s'\n",node->value.element.name);

			if (strcmp(node->value.element.name, pfx(&ctx,"Object")) != 0) {
				a1logd(g_log, 6, "cxf2ti3: skipping non %s node\n","Object");
				goto next; 
			}
 
			/* Get the ObjectType attribute of Object */
			if ((attr = mxmlElementGetAttr(node, "ObjectType")) == NULL) {
				a1logd(g_log, 6, "cxf2ti3: skipping node without ObjectType\n");
				goto next;							/* Skip this one */
			}
	
#ifdef NEVER
			/* Check the attribute value (Should we ?) */
			/* Known values are: Standard, Color */
			/* May be Trial, Target, Substrate, Colorant, ... ? */
			if (strcmp(attr, "Standard") != 0
			 && strcmp(attr, "Color") != 0) {
				a1logd(g_log, 6, "cxf2ti3: skipping node with ObjectType = '%s'\n",attr);
				goto next;							/* Skip this one */
			}
#endif

			if ((ppvals = mxmlFindElement(node, node, pfx(&ctx,"ColorValues"), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
				a1logd(g_log, 4, "cxf2ti3: no reference ColorValues element - skipping\n");
				goto next;
			}

			if (phase == 1) {
				sprintf(buf,"%d", id);
				setel[six++].c = buf;  
			}

			if ((attr = mxmlElementGetAttr(node, "Name")) != NULL) {
				a1logd(g_log, 6, "cxf2ti3: got sample %d name %s\n",id,attr);
				if (!has_loc) {
					if (phase == 1)
						error("Name field found after first ColorValue in '%s'",inname);
					fix_loc = ocg->add_field(ocg, 0, "SAMPLE_LOC", cs_t);
					nsetel++;
					has_loc = 1;	
				} else {
					setel[six++].c = strdup(attr);  
				}
			}

			/* See if there is ColorCIELab */
			if ((pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorCIELab"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL
			 || (pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorSpaceCIELab"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
				char *key[3] = { "L", "A", "B" };

				a1logd(g_log, 4, "cxf2ti3: got ColorCIELab\n");

				/* Check which color specification is being used with Lab */
				if (LabSpecification == NULL) {
					LabSpecification = mxmlElementGetAttr(pvals, "ColorSpecification");
					a1logd(g_log, 4, "cxf2ti3: got LabSpecification '%s'\n",LabSpecification); 
				}

				for (j = 0; j < 3; j++) {
					if ((val = mxmlFindElement(pvals, pvals, pfx(&ctx, key[j]), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
						a1logd(g_log, 6, "cxf2ti3: failed to find ColorCIELab component %s\n",key[j]);
						break;		/* oops */
					}
					Lab[j] = mxmlGetReal(val);
					a1logd(g_log, 6, "cxf2ti3: got ColorCIELab component %s value %f\n",key[j],Lab[j]);
				}
				if (j >= 3) {
					if (!has_lab) {
						if (phase == 1)
							error("Lab field found after first ColorValue in '%s'",inname);
						fix_rgb = ocg->add_field(ocg, 0, "LAB_L", r_t);
						          ocg->add_field(ocg, 0, "LAB_A", r_t);
						          ocg->add_field(ocg, 0, "LAB_B", r_t);
						nsetel += 3;
						has_lab = 1;
					} else {
						setel[six++].d = Lab[0];  
						setel[six++].d = Lab[1];
						setel[six++].d = Lab[2];
					}
				}
			}

			/* See if there is ColorCIEXYZ */
			if ((pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorCIEXYZ"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL
			 || (pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorIEXYZ"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL
			 || (pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorSpaceCIEXYZ"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
				char *key[3] = { "X", "Y", "Z" };

				a1logd(g_log, 4, "cxf2ti3: got ColorCIEXYZ\n");

				/* Check which color specification is being used with XYZ */
				if (XYZSpecification == NULL) {
					XYZSpecification = mxmlElementGetAttr(pvals, "ColorSpecification");
					a1logd(g_log, 4, "cxf2ti3: got XYZSpecification '%s'\n",XYZSpecification); 
				}

				for (j = 0; j < 3; j++) {
					if ((val = mxmlFindElement(pvals, pvals, pfx(&ctx, key[j]), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
						a1logd(g_log, 6, "cxf2ti3: failed to find ColorCIEXYZ component %s\n",key[j]);
						break;		/* oops */
					}
					XYZ[j] = mxmlGetReal(val);
					a1logd(g_log, 6, "cxf2ti3: got ColorCIEXYZ component %s value %f\n",key[j],XYZ[j]);
				}
				if (j >= 3) {
					if (!has_xyz) {
						if (phase == 1)
							error("XYZ field found after first ColorValue in '%s'",inname);
						fix_xyz = ocg->add_field(ocg, 0, "XYZ_X", r_t);
						          ocg->add_field(ocg, 0, "XYZ_Y", r_t);
						          ocg->add_field(ocg, 0, "XYZ_Z", r_t);
						nsetel += 3;
						has_xyz = 1;
					} else {
						setel[six++].d = XYZ[0];  
						setel[six++].d = XYZ[1];
						setel[six++].d = XYZ[2];
					}
				}
			}

			/* Read any spectral values */
			if ((pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ReflectanceSpectrum"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
				const char *elem = NULL;

				a1logd(g_log, 4, "cxf2ti3: got ReflectanceSpectrum\n");

				/* Check which color specification is being used with XYZ */
				if (SpectSpecification == NULL) {
					SpectSpecification = mxmlElementGetAttr(pvals, "ColorSpecification");
					a1logd(g_log, 4, "cxf2ti3: got SpectSpecification '%s'\n",SpectSpecification); 
				}

				if ((elem = mxmlElementGetAttr(pvals, "StartWL")) != NULL) {
					int minwl = 0;
					a1logd(g_log, 4, "cxf2ti3: got StartWL '%s'\n",elem); 
					minwl = atoi(elem);

					if (specmin != 0) {
						if (specmin != minwl)
							error("Inconsistent StartWL %d vs. %d from '%s'",specmin,minwl,inname);
					} else {
						specmin = minwl;
					}
				}

				if ((elem = mxmlElementGetAttr(pvals, "Increment")) != NULL) {
					int inc = 0;
					a1logd(g_log, 4, "cxf2ti3: got Increment '%s'\n",elem); 
					inc = atoi(elem);

					if (specinc != 0) {
						if (specinc != inc)
							error("Inconsistent Increment %d vs. %d from '%s'",specinc,inc,inname);
					} else {
						specinc = inc;
					}
				}

				if (specmin == 0 || specinc == 0)
					error("specmin %d specinc %d from '%s'",specmin,specinc,inname);

				/* mxml stashes multiple elements as children.. */
				if ((val = mxmlGetFirstChild(pvals)) == NULL)
					error("No spectral values in '%s'\n",inname);

				for (j = 0; j < XSPECT_MAX_BANDS; j++) {
					spect.spec[j] = mxmlGetReal(val);
					a1logd(g_log, 6, "cxf2ti3: got spect component %d value %f\n",j,spect.spec[j]);
					val = mxmlGetNextSibling(val);
					if (val == NULL)
						break;
				}
				spect.spec_n = specnum = j+1;
				spect.spec_wl_short = (double)specmin;
				specmax  = specmin + (spect.spec_n-1) * specinc;
				spect.spec_wl_long = (double)specmax;
				spect.norm = 100.0;

				a1logd(g_log, 6, "cxf2ti3: wl num %d short %f long %f\n",spect.spec_n,spect.spec_wl_short,spect.spec_wl_long);

				if (!has_spect) {
					if (phase == 1)
						error("Spectral field found after first ColorValue in '%s'",inname);

					sprintf(buf,"%d", spect.spec_n);
					ocg->add_kword(ocg, 0, "SPECTRAL_BANDS",buf, NULL);
					sprintf(buf,"%f", spect.spec_wl_short);
					ocg->add_kword(ocg, 0, "SPECTRAL_START_NM",buf, NULL);
					sprintf(buf,"%f", spect.spec_wl_long);
					ocg->add_kword(ocg, 0, "SPECTRAL_END_NM",buf, NULL);

					/* Generate fields for spectral values */
					for (i = 0; i < spect.spec_n; i++) {
						int nm, fix;
				
						/* Compute nearest integer wavelength */
						nm = (int)(spect.spec_wl_short + ((double)i/(spect.spec_n-1.0))
						            * (spect.spec_wl_long - spect.spec_wl_short) + 0.5);
						
						sprintf(buf,"SPEC_%03d",nm);
						fix = ocg->add_field(ocg, 0, buf, r_t);
						nsetel++;
						if (i == 0)
							fix_spect = fix;
					}
					has_spect = 1;
				} else {
					for (j = 0; j < spect.spec_n; j++)
						setel[six++].d = spect.spec[j];  
				}
			}
		
			/* Read any device color values */
			if ((ppvals = mxmlFindElement(node, node, pfx(&ctx,"DeviceColorValues"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {

				a1logd(g_log, 4, "cxf2ti3: got DeviceColorValues\n");
			}

			if (ppvals != NULL) {

				// ~~99 could add read of other device type here

				/* See if there is ColorRGB */
				if ((pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorRGB"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL
				 || (pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorSpaceRGB"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
					char *key[3] = { "Red", "Green", "Blue" };

					a1logd(g_log, 4, "cxf2ti3: got ColorRGB\n");

					for (j = 0; j < 3; j++) {
						if ((val = mxmlFindElement(pvals, pvals, pfx(&ctx, key[j]), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
							a1logd(g_log, 6, "cxf2ti3: failed to find ColorRGB component %s\n",key[j]);
							break;		/* oops */
						}
						RGB[j] = mxmlGetReal(val);
						a1logd(g_log, 6, "cxf2ti3: got ColorRGB component %s value %f\n",key[j],RGB[j]);
					}
					if (j >= 3) {
						if (!has_rgb) {
							if (phase == 1)
								error("RGB field found after first ColorValue in '%s'",inname);
							fix_rgb = ocg->add_field(ocg, 0, "RGB_R", r_t);
							          ocg->add_field(ocg, 0, "RGB_G", r_t);
							          ocg->add_field(ocg, 0, "RGB_B", r_t);
							nsetel += 3;
							has_rgb = 1;
						} else {
							setel[six++].d = RGB[0];  
							setel[six++].d = RGB[1];
							setel[six++].d = RGB[2];
						}
					}
				}

				/* See if there is ColorCMYK */
				if ((pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorCMYK"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL
				 || (pvals = mxmlFindElement(ppvals, ppvals, pfx(&ctx,"ColorSpaceCMYK"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
					char *key[4] = { "Cyan", "Magenta", "Yellow", "Black" };

					a1logd(g_log, 4, "cxf2ti3: got ColorCMYK\n");

					for (j = 0; j < 4; j++) {
						if ((val = mxmlFindElement(pvals, pvals, pfx(&ctx, key[j]), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
							a1logd(g_log, 6, "cxf2ti3: failed to find ColorCMYK component %s\n",key[j]);
							break;		/* oops */
						}
						CMYK[j] = mxmlGetReal(val);
						a1logd(g_log, 6, "cxf2ti3: got ColorCMYK component %s value %f\n",key[j],CMYK[j]);
					}
					if (j >= 4) {
						if (!has_cmyk) {
							if (phase == 1)
								error("CMYK field found after first ColorValue in '%s'",inname);
							fix_cmyk = ocg->add_field(ocg, 0, "CMYK_C", r_t);
							           ocg->add_field(ocg, 0, "CMYK_M", r_t);
							           ocg->add_field(ocg, 0, "CMYK_Y", r_t);
							           ocg->add_field(ocg, 0, "CMYK_K", r_t);
							nsetel += 4;
							has_cmyk = 1;
						} else {
							setel[six++].d = CMYK[0];  
							setel[six++].d = CMYK[1];
							setel[six++].d = CMYK[2];
							setel[six++].d = CMYK[3];
						}
					}
				}
			}
			if (phase == 1)
				ocg->add_setarr(ocg, 0, setel);
		
			id++;

		next:;	/* Next color */
			node = mxmlGetNextSibling(node);

			if (phase == 0) {
				/* Decide what sort of file we're creating */
				if ((has_lab || has_xyz || has_spect)
				 && (has_rgb || has_cmyk)) {
					ocg->add_other(ocg, "CTI3"); 	/* our special type is Calibration Target Information 3 */
					ocg->set_table_type(ocg, 0, tt_other, 0);	/* Set it to TI3 */
					if (!got_descriptor)
						ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
					strcat(outname,".ti3");
			
				} else if (has_lab || has_xyz || has_spect) {
					ocg->add_other(ocg, "CGATS.17");
					ocg->set_table_type(ocg, 0, tt_other, 0);	/* Set it to CGATS */
					if (!got_descriptor)
						ocg->add_kword(ocg, 0, "DESCRIPTOR", "Input CIE reference file",NULL);
					strcat(outname,".cie");
				} else {
					error("Can't figure out what sort output file to use for '%s' is",inname);
				}

				break;
			}
		}
		if (setel != NULL)
			free(setel);
	}
	mxmlDelete(tree);

	/* - - - - - - - - - - - - */

	if (!got_orginator) {
		sprintf(buf, "ArgyllCMS translation from CxF file '%s'",inname);
		ocg->add_kword(ocg, 0, "ORIGINATOR", buf, NULL);
	}

	if (!got_created) {
		atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
		ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	}

	/* Guess what scale the spectral data is set to */
	if (has_spect) {
		double maxv = 0.0, val_scale;

		for (i = 0; i < ocg->t[0].nsets; i++) {
			for (j = 0; j < specnum; j++) {
				double vv;
				vv = *((double *)ocg->t[0].fdata[i][fix_spect + j]);
				if (vv > maxv)
					maxv = vv;
			}
		}

		if (maxv < 10.0) {
			val_scale = 100.0/1.0;
			if (verb) printf("Spectral max found = %f, scale by 100.0\n",maxv);
		} else if (maxv > 160.0) {
			val_scale = 100.0/255.0;
			if (verb) printf("Spectral max found = %f, scale by 100/255\n",maxv);
		} else {
			val_scale = 100.0/100.0;
			if (verb) printf("Spectral max found = %f, scale by 1.0\n",maxv);
		}

		for (i = 0; i < ocg->t[0].nsets; i++) {
			for (j = 0; j < specnum; j++) {
				double vv;
				vv = *((double *)ocg->t[0].fdata[i][fix_spect + j]);
				vv *= val_scale;
				*((double *)ocg->t[0].fdata[i][fix_spect + j]) = vv;
			}
		}
	}

	/* Guess what scale the RGB data is set to */
	if (has_rgb) {
		double maxv = 0.0, val_scale;

		for (i = 0; i < ocg->t[0].nsets; i++) {
			for (j = 0; j < 3; j++) {
				double vv;
				vv = *((double *)ocg->t[0].fdata[i][fix_rgb + j]);
				if (vv > maxv)
					maxv = vv;
			}
		}

		if (maxv < 10.0) {
			val_scale = 100.0/1.0;
			if (verb) printf("RGB max found = %f, scale by 100.0\n",maxv);
		} else if (maxv > 160.0) {
			val_scale = 100.0/255.0;
			if (verb) printf("RGB max found = %f, scale by 100/255\n",maxv);
		} else {
			val_scale = 100.0/100.0;
			if (verb) printf("RGB max found = %f, scale by 1.0\n",maxv);
		}

		for (i = 0; i < ocg->t[0].nsets; i++) {
			for (j = 0; j < 3; j++) {
				double vv;
				vv = *((double *)ocg->t[0].fdata[i][fix_rgb + j]);
				vv *= val_scale;
				*((double *)ocg->t[0].fdata[i][fix_rgb + j]) = vv;
			}
		}
	}


#ifdef NEVER
	if (disp) {
		ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);	/* What sort of device this is */

		if (disp == 2)
			ocg->add_kword(ocg, 0, "NORMALIZED_TO_Y_100","NO", NULL);	/* What sort of device this is */
	}else if (inp)
		ocg->add_kword(ocg, 0, "DEVICE_CLASS","INPUT", NULL);	/* What sort of device this is */
	else
		ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */

	/* Note what instrument the chart was read with */
	/* Assume this - could try reading from file INSTRUMENTATION "SpectroScan" ?? */
	ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", inst_name(instSpectrolino) , NULL);

	if (calstd != xcalstd_none)
		ocg->add_kword(ocg, 0, "DEVCALSTD", xcalstd2str(calstd), NULL);

	if (isPolarize)
		ocg->add_kword(ocg, 0, "INSTRUMENT_FILTER", "POLARIZED", NULL);
#endif // NEVER

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	a1logd(g_log, 1, "cxf2ti3: done - %d patches\n",ocg->t[0].nsets);

	ocg->del(ocg);

	return 0;
}
























