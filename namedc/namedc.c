
/* 
 * Argyll Color Correction System
 * Named Color Access Library
 *
 * Author: Graeme W. Gill
 * Date:   3/12/2013
 *
 * Copyright 2013 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * Currently supports CxF2, CxF3 & ICC named color profiles.
 */

/*
 * TTBD:
 *
 * Probably want to add some other formats:
 * (see <http://www.selapa.net/swatchbooker/>
 *  and <http://www.selapa.net/swatches/colors/fileformats.php>  )
 *
 * 	Adobe .aco
 * 	      .acb
 * 	      .acbl
 * 	      .ase
 * 	      .bcf
 * corel  .cpl
 * corel  .xml
 * autocad .acb			RGB only 
 *
 * and lower quality RGB formats:
 *
 * Adobe  .acf		?? May support Lab
 * corel  .acb
 *
 * Paint Shop Pro .pal ???
 *
 * RAL ?
 */

#undef DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <errno.h>
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
#include "namedc.h"
#ifdef STANDALONE_TEST
#include "ui.h"
#endif /* STANDALONE_TEST */

#ifndef STANDALONE_TEST

# define DEB6 6
# define DEB4 4

#ifdef NT       /* You'd think there might be some standards.... */
# ifndef __BORLANDC__
#  define stricmp _stricmp
# endif
#else
# define stricmp strcasecmp
#endif

/* Forward declarations */
static void clear_namedc(namedc *p);
static void clear_nce(nce *p);

static ORD32 do_hash(ORD32 hash, const char *s) {
	int i;
    for (i = 0; s[i] != '\000'; i++) {
        hash += s[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    return hash;
}

static int do_hash2(const char *s1, const char *s2) {
	ORD32 hash;

	hash = do_hash(0, s1);
	hash = do_hash(hash, s2);

	return (int)hash;
}

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

	if (pfxcmp(pname, "ColorValues") == 0
	 && pfxcmp(name, "ReflectanceSpectrum") == 0)
		return MXML_REAL;

	if ((pfxcmp(pname, "ColorCIELab") == 0
	  || pfxcmp(pname, "ColorSpaceCIELab") == 0)
	 && (pfxcmp(name, "L") == 0
	  || pfxcmp(name, "A") == 0
	  || pfxcmp(name, "B") == 0))
		return MXML_REAL;

	if ((pfxcmp(pname, "ColorCIEXYZ") == 0
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
	namedc *p = (namedc *)data;

	if (event == MXML_SAX_ELEMENT_OPEN) {
		const char *pname = mxmlGetElement(node);

		if (pfxcmp(pname, "Resources") == 0) {
//printf("~1 open of Resources\n");
			p->indata = 1;
		}

		if ((p->options & NAMEDC_OP_NODATA) == 0
		 || !p->indata) {
//printf("~1 options 0x%x, indata %d, retaining '%s'\n",p->options,p->indata,pname);
a1logd(p->log, 4, "sax_cb: retaining %s\n",pname);
			mxmlRetain(node);
		}

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
static char *pfx(namedc *p, char *key) {
	if (p->pfx[0] == '\000')
		return key;

	snprintf(p->prefix, NAMEDC_PLEN, "%s:%s", p->pfx, key);

	return p->prefix;
}

/* Return a temporary key with prefix at front and after each '/' */
static char *pfxp(namedc *p, char *path) {
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
static char *sfx(namedc *p, char *key) {
	if (p->pfx[0] == '\000')
		return key;

	snprintf(p->prefix, 200, "%s:%s", key, p->pfx);

	return p->prefix;
}

/* Read in a namedc from a cxf2 or cxf3 file */
/* Return nz on error */
static int read_cxf(namedc *p, const char *filename, int options) {
    FILE *fp;
	char *pfilename;
    mxml_node_t *tree, *cxf, *pnode, *node;
	const char *attr, *name;
	int cxf2 = 0;
	int specmin = 0, specinc = 0;
	int i, j;

	a1logd(p->log, 1, "read_cxf: file '%s' options 0x%x\n",filename,options);
	
	if (p->filename == NULL) {
		if ((pfilename = strdup(filename)) == NULL) {
			snprintf(p->err, NAMEDC_ERRL, "Malloc of filename failed");
			a1logd(p->log, 1, "read_cxf: %s\n",p->err);
			return p->errc = 2;
		}
		clear_namedc(p);
		p->filename = pfilename;
	}
	p->options = options;

    if ((fp = fopen(p->filename, "r")) == NULL) {
		snprintf(p->err, NAMEDC_ERRL, "Opening XML file '%s' failed with %s",
		                                     p->filename,strerror(errno));
		a1logd(p->log, 1, "read_cxf: %s\n",p->err);
		return p->errc = 1;
	}

#ifdef NEVER
    tree = mxmlLoadFile(NULL, fp, type_cb);
#else
	p->indata = 0;
	tree = mxmlSAXLoadFd(NULL, fileno(fp), type_cb, sax_cb, (void *)p);
#endif
    fclose(fp);

	if (tree == NULL) {
		snprintf(p->err, NAMEDC_ERRL, "Parsing XML file '%s' failed",p->filename);
		a1logd(p->log, 1, "read_cxf: %s\n",p->err);
		return p->errc = 1;
	}

	if ((cxf = mxmlFindElement(tree, tree, NULL, NULL, NULL, MXML_DESCEND_FIRST)) == NULL
	 || mxmlGetType(cxf) != MXML_ELEMENT) {
		snprintf(p->err, NAMEDC_ERRL, "Failed to find top element in '%s'",p->filename);
		a1logd(p->log, 1, "read_cxf: %s\n",p->err);
		mxmlDelete(tree);
		return p->errc = 1;
	}
	name = cxf->value.element.name;

	if ((attr = strchr(cxf->value.element.name, ':')) != NULL) {
		int len = attr - name;
		if (len > 99)
			len = 99;
		strncpy(p->pfx, name, len);
		p->pfx[len] = '\000';
	} else {
		p->pfx[0] = '\000';
	}

	a1logd(p->log, 4, "read_cxf: prefix '%s'\n",p->pfx);

	if (strcmp(name, pfx(p, "CxF")) != 0) {
		snprintf(p->err, NAMEDC_ERRL, "Top element not called CxF in '%s'",p->filename);
		a1logd(p->log, 1, "read_cxf: %s\n",p->err);
		mxmlDelete(tree);
		return p->errc = 1;
	}

	/* Check that its CxF3 */
	if ((attr = mxmlElementGetAttr(cxf, sfx(p, "xmlns"))) == NULL) {
		snprintf(p->err, NAMEDC_ERRL, "Failed to find CxF attribute %s in '%s'",sfx(p, "xmlns"), p->filename);
		a1logd(p->log, 1, "read_cxf: %s\n",p->err);
		mxmlDelete(tree);
		return p->errc = 1;
	}

	if (strcmp(attr, "http://colorexchangeformat.com/CxF3-core") != 0) {
		if (strcmp(attr, "http://colorexchangeformat.com/v2") != 0) {
			snprintf(p->err, NAMEDC_ERRL, "File '%s' is not CxF format",p->filename);
			a1logd(p->log, 1, "read_cxf: %s\n",p->err);
			mxmlDelete(tree);
			return p->errc = 1;
		} else {
			cxf2 = 1;
			a1logd(p->log, 1, "read_cxf: This is a CxF2 file\n");
		}
	} else {
		a1logd(p->log, 1, "read_cxf: This is a CxF3 file\n");
	}

	/* Grab the description */
	if (cxf2) {
		if ((node = mxmlFindPathNode(cxf, pfxp(p,"Palette"))) == NULL
		 || (name = mxmlElementGetAttr(node, "PaletteName")) == NULL)
			name = NULL;

		if (name == NULL)
			name = p->filename;
		if ((p->description = strdup(name)) == NULL) {
			snprintf(p->err, NAMEDC_ERRL, "Malloc of description string failed");
			a1logd(p->log, 1, "read_cxf: %s\n",p->err);
			mxmlDelete(tree);
			return p->errc = 2;
		}
		a1logd(p->log, 2, "read_cxf: description '%s'\n",p->description);
	
		p->hash = do_hash2(p->filename, p->description);
	
	} else {	/* else cxf3 */
		if ((node = mxmlFindPathNode(cxf, pfxp(p,"FileInformation/Description"))) == NULL)
			name = NULL;
		else
			name = mxmlGetOpaque(node);
		if (name == NULL)
			name = p->filename;
		if ((p->description = strdup(name)) == NULL) {
			snprintf(p->err, NAMEDC_ERRL, "Malloc of description string failed");
			a1logd(p->log, 1, "read_cxf: %s\n",p->err);
			mxmlDelete(tree);
			return p->errc = 2;
		}
		a1logd(p->log, 2, "read_cxf: description '%s'\n",p->description);
	
		p->hash = do_hash2(p->filename, p->description);
	}

	/* Look through the color specifications and see if there are spectral details */
	pnode = mxmlFindPathNode(cxf, pfxp(p,"Resources/ColorSpecificationCollection/ColorSpecification/WavelengthRange"));
	if (pnode == NULL)
		pnode = mxmlFindPathNode(cxf, pfxp(p,"Resources/ColorSpecificationCollection/ColorSpecification/MeasurementSpec/WavelengthRange"));
	while (pnode != NULL) {
		name = mxmlElementGetAttr(pnode, "StartWL");
		if (name != NULL) {
			specmin = atoi(name);
			a1logd(p->log, 2, "read_cxf: got StartWL %d\n",specmin);
		}
		name = mxmlElementGetAttr(pnode, "Increment");
		if (name != NULL) {
			specinc = atoi(name);
			a1logd(p->log, 2, "read_cxf: got Increment %d\n",specinc);
		}
		pnode = mxmlGetNextSibling(pnode);
	}

	if ((p->options & NAMEDC_OP_NODATA) == 0) {
		char *SampleKey = NULL;
		char *SampleNameKey = NULL;
		char *SampleLabKey = NULL;
		const char *SpectSpecification = NULL;
		char *SampleXYZKey = NULL;
		char *SampleCMYKKey = NULL;
		const char *colorSpecification = NULL;	/* cxf3 ColorSpecification of Lab or XYZ */
		xsp2cie *sp2cie = NULL;

		if (cxf2) {
			/* Locate the Palette/ColorSet node */
			if ((pnode = mxmlFindPathNode(cxf, pfxp(p,"Palette/ColorSet"))) == NULL) {
				snprintf(p->err, NAMEDC_ERRL, "Failed to find Resources/ObjectCollection in '%s'",p->filename);
				a1logd(p->log, 1, "read_cxf: %s\n",p->err);
				mxmlDelete(tree);
				return p->errc = 1;
			}
			SampleKey = "Color";
			SampleNameKey = "ColorName";
			SampleLabKey = "ColorSpaceCIELab";
			SampleXYZKey = "ColorSpaceCIEXYZ";
			SampleCMYKKey = "ColorSpaceCMYK";

		} else {
			/* Locate the Resources/ObjectCollection node */
			if ((pnode = mxmlFindPathNode(cxf, pfxp(p,"Resources/ObjectCollection"))) == NULL) {
				snprintf(p->err, NAMEDC_ERRL, "Failed to find Resources/ObjectCollection in '%s'",p->filename);
				a1logd(p->log, 1, "read_cxf: %s\n",p->err);
				mxmlDelete(tree);
				return p->errc = 1;
			}
			SampleKey = "Object"; 
			SampleNameKey = "Name";
			SampleLabKey = "ColorCIELab";
			SampleXYZKey = "ColorCIEXYZ";
			SampleCMYKKey = "ColorCMYK";
		}

		/* - - - - - - - - - - - - */
		/* Read all the colors. */
		/* Start with first node */
		node = mxmlFindElement(pnode, pnode, pfx(p,SampleKey), NULL, NULL, MXML_DESCEND_FIRST);
		for (i = 0; node != NULL; i++) {
			int j;
		    mxml_node_t *ppvals, *pvals, *val;
			const char *name;
			double Lab[3];
			int Lab_v = 0;
			xspect *sp = NULL;
			double dev[MAX_CHAN];
			int dev_n = 0;
			icColorSpaceSignature devSig = icMaxEnumData;
			
			if (mxmlGetType(node) != MXML_ELEMENT) {
				a1logd(p->log, DEB6, "read_cxf: skipping non element node type %d\n",mxmlGetType(node));
				goto next;
			}
			a1logd(p->log, DEB6, "read_cxf: read node '%s'\n",node->value.element.name);

			if (strcmp(node->value.element.name, pfx(p,SampleKey)) != 0) {
				a1logd(p->log, DEB6, "read_cxf: skipping non %s node\n",SampleKey);
				goto next; 
			}
 
			if (!cxf2) {
				/* Get the ObjectType attribute of Object */
				if ((attr = mxmlElementGetAttr(node, "ObjectType")) == NULL) {
					a1logd(p->log, DEB6, "read_cxf: skipping node without ObjectType\n");
					goto next;							/* Skip this one */
				}
	
#ifdef NEVER
				/* Check the attribute value (Should we ?) */
				/* Known values are: Standard, Color */
				/* May be Trial, Target, Substrate, Colorant, ... ? */
				if (strcmp(attr, "Standard") != 0
				 && strcmp(attr, "Color") != 0) {
					a1logd(p->log, DEB6, "read_cxf: skipping node with ObjectType = '%s'\n",attr);
					goto next;							/* Skip this one */
				}
#endif
			}

			if ((attr = mxmlElementGetAttr(node, SampleNameKey)) == NULL) {
				a1logd(p->log, DEB6, "read_cxf: skipping without %s\n",SampleNameKey);
				goto next;							/* Skip this one */
			}
			name = attr;

			a1logd(p->log, DEB4, "read_cxf: got color %d name '%s'\n",i,name);

			if (!cxf2) {
				if ((ppvals = mxmlFindElement(node, node, pfx(p,"ColorValues"), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
					a1logd(p->log, DEB4, "read_cxf: no reference ColorValuesname  - skipping\n");
					goto next;
				}
			} else {
				ppvals = node;
			}

			/* Read any spectral values */
			if ((p->options & NAMEDC_OP_NOSPEC) == 0
			 && (pvals = mxmlFindElement(ppvals, ppvals, pfx(p,"ReflectanceSpectrum"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
				const char *elem = NULL;

				a1logd(p->log, 4, "read_cxf: got ReflectanceSpectrum\n");

				/* Check which color specification is being used with XYZ */
				if (SpectSpecification == NULL) {
					SpectSpecification = mxmlElementGetAttr(pvals, "ColorSpecification");
					a1logd(p->log, 4, "read_cxf: got SpectSpecification '%s'\n",SpectSpecification); 
				}

				if ((elem = mxmlElementGetAttr(pvals, "StartWL")) != NULL) {
					int minwl = 0;
					a1logd(p->log, 4, "read_cxf: got StartWL '%s'\n",elem); 
					minwl = atoi(elem);

					if (specmin != 0) {
						if (specmin != minwl) {
							a1logd(p->log,1,"Inconsistent StartWL %d vs. %d from '%s'\n",specmin,minwl,p->filename);
							goto skip_spectral;
						}
					} else {
						specmin = minwl;
					}
				}

				if ((elem = mxmlElementGetAttr(pvals, "Increment")) != NULL) {
					int inc = 0;
					a1logd(p->log, 4, "read_cxf: got Increment '%s'\n",elem); 
					inc = atoi(elem);

					if (specinc != 0) {
						if (specinc != inc) {
							a1logd(p->log,1,"Inconsistent Increment %d vs. %d from '%s'\n",specinc,inc,p->filename);
							goto skip_spectral;
						}
					} else {
						specinc = inc;
					}
				}

				if (specmin == 0 || specinc == 0) {
					a1logd(p->log,1,"Missing: specmin %d specinc %d from '%s'\n",specmin,specinc,p->filename);
					goto skip_spectral;
				}

				/* mxml stashes multiple elements as children.. */
				if ((val = mxmlGetFirstChild(pvals)) == NULL) {
					a1logd(p->log,1,"No spectral values in '%s'\n",p->filename);
					goto skip_spectral;
				}

				if ((sp = calloc(1, sizeof(xspect))) == NULL) {
					snprintf(p->err, NAMEDC_ERRL, "Malloc of xspect failed");
					a1logd(p->log, 1, "read_cxf: %s\n",p->err);
					mxmlDelete(tree);
					return p->errc = 2;
				}

				for (j = 0; j < XSPECT_MAX_BANDS; j++) {
					sp->spec[j] = 100.0 * mxmlGetReal(val);
					a1logd(p->log, 6, "read_cxf: got spect component %d value %f\n",j,sp->spec[j]);
					val = mxmlGetNextSibling(val);
					if (val == NULL)
						break;
				}
				sp->spec_n = j+1;
				sp->spec_wl_short = (double)specmin;
				sp->spec_wl_long = (double)(specmin + (sp->spec_n-1) * specinc);
				sp->norm = 100.0;

				a1logd(p->log, 6, "read_cxf: wl num %d short %f long %f\n",sp->spec_n,sp->spec_wl_short,sp->spec_wl_long);
			}
		  skip_spectral:;

			/* See if there is ColorCIELab */
			if ((pvals = mxmlFindElement(ppvals, ppvals, pfx(p,SampleLabKey), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
				char *key[3] = { "L", "A", "B" };

				/* Check which color specification is being used with Lab */
				if (!cxf2 && colorSpecification == NULL) {
					colorSpecification = mxmlElementGetAttr(pvals, "ColorSpecification");
				}

				for (j = 0; j < 3; j++) {
					if ((val = mxmlFindElement(pvals, pvals, pfx(p, key[j]), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
						a1logd(p->log, DEB6, "read_cxf: failed to find ColorCIELab component %s\n",key[j]);
						break;		/* oops */
					}
					Lab[j] = mxmlGetReal(val);
					a1logd(p->log, DEB6, "read_cxf: got ColorCIELab component %s value %f\n",key[j],Lab[j]);
				}
				if (j >= 3)
					Lab_v = 1;
			}

			if (!Lab_v) {

				a1logd(p->log, DEB6, "read_cxf: no valid ColorCIELab value, look for XYZ\n");

				/* See if there is ColorCIEXYZ instead */
				if ((pvals = mxmlFindElement(ppvals, ppvals, pfx(p,SampleXYZKey), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
					char *key[3] = { "X", "Y", "Z" };

					/* Check which color specification is being used with Lab */
					if (!cxf2 && colorSpecification == NULL) {
						colorSpecification = mxmlElementGetAttr(pvals, "ColorSpecification");
					}

					for (j = 0; j < 3; j++) {
						if ((val = mxmlFindElement(pvals, pvals, pfx(p, key[j]), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
							a1logd(p->log, DEB6, "read_cxf: failed to find ColorCIEXYZ component %s\n",key[j]);
							break;		/* oops */
						}
						Lab[j] = mxmlGetReal(val);
						a1logd(p->log, DEB6, "read_cxf: got ColorCIEXYZ component %s value %f\n",key[j],Lab[j]);
					}
					if (j >= 3) {
						icmXYZ2Lab(&icmD50, Lab, Lab);
						Lab_v = 1;
					}
				}
			}

			if (sp == NULL && !Lab_v) {
				a1logd(p->log, DEB6, "read_cxf: no spectral or CIE value found - skipping color\n");
				goto next;
			}
		
			/* Read any device color values */
			if (!cxf2) {
				if ((ppvals = mxmlFindElement(node, node, pfx(p,"DeviceColorValues"), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {

					a1logd(p->log, DEB4, "read_cxf: got DeviceColorValues\n");
				}
			}

			if (ppvals != NULL) {

				// ~~99 could add read of other device values here

				/* See if there is ColorCMYK */
				if ((pvals = mxmlFindElement(ppvals, ppvals, pfx(p,SampleCMYKKey), NULL, NULL, MXML_DESCEND_FIRST)) != NULL) {
					char *key[4] = { "Cyan", "Magenta", "Yellow", "Black" };

					a1logd(p->log, DEB4, "read_cxf: got ColorCMYK\n");
					for (j = 0; j < 4; j++) {
						if ((val = mxmlFindElement(pvals, pvals, pfx(p, key[j]), NULL, NULL, MXML_DESCEND_FIRST)) == NULL) {
							a1logd(p->log, DEB6, "read_cxf: failed to find ColorCMYK component %s\n",key[j]);
							break;		/* oops */
						}
						dev[j] = mxmlGetReal(val);
						a1logd(p->log, DEB6, "read_cxf: got ColorCMYK component %s value %f\n",key[j],dev[j]);
					}
					if (j >= 4) {
						dev_n = j;
						devSig = icSigCmykData;
					}
				}
			}
		
			/* Add an entry */
			if (p->count >= p->count_a) {
				unsigned int count_n;
				count_n = p->count_a + 4096/sizeof(nce);

				a1logd(p->log, 8, "read_cxf: increasing data array size to %d\n",count_n);
				if ((p->data = recalloc(p->data, p->count_a, sizeof(nce), count_n, sizeof(nce))) == NULL) {
					snprintf(p->err, NAMEDC_ERRL, "Malloc of data size %d failed",p->count_a);
					a1logd(p->log, 1, "read_cxf: %s\n",p->err);
					mxmlDelete(tree);
					return p->errc = 2;
				}
				p->count_a = count_n;
			}
			clear_nce(&p->data[p->count]);

			if ((p->data[p->count].name = strdup(name)) == NULL) {
				snprintf(p->err, NAMEDC_ERRL, "Malloc of color name string failed");
				a1logd(p->log, 1, "read_cxf: %s\n",p->err);
				mxmlDelete(tree);
				return p->errc = 2;
			}
			
			if (sp != NULL) {
				p->data[p->count].sp = sp;
			}

			if (Lab_v) {
				p->data[p->count].Lab[0] = Lab[0];
				p->data[p->count].Lab[1] = Lab[1];
				p->data[p->count].Lab[2] = Lab[2];
				p->data[p->count].Lab_v = 1;

			} else {

				/* match() expects lab values */
				if (sp2cie == NULL) {
					if ((sp2cie = new_xsp2cie(p->ill, 0.0, NULL, p->obs, NULL, icSigLabData, 0)) == NULL) {
						snprintf(p->err, NAMEDC_ERRL, "creating spectral conversion failed");
						a1logd(p->log, 1, "read_cxf: %s\n",p->err);
						mxmlDelete(tree);
						return p->errc = 2;
					}
				}
				sp2cie->convert(sp2cie, p->data[p->count].Lab, p->data[p->count].sp); 
				p->data[p->count].Lab_v = 1;
			}

			if (dev_n > 0 && devSig != icMaxEnumData) {
				for (j = 0; j < dev_n; j++)
					p->data[p->count].dev[j] = dev[j];
				p->data[p->count].dev_n = dev_n;
				p->data[p->count].devSig = devSig;
			} else {
				p->data[p->count].dev_n = 0;
				p->data[p->count].devSig = icMaxEnumData;
			}

			a1logd(p->log, 8, "read_cxf: added color %d\n",p->count);

			p->count++;

		next:;	/* Next color */
			node = mxmlGetNextSibling(node);
		}

		/* - - - - - - - - - - - - */
		/* Read extra information */
		if (cxf2) {
			/* Grab the creator */
			if ((node = mxmlFindPathNode(cxf, pfxp(p,"Preamble/Header/Creator"))) == NULL)
				name = NULL;
			else
				name = mxmlGetOpaque(node);
	
			if (name == NULL)
				name = "[Unknown]";
			if ((p->creator = strdup(name)) == NULL) {
				snprintf(p->err, NAMEDC_ERRL, "Malloc of creator string failed");
				a1logd(p->log, 1, "read_cxf: %s\n",p->err);
				mxmlDelete(tree);
				return p->errc = 2;
			}
			a1logd(p->log, 2, "read_cxf: creator '%s'\n",p->creator);

			/* Grab the illuminant type */
			if ((node = mxmlFindPathNode(cxf, pfxp(p,"Palette/ColorSet/CollectionColorSpaceSpecification/ColorSpaceSpecificationSpectrumTristimulus/IlluminationOptions/Illuminant"))) == NULL
			 || (name = mxmlGetOpaque(node)) == NULL) {
				a1logd(p->log, 2, "read_cxf: failed to locate Illuminant - assuming D50\n");
				p->ill = icxIT_D50;
			} else {
				if (strcmp(name, "Illuminant_A") == 0)
					p->ill = icxIT_A;
				else if (strcmp(name, "Illuminant_D50") == 0)
					p->ill = icxIT_D50;
				else if (strcmp(name, "Illuminant_D55") == 0)
					p->ill = icxIT_D55;
				else if (strcmp(name, "Illuminant_D65") == 0)
					p->ill = icxIT_D65;
				else if (strcmp(name, "Illuminant_E") == 0)
					p->ill = icxIT_E;
				else {
					a1logd(p->log, 2, "read_cxf: Illuminant '%s' unrecognised\n",name);
					p->ill = icxIT_D50;
				}
			}
			a1logd(p->log, 2, "read_cxf: illuminant '%s'\n",icm2str(icmIlluminant, p->ill));

			/* Grab the observer type */
			if ((node = mxmlFindPathNode(cxf, pfxp(p,"Palette/ColorSet/CollectionColorSpaceSpecification/ColorSpaceSpecificationSpectrumTristimulus/FieldOfView"))) == NULL
			 || (name = mxmlGetOpaque(node)) == NULL) {
				a1logd(p->log, 2, "read_cxf: failed to locate FieldOfView - assuming 2 degree\n");
				p->obs = icxOT_CIE_1931_2;
			} else {
				if (strcmp(name, "FieldOfView_2_Degree") == 0)
					p->obs = icxOT_CIE_1931_2;
				else if (strcmp(name, "FieldOfView_10_Degree") == 0)
					p->obs = icxOT_CIE_1964_10;
				else {
					a1logd(p->log, 2, "read_cxf: Observer '%s' unrecognised\n",name);
					p->obs = icxOT_CIE_1931_2;
				}
			}
			a1logd(p->log, 2, "read_cxf: observer '%s'\n",icm2str(icmStandardObserver, p->obs));

		} else {
			int found_io = 0;

			/* Grab the creator */
			if ((node = mxmlFindPathNode(cxf, pfxp(p,"FileInformation/Creator"))) == NULL)
				name = NULL;
			else
				name = mxmlGetOpaque(node);
	
			if (name == NULL)
				name = "[Unknown]";
			if ((p->creator = strdup(name)) == NULL) {
				snprintf(p->err, NAMEDC_ERRL, "Malloc of creator string failed");
				a1logd(p->log, 1, "read_cxf: %s\n",p->err);
				mxmlDelete(tree);
				return p->errc = 2;
			}
			a1logd(p->log, 2, "read_cxf: creator '%s'\n",p->creator);

			if (colorSpecification != NULL) {
				a1logd(p->log, 2, "read_cxf: colorSpecification '%s'\n",colorSpecification);
			}

			/* Look through the color specifications and find the one that matches */
			/* the Lab or XYZ color specification */
			pnode = mxmlFindPathNode(cxf, pfxp(p,"Resources/ColorSpecificationCollection/ColorSpecification"));
			while (pnode != NULL) {
				name = mxmlElementGetAttr(pnode, "Id");

				if (colorSpecification == NULL
				 || (name != NULL && strcmp(name, colorSpecification) == 0)) {
		
					a1logd(p->log, 2, "read_cxf: found node matching colorSpecification '%s'\n",colorSpecification);

					/* Grab the illuminant type */
					if ((node = mxmlFindPathNode(pnode, pfxp(p,"TristimulusSpec/Illuminant"))) == NULL
					 || (name = mxmlGetOpaque(node)) == NULL) {
						a1logd(p->log, 2, "read_cxf: failed to locate Illuminant - assuming D50\n");
						p->ill = icxIT_D50;
					} else {
						if (strcmp(name, "A") == 0)
							p->ill = icxIT_A;
						else if (strcmp(name, "D50") == 0)
							p->ill = icxIT_D50;
						else if (strcmp(name, "D55") == 0)
							p->ill = icxIT_D55;
						else if (strcmp(name, "D65") == 0)
							p->ill = icxIT_D65;
						else if (strcmp(name, "E") == 0)
							p->ill = icxIT_E;
						else {
							a1logd(p->log, 2, "read_cxf: Illuminant '%s' unrecognised\n",name);
							p->ill = icxIT_D50;
						}
					}
					a1logd(p->log, 2, "read_cxf: illuminant '%s'\n",icm2str(icmIlluminant, p->ill));
		
					/* Grab the first observer type */
					if ((node = mxmlFindPathNode(cxf, pfxp(p,"Resources/ColorSpecificationCollection/ColorSpecification/TristimulusSpec/Observer"))) == NULL
					 || (name = mxmlGetOpaque(node)) == NULL) {
						a1logd(p->log, 2, "read_cxf: failed to locate Observer - assuming 2 degree\n");
						p->obs = icxOT_CIE_1931_2;
					} else {
						if (strcmp(name, "2_Degree") == 0)
							p->obs = icxOT_CIE_1931_2;
						else if (strcmp(name, "10_Degree") == 0)
							p->obs = icxOT_CIE_1964_10;
						else {
							a1logd(p->log, 2, "read_cxf: Observer '%s' unrecognised\n",name);
							p->obs = icxOT_CIE_1931_2;
						}
					}
					a1logd(p->log, 2, "read_cxf: observer '%s'\n",icm2str(icmStandardObserver, p->obs));
					found_io = 1;
					break;
				}
				pnode = mxmlGetNextSibling(pnode);
			}

			if (!found_io) {
				p->ill = icxIT_D50;
				p->obs = icxOT_CIE_1931_2;
				a1logd(p->log, 2, "read_cxf: failed to locate ColorSpecification - assuming D50 2 degree observer\n");
			}

			/* Other values of interest:
			   Resources/ColorSpecificationCollection/ColorSpecification/ColorSpecification/MeasurementSpec/MeasurementType	ie. "Spectrum_Reflectance"
			   Resources/ColorSpecificationCollection/ColorSpecification/ColorSpecification/MeasurementSpec/CalibrationStandard	ie. "GMDI" for Gretag Macbeth Calibration, "XRGA" for X-Rite Graphic Arts Standard.
			   Resources/ColorSpecificationCollection/ColorSpecification/ColorSpecification/MeasurementSpec/Device/DeviceFilter	ie. "Filter_None" "Filter_UVExcluded", 
			 */
		}
	}

	mxmlDelete(tree);

	a1logd(p->log, 1, "read_cxf: done - %d colors\n",p->count);
	return 0;
}

/* Read in a namedc from a ICC file */
/* Return nz on error */
static int read_icc(namedc *p, const char *filename, int options) {
	char *pfilename;
	icmFile *fp;
	icc *icco;

	a1logd(p->log, 1, "read_icc: file '%s' options 0x%x\n",filename,options);

	if (p->filename == NULL) {
		if ((pfilename = strdup(filename)) == NULL) {
			snprintf(p->err, NAMEDC_ERRL, "Malloc of filename failed");
			a1logd(p->log, 1, "read_icc: %s\n",p->err);
			return p->errc = 2;
		}
		clear_namedc(p);
		p->filename = pfilename;
	}
	p->options = options;
	
	/* Open up the file for reading */
	if ((fp = new_icmFileStd_name(p->filename,"r")) == NULL) {
		snprintf(p->err, NAMEDC_ERRL, "Opening ICC file '%s' failed with %s",
		                                     p->filename,strerror(errno));
		a1logd(p->log, 1, "read_icc: %s\n",p->err);
		return p->errc = 1;
	}

	if ((icco = new_icc()) == NULL) {
		snprintf(p->err, NAMEDC_ERRL, "Creation ICC object failed");
		a1logd(p->log, 1, "read_icc: %s\n",p->err);
		fp->del(fp);
		return p->errc = 1;
	}

	if (icco->read(icco, fp, 0) != 0) {
		snprintf(p->err, NAMEDC_ERRL, "Failed to read '%s'\n",p->filename);
		a1logd(p->log, 1, "read_icc: %s\n",p->err);
		icco->del(icco);
		fp->del(fp);
		return p->errc = 1;
	}

	if (icco->header->deviceClass != icSigNamedColorClass) {
		snprintf(p->err, NAMEDC_ERRL, "Not a named color profile");
		a1logd(p->log, 1, "read_icc: %s\n",p->err);
		icco->del(icco);
		fp->del(fp);
		return p->errc = 1;
	}

	if (icco->header->pcs != icSigXYZData
	 && icco->header->pcs != icSigLabData) {
		snprintf(p->err, NAMEDC_ERRL, "Unrecognised PCS");
		a1logd(p->log, 1, "read_icc: %s\n",p->err);
		icco->del(icco);
		fp->del(fp);
		return p->errc = 1;
	}

	if (icco->find_tag(icco, icSigNamedColor2Tag) != 0) {
		snprintf(p->err, NAMEDC_ERRL, "Can't find ncl2 tag");
		a1logd(p->log, 1, "read_icc: %s\n",p->err);
		icco->del(icco);
		fp->del(fp);
		return p->errc = 1;
	}

	{
		icmTextDescription *txd;

		/* Try and read the tag from the file */
		if ((txd = (icmTextDescription *)icco->read_tag(icco, icSigProfileDescriptionTag)) == NULL) {
			snprintf(p->err, NAMEDC_ERRL, "No description tag");
			a1logd(p->log, 1, "read_icc: %s\n",p->err);
			icco->del(icco);
			fp->del(fp);
			return p->errc = 1;
		}

		if ((p->description = strdup(txd->desc)) == NULL) {
			snprintf(p->err, NAMEDC_ERRL, "Malloc of description string failed");
			a1logd(p->log, 1, "read_icc: %s\n",p->err);
			icco->del(icco);
			fp->del(fp);
			return p->errc = 2;
		}
		a1logd(p->log, 2, "read_icc: description '%s'\n",p->description);
	}

	p->hash = do_hash2(p->filename, p->description);

	if ((p->options & NAMEDC_OP_NODATA) == 0) {
		icmNamedColor *tag;
		char name[3 * 32];
		double Lab[3];
		int Lab_v = 0;
		double dev[MAX_CHAN];
		int dev_n = 0;
		icColorSpaceSignature devSig = icMaxEnumData;
		int i, j;
		
#ifdef NEVER
		/* See if there is a measurementType tag */
		icmMeasurement *meastag;
		if ((meastag = (icmMeasurement *)icco->read_tag(icco, icSigMeasurementTag)) != NULL) {
			p->ill = meastag->illuminant;
			p->obs = meastag->observer;
			a1logd(p->log, 2, "read_cxf: assuming D50 illuminant and 2 degree observer\n");
		} else {
			p->ill = icIlluminantD50;
			p->obs = icStdObs1931TwoDegrees;
		}
		a1logd(p->log, 2, "read_cxf: illuminant '%s'\n",icm2str(icmIlluminant, p->ill));
		a1logd(p->log, 2, "read_cxf: observer '%s'\n",icm2str(icmStandardObserver, p->obs));
#endif

		if ((tag = (icmNamedColor *)icco->read_tag(icco, icSigNamedColor2Tag)) == NULL) {
			snprintf(p->err, NAMEDC_ERRL, "Can't read ncl2 tag");
			a1logd(p->log, 1, "read_icc: %s\n",p->err);
			icco->del(icco);
			fp->del(fp);
			return p->errc = 1;
		}

		for (i = 0; i < tag->count; i++) {

			/* Get the name */
			strcpy(name, tag->prefix);
			strcat(name, tag->data[i].root);
			strcat(name, tag->suffix);

			a1logd(p->log, DEB4, "read_icc: got color %d name '%s'\n",i,name);

			if (icco->header->pcs == icSigXYZData) {
				Lab[0] = tag->data[i].pcsCoords[0];
				Lab[1] = tag->data[i].pcsCoords[1];
				Lab[2] = tag->data[i].pcsCoords[2];
				icmXYZ2Lab(&icmD50, Lab, Lab);

			} else {		/* Lab */
				Lab[0] = tag->data[i].pcsCoords[0];
				Lab[1] = tag->data[i].pcsCoords[1];
				Lab[2] = tag->data[i].pcsCoords[2];
			}
			Lab_v = 1;
			a1logd(p->log, DEB6, "read_icc: got ColorCIELab value %s\n",icmPdv(3, Lab));

			/* Add any device space */
			devSig = icco->header->colorSpace;
			dev_n = icmCSSig2nchan(devSig);
			for (j = 0; j < dev_n; j++)
				dev[j] = tag->data[i].deviceCoords[j];

			a1logd(p->log, DEB6, "read_icc: got %s value %s\n",
			       icm2str(icmColorSpaceSignature, devSig), icmPdv(dev_n, dev));

			/* Add an entry */
			if (p->count >= p->count_a) {
				unsigned int count_n;
				count_n = p->count_a + 4096/sizeof(nce);

				a1logd(p->log, 8, "read_icc: increasing data array size to %d\n",count_n);
				if ((p->data = recalloc(p->data, p->count_a, sizeof(nce), count_n, sizeof(nce))) == NULL) {
					snprintf(p->err, NAMEDC_ERRL, "Malloc of data size %d failed",p->count_a);
					a1logd(p->log, 1, "read_icc: %s\n",p->err);
					icco->del(icco);
					fp->del(fp);
					return p->errc = 2;
				}
				p->count_a = count_n;
			}
			clear_nce(&p->data[p->count]);

			if ((p->data[p->count].name = strdup(name)) == NULL) {
				snprintf(p->err, NAMEDC_ERRL, "Malloc of color name string failed");
				a1logd(p->log, 1, "read_icc: %s\n",p->err);
				icco->del(icco);
				fp->del(fp);
				return p->errc = 2;
			}
			
			if (Lab_v) {
				p->data[p->count].Lab[0] = Lab[0];
				p->data[p->count].Lab[1] = Lab[1];
				p->data[p->count].Lab[2] = Lab[2];
				p->data[p->count].Lab_v = 1;
			}

			if (dev_n > 0 && devSig != icMaxEnumData) {
				for (j = 0; j < dev_n; j++)
					p->data[p->count].dev[j] = dev[j];
				p->data[p->count].dev_n = dev_n;
				p->data[p->count].devSig = devSig;
			} else {
				p->data[p->count].dev_n = 0;
				p->data[p->count].devSig = icMaxEnumData;
			}

			a1logd(p->log, 8, "read_icc: added color %d\n",p->count);

			p->count++;

		next:;	/* Next color */
		}
	}

	icco->del(icco);
	fp->del(fp);

	return 0;
}

/* Read any format named color files */
/* Return nz on error */
static int read_nc(namedc *p, const char *filename, int options) {
	int rv;

	if ((p->format == 0 || p->format == 1)
	 && (rv = read_cxf(p, filename, options)) == 0) {
		p->format = 1;
		return rv;
	} 

	if ((p->format == 0 || p->format == 2)
	 && (rv = read_icc(p, filename, options)) == 0) {
		p->format = 2;
		return rv;
	} 

	/* Try other file types here */
	
	return rv;
}


/* Return the index of the best mataching color, -1 on error. */
/* Lab[] is assumed to be D50, 2 degree standard observer based CIE value, */
/* and the spec value should only be provided if this is a reflective or */
/* transmissive measurement, NULL if emissive. */
/* If named color library is expects other than D50, 2 degree, then */
/* it will use the spectral value if not NULL, or chromatically */
/* adapt the Lab value. */
/* deType == 0 DE76 */
/* deType == 1 DE94 */
/* deType == 2 DE2000 */
/* if de != NULL, return the delta E */
int match(struct _namedc *p, double *de, double *pLab, xspect *rspect, int deType) {
	int i, bix = -1;
	double bde = 1e99;
	double Lab[3];

	if (p->filename == NULL) {		/* We haven't been opened */
		snprintf(p->err, NAMEDC_ERRL, "We haven't been opened");
		a1logd(p->log, 1, "match: %s\n",p->err);
		return -1;
	}

	/* If the colors haven't been read yet, read them now */
	if (p->data == NULL || (p->options & NAMEDC_OP_NODATA)) {
		if (read_nc(p, NULL, (p->options & ~NAMEDC_OP_NODATA))) {
			a1logd(p->log, 1, "match: on demand data load failed with '%s'\n",p->err);
			return -1;
		}
		a1logd(p->log, 1, "match: after loading there are %d colors\n",p->count);
	}

	icmCpy3(Lab, pLab);

	if (p->ill != icIlluminantD50 || p->obs != icStdObs1931TwoDegrees) {
		if (rspect != NULL) {

			if (p->sp2cie == NULL) {
				if ((p->sp2cie = new_xsp2cie(p->ill, 0.0, NULL, p->obs, NULL, icSigLabData, 0)) == NULL) {
					snprintf(p->err, NAMEDC_ERRL, "creating spectral conversion failed");
					a1logd(p->log, 1, "match: %s\n",p->err);
					return -1;
					
				}
			}

			/* Convert spectrum to the Lab we want */
			p->sp2cie->convert(p->sp2cie, Lab, rspect); 

		} else {
			if (p->chrom[0][0] <= -1e38) {
				double wXYZ[3];

				// Special case this for a consistent value with ICC profiles
				if (p->obs == icxOT_CIE_1931_2 && p->ill == icxIT_D65) {
					icmCpy3(wXYZ, icmD65_ary3);
				} else { 
					xsp2cie *tt;
					xspect ts;
					// Get the XYZ of the given white point for the illuminant and observer
					if ((tt = new_xsp2cie(p->ill, 0.0, NULL, p->obs, NULL, icSigXYZData, 0)) == NULL) {
						snprintf(p->err, NAMEDC_ERRL, "creating spectral conversion failed");
						a1logd(p->log, 1, "match: %s\n",p->err);
						return -1;
					}
					if (standardIlluminant(&ts, icxIT_E, 0.0)) {
						snprintf(p->err, NAMEDC_ERRL, "match: creating E type spectrum failed");
						a1logd(p->log, 1, "match: %s\n",p->err);
						return -1;
					} 
					tt->convert(tt, wXYZ, &ts);
					tt->del(tt);
				}
				icmAry2XYZ(p->dXYZ, wXYZ);

				/* Chreate chromatic adapation matrix from D50 to named color */
				icmChromAdaptMatrix(ICM_CAM_BRADFORD, p->dXYZ, icmD50, p->chrom);
			}
			icmLab2XYZ(&icmD50, Lab, pLab);
			icmMulBy3x3(Lab, p->chrom, Lab);
			icmXYZ2Lab(&p->dXYZ, Lab, Lab);
		}
	}
 
	if (deType == 0) {
		for (i = 0; i < p->count; i++) {
			double de = icmLabDEsq(Lab, p->data[i].Lab);
			if (de < bde) {
				bde = de;
				bix = i;
			}
		}

	} else if (deType == 1) {
		for (i = 0; i < p->count; i++) {
			double de = icmCIE94sq(Lab, p->data[i].Lab);
			if (de < bde) {
				bde = de;
				bix = i;
			}
		}

	} else if (deType == 2) {
		for (i = 0; i < p->count; i++) {
			double de = icmCIE2Ksq(Lab, p->data[i].Lab);
			if (de < bde) {
				bde = de;
				bix = i;
			}
		}
	} else {
		snprintf(p->err, NAMEDC_ERRL, "Unnown deType %d",deType);
		a1logd(p->log, 1, "match: %s\n",p->err);
		return -1;
	}

	if (bix < 0) {
		snprintf(p->err, NAMEDC_ERRL, "No colors to match against");
		a1logd(p->log, 1, "match: %s\n",p->err);
		return -1;
	}
	if (de != NULL) {
		*de = sqrt(bde);
	}
	return bix;
}

/* Free an entry */
static void clear_nce(nce *p) {
	if (p != NULL) {
		if (p->name != NULL)
			free(p->name);
		p->name = NULL;
		if (p->sp != NULL)
			free(p->sp);
		p->sp = NULL;
	}
}

/* Free the contents */
static void clear_namedc(namedc *p) {
	if (p != NULL) {
		int i;
		if (p->filename != NULL)
			free(p->filename);
		p->filename = NULL;
		if (p->creator != NULL)
			free(p->creator);
		p->creator = NULL;
		if (p->description != NULL)
			free(p->description);
		p->description = NULL;
		if (p->data != NULL) {
			for (i = 0; i < p->count; i++)
				clear_nce(&p->data[i]);
			free(p->data);
			p->data = NULL;
		}
		p->count = 0;
	}
}

/* Delete it */
static void del_namedc(namedc *p) {
	if (p != NULL) {
		clear_namedc(p);
		p->log = del_a1log(p->log);
		if (p->sp2cie != NULL)
			p->sp2cie->del(p->sp2cie);
		free(p);
	}
}

/* Allocate a new, uninitialised namedc */
/* Note thate black and white points aren't allocated */
namedc *new_namedc(a1log *log) {
	namedc *p;

	a1logd(log, 1, "new_cxf\n");

	if ((p = (namedc *)calloc(1, sizeof(namedc))) == NULL) {
		a1logd(p->log, 1, "new_cxf: calloc failed\n");
		return NULL;
	}

	p->log = new_a1log_d(log);

	/* Init method pointers */
	p->del        = del_namedc;
	p->read_cxf   = read_cxf;
	p->read_icc   = read_icc;
	p->read       = read_nc;
	p->match      = match;

	p->chrom[0][0] = -1e38;

	return p;
}


/* =========================================================================== */

#else /* STANDALONE_TEST */

void usage(char *diag, ...) {
	fprintf(stderr,"Test namedc library\n");
	fprintf(stderr,"Author: Graeme W. Gill\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: namedc [-D level] infile\n");
	fprintf(stderr," -D level               Debug level 1-9\n");
	exit(1);
}

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char inname[MAXNAMEL+1];
	int debug = 0, i;
	namedc *p;

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

			/* Debug level */
			else if (argv[fa][1] == 'D') {
				fa = nfa;
				if (na != NULL) 
					debug = atoi(na);
			}

			else 
				usage("Unknown option '%c'",argv[fa][1]);
		}
		else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage("Missing input filename");
	strncpy(inname,argv[fa++],MAXNAMEL); inname[MAXNAMEL] = '\000';

	g_log->debug = debug;

	if ((p = new_namedc(g_log)) == NULL)
		error("new_namedc failed\n");

	/* Open non-data */
	if (p->read(p, inname, NAMEDC_OP_NODATA)) {
		error("read failed with '%s'\n",p->err);
	}

	printf("Palette is '%s'\n",p->description);
	{
		double Lab[3], de;
		int ix;
		
		Lab[0] = 50.0;
		Lab[1] = 20.0;
		Lab[2] = -10.0;

		if ((ix = p->match(p, &de, Lab, NULL, 0)) < 0)
			error(" match failed with '%s'\n",p->err);
		printf("Matched color '%s' with DE76 %f\n",p->data[ix].name,de);

		if ((ix = p->match(p, &de, Lab, NULL, 1)) < 0)
			error(" match failed with '%s'\n",p->err);
		printf("Matched color '%s' with DE94 %f\n",p->data[ix].name,de);

		if ((ix = p->match(p, &de, Lab, NULL, 2)) < 0)
			error(" match failed with '%s'\n",p->err);
		printf("Matched color '%s' with DE00 %f\n",p->data[ix].name,de);
	}
	
#ifdef NEVER
	printf("Loaded %d colors\n",p->count);
	for (i = 0; i < p->count; i++) {
		printf("Color %d name '%s' = %f %f %f\n",
		       i, p->data[i].name, p->data[i].Lab[0], p->data[i].Lab[1], p->data[i].Lab[2]);
		if (p->data[i].devSig != icMaxEnumData) {
			printf("%s = %s\n", icm2str(icmColorSpaceSignature, p->data[i].devSig),
			             icmPdv(p->data[i].dev_n, p->data[i].dev));
		}
	}
#endif
	p->del(p);

	return 0;
}



#endif /* STANDALONE_TEST */
























