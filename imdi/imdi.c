
/* Integer Multi-Dimensional Interpolation */

/*
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This is the implementation of the run time code.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "imdi.h"
#include "imdi_tab.h"
#include "imdi_k.h"			/* Declaration of all the kernel functions */

#undef VERBOSE
#undef VVERBOSE

static unsigned int imdi_get_check(imdi *im);
static void imdi_reset_check(imdi *im);
static void imdi_info(imdi *s, unsigned long *size, int *gres, int *sres);
static void imdi_del(imdi *im);
static void interp_match(imdi *s, void **outp, int outst, void **inp, int inst,
                         unsigned int npixels);


/* Create a new imdi */
/* Return NULL if request is not supported */
/* Note that we use the high level pixel layout description to locate */
/* a suitable run-time. */
imdi *new_imdi(
	int id,				  /* Number of input dimensions */
	int od,				  /* Number of output lookup dimensions */
	                      /* Number of output channels written = od - no. of oopt skip flags */
	imdi_pixrep in,		  /* Input pixel representation */
	int in_signed,		  /* Bit flag per channel, NZ if treat as signed */
	int *inm,			  /* Input raster channel to callback channel mapping, NULL for none. */
	imdi_iprec iprec,	  /* Internal processing precision */
	imdi_pixrep out,	  /* Output pixel representation */
	int out_signed,		  /* Bit flag per channel, NZ if treat as signed */
	int *outm,			  /* Output raster channel to callback channel mapping, NULL for none. */
	                      /* Mapping must include skipped channels. */
	int res,			  /* Desired table resolution */
	imdi_ooptions oopt,   /* Output per channel options (by callback channel) */
	unsigned int *checkv, /* Output channel check values (by callback channel, NULL == 0) */
	imdi_options opt,	  /* Direction and stride options */

	/* Callbacks to lookup the imdi table values. */
	/* (Skip output channels are looked up) */
	void (*input_curves) (void *cntx, double *out_vals, double *in_vals),
	void (*md_table)     (void *cntx, double *out_vals, double *in_vals),
	void (*output_curves)(void *cntx, double *out_vals, double *in_vals),
	void *cntx		/* Context to callbacks */
) {
	int i;
	int prec;					/* Target internal precision */
	int bk = -1;				/* Best table */
	int bfig = 0x7fffffff;		/* Best tables figure of merit */
	int bstres = 0;				/* Best tables target stres */
	genspec gs;					/* Generation specification */
	tabspec ts;					/* Table specifications */
	genspec bgs;				/* Best gen spec */
	tabspec bts;				/* Best tab spec */
	imdi_conv bcnv = conv_none;	/* Best tables conversion flags */
	imdi_ooptions Ooopt;		/* oopt re-aranged to correspond to output channel index */
	
	imdi *im;

	/* Compute the Output channel index oopt mask */
	if (outm == NULL)
		Ooopt = oopt;
	else {
		int ff;
		Ooopt = 0;
		for (i = 0; i < od; i++) {
			ff = OOPTX(oopt,outm[i]);
			Ooopt |= OOPT(ff,i);
		} 
	}

#ifdef VERBOSE
	printf("new_imdi called with id %d, od %d, res %d Ooopt 0x%x, opt 0x%x\n", id, od, res, Ooopt, opt);
	printf("about to checking %d kernels\n", no_kfuncs);
#endif

	/* Figure out the internal precision requested */
	COMPUTE_IPREC(prec, in, iprec, out)

	/* Zero out the gen and tabspecs, to give diff a place to start */
	memset((void *)&gs, 0, sizeof(genspec));
	memset((void *)&ts, 0, sizeof(tabspec));

	/* The first thing to do is see if there is an available kernel function */
	for (i = 0; i < no_kfuncs; i++) {
		int stres;					/* Computed stres needed */
		imdi_conv cnv = conv_none;	/* Conversions needed for this choice */
		int fig;					/* Figure of merit - smaller is better */
		ktable[i].gentab(&gs, &ts);	/* Udate the kernel functions genspec and tabspec */

#ifdef VERBOSE
	printf("\n");
	printf("kernel %d has id %d, od %d, irep %d orep %d oopt 0x%x, opt 0x%x\n",
	        i, gs.id, gs.od, gs.irep, gs.orep, gs.oopt, gs.opt);
	printf("Input req is id %d, od %d, irep %d orep %d, oopt 0x%x, opt 0x%x\n",
	        id, od, in, out, Ooopt, opt);
#endif
		/* First check mandatory things */
		if (!(
			  id == gs.id			/* Input dimension */
		   && od == gs.od			/* Output dimension */
		)) {
#ifdef VVERBOSE
		printf("  Input or Output dimension mismatch\n");
#endif
			continue;
		}

		/* Check if we have an exact or conversion match for input stride */
		if (!(
		    (opt & opts_istride) == (gs.opt & opts_istride)
		 || (!(opt & opts_istride) && (gs.opt & opts_istride))
		)) {
#ifdef VVERBOSE
		printf("  Input stride mismatch\n");
#endif
			continue;
		}

		/* Check if we have an exact or conversion match for output stride */
		if (!(
		    (opt & opts_ostride) == (gs.opt & opts_ostride)
		 || (!(opt & opts_ostride) && (gs.opt & opts_ostride))
		)) {
#ifdef VVERBOSE
		printf("  Output stride mismatch\n");
#endif
			continue;
		}

		/* Check if we have an exact or conversion match for input representation */
		if (!(
		     in == gs.irep
		 || (in == pixint8 && gs.irep == planeint8 && (gs.opt & opts_istride))	
		 || (in == pixint16 && gs.irep == planeint16 && (gs.opt & opts_istride))
		)) {
#ifdef VVERBOSE
		printf("  Input representation mismatch\n");
#endif
			continue;
		}

		/* Check if we have an exact or conversion match for output representation */
		if (!(
		      out == gs.orep
		  || (out == pixint8 && gs.orep == planeint8 && (gs.opt & opts_ostride))	
		  || (out == pixint16 && gs.orep == planeint16 && (gs.opt & opts_ostride))
		)) {
#ifdef VVERBOSE
		printf("  Output representation mismatch\n");
#endif
			continue;
		}

		/* See if we have the output per channel options needed */
		if (!((Ooopt & gs.oopt) == Ooopt)) {
#ifdef VVERBOSE
		printf("  Output per channel options mismatch\n");
#endif
			continue;
		}

		/* See if we have an exact or conversion match for reverse direction */
		if (!(
		    ((!(opt & opts_bwd) && !(opt & opts_fwd)))		/* Don't care */
		 || (((opt & opts_bwd) && (gs.opt & opts_bwd)))		/* Match */
		 || (((opt & opts_fwd) && (gs.opt & opts_fwd)))		/* Match */
		 || ((gs.opt & opts_istride) && (gs.opt & opts_ostride))	/* Can convert */
		)) {
#ifdef VVERBOSE
		printf("  Direction mismatch\n");
#endif
			continue;
		}

#ifdef VERBOSE
		printf("  found match\n");
#endif
		fig = 0;

		/* Apply penalty if we will have to do a conversion match */
		if ((opt & opts_istride) != (gs.opt & opts_istride)) {
			cnv |= conv_istr;
			fig += 1000;
#ifdef VERBOSE
			printf("  needs input stride conversion though\n");
#endif
		}

		if ((opt & opts_ostride) != (gs.opt & opts_ostride)) {
			cnv |= conv_ostr;
			fig += 1000;
#ifdef VERBOSE
			printf("  needs output stride conversion though\n");
#endif
		}

		if (in != gs.irep) {
			cnv |= conv_irep;
			fig += 5000;
#ifdef VERBOSE
			printf("  needs input representation conversion though\n");
#endif
		}

		if (out != gs.orep) {
			cnv |= conv_orep;
			fig += 5000;
#ifdef VERBOSE
			printf("  needs output representation conversion though\n");
#endif
		}

		/* If we don't need output checking or skipping, but we've got it */
		if ((~Ooopt & gs.oopt) != 0) {
			unsigned int tmp = (~Ooopt & gs.oopt);
			/* Count number of bits set that we don't need */
  			/* (From HACKMEM 169) */
		    tmp = tmp - ((tmp >> 1) & 033333333333) - ((tmp >> 2) & 011111111111);
		    tmp = ((tmp + (tmp >> 3)) & 030707070707) % 63;
			fig += tmp;
#ifdef VERBOSE
			printf("  has output options we don't need though\n");
#endif
		}

		/* If we've got output skipping, we need to skip pointers in interp. */
		if ((Ooopt & OOPTS_SKIP) != 0) {
			cnv |= conv_skip;
		}

		if (((opt & opts_fwd) && (gs.opt & opts_bwd))
		 || ((opt & opts_bwd) && (gs.opt & opts_fwd))) {
			cnv |= conv_rev;
			fig += 1000;
#ifdef VERBOSE
			printf("  needs direction conversion though\n");
#endif
		}

		if (prec != gs.prec) { 	/* Internal precision mismatch */
			fig += 100000;		/* Will have a major effect, so discourage using it */
#ifdef VERBOSE
			printf("  internal precision doesn't match though (want %d, got %d)\n",prec,gs.prec);
#endif
		}

		/* If we've got a choice of algorithm possible */
		if (gs.opt & (opts_splx_sort | opts_sort_splx)) {

			/* If we've got a preference of algorithm, */
			if ((opt & (opts_splx_sort | opts_sort_splx))) {

				/* and there is a mismatch */
				if (((opt & opts_splx_sort) && ts.sort != 0)
				 || ((opt & opts_sort_splx) && ts.sort == 0)) {
					fig += 10000;		/* Will have a great effect, so discourage using it */
#ifdef VERBOSE
					printf("  sort/simplex algorithm mismatch\n");
#endif
				}

			} else {	/* There is no preference chosen at run time, */
				/* and there is a mismatch to the compile time preference */
				if (((gs.opt & opts_splx_sort) && ts.sort != 0)
				 || ((gs.opt & opts_sort_splx) && ts.sort == 0)) {
					fig += 10000;		/* Will have a great effect, so discourage using it */
#ifdef VERBOSE
					printf("  sort/simplex algorithm mismatch\n");
#endif
				}
			}
		}

		if (ts.sort) {
			stres = 0;
#ifdef VERBOSE
		printf("gres = %d, gs.gres = %d\n",res,gs.itres);
#endif
			/* We want one that is equals or exceeds the desired */
			/* resolution, but doesn't exceed it too much, or the code */
			/* will be inefficient. */
			/* If there are no routines that can meet the desired precision, */
			/* then it is ok to use the one closest to the desired precision. */
			if (gs.itres >= res) {
				fig += 10 * (gs.itres - res);
			} else {
				fig += 10000 + 10 * (res - gs.itres);
			}
		} else {
			/* compute the needed stres (Assuming not sort) */
			stres = ((1 << gs.prec)-1 + res-2)/(res-1);

#ifdef VERBOSE
		printf("gres = %d, sres = %d, gs.gres = %d, gs.sres = %d\n",res,stres,gs.itres,gs.stres);
#endif
			/* We want one that is equals or exceeds the desired */
			/* resolution, but doesn't exceed it too much, or the code */
			/* will be inefficient. */
			/* If there are no routines that can meet the desired precision, */
			/* then it is ok to use the one closest to the desired precision. */
			if (gs.itres >= res && gs.stres >= stres) {
				fig += 10 * (gs.itres - res) + (gs.stres - stres);
			} {
				if (gs.itres < res) {
					fig += 10000 + 10 * (res - gs.itres);
				}
				if (gs.stres < stres) {
					fig += 1000 + 10 * (stres - gs.stres);
				}
			}
		}

#ifdef VERBOSE
		printf("  figure of merit %d\n",fig);
#endif
		/* Is this the best one so far ? */
		if (fig < bfig) {
			bfig = fig;
			bk = i;
			bstres = stres;
			bgs = gs;		/* Structure copy */
			bts = ts;		/* Structure copy */
			bcnv = cnv;
#ifdef VERBOSE
			printf("  best so far\n");
#endif
		}
	}

	if (bk < 0) {
#ifdef VERBOSE
		printf("new_imdi failed - dimensionality or representations couldn't be matched\n");
#endif
		return NULL;	/* Nothing matches */
	}

	if ((im = (imdi *)calloc(1, sizeof(imdi))) == NULL) {
#ifdef VERBOSE
		printf("new_imdi malloc imdi failed\n");
#endif
		/* Should we return an error somehow ? */
		return NULL;
	}

	/* We've decided kernel function bk is going to be the best, */
	/* so now setup the appropriate tables to use with it. */
	if (bgs.itres > res)
		bgs.itres = res;		/* Tell table create what the res is */
	if (bgs.stres > bstres)
		bgs.stres = bstres;

	/* Tel table setup how to treat integer input in per channel lookups */
	bgs.in_signed = in_signed;
	bgs.out_signed = out_signed;

#ifdef VERBOSE
	if (!bts.sort) {
		if ((bgs.stres * (bgs.itres-1)) < ((1 << bgs.prec)-1)) {
			printf("Warning: table chosen doesn't reach desired precision!\n");
			printf("Wanted %d, got %d\n", ((1 << bgs.prec)-1), (bgs.stres * (bgs.itres-1)));
		}
	}
#endif

	/* Allocate and initialise the appropriate tables */
	im->impl = (void *)imdi_tab(&bgs, &bts, bcnv, in, out, ktable[bk].interp,
	                            inm, outm, oopt, checkv, input_curves, md_table,
	                            output_curves, cntx);

	if (im->impl == NULL) {
#ifdef VERBOSE
		printf("imdi_tab failed\n");
#endif
		imdi_del(im);
		return NULL;
	}

#ifdef VERBOSE
	if (bcnv != conv_none)
		printf("imdi_tab: using a runtime match, cnv flags 0x%x\n",bcnv);
#endif

	if (bcnv == conv_none)		/* No runtime match conversion needed */
		im->interp  = ktable[bk].interp;	
	else
		im->interp  = interp_match;
	im->get_check   = imdi_get_check;
	im->reset_check = imdi_reset_check;
	im->info        = imdi_info;
	im->del         = imdi_del;

#ifdef VERBOSE
	printf("new_imdi returning 0x%x\n", im);
#endif
	return im;
}

/* Runtime matching adapter */
/* We can emulate many combinations of interpolation function */
/* by converting to a routine that supports stride, or one that */
/* supports plane interleaved and stride. */
/* This also lets us emulate functions that can convert between */
/* pixel and plane interleaved at the same time as doing the */
/* color interpolation. Warp is in components (NOT bytes) */
/* We cope with skipping writes to output channels by */
/* only implementing that support in plane interleaved, */
/* and depending on that being used for pixel interleaved. */
static void interp_match(
imdi *s,
void **outp, int outst,		/* Output pointers and stride */
void **inp, int inst,		/* Input pointers and stride */
unsigned int npixels		/* Number of pixels */
) {
	int j, i;
	void *minp[IXDI];
	void *moutp[IXDI];
	imdi_imp *impl = (imdi_imp *)s->impl;
	imdi_conv cnv = impl->cnv;

	/* Need to supply default stride. */
	/* This is simply the input dimension for pixel interleaved, */
	/* and 1 for plane interleaved. */
	if (cnv & conv_istr) {
		if (impl->cirep == pixint8 || impl->cirep == pixint16)
			inst = impl->id; 
		else
			inst = 1;
	}
	if (cnv & conv_ostr) {
		if (impl->corep == pixint8 || impl->corep == pixint16)
			outst = impl->wod; 
		else
			outst = 1;
	}

	/* Convert from pixel to plane by setting the plane pointers */
	/* to each pixel component, and using the pixel interleaved stride */
	if (cnv & conv_irep) {
		if (impl->cirep == pixint8) {	/* Convert from pix8 to plane8 */
			for (j = 0; j < impl->id; j++)
				minp[j] = (void *)((char *)inp[0] + j);
		} else if (impl->cirep == pixint16) {	/* Convert from pix16 to plane16 */
			for (j = 0; j < impl->id; j++)
				minp[j] = (void *)((char *)inp[0] + 2 * j);
		}
	} else {	/* Copy pointers, because we call with them. */
		if (impl->cirep == pixint8 || impl->cirep == pixint16) {
			minp[0] = inp[0];
		} else {
			for (j = 0; j < impl->id; j++)
				minp[j] = inp[j];
		}
	}
	if (cnv & conv_orep) {
		if (impl->corep == pixint8) {	/* Convert from pix8 to plane8 */
			for (j = i = 0; j < impl->od; j++) {
				if ((impl->skipf & (1 << j)) == 0)	/* If not skipped */
					moutp[j] = (void *)((char *)outp[0] + i++);
				else
					moutp[j] = NULL;
			}
		} else if (impl->corep == pixint16) {	/* Convert from pix16 to plane16 */
			for (j = i = 0; j < impl->od; j++)
				if ((impl->skipf & (1 << j)) == 0)	/* If not skipped */
					moutp[j] = (void *)((char *)outp[0] + 2 * i++);
				else
					moutp[j] = NULL;
		}
	} else {	/* Copy pointers, because we call with them. */
		if (impl->corep == pixint8 || impl->corep == pixint16) {
			moutp[0] = outp[0];
		} else {	/* Plane interleaved */
			for (j = i = 0; j < impl->od; j++) {
				if ((impl->skipf & (1 << j)) == 0)	/* If not skipped */
					moutp[j] = outp[i++];
				else
					moutp[j] = NULL;
			}
		}
	}

	/* Convert fwd function to a backwards one, or visa-versa. */
	/* Move the pointers to the last pixel, and reverse the stride */
	if (cnv & conv_rev) {
		if (impl->firep == pixint8) {
			minp[0] = (void *)((char *)minp[0] + inst * (npixels - 1));
		} else if (impl->firep == pixint16) {
			minp[0] = (void *)((char *)minp[0] + inst * 2 * (npixels - 1));
		} else if (impl->firep == planeint8) {
			for (j = 0; j < impl->id; j++)
				minp[j] = (void *)((char *)minp[j] + inst * (npixels - 1));
		} else if (impl->firep == planeint16) {
			for (j = 0; j < impl->id; j++)
				minp[j] = (void *)((char *)minp[j] + inst * 2 * (npixels - 1));
		}
		inst = -inst;
		if (impl->forep == pixint8) {
			moutp[0] = (void *)((char *)moutp[0] + outst * (npixels - 1));
		} else if (impl->forep == pixint16) {
			moutp[0] = (void *)((char *)moutp[0] + outst * 2 * (npixels - 1));
		} else if (impl->forep == planeint8) {
			for (j = 0; j < impl->od; j++)
				moutp[j] = (void *)((char *)moutp[j] + outst * (npixels - 1));
		} else if (impl->forep == planeint16) {
			for (j = 0; j < impl->od; j++)
				moutp[j] = (void *)((char *)moutp[j] + outst * 2 * (npixels - 1));
		}
		outst = -outst;
	}

	/* Call the conversion routine we have */
	impl->interp(s, moutp, outst, minp, inst, npixels);
}

/* Get the per output channel check flags - bit corresponds to output interpolation channel */
static unsigned int imdi_get_check(imdi *im) {
	imdi_imp *impl = (imdi_imp *)im->impl;
	unsigned int rv;
	int i;

	/* Convert from output channel index to callback index */
	for (rv = 0, i = 0; i < impl->od; i++) {
		unsigned int b;
		b = 1 & (impl->checkf >> i);		/* Output index flag */
		rv |= (b << impl->im_map[i]);		/* to callback index */
	}
	return rv;
}

/* Reset the output check flags (not reset by interp) */
static void imdi_reset_check(imdi *im) {
	imdi_imp *impl = (imdi_imp *)im->impl;

	impl->checkf = 0;
}

/* Return some information about the imdi */
static void imdi_info(imdi *im, unsigned long *psize, int *pgres, int *psres) {
	imdi_imp *impl = (imdi_imp *)im->impl;
	unsigned long size;

	size = sizeof(imdi);

	if (impl != NULL) {
		size += impl->size;
		if (psize != NULL)
			*psize = size;

		if (pgres != NULL)
			*pgres = impl->gres;

		if (psres != NULL)
			*psres = impl->sres;
	}
}


/* Delete the object */
static void imdi_del(imdi *im) {
	/* Free all the allocated tables */
	if (im->impl != NULL)
		imdi_tab_free((imdi_imp *)im->impl);

	/* Free this structure */
	free(im);
}


















