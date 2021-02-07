#ifndef IMDI_H
#define IMDI_H

/* Integer Multi-Dimensional Interpolation */

/*
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This file provides the common definitions for IMDI, and
 * in particular, the runtime conversion object. 
 * Actual runtime details are kept opaque.
 */

#include "imdi_utl.h"
#include "imdi_arch.h"
#include "imdi_gen.h"

/* IMDI Object */
struct _imdi {
	void *impl;			/* Opaque pointer to implementation information (type imdi_imp *) */

	/* Do the interpolation. */

	/* Each pointer corresponds to the colors plane for plane interleaved. */
	/* pointer[0] is used for pixel interleave. */
	
	/* Stride is only obeyed if the appropriate option flag was set */
	/* in new_imdi, and is in pixel components, and effectively defaults */
	/* to 1 for plane interleaved, and id and od (adjusted for skip) for pixel interleave, */
	/* so stride is in color components (NOT bytes) */
	/* Output pointers and data must only reference non-skipped output channels. */

	/* Note that once an imdi is created, multiple can call interp() without */
	/* interfering with each other, allowing parallel execution. */
	void (*interp)(struct _imdi *s, void **outp, int outst,		/* Ouput pointers and stride */
	                                void **inp, int inst,		/* Input pointers and stride */
	                                unsigned int npixels);		/* Number of pixels */

	/* Return some information about the imdi */
	void (*info)(struct _imdi *s, unsigned long *size, int *gres, int *sres);

	/* Get the per output channel check flags (bit is indexed by callback channel) */
	/* Flag gets set if output != checkv */
	unsigned int (*get_check)(struct _imdi *s);
	
	/* Reset the output check flags (flag is not reset by interp) */
	void (*reset_check)(struct _imdi *s);

	/* Delete this object */
	void (*del)(struct _imdi *s);

}; typedef struct _imdi imdi;

/* Create a new imdi. */
/* Return NULL if request is not supported */
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
);

#endif /* IMDI_H */






