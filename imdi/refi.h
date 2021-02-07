
/* Reference floating point interpolator constructed out of rspl's */
/* This provides imdi functionality in floating point */
/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include "../rspl/rspl.h"

/* ------------------------------------------------ */

typedef struct {
	int id, od;		/* Input and output dimensions */
	int inres;		/* Desired input table resolution */
	int clutres;	/* Desired clut table resolution */
	int outres;		/* Desired output table resolution */
	rspl *in[MXDI];
	rspl *clut;
	rspl *out[MXDO];

	void (*input_curves) (void *cntx, double *out_vals, double *in_vals);
	void (*md_table)     (void *cntx, double *out_vals, double *in_vals);
	void (*output_curves)(void *cntx, double *out_vals, double *in_vals);
	void *cntx;		/* Context to callbacks */
	int chan;		/* Current callback channel */
} refi;

refi *new_refi(
	int id,			/* Number of input dimensions */
	int od,			/* Number of output dimensions */
	int inres,		/* Desired input table resolution */
	int clutres,	/* Desired clut table resolution */
	int outres,		/* Desired output table resolution */

	/* Callbacks to lookup the table values */
	void (*input_curves) (void *cntx, double *out_vals, double *in_vals),
	void (*md_table)     (void *cntx, double *out_vals, double *in_vals),
	void (*output_curves)(void *cntx, double *out_vals, double *in_vals),
	void *cntx		/* Context to callbacks */
);

void refi_free(refi *r);


/* Component interpolations */
void refi_input(void *cntx, double *out_vals, double *in_vals);
void refi_clut(void *cntx, double *out_vals, double *in_vals);
void refi_output(void *cntx, double *out_vals, double *in_vals);

