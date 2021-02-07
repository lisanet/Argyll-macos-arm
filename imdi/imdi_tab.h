#ifndef IMDI_TAB_H
#define IMDI_TAB_H

/* Integer Multi-Dimensional Interpolation */

/*
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* 
 * Implementation details needed for table initialisation for a particular
 * kernel. This is private implementation for imdi.[ch]
 *
 * The tabspec structure holds detailed information on the algorithms used
 * by the runtime code, and (implicit in this) the layout of the runtime
 * tables needed to match the algorithm. There are also implicit dependencies on
 * the genspec structure, since this determines the overall features
 * supported by a particular pixel kernel module.
 *
 * This is effectively the product of the genspec, the architechure,
 * and the coding choices made by the code generator
 * (ie. gen_c_kernel() in cgen.c)
 *
 */

/* entries marked with '#' are not currently used by imdi_tab() */
/* NOTE :- if you change this, you need to change the code in cgen.c */
/* labeled !genspec and tabspec delta code! */
struct _tabspec {

	int sort;	/* NZ for explicit sort rather than simplex table lookup */
	int it_xs;	/* NZ if separate interp index and simplex index/Weighting+Offset values  */
	int wo_xs;	/* NZ if separate weighting and vertex offset entries are to be used */

	int it_ix;	/* Non-zero if input value extraction should be done in input table */
	int it_ab;	/* Input table entry size in bits */
	int it_ts;	/* Input table :- total input table entry size in bytes */
				/* Bit packing order is (ms to ls) :
	          	    sort: ix, we, vo
	          	    sort: ix, wo
	          	   !sort: ix, sx
				 */

				/* Interpolation index is always in the input table */
	int ix_ab;	/* # Interpolation index entry size in bits */
	int ix_es;	/* Interpolation index entry size in bytes */
	int ix_eo;	/* Interpolation index entry offset in bytes */

				/* Simplex Index is always in the input table */
	int sx_ab;	/* Simplex Index entry size in bits */
	int sx_es;	/* Simplex Index entry size in bytes */
	int sx_eo;	/* Simplex Index entry offset in bytes */

	int sm_ts;	/* Simplex table entry total size in bytes */
				/* Bit packing order is (ms to ls) : we, vo */

				/* Combined Weighting + Offset may be in input table or Simplex entry */
	int wo_ab;	/* Combined Weighting + Offset entry size in bits */
	int wo_es;	/* Combined Weighting + Offset entry size in bytes */
	int wo_eo;	/* Combined Weighting + Offset entry offset in bytes */

				/* Weighting may be in input table or Simplex entry */
	int we_ab;	/* # Weighting entry size in bits */
	int we_es;	/* Weighting entry size in bytes */
	int we_eo;	/* Weighting entry offset in bytes */

				/* Vertex offset may be in input table or Simplex entry */
	int vo_ab;	/* Vertex Offset entry size in bits */
	int vo_es;	/* Vertex Offset entry size in bytes */
	int vo_eo;	/* Vertex Offset entry offset in bytes */
	int vo_om;	/* Vertex Offset scaling multiplier */

	int im_cd;	/* Non-zero if interpolation table entries are padded with fraction */
	int im_ts;	/* Interp. multidim :- total interp table entry size in bytes */
	int im_oc;	/* # Interp. multidim :- offset scale to apply to index into interp entry */
	int im_fs;	/* Interp. multidim :- full table entry size in bytes */
	int im_fn;	/* Interp. multidim :- number of full entries */
	int im_fv;	/* Interp. multidim :- output values per full entry . */
	int im_ps;	/* Interp. multidim :- partial table entry size in bytes, used & unsused */
	int im_pn;	/* Interp. multidim :- number of partial entries - must be 0 or 1 */
	int im_pv;	/* Interp. multidim :- used output values per partial entry . */

	int ot_ts;	/* Output table :- total entry size in bytes of every table */
	int ot_off[IXDO];	/* Offset for each output value within the output word needed */
	int ot_bits[IXDO];	/* Number of bits for value within the output word needed */

	/* Associated interpolation function */
	void (*interp)(struct _imdi *s, void **inp, void **outp, unsigned int npix); /* At run time */
}; typedef struct _tabspec tabspec;

/* Runtime conversion needed */
typedef enum {
	conv_none  = 0x00,	/* No conversion needed */
	conv_istr  = 0x01,	/* Input stride conversion */
	conv_ostr  = 0x02,	/* Output stride conversion */
	conv_irep  = 0x04,	/* Input representation conversion */
	conv_orep  = 0x08,	/* Output representation conversion */
	conv_rev   = 0x10,	/* Reverse direction conversion */
	conv_skip  = 0x20	/* Skip output channel write conversion */
} imdi_conv;

/* The actual run time table that tabspec describes */
typedef struct {
	/* Runtime setup */
	int id;						/* Number of input dimensions */
	int od;						/* Number of output dimensions (including skip channels) */
	int wod;					/* Number of written output dimensions ( < od if skipf != 0) */
	int it_map[IXDI];			/* Mapping from input raster channels to callback channels. */
	int im_map[IXDO];			/* Mapping from output raster channels to callback channels. */
	imdi_pixrep cirep;			/* High level input pixel representation called with  */
	imdi_pixrep corep;			/* High level output pixel representation called with */
	imdi_pixrep firep;			/* High level input pixel representation of interp func. */
	imdi_pixrep forep;			/* High level output pixel representation of interp func. */
	imdi_conv cnv;				/* Runtime argument conversion needed */
	void (*interp)(struct _imdi *s, void **outp, int outst,	/* Underlying conversion function */
	                                void **inp, int inst,
	                                unsigned int npixels);
	/* Output channel check data */
	unsigned long checkv[IXDO];	/* Output per channel check values. Set flag if != checkv */
	unsigned int checkf;		/* Output per channel check flags (one per bit) */
	unsigned int skipf;			/* Output per channel skip flags (one per bit) */

	/* Table data */
	void *in_tables[IXDI];		/* Input dimension input lookup tables */
    void *sw_table;				/* Simplex weighting lookup table */
    void *im_table;				/* Interpolation Multi-dimensional lookup table */
	void *out_tables[IXDO];		/* Output dimension output lookup tables */
	int nintabs;				/* Number of input tables */
	int nouttabs;				/* Number of output tables */

	/* Extra reporting data */
	unsigned long size;			/* Number of bytes allocated to imdi_imp */
	unsigned int gres, sres;	/* Grid and simplex table resolutions. sres = 0 = sort */
} imdi_imp;

/*
 * The runtime function that knows how to setup an imdi_imp
 * table for for our chosen kernel and the color mapping we
 * want to perform.
 */

imdi_imp *
imdi_tab(
	genspec *gs,		/* Pointer to gen spec */
	tabspec *ts,		/* Pointer to table spec */
	imdi_conv cnv,		/* Runtime argument conversion needed */
	imdi_pixrep irep,	/* High level input pixel representation to match  */
	imdi_pixrep orep,	/* High level output pixel representation to match */
	void (*interp)(struct _imdi *s, void **outp, int outst,	/* Underlying conversion function */
	                                void **inp, int inst,
	                                unsigned int npixels),
	int *inm,			/* Input raster channel to callback channel mapping, NULL for none. */
	int *outm,			/* Output raster channel to callback channel mapping, NULL for none. */
	imdi_ooptions oopt,	  /* Output per channel options (Callback channel, NOT written channel) */
	unsigned int *checkv, /* Output channel check values (Callback channel, NULL for none == 0. */

	/* Callbacks to initialse the imdi table values */
	void (*input_curves) (void *cntx, double *out_vals, double *in_vals),
	void (*md_table)     (void *cntx, double *out_vals, double *in_vals),
	void (*output_curves)(void *cntx, double *out_vals, double *in_vals),
	void *cntx		/* Context of callbacks */
);

void imdi_tab_free(imdi_imp *it);

#endif /* IMDI_TAB_H */
