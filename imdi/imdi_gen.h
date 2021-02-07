#ifndef IMDI_GEN_H
#define IMDI_GEN_H

/* Integer Multi-Dimensional Interpolation */

/*
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * Pixel code generation definitions.
 * 
 * This defines a particular combination of pixel layout,
 * number of channels, and other runtime requestable features.
 * This is used by the code generation to setup anticipated
 * support for these requests, used by the runtime to request
 * the features, and by the runtime to identify the correct
 * generated module to use.
 *
 */


/* -------------------------------------------------- */
/* High level kernel generation description */


/* This is a high level macro desciption of the pixel layout. */
/* It can be expanded by adding a new enumeration, and then */
/* implementing the code in imdi_gen to translate the enumeration */
/* into the exact pixlayout structure details. */

typedef enum {
	invalid_rep = 0x00,
	pixint8     = 0x01,		/* 8 Bits per value, pixel interleaved, no padding */
	planeint8   = 0x02,		/* 8 bits per value, plane interleaved */
	pixint16    = 0x03,		/* 16 Bits per value, pixel interleaved, no padding */
	planeint16  = 0x04		/* 16 bits per value, plane interleaved */
} imdi_pixrep;

/* The internal processing precision */
typedef enum {
	prec_min   = 0,		/* Minimum of input and output precision */
	prec_max   = 1,		/* Maximum of input and output precision */
	prec_in    = 2,		/* Same as input */
	prec_out   = 3,		/* Same as output */
	prec_p8    = 4,		/* 8 Bits precision */
	prec_p16   = 5		/* 16 Bits precision */
} imdi_iprec;

/* Output per channel processing options. For the compiled */
/* kernels, this is indexed by physical output channel, */
/* (not to be confused with the runtime oopt which is indexed */
/* by callback channel). */
/* - shift left by channel no * 2 to access option. */
/* Note that oopts_skip valid only for plane interleaved, */
/* and doesn't change the number output channels looked up, */
/* it changes the number of output channels written, hence */
/* must be matched by the arguments to interp() (number of */
/* output pointers, stride etc.) */ 
typedef enum {
	oopts_none     = 0x0,	/* No extra options */
	oopts_check    = 0x1,	/* Test channel value against trigger */
	oopts_skip     = 0x2,	/* Don't write channel value */
	oopts_chskp    = 0x3	/* oopts_check + oopts_skip */
} imdi_ooptions;

/* Macro to create an instance of imdi_ooptions for a particular channel */
#define OOPT(flag,chan) ((flag) << (chan * 2))

/* Macro to extract the flags particular channel */
#define OOPTX(val,chan) (0x3 & ((val) >> (chan * 2)))

/* Value of no oopt flags set */
#define OOPTS_NONE 0x00000000

/* Mask to detect any oopts_check flags that are set */
#define OOPTS_CHECK 0x55555555

/* Mask to ignore any oopts_check flags that are set */
#define OOPTS_NOT_CHECK 0xAAAAAAAA

/* Mask to detect any oopts_skip flags that are set */
#define OOPTS_SKIP 0xAAAAAAAA

/* Mask to ignore any oopts_skip flags that are set */
#define OOPTS_NOT_SKIP 0x55555555

/* Processing options */
/* Note that there will be a choice between simplex lookup or sort algoritm for */
/* only a subset of kernels (8 bits processing precision, <= 4 input channels) */
typedef enum {
	opts_none      = 0x00,	/* Forward direction, no stride */
	opts_fwd       = 0x01,	/* Forwards direction (default is gen fwd/run don't care) */
	opts_bwd       = 0x02,	/* Backwards direction (default is gen fwd/run don't care) */
	opts_istride   = 0x04,	/* Stride on input */
	opts_ostride   = 0x08,	/* Stride on output */

	opts_splx_sort = 0x10,	/* Prefer simplex over sort algorithm (generate both)  */
	opts_sort_splx = 0x20,	/* Force sort algorithm, rather than simplex table (generate both)  */
	opts_splx      = 0x40,	/* Generate simplex only (when possible), default is sort only. */

	opts_end       = 0x80	/* End marker */
} imdi_options;

/* This sructure allows a series of related kernels to be generated */
#define MX_GDCS 20
									  /* * means multiplies combination */
									  /* + means lockstep with previous line */
typedef struct {
	int idcombs[MX_GDCS];			  /* * Input dimension combinations (0 at end) */
	int itres[MX_GDCS];				  /* + Interpolation table resolutions for */
	int stres[MX_GDCS];				  /* + Simplex table resolutions */

	int odcombs[MX_GDCS];			  /* * Output dimensions combinations (0 at end) */
	imdi_ooptions ooptcombs[MX_GDCS]; /* + output per channel options */

	imdi_pixrep incombs[MX_GDCS];	  /* * Input pixel representation */
	imdi_iprec  iprecs[MX_GDCS];	  /* + Internal precision */
	imdi_pixrep outcombs[MX_GDCS];	  /* + Output pixel representation */

	imdi_options optcombs[MX_GDCS];	  /* * Direction and stride options */
} gendesc;

/* -------------------------------------------------- */
/* Detailed level of generation specification */

/* Pixel layout: */
/* Each channel value is assumed to be read with a 		              */
/* native machine single read of size bpch, 			              */
/* and bov and bpv being bit indexes into that value.	              */
/*                                        				              */
/* If pint == 0, then each read will be of size bpch[], and will      */
/* have its own pointer, will be incremented by chi[].	              */
/* If pint != 0, then the reads will be size bpch[0], from successive */
/* locations, chi[] apart. 								              */
/*                                        				              */
/* If packed == 0, then separate reads are needed for each input      */
/* channel.												              */
/* If packed != 0, then the channel values will be extracted          */
/* from a single read of size bpch[0]					              */
/*                                        				              */
/* Note that at all times the bit offset and size values              */
/* will be obeyed for each input value. 				              */

/* NOTE :- if you change this, you need to change the code in cgen.c */
/* labeled !genspec and tabspec delta code! */
typedef struct {
	int bpch[IXDIDO];	/* Bits per channel read (must be divisible by 8) */
	int chi[IXDIDO];	/* channel increment in multiples of bpch[] (0 == dimensionality) */
	int bov[IXDIDO];	/* Bit offset to value within channel */
	int	bpv[IXDIDO];	/* Bits per value within channel */
	int pint;			/* Flag - nonz if pixel interleaved (ie. reads from successice locations) */
	int packed;			/* Flag - nonz if all channels are packed into a single read */
} pixlayout;

/* Structure that specifies the configuration of a generated interpolation kernel. */
/* NOTE :- if you change this, you need to change the code in cgen.c */
/* labeled !genspec and tabspec delta code! */
typedef struct {
	/* Input to code generator */
	int prec;			/* Internal precision:- either 8 or 16 bits */ 
	int id;				/* Number of input dimensions used for interpolation */
	int od;				/* Number of output dimensions created by interpolation */
	imdi_pixrep irep;	/* High level input pixel representation */
	imdi_pixrep orep;	/* High level output pixel representation */
	int in_signed;		/* Bit flag per channel, NZ if treat as signed (runtime setup) */
	int out_signed;		/* Bit flag per channel, NZ if treat as signed (runtime setup) */
	pixlayout in;		/* Input pixel layout */
	pixlayout out;		/* Output pixel layout */
	imdi_ooptions oopt;	/* output per channel options */
	imdi_options opt;	/* Direction and stride options */
	int itres;			/* Interpolation table resolution */
	int stres;			/* Simplex table resolution */

	/* Returned value */
	char kkeys[100];		/* Kernel keys */
	char kdesc[100];		/* At genspec time */
	char kname[100];		/* At generation time */
} genspec;

/* - - - - - - - - - - - - - - - - - - - - - - - */

/* Translate between high level and low level generation specification */
int set_genspec(genspec *gs, gendesc *gd, int comb, mach_arch *a);

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Convert the input and output pixel representation and the internal */
/* precision choice into a precision in bits */
/* (Used at code gen and code execution, hence defined as a macro) */
#define COMPUTE_IPREC(PRECOUT, INREP, IPREC, OUTREP)		\
{																	\
	int _iprec, _oprec;												\
																	\
	switch (INREP) {												\
		default:													\
		case pixint8:												\
		case planeint8:												\
			_iprec = 8;												\
			break;													\
		case pixint16:												\
		case planeint16:											\
			_iprec = 16;											\
			break;													\
	}																\
	switch (OUTREP) {												\
		default:													\
		case pixint8:												\
		case planeint8:												\
			_oprec = 8;												\
			break;													\
		case pixint16:												\
		case planeint16:											\
			_oprec = 16;											\
			break;													\
	}																\
	switch (IPREC) {												\
		default:													\
		case prec_min:												\
			PRECOUT = _iprec < _oprec ? _iprec : _oprec;			\
			break;													\
		case prec_max:												\
			PRECOUT = _iprec > _oprec ? _iprec : _oprec;			\
			break;													\
		case prec_in:												\
			PRECOUT = _iprec;										\
			break;													\
		case prec_out:												\
			PRECOUT = _oprec;										\
			break;													\
		case prec_p8:												\
			PRECOUT = 8;											\
			break;													\
		case prec_p16:												\
			PRECOUT = 16;											\
			break;													\
	}																\
}																	\


/* - - - - - - - - - - - - - - - - - - - - - - - */

struct _tabspec;

/* Supported code generators: */

/* The 'C' code generator */
/* Generate a source file to implement the specified */
/* interpolation kernel. Fill in return values and return 0 if OK. */
/* Return 1 if this kernel should be skipped (ie. force sort and sort not forced) */
/* and some other anon-zero on error. */
int gen_c_kernel(genspec *g, struct _tabspec *t, mach_arch *a,
                 FILE *fp, int index, genspec *og, struct _tabspec *ot);

/* asm, MMX, etc. generators declarations go here ! */

#endif /* IMDI_GEN_H */











