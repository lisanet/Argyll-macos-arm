
/* Integer Multi-Dimensional Interpolation */

/*
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * Code generation support that translates between the higher level
 * gendesc, and the low level genspec.
 */

#undef VERBOSE

/*
	Ideal grid resolutions for 8/16 bit precision calculations:

	There are a limited number of grid resolution and grid cells
	sub resolution values that work perfectly, and allow every possible
	input value to map to a precise interpolation value. For any
	system that deals with Lab data, it is also better if there is
	a grid point near a,b == 0, meaning that an odd grid size is
	preferable for 3 dimensional Lab input conversions. Reasonable memory
	usage is usually somewhere between 100K and 1M entries for
	the grid table, limiting the maximum grid resolution.
	The following lists the some of the possible grid and sub grid
	resolutions in order of best to worst. (Fewer jumps is better).


	Grid	Sub8	Round8	Sub16	Round16	Jumps	
	4		85		85		21845	21845	0
	6		51		51		13107	13107	0
	16		17		17		4369	4369	0
	18		15		15		3855	3855	0
	52		5		5		1285	1285	0
	86		3		3		771		771		0
	256		1		1		257		257		0
	258						255		255		0
	772						85		85		0
	1286					51		51		0
	3856					17		17		0
	4370					15		15		0
	13108					5		5		0
	21846					3		3		0
	65536					1		1		0
	3		127.5	128						1
	5		63.75	64						1
	9		31.875	32						1
	17		15.937	16						1
	33		7.968	8						1
	65		3.984	4						1
	128		2.007	2						1
	129		1.992	2						1
	255		1.003	1						1
	12		23.188	23						2
	24		11.086	11						2
	254		1.007	1						2
	7		42.5	43						3
	8		36.428	36						3
	10		28.333	28						3
	13		21.25	21						3
	15		18.214	18						3
	19		14.166	14						3
	22		12.142	12						3
	29		9.107	9						3
	37		7.083	7						3
	43		6.071	6						3
	44		5.930	6						3
	64		4.047	4						3
	85		3.035	3						3
	87		2.965	3						3
	127		2.023	2						3
	130		1.976	2						3
	253		1.011	1						3

	[8  bit: sub = 255/(grid-1), jumps = abs((grid-1)*round(sub)-255)]
	[16 bit: sub = 65535/(grid-1), jumps = abs((grid-1)*round(sub)-65535)]

	The above takes into consideration the mapping of the sub-cell or
	simplex table resolution, but doesn't take into account the quantizing
	of the sub-cell weighting values into the range 0..256 or 0..65536.

	This will be best when round(sub) divides evenly into 256 or 65536,
	so grid resolutions of 3, 5, 9, 17, 33, 65, 128, 129, 255 may be the
	best choice for sort algorithm grid resolutions.

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "copyright.h"
#include "aconfig.h"

#include "imdi_utl.h"
#include "imdi_arch.h"
#include "imdi_gen.h"

static void translate_pixrep(pixlayout *pl, char **desc, int *prec, imdi_pixrep rep, int dim, mach_arch *a);

/* Translate between a gendesc and genspec */
int 				/* Return number of combinations possible */
set_genspec(
genspec *gs,		/* Return specification */
gendesc *gd,		/* Input description */
int comb,			/* Combination index */
mach_arch *a		/* Machine architecture */
) {
	int nc = 1;			/* Number of combinations */
	int nidc, id;		/* Number of input dimension combinations, current index */
	int nodc, od;		/* Number of output dimension combinations, current index */
	int nirc, ir;		/* Number of input, prec. & out representations combinations, cnt index */
	int pc, or;			/* Precision combination index, Output combination index */
	int ndc, opt;		/* Number of direction stride ombinations, current index */
	int tt;

	/* Figure out how many combinations there are */
	for (nidc = 0; gd->idcombs[nidc] != 0; nidc++)	/* Input dimension */
		;
	nc *= nidc;
	for (nodc = 0; gd->odcombs[nodc] != 0; nodc++)	/* Output dimension */
		;
	nc *= nodc;
	for (nirc = 0; gd->incombs[nirc] != 0; nirc++)	/* Input/Precision/Output representation */
		;
	nc *= nirc;
	for (tt = 0; gd->iprecs[tt] != 0; tt++)	/* Internal precision representation */
		;
	if (nirc != tt) {
		fprintf(stderr,"imdi_gen: Must be equal numberof input and precision representations\n");
		exit(-1);
	}
	for (tt = 0; gd->outcombs[tt] != 0; tt++)	/* Output representation */
		;
	if (nirc != tt) {
		fprintf(stderr,"imdi_gen: Must be equal numberof input and output representations\n");
		exit(-1);
	}
		
	for (ndc = 0; gd->optcombs[ndc] != opts_end; ndc++)	/* Direction and stride options */
		;
	nc *= ndc;

	if (nc <= 0) {
		fprintf(stderr,"imdi_gen: no valid combinations specified for kernel generation\n");
		exit(-1);
	}
	if (comb < nc) {	/* If we are within a legal combination */
		int iprec, oprec;
		char *idesc, *odesc;	/* Representation descriptions */
		char *ddesc;			/* Direction description */

		id = comb % nidc;
		comb /= nidc;
		od = comb % nodc;
		comb /= nodc;
		ir = comb % nirc;
		comb /= nirc;
		pc = ir;			/* In and precision combs track together */
		or = ir;			/* In and out combs track together */
		opt = comb % ndc;
		comb /= ndc;
#ifdef VERBOSE
		printf("Combination id = %d, od = %d, ir = %d, pc = %d, or = %d, opt = %d\n",
		       id,od,ir,pc,or,opt);
#endif /* VERBOSE */

		gs->id = gd->idcombs[id];		/* Input dimensions */
		gs->itres = gd->itres[id];		/* Interpolation table resolution */
		gs->stres = gd->stres[id];		/* Simplex table resolution */
		gs->od = gd->odcombs[od];		/* Output dimensions */
		gs->oopt = gd->ooptcombs[od];	/* Output per channel options */

		if (gs->id > IXDI) {
			fprintf(stderr,"imdi_gen: Input dimension %d exceeds limit %d\n",gs->id,IXDI);
			exit(-1);
		}
		if (gs->od > IXDO) {
			fprintf(stderr,"imdi_gen: Output dimension %d exceeds limit %d\n",gs->od,IXDO);
			exit(-1);
		}
		/* Input representation */
		gs->irep = gd->incombs[ir];		/* Keep a copy of this */
		translate_pixrep(&gs->in,  &idesc, &iprec, gd->incombs[ir], gs->id, a);
		gs->in_signed = 0x0;			/* Not used during generation, used at runtime setup */

		/* Output representation */
		gs->orep = gd->outcombs[or];		/* Keep a copy of this */
		translate_pixrep(&gs->out, &odesc, &oprec, gd->outcombs[or], gs->od, a);
		gs->out_signed = 0x0;				/* Not used during generation, used at runtime setup */

		COMPUTE_IPREC(gs->prec, gs->irep, gd->iprecs[pc], gs->orep)

		gs->opt = gd->optcombs[opt];	/* Direction and stride options */

		if (((gs->opt & opts_splx_sort) && (gs->opt & opts_sort_splx))
		 || ((gs->opt & opts_splx_sort) && (gs->opt & opts_splx))
		 || ((gs->opt & opts_sort_splx) && (gs->opt & opts_splx))) {
			fprintf(stderr,"imdi_gen: Conflict in simplex/sort preferences\n");
			exit(-1);
		}

		ddesc = gs->opt & opts_bwd ? "b" : "f";	/* Direction description */

#ifdef VERBOSE
		printf("translates to prec = %d, id = %d, od = %d, itres %d, stdres %d\n",
		       gs->prec,gs->id,gs->od,gs->itres,gs->stres);
#endif /* VERBOSE */

		/* Create a concise description string */
		sprintf(gs->kdesc,"%d_%d_%s_%s_%s", gs->id, gs->od, idesc, odesc, ddesc);
	}
	return nc;
}

/* Convert the high level pixrep into the lower level pixel layout */
static void
translate_pixrep(
pixlayout *pl,		/* pixlayout to initialise */
char **desc,		/* Return description identifier (may be NULL) */
int *prec,			/* Return basic precision specifier (may be NULL) */
imdi_pixrep rep,	/* Representation to be translated */
int dim,			/* Number of dimensions (values/pixel) */
mach_arch *a		/* Machine architecture */
) {
	switch (rep) {

		case pixint8: {		/* 8 Bits per value, pixel interleaved, no padding */
			int i;

			/* Could optimise this to packed for dim == 4 ~~~~ */

			pl->pint = 1;		/* pixel interleaved */
			pl->packed = 0;		/* Not packed */
	
			for (i = 0; i < dim; i++) {
				pl->bpch[i] = 8;	/* Bits per channel */
				pl->chi[i] = dim;	/* Channel increment */
				pl->bov[i] = 0;		/* Bit offset to value within channel */
				pl->bpv[i] = 8;		/* Bits per value within channel */
			}
			
			if (prec != NULL)
				*prec = 8;			/* Basic 8 bit precision */
			if (desc != NULL)
				*desc = "i8";		/* Interleaved 8 bit */
		} break;
	
		case planeint8: {		/* 8 bits per value, plane interleaved */
			int i;

			pl->pint = 0;		/* Not pixel interleaved */
			pl->packed = 0;		/* Not packed */
	
			for (i = 0; i < dim; i++) {
				pl->bpch[i] = 8;	/* Bits per channel */
				pl->chi[i] = 1;		/* Channel increment */
				pl->bov[i] = 0;		/* Bit offset to value within channel */
				pl->bpv[i] = 8;		/* Bits per value within channel */
			}
			
			if (prec != NULL)
				*prec = 8;			/* Basic 8 bit precision */
			if (desc != NULL)
				*desc = "p8";		/* Planar 8 bit */

		} break;

		case pixint16: {		/* 16 Bits per value, pixel interleaved, no padding */
			int i;

			/* Could optimise this to packed for dim == 4 ~~~~ */

			pl->pint = 1;		/* pixel interleaved */
			pl->packed = 0;		/* Not packed */
	
			for (i = 0; i < dim; i++) {
				pl->bpch[i] = 16;	/* Bits per channel */
				pl->chi[i] = dim;	/* Channel increment */
				pl->bov[i] = 0;		/* Bit offset to value within channel */
				pl->bpv[i] = 16;	/* Bits per value within channel */
			}
			
			if (prec != NULL)
				*prec = 16;			/* Basic 8 bit precision */
			if (desc != NULL)
				*desc = "i16";		/* Interleaved 16 bit */
		} break;
	
		case planeint16: {		/* 16 bits per value, plane interleaved */
			int i;

			pl->pint = 0;		/* Not pixel interleaved */
			pl->packed = 0;		/* Not packed */
	
			for (i = 0; i < dim; i++) {
				pl->bpch[i] = 16;	/* Bits per channel */
				pl->chi[i] = 1;		/* Channel increment */
				pl->bov[i] = 0;		/* Bit offset to value within channel */
				pl->bpv[i] = 16;	/* Bits per value within channel */
			}
			
			if (prec != NULL)
				*prec = 16;			/* Basic 8 bit precision */

			if (desc != NULL)
				*desc = "p16";		/* Planar 8 bit */
		} break;

		default: {
			fprintf(stderr,"Warning: Unknown pixel representation %d\n",rep);
		} break;
	}
}

