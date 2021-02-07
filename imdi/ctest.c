
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

/* Test frame for cgen.c IMDI generation code */
/*
 * Copyright 2000 - 2006 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include "imdi.h"
#include "imdi_tab.h"


int
main(void) {
	int rv;
	genspec gs, ogs;
	tabspec ts, ots;
	mach_arch  ar;
	FILE *kcode;

	/* Zero out the gen and tabspecs, to give diff a place to start */
	memset((void *)&ogs, 0, sizeof(genspec));
	memset((void *)&gs, 0, sizeof(genspec));
	memset((void *)&ots, 0, sizeof(tabspec));
	memset((void *)&ts, 0, sizeof(tabspec));

	printf("Testing gen_c_kernel\n");

	if ((kcode = fopen("imdi_k99.c","w")) == NULL) {
		printf("Couldn't open 'imdi_k99.c'\n");
		return -1;
	}

	/* Setup an interpolation kernel specification */
	gs.id = 4;			/* Number of input dimensions */

	gs.in.pint = 1;		/* Flag - nonz if pixel interleaved */
	gs.in.packed = 0;	/* Flag - nonz if channels packed into one read */

	gs.in.bpch[0] = 8;	/* Bits per channel */
	gs.in.chi[0] = 4;	/* Channel increment */
	gs.in.bov[0] = 0;	/* Bit offset to value within channel */
	gs.in.bpv[0] = 8;	/* Bits per value within channel */

	gs.in.bpch[1] = 8;	/* Bits per channel */
	gs.in.chi[1] = 4;	/* Channel increment */
	gs.in.bov[1] = 0;	/* Bit offset to value within channel */
	gs.in.bpv[1] = 8;	/* Bits per value within channel */

	gs.in.bpch[2] = 8;	/* Bits per channel */
	gs.in.chi[2] = 4;	/* Channel increment */
	gs.in.bov[2] = 0;	/* Bit offset to value within channel */
	gs.in.bpv[2] = 8;	/* Bits per value within channel */

	gs.in.bpch[3] = 8;	/* Bits per channel */
	gs.in.chi[3] = 4;	/* Channel increment */
	gs.in.bov[3] = 0;	/* Bit offset to value within channel */
	gs.in.bpv[3] = 8;	/* Bits per value within channel */

	gs.od = 1;			/* Number of output dimensions */

	gs.out.pint = 1;	/* Flag - nonz if pixel interleaved */
	gs.out.packed = 0;	/* Flag - nonz if channels packed into one write */

	gs.out.bpch[0] = 8;	/* Bits per channel */
	gs.out.chi[0] = 1;	/* Channel increment */
	gs.out.bov[0] = 0;	/* Bit offset to value within channel */
	gs.out.bpv[0] = 8;	/* Bits per value within channel */

	gs.out.bpch[1] = 8;	/* Bits per channel */
	gs.out.chi[1] = 4;	/* Channel increment */
	gs.out.bov[1] = 0;	/* Bit offset to value within channel */
	gs.out.bpv[1] = 8;	/* Bits per value within channel */

	gs.out.bpch[2] = 8;	/* Bits per channel */
	gs.out.chi[2] = 4;	/* Channel increment */
	gs.out.bov[2] = 0;	/* Bit offset to value within channel */
	gs.out.bpv[2] = 8;	/* Bits per value within channel */

	gs.out.bpch[3] = 8;	/* Bits per channel */
	gs.out.chi[3] = 4;	/* Channel increment */
	gs.out.bov[3] = 0;	/* Bit offset to value within channel */
	gs.out.bpv[3] = 8;	/* Bits per value within channel */

	gs.out.bpch[4] = 8;	/* Bits per channel */
	gs.out.chi[4] = 4;	/* Channel increment */
	gs.out.bov[4] = 0;	/* Bit offset to value within channel */
	gs.out.bpv[4] = 8;	/* Bits per value within channel */

	gs.out.bpch[5] = 8;	/* Bits per channel */
	gs.out.chi[5] = 4;	/* Channel increment */
	gs.out.bov[5] = 0;	/* Bit offset to value within channel */
	gs.out.bpv[5] = 8;	/* Bits per value within channel */

	gs.opt = opts_none;		/* Direction and stride options */

	gs.prec  = 8;		/* Precsision needed */
	gs.itres = 16;		/* Interpolation table resolution */
	gs.stres = 17;		/* Simplex table resolution */

	/* Setup a machine architecture */

	ar.bigend = 0;		/* Non-zero if this is a bigendian architecture */
	ar.uwa = 0;			/* Use wide memory access */
	ar.shfm = 0;		/* Use shifts to mask values */
	ar.oscale = 8;		/* Has scaled indexing up to * 8 */
	ar.smmul = 0;		/* Doesn't have fast small multiply for index scaling */

	ar.pbits = 32;		/* Number of bits in a pointer */

	ar.nords = 3;		/* Number of ord types */
	ar.ords[0].bits = 8;
	ar.ords[0].name = "unsigned char";
	ar.ords[0].align = 1;
	ar.ords[1].bits = 16;
	ar.ords[1].name = "unsigned short";
	ar.ords[1].align = 1;
	ar.ords[2].bits = 32;
	ar.ords[2].name = "unsigned int";
	ar.ords[2].align = 1;
#ifdef ALLOW64
	ar.ords[3].bits = 64;
	ar.ords[3].name = "unsigned longlong";
	ar.ords[3].align = 0;
#endif /* ALLOW64 */
	ar.natord = 2;		/* Most natural type */

	ar.nints = 3;		/* Number of int types */
	ar.ints[0].bits = 8;
	ar.ints[0].name = "signed char";
	ar.ints[0].align = 1;
	ar.ints[1].bits = 16;
	ar.ints[1].name = "short";
	ar.ints[1].align = 1;
	ar.ints[2].bits = 32;
	ar.ints[2].name = "int";
	ar.ints[2].align = 1;
#ifdef ALLOW64
	ar.ints[3].bits = 64;
	ar.ints[3].name = "longlong";
	ar.ints[3].align = 0;
#endif /* ALLOW64 */
	ar.natint = 2;		/* Most natural type */

	rv = gen_c_kernel(&gs, &ts, &ar, kcode, 99, &ogs, &ots);

	fclose(kcode);

	return 0;
}
