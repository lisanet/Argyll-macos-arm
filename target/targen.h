
#ifndef TARGEN_H
/* 
 * Argyll Color Correction System
 * Test target chart Generator.
 *
 * Author: Graeme W. Gill
 * Date:   7/4/2002
 *
 * Copyright 2002, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/* Targen generation common defines */

#define MXTD ICX_MXINKS		/* From xicc/xcolorants.h */

/* A fixed point position */
struct _fxpos {
	double p[MXTD];		/* Device coordinate position */
	double v[MXTD];		/* Room for perceptual value */
	int eloc;			/* if >= 0, cgats index, to even location */
}; typedef struct _fxpos fxpos;


#define TARGEN_H
#endif /* TARGEN_H */
