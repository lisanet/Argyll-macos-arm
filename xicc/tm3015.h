
#ifndef TM3015_H
#define TM3015_H

/* 
 * Author:  Graeme W. Gill
 * Date:    4/1/19
 * Version: 1.00
 *
 * Copyright 2019 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 2 or later :-
 * see the License2.txt file for licencing details.
 */

/*
 * These functions compute IES TM-30-15,
 * "IES Method for Evaluating Light Source Color Rendition"
 */

/*
 * TTBD:
 *
 */

#ifdef __cplusplus
	extern "C" {
#endif

#define IES_TM_30_15_BINS 16
#define IES_TM_30_15_ESAMPLES 99

/* Compute the EBU TLCI-2012 Rf & Rg */
/* Return 1 on invalid, 2 on error */
/* Invalid is when sample is not white enough. */
int icx_IES_TM_30_15(
double *pRf,			/* Return Rf */
double *pRg,			/* Return Rg */
double *pcct,			/* Return correlated color temperature */
double *pdc,			/* Return signed Du'v' to locus */
double pbins[IES_TM_30_15_BINS][2][3],		/* If not NULL, return ref & tsamp Jab */
xspect *tsamp			/* Illuminant sample to compute TLCI of */
);


#ifdef __cplusplus
	}
#endif

#endif /* TM3015_H */






































