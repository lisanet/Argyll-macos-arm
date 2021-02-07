
#ifndef NUMLIB_H
#define NUMLIB_H

/* Numerical routine library interface declarations */

#include "numsup.h"		/* Support routines, macros */
#include "dnsq.h"		/* Non-linear equation solver */
#include "gnewt.h"		/* Global Newton non-linear equation solver */
#include "powell.h"		/* Powell multi dimentional minimiser */
#include "dhsx.h"		/* Downhill simplex multi dimentional minimiser */
#include "varmet.h"		/* Variable Metric multi dimentional minimiser */
#include "ludecomp.h"	/* LU decomposition matrix solver */
#include "svd.h"		/* Singular Value decomposition matrix solver */
#include "zbrent.h"		/* 1 dimentional brent root search */
#include "rand.h"		/* Random number generators */
#include "sobol.h"		/* Sub-random vector generators */
#include "aatree.h"		/* Anderson balanced binary tree */
#include "quadprog.h"	/* Quadradic Programming solution */
#include "roots.h"		/* Quadratic, Cubic and Quartic root solving */

#endif /* NUMLIB_H */
