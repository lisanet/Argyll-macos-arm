#ifndef SVD_H
#define SVD_H

/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifdef __cplusplus
	extern "C" {
#endif

/*

	U[] decomposes A[]'s columns into orthogonal, singular vectors.
	U[]'s columns are vectors that sum to 1.0, i.e. they leave a vectors normal unchanged.
	The inverse of U[] is its transpose.
	U[]'s columns corresponding to non-zero W[] are the orthonormal vectors that span
	the (input) RANGE space. Columns corresponding to zero W[] are zero.

	W[] will return singlular values, i.e. the weighting of the singular vectors.
	It's inverse is is the reciprical of its elements.

	V[] composes the singular vectors back into A[]'s rows.
	V[]'s columns and rows are orthonormal.
	V[]'s columns corresponding to non-zero W[] are the orthonormal vectors that span
	the (output) RANGE space. 
	V[]'s columns corresponding to zero W[] are the (output) othonormal vectors that span
	the NULL space.
	The inverse of V[] is its transpose.

	To re-form, A = U.W.Vt, i.e. multiply by transpose of V.

	i.e. mult. input vector[m] by U[] converts to [n] compact, orthogonal
	basis vectors. W then scales them appropiately, setting null space
	vectors to 0. V[] then transforms from the orthogonal basis vectors
	to A[]'s output space.

	To reveal NULL space, make sure n >= m, since U[] vectors corrsponding
	to zero's are set to zero.

	Inverse of A = V.1/W.Ut ?
*/

/* Compute Singular Value Decomposition of A = U.W.Vt */
/* Return status value: */
/* 0 = no error */
/* 1 - Too many itterations */
int svdecomp(
double **a,		/* A[0..m-1][0..n-1], return U[0..m-1][0..n-1] */
double  *w,		/* return W[0..n-1] = singular values */
double **v,		/* return V[0..n-1][0..n-1] (not transpose!) */
int      m,		/* Number of equations */
int      n		/* Number of unknowns */
);

/* Threshold the singular values W[] */ 
void svdthresh(
double w[],		/* Singular values */
int      n		/* Number of unknowns */
);

/* Threshold the singular values W[] to give a specific dof */ 
void svdsetthresh(
double w[],		/* Singular values */
int      n,		/* Number of unknowns */
int      dof	/* Expected degree of freedom */
);

/* Use output of svdcmp() to solve overspecified and/or */
/* singular equation A.x = b */
int svdbacksub(
double **u,		/* U[0..m-1][0..n-1] U, W, V SVD decomposition of A[][] */
double  *w,		/* W[0..n-1] */
double **v,		/* V[0..n-1][0..n-1] (not transpose!) */
double b[],		/* B[0..m-1]  Right hand side of equation */
double x[],		/* X[0..n-1]  Return solution. (May be the same as b[]) */
int      m,		/* Number of equations */
int      n		/* Number of unknowns */
);

/* Solve the equation A.x = b using SVD */
/* (The w[] values are thresholded for best accuracy) */
/* Return non-zero if no solution found */
/* !!! Note that A[][] will be changed !!! */
int svdsolve(
double **a,		/* A[0..m-1][0..n-1] input A[][], will return U[][] */
double b[],		/* B[0..m-1]  Right hand side of equation, return solution */
int      m,		/* Number of equations */
int      n		/* Number of unknowns */
);

/* Solve the equation A.x = b using SVD */
/* The top s out of n singular values will be used */
/* Return non-zero if no solution found */
/* !!! Note that A[][] will be changed !!! */
int svdsolve_s(
double **a,		/* A[0..m-1][0..n-1] input A[][], will return U[][] */
double b[],		/* B[0..m-1]  Right hand side of equation, return solution */
int      m,		/* Number of equations */
int      n,		/* Number of unknowns */
int      s		/* Number of unknowns */
);

/* Solve the equation A.x = b using Direct calculation, LU or SVD as appropriate */
/* Return non-zero if no solution found */
/* !!! Note that A[][] will be changed !!! */
int gen_solve_se(
double **a,		/* A[0..m-1][0..n-1] input A[][], will return U[][] */
double b[],		/* B[0..m-1]  Right hand side of equation, return solution */
int      m,		/* Number of equations */
int      n		/* Number of unknowns */
);

/* Compute the inverse matrix Ai[[0..n-1][0..m-1] from SVD components */
static void svdinverse(
double **u,		/* U[0..m-1][0..n-1] U, W, V SVD decomposition of A[][] */
double  *w,		/* W[0..n-1] */
double **v,		/* V[0..n-1][0..n-1] (not transpose!) */
double **ia,	/* iA[0..n-1][0..m-1] return inverse of A */
int      m,		/* Number of equations */
int      n		/* Number of unknowns */
);

#ifdef __cplusplus
	}
#endif

#endif /* SVD_H */
