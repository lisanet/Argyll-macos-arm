
#ifndef LUDECOMP_H
#define LUDECOMP_H

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

/* NOTE:- lu decomp rearanges the rows of the matrix */
/* by swapping pointers rather than exchanging data, */
/* so the matrix must be addressed by the **pointer */
/* if it is re-used after an ludecomp!!! */

/* Solve the simultaneous linear equations A.X = B */
/* Return 1 if the matrix is singular, 0 if OK */
int
solve_se(
double **a,	/* A[][] input matrix, returns LU decimposition of A */
double  *b,	/* B[]   input array, returns solution X[] */
int      n	/* Dimensionality */
);

/* Solve the simultaneous linear equations A.X = B, with polishing */
/* Return 1 if the matrix is singular, 0 if OK */
int
polished_solve_se(
double **a,	/* A[][] input matrix, returns LU decimposition of A */
double  *b,	/* B[]   input array, returns solution X[] */
int      n	/* Dimensionality */
);

/* Decompose the square matrix A[][] into lower and upper triangles */
/* NOTE that rows get swaped by swapping matrix pointers! */
/* Return 1 if the matrix is singular. */
int
lu_decomp(
double **a,		/* A input array, output upper and lower triangles. */
int      n,		/* Dimensionality */
int     *pivx,	/* Return pivoting row permutations record */
double  *rip	/* Row interchange parity, +/- 1.0, used for determinant */
);

/* Solve a set of simultaneous equations from the */
/* LU decomposition, by back substitution. */
void
lu_backsub(
double **a,		/* A[][] LU decomposed matrix */
int      n,		/* Dimensionality */
int     *pivx,	/* Pivoting row permutations record */
double  *b		/* Input B[] vecor, return X[] */
);

/* Polish a solution for equations */
void
lu_polish(
double **a,			/* Original A[][] matrix */
double **lua,		/* LU decomposition of A[][] */
int      n,			/* Dimensionality */
double  *b,			/* B[] vector of equation */
double  *x,			/* X[] solution to be polished */
int     *pivx		/* Pivoting row permutations record */
);

/* Invert a matrix A using lu decomposition */
/* NOTE that it returns transposed inverse by normal convention. */
/* Use sym_matrix_trans() to fix this, or use matrix_trans_mult() */ 
/* Return 1 if the matrix is singular, 0 if OK */
int
lu_invert(
double **a,	/* A[][] input matrix, returns inversion of A transposed */
int      n	/* Dimensionality */
);

/* Invert a matrix A using lu decomposition */
/* The normal convention inverse is returned */
/* Return 1 if the matrix is singular, 0 if OK */
int
lu_invert_normal(
double **a,	/* A[][] input matrix, returns inversion of A */
int      n	/* Dimensionality */
);

/* Invert a matrix A using lu decomposition, and polish it. */
/* NOTE that it returns transposed inverse by normal convention. */
/* Use sym_matrix_trans() to fix this, or use matrix_trans_mult() */
/* Return 1 if the matrix is singular, 0 if OK */
int
lu_polished_invert(
double **a,	/* A[][] input matrix, returns inversion of A */
int      n	/* Dimensionality */
);

/* Pseudo-Invert matrix A using lu decomposition */
/* Return 1 if the matrix is singular, 0 if OK */
int
lu_psinvert(
double **out,	/* Output[0..N-1][0..M-1] */
double **in,	/*  Input[0..M-1][0..N-1] input matrix */
int      m,		/* In Rows */
int      n		/* In Columns */
);

/* - - - - - - - - - - - - - - - - - - - */

/* Use Cholesky decomposition on a symetric positive-definite matrix. */
/* Only the upper triangle of the matrix A is accessed. */
/* L returns the decomposition */
/* Return nz if A is not positive-definite */
int llt_decomp(double **L, double **A, int n);


/* Solve a set of simultaneous equations A.x = b from the */
/* LLt decomposition, by back substitution. */
void llt_backsub(
double **L,			/* A[][] LLt decomposition in lower triangle */
int n,				/* Dimensionality */
double *b,			/* Input B[] */
double *x			/* Return X[] (may be same as B[]) */
);


#ifdef __cplusplus
	}
#endif

#endif /* LUDECOMP_H */
