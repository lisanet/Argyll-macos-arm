
/* Test the LU decomposition code */
/*
 * Copyright 1999 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include "numlib.h"

int test(int n, double **a, double *b);

int main() {
	double **a;
	double *b;
	int n = 4;
	int rv;

	a = dmatrix(0, n-1, 0, n-1);
	b = dvector(0, n-1);

	a[0][0] = 1.0;
	a[0][1] = 2.0;
	a[0][2] = 3.1415926;
	a[0][3] = 5.1415926;

	a[1][0] = 3.0;
	a[1][1] = 2.0;
	a[1][2] = 4.0;
	a[1][3] = -0.1;

	a[2][0] = -1.0;
	a[2][1] = 0.5;
	a[2][2] = 1.0;
	a[2][3] = 1.5;

	a[3][0] = 11.0;
	a[3][1] = 9.0;
	a[3][2] = 15.0;
	a[3][3] = 0.9;

	b[0] = 6.0;
	b[1] = 4.0;
	b[2] = 4.5;
	b[3] = -10.0;

	if ((rv = test(n, a, b)) != 0) {
		if (rv == 1)
			printf("LU test failed due to singularity\n");
		else {
			printf("LU test failed to verify\n");
			printf("Got solution %f %f %f %f\n",b[0],b[1],b[2],b[3]);
		}
	} else {
		printf("Got verified solution %f %f %f %f\n",b[0],b[1],b[2],b[3]);
	}
	return 0;
}


int
test(
int      n,	/* Dimensionality */
double **a,	/* A[][] input matrix, returns LU decimposition of A */
double  *b	/* B[]   input array, returns solution X[] */
) {
	int i, j;
	double rip;		/* Row interchange parity */
	int *pivx;
	int rv = 0;

	double **sa;		/* save input matrix values */
	double *sb;			/* save input vector values */

	pivx = ivector(0, n-1);
	sa = dmatrix(0, n-1, 0, n-1);
	sb = dvector(0, n-1);

	/* Copy input matrix and vector values */
	for (i = 0; i < n; i++) {
		sb[i] = b[i];
		for (j = 0; j < n; j++)
			sa[i][j] = a[i][j];
	}

	if (lu_decomp(a, n, pivx, &rip)) {
		free_dvector(sb, 0, n-1);
		free_dmatrix(sa, 0, n-1, 0, n-1);
		free_ivector(pivx, 0, n-1);
		return 1;
	}

	lu_backsub(a, n, pivx, b);

	/* Check that the solution is correct */
	for (i = 0; i < n; i++) {
		double sum, temp;
		sum = 0.0;
		for (j = 0; j < n; j++)
			sum += sa[i][j] * b[j];
//printf("~~ check %d = %f, against %f\n",i,sum,sb[i]);
		temp = fabs(sum - sb[i]);
		if (temp > 1e-6)
			rv = 2;
	}
	free_dvector(sb, 0, n-1);
	free_dmatrix(sa, 0, n-1, 0, n-1);
	free_ivector(pivx, 0, n-1);
	return rv;
}

