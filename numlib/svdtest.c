
/* SVD test */
/* Verify that the SVD solver does what it is supposed to. */
/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* We assume two device dependent variables, and one objective function value */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#if defined(__IBMC__)
#include <float.h>
#endif

#include "numlib.h"

int main(void) {
	int its;
	int i,j,x;
	double **a;		/* A[0..M-1][0..N-1] input */
	double **u;		/* U[0..M-1][0..N-1] output */
	double *w;		/* W[0..N-1]	     output */
	double **v;		/* V[0..N-1][0..N-1] output */
	double **a2;	/* A[0..M-1][0..N-1] check on a */
	double **t;		/* A[0..M-1][0..N-1] temp */
	int m,n;

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	printf("Test SVD\n");

	for (its = 400; its > 0; its--) {
		int bad;
		int bad0;

		m  = i_rand(1,30);		/* Number of equations */
		n  = i_rand(1,30);		/* Number of unknowns */

		a  = dmatrix(0,m-1, 0,n-1);
		t  = dmatrix(0,m-1, 0,n-1);
		a2 = dmatrix(0,m-1, 0,n-1);
		u  = dmatrix(0,m-1, 0,n-1);
		w  = dvector(0,n-1);
		v  = dmatrix(0,n-1, 0,n-1);

		printf("Testing %d by %d\n",m,n);

		/* Create A matrix */
		for (j = 0; j < m; j++)
			for (i = 0; i < n; i++)
				a[j][i] = d_rand(-10.0, 10.0); 

		/* Setup u */
		for (j = 0; j < m; j++)
			for (i = 0; i < n; i++)
				u[j][i] = a[j][i];

		/* decompose A into U, W and V */
		svdecomp(u, w, v, m, n);

		/* Check results by computing a2 = U.W.Vt */
		for (j = 0; j < m; j++) {	/* U.W */
			for (i = 0; i < n; i++) {
				t[j][i] = 0.0;
				for (x = 0; x < n; x++) {
					if (x == i)
						t[j][i] += u[j][x] * w[x];
				}
			}
		}
		for (j = 0; j < m; j++) {	/* .Vt */
			for (i = 0; i < n; i++) {
				a2[j][i] = 0.0;
				for (x = 0; x < n; x++) {
					a2[j][i] += t[j][x] * v[i][x];
				}
			}
		}
		
		/* Now check */
		bad = 0;
		for (j = 0; j < m; j++)
			for (i = 0; i < n; i++) {
				double tt;
				tt = a2[j][i] - a[j][i];
				tt = fabs(tt);
				if (tt > 0.0000001) {
					bad = 1;
				}
		}
		if (bad)
			printf("A == U.W.Vt Check failed!\n");


		/* Check that U and Ut are inverses */
		bad = bad0 = 0;
		for (j = 0; j < n; j++) {
			for (i = 0; i < n; i++) {
				double t2, tt = 0.0;
				for (x = 0; x < m; x++) {
					tt += u[x][j] * u[x][i];
				}
				t2 = tt;
				if (i == j)
					tt -= 1.0;
				tt = fabs(tt);
				if (tt > 0.0000001) {
					if (i == j && fabs(t2) < 0.0000001)
						bad0++;				/* Unexpected zero diagonal */
					else {
						bad = 1;
						printf("Possible U error at %d %d = %f \n",j,i,tt);
					}
				}
			}
		}
		/* Expect n-m diagnals to be 0 instead of 1 if m < n */
		if (bad || (m >= n && bad0) || (m < n && bad0 != n-m))
			printf("U,Ut == 1 Check failed!\n");

		/* Check that V and Vt are inverses */
		bad = 0;
		for (j = 0; j < n; j++) {
			for (i = 0; i < n; i++) {
				double tt = 0.0;
				for (x = 0; x < n; x++) {
					tt += v[j][x] * v[i][x];
				}
				if (i == j)
					tt -= 1.0;
				tt = fabs(tt);
				if (tt > 0.0000001) {
					bad = 1;
					printf("V Error at %d %d = %f \n",j,i,tt);
				}
			}
		if (bad)
			printf("V,Vt == 1 Check failed!\n");
		}

#ifdef NEVER
		printf("A = %f %f %f\n    %f %f %f\n    %f %f %f\n\n",
		       a[0][0], a[0][1], a[0][2], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2]);

		printf("u = %f %f %f\n    %f %f %f\n    %f %f %f\n\n",
		       u[0][0], u[0][1], u[0][2], u[1][0], u[1][1], u[1][2], u[2][0], u[2][1], u[2][2]);

		printf("w = %f %f %f\n\n", w[0],w[1],w[2]);

		printf("V = %f %f %f\n    %f %f %f\n   %f %f %f\n\n",
	          v[0][0], v[0][1], v[0][2], v[1][0], v[1][1], v[1][2], v[2][0], v[2][1], v[2][2]);;

		printf("A2 = %f %f %f\n    %f %f %f\n    %f %f %f\n\n",
		       a2[0][0], a2[0][1], a2[0][2], a2[1][0], a2[1][1],
		       a2[1][2], a2[2][0], a2[2][1], a2[2][2]);
#endif

		free_dmatrix(a,  0,m-1, 0,n-1);
		free_dmatrix(t,  0,m-1, 0,n-1);
		free_dmatrix(a2, 0,m-1, 0,n-1);
		free_dmatrix(u,  0,m-1, 0,n-1);
		free_dvector(w,  0,n-1);
		free_dmatrix(v,  0,n-1, 0,n-1);
	}
	return 0;
}





































