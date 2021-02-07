/*

 This file contains just an example on how to set-up the matrices for using with
 the quadprog() function.
 
 The test problem is the following:
 
 Given:
 G =  4 -2   g0^T = [6 0]
     -2  4       
 
 Solve:
 min f(x) = 1/2 x G x + g0 x
 s.t.
   x_1 + x_2 = 3
   x_1 >= 0
   x_2 >= 0
   x_1 + x_2 >= 2
 
 The solution is x^T = [1 2] and f(x) = 12
 
*/

#include "stdio.h"
#include "numlib.h"
#include "quadprog.h"

int main() {
	double **GG, **G, **CE, **CI;
	double *g0, *ce0, *ci0, *x;
	int n, m, p;
	int i, j;
	double f, sum = 0.0;
	
	n = 2;

	/* Orginal, unchanged G */
	GG = dmatrix(0, n-1, 0, n-1);
	GG[0][0] = 4.0;
	GG[0][1] = -2.0;
	GG[1][0] = -2.0;
	GG[1][1] = 4.0;

	/* Working copy - gets modified */
	G = dmatrix(0, n-1, 0, n-1);
	copy_dmatrix(G, GG, 0, n-1, 0, n-1);

	g0 = dvector(0, n-1);
	g0[0] = 6.0;
	g0[1] = 0.0;

	p = 1;

	CE = dmatrix(0, n-1, 0, p-1);
	CE[0][0] = 1.0;
	CE[1][0] = 1.0;

	ce0 = dvector(0, p-1);
	ce0[0] = -3.0;
	
	m = 3;
	CI = dmatrix(0, n-1, 0, m-1);
	CI[0][0] = 1.0;
	CI[0][1] = 0.0;
	CI[0][2] = 1.0;
	CI[1][0] = 0.0;
	CI[1][1] = 1.0;
	CI[1][2] = 1.0;
	
	ci0 = dvector(0, m-1);
	ci0[0] = 0.0;
	ci0[1] = 0.0;
	ci0[2] = -2.0;

	x = dvector(0, n-1);

	f = quadprog(x, G, g0, CE, ce0, CI, ci0, n, p, m);

	if (f == QP_INFEASIBLE)
		printf("Failed to find solution\n");
	else {
		printf("Function value at %f %f = %f\n",x[0],x[1],f);
	
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				sum += x[i] * GG[i][j] * x[j];
		sum *= 0.5;	
		
		for (i = 0; i < n; i++)
			sum += g0[i] * x[i];
	
		printf("Double check funtion value %f\n",sum);
	}

	free_dmatrix(GG, 0, n-1, 0, n-1);
	free_dmatrix(G, 0, n-1, 0, n-1);
	free_dvector(g0, 0, n-1);
	free_dmatrix(CE, 0, n-1, 0, p-1);
	free_dvector(ce0, 0, p-1);
	free_dmatrix(CI, 0, n-1, 0, m-1);
	free_dvector(ci0, 0, m-1);
	free_dvector(x, 0, n-1);

	return 0;
}
