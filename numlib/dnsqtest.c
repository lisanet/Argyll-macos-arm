
/*
 * Copyright 1999 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* Example use of dnsqe()                                                  */
/*                                                                         */
/*       The problem is to determine the values of X(1), X(2), ..., X(9),  */
/*       which solve the system of tridiagonal equations                   */
/*                                                                         */
/*       (3-2*X(1))*X(1)           -2*X(2)                   = -1          */
/*               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8   */
/*                                   -X(8) + (3-2*X(9))*X(9) = -1          */
/*                                                                         */
/*       Final approximate solution:                                       */
/*                                                                         */
/*       -0.5706545E+00                                                    */
/*       -0.6816283E+00                                                    */
/*       -0.7017325E+00                                                    */
/*       -0.7042129E+00                                                    */
/*       -0.7013690E+00                                                    */
/*       -0.6918656E+00                                                    */
/*       -0.6657920E+00                                                    */
/*       -0.5960342E+00                                                    */
/*       -0.4164121E+00                                                    */

#include "numlib.h"

/* Compute norm of a vector */
static double denorm(int n, double *x);

int fcn(void *fdata, int n, double *x, double *fvec, int iflag);

double expect[9] = {
       -0.5706545E+00,
       -0.6816283E+00,
       -0.7017325E+00,
       -0.7042129E+00,
       -0.7013690E+00,
       -0.6918656E+00,
       -0.6657920E+00,
       -0.5960342E+00,
       -0.4164121E+00 };

int main(void)
{
	int n = 9 /* 9 */;			/* Problem vector size */
	double x[9];		/* Function input values */
	double fvec[9];		/* Function output values */
	double ss;			/* Search area */
	int info, j;
	double fnorm;
	int nprint = 0;		/* Itteration debugging print = off */
	double tol;

	error_program = "dnsqtest";	/* Set global error reporting string */

	/*	 Driver for dnsqe example. */
	/* Not supplying Jacobian, use approximation */

	/*	 The following starting values provide a rough solution. */
	for (j = 1; j <= 9; ++j) {
		x[j - 1] = -1.f;
	}
	ss = 0.1;

	nprint = 0;

	/*	 Set tol to the square root of the machine precision. */
	/*	 Unless high precision solutions are required, */
	/*	 this is the recommended setting. */

	tol = M_SQRT_DIVER;

	info = dnsqe(NULL, fcn, NULL, n, x, ss, fvec, tol, tol, 0, nprint);
	fnorm = denorm(n, fvec);
	fprintf(stdout,"Final L2 norm of the residuals = %e\n",fnorm);
	fprintf(stdout,"Exit return value = %d (1 = sucess)\n",info);
	fprintf(stdout,"Final approximate solution:\n");
	for (j = 0; j < n; j++) {
		fprintf(stdout,"x[%d] = %f, expect %f\n",j,x[j], expect[j]);
	}

	return 0;
} /* main() */

/* Function being solved */
int fcn(void *fdata, int n, double *x, double *fvec, int iflag)
{
	double temp, temp1, temp2;
	int k;

	/* Function Body */
	for (k = 0; k < n; ++k) {
		temp = (3.0 - 2.0 * x[k]) * x[k];
		temp1 = 0.0;
		if (k != 0) {
			temp1 = x[k-1];
		}
		temp2 = 0.0;
		if (k != ((n)-1))
			temp2 = x[k+1];
		fvec[k] = temp - temp1 - 2.0 * temp2 + 1.0;
		if (iflag == 0)
			printf("x[%d] = %f, fvec[%d] + %f\n",k,x[k],k,fvec[k]);

#ifdef DEBUG
printf("~~ x[%d] = %f, fvec[%d] + %f\n",k,x[k],k,fvec[k]);
#endif /* DEBUG */
	}
	/* Return < 0 to abort */
	return 0;
}

/* - - - - - - - - - - - - - - - - - - - */

static double denorm(
	int n,			/* Size of x[] */
	double x[])		/* Input vector */
{
	/* Initialized data */
	static double rdwarf = 3.834e-20;
	static double rgiant = 1.304e19;

	/* Local variables */
	static double xabs, x1max, x3max;
	static int i;
	static double s1, s2, s3, agiant, floatn;
	double ret_val, td;

	s1 = 0.0;	/* Large component */
	s2 = 0.0;	/* Intermedate component */
	s3 = 0.0;	/* Small component */
	x1max = 0.0;
	x3max = 0.0;
	floatn = (double) (n + 1);
	agiant = rgiant / floatn;
	for (i = 0; i < n; i++) {
		xabs = (td = x[i], fabs(td));

		/* Sum for intermediate components. */
		if (xabs > rdwarf && xabs < agiant) {
			td = xabs;				 	/* Computing 2nd power */
			s2 += td * td;

		/* Sum for small components. */
		} else if (xabs <= rdwarf) {
			if (xabs <= x3max) {
				if (xabs != 0.0) {		/* Computing 2nd power */
				td = xabs / x3max;
				s3 += td * td;
				}
			} else { /* Computing 2nd power */
				td = x3max / xabs;
				s3 = 1.0 + s3 * (td * td);
				x3max = xabs;
			}

		/* Sum for large components. */
		} else {
			if (xabs <= x1max) {		/* Computing 2nd power */
				td = xabs / x1max;
				s1 += td * td;
			} else {					/* Computing 2nd power */
				td = x1max / xabs;
				s1 = 1.0 + s1 * (td * td);
				x1max = xabs;
			}
		}
	}

	/* Calculation of norm. */
	if (s1 != 0.0) {		/* Large is present */
		ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
	} else {				/* Medium and small are present */
		if (s2 == 0.0) {
			ret_val = x3max * sqrt(s3);		/* Small only */
		} else {
			if (s2 >= x3max) {		/* Medium larger than small */
				ret_val = sqrt(s2 * (1.0 + x3max / s2 * (x3max * s3)));
			} else {				/* Small large than medium */
				ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
			}
		}
	}
	return ret_val;
}

