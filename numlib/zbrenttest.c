/* Do a test of the brent 1d root finder */

/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* Solves (3 - 2*X) * X = -1  */

#include "numlib.h"

int its = 0;
double fcn(void *fdata, double tp);

/* Expected solution */
double expect = -2.80776403408306e-01;

int main(void)
{
	int rv;
	double x1, x2, soln;

	/* Attempt bracket the solution */
	x1 = -0.1; 
	x2 =  0.1;
	rv = zbrac(&x1, &x2, fcn, NULL);
	if (rv != 0) {
		error("zbrack failed with rv = %d\n",rv);
	}
	printf("Got bracket %f, %f in %d itterations\n",x1,x2,its);

	/* Solve the equation */
	its = 0;
	rv = zbrent(&soln, x1, x2, 1e-6, fcn, NULL);
	if (rv != 0) {
		error("xbrent failed with rv = %d\n",rv);
	}

	printf("Solution = %f, expected = %f in %d itterations\n",soln, expect, its);
	return 0;

} /* main() */

/* Function being solved */
double fcn(void *fdata, double tp)
{
	double temp;
	its++;
	temp = ((3.0 - 2.0 * tp) * tp) + 1.0;
/* printf("~~ %f returns %f\n",tp,temp); */
	return temp;
}

