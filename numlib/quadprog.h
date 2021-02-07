
#ifndef QUADPROG_H
#define QUADPROG_H

#ifdef __cplusplus
	extern "C" {
#endif

/* 
	Quadradic Programming solution.

	Based on Luca Di Gaspero's QuadProgpp++, made available under the MIT license.

	Translation to C by Graeme W. Gill, Copyright 2017, also licensed under the
	MIT license.
*/

/* 

 The quadprog_solve() function implements the algorithm of Goldfarb and Idnani 
 for the solution of a (convex) Quadratic Programming problem
 by means of an active-set dual method.
	 
The problem is in the form:

min 0.5 * x G x + g0 x
s.t.
    CE^t x + ce0 = 0
    CI^t x + ci0 >= 0
	 
 The matrix and vectors dimensions are as follows:
      G: n * n
	 g0: n
				
	 CE: n * p
	ce0: p
				
	 CI: n * m
    ci0: m

      x: n
 
 References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
             strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

 NOTE that this function doesn't work if the objective function (i.e. G + g0) has lower dimension
 than n. Using dummy terms doesn't seem to work around this.

 It only works with convex problems - i.e. those where it is possible to reach
 the optimum starting anywhere on the constraint surface.

*/

#define QP_INFEASIBLE 1.0E300

double quadprog(	/* Return solution cost, QP_INFEASIBLE if infeasible/error */
    double *x,		/* Return x[n] value */
	double **G,		/* G[n][n]	Quadratic combination matrix - modified */
	double *g0,		/* g0[n]    Direct vector */
    double **CE,	/* CE[n][p]	Equality constraint matrix */
	double *ce0, 	/* ce0[p]   Equality constraing constants */ 
	double **CI,	/* CI[n][m]	Constraint matrix */
    double *ci0,	/* cie[m]	Constraint constants */ 
	int n,			/* Number of variables */ 
	int p,			/* Number of equalities */ 
	int m			/* Number of constraints */ 
);

#ifdef __cplusplus
	}
#endif

#endif /*define QUADPROG_H */
