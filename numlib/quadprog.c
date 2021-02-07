
/* 
	Quadradic Programming soliution.

	Based on Luca Di Gaspero's QuadProgpp++, made available with
	the following license:

		MIT License
		
		Copyright (c) 2007-2016 Luca Di Gaspero
		
		Permission is hereby granted, free of charge, to any person obtaining a copy
		of this software and associated documentation files (the "Software"), to deal
		in the Software without restriction, including without limitation the rights
		to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
		copies of the Software, and to permit persons to whom the Software is
		furnished to do so, subject to the following conditions:
		
		The above copyright notice and this permission notice shall be included in all
		copies or substantial portions of the Software.
		
		THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
		IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
		FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
		AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
		LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
		OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
		SOFTWARE.

	Translation to C by Graeme W. Gill, Copyright 2017, also licensed under the
	above license.
*/

/* 

 The quadprog_solve() function implements the algorithm of Goldfarb and Idnani 
 for the solution of a (convex) Quadratic Programming problem
 by means of an active-set dual method.
	 
The problem is in the form:

min 0.5 * x G x + g0 x
s.t.
    CE^T x + ce0 = 0
    CI^T x + ci0 >= 0
	 
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

*/

#include "numsup.h"
#include "quadprog.h"

#undef TRACE_SOLVER		/* Print progress */

/* Utility functions for updating some data needed by the solution method  */
INLINE static void compute_d(double *d, double **J, double *np, int n);
INLINE static void update_z(double *z, double **J, double *d, int iq, int n);
INLINE static void update_r(double **R, double *r, double *d, int iq, int n);
static int add_constraint(double **R, double **J, double *d, int *piq, double *prnorm, int n);
INLINE static void delete_constraint(double **R, double **J, int *A, double *u,
	                   int n, int p, int *piq, int l);

/* Utility functions for computing the Cholesky decomposition and solving */
/* linear systems */
static void cholesky_decomposition(double **A, int n);
static void cholesky_solve(double **L, double *x,  double *b, int n);
INLINE static void forward_elimination(double **L, double *y, double *b, int n);
INLINE static void backward_elimination(double **U, double *x, double *y, int n);

/* Utility functions for computing the scalar product and the euclidean  */
/* distance between two numbers */
INLINE static double scalar_product(double *x, double *y, int n);
INLINE static double distance(double a, double b);

/* Utility functions for printing vectors and matrices */
static void print_matrix(char* name, double **A, int n, int m);

static void print_vector(char* name, double *v, int n);
static void print_ivector(char* name, int *v, int n);

#define FMIN(A,B) ((A) < (B) ? (A) : (B))

#define FMAX(A,B) ((A) > (B) ? (A) : (B))

#define ASZ 10
#define VSZ 20

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
) {
	int i, j, k, l;	/* indices */
	int ip;					 	/* index of the constraint to be added to the active set */
	double f_value = QP_INFEASIBLE, psi, c1, c2, sum, ss, R_norm;
	double inf;
	double t, t1, t2;		/* t is the step length, which is the minimum of the partial */
							/* step length t1 and the full step length t2 */
	int q, iq, iter = 0;

	double **R, **J;
	double *s, *z, *r, *d, *np, *u, *x_old, *u_old;
	int *A, *A_old, *iai;
	short *iaexcl;
	
	double *_R[ASZ], __R[ASZ][ASZ], *_J[ASZ], __J[ASZ][ASZ];
	double _s[VSZ], _z[VSZ], _r[VSZ], _d[VSZ], _np[VSZ], _u[VSZ], _x_old[VSZ], _u_old[VSZ];
	int _A[VSZ], _A_old[VSZ], _iai[VSZ];
	short _iaexcl[VSZ];
	
	/* Avoid malloc's if we can.. */
	if (n <= ASZ) {
		R = _R;
		J = _J;
		for (i = 0; i < n; i++) {
			_R[i] = __R[i];
			_J[i] = __J[i];
		}
	} else {
		R = dmatrix(0, n-1, 0, n-1);
		J = dmatrix(0, n-1, 0, n-1);
	}

	if (n <= VSZ) {
		x_old = _x_old;
		z = _z;
		d = _d;
		np = _np;
	} else {
		x_old = dvector(0, n-1);
		z = dvector(0, n-1);
		d = dvector(0, n-1);
		np = dvector(0, n-1);
	}

	if ((m + p) <= VSZ) {
		s = _s;
		r = _r;
		u = _u;
		u_old = _u_old;
		A = _A;
		A_old = _A_old;
		iai = _iai;
		iaexcl = _iaexcl;
	} else {
		s = dvector(0, m + p-1);
		r = dvector(0, m + p-1);
		u = dvector(0, m + p-1);
		u_old = dvector(0, m + p-1);
		A = ivector(0, m + p-1);
		A_old = ivector(0, m + p-1);
		iai = ivector(0, m + p-1);
		iaexcl = svector(0, m + p-1);
	}

	q = 0;	/* size of the active set A (containing the indices of the active constraints) */

#ifdef TRACE_SOLVER
	printf("\nStarting solve_quadprog\n");
	print_matrix("G", G, n, n);
	print_vector("g0", g0, n);
	print_matrix("CE", CE, n, p);
	print_vector("ce0", ce0, p);
	print_matrix("CI", CI, n, m);
	print_vector("ci0", ci0, m);
#endif	
	
	/*
	 * Preprocessing phase
	 */
	
	/* compute the trace of the original matrix G */
	c1 = 0.0;
	for (i = 0; i < n; i++)
		c1 += G[i][i];
	/* decompose the matrix G in the form L^T L */
	cholesky_decomposition(G, n);
#ifdef TRACE_SOLVER
	print_matrix("G", G, n, n);
#endif

	/* initialize the matrix R */
	for (i = 0; i < n; i++) {
		d[i] = 0.0;
		for (j = 0; j < n; j++)
			R[i][j] = 0.0;
	}
	R_norm = 1.0;	/* this variable will hold the norm of the matrix R */
	
	/* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
	c2 = 0.0;
	for (i = 0; i < n; i++) {
		d[i] = 1.0;
		forward_elimination(G, z, d, n);
		for (j = 0; j < n; j++)
			J[i][j] = z[j];
		c2 += z[i];
		d[i] = 0.0;
	}
#ifdef TRACE_SOLVER
	print_matrix("J", J, n, n);
#endif
	
	/* c1 * c2 is an estimate for cond(G) */
	
	/* Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x  */
	/* this is a feasible point in the dual space */
	/* x = G^-1 * g0 */
	cholesky_solve(G, x, g0, n);
	for (i = 0; i < n; i++)
		x[i] = -x[i];
	/* and compute the current solution value */ 
	f_value = 0.5 * scalar_product(g0, x, n);
#ifdef TRACE_SOLVER
	printf("Unconstrained solution: f_value %f\n",f_value);
	print_vector("x", x, n);
#endif
	
	/* Add equality constraints to the working set A */
	iq = 0;
	for (i = 0; i < p; i++) {
		for (j = 0; j < n; j++)
			np[j] = CE[j][i];
		compute_d(d, J, np, n);
		update_z(z, J, d, iq, n);
		update_r(R, r, d, iq, n);
#ifdef TRACE_SOLVER
		print_matrix("R", R, n, n);
		print_vector("z", z, n);
		print_vector("r", r, iq);
		print_vector("d", d, n);
#endif
		
		/* compute full step length t2: i.e., the minimum step in primal space s.t. */
		/* the contraint becomes feasible */
		t2 = 0.0;
		if (fabs(scalar_product(z, z, n)) > DBL_EPSILON) /* i.e. z != 0 */
			t2 = (-scalar_product(np, x, n) - ce0[i]) / scalar_product(z, np, n);
		
		/* set x = x + t2 * z */
		for (k = 0; k < n; k++)
			x[k] += t2 * z[k];
		
		/* set u = u+ */
		u[iq] = t2;
		for (k = 0; k < iq; k++)
			u[k] -= t2 * r[k];
		
		/* compute the new solution value */
		f_value += 0.5 * (t2 * t2) * scalar_product(z, np, n);
		A[i] = -i - 1;
		
		if (!add_constraint(R, J, d, &iq, &R_norm, n)) {		
			/* Equality constraints are linearly dependent */
#ifdef TRACE_SOLVER
			printf("Constraints are linearly dependent\n");
#endif
			goto done;
		}
	}
	
	/* set iai = K \ A */
	for (i = 0; i < m; i++)
		iai[i] = i;

	/* Hmm. Use of crossing goto's is a little Fortran like... */

l1:	iter++;
#ifdef TRACE_SOLVER
	print_vector("x", x, n);
#endif
	/* step 1: choose a violated constraint */
	for (i = p; i < iq; i++) {
		ip = A[i];
		iai[ip] = -1;
	}
	
	/* compute s[x] = ci^T * x + ci0 for all elements of K \ A */
	ss = 0.0;
	psi = 0.0;	/* this value will contain the sum of all infeasibilities */
	ip = 0;		/* ip will be the index of the chosen violated constraint */
	for (i = 0; i < m; i++) {
		iaexcl[i] = 1;
		sum = 0.0;
		for (j = 0; j < n; j++)
			sum += CI[j][i] * x[j];
		sum += ci0[i];
		s[i] = sum;
		psi += FMIN(0.0, sum);
	}
#ifdef TRACE_SOLVER
	print_vector("s", s, m);
#endif
	
	if (fabs(psi) <= m * DBL_EPSILON * c1 * c2* 100.0) {
		/* numerically there are no infeasibilities anymore */
#ifdef TRACE_SOLVER
		printf("numerically there are no infeasibilities anymore\n");
#endif
		q = iq;
		goto done;
	}
	
	/* save old values for u and A */
	for (i = 0; i < iq; i++) {
		u_old[i] = u[i];
		A_old[i] = A[i];
	}
	
	/* and for x */
	for (i = 0; i < n; i++)
		x_old[i] = x[i];
	
l2:	/* Step 2: check for feasibility and determine a new S-pair */
	for (i = 0; i < m; i++) {
		if (s[i] < ss && iai[i] != -1 && iaexcl[i]) {
			ss = s[i];
			ip = i;
		}
	}
	if (ss >= 0.0) {
#ifdef TRACE_SOLVER
		printf(" ss %f >= 0.0\n",ss);
#endif
		q = iq;
		goto done;
	}
	
	/* set np = n[ip] */
	for (i = 0; i < n; i++)
		np[i] = CI[i][ip];
	/* set u = [u 0]^T */
	u[iq] = 0.0;
	/* add ip to the active set A */
	A[iq] = ip;
	
#ifdef TRACE_SOLVER
	printf("Trying with constraint %d\n",ip);
	print_vector("np", np, n);
#endif
	
l2a:
	/* Step 2a: determine step direction */
	/* compute z = H np: the step direction in the primal space (through J, see the paper) */
	compute_d(d, J, np, n);
	update_z(z, J, d, iq, n);
	/* compute N* np (if q > 0): the negative of the step direction in the dual space */
	update_r(R, r, d, iq, n);
#ifdef TRACE_SOLVER
	printf("Step direction z\n");
	print_vector("z", z, n);
	print_vector("r", r, iq + 1);
	print_vector("u", u, iq + 1);
	print_vector("d", d, n);
	print_ivector("A", A, iq + 1);
#endif
	
	/* Step 2b: compute step length */
	l = 0;
	/* Compute t1: partial step length (maximum step in dual space */
	/* without violating dual feasibility */
	t1 = QP_INFEASIBLE; /* +inf */
	/* find the index l s.t. it reaches the minimum of u+[x] / r */
	for (k = p; k < iq; k++) {
		if (r[k] > 0.0) {
			if (u[k] / r[k] < t1) {
				t1 = u[k] / r[k];
				l = A[k];
			}
		}
	}
	/* Compute t2: full step length (minimum step in primal space */
	/* such that the constraint ip becomes feasible */
	if (fabs(scalar_product(z, z, n)) > DBL_EPSILON) { /* i.e. z != 0 */
		t2 = -s[ip] / scalar_product(z, np, n);
		if (t2 < 0) /* patch suggested by Takano Akio for handling numerical inconsistencies */
			t2 = QP_INFEASIBLE;
	} else
		t2 = QP_INFEASIBLE; /* +inf */
	
	/* the step is chosen as the minimum of t1 and t2 */
	t = FMIN(t1, t2);
#ifdef TRACE_SOLVER
	printf("Step size: %e (t1 = %e, t2 = %e\n",t,t1,t2);
#endif
	
	/* Step 2c: determine new S-pair and take step: */
	
	/* case (i): no step in primal or dual space */
	if (t >= QP_INFEASIBLE) {
		/* QPP is infeasible */
		// FIXME: unbounded to raise
#ifdef TRACE_SOLVER
		printf(" QPP is infeasible\n");
#endif
		q = iq;
		f_value = QP_INFEASIBLE;
		goto done;
	}
	/* case (ii): step in dual space */
	if (t2 >= QP_INFEASIBLE) {
		/* set u = u +	t * [-r 1] and drop constraint l from the active set A */
		for (k = 0; k < iq; k++)
			u[k] -= t * r[k];
		u[iq] += t;
		iai[l] = l;
		delete_constraint(R, J, A, u, n, p, &iq, l);
#ifdef TRACE_SOLVER
		printf(" in dual space: %f\n",f_value);
		print_vector("x", x, n);
		print_vector("z", z, n);
		print_ivector("A", A, iq + 1);
#endif
		goto l2a;
	}
	
	/* case (iii): step in primal and dual space */
	
	/* set x = x + t * z */
	for (k = 0; k < n; k++)
		x[k] += t * z[k];
	/* update the solution value */
	f_value += t * scalar_product(z, np, n) * (0.5 * t + u[iq]);
	/* u = u + t * [-r 1] */
	for (k = 0; k < iq; k++)
		u[k] -= t * r[k];
	u[iq] += t;
#ifdef TRACE_SOLVER
	printf(" in both spaces: %f\n",f_value);
	print_vector("x", x, n);
	print_vector("u", u, iq + 1);
	print_vector("r", r, iq + 1);
	print_ivector("A", A, iq + 1);
#endif
	
	if (fabs(t - t2) < DBL_EPSILON) {
#ifdef TRACE_SOLVER
		printf("Full step has been taken %f\n",t);
		print_vector("x", x, n);
#endif
		/* full step has been taken */
		/* Add constraint ip to the active set*/
		if (!add_constraint(R, J, d, &iq, &R_norm, n)) {
#ifdef TRACE_SOLVER
			printf("add_constraint failed - removing constraint ant step\n",t);
#endif
			iaexcl[ip] = 0;
			delete_constraint(R, J, A, u, n, p, &iq, ip);
#ifdef TRACE_SOLVER
			print_matrix("R", R, n, n);
			print_ivector("A", A, iq);
			print_ivector("iai", iai, m+p);
#endif
			for (i = 0; i < m; i++)
				iai[i] = i;
			for (i = p; i < iq; i++) {
				A[i] = A_old[i];
				u[i] = u_old[i];
				iai[A[i]] = -1;
			}
			for (i = 0; i < n; i++)
				x[i] = x_old[i];
			goto l2; /* go to step 2 */
		} else
			iai[ip] = -1;
#ifdef TRACE_SOLVER
		print_matrix("R", R, n, n);
		print_ivector("A", A, iq);
		print_ivector("iai", iai, m + p);
#endif
		goto l1;
	}
	
	/* a patial step has taken */
#ifdef TRACE_SOLVER
	printf("Partial step has taken %f\n",t);
	print_vector("x", x, n);
#endif
	/* drop constraint l */
	iai[l] = l;
	delete_constraint(R, J, A, u, n, p, &iq, l);
#ifdef TRACE_SOLVER
	print_matrix("R", R, n, n);
	print_ivector("A", A, iq);
#endif
	
	/* update s[ip] = CI * x + ci0 */
	sum = 0.0;
	for (k = 0; k < n; k++)
		sum += CI[k][ip] * x[k];
	s[ip] = sum + ci0[ip];
	
#ifdef TRACE_SOLVER
	print_vector("s", s, m);
#endif
	goto l2a;

  done:;

	/* Cleanup and return f value */
	if (n > ASZ) {
		free_dmatrix(R, 0, n-1, 0, n-1);
		free_dmatrix(J, 0, n-1, 0, n-1);
	}

	if (n > VSZ) {
		free_dvector(x_old, 0, n-1);
		free_dvector(z, 0, n-1);
		free_dvector(d, 0, n-1);
		free_dvector(np, 0, n-1);
	}

	if ((m + p) > VSZ) {
		free_dvector(s, 0, m + p-1);
		free_dvector(r, 0, m + p-1);
		free_dvector(u, 0, m + p-1);
		free_dvector(u_old, 0, m + p-1);
		free_ivector(A, 0, m + p-1);
		free_ivector(A_old, 0, m + p-1);
		free_ivector(iai, 0, m + p-1);
		free_svector(iaexcl, 0, m + p-1);
	}

#ifdef TRACE_SOLVER
	printf("Returning f_value %e\n",f_value);
	print_vector("Returning x", x, n);
#endif

	return f_value;
}

INLINE static void compute_d(double *d, double **J, double *np, int n) {
	int i, j;
	double sum;
	
	/* compute d = H^T * np */
	for (i = 0; i < n; i++) {
		sum = 0.0;
		for (j = 0; j < n; j++)
			sum += J[j][i] * np[j];
		d[i] = sum;
	}
}

INLINE static void update_z(double *z, double **J, double *d, int iq, int n) {
	int i, j;
	
	/* setting of z = H * d */
	for (i = 0; i < n; i++) {
		z[i] = 0.0;
		for (j = iq; j < n; j++)
			z[i] += J[i][j] * d[j];
	}
}

INLINE static void update_r(double **R, double *r, double *d, int iq, int n) {
	int i, j;
	double sum;
	
	/* setting of r = R^-1 d */
	for (i = iq - 1; i >= 0; i--) {
		sum = 0.0;
		for (j = i + 1; j < iq; j++)
			sum += R[i][j] * r[j];
		r[i] = (d[i] - sum) / R[i][i];
	}
}

static int add_constraint(double **R, double **J, double *d, int *piq, double *prnorm, int n) {
	int iq = *piq;
	int i, j, k;
	double cc, ss, h, t1, t2, xny;
	
#ifdef TRACE_SOLVER
	printf("Add constraint %d",iq);
#endif

	/* we have to find the Givens rotation which will reduce the element */
	/* d[j] to zero. if it is already zero we don't have to do anything, except */
	/* of decreasing j */	
	for (j = n - 1; j >= iq + 1; j--) {
		/* The Givens rotation is done with the matrix (cc cs, cs -cc). */
		/* If cc is one, then element (j) of d is zero compared with element */
		/* (j - 1). Hence we don't have to do anything.  */
		/* If cc is zero, then we just have to switch column (j) and column (j - 1)  */
		/* of J. Since we only switch columns in J, we have to be careful how we */
		/* update d depending on the sign of gs. */
		/* Otherwise we have to apply the Givens rotation to these columns. */
		/* The i - 1 element of d has to be updated to h. */
		cc = d[j - 1];
		ss = d[j];
		h = distance(cc, ss);
		if (fabs(h) < DBL_EPSILON)		/* h == 0 */
			continue;
		d[j] = 0.0;
		ss = ss / h;
		cc = cc / h;
		if (cc < 0.0) {
			cc = -cc;
			ss = -ss;
			d[j - 1] = -h;
		} else
			d[j - 1] = h;
		xny = ss / (1.0 + cc);
		for (k = 0; k < n; k++) {
			t1 = J[k][j - 1];
			t2 = J[k][j];
			J[k][j - 1] = t1 * cc + t2 * ss;
			J[k][j] = xny * (t1 + J[k][j - 1]) - t2;
		}
	}
	/* update the number of constraints added*/
	*piq = ++iq;
	/* To update R we have to put the iq components of the d vector into column iq - 1 of R */
	for (i = 0; i < iq; i++)
		R[i][iq - 1] = d[i];
#ifdef TRACE_SOLVER
	printf("/%d\n", iq);
	print_matrix("R", R, iq, iq);
	print_matrix("J", J, n, n);
	print_vector("d", d, iq);
#endif
	
	if (fabs(d[iq - 1]) <= DBL_EPSILON * *prnorm) {
		/* problem degenerate */
		return 0;
	}
	*prnorm = (double)FMAX(*prnorm, fabs(d[iq - 1]));

	return 1;
}

static void delete_constraint(double **R, double **J, int *A, double *u,
	                   int n, int p, int *piq, int l) {
	int iq = *piq;
	int i, j, k, qq = -1; // just to prevent warnings from smart compilers
	double cc, ss, h, xny, t1, t2;
	
#ifdef TRACE_SOLVER
	printf("Delete constraint %d %d",l,iq);
#endif

	/* Find the index qq for active constraint l to be removed */
	for (i = p; i < iq; i++) {
		if (A[i] == l) {
			qq = i;
			break;
		}
	}
			
	/* remove the constraint from the active set and the duals */
	for (i = qq; i < iq - 1; i++) {
		A[i] = A[i + 1];
		u[i] = u[i + 1];
		for (j = 0; j < n; j++)
			R[j][i] = R[j][i + 1];
	}
			
	A[iq - 1] = A[iq];
	u[iq - 1] = u[iq];
	A[iq] = 0; 
	u[iq] = 0.0;
	for (j = 0; j < iq; j++)
		R[j][iq - 1] = 0.0;

	/* constraint has been fully removed */
	*piq = --iq;

#ifdef TRACE_SOLVER
	printf("/%d\n",iq);
#endif 
	
	if (iq == 0)
		return;
	
	for (j = qq; j < iq; j++) {
		cc = R[j][j];
		ss = R[j + 1][j];
		h = distance(cc, ss);
		if (fabs(h) < DBL_EPSILON) // h == 0
			continue;
		cc = cc / h;
		ss = ss / h;
		R[j + 1][j] = 0.0;
		if (cc < 0.0) {
			R[j][j] = -h;
			cc = -cc;
			ss = -ss;
		} else
			R[j][j] = h;
		
		xny = ss / (1.0 + cc);
		for (k = j + 1; k < iq; k++) {
			t1 = R[j][k];
			t2 = R[j + 1][k];
			R[j][k] = t1 * cc + t2 * ss;
			R[j + 1][k] = xny * (t1 + R[j][k]) - t2;
		}
		for (k = 0; k < n; k++) {
			t1 = J[k][j];
			t2 = J[k][j + 1];
			J[k][j] = t1 * cc + t2 * ss;
			J[k][j + 1] = xny * (J[k][j] + t1) - t2;
		}
	}
}

INLINE static double distance(double a, double b) {
	double a1, b1, t;
	a1 = fabs(a);
	b1 = fabs(b);
	if (a1 > b1) {
		t = (b1 / a1);
		return a1 * sqrt(1.0 + t * t);
	} else if (b1 > a1) {
		t = (a1 / b1);
		return b1 * sqrt(1.0 + t * t);
	}
	return a1 * sqrt(2.0);
}


INLINE static double scalar_product(double *x, double *y, int n) {
	int i;
	double sum;
	
	sum = 0.0;
	for (i = 0; i < n; i++)
		sum += x[i] * y[i];
	return sum;			
}

// !!! this doesn't work for semi-definite matricies, ie.
// they will have diagonal sum == 0.0 !!!
static void cholesky_decomposition(double **A, int n) {
	int i, j, k;
	double sum;
	
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			sum = A[i][j];
			for (k = i - 1; k >= 0; k--)
				sum -= A[i][k] * A[j][k];
#ifdef TRACE_SOLVER
			printf("sum[%d][%d] = %f\n",i,j,sum);
#endif
			if (i == j) {
				if (sum <= 0.0)
				{
					// raise error
					print_matrix("A", A, n, n);
					error("QuadProg:cholesky decomposition, matrix is not postive definite, sum: %e",sum);
				}
					A[i][i] = sqrt(sum);
			} else
				A[j][i] = sum / A[i][i];
		}
		for (k = i + 1; k < n; k++)
			A[i][k] = A[k][i];
	} 
}

static void cholesky_solve(double **L, double *x,  double *b, int n) {
	double *y, _y[VSZ];
	
	if (n <= VSZ)
		y = _y;
	else
		y = dvector(0, n-1);

	/* Solve L * y = b */
	forward_elimination(L, y, b, n);

	/* Solve L^T * x = y */
	backward_elimination(L, x, y, n);

	if (n > VSZ)
		free_dvector(y, 0, n-1);
}

INLINE static void forward_elimination(double **L, double *y, double *b, int n) {
	int i, j;
	
	y[0] = b[0] / L[0][0];
	for (i = 1; i < n; i++) {
		y[i] = b[i];
		for (j = 0; j < i; j++)
			y[i] -= L[i][j] * y[j];
		y[i] = y[i] / L[i][i];
	}
}

INLINE static void backward_elimination(double **U, double *x, double *y, int n) {
	int i, j;
	
	x[n - 1] = y[n - 1] / U[n - 1][n - 1];
	for (i = n - 2; i >= 0; i--) {
		x[i] = y[i];
		for (j = i + 1; j < n; j++)
			x[i] -= U[i][j] * x[j];
		x[i] = x[i] / U[i][i];
	}
}

/* --------------------------------------------------- */

static void print_matrix(char *name, double **A, int n, int m) {
	int i, j;

	printf("%s: \n",name);
	for (i = 0; i < n; i++) {
		printf(" ");
		for (j = 0; j < m; j++)
			printf("%f%s", A[i][j], j < (m-1) ? ", " : "");
		printf("\n");
	}
	printf("\n");
}

static void print_vector(char *name, double *v, int n) {
	int i;

	printf("%s: \n",name);
	printf(" ");
	for (i = 0; i < n; i++)
		printf("%f%s", v[i], i < (n-1) ? ", " : "");
	printf("\n\n");
}

static void print_ivector(char *name, int *v, int n) {
	int i;

	printf("%s: \n",name);
	printf(" ");
	for (i = 0; i < n; i++)
		printf("%d%s", v[i], i < (n-1) ? ", " : "");
	printf("\n\n");
}

