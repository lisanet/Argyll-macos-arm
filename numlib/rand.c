/* Integer and floating point random number generator routines */
/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include "numsup.h"
#include "rand.h"

/* 32 bit pseudo random sequencer based on XOR feedback */
/* generates number between 1 and 4294967295 */
#define PSRAND32(S) (((S) & 0x80000000) ? (((S) << 1) ^ 0xa398655d) : ((S) << 1))

/* 32 bit linear congruent generator */
/* generates number between 0 and 4294967295 */
/* (From Knuth & H.W.Lewis) */
#define PSRAND32L(S) ((S) * 1664525L + 1013904223L)

/* - - - - - - - - - - - - - - - */
/* Global state random generator */

static rand_state g_rand = { 0 }; /* Use default seed for global state generator */

/* Return a 32 bit number between 0 and 4294967295 */
/* Use Knuth shuffle to improve PSRAND32 sequence */
unsigned int
rand32(					/* Return 32 bit random number */
unsigned int seed		/* Optional seed. Non-zero re-initialized with that seed */
) {
	return rand32_th(NULL, seed);
}

/* return a random number between 0.0 and 1.0 */
/* based on rand32 */
double ranno(void) {
	return rand32_th(NULL, 0) / 4294967295.0;
}

/* Return a uniform random double in the range min to max */
double
d_rand(double min, double max) {
	return d_rand_th(NULL, min, max);
}

/* Return a squared distribution random double in the range min to max */
double
d2_rand(double min, double max) {
	return d2_rand_th(NULL, min, max);
}

/* Return a random integer in the range min to max inclusive */
int
i_rand(int min, int max) {
	return i_rand_th(NULL, min, max);
}

/* Return a random floating point number with a gausian/normal */
/* distribution, centered about 0.0, with standard deviation 1.0 */
/* This uses the Box-Muller transformation */
double norm_rand(void) {
	return norm_rand_th(NULL);
}

/* - - - - - - - - - - - - - - - */
/* Explicit state random generator */

/* Init rand_state to default */
void rand_init(rand_state *p) {
	if (p == NULL)
		p = &g_rand;
	memset((void *)p, 0, sizeof(rand_state));
}

/* Return a 32 bit number between 0 and 4294967295 */
/* Use Knuth shuffle to improve PSRAND32 sequence */
unsigned int
rand32_th(rand_state *p,
unsigned int seed		/* Optional seed. Non-zero re-initialized with that seed */
) {
	int i;

	if (p == NULL)
		p = &g_rand;

	if (seed != 0) {
//printf("~1 rand 0x%x seed 0x%x\n",p,seed);
		rand_init(p);
		p->ran = seed;
	}

	/* Init random storage locations */
	if (p->pvs_inited == 0) {
		if (p->ran == 0)
			p->ran = RAND_SEED;
		for (i = 0; i < RAND_TSIZE; i++)
  			p->pvs[i] = p->ran = PSRAND32(p->ran);
		p->last = p->ran;
		p->pvs_inited = 1;
	}
	i = p->last % RAND_TSIZE;				/* New location */
	p->last = p->pvs[i];					/* Value generated */
  	p->pvs[i] = p->ran = PSRAND32(p->ran);	/* New value */

//printf("~1 rand 0x%x ret 0x%x\n",p,p->last-1);
	return p->last-1;
}

/* return a random number between 0.0 and 1.0 */
/* based on rand32 */
double ranno_th(rand_state *p) {
	return rand32_th(p, 0) / 4294967295.0;
}

/* Return a uniform random double in the range min to max */
double
d_rand_th(rand_state *p, double min, double max) {
	return min + (max - min) * ranno_th(p);
}

/* Return a squared distribution random double in the range min to max */
double
d2_rand_th(rand_state *p, double min, double max) {
	double val = ranno_th(p);
	return min + (max - min) * val * val;
}

/* Return a random integer in the range min to max inclusive */
int
i_rand_th(rand_state *p, int min, int max) {
	return min + (int)floor(0.5 + ((double)(max - min)) * ranno_th(p));
}

/* Return a random floating point number with a gausian/normal */
/* distribution, centered about 0.0, with standard deviation 1.0 */
/* This uses the Box-Muller transformation */
double norm_rand_th(rand_state *p) {
	if (p == NULL)
		p = &g_rand;

	if (p->r2 == 0) {				/* No previously calculated number */
		double v1, v2, t1, t2, r1;
		do {
			v1 = d_rand_th(p, -1.0, 1.0);
			v2 = d_rand_th(p, -1.0, 1.0);
			t1 = v1 * v1 + v2 * v2;
		} while (t1 == 0.0 || t1 >= 1.0);
		t2 = sqrt(-2.0 * log(t1)/t1);
		p->nr2 = v2 * t2;			/* One for next time */
		p->r2 = 1;
		r1 = v1 * t2;
		return r1;
	} else {						/* Return previously calculated number */
		p->r2 = 0;
		return p->nr2;
	}
}







