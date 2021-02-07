#ifndef RAND_H
#define RAND_H

/*
 * Copyright 1998 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifdef __cplusplus
	extern "C" {
#endif

/* - - - - - - - - - - - - - - - */
/* Global state random generator */

/* Return a random number between 0 and 4294967294 */
unsigned int
rand32(						/* Return 32 bit random number */
unsigned int seed);			/* Optional seed. Non-zero re-initialized with that seed */

/* Return a random integer in the range min to max inclusive */
int i_rand(int min, int max);

/* Return a uniform random double in the range min to max inclusive */
double d_rand(double min, double max);

/* Return a squared distribution random double in the range min to max inclusive */
double d2_rand(double min, double max);

/* Return a random floating point number with a gausian/normal */
/* distribution, centered about 0.0, with standard deviation 1.0 */
/* and an average deviation of 0.564 */
double norm_rand(void);

/* Scale normal value by this to give it a mean absolute deviation of 1.0 */
#define NORM_RAND_ABS_SCALE 0.62665706865775

/* - - - - - - - - - - - - - - - */
/* Explicit state random generator */
/* Use NULL for global state */

#define RAND_TSIZE 2843      	/* Prime */
#define RAND_SEED 0x12345678    /* Default seed */

/* Should set to 0 to default intialize */
typedef struct {
	int pvs_inited;
	unsigned int ran, last;
	unsigned int pvs[RAND_TSIZE];

	/* normal distribution 2nd value */
	int r2;			/* 2nd value available */
	double nr2;		/* 2nd value */
} rand_state;

/* Init rand_state to default */
void rand_init(rand_state *p);

/* Return a random number between 0 and 4294967294 */
unsigned int
rand32_th(rand_state *p,
unsigned int seed);			/* Optional seed. Non-zero re-initialized with that seed */

/* Return a random integer in the range min to max inclusive */
int i_rand_th(rand_state *p, int min, int max);

/* Return a uniform random double in the range min to max inclusive */
double d_rand_th(rand_state *p, double min, double max);

/* Return a squared distribution random double in the range min to max inclusive */
double d2_rand_th(rand_state *p, double min, double max);

/* Return a random floating point number with a gausian/normal */
/* distribution, centered about 0.0, with standard deviation 1.0 */
/* and an average deviation of 0.564 */
double norm_rand_th(rand_state *p);

#ifdef __cplusplus
	}
#endif

#endif /* RAND_H */
