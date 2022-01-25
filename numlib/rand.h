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
/* and an average deviation of sqrt(2/pi) */
/* 0.79788456080286535587989211986876 = */
double norm_rand(void);

/* Scale normal vector components by the following to give the vector */
/* a mean absolute deviation of 1.0 for the corresponding dimension. */
/* [ Formula is 1.0/(sqrt(2.0) * gamma_func((d+1.0)/2.0)/gamma_func(d/2.0)) ] */

#define NORM_RAND_ABS_SCALE_1 1.25331413731550320
#define NORM_RAND_ABS_SCALE_2 0.79788456080287129
#define NORM_RAND_ABS_SCALE_3 0.62665706865775450
#define NORM_RAND_ABS_SCALE_4 0.53192304053525241
#define NORM_RAND_ABS_SCALE_5 0.46999280149331241
#define NORM_RAND_ABS_SCALE_6 0.42553843242820782
#define NORM_RAND_ABS_SCALE_7 0.39166066791109894
#define NORM_RAND_ABS_SCALE_8 0.36474722779559393
#define NORM_RAND_ABS_SCALE_9 0.34270308442220981
#define NORM_RAND_ABS_SCALE_10 0.32421975804054332
#define NORM_RAND_ABS_SCALE_11 0.30843277597999841
#define NORM_RAND_ABS_SCALE_12 0.29474523458229918
#define NORM_RAND_ABS_SCALE_13 0.28273004464831558
#define NORM_RAND_ABS_SCALE_14 0.27207252422983330
#define NORM_RAND_ABS_SCALE_15 0.26253504145918249
#define NORM_RAND_ABS_SCALE_16 0.25393435594782049
#define NORM_RAND_ABS_SCALE_17 0.24612660136794409
#define NORM_RAND_ABS_SCALE_18 0.23899704089211726
#define NORM_RAND_ABS_SCALE_19 0.23245290129195134

/* Table of the above values indexed by dimension */
#define NORM_RAND_ABS_SCALE_MAXD 19
extern double NORM_RAND_ABS_SCALE[NORM_RAND_ABS_SCALE_MAXD+1];

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

/* Set random dvector */
void vect_rand(double *d, double min, double max, int len);

#ifdef __cplusplus
	}
#endif

#endif /* RAND_H */
