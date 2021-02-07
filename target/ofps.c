
/* 
 * ArgyllCMS Color Correction System
 *
 * Optimised Farthest Point Sampling - NN
 *
 * Author: Graeme W. Gill
 * Date:   6/9/2004
 *
 * Copyright 2004, 2009 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* Latest version using vertex nets to reduce internal accounting overhead, */
/* in an attempt to improve performance scaling with larger numbers of points. */

/* TTBD:

	This code shouldn't exit on an error - this causes an unnecessary failure
	when ofps is used to evaluate the point distribution of other
	distribution algorithms.

	There is a bug when the ink limit == dimensions-1 (200% for CMY), and
	the number of bit mask then exceeds > 32. This is not so +/- 0.2% either side
	of 200%.
	(see "Hack to workaround pathalogical")

	There is a bug for CMYK when the ink limit == 100%
	(see "Hack to workaround pathalogical")

	One way of addressing the performance issues would be to use multiple
	threads to call dnsq. A pool could be setup, one for each CPU.

	Some profiles are too rough, and slow/stall vertex placement.
	Reducing the cache grid and/or smoothing the rspl values
	may mitigate this to some degree, and make this more robust ??

 */

/*
	Description:

		We create a function that estimates the sample positioning error at
	any location based on a weighted combination of perceptual and device distance,
	and perceptual function curvature.	

	We then initially add sampling points at the largest estimated error
	verticies of the veronoi natural neighbourhood.
	This gives us an optimal distribution measuring in estimated position
	error within a tollerance of 2:1
	
	We then iteratively improve the distribution of point nodes by
	moving them in the direction of the adjacent vertex with the
	largest estimated sampling error, so that the errors are equally
	distributed. 

	To ensure that there is a good distribution of sampling
 	points on the edges and faces of the gamut, the initial
	points are given a slighte weighting towards these
	elements, and then fastened to them. The points on
	each lower dimensional element (edge, face) is then
	optimized only amongst themselves, while higher
	dimension points are aware of the lower dimension
	ones. In this way the distribution of points on
	lower dimensional surfaces is well spread, while
	the higher dimension points take their positions into account.
 */

/*
	Failings:

	The distribution near the gamut surfaces has a characteristic
	"buffer zone" layer that is not very nice. This is because
	the surface concentrate the sufrace points forming a "force field".
	It would be good to add a tweak factor to reduce this surface "gang effect".


	The initial allocation of points to lower dimension surfaces
	is a bit haphazard. It would be nice to have some mechanism
	to add or subtract points to/from lower dimensional surfaces
	if they were over or under sampled compared to everything else.


	While the current algorithm meets many goals, such as minimizing the	
	maximum estimated error from any point in the space to the nearest
	node, and placing nodes on sub dimensional surfaces with distributions
	optimal within that sub dimensions, it has one obvious failing, and
	that is that it fails to stratify the samples. 

	So if a device is dominantly channel indepenedent, it doesn't
	take advantage of the number of samples used to fully explore
	all possible device channel values. This also applies at
	higher dimensions (ie. the CMYK values exploring response
	to different K values doesn't spread the CMY values
	evenly apart.)

	Stratification seems to be somewhat at odds with the primary goal
	of minimizing the maximum estimated error from any point in the
	space to the nearest node, but it would be good if there were
	some way of setting this balance.

	How could stratification be added to the current approach ?

	In general there are many sub-dimensions views, not all of
	which would probably be regarded as important.

	To measure spread, independent voronoi tessellations of
	these sub dimensions would be needed, and they could be
	used partly driver optimization (??).

	For 1D device channels this wouldn't be so hard to
	add, although it's hard to know how effective it would
	be, or whether it would wreck the ND optimization. It
	might also be unecessary if per channel calibration
	has been applied.

	For CMY this would need a 3D shadow veronoi.

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include "aconfig.h"
#include "numlib.h"
#include "sort.h"
#include "counters.h"
#include "icc.h"
#include "xicc.h"
#include "xcolorants.h"
#include "targen.h"
#include "rspl.h"
#include "conv.h"
#include "ofps.h"

//#include <iperf.h>

#undef DEBUG
#undef WARNINGS		/* Print warnings within DEBUG */
#undef STATS		/* Show function stats */

	/* Optimal fully adapted weightings : */
#define ADAPT_PERCWGHT 0.65			/* Degree of perceptual adaptation */
#define ADAPT_CURVWGHT 1.0			/* Degree of curvature */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifndef STANDALONE_TEST	// targen settings

# define DOOPT				/* Do optimization */
# define INDEP_SURFACE		/* Make surface point distribution and optimization independent */
# undef MAXINDEP_2D			/* Limit independent surfaces to 2D */
							/* Seems to be best for ink limited devices to #undef ? */
//# define GAMUT_EDGE_FUDGE 1.5		/* Fudge factor to counteract gamut suface barrier effect */
							// This increases edge point density as a side effect ??

# define KEEP_SURFACE		/* Keep surface points on the surface during opt. */
# define INITIAL_SURFACE_PREF 1.50	/* Extra weighting for surface points at start of seeding */
# define FINAL_SURFACE_PREF 0.80	/* Extra weighting for surface points by end of seeding */

# define SURFTOL 0.0001		/* Proportion of average spacing to force to gamut boundary */
# define RANDOM_PERTERB		/* Perpterb initial placement to break up patterns */

/* Good mode */
# define PERTERB_AMOUNT 0.5	/* and to aid surface point distribution with INDEP_SURFACE */
# define OPT_MAXITS 20			/* Number of optimisation itterations (0 to disable optimisation) */
# define OPT_TRANS_ITTERS 18		/* Numbers of itterations to transition overshoot and sepw */
# define OPT_TRANS_POW 1.6		/* Power curve to blend along */
# define OPT_INITIAL_OVERSHOOT 1.9   /* Optimisation movement initial overshoot */
# define OPT_FINAL_OVERSHOOT 0.1   /* Optimisation movement final overshoot */
# define OPT_INITIAL_SEP_WEIGHT 0.7		/* Weight to give separation of nodes during opt */
# define OPT_FINAL_SEP_WEIGHT 0.3		/* Weight to give separation of nodes during opt */
# define OPT_STOP_TOL 0.0005    	/* Stopping tollerance */

/* Fast mode */
# define PERTERB_AMOUNT_2 0.1	/* and to aid surface point distribution with INDEP_SURFACE */
# define OPT_MAXITS_2 6			/* Number of optimisation itterations (0 to disable optimisation) */
# define OPT_TRANS_ITTERS_2 5		/* Numbers of itterations to transition overshoot and sepw */
# define OPT_TRANS_POW_2 1.7		/* Power curve to blend along */
# define OPT_INITIAL_OVERSHOOT_2 1.6   /* Optimisation movement initial overshoot */
# define OPT_FINAL_OVERSHOOT_2 0.05   /* Optimisation movement final overshoot */
# define OPT_INITIAL_SEP_WEIGHT_2 0.8		/* Weight to give separation of nodes during opt */
# define OPT_FINAL_SEP_WEIGHT_2 0.3		/* Weight to give separation of nodes during opt */
# define OPT_STOP_TOL_2 0.001    	/* Stopping tollerance */

/* Diagnostic settings */
# undef DUMP_STRUCTURE		/* Dump internal node & vertex structure */
# undef DUMP_PLOT_SEED	/* Show on screen plot for each initial seed point */
# undef DUMP_PLOT		/* Show on screen plot after each itteration */
# define DUMP_VTX 1		/* Display the vertex locations too */
# define DUMP_PLA 1		/* Display the node planes too */
# define PERC_PLOT 0	/* Emit perceptive space plots */
# define DO_WAIT 1 		/* Wait for user key after each plot */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#else	/* ofps standalone test settings */

# define DOOPT			/* Do optimization */
# define INDEP_SURFACE		/* Make surface point distribution and optimization independent */
# define MAXINDEP_2D		/* Limit independent surfaces to 2D */

//# define GAMUT_EDGE_FUDGE 1.5		/* Fudge factor to counteract gamut suface barrier effect */

# define KEEP_SURFACE		/* Keep surface points on the surface during opt. */
# define INITIAL_SURFACE_PREF 1.60	/* Extra weighting for surface points at start of seeding */
# define FINAL_SURFACE_PREF 0.80	/* Extra weighting for surface points by end of seeding */

# define SURFTOL 0.0001		/* Proportion of averag spacing to force to gamut boundary */
# define RANDOM_PERTERB		/* Perpterb initial placement to break up patterns, */
# define PERTERB_AMOUNT 0.5
 
# define OPT_MAXITS 20		/* Number of optimisation itterations (0 to disable optimisation) */
# define OPT_TRANS_ITTERS 18		/* Numbers of itterations to transition overshoot and sepw */
# define OPT_TRANS_POW 2.5		/* Power curve to blend along */
# define OPT_INITIAL_OVERSHOOT 0.8   /* Optimisation movement initial overshoot */
# define OPT_FINAL_OVERSHOOT 0.1   /* Optimisation movement final overshoot */
# define OPT_INITIAL_SEP_WEIGHT 0.9		/* Weight to give separation of nodes during opt */
# define OPT_FINAL_SEP_WEIGHT 0.3		/* Weight to give separation of nodes during opt */
# define OPT_STOP_TOL 0.0005    /* Device stopping tollerance */

/* Diagnostic settings */
# undef DUMP_STRUCTURE		/* Dump internal node & vertex structure */
# undef DUMP_PLOT_SEED	/* Show on screen plot for each initial seed point */
# undef DUMP_PLOT_NCOL	/* Show on screen plot after adding neighbours, before collecting */
# define DUMP_PLOT		/* Show on screen plot after each itteration */
# undef DUMP_PLOT_RESEED	/* Show on screen plot for each re-seed point */
# undef DUMP_OPT_PLOT	/* Show on screen plot for each optimization pass */
# undef DUMP_PLOT_BEFORFIXUP	/* Show plot after reposition but before fixups */
# undef DUMP_PLOT_EACHFIXUP		/* Show each node fixup */
# define DUMP_VTX 1		/* Display the vertex locations too */
# define DUMP_PLA 1		/* Display the node planes too */
# define PERC_PLOT 0	/* Emit perceptive space plots */
# define DO_WAIT 1 		/* Wait for user key after each plot */
# undef DUMP_EPERR		/* Create .tiff of eperr */
# undef DUMP_FERR /* 10000 */	/* Create .tiff of function error >= 20 and stop. */
//# define SA_ADAPT 0.001	/* Standalone test, adaptation level */

# define SA_ADAPT -1.0		/* Standalone test, adaptation level (< 0.0, use individual) */
# define SA_DEVD_MULT 1.0	/* Delta E for each percent of device space distance */
# define SA_PERC_MULT 0.0	/* Delta E for each delta E of perceptual space distance */
# define SA_INTERP_MULT 0.0	/* Delta E for each delta E of estimated interpolation error */

#endif /* NEVER */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Overall algorithm */

#define NINSERTTRIES 100	/* Number of seedin insert tries befor failing with error() */

#define NUMTOL 1e-16	/* Numerical tollerance */
#define FTOL 1e-8		/* dnsqe function tollerance */
#define FGPMUL 5.0		/* Weighting of gamut plane error into dnsqe function */
#define COINTOL 1e-8	/* Tollerance for point cooincidence */
#define ILIMITEPS 1e-6	/* imin, imax and ilimit clip test margine */

#define FASTREJMULT1 20.5 	 /* Fast reject layer distance fudge factor */
#define FASTREJECTMULT 0.08  /* Fast reject cell skip threshold fudge factor */

#define CELLMAXEPERRFF 2.2  /* Cell worst case eperr from center fudge factor */

#define TNPAGRID 0.8		/* Target nodes per accelleration & peceptual cache grid cell */
#define TNPAGRIDMINRES 7	/* Perceptual cache grid min resolution */
#define TNPAGRIDMAXRES 33	/* Perceptual cache grid max resolution */
#undef FORCE_INCREMENTAL	/* Force incremental update after itteration */
#define FORCE_RESEED		/* Force reseed after itteration */
#define MAXTRIES 41		/* Maximum dnsq tries before giving up */
#define CACHE_PERCEPTUAL		/* Cache the perceptual lookup function */
#define USE_DISJOINT_SETMASKS		/* Reduce INDEP_SURFACE setmask size */ 

/* Sanity checks (slow) */
#undef SANITY_CHECK_SEED	/* Sanity check the selection of the bigest eperr seed vertex */
#undef SANITY_CHECK_HIT	/* Sanity check the hit detection */
#undef SANITY_CHECK_HIT_FATAL		/* throw fatal error in sanity check */
#undef SANITY_CHECK_FIXUP			/* Check that fixup was sucessful */
#undef SANITY_CHECK_FIXUP_FATAL	/* throw fatal if it wasn't */
#undef SANITY_CHECK_CLOSEST	/* Check that ofps_findclosest_xx() returns correct result */
#undef SANITY_CHECK_CLOSEST_FATAL	/* throw fatal error on sanity fail */
#undef SANITY_CHECK_CONSISTENCY		/* Check internal consistency at each itteration */
#undef SANITY_CHECK_CONSISTENCY_FATAL	/* Throw a fatal if it wasn't */

#undef SANITY_RESEED_AFTER_FIXUPS	/* Re-create veronoi from scratch after incremental update. */
#undef SANITY_CHECK_EXAUSTIVE_SEARCH_FOR_VERTEXES	/* Do very, vert slow search for all vertexes */

#define ALWAYS
#undef NEVER

#ifdef STATS
# include "conv.h"		/* System dependent convenience functions */
#endif

#if defined(DEBUG) || defined(DUMP_PLOT_SEED) || defined(DUMP_PLOT)
# include "plot.h"
# include "ui.h"
#endif

#if defined(DUMP_EPERR) || defined(DUMP_FERR)
#include "tiffio.h"
struct _vopt_cx;
static void dump_dnsqe(ofps *s, char *fname, int *nix, struct _vopt_cx *cx);
#endif

#if defined(DEBUG) || defined(DUMP_PLOT_SEED) || defined(DUMP_PLOT)
static void dump_image(ofps *s, int pcp, int dwt, int vtx, int dpla, int ferr, int noi);
#endif
#if defined(DEBUG) || defined (SANITY_CHECK_CONSISTENCY) 
static void sanity_check(ofps *s, int check_nodelists);
#endif
#if defined(DEBUG) || defined(DUMP_STRUCTURE)
static void dump_node_vtxs(ofps *s, int check_nodelists);
//static void dump_node_vtxs2(ofps *s, char *cc);
#endif

static void ofps_binit(ofps *s);
static void ofps_stats(ofps *s);
static int ofps_point2cell(ofps *s, double *v, double *p);
static void ofps_gridcoords(ofps *s, int *c, double *v, double *p);
static void ofps_add_nacc(ofps *s, node *n);
static void ofps_rem_nacc(ofps *s, node *n);
static void ofps_add_vacc(ofps *s, vtx *vx);
static void ofps_rem_vacc(ofps *s, vtx *vx);
static void ofps_add_vseed(ofps *s, vtx *vx);
static void ofps_rem_vseed(ofps *s, vtx *vx);
static void ofps_re_create_node_node_vtx_lists(ofps *s);
static void do_batch_update1(ofps *s, int fixup);
static void do_batch_update2(ofps *s, int fixup);

static node *ofps_findclosest_node(ofps *s, double *ceperr, vtx *vx);
//static vtx *ofps_findclosest_vtx(ofps *s, double *ceperr, node *nn);
static int ofps_findhit_vtxs(ofps *s, node *nn);

static char *pco(int di, int *co);
static char *ppos(int di, double *p);
static char *pcomb(int di, int *n);
static char *peperr(double eperr);
static char *psm(ofps *s, setmask *sm);

/* Check the incremental vertexes against the re-seeded vertexes */
static void save_ivertexes(ofps *s);
static int check_vertexes(ofps *s);

/* Check that no node is closer to a vertex than its parent */
static int check_vertex_closest_node(ofps *s);

/* Do an exaustive check for missing vertexes */
static void check_for_missing_vertexes(ofps *s);

/* --------------------------------------------------- */
/* Setmask manipulation functions */

#ifdef USE_DISJOINT_SETMASKS
	/* We assume the number of words is <= 1, */
	/* and we can use macros */

/* Signal this is a single word mask by using -ve no. of bits */
#define sm_init(s, nbits)  _sm_init(s, -(nbits))

#define sm_cp(s, sm_B, sm_A) \
	((sm_B)->m[0] = (sm_A)->m[0])

#define sm_or(s, sm_C, sm_A, sm_B) \
	((sm_C)->m[0] = (sm_A)->m[0] | (sm_B)->m[0])

#define sm_orand(s, sm_D, sm_A, sm_B, sm_C) \
	((sm_D)->m[0] = (sm_A)->m[0] | ((sm_B)->m[0] & (sm_C)->m[0]))

#define sm_and(s, sm_C, sm_A, sm_B) \
	((sm_C)->m[0] = (sm_A)->m[0] & (sm_B)->m[0])

#define sm_andnot(s, sm_C, sm_A, sm_B) \
	((sm_C)->m[0] = (sm_A)->m[0] & (s->lwmask ^ (sm_B)->m[0]))

#define sm_andand(s, sm_D, sm_A, sm_B, sm_C) \
	((sm_D)->m[0] = (sm_A)->m[0] & (sm_B)->m[0] & (sm_C)->m[0])

#define sm_test(s, sm_A) \
	((sm_A)->m[0] & s->lwmask)

#define sm_equal(s, sm_A, sm_B) \
	(((sm_A)->m[0] & s->lwmask) == ((sm_B)->m[0] & s->lwmask))

#define sm_andtest(s, sm_A, sm_B) \
	((sm_A)->m[0] & (sm_B)->m[0])

#define sm_andnottest(s, sm_A, sm_B) \
	((sm_A)->m[0] & (s->lwmask ^ (sm_B)->m[0]))

#define sm_vtx_vtx(s, v1, v2) \
	((v1)->vm.m[0] & (v2)->vm.m[0] & s->sc[(v1)->cmask & (v2)->cmask].a_sm.m[0])

#define sm_vtx_node(s, vx, nn) \
	((vx)->vm.m[0] & s->sc[(nn)->pmask].a_sm.m[0] & s->sc[(vx)->cmask & (nn)->pmask].a_sm.m[0])

#else

#define sm_init(s, nbits) _sm_init(s, nbits)
#define sm_cp(s, sm_B, sm_A) _sm_cp(s, sm_B, sm_A)
#define sm_or(s, sm_C, sm_A, sm_B) _sm_or(s, sm_C, sm_A, sm_B)
#define sm_orand(s, sm_D, sm_A, sm_B, sm_C) _sm_orand(s, sm_D, sm_A, sm_B, sm_C)
#define sm_and(s, sm_C, sm_A, sm_B) _sm_and(s, sm_C, sm_A, sm_B)
#define sm_andnot(s, sm_C, sm_A, sm_B) _sm_andnot(s, sm_C, sm_A, sm_B)
#define sm_andand(s, sm_D, sm_A, sm_B, sm_C) _sm_andand(s, sm_D, sm_A, sm_B, sm_C)
#define sm_test(s, sm_A) _sm_test(s, sm_A)
#define sm_equal(s, sm_A, sm_B) _sm_equal(s, sm_A, sm_B)
#define sm_andtest(s, sm_A, sm_B) _sm_andtest(s, sm_A, sm_B)
#define sm_andnottest(s, sm_A, sm_B) _sm_andnottest(s, sm_A, sm_B)
#define sm_vtx_node(s, vx, nn) _sm_vtx_node(s, vx, nn)
#define sm_vtx_vtx(s, v1, v2) _sm_vtx_vtx(s, v1, v2)

#endif

/* Compute set mask parameters */
static void _sm_init(ofps *s, int nbits) {

//printf("~1 _sm_init with %d bits\n",nbits);
	s->bpsmw = sizeof(unsigned int) * 8;

	if (nbits < 0) {	/* Macro initialisation */
#ifdef DEBUG
		printf("Disjoint sets being used\n");
#endif
		nbits = -nbits;
		if (nbits > s->bpsmw)
			error("Attempt to use macro setmasks when nbits %d > a words bits %d",nbits,s->bpsmw);
	}

	s->smbits = nbits;
	s->nsmw = (s->smbits + s->bpsmw - 1)/s->bpsmw;
	s->lwmask = ~0;
	s->lwmask >>= s->nsmw * s->bpsmw - s->smbits;			/* Number of unused bits */
	if (s->nsmw > MXSMASKW)
		error("Not enough words for %d setmask bits, got %d need %d\n",s->smbits,MXSMASKW,s->nsmw);
}
	
/* Copy a setmask */
static void _sm_cp(ofps *s, setmask *sm_B, setmask *sm_A) {
	int i;

	for (i = 0; i < s->nsmw; i++)
		sm_B->m[i] = sm_A->m[i];
}

/* Set the whole mask to zero or one */
static void sm_set(ofps *s, setmask *sm, int val) {
	int i;
	unsigned int vv = 0;

	if (val & 1)
		vv = ~0;
	for (i = 0; i < s->nsmw; i++)
		sm->m[i] = vv;
	sm->m[i-1] &= s->lwmask;
}

/* Set the given bit to zero or one */
static void sm_setbit(ofps *s, setmask *sm, int bit, int val) {
	int i;
	unsigned int vv = 0;

	if (bit > s->smbits)
		error("assert, trying to set bit %d outside setmask size %d",bit,s->smbits);
	i = bit / s->bpsmw;
	vv = 1 << bit % s->bpsmw;
	if (val & 1)
		sm->m[i] |= vv;
	else
		sm->m[i] &= ~vv;
}

/* C = A | B */
static void _sm_or(ofps *s, setmask *sm_C, setmask *sm_A, setmask *sm_B) {
	int i;

	for (i = 0; i < s->nsmw; i++)
		sm_C->m[i] = sm_A->m[i] | sm_B->m[i];
}

/* D = A | (B & C) */
static void _sm_orand(ofps *s, setmask *sm_D, setmask *sm_A, setmask *sm_B, setmask *sm_C) {
	int i;

	for (i = 0; i < s->nsmw; i++)
		sm_D->m[i] = sm_A->m[i] | (sm_B->m[i] & sm_C->m[i]);
}

/* C = A & B */
/* Return zero if result is zero */
static unsigned int _sm_and(ofps *s, setmask *sm_C, setmask *sm_A, setmask *sm_B) {
	unsigned int vv = 0;
	int i;

	for (i = 0; i < s->nsmw; i++)
		vv |= sm_C->m[i] = sm_A->m[i] & sm_B->m[i];
	return vv;
}

/* C = A & ~B */
/* Return zero if result is zero */
static unsigned int _sm_andnot(ofps *s, setmask *sm_C, setmask *sm_A, setmask *sm_B) {
	unsigned int vv = 0;
	int i;

	for (i = 0; i < s->nsmw; i++) {
		if (i < (s->nsmw-1))
			vv |= sm_C->m[i] = sm_A->m[i] & ~sm_B->m[i];
		else
			vv |= sm_C->m[i] = sm_A->m[i] & (s->lwmask ^ sm_B->m[i]);
	}
	return vv;
}

/* D = A & B & C */
/* Return zero if result is zero */
static unsigned int _sm_andand(ofps *s, setmask *sm_D, setmask *sm_A, setmask *sm_B, setmask *sm_C) {
	unsigned int vv = 0;
	int i;

	for (i = 0; i < s->nsmw; i++)
		vv |= sm_D->m[i] = sm_A->m[i] & sm_B->m[i] & sm_C->m[i];
	return vv;
}

/* Return zero if result is zero */
static unsigned int _sm_test(ofps *s, setmask *sm_A) {
	unsigned int vv = 0;
	int i;

	for (i = 0; i < s->nsmw; i++) {
		if (i < (s->nsmw-1))
			vv |= sm_A->m[i];
		else
			vv |= sm_A->m[i] & s->lwmask;
	}
	return vv;
}

/* Return nz if the two are equal */
static unsigned int _sm_equal(ofps *s, setmask *sm_A, setmask *sm_B) {
	int i;

	for (i = 0; i < s->nsmw; i++) {
		if (i < (s->nsmw-1)) {
			if (sm_A->m[i] != sm_B->m[i])
				return 0;
		} else {
			if ((sm_A->m[i] & s->lwmask) != (sm_B->m[i] & s->lwmask))
				return 0;
		}
	}
	return 1;
}

/* A & B and return zero if the result was zero. */
static unsigned int _sm_andtest(ofps *s, setmask *sm_A, setmask *sm_B) {
	unsigned int vv = 0;
	int i;

	for (i = 0; i < s->nsmw; i++)
		vv |= sm_A->m[i] & sm_B->m[i];

	return vv;
}

/* A & ~B and return zero if result is zero */
static unsigned int _sm_andnottest(ofps *s, setmask *sm_A, setmask *sm_B) {
	unsigned int vv = 0;
	int i;

	for (i = 0; i < s->nsmw; i++) {
		if (i < (s->nsmw-1))
			vv |= sm_A->m[i] & ~sm_B->m[i];
		else
			vv |= sm_A->m[i] & (s->lwmask ^ sm_B->m[i]);
	}
	return vv;
}

/* Test if two vertexes interact */
/* return nz if they do */
static unsigned int _sm_vtx_vtx(ofps *s, vtx *v1, vtx *v2) {
	unsigned int vv = 0;
	int i;

#ifdef USE_DISJOINT_SETMASKS
	/* Because the mask bits are re-used across disjoint sets, */
	/* we have to discount any intersection that occurs where */
	/* the two items are disjoint, with the exception of the full-d set. */
	for (i = 0; i < s->nsmw; i++)
		vv |= v1->vm.m[i] & v2->vm.m[i] & s->sc[v1->cmask & v2->cmask].a_sm.m[i];
#else
	for (i = 0; i < s->nsmw; i++)
		vv |= v1->vm.m[i] & v2->vm.m[i];
# endif
	return vv;
}

/* Test if a vertex and node interact */
static unsigned int _sm_vtx_node(ofps *s, vtx *vx, node *nn) {
	unsigned int vv = 0;
	int i;

	/* Because the mask bits are re-used across disjoint sets, */
	/* we have to discount any intersection that occurs where */
	/* the two items are disjoint, with the exception of the full-d set. */
#ifdef USE_DISJOINT_SETMASKS
	for (i = 0; i < s->nsmw; i++)
		vv |= vx->vm.m[i] & s->sc[nn->pmask].a_sm.m[i] & s->sc[vx->cmask & nn->pmask].a_sm.m[i];
#else
	for (i = 0; i < s->nsmw; i++)
		vv |= vx->vm.m[i] & s->sc[nn->pmask].a_sm.m[i];
# endif
	return vv;
}

/* Utility - return a string containing the mask in hex */
static char *psm(ofps *s, setmask *sm) {
	static char buf[5][200];
	static int ix = 0;
	int e, f;
	char *bp;

	if (++ix >= 5)
		ix = 0;
	bp = buf[ix];

	sprintf(bp, "0x"); bp += strlen(bp);
	for (f = 0, e = s->nsmw-1; e >= 0; e--) {
		if (f || e == 0 || sm->m[e] != 0) {
			if (f) {
				sprintf(bp, "%08x", sm->m[e]); bp += strlen(bp);
			} else {
				sprintf(bp, "%x", sm->m[e]); bp += strlen(bp);
				f = 1;
			}
		}
	}
	return buf[ix];
}

/* --------------------------------------------------- */
/* Swap the location of a node in s->n[]. This is assumed to */
/* be done _before_ a node is added to the veroni */
static void swap_nodes(ofps *s, int i, int j) {
	node *n;
	int xx;

	n = s->n[i];
	s->n[i] = s->n[j];
	s->n[j] = n;
	
	/* fix index number */
	xx = s->n[i]->ix;
	s->n[i]->ix = s->n[j]->ix;
	s->n[j]->ix = xx;

	xx = s->n[i]->ixm;
	s->n[i]->ixm = s->n[j]->ixm;
	s->n[j]->ixm = xx;
}

/* Shuffle all the nodes in the list along to put */
/* the given node at the start. */
static void move_node_to_front(ofps *s, int i) {
	node *n;
	int j;

	n = s->n[i];

	for (j = 1; j <= i; j++) 
		s->n[j] = s->n[j-1];

	s->n[0] = n;
	
	/* Fix ->ix and ixm */
	for (j = 0; j <= i; j++) {
		int bitp;
		s->n[j]->ix = i;

		bitp = 31 & (j + (j >> 4) + (j >> 8) + (j >> 12));
		s->n[j]->ixm = (1 << bitp);
	}
}

/* Randomly shuffle all the nodes */
static void shuffle_node_order(ofps *s) {
	int i;

	for (i = 0; i < s->tinp; i++) {
		swap_nodes(s, i, i_rand(0, s->tinp-1));
	}
}

/* Reverse the nodes order */
static void reverse_node_order(ofps *s) {
	int i, j;

	for (i = 0, j = s->tinp-1; i < j; i++, j--) {
		swap_nodes(s, i, j);
	}
}

/* --------------------------------------------------- */
/* Default convert the nodes device coordinates into approximate perceptual coordinates */
/* (usually overriden by caller supplied function) */
static void
default_ofps_to_percept(void *od, double *p, double *d) {
	ofps *s = (ofps *)od;
	int e;

	/* Default Do nothing - copy device to perceptual. */
	for (e = 0; e < s->di; e++) {
		double tt = d[e];
		p[e] = tt * 100.0;
	}
}

/* Filtered perceptual lookup, used for setting up rspl cache values. */
/* Input is device values, output L*a*b* like perceptual values. */
static void
filtered_ofps_to_percept(void *ss, double *p, double *d) {
	ofps *s = (ofps *)ss;
	double rad = 1.0/s->pcache_res;		/* Filter radius = grid res. */
	int fres = 2;						/* +/- 2 around center */
	int e, f;
	double pw;				/* Accumulated Weight */
	double off[MXPD];		/* Offset value */
	double out[MXPD];		/* Offset output value */
	double roff[MXPD];		/* Reflection offset value (for clip case) */
	double rout[MXPD];		/* Reflected offset output value */
	DCOUNT(co, MXPD, s->di, -fres, -fres, fres+1);

	if (rad > 0.05)			/* Don't loose too much detail */
		rad = 0.05;

//printf("filtered called with %s\n",debPdv(s->di,d)); 

	for (f = 0; f < s->di; f++) 
		p[f] = 0.0;			/* Accumulated value */
	pw = 0.0;

	DC_INIT(co);

	while (!DC_DONE(co)) {
		double tw = 1.0;
		int clip = 0;
		
//printf(" sub samp at %s\n", debPiv(s->di, co));
		for (e = 0; e < s->di; e++) {
			double w, ov;
			ov = ((double)co[e])/(fres+1) * rad;
			roff[e] = off[e] = d[e] + ov;
			if (off[e] < 0.0) {
				off[e] = 0.0;
				roff[e] = 0.0 - ov;
				clip = 1;
			}
			else if (off[e] > 1.0) {
				off[e] = 1.0;
				roff[e] = 1.0 - ov;
				clip = 1;
			}
//printf(" w[%d] = %f\n",e,(fres+1 - fabs((double)co[e]))/(fres+1)); 
			w = (fres+1 - fabs((double)co[e]))/(fres+1);
			tw *= w;
		}
//printf(" off %s wt %f\n", debPdv(s->di,off),tw);
		s->percept(s->od, out, off);

		/* For clipped case, use reflected value from reflected location */
		if (clip) {
			s->percept(s->od, rout, roff);
//printf(" roff %s\n", debPdv(s->di, roff));
//printf(" out %s\n", debPdv(s->di, out));
//printf(" rout %s\n", debPdv(s->di, rout));

			for (f = 0; f < s->di; f++)
				out[f] = 2 * out[f] - rout[f];

//printf(" eout %s\n", debPdv(s->di, out));
		}
		for (f = 0; f < s->di; f++)
			p[f] += tw * out[f];
		pw += tw;

		DC_INC(co);
	}
	for (f = 0; f < s->di; f++)
		p[f] /= pw;

//s->percept(s->od, out, d);
//printf(" u out %s\n", debPdv(s->di, out));
//printf(" f out %s\n", debPdv(s->di, p));
//printf("\n");
}

/* Cached perceptual lookup */
static void
ofps_cache_percept(void *od, double *p, double *d) {
	int e;
	co tp;
	rspl *pc = (rspl *)od;

	for (e = 0; e < pc->di; e++)
		tp.p[e] = d[e];
	pc->interp(pc, &tp);
	for (e = 0; e < pc->fdi; e++)
		p[e] = tp.v[e];
}

/* Return the distance of the device value from the device gamut */
/* This will be -ve if the point is outside */
/* If bvp is non-null, the index of the closest dim times 2 */
/* will be returned for the 0.0 boundary, dim * 2 + 1 for the 1.0 */
/* boundary, and di * 2 for the ink limit boundary. */
static double
ofps_in_dev_gamut(ofps *s, double *d, int *bvp) {
	int e, di = s->di;
	double tt;
	double dd = 100.0;				/* Worst distance outside */
	double ss = 0.0;				/* Sum of values */
	int bv = di;
	for (e = 0; e < di; e++) {
		tt = d[e] - s->imin[e];
		if (tt < dd) {
			dd = tt;
			bv = e * 2;
		}
		tt = s->imax[e] - d[e];
		if (tt < dd) {
			dd = tt;
			bv = e * 2 + 1;
		}
		ss += d[e];					/* Track sum */
	}
	ss = (s->ilimit - ss)/di;	/* Axis aligned distance to ink limit */
	tt = sqrt((double)di) * ss;	/* Diagonal distance to ink limit */
	if (tt < dd) {
		dd = tt;
		bv = di * 2;
	}
	if (bvp != NULL)
		*bvp = bv;
	return dd;
}

#ifdef NEVER	/* Allow performance trace on ofps_clip_point usage */
static int ofps_clip_point(ofps *s, double *cd, double *d);

static int ofps_clip_point1(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point2(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point3(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point4(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point5(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point6(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point7(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point8(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point9(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }
static int ofps_clip_point10(ofps *s, double *cd, double *d) {
	return ofps_clip_point(s, cd, d); }

#else	/* Production code */
#define ofps_clip_point1 ofps_clip_point
#define ofps_clip_point2 ofps_clip_point
#define ofps_clip_point3 ofps_clip_point
#define ofps_clip_point4 ofps_clip_point
#define ofps_clip_point5 ofps_clip_point
#define ofps_clip_point6 ofps_clip_point
#define ofps_clip_point7 ofps_clip_point
#define ofps_clip_point8 ofps_clip_point
#define ofps_clip_point9 ofps_clip_point
#define ofps_clip_point10 ofps_clip_point
#endif

/* Given the new intended device coordinates, */
/* clip the new position to the device gamut edge */
/* return non-zero if the point was clipped */
static int
ofps_clip_point(ofps *s, double *cd, double *d) {
	int di = s->di;
	double ss = 0.0;
	int rv = 0;

#define STEP(IX)								\
		cd[IX] = d[IX]; 						\
		if (cd[IX] < s->imin[IX]) {				\
			cd[IX] = s->imin[IX];				\
			if (cd[IX] < (s->imin[IX] - ILIMITEPS)) 	\
				rv |= 1;						\
		} else if (cd[IX] > s->imax[IX]) {		\
			cd[IX] = s->imax[IX];				\
			if (cd[IX] > (s->imax[IX] + ILIMITEPS)) 	\
				rv |= 1;						\
		}										\
		ss += cd[IX];						

	switch (di) {
		case 4:
			STEP(3)
		case 3:
			STEP(2)
		case 2:
			STEP(1)
		case 1:
			STEP(0)
	}
#undef STEP
	if (ss > s->ilimit) {
		if (ss > (s->ilimit + ILIMITEPS))
			rv |= 1;
		ss = (ss - s->ilimit)/s->di;
		switch (di) {
			case 4:
				cd[3] -= ss;
			case 3:
				cd[2] -= ss;
			case 2:
				cd[1] -= ss;
			case 1:
				cd[0] -= ss;
		}
	}
	return rv;
}

/* Given the new intended device coordinates, */
/* return non-zero if the point would be clipped. */
static int
ofps_would_clip_point(ofps *s, double *d) {
	int e;
	double ss;
	for (ss = 0.0, e = 0; e < s->di; e++) {
		if (d[e] < (s->imin[e] - ILIMITEPS))
			return 1;
		else if (d[e] > (s->imax[e] + ILIMITEPS))
			return 1;
		ss += d[e];
	}
	if (ss > (s->ilimit + ILIMITEPS))
		return 1;
	return 0;
}

/* Return a out of gamut value. */
/* 0.0 is returned if the posistion is in gamut */
static double ofps_oog(ofps *s, double *p) {
	int e, di = s->di;
	double ss, oog = 0.0;

	for (ss = 0.0, e = 0; e < di; e++) {
		if (p[e] < (s->imin[e])) {
			double tt = s->imin[e] - p[e];
			if (tt > oog) oog = tt; 
		} else if (p[e] > (s->imax[e])) {
			double tt = p[e] - s->imax[e];
			if (tt > oog) oog = tt; 
		}
		ss += p[e];
	}
	if (ss > s->ilimit) {
		double tt;
		ss = (ss - s->ilimit)/di;	/* Axis aligned distance to ink limit */
		tt = sqrt((double)di) * ss;	/* Diagonal distance to ink limit */
		if (tt > oog)
			oog = tt; 
	}
	return oog;
}

/* Unbounded perceptual lookup. */
/* return nz if it was actually clipped and extended */
static int ofps_cc_percept(ofps *s, double *v, double *p) {
	co cp;
	int clip;

	clip = ofps_clip_point(s, cp.p, p);

	if (s->pcache) {		/* In line this for speed */
		int e, di = s->di;
		
		s->pcache->interp(s->pcache, &cp);
		for (e = 0; e < di; e++)
			v[e] = cp.v[e];

	} else {
		s->percept(s->od, v, cp.p); 
	}

	/* Extend perceptual value using matrix model */
	if (clip) {
		int e, di = s->di;
		double mcv[MXPD], zv[MXPD];

#ifdef DEBUG
		if (s->pmod_init == 0)
			error("ofps_cc_percept() called before pmod has been inited");
#endif
		/* Lookup matrix mode of perceptual at clipped device */
		icxCubeInterp(s->pmod, di, di, mcv, cp.p);

		/* Compute a correction factor to add to the matrix model to */
		/* give the actual perceptual value at the clipped location */
		for (e = 0; e < di; e++)
			zv[e] = v[e] - mcv[e]; 

		/* Compute the unclipped matrix model perceptual value */
		icxCubeInterp(s->pmod, di, di, v, p);

		/* Add the correction value to it */
		for (e = 0; e < di; e++)
			v[e] += zv[e]; 
	}
	return clip;
}

/* --------------------------------------------------- */
/* Vertex alloc/free support */

/* Check if a vertex is in the cache index, */
/* and return it if it is. Return NULL otherwise */
static vtx *vtx_cache_get(ofps *s, int *nix) {
	int e, di = s->di;
	unsigned int hash;
	vtx *vx;

	hash = (unsigned int)nix[MXPD+1];		/* We assume it was put there by sort */

	for (vx = s->vch[hash]; vx != NULL; vx = vx->chn) {
		for (e = 0; e <= di; e++) { /* See if it is a match */
			if (nix[e] != vx->nix[e])
				break;
		}
		if (e > di) {	/* It is */
			return vx;
		}
	}
	return vx;
}

/* Add a vertex to the cache index */
static void vtx_cache_add(ofps *s, vtx *vv) {
	int e, di = s->di;
	unsigned int hash;

	hash = (unsigned int)vv->nix[MXPD+1];

	/* Add it to the list */
	vv->chn = s->vch[hash];
	if (s->vch[hash] != NULL)
		s->vch[hash]->pchn = &vv->chn;
	s->vch[hash] = vv;
	vv->pchn = &s->vch[hash];
}

/* Remove a vertex from the cache index */
static void vtx_cache_rem(ofps *s, vtx *vv) {
	int e, di = s->di;
	unsigned int hash;
	vtx *vx;

	hash = (unsigned int)vv->nix[MXPD+1];

	for (vx = s->vch[hash]; vx != NULL; vx = vx->chn) {
		if (vx == vv) {
			if (vx->pchn != NULL) {
				*vx->pchn = vx->chn;
				if (vx->chn != NULL)
					vx->chn->pchn = vx->pchn;
			}		
			return;
		}
	}
	/* Hmm. not in cache */
}

/* Each vtx returned gets a unique serial number */
static vtx *new_vtx(ofps *s) {
	vtx *vv;

	if (s->fvtx != NULL) {	/* re-use one we've got */
		vv = s->fvtx;
		s->fvtx = vv->link;
		memset((void *)vv, 0, sizeof(vtx));

	} else {
		if ((vv = (vtx *)calloc(sizeof(vtx), 1)) == NULL)
			error("ofps: malloc failed on new vertex");
	}

	/* Link vertex to currently used list */
	vv->link = s->uvtx;
	if (s->uvtx != NULL)
		s->uvtx->plp = &vv->link;
	s->uvtx = vv;
	vv->plp = &s->uvtx;
	vv->no = s->nxvno++;
	s->nv++;

	s->nvtxcreated++;

	vv->fuptol = NUMTOL;

	return vv;
}

/* Remove a vertx from the used list, and put it on the hidden list. */
/* (Used for making inside and outside vertexes unavailabe) */
static void remu_vtx(ofps *s, vtx *v) {
//printf("~1 remu_vtx called on no %d\n",v->no);

	/* Remove it from the used list */
	if (v->plp != NULL) {		/* If is on used list, remove it */
		*v->plp = v->link;
		if (v->link != NULL)
			v->link->plp = v->plp;
	}
	v->plp = NULL;
	v->link = NULL;

	/* Add it to the hidden list */
	v->link = s->hvtx;
	if (s->hvtx != NULL)
		s->hvtx->plp = &v->link;
	s->hvtx = v;
	v->plp = &s->hvtx;

	s->nv--;			/* Don't count hidden verts */
}

/* Remove a vertex from the cache and spatial accelleration grid, */
/* and and the used list. */
static void del_vtx1(ofps *s, vtx *vx) {
	node *nn, *nnn;

	if (vx->plp != NULL) {		/* If is on used list, remove it */
		*vx->plp = vx->link;
		if (vx->link != NULL)
			vx->link->plp = vx->plp;
		s->nv--;
	}
	
	if (vx->pfchl != NULL) {		/* If is on fixup check list, remove it */
		*vx->pfchl = vx->fchl;
		if (vx->fchl != NULL)
			vx->fchl->pfchl = vx->pfchl;
	}
	vx->pfchl = NULL;
	vx->fchl = NULL;

	if (vx->ofake == 0) {

		/* Remove it from cache */
		vtx_cache_rem(s, vx);

		/* Remove it from spatial accelleration grid */ 
		ofps_rem_vacc(s, vx);

		/* Remove it from seeding group */
		ofps_rem_vseed(s, vx);
	}

	/* Remove vertex from the s->svtxs[] list */
	/* so that it doesn't get used in fixups. */
	if (vx->psvtxs != NULL) {
		*vx->psvtxs = NULL;
		vx->psvtxs = NULL;
	}
}

/* Delete a vertex by removing it from the cache and spatial accelleration grid, */
/* and then moving it to the free list */
/* (It's assumed that it's been removed from all other */
/*  structures by the caller) */
static void del_vtx(ofps *s, vtx *vx) {

//printf("~1 del_vtx called on no %d\n",vx->no);

	/* Remove it from various lists */
	del_vtx1(s, vx);

	/* Free vertex net neighbours list */
	if (vx->nv != NULL) {
/* !!!! if we need this, we're referencing deleted vertexes !!!! */
/*		vx->nnv = vx->_nnv = 0; */
		free(vx->nv);
/*		vx->nv = NULL; */
	}

	/* Add to free list */
	vx->link = s->fvtx;
	vx->plp = NULL;
	s->fvtx = vx;

	s->nvtxdeleted++;
}

/* Add a vertex to a vertex's net */
static void vtx_add_vertex(ofps *s, vtx *vv, vtx *vx) {
	int i;

//printf("~1 Adding vertex no %d comb %s vm %s to vtx %d\n",vx->no,pcomb(s->di,vx->nix),psm(s,&vx->vm),vv->no);
//printf("~1 Adding vertex no %d to vtx no %d\n",vx->no,vv->no);
//printf("~1 Before add, no %d list is :",vv->no); for (i = 0; i < vv->nnv; i++) printf(" %d",vv->nv[i]->no); printf("\n");
	if (vv->_nnv == 0) {
		vv->_nnv = 4;		/* Initial allocation */
		if ((vv->nv = (vtx **)malloc(sizeof(vtx *) * vv->_nnv)) == NULL)
			error("ofps: malloc failed on node vertex pointers");
	} else if (vv->nnv >= vv->_nnv) {
		vv->_nnv *= 2;		/* Double allocation */
		if ((vv->nv = (vtx **)realloc(vv->nv, sizeof(vtx *) * vv->_nnv)) == NULL)
			error("ofps: realloc failed on node vertex pointers");
	}
#ifdef DEBUG
{
	int i;

	/* Check that we're not adding ourself */
	if (vx == vv) {
		printf("Adding vtx no %d comb %s to itself!\n",vx->no,pcomb(s->di,vx->nix)); fflush(stdout);
		error("Adding vtx no %d comb %s to itself!\n",vx->no,pcomb(s->di,vx->nix));
	}
	
	/* Check that the vertex is not already here */
	for (i = 0; i < vv->nnv; i++) {
		if (vx == vv->nv[i]) {
			printf("Adding vtx no %d comb %s to vtx %d when already there!\n",vx->no,pcomb(s->di,vx->nix),vv->no); fflush(stdout);
			fprintf(stderr,"Adding vtx no %d comb %s to vtx %d when already there!\n",vx->no,pcomb(s->di,vx->nix),vv->no); fflush(stdout);
//*((char *)0) = 55;
			return;
		}
	}
}
#endif /* DEBUG */

	vv->nv[vv->nnv++] = vx;

//printf("~1 After add, no %d list is :",vv->no); for (i = 0; i < vv->nnv; i++) printf(" %d",vv->nv[i]->no); printf("\n");
}

/* Delete a vertex to a vertex's net */
static void vtx_rem_vertex(ofps *s, vtx *vv, vtx *vx) {
	int i, j;

//printf("~1 Removing vertex no %d comb %s vm %s from vtx %d\n",vx->no,pcomb(s->di,vx->nix),psm(s,&vx->vm),vv->no);
//printf("~1 Before delete, no %d list is :",vv->no); for (i = 0; i < vv->nnv; i++) printf(" %d",vv->nv[i]->no); printf("\n");

	for (i = j = 0; i < vv->nnv; i++) {
		if (vv->nv[i] != vx) {
			vv->nv[j] = vv->nv[i];
			j++;
		}
	}
	vv->nnv = j;
//printf("~1 After delete, no %d list is :",vv->no); for (i = 0; i < vv->nnv; i++) printf(" %d",vv->nv[i]->no); printf("\n");
}

/* Given two vertexes, check if they are neighbours, and if they are, */
/* add them to each others neighbourhood. */
/* If fixup is set, first check that they arn't already neighbours */
/* Return nz if ther were added to each other */
static int vtx_cnd_biadd_vtx(ofps *s, vtx *vx1, vtx *vx2, int fixup) {
	int f, ff, e, di = s->di;
	int aa, bb, cc;		/* Probable hit check */
	int nnm, nmix;

	if (vx1 == vx2)
		return 0;

#ifdef NEVER	/* vertex net needs all neighbours ? */
#ifdef INDEP_SURFACE
	/* Can only have a net between them if they are visible to each other */
    if (sm_vtx_vtx(s, vx1, vx2) == 0)
		return 0;
#endif
#endif

	/* Use the nixm to quickly check if all but one parent node matches */
	aa = vx1->nix[MXPD+2];	/* nixm */
	bb = vx2->nix[MXPD+2];	/* nixm */
	if ((aa & bb) == 0 || (cc = aa & ~bb, (cc & (cc-1)) != 0)) {
		return 0;			/* It's certainly not */
	}

	/* Do an exact check of all except one node match */
	for (nnm = ff = e = 0; e <= di; e++) {
		for (f = ff; f <= di; f++) {
			if (vx1->nix[e] == vx2->nix[f]) {
				ff = f;			/* Start from here next time */
				break;
			}
			if (vx1->nix[e] > vx2->nix[f])	/* No point in looking further */
				f = di;
		}
		if (f > di) {	/* Didn't match */
			if (++nnm > 1)
				break;
			nmix = e;
		}
	}
	if (e <= di) {
		return 0;			/* No match */
	}
	
	if (nnm == 0) {
		fflush(stdout);
		error("ofps: two vertexes have the same nodes !\n"
			  "no %d at %s nix %s\nno %d at %s nix %s",
		vx1->no,ppos(di,vx1->p),pcomb(di,vx1->nix),
		vx2->no,ppos(di,vx2->p),pcomb(di,vx2->nix));
	}

	/* If fixup or not INDEP_SURFACE, check that the vertex */
	/* is not already here */
#ifndef INDEP_SURFACE
	if (fixup)
#endif
	{
		int i;
		vtx *va = vx1, *vb = vx2;

		if (vx1->nnv > vx2->nnv) {	/* Search the shortest list */
			va = vx2;
			vb = vx1;
		}
		for (i = 0; i < va->nnv; i++) {
			if (vb == va->nv[i]) {
				return 0;
			}
		}
	}

//printf("~1 Adding net between vtx no %d and no %d\n",vx1->no,vx2->no);
	/* vx2 is a neighbour, so add it to the vtx net */
	vtx_add_vertex(s, vx1, vx2);

	/* The reverse must apply too */
	vtx_add_vertex(s, vx2, vx1);

	return 1;
}

#ifdef NEVER	/* Not used */

/* Clear any veroinoi content of a vertex, but not the vertex itself. */
static void vtx_clear(ofps *s, vtx *v) {

	/* Clear the list of vertex net neighbours */
	v->nnv  = 0;
}

/* Free any allocated content of a vertex, but not the vertex itself. */
static void vtx_free(ofps *s, vtx *v) {

//printf("~1 freeing node ix %d and all contents\n",v->ix);

	/* Free up list of Vertex net neighbours */
	if (v->nv != NULL) {
		v->nnv = v->_nnv = 0;
		free(v->nv);
		v->nv = NULL;
	}
}
#endif	/* NEVER */

/* vertex binary tree support */

static int vtx_aat_cmp_eperr(const void *p1, const void *p2) {
	return ((vtx *)p1)->eperr == ((vtx *)p2)->eperr ? 0 :
	                           (((vtx *)p1)->eperr < ((vtx *)p2)->eperr ? -1 : 1);
}

static int vtx_aat_cmp_eserr(const void *p1, const void *p2) {
	return ((vtx *)p1)->eserr == ((vtx *)p2)->eserr ? 0 :
	                           (((vtx *)p1)->eserr < ((vtx *)p2)->eserr ? -1 : 1);
}


/* --------------------------------------------------- */
/* Midpoint alloc/free support */

/* Each mid returned gets a unique serial number */
/* and a refc or 0 */
static mid *new_mid(ofps *s) {
	mid *p;

	if (s->fmid != NULL) {	/* re-use one we've got */
		p = s->fmid;
		s->fmid = p->link;
		memset((void *)p, 0, sizeof(mid));

	} else {
		if ((p = (mid *)calloc(sizeof(mid), 1)) == NULL)
			error("ofps: malloc failed on new midpoint");
	}

	/* Link midpoint to currently used list */
	p->link = s->umid;
	if (s->umid != NULL)
		s->umid->plp = &p->link;
	s->umid = p;
	p->plp = &s->umid;
	p->no = s->nxmno++;

	return p;
}

/* Decrement reference count, and midpoint to the free list */
static void del_mid(ofps *s, mid *p) {
//printf("~1 del_mid called on no %d, refc = %d\n",p->no,p->refc);
	if (--p->refc <= 0) {

		if (p->plp != NULL) {		/* If is on used list, remove it */
			*p->plp = p->link;
			if (p->link != NULL)
				p->link->plp = p->plp;
		}

		p->link = s->fmid;		/* Add to free list */
		p->plp = NULL;
		s->fmid = p;
		p->refc = 0;
	}
}

/* --------------------------------------------------- */
/* Node basic support functions */

/* Clear any veroinoi content of a node, but not the node itself. */
static void node_clear(ofps *s, node *p) {

	/* Clear the list of Voronoi verticies */
	p->nvv = 0;

	/* Clear any midpoints and nodes */
	while (p->nvn > 0) {
		if (p->mm[--p->nvn] != NULL)
			del_mid(s, p->mm[p->nvn]);
	}
}

/* Free any allocated content of a node, but not the node itself. */
static void node_free(ofps *s, node *p) {

//printf("~1 freeing node ix %d and all contents\n",p->ix);

	/* Free up list of Voronoi verticies */
	if (p->vv != NULL) {
		free(p->vv);
		p->vv = NULL;
		p->nvv = p->_nvv = 0;
	}

	/* Free up list of voronoi node indexes */
	if (p->vn != NULL) {
		while (p->nvn > 0) {
			if (p->mm[--p->nvn] != NULL)
				del_mid(s, p->mm[p->nvn]);
		}
		p->nvn = p->_nvn = 0;
		free(p->vn);
		free(p->mm);
		p->vn = NULL;
		p->mm = NULL;
	}

	p->nsp = 0;		/* No list of surface planes */
}

/* Add a vertex to the node vertex list */
static void node_add_vertex(ofps *s, node *pp, vtx *vx) {

//printf("~1 Adding vertex no %d comb %s vm %s to node %d\n",vx->no,pcomb(s->di,vx->nix),psm(s,&vx->vm),pp->ix);
#ifdef DEBUG
	if (vx->del)
		warning("!!!!! adding vertex no %d with delete flag set!!!",vx->no);
#endif

	if (pp->_nvv == 0) {
		pp->_nvv = 4;		/* Initial allocation */
		if ((pp->vv = (vtx **)malloc(sizeof(vtx *) * pp->_nvv)) == NULL)
			error("ofps: malloc failed on node vertex pointers");
	} else if (pp->nvv >= pp->_nvv) {
		pp->_nvv *= 2;		/* Double allocation */
		if ((pp->vv = (vtx **)realloc(pp->vv, sizeof(vtx *) * pp->_nvv)) == NULL)
			error("ofps: realloc failed on node vertex pointers");
	}

#ifdef DEBUG
{
	int i;

	/* Check that the vertex is not already here */
	for (i = 0; i < pp->nvv; i++) {
		if (vx == pp->vv[i]) {
			printf("Adding vtx no %d comb %s when already there!\n",vx->no,pcomb(s->di,vx->nix)); fflush(stdout);
			error("Adding vtx no %d comb %s when already there!",vx->no,pcomb(s->di,vx->nix)); fflush(stdout);
		}
	}
}
#endif /* DEBUG */

	pp->vv[pp->nvv++] = vx;

#ifdef NEVER
#ifdef DEBUG
	{
		int e, di = s->di;
		printf("~1 +++ Node ix %d add vtx no %d pos %s err %f @ %d:",pp->ix,vx->no,ppos(di,vx->p),vx->eperr,pp->nvv);
		for (e = 0; e < pp->nvv; e++)
			printf("%d ",pp->vv[e]->no);
		printf("\n");
	}
#endif
#endif
}

/* Remove a vertex from the node vertex list */
static void node_rem_vertex(ofps *s, node *pp, vtx *vx) {
	int i, j;

//printf("~1 Removing vertex no %d comb %s vm %s from node ix %d\n",vx->no,pcomb(s->di,vx->nix),psm(s,&vx->vm),pp->ix);
//printf("~1 Before delete, no %d list is :",vv->no); for (i = 0; i < vv->nnv; i++) printf(" %d",vv->nv[i]->no); printf("\n");

	for (i = j = 0; i < pp->nvv; i++) {
		if (pp->vv[i] != vx) {
			pp->vv[j] = pp->vv[i];
			j++;
		}
	}
	pp->nvv = j;
}

/* Add a node index to the node */
static void node_add_nix(ofps *s, node *pp, int ix) {

	if (pp->_nvn == 0) {
		pp->_nvn = 4;		/* Initial allocation */
		if ((pp->vn = (int *)malloc(sizeof(int) * pp->_nvn)) == NULL)
			error("ofps: malloc failed on node index list");
		if ((pp->mm = (mid **)malloc(sizeof(mid *) * pp->_nvn)) == NULL)
			error("ofps: malloc failed on midpoint pointer list");
	} else if (pp->nvn >= pp->_nvn) {
		pp->_nvn *= 2;		/* Double allocation */
		if ((pp->vn = (int *)realloc(pp->vn, sizeof(int) * pp->_nvn)) == NULL)
			error("ofps: realloc failed on node index list");
		if ((pp->mm = (mid **)realloc(pp->mm, sizeof(mid *) * pp->_nvn)) == NULL)
			error("ofps: realloc failed on midpoint pointer list");
	}
	pp->vn[pp->nvn] = ix;
	pp->mm[pp->nvn++] = NULL;

#ifdef NEVER
#ifdef DEBUG
	{
		int e, di = s->di;
		printf("~1 +++ Node ix %d add node ix %d at %s @ %d: ",pp->ix,ix,ppos(di,pp->p),pp->nvn);
		for (e = 0; e < pp->nvn; e++)
			printf("%d ",pp->vn[e]);
		printf("\n");
	}
#endif
#endif
}

/* Recompute a nodes neighborhood nodes. */
/* (Invalidates and deletes any midpoints) */
static void node_recomp_nvn(
	ofps *s,
	node *pp
) {
	int e, di = s->di;
	int i, j, k;

#ifdef DEBUG
	printf("node_recomp_nvn for node ix %d\n",pp->ix);
#endif
	/* Clear any midpoints and nodes */
	while (pp->nvn > 0) {
		if (pp->mm[--pp->nvn] != NULL)
			del_mid(s, pp->mm[pp->nvn]);
	}
	s->nvnflag++;					/* Make sure each node is only added once */
	pp->nvnflag = s->nvnflag;		/* Don't put self in list */

	for (i = 0; i < pp->nvv; i++) {		/* For each vertex */
		double rads;
		vtx *vv = pp->vv[i];

		for (j = 0; j <= di; j++) {		/* For each node in vertex */
			int ix = vv->nix[j];
			node *ap = s->n[ix];

			if (ap->nvnflag == s->nvnflag)
				continue;				/* Already done that node */

			node_add_nix(s, pp, ix);
			ap->nvnflag = s->nvnflag;	/* Don't worry about it again */
		}
	}
}

/* Sort a vertex node index array of di+1 nodes, */
/* and add a hash at MXDP+1, and nixm at MXDP+2 */
/* This is to speed up searching for a match */
/* Sort largest to smallest (so fake gamut nodes are last) */
static void sort_nix(ofps *s, int *nix) {
	int i, j, t;
	int di = s->di;		/* There are di+1 nodes */
	unsigned int hash = 0;
	int nixm = 0;

	/* Do a really simple exchange sort */
	for (i = 0; i < di; i++) {
		for (j = i+1; j <= di; j++) {
			if (nix[i] < nix[j]) {
				t = nix[j]; 
				nix[j] = nix[i];
				nix[i] = t;
			}
		}
	}

	/* And then compute the hash and nixm */
	for (i = 0; i <= di; i++) {
		int bitp, ix = nix[i];
		hash = hash * 17 + nix[i];
		bitp = 31 & (ix + (ix >> 4) + (ix >> 8) + (ix >> 12));
		nixm |= (1 << bitp);
	}
	hash %= VTXCHSIZE;

	nix[MXPD+1] = (int)hash;
	nix[MXPD+2] = nixm;
}

/* Check if the given locate is on the gamut boundary surface, */
/* and return the corresponding plane mask */
static unsigned int check_pos_gsurf(ofps *s, double *p) {
	int i, e, di = s->di;
	unsigned int pmask = 0;

	/* For all the gamut boundary planes */
	for (i = 0; i < s->nbp; i++) {
		pleq *vp = &s->gpeqs[i];
		double v;

		for (v = vp->pe[di], e = 0; e < di; e++) 	/* Compute relation to plane equation */
			v += vp->pe[e] * p[e];
		/* See if this location close to, or outside boundary plane */
		if (v > -s->surftol)
			pmask |= (1 << i); 
	}
#ifdef MAXINDEP_2D
	if (s->sc[pmask].valid == 0)
		pmask = 0;
#endif 
	return pmask;
}

/* Check if the given node is on the gamut boundary surface, */
/* and record the number of surfaces it is on. */
/* Return nz if the node is on one or more gamut boundaries. */
/* Set the state of the node clip flag too. */
static int det_node_gsurf(ofps *s, node *n, double *p) {
	int i, e, di = s->di;
	double ss;

	n->nsp = 0;
	n->pmask = 0;

	/* For all the gamut boundary planes */
	for (i = 0; i < s->nbp; i++) {
		pleq *vp = &s->gpeqs[i];
		double v;

		for (v = vp->pe[di], e = 0; e < di; e++) 	/* Compute relation to plane equation */
			v += vp->pe[e] * p[e];

		/* See if this location close to, or outside boundary plane */
		if (v > -s->surftol) {
			/* Add pointer to plane it falls on */
			n->sp[n->nsp++] = vp;
			n->pmask |= (1 << i); 
			if (n->nsp > MXPD+1)
				error("Assert in ofps det_node_gsurf : nsp %d > MXPD +1 %d",n->nsp, MXPD+1);
		}
	}

#ifdef MAXINDEP_2D
	if (s->sc[n->pmask].valid == 0) {
		n->pmask = 0;
		n->nsp = 0;
	}
#endif 
//printf("~1 node pmask = 0x%x\n",n->pmask);

	return (n->nsp > 0);
}

/* Compute a cmask from an nix */
static unsigned int comp_cmask(ofps *s, int *nix) {
	unsigned int smask = 0, cmask = ~0;
	int i, e, di = s->di;

	/* The composition mask indicates all the common surface planes */
	/* that a vertexes parent nodes lie on. Given this, one expects */
	/* the resulting location to be the same of within this. */
	for (e = 0; e <= di; e++) {
		int ix = nix[e];
		if (ix < 0 && ix >= (-s->nbp)) {	/* If fake surface node */
			smask |= 1 << (-ix-1);
		} else if (ix >= 0) {			/* If real node */
			cmask &= s->n[ix]->pmask;
		}
	}
	if (smask != 0)
		cmask &= smask;

#ifdef MAXINDEP_2D
	if (s->sc[cmask].valid == 0)
		cmask = 0;
#endif 

	return cmask;
}

/* Check if the given vertex is on the gamut boundary surface, */
/* and record the number of surfaces it is on. */
/* Also compute its cmask based on its parent nodes. */
/* Set the state of the clip flag too. */
static void det_vtx_gsurf(ofps *s, vtx *vx) {
	int i, e, di = s->di;
	unsigned int smask = 0;

	vx->nsp = 0;

	/* The composition mask indicates all the common surface planes */
	/* that a vertexes parent nodes lie on. Given this, one expects */
	/* the resulting location to be the same of within this. */
	vx->cmask = ~0;
	for (e = 0; e <= di; e++) {
		int ix = vx->nix[e];
		if (ix < 0 && ix >= -s->nbp) {	/* If fake surface node */
			smask |= 1 << (-ix-1);
		} else if (ix >= 0) {
			vx->cmask &= s->n[ix]->pmask;
		}
	}
	if (smask != 0)
		vx->cmask &= smask;

#ifdef MAXINDEP_2D
	if (s->sc[vx->cmask].valid == 0)
		vx->cmask = 0;
#endif 

	/* For all the gamut boundary planes */
	for (i = 0; i < s->nbp; i++) {
		pleq *vp = &s->gpeqs[i];
		double v;

		for (v = vp->pe[di], e = 0; e < di; e++) 	/* Compute relation to plane equation */
			v += vp->pe[e] * vx->p[e];

		/* See if this location close to, or outside boundary plane */
		if (v > -s->surftol) {
			vx->sp[vx->nsp++] = vp;
			vx->pmask |= (1 << i); 
			if (vx->nsp > di+1)
				error("Assert in ofps det_vtx_gsurf : nsp %d > di+1 %d",vx->nsp, di+1);
		}
	}
#ifdef MAXINDEP_2D
	if (s->sc[vx->pmask].valid == 0) {
		vx->pmask = 0;
		vx->nsp = 0;
	}
#endif 

//printf("~1 vertex pmask = 0x%x, cmask = 0x%x\n",vx->pmask,vx->cmask);
}

/* Given a device position and a list of surface planes, */
/* move the position to lie on the closest location */
/* on those planes. */
static void confineto_gsurf(ofps *s, double *p, pleq **psp, int nsp) {

	if (nsp > 0) {		/* It's a surface point, so keep it on the surface */
		int i, j, e, di = s->di;
		double nn, np, q;

		/* Special case the common situation for speed. */
		if (nsp == 1) {
			pleq *sp = psp[0];

			/* Compute the dot product of the plane equation normal */
			for (nn = 0.0, e = 0; e < di; e++)
				nn += sp->pe[e] * sp->pe[e]; 

			/* Compute the dot product of the plane equation and the point */
			for (np = 0.0, e = 0; e < di; e++)
				np += sp->pe[e] * p[e]; 

			/* Compute the parameter */
			q = (sp->pe[di] + np)/nn;

			/* Compute the closest point */
			for (e = 0; e < di; e++)
				p[e] -= q * sp->pe[e];

		/* General case using matrix solution. */
		/* We compute the proportion of each plane normal vector to add to point */
		/* to map it onto all planes simultaniously (ie. to map the point to */
		/* the intersection of all the planes). */
		} else if (nsp > 1) {
			double **ta, *TTA[MXPD + 1], TA[MXPD+1][MXPD + 1];
			double *tb, TB[MXPD + 1];

			for (e = 0; e < nsp; e++)
				TTA[e] = TA[e];
			ta = TTA;
			tb = TB;

			/* For each combination of planes */
			for (i = 0; i < nsp; i++) {
				pleq *spi = psp[i];
				for (j = i; j < nsp; j++) {
					pleq *spj = psp[j];
					double vv;

					/* Compute dot product of the two normals */
					for (vv = 0.0, e = 0; e < di; e++)
						vv += spi->pe[e] * spj->pe[e];
					ta[j][i] = ta[i][j] = vv;		/* Use symetry too */
				}

				/* Compute right hand side */
				for (tb[i] = 0.0, e = 0; e < di; e++)
					tb[i] += spi->pe[e] * p[e];		/* Dot prod of plane normal and point */
				tb[i] += spi->pe[di];				/* plus plane constant */
			}
			/* Solve the simultaneous linear equations A.x = B */
			/* Return 1 if the matrix is singular, 0 if OK */
			if (solve_se(ta, tb, nsp) == 0) { 
				/* Compute the closest point */
				for (i = 0; i < nsp; i++) {
					pleq *spi = psp[i];
					for (e = 0; e < di; e++)
						p[e] -= tb[i] * spi->pe[e];
				}
			}
		}
		/* The mapping may leave it out of gamut */
		ofps_clip_point2(s, p, p);
	}
}

/* Given a device position and a list of surface planes, */
/* check that the point lies on all the planes. */
/* Return NZ if it does, Z if it doesn't */
static int checkon_gsurf(ofps *s, double *p, pleq **psp, int nsp) {
	int i, e, di = s->di;
	double nn, np, q;

	if (nsp == 0) 
		return 1;

	for (i = 0; i < nsp; i++) {
		pleq *vp = psp[i];
		double v;

		for (v = vp->pe[di], e = 0; e < di; e++) 	/* Compute relation to plane equation */
			v += vp->pe[e] * p[e];

		/* See if this location close to, or outside boundary plane */
		if (fabs(v) > s->surftol) {
			return 0;
		}
	}
	return 1;
}

/* Compute the estimated positioning error given two locations. */
/* [ This seems to be the critical inner loop in regard to */
/*   overall speed. The dominant callers are dnsq_solver() 10%, */
/*   followed by add_node2voronoi() 5%, others <= 1% ] */

#ifdef NEVER	/* Allow performance trace on eperr usage */
static double ofps_comp_eperr(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp);
static double ofps_comp_eperr1(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); }
static double ofps_comp_eperr2(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); }
static double ofps_comp_eperr3(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); }
static double ofps_comp_eperr4(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); }
static double ofps_comp_eperr5(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); }
static double ofps_comp_eperr6(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); }
static double ofps_comp_eperr7(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); }
static double ofps_comp_eperr8(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); } 
static double ofps_comp_eperr9(ofps *s, double *pddist, double *v, double *p, double *nv, double *np, int nsp) {
	return ofps_comp_eperr(s, pddist, v, p, nv, np, nsp); }
#else	/* Production code */
#define ofps_comp_eperr1 ofps_comp_eperr
#define ofps_comp_eperr2 ofps_comp_eperr
#define ofps_comp_eperr3 ofps_comp_eperr
#define ofps_comp_eperr4 ofps_comp_eperr
#define ofps_comp_eperr5 ofps_comp_eperr
#define ofps_comp_eperr6 ofps_comp_eperr
#define ofps_comp_eperr7 ofps_comp_eperr
#define ofps_comp_eperr8 ofps_comp_eperr
#define ofps_comp_eperr9 ofps_comp_eperr
#endif

static double ofps_comp_eperr(
	ofps *s,
	double *pddist,	/* If not NULL, return the device distance */
	double *v,		/* Device perceptual value */ 
	double *p,		/* Device sample location to be evaluated */
	double *nv,		/* Other perceptual value */
	double *np,		/* Other sample location value */
	int nsp			/* Number of surface planes */
) {
	int ii, e, f, di = s->di;
	int isc;
	double tt, ddist, pdist;
	double eperr;

	/* Uncertaintly error computed from device and perceptual distance */

#ifndef NEVER	/* unrolled code */

#if MXPD > 4
# error "ofps.c: Need to expand switch  code for MXPD > 4"
#endif
	/* Unrole the loop */
	pdist = ddist = 0.0;
	switch (di) {
		case 4:
			tt = (p[3] - np[3]); /* Device distance */
			ddist += tt * tt;
			tt = (v[3] - nv[3]); /* Perceptual distance */
			pdist += tt * tt;
		case 3:
			tt = (p[2] - np[2]); /* Device distance */
			ddist += tt * tt;
			tt = (v[2] - nv[2]); /* Perceptual distance */
			pdist += tt * tt;
		case 2:
			tt = (p[1] - np[1]); /* Device distance */
			ddist += tt * tt;
			tt = (v[1] - nv[1]); /* Perceptual distance */
			pdist += tt * tt;
		case 1:
			tt = (p[0] - np[0]); /* Device distance */
			ddist += tt * tt;
			tt = (v[0] - nv[0]); /* Perceptual distance */
			pdist += tt * tt;
	}
#else
	/* General code */
	for (pdist = ddist = 0.0, e = 0; e < di; e++) {
		
		/* Compute the device distance */
		tt = (p[e] - np[e]);
		ddist += tt * tt;

		/* Compute the perceptual distance */
		tt = (v[e] - nv[e]);
		pdist += tt * tt;
	}
#endif

	if (pddist != NULL)
		*pddist = ddist;

	ddist *= 100.0 * 100.0;

//printf("~1 Device distance = %f, dev error = %f\n",ddist,s->devd_wght * ddist);

	ddist = sqrt(ddist);
	pdist = sqrt(pdist);

	eperr = s->devd_wght * ddist + s->perc_wght * pdist;

#ifdef GAMUT_EDGE_FUDGE
	/* Fudge factor to prevent gap at gamut boundaries */
	if (nsp > 0) {
//		int nn;
//		for (nn = 0; nn < nsp; nn++)
			eperr *= GAMUT_EDGE_FUDGE;
	}
#endif /* GAMUT_EDGE_FUDGE */

//printf("~1 Percept distance = %f, perc error = %f\n",pdist,s->perc_wght * pdist);
	return eperr;
}

/* Compute the per node estimated position and interpolation errors */
/* given a location and a list of up to di+1 neighborhood measurement nodes. */
static void ofps_pn_eperr(
	ofps *s,
	double *ce,		/* return the curvature/interpolation error for each node (may be NULL) */
	double *ee,		/* return the uncertaintly error for each node */
	double *sv,		/* Perceptual value if known, othewise NULL */
	double *sp,		/* Device sample location to be evaluated */
	node **nds,		/* Array of pointers to measurement nodes */
	int nnds		/* Number of measurement nodes (>= 1) */
) {
	int ii, e, di = s->di;
	node *np;
	double _sv[MXPD];			/* Sample perceptual value */
	double iv[MXPD];			/* Interpolated perceptual value */

	/* Lookup perceptual value at sample point location */
	if (sv == NULL) {
		sv = _sv;
		ofps_cc_percept(s, sv, sp);
	}

	/* Uncertaintly error computed from device and perceptual distance */
	for (ii = 0; ii < nnds; ii++)
		ee[ii] = ofps_comp_eperr1(s, NULL, sv, sp, nds[ii]->v, nds[ii]->p, nds[ii]->nsp);

	if (ce == NULL)
		return;

	/* This could be made more efficient by only computing it for every */
	/* vertex during the initial seeding, and then on subsequent */
	/* passes only computing it once the re-seed/fixups are done. */
	if (s->curv_wght != 0.0) {	/* Don't waste the time unless it's used */

		/* Compute an error estimate that's related to curvature */
		for (ii = 0; ii < nnds; ii++) {
			double cp[MXPD];		/* Midway points location */
			double civ[MXPD];		/* Midway points interpolated perceptual value */
			double cv[MXPD];		/* Midway points actual perceptual value */

			/* Compute a point midway between the sample location and the node */
			for (e = 0; e < di; e++) {
				cp[e] = 0.5 * (sp[e] + nds[ii]->p[e]);
				civ[e] = 0.5 * (sv[e] + nds[ii]->v[e]);
			}

			/* Look the actual perceptual value */
			/* (Computing this for each vertex consumes about 8% of execution time) */
			ofps_cc_percept(s, cv, cp);

			/* Compute the difference between the interpolated and actual perceptual values */
			for (ce[ii] = 0.0, e = 0; e < di; e++) {
				double tt;
				tt = civ[e] - cv[e];
				ce[ii] += tt * tt;
			}
			ce[ii] = s->curv_wght * sqrt(ce[ii]);
#ifdef GAMUT_EDGE_FUDGE
			/* Fudge factor to prevent gap at gamut boundaries */
			if (nds[ii]->nsp > 0) {
//				int nn;
//				for (nn = 0; nn < nds[ii]->nsp; nn++)
					ce[ii] *= GAMUT_EDGE_FUDGE;
			}
#endif /* GAMUT_EDGE_FUDGE */

		}
	} else {
		for (ii = 0; ii < nnds; ii++)
			ce[ii] = 0.0;
	}
}

/* Compute the estimated position error given the results */
/* of ofps_pn_eperr() */ 
static double ofps_eperr2(
	double *ee,		/* Uncertainty error for each node */
	int nnds		/* Number of measurement nodes */
) {
	double eperr;
	int ii;

	for (eperr = 1e80, ii = 0; ii < nnds; ii++) {
		if (ee[ii] < eperr)
			eperr = ee[ii];
	}

//printf("~1 ofps_eperr returning %f\n",eperr);
	return eperr;
}

/* Compute the estimated sampling error of a location given */
/* the results of ofps_pn_eperr() */ 
static double ofps_eserr2(
	double *ce,		/* The estimated curvature/interpolation error for each node */
	double *ee,		/* Uncertainty error for each node */
	int nnds		/* Number of measurement nodes */
) {
	int ii;
	double eserr;
	double mxce;

#ifdef NEVER
	/* We assume errors are inverse probabilities, */
	/* and sum inverse squares. */
	for (mxce = 0.0, eserr = 0.0, ii = 0; ii < nnds; ii++) {
		double tt;

		tt = ee[ii] * ee[ii];
		if (tt > NUMTOL)
			eserr += 1.0/tt;
		else {					/* One error is close to zero */
			eserr = 0.0;
			break;
		}
		if (ce[ii] > mxce)
			mxce = ce[ii];
	}
	if (ii >= nnds)
		eserr = 1.0/sqrt(eserr);
	eserr += mxce;

#else	/* This seems best ? */
	/* Return nearest neighbor error metric for the moment, */
	/* with the ce being the maximum of the surrounders. */
	for (mxce = 0.0, eserr = 1e80, ii = 0; ii < nnds; ii++) {
		if (ee[ii] < eserr)
			eserr = ee[ii];
		if (ce[ii] > mxce)
			mxce = ce[ii];
	}
	eserr += mxce;
#endif

//printf("~1 ofps_eserr returning %f\n",eserr);
	return eserr;
}

/* - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - */
/* Finding vertex location code using dnsqe() */

typedef struct {
	double p[MXPD];
} loc;

/* Context for callback */
struct _vopt_cx {
	ofps *s;
	node *nds[MXPD+1];	/* List of real nodes */
	int nn;				/* Number of real nodes */
	int on;				/* Index of odd node */ 

	pleq *sp[MXPD+1];	/* List of touched gamut surface planes */
	int nsp;			/* Number of touched gamut surface planes */

	double srad;		/* Search radius used */
	double stp[MXPD];	/* Starting point used */

#ifdef DUMP_FERR
	/* Debug: */
	int debug;			/* nz to trace search path */
	loc *clist;			/* List of points sampled */
	int _nl;			/* Allocated size */
	int nl;				/* Number of points */
#endif

}; typedef struct _vopt_cx vopt_cx;

/* calculate the functions at x[] */
int dnsq_solver(	/* Return < 0 on abort */
	void *fdata,	/* Opaque data pointer */
	int n,			/* Dimenstionality */
	double *x,		/* Multivariate input values */
	double *fvec,	/* Multivariate output values */
	int iflag		/* Flag set to 0 to trigger debug output */
) {
	vopt_cx *cx = (vopt_cx *)fdata;
	ofps *s = cx->s;
	int k, e, di = s->di;
	int nn_1 = cx->nn-1;
	double sv[MXPD];
	double cee[MXPD+1], teperr;

#ifdef DUMP_FERR
	/* record the points we visited */
	if (cx->debug) {
		if (cx->nl >= cx->_nl) {
			cx->_nl = 2 * cx->_nl + 5;
			if ((cx->clist = (loc *)realloc(cx->clist, sizeof(loc) * cx->_nl)) == NULL)
				error("ofps: malloc failed on debug location array %d", cx->_nl);
		}
		for (e = 0; e < di; e++)
			cx->clist[cx->nl].p[e] = x[e];
		cx->nl++;
	}
#endif
//printf("~1 dnsq_solver got %d nodes and %d planes\n",cx->nn,cx->nsp);

	/* Get eperr at each real node */
	ofps_cc_percept(s, sv, x);	/* We have to compute it */
	for (k = 0; k < cx->nn; k++) {

		cee[k] = ofps_comp_eperr2(s, NULL, sv, x, cx->nds[k]->v, cx->nds[k]->p, cx->nds[k]->nsp);
	}

//fprintf(stderr,"~1 maxeperr = %f\n",cmax);

//printf("~1 error =");
//for (k = 0; k < cx->nn; k++)
//	printf(" %f",cee[k]);
//printf("\n");

	/* We need to create nn-1 output values from nn eperr's */
	/* cx->on is the odd one out. */
	/* Difference to average (best) */
	for (teperr = 0.0, k = 0; k < cx->nn; k++)
		teperr += cee[k];	
	teperr /= (double)cx->nn;

	for (k = e = 0; k < cx->nn; k++) {
		if (k != cx->on) {
			fvec[e++] = (teperr - cee[k]);
		}
	}

	/* Compute plane errors */
	for (k = 0; k < cx->nsp; k++) {
		pleq *pl = cx->sp[k];
		double v;

		for (v = pl->pe[di], e = 0; e < di; e++) 	/* Compute relation to plane equation */
			v += pl->pe[e] * x[e];
		fvec[nn_1 + k] = FGPMUL * v;
	}

	s->funccount++;

//for (k = 0; k < nn_1; k++)
//printf("~1 fvec[%d] = %f\n",k,fvec[k]);
//printf("dnsq_solver returning %s from %s\n",ppos(di,fvec),ppos(di,x));
//fprintf(stderr,"dnsq_solver returning %s from %s\n",ppos(di,fvec),ppos(di,x));

	return 0;
}

/* Locate a vertex position that has the eperr from all the real nodes */
/* being equal. Set eperr, eserr and subjective value v[] too. */ 
/* vv->ceperr contains the current eperr that must be bettered. */
/* Return 0 if suceeded, 1 if best result is out of tollerance, 2 if failed. */
static int position_vtx(
	ofps *s,
	nodecomb *vv,		/* Return the location and its error */
	int startex,		/* nz if current position is to be used as initial start point */
	int repos,			/* nz after an itteration and we expect out of gamut */
	int fixup			/* nz if doing fixups after itteration and expect out of gamut ??? */
) {
	int e, di = s->di;
	int k, ii;
	double tw;
	double atp[MXPD], mct[MXPD];	/* Average node position, middle of closest 2 nodes */	
	double osp[MXPD];				/* Original start position, start position */	
	double cdist, fdist;			/* Closest/furthest two points distance apart */
	double bsrad;					/* Basic search radius, search radius */
	double ftol = FTOL;				/* Final eperr tolerance */
	int tfev = 0, tcalls = 0;		/* Track average successful fevs */
	int maxfev = 50;				/* Maximum function evaluations */
	int tries, notries = MAXTRIES;	/* Point to give up on. Giving up will error. */
	vopt_cx cx, pcx;				/* dnsq context + previous context */

#ifdef DEBUG
	printf("Position_vtx called for comb %s\n",pcomb(di,vv->nix));
#endif

	s->positions++;
	s->sob->reset(s->sob);

#ifdef DUMP_FERR
	cx.debug = 0;
#endif

	/* Setup for dnsq to optimize for equal eperr */
	cx.s = s;

	/* Pointers to real nodes. Although we allow for the */
	/* fake inner/outer nodes, eperr() will fail them later. */
	for (ii = e = 0; e <= di; e++) {
		if (vv->nix[e] >= 0 || vv->nix[e] < -s->nbp)
			cx.nds[ii++] = s->n[vv->nix[e]];
	}
	cx.nn = ii;

	if (ii == 0) {
		fflush(stdout);
		error("ofps: unexpectedely got no real nodes in vertex position %s",pcomb(di,vv->nix));
	}

#ifdef DEBUG
		printf("%d real nodes\n",ii);
		for (k = 0; k < ii; k++)
			printf("Node ix %d at %s\n",cx.nds[k]->ix,ppos(di,cx.nds[k]->p));
#endif

	/* Setup gamut suface planes */
	cx.nsp = 0;
	if (ii < (di+1)) {
		/* Create list of pointers to the gamut surface planes involved */
		for (e = 0, k = ii; k <= di; e++, k++) {
#ifdef DEBUG
			printf("Adding plane for node ix %d\n",vv->nix[k]);
#endif
			cx.sp[e] = &s->gpeqs[-1 - vv->nix[k]];
		}
		cx.nsp = e;
	}

	/* If there is only one node, map it to the planes */
	/* and we're done. */
	if (ii == 1) {
		double ee[MXPD+1];

		for (e = 0; e < di; e++)
			vv->p[e] = cx.nds[0]->p[e];
		
		confineto_gsurf(s, vv->p, cx.sp, cx.nsp);

		if (checkon_gsurf(s, vv->p, cx.sp, cx.nsp) == 0) { 
#ifdef DEBUG
			printf("Single node comb %s failed to confine to gamut surface\n",pcomb(di,vv->nix));
#endif
			return 2;
		}

		/* Compute perceptual (can't clip because of confine) */
		s->percept(s->od, vv->v, vv->p);

		/* Compute the eperr's for each node. */
		ofps_pn_eperr(s, vv->ce, ee, vv->v, vv->p, cx.nds, cx.nn);

		/* Compute errors at returned location */
		vv->eperr = ofps_eperr2(ee, cx.nn);
		vv->eserr = ofps_eserr2(vv->ce, ee, cx.nn);

#ifdef DEBUG
		printf("Single node, returning comb %s opt pos = %s, eperr = %f, eserr = %f\n",pcomb(di,vv->nix),ppos(di,vv->p),vv->eperr,vv->eserr);
#endif
		return 0;
	}
	{
		/* Compute average of real nodes */
		for (e = 0; e < di; e++)
			atp[e] = 0.0;
		for (tw = 0.0, k = 0; k < ii; k++) {
			double w = 1.0;

			for (e = 0; e < di; e++)
				atp[e] += cx.nds[k]->p[e];
			tw += w;
		}
		for (e = 0; e < di; e++)
			atp[e] /= tw;

#ifdef DEBUG
		printf("Average of %d real node start pos = %s\n",ii,ppos(di,atp));
#endif
	}

	/* Locate the closest and furthest two nodes */
	{
		double ceperr = 1e200;
		int i, j, bi = 0, bj = 0;

		/* Find the two vectors that have the closest eperr. Brute force search */
		/* and track the device position for the two points involved. */
		/* Also locate the smallest device distance. */
		cdist = 1e200;
		fdist = -1.0;
		for (i = 0; i < (ii-1); i++) {
			for (j = i+1; j < ii; j++) {
				double dist;
				dist = ofps_comp_eperr3(s, NULL, cx.nds[i]->v, cx.nds[i]->p, cx.nds[j]->v, cx.nds[j]->p, cx.nds[i]->nsp) ;
				if (dist < ceperr) {
					ceperr = dist;
					bi = i;
					bj = j;
				}
				for (dist = 0.0, e = 0; e < di; e++) {
					double tt = cx.nds[i]->p[e] - cx.nds[j]->p[e];
					dist += tt * tt;
				}
				if (dist < cdist)
					cdist = dist;
				if (dist > fdist)
					fdist = dist;
			}
		}
		
		fdist = sqrt(fdist);
		cdist = sqrt(cdist);

		/* Compute the middle of the two closest eperr nodes */
		for (e = 0; e < di; e++)
			mct[e] = 0.5 * (cx.nds[bi]->p[e] + cx.nds[bj]->p[e]);

		/* Set a step/search radius based on the distance */
		/* between the two closest device distance nodes. */
		if (cdist < COINTOL) {
#ifdef DEBUG
			printf("Two nodes are cooincident! - dnsq will fail!\n");
#endif
			if (s->verb > 1)
				warning("Two nodes are cooincident! ix %d, pos %s and ix %d pos %s",cx.nds[bi]->ix,ppos(di,cx.nds[bi]->p),cx.nds[bj]->ix,ppos(di,cx.nds[bj]->p));
		}
		bsrad = 0.2 * cdist;
		if (bsrad < 1e-5)
			bsrad = 1e-5;
	}

	/* Set initial starting position */
//	if (startex && ! ofps_would_clip_point(s, vv->p)) { }
	if (startex) {

		for (e = 0; e < di; e++)
			osp[e] = vv->p[e];

//		ofps_clip_point(s, vv->p, vv->p);
#ifdef DEBUG
		printf("Startex startposition = %s\n",ppos(di,atp));
#endif
	} else {
		double mwt = 0.3;
		/* Start at equalateral point between two closest */
		/* nodes towards average. */
		for (e = 0; e < di; e++) {
//			osp[e] = mct[e];		/* Best for 2D ? */
//			osp[e] = atp[e];		/* best for 3D/4D ? */
			osp[e] = mwt * mct[e] + (1.0 - mwt) * atp[e];	/* Good compromize */
		}
	}

	/* Try our computed starting position first, and if that fails, */
	/* retry with a random offset starting location. */
	cx.srad = bsrad;
	for (e = 0; e < di; e++)
		cx.stp[e] = osp[e];
	cx.on = 0;

	for (tries = 0; tries < notries; tries++) {
		int rv;
		double fvec[MXPD];			/* Return function value at solution */
		int cfunccount;

		if (tries > 0) {	/* Determine a starting point */

			/* Try all possible odd one outs */
			cx.on++;

			/* On carry, use a random start offset */
			/* (Tried culling random starts and odd one outs */
			/* by picking one with a low norm, but */
			/* while this reduced the number of small */
			/* retries, it worsened the number of failures */
			/* and sucesses after a large number of retries.) */
			if (cx.on >= cx.nn) {
				double rscale = 1.0;	/* Random scale */
				double fval[MXPD];
				int nc;

				s->sob->next(s->sob, cx.stp);

				/* Scale random value around original starting point */
				for (e = 0; e < di; e++) {
					cx.stp[e] = cx.stp[e] * 2.0 - 1.0;		/* Make -1.0 to 1.0 range */
					if (cx.stp[e] < 0.0) {
						cx.stp[e] *= rscale * (osp[e] - s->imin[e]);
					} else {
						cx.stp[e] *= rscale * (s->imax[e] - osp[e]);
					}
					cx.stp[e] += atp[e];
				}
				ofps_clip_point4(s, cx.stp, cx.stp);
				cx.on = 0;
			}
		}

		/* Set start position */
		for (e = 0; e < di; e++)
			vv->p[e] = cx.stp[e];

#ifdef DEBUG
		printf("Starting location = %s, srad = %f, on = %d\n",ppos(di,cx.stp),cx.srad,cx.on);
#endif

//printf("\nStarting location = %s, srad = %f\n",ppos(di,cx.stp),cx.srad);
		/* Locate vertex */
		cfunccount = s->funccount;
		s->dnsqs++;
		if (tcalls == 0)
			maxfev = 500;
		else
			maxfev = 2 * tfev/tcalls; 
		rv = dnsqe((void *)&cx, dnsq_solver, NULL, di, vv->p, cx.srad, fvec, 0.0, ftol, maxfev, 0);
		if ((s->funccount - cfunccount) > 20) {
//printf("More than 20: %d\n",s->funccount - cfunccount);
		}
		if ((s->funccount - cfunccount) > s->maxfunc) {
			s->maxfunc = (s->funccount - cfunccount);
//printf("New maximum %d\n",s->maxfunc);
		}

		if (rv != 1 && rv != 3) {
			/* Fail to converge */
#ifdef DEBUG
			printf("dnsqe fail to converge, retuned %d\n",rv);
#endif
		} else {
			double ee[MXPD+1];
			double max, min;
			double ple = 0.0;		/* Gamut plane error */ 

			/* Evaluate the result */

			/* Update average function evaluations */
			tcalls++;
			tfev += s->funccount - cfunccount;
			
#ifdef DEBUG
			printf("dnsq pos %s\n",ppos(di,vv->p));
			if (rv == 3)
				printf("dnsq returned 3 - dtol too small\n");
#endif

			/* Compute perceptual */
			ofps_cc_percept(s, vv->v, vv->p);

			/* Compute the eperr's for each node. */
			ofps_pn_eperr(s, vv->ce, ee, vv->v, vv->p, cx.nds, cx.nn);

			min = 1e80, max = -1e80;
			for (e = 0; e < cx.nn; e++) {
				if (min > ee[e])
					min = ee[e];
				if (max < ee[e])
					max = ee[e];
			}

			/* Compute the worst plane equation error */
			for (k = 0; k < cx.nsp; k++) {
				double v;
				for (v = cx.sp[k]->pe[di], e = 0; e < di; e++) 	/* Compute relation to plane equation */
					v += cx.sp[k]->pe[e] * vv->p[e];
				v = fabs(v * FGPMUL);
//printf("~1 gamut plane %d err = %f\n",k,v);
				if (v > ple)
					ple = v;
			}

#ifdef DEBUG
			printf("new vertex pos %s has eperrs match by %f & gamut plane %f\n",ppos(di,vv->p),max-min,ple);
#endif
			/* If not dtol too large, Check that the balance is acceptable */
			if (/* rv != 3 && */
			     (((cx.nn > 1) && (max - min) > (ftol * 2.0))
			   || ((cx.nsp > 0) && ple > (ftol * 2.0))))	{
				/* Don't use this */
#ifdef DEBUG
				printf("new vertex pos %s doesn't have sufficient matching eperrs and on gamut plane\n",ppos(di,vv->p),max-min,ple);
#endif
			} else {
				/* eperr balance is acceptable, so further */
				/* evaluate the location found. */
				double ss;

				s->sucfunc += (s->funccount - cfunccount);
				s->sucdnsq++;

				/* Compute how much the result is out of gamut */
				vv->oog = ofps_oog(s, vv->p);
#ifdef DEBUG
				if (vv->oog > 0.01)
					printf("dnsq returned out of gamut result by %e\n", vv->oog);
#endif

				/* Compute errors at returned location */
				vv->eperr = ofps_eperr2(ee, cx.nn);
				vv->eserr = ofps_eserr2(vv->ce, ee, cx.nn);

				/* Decide whether a vertex location is acceptable */
				/* We accept a point that has an acceptable error balance */
				/* and improves the eperr, and is in gamut if this is not a repos. */
				/* (There's some mystery stuff in here for fixups) */
				if ((cx.nn <= 1) || ((max - min) <= (ftol * 2.0)
				 &&  ((!repos && vv->oog <= 0.01 && vv->eperr < (vv->ceperr + 0.1))
				   || ( repos && vv->oog <= 0.01 && vv->eperr < (5.0 * vv->ceperr + 20.0))
				   || ( repos && vv->oog > 0.0 && vv->eperr < 1000.0) 
				   || ( fixup && vv->oog < 20.0 && vv->eperr < (vv->ceperr + 0.01))
				))) {

					if (tries > s->maxretries)
						s->maxretries = tries;
#ifdef DEBUG
					printf(" - comb %s suceeded on retry %d (max %d)\n",pcomb(di,vv->nix),tries,s->maxretries);
					printf("       oog = %f, eperr = %f, ceperr = %f\n",vv->oog,vv->eperr,vv->ceperr);
#endif
//if (tries > 10)
//	printf(" - comb %s suceeded on retry %d (max %d)\n",pcomb(di,vv->nix),tries,s->maxretries);
// 
//printf("Solution for comb %s has eperr %f < ceperr %f and not out of gamut by %f, retry %d\n",pcomb(di,vv->nix),vv->eperr,vv->ceperr,vv->oog,tries+1);
//printf("Solution is at %s (%s)\n",ppos(di,vv->p),ppos(di,vv->v));

// 
//if (repos) printf("~1 vtx no %d dtav = %f, fdist = %f\n",vv->dtav,fdist);
//if (repos && vv->vv->no == 889) printf("~1 vtx no %d dtav = %f, fdist = %f\n",vv->vv->no,vv->dtav,fdist);

#ifdef DUMP_FERR 		/* Create .tiff of dnsq function error */
					if (tries >= DUMP_FERR) {
						printf("Suceeded on retry %d, dumping debug rasters\n",tries);

						/* Re-run the last unsucessful dnsq, to trace the path */
						pcx.debug = 1;
						pcx.clist = NULL;
						pcx._nl = 0;
						pcx.nl = 0;
						dnsqe((void *)&pcx, dnsq_solver, NULL, di, pcx.stp, pcx.srad, fvec, 0.0, ftol, maxfev, 0);
						pcx.debug = 0;
						dump_dnsqe(s, "dnsq_fail2.tif", vv->nix, &pcx);
						free(pcx.clist);

						/* Re-run the first unsucessful dnsq, to trace the path */
						pcx.debug = 1;
						pcx.clist = NULL;
						pcx._nl = 0;
						pcx.nl = 0;
						pcx.on = 0;					/* First odd one out */
						for (e = 0; e < di; e++)
							pcx.stp[e] = atp[e];		/* best start ? */
						dnsqe((void *)&pcx, dnsq_solver, NULL, di, pcx.stp, pcx.srad, fvec, 0.0, ftol, maxfev, 0);
						pcx.debug = 0;
						dump_dnsqe(s, "dnsq_fail1.tif", vv->nix, &pcx);
						free(pcx.clist);

						/* Re-run the sucessful dnsq, to trace the path */
						cx.debug = 1;
						cx.clist = NULL;
						cx._nl = 0;
						cx.nl = 0;
						dnsqe((void *)&cx, dnsq_solver, NULL, di, cx.stp, cx.srad, fvec, 0.0, ftol, maxfev, 0);
						cx.debug = 0;
						dump_dnsqe(s, "dnsq_suc.tif", vv->nix, &cx);
						free(cx.clist);
						exit(0);
					}
#endif

					break;				/* Use the result now */
				}
#ifdef DEBUG
				printf("Solution for comb %s has eperr %f > ceperr %f or out of gamut by %f, retry %d\n",pcomb(di,vv->nix),vv->eperr,vv->ceperr,vv->oog,tries+1);
#endif
			}
		}
		pcx = cx;		/* Save unsucessful context for debug */
	}	/* Retry */

	/* If we've run out of tries, return the best solution we found */
	if (tries >= notries) {

		/* Show up if this ever gets used */
		for (e = 0; e < di; e++)
			vv->p[e] = -0.1;
		ofps_cc_percept(s, vv->v, vv->p);
#ifdef DEBUG
		printf("vertex location solving failed after %d tries\n",tries);
#endif
		if (s->verb > 1)
			warning("vertex location solving failed after %d tries",tries);
		return 2;		/* Don't use this vertex */
	}

#ifdef DEBUG
	printf("Returning comb %s opt pos = %s val = %s, eperr = %f, eserr = %f\n",pcomb(di,vv->nix),ppos(di,vv->p),ppos(di,vv->v),vv->eperr,vv->eserr);
#endif

	return 0;
}

/* --------------------------------------------------------- */
/* Deal with creating a dummy vertex to represent one that */
/* can't be positioned. We simply locate the best point we can. */

/* calculate the functions at x[] */
double powell_solver(	/* Return < 0 on abort */
	void *fdata,	/* Opaque data pointer */
	double *x		/* Multivariate input values */
) {
	vopt_cx *cx = (vopt_cx *)fdata;
	ofps *s = cx->s;
	int k, e, di = s->di;
	int nn_1 = cx->nn-1;
	double sv[MXPD];
	double cee[MXPD+1], teperr;
	double ss, oog;
	double rv = 0.0;

//printf("~1 powell_solver got %d nodes and %d planes\n",cx->nn,cx->nsp);

	/* Get eperr at each real node */
	ofps_cc_percept(s, sv, x);	/* We have to compute it */
	for (k = 0; k < cx->nn; k++)
		cee[k] = ofps_comp_eperr2(s, NULL, sv, x, cx->nds[k]->v, cx->nds[k]->p, cx->nds[k]->nsp);

//fprintf(stderr,"~1 maxeperr = %f\n",cmax);

//printf("~1 eprror =");
//for (k = 0; k < cx->nn; k++)
//	printf(" %f",cee[k]);
//printf("\n");

	/* The error is zero if the input value */
	/* is within gamut, all the real node eperr's are */
	/* the same, and the ditance to gamut planes is zero. */

	/* Compute average eperr */
	for (teperr = 0.0, k = 0; k < cx->nn; k++)
		teperr += cee[k];	
	teperr /= (double)cx->nn;

//printf("~1 average = %f\n", teperr);

	/* Add diference to average */
	for (k = e = 0; k < cx->nn; k++) {
		if (k != cx->on) {
			double tt; 
			tt = teperr - cee[k];
			rv += tt * tt;
		}
	}
//printf("~1 after diff to avg rv = %f\n", rv);

	/* Compute distances to planes */
	for (k = 0; k < cx->nsp; k++) {
		pleq *pl = cx->sp[k];
		double v;

		for (v = pl->pe[di], e = 0; e < di; e++) 	/* Compute relation to plane equation */
			v += pl->pe[e] * x[e];
		v *= FGPMUL;
		rv += v * v;
	}
//printf("~1 after diff to planes rv = %f\n", rv);

	/* Compute distance out of gamut */

	for (ss = oog = 0.0, e = 0; e < di; e++) {
		if (x[e] < (s->imin[e])) {
			double tt = s->imin[e] - x[e];
			if (tt > oog) oog = tt; 
		} else if (x[e] > (s->imax[e])) {
			double tt = x[e] - s->imax[e];
			if (tt > oog) oog = tt; 
		}
		ss += x[e];
	}
	if (ss > s->ilimit) {
		double tt;
		ss = (ss - s->ilimit)/di;	/* Axis aligned distance to ink limit */
		tt = sqrt((double)di) * ss;	/* Diagonal distance to ink limit */
		if (tt > oog)
			oog = tt; 
	}

	rv += 1000.0 * oog * oog;
//printf("~1 after oog rv = %f\n", rv);

//printf("powell_solver returning %f from %s\n",rv, ppos(di,x));
//fprintf(stderr,"powell_solver returning %f from %s\n",rv, ppos(di,x));

	return 0;
}

/* Fake up a vertex position when position_vtx() has failed. */
static void dummy_vtx_position(
	ofps *s,
	vtx *ev1, vtx *ev2,	/* Deleted and non-deleted vertexes */
	nodecomb *vv		/* Return the location and its error */
) {
	int e, di = s->di;
	int k, ii;
	vopt_cx cx;						/* dnsq context */
	double ss[MXPD];
	double bl;						/* Location on path between del and !del vertexes */
	double ee[MXPD+1];

#ifdef DEBUG
	printf("dummy_vtx_position called for comb %s\n",pcomb(di,vv->nix));
#endif

#ifdef DUMP_FERR
	cx.debug = 0;
#endif

	/* Setup for dnsq to optimize for equal eperr */
	cx.s = s;

	/* Pointers to real nodes. Although we allow for the */
	/* fake inner/outer nodes, eperr() will fail them later. */
	for (ii = e = 0; e <= di; e++) {
		if (vv->nix[e] >= 0 || vv->nix[e] < -s->nbp)
			cx.nds[ii++] = s->n[vv->nix[e]];
	}
	cx.nn = ii;

	if (ii == 0) {
		fflush(stdout);
		error("ofps: unexpectedely got no real nodes in vertex position %s",pcomb(di,vv->nix));
	}

#ifdef DEBUG
		printf("%d real nodes\n",ii);
		for (k = 0; k < ii; k++)
			printf("Node ix %d at %s (%s)\n",cx.nds[k]->ix,ppos(di,cx.nds[k]->p),ppos(di,cx.nds[k]->v));
#endif

	/* Setup gamut suface planes */
	cx.nsp = 0;
	if (ii < (di+1)) {
		/* Create list of pointers to the gamut surface planes involved */
		for (e = 0, k = ii; k <= di; e++, k++) {
#ifdef DEBUG
			printf("Adding plane for node ix %d\n",vv->nix[k]);
#endif
			cx.sp[e] = &s->gpeqs[-1 - vv->nix[k]];
		}
		cx.nsp = e;
	}

	/* Set search area */
	for (e = 0; e < di; e++)
		ss[e] = 0.001;

	/* Compute a position on the locus between the del and !del vertexes */
	bl = (ev1->nba_eperr - ev2->eperr)/(ev1->eperr - ev2->eperr);
	if (bl < 0.0)
		bl = 0.0;
	else if (bl > 1.0)
		bl = 1.0;
	for (e = 0; e < di; e++) { 
		vv->p[e] = bl * vv->v1[0]->p[e] + (1.0 - bl) * vv->v2[0]->p[e];
	}

	ofps_clip_point5(s, vv->p, vv->p);

#ifndef NEVER
	/* Seem to often fail due to pathalogical condition for max type eperr() */
	if (powell(NULL, di, vv->p, ss, 1e-5, 1000, powell_solver, &cx, NULL, NULL)) {
		warning("dummy_vtx_position powell failed");
	}
#endif

	ofps_clip_point5(s, vv->p, vv->p);

	/* Compute perceptual (was clipped above) */
	s->percept(s->od, vv->v, vv->p);

	/* Compute the eperr's for each node. */
	ofps_pn_eperr(s, vv->ce, ee, vv->v, vv->p, cx.nds, cx.nn);

	/* Compute errors at returned location */
	vv->eperr = ofps_eperr2(ee, cx.nn);
	vv->eserr = ofps_eserr2(vv->ce, ee, cx.nn);

#ifdef DEBUG
	printf("Returning comb %s opt pos = %s val = %s, eperr = %f, eserr = %f\n",pcomb(di,vv->nix),ppos(di,vv->p),ppos(di,vv->v),vv->eperr,vv->eserr);
#endif
}

/* --------------------------------------------------- */
/* Vertex add routines */

/* Comlete adding a node to a Voronoi surface. */
/* It's assumed that the hit nodes have been added to the s->nxh list */ 
/* and the s->nvcheckhits set to the number of hit vertexes. */
/* The nodes may be fake gamut boundary nodes, but must have */
/* real vertexes. Vertexes that are marked for deletion or new ones */
/* will be added to the batch update list for later execution. */
/* Nodes that have had vertexes added to them will added to the node 'to be updated' list */
/* and then a batch update will be executed. */
/* Return 0 if it wasn't added, */
/* Return 1 if it was added */
static int add_to_vsurf(
ofps *s,
node *nn,		/* Node to add */
int fixup,		/* 0 = seed, 1 = fixup ?? */
int abortonfail	/* 0 = ignore position failures, 1 = abort add if there are any failures */
) {
	int e, ff, f, di = s->di;
	int i, j, k, ndi;
	vtx *tev;
	vtx *ev1, *ev2;	/* Deleted and non-deleted vertexes */
	int ndelvtx;	/* Number of vertexes to delete */ 
	int nncombs;	/* Number of node combinations generated, allocated. */

#ifdef DEBUG
	printf("\nAdd_to_vsurf node ix %d (p %s), i_sm %s, a_sm %s\n",nn->ix, ppos(di,nn->p),psm(s,&s->sc[nn->pmask].i_sm),psm(s,&s->sc[nn->pmask].a_sm));
#endif

#ifdef DEBUG
	if (nn->ix < -s->nbp) {
		printf("Fake node involved\n");
	}
#endif

	/* Update stats for the hit vertexes */
	s->nsurfadds++;
	s->nhitv += s->nvcheckhits;
	if (s->nvcheckhits > s->maxhitv)
		s->maxhitv = s->nvcheckhits;

	if (s->nvcheckhits == 0) {	/* Node doesn't improve any vertex eperrs */
#ifdef DEBUG
		printf("Add_to_vsurf done - not better, not added\n");
#endif
		return 0;
	}

#ifdef DEBUG
	printf("There are %d vertexes to replace, and %d potential replacement vertexes\n",s->nvcheckhits, di * s->nvcheckhits);
	printf("Vertexes marked for deletion are:\n");
	for (ev1 = s->nxh; ev1 != NULL; ev1 = ev1->nxh)
		printf("  vtx no %d nix %s\n",ev1->no,pcomb(di,ev1->nix));
#endif

	/* Generate all the potential new node combinations/replacement verticies. */
	/* We check each deleted vertex against its non-deleted neighbours. */
	/* The same replacement combination may be generated more than once. */
	for (nncombs = 0, ev1 = s->nxh; ev1 != NULL; ev1 = ev1->nxh) {

		for (ndi = 0; ndi < ev1->nnv; ndi++) {
			int nix[MXNIX];
#ifdef INDEP_SURFACE
			setmask cvm;		/* visibility setmask for each vertex combination */
#endif
			ev2 = ev1->nv[ndi];

			if (ev2->del != 0)
				continue;		/* Can't pair with another deleted vertex */

#ifdef DEBUG
			printf("\nDealing with vertex pair del no %d nix %s, and !del no %d nix %s\n",ev1->no,pcomb(di,ev1->nix),ev2->no,pcomb(di,ev2->nix));
#endif

#ifdef DEBUG
		{
			int aa, bb, cc;		/* Probable hit check */
			int nnm, nmix;

			/* Use the nixm to quickly check if all but one parent node matches */
			aa = ev1->nix[MXPD+2];	/* nixm */
			bb = ev2->nix[MXPD+2];	/* nixm */
			if ((aa & bb) == 0 || (cc = aa & ~bb, (cc & (cc-1)) != 0)) {
				error("Vertexes %d comb %s and %d comb %s are vn neighbours that shouldn't be!", ev1->no,pcomb(di,ev1->nix),ev2->no,pcomb(di,ev2->nix));
			}

			/* Do an exact check of all except one node match */
			for (nnm = ff = e = 0; e <= di; e++) {
				for (f = ff; f <= di; f++) {
					if (ev1->nix[e] == ev2->nix[f]) {
						ff = f;			/* Start from here next time */
						break;
					}
					if (ev1->nix[e] > ev2->nix[f])	/* No point in looking further */
						f = di;
				}
				if (f > di) {	/* Didn't match */
					if (++nnm > 1)
						break;
					nmix = e;
				}
			}
			if (e <= di) {
				error("Vertexes %d comb %s and %d comb %s are vn neighbours that shouldn't be!", ev1->no,pcomb(di,ev1->nix),ev2->no,pcomb(di,ev2->nix));
			}
		}
#endif	/* DEBUG */
			/* Create the node combination */
			for (e = 0; e <= di; e++) {
				nix[e] = ev1->nix[e];
				for (f = 0; f <= di; f++) {
					if (nix[e] == ev2->nix[f])
						break;
				}
				if (f > di) 	/* Found one different */
					nix[e] = nn->ix;
			}
			sort_nix(s, nix);

			/* Check that the same node doesn't appear twice */
			/* (~~99 Why do we need this - does it ever happen ??) */
			for (e = 0; e < di; e++) {
				for (k = e+1; k <= di; k++) {
					if (nix[e] == nix[k]) {
#ifdef DEBUG
						printf("New vertex with duplicate nodes %s from vertexes %d comb %s and %d comb %s ignored\n",pcomb(di,nix),ev1->no,pcomb(di,ev1->nix),ev2->no,pcomb(di,ev2->nix));
						if (s->verb > 1)
							warning("New vertex with duplicate nodes %s from vertexes %d comb %s and %d comb %s ignored",pcomb(di,nix),ev1->no,pcomb(di,ev1->nix),ev2->no,pcomb(di,ev2->nix));
#endif
						break;
					}
				}
				if (k <= di)
					break;
			}
			if (e < di)	
				continue;

			/* See if the combination has at least one real node, */
			/* and none of the inner or outer fake nodes. */
			/* (~~99 Do we need this - does it ever happen ??) */
			k = 0;
			for (e = 0; e <= di; e++) {
				if (nix[e] >= 0)
					k |= 1;			/* Found one real node */
				else if (nix[e] < -s->nbp)
					break;				/* There's a fake inner or outer node though */
			}
			if (e <= di || k == 0) {
#ifdef DEBUG
				printf("Combination ix: %s, skipped because it has no real nodes\n",pcomb(di,nix));
				if (s->verb > 1)
					warning("Combination ix: %s, skipped because it has no real nodes",pcomb(di,nix));
#endif
				continue;
			}

#ifdef INDEP_SURFACE
			/* Compute the pertinent visibility mask for this vertex creation. */
			/* Note that we keep pairs of vertexes that aren't visible to each other */
			/* so that we can add them to the vertex net. */
			sm_andand(s, &cvm, &ev1->vm, &ev2->vm, &s->sc[comp_cmask(s, nix)].a_sm);
#ifdef DEBUG
			printf("Combination ix: %s, vm %s, eperrs %f to %f\n",pcomb(di,nix),psm(s,&cvm),ev1->eperr,ev2->eperr);
#endif
#else 	/* !INDEP_SURFACE */
#ifdef DEBUG
			printf("Combination ix: %s, eperrs %f to %f\n",pcomb(di,nix),ev1->eperr,ev2->eperr);
#endif
#endif	/* !INDEP_SURFACE */

			/* See if this combination is already in the list */
			/* due to a pair having the same common nodes. */
			for (k = 0; k < nncombs; k++) {

				if (s->combs[k].nix[MXPD+1] != nix[MXPD+1]) /* Hashes don't match */
					continue;

				for (e = 0; e <= di; e++) {					/* Do full check */
					if (s->combs[k].nix[e] != nix[e]) {
						break;	/* No match */
					}
				}
				if (e > di) {	/* Match */
					double ceperr;
					if (s->combs[k].count >= s->combs[k]._count) {
						s->combs[k]._count = 2 * s->combs[k]._count + 5;
						if ((s->combs[k].v1 = (vtx **)realloc(s->combs[k].v1,
						                              sizeof(vtx *) * s->combs[k]._count)) == NULL)
							error ("ofps: malloc failed on node combination vertex list %d",
							                                                    s->combs[k]._count);
						if ((s->combs[k].v2 = (vtx **)realloc(s->combs[k].v2,
						                              sizeof(vtx *) * s->combs[k]._count)) == NULL)
							error ("ofps: malloc failed on node combination vertex list %d",
							                                                    s->combs[k]._count);
					}
					s->combs[k].v1[s->combs[k].count] = ev1;
					s->combs[k].v2[s->combs[k].count] = ev2;
					s->combs[k].count++;

					/* Update ceperr if this is higher */
					ceperr = ev1->eperr > ev2->eperr ? ev1->eperr : ev2->eperr;
					if (ceperr > s->combs[k].ceperr)
						s->combs[k].ceperr = ceperr;
#ifdef DEBUG
					printf("Vertex generation count now %d with ceperr %f\n",s->combs[k].count,s->combs[k].ceperr);
#endif
#ifdef INDEP_SURFACE
					sm_or(s, &s->combs[k].vm, &s->combs[k].vm, &cvm);
#ifdef DEBUG
					printf("Vertex combination vm now %s\n",psm(s,&s->combs[k].vm));
#endif
#endif
					break;
				}
			}
			if (k <  nncombs) 	
				continue;		/* Already on list */

			/* Add this combination to the list as a new entry */
			if (nncombs >= s->_ncombs) {
				int o_ncombs = s->_ncombs;
				s->_ncombs = 2 * s->_ncombs + 5;
				if ((s->combs = (nodecomb *)realloc(s->combs, sizeof(nodecomb) * s->_ncombs)) == NULL)
					error ("ofps: malloc failed on node combination array length %d", s->_ncombs);
				memset((void *)(s->combs + o_ncombs), 0,
				                                      (s->_ncombs - o_ncombs) * sizeof(nodecomb));
			}

			if (1 >= s->combs[nncombs]._count) {
				s->combs[nncombs]._count = 2 * s->combs[nncombs]._count + 5;
				if ((s->combs[nncombs].v1 = (vtx **)realloc(s->combs[nncombs].v1,
				                              sizeof(vtx *) * s->combs[nncombs]._count)) == NULL)
					error ("ofps: malloc failed on node combination vertex list %d",
					                                                    s->combs[nncombs]._count);
				if ((s->combs[nncombs].v2 = (vtx **)realloc(s->combs[nncombs].v2,
				                              sizeof(vtx *) * s->combs[nncombs]._count)) == NULL)
					error ("ofps: malloc failed on node combination vertex list %d",
					                                                    s->combs[nncombs]._count);
			}
			s->combs[nncombs].v1[0] = ev1;
			s->combs[nncombs].v2[0] = ev2;
			s->combs[nncombs].count = 1;
			for (e = 0; e <= di; e++)
				s->combs[nncombs].nix[e] = nix[e];
			s->combs[nncombs].nix[MXPD+1] = nix[MXPD+1];	/* Copy Hash */
			s->combs[nncombs].nix[MXPD+2] = nix[MXPD+2];	/* Copy nixm */
			s->combs[nncombs].ceperr = ev1->eperr > ev2->eperr ? ev1->eperr : ev2->eperr;
			s->combs[nncombs].startex = 0;
			s->combs[nncombs].pvalid = 0;
			s->combs[nncombs].vv = NULL;
#ifdef INDEP_SURFACE
			sm_cp(s, &s->combs[nncombs].vm, &cvm);
#else	/* !INDEP_SURFACE */
			sm_set(s, &s->combs[nncombs].vm, 0);		/* Not used */
#endif	/* !INDEP_SURFACE */

			nncombs++;
			
#ifdef DEBUG
			printf("Adding combination to list with ceperr %f, vm %s, list size %d\n",s->combs[nncombs-1].ceperr,psm(s,&s->combs[nncombs-1].vm), nncombs);
#endif
		}
	}

#ifdef DEBUG
	printf("\nThere are %d unique node combinations in list, locating combs. in list:\n",nncombs);
#endif

	/* Locate the replacement vertex positions */
	for (i = 0; i < nncombs; i++) { 

		ev1 = s->combs[i].v1[0];
		ev2 = s->combs[i].v2[0];

#ifdef INDEP_SURFACE
		/* Ignore pairs of vertexes that don't form a visible new combination */
	    if (sm_test(s, &s->combs[i].vm) == 0) {
#ifdef DEBUG
			printf("Combination ix: %s, skipped because vm %s == 0x0\n",pcomb(di,s->combs[i].nix),psm(s,&s->combs[i].vm));
#endif
			continue;
		}
#endif	/* INDEP_SURFACE */

#ifdef DEBUG
		printf("\nNode combination ix: %s\n",pcomb(di,s->combs[i].nix));
#endif
		/* Try and locate existing vertex that is due to the same nodes */
		if ((s->combs[i].vv = vtx_cache_get(s, s->combs[i].nix)) != NULL) {
#ifdef DEBUG
			printf("Vertex is same as existing no %d\n",s->combs[i].vv->no);
#endif

			s->combs[i].vv->add = 2;		/* Updated vertex */

			/* If a new vertex is not the same as the two nodes it's being created */
			/* from yet has been marked for deletion, reprieve it. */
			if (s->combs[i].vv->del) {

				s->combs[i].vv->del = 0;
#ifdef DEBUG
				printf("New existing vertex no %d is deleted vertex - reprieve it\n",s->combs[i].vv->no);
#endif
			}
		}

		/* We need to create a replacement vertex, locate position for it */
		if (s->combs[i].vv == NULL) {
#ifdef DEBUG
			printf("About to locate comb ix: %s, ceperr %f\n",pcomb(di,s->combs[i].nix),s->combs[i].ceperr);
#endif
//printf("~1 About to locate comb ix: %s, ceperr %f\n",pcomb(di,s->combs[i].nix),s->combs[i].ceperr);
			/* Compute a starting position between the deleted/not deleted pair */
			/* This seems very slightly better than the default mct[] + atp[]  scheme. */
			if (nn->ix >= 0) {	/* If not boundary */
				double bl;
				bl = (ev1->nba_eperr - ev2->eperr)/(ev1->eperr - ev2->eperr);
				if (bl < 0.0)
					bl = 0.0;
				else if (bl > 1.0)
					bl = 1.0;
				for (e = 0; e < di; e++) { 
					s->combs[i].p[e] = bl * s->combs[i].v1[0]->p[e] + (1.0 - bl) * s->combs[i].v2[0]->p[e];
				}
				ofps_clip_point5(s, s->combs[i].p, s->combs[i].p);
//printf("Startex is %s\n",ppos(di,s->combs[i].p));
				s->combs[i].startex = 1;
			}
			/* find vertex position of max eperr */
			if (position_vtx(s, &s->combs[i], s->combs[i].startex, 0, fixup) != 0) {
				if (s->verb > 1)
					warning("Unable to locate vertex at node comb %s\n",pcomb(di,s->combs[i].nix));
				s->posfails++;
				s->posfailstp++;
				if (abortonfail)
					break;

			} else {
				s->combs[i].pvalid = 1;
			}
		}
	}	/* Next replacement vertex */

	/* If we aborted because abortonfail is set and we failed to place a new node, */
	/* erase our tracks and return failure. */
	if (i < nncombs) {
		for (i = 0; i < nncombs; i++) { 
			if (s->combs[i].vv != NULL) {
				s->combs[i].vv->add = 0;
				s->combs[i].vv->del = 0;
			}
		}
		return 0;
	}

#ifdef DEBUG
	printf("\nNow converting positioned combinations to vertexes\n");
#endif
	/* Convert from computed position to vertexes */
	for (i = 0; i < nncombs; i++) { 

#ifdef INDEP_SURFACE
		/* Ignore combo that doesn't form a visible new vertex */
		if (sm_test(s, &s->combs[i].vm) == 0) 
			continue;
#endif	/* INDEP_SURFACE */

		if (s->combs[i].vv == NULL) {		/* Not an existing vertex */

			if (s->combs[i].pvalid == 0) {	/* No valid vertex found */

				/* [ If a valid vertex location was not found we tried */
				/*   using the deleted vertex's location instead, and */
				/*   relying on it getting deleted at some later stage. */
				/*   This seems to stuff the incremental fixup code up */
				/*   completely, so we create a summy vertex position instead. ] */

				/* Fake up a vertex position when position_vtx() has failed. */
				dummy_vtx_position(s, ev1, ev2, &s->combs[i]);
				goto lnew_vtx;

			} else {

 lnew_vtx:;
				/* Allocate space for new vertex, and create it from location */
				s->combs[i].vv = new_vtx(s);

#ifdef DEBUG
				printf("Converting comb %s to from del %d !del %d to vtx no %d vm %s\n",pcomb(di,s->combs[i].nix),s->combs[i].v1[0]->no, s->combs[i].v2[0]->no, s->combs[i].vv->no,psm(s,&s->combs[i].vm));
#endif

				for (e = 0; e < di; e++) {
					s->combs[i].vv->nix[e] = s->combs[i].nix[e];
					s->combs[i].vv->ce[e] = s->combs[i].ce[e];
					s->combs[i].vv->p[e] = s->combs[i].p[e];
					s->combs[i].vv->v[e] = s->combs[i].v[e];
				}
				s->combs[i].vv->nix[e] = s->combs[i].nix[e];
				s->combs[i].vv->nix[MXPD+1] = s->combs[i].nix[MXPD+1];	/* Copy Hash */
				s->combs[i].vv->nix[MXPD+2] = s->combs[i].nix[MXPD+2];	/* Copy nixm */

				s->combs[i].vv->eperr = s->combs[i].eperr;
				s->combs[i].vv->eserr = s->combs[i].eserr;

				/* Count the number of gamut surfaces the vertex falls on */
				det_vtx_gsurf(s, s->combs[i].vv);

				/* Check if the node and vertex cooincide, and aren't going to move */
				if (nn->nsp == di && s->combs[i].vv->nsp == di
				 && nn->pmask == s->combs[i].vv->pmask) {
//printf("~1 Trapped node and vertex coincide - mark vertex as ghost\n");
						s->combs[i].vv->ghost = 1;
				}
				s->combs[i].vv->del = 0;
				s->combs[i].vv->add = 1;	/* New vertex */
			}
		}

		if (s->combs[i].vv != NULL) {		/* There is a new or updated vertex */
#ifdef DEBUG
			printf("Vertex no %d pmask 0x%x cmask 0x%x vm %s at %s being added to batch list\n",s->combs[i].vv->no, s->combs[i].vv->pmask, s->combs[i].vv->cmask, psm(s,&s->combs[i].vm),ppos(di,s->combs[i].vv->p));
#endif
			/* Add to batch update list if it is not already there */
			if (s->combs[i].vv->bch == 0) {
//printf("~1 adding vtx 0x%x no %d to batch list\n",s->combs[i].vv,s->combs[i].vv->no);
				s->combs[i].vv->batch = s->batch;
				s->batch = s->combs[i].vv;
				s->combs[i].vv->bch = 1;
			}

#ifdef INDEP_SURFACE
			/* Set or update the vm */ 
			sm_or(s, &s->combs[i].vv->buvm, &s->combs[i].vv->buvm, &s->combs[i].vm);
//printf("~1 Vertex no %d buvm set/update to %s\n",s->combs[i].vv->no,psm(s,&s->combs[i].vv->buvm));
#endif	/* INDEP_SURFACE */
		}
	}

	/* Add all vtx marked for deletion to the batch update list. */ 
	for (ev1 = s->nxh; ev1 != NULL; ev1 = ev1->nxh) {
			
#ifdef DEBUG
		printf("Vertex no %d being added to pending delete batch list, bdvm %s\n",ev1->no,psm(s,&ev1->bdvm));
#endif

		if (ev1->bch == 0) {
//printf("~1 adding vtx 0x%x no %d to batch list\n",ev1,ev1->no);
			ev1->batch = s->batch;
			s->batch = ev1;
			ev1->bch = 1;
		}

#ifdef INDEP_SURFACE
		/* Add node setmask to those that will be removed delete vertex visibility */ 
# ifdef USE_DISJOINT_SETMASKS
		sm_orand(s, &ev1->bdvm, &ev1->bdvm, &s->sc[nn->pmask].a_sm, &s->sc[ev1->cmask & nn->pmask].a_sm);
# else
		sm_or(s, &ev1->bdvm, &ev1->bdvm, &s->sc[nn->pmask].a_sm);
# endif
#endif	/* INDEP_SURFACE */
	}

	/* Do first part of batch update. */
	/* This will reset ->del on nodes that will be retained */
	do_batch_update1(s, fixup);

	/* Remove deleted vertex's from the vertex net, and their */
	/* parent nodes. */
	{
		vtx *vx1, *vx2;
		for (vx2 = s->nxh; vx2 != NULL; vx2 = vx2->nxh) {
			int aa, bb, cc;		/* Probable hit check */
			int nnm, nmix;
	
//printf("~1 Removing deleted vertex no %d from net\n",vx2->no);
			if (vx2->del == 0) {	/* It's not really being deleted */
//printf("~1 vtx no %d is being retained\n",vx2->no);
				continue;
			}

			/* Remove from vertex net */
			for (j = 0; j < vx2->nnv; j++) {
				vx1 = vx2->nv[j];
	
//printf("~1 Removing vtx no %d from vtx no %d\n",vx2->no, vx1->no);
//				if (vx1->del == 0) { }		/* Speed optimization */
				{
					vtx_rem_vertex(s, vx1, vx2);
				}
//else printf("~1 Not removing from vtx no %d because it will be deleted anyway\n",vx1->no);
			}
			vx2->nnv = 0;

			/* Remove from parent nodes */
			for (e = 0; e <= di; e++) {
				int ix = vx2->nix[e];
				node *pp = s->n[ix];
				node_rem_vertex(s, pp, vx2);
			}
		}
	}

	/* Create/modify the vertex neighbour net lists for all the new vertexes. */
	/* Use a brute force search of local nodes to create the vertex net. */
	{
		vtx *vx1, *vx2;

		/* For each new vertx */
		for (i = 0; i < nncombs; i++) {
			vx1 = s->combs[i].vv;

			if (vx1 == NULL || vx1->del)
				continue;

			/* Possibly add other new vertexes as neighbours */
			for (j = i+1; j < nncombs; j++) { 

				vx2 = s->combs[j].vv;

				if (vx2 != NULL && vx2->del == 0)
					vtx_cnd_biadd_vtx(s, vx1, vx2, fixup);
			}

			/* Possibly add deleted and non-deleted vertexes */
			for (k = 0; k < s->combs[i].count; k++) {

#ifdef INDEP_SURFACE
				/* Add deleted vertex if it isn't going to be deleted */
				if (s->combs[i].v1[k]->del == 0) {
					vtx_cnd_biadd_vtx(s, vx1, s->combs[i].v1[k], fixup);
				}
#endif /* INDEP_SURFACE */

				/* Add non-deleted vertex */
				if (s->combs[i].v2[k]->del == 0)
					vtx_cnd_biadd_vtx(s, vx1, s->combs[i].v2[k], fixup);
			}

			/* Add any existing vertexes of the node we're re-adding */ 
			if (fixup) {
				for (j = 0; j < nn->nvv; j++) {
					if (nn->vv[j]->del)
						continue;
					vtx_cnd_biadd_vtx(s, vx1, nn->vv[j], fixup);
				}
			}
		}
	}

	/* Do second part of batch update */
	do_batch_update2(s, fixup);

#ifdef DEBUG
	printf("Add_to_vsurf done - added node %d\n",nn->ix);
#endif

	/* If we want intermediate fixup state: */
	/* dump_node_vtxs(s, 0); */
	/* sanity_check(s, 0); */

#ifdef NEVER
{
	vtx *vx;

	/* Dump vertex and associated vertex information */
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		printf("Vertex no %d has Vtx net:",vx->no);
		for (j = 0; j < vx->nnv; j++) {
			vtx *vx2 = vx->nv[j];
			printf(" %d",vx2->no);
		}
		printf("\n");
	}
	printf("\n");
	fflush(stdout);
}
#endif /* NEVER */

	return 1;
}

/* - - - - - - - - - - - */
/* Deal with verticies marked for deletion or addition, */
/* as well as updating the nodes consequently affects. */
/* If fixup is set, add any new or updates vertexes to the s->fchl */

/* Do the first part of the batch update */
static void do_batch_update1(ofps *s, int fixup) {
	int e, di = s->di;
	vtx *vv, *nvv;
	node *pp;

#ifdef DEBUG
	printf("Doing batch update to add/delete vertexes - 1\n");

#endif
	/* Update a vertexes vm, and decide whether it is going */
	/* to be deleted or just hidden. */
	for (vv = s->batch; vv != NULL; vv = vv->batch) {

#ifdef DEBUG
		printf("Pending vtx no %d del %d, add %d, vm %s |= %s &= %s\n",vv->no,vv->del,vv->add,psm(s,&vv->vm),psm(s,&vv->buvm),psm(s,&vv->bdvm));
		if (vv->ofake)
			error("An ofake vertex no %d was hit!\n",vv->ofake);
#endif

		if (vv->add == 1) {		/* New node */

#ifdef INDEP_SURFACE
			sm_or(s, &vv->vm, &vv->vm, &vv->buvm);
#ifdef DEBUG
			printf("Set vertex no %d vm to %s\n",vv->no,psm(s,&vv->vm));
#endif
#endif	/* INDEP_SURFACE */

			/* Add it to the cache */
			vtx_cache_add(s, vv);
	
			/* Add it to the spatial accelleration grid */
			ofps_add_vacc(s, vv);
	
			/* Add to seeding lists */
			ofps_add_vseed(s, vv);

		} else if (vv->add == 2) { /* Update the visibility setmask */
			int was_inseed = 0, is_inseed = 0;

#ifdef INDEP_SURFACE
		    if (sm_andtest(s, &s->sc[0].a_sm, &vv->vm) != 0)
				was_inseed = 1;

			sm_or(s, &vv->vm, &vv->vm, &vv->buvm);
#ifdef DEBUG
			printf("Updated vertex no %d vm to %s\n",vv->no,psm(s,&vv->vm));
#endif
		    if (sm_andtest(s, &s->sc[0].a_sm, &vv->vm) != 0)
				is_inseed = 1;

			if (vv->used == 0) {
				/* Adjust presense in eserr tree if visibility has changed */
				if (was_inseed && !is_inseed) {
//printf("Removing (1) vtx no %d, used %d, eserr %f, vm %s nsp %d\n",vv->no,vv->used,vv->eserr,psm(s,&vv->vm),vv->nsp);
					if ((aat_aerase(s->vtrees[vv->nsp], (void *)vv)) == 0)
						error("aat_aerase vertex failed to find vertex no %d (1)", vv->no);
				} else if (!was_inseed && is_inseed) {
//printf("Adding (1) vtx no %d, used %d, eserr %f, vm %s nsp %d\n",vv->no,vv->used,vv->eserr,psm(s,&vv->vm),vv->nsp);
					if ((aat_ainsert(s->vtrees[vv->nsp], (void *)vv)) == 0)
						error("aat_ainsert vertex malloc failed");
				}
			}
#endif	/* INDEP_SURFACE */

		} else if (vv->del != 0) {

#ifdef INDEP_SURFACE
			int was_inseed = 0, is_inseed = 0;

		    if (sm_andtest(s, &s->sc[0].a_sm, &vv->vm) != 0)
				was_inseed = 1;
//printf("Checked was_inseed %d for vtx no %d, used %d, eserr %f, vm %s nsp %d\n",was_inseed,vv->no,vv->used,vv->eserr,psm(s,&vv->vm),vv->nsp);

			/* Remove visibility due to any deletes */ 
			sm_andnot(s, &vv->vm, &vv->vm, &vv->bdvm);

		    if (sm_andtest(s, &s->sc[0].a_sm, &vv->vm) != 0)
				is_inseed = 1;
//printf("Checking is_inseed %d for vtx no %d, used %d, eserr %f, vm %s nsp %d\n",is_inseed,vv->no,vv->used,vv->eserr,psm(s,&vv->vm),vv->nsp);

			/* Adjust presense in eserr tree if visibility has changed */
			if (vv->used == 0 && was_inseed && !is_inseed) {
//printf("Removing (2) vtx no %d, used %d, eserr %f, vm %s nsp %d\n",vv->no,vv->used,vv->eserr,psm(s,&vv->vm),vv->nsp);
				if ((aat_aerase(s->vtrees[vv->nsp], (void *)vv)) == 0)
					error("aat_aerase vertex failed to find vertex no %d (2)", vv->no);
			}
#ifdef DEBUG
			printf("Delete vertex no %d vm to %s\n",vv->no,psm(s,&vv->vm));
#endif
			/* Don't delete vertex if it remains visible to some sub-surfaces. */
			if (sm_test(s, &vv->vm) != 0) {
				vv->del = 0;
				vv->add = 2;		/* Update it instead */
#ifdef DEBUG
				printf("Retaining vtx no %d marked for deletion because vm is %s\n",vv->no,psm(s,&vv->vm));
#endif
			}
#endif	/* INDEP_SURFACE */
		}
	}
}

/* Do the second part of the batch update */
static void do_batch_update2(ofps *s, int fixup) {
	int e, di = s->di;
	vtx *vv, *nvv;
	node *pp;

#ifdef DEBUG
	printf("Doing batch update to add/delete vertexes - 2\n");

#endif
	/* Add or delete a vertex */
	for (vv = s->batch; vv != NULL; vv = nvv) {

		nvv = vv->batch;
		vv->batch = NULL;		/* Ready for next time */
		vv->bch = 0;

		/* Setup vertex ready for another round */
		sm_set(s, &vv->buvm, 0);		/* Ready to OR in new visibility next time */
		sm_set(s, &vv->bdvm, 0);		/* Ready for OR in visibility to be remove next time */

		if (vv->del) {				/* delete vertex */

#ifdef DEBUG
			printf("Deleting vertex no %d\n",vv->no); fflush(stdout);
#endif
			/* Add all the parent nodes of this vertex to the update list */
			for (e = 0; e <= di; e++) {
				int ix = vv->nix[e];
				node *pp = s->n[ix];
	
				if (pp->upflag != s->flag) {
					pp->nup = s->nup;
					s->nup = pp;
					pp->upflag = s->flag;
				}
				/* During fixups, maintain nodes vertexes lists */
				/* (During re-seeding we update it as a batch) */
				if (fixup)
					node_rem_vertex(s, pp, vv);
			}
			del_vtx(s, vv);

#ifdef DEBUG
	{
		vtx *vx;
		int i, k;
	
		printf("~1 checking that no references to vertex remain after delete:\n");
		/* Check vertexes references to vertexes */
		for (vx = s->uvtx; vx != NULL; vx = vx->link) {
			for (k = 0; k < vx->nnv; k++) {
				if (vx->nv[k] == vv) {
					printf("Vertex 0x%x no %d still in no %d nv list after deletion\n",vv,vv->no,vx->no);
#ifdef WARNINGS
					warning("Vertex 0x%x no %d still in no %d nv list after deletion",vv,vv->no,vx->no);
#endif
				}
			}
		}
		/* Check nodes references to vertexes */
		for (i = 0; i < s->np; i++) {		/* For all nodes */
			node *p1 = s->n[i];
			for (k = 0; k < p1->nvv; k++) {	/* For all its vertexes */
				if (p1->vv[k] == vv) {
					printf("Vertex 0x%x no %d still in ix %d vv list after deletion\n",vv,vv->no,p1->ix);
#ifdef WARNINGS
					warning("Vertex 0x%x no %d still in ix %d vv list after deletion",vv,vv->no,p1->ix);
#endif
				}
			}
		}
		/* Fixup sorted list reference to vertex */
		for (i = 0; i < s->nsvtxs; i++) {
			if (s->svtxs[i] == vv) {
				printf("Vertex 0x%x no %d still in svtxs[%d] after deletion\n",vv,vv->no,i);
#ifdef WARNINGS
				warning("Vertex 0x%x no %d still in svtxs[%d] after deletion",vv,vv->no,i);
#endif
			} 
		}
	}
#endif /* DEBUG */
	
		} else if (vv->add != 0) {		/* New or updated vertex, set updates & checks */

#ifdef DEBUG
			printf("Adding vertex no %d\n",vv->no);
#endif

			/* Add all the parent nodes of this vertex to the update list */
			for (e = 0; e <= di; e++) {
				node *pp = s->n[vv->nix[e]];

				if (pp->upflag != s->flag) {
					/* Add node to update list */
					pp->nup = s->nup;
					s->nup = pp;
					pp->upflag = s->flag;
				}

				/* During fixups, maintain nodes vertexes lists. */
				/* (During re-seeding we update it as a batch) */
				if (fixup && vv->add == 1)
					node_add_vertex(s, pp, vv);
			}

			/* If this is a fixup and the vertex hasn't been added */
			/* to the "check" list, do so */
			if (fixup && vv->fflag != s->fflag) {
				vv->fchl = s->fchl; /* Add vertex to the "to be checked" list */
				if (s->fchl != NULL)
					s->fchl->pfchl = &vv->fchl;
				s->fchl = vv;
				vv->pfchl = &s->fchl;
				vv->fflag = s->fflag; 
#ifdef DEBUG
				printf("Adding vtx no %d to check list due to addition\n",vv->no);
#endif
			}
		}
	}

	s->batch = NULL;	/* Nothing in pending delete list */
	s->nup = NULL;		/* Nothing in nodes to be updated list */
}

/* ------------------------------------------------------------------------------- */

/* Do a quick count of the number of verticies hit by their */
/* neighbour nodes. This is used during itteration to decide */
/* whether to reseed or fixup. */
/* Return the number of vertexes hit */
static int
ofps_quick_check_hits(ofps *s) {
	int i, j, k, e, di = s->di;
	int nvxhits = 0;

	/* For all nodes */
	for (i = -s->nbp; i < s->np; i++) {
		node *nn = s->n[i];		/* Node being considered */

		s->flag++;				/* Marker flag for testing this node */
		nn->flag = s->flag;

		/* Check all the neighbors nodes */
		for (j = 0; j < nn->nvn; j++) {
			node *pp = s->n[nn->vn[j]];

			/* Test nn against all of pp's vertexes */ 
			for (k = 0; k < pp->nvv; k++) {
				vtx *vx = pp->vv[k];
				
				if (vx->cflag == s->flag)
					continue;			/* Don't test same node twice */
				vx->cflag = s->flag;

				/* If node that we're testing against is in vertex */
				/* ignore it, we expect them to hit. */
				for (e = 0; e <= di; e++) {
					if (nn->ix == vx->nix[e])
						break;
				}
				if (e <= di) {
					continue;
				}

#ifdef INDEP_SURFACE
				/* Check if this vertex is visible to this node */
				if (sm_vtx_node(s, vx, nn) == 0) {
					continue;	/* It's hidden */
				}
#endif	/* INDEP_SURFACE */

				if (nn->ix < 0) {
					pleq *vp = &s->gpeqs[-1 - nn->ix];
					double v = 0.0;
					
					/* See if the vertex is on the wrong side of the plane */
					for (v = vp->pe[di], e = 0; e < di; e++)
						v += vp->pe[e] * vx->p[e];

					if (v > 0.0) {
						nvxhits++;
					}
				} else {
					double eperr = ofps_comp_eperr7(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);

					/* See if the vertex eperr will be improved */
					if (eperr < (vx->eperr + 0.0)) {
						nvxhits++;
					}
				}
			}
		}
	}
	return nvxhits;
}

/* ------------------------------------------------------------------------------- */

/* Recursive hit search routine: */

/* Test a vertex for a possible hit from a node. */
/* Recursively search dorec recursions until there is at least */
/* one hit, and stop after beyhit recursions beyond the hit region. */
/* s->flag is assumed to be relevant for the given node. */
/* Any hit vertexes are added to the s->nxh list. */
/* s->vvchecks and s->nvcheckhits will be updated. */
/* Return nz if vertex was hit */

/* Breadth first search. */
/* (Assumes s->flag has been incremented) */
/* [We're assuming that there is a single connected hit region, which */
/*  is not valid for SUBD. ] */
static int
ofps_check_vtx(ofps *s, node *nn, vtx *vx, int dorec, int beyhit) {
	int i, j, e, di = s->di;
	vtx *slist = NULL;			/* Next to search list */
	int dist;					/* Distance from initial vertex */
	int hit = 0;
	double tol = 0.0;			/* Tollerance */

#ifdef DEBUG
	printf("ofps_check_vtx() for node ix %d starting at vertex no %d, dorec %d, beyhit %d\n",nn->ix,vx->no,dorec,beyhit);
#endif

#ifdef SANITY_CHECK_HIT
	if (nn->ix < -s->nbp)
		error("Calling ofps_check_node on fake outside node");
#endif

	if (vx->cflag == s->flag)
		return vx->del;		/* Already been checked */

	/* Put the starting node on the search list */
	vx->slist = slist;
	slist = vx;
	vx->sflag = s->flag;	/* Mark as done for pre-hit search */

	/* until we run out of vertexes, or we are done */
	for (dist = 0; slist != NULL && dist <= dorec; dist++) {
		vtx *nvx;

		/* For each vertex in the search list, check it and recursion. */
		for (vx = slist, slist = NULL; vx != NULL; vx = nvx) {
			nvx = vx->slist;
			vx->opqsq = 0;		/* Not on the list anymore */

			if (vx->ofake)
				continue;		/* ofake vertexes can't be hit */
#ifdef DEBUG
			printf("%d: Checking vtx no %d %s\n",dist, vx->no,hit ? "Post-Hit" : "Pre-Hit");
#endif
#ifdef INDEP_SURFACE
			/* Only check for hit if the vertex is visible to the node */
			if (sm_vtx_node(s, vx, nn) == 0) {
# ifdef DEBUG
				printf("%d: Vertex no %d xmask 0x%x vm %s isn't visible to ix %d pmask 0x%x a_sm %s\n",dist,vx->no,vx->cmask,psm(s,&vx->vm),nn->ix,nn->pmask,psm(s,&s->sc[nn->pmask].a_sm));
# endif	/* DEBUG */
				continue;
			}
#endif	/* INDEP_SURFACE */

			/* If the vertex hasn't been checked yet: */
			if (vx->cflag != s->flag) {
				vx->add = 0;
				vx->del = 0;
				vx->par = 0;

				s->vvchecks++;			/* Checking a vertex */

				/* Check if node is already parent to this vertex. */
				/* This only happens during fixups if the reposition fails and we */
				/* retain the vertex with the deleted vertex location (not currently */
				/* done), or by slim numerical margine, so ignore such hits. */
				/* We treat a parent as a hit node for the purposes of recursion, */
				/* and add it to a special list used to complete the vertex net. */
				if (nn->ixm & vx->nix[MXPD+2]) {	/* Is in nixm */
					for (e = 0; e <= di; e++) {		/* Do exact check */
						if (nn->ix == vx->nix[e])
							break;
					}
					if (e <= di) {
#ifdef DEBUG
						printf("Vertex no %d has already got node ix %d\n",vx->no,nn->ix);
#endif
						vx->par = 1;
					}
				}

				if (nn->ix < 0) {
					pleq *vp = &s->gpeqs[-1 - nn->ix];
					double v = 0.0;
					
					/* See if the vertex is on the wrong side of the plane */
					for (v = vp->pe[di], e = 0; e < di; e++)
						v += vp->pe[e] * vx->p[e];

					if (!vx->par && v > tol) {
						s->nvcheckhits++;
						if (!hit)
							slist = nvx = NULL;		/* Abort pre-hit search */
						hit = 1;
						vx->del = 1;		/* Mark for deletion */
						vx->nxh = s->nxh;	/* Add vertex to list */
						s->nxh = vx;
						vx->disth = 0;		/* This is a hit vertex */
						vx->hflag = s->flag;
#ifdef DEBUG
						printf("%d: Gamut surface boundary plain hit by %f\n",dist,v);
#endif
					}
#ifdef DEBUG
					  else {			/* If worse */
						printf("%d: Gamut surface boundary plain miss by %f\n",dist,v);
					}
#endif
				} else {	/* Node rather than boundary plane */

					/* nba_eperr is assumed to be valid if vx->cflag == s->flag */
					vx->nba_eperr = ofps_comp_eperr7(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);
#ifdef DEBUG
					printf("%d: Computing nba_eperr of %f for vtx no %d\n",dist, vx->nba_eperr, vx->no);
#endif
					/* See if the vertex eperr will be improved */
					if (!vx->par && (vx->eperr - vx->nba_eperr) > tol) {
						s->nvcheckhits++;
						if (!hit)
							slist = nvx = NULL;		/* Abort pre-hit search */
						hit = 1;
						vx->del = 1;		/* Mark for deletion */
						vx->nxh = s->nxh;	/* Add vertex to list */
						s->nxh = vx;
						vx->disth = 0;		/* This is a hit vertex */
						vx->hflag = s->flag;
#ifdef DEBUG
						printf("%d: Vertex error improvement hit by %f (%f < %f)\n",dist, vx->eperr-vx->nba_eperr,vx->nba_eperr,vx->eperr);

						if (vx->par) {
							printf("Vertex no %d hit by its own parent ix %d\n",vx->no, nn->ix);
#ifdef WARNINGS
							warning("Vertex no %d hit by its own parent ix %d",vx->no, nn->ix);
#endif
						}
#endif
					}
#ifdef DEBUG
					  else {			/* If worse */
						printf("%d: Vertex error not hit by %f (%f < %f)\n",dist, vx->eperr-vx->nba_eperr,vx->nba_eperr,vx->eperr);
					}
#endif
				}
				vx->cflag = s->flag;
// ~~777
//if (vx->del && i_rand(1,1000) == 15) {
//	printf("~1 failing to check vertex no %d\n",vx->no);
//	vx->cflag = s->flag -1;
//	vx->del = 0;
//}
			}

			/* Decide whether to recurse by adding vertexes to the new list */
			if (!hit) {

				/* Pre-hit recursion */
				if (dist < dorec) {		/* Still within search radius */

					/* Add all the unsearched vertexes neighbors to the next search list */
					for (j = 0; j < vx->nnv; j++) {
						vtx *vx2 = vx->nv[j];
		
						if (vx2->sflag == s->flag)	/* Already been pre-hit searched */
							continue;
#ifdef DEBUG
						printf("%d: Adding vtx no %d to next pre-hit search list\n",dist, vx2->no);
#endif
						/* Put the neighbour node on the search list */
						vx2->slist = slist;
						slist = vx2;
						vx2->sflag = s->flag;
					}
				}

			} else {

				/* Post hit recursion */
				if (vx->disth <= beyhit) {		/* Still within post-hit search radius */
					int disth = vx->disth + 1;	/* Neighbours distance */

					/* Add all the unsearched vertexes neighbors to the next search list */
					for (j = 0; j < vx->nnv; j++) {
						vtx *vx2 = vx->nv[j];
		
						if (vx2->hflag == s->flag) {	/* Already been post-hit searched */
							if (disth >= vx2->disth)
								continue;				/* So skip it */

							/* The already post-hit searched neighbour has an improved distance */ 
#ifdef DEBUG
							printf("%d: Improving ph searched vtx %d disth from %d to %d\n",dist, vx2->no,vx2->disth,disth);
#endif
							vx2->disth = disth;		/* Improved distance to hit though */
							if (vx2->disth > beyhit || vx2->opqsq)
								continue;		/* But it's still too far, or already on the list */
							/* Search this neighbour again now that it is within radius */
						} else {
							vx2->disth = disth;				/* Set hit distance */
						}
#ifdef DEBUG
						printf("%d: Adding vtx no %d to next post-hit search list (disth = %d)\n",dist, vx2->no,vx2->disth);
#endif
						/* Put the neighbour node on the search list */
						vx2->slist = slist;
						slist = vx2;
						vx2->sflag = vx2->hflag = s->flag;
						vx2->opqsq = 1;			/* On the list */
					}
				}
#ifdef DEBUG
				  else
					printf("%d: Vertex %d disth %d is > beyhit %d so not recursing\n",dist, vx->no,vx->disth,beyhit);
#endif

			}
		}	/* Next vertex in current list */
#ifdef DEBUG
		printf("Finished inner loop because vx 0x%x = NULL\n",vx);
#endif
	}	/* Next list */
#ifdef DEBUG
	printf("Finished outer loop because slist 0x%x = NULL, || dist %d > dorec %d\n",slist,dist,dorec);
#endif

	return hit;
}

/* Non-recursive version of above used for sanity checking */
/* that doesn't set any flags on vx */
static int
ofps_check_vtx_sanity(ofps *s, node *nn, vtx *vx, int fixit) {
	int i, j, e, di = s->di;
	vtx *slist = NULL;			/* Next to search list */
	int dist;					/* Distance from initial vertex */
	double tol = 1e-6;
	int hit = 0;
	int par = 0;

#ifdef DEBUG
	printf("ofps_check_vtx_sanity() for node ix %d and vertex no %d\n",nn->ix,vx->no);
#endif
	if (vx->cflag == s->flag) {
#ifdef DEBUG
		printf("Returning alread calculated del = %d\n",vx->del);
#endif
		return vx->del;		/* Already been checked */
	}

	if (vx->ofake) {	/* ofake nodes can't be hit */
#ifdef DEBUG
		printf("Returning ofake del = 0\n");
#endif
		return 0;
	}

#ifdef INDEP_SURFACE
	/* Only check for hit if the vertex is visible to the node */
	if (sm_vtx_node(s, vx, nn) == 0) {
#ifdef DEBUG
		printf("Returning non-visible del = 0\n");
#endif
		return 0;
	}
#endif	/* INDEP_SURFACE */

	/* Check if node is already parent to this vertex. */
	/* This only happens during fixups if the reposition fails and we */
	/* retain the vertex with the deleted vertex location (not currently */
	/* done), or by slim numerical margine, so ignore such hits. */
	/* We treat a parent as a hit node for the purposes of recursion, */
	/* and add it to a special list used to complete the vertex net. */
	if (nn->ixm & vx->nix[MXPD+2]) {	/* Is in nixm */
		for (e = 0; e <= di; e++) {		/* Do exact check */
			if (nn->ix == vx->nix[e])
				break;
		}
		if (e <= di)
			par = 1;
	}

	if (nn->ix < 0) {
		pleq *vp = &s->gpeqs[-1 - nn->ix];
		double v = 0.0;
		
		/* See if the vertex is on the wrong side of the plane */
		for (v = vp->pe[di], e = 0; e < di; e++)
			v += vp->pe[e] * vx->p[e];

		if (!par && v > tol) {
			hit = 1;

			if (fixit) {
				vx->slist = slist;
				slist = vx;
				vx->sflag = s->flag;
			}
		}
	} else {	/* Node rather than boundary plane */
		double nba_eperr;

		nba_eperr = ofps_comp_eperr7(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);

		/* See if the vertex eperr will be improved */
		if (!par && (vx->eperr - nba_eperr) > tol) {
			hit = 1;
			if (fixit) {
				vx->slist = slist;
				slist = vx;
				vx->sflag = s->flag;
			}
		}
	}

#ifdef DEBUG
	printf("Returning computed del = %d\n",hit);
#endif
	return hit;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

/* Add a node to the currnent voronoi. */
/* Return nz if the addition fails due to there being no vetex hits or a cooincince. */
/* Return nz if abortonfail is set and we fail to position the node. */
/* (This theoretically shouldn't happen, but does, due to the perceptual */
/*  geometry ?) */ 
static int add_node2voronoi(
ofps *s,
int poi,		/* Index of sample point to update/create Voronoi surface */
int abortonfail	/* 0 = ignore position failures, 1 = abort add if there are any failures */
) {
	node *nn = s->n[poi];		/* Node in question */
	int e, di = s->di;
	int i, j;
	vtx *vx = NULL;			/* Closest vertex */

#ifdef DEBUG
	printf("\nAdding Node ix %d pmask 0x%x at %s (perc %s) to Voronoi surface\n",poi,nn->pmask,ppos(di,nn->p),ppos(di,nn->v));
#endif

	if (poi < 0)
		error("Attempt to add fake point to veronoi surface");

	if (nn->nvv > 0)
		error("ofps: assert, node vertex info should be empty on add_node2voronoi() entry");

	for (i = 0; i < 20; i++) {
		int pci;		/* Point cell list index */
		acell *cp;		/* Acceleration cell */
		node *pp;

		/* Check if by some misfortune, this node colides with an existing node. */
		pci = ofps_point2cell(s, nn->v, nn->p);	/* Grid index of cell of interest */
		cp = &s->grid[pci];
		for (pp = cp->head; pp != NULL; pp = pp->n) {
			for (e = 0; e < di; e++) {
				if (fabs(nn->v[e] - pp->v[e]) > (COINTOL * 100.0))
					break;		/* Not cooincident */
			}
			if (e >= di) { 	/* Cooincident */
#ifdef DEBUG
				printf("Node oint collides with existing - joggling it\n");
//				warning("Node oint collides with existing - joggling it");
#endif
				/* Joggle it's position */
				for (e = 0; e < di; e++) {
					if (nn->p[e] < 0.5)
						nn->p[e] += d_rand(0.0, 1e-4);
					else
						nn->p[e] -= d_rand(0.0, 1e-4);
				}
				/* Ignore confine planes. Next itter should fix it anyway ? */
				ofps_clip_point6(s, nn->p, nn->p);
				s->percept(s->od, nn->v, nn->p);
				break;
			}
		}
		if (pp == NULL)
			break;
	}
	if (i >= 20) {
		if (s->verb > 1)
			warning("add_node2voronoi: Assert, was unable to joggle cooincindent point");
		return 1;
	}

#ifdef DEBUG
	printf("Locating all the hit vertexs\n");
#endif

	s->nvcheckhits = 0;	/* Count number of vertexes hit by recursive check. */
	s->batch = NULL;	/* Nothing in pending delete list */
	s->nup = NULL;		/* Nothing in nodes to be updated list */
	s->flag++;			/* Marker flag for adding this node */
	s->nxh = NULL;		/* Nothing in nodes hit list */

#ifdef DEBUG
	printf("Done check of vertexes for hits\n");
#endif

	/* Number of nodes that would be checked by exaustive search */
	s->vvpchecks += s->nv;

	ofps_findhit_vtxs(s, nn);

#ifdef SANITY_CHECK_HIT
#ifdef DEBUG
	printf("Doing sanity check of hits\n");
#endif
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {		/* Check all vertexes */

		if (vx->cflag != s->flag && ofps_check_vtx_sanity(s, nn, vx, 0)) {
			warning("!!!!!! Sanity: Add hit missed vertex no %d at %s !!!!!!",vx->no,ppos(di,vx->p));
			printf("!!!!!! Sanity: Add hit missed vertex no %d at %s !!!!!!\n",vx->no,ppos(di,vx->p));
			/* Don't stop for out of gamut vertexes that would have been hit */
			if (ofps_would_clip_point(s, vx->p))
				continue;

			/* Check if any of it's neighbours have been checked. */
			for (j = 0; j < vx->nnv; j++) {
				vtx *vx2 = vx->nv[j];

				if (vx2->cflag == s->flag)
					break;				/* Yes */
			}
			if (j >= vx->nnv) {
				warning("!!!!!! Sanity: Missed vertex was in isolated region");
			} else {
				warning("!!!!!! Sanity: Missed vertex was adjacent to no %d", vx->nv[j]->no);
			}
#ifdef SANITY_CHECK_HIT_FATAL
			error("Failed to locate all hit vertexes");
#endif /* SANITY_CHECK_HIT_FATAL */
		}
	}
#endif	/* SANITY_CHECK_HIT */
			
#ifdef DEBUG
	printf("There were %d vertexes that will be hit by adding node\n",s->nvcheckhits);
#endif

	if (s->nvcheckhits == 0) {
			if (s->verb > 1)
			warning("Failed to get any vertex hits when adding a new node ix %d at %s",nn->ix,ppos(di,nn->p));
		return 1;
	}

	/* Now turn all the hit vertexes into new vertexes. */
	if (add_to_vsurf(s, nn, 0, abortonfail) > 0) {
		s->add_hit++;
	} else {
		if (abortonfail)
			return 1;
		s->add_mis++;
	}

	ofps_add_nacc(s, nn);			/* Add to spatial accelleration grid */

	s->np++;

#ifdef DEBUG
	printf("Done add_node2voronoi()\n");
#endif

	return 0;
}

/* ------------------------------------------------------------------------------- */

/* Given a list of di plane equations, */
/* compute the intersection point. */
/* return nz if there is no intersection */
static int comp_vtx(ofps *s, double *p, pleq **peqs) {
	int i, e, di = s->di;
	double **ta, *TTA[MXPD], TA[MXPD][MXPD];

	for (e = 0; e < di; e++)
		TTA[e] = TA[e];
	ta = TTA;

	for (i = 0; i < di; i++) {
		for (e = 0; e < di; e++)
			ta[i][e] = peqs[i]->pe[e];	/* Plane normal becomes row of matrix */
		p[i] = -peqs[i]->pe[di];		/* Plane constant becomes target */
	}
	/* Solve the simultaneous linear equations A.x = B */
	/* Return 1 if the matrix is singular, 0 if OK */
	if (polished_solve_se(ta, p, di))
		return 1;

	return 0;
}

/* --------------------------------------------------- */

/* Use a brute force search to (re-)create the vertex net. */
/* This is used in initialization. */
static void create_vtx_net(ofps *s) {
	int ff, f, e, di = s->di;
	vtx *vx1, *vx2;
	
#ifdef DEBUG
	printf("Doing create_vtx_net\n");
#endif
	
	/* For each vertx */
	for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link) {

		vx1->nnv = 0;		/* Clear the current list */
#ifdef DEBUG
		printf("Creating neighbourhood net for vtx no %d\n",vx1->no);
#endif

		/* Search all other vertexes for neighbours */
		for (vx2 = s->uvtx; vx2 != NULL; vx2 = vx2->link) {
			int aa, bb, cc;		/* Probable hit check */
			int nnm, nmix;

//printf("~1 checking against vtx %d\n",vx2->no);
			if (vx1 == vx2) {
//printf("~1 skip because it's the same\n");
				continue;
			}


			/* Use the nixm to quickly check if all but one parent node matches */
			aa = vx1->nix[MXPD+2];	/* nixm */
			bb = vx2->nix[MXPD+2];	/* nixm */
			if ((aa & bb) == 0 || (cc = aa & ~bb, (cc & (cc-1)) != 0)) {
//printf("~1 skip because nixm 0x%x and 0x%x don't match\n",aa,bb);
				continue;		/* It's certainly not */
			}

			/* Do an exact check of all except one node match */
			for (nnm = ff = e = 0; e <= di; e++) {
				for (f = ff; f <= di; f++) {
					if (vx1->nix[e] == vx2->nix[f]) {
						ff = f;			/* Start from here next time */
						break;
					}
					if (vx1->nix[e] > vx2->nix[f])	/* No point in looking further */
						f = di;
				}
				if (f > di) {	/* Didn't match */
					if (++nnm > 1)
						break;
					nmix = e;
				}
			}
			if (e <= di) {
//printf("~1 skip because nix %s and %s aren't one different\n",pcomb(di,vx1->nix),pcomb(di,vx2->nix));
				continue;			/* No match */
			}
			
			if (nnm == 0) {
				error("ofps: two vertexes have the same nodes !\n"
					  "no %d at %s nix %s\nno %d at %s nix %s",
				vx1->no,ppos(di,vx1->p),pcomb(di,vx1->nix),
				vx2->no,ppos(di,vx2->p),pcomb(di,vx2->nix));
			}

			/* vx2 is a neighbour, so add it to the vtx net */
			vtx_add_vertex(s, vx1, vx2);
//printf("~1 brute force: adding vtx %d as neighbour to %d\n",vx2->no,vx1->no);
		}
	}
}

/* --------------------------------------------------- */

/* Use a brute force search to discover all the valid */
/* sub-surface combinations. */
static void discover_subsuf(ofps *s) {
	int co;
	double p[MXPD];
	pleq *peqs[MXPD];
	int i, j, k, e, di = s->di;
	setmask acm;		/* Accumulated mask for special last entry */

	if (s->sminit)
		return;		/* Do this once */

#ifdef DEBUG
	printf("Computing subd face combinations\n");
#endif

	if ((s->sc = (surfcomb *)calloc(sizeof(surfcomb), (1 << s->nbp))) == NULL)
		error ("ofps: malloc failed on sufcomb array");

	for (co = 0; co < (1 << s->nbp); co++) {
		
		s->sc[co].co = co;

		/* Count number of planes */
		for (i = e = 0; e < s->nbp; e++) {
			if (co & (1 << e))
				i++;
			/* Skip combo if odd and even dimension planes are set */
			if ((e & 1) == 0 && e < (2 * di) && (co & (1 << e)) && (co & (1 << (e+1))))
				break;
		}
		s->sc[co].nos = i;
		if (i > di || e < s->nbp) {
			s->sc[co].valid = 0;
			continue;
		}

		/* Check that the combination results in a valid */
		if (i == di) {
			for ( j = e = 0; e < s->nbp; e++) {
				if (co & (1 << e))
					peqs[j++] = &s->gpeqs[e];
			}
			if (comp_vtx(s, p, peqs) != 0 || ofps_would_clip_point(s, p)) {
				s->sc[co].valid = 0;
			} else {
				s->sc[co].valid = 1;
			}
		} else {
			if (co == 0 || i == 1) {
				s->sc[co].valid = 1;
			} else {
				s->sc[co].valid = -1;
			}
		}
//printf("~1 val %s sc[%d].valid = %d\n",icmPdv(di, p), co,s->sc[co].valid);
	}
	/* Go through the unknown combinations, and see if there */
	/* is a valid lower dimensional combination that is valid. */
	for (co = 0; co < (1 << s->nbp); co++) {
		if (s->sc[co].valid == -1) {
			for (i = co+1; i < (1 << s->nbp); i++) {
				if ((i & co) == co && s->sc[i].valid == 1) {
					s->sc[co].valid = 1;
					break;
				}
			}
			if (i >= (1 << s->nbp))		/* Failed to find a valid combination */
				s->sc[co].valid = 0;
		}
	}
#ifdef USE_DISJOINT_SETMASKS
	/* We can reduce the number of setmask bits by figuring out which */
	/* combinations are disjoint, and using the same setmask bits for disjoint */
	/* combinations. For CMYK, this reduces the setmask from 80-100 to less than 32 bits, */
	/* permiting faster mask manipulation. */ 
	{
		surfcomb *scp, *zd = NULL;	/* Zero Dimension combinations */
		surfcomb *sets = NULL;		/* Sets at a given nos */
		int nsets = 0;				/* Current number of sets */
		int _nsets = 0;				/* Allocated array size */
		int nos;					/* Number of surfaces */

		/* init the circular lists, and add the 0D points to their list */
		for (k = co = 0; co < (1 << s->nbp); co++) {
			s->sc[co].ds = &s->sc[co];	/* Init circular list to itself */
			if (s->sc[co].valid == 0)
				continue;
			/* Create a list of 0D points and count them */
			if (s->sc[co].nos == di)  {
				k++;
				if (zd == NULL)
					zd = &s->sc[co];
				else {
					s->sc[co].ds = zd->ds;
					zd->ds = &s->sc[co];
				}
			}
		}

		if (zd == NULL)
			error("No zero-dim surface combinations (s->nbp = %d)",s->nbp);

//printf("~1 total 0D points = %d\n",k);

		/* Temporarily use the setmask to track 0D hits */
		sm_init(s, k);
	
		k = 2;		/* Count total disjoint sets, including 2 for di D and 0 D */

		/* Locates sets for each dimension level */
		for (nos = 1; nos < di; nos++) {
			nsets = 0;

//printf("~1 doing nos = %d\n",nos);
			/* Add the next combination to the sets */
			for (co = 0; co < (1 << s->nbp); co++) {
				if (s->sc[co].valid == 0 || s->sc[co].nos != nos)
					continue;

//printf("~1 checking combo 0x%x\n",co);
				/* Figure out 0D hits on this combo */
				i = 0;
				scp = zd;
				do {
					if ((co & scp->co) == co)
						sm_setbit(s, &s->sc[co].i_sm, i, 1);
					i++;
					scp = scp->ds;
				} while(scp != zd);
//printf("~1 combo 0x%x has hits %s\n",co,psm(s,&s->sc[co].i_sm));

				/* Search through the existing sets, and see */
				/* if this combo is disjoint */
				for (j = 0; j < nsets; j++) {
					setmask tsm;
					
					if (sm_and(s, &tsm, &sets[j].i_sm, &s->sc[co].i_sm) == 0) {
						/* Add this combo to the existing set */

//printf("~1 adding to set %d\n",j);
						s->sc[co].ds = sets[j].ds->ds;
						sets[j].ds->ds = &s->sc[co];
						sm_or(s, &sets[j].i_sm, &sets[j].i_sm, &s->sc[co].i_sm);
						break;
					}
//else printf("Miss on set %d hits %s, AND %s\n",j,psm(s,&sets[j].i_sm),psm(s,&tsm));
				}
				/* If we can't use an existing set, create a new one */
				if (j >= nsets) {
					if (nsets >= _nsets) {
						_nsets = 2 * _nsets + 5;
						if ((sets = (surfcomb *)realloc(sets, sizeof(surfcomb) * _nsets)) == NULL)
							error("malloc failed on disjoint sets size %d", _nsets);
					}
					sm_cp(s, &sets[j].i_sm, &s->sc[co].i_sm);	/* Hits to this set */
					sets[j].ds = &s->sc[co];	/* Only entry in circular list */
//printf("New set %d hits %s\n",j,psm(s,&sets[j].i_sm));
					nsets++;
					k++;
				}
			}
		}

		if (sets != NULL)
			free(sets);

#ifdef DEBUG
		printf("Total number of setmask disjoint sets = %d\n",k);
#endif

		/* Setup the setmask params */
		sm_init(s, k);
	
		/* Assign the individual setmask bits */
		for (i = co = 0; co < (1 << s->nbp); co++) {
			if (s->sc[co].valid == 0 || s->sc[co].smset == 1)
				continue;

//printf("~1 setting mask bit on comb 0x%x and its set\n",co);
			/* Assign setmask to all in this set */
			scp = &s->sc[co];
			do {
//printf("~1 setting mask bit %d on comb 0x%x\n",i,scp->co);
				sm_set(s, &scp->i_sm, 0);			/* Clear temporary hit mask */
				sm_setbit(s, &scp->i_sm, i, 1);
				scp->smset = 1;
				scp = scp->ds;
			} while(scp != &s->sc[co]);
			i++;
		}

	}
#else /* !USE_DISJOINT_SETMASKS */

	/* Count the number of valid combinations */
	for (i = co = 0; co < (1 << s->nbp); co++) {
		if (s->sc[co].valid == 0)
			continue;
		i++;
	}
#ifdef DEBUG
	printf("Total number of setmask sets = %d\n",i);
#endif

	/* Setup the setmask params */
	sm_init(s, i);

	/* Assign the individual setmask bits */
	for (i = co = 0; co < (1 << s->nbp); co++) {
		if (s->sc[co].valid == 0)
			continue;
		sm_setbit(s, &s->sc[co].i_sm, i, 1);
		i++;
	}
#endif /* !USE_DISJOINT_SETMASKS */

	sm_set(s, &acm, 0);		/* Init overall accumulated mask */

	/* Compute the accumulated setmask bits */
	for (i = 0; i < (1 << s->nbp); i++) {
		if (s->sc[i].valid == 0)
			continue;
		for (j = 0; j < (1 << s->nbp); j++) {
			if ((i & j) != j || s->sc[j].valid == 0)
				continue;
			sm_or(s, &s->sc[i].a_sm, &s->sc[i].a_sm, &s->sc[j].i_sm);
		}
		sm_or(s, &acm, &acm, &s->sc[i].i_sm);
	}

	/* Set special "all planes, all valid" combination as the last */
	/* entry for use by fake surface nodes. */
	s->sc[(1 << s->nbp)-1].valid = 1;
	s->sc[(1 << s->nbp)-1].nos = s->nbp; 
	sm_cp(s, &s->sc[(1 << s->nbp)-1].i_sm, &acm);
	sm_cp(s, &s->sc[(1 << s->nbp)-1].a_sm, &acm);
	s->sc[(1 << s->nbp)-1].smset = 1;
	s->sc[(1 << s->nbp)-1].ds = NULL;

#ifdef MAXINDEP_2D
	/* Go through the combinations and invalidate any */
	/* that are not full-d or more than 2D */
	for (co = 0; co < (1 << s->nbp); co++) {
		if (s->sc[co].valid) {
//			if (s->sc[co].nos != 0 && (di - s->sc[co].nos) > 1) // test in 3D
			if (s->sc[co].nos != 0 && (di - s->sc[co].nos) > 2)
				s->sc[co].valid = 0;
		}
	}
#endif /* MAXINDEP_2D */

#ifdef DEBUG
	/* Print diagnostics */
	for (i = 0; i < (1 << s->nbp); i++) {
		if (s->sc[i].valid == 0)
			continue;
		printf(" Mask 0x%x, setmasks i = %s, a = %s\n",i,psm(s,&s->sc[i].i_sm),psm(s,&s->sc[i].a_sm));
	}
#endif

	s->sminit = 1;
}
/* --------------------------------------------------- */

/* Compute a simple but unbounded model of the */
/* perceptual function. We use the current vertex values */
/* to setup the model */
/* (It would be faster to do the optimization per output channel!) */

/* Matrix optimisation function handed to powell() */
static double xfitfunc(void *edata, double *x) {
	ofps *s = (ofps *)edata;
	int e, di = s->di;
	double rv = 0.0;
	vtx *vx;

	/* For all the vertexes */
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		double v[MXPD], ev;

		/* Apply matrix cube interpolation */
		icxCubeInterp(x, di, di, v, vx->p);

		/* Evaluate the error */
		for (ev = 0.0, e = 0; e < di; e++) {
			double tt;
			tt = vx->v[e] - v[e];
			ev += tt * tt;
		}
		rv += ev;
	}

//	printf("~1 rv = %f\n",rv);

	return rv;
}

/* Fit the unbounded perceptual model to just the inside vertexes */
static void init_pmod(ofps *s) {
	int e, di = s->di;
	double sa[MXPD * (1 << MXPD)];
	double rerr;

	/* Setup matrix to be closest values initially */
	for (e = 0; e < (1 << di); e++) {	/* For each colorant combination */
		int j, f;
		double bdif = 1e6;
		double ov[MXPD];
		vtx *vx, *bvx = NULL;
	
		/* Search the vertex list to find the one closest to this input combination */
		for (vx = s->uvtx; vx != NULL; vx = vx->link) {
			double dif = 0.0;

			if (vx->ofake)
				continue;		/* Ignore outside vertexes */

			for (j = 0; j < di; j++) {
				double tt;
				if (e & (1 << j))
					tt = s->imax[j] - vx->p[j];
				else
					tt = s->imin[j] - vx->p[j];
				dif += tt * tt;
			}
			if (dif < bdif) {		/* best so far */
				bdif = dif;
				bvx = vx;
				if (dif < 0.001)
					break;			/* Don't bother looking further */
			}
		}
		for (f = 0; f < di; f++)
			 s->pmod[f * (1 << di) + e] = bvx->v[f];
	}

	for (e = 0; e < (di * (1 << di)); e++)
		sa[e] = 10.0;

	if (powell(&rerr, di * (1 << di), s->pmod, sa, 0.001, 1000,
	                                    xfitfunc, (void *)s, NULL, NULL) != 0) {
		if (s->verb > 1)
			warning("Powell failed to converge, residual error = %f",rerr);
	}

#ifdef DEBUG
	printf("Perceptual model fit residual = %f\n",sqrt(rerr));
#endif
	s->pmod_init = 1;
}

/* --------------------------------------------------- */
/* Init fake node contents, and setup the initial */
/* voronoi surface with the first node. */
static void ofps_binit(ofps *s) {
	int e, di = s->di;
	int doink = 0;
	int i, j;
	DCOUNT(co, MXPD, di, 0, 0, 2);	/* Count through corner verticies */
	int iix = -2 * di - 2;  /* Fake inside node index */
	node *inp = s->n[iix];	/* Fake inside node */
	int oix = -2 * di - 3;	/* Fake outside node index */
	node *onp = s->n[oix];	/* Fake outside node */
	double ivtx_whts[] = { /* Initial vertex weightings */
		3.7144267283692024e+165,
		1.3997102851752585e-152,
		6.1677886722367450e+223,
		1.7281009363426126e+097,
		2.0087766625640005e-139,
		4.9406564584124654e-323,
		7.7791723264315535e-260,
		8.5733372291341995e+170,
		6.0046007797559735e-067,
		2.8214561724952793e+243,
		5.0132438738338732e+262,
		1.6259745436952323e-260,
		7.9968034958246946e+001
	};
	double vtxwt;			/* Combined initial vertexnode weighting */
	unsigned int fullmask = 0;

#ifdef DEBUG
	printf("Binit called\n");
#endif
	if (s->ilimit < (double)di) 	/* Ink limit is active */
		doink = 1;

	/* Init fake inside and outside node */
	inp->ix = iix;
	inp->fx = 1;
	inp->pmask = 0;
	onp->ix = oix;
	onp->fx = 1;
	onp->pmask = 0;

	for (i = 0; i < (2 * di); i++)
		fullmask |= 1 << i;
	if (doink)
		fullmask |= 1 << i;

	/* Init the axis aligned gamut surface plane equations */
	/* and also setup nodes that are indexes by the fake indexes */
	/* with just the information that will be used. */ 
	for (i = 0; i < (2 * di); i++) {	/* unit cell at 0 */
		int ii = i >> 1;				/* Dimension */
		int ix;					/* Surface "node" index */
		pleq *vp;				/* plane being initialized */
		node *np;

		ix = -i-1;						/* -1 to -2di fake other nodes */
		vp = &s->gpeqs[-1-ix];			/* Pointer to plane associated with fake node */
		vp->ix = ix;

		for (e = 0; e < di; e++)
			vp->pe[e] = 0.0;
		vp->pe[ii] = i & 1 ? 1.0 : -1.0;					/* Normal */
		vp->pe[di] = i & 1 ? -s->imax[ii] : s->imin[ii];	/* Constant */

		np = s->n[ix];
		np->ix = ix;
		np->fx = 1;		/* They don't move */
		np->nsp = 1;
		np->sp[0] = vp;
		np->pmask = fullmask; 
//		np->pmask = 1 << i;		/* fake surface node pmask is itself for cmask ?? */
//		onp->pmask = inp->pmask |= np->pmask;
	}
	s->nbp = 2 * di;		/* Number of boundary planes */

	/* Add ink limit surface plane and its fake node */
	if (doink) {	/* Ink limit plane is orthogonal to diagonal */
		int ix;					/* Surface "node" index */
		pleq *vp;				/* plane being initialized */
		node *np;
		double len;

		ix = -i-1;						/* -1 to -2di fake other nodes */
		vp = &s->gpeqs[-1-ix];			/* Pointer to plane associated with fake node */
		vp->ix = ix;
		len = 1.0/sqrt((double)di);		/* Normalised length */
		for (e = 0; e < di; e++)
			vp->pe[e] = len;
		vp->pe[di] = -s->ilimit * len;

		np = s->n[ix];
		np->ix = ix;
		np->fx = 1;		/* They don't move */
		np->nsp = 1;
		np->sp[0] = vp;
		np->pmask = fullmask; 
//		np->pmask = 1 << i;		/* fake surface node pmask is itself for cmask ?? */
//		onp->pmask = inp->pmask |= np->pmask;
		s->nbp++;		/* Number of boundary planes */
	} else {
		s->n[-i-1]->ix = -i-1;		/* Label unused node */
	}

#ifdef DEBUG
	printf("Number of boundary planes = %d\nDiscovering all valid veronoi sub-surfaces\n",s->nbp);
#endif

	discover_subsuf(s);

#ifdef DEBUG
	printf("Creating rectangular initial vertexes\n");
#endif

	/* Compute initial node weighting */
	for (vtxwt = 0.0, i = 0; i < (sizeof(ivtx_whts)/sizeof(double)-1); i++)
		vtxwt += log(ivtx_whts[i]);
	vtxwt += ivtx_whts[i];

	/* Create initial verticies, one for each di combination of planes, */
	/* and keep the ones that are in gamut. */
	DC_INIT(co);
	while(!DC_DONE(co)) {
		double p[MXPD];
		vtx *vi, *vo;			/* Inside and outside vertex */

		/* Compute vertex location */
		for (e = 0; e < di; e++) {
			if (co[e] != 0)
				p[e] = s->imax[e];
			else
				p[e] = s->imin[e];
		}

		if (ofps_would_clip_point(s, p)) {
#ifdef DEBUG
			printf("Position %s rejected, out of gamut\n",ppos(di,p));
#endif
			goto next_co;
		}
#ifdef DEBUG
		printf("Position %s accepted\n",ppos(di,p));
#endif

		vi = new_vtx(s);
		vo = new_vtx(s);

		for (e = 0; e < di; e++)
			vi->p[e] = vtxwt * p[e];
		ofps_cc_percept(s, vi->v, vi->p);

		for (e = 0; e < di; e++)
			vo->p[e] = (10.0 * (p[e] - 0.5)) + 0.5;
//		ofps_cc_percept(s, vo->v, vo->p);

		/* Compute nodes involved */
		for (e = 0; e < di; e++) {
			if (co[e] == 0) {
				vo->nix[e] = vi->nix[e] = -1 - (2 * e + 0);
			} else {
				vo->nix[e] = vi->nix[e] = -1 - (2 * e + 1);
			}
		}

		vi->nix[di] = iix;		/* First nodee */
#ifdef DEBUG
		printf("ivertex nix %s\n",pcomb(di,vi->nix));
#endif
		sort_nix(s, vi->nix);
		vi->eperr = 10000.0;		/* Very bad, so they get chosen first */
		vi->eserr = 10000.0;

		det_vtx_gsurf(s, vi);		/* Set pmask & cmask */
		sm_cp(s, &vi->vm, &s->sc[vi->cmask].a_sm);	/* Set visibility */

		vi->ifake = 1;			/* Inside fake */

		vtx_cache_add(s, vi);	 /* Add it to the vertex cache and spatial accelleration grid */
		ofps_add_vacc(s, vi);
		ofps_add_vseed(s, vi);

		vo->nix[di] = oix;		/* Fake outside node */
#ifdef DEBUG
		printf("overtex nix %s\n",pcomb(di,vo->nix));
#endif
		sort_nix(s, vo->nix);
		vo->eperr = vtxwt * -9.0;		/* Better than zero error */
		vo->eserr = vtxwt * -9.0;
		/* Leave pmask,cmask = 0 */
		vo->pmask = vi->pmask;	/* Copy from inner vertexes */
		vo->cmask = vi->cmask;
		sm_cp(s, &vo->vm, &s->sc[vo->cmask].a_sm);	/* Set visibility */

		vo->ofake = 1;			/* Outside fake - don't plot vnets and don't use */
								/* for perceptual function extension. */
		vo->used = 1;			/* Not a candidate for seeding */
	
	next_co:;
		DC_INC(co);
	}

	/* Add ink limit  vertexes */
	if (doink) {	/* Ink limit plane is orthogonal to diagonal */
		COMBO(nco, MXPD, di-1, s->nbp-1);		/* di-1 out of neighbor nodes combination counter */

#ifdef DEBUG
		printf("Creating ink limit vertexes\n");
#endif
		/* Intersect the ink limit plane with each combination of */
		/* it and and di-1 of the existing planes, to generate */
		/* potential vertexes, and keep the ones that are in gamut. */
		CB_INIT(nco);
		while (!CB_DONE(nco)) {
			pleq *peqs[MXPD];
			double p[MXPD];

			for (e = 0; e < (di-1); e++) {
				peqs[e] = &s->gpeqs[nco[e]];
			}
			peqs[e] = &s->gpeqs[2 * di];

			/* Compute device location of intersection */
			if (comp_vtx(s, p, peqs) == 0) {
				vtx *vi, *vo;			/* Inside and outside vertex */

				if (ofps_would_clip_point(s, p)) {
#ifdef DEBUG
					printf("Position %s rejected, out of gamut\n",ppos(di,p));
#endif
					goto next_nco;
				}
#ifdef DEBUG
				printf("Position %s accepted\n",ppos(di,p));
#endif

				vi = new_vtx(s);
				vo = new_vtx(s);

				/* Device and perceptual */
				for (e = 0; e < di; e++)
					vi->p[e] = vtxwt * p[e];
				ofps_cc_percept(s, vi->v, vi->p);

				for (e = 0; e < di; e++)
					vo->p[e] = (10.0 * (p[e] - 0.5)) + 0.5;
//				ofps_cc_percept(s, vo->v, vo->p);

				for (e = 0; e < (di-1); e++)
					vo->nix[e] = vi->nix[e] = -1-nco[e];	/* Fake gamut surface plane nodes */
				vo->nix[e] = vi->nix[e] = -2 * di -1;		/* Fake ink limit node */

				vi->nix[di] = iix;		/* First node */
#ifdef DEBUG
				printf("ivertex nix %s\n",pcomb(di,vi->nix));
#endif
				sort_nix(s, vi->nix);
				vi->eperr = 10000.0;		/* Very bad */
				vi->eserr = 10000.0;

				det_vtx_gsurf(s, vi);	/* Set pmask & cmask */
				sm_cp(s, &vi->vm, &s->sc[vi->cmask].a_sm);	/* Set visibility */
		
				vi->ifake = 1;			/* Inside fake */

				vtx_cache_add(s, vi);	 /* Add to vertex cache and spatial accelleration grid */
				ofps_add_vacc(s, vi);
				ofps_add_vseed(s, vi);

				vo->nix[di] = oix;		/* Fake outside node */
#ifdef DEBUG
				printf("overtex nix %s\n",pcomb(di,vo->nix));
#endif
				sort_nix(s, vo->nix);
				vo->eperr = vtxwt * -9.0;		/* Better than zero error */
				vo->eserr = vtxwt * -9.0;
				/* Leave pmask,cmask = 0 */
				vo->pmask = vi->pmask;	/* Copy from inner vertexes */
				vo->cmask = vi->cmask;
				sm_cp(s, &vo->vm, &s->sc[vo->cmask].a_sm);	/* Set visibility */

				vo->ofake = 1;			/* Outside fake - don't plot vnets and don't use */
										/* for perceptual function extension. */
				vo->used = 1;			/* Not a candidate for seeding */
			}
		next_nco:;
			CB_INC(nco);
		}			/* Next combination */
	}

	/* Create an initial vertex network */
	create_vtx_net(s);

	/* Fit the unbounded perceptual model to just the inside vertexes */
	if (s->pmod_init == 0)
		init_pmod(s);

	/* Compute the nodes node and vertex lists */
	ofps_re_create_node_node_vtx_lists(s);

#ifdef DUMP_STRUCTURE
	printf("Done binit\n");
	dump_node_vtxs(s, 1);
//	dump_node_vtxs2(s, "Done binit");
	printf("=========================================================================\n");
#endif
#ifdef DUMP_PLOT_SEED
	dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 0, -1);		/* Device, No wait, verticies */
#endif /* DUMP_PLOT_SEED */
}

/* --------------------------------------------------- */
/* Setup the perceptual lookup cache */
static void
ofps_init_pcache(ofps *s) {
	int i, e;
	int di = s->di;
	int gr, gres[MXPD];
	int tinp = s->tinp;

#ifdef DEBUG
	printf("Initializing perceptual lookup cache\n");
#endif

	/* Choose a grid resolution that aims for aproximately TNPAGRID nodes per grid */
	if (tinp > 10000)
		tinp = 10000;
	gr = (int)(pow(tinp/TNPAGRID, 1.0/di) + 0.5);
	gr |= 1;			/* make it odd */

	if (gr < TNPAGRIDMINRES)
		gr = TNPAGRIDMINRES;
	if (gr > TNPAGRIDMAXRES)
		gr = TNPAGRIDMAXRES;

#ifndef DEBUG
	if (s->verb)
#endif
	{
		printf("Perceptual cache resolution = %d\n",gr);
		printf("Seeding cache..."); fflush(stdout);
	}

	/* Create a rspl to cache the perceptual lookup */
	if ((s->pcache = new_rspl(RSPL_NOFLAGS, s->di, s->di)) == NULL)
		error("new_rspl failed");

	for (e = 0; e < di; e++)
		gres[e] = gr;
	s->pcache_res = gr;

//	s->pcache->set_rspl(s->pcache, RSPL_SET_APXLS, s->od, s->percept, NULL, NULL, gres, NULL, NULL);

	/* Filtering seems to make this more robust for some profiles, less for others. */
	if (s->percept != default_ofps_to_percept)
		s->pcache->set_rspl(s->pcache, RSPL_NOFLAGS, s, filtered_ofps_to_percept,
		                                            NULL, NULL, gres, NULL, NULL);
	else
		s->pcache->set_rspl(s->pcache, RSPL_NOFLAGS, s->od, s->percept,
		                                  NULL, NULL, gres, NULL, NULL);

	/* Hmm. Should we store the underlying ->percept & ->od somewhere before we overwrite it ? */
	s->percept = ofps_cache_percept;
	s->od = s->pcache;

#ifndef DEBUG
	if (s->verb)
#endif
		printf("done\n");

}

/* --------------------------------------------------- */
/* Setup the acceleration grid structure and perceptual cache. */
/* The grid is in device space, although it is used to find the point */
/* with the smallest eperr. */
/* (Note that ofps_cc_percept() can't be called on clipped values yet) */
static void
ofps_init_acc1(ofps *s) {
	int i, e;
	int di = s->di;
	int gres[MXPD];
	int tinp = s->tinp;

#ifdef DEBUG
	printf("Initializing accelleration array (1)\n");
#endif

	/* Create acceleration grid array */

	/* Choose a grid resolution that aims for aproximately TNPAGRID nodes per grid */
	if (tinp > 10000)
		tinp = 10000;
	s->agres = (int)(pow(tinp/TNPAGRID, 1.0/di) + 0.5);
	if (s->agres < 1)
		s->agres = 1;

	if (s->verb)
		printf("Acceleration grid res = %d\n",s->agres);

#ifdef DEBUG
	printf("Acceleration grid res = %d\n",s->agres);
#endif

	/* Cell width in grid units */
	s->gw = 1.0/s->agres;

	/* Compute grid index multipliers */
	/* (We allocate an two extra rows for boundary cells to be looked up.) */
	for (s->gim[0] = 1, e = 1; e < di; s->gim[e] = s->gim[e-1] * (s->agres+2), e++)
		;
	
	/* Compute cell diagonal distance */
	s->gcd = sqrt((double)di * s->gw * s->gw);

	/* Compute number of cells in grid (with two extra rows) */
	for (s->nig = 1, e = 0; e < di; e++)
		s->nig *= (s->agres+2);

	/* Allocate grid (with two extra rows) */
	if ((s->_grid = (acell *)malloc(sizeof(acell) * s->nig)) == NULL)
		error ("ofps: malloc failed for acceleration grid");

	/* Set pointer to base of grid without extra row */
	for (s->grid = s->_grid, e = 0; e < di; e++)
		s->grid += s->gim[e];

	/* Initialise grid (including extra gruard rows) */
	{
		DCOUNT(co, MXPD, di, -1, -1, (s->agres+1));

		i = 0;
		DC_INIT(co);
		while (!DC_DONE(co)) {
			acell *cp = &s->_grid[i];
			unsigned int gflag = 0;

			for (e = 0; e < di; e++) {
				if (co[e] < 0 || co[e] >= s->agres)
					gflag = BOUND_GFLAG;
				cp->co[e] = co[e];			/* Grid coordinate of base of cell */
				cp->p[e] = co[e] * s->gw;	/* Device coord of base of cell */
				cp->cp[e] = (co[e] + 0.5) * s->gw;	/* Device coord of center of cell */
			}

			cp->gflag = gflag;
			cp->head = NULL;
			cp->vhead = NULL;

			DC_INC(co);
			i++;
		}
	}
	s->gflag = 0;

	/* Create the neighbour offset list */

	/* There are 3^di -1 neighbours  for each cell */
	for (s->nacnl = 1, e = 0; e < di; s->nacnl *= 3, e++)
		;
	s->nacnl--;

	if ((s->acnl = (int *)malloc(sizeof(int) * s->nacnl)) == NULL)
		error ("ofps: malloc failed on acnl list");

	/* Initialise list from cube */
	{
		DCOUNT(co, MXPD, di, -1, -1, 2);
		
		i = 0;
		DC_INIT(co);
		while (!DC_DONE(co)) {

			/* check we're not at the center cell */
			for (e = 0; e < di; e++) {
				if (co[e] != 0)
					break;
			}
			if (e < di) {		/* Not center cell */
				/* Compute offset */
				for (s->acnl[i] = 0, e = 0; e < di; e++) {
					s->acnl[i] += co[e] * s->gim[e]; 
				}
//printf("~1 acnl[%d] for co %s = %d\n",i,pco(di,co),s->acnl[i]);
				i++;
			}
			DC_INC(co);
		}
	}
}

/* Init the grid location p[] and v[] values */
static void
ofps_init_acc2(ofps *s) {
	int i, e, di = s->di;
	int k;
	DCOUNT(co, MXPD, di, 0, 0, 2);
	double maxratio, avgratio, noratio;
	double aitters = 0.0;

#ifdef DEBUG
	printf("Initializing accelleration array (2)\n");
#endif

	for (i = 0; i < s->nig; i++) {
		acell *cp = &s->_grid[i];

		/* Lookup perceptual base and center values */
		ofps_cc_percept(s, cp->v, cp->p);
		ofps_cc_percept(s, cp->cv, cp->cp);
	}

	/* Compute the worst case eperr from a corner to the center */ 
	maxratio = -1.0;
	avgratio = noratio = 0.0;
	for (i = 0; i < s->nig; i++) {
		acell *cp = &s->_grid[i];
		double ratio;
		double mov[MXPD];
	
		if (cp->gflag == BOUND_GFLAG)
			continue;
		
#define ACELITERS 20
		for (k = 0; ; k++) {
			double eperr_avg, eperr_min, eperr_max, no;
			cp->eperr = 0.0;
	
			DC_INIT(co);
			eperr_avg = no = 0.0;
			eperr_min = 1e300;
			eperr_max = -1.0;
			while (!DC_DONE(co)) {
				acell *np = &s->_grid[i];
				double eperr;
				int j;
	
				/* Locate cell corner */
				for (j = 0, e = 0; e < di; e++)
					j += co[e] * s->gim[e]; 
				np = &s->_grid[i + j];
	
				/* eperr from that corner to center of this cell */
				eperr = ofps_comp_eperr(s, NULL, cp->cv, cp->cp, np->v, np->p, 0);
				eperr_avg += eperr;
				if (eperr > eperr_max)
					eperr_max = eperr;
				if (eperr < eperr_min)
					eperr_min = eperr;

				no++;

				if (eperr > cp->eperr)
					cp->eperr = eperr;
	
				DC_INC(co);
			}
			eperr_avg /= no;

			ratio = eperr_max/eperr_min;

			if (k >= ACELITERS || ratio < 1.2) {
				avgratio += ratio;
				noratio++;
				if (ratio > maxratio)
					maxratio = ratio;

				break;

			} else {

				/* Adjust the center position to minimuze range of eperr's */
				for (e = 0; e < di; e++)
					mov[e] = 0.0;

//printf("~1 cp was at %s\n",ppos(di,cp->cp));
				DC_INIT(co);
				while (!DC_DONE(co)) {
					acell *np = &s->_grid[i];
					double eperr, wf;
					int j;
		
					/* Locate cell corner */
					for (j = 0, e = 0; e < di; e++)
						j += co[e] * s->gim[e]; 
					np = &s->_grid[i + j];
		
					/* Compose new center point from weighted corner points. */
					/* Weighting is proportional to eperr value */
					eperr = ofps_comp_eperr(s, NULL, cp->cv, cp->cp, np->v, np->p, 0);

					if (eperr < eperr_avg) {
						/* Move away from corner */
						wf = (eperr_avg - eperr)/eperr_avg;
//printf("~1 eperr %f, avg %f, min %f, wf %f\n",eperr,eperr_avg,eperr_min,wf);
						for (e = 0; e < di; e++)
							mov[e] += wf * (cp->cp[e] - np->p[e]); 
					} else {
						/* Move towards corner */
						wf = (eperr - eperr_avg)/eperr_avg;
//printf("~1 eperr %f, avg %f, max %f, wf %f\n",eperr,eperr_avg,eperr_max,wf);
						for (e = 0; e < di; e++)
							mov[e] += wf * (np->p[e] - cp->cp[e]); 
					}

					DC_INC(co);
				}
				for (e = 0; e < di; e++) {
					mov[e] = 1.2 * mov[e] / no;
					cp->cp[e] += mov[e];
				}
				ofps_cc_percept(s, cp->cv, cp->cp);
//printf("~1 moving by %s to %s\n",ppos(di,mov),ppos(di,cp->cp));
			}
		}
		aitters += k;

		cp->eperr *= CELLMAXEPERRFF;		/* Times the fudge factor */
	}
	aitters /= s->nig;

	avgratio /= noratio;
#ifdef DEBUG
	printf("Average acell eperr ratio = %f, maximum = %f, avg itters %f\n",avgratio,maxratio,aitters);

#endif

	s->agrid_init = 1;
}

/* Convert a location into an acceleration cell index */
static int
ofps_point2cell(ofps *s, double *v, double *p) {
	int i, e, di = s->di;
	int agres = s->agres;
	double pp[MXPD];

	ofps_clip_point(s, pp, p);

	for (i = e = 0; e < di; e++) {
		int t;
		t = (int)floor(agres * pp[e]);
		if (t < 0)
			t = 0;
		else if (t >= agres)
			t = (agres-1);
		i += s->gim[e] * t;
	}
	return i;
}

/* Diagnostic: Return the grid coordinates */
static void ofps_gridcoords(ofps *s, int *c, double *v, double *p) {
	int i, e, di = s->di;
	int agres = s->agres;
	double pp[MXPD];

	ofps_clip_point(s, pp, p);

	for (i = e = 0; e < di; e++) {
		int t;
		c[e] = (int)floor(agres * pp[e]);
		if (c[e] < 0)
			c[e] = 0;
		else if (c[e] >= agres)
			c[e] = (agres-1);
	}
}

/* Add a node to the spatial acceleration grid */
/* Note that little more than the node perceptual value */
/* may be valid when this is called. */
static void
ofps_add_nacc(ofps *s, node *n) {
	int pci;
	acell *cp;
	
	if (n->ix < 0)
		return;

	pci = ofps_point2cell(s, n->v, n->p);
	cp = &s->grid[pci];
	n->n = cp->head;
	if (cp->head != NULL)
		cp->head->pn = &n->n;
	cp->head = n;
	n->pn = &cp->head;
	n->pci = pci; 
	n->cell = cp;

#ifdef SANITY_CHECK_CLOSEST
	if (s->agrid_init) {
		double eperr;
		/* Check that the eperr to the center of the cell */
		/* is less than the worst case for that cell */
		eperr = ofps_comp_eperr(s, NULL, cp->cv, cp->cp, n->v, n->p, n->nsp);
		if (eperr > cp->eperr) {
			warning("Sanity check ofps_add_nacc() node ix %d eperr %f > cell eperr %f",n->ix,eperr,cp->eperr);
			printf("Sanity check ofps_add_nacc() node ix %d eperr %f > cell eperr %f\n",n->ix,eperr,cp->eperr);
#ifdef SANITY_CHECK_CLOSEST_FATAL
			error("ofps_add_nacc cell eperr failed");
#endif
		}
	}
#endif	/* SANITY_CHECK_CLOSEST */
}

/* Remove a node from the spatial acceleration grid */
static void
ofps_rem_nacc(ofps *s, node *n) {
	if (n->ix < 0)
		return;
	if (n->pn != NULL) {		/* If is on acceleration list, remove it */
		*n->pn = n->n;
		if (n->n != NULL)
			n->n->pn = n->pn;
	}
	n->pn = NULL;
	n->n = NULL;
}

/* Add a vertex to the spatial acceleration grid */
static void
ofps_add_vacc(ofps *s, vtx *vx) {
	int pci;
	acell *cp;

	/* Normal spatial acceleration grid */
	pci = ofps_point2cell(s, vx->v, vx->p);
	cp = &s->grid[pci];
	vx->n = cp->vhead;
	if (cp->vhead != NULL)
		cp->vhead->pn = &vx->n;
	cp->vhead = vx;
	vx->pn = &cp->vhead;
	vx->pci = pci;

#ifdef DEBUG
	printf("Adding vertex no %d to spatial accelleration grid in cell %d\n",vx->no,pci);
#endif

#ifdef SANITY_CHECK_CLOSEST
	if (s->agrid_init) {
		int e, di = s->di;
		double eperr;
		double p[MXPD], v[MXPD];

		/* Check that the eperr to the center of the cell */
		/* is less than the worst case for that cell */

		/* Clip point in case it lies outside the grid, */
		/* and would give an excessive eperr */
		for (e = 0; e < di; e++) {
			p[e] = vx->p[e];
			if (p[e] < 0.0)
				p[e] = 0.0;
			else if (p[e] > 1.0)
				p[e] = 1.0;
		}
		ofps_cc_percept(s, v, p);
		eperr = ofps_comp_eperr(s, NULL, cp->cv, cp->cp, v, p, 0);

		if (eperr > cp->eperr) {

//printf("~1 Cell ix %d co %s center %s (%s), vtx at %s (%s) clipped to %s (%s)\n",pci,pco(s->di,cp->co),ppos(di,cp->cp),ppos(di,cp->cv),ppos(di,vx->p),ppos(di,vx->v),ppos(di,p),ppos(di,v));
			warning("Sanity check ofps_add_vacc() vtx no %d eperr %f > cell eperr %f",vx->no,eperr,cp->eperr);
			printf("Sanity check ofps_add_vacc() vtx no %d eperr %f > cell eperr %f\n",vx->no,eperr,cp->eperr);

#ifdef SANITY_CHECK_CLOSEST_FATAL
			error("ofps_add_vacc cell eperr failed");
#endif
		}
	}
#endif	/* SANITY_CHECK_CLOSEST */
}

/* Add a vertex to the seeding groups */
static void
ofps_add_vseed(ofps *s, vtx *vx) {
	double oog;
	int pci;
	acell *cp;

#ifdef DEBUG
	printf("Adding vertex no %d to sorted binary tree\n",vx->no);
#endif

	/* Add the vertex to the sorted binary trees */
	if ((aat_ainsert(s->vtreep, (void *)vx)) == 0)
		error("aat_ainsert vertex malloc failed");

	/* Out of gamut vertexes are not candidates for seeds */
	if ((oog = ofps_oog(s, vx->p)) > COINTOL) {
		vx->used = 1;
//printf("Setting used on vtx no %d, used %d, eserr %f, vm %s nsp %d oog by %e\n",vx->no,vx->used,vx->eserr,psm(s,&vx->vm),vx->nsp, oog);
	}

	if (vx->used == 0) {
#ifdef INDEP_SURFACE
		/* Only pick full dimensional visible vertexes for seeding group, */
		/* since only they have a full-d error value. */
		if (sm_andtest(s, &s->sc[0].a_sm, &vx->vm) != 0) {
#endif
//printf("Adding (3) vtx no %d, used %d, eserr %f, vm %s nsp %d\n",vx->no,vx->used,vx->eserr,psm(s,&vx->vm),vx->nsp);
			if ((aat_ainsert(s->vtrees[vx->nsp], (void *)vx)) == 0)
				error("aat_ainsert vertex malloc failed");
#ifdef INDEP_SURFACE
		} else {
//printf("Not adding (2) vtx no %d, used %d, eserr %f, vm %s nsp %d\n",vx->no,vx->used,vx->eserr,psm(s,&vx->vm),vx->nsp);
		}
#endif
	}
//else printf("Not adding (3) vtx no %d, used %d, eserr %f, vm %s nsp %d\n",vx->no,vx->used,vx->eserr,psm(s,&vx->vm),vx->nsp);
}

/* Remove a vertex from the seeding groups */
static void
ofps_rem_vseed(ofps *s, vtx *vx) {

#ifdef DEBUG
	printf("Removing vertex no %d  from sorted binary tree\n",vx->no);
#endif

	/* Remove the vertex from the sorted binary tree */
	if ((aat_aerase(s->vtreep, (void *)vx)) == 0)
		error("aat_aerase vertex failed to find vertex no %d (3)", vx->no);

	if (vx->used == 0) {
#ifdef INDEP_SURFACE
	/* Only pick full dimensional visible vertexes for seeding group, */
	/* since only they have a full-d error value. */
    if (sm_andtest(s, &s->sc[0].a_sm, &vx->vm) != 0) {
#endif
//printf("Removing (4) vtx no %d,  0x%x, used %d, eserr %f, vm %s nsp %d\n",vx->no,vx, vx->used,vx->eserr,psm(s,&vx->vm),vx->nsp);
		if ((aat_aerase(s->vtrees[vx->nsp], (void *)vx)) == 0)
			error("aat_aerase vertex failed to find vertex no %d (4)", vx->no);
#ifdef INDEP_SURFACE
		} else {
//printf("Not removing (1) vtx no %d, 0x%x, used %d, eserr %f, vm %s nsp %d\n",vx->no,vx,vx->used,vx->eserr,psm(s,&vx->vm),vx->nsp);
		}
#endif
	}
// else printf("Not removing (2) vtx no %d, 0x%x, used %d, eserr %f, vm %s nsp %d\n",vx->no,vx,vx->used,vx->eserr,psm(s,&vx->vm),vx->nsp);
}

/* Remove a vertex from the spatial acceleration grid */
static void
ofps_rem_vacc(ofps *s, vtx *vx) {

	/* Remove from spatial acceleration grid */
	if (vx->pn != NULL) {		/* If is on acceleration list, remove it */
		*vx->pn = vx->n;
		if (vx->n != NULL)
			vx->n->pn = vx->pn;
	}
	vx->pn = NULL;
	vx->n = NULL;
}

/* Clear the spatial acceleration grid */
static void
ofps_reset_acc(ofps *s) {
	int i;

	for (i = 0; i < s->nig; i++) {
		acell *cp = &s->_grid[i];
		if (cp->gflag != BOUND_GFLAG)
			cp->gflag = 0;
		cp->head = NULL;
		cp->vhead = NULL;
	}
	s->gflag = 0;
}

/* --------------------------------------------------- */

/* Creat a randomized order list of pointers to the fixed points */
static void
ofps_setup_fixed(
ofps *s,
fxpos *fxlist,			/* List of existing fixed points */
int fxno				/* Number in fixed list */
) {
	int e, di = s->di;
	int i, j;

	s->fnp = 0;
	if (fxno == 0)
		return;

	/* Allocate a list of pointers sufficient for all the fixed points */
	if ((s->ufx = (fxpos **)calloc(sizeof(fxpos *), fxno)) == NULL)
		error ("ofps: malloc failed on pointers to fixed list");

	/* Add each fixed point to the list */
	for (i = 0; i < fxno; i++) {

		/* Clip the fixed point */
		ofps_clip_point7(s, fxlist[i].p, fxlist[i].p);

		/* Comute perceptual attributes */
		s->percept(s->od, fxlist[i].v, fxlist[i].p);

		/* Skip any duplicate points, or Voronoi will get confused.. */
		for (j = 0; j < s->fnp; j++) {
			for (e = 0; e < di; e++) {
				if (fabs(s->ufx[j]->p[e] - fxlist[i].p[e]) > 1e-5)
					break;		/* Not a match */
			}
			if (e >= di)
				break;			/* Is a match */
		}
		if (j < s->fnp)
			continue;			/* Skip adding this point */

		s->ufx[s->fnp++] = &fxlist[i];
	}

	/* Randomly shuffle the fixed points */
	for (i = 0; i < s->fnp; i++) {
		fxpos *tp;

		j = i_rand(0, s->fnp-1);

		/* Swap the pointers */
		tp = s->ufx[i];
		s->ufx[i] = s->ufx[j];
		s->ufx[j] = tp;
	}
	s->tinp -= (fxno - s->fnp);

}

/* Seed the object with any fixed points */
/* (I think this is only used if ofps is used to check the stats */
/*  on all the points. ) */
/* Return NZ on failure */
static int
ofps_add_fixed(
ofps *s
) {
	int e, di = s->di;
	int i, j, ii;

	/* Add fixed points if there are any */
	if (s->fnp == 0)
		return 0;

	if (s->verb)
		printf("Adding %d unique fixed points\n",s->fnp);

	for (i = 0; i < s->fnp; i++) {
		node *p = s->n[i];	/* Destination for point */

		/* Make sure that fixed point is within our gamut */
		ofps_clip_point(s, s->ufx[i]->p, s->ufx[i]->p);

		for (e = 0; e < di; e++) {	/* copy device and perceptual coords */
			p->op[e] = p->p[e] = s->ufx[i]->p[e];
			p->v[e] = s->ufx[i]->v[e];
		}

		/* Count gamut surface planes it lies on */
		det_node_gsurf(s, p, p->p);

		p->fx = 1;							/* is a fixed point */

		/* Compute the Voronoi for it, and inc s->np */
		if (add_node2voronoi(s, s->np, 1)) {
			/* In theory we could try adding points in a different order, */
			/* by resetting the voronoi, shuffling all the fixedpoints */
			/* and re-adding them again. */
			warning("Adding a fixed point failed to hit any vertexes, and no points to swap with!");	
			return 1;
		}

		if (s->verb)
			printf("%cAdded fixed %d/%d",cr_char,i,s->fnp); fflush(stdout);

#ifdef DUMP_STRUCTURE
		printf("Done node %d\n",s->np);
		dump_node_vtxs(s, 0);
		printf("=========================================================================\n");
#endif
#ifdef DUMP_PLOT_SEED
		dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 0, -1);		/* Device, No wait, verticies */
#endif /* DUMP_PLOT_SEED */
	}

	return 0;
}

/* Seed the object with any fixed and movable incremental farthest points. */
/* (We are only called if there is at leaset one movable point) */
static void
ofps_seed(ofps *s) {
	int e, di = s->di;
	int ii, i, j, k, fc;
	double rerr;
	int needfirst = 0;	/* Need a special first seed point (seems better without this ?) */
	int dofixed = 0;	/* Do a fixed point next */
	int abortonfail = 0;	/* Abort on failing to add fixed points (isn't always good ?) */
	int nsp = 0;		/* Number of surface points */
	aat_atrav_t *aat_tr;

	if (s->verb)
		printf("\n");

	if ((aat_tr = aat_atnew()) == NULL)
		error("aat_atnew returned NULL");

	if (s->verb) {
		printf("There are %d unique fixed points to add (%d total fixed points)\n",s->fnp, s->fxno);
		printf("There are %d far spread points to add\n",s->tinp - s->fnp);
	}

	if (!needfirst && s->fnp > 1) {	/* There are fixed points to add */
		dofixed = s->fnp > 2 ? 2 : 1;
	}

	/* Seed all the points. */
	/* (i is the node we're creating, j is the verbose interval count, */
	/*  fc is the count of the fixed points added, ii is the movable count) */
	for (fc = j = i = ii = 0; i < s->tinp; i++, j++) {
		node *p = s->n[i];		/* New node */
		double spref_mult;

		/* Compute current surface preference weighting */
		if (s->tinp - s->fnp <= 1)	/* Prevent divide by zero */
			spref_mult = 0.0;
		else
			spref_mult = ii/(s->tinp - s->fnp - 1.0);
		spref_mult = (1.0 - spref_mult) * s->ssurfpref + spref_mult * s->esurfpref;

		if (needfirst) {			/* No initial fixed points, so seed the first */
									/* point as a special. */
			double min[MXPD], max[MXPD];

			p->fx = 0;

			/* If there are no fixed points in the bulk, make the */
			/* first point such a point, to avoid pathology. */
			for (e = 0; e < di; e++)
				p->p[e] = s->imin[e] + (s->imax[e] - s->imin[e]) * 1.0/4.141592654;

			/* Clip the new location */
			ofps_clip_point8(s, p->p, p->p);

			s->percept(s->od, p->v, p->p);
#ifdef DEBUG
			printf("Creating first seed point (moveable %d out of %d)\n",i+1,s->tinp);
#endif
		} else if (dofixed) {	/* Setup to add a fixed point */

#ifdef DEBUG
			printf("Adding fixed point %d out of %d\n",fc+1,s->fnp);
#endif

			/* (Fixed points are already clipped and have perceptual value) */
			for (e = 0; e < di; e++) {	/* copy device and perceptual coords */
				p->p[e] = s->ufx[fc]->p[e];
				p->v[e] = s->ufx[fc]->v[e];
			}

			p->fx = 1;							/* is a fixed point */

			/* Count gamut surface planes it lies on */
			if (det_node_gsurf(s, p, p->p) != 0)
				nsp++;

		} else {	/* Setup to add a movable point */
			int k;
			int sf = 0;

			double spref_mult;
			double mx;
			vtx *vx, *bvx;
			double spweight[MXPD+1];		/* Surface preference weight table */
			double bspweight;				/* Biggest weight */

#ifdef DEBUG
			printf("Adding movable point %d out of %d\n",i+1,s->tinp);
#endif

			p->fx = 0;

			/* Compute current surface preference weighting */
			if (s->tinp - s->fnp <= 1)	/* Prevent divide by zero */
				spref_mult = 0.0;
			else
				spref_mult = ii/(s->tinp - s->fnp - 1.0);
			spref_mult = (1.0 - spref_mult) * s->ssurfpref + spref_mult * s->esurfpref;
	
			/* Until we use the next vertex, keep looking for a movable point */
			for (;;) {

#ifdef NEVER	/* DEBUG: Show the contents of each list */
				printf("All vertex list:\n");
				for (vx = s->uvtx; vx != NULL; vx = vx->link) {
					printf("   Vtx %d, used %d, eserr %f, weserr %f, vm %s\n",vx->no,vx->used,vx->eserr,vx->eserr * spweight[vx->nsp],psm(s,&vx->vm));
				}

				for (e = 0; e <= (di+1); e++) {
					printf("Sorted tree vertex list for nsp %d :\n",e);

					for (vx = aat_atlast(aat_tr, s->vtrees[e]); vx != NULL; vx = aat_atprev(aat_tr)) {
						printf("   Vtx %d, used %d, eserr %f, weserr %f, vm %s\n",vx->no,vx->used,vx->eserr,vx->eserr * spweight[vx->nsp],psm(s,&vx->vm));
					}
				}
#endif
				/* Compute the surface weighting multiplier for */
				/* each possible nsp + 1 */
				spweight[0] = 1.0;

				for (e = 1; e <= di; e++)
					spweight[e] = spweight[e-1] * spref_mult;

				/* Locate the Voronoi vertex with the greatest distance to a sampling points */
				for (mx = -1.0, bvx = NULL, e = 0; e <= (di+1); e++) {
					double weserr;

					/* Get largest eserr vertex for this nsp */
					if ((vx = aat_atlast(aat_tr, s->vtrees[e])) == NULL)
						continue;

					weserr = vx->eserr * spweight[vx->nsp];

//printf("~1 considering vertex no %d, eserr %f, weserr %f\n",vx->no,vx->eserr,weserr);
					if (weserr > mx) {
						mx = weserr;
						bvx = vx;
					}
				}
//if (bvx != NULL) printf("~1 got vertex no %d, eserr %f, weserr %f\n",bvx->no,bvx->eserr,mx);

#ifdef SANITY_CHECK_SEED
				/* Do exaustive search of candidate vertexes */
				for (vx = s->uvtx; vx != NULL; vx = vx->link) {
					double tweserr;

					tweserr = vx->eserr * spweight[vx->nsp];

					if (vx->used == 0 && tweserr > (mx + 10.0 * NUMTOL)
#ifdef INDEP_SURFACE
					/* Only pick full dimensional visible vertexes, */
					/* since only they have a full-d error value. */
					    && sm_andtest(s, &s->sc[0].a_sm, &vx->vm) != 0
#endif
					 ) {
						warning("!!!!!! Sanity: Didn't pick largest eperr vtx no %d %f, picked %d %f instead !!!!!!",vx->no,tweserr,bvx->no,mx);
						printf("!!!!!! Sanity: Didn't pick largest eperr vtx no %d %f, picked %d %f instead !!!!!!\n",vx->no,tweserr,bvx->no,mx);
						mx = tweserr;
						bvx = vx;
					}
				}
#endif /* SANITY_CHECK_SEED */

				if (bvx == NULL) {		/* We've failed to find a movable point */
										/* This could be because there are fixed points */
										/* at all candidate locations. */

					if (fc < s->fnp) {	/* Use a fixed point */
						/* (Fixed points are already clipped and have perceptual value) */
						for (e = 0; e < di; e++) {	/* copy device and perceptual coords */
							p->p[e] = s->ufx[fc]->p[e];
							p->v[e] = s->ufx[fc]->v[e];
						}
			
						p->fx = 1;							/* is a fixed point */
		
						/* Count gamut surface planes it lies on */
						if (det_node_gsurf(s, p, p->p) != 0)
							nsp++;
						break;		/* Go and use this point */
					}
					error("ofps: assert, there are no vertexes to choose in initial seed\n");
				}

#ifdef DEBUG
				printf("Picking vertex no %d at %s with weserr %f, mask 0x%x\n",bvx->no,ppos(di,bvx->p),mx,bvx->pmask);
#endif

				/* Don't pick the vertex again */
				bvx->used = 1;
//printf("Removing (4) vtx no %d, used %d, eserr %f, vm %s nsp %d\n",bvx->no,bvx->used,bvx->eserr,psm(s,&bvx->vm),bvx->nsp);
				if ((aat_aerase(s->vtrees[bvx->nsp], (void *)bvx)) == 0)
					error("aat_aerase vertex failed to find vertex no %d (5)", bvx->no);

				/* Add the new node */
				for (e = 0; e < di; e++)
					p->p[e] = bvx->p[e];

				/* Count gamut surface planes it lies on */
				if (det_node_gsurf(s, p, p->p) != 0)
					nsp++;

#ifdef RANDOM_PERTERB 
				/* Compute radius of closest real node to vertex */
				rerr = 1.0;
				for (k = 0; k <= di; k++) {
					double rads;
					int ix = bvx->nix[k]; 
					node *np;
					if (ix < 0)
						break;				/* Done */
					np = s->n[ix];
					for (rads = 0.0, e = 0; e < di; e++) {
						double tt = bvx->p[e] - np->p[e];
						rads += tt * tt;
					}
					if (rads < rerr)
						rerr = rads;
				}
				rerr = s->lperterb * sqrt(rerr);
				if (rerr < 0.001)
					rerr = 0.001;

				/* Add a random offset to the position, and retry */
				/* if it collides with an existing node or future fixed point */
				for (k = 0; k < 20; k++) {
					int pci;		/* Point list index */
					acell *cp;				/* Acceleration cell */
					node *p1;
				
					for (e = 0; e < di; e++)
						p->p[e] = bvx->p[e] + d_rand(-rerr, rerr);

					/* Confine node to planes vertex was on, */
					/* bit not if we're having trouble avoiding collissions */
					if (bvx->nsp > 0)
						confineto_gsurf(s, p->p, p->sp, p->nsp);

					ofps_clip_point9(s, p->p, p->p);

					s->percept(s->od, p->v, p->p);

					pci = ofps_point2cell(s, p->v, p->p);	/* Grid index of cell of interest */
					if ((cp = &s->grid[pci]) == NULL)
						break;						/* Nothing in cell */
					for (p1 = cp->head; p1 != NULL; p1 = p1->n) {
						for (e = 0; e < di; e++) {
							if (fabs(p->v[e] - p1->v[e]) > COINTOL)
								break;		/* Not cooincident */
						}
						if (e >= di) { 	/* Cooincident */
#ifdef DEBUG
							printf("Random offset node ix %d at %s collides with ix %d at %s - retry random %d\n",p->ix,ppos(di,p->p),p1->ix,ppos(di,p1->p),k);
#endif
							break;		/* Retry */
						}
					}
					if (p1 != NULL)		/* Coincident */
						continue;

					/* Check movable point against fixed points that are yet to be added */
					/* (~~~ Ideally we should use an accelleration structure to check */
					/*      for cooincidence rather than doing an exaustive search. ~~~) */
					if (fc < s->fnp) {	/* There are more fixed points to add */
						int f;

						for (f = fc; f < s->fnp; f++) {
							for (e = 0; e < di; e++) {
								if (fabs(p->p[e] - s->ufx[f]->p[e]) > COINTOL)
									break;		/* Not cooincident */
							}
							if (e >= di) { 	/* Cooincident */
#ifdef DEBUG
								printf("Movable node ix %d at %s collides with fixed point %d at %s - retry movable %d\n",p->ix,ppos(di,p->p),f,ppos(di,s->ufx[f]->p),k);
#endif
								break;
							}
						}
						if (f >= s->fnp)
							break;		/* movable point is not cooincident */
					} else {
						break;			/* Not cooincident, so OK */
					}
				}
				if (k >= 20) {
					/* This can happen if we didn't pick the absolute largest weserr, */
					/* and the vertex we ended up with is being confined to the same */
					/* location as an existing vertex by the planes is on. */
					/* (Why does it have an weperr > 0.0 then ????) */
					/* Give up on this point and chose another one. */
					continue;

//					error("ofps_seed: Assert, was unable to joggle cooincindent point");
				}
#else /* !RANDOM_PERTERB */
				/* Confine node to planes vertex was on */
				if (bvx->nsp > 0)
					confineto_gsurf(s, p->p, p->sp, p->nsp);

				ofps_clip_point9(s, p->p, p->p);

				s->percept(s->od, p->v, p->p);
#endif /* !RANDOM_PERTERB */

				/* Added this movable point, so chosen the next point */
				break;
			}	/* keep looking for a movable point */
		}

		/* We now have a first/fixed/moevable point to add */

/* hack test */
//p->p[0] = d_rand(0.0, 1.0);
//p->p[1] = d_rand(0.0, 1.0);
//ofps_cc_percept(s, p->v, p->p);

		/* Establish original position */
		for (e = 0; e < di; e++)
			p->op[e] = p->p[e];

		/* Compute the Voronoi for it, and inc s->np */
		/* Fail if we get a position fail */
		if (add_node2voronoi(s, i, dofixed && abortonfail)) {
			if (dofixed) {
				/* Pospone adding this vertex */
				if ((s->fnp - fc) >= (s->tinp - i - 1))	{	/* No room for moveable points */
//					error("Adding fixed point failed to hit any vertexes or posn. failed");	
					abortonfail = 0;
				} else
					dofixed = 0;
			}
			if (needfirst) {
				/* Hmm. The first seed point has failed. What should we do ? */
				error("Adding first seed point failed to hit any vertexes or posn. failed");	
			}
	
			/* Skip this point */
			--i;
			--j;
			continue;
		}

		/* Suceeded in adding the point */
		if (p->fx) {	/* Fixed point was added */
			fc++;
			if (dofixed > 0)		/* May not have been triggered by dofixed */
				dofixed--;
			if ((s->fnp - fc) >= (s->tinp - i - 1))	{	/* No room for moveable points */
				dofixed = s->fnp - fc;					/* Do all the fixed */
			}
		} else {		/* Movable point */
			ii++;
			if (fc < s->fnp) {	/* There are more fixed points to add */
				dofixed = s->fnp - fc;
				/* Add fixed 2 at a time to try and minimize the disruption */
				/* of the movable point edge priority */
				if (dofixed > 2)
					dofixed = 2;
			}
		}

		if (s->verb && (j == 11 || i == (s->tinp-1))) {
			printf("%cAdded %d/%d",cr_char,s->np,s->tinp); fflush(stdout);
			j = 0;
		}
#ifdef DUMP_STRUCTURE
		printf("Done node %d\n",i);
		dump_node_vtxs(s, 0);
		printf("=========================================================================\n");
#endif
#ifdef DUMP_PLOT_SEED
		dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 0, -1);		/* Device, No wait, verticies */
#endif /* DUMP_PLOT_SEED */

		needfirst = 0;		/* Must have done first */
	}
//printf("Number of gamut surface points = %d\n",nsp);

	aat_atdelete(aat_tr);

	if (s->verb)
		printf("\n");
}

/* Recreate the Voronoi diagram with the current point positions */
static void
ofps_redo_voronoi(
ofps *s
) {
	vtx *vx, *nvx;
	int i, j, k, e, di = s->di;

	/* Retry if we get a failure to add a point */
	for (k = 0; k < NINSERTTRIES; k++) {

		/* (~9 should think about smoothing the pre-conditioning lookup */
		/*   if the number of tries is high. Add this to rspl.) */
 
		/* Clear the voronoi nodes */
		node_clear(s, s->n[-2 * di - 2]);
		for (j = -s->gnp; j < s->np; j++)
			node_clear(s, s->n[j]);

		/* Delete the voronoi verticies */
		for (vx = s->uvtx; vx != NULL; vx = nvx) {
			nvx = vx->link;
			del_vtx(s, vx);
		}
		s->uvtx = NULL;

		/* Clear out the spatial acceleration grid */
		ofps_reset_acc(s);

		if (s->nv != 0)
			warning("ofps: Assert, clear didn't leave us with 0 vertexes");

		if (s->umid != NULL)
			warning("ofps: Assert, clear didn't empty used midpoint list");

		if (aat_asize(s->vtreep) != 0)
			warning("ofps: Assert, clear didn't empty vertex tree");
		for (e = 0; e <= (di+1); e++) {
			if (aat_asize(s->vtrees[e]) != 0)
				warning("ofps: Assert, clear didn't empty vertex tree");
		}

		/* Set number of points in voronoi to zero */
		s->np = 0;

		/* Initialse the empty veronoi etc. */ 
		ofps_binit(s);

		s->posfailstp = 0;

		/* Add all points in again. */
		for (i = 0 ;i < s->tinp; i++) {	/* Same order as before */

			/* Compute the Voronoi for it (will add it to spatial accelleration grid) */
			/* and increment s->np */
			if (add_node2voronoi(s, i, 0)) {

				/* Hmm. Shuffle and retry the whole thing. */
				shuffle_node_order(s);
				break;
			}

			/* If it's not going well, re-shuffle and abort too */
			if (i > 10 && s->posfailstp/(1.0+i) > 0.2) {
//printf("~1 after node %d, posfailes = %d, prop %f\n",i,s->posfailstp, s->posfailstp/(1.0+i));
				/* Hmm. Shuffle and and retry the whole thing. */
				if (s->verb > 1)
					warning("Too many nodes are failing to be inserted - reshuffling and re-starting\n");
				shuffle_node_order(s);
				break;
			}

#ifdef DUMP_STRUCTURE
			printf("Done node %d\n",i);
			dump_node_vtxs(s, 0);
//		ofps_re_create_node_node_vtx_lists(s);
//	if ((s->optit+1) >= 4)
//		{ char buf[200]; sprintf(buf, "Itteration %d node ix %d",s->optit+1,s->np-1); dump_node_vtxs2(s, buf); }
			printf("=========================================================================\n");
#endif
#ifdef DUMP_PLOT_RESEED
			dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 0, -1);		/* Device, No wait, verticies */
#endif /* DUMP_PLOT_RESEED */
		}
		if (i >= s->tinp) {
#ifdef DEBUG
			if (k > 1) printf("Took %d retries\n",k-1);
#endif /* DEBUG */
			break;
		}
		/* Retry the whole thing */
	}
	if (k >= NINSERTTRIES)
		error("Failed to re-seed the veronoi after %d tries - too many node insertion failures ?",NINSERTTRIES);
}

/* ----------------------------------------------------------- */
/* Ideas for improving the accelleration:

	When there is no SUBD, then it is possible that the node
	neighbourhood net (if it is kept up to date during seeding)
	could be used to locate the closest node and then vertex.
	(It can't be used for fixup, because the voronoi properties
	aren't true during fixup.)
	Starting at at the first node found using the spiral structure,
	check all it's neigbours and if it's neighbour is closer to
	the target, switch to it. If no neighbour is closer,
	then that is the closest node. 
	The closest vertex is then connected to the closest node ?
	
*/

#undef DEBUG_FCLOSE

/* Given a node, locate all vertexes that it hits. */
/* s->flag is assumed to be relevant for the given node. */
/* Any hit vertexes are added to the s->nxh list. */
/* s->vvchecks and s->nvcheckhits will be updated. */
/* Return nz if vertexs were hit */
/* (This only returns visible vertexes.) */
static int ofps_findhit_vtxs(ofps *s, node *nn) {
	int e, di = s->di;
	int i, j;
	int pci;				/* Point cell index */
	acell *cp;
	vtx *vx;
	double beperr, eperr;
	acell *slist = NULL, *sliste = NULL;		/* Next to search list */
	int hit = 0;

	if (nn->ix < 0)
		error("ofps_findhit_vtxs given gamut boudary node ix %d",nn->ix);

#ifdef DEBUG
	if (s->agrid_init == 0)
		error("ofps_findhit_vtxs() called before agrid_init");
#endif

#ifdef DEBUG_FCLOSE
	printf("\nLocating a hit vtx to node at p = %s, v = %s\n", ppos(di,nn->p), ppos(di,nn->v));
#endif

	/* Determine the largest eperr of any vertex */
	{
		aat_atrav_t *aat_tr;

		beperr = 1e300;

		if ((aat_tr = aat_atnew()) == NULL)
			error("aat_atnew returned NULL");

		/* Find the largest vertex eperr visible to the node */
		for (vx = aat_atlast(aat_tr, s->vtreep); vx != NULL; vx = aat_atprev(aat_tr)) {
#ifdef INDEP_SURFACE
			if (sm_vtx_node(s, vx, nn) == 0)
				continue;
#endif	/* INDEP_SURFACE */
			beperr = vx->eperr;
			break;
		}
		aat_atdelete(aat_tr);

#ifdef DEBUG_FCLOSE
		printf("Largest eperr of any vertex = %f\n", beperr);
//		fprintf(stderr,"Largest eperr of any vertex = %f\n", beperr);
#endif
	}

	s->nvfschd += s->nv;		/* Number of vertexes in a full search */
	s->naccsrch++;				/* Number of searches */

	/* Do a breadth first seed search for any hit vertexes, or until */
	/* we run out of cells that could possibly be hits. */

	/* Locate a starting cell using the grid */
	pci = ofps_point2cell(s, nn->v, nn->p);	/* Grid index of cell of interest */
	cp = &s->grid[pci];

	s->gflag++;					/* cell touched flag */ 

	/* Put the starting cell on the search list */
#ifdef DEBUG_FCLOSE
	printf("Adding cell ix %d co %s to slist\n",cp - s->grid, pco(di,cp->co));
#endif
#ifdef NEVER
	if (sliste == NULL) {		/* First in empty list */
		slist = cp;
	} else {
		sliste->slist = cp; 	/* Add to end of list */
	}
	sliste = cp;
	cp->slist = NULL;
#else
	cp->slist = slist;
	slist = cp;
#endif
	cp->gflag = s->gflag;	/* Cell is on list to be searched */

	/* until we run out of cells to search */
	for (;slist != NULL;) {
		acell *ncp;

		/* For each cell in the search list, check it and recursion. */
		for (cp = slist, slist = sliste = NULL; cp != NULL; cp = ncp) {
			double ceperr;
			ncp = cp->slist;

#ifdef DEBUG_FCLOSE
			printf("Checking cell ix %d co %s\n",cp - s->grid,pco(di,cp->co));
#endif

			/* Compute the smallest eperr possible in this cell, by computing the */
			/* eperr of the cell center to the node minus the estimated */ 
			/* largest eperr of any point within the cell to the center. */
			ceperr = ofps_comp_eperr(s, NULL, cp->v, cp->p, nn->v, nn->p, nn->nsp);
			eperr = ceperr - cp->eperr;
			
//printf("~1 ceperr %f, cp->eperr %f, eperr %f, beperr %f\n",ceperr,cp->eperr,eperr,beperr);
			/* If smallest possible eperr is larger than largest vertexe eperr */
			if (eperr > beperr) {
//printf("~1 skipping cell\n");
#ifdef SANITY_CHECK_CLOSEST
				/* Check all nodees in the cell anyway */
				for (vx = cp->vhead; vx != NULL; vx = vx->n) {
					int par = 0;

					if (vx->cflag == s->flag)
						continue;
#ifdef INDEP_SURFACE
					if (sm_vtx_node(s, vx, nn) == 0)
						continue;
#endif	/* INDEP_SURFACE */
		
					if (nn->ixm & vx->nix[MXPD+2]) {	/* Is in nixm */
						for (e = 0; e <= di; e++) {		/* Do exact check */
							if (nn->ix == vx->nix[e])
								break;
						}
						if (e <= di)
							par = 1;
					}

					eperr = ofps_comp_eperr7(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);

					if (!par && (vx->eperr - eperr) > 0.0) {
//printf("~1 Node ix %d at %s (%s)\n   Cell ix %d co %s center %s (%s),\n   vtx no %d at %s (%s)\n",nn->ix, ppos(di,nn->p),ppos(di,nn->v),cp - s->grid,pco(s->di,cp->co),ppos(di,cp->cp),ppos(di,cp->cv),vx->no, ppos(di,vx->p),ppos(di,vx->v));
						warning("Sanity check ofps_findhit_vtxs() cell skip failed, hit on vtx no %d, eperr %f < vx->eperr %f, cell ix %d eperr %f, est min eperr %f",vx->no,eperr,vx->eperr,cp - s->grid,ceperr,ceperr - cp->eperr);
						printf("Sanity check ofps_findhit_vtxs() cell skip failed, hit on vtx no %d, eperr %f < vx->eperr %f, cell ix %d eperr %f, est min eperr %f\n",vx->no,eperr,vx->eperr,cp - s->grid,ceperr,ceperr - cp->eperr);
#ifdef SANITY_CHECK_CLOSEST_FATAL
						error("findclosest node cell skip failed");
#endif
					}
				}
#endif /* SANITY_CHECK_CLOSEST */
				continue;			/* Cell is not worth searching */	
			}

			/* Search the cell */
			s->ncellssch++;

			/* For vertexes in this cell */
			for (vx = cp->vhead; vx != NULL; vx = vx->n) {
#ifdef DEBUG
				printf("Checking vtx no %d\n",vx->no);
#endif
				/* If the vertex has already been checked */
				if (vx->cflag == s->flag)
					continue;

				if (vx->ofake)		/* ofake vertexes can't be hit */
					continue;

#ifdef INDEP_SURFACE
				/* Only check for hit if the vertex is visible to the node */
				if (sm_vtx_node(s, vx, nn) == 0) {
# ifdef DEBUG
					printf("Vertex no %d xmask 0x%x vm %s isn't visible to ix %d pmask 0x%x a_sm %s\n",vx->no,vx->cmask,psm(s,&vx->vm),nn->ix,nn->pmask,psm(s,&s->sc[nn->pmask].a_sm));
# endif	/* DEBUG */
					continue;
				}
#endif	/* INDEP_SURFACE */
	
				vx->add = 0;
				vx->del = 0;
				vx->par = 0;

				s->vvchecks++;			/* Checking a vertex */

				/* Check if node is already parent to this vertex. */
				/* This only happens during fixups if the reposition fails and we */
				/* retain the vertex with the deleted vertex location (not currently */
				/* done), or by slim numerical margine, so ignore such hits. */
				/* We treat a parent as a hit node for the purposes of recursion, */
				/* and add it to a special list used to complete the vertex net. */
				if (nn->ixm & vx->nix[MXPD+2]) {	/* Is in nixm */
					for (e = 0; e <= di; e++) {		/* Do exact check */
						if (nn->ix == vx->nix[e])
							break;
					}
					if (e <= di) {
#ifdef DEBUG
						printf("Vertex no %d has already got node ix %d\n",vx->no,nn->ix);
#endif
						vx->par = 1;
					}
				}

				/* nba_eperr is assumed to be valid if vx->cflag == s->flag */
				vx->nba_eperr = ofps_comp_eperr7(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);
#ifdef DEBUG
				printf("Computing nba_eperr of %f for vtx no %d\n",vx->nba_eperr, vx->no);
#endif
				/* See if the vertex eperr will be improved */
				if (!vx->par && (vx->eperr - vx->nba_eperr) > 0.0) {
					s->nvcheckhits++;
					hit = 1;
					vx->del = 1;		/* Mark for deletion */
					vx->nxh = s->nxh;	/* Add vertex to list */
					s->nxh = vx;
					vx->hflag = s->flag;
#ifdef DEBUG
					printf("Vertex error improvement hit by %f (%f < %f)\n",vx->eperr-vx->nba_eperr,vx->nba_eperr,vx->eperr);

					if (vx->par) {
						printf("Vertex no %d hit by its own parent ix %d\n",vx->no, nn->ix);
						warning("Vertex no %d hit by its own parent ix %d",vx->no, nn->ix);
					}
#endif
				}
#ifdef DEBUG
				  else {			/* If worse */
					printf("Vertex error not hit by %f (%f < %f)\n",vx->eperr-vx->nba_eperr,vx->nba_eperr,vx->eperr);
				}
#endif
				vx->cflag = s->flag;

			}	/* Next vertex in cell */

			/* Put all this cells neighbours on the search list */
			/* (This is probably the critical inner loop. If ->acnl was */
			/* scaled by sizeof(acell), then the implicit multiply could */
			/* be avoided) */
			for (j = 0; j < s->nacnl; j++) {
				acell *nc = cp + s->acnl[j];
		
				if (nc->gflag >= s->gflag)
					continue;

#ifdef DEBUG_FCLOSE
				printf("Adding cell ix %d co %s to slist\n",nc - s->grid, pco(di,nc->co));
#endif
#ifdef NEVER
				if (sliste == NULL) {		/* First in empty list */
					slist = nc;
				} else {
					sliste->slist = nc; 	/* Add to end of list */
				}
				sliste = nc;
				nc->slist = NULL;
#else
				nc->slist = slist;
				slist = nc;
#endif
				nc->gflag = s->gflag;	/* Cell is on list to be searched */
			}
		}	/* Next cell in current list */
//printf("~1 don that search list\n");
	}	/* Next list */
//printf("~1 no more search lists\n");

	return hit;
}

#ifdef DEBUG_FCLOSE
#undef DEBUG_FCLOSE
#endif


#undef DEBUG_FCLOSE

/* Given a vertex, locate the smallest eperr node. */
/* Return NULL if none, and if ceperr is not NULL, set it to the */
/* eperr to the returned node. */ 
/* (This only returns visible nodes.) */
static node *ofps_findclosest_node(ofps *s, double *ceperr, vtx *vx) {
	int e, di = s->di;
	int i, j;
	int pci;				/* Point cell index */
	acell *cp;
	double eperr, beperr = 1e300;	/* eperr of closest node */
	node *bno = NULL;			/* Closest node */
	acell *slist = NULL, *sliste = NULL;		/* Next to search list */

#ifdef DEBUG
	if (s->agrid_init == 0)
		error("ofps_findclosest_node() called befor agrid_init");
#endif

#ifdef DEBUG_FCLOSE
	printf("\nLocating closest node to vtx at p = %s, v = %s\n", ppos(di,vx->p), ppos(di,vx->v));
#endif

	s->nnfschd += s->np;		/* Number of nodes in a full search */
	s->naccsrch++;				/* Number of searches */

	/* Do a breadth first seed search for any better nodees, or until */
	/* we run out of cells that could improve on the current best. */

	/* Locate a starting cell using the grid */
	pci = ofps_point2cell(s, vx->v, vx->p);	/* Grid index of cell of interest */
	cp = &s->grid[pci];

	s->gflag++;					/* cell touched flag */ 

	/* Put the starting cell on the search list */
#ifdef DEBUG_FCLOSE
	printf("Adding cell ix %d co %s to slist\n",cp - s->grid, pco(di,cp->co));
#endif
#ifdef NEVER
	if (sliste == NULL) {		/* First in empty list */
		slist = cp;
	} else {
		sliste->slist = cp; 	/* Add to end of list */
	}
	sliste = cp;
	cp->slist = NULL;
#else
	cp->slist = slist;			/* Add it to start of list */
	slist = cp;	
#endif
	cp->gflag = s->gflag;	/* Cell is on list to be searched */

	/* until we run out of cells to search */
	for (;slist != NULL;) {
		acell *ncp;

		/* For each cell in the search list, check it and recursion. */
		for (cp = slist, slist = sliste = NULL; cp != NULL; cp = ncp) {
			double ceperr;
			ncp = cp->slist;

#ifdef DEBUG_FCLOSE
			printf("Checking cell ix %d co %s\n",cp - s->grid,pco(di,cp->co));
#endif

			/* Compute the eperr of the cell center to the vtx minus the estimated */ 
			/* largest eperr of any point within the cell to the center. */
			ceperr = ofps_comp_eperr(s, NULL, cp->v, cp->p, vx->v, vx->p, vx->nsp);
			eperr = ceperr - cp->eperr;
			
			/* If the cell is worth searching */
			if (eperr < beperr) {
				node *no;

				/* Search the cell */
				s->ncellssch++;

				for (no = cp->head; no != NULL; no = no->n) {

#ifdef INDEP_SURFACE
					/* Check if this node is visible to this vtx */
					if (sm_vtx_node(s, vx, no) == 0) {
						continue;	/* It's hidden */
					}
#endif	/* INDEP_SURFACE */

					/* Compute the eperr between the node to the new vtx */
					eperr = ofps_comp_eperr(s, NULL, no->v, no->p, vx->v, vx->p, vx->nsp);
					if (eperr < beperr) {
						bno = no;
						beperr = eperr;
#ifdef DEBUG_FCLOSE
						printf("Improved to node ix %d eperr\n",bno->ix,beperr);
#endif
					}
				}

				/* Put all this cells neighbours on the search list */
				/* (This is probably the critical ivxer loop. If ->acnl was */
				/* scaled by sizeof(acell), then the implicit multiply could */
				/* be avoided) */
				for (j = 0; j < s->nacnl; j++) {
					acell *nc = cp + s->acnl[j];
			
					if (nc->gflag >= s->gflag)
						continue;

#ifdef DEBUG_FCLOSE
					printf("Adding cell ix %d co %s to slist\n",nc - s->grid, pco(di,nc->co));
#endif
#ifdef NEVER
					if (sliste == NULL) {		/* First in empty list */
						slist = nc;
					} else {
						sliste->slist = nc; 	/* Add to end of list */
					}
					sliste = nc;
					nc->slist = NULL;
#else
					nc->slist = slist;			/* Add it to start of list */
					slist = nc;
#endif
					nc->gflag = s->gflag;	/* Cell is on list to be searched */
				}
			}
#ifdef SANITY_CHECK_CLOSEST
			/* Check all nodees in the cell anyway */
			  else {
				double teperr;
				node *no;

				for (no = cp->head; no != NULL; no = no->n) {

#ifdef INDEP_SURFACE
					/* Check if this node is visible to this vtx */
					if (sm_vtx_node(s, vx, no) == 0) {
						continue;	/* It's hidden */
					}
#endif	/* INDEP_SURFACE */

					/* Compute the eperr between the node to the new vtx */
					teperr = ofps_comp_eperr(s, NULL, no->v, no->p, vx->v, vx->p, vx->nsp);
					if (teperr < beperr) {
						warning("Sanity check ofps_findclosest_node() cell skip failed, estimated %f from cellc eperr %f - cell eperr %f, found %f from node ix %d",eperr,ceperr,cp->eperr,teperr,no->ix);
						printf("Sanity check ofps_findclosest_node() cell skip failed, estimated %f from cellc eperr %f - cell eperr %f, found %f from node ix %d\n",eperr,ceperr,cp->eperr,teperr,no->ix);
#ifdef SANITY_CHECK_CLOSEST_FATAL
						error("findclosest node cell skip failed");
#endif
					}
				}
			}
#endif	/* SANITY_CHECK_CLOSEST */

		}	/* Next cell in current list */
#ifdef DEBUG_FCLOSE
		printf("Finished ivxer loop because p 0x%x = NULL\n",cp);
#endif
	}	/* Next list */
#ifdef DEBUG_FCLOSE
	printf("Finished outer loop because slist 0x%x = NULL\n",slist);
#endif

#ifdef DEBUG_FCLOSE
	if (bno == NULL)
		printf("Failed to find a closest node");
	else
		printf("Returning best node ix %d, eperr %f\n",bno->ix,beperr);
#endif

#ifdef SANITY_CHECK_CLOSEST
	/* Use exaustive search */
	{
		double ch_beperr = 1e300;	/* Device distance squared of closest vertex */
		node *ch_bno = NULL;
		for (i = 0; i < (s->np-1); i++) {
			node *nn = s->n[i];
			double eperr;

#ifdef INDEP_SURFACE
			/* Check if this vertex is visible to this node */
			if (sm_vtx_node(s, vx, nn) == 0) {
				continue;	/* It's hidden */
			}
#endif	/* INDEP_SURFACE */

			/* Compute the eperr between the node and the vertex */
			eperr = ofps_comp_eperr(s, NULL, nn->v, nn->p, vx->v, vx->p, vx->nsp);
			if (eperr < ch_beperr) {
				ch_bno = nn;
				ch_beperr = eperr;
			}
		}

		if (ch_bno != NULL && ch_beperr + 1e-3 < beperr) {
			if (bno == NULL) {
				warning("Sanity check ofps_findclosest_node() failed,\n   found none, should be ix %d dist %f",ch_bno->ix,ch_beperr); 
				printf("Sanity check ofps_findclosest_node() failed,\n   found none, should be ix %d dist %f\n",ch_bno->ix,ch_beperr); 
			} else {
				warning("Sanity check ofps_findclosest_node() failed,\n   found ix %d dist %f, should be ix %d dist %f",bno->ix,beperr,ch_bno->ix,ch_beperr); 
				printf("Sanity check ofps_findclosest_node() failed,\n   found ix %d dist %f, should be ix %d dist %f\n",bno->ix,beperr,ch_bno->ix,ch_beperr); 
			}
#ifdef SANITY_CHECK_CLOSEST_FATAL
			error("findclosest node failed");
#endif
		}
	}
#endif

	if (bno != NULL && ceperr != NULL)
		*ceperr = beperr;

	return bno;
}

/* ----------------------------------------------------------- */

#ifdef NEVER		/* No longer used */

/* Given a node, locate the smallest eperr vertex. */
/* Return NULL if none, and if ceperr is not NULL, set it to the */
/* eperr to the returned vertex. */ 
/* (This only returns visible vertexes.) */
static vtx *ofps_findclosest_vtx(ofps *s, double *ceperr, node *nn) {
	int e, di = s->di;
	int i, j;
	int pci;				/* Point cell index */
	acell *cp;
	double eperr, beperr = 1e300;	/* eperr of closest vertex */
	vtx *bvx = NULL;			/* Closest vertex */
	acell *slist = NULL, *sliste = NULL;		/* Next to search list */

#ifdef DEBUG
	if (s->agrid_init == 0)
		error("ofps_findclosest_vtx() called befor agrid_init");
#endif

#ifdef DEBUG_FCLOSE
	printf("\nLocating closest vtx to node at p = %s, v = %s\n", ppos(di,nn->p), ppos(di,nn->v));
#endif

	s->nvfschd += s->nv;		/* Number of vertexes in a full search */
	s->naccsrch++;				/* Number of searches */

	/* Do a breadth first seed search for any better vertexes, or until */
	/* we run out of cells that could improve on the current best. */

	/* Locate a starting cell using the grid */
	pci = ofps_point2cell(s, nn->v, nn->p);	/* Grid index of cell of interest */
	cp = &s->grid[pci];

	s->gflag++;					/* cell touched flag */ 

	/* Put the starting cell on the search list */
#ifdef DEBUG_FCLOSE
	printf("Adding cell ix %d co %s to slist\n",cp - s->grid, pco(di,cp->co));
#endif
	if (sliste == NULL) {		/* First in empty list */
		slist = cp;
	} else {
		sliste->slist = cp; 	/* Add to end of list */
	}
	sliste = cp;
	cp->slist = NULL;
	cp->gflag = s->gflag;	/* Cell is on list to be searched */

	/* until we run out of cells to search */
	for (;slist != NULL;) {
		acell *ncp;

		/* For each cell in the search list, check it and recursion. */
		for (cp = slist, slist = sliste = NULL; cp != NULL; cp = ncp) {
			double ceperr;
			ncp = cp->slist;

#ifdef DEBUG_FCLOSE
			printf("Checking cell ix %d co %s\n",cp - s->grid,pco(di,cp->co));
#endif

			/* Compute the eperr of the cell center to the node minus the estimated */ 
			/* largest eperr of any point within the cell to the center. */
			ceperr = ofps_comp_eperr(s, NULL, cp->v, cp->p, nn->v, nn->p, nn->nsp);
			eperr = ceperr - cp->eperr;
			
			/* If the cell is worth searching */
			if (eperr < beperr) {
				vtx *vx;

				/* Search the cell */
				s->ncellssch++;

				for (vx = cp->vhead; vx != NULL; vx = vx->n) {

#ifdef INDEP_SURFACE
					/* Check if this vertex is visible to this node */
					if (sm_vtx_node(s, vx, nn) == 0) {
						continue;	/* It's hidden */
					}
#endif	/* INDEP_SURFACE */

					/* Compute the eperr between the vertex to the new node */
					eperr = ofps_comp_eperr(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);
					if (eperr < beperr) {
						bvx = vx;
						beperr = eperr;
#ifdef DEBUG_FCLOSE
						printf("Improved to vtx no %d eperr\n",bvx->no,beperr);
#endif
					}
				}

				/* Put all this cells neighbours on the search list */
				/* (This is probably the critical inner loop. If ->acnl was */
				/* scaled by sizeof(acell), then the implicit multiply could */
				/* be avoided) */
				for (j = 0; j < s->nacnl; j++) {
					acell *nc = cp + s->acnl[j];
			
					if (nc->gflag >= s->gflag)
						continue;

#ifdef DEBUG_FCLOSE
					printf("Adding cell ix %d co %s to slist\n",nc - s->grid, pco(di,nc->co));
#endif
					if (sliste == NULL) {		/* First in empty list */
						slist = nc;
					} else {
						sliste->slist = nc; 	/* Add to end of list */
					}
					sliste = nc;
					nc->slist = NULL;
					nc->gflag = s->gflag;	/* Cell is on list to be searched */
				}
			}
#ifdef SANITY_CHECK_CLOSEST
			/* Check all vertexes in the cell anyway */
			  else {
				double teperr;
				vtx *vx;

				for (vx = cp->vhead; vx != NULL; vx = vx->n) {

#ifdef INDEP_SURFACE
					/* Check if this vertex is visible to this node */
					if (sm_vtx_node(s, vx, nn) == 0) {
						continue;	/* It's hidden */
					}
#endif	/* INDEP_SURFACE */

					/* Compute the eperr between the vertex to the new node */
					teperr = ofps_comp_eperr(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);
					if (teperr < beperr) {
						warning("Sanity check ofps_findclosest_vtx() cell skip failed, estimated %f from cellc eperr %f - cell eperr %f, found %f from vtx no %d",eperr,ceperr,cp->eperr,teperr,vx->no);
						printf("Sanity check ofps_findclosest_vtx() cell skip failed, estimated %f from cellc eperr %f - cell eperr %f, found %f from vtx no %d\n",eperr,ceperr,cp->eperr,teperr,vx->no);
#ifdef SANITY_CHECK_CLOSEST_FATAL
						error("findclosest vertex cell skip failed");
#endif
					}
				}
			}
#endif	/* SANITY_CHECK_CLOSEST */

		}	/* Next cell in current list */
#ifdef DEBUG_FCLOSE
		printf("Finished inner loop because p 0x%x = NULL\n",cp);
#endif
	}	/* Next list */
#ifdef DEBUG_FCLOSE
	printf("Finished outer loop because slist 0x%x = NULL\n",slist);
#endif

#ifdef DEBUG_FCLOSE
	if (bvx == NULL)
		printf("Failed to find a closest vertex");
	else
		printf("Returning best vtx no %d, eperr %f\n",bvx->no,beperr);
#endif

#ifdef SANITY_CHECK_CLOSEST
	/* Use exaustive search */
	{
		double ch_beperr = 1e300;	/* Device distance squared of closest vertex */
		vtx *vx, *ch_bvx = NULL;
		for (vx = s->uvtx; vx != NULL; vx = vx->link) {		/* Check all vertexes */
			double eperr;

#ifdef INDEP_SURFACE
			/* Check if this vertex is visible to this node */
			if (sm_vtx_node(s, vx, nn) == 0) {
				continue;	/* It's hidden */
			}
#endif	/* INDEP_SURFACE */

			/* Compute the eperr between the vertex to the new node */
			eperr = ofps_comp_eperr(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);
			if (eperr < ch_beperr) {
				ch_bvx = vx;
				ch_beperr = eperr;
			}
		}

		if (ch_bvx != NULL && ch_beperr + 1e-3 < beperr) {
			if (bvx == NULL) {
				warning("Sanity check ofps_findclosest_vtx() failed,\n   found none, should be no %d dist %f",ch_bvx->no,ch_beperr); 
				printf("Sanity check ofps_findclosest_vtx() failed,\n   found none, should be no %d dist %f\n",ch_bvx->no,ch_beperr); 
			} else {
				warning("Sanity check ofps_findclosest_vtx() failed,\n   found no %d dist %f, should be no %d dist %f",bvx->no,beperr,ch_bvx->no,ch_beperr); 
				printf("Sanity check ofps_findclosest_vtx() failed,\n   found no %d dist %f, should be no %d dist %f\n",bvx->no,beperr,ch_bvx->no,ch_beperr); 
			}
#ifdef SANITY_CHECK_CLOSEST_FATAL
			error("findclosest vertex failed");
#endif
		}
	}
#endif

	if (bvx != NULL && ceperr != NULL)
		*ceperr = beperr;

	return bvx;
}

/* Given a node, locate a vertex that it hits. */
/* Return NULL if none, and if ceperr is not NULL, set it to the */
/* eperr to the returned vertex. */ 
/* (This only returns visible vertexes.) */
static vtx *ofps_findhit_vtx(ofps *s, double *ceperr, node *nn) {
	int e, di = s->di;
	int i, j;
	int pci;				/* Point cell index */
	acell *cp;
	double eperr, beperr = 1e300;	/* eperr of closest vertex */
	vtx *bvx = NULL;			/* Closest vertex */
	acell *slist = NULL, *sliste = NULL;		/* Next to search list */

#ifdef DEBUG
	if (s->agrid_init == 0)
		error("ofps_findhit_vtx() called befor agrid_init");
#endif

#ifdef DEBUG_FCLOSE
	printf("\nLocating a hit vtx to node at p = %s, v = %s\n", ppos(di,nn->p), ppos(di,nn->v));
#endif

	s->nvfschd += s->nv;		/* Number of vertexes in a full search */
	s->naccsrch++;				/* Number of searches */

	/* Do a breadth first seed search for any hit vertexes, or until */
	/* we run out of cells that could improve on the current best. */

	/* Locate a starting cell using the grid */
	pci = ofps_point2cell(s, nn->v, nn->p);	/* Grid index of cell of interest */
	cp = &s->grid[pci];

	s->gflag++;					/* cell touched flag */ 

	/* Put the starting cell on the search list */
	for (j = 0; j < s->nacnl; j++) {
#ifdef DEBUG_FCLOSE
		printf("Adding cell ix %d co %s to slist\n",cp - s->grid, pco(di,cp->co));
#endif
		if (sliste == NULL) {		/* First in empty list */
			slist = cp;
		} else {
			sliste->slist = cp; 	/* Add to end of list */
		}
		sliste = cp;
		cp->slist = NULL;
		cp->gflag = s->gflag;	/* Cell is on list to be searched */
	}

	/* until we run out of cells to search */
	for (;slist != NULL;) {
		acell *ncp;

		/* For each cell in the search list, check it and recursion. */
		for (cp = slist, slist = sliste = NULL; cp != NULL; cp = ncp) {
			ncp = cp->slist;

#ifdef DEBUG_FCLOSE
			printf("Checking cell ix %d co %s\n",cp - s->grid,pco(di,cp->co));
#endif

			/* If the cell is worth searching */
			if (1) {
				vtx *vx;

				/* Search the cell */
				s->ncellssch++;

				for (vx = cp->vhead; vx != NULL; vx = vx->n) {

#ifdef INDEP_SURFACE
					/* Check if this vertex is visible to this node */
					if (sm_vtx_node(s, vx, nn) == 0) {
						continue;	/* It's hidden */
					}
#endif	/* INDEP_SURFACE */

					/* Compute the eperr between the vertex to the new node */
					eperr = ofps_comp_eperr(s, NULL, vx->v, vx->p, nn->v, nn->p, nn->nsp);
					if (eperr < vx->eperr) {
						bvx = vx;
						beperr = eperr;
#ifdef DEBUG_FCLOSE
						printf("Found hit vtx no %d eperr\n",bvx->no,beperr);
#endif
						break;
					}
				}
				if (vx != NULL)
					break;

				/* Put all this cells neighbours on the search list */
				/* (This is probably the critical inner loop. If ->acnl was */
				/* scaled by sizeof(acell), then the implicit multiply could */
				/* be avoided) */
				for (j = 0; j < s->nacnl; j++) {
					acell *nc = cp + s->acnl[j];
			
					if (nc->gflag >= s->gflag)
						continue;

#ifdef DEBUG_FCLOSE
					printf("Adding cell ix %d co %s to slist\n",nc - s->grid, pco(di,nc->co));
#endif
					if (sliste == NULL) {		/* First in empty list */
						slist = nc;
					} else {
						sliste->slist = nc; 	/* Add to end of list */
					}
					sliste = nc;
					nc->slist = NULL;
					nc->gflag = s->gflag;	/* Cell is on list to be searched */
				}
			}
		}	/* Next cell in current list */
	}	/* Next list */

#ifdef DEBUG_FCLOSE
	if (bvx == NULL)
		printf("Failed to find a hit vertex");
	else
		printf("Returning hit vtx no %d, eperr %f\n",bvx->no,beperr);
#endif

	if (bvx != NULL && ceperr != NULL)
		*ceperr = beperr;

	return bvx;
}

#ifdef DEBUG_FCLOSE
#undef DEBUG_FCLOSE
#endif

#endif /* NEVER */

/* ----------------------------------------------------------- */

/* Re-position the vertexes given the current point positions, */
/* and fixup the veronoi. */
static void
ofps_repos_and_fix_voronoi(
ofps *s
) {
	int e, di = s->di;
	int i, j, k;
	node *nds[MXPD+1];	/* Real nodes of vertex */
	int ii;				/* Number of real nodes */
	double ee[MXPD+1];	/* Per node estimated error */
	vtx *vx;
	node *nn, *pp;
	int nfuxups, l_nfuxups;		/* Count of fixups */
	int csllow;					/* Count since last low */
	int mxcsllow = 5;			/* Threshold to give up */

#ifdef DEBUG
	printf("Repositioning vertexes\n");
#endif

	/* Re-position the vertexes to match optimized node positions */
	s->fchl = NULL;
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		nodecomb nc;

		if (vx->ifake || vx->ofake)
			continue;

		vx->p_eperr = vx->eperr;

		/* Pointers to real nodes. */
		for (ii = e = 0; e <= di; e++) {
			if (vx->nix[e] >= 0)
				nds[ii++] = s->n[vx->nix[e]];
			if (vx->nix[e] < -s->nbp)
				error("ofps_repos_and_fix_voronoi() got fake node no %d comb %s fake %d",vx->no,pcomb(di,vx->nix),vx->ofake);
		}

		/* Compute the current eperr at the vertex given the repositioned nodes, */
		/* to set acceptance threshold for repositioned vertex. */
		ofps_pn_eperr(s, NULL, ee, vx->v, vx->p, nds, ii);
		nc.ceperr = ofps_eperr2(ee, ii);
 
		/* Setup to re-position the vertex */
		memset((void *)&nc, 0, sizeof(nodecomb));
		for (e = 0; e < di; e++) {
			nc.nix[e] = vx->nix[e];
			nc.p[e] = vx->p[e]; 
			nc.v[e] = vx->v[e]; 
		}
		nc.nix[e] = vx->nix[e];

#ifdef DEBUG
		printf("Repositioning vertex no %d nodes %s at %s, ceperr %f\n",vx->no,pcomb(di,vx->nix),ppos(di,vx->p),nc.ceperr);
#endif

		/* We're about to change the position and eperr: */
		ofps_rem_vacc(s, vx);
		ofps_rem_vseed(s, vx);

		if (position_vtx(s, &nc, 1, 1, 0) == 2) {
			/* Just leave it where it was. Perhaps fixups will delete it */
			if (s->verb > 1)
				warning("re_position_vtx failed for vtx no %d at %s",vx->no,ppos(di,vx->p));
		} else {
//printf("~1 moved from %s to %s\n",ppos(di,vx->p),ppos(di,nc.p));

			for (e = 0; e < di; e++) {
				vx->p[e] = nc.p[e];
				vx->v[e] = nc.v[e];
			}
			vx->eperr = nc.eperr;
			vx->eserr = nc.eserr;
		}

		/* Count the number of gamut surfaces the vertex falls on */
		det_vtx_gsurf(s, vx);

		/* We've changed the position and eperr: */
		ofps_add_vacc(s, vx);
		ofps_add_vseed(s, vx);

		/* Add all vertexes to the "to be checked" list */
		vx->fchl = s->fchl; /* Add vertex to the "to be checked" list */
		if (s->fchl != NULL)
			s->fchl->pfchl = &vx->fchl;
		s->fchl = vx;
		vx->pfchl = &s->fchl;
		vx->fflag = s->fflag; 
		vx->fupcount = 0;
		vx->fuptol = NUMTOL;
	}

#ifdef DUMP_PLOT_BEFORFIXUP
	printf("Before applying fixups:\n");
	dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 0, -1);		/* Device, No wait, verticies */
#endif /* DUMP_PLOT_BEFORFIXUP */

	/* Now fixup the veroni. */
#ifdef DEBUG
	printf("Doing fixups:\n");
#endif

	/* We loop until the check list is empty */
	l_nfuxups = 1e9;
	csllow = 0;
	while (s->fchl != NULL && csllow < mxcsllow) {
		vtx *nvx;

		s->fflag++;			/* Fixup round flag */
		s->nsvtxs = 0;
		nfuxups = 0;

#ifdef DEBUG
		printf("\nFixup round %d%s\n",s->fflag, s->fchl == NULL ? "" : " fchl != NULL");
#endif

		/* out of gamut, or whether the closest node to it */
		/* is not one of its parent nodes. */ 
		for (vx = s->fchl; vx != NULL; vx = nvx) {
			double ceperr;		/* eperr to closest node */
			int hit = 0;

		/* For each vertex on the check list, check if it is */
			nvx = vx->fchl;
#ifdef DEBUG
			printf("Checking vtx no %d, fuptol %e\n",vx->no,vx->fuptol);
#endif

			vx->hnode = NULL;
			nn = NULL;
			/* Check if the vertex position is clipped, */
			/* and add fake boundary node if it is */
			/* For all the gamut boundary planes: */
			for (i = 0; i < s->nbp; i++) {
				pleq *vp = &s->gpeqs[i];
				double v;
		
				nn = s->n[-1-i];

#ifdef INDEP_SURFACE
				/* Check if this vertex is visible to this node */
				if (sm_vtx_node(s, vx, nn) == 0) {
					continue;
				}
#endif	/* INDEP_SURFACE */

				for (v = vp->pe[di], e = 0; e < di; e++)
					v += vp->pe[e] * vx->p[e];
				if (v > vx->fuptol) {

#ifdef NEVER
					/* This is an optimization: */
					/* Check whether nn is already a parent of the node */
					for (e = 0; e <= di; e++) {
						if (nn->ix == vx->nix[e])
							break;
					}
					if (e <= di) {
						continue;			/* It is */
					}
#endif

					/* Add all the vertexes parent nodes to the nearest nodes "add" list */
#ifdef DEBUG
					printf("Vertex no %d hit by boundary node ix %d by %e\n",vx->no,nn->ix,v);
#endif
					hit = 1;
					if (vx->hnode == NULL) { 
						vx->hnode = nn;
						vx->hitmarg = 50.0 * v;

						if (s->nsvtxs >= s->_nsvtxs) {
							s->_nsvtxs = 2 * s->_nsvtxs + 5;
							if ((s->svtxs = (vtx **)realloc(s->svtxs, sizeof(vtx *) * s->_nsvtxs)) == NULL)
								error("ofps: malloc failed on svtxs%d", s->_nsvtxs);
						}
						s->svtxs[s->nsvtxs]= vx;
						vx->psvtxs = &s->svtxs[s->nsvtxs];
						s->nsvtxs++;

					} else if (50.0 * v > vx->hitmarg) {
						vx->hnode = nn;
						vx->hitmarg = 50.0 * v;
					}
#ifdef DEBUG
					printf("Added vtx no %d to node %d for fixup\n",vx->no,nn->ix);
#endif
				}
			}

			/* Or locate the nearest node to the vertex. */
			/* (This only returns visible nodes) */
			if ((nn = ofps_findclosest_node(s, &ceperr, vx)) != NULL) {
				double errimp = vx->eperr - ceperr;

				/* See if it is closer than the parent nodes */
				if (errimp >= vx->fuptol) {		/* It is */

					/* Add the vertexe to the "to be fixed" list */
#ifdef DEBUG
					printf("Vertex no %d hit by node ix %d by %e\n",vx->no,nn->ix,errimp);
#endif
					hit = 1;

					if (vx->hnode == NULL) { 
						vx->hnode = nn;
						vx->hitmarg = errimp;

						if (s->nsvtxs >= s->_nsvtxs) {
							s->_nsvtxs = 2 * s->_nsvtxs + 5;
							if ((s->svtxs = (vtx **)realloc(s->svtxs, sizeof(vtx *) * s->_nsvtxs)) == NULL)
								error("ofps: malloc failed on svtxs%d", s->_nsvtxs);
						}
						s->svtxs[s->nsvtxs]= vx;
						vx->psvtxs = &s->svtxs[s->nsvtxs];
						s->nsvtxs++;

					} else if (errimp > vx->hitmarg) {
						vx->hnode = nn;
						vx->hitmarg = errimp;
					}
#ifdef DEBUG
					printf("Added node %d to vtx no %d fixup\n",nn->ix, vx->no);
#endif
				}
			}

		next_vtx:;
			if (hit) {
				vx->fupcount++;
				vx->fuptol *= 2.0;
				nfuxups++;
			}

			/* Remove this vertex from the check list */
			if (vx->pfchl != NULL) {		/* If is on fixup check list, remove it */
				*vx->pfchl = vx->fchl;
				if (vx->fchl != NULL)
					vx->fchl->pfchl = vx->pfchl;
			}
			vx->pfchl = NULL;
			vx->fchl = NULL;
		}
		if (s->fchl != NULL)
			error("Check list should be empty!");

		if (nfuxups < l_nfuxups) {
			l_nfuxups = nfuxups;
			csllow = 0;
		} else {
			csllow++;
		}

#ifdef DEBUG
		printf("\nAbout to fixup %d marked vertexes\n",nfuxups);
#endif
		/* Smallest error to largest seems best, */
		/* probably because the closer nodes cut off the */
		/* further ones, reducing the number of redundant create/deletes */
#define HEAP_COMPARE(A,B) ((A)->hitmarg < (B)->hitmarg)
		HEAPSORT(vtx *, s->svtxs, s->nsvtxs);
#undef HEAP_COMPARE

		/* Fixup the back references after the sort */
		for (i = 0; i < s->nsvtxs; i++) {
			vx = s->svtxs[i];
			vx->psvtxs = &s->svtxs[i];
		}

		/* For each vertex on the "to be fixed" list, */
		/* search for hits by the node starting at that vertex, */
		/* and recursively locate all the hit vertexes. */
		for (i = 0; i < s->nsvtxs; i++) {

			if ((vx = s->svtxs[i]) == NULL) {
				continue;			/* Vertex got deleted by a previous fix */
			}
			nn = vx->hnode; 
			s->svtxs[i] = NULL;
			vx->psvtxs = NULL;

			s->nvcheckhits = 0;	/* Count number of vertexes hit by recursive check. */
			s->batch = NULL;	/* Nothing in pending delete list */
			s->nup = NULL;		/* Nothing in nodes to be updated list */
			s->flag++;			/* Marker flag for adding this node */
			s->nxh = NULL;		/* Nothing in nodes hit list */

#ifdef DEBUG
			printf("\nFixing up node ix %d starting at vx no %d\n",nn->ix,vx->no);
//			fprintf(stderr,"Fixing up node ix %d starting at vx no %d\n",nn->ix,vx->no);
#endif
			/* Recursively search for all vertexes hit by the new node */
			/* Note that we don't care that this only finds connected hits, */
			/* since there should be a separate s->svtxs[] entry for a hit by this */
			/* node on a disconnected region. */
			ofps_check_vtx(s, nn, vx, 100000, 0);

#ifdef DEBUG
			printf("Fixing up node ix %d, %d vertexes hit by it\n",nn->ix,s->nvcheckhits);
#endif

			/* Number of nodes that would be checked by exaustive search */
			s->vvpchecks += s->nv;
			
			/* Now re-add the node to the veronoi */
			if (add_to_vsurf(s, nn, 1, 0) > 0) {
				s->add_hit++;
#ifdef DUMP_PLOT_EACHFIXUP
				printf("After adding node ix %d at %s to vurf\n",nn->ix,ppos(di,nn->p));
				dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 0, -1);		/* Device, No wait, verticies */
#endif /* DUMP_PLOT_EACHFIXUP */
			} else {
#ifdef DUMP_PLOT_EACHFIXUP
				printf("Adding node ix %d at %s to vurf was miss\n",nn->ix,ppos(di,nn->p));
#endif /* DUMP_PLOT_EACHFIXUP */
				s->fadd_mis++;
			}
		}
		s->nsvtxs = 0;
	}	/* Loop until there are no more vertexes to check */

#ifdef DEBUG
	printf("Done fixups s->fchl 0x%x == NULL or csllow %d >= %d\n",s->fchl,csllow,mxcsllow);
#endif

#ifdef SANITY_CHECK_FIXUP 
	/* Check that no node other than a parent is closer to any vertex */
	if (check_vertex_closest_node(s)) {
#ifdef SANITY_CHECK_FIXUP_FATAL
		error("!!!!!! Sanity: Fixup didn't work");
#endif /* SANITY_CHECK_FIXUP_FATAL */
	}
#endif /* SANITY_CHECK_FIXUP */

#ifdef DEBUG
		printf("Applied fixups\n");
#endif
}

/* --------------------------------------------------- */
/* After seeding or re-positioning, create the node */
/* neighbour node and vertex lists. */
/* (Invalidates and deletes any midpoints) */
static void ofps_re_create_node_node_vtx_lists(ofps *s) {
	int i, e, di = s->di;
	vtx *vx;

	/* for each node, clear its vertex list */
	for (i = -s->gnp; i < s->np; i++) {
		node *p = s->n[i];
		p->nvv = 0;
	}

	/* For each vertex, add it to each of its parent nodes */
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		for (e = 0; e <= di; e++) {
			node *p = s->n[vx->nix[e]];
			node_add_vertex(s, p, vx);
		}
	}

	/* For each node, recompute its neighbourhood nodes */
	for (i = -s->gnp; i < s->np; i++) {
		node *p = s->n[i];
		node_recomp_nvn(s, p);	/* Recompute the nodes associated vertex nodes */
	}
}

/* --------------------------------------------------- */
/* Midpoints */

/* Finding midpoint location code using dnsqe() */

/* Context for callback */
typedef struct {
	ofps *s;
	node *nds[2];	/* List of nodes */
} mopt_cx;

/* calculate the functions at x[] */
int dnsq_mid_solver(	/* Return < 0 on abort */
	void *fdata,	/* Opaque data pointer */
	int n,			/* Dimenstionality */
	double *x,		/* Multivariate input values */
	double *fvec,	/* Multivariate output values */
	int iflag		/* Flag set to 0 to trigger debug output */
) {
	mopt_cx *cx = (mopt_cx *)fdata;
	ofps *s = cx->s;
	int e, di = s->di;
	double pos[MXPD], sv[MXPD];
	double cee[2], teperr;

//printf("~1 dnsq_solver got %d nodes and %d planes\n",cx->nn,cx->nsp);

	/* Compute pos as interpolation between node 0 and 1 */
	for (e = 0; e < di; e++)
		pos[e] = cx->nds[0]->p[e] * (1.0 - x[0]) + cx->nds[1]->p[e] * x[0];

	ofps_cc_percept(s, sv, pos);

	/* Get eperr */
	cee[0] = ofps_comp_eperr8(s, NULL, sv, pos, cx->nds[0]->v, cx->nds[0]->p, cx->nds[0]->nsp);
	cee[1] = ofps_comp_eperr8(s, NULL, sv, pos, cx->nds[1]->v, cx->nds[1]->p, cx->nds[1]->nsp);

//printf("~1 error = %f, %f", cee[0], cee[1]);

	teperr = 0.5 * (cee[0] + cee[1]);

	fvec[0] = teperr - cee[0];

//	printf("dnsq_mid_solver returning %f from %f\n",fvec[0],x[0]);

	return 0;
}

/* Create or re-create all the midpoints, given the vertexes are done. */
static void
ofps_create_mids(ofps *s) {
	int e, di = s->di;
	int i, j, k;
	double rerr;
	int nsp = 0;		/* Number of surface points */
	double dnsqtol = 1e-6;		/* Solution tollerance to aim for */
	vopt_cx cx;
	double fvec[1];
	int rv;

	cx.s = s;

//printf("~1 creating mid points\n");
	/* Clear any existing midpoints */
	for (i = 0; i < s->tinp; i++) {
		node *p = s->n[i];

		if (p->ix < 0)
			break;		/* Done when we get to gamut boundary nodes */

		for (j = 0; j < p->nvn; j++) {
			if (p->mm[j] != NULL) {
				del_mid(s, p->mm[j]);
				p->mm[j] = NULL;
			}
		}
	}

	/* For each node, make sure it and each neighbor node have a shared midpoint */ 
	for (i = 0; i < s->tinp; i++) {
		node *p = s->n[i];

		if (p->ix < 0)
			break;		/* Done when we get to gamut boundary nodes */

		/* For each neighbor node, create midpoint */
		for (j = 0; j < p->nvn; j++) {
			mid *mp;
			node *p2;
			double ee[2];

			if (p->vn[j] < 0 || p->mm[j] != NULL)
				continue;	/* Gamut boundary or already got a midpoint */

			/* Create a midpoint between node p->ix and p->vn[j] */
			p2 = s->n[p->vn[j]];
			mp = new_mid(s);
			mp->refc++;

			p->mm[j] = mp;
			for (k = 0; k < p2->nvn; k++) {
				if (p2->vn[k] == p->ix) {
					p2->mm[k] = mp;
					mp->refc++;
					break;
				}
			}

			mp->nix[0] = p->ix;
			mp->nix[1] = p->vn[j];
//printf("~1 creating midpoint %d between nodes %d %d\n",mp->no,p->ix,p->vn[j]);
				
			cx.nds[0] = p;
			cx.nds[1] = p2;
			mp->np = 0.5;

			/* Locate mid point */
			if ((rv = dnsqe((void *)&cx, dnsq_mid_solver, NULL, 1, &mp->np,
			                0.2, fvec, 0.0, dnsqtol, 0, 0)) != 1 && rv != 3) {
				error("ofps: Locating midpoint failed with %d",rv);
			}

			for (e = 0; e < di; e++)
				mp->p[e] = p->p[e] * (1.0 - mp->np) + p2->p[e] * mp->np;
			ofps_cc_percept(s, mp->v, mp->p);

			/* Compute the eperr's for midpoint */
			ofps_pn_eperr(s, mp->ce, ee, mp->v, mp->p, cx.nds, 2);
			mp->eperr = ofps_eperr2(ee, 2);
			mp->eserr = ofps_eserr2(mp->ce, ee, 2);
//printf("~1 location %s (%s) eperr %f\n",ppos(di,mp->p), ppos(di,mp->v), mp->eperr);
		}
	}
}

/* --------------------------------------------------- */
/* Statistics: Compute serr stats. */

static void ofps_stats(ofps *s) {
	int e, di = s->di;
	int i, j;
	double acnt;
	vtx *vx;
	mid *mp;

//printf("~1 stats called\n");
	s->mn = 1e80;
	s->mx = -1e80;
	s->av = 0.0;
	acnt = 0.0;

	/* Vertex stats */
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		double es;

		if (vx->ghost)		/* Skip a ghost (coincident) vertex */
			continue;

#ifdef INDEP_SURFACE
		/* Ignore vertexes that aren't full dimension. */
	    if (sm_andtest(s, &s->sc[0].a_sm, &vx->vm) == 0)
			continue;
#endif
		es = vx->eserr;
		if (es >= 0.0 && es < s->mn)
			s->mn = es;
		if (es > s->mx)
			s->mx = es;
		s->av += es;
		acnt++;
	}

	s->av /= acnt;

	/* Midpoint/node stats */
	for (s->smns = 1e6, mp = s->umid; mp != NULL; mp = mp->link) {

		if (mp->nix[0] < 0 || mp->nix[1] < 0
		 || mp->eserr < 0.0)
			continue;		/* Skip fake points */

		if (mp->eserr < s->smns) {
			s->smns = mp->eserr;
		}
	}
	s->smns *= 2.0;		/* Error distance between nodes is double error to midpoint */
}

/* --------------------------------------------------- */
/* Support accessing the list of generated sample points */

/* Reset the read index */
static void
ofps_reset(ofps *s) {
	s->rix = 0;
}

/* Read the next non-fixed point value */
/* Return nz if no more */
static int
ofps_read(ofps *s, double *p, double *v) {
	int e;

	/* Advance to next non-fixed point */
	while(s->rix < s->np && s->n[s->rix]->fx)
		s->rix++;
	
	if (s->rix >= s->np)
		return 1;

	/* Return point info to caller */
	for (e = 0; e < s->di; e++) {
		if (p != NULL)
			p[e] = s->n[s->rix]->p[e];
		if (v != NULL)
			v[e] = s->n[s->rix]->v[e];
	}
	s->rix++;

	return 0;
}

/* --------------------------------------------------- */

/* Compute more optimum location for node amongst the surrounding */
/* verticies. The result is put in ->np[] and ->nv[]. */
/* The main aim is to minimize the maximum eserr of any vertex, */
/* but moving away from low midpoint eserr's improves the */
/* convergence rate and improves the eveness of the result. */
static void comp_opt(ofps *s, int poi, double oshoot, double sep_weight) {
	node *pp;						/* Node in question */
	int e, di = s->di;
	double radsq = -1.0;			/* Span/radius squared */
	double rad;
	double sum;
	int i;
	int bi = 0, bj = 0;

	pp = s->n[poi];		/* Node in question */

	/* Move towards vertex with highest eserr approach */
	if (pp->nvv > 0) {
		double aerr1, werr1, berr1, cnt1;	/* Average, worst, best eserr from vertexes */
		int weix1, beix1;			/* Worst, best error vertex index */
		double ov1[MXPD];			/* Optimization vector towards largest vertex error */
		double aerr2, werr2, berr2, cnt2;	/* Average, worst, best eserr from midpoints */
		int weix2, beix2;			/* Worst, best error midpoint index */
		double ov2[MXPD];			/* Optimization vector away from smallest midpoint error */


//printf("\n --------------------------------\n");
//printf("~1 Optimizing ix %d, %f %f\n",poi,pp->p[0],pp->p[1]);

		/* Compute the average and locate the largest error from verticies */
		for (aerr1 = cnt1 = 0.0, werr1 = -1.0, berr1 = 1e80, i = 0; i < pp->nvv; i++) {
			vtx *vp = pp->vv[i];

#ifdef INDEP_SURFACE
			/* Ingnore vertexes that are not visible to this node. */
			if (sm_vtx_node(s, vp, pp) == 0) {
				continue;
			}
#endif
			aerr1 += vp->eserr;
			cnt1++;

//printf("~1 Vertex no %d at %f %f serr = %f\n",vp->no,vp->p[0],vp->p[1],vp->eserr);
			if (vp->eserr > werr1) {
				werr1 = vp->eserr;
				weix1 = i;
			}
			if (vp->eserr < berr1) {
				berr1 = vp->eserr;
				beix1 = i;
			}
		}

		if (cnt1 > 0.0 && werr1 > NUMTOL && berr1 > NUMTOL) {
			double wbf, bbf;
			double towards = 0.8;	/* Amount to weight vector towards from closest */

			/* Compute a blend factor that takes the current */
			/* location towards the worst vertex error and */
			/* away from the best */
			aerr1 /= cnt1;
			wbf = towards         * (werr1 - aerr1)/werr1;
			bbf = (1.0 - towards) * (aerr1 - berr1)/berr1;
//			wbf = towards         * (werr1 - aerr1)/aerr1;
//			bbf = (1.0 - towards) * (aerr1 - berr1)/aerr1;

			for (e = 0; e < di; e++)
				ov1[e] = wbf * (pp->vv[weix1]->p[e] - pp->p[e])
				       + bbf * (pp->p[e] - pp->vv[beix1]->p[e]);
//printf("~1 moved %f %f towards vtx no %d at %f %f\n",ov1[0],ov1[2],pp->vv[weix1]->no,pp->vv[weix1]->p[0],pp->vv[weix1]->p[1]);
		} else {
			for (e = 0; e < di; e++)
				ov1[e] = 0.0;
		}

		/* Compute the average and locate the smallest error from midpoints */
		for (aerr2 = cnt2 = 0.0, werr2 = 1e80, berr2 = -1.0, i = 0; i < pp->nvn; i++) {
			mid *mp = pp->mm[i];
			node *on = s->n[pp->vn[i]];		/* Other node involved */

			if (mp == NULL || mp->nix[0] < 0 || mp->nix[1] < 0)
				continue;			/* Must be a fake gamut boundary node */

#ifdef INDEP_SURFACE
			/* Ingnore nodes of higher dimension */
			if ((pp->pmask & on->pmask) != pp->pmask) {
				continue;
			}
#endif
			aerr2 += mp->eserr;
			cnt2++;
//printf("~1 plane no %d from node ix %d serr = %f\n",mp->no,on->ix,mp->eserr);
			if (mp->eserr < werr2) {
				werr2 = mp->eserr;
				weix2 = i;
			}
			if (mp->eserr > berr2) {
				berr2 = mp->eserr;
				beix2 = i;
			}
		}

		if (cnt2 > 0.0 && werr2 > NUMTOL && berr2 > NUMTOL) {
			double wbf, bbf;
			double away = 0.8;		/* Amount to weight vector away from closest */

			/* Compute a blend factor that takes the current */
			/* location away from the worst plane error */
			aerr2 /= cnt2;
			wbf = away         * (aerr2 - werr2)/werr2;
			bbf = (1.0 - away) * (berr2 - aerr2)/berr2;
//			wbf = away         * (aerr2 - werr2)/aerr2;
//			bbf = (1.0 - away) * (berr2 - aerr2)/aerr2;

			for (e = 0; e < di; e++)
				ov2[e] = wbf * (pp->p[e] - pp->mm[weix2]->p[e])
				       + bbf * (pp->mm[beix2]->p[e] - pp->p[e]);
//printf("~1 moved %f %f away from node ix %d at %f %f\n",ov2[0],ov2[1],pp->vn[weix2],pp->mm[weix2]->p[0],pp->mm[weix2]->p[1]);
		} else {
			for (e = 0; e < di; e++)
				ov2[e] = 0.0;
		}
//printf("~1 ov1 = %f %f, ov2 = %f %f, sep weight %f\n",ov1[0], ov1[1], ov2[0], ov2[1], sep_weight);

		/* Move the node by the sum of the two vectors */ 
		for (e = 0; e < di; e++)
			pp->np[e] = pp->p[e] + (1.0 - sep_weight) * ov1[e] + sep_weight * ov2[e];
//printf("~1 moved node %d by %f %f\n",pp->ix, (1.0 - sep_weight) * ov1[0] + sep_weight * ov2[0],(1.0 - sep_weight) * ov1[1] + sep_weight * ov2[1]);
	}

//printf("~1 check moved by %f %f\n",pp->np[0] - pp->p[0], pp->np[1] - pp->p[1]);
	/* Apply overshoot/damping */
	for (e = 0; e < di; e++)
		pp->np[e] = pp->p[e] + (pp->np[e] - pp->p[e]) * oshoot;
//printf("~1 after overshoot of %f got %f %f\n",oshoot,pp->np[0],pp->np[1]);

	/* Clip the new location */
	ofps_clip_point10(s, pp->np, pp->np);

#if defined(KEEP_SURFACE) || defined(INDEP_SURFACE)
	if (pp->nsp > 0) {
		confineto_gsurf(s, pp->np, pp->sp, pp->nsp);
	}
#endif
	/* Update perceptual */
	s->percept(s->od, pp->nv, pp->np);	/* Was clipped above */

	/* Compute how far the point has moved */
	/* (?? maybe should change this to change in average or max eserr ??) */
	for (sum = 0.0, e = 0; e < di; e++) {
		double tt = pp->np[e] - pp->p[e];
		sum += tt * tt;
//printf("~1 total motion = %f\n",sqrt(sum));
	}
	if (sum > s->mxmvsq)	/* Track maximum movement */
		s->mxmvsq = sum;
}

static void
ofps_optimize(
ofps *s
) {
	int maxits;
	int transitters;
	double transpow;
	double oshoot, ioshoot, foshoot;
	double sepw, isepw, fsepw;
	double stoptol;
	int e, di = s->di;
	int i, j;

	/* Default is "Good" */
	maxits = OPT_MAXITS;
	transitters = OPT_TRANS_ITTERS;
	transpow = OPT_TRANS_POW;
	ioshoot = OPT_INITIAL_OVERSHOOT;
	foshoot = OPT_FINAL_OVERSHOOT;
	isepw = OPT_INITIAL_SEP_WEIGHT;
	fsepw = OPT_FINAL_SEP_WEIGHT;
	stoptol = OPT_STOP_TOL;

#ifdef OPT_MAXITS_2
	/* Option is "Fast" */
	if (s->good == 0) {
		maxits = OPT_MAXITS_2;
		transitters = OPT_TRANS_ITTERS_2;
		transpow = OPT_TRANS_POW_2;
		ioshoot = OPT_INITIAL_OVERSHOOT_2;
		foshoot = OPT_FINAL_OVERSHOOT_2;
		isepw = OPT_INITIAL_SEP_WEIGHT_2;
		fsepw = OPT_FINAL_SEP_WEIGHT_2;
		stoptol = OPT_STOP_TOL_2;
	} 
#endif /* OPT_MAXITS_2 */

	oshoot = ioshoot;
	for (s->optit = 0; s->optit < maxits; s->optit++) {	/* Up to maximum number of itterations */
		vtx *vx;
		double bf = 1.0;
		int nvxhits;
		double hratio, thresh;
		int doinc = 0;

		s->mxmvsq = 0.0;
		if (s->optit < transitters)
			bf = s->optit/(double)transitters;
		bf = pow(bf, transpow);
		oshoot = (1.0 - bf) * ioshoot + bf * foshoot;
		sepw = (1.0 - bf) * isepw  + bf * fsepw;

		/* Compute optimized node positions */
		for (i = 0; i < s->tinp; i++) {

			if (s->n[i]->fx)
				continue;		/* Ignore fixed points */ 

			comp_opt(s, i, oshoot, sepw);
		}

		/* Then update their positions to the optimized ones */
		for (i = 0; i < s->tinp; i++) {
			node *pp = s->n[i];

			if (pp->fx)
				continue;		/* Ignore fixed points */ 

			ofps_rem_nacc(s, pp);			/* Remove from spatial accelleration grid */
			
			for (e = 0; e < di; e++) {
				pp->op[e] = pp->p[e];		/* Record previous position */
				pp->p[e] = pp->np[e];		/* Move to optimized location */
				pp->v[e] = pp->nv[e];
			}
			ofps_add_nacc(s, pp);			/* Add to spatial acceleration grid */
		}

		/* Make sure that the optimized nodes don't accidentaly collide */
		for (i = 0; i < s->tinp; i++) {
			node *pp = s->n[i];

			if (pp->fx)
				continue;		/* Ignore fixed points */ 

			for (j = 0; j < 20; j++) {	/* Retry until not cooincident */
				int pci;				/* Point list index */
				acell *cp;				/* Acceleration cell */
				node *p1;
				pci = ofps_point2cell(s, pp->v, pp->p);		/* Grid index of cell of interest */

				cp = &s->grid[pci];
				for (p1 = cp->head; p1 != NULL; p1 = p1->n) {
					if (p1 == pp)
						continue;
					for (e = 0; e < di; e++) {
						if (fabs(pp->p[e] - p1->p[e]) > COINTOL)
							break;		/* Not cooincident */
					}
					if (e >= di) { 	/* Cooincident */
#ifdef DEBUG
						printf("Optimized node ix %d at %s collides with ix %d at %s - joggling it %d\n",pp->ix,ppos(di,pp->p),p1->ix,ppos(di,p1->p),i);
						warning("Optimized node ix %d at %s collides with ix %d at %s - joggling it %d",pp->ix,ppos(di,pp->p),p1->ix,ppos(di,p1->p),i);
#endif
						ofps_rem_nacc(s, pp);			/* Remove from spatial accelleration grid */

						/* Joggle it's position */
						for (e = 0; e < di; e++) {
							if (pp->p[e] < 0.5)
								pp->p[e] += d_rand(0.0, 1e-4);
							else
								pp->p[e] -= d_rand(0.0, 1e-4);
						}
						/* Ignore confine planes. Next itter should fix it anyway ? */
						ofps_clip_point10(s, pp->p, pp->p);

						/* Update perceptual (was clipped above) */
						s->percept(s->od, pp->v, pp->p);

						break;
					}
				}
				if (p1 == NULL)
					break;
			}
			if (j >= 20)
				error("ofps_optimize: Assert, was unable to joggle cooincindent point");
		}

		/* Ideally the fixup method should create and delete fewer vertexes */
		/* than reseeding, hence always be faster, but in practice this doesn't */
		/* seem to be so. Perhaps this is because the fixups are being */
		/* done in a far from optimal order ? What this means is that often */
		/* for big movements reseeding will be faster. To get the best of both, */
		/* we try and estimate when the fixup method will break even with */
		/* reseeding, and switch over. */

		/* Estimate how many vertexes will be hit by the move */
		nvxhits = ofps_quick_check_hits(s);

		/* Decide which way to go */
		thresh = 1.0/(di * di);
		hratio = nvxhits/(double)s->nv;
//printf("~1 quick check of vertex hits = %d, ratio %f, threshold %f\n",nvxhits,hratio,thresh);

		/* Hmm. Re-seed seems to sometimes be slower than expected for > 3D, */
		/* so don't use it. */
		if (hratio < thresh && di < 4) {
			doinc = 1;
		}

#ifdef FORCE_RESEED		/* Force reseed after itteration */
		doinc = 0;
#else
# ifdef FORCE_INCREMENTAL	/* Force incremental update after itteration */
		doinc = 1;
# endif
#endif
		/* Incrementally update veronoi */
		if (doinc) {

			if (s->verb)
				printf("Fixing up veronoi\n");

			/* Re-position the vertexes, and fixup the veronoi */
			ofps_repos_and_fix_voronoi(s);

		/* Reseed the veronoi */
		} else {

			if (s->verb)
				printf("Re-seeding\n");

			/* remove nodes from the spatial acceleration grid. */
			for (i = 0; i < s->tinp; i++) {
				node *pp = s->n[i];
				ofps_rem_nacc(s, pp);			/* Remove from spatial accelleration grid */
			}
	
			/* And recompute veronoi, and add to spatial accelleration grid. */
			ofps_redo_voronoi(s);
		}
		ofps_re_create_node_node_vtx_lists(s);
		ofps_create_mids(s);

		ofps_stats(s);
		if (s->verb) {
			printf("It %d: Maxmv = %f, MinPoint = %.3f, Min = %.3f, Avg. = %.3f, Max = %.3f, %.1f secs.\n",s->optit+1,sqrt(s->mxmvsq),s->smns,s->mn,s->av,s->mx,(msec_time() - s->l_mstime) / 1000.0);
#ifdef STATS
			printf("Current vtx %d, created %d, deleted %d, positioned %d\n", s->nv,s->nvtxcreated - s->l_nvtxcreated,s->nvtxdeleted - s->l_nvtxdeleted, s->positions - s->l_positions);
			s->l_positions = s->positions;
			s->l_nvtxcreated = s->nvtxcreated;
			s->l_nvtxdeleted = s->nvtxdeleted;
#endif
			s->l_mstime = msec_time();
		}

#ifdef DUMP_STRUCTURE
		dump_node_vtxs(s, 1);
//		{ char buf[200]; sprintf(buf, "After itteration %d",s->optit+1); dump_node_vtxs2(s, buf); }
		printf("=========================================================================\n");
#endif
#ifdef DUMP_PLOT
		dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 1, -1);		/* Device, wait, verticies */
#endif /* DUMP_PLOT */

#ifdef SANITY_RESEED_AFTER_FIXUPS
		/* For debugging, replace the incremental fixed up veronoi with */
		/* a from scratch one.  */

		if (s->verb)
			printf("Re-seeding after fixup:\n");

		/* Save the current incremental vertexes */
		save_ivertexes(s);

		ofps_redo_voronoi(s);
		ofps_re_create_node_node_vtx_lists(s);
		ofps_create_mids(s);

		ofps_stats(s);
		if (s->verb) {
			printf("It %d: Maxmv = %f, MinPoint = %.3f, Min = %.3f, Avg. = %.3f, Max = %.3f, %.1f secs.\n",s->optit+1,sqrt(s->mxmvsq),s->smns,s->mn,s->av,s->mx,(msec_time() - s->l_mstime) / 1000.0);
#ifdef STATS
			printf("Current vtx %d, created %d, deleted %d, positioned %d\n", s->nvtxcreated - s->l_nvtxcreated,s->nvtxdeleted - s->l_nvtxdeleted, s->positions - s->l_positions);
			s->l_positions = s->positions;
			s->l_nvtxcreated = s->nvtxcreated;
			s->l_nvtxdeleted = s->nvtxdeleted;
#endif
			s->l_mstime = msec_time();
		}

		/* Check that no node other than a parent is closer to any vertex */
		if (check_vertex_closest_node(s)) {
			warning("Verify that re-seed leaves only parents closest to vertexes failed");
		}

		/* Check the incremental vertexes against the re-seeded vertexes */
		if (check_vertexes(s)) {
			warning("Verify of incremental vertexes failed!");
			printf("Verify of incremental vertexes failed!\n");
		} else {
			warning("Verify of incremental vertexes suceeded!");
		}
#ifdef DUMP_STRUCTURE
		dump_node_vtxs(s, 1);
#endif
#ifdef DUMP_PLOT
		dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 1, -1);		/* Device, wait, verticies */
#endif /* DUMP_PLOT */
#endif /* SANITY_RESEED_AFTER_FIXUPS */

		if (sqrt(s->mxmvsq) < stoptol)
			break;
	}
}

/* ------------------------------------------------------------------------ */
/* Main object creation/destruction */

/* Destroy ourselves */
static void
ofps_del(ofps *s) {
	int i, e, di = s->di;

	if (s->ufx != NULL)
		free(s->ufx);

	/* Free our nodes */
	for (i = 0; i < s->np; i++) {
		node_free(s, s->n[i]);
	} 
	s->n -= s->gnp;		/* Fixup offset */
	free(s->n);
	free(s->_n);

	/* Any free vertexes */
	while (s->fvtx != NULL) {
		vtx *p = s->fvtx;
		s->fvtx = p->link;
		free(p);
	}

	/* Any other allocations */
	s->sob->del(s->sob);
	if (s->combs != NULL) {
		for (i = 0; i < s->_ncombs; i++) {
			if (s->combs[i].v1 != NULL)
				free(s->combs[i].v1);
			if (s->combs[i].v2 != NULL)
				free(s->combs[i].v2);
		}
		free(s->combs);
	}
	if (s->sc)
		free(s->sc);

	if (s->svtxs != NULL)
		free(s->svtxs);

	if (s->_grid != NULL)
		free(s->_grid);

	if (s->acnl != NULL)
		free(s->acnl);

	if (s->vtreep != NULL)
		aat_adelete(s->vtreep);

	for (e = 0; e <= (di+1); e++) {
		if (s->vtrees[e] != NULL)
			aat_adelete(s->vtrees[e]);
	}

	if (s->pcache != NULL)
		s->pcache->del(s->pcache);

	free(s);
}

/* Constructor */
ofps *new_ofps(
int verb,				/* Verbosity level, 1 = progress, 2 = warnings */
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
int tinp,				/* Total number of points to generate, including fixed */
int good,				/* 0 = fast, 1 = good */
double dadaptation,		/* Degree of adaptation to device characteristic 0.0 - 1.0 */
double devd_wght,		/* Device space weighting (if dad < 0) */
double perc_wght,		/* Perceptual space weighting (if dad < 0) */
double curv_wght,		/* Curvature weighting (if dad < 0) */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
) {
	return new_ofps_ex(verb, di, ilimit, NULL, NULL, tinp, good,
	                   dadaptation, devd_wght, perc_wght, curv_wght,
	                   fxlist, fxno, percept, od, 0, -1);
}

/* Extended constructor */
ofps *new_ofps_ex(
int verb,				/* Verbosity level, 1 = progress, 2 = warnings */
int di,					/* Dimensionality of device space */
double ilimit,			/* Total ink limit (sum of device coords max) */
double *imin,			/* Ink limit - limit on min of p[], usually >= 0.0 (may be NULL) */
double *imax,			/* Ink limit - limit on min of p[], usually <= 1.0 (may be NULL) */
int tinp,				/* Total number of points to generate, including fixed */
int good,				/* 0 = fast, 1 = good */
double dadaptation,		/* Degree of adaptation to device characteristic 0.0 - 1.0 */
double devd_wght,		/* Device space weighting (if dad < 0) */
double perc_wght,		/* Perceptual space weighting (if dad < 0) */
double curv_wght,		/* Curvature weighting (if dad < 0) */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od,				/* context for Perceptual function */
int ntostop,			/* Debug - number of points until diagnostic stop */
int nopstop				/* Debug - number of optimizations until diagnostic stop, -1 = not */
) {
	int i, e;
	ofps *s;
	long stime,ttime;

	stime = clock();

	if ((s = (ofps *)calloc(sizeof(ofps), 1)) == NULL)
		error ("ofps: malloc failed on new ofps");

	if (di > MXPD)
		error ("ofps: Can't handle di %d",di);

	s->verb = verb;
	s->ntostop = ntostop;
	s->nopstop = nopstop;

	if ((s->sob = new_sobol(di)) == NULL)
		error ("ofps: new_sobol %d failed", di);
	
	if (s->verb)
		printf("Degree of adaptation: %.3f\n", dadaptation);

	/* Set internal values explicitly */
	if (dadaptation < 0.0) {
		s->devd_wght = devd_wght;
		s->perc_wght = perc_wght;
		s->curv_wght = curv_wght;

	/* Set values implicitly with adapation level */
	} else {
		if (dadaptation > 1.0)
			dadaptation = 1.0;

		/* Convert to internal numbers */
		s->perc_wght = ADAPT_PERCWGHT * dadaptation;
		s->curv_wght = ADAPT_CURVWGHT * dadaptation * dadaptation;
		s->devd_wght = 1.0 - s->perc_wght;
	}
	if (s->verb)
		printf("Adaptation weights: Device = %.3f, Perceptual = %.3f, Curvature = %.3f\n",
		                                             s->devd_wght,s->perc_wght,s->curv_wght);

	s->di = di;

	if (tinp < fxno)	/* Make sure we return at least the fixed points */
		tinp = fxno;

	s->fxno = fxno;		/* Number of fixed points provided */
	s->tinp = tinp;		/* Target total number of points */

	/* Hack to workaround pathalogical case. At ilimit == di-2.0, we get > 32 bits */
	/* of mask for CMYK */
	if (di >= 3
	 && ilimit >= (di-2.0 - 2 * ILIMITEPS)
	 && ilimit <= (di-2.0 + 2 * ILIMITEPS))
		ilimit = di-2.0 - 2 * ILIMITEPS;

	/* Hack to workaround pathalogical case. At ilimit == 100% we get a failure */
	/* to add any variable steps */
	if (ilimit > 0.9999 && ilimit < 1.0001)
		ilimit = 0.9999;

	s->ilimit = ilimit;

	for (e = 0; e < di; e++) {
		if (imin != NULL)
			s->imin[e] = imin[e];
		else
			s->imin[e] = 0.0;
			
		if (imax != NULL)
			s->imax[e] = imax[e];
		else
			s->imax[e] = 1.0;
	}

	/* Compute an approximate half expected sample point spacing, */
	/* and setup seeding acceleration grid. */
	{
		double vol = 1.0;
		double eprange;

		for (e = 0; e < di; e++)
			vol *= s->imax[e] - s->imin[e];

		vol /= tinp;				/* Approx vol per point */
		vol = pow(vol, 1.0/di);		/* Distance per point */

		s->surftol = SURFTOL * vol;
//printf("~1 surftol = %f\n",s->surftol);
	}

#ifdef STANDALONE_TEST
	/* If no perceptual function given, use default */
	if (percept == NULL) {
		s->percept = default_ofps_to_percept;
		s->od = s;
	} else {
		s->percept = percept;
		s->od = od;
	}
#else
	/* If no perceptual function given, use default */
//warning("~1 new_ofps_ex() forcing default perceptual function");
	if (percept == NULL) {
		s->percept = default_ofps_to_percept;
		s->od = s;
	} else {
		s->percept = percept;
		s->od = od;
	}
#endif

	s->good = good;		/* Fast/Good flag */
	s->lperterb = PERTERB_AMOUNT;
#ifdef OPT_MAXITS_2
	if (s->good == 0)
		s->lperterb = PERTERB_AMOUNT_2;
#endif
	s->ssurfpref = INITIAL_SURFACE_PREF;
	s->esurfpref = FINAL_SURFACE_PREF;

	/* Init method pointers */
	s->reset = ofps_reset;
	s->read  = ofps_read;
	s->stats = ofps_stats;
	s->del   = ofps_del;

	s->gnp = 2 * di + 1 + 2;	/* Gamut boundary + inside/outside fake points */
								/* -1 to -2di-1 are fake boundary nodes indexes, */
								/* with -2di-1 being the ink limit boundary. */
								/* -2di-2 is the fake inside node. */
								/* -2di-3 is the fake outside node. */

	/* Allocate the space for the target number of points */
	if ((s->_n = (node *)calloc(sizeof(node), s->gnp + s->tinp)) == NULL)
		error ("ofps: malloc failed on sample nodes");
	if ((s->n = (node **)calloc(sizeof(node *), s->gnp + s->tinp)) == NULL)
		error ("ofps: malloc failed on sample nodes");
	s->n += s->gnp;		/* Allow -ve index for fake points */
	for (i = -s->gnp; i < s->tinp; i++) {
		int bitp;
		s->n[i] = &s->_n[i + s->gnp];
		s->n[i]->ix = i;

		bitp = 31 & (i + (i >> 4) + (i >> 8) + (i >> 12));
		s->n[i]->ixm = (1 << bitp);
	}

	s->np = s->fnp = 0;

#ifdef STATS
	/* Save current counts to report stats after a pass */
	s->l_positions = s->positions;
	s->l_nvtxcreated = s->nvtxcreated;
	s->l_nvtxdeleted = s->nvtxdeleted;
#endif
	s->l_mstime = msec_time();

	/* Setup the eperr sorted trees */
	if ((s->vtreep = aat_anew(vtx_aat_cmp_eperr)) == NULL)
		error("Allocating aat tree failed");

	/* One sorted tree per number of surface planes */
	for (e = 0; e <= (di+1); e++) {
		if ((s->vtrees[e] = aat_anew(vtx_aat_cmp_eserr)) == NULL)
			error("Allocating aat tree failed");
	}

#ifdef CACHE_PERCEPTUAL
	ofps_init_pcache(s);
# endif /* CACHE_PERCEPTUAL */
	
	/* Setup spatial acceleration grid */
	ofps_init_acc1(s);

	/* Initialse the empty veronoi etc. */ 
	ofps_binit(s);

	/* Setup spatial acceleration grid (2) */
	ofps_init_acc2(s);

	/* Setup the fixed points */
	ofps_setup_fixed(s, fxlist, fxno);

	if (fxno > 0 && tinp <= fxno) {		/* There are no moveable points to create */

		/* Add the fixed points */
		if (ofps_add_fixed(s)) {
			s->del(s);
			return NULL;
		}

		if (s->verb && fxno > 0) {
			ofps_stats(s);
			printf("After fixed points: MinPoint = %.3f, Min = %.3f, Avg. = %.3f, Max = %.3f\n",s->smns,s->mn,s->av,s->mx);
		}
	}

	if (tinp > fxno) {		/* There are movable points to create */

		/* Add the fixed points and create the moveable points */
		ofps_seed(s);
		ofps_re_create_node_node_vtx_lists(s);
		ofps_create_mids(s);

		ofps_stats(s);
		if (s->verb) {
			printf("After seeding points: MinPoint = %.3f, Min = %.3f, Avg. = %.3f, Max = %.3f, %.1f secs\n",s->smns,s->mn,s->av,s->mx,(msec_time() - s->l_mstime) / 1000.0);
	
#ifdef STATS
			printf("Current vtx %d, created %d, deleted %d, positioned %d\n", s->nv,s->nvtxcreated - s->l_nvtxcreated,s->nvtxdeleted - s->l_nvtxdeleted, s->positions - s->l_positions);
			s->l_positions = s->positions;
			s->l_nvtxcreated = s->nvtxcreated;
			s->l_nvtxdeleted = s->nvtxdeleted;
#endif
			s->l_mstime = msec_time();
		}
# ifdef DUMP_STRUCTURE
		printf("After seeding:\n");
		dump_node_vtxs(s, 1);
//		dump_node_vtxs2(s, "After seeding");
#else /* !DUMP_STRUCTURE */
#ifdef SANITY_CHECK_CONSISTENCY
		sanity_check(s, 1);
#endif
#endif /* !DUMP_STRUCTURE */
#ifdef DUMP_PLOT
		dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 1, -1);		/* Device, No wait, no verticies */
#endif /* DUMP_PLOT */

#ifdef DOOPT
		/* Do the optimization */
		ofps_optimize(s);
#endif /* DOOPT */

# ifdef DUMP_STRUCTURE
		printf("After optimization:\n");
		dump_node_vtxs(s, 1);
//		dump_node_vtxs2(s, "After optimization");
#else /* !DUMP_STRUCTURE */
#ifdef SANITY_CHECK_CONSISTENCY
		sanity_check(s, 1);
#endif
#endif /* !DUMP_STRUCTURE */
		ofps_stats(s);
		if (s->verb)
			printf("After optimization: MinPoint = %.3f, Min = %.3f, Avg. = %.3f, Max = %.3f\n",s->smns, s->mn,s->av,s->mx);
#ifdef DUMP_PLOT
		dump_image(s, PERC_PLOT, DO_WAIT, DUMP_VTX, DUMP_PLA, 1, -1);		/* Device, wait, verticies */
#endif /* DUMP_PLOT */
	}

	ofps_reset(s);		/* Reset read index */

#if defined(DEBUG) || defined(STATS)
	{
		vtx *vx;
		int novtx = 0;
		int totvtxverts = 0; 
		int maxvtxverts = 0;
		vtx **svtxs;		/* Sorted vertexes by number of vertexes */

		ttime = clock() - stime;
		printf("Execution time = %f seconds\n",ttime/(double)CLOCKS_PER_SEC);

		/* Look at the vertexes */
		for (novtx = 0, vx = s->uvtx; vx != NULL; vx = vx->link, novtx++)
			;

		if ((svtxs = (vtx **)malloc(sizeof(vtx *) * novtx)) == NULL)
			error ("ofps: malloc failed on vertex pointer list");

		/* Look at the vertexes */
		for (novtx = 0, vx = s->uvtx; vx != NULL; vx = vx->link, novtx++) {

			svtxs[novtx] = vx;

			totvtxverts += vx->nnv; 
			if (vx->nnv > maxvtxverts)
				maxvtxverts = vx->nnv;
		}

#define HEAP_COMPARE(A,B) ((A)->nnv > (B)->nnv)
		HEAPSORT(vtx *, svtxs, novtx);
#undef HEAP_COMPARE

//		printf("Top 20 vertexes per vertex:\n");
//		for (i = 0; i < 20 && i < novtx; i++) {
//			printf("  Vtx no %d, no vtxs = %d\n",svtxs[i]->no,svtxs[i]->nnv);
//		}

		fprintf(stderr,"Average vertexes per vertex %.1f, max %d\n",totvtxverts/(double)novtx,maxvtxverts);
		fprintf(stderr,"Average hit vertexes per add %.1f\n",s->nhitv/(double)s->nsurfadds,s->maxhitv);
		fprintf(stderr,"Total number of vertex = %d\n",novtx); 
		fprintf(stderr,"Total vertex positions = %d\n",s->positions); 
		fprintf(stderr,"Total dnsqs = %d\n",s->dnsqs); 
		fprintf(stderr,"Total function calls = %d\n",s->funccount); 
		fprintf(stderr,"Average dnsqs/position = %.2f\n",s->dnsqs/(double)s->positions); 
		fprintf(stderr,"Average function calls/dnsq = %.1f\n",s->funccount/(double)s->dnsqs); 
		fprintf(stderr,"Maximum function calls/dnsq = %d\n",s->maxfunc); 
		fprintf(stderr,"Average function calls/sucessful dnsq = %.2f\n",s->sucfunc/(double)s->sucdnsq); 
		fprintf(stderr,"Average function calls/position = %.1f\n",s->funccount/(double)s->positions); 
		fprintf(stderr,"Maximum tries for dnsq sucess %d\n",s->maxretries); 
		fprintf(stderr,"Number of position_vtx failures %d\n",s->posfails); 
		fprintf(stderr,"Vertex hit check efficiency = %.1f%%\n",100.0 * (1.0 - s->vvchecks/(double)s->vvpchecks));
		fprintf(stderr,"Average accell cells searched = %.2f\n",s->ncellssch/(double)s->naccsrch);
		fprintf(stderr,"add_to_vsurf hit rate = %.1f%%\n",100.0 *  s->add_hit/(s->add_hit + s->add_mis));
#ifdef  DOOPT
		fprintf(stderr,"fixup add_to_vsurf hit rate = %.1f%%\n",100.0 *  s->fadd_hit/(s->fadd_hit + s->fadd_mis));
		fprintf(stderr,"Vertex closest search efficiency = %.1f%%\n",100.0 * (1.0 - s->nvschd/(double)s->nvfschd));
		fprintf(stderr,"Node closest search efficiency = %.1f%%\n",100.0 * (1.0 - s->nnschd/(double)s->nnfschd));
#endif

		free(svtxs);
	}
#endif

	return s;
}

/* =================================================== */

#ifdef STANDALONE_TEST

/* Graphics Gems curve */
static double gcurve(double vv, double g) {
	if (g >= 0.0) {
		vv = vv/(g - g * vv + 1.0);
	} else {
		vv = (vv - g * vv)/(1.0 - g * vv);
	}
	return vv;
}


static void sa_percept(void *od, double *p, double *d) {
	double dd[2];

	/* Default linear */
	p[0] = 100.0 * (dd[0] = d[0]);
	p[1] = 100.0 * (dd[1] = d[1]);

	/* Normal non-linear test */
//	p[0] = 100.0 * gcurve(dd[0], -8.0);
//	p[1] = 100.0 * gcurve(dd[1], 4.0);

	/* More extreme non-linear test */
	p[0] = 100.0 * gcurve(dd[0], -16.0);
	p[1] = 100.0 * gcurve(dd[1], 8.0);

	/* An X break point to test curvature weighting */
//	if (dd[0] < 0.5)
//		p[0] = 100.0 * 0.6 * dd[0];
//	else
//		p[0] = 100.0 * (0.3 + 1.4 * (dd[0] - 0.5));
//	p[1] = 100.0 * dd[1];

//	if (dd[0] < 0.0)
//		dd[0] = 0.0;
//	if (dd[1] < 0.0)
//		dd[1] = 0.0;
//	p[0] = 100.0 * pow(dd[0], 0.5);
//	p[1] = 100.0 * pow(dd[1], 1.0);
//	p[1] = 0.8 * p[1] + 0.2 * p[0];

	/* One that causes dnsq failures due to ACCELL failure */
//	p[0] = gcurve(dd[0], -4.0);
//	p[1] = gcurve(dd[1], 2.0);
//	p[0] = 100.0 * gcurve(0.6 * p[0] + 0.4 * p[1], 2.0);
//	p[1] = 100.0 * gcurve(0.1 * p[1] + 0.9 * p[1], -4.0);

//	p[0] = 100.0 * dd[0] * dd[0];
//	p[1] = 100.0 * dd[1] * dd[1];
}

int
main(argc,argv)
int argc;
char *argv[];
{
	int npoints = 55;
	int ntostop = 0;
	int nopstop = 0;
	ofps *s;
	fxpos fx[4];		/* Any fixed points */
	int nfx = 0;

	error_program = argv[0];

	printf("Standalone test of ofps, args are: no. of points, default %d, points to skip before diag. plots, optim passes to skip\n",npoints);

	if (argc > 1)
		npoints = atoi(argv[1]);

	if (argc > 2)
		ntostop = atoi(argv[2]);

	if (argc > 3)
		nopstop = atoi(argv[3]);

	fx[0].p[0] = 0.5;
	fx[0].p[1] = 0.5;

	fx[1].p[0] = 0.145722;
	fx[1].p[1] = 0.0;

	fx[2].p[0] = 1.0;
	fx[2].p[1] = 0.104414;

	nfx = 0;

	/* Create the required points */
	s = new_ofps_ex(1, 2, 1.5, NULL, NULL, npoints, 1,
//	s = new_ofps_ex(1, 2, 2.5, NULL, NULL, npoints, 1,
	             SA_ADAPT, SA_DEVD_MULT, SA_PERC_MULT, SA_INTERP_MULT,
	             fx, nfx, sa_percept, (void *)NULL, ntostop, nopstop);

#ifdef DUMP_PLOT
	printf("Device plot (with verts):\n");
	dump_image(s, 0, DO_WAIT, 1, DUMP_PLA, 1, -1);
	printf("Device plot:\n");
	dump_image(s, 0, DO_WAIT, 0, 0, 1, -1);
	printf("Perceptual plot (with verts):\n");
	dump_image(s, 1, DO_WAIT, 1, DUMP_PLA, 1, -1);
	printf("Perceptual plot:\n");
	dump_image(s, 1, DO_WAIT, 0, 0, 1, -1);
#endif /* DUMP_PLOT */

	s->del(s);

	return 0;
}

#endif /* STANDALONE_TEST */

#define WIDTH 400			/* Raster size for debug plots */
#define HEIGHT 400

/* Utility - return a string containing the di coord */
static char *pco(int di, int *co) {
	static char buf[5][200];
	static int ix = 0;
	int e;
	char *bp;

	if (++ix >= 5)
		ix = 0;
	bp = buf[ix];

	for (e = 0; e < di; e++) {
		if (e > 0)
			*bp++ = ' ';
		sprintf(bp, "%d", co[e]); bp += strlen(bp);
	}
	return buf[ix];
}

/* Utility - return a string containing the di vector */
static char *ppos(int di, double *p) {
	static char buf[5][200];
	static int ix = 0;
	int e;
	char *bp;

	if (++ix >= 5)
		ix = 0;
	bp = buf[ix];

	for (e = 0; e < di; e++) {
		double val = p[e];
		/* Make -0.00000000 turn into 0.000 for cosmetics */
		if (val < 0.0 && val >-1e-9)
			val = 0.0;
		if (e > 0)
			*bp++ = ' ';
		sprintf(bp, "%f", val); bp += strlen(bp);
	}
	return buf[ix];
}

/* Utility - return a string containing the di+1 combination */
static char *pcomb(int di, int *n) {
	static char buf[5][200];
	static int ix = 0;
	int e;
	char *bp;

	if (++ix >= 5)
		ix = 0;
	bp = buf[ix];

	for (e = 0; e <= di; e++) {
		if (e > 0)
			*bp++ = ' ';
		sprintf(bp, "%d", n[e]); bp += strlen(bp);
	}
	return buf[ix];
}

/* Utility - return a string containing the eperr/eserr value */
static char *peperr(double eperr) {
	static char buf[5][200];
	static int ix = 0;
	int e;
	char *bp;

	if (++ix >= 5)
		ix = 0;
	bp = buf[ix];

	if (eperr >= 1e50)
		sprintf(bp,"%s", "Big");
	else
		sprintf(bp,"%f",eperr);
	return buf[ix];
}

/* --------------------------------------------------------------- */
#if defined(DEBUG) || defined(DUMP_PLOT_SEED) || defined(DUMP_PLOT)

/* Dump the current point positions to a plot window file */
static void
dump_image(
	ofps *s,
	int pcp,			/* Do perceptual plot */
	int dwt,			/* Do wait for a key */
	int dvx,			/* Dump voronoi verticies and mid points */
	int dpla,			/* Dump node planes */ 
	int ferr,			/* Show final error rather than seeding error */
	int noi				/* -1 for general state, node of interest for particular */
) {
	int i, j, k, e, di = s->di;
	double minx, miny, maxx, maxy;
	static double *x1a = NULL;		/* Previous sample locations */
	static double *y1a = NULL;
	static double *x2a = NULL;		/* Current sample locations */
	static double *y2a = NULL;
	static char *_ntext, **ntext;
	static int _n3 = 0;				/* Current Voronoi verticies */
	static double *x3a = NULL;
	static double *y3a = NULL;
	static plot_col *mcols = NULL;
	static char *_mtext, **mtext;
	int n3;
	static double *x4a = NULL;		/* plane vectors */
	static double *y4a = NULL;
	static double *x5a = NULL;
	static double *y5a = NULL;
	static plot_col *ocols = NULL;
	static int _o4 = 0;
	int o4;

	if (pcp != 0) {	/* Perceptual range */
		vtx *vx;
		minx = miny = 1e60;
		maxx = maxy = -1e60;
		for (vx = s->uvtx; vx != NULL; vx = vx->link) {
			double v[MXPD];

			if (vx->v[0] < minx)
				minx = vx->v[0];
			if (vx->v[1] < miny)
				miny = vx->v[1];
			if (vx->v[0] > maxx)
				maxx = vx->v[0];
			if (vx->v[1] > maxy)
				maxy = vx->v[1];
		}
	} else {
		minx = 0.0;	/* Assume */
		miny = 0.0;
		maxx = 1.0;
		maxy = 1.0;
	}

#ifdef NEVER
	/* Expand the range a little */
	minx -= 0.1 * (maxx - minx);
	maxx += 0.1/1.1 * (maxx - minx);
	miny -= 0.1 * (maxy - miny);
	maxy += 0.1/1.1 * (maxy - miny);
#endif
	
	if (x1a == NULL) {
		if ((x1a = (double *)malloc(s->tinp * sizeof(double))) == NULL)
			error ("ofps: malloc failed x1a");
		if ((y1a = (double *)malloc(s->tinp * sizeof(double))) == NULL)
			error ("ofps: malloc failed ya1");
		if ((x2a = (double *)malloc(s->tinp * sizeof(double))) == NULL)
			error ("ofps: malloc failed x2a");
		if ((y2a = (double *)malloc(s->tinp * sizeof(double))) == NULL)
			error ("ofps: malloc failed y2a");
		if ((_ntext = (char *)malloc(s->tinp * 10 * sizeof(char))) == NULL)
			error ("ofps: malloc failed _ntext");
		if ((ntext = (char **)malloc(s->tinp * sizeof(char *))) == NULL)
			error ("ofps: malloc failed ntext");
		for (i = 0; i < s->tinp; i++)
			ntext[i] = _ntext + i * 10;
	}

	/* Add sample node location */
	for (i = 0; i < s->np; i++) {
		node *p = s->n[i];
		
		if (pcp != 0) {
			double ov[MXPD];
			ofps_cc_percept(s, ov, p->op);
			x1a[i] = ov[0];
			y1a[i] = ov[1];
			x2a[i] = p->v[0];
			y2a[i] = p->v[1];
		} else {
			x1a[i] = p->op[0];
			y1a[i] = p->op[1];
			x2a[i] = p->p[0];
			y2a[i] = p->p[1];
		}
		sprintf(ntext[i],"%d",p->ix);
//		sprintf(ntext[i],"",p->ix);
	}

	if (dvx) {
		vtx *vx;
		mid *mp;
		node *p = NULL;
//		double rgb0[3] = { 0.0, 0.5, 0.5 };		/* "cool" */
//		double rgb1[3] = { 1.0, 0.5, 0.0 };		/* "warm" */
		double rgb0[3] = { 0.0, 1.0, 0.0 };		/* "cool" */
		double rgb1[3] = { 1.0, 0.0, 0.5 };		/* "warm" */
		double mine, maxe;				/* Min and max vertex eserr */

		if (noi >= 0)
			p = s->n[noi];

		if (x3a == NULL) {		/* Initial allocation */
			_n3 = s->np * 4;
			if ((x3a = (double *)malloc(_n3 * sizeof(double))) == NULL)
				error ("ofps: malloc failed x3a");
			if ((y3a = (double *)malloc(_n3 * sizeof(double))) == NULL)
				error ("ofps: malloc failed y3a");
			if ((mcols = (plot_col *)malloc(_n3 * sizeof(plot_col))) == NULL)
				error ("ofps: malloc failed mcols");
			if ((_mtext = (char *)malloc(_n3 * 10 * sizeof(char))) == NULL)
				error ("ofps: malloc failed _mtext");
			if ((mtext = (char **)malloc(_n3 * sizeof(char *))) == NULL)
				error ("ofps: malloc failed mtext");
			for (i = 0; i < _n3; i++)
				mtext[i] = _mtext + i * 10;
		}

		/* Compute min & max serr for each vertex */
		mine = 1e6;
		maxe = -1e6;
		for (vx = s->uvtx; vx != NULL; vx = vx->link) {
			if (vx->ghost)
				continue;
			if (vx->eserr > maxe)
				maxe = vx->eserr;
			if (vx->eserr > NUMTOL && vx->eserr < mine)
				mine = vx->eserr;
		}
		if ((maxe - mine) < 10.0)
			maxe = mine + 1.0;

		/* Add mid points */
		for (n3 = 0, mp = s->umid; mp != NULL; mp = mp->link, n3++) {

			if (n3 >= _n3) {		/* need more space */
				_n3 = 2 * _n3 + 5;
				if ((x3a = (double *)realloc(x3a, _n3 * sizeof(double))) == NULL)
					error ("ofps: realloc failed x3a %d",_n3);
				if ((y3a = (double *)realloc(y3a, _n3 * sizeof(double))) == NULL)
					error ("ofps: realloc failed y3a");
				if ((mcols = (plot_col *)realloc(mcols, _n3 * sizeof(plot_col))) == NULL)
					error ("ofps: realloc failed mcols");
				if ((_mtext = (char *)realloc(_mtext, _n3 * 10 * sizeof(char))) == NULL)
					error ("ofps: realloc failed _mtext");
				if ((mtext = (char **)realloc(mtext, _n3 * sizeof(char *))) == NULL)
					error ("ofps: realloc failed mtest");
				for (i = 0; i < _n3; i++)
					mtext[i] = _mtext + i * 10;
			}
			if (pcp != 0) {
				x3a[n3] = mp->v[0];
				y3a[n3] = mp->v[1];
			} else {
				x3a[n3] = mp->p[0];
				y3a[n3] = mp->p[1];
			}

			/* Show mid points in grey */
			mcols[n3].rgb[0] = 0.85;
			mcols[n3].rgb[1] = 0.85;
			mcols[n3].rgb[2] = 0.85;

			sprintf(mtext[n3],"%s","");
			sprintf(mtext[n3],"%d",mp->no);
//			sprintf(mtext[n3],"%d",(int)(mp->eserr + 0.5));
		}

		/* Add Voronoi verticies */
		for (vx = s->uvtx; vx != NULL; vx = vx->link, n3++) {

			if (n3 >= _n3) {		/* need more space */
				_n3 = _n3 * 2 + 5;
				if ((x3a = (double *)realloc(x3a, _n3 * sizeof(double))) == NULL)
					error ("ofps: realloc failed x3a %d",_n3);
				if ((y3a = (double *)realloc(y3a, _n3 * sizeof(double))) == NULL)
					error ("ofps: realloc failed y3a");
				if ((mcols = (plot_col *)realloc(mcols, _n3 * sizeof(plot_col))) == NULL)
					error ("ofps: realloc failed mcols");
				if ((_mtext = (char *)realloc(_mtext, _n3 * 10 * sizeof(char))) == NULL)
					error ("ofps: realloc failed _mtext");
				if ((mtext = (char **)realloc(mtext, _n3 * sizeof(char *))) == NULL)
					error ("ofps: realloc failed mtext");
				for (i = 0; i < _n3; i++)
					mtext[i] = _mtext + i * 10;
			}
			if (pcp != 0) {
				x3a[n3] = vx->v[0];
				y3a[n3] = vx->v[1];
			} else {
				x3a[n3] = vx->p[0];
				y3a[n3] = vx->p[1];
			}

			/* Show the vertexes as warm to cold, depending on their eserr */
			if (p == NULL) {
				double bf;

				bf = (vx->eserr - mine)/(maxe - mine);
				if (bf < 0.0)
					bf = 0.0;
				if (bf > 1.0)
					bf = 1.0;
				
				for (e = 0; e < 3; e++)
					mcols[n3].rgb[e] = bf * rgb1[e] + (1.0 - bf) * rgb0[e];

//printf("~1 serr = %f, color = %f %f %f\n",vx->eserr, mcols[n3].rgb[0], mcols[n3].rgb[1], mcols[n3].rgb[2]);
//				sprintf(mtext[n3],"");
//				sprintf(mtext[n3],"%d",(int)(vx->eserr + 0.5));

#ifndef NEVER /* Vertex no */
				sprintf(mtext[n3],"%d",vx->no);
#endif

#ifdef NEVER /* Vertex no and eserr */
				if (vx->eserr >= 1e50)
					sprintf(mtext[n3],"%d:Big",vx->no);
				else
					sprintf(mtext[n3],"%d:%d",vx->no,(int)(vx->eserr + 0.5));
#endif

#ifdef NEVER /* eserr */
				if (vx->eserr >= 1e50)
					sprintf(mtext[n3],"Big");
				else
					sprintf(mtext[n3],"%d",(int)(vx->eserr + 0.5));
#endif

			/* Highlight the vertcies of interest */
			} else {
				for (j = 0; j < p->nvv; j++) {
					if (p->vv[j] == vx)
						break;
				}
				if (j < p->nvv) {			/* Vertex associated with node of interest */
					mcols[n3].rgb[0] = 0.1;
					mcols[n3].rgb[1] = 0.9;
					mcols[n3].rgb[2] = 0.9;

					sprintf(mtext[n3],"%d",(int)(vx->eserr + 0.5));
			
				} else {
					mcols[n3].rgb[0] = 0.82;		/* default color */
					mcols[n3].rgb[1] = 0.59;
					mcols[n3].rgb[2] = 0.0;
	
					sprintf(mtext[n3],"%s","");
				}
			}
		}
#ifdef DUMP_EPERR		/* Create .tiff of eperr */
		if (s->np >= s->ntostop) {

			unsigned char pa[WIDTH * 3];
			char *name = "ofps.tif";
			int width = WIDTH;
			int height = HEIGHT;
			int x, y;
			TIFF *tif;
			double pos[MXPD], vpos[MXPD];
			double rgb_low[3] = { 0.0, 1.0, 0.0 };		/* "low error" */
			double rgb_high[3] = { 1.0, 0.0, 0.0 };		/* "high error" */

			if ((tif = TIFFOpen(name, "w")) == NULL) {
				fprintf(stderr,"Failed to open output TIFF file '%s'\n",name);
				exit (-1);
			}

			TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,  width);
			TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
			TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
			TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
			TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
			TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
			TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
			TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

			mine = 0.0;

			for (y = 0; y < height; y++) {
				pos[1] = 1.0 - y/(height-1.0);

				/* Fill in pa[] with colors for this line */
				for (x = 0; x < width; x++) {
					double ss;
					unsigned char *dp;
					double bf;
					double beserr, eserr;

					dp = pa + x * 3;
					pos[0] = x/(width-1.0);
					dp[0] = dp[1] = dp[2] = 0;
//printf("~1 doing %d %d pos %f %f\n",x,y,pos[0],pos[1]);

					/* Lookup perceptual value at sample point location */
					ofps_cc_percept(s, vpos, pos);
	
					/* See if the sample is in gamut */
					for (ss = 0.0, e = 0; e < s->di; e++) {
						if (pos[e] < s->imin[e]
						 || pos[e] > s->imax[e])
							break;
						ss += pos[e];
					}
					if (e < s->di || ss > (s->ilimit + ILIMITEPS)) {
//printf("~1 out of gamut\n");
						continue;
					}

					/* We determine the eserr by evaluating eserr for */
					/* every node, and keeping the smallest. */
					/* (This could be speeded up by using nearest search function) */
					beserr = 1e80;
					for (i = 0; i < s->np; i++) {
						node *np = s->n[i];
					
						eserr = ofps_comp_eperr9(s, NULL, vpos, pos, np->v, np->p, np->nsp);
						if (eserr < beserr)
							beserr = eserr;
					}
					bf = (beserr - mine)/(maxe - mine);
//printf("~1 beserr = %f, bf = %f\n",beserr,bf);
					if (bf < 0.0)
						bf = 0.0;
					if (bf > 1.0)
						bf = 1.0;
					
					for	(e = 0; e < 3; e++)
						dp[e] = (int)(255.0 * (bf * rgb_high[e] + (1.0 - bf) * rgb_low[e]) + 0.5);

				}
				if (TIFFWriteScanline(tif, (tdata_t)pa, y, 0) < 0) {
					fprintf(stderr,"WriteScanline Failed at line %d\n",y);
					exit (-1);
				}
			}
			(void) TIFFClose(tif);
		}
#endif /* DUMP_EPERR */
	}
	
	/* Show veronoi planes by plotting the vertex network */
	if (dpla) {
		vtx *vx1, *vx2;

		if (x4a == NULL) {
			_o4 = s->tinp;
			if ((x4a = (double *)malloc(_o4 * sizeof(double))) == NULL)
				error ("ofps: malloc %d failed",_o4);
			if ((y4a = (double *)malloc(_o4 * sizeof(double))) == NULL)
				error ("ofps: malloc %d failed",_o4);
			if ((x5a = (double *)malloc(_o4 * sizeof(double))) == NULL)
				error ("ofps: malloc %d failed",_o4);
			if ((y5a = (double *)malloc(_o4 * sizeof(double))) == NULL)
				error ("ofps: malloc %d failed",_o4);
			if ((ocols = (plot_col *)malloc(_o4 * sizeof(plot_col))) == NULL)
				error ("ofps: malloc %d failed",_o4);
		}

		/* Add normal planes then subd planes, so that subd are always on top */
			o4 = 0;
#ifdef INDEP_SURFACE
		for (k = 0; k < 2; k++)		/* Do two passes */
#else
		for (k = 0; k < 1; k++)
#endif
		{
			/* Add node planes */
			for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link) {

				/* Don't plot faces involving the fake inside or outside node */
				for (e = 0; e <= di; e++) {
					if (vx1->nix[e] < -s->nbp)
						break;
				}
				if (e <= di)
					continue;

				for (j = 0; j < vx1->nnv; j++) {
					vx2 = vx1->nv[j];

					/* Don't plot faces involving the fake inside or outside node */
					for (e = 0; e <= di; e++) {
						if (vx2->nix[e] < -s->nbp)
							break;
					}
					if (e <= di)
						continue;

#ifdef INDEP_SURFACE
					if (sm_andtest(s, &vx1->vm, &s->sc[0].a_sm) == 0
					 || sm_andtest(s, &vx2->vm, &s->sc[0].a_sm) == 0) {	/* Subd plane */
						if (k == 0)
							continue;		/* Doing non-zubd pass */				
					} else {
						if (k == 1)
							continue;		/* Doing subd pass */
					}
#endif

					if (o4 >= _o4) {		/* need more space */
						_o4 *= 2;
						if ((x4a = (double *)realloc(x4a, _o4 * sizeof(double))) == NULL)
							error ("ofps: realloc x4a %d failed", _o4);
						if ((y4a = (double *)realloc(y4a, _o4 * sizeof(double))) == NULL)
							error ("ofps: realloc y4a %d failed", _o4);
						if ((x5a = (double *)realloc(x5a, _o4 * sizeof(double))) == NULL)
							error ("ofps: realloc x5a %d failed", _o4);
						if ((y5a = (double *)realloc(y5a, _o4 * sizeof(double))) == NULL)
							error ("ofps: realloc y5a %d failed", _o4);
						if ((ocols = (plot_col *)realloc(ocols, _o4 * sizeof(plot_col))) == NULL)
							error ("ofps: realloc y5a %d failed", _o4);
					}

					if (pcp != 0) {
						x4a[o4] = vx1->v[0];
						y4a[o4] = vx1->v[1];
						x5a[o4] = vx2->v[0];
						y5a[o4] = vx2->v[1];
					} else {
						x4a[o4] = vx1->p[0];
						y4a[o4] = vx1->p[1];
						x5a[o4] = vx2->p[0];
						y5a[o4] = vx2->p[1];
					}

#ifdef INDEP_SURFACE
					/* Show the sub dimension outline in apricot */
					if (k == 1) {
						ocols[o4].rgb[0] = 1.0;		/* Apricot */
						ocols[o4].rgb[1] = 0.52;
						ocols[o4].rgb[2] = 0.57;
					} else
#endif
					{
						ocols[o4].rgb[0] = 0.5;		/* Light Blue */
						ocols[o4].rgb[1] = 0.9;
						ocols[o4].rgb[2] = 0.9;
					}
					o4++;
				}
			}
		}
	}

	if ((s->nopstop >= 0 && s->optit < s->nopstop) || s->np < s->ntostop)
		dwt = 0;

	/* Plot the vectors */
	do_plot_vec2(minx, maxx, miny, maxy, 
				x1a, y1a, x2a, y2a, ntext, s->np, dwt,
	            x3a, y3a, mcols, mtext, dvx ? n3 : 0,
	            x4a, y4a, x5a, y5a, ocols, dpla ? o4 : 0);

}

#endif /* DEBUG || DUMP_PLOT */

/* ------------------------------------------------------------------- */
#ifdef SANITY_RESEED_AFTER_FIXUPS

/* Save the current used vertexes to the i_uvtx list, */
/* so that they can be verified against the re-seeded vertexes */
static void save_ivertexes(ofps *s) {
	vtx *vx, *nvx;

	s->i_uvtx = NULL;

	for (vx = s->uvtx; vx != NULL; vx = nvx) {
		nvx = vx->link;

		/* Remove the vertex from used and other lists */
		del_vtx1(s, vx);

		/* Add it to the i_uvtx list */
		vx->link = s->i_uvtx;
		s->i_uvtx = vx;
	}
}

/* Check the incremental vertexes against the re-seeded vertexes */
static int check_vertexes(ofps *s) {
	int i, j, e, k, di = s->di;
	vtx *v1, *v2;
	int fail = 0;

	printf("Verifying incremental vertexes against re-seeded:\n");

	/* For each reference (re-seeded) vertex */
	for (v1 = s->uvtx; v1 != NULL; v1 = v1->link) {

		/* Locate the equivalent incremental vertex */
		for (v2 = s->i_uvtx; v2 != NULL; v2 = v2->link) {
			for (e = 0; e <= di; e++) {
				if (v1->nix[e] != v2->nix[e])
					break;
			}
			if (e > di)
				break;		/* Found it */
		}
		if (v2 == NULL) {
			printf("Missing vertex no %d comb %s\n",v1->no,pcomb(di,v1->nix));
			fail = 1;
			continue;
		}

		/* Check the vertex location */
		for (e = 0; e < di; e++) {
			if (fabs(v1->p[e] - v2->p[e]) > 1e-5) {
				break;
			}
		}
		if (e < di) {
			printf("Vertex no %d (%d) comb %s in different location %s, should be %s\n",v1->no,v2->no,pcomb(di,v1->nix),ppos(di,v2->p),ppos(di,v1->p));
			fail = 1;
		}
		/* Check the eserr */
		if (fabs(v1->eserr - v2->eserr) > 1e-3) {
			printf("Vertex no %d (%d) comb %s has different eserr %f, should be %f\n",v1->no,v2->no,pcomb(di,v1->nix),v2->eserr,v1->eserr);
			fail = 1;
		}

		/* Check setmask */
		if (!_sm_equal(s, &v1->vm, &v2->vm)) {
			printf("Vertex no %d (%d) comb %s has different vm %s, should be %s\n",v1->no,v2->no,pcomb(di,v1->nix),psm(s,&v2->vm),psm(s,&v1->vm));
			fail = 1;
		}

		/* Check that the vertex nets are the same */
		for (i = 0; i < v1->nnv; i++) {
			vtx *vv1 = v1->nv[i];

			for (j = 0; j < v2->nnv; j++) {
				vtx *vv2 = v2->nv[j];

				for (e = 0; e <= di; e++) {
					if (vv1->nix[e] != vv2->nix[e])
						break;
				}
				if (e > di)
					break;		/* Found it */
			}
			if (j >= v2->nnv) {
				printf("Vertex no %d comb %s, i_ missing neighbour no %d comb %s\n",v1->no,pcomb(di,v1->nix),vv1->no,pcomb(di,vv1->nix));
				fail = 1;
			}
		}
		for (j = 0; j < v2->nnv; j++) {
			vtx *vv2 = v2->nv[j];

			for (i = 0; i < v1->nnv; i++) {
				vtx *vv1 = v1->nv[i];

				for (e = 0; e <= di; e++) {
					if (vv1->nix[e] != vv2->nix[e])
						break;
				}
				if (e > di)
					break;		/* Found it */
			}
			if (i >= v1->nnv) {
				printf("Vertex no %d comb %s, i_ extra neighbour no (%d) comb %s\n",v1->no,pcomb(di,v1->nix),vv2->no,pcomb(di,vv2->nix));
				fail = 1;
			}
		}
	}

	/* For each incremental vertex, check that there is a corresponding re-seeded vertex */
	for (v2 = s->i_uvtx; v2 != NULL; v2 = v2->link) {

		for (v1 = s->uvtx; v1 != NULL; v1 = v1->link) {
			for (e = 0; e <= di; e++) {
				if (v1->nix[e] != v2->nix[e])
					break;
			}
			if (e > di)
				break;		/* Found it */
		}
		if (v1 == NULL) {
			printf("Extra vertex no (%d) comb %s\n",v2->no,pcomb(di,v2->nix));
			fail = 1;
		}
	}

	if (fail)
		printf("Failed to verify incremental vertexes against re-seeded:\n");
	else
		printf("Successfully verified incremental vertexes against re-seeded\n");

	return fail;
}

#endif /* SANITY_RESEED_AFTER_FIXUPS */

/* ------------------------------------------------------------------- */
/* Do an exaustive, very slow check for missing vertexes */
/*
	This may be really, really, really slow. 

	For every possible combination of di+1 nodes,
	locate the corresponding vertex. If it is 
	locatable, check that no other node is closer to it.
	If it meets these conditions, then check that it is in the veronoi surface.
 */
static void check_for_missing_vertexes(ofps *s) {
	int e, di = s->di;
	vtx *vx;
	COMBO(co, MXPD+1, di+1, s->np + s->nbp);		/* di-1 out of neighbor nodes combination counter */
	nodecomb vv;
	int lsc = -100;
	int isok = 1;

	printf("Doing exaustive check for missing vertexes:\n");

	/* Mark all the vertexes so that we can tell if any are missed. */
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		vx->sch = 0;
	}

	CB_INIT(co);
	while (!CB_DONE(co)) {
		int rl = 0;

		memset((void *)&vv, 0, sizeof(nodecomb));
		for (e = 0; e <= di; e++) {
			vv.nix[e] = co[e] - s->nbp;
			if (vv.nix[e] >= 0)
				rl = 1;
		}
		if (rl == 0)
			goto next_comb;		/* No real nodes */

		vv.vv = NULL;
		vv.ceperr = 1e100;

		sort_nix(s, vv.nix);

		printf("Comb %s\n",pcomb(di, vv.nix));
		if (lsc != vv.nix[di]) {
			fprintf(stderr,"digit %d\n",vv.nix[di]);
			lsc = vv.nix[di];
		}

		if (position_vtx(s, &vv, 0, 0, 0) == 0) {
			int ix;
			double eperr;
			node *nn;

			printf(" Located at %s (%s), eperr %f\n",ppos(di,vv.p),ppos(di,vv.v),vv.eperr);

			/* Check that the point is not out of gamut */
			if (ofps_in_dev_gamut(s, vv.p, NULL) < -s->surftol) {
				printf(" vertex is out of gamut\n");
				goto not_valid;
			}

			/* Check that no other vertex is closer */
			for (ix = 0; ix < s->np; ix++) {
				for (e = 0; e <= di; e++) {
					if (vv.nix[e] == ix)
						break;
				}
				if (e <= di)
					continue;	/* Is a parent */

				nn = s->n[ix];
				eperr = ofps_comp_eperr(s, NULL, nn->v, nn->p, vv.v, vv.p, 0);

				printf(" eperr to ix %d is %f\n",nn->ix,eperr);
				if (eperr < vv.eperr) {
					printf("vertex is closer to node ix %d\n",nn->ix);
					break;
				} 
			}
			if (ix >= s->np) {

				printf("Point %s is valid\n",pcomb(di,vv.nix));

				/* see if we've created it */
				for (vx = s->uvtx; vx != NULL; vx = vx->link) {

					for (e = 0; e <= di; e++) {
						if (vx->nix[e] != vv.nix[e])
							break;
					}
					if (e > di)
						break;		/* Found it */
				}
				if (vx == NULL) {
					printf("Can't find vertex %s at %s (%s)\n",pcomb(di,vv.nix),ppos(di,vv.p),ppos(di,vv.v));
					fprintf(stderr,"Can't find vertex %s at %s (%s)\n",pcomb(di,vv.nix),ppos(di,vv.p),ppos(di,vv.v));
					isok = 0;
				} else {
					vx->sch = 1;
					printf("Found vertex no %d nix %s at %s (%s) OK\n",vx->no,pcomb(di,vv.nix),ppos(di,vv.p),ppos(di,vv.v));
					fprintf(stderr,"Found vertex no %d nix %s at %s (%s) OK\n",vx->no,pcomb(di,vv.nix),ppos(di,vv.p),ppos(di,vv.v));
				}
			}
	not_valid:;
		} else {
			printf(" Failed to locate %s\n",pcomb(di,vv.nix));
		}
	next_comb:;
		CB_INC(co);
	}
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		if (vx->sch == 0) {
			for (e = 0; e <= di; e++) {
				if (vx->nix[e] < -s->nbp)	/* involves inside or outside fake point */
					break;
			}
			if (e <= di)
				continue;	/* Ignore */
			printf("Extra vertex no %d nix %s at %s (%s) OK\n",vx->no,pcomb(di,vx->nix),ppos(di,vx->p),ppos(di,vx->v));
			fprintf(stderr,"Extra vertex no %d nix %s at %s (%s) OK\n",vx->no,pcomb(di,vx->nix),ppos(di,vx->p),ppos(di,vx->v));
			isok = 0;
		}
	}
	if (isok) {
		printf("Check for missing veftexes is OK\n");
		fprintf(stderr,"Check for missing veftexes is OK\n");
	} else {
		printf("Check for missing veftexes FAILED\n");
		fprintf(stderr,"Check for missing veftexes FAILED\n");
	}
}

/* ------------------------------------------------------------------- */
/* Check the veronoi to check that no node other than the parent */
/* node is closer to any vertex. */
/* return nz if there is a problem */
static int check_vertex_closest_node(ofps *s) {
	int i, e, di = s->di;
	node *nn, *pp;
	vtx *vx;

	/* Check that no node other than a parent is closer to any vertex */
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		double ceperr;

		if (vx->ofake)
			continue;

		/* Check if the vertex position is clipped by a gamut boundary. */
		for (i = 0; i < s->nbp; i++) {
			pleq *vp = &s->gpeqs[i];
			double v;
	
			pp = s->n[-1-i];
#ifdef INDEP_SURFACE
			/* Check if this vertex is visible to this node */
			if (sm_vtx_node(s, vx, pp) == 0) {
				continue;	/* It's hidden */
			}
#endif	/* INDEP_SURFACE */

			for (v = vp->pe[di], e = 0; e < di; e++)
				v += vp->pe[e] * vx->p[e];
			if (v > 2.0 * NUMTOL) {
				/* Check whether pp is already a parent of the node */
				for (e = 0; e <= di; e++) {
					if (pp->ix == vx->nix[e])
						break;
				}
				if (e <= di)
					continue;	/* It is */

#ifdef DEBUG
				printf("Vertex %d parents %s is clipped by boundary node %d by %e\n", vx->no,pcomb(di,vx->nix),pp->ix,v);
#endif
				warning("Vertex %d parents %s is clipped by boundary node %d by %e", vx->no,pcomb(di,vx->nix),pp->ix,v);
				return 1;
			}
		}

		/* locate the nearest node to the vertex */
		if ((nn = ofps_findclosest_node(s, &ceperr, vx)) == NULL)	
			continue;

		/* See if it is closer than the parent nodes */
		if ((vx->eperr - ceperr) < 2.0 * NUMTOL)
			continue;		/* No it's not */

		/* Check whether nn is already a parent of the node */
		for (e = 0; e <= di; e++) {
			if (nn->ix == vx->nix[e])
				break;
		}
		if (e <= di)
			continue;		/* A parent */

#ifdef DEBUG
		printf("Vertex %d is closer to %d (%f) than parent nodes %s (%f) by %e\n",vx->no,nn->ix,ceperr,pcomb(di,vx->nix),vx->eperr, ceperr - vx->eperr);
#endif
		warning("Vertex %d is closer to %d (%f) than parent nodes %s (%f) by %e",vx->no,nn->ix,ceperr,pcomb(di,vx->nix),vx->eperr, ceperr - vx->eperr);
		return 1;
	}
	fflush(stdout);

	return 0;
}

/* ------------------------------------------------------------------- */

#if defined(DEBUG) || defined(DUMP_PLOT) || defined (SANITY_CHECK_CONSISTENCY) || defined(DUMP_STRUCTURE)
/* Do some sanity checking on the points */
static void
sanity_check(
	ofps *s,
	int check_nodelists			/* nz to check node lists */
) {
	int i, j, k, e, di = s->di;
	vtx *vx1, *vx2;
	int fail = 0;		/* 0 = pass, 1 = soft fail, 2 = hard fail */

#ifdef DEBUG
	printf("Running sanity check...\n");
#endif

	/* See if any of the sample nodes are near the same location */
	for (i = 0; i < (s->np-1); i++) {
		node *p1 = s->n[i];
		for (j = i+1; j < s->np; j++) {
			node *p2 = s->n[j];
			double rad;
			for (rad = 0.0, e = 0; e < di; e++) {
				double tt = p1->p[e] - p2->p[e];
				rad += tt * tt;
			}
			rad = sqrt(rad);
			if (rad < 1e-5) {
#ifdef DEBUG
				printf("Nodes ix %d and ix %d are at %s and %s\n", i,j,ppos(di,p1->p),ppos(di,p2->p));
#endif
				fail = 2;
			}
		}
	}

	/* See if any of the vertexes have the same node combinations */
	for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link)
		vx1->sch = 0;

	for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link) {
		if (vx1->sch)
			continue;
		for (vx2 = vx1->link; vx2 != NULL; vx2 = vx2->link) {
			if (vx2->sch)
				continue;
			for (e = 0; e <= di; e++) {
				if (vx1->nix[e] != vx2->nix[e])
					break;
			}
			if (e > di) {
				vx1->sch = vx2->sch = 1;		/* Don't do these again */
#ifdef DEBUG
				printf("Vertex ix %d and ix %d have same nix %s\n", vx1->no,vx2->no,pcomb(di,vx1->nix));
#endif
				fail = 2;
			}
		}
	}

	/* See if any of the vertexes are at the same location */
	for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link)
		vx1->sch = 0;

	for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link) {
		if (vx1->sch)
			continue;
		for (vx2 = vx1->link; vx2 != NULL; vx2 = vx2->link) {
			double rad;
			if (vx2->sch)
				continue;
			for (rad = 0.0, e = 0; e < di; e++) {
				double tt = vx1->p[e] - vx2->p[e];
				rad += tt * tt;
			}
			rad = sqrt(rad);
			if (rad < 1e-10) {
				vx1->sch = vx2->sch = 1;		/* Don't do these again */
#ifdef DEBUG
				printf("Vertex no %d nix %s vm %s and no %d nix %s vm %s are at %s and %s", vx1->no,pcomb(di,vx1->nix),psm(s,&vx1->vm),vx2->no,pcomb(di,vx2->nix),psm(s,&vx2->vm),ppos(di,vx1->p),ppos(di,vx2->p));
				if (fabs(vx1->eperr - vx2->eperr) > 1e-5)
					printf(" and errs %f %f\n",vx1->eperr,vx2->eperr);
				else
					printf("\n");
#endif
				/* See if the two vertexes are both visible to each other */
				if (sm_vtx_vtx(s, vx1, vx2) != 0) {
					fail = 2;
				}
			}
		}
	}

	/* See if any of the nodes and vertexes are near the same location */
	for (i = 0; i < s->np; i++) {
		node *p1 = s->n[i];
		for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link) {
			double rad;
			for (rad = 0.0, e = 0; e < di; e++) {
				double tt = p1->p[e] - vx1->p[e];
				rad += tt * tt;
			}
			rad = sqrt(rad);
			if (rad < 1e-5) {
#ifdef DEBUG
				printf("Node ix %d and Vertex no %d are at %s and %s%s", i,vx1->no,ppos(di,p1->p),ppos(di,vx1->p),vx1->ghost ? " (ghost)" : "");
				if (vx1->eperr > 1e-5)
					printf(" and err %f\n",vx1->eperr);
				else
					printf("\n");
#endif
				if (vx1->ghost == 0)
					fail = 1;
			}
		}
	}

	/* Check every node appears in at least one vertex */
	for (i = 0; i < s->np; i++) {		/* For all nodes */
		node *p1 = s->n[i];
		vtx *vx;
		for (vx = s->uvtx; vx != NULL; vx = vx->link) {
			for (e = 0; e <= di; e++) {
				if (vx->nix[e] == p1->ix)
					break;		/* yes */
			}
			if (e <= di)
				break;		/* yes */
		}
		if (vx == NULL) {
#ifdef DEBUG
			printf("Node ix %d has no vertexes that refer to it\n", p1->ix);
#endif
			fail = 2;
		}
	}

	if (check_nodelists) {
		/* See if any vertexes do not appear in their constituent nodes */
		/* vertex list, or whether verexes nodes don't appear in neighbour list. */
		for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link) {
			for (e = 0; e <= di; e++) {
				int ix = vx1->nix[e];
				node *pp = s->n[ix];
	
				for (j = 0; j < pp->nvv; j++) {
					if (pp->vv[j] == vx1)
						break;
				} 
				if (j >= pp->nvv) {
#ifdef DEBUG
					printf("Vertex no %d nix %s doesn't appear in node ix %d\n", vx1->no,pcomb(di,vx1->nix),pp->ix);
#endif
					fail = 2;
				}
			}
		}
	}

	/* Check that every vertex of a node contains that node. */
	for (i = 0; i < s->np; i++) {		/* For all nodes */
		node *p1 = s->n[i];
		for (j = 0; j < p1->nvv; j++) {	/* For all its vertexes */
			vtx *vx = p1->vv[j];

			for (e = 0; e <= di; e++) {	/* All vertexes parent nodes */
				int ix = vx->nix[e];
				node *pp = s->n[ix];

				if (ix == p1->ix)
					break;
			}
			if (e > di) {
#ifdef DEBUG
				printf("Node ix %d has vtx no %d nix %s that doesn't contain node\n", p1->ix, vx->no,pcomb(di,vx->nix));
#endif
				fail = 2;
			}
		}
	}

	if (check_nodelists) {
		/* Check that a node contains as neighbours all the parent */
		/* nodes of its vertexes */
		for (i = 0; i < s->np; i++) {		/* For all nodes */
			node *p1 = s->n[i];
			for (j = 0; j < p1->nvv; j++) {	/* For all its vertexes */
				vtx *vx = p1->vv[j];

				for (e = 0; e <= di; e++) {
					int ix = vx->nix[e];
					node *pp = s->n[ix];

					if (ix == p1->ix)
						continue;			/* Neighbours don't include self */
					for (k = 0; k < p1->nvn; k++) {
						if (p1->vn[k] == ix)
							break;
					}
					if (k >= p1->nvn) {
#ifdef DEBUG
						printf("Node ix %d has vtx no %d nix %s where neighbour ix %d is missing\n", p1->ix, vx->no,pcomb(di,vx->nix),ix);
#endif
						fail = 2;
					}
				}
			}
		}
	}

	/* Check that the vertex net is correct */
	{
		int ff, f, e, di = s->di;
		vtx *vx1, *vx2;
		int nnv = 0;
		int _nnv = 0;
		struct _vtx **nv = NULL;
		
		/* Do a brute force search to locate all this vertexes net neighbours */
		for (vx1 = s->uvtx; vx1 != NULL; vx1 = vx1->link) {

			nnv = 0;		/* Clear the current list */

			/* Search all other vertexes for neighbours */
			for (vx2 = s->uvtx; vx2 != NULL; vx2 = vx2->link) {
				int aa, bb, cc;		/* Probable hit check */
				int nnm, nmix;

				if (vx1 == vx2)
					continue;

#ifdef NEVER	/* vertex net needs all neighbours ? */
#ifdef INDEP_SURFACE
			    if (sm_vtx_vtx(s, vx1, vx2) == 0)
					continue;
#endif /* INDEP_SURFACE */
#endif

				/* Use the nixm to quickly check if all but one parent node matches */
				aa = vx1->nix[MXPD+2];	/* nixm */
				bb = vx2->nix[MXPD+2];	/* nixm */
				if ((aa & bb) == 0 || (cc = aa & ~bb, (cc & (cc-1)) != 0))
					continue;		/* It's certainly not */

				/* Do an exact check of all except one node match */
				for (nnm = ff = e = 0; e <= di; e++) {
					for (f = ff; f <= di; f++) {
						if (vx1->nix[e] == vx2->nix[f]) {
							ff = f;			/* Start from here next time */
							break;
						}
						if (vx1->nix[e] > vx2->nix[f])	/* No point in looking further */
							f = di;
					}
					if (f > di) {	/* Didn't match */
						if (++nnm > 1)
							break;
						nmix = e;
					}
				}
				if (e <= di)
					continue;			/* No match */
				
				if (nnm == 0) {
					error("ofps: two vertexes have the same nodes !\n"
						  "no %d at %s nix %s\nno %d at %s nix %s",
					vx1->no,ppos(di,vx1->p),pcomb(di,vx1->nix),
					vx2->no,ppos(di,vx2->p),pcomb(di,vx2->nix));
				}
				if (nnv >= _nnv) {
					_nnv = 2 * _nnv + 1;
					if ((nv = (vtx **)realloc(nv, sizeof(vtx *) * _nnv)) == NULL)
						error("ofps: realloc failed on node vertex pointers");
				}
				nv[nnv++] = vx2;
			}

			/* Now check that the vertex nets match */
			for (i = 0; i < nnv; i++) {
				for (j = 0; j < vx1->nnv; j++) {
					if (nv[i] == vx1->nv[j])
						break;
				}
				if (j >= vx1->nnv) {
					printf("Vtx no %d is missing vtx no %d from net\n",vx1->no,nv[i]->no);
					fail = 2;
				}
			}
			for (j = 0; j < vx1->nnv; j++) {
				for (i = 0; i < nnv; i++) {
					if (nv[i] == vx1->nv[j])
						break;
				}
				if (i >= nnv) {
					printf("Vtx no %d has extra vtx no %d in net\n",vx1->no,vx1->nv[j]->no);
					fail = 2;
				}
			}
		}
	}
	if (fail) {
		if (fail == 1)
			warning("Internal consistency check failed (soft)");
		else
			warning("Internal consistency check failed");
#ifdef DEBUG
		if (fail == 1)
			printf("Internal consistency check failed (soft)\n");
		else
			printf("Internal consistency check failed\n");
		fflush(stdout);
#endif
#ifdef SANITY_CHECK_CONSISTENCY_FATAL
		if (fail == 2)
			error("Internal consistency check failed");
#endif
	}

#ifdef SANITY_CHECK_EXAUSTIVE_SEARCH_FOR_VERTEXES
	check_for_missing_vertexes(s);
#endif
}
#endif /* SANITY_CHECK_CONSISTENCY */

#if defined(DEBUG) || defined(DUMP_STRUCTURE)

/* ------------------------------------------------------------------- */

/* Dump the node & vertex relationship */
static void
dump_node_vtxs(
	ofps *s,
	int check_nodelists
) {
	int i, j, e, di = s->di;
	vtx *vx;

	printf("\n");
	printf("Dumping current state...\n");

	/* Dump node information */
	for (i = -s->gnp; i < s->np; i++) {
		node *p1 = s->n[i];
		printf("Node ix %d, pos %s, mask 0x%x, asm %s\n",p1->ix,ppos(di, p1->p),p1->pmask,psm(s,&s->sc[p1->pmask].a_sm));
	}
	printf("\n");

	/* Dump vertex information */
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		if (vx->ofake == 0)
			printf("Vertex no %d, pmask 0x%x, cmask 0x%x, vm %s\n      pos %s, nix %s, eperr = %s, eserr = %s\n",vx->no,vx->pmask,vx->cmask,psm(s,&vx->vm), ppos(di, vx->p), pcomb(di,vx->nix), peperr(vx->eperr), peperr(vx->eserr));
		else
			printf("Vertex no %d, pmask 0x%x, cmask 0x%x, vm %s, nix %s, eperr = %s, eserr = %s (ofake)\n",vx->no,vx->pmask,vx->cmask,psm(s,&vx->vm), pcomb(di,vx->nix), peperr(vx->eperr), peperr(vx->eserr));
	}
	printf("\n");

	/* Dump vertex and associated vertex information */
	for (vx = s->uvtx; vx != NULL; vx = vx->link) {
		printf("Vertex no %d has Vtx net:",vx->no);
		for (j = 0; j < vx->nnv; j++) {
			vtx *vx2 = vx->nv[j];
			printf(" %d",vx2->no);
		}
		printf("\n");
	}
	printf("\n");

	/* Dump node and associated vertex information */
	for (i = -s->nbp; i < s->np; i++) {
		node *p1 = s->n[i];
		printf("Node ix %d, pos %s, mask 0x%x, a_sm %s:\n",p1->ix,ppos(di, p1->p),p1->pmask,psm(s,&s->sc[p1->pmask].a_sm));
		for (j = 0; j < p1->nvv; j++) {
			vtx *vx = p1->vv[j];
			if (vx->ofake == 0)
				printf("  Vtx no %d pmask 0x%x cmask 0x%x vm %s pos %s nix %s eserr %s\n",vx->no,vx->pmask,vx->cmask,psm(s,&vx->vm),ppos(di, vx->p), pcomb(di, vx->nix), peperr(vx->eserr));
			else
				printf("  Vtx no %d pmask 0x%x cmask 0x%x vm %s nix %s eserr %s (ofake)\n",vx->no,vx->pmask,vx->cmask,psm(s,&vx->vm),pcomb(di,vx->nix), peperr(vx->eserr));
		}
		for (j = 0; j < p1->nvn; j++) {
			int ix = p1->vn[j];
			if (ix >= 0) {
				node *n1 = s->n[ix];
				printf("  Assoc. node ix %d pos %s\n",ix,ppos(di, n1->p));
			} else {
				printf("  Assoc. node ix %d\n",ix);
			}
		}
	}
	printf("\n");

	sanity_check(s, check_nodelists);

	fflush(stdout);
}

#ifdef NEVER
/* Special dump the node & vertex relationship, */
/* for comparing with "good" output. */
/* Deal with vertexe order and numbering. */
/* Note that the "new" "bad" code needs ofake vertexes */
/* to work, wheras the "old" "good" code doesn't, so */
/* skip reporting ofake vetexes. */
static void
dump_node_vtxs2(
	ofps *s,
	char *com
) {
	int i, j, e, di = s->di;
	vtx *vx;
	FILE *fp;
	static int cc = 0;
	int showofake = 1;
	vtx **vlist;

	if (cc == 0) {
		if ((fp = fopen("bad.log","w")) == NULL)
			error("Unable to open file '%s'\n","bad.log");
		cc = 1;
	} else {
		if ((fp = fopen("bad.log","a")) == NULL)
			error("Unable to open file '%s'\n","bad.log");
	}

	fprintf(fp,"\n");
	fprintf(fp,"Dumping current state (%s) ...\n",com);

	/* Dump node information */
	for (i = -s->gnp; i < s->np; i++) {
		node *p1 = s->n[i];
		fprintf(fp,"Node ix %d, pos %s, mask 0x%x, asm %s\n",p1->ix,ppos(di, p1->p),p1->pmask,psm(s,&s->sc[p1->pmask].a_sm));
	}
	fprintf(fp,"\n");

	/* Sort the vertexes by their nix */
	{
		int scl = s->gnp + s->np;
		int nv;

		int nused;
		for (nused = 0, vx = s->uvtx; vx != NULL; vx = vx->link)
			nused++;
if (nused != s->nv) error("s->nv %d doesn't match uvtx list %d",s->nv,nused);

//printf("~1 number of vertexes = %d\n",s->nv);
		if ((vlist = (vtx **)malloc(sizeof(vtx *) * s->nv)) == NULL)
			error ("ofps: malloc failed on sorted vertex list");
		for (i = nv = 0, vx = s->uvtx; vx != NULL; vx = vx->link, i++) {
			if (!showofake && vx->ofake)
				continue;
			vlist[nv++] = vx;
	
			/* Convert nix into sort index */
			vx->sch = 0;
			for (e = 0; e <= di; e++) {
				vx->sch = scl * vx->sch + (vx->nix[e] + s->gnp); 
			}
//printf("~1 nix %s ix %d\n",pcomb(di,vx->nix),vx->sch);
//fflush(stdout);
		}

		/* Sort */
#define HEAP_COMPARE(A,B) ((A)->sch < (B)->sch)
		HEAPSORT(vtx *, vlist, nv);
#undef HEAP_COMPARE

		for (i = 0; i < nv; i++) {
			vx = vlist[i];
			vx->sch = i;		/* Sorted index */
		}

		/* Dump vertex information */
		for (i = 0; i < nv; i++) {
			vx = vlist[i];
			fprintf(fp,"Vertex no %d, pmask 0x%x, cmask 0x%x, vm %s\n      pos %s, nix %s, eserr = %s\n",vx->sch,vx->pmask,vx->cmask,psm(s,&vx->vm), ppos(di, vx->p), pcomb(di,vx->nix), peperr(vx->eserr));
		}
		fprintf(fp,"\n");
		free(vlist);
	}

	/* Dump node and associated vertex information */
	for (i = -s->nbp; i < s->np; i++) {
		node *p1 = s->n[i];
		vtx **vv;
		int nvv;
		int *vn;

		fprintf(fp,"Node ix %d, pos %s, mask 0x%x, a_sm %s\n",p1->ix,ppos(di, p1->p),p1->pmask,psm(s,&s->sc[p1->pmask].a_sm));

		/* Display the vertexes in order */
		if ((vv = (vtx **)malloc(sizeof(vtx *) * p1->nvv)) == NULL)
			error ("ofps: malloc failed on sorted vertex list");
		for (nvv = j = 0; j < p1->nvv; j++) {
			if (!showofake && p1->vv[j]->ofake)
				continue;
			vv[nvv++] = p1->vv[j];
		}
#define HEAP_COMPARE(A,B) ((A)->sch < (B)->sch)
		HEAPSORT(vtx *, vv, nvv);
#undef HEAP_COMPARE
		for (j = 0; j < nvv; j++) {
			vtx *vx = vv[j];
			fprintf(fp,"  Vtx no %d pmask 0x%x cmask 0x%x vm %s pos %s nix %s eserr %s\n",vx->sch,vx->pmask,vx->cmask,psm(s,&vx->vm),ppos(di, vx->p), pcomb(di, vx->nix), peperr(vx->eserr));
		}
		free(vv);

		/* Sort the nodes to be in order */
		if ((vn = (int *)malloc(sizeof(int) * p1->nvn)) == NULL)
			error ("ofps: malloc failed on sorted vertex list");
		for (j = 0; j < p1->nvn; j++)
			vn[j] = p1->vn[j];
#define HEAP_COMPARE(A,B) ((A) < (B))
		HEAPSORT(int, vn, p1->nvn);
#undef HEAP_COMPARE
		for (j = 0; j < p1->nvn; j++) {
			int ix = vn[j];
			if (ix >= 0) {
				node *n1 = s->n[ix];
				fprintf(fp,"  Assoc. node ix %d pos %s\n",ix,ppos(di, n1->p));
			} else {
				fprintf(fp,"  Assoc. node ix %d\n",ix);
			}
		}
		free(vn);
	}
	printf("\n");
	fflush(fp);
	fclose(fp);
}
#endif /* NEVER */
#endif /* DEBUG || DUMP_PLOT || DUMP_STRUCTURE */


/* --------------------------------------------------------------- */
#ifdef DUMP_FERR		/* Create .tiff of dnsq function error */

/* Draw a line in the output diagnostic raster */
static int
show_line(
ofps *s,							/* ofps object */
int x1, int y1, int x2, int y2,		/* line start and end points */
unsigned char rgb[3],						/* Color */
unsigned char *base,				/* Raster base of line */
int pitch,
int width,
int height
) {
	unsigned char *pp;
	int ow = width, oh = height;	/* width and height of raster for clipping */
	int dx, dy;			/* Line deltas */
	int adx, ady;		/* Absolute deltas */

	int e, k1, k2;		/* Error and axial/diagonal error change values */
	int m1,m2;		/* axial/diagonal coordinate change values */

	int ll;				/* Line length */

	/* Do a crude clip */
	if (x1 < 0)
		x1 = 0;
	if (x1 >= ow)
		x1 = ow-1;
	if (x2 < 0)
		x2 = 0;
	if (x2 >= ow)
		x2 = ow-1;
	if (y1 < 0)
		y1 = 0;
	if (y1 >= oh)
		y1 = oh-1;
	if (y2 < 0)
		y2 = 0;
	if (y2 >= oh)
		y2 = oh-1;

	/* calculate the standard constants */
	dx = x2 - x1;
	dy = y2 - y1;

	if(dx  < 0) {
		m1 = -3;		/* x is going backwards */
		adx = -dx;		/* make this absolute */
	} else {
		m1 = 3;			/* x is going forwards */
		adx = dx;
	}
	
	e = 0;
	if(dy < 0) {
		m2 = -pitch;	/* y is going upwards (decreasing) */
		ady = -dy;		/* make this absolute */
		e = -1;			/* make lines retraceable */
	} else {
		m2 = pitch;		/* y is going downwards (increasing) */
		ady = dy;
	}

	/* m1 has been set to x increment, m2 to y increment */

	m2 += m1;			/* make m2 the diagonal address increment */
						/* and m1 the x axial inrement */
	if(adx > ady) {			/* x is driven */
		ll = adx;
		k1 = 2 * ady;
		k2 = 2 * (ady - adx);
		e += k1 - adx;
	} else {
		ll = ady;
		k1 = 2 * adx;
		k2 = 2 * (adx - ady);
		e += k1 - ady;
		m1 = m2 - m1;		/* Make m1 the y increment */
	}

	/* Start pixel of line */
	pp = base + y1 * pitch + 3 * x1;

	ll++;	/* Draw start and end point */

	while( ll > 0) {
		while(e < 0 && ll > 0) {
			pp[0] = rgb[0];
			pp[1] = rgb[1];
			pp[2] = rgb[2];
			pp += m1;
			e += k1;
			ll--;
		}
		while(e >= 0 && ll > 0) {
			pp[0] = rgb[0];
			pp[1] = rgb[1];
			pp[2] = rgb[2];
			pp += m2;
			e += k2;
			ll--;
		}
	}
	return 0;
}

/* Dump a TIFF of the dnsq function values for a given point/plane set. */
static void
dump_dnsqe(
	ofps *s,
	char *fname,
	int *nix,
	vopt_cx *cx
) {
	int i, j, e, di = s->di;
	unsigned char *base, *pa, col[2][3];
	int width = WIDTH;
	int height = HEIGHT;
	int pitch = width * 3; 
	int x, y;
	TIFF *tif;
	double pos[MXPD], fval[MXPD];
	double angle, mag;

	printf("Dumping dnsqe error for combination %s\n",pcomb(di,nix));

	if ((tif = TIFFOpen(fname, "w")) == NULL) {
		fprintf(stderr,"Failed to open output TIFF file '%s'\n",fname);
		exit (-1);
	}

	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,  width);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

	/* allocate a raster */
	if ((base = (unsigned char *)malloc(sizeof(unsigned char) * height * pitch)) == NULL)
		error ("ofps: malloc failed on diagnostic raster");
	
	for (y = 0; y < height; y++) {
		pos[1] = 1.0 - y/(height-1.0);
		pos[1] = 1.4 * pos[1] - 0.2;
		pa = base + y * pitch;

		/* Fill in pa[] with colors for this line */
		for (x = 0; x < width; x++) {
			double ss;
			unsigned char *dp;
			double beserr, eserr;
			double bf, rgb[3];
			double oog = 1.0;
			double escale = 10.0;		/* Error value scaling */

			dp = pa + x * 3;
			pos[0] = x/(width-1.0);
			pos[0] = 1.4 * pos[0] - 0.2;
			dp[0] = dp[1] = dp[2] = 255;
//printf("~1 doing %d %d pos %f %f\n",x,y,pos[0],pos[1]);

			/* Se if the sample is in gamut */
			for (ss = 0.0, e = 0; e < s->di; e++) {
				if (pos[e] < s->imin[e]
				 || pos[e] > s->imax[e])
					break;
				ss += pos[e];
			}
			if (e < s->di || ss > (s->ilimit + ILIMITEPS)) {
				oog = 0.7;			/* Show gamut boundary */
			}

#ifdef NEVER		/* Test colors out */
			fval[0] = pos[0] * 2.0 * escale - escale;
			fval[1] = pos[1] * 2.0 * escale - escale;
#else
			/* Lookup the function value here */
			dnsq_solver(cx, di, pos, fval, 0);  
#endif

			/* Turn the two values into colors. */
			for (e = 0; e < di; e++) {
				fval[e] = (fval[e] / escale);
				if (fval[e] >= 0.0)
					fval[e] = pow(fval[e], 0.5);
				else
					fval[e] = -pow(-fval[e], 0.5);
			}

			/* Convert to angle and magnitude */
			angle = 180.0/3.1415926 * atan2(fval[0], fval[1]);
			if (angle < 0.0)
				angle += 360.0;
			else if (angle > 360.0)
				angle -= 360.0;
			mag = sqrt(fval[0] * fval[0] + fval[1] * fval[1]);
			if (mag > 1.0)
				mag = 1.0;

			rgb[0] = rgb[1] = rgb[1] = 0.0;
			if (angle < 120.0) {				/* red to green */
				bf = angle / 120.0;
				rgb[0] = 1.0 - bf;
				rgb[1] = bf;
				rgb[2] = 0.0;
			} else if (angle < 240.0) {		/* green to blue */
				bf = (angle - 120.0) / 120.0;
				rgb[0] = 0.0;
				rgb[1] = 1.0 - bf;
				rgb[2] = bf;
			} else {						/* blue to red */
				bf = (angle - 240.0) / 120.0;
				rgb[0] = bf;
				rgb[1] = 0.0;
				rgb[2] = 1.0 - bf;
			}

			/* Scale to black with magnitude */
			for	(e = 0; e < 3; e++) {
				rgb[e] = 1.0 - rgb[e];
				rgb[e] = (1.0 - mag) * 0.0 + mag * rgb[e];
			}

			for	(e = 0; e < 3; e++)
				dp[e] = (int)(255.0 * oog * rgb[e] + 0.5);
		}
	}


	/* Show the path the dnsq sampled */
	col[0][0] = col[0][1] = 255, col[0][2] = 128;
	col[1][0] = 128, col[1][1] = col[1][2] = 255;

	for (i = 0; i < (cx->nl-1); i++) {
//printf("~1 line %d: %f %f -> %f %f\n",i, cx->clist[i].p[0], cx->clist[i].p[1], cx->clist[i+1].p[0], cx->clist[i+1].p[1]);
		show_line(s, 
			(int)(((cx->clist[i].p[0] + 0.2) / 1.4) * (width - 1.0) + 0.5),
			(int)((1.0 - ((cx->clist[i].p[1] + 0.2) / 1.4)) * (height - 1.0) + 0.5),
			(int)(((cx->clist[i+1].p[0] + 0.2) / 1.4) * (width - 1.0) + 0.5),
			(int)((1.0 - ((cx->clist[i+1].p[1] + 0.2) / 1.4)) * (height - 1.0) + 0.5),
			col[i & 1], base, pitch, width, height); 
	}

	/* Write the raster out */
	for (y = 0; y < height; y++) {
		pa = base + y * pitch;

		if (TIFFWriteScanline(tif, (tdata_t)pa, y, 0) < 0) {
			fprintf(stderr,"WriteScanline Failed at line %d\n",y);
			exit (-1);
		}
	}
	(void) TIFFClose(tif);
	free(base);
}
#endif /* DUMP_FERR */

/* --------------------------------------------------------------- */

#ifdef NEVER

		/* Compute an aproximate bounding shere, and use */
		/* the center of it as the start point. */

		double radsq = -1.0;			/* Span/radius squared */
		double rad;
		double sum;
		int i, j;
		int bi = 0, bj = 0;

		/* Find the two vectors that are farthest apart. Brute force search */
		/* Also track the device position for the points used to define the shere */
		for (i = 0; i < (ii-1); i++) {
			for (j = i+1; j < ii; j++) {
				for (sum = 0.0, e = 0; e < di; e++) {
					double tt = cx.nds[i]->p[e] - cx.nds[j]->p[e];
					sum += tt * tt;
				}
				if (sum > radsq) {
					radsq = sum;
					bi = i;
					bj = j;
				}
			}
		}
		
		/* Set initial bounding sphere */
		for (e = 0; e < di; e++)
			atp[e] = 0.5 * (cx.nds[bi]->p[e] + cx.nds[bj]->p[e]);
		radsq /= 4.0;			/* diam^2 -> rad^2 */
		rad = sqrt(radsq);

		/* Go though all the points again, expanding sphere if necessary */
		for (i = 0; i < ii; i++) {

			if (i == bi || i == bj)
				continue;

			/* Compute distance squared of vertex to bounding sphere center */
			for (sum = 0.0, e = 0; e < di; e++) {
				double tt = cx.nds[i]->p[e] - atp[e];
				sum += tt * tt;
			}
			if (sum > radsq) {
				double tt;

				sum = sqrt(sum) + 1e-10;			/* Radius to point */
				rad = 0.5 * (rad + sum);
				radsq = rad * rad;
				tt = sum - rad;
				for (e = 0; e < di; e++)
					atp[e] = (rad * atp[e] + tt * cx.nds[i]->p[e])/sum;
			}
		}

/* Given two sample point indexes, compute the plane between them. */
/* (This will fail with a divide by zero error if two points are coincident) */
static void comp_pleq(ofps *s, pleq *vp, int ix1, int ix2) {
	node *p0 = s->n[ix1], *p1 = s->n[ix2];
	int e, di = s->di;
	double cp[MXPD];
	double sum = 0.0;

	/* Compute plane normal from ix1 to ix2 */
	for (e = 0; e < di; e++) {
		double tt = p1->p[e] - p0->p[e];
		vp->pe[e] = tt;
		sum += tt * tt;
	}
	sum = sqrt(sum);

	/* Normalise it */
	for (e = 0; e < di; e++)
		vp->pe[e] /= sum;

	/* Compute mid point */
	for (e = 0; e < di; e++) 
		cp[e] = 0.5 * (p1->p[e] + p0->p[e]);

	/* Compute the plane equation constant */
	for (vp->pe[di] = 0.0, e = 0; e < di; e++) 
		vp->pe[di] -= vp->pe[e] * cp[e];
}

#endif // NEVER

