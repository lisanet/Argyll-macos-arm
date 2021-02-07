
/* First cut at lchw weighted. Problems with list size, memory use and */
/* performance. Version uses direct bwd cell nnrev[] creation */ 

/* 
 * Argyll Color Correction System
 * Multi-dimensional regularized spline data structure
 *
 * Reverse interpolation support code.
 *
 * Author: Graeme W. Gill
 * Date:   30/1/00
 *
 * Copyright 1999 - 2008 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Latest simplex/linear equation version.
 */

/* TTBD:

	Add option/function to return a gamut surface triangle list
	based on the rev setup thinned vertex list.
	Need to add code to convert over ink edges to triangles
	and then shadow test them though.

	XYZ PCS doesn't work with a LCh weighting, although this is
	no an issue when xicc uses separate Jab rspl for clip case (CAM CLIP).

	Allow function callback to set auxiliary values for 
	flag RSPL_AUXLOCUS. 
	How to pass enough info back to aux_compute() ?

	Should auxil return multiple solutions if it finds them ???

	Sometimes slivers remain in the surface in the exact
	direction of the focal point. See test/HarveyMiller colprof -qu
	with #define REVVRML. Probably not actually a problem, just not 100%
	correct gamut surface.

 */

/* TTBD:

	Get rid of error() calls - return status instead

	Need to add a hefty overview and explanation of
	how all this works, before I forget it !

	ie:

	  Basic function requirements:  exact, auxil, locus, clip

	  Fwd cell - fxcell list lookup

	  Basic layout di -> fdi + auxils + ink limit

	  Basic search strategy

	  Sub Simplex decomposition & properties

	  How each type of function finds solutions
		Sub-simplex dimensionality & dof + target dim & dof
		Linear algebra choices.
		
	  How final solutions are chosen

 */

/* PROBLEMS:

	Sometimes the aux locus doesn't correspond exactly to
	the inversion :- ie. one locus segment is returned,
	yet the inversion can't return a solution with
	a particular aux target that lies within that segment.
	(1150 near black, k ~= 0.4).

 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <memory.h>
#include <time.h>

#ifdef NT 
# ifdef WINVER
#  undef WINVER
# endif
# define WINVER 0x0500		/* We need 2k features */
# include <windows.h>
#else
# include <unistd.h>
# ifdef __APPLE__
#  include <fcntl.h>
#  include <sys/types.h>
#  include <sys/sysctl.h>
# endif
#endif

#define INKSCALE 5000.0	/* For ink limit weighting to fudge SVD least squares solution */

#include "rspl_imp.h"
#include "numlib.h"
#include "sort.h"		/* Heap sort */
#include "counters.h"	/* Counter macros */

//#define DMALLOC_GLOBALS
//#include "dmalloc.h"
//#undef DMALLOC_GLOBALS

#define DOSORT			/* [def] Cell sort for better speed */
#undef EN_UNTWIST		/* [und] Force attempt to try and untwist gamut surface */
						/* - Seems to improve some, make some worse ?? (i.e Bonet) */
						/* By default this is controlled using ARGYLL_UNTWIST_GAMUT_SURFACE */
						/* environment variable. */

#undef REVTABLESTATS	/* [und] Reverse table stats */
#undef REVVRML			/* [und] Reverse table 3D plots */

#undef DEBUG1			/* [und] Higher level code */
#undef DEBUG2			/* [und] Lower level code */

/* Debug memory usage accounting */
#ifdef NEVER
#ifdef NEVER
int thissz, lastsz = -1;
#define INCSZ(s, bbb) {						\
					(s)->rev.sz += (bbb);	\
					(s)->rev.thissz = (s)->rev.sz/1000000;		\
                    if ((s)->rev.thissz != (s)->rev.lastsz) fprintf(stderr,"~1 0x%x: %s, %d: rev size = %d Mbytes, delta %d, limit %d\n",((int)(s) >> 8) & 0xf, __FILE__, __LINE__,(s)->rev.thissz,(bbb),(s)->rev.max_sz/1000000);	\
					(s)->rev.lastsz = (s)->rev.thissz;	\
					}
#define DECSZ(s, bbb) {						\
					(s)->rev.sz -= (bbb);	\
					(s)->rev.thissz = (s)->rev.sz/1000000;		\
                    if ((s)->rev.thissz != (s)->rev.lastsz) fprintf(stderr,"~1 0x%x: %s, %d: rev size = %d Mbytes, delta %d, limit %d\n",((int)(s) >> 8) & 0xf, __FILE__, __LINE__,(s)->rev.thissz,-(bbb),(s)->rev.max_sz/1000000);	\
					(s)->rev.lastsz = (s)->rev.thissz;	\
					}
#else
#define INCSZ(s, bbb) (s)->rev.sz += (bbb);	\
                     fprintf(stderr,"%s, %d: rev.sz += %d\n",__FILE__, __LINE__, bbb)
#define DECSZ(s, bbb) (s)->rev.sz -= (bbb);	\
                     fprintf(stderr,"%s, %d: rev.sz -= %d\n",__FILE__, __LINE__, bbb)
#endif
#else
#define INCSZ(s, bbb) (s)->rev.sz += (bbb)
#define DECSZ(s, bbb) (s)->rev.sz -= (bbb)
#endif

/* Set STATS in rev.h */

/* Print a vectors value */
#define DBGVI(text, dim, out, vec, end)			\
{	int pveci;									\
	printf("%s",text);							\
	for (pveci = 0 ; pveci < (dim); pveci++)		\
		printf(out,(vec)[pveci]);				\
	printf(end);								\
}

/* Print a matrix value */
#define DBGMI(text, rows, cols, out, mat, end)		\
{	int pveci, pvecr;								\
	printf("%s",text);								\
	for (pvecr = 0 ; pvecr < (rows); pvecr++) {		\
		for (pveci = 0 ; pveci < (cols); pveci++)		\
			printf(out,(mat)[pvecr][pveci]);		\
	if ((pvecr+1) < (rows))							\
		printf("\n");								\
	}												\
	printf(end);									\
}

#if defined(DEBUG1) || defined(DEBUG2)
# define REVTABLESTATS	/* [und] Reverse table stats */
#endif

#ifdef REVTABLESTATS
# pragma message("!!!!!!!!! REVTABLESTATS set in rspl/rev.c !!!!!!!!!!!")
#endif

#ifdef REVVRML
# pragma message("!!!!!!!!! REVVRML set in rspl/rev.c !!!!!!!!!!!")
# include "vrml.h"
#endif

#ifdef CHECK_NNLU
# pragma message("!!!!!!!!! CHECK_NNLU set in rspl/rspl.h !!!!!!!!!!!")
#endif

/* Do an arbitrary printf */
#define DBGI(text) printf text ;

#undef DEBUG
#undef DBG
#undef DBGV
#undef DBGM

#undef NEVER
#define ALWAYS

#ifdef DEBUG1
#undef DBGS
#undef DBG
#undef DBGV
#undef DBGM
#define DEBUG
#define DBGS(xxx) xxx
#define DBG(xxx) DBGI(xxx)
#define DBGV(xxx) DBGVI xxx
#define DBGM(xxx) DBGMI xxx
#else
#undef DEBUG
#undef DBGS
#undef DBG
#undef DBGV
#undef DBGM
#define DBGS(xxx) 
#define DBG(xxx) 
#define DBGV(xxx) 
#define DBGM(xxx) 
#endif

/* Debug string routines */
static char *pcellorange(fxcell *c);

/* Convention is to use:
   i to index grid points u.a
   n to index data points d.a
   e to index position dimension di
   f to index output function dimension fdi
   j misc and cube corners
   k misc
 */

#define	EPS (2e-6)			/* 2e-6 Allowance for numeric error */

static void make_rev(rspl *s);
static void init_revaccell(rspl *s);

static fxcell *get_fxcell(schbase *b, int ix, int force);
static void uncache_fxcell(revcache *r, fxcell *cp);
#define unget_fxcell(r, cp) uncache_fxcell(r, cp)		/* These are the same */
static void invalidate_revaccell(rspl *s);
static int decrease_revcache(revcache *rc);

/* ====================================================== */

static schbase *init_search(rspl *s, int flags, double *av, int *auxm,
                        double *v, double *cdir, co *cpp, int mxsoln, enum ops op);
static void adjust_search(rspl *s, int flags, double *av, enum ops op);
static schbase *set_search_limit(rspl *s, double (*limit)(void *vcntx, double *in),
                                 void *lcntx, double limitv);
static void set_lsearch(rspl *s, int e);
static void free_search(schbase *b);

static int *calc_fwd_cell_list(rspl *s, double *v);

static int *calc_fwd_nn_cell_list(rspl *s, double *v);

static void init_line_eq(schbase *b, double st[MXRO], double de[MXRO]);
static int *init_line(rspl *s, line *l, double st[MXRO], double de[MXRO]);
static int *next_line_cell(line *l);

static void search_list(schbase *b, int *rip, unsigned int tcount);

static void clear_limitv(rspl *s);

static double get_limitv(schbase *b, int ix,	float *fcb, double *p);

#ifdef STATS
static char *opnames[6] = { "exact", "clipv", "clipn", "auxil", "locus" };
#endif /* STATS */

#define INF_DIST 1e38		/* Stands for infinite "current best" distance */

/* ====================================================== */
/* Globals that track overall usage of reverse cache to aportion memory */
/* This is incremented for rspl with di > 1 when rev.rev_valid != 0 */
size_t g_avail_ram = 0;			/* Total maximum memory to be used */
size_t g_test_ram = 0;			/* Amount of memory that has been tested to be allocatable */
int g_no_rev_cache_instances = 0;
rev_struct *g_rev_instances = NULL;

/* ------------------------------------------------------ */
/* Retry allocation routines - if the malloc fails,       */
/* try reducing the cache size and trying again */
/* (This won't catch the problem if it occurs in a malloc outside rev) */

/* When a malloc fails, reduce the maximum cache to */
/* it's current allocation minus the given size. */
static void rev_reduce_cache(size_t size) {
	rev_struct *rsi;
	size_t ram;

	/* Compute how much ram is currently allocated */
	for (ram = 0, rsi = g_rev_instances; rsi != NULL; rsi = rsi->next)
		ram += rsi->sz;

	if (size > ram)
		error("rev_reduce_cache: run out of rev virtual memory! (want %d, got %d)",size,ram);

//printf("~1 size = %d, g_test_ram = %d\n",size,g_test_ram);
//printf("~1 rev: Reducing cache because alloc of %d bytes failed. Reduced from %d to %d MB\n",
//size, g_avail_ram/1000000, (ram - size)/1000000);
	ram = g_avail_ram = ram - size;

	/* Aportion the memory, and reduce the cache allocation to match */
	ram /= g_no_rev_cache_instances; 
	for (rsi = g_rev_instances; rsi != NULL; rsi = rsi->next) {
		revcache *rc = rsi->cache;

		rsi->max_sz = ram;
		while (rc->nunlocked > 0 && rsi->sz > rsi->max_sz) {
			if (decrease_revcache(rc) == 0)
				break;
		}
//printf("~1 rev instance ram = %d MB\n",rsi->sz/1000000);
	}
	if (g_rev_instances != NULL && g_rev_instances->sb->s->verbose)
		printf("%cThere %s %d rev cache instance%s with %lu Mbytes limit\n",
              cr_char,
				g_no_rev_cache_instances > 1 ? "are" : "is",
                   g_no_rev_cache_instances,
				g_no_rev_cache_instances > 1 ? "s" : "",
                   (unsigned long)ram/1000000);
}

/* Check that the requested allocation plus 20 M Bytes */
/* can be allocated, and if not, reduce the rev-cache limit. */
/* This is so as to detect running out of VM before */
/* we actually run out and (on OS X) avoid emitting a warning. */
static void rev_test_vram(size_t size) {
	char *a1;
#ifdef __APPLE__
	int old_stderr, new_stderr;

	/* OS X malloc() blabs about a malloc failure. This */
	/* will confuse users, so we temporarily redirect stdout */
	fflush(stderr);
	old_stderr = dup(fileno(stderr));
	new_stderr = open("/dev/null", O_WRONLY | O_APPEND);
	dup2(new_stderr, fileno(stderr));
#endif
	size += 20 * 1024 * 1024;	/* This depends on the VM region allocation size */
	if ((a1 = malloc(size)) == NULL) {
		rev_reduce_cache(size);
	} else {
		free(a1);
	}
	g_test_ram = size/2;		/* Allow for twice as much VM to be used for each allocation */
#ifdef __APPLE__
	fflush(stderr);
	dup2(old_stderr, fileno(stderr));	/* Restore stderr */
	close(new_stderr);
	close(old_stderr);
#endif
}

static void *rev_malloc(rspl *s, size_t size) {
	void *rv;

	if ((size + 1 * 1024 * 1024) > g_test_ram)
		rev_test_vram(size);
	if ((rv = malloc(size)) == NULL) {
		rev_reduce_cache(size);
		rv = malloc(size);
	}
	if (rv != NULL)
		g_test_ram -= size;

	return rv;
}

static void *rev_calloc(rspl *s, size_t num, size_t size) {
	void *rv;

	if (((num * size) + 1 * 1024 * 1024) > g_test_ram)
		rev_test_vram(size);
	if ((rv = calloc(num, size)) == NULL) {
		rev_reduce_cache(num * size);
		rv = calloc(num, size);
	}
	if (rv != NULL)
		g_test_ram -= size;

	return rv;
}

static void *rev_realloc(rspl *s, void *ptr, size_t size) {
	void *rv;

	if ((size + 1 * 1024 * 1024) > g_test_ram)
		rev_test_vram(size);
	if ((rv = realloc(ptr, size)) == NULL) {
		rev_reduce_cache(size);		/* approximation */
		rv = realloc(ptr, size);
	}
	if (rv != NULL)
		g_test_ram -= size;

	return rv;
}


/* ====================================================== */
/* Set the ink limit information for any reverse interpolation. */
/* Calling this will clear the reverse interpolaton cache and acceleration structures. */
static void
rev_set_limit_rspl(
	rspl *s,		/* this */
	double (*limit)(void *lcntx, double *in),	/* Optional input space limit function. Function */
					/* should evaluate in[0..di-1], and return number that is not to exceed */
					/* limitv. NULL if not used */
	void *lcntx,	/* Context passed to limit() */
	double limitv	/* Value that limit() is not to exceed */
) {
	schbase *b;

	DBG(("rev: setting ink limit function %p and limit %f\n",limit,limitv));

	/* This is a restricted size function */
	if (s->di > MXRI)
		error("rspl: rev_set_limit can't handle di = %d",s->di);
	if (s->fdi > MXRO)
		error("rspl: rev_set_limit can't handle fdi = %d",s->fdi);

	b = set_search_limit(s, limit, lcntx, limitv);	/* Init and set limit info */

	if (s->rev.inited) {		/* If cache and acceleration has been allocated */
		invalidate_revaccell(s);		/* Invalidate the reverse cache */
	}

	/* Invalidate any ink limit values cached with the fwd grid data */
	clear_limitv(s);
}

/* Get the ink limit information for any reverse interpolation. */
static void
rev_get_limit_rspl(
	rspl *s,		/* this */
	double (**limitf)(void *lcntx, double *in),	/* Return pointer to function of NULL if not set */
	void **lcntx,	/* return context pointer */
	double *limitv	/* Return limit value */
) {
	schbase *b = s->rev.sb;

	/* This is a restricted size function */
	if (s->di > MXRI)
		error("rspl: rev_get_limit can't handle di = %d",s->di);
	if (s->fdi > MXRO)
		error("rspl: rev_get_limit can't handle fdi = %d",s->fdi);

	if (b == NULL) {
		*limitf = NULL;
		*lcntx = NULL;
		*limitv = 0.0;
	} else {
		*limitf = s->limitf;
		*lcntx = s->lcntx;
		*limitv = s->limitv/INKSCALE;
	}
}

/* Set the RSPL_NEARCLIP LCh weightings. */
/* Will only work with L*a*b* like output spaces. */
/* Calling this will clear the reverse interpolaton cache. */
static void rev_set_lchw(
	struct _rspl *s,	/* this */
	double lchw[MXRO]	/* Weighting */
) {
	int f;

	DBG(("rev: setting LCH weightings %f %f %f \n",lchw[0], lchw[1], lchw[2]));

	/* This is a restricted size function */
	if (s->di > MXRI)
		error("rspl: rev_set_lchw can't handle di = %d",s->di);
	if (s->fdi > MXRO || s->fdi != 3)
		error("rspl: rev_set_lchw can't handle fdi = %d",s->fdi);

	s->rev.lchweighted = 1;
	for (f = 0; f < s->fdi; f++) {
		s->rev.lchw[f] = lchw[f];
		s->rev.lchw_sq[f] = s->rev.lchw[f] * s->rev.lchw[f];
	}
	s->rev.lchw_chsq = s->rev.lchw_sq[1] - s->rev.lchw_sq[2];	/* C - H squared weight */

	if (s->rev.inited) {				/* If cache and acceleration has been allocated */
		invalidate_revaccell(s);		/* Invalidate the reverse cache */
	}
}

#define RSPL_CERTAIN 0x80000000 						/* WILLCLIP hint is certain */
#define RSPL_WILLCLIP2 (RSPL_CERTAIN | RSPL_WILLCLIP)	/* Clipping will certainly be needed */

#ifdef CHECK_NNLU
static void check_nn(rspl *s, double *oval, co *cpp);
static void print_nnck(rspl *s);
#endif

/* Do reverse interpolation given target output values and (optional) auxiliary target */
/* input values. Return number of results and clipping flag. If return value == mxsoln, */
/* then there might be more results. The target values returned will correspond to the */
/* actual (posssibly clipped) point. The return value is the number of solutions + */
/* a clipped flag. Properly set hint flags improve performance, but a correct result should */
/* be returned if the RSPL_NEARCLIP is set, even if they are not set correctly. */
/* If RSPL_NONNSETUP is set, then rev.fastsetup will be set for this call, avoiding */
/* initialization of the nngrid if RSPL_NEARCLIP hasn't been used before. */ 
static int
rev_interp_rspl(
	rspl *s,		/* this */
	int flags,		/* Hint flag */
	int mxsoln,		/* Maximum number of solutions allowed for */
	int *auxm,		/* Array of di mask flags, !=0 for valid auxliaries (NULL if no auxiliaries) */
	double cdir[MXRO],	/* Clip vector direction and length - NULL if not used */
	co *cpp			/* Given target output space value in cpp[0].v[] +  */
					/* target input space auxiliaries in cpp[0].p[], return */
					/* input space solutions in cpp[0..retval-1].p[], and */
) {
	int e, di = s->di;
	int fdi = s->fdi;
	int i, *rip = NULL;
	schbase *b = NULL;		/* Base search information */
	double auxv[MXRI];		/* Locus proportional auxiliary values */
	int didclip = 0;		/* flag - set if we clipped the target */
	int fastsetup;			/* fastsetup on entry */
	
	DBGV(("\nrev interp called with out targets", fdi, " %f", cpp[0].v, "\n"));

	/* This is a restricted size function */
	if (di > MXRI)
		error("rspl: rev_interp can't handle di = %d",di);
	if (fdi > MXRO)
		error("rspl: rev_interp can't handle fdi = %d",fdi);

	if (auxm != NULL) {
		double ax[MXRI];
		for (i = 0; i < di; i++) {
			if (auxm[i] != 0)
				ax[i] = cpp[0].p[i];
			else
				ax[i] = 0.0;
		}
		DBGV(("                  auxiliaries mask", di, " %d", auxm, "\n"));
		DBGV(("                auxiliaries values", di, " %f", ax, "\n"));
	}
	DBG(("di = %d, fdi = %d\n",di, fdi));
	DBG(("flags = 0x%x\n",flags));

	fastsetup = s->rev.fastsetup;	/* fastsetup on entry */
	if (flags & RSPL_NONNSETUP)		/* Avoid triggering nnsetup on this call */
		s->rev.fastsetup = 1;

	mxsoln &= RSPL_NOSOLNS;			/* Prevent silliness */

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Auxiliary is proportion of locus, so we need to find locus extent */	
	if (flags & RSPL_AUXLOCUS) {
		DBG(("rev interp has aux targets as proportion of locus\n"));

		flags &= ~RSPL_WILLCLIP;		/* Reset hint flag, as we will figure it out */

		/* For each valid auxiliary */
		for (e = 0; e < di; e++) {
			if (auxm[e] == 0)
				continue;			/* Skip unsused auxiliaries */
	
			/* Do search for min and max */
			DBG(("rev locus searching for aux %d min/max\n", e));
			if (b == NULL) {
				b = init_search(s, flags, cpp[0].p, auxm, cpp[0].v, cdir, cpp, mxsoln, locus);
#ifdef STATS
				s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
			} else
				set_lsearch(s, e);		/* Reset locus search for next auxiliary */

			if (rip == NULL) {		/* Not done this yet */
				rip = calc_fwd_cell_list(s, cpp[0].v); /* Reverse grid index for out target */
				if (rip == NULL) {
					DBG(("Got NULL list (point outside range) for auxiliary locus search\n"));
					flags |= RSPL_WILLCLIP2;
					break;
				}
			}
	
			search_list(b, rip, s->get_next_touch(s)); /* Setup, sort and search the list */
	
			if (b->min > b->max) {			/* Failed to find locus */
				DBG(("rev interp failed to find locus for aux %d, so expect clip\n",e));
				flags |= RSPL_WILLCLIP2;
				break;
			}
			auxv[e] = (cpp[0].p[e] * (b->max - b->min)) + b->min;
		}

		DBG(("rev interp got all locuses, so expect exact result\n",e));
		if (!(flags & RSPL_WILLCLIP)) {
			flags |= RSPL_EXACTAUX;				/* Got locuses, so expect exact result */
		}
	}

	/* Init the search information */
	if (b == NULL)
		b = init_search(s, flags, cpp[0].p, auxm, cpp[0].v, cdir, cpp, mxsoln, exact);
	else
		adjust_search(s, flags, auxv, exact);		/* Using proportion of locus aux */

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If hinted that we will not need to clip, look for exact solution. */
	if (!(flags & RSPL_WILLCLIP)) {
		DBG(("Hint we won't clip, so trying exact search\n"));

		/* First do an exact search (init will select auxil if requested) */
		adjust_search(s, flags, NULL, exact);
	
		/* Figure out the reverse grid index appropriate for this request */
		if (rip == NULL)	/* Not done this yet */
			rip = calc_fwd_cell_list(s, cpp[0].v);
	
#ifdef STATS
			s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
		if (rip != NULL) {
			/* Setup, sort and search the list */
			search_list(b, rip, s->get_next_touch(s));
		} else {
			DBG(("Got NULL list (point outside range) for first exact fxcell\n"));
		}
	
		/* If we selected exact aux, but failed to find a solution, relax expectation */
		if (b->nsoln == 0 && b->naux > 0 && (flags & RSPL_EXACTAUX)) {
//printf("~1 relaxing notclip expactation when nsoln == %d, naux = %d, falgs & RSPL_EXACTAUX = 0x%x\n", b->nsoln,b->naux,flags & RSPL_EXACTAUX);
			DBG(("Searching for exact match to auxiliary target failed, so try again\n"));
			adjust_search(s, flags & ~RSPL_EXACTAUX, NULL, exact);

#ifdef STATS
			s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
			/* Candidate cell list should be the same */
			if (rip != NULL) {
				/* Setup, sort and search the list */
				search_list(b, rip, s->get_next_touch(s));
			} else {
				DBG(("Got NULL list (point outside range) for nearest search fxcell\n"));
			}
		}
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If the exact search failed, and we should look for a nearest solution */
	if (b->nsoln == 0 && (flags & RSPL_NEARCLIP)) {
#ifdef CHECK_NNLU
		int f, fdi = s->fdi;
		double oval[MXRO];			/* Save the input target value for check_nn() */

		for (f = 0; f < fdi; f++)
			oval[f] = cpp[0].v[f];
#endif
		DBG(("Trying nearest search\n"));

#ifdef STATS
		s->rev.st[b->op].searchcalls++;
#endif	/* STATS */

		/* We get returned a list of cube base indexes of all cubes that have */
		/* the closest valid vertex value to the target value. */
		adjust_search(s, flags, NULL, clipn);

		/* Get list of cells enclosing nearest vertex */
		if ((rip = calc_fwd_nn_cell_list(s, cpp[0].v)) != NULL) {
			search_list(b, rip, s->get_next_touch(s)); /* Setup, sort and search the list */
		} else {
			DBG(("Got NULL list! (point inside gamut \?\?) for nearest search\n"));
		}

		if (b->nsoln > 0) {
			didclip = RSPL_DIDCLIP;
#ifdef CHECK_NNLU
			check_nn(s, oval, cpp);		/* Run diagnostic to check sanity of result */
#endif
		}
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If we still don't have a solution, do a vector direction clip */
	if (b->nsoln == 0 && b->canvecclip) {
		/* Find clipping solution in vector direction */
		line ln;				/* Structure to hold line context */
		unsigned int tcount;	/* grid touch count for this operation */

		DBG(("Starting a clipping vector search now!!\n"));

		adjust_search(s, flags, NULL, clipv);

		tcount = s->get_next_touch(s);		/* Get next grid touched generation count */

#ifdef STATS
		s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
		init_line_eq(b, b->v, cdir);				/* Init the implicit line equation */
		rip = init_line(s, &ln, cpp[0].v, cdir);	/* Init the line cell dda */
//~~1 HACK!!! should be <= 1.0 !!!
		for (; ln.t <= 2.0; rip = next_line_cell(&ln)) {
			if (rip == NULL) {
				DBG(("Got NULL list for this fxcell\n"));
				continue;
			}

			/* Setup, sort and search the list */
			search_list(b, rip, tcount);

			/* If we have found a solution, then abort the search - */
			/* this line will be taking us away from the best solution. */
			if (b->nsoln > 0)
				break;
		}
		if (b->nsoln > 0)
			didclip = RSPL_DIDCLIP;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If the clipped solution seems to have been jumping to conclusions, */
	/* search for an exact solution. */
	if (didclip && (flags & RSPL_WILLCLIP && !(flags & RSPL_CERTAIN))
	 && (b->cdist/s->get_out_scale(s)) < 0.002) {
		co c_cpp       = b->cpp[0];	/* Save clip solution in case we want it */
		double c_idist = b->idist;	
		int c_iabove   = b->iabove;	
		int c_nsoln    = b->nsoln;
		int c_pauxcell = b->pauxcell;
		double c_cdist = b->cdist;
		int c_iclip    = b->iclip;

		DBG(("Trying exact search again\n"));

		/* Do an exact search (init will select auxil if requested) */
		adjust_search(s, flags & ~RSPL_WILLCLIP, NULL, exact);
	
		/* Figure out the reverse grid index appropriate for this request */
		rip = calc_fwd_cell_list(s, cpp[0].v);
	
#ifdef STATS
		s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
		if (rip != NULL) {
			/* Setup, sort and search the list */
			search_list(b, rip, s->get_next_touch(s));
		} else {
			DBG(("Got NULL list (point outside range) for first exact fxcell\n"));
		}
	
		/* If we selected exact aux, but failed to find a solution, relax expectation */
		if (b->nsoln == 0 && b->naux > 0 && (flags & RSPL_EXACTAUX)) {
			DBG(("Searching for exact match to auxiliary target failed, so try again\n"));
//printf("~1 relaxing didclip expactation when nsoln == %d, naux = %d, flags & RSPL_EXACTAUX = 0x%x\n", b->nsoln,b->naux,flags & RSPL_EXACTAUX);
			adjust_search(s, flags & ~RSPL_EXACTAUX, NULL, exact);

#ifdef STATS
			s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
			/* Candidate cell list should be the same */
			if (rip != NULL) {
				/* Setup, sort and search the list */
				search_list(b, rip, s->get_next_touch(s));
			} else {
				DBG(("Got NULL list (point outside range) for nearest search fxcell\n"));
			}
		}

		/* If we did get an exact solution */
		if (b->nsoln > 0) {
			DBG(("Deciding to return exact solution after finding clipped\n"));
			didclip = 0;		/* Reset did-clip and return exact solution */

		} else {
			DBG(("keeping clipped solution\n"));
			/* Restore the clipped solution */
			b->cpp[0] = c_cpp;
			b->idist = c_idist;	
			b->iabove = c_iabove;	
			b->nsoln = c_nsoln;
			b->pauxcell = c_pauxcell;
			b->cdist = c_cdist;
			b->iclip = c_iclip;
		}
	}

	if (b->nsoln > 0) {
		DBGV(("rev interp returning 1st soln: ",di," %f", cpp[0].p, "\n"));
	}
	DBG(("rev interp returning %d solutions%s\n",b->nsoln, didclip ? " [clip]" : ""));

	s->rev.fastsetup = fastsetup;	/* retore fastsetup state */

	return b->nsoln | didclip;
}

/* ------------------------------------------------------------------------------------ */
/* Do reverse search for the auxiliary min/max ranges of the solution locus for the */
/* given target output values. */
/* Return number of locus segments found, up to mxsoln. 0 will be returned if no solutions */
/* are found. */

static int
rev_locus_segs_rspl (
	rspl *s,		/* this */
	int *auxm,		/* Array of di mask flags, !=0 for valid auxliaries (NULL if no auxiliaries) */
	co *cpp,		/* Input value in cpp[0].v[] */
	int mxsoln,		/* Maximum number of solutions allowed for */
	double min[][MXRI],	/* Array of min[MXRI] to hold return segment minimum values. */
	double max[][MXRI]	/* Array of max[MXRI] to hold return segment maximum values. */
) {
	int e, di = s->di;
	int f, fdi = s->fdi;
	int six;		/* solution index */
	int *rip = NULL;
	int rv = 1;				/* Return value */
	schbase *b = NULL;		/* Base search information */
	
	DBGV(("rev locus called with out targets", fdi, " %f", cpp[0].v, "\n"));
	
	/* This is a restricted size function */
	if (di > MXRI)
		error("rspl: rev_locus_segs can't handle di = %d",di);
	if (fdi > MXRO)
		error("rspl: rev_locus_segs can't handle fdi = %d",fdi);

	if (mxsoln < 1) {
		return 0;			/* Guard against silliness */
	}

	if (auxm != NULL) {
		int i;
		double ax[MXRI];
		for (i = 0; i < di; i++) {
			if (auxm[i] != 0)
				ax[i] = cpp[0].p[i];
			else
				ax[i] = 0.0;
		}
		DBGV(("                  auxiliaries mask", di, " %d", auxm, "\n"));
		DBGV(("                auxiliaries values", di, " %f", ax, "\n"));
	}

	/* Init default return values */
	for (six = 0; six < mxsoln; six++) {
		for (e = 0; e < di; e++) {
			if (auxm[e] == 0) {
				min[six][e] = max[six][e] = 0;	/* Return 0 for unused auxiliaries */
			} else {
				min[six][e] = 1.0;			/* max < min indicates invalid range */
				max[six][e] = 0.0;
			}
		}
	}

	/* For each valid auxiliary */
	for (e = 0; e < di; e++) {
		if (auxm[e] == 0)
			continue;			/* Skip unsused auxiliaries */

		/* Do search for min and max */
		DBG(("rev locus searching for aux %d min/max\n", e));
		if (b == NULL)
			b = init_search(s, 0, cpp[0].p, auxm, cpp[0].v, NULL, cpp, mxsoln, locus);
		else
			set_lsearch(s, e);		/* Reset locus search for next auxiliary */

		if (rip == NULL) {		/* Not done this yet */
			rip = calc_fwd_cell_list(s, cpp[0].v); /* Reverse grid index for this request */
			if (rip == NULL) {
				DBG(("Got NULL list (point outside range) for auxiliary locus search\n"));
				rv = 0;
				break;
			}
		}

		search_list(b, rip, s->get_next_touch(s)); /* Setup, sort and search the list */

		if (b->min > b->max) {
			rv = 0;				/* Failed to find a result */
			break;
		}

		if (b->asegs == 0) {		/* Overall min max only */

			min[0][e] = b->min;		/* Save single result */
			max[0][e] = b->max;

		} else {				/* Tracking auxiliary segments */
			int si;					/* Start i */
			int i, j, ff;

			/* Sort the segment list */
#define 	HEAP_COMPARE(A,B) (A.xval < B.xval)
			HEAPSORT(axisec, b->axisl, b->axisln)
#undef 		HEAP_COMPARE

#ifdef NEVER
for (i = 0; i < b->axisln; i++) {
printf("~2 xval = %f, verts = ",b->axisl[i].xval);
for (f = 0; f < b->axisl[i].nv; f++)
printf(" %d", b->axisl[i].vix[f]);
printf("\n");
}
#endif
			/* Find the segments by finding common verticies */
			six = si = i = 0;

			min[six][e] = b->axisl[i].xval;

			for (i++; i < (b->axisln-1); i++) {
				/* Check if any i and i-1 to j are connected */
				for (j = i-1; j >= si; j--) {
					for (f = 0; f < b->axisl[j].nv; f++) {
						for (ff = 0; ff < b->axisl[i].nv; ff++) {
							if (b->axisl[j].vix[f] == b->axisl[i].vix[ff])
								break;		/* Found a link */
						}
						if (ff < b->axisl[i].nv)
							break;
					}
					if (f < b->axisl[j].nv)
						break;
				}
				if (j < si) {	/* Wasn't linked */
					int ii, jj;
					/* Think we found a break. Check that all the rest of */
					/* the entries don't have any links to the previous group */

					/* This could be rather a slow way of checking ! (On^2) */
					for (ii = i+1; ii < (b->axisln); ii++) {
						for (jj = i-1; jj >= si; jj--) {
							for (f = 0; f < b->axisl[jj].nv; f++) {
								for (ff = 0; ff < b->axisl[ii].nv; ff++) {
									if (b->axisl[jj].vix[f] == b->axisl[ii].vix[ff])
										break;		/* Found a link */
								}
								if (ff < b->axisl[ii].nv)
									break;
							}
							if (f < b->axisl[jj].nv)
								break;
						}
						if (jj >= si)
							break;
					}
					if (ii >= b->axisln) {	/* Wasn't forward linked */
						/* Nothing ahead links to last group */
						max[six][e] = b->axisl[i-1].xval;

						/* If we run out of solution space */
						/* merge the last segments */
						if ((six+1) < mxsoln) {
							six++;
							min[six][e] = b->axisl[i].xval;
						}
					}
				}
			}
			max[six++][e] = b->axisl[i].xval;

			if (six > rv)
				rv = six;
		}
	}

#ifdef STATS
	s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
	if (rv) {
		for (six = 0; six < rv; six++) {
			DBG(("rev locus returning:\n"));
			DBGV(("     min", di, " %f", min[six], "\n"));
			DBGV(("     max", di, " %f", max[six], "\n"));
		}
	}

	DBG(("rev locus returning status %d\n",rv));
	return rv;
}

/* ------------------------------------------------------------------------------------ */
typedef double mxdi_ary[MXRI];

/* Do reverse search for the locus of the auxiliary input values given a target output. */
/* Return 1 on finding a valid solution, and 0 if no solutions are found. */
static int
rev_locus_rspl(
	rspl *s,		/* this */
	int *auxm,		/* Array of di mask flags, !=0 for valid auxliaries (NULL if no auxiliaries) */
	co *cpp,		/* Input value in cpp[0].v[] */
	double min[MXRI],/* Return minimum auxiliary values */
	double max[MXRI] /* Return maximum auxiliary values */
) {

	/* Use segment routine to compute oveall locus */
	return rev_locus_segs_rspl (s, auxm, cpp, 1, (mxdi_ary *)min, (mxdi_ary *)max);
}

/* ------------------------------------------------------------------------------------ */

#ifdef DEBUG2
#define DEBUG
#undef DBG
#undef DBGV
#undef DBGM
#define DBG(xxx) DBGI(xxx)
#define DBGV(xxx) DBGVI xxx
#define DBGM(xxx) DBGMI xxx
#else
#undef DEBUG
#undef DBG
#undef DBGV
#undef DBGM
#define DBG(xxx) 
#define DBGV(xxx) 
#define DBGM(xxx) 
#endif

/* ------------------------------------------------ */
/* subroutines of top level reverse lookup routine */

static int exact_setsort(schbase *b, fxcell *c);
static int exact_compute(schbase *b, simplex *x);

static int auxil_setsort(schbase *b, fxcell *c);
static int auxil_check(schbase *b, fxcell *c);
static int auxil_compute(schbase *b, simplex *x);

static int locus_setsort(schbase *b, fxcell *c);
static int locus_check(schbase *b, fxcell *c);
static int locus_compute(schbase *b, simplex *x);

static int clipv_setsort(schbase *b, fxcell *c);
static int clipv_check(schbase *b, fxcell *c);
static int clipv_compute(schbase *b, simplex *x);

static int clipn_setsort(schbase *b, fxcell *c);
static int clipn_check(schbase *b, fxcell *c);
static int clipn_compute(schbase *b, simplex *x);

/* Allocate the search base structure */
static schbase *
alloc_sb(rspl *s) {
	schbase *b;
	if ((b = s->rev.sb = (schbase *)rev_calloc(s, 1, sizeof(schbase))) == NULL)
		error("rspl malloc failed - rev.sb structure");
	INCSZ(s, sizeof(schbase));

	b->s = s;				/* rsp */
	b->pauxcell =			/* Previous solution cell indexes */
	b->plmaxcell = 
	b->plmincell = -1;

	return b;
}

/* Free the search base structure */
static void
free_sb(schbase *b) {
	DECSZ(b->s, sizeof(schbase));
	free(b);
}

/* Do the basic search type independent initialization */
static schbase *	/* Return pointer to base search information */
init_search(
	rspl *s,		/* rsp; */
	int flags,		/* Hint flag */

	double *av,		/* Auxiliary input values - may be NULL */
	int *auxm,		/* Array of di mask flags, !=0 for valid auxliaries (NULL if no auxiliaries) */
					/* Locus search will search for max/min of first valid auxlilary */
	double *v,		/* Output value target, NULL if none */
	double *cdir,	/* Clip vector direction/LCh weighting, NULL if none */
	co *cpp,		/* Array that hold solutions, NULL if none. */
	int mxsoln,		/* Maximum number of solutions allowed for */
	enum ops op		/* Type of reverse search operation requested */
) {
	schbase *b = NULL;		/* Pointer to search base information structure */
	int e, di = s->di;
	int f, fdi = s->fdi;

	DBG(("Initializing search di %d fdi %d\n",s->di,s->fdi));

	if (s->rev.inited == 0) 	/* Compute reverse info if it doesn't exist */
		make_rev(s);

	/* If first time initialisation (Fourth section init) */
	if ((b = s->rev.sb) == NULL)
		b = alloc_sb(s);

	/* Init some basic search info */
	b->op    = op;				/* operation */
	b->flags = flags;			/* hint flags */
	b->canvecclip = 0;			/* Assume invalid clip direction */

	b->ixc = (1<<di)-1;			/* Cube index of corner that holds maximum input values */

	/* Figure out if auxiliaries have been requested */
	b->naux = 0;
	b->auxbm = 0;
	if (auxm != NULL) {
		unsigned bm;

		if (mxsoln > 1)
			b->asegs = 1;		/* Find all segments */
		else
			b->asegs = 0;		/* Find only overall aux locus range */

		for (e = di-1, bm = 1 << e; e >= 0; e--, bm >>= 1) {	/* Record auxiliary mask bits */
			if (av != NULL)
				b->av[e] = av[e];	/* Auxiliary target values */
			b->auxm[e] = auxm[e];	/* Auxiliary mask */
			if (auxm[e] != 0) {
				b->auxbm |= bm;			/* Auxiliary bit mask */
				b->auxi[b->naux++] = e;	/* Index of next auxiliary input to be used */
				/* Auxiliary locus extent */
				b->lxi = e;			/* Assume first one */
				b->max = -INF_DIST;	/* In case searching for max */
				b->min =  INF_DIST;	/* In case searching for minimum */
				b->axisln = 0;		/* No intersects in list */
			}
		}
	}

	/* Figure out if the clip direction is meaningfull */
	/* Check that the clip vector makes sense */
	if (!(flags & RSPL_NEARCLIP)  && cdir != NULL) {	/* Clip vector is specified */
		double ss;
		for (ss = 0.0, f = 0; f < fdi; f++) {
			double tt = cdir[f];
			b->cdir[f] = tt;
			ss += tt * tt;
		}

		if (ss > 1e-6) {
			b->canvecclip = 1;	/* It has a non-zero length */
			ss = sqrt(ss);
			/* Compute normalised clip vector direction */
			for (f = 0; f < fdi; f++) {
				b->ncdir[f] = b->cdir[f]/ss;
			}
		}
	}

	if (di <= fdi)		/* Only allow auxiliaries if di > fdi */
		b->naux = 0;

	/* Switch to appropriate operation */
	if (b->op == exact && (b->naux > 0 || di != fdi)) {
		b->op = auxil;
	} else if (b->op == auxil && b->naux == 0 && di == fdi) {
		b->op = exact;
	}

	/* Set appropriate functions for type of operation */
	switch (b->op) {
		case exact:
			b->setsort = exact_setsort;
			b->check   = NULL;
			b->compute = exact_compute;
			b->snsdi = b->ensdi = di;	/* Search full dimension simplex, expect point soln. */
			break;
		case auxil:
			b->setsort = auxil_setsort;
			b->check   = auxil_check;
			b->compute = auxil_compute;
			b->snsdi = di;				/* Start here DOF = di-fdi locus solutions */
			b->ensdi = fdi;				/* End with DOF = 0 for point solutions */
			break;
		case locus:
			b->setsort = locus_setsort;
			b->check   = locus_check;
			b->compute = locus_compute;
			b->snsdi = b->ensdi = fdi;	/* Search for point solutions */
			break;
		case clipv:
			b->setsort = clipv_setsort;
			b->check   = clipv_check;
			b->compute = clipv_compute;
											/* Clip vector 1 dimension in output space, */
			b->snsdi = b->ensdi = fdi-1;	/* search planes for combined point solution */
			break;
		case clipn:
			b->setsort = clipn_setsort;
			b->check   = clipn_check;
			b->compute = clipn_compute;
			b->snsdi = 0;				/* Start with DOF = 0 for point solutions */
			b->ensdi = fdi-1;			/* End on DOF = di-fdi-1 on surfaces of simplexes */
			break;
		default:
			error("init_search: Unknown operation %d\n",b->op);
	}

	if (v != NULL) {
		for (f = 0; f < fdi; f++)	/* Record target output values */
			b->v[f] = v[f];
		b->v[fdi] = s->limitv;		/* Limitvalue is output target for limit clip subsimplexes */
	}

	b->mxsoln = mxsoln;				/* Allow solutions to be returned */
	b->cpp    = cpp;				/* Put solutions here */
	b->nsoln = 0;					/* No solutions at present */
	b->iclip = 0;					/* Default solution isn't above ink limit */

	if (flags & RSPL_EXACTAUX)		/* Expect to be able to match auxiliary target exactly */
		b->idist = 2.0 * EPS;		/* Best input distance to beat - helps sort/triage */
	else
		b->idist = INF_DIST;		/* Best input distance to beat. */
	b->iabove = 0;					/* Best isn't known to be above (yet) */

	b->cdist = INF_DIST;			/* Best clip distance to beat. */

	DBG(("Search initialized\n"));

	return b;
}

/* Adjust the search */
static void
adjust_search(
	rspl *s,		/* rsp; */
	int flags,		/* Hint flag */
	double *av,		/* Auxiliary input values - may be NULL */
	enum ops op		/* Type of reverse search operation requested */
) {
	schbase *b = s->rev.sb;		/* Pointer to search base information structure */
	int e, di = s->di;
	int fdi = s->fdi;

	DBG(("Adjusting search\n"));

	b->op    = op;				/* operation */
	b->flags = flags;			/* hint flags */

	/* Switch from exact to aux if we need to */
	if (b->op == exact && (b->naux > 0 || di != fdi)) {
		b->op = auxil;
	} else if (b->op == auxil && b->naux == 0 && di == fdi) {
		b->op = exact;
	}

	/* Update auxiliary target values */
	if (av != NULL) {
		for (e = 0; e < b->naux; e++) {
			int ee = b->auxi[e];
			b->av[ee] = av[ee];
		}
	}

	/* Set appropriate functions for type of operation */
	switch (b->op) {
		case exact:
			b->setsort = exact_setsort;
			b->check   = NULL;
			b->compute = exact_compute;
			b->snsdi = b->ensdi = di;		/* Expect point solution */
			break;
		case auxil:
			b->setsort = auxil_setsort;
			b->check   = auxil_check;
			b->compute = auxil_compute;
			b->snsdi = di;				/* Start here DOF = di-fdi locus solutions */
			b->ensdi = fdi;				/* End with DOF = 0 for point solutions, */
			break;						/* will early exit DOF if good soln found. */
		case locus:
			b->setsort = locus_setsort;
			b->check   = locus_check;
			b->compute = locus_compute;
			b->snsdi = b->ensdi = fdi;	/* Search for point solutions */
			break;
		case clipv:
			b->setsort = clipv_setsort;
			b->check   = clipv_check;
			b->compute = clipv_compute;
											/* Clip vector 1 dimension in output space, */
			b->snsdi = b->ensdi = fdi-1;	/* so the intersection with the simplex is a point. */
			break;
		case clipn:
			b->setsort = clipn_setsort;
			b->check   = clipn_check;
			b->compute = clipn_compute;
			b->snsdi = 0;				/* Start with DOF = 0 for point solutions */
			b->ensdi = fdi-1;			/* End on DOF = di-fdi-1 on surfaces of simplexes */
			break;						/* Will go through all DOF */
		default:
			error("init_search: Unknown operation %d\n",b->op);
	}

	b->nsoln = 0;					/* No solutions at present */

	if (flags & RSPL_EXACTAUX)		/* Expect to be able to match auxiliary target exactly */
		b->idist = 2.0 * EPS;		/* Best input distance to beat - helps sort/triage */
	else
		b->idist = INF_DIST;		/* Best input distance to beat. */
	b->iabove = 0;					/* Best isn't known to be above (yet) */

	b->cdist = INF_DIST;			/* Best clip distance to beat. */

	DBG(("Search adjusted\n"));
}

/* Adjust existing locus search for a different auxiliary */
static void
set_lsearch(
rspl *s,
int e			/* Next auxiliary */
) {
	schbase *b = s->rev.sb;		/* Pointer to search base information structure */

	b->lxi = e;			/* Assume first one */
	b->max = -INF_DIST;	/* In case searching for max */
	b->min =  INF_DIST;	/* In case searching for minimum */
	b->axisln = 0;		/* No intersects in list */
}

/* Set the limit search information */
/* Note this doesn't create or init the main rev information. */
static schbase *	/* Return pointer to base search information */
set_search_limit(
	rspl *s,		/* rsp; */
	double (*limitf)(void *vcntx, double *in),	/* Optional input space limit function. Function */
					/* should evaluate in[0..di-1], and return number that is not to exceed */
					/* limitv. NULL if not used */
	void *lcntx,	/* Context passed to limit() */
	double limitv	/* Value that limitf() is not to exceed */
) {
	schbase *b = NULL;		/* Pointer to search base information structure */

	/* If sb info needs initialising (Fourth section init) */
	if ((b = s->rev.sb) == NULL) {
		b = alloc_sb(s);
	}

	s->limitf = limitf;				/* Input limit function */
	s->lcntx  = lcntx; 				/* Context passed to limit() */
	s->limitv = INKSCALE * limitv; 	/* Context passed to values not to be exceeded by limit() */
	if (limitf != NULL) {
		s->limiten = 1;				/* enable limiting by default */
	} else
		s->limiten = 0;				/* No limit function, so limiting not enabled. */

	return b;
}

/* Free any search specific data, plus the search base. */
static void
free_search(
schbase *b	/* Base search information */
) {
	DBG(("Freeing search\n"));

	/* Clip line implicit equation (incuding space for ink target) */
	if (b->cla != NULL) {
		int fdi = b->s->fdi;
		free_dmatrix(b->cla, 0, fdi-1, 0, fdi);
		b->cla = NULL;
	}

	/* Auxiliary segment list */
	if (b->axislz > 0) {
		free(b->axisl);
		DECSZ(b->s, b->axislz * sizeof(axisec));
		b->axisl = NULL;
		b->axislz = 0;
		b->axisln = 0;
	}

	/* Sorted cell list */
	if (b->lclistz > 0) {
		free(b->lclist);
		DECSZ(b->s, b->lclistz * sizeof(fxcell *));
		b->lclist = NULL;
		b->lclistz = 0;
	}

	/* Simplex filter list */
	if (b->lsxfilt > 0) {
		free(b->sxfilt);
		DECSZ(b->s, b->lsxfilt * sizeof(char));
		b->sxfilt = NULL;
		b->lsxfilt = 0;
	}

	free_sb(b);
}

/* Return the pointer to the list of fwd cells given */
/* the target output values. The pointer will be to the first */
/* index in the list (ie. list address + 3) */
/* Return NULL if none in list (out of gamut). */
static int *
calc_fwd_cell_list(
	rspl *s,		/* this */
	double *v		/* Output values */
) {
	int f, fdi = s->fdi;
	int **rpp;
	int rgres_1 = s->rev.res - 1;

	if (s->rev.rev_valid == 0)
		init_revaccell(s);
		
	for (rpp = s->rev.rev, f = 0; f < fdi; f++) {
		int mi;
		double t = (v[f] - s->rev.gl[f])/s->rev.gw[f];
		mi = (int)floor(t);				/* Grid coordinate */
		if (mi < 0 || mi > rgres_1) { 	/* If outside valid reverse range */
			return NULL;
		}
		rpp += mi * s->rev.coi[f];	/* Accumulate reverse grid pointer */
	}
	s->rev.sb->rix = rpp - s->rev.rev;	/* Set diagnostic value */

	if (*rpp == NULL)
		return NULL;
	return (*rpp) + 3;
}

void alloc_simplexes(fxcell *c, int nsdi);

/* Given a pointer to a list of fwd cells, cull cells that */
/* cannot contain or improve the solution, sort the list, */
/* and then compute the final best solution. */
static void
search_list(
schbase *b,				/* Base search information */
int     *rip,			/* Pointer to first index in cell list */
unsigned int tcount		/* grid touch count for this operation */
) {
	rspl *s = b->s;
	int nsdi;
	int i;
	int nilist;				/* Number in cell list */
	unsigned int stouch;	/* Simplex touch count */
	
	DBG(("search_list called\n"));

	/* (rip[-3] contains allocation for fwd cells in the list) */
	/* (rip[-2] contains the index of the next free entry in the list) */
	/* (rip[-1] contains the reference count for the list) */
	if (b->lclistz < rip[-3]) {	/* Allocate more space if needed */

		if (b->lclistz > 0) {	/* Free old space before allocating new */
			free(b->lclist);
			DECSZ(b->s, b->lclistz * sizeof(fxcell *));
		}
		b->lclistz = 0;
		/* Allocate enough space for all the candidate cells */
		if ((b->lclist = (fxcell **)rev_malloc(s, rip[-3] * sizeof(fxcell *))) == NULL)
			error("rev: malloc failed - candidate cell list, count %d",rip[-3]);
		b->lclistz = rip[-3];	/* Current allocated space */
		INCSZ(b->s, b->lclistz * sizeof(fxcell *));
	}
		
	/* Get the next simplex touch count, so that we don't search shared */
	/* face simplexes more than once in this pass through the cells. */
	if ((stouch = ++s->rev.stouch) == 0) {	/* If touch count rolls over */
		fxcell *cp;
		stouch = s->rev.stouch = 1;

		/* For all of the cells */
		DBG(("touch has rolled over, resetting it\n"));
		for (cp = s->rev.cache->mrubot; cp != NULL; cp = cp->mruup) {
			int nsdi;
	
			if (cp->s == NULL)	/* Cell has never been used */
				continue;

			/* For all the simplexes in the fxcell */
			for (nsdi = 0; nsdi <= s->di; nsdi++) {
				if (cp->sx[nsdi] != NULL) {
					int si;

					for (si = 0; si < cp->sxno[nsdi]; si++) {
						cp->sx[nsdi][si]->touch = 0;
					}
				}
			}
		}
	}

	/* For each chunk of the list that we can fit in the rcache: */
	for (; *rip != -1;)  {

		/* Go through all the candidate fwd cells, and build up the list of search cells */
		for (nilist = 0; *rip != -1; rip++)  {
			int ix = *rip;				/* Fwd cell index */
			float *fcb = s->g.a + ix * s->g.pss;	/* Pointer to base float of fwd cell */
			fxcell *c;

			if (TOUCHF(fcb) >= tcount) {	/* If we have visited this cell before */
				DBG((" Already touched cell index %d\n",ix));
				continue;
			}
			/* Get pointers to cells from cache, and lock it in the cache */
			if ((c = get_fxcell(b, ix, nilist == 0 ? 1 : 0)) == NULL) {
				static int warned = 0;
				if (!warned) {
					warning("%cWarning - Reverse Cell Cache exausted, processing in chunks",cr_char);
					warned = 1;
				}
				DBG(("revcache is exausted, do search in chunks\n"));
				if (nilist == 0) {
					/* This should never happen, because nz force should prevent it */
					revcache *rc = s->rev.cache;
					fxcell *cp;
					int nunlk = 0;
					/* Double check that there are no unlocked cells */
					for (cp = rc->mrubot; cp != NULL && cp->refcount > 0; cp = cp->mruup) {
						if (cp->refcount == 0)
							nunlk++;
					}
					fprintf(stdout,"Diagnostic: rev.sz = %lu, rev.max_sz = %lu, numlocked = %d, nunlk = %d\n",
					               (unsigned long)rc->s->rev.sz, (unsigned long)rc->s->rev.max_sz,
					               rc->nunlocked, nunlk);
					error("Not enough memory to process in chunks");
				}
				break;		/* cache has run out of room, so abandon, and do it next time */
			}

			DBG(("checking out cell %d range %s\n",ix,pcellorange(c)));
			TOUCHF(fcb) = tcount;			/* Touch it */

			/* Check mandatory conditions, and compute search key */
			if (!b->setsort(b, c)) {
				DBG(("cell %d rejected from list\n",ix));
				unget_fxcell(s->rev.cache, c);
				continue;
			}
			DBG(("cell %d accepted into list\n",ix));

			b->lclist[nilist++] = c; /* Cell is accepted as recursion candidate */
		}

		if (nilist == 0) {
			DBG(("List was empty\n"));
		}

#ifdef DOSORT
		/* If appropriate, sort child cells into best order */
		/* == sort key smallest to largest */
		switch (b->op) {
			case locus:
				{	/* Special case, adjust sort values */
					double min = INF_DIST, max = -INF_DIST;
					for (i = 0; i < nilist; i++) {
						fxcell *c = b->lclist[i];
						if (c->sort < min)
							min = c->sort;
						if (c->sort > max)
							max = c->sort;
					}
					max = min + max;	/* Total of min/max */
					min = 0.5 * max;	/* Average sort value */
					for (i = 0; i < nilist; i++) {
						fxcell *c = b->lclist[i];
						if (c->ix == b->plmincell || c->ix == b->plmaxcell) {
							c->sort = -1.0;		/* Put previous solution cells at head of list */
						} else if (c->sort > min) {
							c->sort = max - c->sort;	/* Reflect about average */
						}
					}
				}
				/* Fall through to sort */
			case auxil:
			case clipv:
			case clipn:
#define 	HEAP_COMPARE(A,B) (A->sort < B->sort)
				HEAPSORT(fxcell *,b->lclist, nilist)
#undef 		HEAP_COMPARE
				break;
			default:
				break;
		}
#endif /* DOSORT */

		DBG(("List sorted, about to search\n"));
#ifdef NEVER
		printf("\n~1 Op = %s, Cell sort\n",opnames[b->op]);
		for (i = 0; i < nilist; i++) {
			printf("~1 List %d, cell %d, sort = %f\n",i,b->lclist[i]->ix,b->lclist[i]->sort);
		}
#endif /* NEVER */

		/* 
			Tried reversing the "for each cell" and "for each level" loops,
			but it made a negligible difference to the performance.
			We choose to have cell on the outer so that we can unlock
			them as we go, so that they may be freed, even though
			this is a couple of percent slower (?).
		 */

		/* For each cell in the list */
		for (i = 0; i < nilist; i++) {
			fxcell *c = b->lclist[i];

#ifdef STATS
			s->rev.st[b->op].csearched++;
#endif /* STATS */

			/* For each dimensionality of sub-simplexes, in given order */
			DBG(("Searching from level %d to level %d\n",b->snsdi, b->ensdi));
			for (nsdi = b->snsdi;;) {
				int j, nospx;					/* Number of simplexes in cell */

				DBG(("\n******************\n"));
				DBG(("Searching level %d\n",nsdi));

				/* For those searches that have an optimisation goal, */
				/* re-check the cell to see if the goal can still improve on. */
				if (b->check != NULL && !b->check(b, c))
					break;

				if (c->sx[nsdi] == NULL) {
					alloc_simplexes(c, nsdi);	/* Do level 1 initialisation for nsdi */
				}

				/* For each simplex in a cell */
				nospx = c->sxno[nsdi];			/* Number of nsdi simplexes */
				for (j = 0; j < nospx; j++) {
					simplex *x = c->sx[nsdi][j];

					if (x->touch >= stouch) {
						continue;						/* We've already seen this one */
					}

					if (s->limiten == 0) {
						if (x->flags & SPLX_CLIPSX)		/* If limiting is disabled, we're */
							continue;					/* not interested in clip plane simplexes */
					}
#ifdef STATS
					s->rev.st[b->op].ssearched++;
#endif /* STATS */
					if (b->compute(b, x)) {
						DBG(("search aborted by compute\n"));
						break;					/* Found enough solutions */
					}
					x->touch = stouch;			/* Don't look at it again */

				}	/* Next Simplex */

				if (nsdi == b->ensdi)
					break;						/* We're done with levels */

				/* Next Simplex dimensionality */
				if (b->ensdi < b->snsdi) {
					if (nsdi == b->snsdi && b->nsoln > 0
					 && (b->op != auxil || b->idist <= 2.0 * EPS))
						break; 		/* Don't continue though decreasing */
									/* sub-simplex dimensions if we found a solution at */
									/* the highest dimension level. */
					nsdi--;
				} else if (b->ensdi > b->snsdi) {
					nsdi++;				/* Continue through increasing sub-simplex dimenionality */
				}						/* until we get to the top. */
			}
			/* Unlock the fxcell now that we're done with it */
			unget_fxcell(s->rev.cache, b->lclist[i]);
		}	/* Next cell */

	}	/* Next chunk */

	DBG(("search_list complete\n"));
	return;
}

/* ------------------------------------- */
/* Vector search in output space support */

/* Setup the line, and fetch the first cell */
/* Return the pointer to the list of fwd cells, NULL if none in list. */
static int *
init_line(
	rspl *s,			/* this */
	line *l,			/* line structure */
	double st[MXRO],	/* start of line */
	double de[MXRO]		/* line direction and length */
) {
	int f, fdi = s->fdi;
	int **rpp;
	int rgres_1 = s->rev.res - 1;
	int nvalid = 0;		/* Flag set if outside reverse grid range */

	DBGV(("Line from ", fdi, " %f", st, "\n"));
	DBGV(("In dir    ", fdi, " %f", de, "\n"));
	DBGV(("gl        ", fdi, " %f", s->rev.gl, "\n"));
	DBGV(("gh        ", fdi, " %f", s->rev.gh, "\n"));
	DBGV(("gw        ", fdi, " %f", s->rev.gw, "\n"));
	
	/* Init */
	l->s = s;
	for (f = 0; f < fdi; f++) {
		l->st[f] = st[f] - s->rev.gl[f];
		l->de[f] = de[f];
		if (de[f] > 0.0)
			l->di[f] = 1;	/* Axis increments */
		else if (de[f] < 0.0)
			l->di[f] = -1;
		else
			l->di[f] = 0;
	}
	l->t = 0.0;
	DBGV(("increments =", fdi, " %d", l->di, "\n"));

	/* Figure out the starting cell */
	for (rpp = s->rev.rev, f = 0; f < fdi; f++) {
		double t = l->st[f]/s->rev.gw[f];
		l->ci[f] = (int)floor(t);					/* Grid coordinate */
		if (l->ci[f] < 0 || l->ci[f] > rgres_1) 	/* If outside valid reverse range */
			nvalid = 1;
		rpp += l->ci[f] * s->rev.coi[f];	/* Accumulate reverse grid pointer */
	}
	DBGV(("current line cell = ", fdi, " %d", l->ci, "")); DBG((",  t = %f, nvalid = %d\n",l->t,nvalid));
#ifdef DEBUG
	{
		int ii;
		double tt;
		printf("Current cell = ");
		for (ii = 0; ii < fdi; ii++) {
			tt = l->ci[ii] * s->rev.gw[ii] + s->rev.gl[ii];
			printf(" %f - %f",tt,tt+s->rev.gw[ii]);
		}
		printf("\n");
	}
#endif	/* DEBUG */
	if (nvalid)
		return NULL;
	if (*rpp == NULL)
		return NULL;
	return *rpp + 3;
}

/* Get the next cell on the line. */
/* Return the pointer to the list of fwd cells, NULL if none in list. */
static int *
next_line_cell(
	line *l		/* line structure */
) {
	rspl *s = l->s;
	int bf = 0, f, fdi = s->fdi;
	int **rpp;
	int rgres_1 = s->rev.res - 1;
	double bt = 100.0;	/* Best (smalest +ve) parameter value to move */

	/* See which axis cell crossing we will hit next */
	for (f = 0; f < fdi; f++) {
		double t;
		if (l->de[f] != 0) {
			t = ((l->ci[f] + l->di[f]) * s->rev.gw[f] - l->st[f])/l->de[f];
			DBG(("t for dim %d = %f\n",f,t));
			if (t < bt) {
				bt = t;
				bf = f;		/* Best direction to move */
			}
		}
	}

	/* Move to the next reverse grid coordinate */
	l->ci[bf] += l->di[bf];
	l->t = bt;

	DBGV(("current line cell =", fdi, " %d", l->ci, "")); DBG((",  t = %f\n",l->t));

#ifdef DEBUG
	{
		int ii;
		double tt;
		printf("Current cell = ");
		for (ii = 0; ii < fdi; ii++) {
			tt = l->ci[ii] * s->rev.gw[ii] + s->rev.gl[ii];
			printf(" %f - %f",tt,tt+s->rev.gw[ii]);
		}
		printf("\n");
	}
#endif	/* DEBUG */

	/* Compute fxcell index */
	for (rpp = s->rev.rev, f = 0; f < fdi; f++) {
		if (l->ci[f] < 0 || l->ci[f] > rgres_1) { 	/* If outside valid reverse range */
			DBG(("Outside list on dim %d, 0 <= %d <= %d\n", f, l->ci[f],rgres_1));
			return NULL;
		}
		rpp += l->ci[f] * s->rev.coi[f];	/* Accumulate reverse grid pointer */
	}
	if (*rpp == NULL)
		return NULL;
	return *rpp + 3;
}

/* ------------------------------------- */
/* Clip nearest support. */

/* Weighted distance function macro: */

#define LCHW_SQ(fname, arg2type) 		\
									\
static double fname(rspl *s, double in1[MXDO], arg2type in2[MXDO]) {	\
	int f, fdi = s->fdi;	\
	double tt, rr = 0.0;	\
	\
	/* Fall back */	\
	if (!s->rev.lchweighted || fdi < 3) {	\
		for (f = 0; f < fdi; f++) {	\
			tt = in1[f] - (double)in2[f];	\
			rr += tt * tt;	\
		}	\
		return rr;	\
	}	\
	\
	{	\
		double dxsq = 0.0, dchsq;	\
		double dlsq, dcsq, dhsq;	\
		double dc, c1, c2;	\
		\
		/* Compute delta L squared and delta E squared */	\
		{	\
			double dl, da, db;	\
			dl = in1[0] - (double)in2[0];	\
			da = in1[1] - (double)in2[1];	\
			db = in1[2] - (double)in2[2];	\
		\
			dlsq = dl * dl;		/* dl squared */	\
			dchsq = da * da + db * db;	\
		}	\
		\
		/* Add any extra dims */	\
		for (f = 3; f < fdi; f++) {	\
			tt = in1[f] - (double)in2[f];	\
			dxsq += tt * tt;	\
		}	\
		\
		/* compute delta chromanance squared */	\
		{	\
			/* Compute chromanance for the two colors */	\
			c1 = sqrt(in1[1] * in1[1] + in1[2] * in1[2]);	\
			c2 = sqrt((double)in2[1] * (double)in2[1] + (double)in2[2] * (double)in2[2]);	\
		\
			dc = c1 - c2;	\
			dcsq = dc * dc;	\
		}	\
		\
		/* Compute delta hue squared */	\
		/* (Hue is simply the orthogonal delta to chromanance in the a*b* plane) */	\
		if ((dhsq = dchsq - dcsq) < 0.0)	\
			dhsq = 0.0;	\
		\
		/* Compute weighted error squared */	\
		rr = dxsq + s->rev.lchw_sq[0] * dlsq + s->rev.lchw_sq[1] * dcsq + s->rev.lchw_sq[2] * dhsq;	\
		\
		return rr;	\
	}	\
}

/* Compute weighted LCh output distance squared. */
/* Weighting is to L,C,h, delta's squared - double[], double[] version */
LCHW_SQ(lchw_sq, double)

/* Weighting is to L,C,h, delta's squared - double[], float[] version */
LCHW_SQ(lchw_sq_f, float)

/*	Notes:

	Estimation accuracy is hobbled by 100% at HWEIGHT 1.0 
	compare to pure euclidean estimate, due to the conservative
	maxDlc maxDh of points in group, but this reduces at larger
	HWEIGHT's. The handicap also decreases quickly with tighter
	group size, since C variation is diminished.

	The handicap limits filtering efficiency for large group to group,
	so ideally group size shouldn't be larger than about 10 DE in diameter.

	It's not clear if any better approach is possible.
*/

#define       NN_GCMIN (1e-6)

/* Create a nn group. */
/* If G != NULL, use it as group center rather than computing from members. */
static void nn_grpinit(rspl *s, nn_grp *p, double **pnts, int npnts, double *G) {
	int f, ee, ff, fdi = s->fdi;
	int i;
	double *min[MXRO], *max[MXRO];	/* Pointers to points with min/max values */
	double rad, radsq = -1.0;		/* Span/radius squared */
	int spf;
	double dxsq = 0.0, desq, dchsq, dlcsq;	
	double dlsq, dcsq, dhsq;	
	double dc, c1, c2;	
	double c, minc = 1e200, maxc = -1.0;

	if (G != NULL) {
		for (f = 0; f < fdi; f++)
			p->bcent[f] = G[f];

		if (fdi >= 3) {
			/* Track minimum and maximum member C squared */
			for (i = 0; i < npnts; i++) {
				c = pnts[i][1] * pnts[i][1] + pnts[i][2] * pnts[i][2];
				if (c < minc)
					minc = c;
				if (c > maxc)
					maxc = c;
			}
		}

	} else if (npnts <= 2) {

		/* Compute center as simple average */
		for (f = 0; f < fdi; f++)
			p->bcent[f] = 0.0;
	
		for (i = 0; i < npnts; i++) {
			for (f = 0; f < fdi; f++)
				p->bcent[f] += pnts[i][f];
	
			if (fdi >= 3) {
				/* Track minimum and maximum member C squared */
				c = pnts[i][1] * pnts[i][1] + pnts[i][2] * pnts[i][2];
				if (c < minc)
					minc = c;
				if (c > maxc)
					maxc = c;
			}
		}
		for (f = 0; f < fdi; f++)
			p->bcent[f] *= 1.0/(double)npnts;

	} else {
		/* We establish a center point in un-weighted space, because this is */
		/* what's needed for in-gamut work, and is computationally faster */
		/* and easier than attempting it using weighted space. */
	
		/* Find verticies of cell that have min and max values in output space */
		for (f = 0; f < fdi; f++)
			min[f] = max[f] = NULL;
	
		for (ee = 0; ee < npnts; ee++) {
			double *vp = pnts[ee];
			for (f = 0; f < fdi; f++) {
				if (min[f] == NULL || min[f][f] > vp[f])
					min[f] = vp;
				if (max[f] == NULL || max[f][f] < vp[f])
					max[f] = vp;
			}
		}
	
		/* Find the pair of points with the largest span (diameter) in output space */
		for (ff = 0; ff < fdi; ff++) {
			double ss;
			for (ss = 0.0, f = 0; f < fdi; f++) {
				double tt;
				tt = max[ff][f] - min[ff][f];
				ss += tt * tt;
			}
			if (ss > radsq) {
				radsq = ss;
				spf = ff;		/* Output dimension max was in */
			}
		}
	
		/* Set initial bounding sphere */
		for (f = 0; f < fdi; f++)
			p->bcent[f] = (max[spf][f] + min[spf][f])/2.0;
		radsq /= 4.0;			/* diam^2 -> rad^2 */
		rad = sqrt(radsq);
		
		/* Go though all the points again, expanding sphere if necessary */
		for (ee = 0; ee < npnts; ee++) {
			double ss;
			double *vp = pnts[ee];
	
			/* Compute distance squared of point to bounding shere */
			for (ss = 0.0, f = 0; f < fdi; f++) {
				double tt = vp[f] - p->bcent[f];
				ss += tt * tt;
			}
			if (ss > radsq) {
				double tt;
				/* DBG(("Expanding bounding sphere by %f\n",sqrt(ss) - rad)); */
	
				ss = sqrt(ss) + EPS;			/* Radius to point */
				rad = (rad + ss)/2.0;
				radsq = rad * rad;
				tt = ss - rad;
				for (f = 0; f < fdi; f++)
					p->bcent[f] = (rad * p->bcent[f] + tt * vp[f])/ss;
			} else {
				/* DBG(("Bounding sphere encloses by %f\n",rad - sqrt(ss))); */
			}
		}
		if (fdi >= 3) {
			/* Establish the minimum and maximum member C squared */
			for (ee = 0; ee < npnts; ee++) {
				c = pnts[ee][1] * pnts[ee][1] + pnts[ee][2] * pnts[ee][2];
				if (c < minc)
					minc = c;
				if (c > maxc)
					maxc = c;
			}
		}
	}

	p->brad = p->bradsq = -1.0;
	p->maxDlc = -1.0;
	p->maxDh = p->maxDh_ = -1.0;
	p->sratio = 1.0;
	p->Wsratio = s->rev.lchw_sq[2]; 
	p->bratio = 1.0;
	p->Wbratio = s->rev.lchw_sq[2]; 
	p->Gc = p->Gc_ = NN_GCMIN;

	/* No weighting */
	if (!s->rev.lchweighted || fdi < 3) {

		for (i = 0; i < npnts; i++) {
			desq = 0.0;
			for (f = 0; f < fdi; f++) {
				double tt = p->bcent[f] - pnts[i][f];
				desq += tt * tt;
			}
			/* Track maximum euclidean  distance */
			if (desq > p->bradsq)
				p->bradsq = desq;
		}
		p->brad = sqrt(p->bradsq);	/* Distance rather than squared */

	/* Weighted */
	} else {
		double maxde = -1.0;

		/* Locate member maximum deltaLC and deltaH */
		for (i = 0; i < npnts; i++) {
	
			/* Compute delta L squared and delta E squared */
			{	
				double dl, dasq, dbsq;	
				dl = p->bcent[0] - pnts[i][0];	
				dlsq = dl * dl;		/* dl squared */	
				dasq = p->bcent[1] - pnts[i][1];	
				dasq *= dasq;
				dbsq = p->bcent[2] - pnts[i][2];	
				dbsq *= dbsq;
			
				dchsq = dasq + dbsq;	
				desq = dlsq + dchsq;	
			}	

			/* Add any extra dims */
			for (f = 3; f < fdi; f++) {
				double tt = p->bcent[f] - pnts[i][f];
				dxsq += tt * tt;
			}
			desq += dxsq;

			/* Track maximum euclidean distance too */
			if (desq > p->bradsq)
				p->bradsq = desq;
		
			/* compute delta chromanance squared */	
			{	
				/* Compute chromanance of member to group center */
				c1 = sqrt(p->bcent[1] * p->bcent[1] + p->bcent[2] * p->bcent[2]);	
				c2 = sqrt(pnts[i][1] * pnts[i][1] + pnts[i][2] * pnts[i][2]);	
			
				dc = c1 - c2;	
				dcsq = dc * dc;	
			}	
			
			/* Compute delta hue squared */	
			/* (Hue is simply the orthogonal delta to chromanance in the a*b* plane) */	
			if ((dhsq = dchsq - dcsq) < 0.0)	
				dhsq = 0.0;	
		
			/* Weighted delta extra + luminance + chromanance squared */
			dlcsq = dxsq + s->rev.lchw_sq[0] * dlsq + s->rev.lchw_sq[1] * dcsq;
	
			/* Using maxDlc & maxDh is an absolute worst case, but */
			/* using a more exact approximation to the worst point */
			/* for a given hue correction factor, doesn't seem to help */
			/* for HWEIGHT > 1.5 */ 

			/* Track maximum weighted deltaLC squared */
			if (dlcsq > p->maxDlc)
				p->maxDlc = dlcsq;
	
			/* Track maximum deltaH squared */
			if (dhsq > p->maxDh)
				p->maxDh = dhsq;
		}
		p->brad = sqrt(p->bradsq);	/* Euclidean distance rather than squared */
		p->maxDh_ = sqrt(p->maxDh);
	
		/* Pre-calculate center C squared */
		p->Gc = p->bcent[1] * p->bcent[1] + p->bcent[2] * p->bcent[2];
		if (p->Gc < NN_GCMIN)
			p->Gc = NN_GCMIN;
		p->Gc_ = sqrt(p->Gc);
	
		/* Calculate hue scale down factor for Group center to smallest member C */
		/* (This is used to scale point/center to center distance) */
		if (minc < p->Gc) {
			p->sratio = sqrt(minc/p->Gc);
			if (s->rev.lchw_sq[2] > 1.0)			/* Slightly improves filter ratio */
				p->Wsratio = (s->rev.lchw_sq[2] - 1.0) * p->sratio + 1.0;
			else
				p->Wsratio = s->rev.lchw_sq[2] * p->sratio;
		}
	
		/* Calculate hue scale up factor for Group center to largest member C */
		/* (This is used to scale point/center to center distance) */
		/* (For group target, multiply group ->bratio values ??) */
		if (maxc > p->Gc) {
			p->bratio = sqrt(maxc/p->Gc);
			if (s->rev.lchw_sq[2] > 1.0)			/* Slightly improves filter ratio */
				p->Wbratio = (s->rev.lchw_sq[2] - 1.0) * p->bratio + 1.0;
			else
				p->Wbratio = s->rev.lchw_sq[2] * p->bratio;
		}
	}
} 

/* Return nz if point is within euclidean bounding sphere. */
/* Also return distance squared in *dist if non-NULL */
static int nn_insphere(rspl *s, double *dist, nn_grp *p, double *src) {
	int f, fdi = s->fdi;
	double desq = 0.0;

	for (f = 0; f < fdi; f++) {
		double tt = p->bcent[f] - src[f];	
		desq += tt * tt;
	}

	if (dist != NULL)
		*dist = desq;

	return desq <= p->bradsq;
}

/* Estimate possible smallest weighted distance of point to group. */
/* If lgst != NULL, also return the estimated largest possible distance. */
static double nn_pntgrp_est(rspl *s, double *lgst, nn_grp *p, double *src) {
	int f, fdi = s->fdi;
	double dxsq = 0.0, desq, dchsq;
	double dlsq, dcsq, dhsq;
	double dc, c1, c2;	
	double Tc;			/* Point chromanance squared */
	double sGrr;		/* Min Point to group center diatance squared */
	double bGrr;		/* Max Point to group center diatance squared */
	double rr;			/* Largest member distance squared */
	double sdist;		/* Min. estimated distance squared */
	double bdist;		/* Max.. estimated distance squared */
	double aratio = 1.0;

	/* If not using LCh weighted distances */
	if (!s->rev.lchweighted || fdi < 3) {

		desq = 0.0;
		for (f = 0; f < fdi; f++) {
			double tt = p->bcent[f] - src[f];	
			desq += tt * tt;
		}

		/* Return largest possible distance */
		if (lgst != NULL) {
			bdist = sqrt(desq) + p->brad + EPS;
			*lgst = bdist;
		}

		/* Return min possible distance */
		sdist = sqrt(desq) - p->brad - EPS;
		if (sdist < 0.0)
			sdist = 0.0;
		return sdist;

	/* We're using LCh weighting, so we need to do some adjustments */
	} else {
		/* Compute components of weighted distance of point */
		/* to group center. */
		{	
			double dl, dasq, dbsq;	
			dl = p->bcent[0] - src[0];	
			dlsq = dl * dl;				/* dl squared */	
			dasq = p->bcent[1] - src[1];	
			dasq *= dasq;
			dbsq = p->bcent[2] - src[2];	
			dbsq *= dbsq;

			dchsq = dasq + dbsq;
		}	

		/* Compute any extra dims */
		for (f = 3; f < fdi; f++) {
			double tt = p->bcent[f] - src[f];
			dxsq += tt * tt;
		}

		/* compute delta chromanance squared of target to group center */	
		{	
			/* Compute delta chromanance between target point and group center */
			c1 = p->Gc_;
			c2 = Tc = src[1] * src[1] + src[2] * src[2];	
			c2 = sqrt(c2);
			dc = c1 - c2;	
			dcsq = dc * dc;	
		}	
		
		/* Compute delta hue squared of target point to group center */	
		/* (Hue is simply the orthogonal delta to chromanance in the a*b* plane) */	
		if ((dhsq = dchsq - dcsq) < 0.0)	
			dhsq = 0.0;	
	
		/* Weighted values of L and C delta's */
		dlsq *= s->rev.lchw_sq[0];
		dcsq *= s->rev.lchw_sq[1];

		/* Most distant member hue delta adjustment factor */
		aratio = s->rev.lchw_sq[2];
		if (Tc > p->Gc) {
			aratio = sqrt(Tc/p->Gc);
			if (s->rev.lchw_sq[2] > 1.0)			/* Slightly improves filter ratio */
				aratio = (s->rev.lchw_sq[2] - 1.0) * aratio + 1.0;
			else
				aratio = s->rev.lchw_sq[2] * aratio;
		}

		/* Adjusted maximum member distance to group center */
		rr = sqrt(p->maxDlc + aratio * p->maxDh);

		/* Return max. possible distance squared */
		if (lgst != NULL) {

			/* Adjusted weighted max. distance squared of target to group center */
			bGrr = dxsq + dlsq + dcsq + dhsq * p->Wbratio;

			/* max. possible distance of target to most distant member */
			bdist = sqrt(bGrr) + rr + EPS;
			*lgst = bdist;
		}

		/* Adjusted weighted min. distance squared of target to group center */
		sGrr = dxsq + dlsq + dcsq + dhsq * p->Wsratio;

		/* min. possible distance of target to most distant member */
		sdist = sqrt(sGrr) - rr - EPS;
		if (sdist < 0.0)
			sdist = 0.0;

		return sdist;
	}
}

/* Estimate possible smallest weighted distance of group to group. */
/* If lgst != NULL, also return the estimated largest possible distance. */
static double nn_grpgrp_est(rspl *s, double *lgst, nn_grp *p1, nn_grp *p2) {
	int f, fdi = s->fdi;
	double dxsq = 0.0, desq, dchsq;
	double dlsq, dcsq, dhsq;
	double dc, c1, c2;	
	double sGrr;		/* Min Point to group center diatance squared */
	double bGrr;		/* Max Point to group center diatance squared */
	double rr1, rr2;		/* Largest member distance squared */
	double sdist;			/* Min. estimated distance squared */
	double bdist;			/* Max.. estimated distance squared */
	double aratio1 = 1.0, aratio2 = 1.0;

	/* If not using LCh weighted distances */
	if (!s->rev.lchweighted || fdi < 3) {

		desq = 0.0;
		for (f = 0; f < fdi; f++) {
			double tt = p1->bcent[f] - p2->bcent[f];	
			desq += tt * tt;
		}

		/* Return largest possible distance */
		if (lgst != NULL) {
			bdist = sqrt(desq) + p1->brad + p2->brad + EPS;
			*lgst = bdist;
		}

		/* Return min possible distance */
		sdist = sqrt(desq) - p1->brad - p2->brad - EPS;
		if (sdist < 0.0)
			sdist = 0.0;
		return sdist;

	/* We're using LCh weighting, so we need to do some adjustments */
	} else {
		double Wratio;

		/* Compute components of weighted distance of group center */
		/* to group center. */
		{	
			double dl, dasq, dbsq;	
			dl = p1->bcent[0] - p2->bcent[0];	
			dlsq = dl * dl;				/* dl squared */	
			dasq = p1->bcent[1] - p2->bcent[1];	
			dasq *= dasq;
			dbsq = p1->bcent[2] - p2->bcent[2];	
			dbsq *= dbsq;
			
			dchsq = dasq + dbsq;
		}	

		/* Compute any extra dims */
		for (f = 3; f < fdi; f++) {
			double tt = p1->bcent[f] - p2->bcent[f];
			dxsq += tt * tt;
		}

		/* compute delta chromanance squared of point to group center */	
		{	
			/* Compute delta chromanance group centers */
			c1 = p1->Gc_;
			c2 = p2->Gc_;
			dc = c1 - c2;	
			dcsq = dc * dc;	
		}	
		
		/* Compute delta hue squared of group centers */	
		/* (Hue is simply the orthogonal delta to chromanance in the a*b* plane) */	
		if ((dhsq = dchsq -  dcsq) < 0.0)	
			dhsq = 0.0;	
	
		/* Weighted values of L and C delta's */
		dlsq *= s->rev.lchw_sq[0];
		dcsq *= s->rev.lchw_sq[1];

		/* Most distant member hue delta adjustment factor */
		aratio1 = aratio2 = s->rev.lchw_sq[2];

		if ((p1->Gc_ + p1->maxDh) > p2->Gc_) {
			aratio2 = (p1->Gc_ + p1->maxDh)/p2->Gc_;
			if (s->rev.lchw_sq[2] > 1.0)			/* Slightly improves filter ratio */
				aratio2 = (s->rev.lchw_sq[2] - 1.0) * aratio2 + 1.0;
			else
				aratio2 = s->rev.lchw_sq[2] * aratio2;

		}
		if ((p2->Gc_ + p2->maxDh) > p1->Gc_) {
			aratio1 = (p2->Gc_ + p2->maxDh)/p1->Gc_;
			if (s->rev.lchw_sq[2] > 1.0)			/* Slightly improves filter ratio */
				aratio1 = (s->rev.lchw_sq[2] - 1.0) * aratio1 + 1.0;
			else
				aratio1 = s->rev.lchw_sq[2] * aratio1;
		}

		/* Adjusted maximum member distance to group center */
		rr1 = sqrt(p1->maxDlc + aratio1 * p1->maxDh);
		rr2 = sqrt(p2->maxDlc + aratio2 * p2->maxDh);

		/* Returne max. possible distance squared */
		if (lgst != NULL) {
			if (s->rev.lchw_sq[2] > 1.0)			/* Slightly improves filter ratio */
				Wratio = (s->rev.lchw_sq[2] - 1.0) * p1->bratio * p2->bratio + 1.0;
			else
				Wratio = s->rev.lchw_sq[2] * p1->bratio * p2->bratio;

			/* Adjusted weighted max. distance squared of group centers */
			bGrr = dxsq + dlsq + dcsq + dhsq * Wratio;

			/* max. possible distance of target to most distant member */
			bdist = sqrt(bGrr) + rr1 + rr2 + EPS;
			*lgst = bdist;
		}

		if (s->rev.lchw_sq[2] > 1.0)			/* Slightly improves filter ratio */
			Wratio = (s->rev.lchw_sq[2] - 1.0) * p1->sratio * p2->sratio + 1.0;
		else
			Wratio = s->rev.lchw_sq[2] * p1->sratio * p2->sratio;

		/* Adjusted weighted min. distance squared of group centers */
		sGrr = dxsq + dlsq + dcsq + dhsq * Wratio;

		/* min. possible distance of target to most distant member */
		sdist = sqrt(sGrr) - rr1 - rr2 - EPS;
		if (sdist < 0.0)
			sdist = 0.0;

		return sdist;
	}
}

/* ------------------------------------------------------------ */
static void fill_nncell(rspl *s, int *co, int ix);

/* Return the pointer to the list of nearest fwd cells given */
/* the target output values. The pointer will be to the first */
/* index in the list (ie. list address + 3) */
/* Return NULL if none in list (out of gamut). */
static int *
calc_fwd_nn_cell_list(
	rspl *s,		/* this */
	double *v		/* Output values */
) {
	int f, fdi = s->fdi, ix;
	int **rpp;
	int rgres_1 = s->rev.res - 1;
	int mi[MXDO];

	if (s->rev.rev_valid == 0)
		init_revaccell(s);

	for (ix = 0, f = 0; f < fdi; f++) {
		double t = (v[f] - s->rev.gl[f])/s->rev.gw[f];
		mi[f] = (int)floor(t);			/* Grid coordinate */
		if (mi[f] < 0)	 				/* Clip to reverse range, so we always return a result  */
			mi[f] = 0;
		else if (mi[f] > rgres_1)
			mi[f] = rgres_1;
		ix += mi[f] * s->rev.coi[f];	/* Accumulate reverse grid index */
	}
	s->rev.sb->rix = ix;				/* Set diagnostic value */

	rpp = s->rev.nnrev + ix;
	if (*rpp == NULL) {
		if (s->rev.fastsetup)
			fill_nncell(s, mi, ix);		/* Fill on-demand */
		if (*rpp == NULL)
			rpp = s->rev.rev + ix;		/* fall back to in-gamut lookup */ 
	}
	if (*rpp == NULL) {
#ifdef CHECK_NNLU
		printf("Got NULL list for nearest search, targ %s,\n coord %s, rix %d\n", debPdv(fdi,v),debPiv(fdi,mi),ix);
		if (ix < 0 || ix >= s->rev.no)
			printf("Index is outside range 0 .. %d\n",s->rev.no-1);
		else {
			if (s->rev.nnrev[ix] == NULL)
				printf(" nnrev = NULL\n");
			else
				printf(" nnrev length = %d\n",s->rev.nnrev[ix][1]-3);
			if (s->rev.rev[ix] == NULL)
				printf(" rev = NULL\n");
			else
				printf(" rev = length = %d\n",s->rev.rev[ix][1]-3);
		}
#endif
		return NULL;
	}
	return (*rpp) + 3;
}

/* =================================================== */
/* The cell and simplex solver top level routines */

static int add_lu_svd(simplex *x);
static int add_locus(schbase *b, simplex *x);
static int add_auxil_lu_svd(schbase *b, simplex *x);
static int within_simplex(simplex *x, double *p);
static int within_simplex_limit(simplex *x, double *p);
static void simplex_to_abs(simplex *x, double *in, double *out);

static int auxil_solve(schbase *b, simplex *x, double *xp);

/* ---------------------- */
/* Exact search functions */
/* Return non-zero if cell is acceptable */
static int exact_setsort(schbase *b, fxcell *c) {
	rspl *s = b->s;
	int f, fdi = s->fdi;
	double ss;

	DBG(("Reverse exact search, evaluate and set sort key on cell\n"));

	/* Check that the target lies within the cell bounding sphere */
	for (ss = 0.0, f = 0; f < fdi; f++) {
		double tt = c->g.bcent[f] - b->v[f];
		ss += tt * tt;
	}
	if (ss > c->g.bradsq) {
		DBG(("Cell rejected - %s outside sphere c %s rad %f\n",debPdv(fdi,b->v),debPdv(fdi,c->g.bcent),sqrt(c->g.bradsq)));
		return 0;
	}

	if (s->limiten != 0 && c->limmin > s->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,s->limitv));
		return 0;
	}

	/* Sort can't be used, because we return all solutions */
	c->sort = 0.0;

	DBG(("Cell is accepted\n"));

	return 1;
}

/* Compute a solution for a given sub-simplex (if there is one) */
/* Return 1 if search should be aborted */
static int exact_compute(schbase *b, simplex *x) {
	rspl *s     = b->s;
	int e, di = s->di, sdi  = x->sdi;
	int f, fdi  = s->fdi;
	int i;
	datai xp;	/* solution in simplex relative coord order */
	datai p;	/* absolute solution */
	int wsrv;	/* Within simplex return value */

	DBG(("\nExact: computing possible solution\n"));

#ifdef DEBUG
	/* Sanity check */
	if (sdi != fdi || sdi != di || x->efdi != fdi) {
		printf("di = %d, fdi = %d\n",di,fdi);
		printf("sdi = %d, efdi = %d\n",sdi,x->efdi);
		error("rspl exact reverse interp called with sdi != fdi, sdi != di, efdi != fdi");
		/* !!! could switch to SVD solution if di != fdi ?? !!! */
	}
#endif

	/* This may not be worth it here since it may not filter out */
	/* many more simplexes than the cube check did. */
	/* This is due to full dimension simplexes all sharing the main */
	/* diagonal axis. */

	/* Check that the target lies within the simplex bounding cube */
	for (f = 0; f < fdi; f++) {
		if (b->v[f] < x->min[f] || b->v[f] > x->max[f]) {
			DBG(("Simplex is rejected - bounding cube\n"));
			return 0;
		}
	}

	/* Create the LU decomp needed to exactly solve */
	if (add_lu_svd(x)) {
		DBG(("LU decomp was singular, skip simplex\n"));
		return 0;
	}

	/* Init the RHS B[] vector (note di == fdi) */
	for (f = 0; f < fdi; f++) {
		xp[f] = b->v[f] - x->v[di][f];
	}

	/* Compute the solution (in simplex space) */
	lu_backsub(x->d_u, sdi, (int *)x->d_w, xp);

	/* Check that the solution is within the simplex & meets ink limit */
	if ((wsrv = within_simplex(x, xp)) == 0) {
		DBG(("Solution rejected because not in simplex\n"));
		return 0;
	}

	/* Convert solution from simplex relative to absolute space */
	simplex_to_abs(x, p, xp);

	/* Check if a very similiar input solution has been found before */
	for (i = 0; i < b->nsoln; i++) {
		double tt;
		for (e = 0; e < di; e++) {
			tt = b->cpp[i].p[e] - p[e];
			if (fabs(tt) > (2 * EPS))
				break;	/* Mismatch */
		}
		if (e >= di)	/* Found good match */
			break;
	}

	/* Probably alias caused by solution lying close to a simplex boundary */
	if (i < b->nsoln) {
		DBG(("Another solution has been found before - index %d\n",i));
		return 0;		/* Skip this, since betters been found before */
	}

	/* Check we haven't overflowed space */
	if (i >= b->mxsoln) {
		DBG(("Run out of space for new solution\n"));
		return 1;		/* Abort */
	}

	DBG(("######## Accepting new solution\n"));

	/* Put solution in place */
	for (e = 0; e < di; e++)
		b->cpp[i].p[e] = p[e];
	for (f = 0; f < fdi; f++)
		b->cpp[i].v[f] = b->v[f];	/* Assumed to be an exact solution */
	if (i == b->nsoln)
		b->nsoln++;
	if (wsrv == 2)					/* Is above (disabled) ink limit */
		b->iclip = 1;
	return 0;
}

/* -------------------------- */
/* Auxiliary search functions */
static int auxil_setsort(schbase *b, fxcell *c) {
	rspl *s = b->s;
	int f, fdi  = b->s->fdi;
	int ee, ixc = b->ixc;
	double ss, sort, nabove;

	DBG(("Reverse auxiliary search, evaluate and set sort key on cell\n"));

	if (b->s->di <= fdi) {	/* Assert */
		error("rspl auxiliary reverse interp called with di <= fdi (%d %d)", b->s->di, fdi);
	}

	/* Check that the target lies within the cell bounding sphere */
	for (ss = 0.0, f = 0; f < fdi; f++) {
		double tt = c->g.bcent[f] - b->v[f];
		ss += tt * tt;
	}
	if (ss > c->g.bradsq) {
		DBG(("Cell rejected - %s outside sphere c %s rad %f\n",debPdv(fdi,b->v),debPdv(fdi,c->g.bcent),sqrt(c->g.bradsq)));
		return 0;
	}

	if (s->limiten != 0 && c->limmin > s->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,s->limitv));
		return 0;
	}

	/* Check if this cell could possible improve b->idist */
	/* and compute sort key as the distance to auxilliary target */
	/* (We may have a non INF_DIST idist before commencing the */
	/* search if we already know that the auxiliary target is */
	/* within gamut - the usual usage case!) */
	for (sort = 0.0, nabove = ee = 0; ee < b->naux; ee++) {
		int ei = b->auxi[ee];
		double tt = (c->p[0][ei] + c->p[ixc][ei]) - b->av[ei];
		sort += tt * tt;
		if (c->p[ixc][ei] >= (b->av[ei] - EPS))		/* Could be above */
			nabove++;
	}

	if (b->flags & RSPL_MAXAUX && nabove < b->iabove) {
		DBG(("Doesn't contain solution that has as many aux above auxiliary goal\n"));
		return 0;
	}
	if (!(b->flags & RSPL_MAXAUX) || nabove == b->iabove) {
		for (ee = 0; ee < b->naux; ee++) {
			int ei = b->auxi[ee];
			if (c->p[0][ei]   >= (b->av[ei] + b->idist)
			 || c->p[ixc][ei] <= (b->av[ei] - b->idist)) {
				DBG(("Doesn't contain solution that will be closer to auxiliary goal\n"));
				return 0;
			}
		}
	}
	c->sort = sort + 0.01 * ss;

	if (c->ix == b->pauxcell)
		c->sort = -1.0;			/* Put previous calls solution cell at top of sort list */

	DBG(("Cell is accepted\n"));
	return 1;
}

/* Re-check whether it's worth searching cell */
static int auxil_check(schbase *b, fxcell *c) {
	int ee, ixc = b->ixc, nabove;

	DBG(("Reverse auxiliary search, re-check cell\n"));

	/* Check if this cell could possible improve b->idist */
	/* and compute sort key as the distance to auxilliary target */

	for (nabove = ee = 0; ee < b->naux; ee++) {
		int ei = b->auxi[ee];
		if (c->p[ixc][ei] >= (b->av[ei] - EPS))		/* Could be above */
			nabove++;
	}

	if (b->flags & RSPL_MAXAUX && nabove < b->iabove) {
		DBG(("Doesn't contain solution that has as many aux above auxiliary goal\n"));
		return 0;
	}
	if (!(b->flags & RSPL_MAXAUX) || nabove == b->iabove) {
		for (ee = 0; ee < b->naux; ee++) {
			int ei = b->auxi[ee];
			if (c->p[0][ei]   >= (b->av[ei] + b->idist)
			 || c->p[ixc][ei] <= (b->av[ei] - b->idist)) {
				DBG(("Doesn't contain solution that will be closer to auxiliary goal\n"));
				return 0;
			}
		}
	}
	DBG(("Cell is still ok\n"));
	return 1;
}

/* Compute a solution for a given simplex (if there is one) */
/* Return 1 if search should be aborted */
static int auxil_compute(schbase *b, simplex *x) {
	rspl *s     = b->s;
	int e, di   = s->di;
	int f, fdi  = s->fdi;
	datai xp;		/* solution in simplex relative coord order */
	datai p;		/* absolute solution */
	double idist;	/* Auxiliary input distance */
	int wsrv;		/* Within simplex return value */
	int nabove;		/* Number above aux target */

	DBG(("\nAuxil: computing possible solution\n"));

#ifdef DEBUG
	{
	unsigned int sum = 0;
	for (f = 0; f <= x->sdi; f++)
		sum += x->vix[f];
	printf("Simplex of cell ix %d, sum 0x%x, sdi = %d, efdi = %d\n",x->ix, sum, x->sdi, x->efdi);
	printf("Target val %s\n",debPdv(fdi, b->v));
	for (f = 0; f <= x->sdi; f++) {
		int ix = x->vix[f], i;
		float *fcb = s->g.a + ix * s->g.pss;	/* Pointer to base float of fwd cell */
		printf("Simplex vtx %d [cell ix %d] val %s\n",f,ix,debPfv(fdi, fcb));
	}
	}
#endif

	/* Check that the target lies within the simplex bounding cube */
	for (f = 0; f < fdi; f++) {
		if (b->v[f] < x->min[f] || b->v[f] > x->max[f]) {
			DBG(("Simplex is rejected - bounding cube\n"));
			return 0;
		}
	}

	/* Check if this cell could possible improve b->idist */
	for (nabove = e = 0; e < b->naux; e++) {
		int ei = b->auxi[e];					/* pmin/max[] is indexed in input space */
		if (x->pmax[ei] >= (b->av[ei] - EPS))	/* Could be above */
			nabove++;
	}
	if ((b->flags & RSPL_MAXAUX) && nabove < b->iabove) {
		DBG(("Simplex doesn't contain solution that has as many aux above auxiliary goal\n"));
		return 0;
	}
	if (!(b->flags & RSPL_MAXAUX) || nabove == b->iabove) {
		for (nabove = e = 0; e < b->naux; e++) {
			int ei = b->auxi[e];					/* pmin/max[] is indexed in input space */
			if (x->pmin[ei] >= (b->av[ei] + b->idist)
			 || x->pmax[ei] <= (b->av[ei] - b->idist)) {
				DBG(("Simplex doesn't contain solution that will be closer to auxiliary goal\n"));
				return 0;
			}
		}
	}

//printf("~~ About to create svd decomp\n");
	/* Create the SVD or LU decomp needed to compute solution or locus */
	if (add_lu_svd(x)) {
		DBG(("SVD decomp failed, skip simplex\n"));
		return 0;
	}

//printf("~~ About to solve locus for aux target\n");
	/* Now solve for locus parameter that minimises */
	/* distance to auxliary target. */
	if ((wsrv = auxil_solve(b, x, xp)) == 0) {
		DBG(("Target auxiliary along locus is outside simplex,\n"));
		DBG(("or computation failed, skip simplex\n"));
		return 0;
	}

//printf("~~ About to convert solution to absolute space\n");
	/* Convert solution from simplex relative to absolute space */
	simplex_to_abs(x, p, xp);

	DBG(("Got solution at %s\n", debPdv(di,p)));

//printf("~~ soln = %f %f %f %f\n",p[0],p[1],p[2],p[3]);
//printf("~~ About to compute auxil distance\n");
	/* Compute distance to auxiliary target */
	for (idist = 0.0, nabove = e = 0; e < b->naux; e++) {
		int ei = b->auxi[e];
		double tt = b->av[ei] - p[ei];
		idist += tt * tt;
		if (p[ei] >= (b->av[ei] - EPS))
			nabove++;
	}
	idist = sqrt(idist);
//printf("~1 idist %f, nabove %d\n",idist, nabove);
//printf("~1 best idist %f, best iabove %d\n",b->idist, b->iabove);

	/* We want the smallest error from auxiliary target */
	if (b->flags & RSPL_MAXAUX) {
		if (nabove < b->iabove || (nabove == b->iabove && idist >= b->idist)) {
			DBG(("nsoln %d, nabove %d, iabove %d, idist = %f, better solution has been found before\n",b->nsoln, nabove, b->iabove, idist));
			return 0;
		}
	} else {
		if (idist >= b->idist) {	/* Equal or worse auxiliary solution */
			DBG(("nsoln %d, idist = %f, better solution has been found before\n",b->nsoln,idist));
			return 0;
		}
	}

	/* Solution is accepted */
	DBG(("######## Accepting new solution with nabove %d <= iabove %d and idist %f <= %f\n",nabove,b->iabove,idist,b->idist));
	for (e = 0; e < di; e++)
		b->cpp[0].p[e] = p[e];
	for (f = 0; f < fdi; f++)
		b->cpp[0].v[f] = b->v[f];	/* Assumed to be an exact solution */
	b->idist = idist;
	b->iabove = nabove;
	b->nsoln = 1;
	b->pauxcell = x->ix;
	if (wsrv == 2)					/* Is above (disabled) ink limit */
		b->iclip = 1;

	return 0;
}

/* ------------------------------------ */
/* Locus range search functions */

static int locus_setsort(schbase *b, fxcell *c) {
	rspl *s = b->s;
	int f, fdi  = s->fdi;
	int lxi = b->lxi;	/* Auxiliary we are finding min/max of */
	int ixc = b->ixc;
	double sort, ss;

	DBG(("Reverse locus evaluate and set sort key on cell\n"));

#ifdef DEBUG
	if (b->s->di <= fdi) {	/* Assert ~1 */
		error("rspl auxiliary locus interp called with di <= fdi");
	}
#endif /* DEBUG */

	/* Check that the target lies within the cell bounding sphere */
	for (ss = 0.0, f = 0; f < fdi; f++) {
		double tt = c->g.bcent[f] - b->v[f];
		ss += tt * tt;
	}
	if (ss > c->g.bradsq) {
		DBG(("Cell rejected - %s outside sphere c %s rad %f\n",debPdv(fdi,b->v),debPdv(fdi,c->g.bcent),sqrt(c->g.bradsq)));
		return 0;
	}

	if (s->limiten != 0 && c->limmin > s->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,s->limitv));
		return 0;
	}

	/* Check if this cell could possible improve the locus min/max */
	if (b->asegs == 0) {	/* If we aren't find all segments of the locus */
		if (c->p[0][lxi] >= b->min && c->p[ixc][lxi] <= b->max ) {
			DBG(("Doesn't contain solution that will expand the locus\n"));
			return 0;
		}
	}

	/* Compute sort index from average of auxiliary values */
	sort = (c->p[0][b->lxi] + c->p[ixc][b->lxi]);
	
	c->sort = sort + 0.01 * ss;

	DBG(("Cell is accepted\n"));
	return 1;
}

/* Re-check whether it's worth searching simplexes */
static int locus_check(schbase *b, fxcell *c) {
	int lxi = b->lxi;	/* Auxiliary we are finding min/max of */
	int ixc = b->ixc;

	DBG(("Reverse locus re-check\n"));

	/* Check if this cell could possible improve the locus min/max */
	if (b->asegs == 0) {	/* If we aren't find all segments of the locus */
		if (c->p[0][lxi] >= b->min && c->p[ixc][lxi] <= b->max ) {
			DBG(("Doesn't contain solution that will expand the locus\n"));
			return 0;
		}
	}

	DBG(("Cell is still ok\n"));
	return 1;
}

static int auxil_locus(schbase *b, simplex *x);

/* We expect to be given a sub-simplex with no DOF, to give an exact solution */
static int locus_compute(schbase *b, simplex *x) {
	rspl *s  = b->s;
	int f, fdi  = s->fdi;
	int lxi  = b->lxi;	/* Auxiliary we are finding min/max of */

	DBG(("\nLocus: computing possible solution\n"));

#ifdef DEBUG
	{
	unsigned int sum = 0;
	for (f = 0; f <= x->sdi; f++)
		sum += x->vix[f];
	printf("Simplex of cell ix %d, sum 0x%x, sdi = %d, efdi = %d\n",x->ix, sum, x->sdi, x->efdi);
	printf("Target val %s\n",debPdv(fdi, b->v));
	for (f = 0; f <= x->sdi; f++) {
		int ix = x->vix[f], i;
		float *fcb = s->g.a + ix * s->g.pss;	/* Pointer to base float of fwd cell */
		double v[MXDO];
		printf("Simplex vtx %d [cell ix %d] val %s\n",f,ix,debPfv(fdi, fcb));
	}
	}
#endif

	/* Check that the target lies within the simplex bounding cube */
	for (f = 0; f < fdi; f++) {
		if (b->v[f] < x->min[f] || b->v[f] > x->max[f]) {
			DBG(("Simplex is rejected - bounding cube\n"));
			return 0;
		}
	}

	/* Check if simplex could possible improve the locus min/max */
	if (b->asegs == 0) {	/* If we aren't find all segments of the locus */
		if (x->pmin[lxi] >= b->min && x->pmax[lxi] <= b->max ) {
			DBG(("Simplex doesn't contain solution that will expand the locus\n"));
			return 0;
		}
	}

//printf("~~ About to create svd decomp\n");
	/* Create the SVD decomp needed to compute solution extreme points */
	if (add_lu_svd(x)) {
		DBG(("SVD decomp failed, skip simplex\n"));
		return 0;
	}

//printf("~~ About to solve locus for aux extremes\n");
	/* Now solve for locus parameter that are at the extremes */
	/* of the axiliary we are interested in. */
	if (!auxil_locus(b, x)) {
		DBG(("Target auxiliary is outside simplex,\n"));
		DBG(("or computation failed, skip simplex\n"));
		return 0;
	}

	return 0;
}

/* ------------------- */
/* Vector clipping search functions */
static int clipv_setsort(schbase *b, fxcell *c) {
	rspl *s = b->s;
	int f, fdi  = s->fdi;
	double ss, dp;

	DBG(("Reverse clipping search evaluate cell\n"));

//printf("~~sphere center = %f %f %f, radius %f\n",c->bcent[0],c->bcent[1],c->bcent[2],sqrt(c->bradsq));
	/* Check if the clipping line intersects the bounding sphere */
	/* First compute dot product cdir . (bcent - v) */
	/* == distance to center of sphere in direction of clip vector */
	for (dp = 0.0, f = 0; f < fdi; f++) {
		dp += b->ncdir[f] * (c->g.bcent[f] - b->v[f]);
	}

	if (s->limiten != 0 && c->limmin > s->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,s->limitv));
		return 0;
	}

//printf("~~ dot product = %f\n",dp);
	/* Now compute closest distance to sphere center */
	for (ss = 0.0, f = 0; f < fdi; f++) {
		double tt = b->v[f] + dp * b->ncdir[f] - c->g.bcent[f];
		ss += tt * tt;
	}

//printf("~~ distance to sphere center = %f\n",sqrt(ss));
	if (ss > c->g.bradsq) {
		DBG(("Cell is rejected - wrong direction or bounding sphere\n"));
		return 0;
	}
	c->sort = dp;		/* May be -ve if beyond clip target point ? */

	DBG(("Cell is accepted\n"));
	return 1;
}

/* Clipping check functions */
/* Note that we don't bother with this check in setsort(), */
/* because we assume that nothing will set a small cdist */
/* before the search commences (unlike auxil). */
/* Note that line search loop exits on finding any solution. */
static int clipv_check(schbase *b, fxcell *c) {

	DBG(("Reverse clipping re-check\n"));

	if (b->cdist < INF_DIST) {	/* If some clip solution has been found */
		int f, fdi = b->s->fdi;
		double dist;
		/* Compute a conservative "best possible solution clip distance" */
		for (dist = 0.0, f = 0; f < fdi ; f++) {
			double tt = (c->g.bcent[f] - b->v[f]);
			dist += tt * tt;
		}
		dist = sqrt(dist); /* Target distance to bounding */

		if (dist >= (c->g.brad + b->cdist)) {	/* Equal or worse clip solution */
			DBG(("Cell best possible solution worse than current\n"));
			return 0;
		}
	}

	DBG(("Cell is still ok\n"));
	return 1;
}

static int vnearest_clip_solve(schbase *b, simplex *x, double *xp, double *xv, double *err);

/* Compute a clip solution */
static int clipv_compute(schbase *b, simplex *x) {
	rspl   *s  = b->s;
	int f, fdi = s->fdi;
	datai p;				/* Input space solution */
	datao v;				/* Output space solution */
	double err;				/* output error of solution */
	int wsrv;	/* Within simplex return value */

	DBG(("Clips: computing possible solution\n"));

	/* Compute a solution value */
	if ((wsrv = vnearest_clip_solve(b, x, p, v, &err)) == 0) {
		DBG(("Doesn't contain a solution\n"));
		return 0;
	}

	/* We want the smallest clip error */
	/* (Should we reject points in -ve vector direction ??) */
	if (err >= b->cdist) {	/* Equal or worse clip solution */
		DBG(("better solution has been found before\n"));
		return 0;
	}

	simplex_to_abs(x, b->cpp[0].p, p);	/* Convert to abs. space & copy */

	DBG(("######## Accepting new clipv solution with error %f\n",err));
#ifdef DEBUG
	if (s->limiten != 0) {
		DBG(("######## Ink value = %f, limit %f\n",get_limitv(b, x->ix, NULL, b->cpp[0].p), s->limitv));
	}
#endif

	/* Put solution in place */
	for (f = 0; f < fdi; f++)
		b->cpp[0].v[f] = v[f];
	b->cdist = err;
	b->nsoln = 1;
	if (wsrv == 2)					/* Is above (disabled) ink limit */
		b->iclip = 1;

	return 0;
}

/* ------------------- */
/* Nearest clipping search functions. */
/* We use weighted distances if lchweighted. */
static int clipn_setsort(schbase *b, fxcell *c) {
	rspl *s = b->s;
	int f, fdi  = s->fdi;
	double ss;

	DBG(("Reverse nearest clipping search evaluate fwd cell ix %d\n",c->ix));
//if (b->rix == 7135) printf("Reverse nearest clipping search evaluate fwd cell ix %d\n",c->ix);

	/* Compute an estimated weighted clip distance from target point to this fxcell */
	ss = nn_pntgrp_est(s, NULL, &c->g, b->v);

	/* Check that the cell could possibly improve the solution */
	if (b->cdist < INF_DIST) {	/* If some clip solution has been found */
		if (ss >= b->cdist) {	/* Equal or worse clip solution */
			DBG(("Cell best possible solution worse than current\n"));

//if (b->rix == 7135) {
//	printf("Cell best possible solution worse than current\n");
//	printf("current dist %f, best to fwd %f\n",b->cdist,ss); 
//}
			return 0;
		}
	}

	if (s->limiten != 0 && c->limmin > s->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,s->limitv));
//if (b->rix == 7135) printf("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,s->limitv);
		return 0;
	}

	c->sort = ss;

	DBG(("Cell is accepted (%f < %f)\n",ss,b->cdist));
//if (b->rix == 7135) printf("Cell is accepted (%f < %f)\n",ss,b->cdist);
	return 1;
}

/* Clipping check functions */
static int clipn_check(schbase *b, fxcell *c) {

	DBG(("Reverse nearest clipping re-check fwd cell ix %d\n",c->ix));
//if (b->rix == 7135) printf("Reverse nearest clipping re-check fwd cell ix %d\n",c->ix);

	if (b->cdist < INF_DIST) {	/* If some clip solution has been found */
		/* re-use sort value, best possible distance to solution */
		if (c->sort >= b->cdist) {	/* Equal or worse clip solution */
			DBG(("Cell best possible solution now worse than current\n"));
//if (b->rix == 7135) {
//	printf("Cell best possible solution now worse than current\n");
//	printf("current dist %f, best to fwd %f\n",b->cdist,c->sort); 
//}
			return 0;
		}
	}

	DBG(("Cell is still ok\n"));
//if (b->rix == 7135) printf("Cell is still ok\n");
	return 1;
}

static int lchw_nnearest_clip_solve(schbase *b, simplex *x, double *xp, double *xv, double *err);
static int nnearest_clip_solve(schbase *b, simplex *x, double *xp, double *xv, double *err);

/* Compute a clip solution */
static int clipn_compute(schbase *b, simplex *x) {
	rspl   *s  = b->s;
	int f, fdi = s->fdi;
	datai p;			/* Simplex input space solution */
	datao v;			/* Output space solution */
	double err;			/* output error of solution */
	int wsrv;			/* Within simplex return value */

	DBG(("Clipn: computing possible solution cell %d, simplex %d, sdi = %d, efdi = %d\n",x->ix,x->si,x->sdi,x->efdi));
//if (b->rix == 7135) printf("Clipn: computing possible solution cell %d, simplex %d, sdi = %d, efdi = %d\n",x->ix,x->si,x->sdi,x->efdi);

	/* Compute a solution value */
	if (s->rev.lchweighted) {
		if ((wsrv = lchw_nnearest_clip_solve(b, x, p, v, &err)) == 0) {
			DBG(("Doesn't contain a solution\n"));
//if (b->rix == 7135) printf("Doesn't contain a solution\n");
			return 0;
		}
	} else {
		if ((wsrv = nnearest_clip_solve(b, x, p, v, &err)) == 0) {
			DBG(("Doesn't contain a solution\n"));
//if (b->rix == 7135) printf("Doesn't contain a solution\n");
			return 0;
		}
	}

	/* We want the smallest clip error */
	if (err >= b->cdist) {	/* Equal or worse clip solution */
		DBG(("better solution has been found before (%f < %f)\n",b->cdist,err));
//if (b->rix == 7135) printf("better solution has been found before (%f < %f)\n",b->cdist,err);
		return 0;
	}

	DBG(("######## Accepting new clipn solution with error %f (replaces %f)\n",err,b->cdist));
//if (b->rix == 7135) printf("######## Accepting new clipn solution with error %f (replaces %f)\n",err,b->cdist);

	simplex_to_abs(x, b->cpp[0].p, p);	/* Convert to abs. space & copy */

	/* Put solution in place */
	for (f = 0; f < fdi; f++)
		b->cpp[0].v[f] = v[f];
	b->cdist = err;
	b->nsoln = 1;
	if (wsrv == 2)					/* Is above (disabled) ink limit */
		b->iclip = 1;

	return 0;
}

/* -------------------------------------------------------- */
/* Cell/simplex solver middle level code */

/* Find the point on this sub-simplexes solution locus that is */
/* closest to the target auxiliary values, and return it in xp[] */
/* Return zero if this point canot be calculated, */
/* or it lies outside the simplex. */
/* Return 1 normally, and 2 if the solution would be over the ink limit */
static int
auxil_solve(
schbase *b,
simplex *x,
double *xp		/* Return solution xp[sdi] */
) {
	rspl *s = b->s;
	int ee, e, di = s->di, sdi = x->sdi; 
	int f, efdi = x->efdi; 
	int dof = sdi-efdi;			 /* Degree of freedom of simplex locus */
	int *icomb = x->psxi->icomb; /* abs -> simplex coordinate translation */
	double auxt[MXRI];			 /* Simplex relative auxiliary targets */
	double bb[MXRI];
	int wsrv;	/* Within simplex return value */

	DBG(("axuil_solve called\n"));

	if (dof < 0)
		error("Error - auxil_solve got sdi < efdi (%d < %d) - don't know how to handle this",sdi, efdi);

	/* If there is no locus, compute an exact solution */
	if (dof == 0) {
		DBG(("axuil_solve dof = zero\n"));

		/* Init the RHS B[] vector (note sdi == efdi) */
		for (f = 0; f < efdi; f++) {
			xp[f] = b->v[f] - x->v[sdi][f];
		}

		/* Compute the solution (in simplex space) */
		lu_backsub(x->d_u, sdi, (int *)x->d_w, xp);

		/* Check that the solution is within the simplex & meets ink limit */
		if ((wsrv = within_simplex(x, xp)) != 0) {
			DBG(("Got solution at %s\n", debPdv(sdi,xp)));
			return wsrv;				/* OK, got solution */
		}

		DBG(("No solution (not within simplex)\n"));
		return 0;
	}

	/* There is a locus, so find solution nearest auxiliaries */

	/* Compute locus for target function values (if sdi > efdi) */
	if (add_locus(b, x)) {
		DBG(("Locus computation failed, skip simplex\n"));
		return 0;
	}

	/* Convert aux targets from absolute space to simplex relative */
	for (e = 0; e < di; e++) {	/* For abs coords */
		int ei = icomb[e];		/* Simplex coord */

		if (ei >= 0 &&  b->auxm[e] != 0) {
			auxt[ei] = (b->av[e] - x->p0[e])/s->g.w[e];	/* Only sets those needed */
		}
	}

	if (dof == 1 && b->naux == 1) {		/* Special case, because it's common and easy! */
		int ei = icomb[b->auxi[0]];		/* Simplex relative auxiliary index */
		double tt;

		DBG(("axuil_solve dof = naux = 1\n"));
		if (ei < 0)
			return 0;					/* Not going to find solution */
		if ((tt = x->lo_l[ei][0]) == 0.0)
			return 0;
		tt = (auxt[ei] - x->lo_bd[ei])/tt;	/* Parameter solution for target auxiliary */

		/* Back substitute parameter */
		for (e = 0; e < sdi; e++) {
			xp[e] = x->lo_bd[e] + tt * x->lo_l[e][0];
		}
		/* Check that the solution is within the simplex & meets ink limit */
		if ((wsrv = within_simplex(x, xp)) != 0) {
			DBG(("Got solution %s\n",debPdv(di,xp)));
			return wsrv;				/* OK, got solution */
		}
		DBG(("No solution (not within simplex)\n"));
		return 0;
	}

	/* Compute the locus decompositions needed (info #5) */
	if (add_auxil_lu_svd(b, x)) {	/* Will set x->naux */
		DBG(("LU/SVD decomp failed\n"));
		return 0;
	}

	/* Setup B[], equation RHS  */
	for (e = ee = 0; ee < b->naux; ee++) {
		int ei = icomb[b->auxi[ee]];		/* Simplex relative auxiliary index */
		if (ei >= 0)						/* Usable auxiliary on this sub simplex */ 
			bb[e++] = auxt[ei] - x->lo_bd[ei];
	}
	if (e != x->naux)	/* Assert */
		error("Internal error - auxil_solve got mismatching number of auxiliaries");

	if (x->naux == dof) {			/* Use LU decomp to solve */
		DBG(("axuil_solve using LU\n"));
		lu_backsub(x->ax_u, dof, (int *)x->ax_w, bb);

	} else if (x->naux > 0) {	/* Use SVD to solve least squares */
		DBG(("axuil_solve using SVD\n"));
		svdbacksub(x->ax_u, x->ax_w, x->ax_v, bb, bb, x->naux, dof);

	} else {	/* x->naux == 0 */
		DBG(("axuil_solve  naux = 0\n"));
		for (f = 0; f < dof; f++)
			bb[f] = 0.0;		/* Use base solution ?? */
	}

	/* Now back substitute the locus parameters */
	/* to calculate the solution point (in simplex space) */
	for (e = 0; e < sdi; e++) {
		double tt;
		for (tt = 0.0, f = 0; f < dof; f++) {
			tt += bb[f] * x->lo_l[e][f];
		}
		xp[e] = x->lo_bd[e] + tt;
	}

	/* Check that the solution is within the simplex & meets ink limit */
	if ((wsrv = within_simplex(x, xp)) != 0) {
		DBG(("Got solution %s\n",debPdv(di,xp)));
		return wsrv;				/* OK, got solution */
	}
	DBG(("No solution (not within simplex)\n"));
	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Compute the min/max values for the current auxiliary of interest. */
/* Return zero if this point canot be calculated, */
/* or it lies outside the simplex. */
/* Return 1 normally, 2 if it would be outside the simplex if limting was enabled */
/* We expect to get a sub-simplex that will give an exact solution. */
static int
auxil_locus(
schbase *b,
simplex *x
) {
	rspl *s = b->s;
	int sdi = x->sdi; 
	int f, efdi = x->efdi; 
	double pp[MXRI];
	int wsrv;	/* Within simplex return value */

	DBG(("axuil_locus called\n"));

	if (sdi != efdi)
		warning("Internal error - auxil_locus got sdi != efdi (%d < %d)",sdi, efdi);

	/* Init the RHS B[] vector (note sdi == efdi) */
	for (f = 0; f < efdi; f++) {
		pp[f] = b->v[f] - x->v[sdi][f];
	}

	/* Compute the solution (in simplex space) */
	lu_backsub(x->d_u, sdi, (int *)x->d_w, pp);

	/* Check that the solution is within the simplex & meets ink limit */
	if ((wsrv = within_simplex(x, pp)) != 0) {
		double xval;
		int lxi = b->lxi;	/* Auxiliary we are finding min/max of (Abs space) */
		int xlxi = x->psxi->icomb[lxi];	/* Auxiliary we are finding min/max of (simplex space) */

		DBG(("Got locus solution within simplex\n"));

		/* Compute auxiliary value for this solution (absolute space) */
		xval = x->p0[lxi];
		if (xlxi >= 0)				/* Simplex param value */
			xval += s->g.w[lxi] * pp[xlxi];
		else if (xlxi == -2)		/* 1 value */
			xval += s->g.w[lxi];
									/* Else 0 value */

		if (b->asegs != 0) {		/* Tracking auxiliary segments */
			if (b->axisln >= b->axislz) {	/* Need some more space in list */
				if (b->axislz == 0) {
					b->axislz = 10;
					if ((b->axisl = (axisec *)rev_malloc(s, b->axislz * sizeof(axisec))) == NULL)
						error("rev: malloc failed - Auxiliary intersect list size %d",b->axislz);
					INCSZ(b->s, b->axislz * sizeof(axisec));
				} else {
					INCSZ(b->s, b->axislz * sizeof(axisec));
					b->axislz *= 2;
					if ((b->axisl = (axisec *)rev_realloc(s, b->axisl, b->axislz * sizeof(axisec)))
					    == NULL)
						error("rev: realloc failed - Auxiliary intersect list size %d",b->axislz);
				}
			}
			b->axisl[b->axisln].xval = xval;
			b->axisl[b->axisln].nv = x->sdi + 1;
			for (f = 0; f <= x->sdi; f++) {
				b->axisl[b->axisln].vix[f] = x->vix[f];
			}
			b->axisln++;
		}

#ifdef DEBUG
		if (xval >= b->min && xval <= b->max)
			DBG(("auxil_locus: solution %f doesn't improve on min %f, max %f\n",xval,b->min,b->max));
#endif
		/* If this solution is expands the min or max, save it */
		if (xval < b->min) {
			DBG(("######## Improving minimum to %f\n",xval));
			b->min = xval;
			b->plmincell = x->ix;
		}
		if (xval > b->max) {
			DBG(("######## Improving maximum to %f\n",xval));
			b->max = xval;
			b->plmaxcell = x->ix;
		}
	} else {
		DBG(("Solution wasn't within the simplex\n"));
		return 0;
	}

	return wsrv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Find the point on the clip line locus and simplexes */
/* valid surface, that is closest to the target output value. */
/* We expect to be given a sub simplex with sdi = fdi-1, and efdi = fdi */
/* or a limit sub-simplex with sdi = fdi, and efdi = fdi+1 */
/* Return zero if solution canot be calculated, */
/* return 1 normally, 2 if solution would be above the (disabled) ink limit */
static int
vnearest_clip_solve(
schbase *b,
simplex *x,
double *xp,		/* Return solution (simplex parameter space) */
double *xv,		/* Return solution (output space) */
double *err		/* Output error distance at solution point */
) {
	rspl *s = b->s;
	int e, sdi = x->sdi; 
	int f, fdi = s->fdi, efdi = x->efdi; 
	int g;
	int wsrv;	/* Within simplex return value */

	double *ta[MXRO], TA[MXRO][MXRO];
	double tb[MXRO];

	DBG(("Vector nearest clip solution called, cell %d, splx %d\n", x->ix, x->si));

	/* Setup temporary matricies */
	for (f = 0; f < sdi; f++) {
		ta[f] = TA[f];
	}

	/* Substitute simplex equation for output values V */
	/* in terms of sub-simplex parameters P, */
	/* into  clip line implicit equation in V, to give */
	/* clip line simplex implicit equation in terms of P (simplex input space) */
	/* If this is a limit sub-simlex, the ink limit part of the clip vector */
	/* equations will be used. */

	/* LHS: ta[sdi][sdi] = cla[sdi][efdi] * vv[efdi][sdi] */
	/* RHS: tb[sdi] = clb[sdi] - cla[sdi][efdi] * vv_di[efdi] */
	for (f = 0; f < sdi; f++) {
		double tt;
		for (e = 0; e < sdi; e++) {
			for (tt = 0.0, g = 0; g < efdi; g++)
				tt += b->cla[f][g] * (x->v[e][g] - x->v[e+1][g]);
			ta[f][e] = tt;
		}
		for (tt = 0.0, g = 0; g < efdi; g++)
			tt += b->cla[f][g] * x->v[sdi][g];
		tb[f] = b->clb[f] - tt;
	}

	/* Compute the solution */
	if (gen_solve_se(ta, tb, sdi, sdi)) {
		DBG(("Equation solution failed!\n"));
		return 0;		/* No solution */
	}

	/* Check that the solution is within the simplex & meets ink limit */
	if ((wsrv = within_simplex(x, tb)) != 0) {
		double dist;				/* distance to clip target */

		DBG(("Got solution within simplex %s\n", debPdv(sdi,tb)));

		/* Compute the output space solution point */
		for (f = 0; f < fdi; f++) {
			double tt = 0.0;
			for (e = 0; e < sdi; e++) {
				tt += (x->v[e][f] - x->v[e+1][f]) * tb[e];
			}
			xv[f] = tt + x->v[sdi][f];
		}

		/* Copy to return array */
		for (e = 0; e < sdi; e++)
			xp[e] = tb[e];

		// ~~~ are we properly checking if the intersection is 
		// ~~~ backwards rather than forwards in the line direction ?

		/* Compute distance to clip target */
		for (dist = 0.0, f = 0; f < fdi ; f++) {
			double tt = (b->v[f] - xv[f]);
			dist += tt * tt;
		}
		DBGV(("Vector clip output soln: ",fdi," %f", xv, "\n"));

		/* Return the solution in xp[]m xv[] and *err */
		*err = sqrt(dist);

		DBG(("Vector clip returning a solution with error %f\n",*err));
		return wsrv;
	}

	DBG(("Vector clip solution not in simplex\n"));
	return 0;		/* No solution */
}

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Find the point on the simplexes valid surface, that is closest */
/* to the target output value, for the linear (unweighted) case. */
/* We expect to be given a sub simplex with sdi = fdi-1, and efdi = fdi */
/* or a limit sub-simplex with sdi = fdi, and efdi = fdi+1 */
/* Return zero if solution canot be calculated, */
/* return 1 normally, 2 if solution would be above the (disabled) ink limit */
static int
nnearest_clip_solve(
schbase *b,
simplex *x,
double *xp,		/* Return solution (simplex parameter space) */
double *xv,		/* Return solution (output space) */
double *err		/* Output error (weighted) distance at solution point */
) {
	rspl *s = b->s;
	int e, sdi = x->sdi; 
	int f, fdi = s->fdi, efdi = x->efdi; 
	double tb[MXRO];		/* RHS & Parameter solution */
	double dist;			/* distance to clip target */
	int wsrv = 0;			/* Within simplex return value */

	DBG(("Nearest clip solution called, cell %d, splx %d\n", x->ix, x->si));

	if (sdi == 0) {		/* Solution is vertex */
		wsrv = 1;
		for (f = 0; f < efdi; f++)
			xv[f] = x->v[0][f]; 		/* Copy vertex value */
		if (x->v[0][fdi] > s->limitv) {
			if (s->limiten) 			/* Needed when limiten == 0 */
				return 0;				/* Over ink limit - no good */
			wsrv = 2;					/* Would be over */
		}
		DBG(("Got assumed vertex solution (vtx ix %d)\n",x->vix[0]));

	/* General linear nearest solver */
	} else {

#ifdef NEVER	/* Don't specialise ink limit version - use INKSCALE fudge instead */
		if (!(x->flags & SPLX_CLIPSX)) {	/* Not an ink limited plane simplex */
		
#endif
			/* Create the SVD decomp needed for least squares solution */
			if (add_lu_svd(x)) {
				DBG(("SVD decomp failed, skip simplex\n"));
				return 0;
			}
		
			/* Setup RHS to solve */
			for (f = 0; f < efdi; f++)
				tb[f] = b->v[f] - x->v[sdi][f]; 

			/* Find least squares solution */
			svdbacksub(x->d_u, x->d_w, x->d_v, tb, tb, efdi, sdi);
	
			/* Check that the solution is within the simplex & meets ink limit */
			if ((wsrv = within_simplex(x, tb)) == 0) {
				DBG(("Nearest clip solution not in simplex\n"));
				return 0;		/* No solution */
			}
	
			DBG(("Got solution within simplex %s\n",debPdv(sdi,tb)));
//if (b->rix == 7135) printf("Got solution within simplex params %s\n",debPdv(sdi,tb));
//if (b->rix == 7135) printf(" verticies ix %s\n",debPiv(sdi+1,x->vix));
	
			/* Compute the output space solution point */
			for (f = 0; f < fdi; f++) {
				double tt = 0.0;
				for (e = 0; e < sdi; e++)
					tt += (x->v[e][f] - x->v[e+1][f]) * tb[e];
				xv[f] = tt + x->v[sdi][f];
			}
//if (b->rix == 7135) printf("Computed Got simplex solution %s\n",debPdv(fdi,xv));
#ifdef NEVER /* ~~1 Haven't figured out equations to make this a special case. */
			 /* Content to use INKSCALE fudge and rely on SVD least squares. */
		} else {
			/* We can't use the given equations, because we want the solution */
			/* to lie exactly on the ink limit plane, and be least squares to the */
			/* other target parameters. */
			/* Extract the ink limit parameters, and transform them into */
			/* a parameterised surface for this simplex. */
			/* Substitute the ink plane equation into the remaining target */
			/* parameter equations, and solve for least squares. */

		}
#endif

		/* Copy to return array */
		for (e = 0; e < sdi; e++)
			xp[e] = tb[e];
	}

	/* Compute weighted distance to clip target */
	dist = sqrt(lchw_sq(s, b->v, xv));

//if (b->rix == 7135 && dist < b->cdist) {
//	printf("Got dist %f from %s -> %s with weight %d, %s\n", dist,debPdv(fdi,b->v),debPdv(fdi,xv),s->rev.lchweighted,debPdv(fdi,s->rev.lchw)); }

	DBGV(("Nearest clip output soln: ",fdi," %f", xv, "\n"));

	/* Return the solution in xp[], xv[] and *err */
	*err = dist;

	DBG(("Nearest clip returning a solution with error %f\n",*err));
	return wsrv;
}


#ifdef NEVER
/* Utility to convert an implicit ink limit plane equation */
/* held at the end of the simplex output value equations), */
/* into a parameterized surface equation. */
static void
compute_param_limit_surface(
schbase *b,
simplex *x
) {
	rspl *s = b->s;
	int ff, f, fdi = s->fdi;
	int i, p;
	double lgst;

double st[MXRO],	/* Start point */
double de[MXRO]		/* Delta */
	DBG(("Computing clipping line implicit equation, dim = %d\n", fdi));
	
	/* Pick a pivot element - the smallest */
	for (lgst = -1.0, p = -1, f = 0; f < fdi; f++) {
		double tt = de[f];
		b->cdir[f] = tt;		/* Stash this away */
		tt = fabs(tt);
		if (tt > lgst) {
			lgst = tt;
			p = f;
		}
	}
	if (p < 0)	/* Shouldn't happen */
		error("rspl rev, internal, trying to cope with zero length clip line\n");
	
	if (b->cla == NULL)
		b->cla = dmatrix(0, fdi-1, 0, fdi);	/* Allow for ink limit supliment */

	for (i = ff = 0;  ff < fdi; ff++) {	/* For the input rows */
		if (ff == p) {
			continue;					/* Skip pivot row */
		}
		for (f = 0; f < fdi; f++) {		/* For input & output columns */
			if (f == p) {
				b->cla[i][f] = -de[ff];	/* Last column is -ve delta value */
			} else if (f == ff) {
				b->cla[i][f] = de[p];	/* Diagonal is pivot value */
			} else {
				b->cla[i][f] = 0.0;		/* Else zero */
			}
		}
		b->clb[i] = de[p] * st[ff] - de[ff] * st[p];
		i++;
	}

	/* Add ink limit target equation - */
	/* interpolated ink value == target */
	if (s->limitf != NULL) {
		for (i = 0;  i < (fdi-1); i++)
			b->cla[i][fdi] = 0.0;

		for (f = 0; f < fdi; f++) 
			b->cla[fdi-1][f] = 0.0;
		
		b->cla[fdi-1][fdi] = 1.0;
		b->clb[fdi-1] = s->limitv;
	}

#ifdef NEVER
/* Verify that the implicit equation is correct */
{
	double pnt[MXRO], v[MXRO];
	double pa;	/* Parameter */
	for (pa = 0.0; pa <= 1.0; pa += 0.125) {
		for (f = 0; f < fdi; f++) {
			pnt[f] = st[f] + pa * de[f];
		}

		/* Verify the implicit equation */
		for (ff = 0; ff < (fdi-1); ff++) {
			v[ff] = 0.0;
			for (f = 0; f < fdi; f++) {
				v[ff] += b->cla[ff][f] * pnt[f];
			}
			v[ff] -= b->clb[ff];
			if (v[ff] < 0.0)
				v[ff] = -v[ff];
			if (v[ff] > 0.000001) {
				printf("Point on clip line = %f %f %f\n",pnt[0],pnt[1],pnt[2]);
				printf("Implicit %d error of = %f\n",ff, v[ff]);
			}
		}
	}
}
#endif /* NEVER */

}

#endif


/* -------------------------------------------------------- */
static int lchw_edge_solve(rspl *s, double *vv, double *p, double *vt, double v[MXRI+1][MXRO+1]);
static int lchw_tri_solve(rspl *s, double *vv, double *p, double *vt, double v[MXRI+1][MXRO+1]); 

/* Find the point on the simplexes valid surface, that is closest */
/* to the target output value, for the LCh weighted case. */
/* We use Newton itteration to solve this for the 1D (line) and 2D (triangle) */
/* cases, and explicitly decode the ink limit surfaces back to point, line */
/* and triangled cases. */ 
/* We expect to be given a sub simplex with sdi = 0..2, and efdi = fdi */
/* or a limit sub-simplex with sdi = 1..3, and efdi = fdi+1 */
/* We bail with an assert if we get more than 2D to solve. */
/* Return zero if solution canot be calculated, */
/* return 1 normally, 2 if solution would be above the (disabled) ink limit */
static int
lchw_nnearest_clip_solve(
schbase *b,
simplex *x,
double *xp,		/* Return solution (simplex parameter space) */
double *xv,		/* Return solution (output space) */
double *err		/* Output error (weighted) distance at solution point */
) {
	rspl *s = b->s;
	int e, ee, sdi = x->sdi; 
	int f, fdi = s->fdi, efdi = x->efdi; 
	double tb[MXRO];		/* RHS & Parameter solution */
	double dist;			/* distance to clip target */
	int wsrv = 0;			/* Within simplex return value */

	DBG(("LChw nearest clip solution called, cell %d, splx %d\n", x->ix, x->si));

	/* - - - - - - - */
	if (sdi == 0) {		/* Solution is vertex */
		wsrv = 1;
		for (f = 0; f < efdi; f++)
			xv[f] = x->v[0][f]; 		/* Copy vertex value */
		if (x->v[0][fdi] > s->limitv) {
			if (s->limiten) 			/* Needed when limiten == 0 */
				return 0;				/* Over ink limit - no good */
			wsrv = 2;					/* Would be over */
		}
		DBG(("Got assumed vertex solution (vtx ix %d)\n",x->vix[0]));

	/* - - - - - - - */
	/* Ink limit simplex case */
	} else if (efdi == (fdi+1)) {

		/* Convert line into vertex and return it */
		if (sdi == 1) {
			wsrv = 1;

			/* Ink limit plane point along line */
			xp[0] = (s->limitv - x->v[1][fdi])/(x->v[0][fdi] - x->v[1][fdi]);

			/* Output value at that point */
			for (f = 0; f < fdi; f++)
				xv[f] = (x->v[0][f] - x->v[1][f]) * xp[0] + x->v[1][f];

			DBG(("Got ink limit point on edge\n"));

		/* Turn triangle into line and solve line. */
		} else if (sdi == 2) {
			int pos = 0, neg = 0;
			int ix[MXRI+1];		/* Odd index and the two other indexes */
			double p[MXRI+1], pp[MXRI+1];
			double v[MXRI+1][MXRO+1];

			/* Count ink limit signs of vertexes */
			for (e = 0; e <= sdi; e++) {
				ix[e] = e;
				if (x->v[e][fdi] > s->limitv)
					pos++;
				else
					neg++;
			}

			/* We expect one vertex to be on the other side of the */
			/* ink limit to the two others. */
			if (pos == 0 || neg == 0)
				error("Ink limit tri doesn't have one opposite sign");

			/* Make the first ix be the odd one */
			if (pos == 1) {
				if (x->v[0][fdi] <= s->limitv) {
					if (x->v[1][fdi] > s->limitv) {
						ix[0] = 1;
						ix[1] = 0;
					} else {
						ix[0] = 2;
						ix[2] = 0;
					}
				}
			} else {
				if (x->v[0][fdi] > s->limitv) {
					if (x->v[1][fdi] <= s->limitv) {
						ix[0] = 1;
						ix[1] = 0;
					} else {
						ix[0] = 2;
						ix[2] = 0;
					}
				}
			}
			
			/* Compute the points on the two edges that cross the ink limit. */
			/* i.e. for edges ix 0..1 & 0..2 */
			pp[0] = (s->limitv - x->v[ix[1]][fdi])/(x->v[ix[0]][fdi] - x->v[ix[1]][fdi]);
			pp[1] = (s->limitv - x->v[ix[2]][fdi])/(x->v[ix[0]][fdi] - x->v[ix[2]][fdi]);
			for (f = 0; f < fdi; f++) {
				v[0][f] = (x->v[ix[0]][f] - x->v[ix[1]][f]) * pp[0] + x->v[ix[1]][f];
				v[1][f] = (x->v[ix[0]][f] - x->v[ix[2]][f]) * pp[1] + x->v[ix[2]][f];
			}
			
			/* Solve it */
			if ((wsrv = lchw_edge_solve(s, xv, p, b->v, v)) != 0) {

				/* Figure out the solution simplex coords */
				/* (p is weighting of lower indexes vertex) */

				/* Convert solution simplex coords into baricentric weighting */ 
				p[1] = 1.0 - p[0];

				/* Sum baricentric weightings for each vertex */
				for (e = 0; e <= sdi; e++)
					xp[e] = 0.0;

				xp[ix[0]] += pp[0] * p[0];
				xp[ix[1]] += (1.0 - pp[0]) * p[0];
				xp[ix[0]] += pp[1] * p[1];
				xp[ix[2]] += (1.0 - pp[1]) * p[1];

				/* Convert back to simplex coords */
				xp[1] = 1.0 - xp[2];
				xp[0] = xp[0];

				DBG(("Got ink limit edge in triangle\n"));
			}

		/* Turn tetrahedron into one or two triangles */
		/* and solve triangles. */
		} else if (sdi == 3) {
			int pos = 0, neg = 0;
			int ix[MXRI+1];		/* Odd index and the three other indexes or 2 + 2 */
			double p[MXRI+1], pp[MXRI+1];
			double v[MXRI+1][MXRO+1];

			/* Count ink limit signs of vertexes */
			for (e = 0; e <= sdi; e++) {
				ix[e] = e;
				if (x->v[e][fdi] > s->limitv)
					pos++;
				else
					neg++;
			}

			/* We expect one or two vertexes t be on the other side of the */
			/* ink limit to the two others. */
			if (pos == 0 || neg == 0)
				error("Ink limit tetrahedron doesn't have one opposite sign");

			/* If we can decompose this into a single triangle */
			if (pos == 1 || neg == 1) {

				/* Make the first ix be the odd one */
				for (e = 0; e <= sdi; e++) {
					if ((pos == 1 && x->v[e][fdi] > s->limitv)
					 || (neg == 1 && x->v[e][fdi] <= s->limitv)) {
						int tt = ix[0];
						ix[0] = e;
						ix[e] = tt;
						break;
					}
				}

				/* Compute the points on the three edges that cross the ink limit. */
				/* i.e. for edges ix 0..1, 0..2 & 0..3 */
				pp[0] = (s->limitv - x->v[ix[1]][fdi])/(x->v[ix[0]][fdi] - x->v[ix[1]][fdi]);
				pp[1] = (s->limitv - x->v[ix[2]][fdi])/(x->v[ix[0]][fdi] - x->v[ix[2]][fdi]);
				pp[2] = (s->limitv - x->v[ix[3]][fdi])/(x->v[ix[0]][fdi] - x->v[ix[3]][fdi]);
				for (f = 0; f < fdi; f++) {
					v[0][f] = (x->v[ix[0]][f] - x->v[ix[1]][f]) * pp[0] + x->v[ix[1]][f];
					v[1][f] = (x->v[ix[0]][f] - x->v[ix[2]][f]) * pp[1] + x->v[ix[2]][f];
					v[2][f] = (x->v[ix[0]][f] - x->v[ix[3]][f]) * pp[2] + x->v[ix[3]][f];
				}
				
				/* Solve it */
				if ((wsrv = lchw_tri_solve(s, xv, p, b->v, v)) != 0) {

					/* Figure out the solution simplex coords */
					/* (p is weighting of lower indexes vertex) */

					/* Convert solution simplex coords into baricentric weighting */ 
					p[2] = 1.0 - p[1];
					p[1] = p[1] - p[0];
					p[0] = p[0];

					/* Sum baricentric weightings for each vertex */
					for (e = 0; e <= sdi; e++)
						xp[e] = 0.0;

					xp[ix[0]] += pp[0] * p[0];
					xp[ix[1]] += (1.0 - pp[0]) * p[0];
					xp[ix[0]] += pp[1] * p[1];
					xp[ix[2]] += (1.0 - pp[1]) * p[1];
					xp[ix[0]] += pp[2] * p[2];
					xp[ix[3]] += (1.0 - pp[2]) * p[2];

					/* Convert back to simplex coords */
					xp[2] = 1.0 - xp[3];
					xp[1] = xp[1] + xp[0];
					xp[0] = xp[0];

					DBG(("Got single ink limit triangle in tetrahedron\n"));
				}

			/* We need to decompose this into two triangles */
			} else {
				int wsrv2 = 0;
				double dist2;
				double xv2[MXRO];		/* 2nd triangle solution */

				/* Make the first two ix's be the same, leaving second two the same. */
				for (e = 1; e <= sdi; e++) {
					if (x->v[0][fdi] > s->limitv && x->v[e][fdi] > s->limitv) {
						int tt = ix[1];
						ix[1] = e;
						ix[e] = tt;
						break;
					}
				}

				/* We choose disjoint vertex pairs as the common edge of the two */
				/* triangles, and then use each of the remaining pairs to form */
				/* the other edges. */
				/* i.e. common edge 0..2 + 1..3, then add 0..3 then 1..2 */
				pp[0] = (s->limitv - x->v[ix[2]][fdi])/(x->v[ix[0]][fdi] - x->v[ix[2]][fdi]);
				pp[1] = (s->limitv - x->v[ix[3]][fdi])/(x->v[ix[1]][fdi] - x->v[ix[3]][fdi]);
				pp[2] = (s->limitv - x->v[ix[3]][fdi])/(x->v[ix[0]][fdi] - x->v[ix[3]][fdi]);
				for (f = 0; f < fdi; f++) {
					v[0][f] = (x->v[ix[0]][f] - x->v[ix[2]][f]) * pp[0] + x->v[ix[2]][f];
					v[1][f] = (x->v[ix[1]][f] - x->v[ix[3]][f]) * pp[1] + x->v[ix[3]][f];
					v[2][f] = (x->v[ix[0]][f] - x->v[ix[3]][f]) * pp[2] + x->v[ix[3]][f];
				}
				
				/* Solve first one */
				if ((wsrv = lchw_tri_solve(s, xv, p, b->v, v)) != 0) {

					dist = sqrt(lchw_sq(s, b->v, xv));

					/* Figure out the solution simplex coords */
					/* (p is weighting of lower indexes vertex) */

					/* Convert solution simplex coords into baricentric weighting */ 
					p[2] = 1.0 - p[1];
					p[1] = p[1] - p[0];
					p[0] = p[0];

					/* Sum baricentric weightings for each vertex */
					for (e = 0; e <= sdi; e++)
						xp[e] = 0.0;

					xp[ix[0]] += pp[0] * p[0];
					xp[ix[2]] += (1.0 - pp[0]) * p[0];
					xp[ix[1]] += pp[1] * p[1];
					xp[ix[3]] += (1.0 - pp[1]) * p[1];
					xp[ix[0]] += pp[2] * p[2];
					xp[ix[3]] += (1.0 - pp[2]) * p[2];

					/* Convert back to simplex coords */
					xp[2] = 1.0 - xp[3];
					xp[1] = xp[1] + xp[0];
					xp[0] = xp[0];
				}

				/* Setup other triangle, 0..2 + 1..3, with 1..2 */
				pp[0] = (s->limitv - x->v[ix[2]][fdi])/(x->v[ix[0]][fdi] - x->v[ix[2]][fdi]);
				pp[1] = (s->limitv - x->v[ix[3]][fdi])/(x->v[ix[1]][fdi] - x->v[ix[3]][fdi]);
				pp[2] = (s->limitv - x->v[ix[2]][fdi])/(x->v[ix[1]][fdi] - x->v[ix[2]][fdi]);
				for (f = 0; f < fdi; f++) {
					v[0][f] = (x->v[ix[0]][f] - x->v[ix[2]][f]) * pp[0] + x->v[ix[2]][f];
					v[1][f] = (x->v[ix[1]][f] - x->v[ix[3]][f]) * pp[1] + x->v[ix[3]][f];
					v[2][f] = (x->v[ix[1]][f] - x->v[ix[2]][f]) * pp[2] + x->v[ix[2]][f];
				}
				
				/* Solve second triangle */
				if ((wsrv2 = lchw_tri_solve(s, xv2, p, b->v, v)) != 0) {

					dist2 = sqrt(lchw_sq(s, b->v, xv2));

					/* Use this second solution */
					if (wsrv == 0 || dist2 < dist) {

						dist = dist2;

						/* Figure out the solution simplex coords */
						/* (p is weighting of lower indexes vertex) */
	
						/* Convert solution simplex coords into baricentric weighting */ 
						p[2] = 1.0 - p[1];
						p[1] = p[1] - p[0];
						p[0] = p[0];
	
						/* Sum baricentric weightings for each vertex */
						for (e = 0; e <= sdi; e++)
							xp[e] = 0.0;
	
						xp[ix[0]] += pp[0] * p[0];
						xp[ix[2]] += (1.0 - pp[0]) * p[0];
						xp[ix[1]] += pp[1] * p[1];
						xp[ix[3]] += (1.0 - pp[1]) * p[1];
						xp[ix[1]] += pp[2] * p[2];
						xp[ix[2]] += (1.0 - pp[2]) * p[2];
	
						/* Convert back to simplex coords */
						xp[2] = 1.0 - xp[3];
						xp[1] = xp[1] + xp[0];
						xp[0] = xp[0];

						for (f = 0; f < fdi; f++)
							xv[f] = xv2[f];

					} else {
						wsrv2 = 0;
					}
				}

#ifdef DEBUG
				if (wsrv2) {
					DBG(("Got second ink limit triangle in tetrahedron\n"));
				} else if (wsrv) {
					DBG(("Got first ink limit triangle in tetrahedron\n"));
				}
#endif
				*err = dist;
				return wsrv;
			}
		} else {
			error("rev: lchw_nnearest_clip_solve sdi = %d\n",sdi);
		}

		/* All solutions computed on the ink limit surface */
		/* are assumed to be valid */

	/* - - - - - - - */
	/* Non-ink limit simplex case */
	} else {

		/* Line */
		if (sdi == 1) {
			wsrv = lchw_edge_solve(s, xv, xp, b->v, x->v);

			DBG(("Got line solution\n"));

		/* Triangle */
		} else if (sdi == 2) {
			wsrv = lchw_tri_solve(s, xv, xp, b->v, x->v);

			DBG(("Got triangle solution\n"));

		/* Oops */
		} else {
			error("rev: lchw_nnearest_clip_solve sdi = %d\n",sdi);
		}

		/* Check that the result is within the ink limit */
		if (wsrv != 0)
			wsrv = within_simplex_limit(x, xp);
	}

	if (wsrv == 0)
		return wsrv;

	/* Compute weighted distance to clip target */
	dist = sqrt(lchw_sq(s, b->v, xv));

	DBGV(("LChw nearest clip output soln: ",fdi," %f", xv, "\n"));

	/* Return the solution in xp[], xv[] and *err */
	*err = dist;

	DBG(("LChw nearest clip returning a solution with error %f\n",*err));

#ifdef NEVER
	{
		double chxv[MXRO]; 
	
		printf("LChw nearest clip returning a solution with error %f\n",dist);
	
		printf("Solution (sx in) %s -> out %s\n", debPdv(sdi, xp), debPdv(fdi, xv));
	
		if (dist < b->cdist) {	/* Equal or worse clip solution */
			printf("Will be new best solution\n");
		}
	
		/* Check the output space solution point */
		for (f = 0; f < fdi; f++) {
			double tt = 0.0;
			for (e = 0; e < sdi; e++)
				tt += (x->v[e][f] - x->v[e+1][f]) * xp[e];
			chxv[f] = tt + x->v[sdi][f];
		}
		for (f = 0; f < fdi; f++) {
			if (fabs(chxv[f] - xv[f]) > 1e-3) {
				break;
			}
		}
		if (f < fdi)
			printf(" ###### Check of out failed: %s\n", debPdv(fdi, chxv));
	}
#endif

	return wsrv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Edge lchw Newton itteration code */ 

#ifdef NEVER		/* Not actually used here */
/* return weighted delta squared for target to edge at param value p */
static double lchw_edge_sq(rspl *s, double *vt, double v[MXRI+1][MXRO+1], double p) {
	int f, fdi = s->fdi;
	double vv[MXRO];			/* Point at parameter location */
	double dlsq;				/* Delta L squared */
	double da, db, dchsq;		/* Delta CH squared */
	double ct, cv, dc, dcsq;	/* Delta C squared */
	double lcomp, chcomp, ccomp;
	double de;

	/* Compute point at parameter location */
	for (f = 0; f < fdi; f++)
		vv[f] = (v[0][f] - v[1][f]) * p + v[1][f];

	/* Delta L component */
	dlsq = vv[0] - vt[0];
	dlsq = dlsq * dlsq;
	lcomp = s->rev.lchw_sq[0] * dlsq;

	/* Delta CH component */
	da = vv[1] - vt[1];
	db = vv[2] - vt[2];
	dchsq = da * da + db * db;
	chcomp = s->rev.lchw_sq[2] * dchsq;

	/* Compute chromanance for the two colors */	
	ct = sqrt(vt[1] * vt[1] + vt[2] * vt[2]);	
	cv = sqrt(vv[1] * vv[1] + vv[2] * vv[2]);	
	dc = ct - cv;	
	dcsq = dc * dc;	
	
	ccomp = s->rev.lchw_chsq * dcsq;		/* w = cw - hw because dh = dch - dc */

	de = lcomp + chcomp + ccomp;

	return de;
}
#endif /* NEVER */

/* return weighted 1st derivativ of delta squared for target to edge at param value p */
static double lchw_edge_Dp_sq(rspl *s, double *vt, double v[MXRI+1][MXRO+1], double p) {
	int f, fdi = s->fdi;
	double vv[MXRO];			/* Point at parameter location */
	double Dvv[MXRO];			/* Derivative wrt p of vv */
	double dl, Ddlsq;			/* Delta L squared */
	double da, Ddasq, db, Ddbsq, Ddchsq;		/* Delta CH squared */
	double ct, cv, Dcv, dc, Ddc, Dvv1sq, Dvv2sq, Ddcsq;	/* Delta C squared */
	double Dlcomp, Dchcomp, Dccomp;
	double Dde;

	/* Compute point at parameter location */
	for (f = 0; f < fdi; f++) {
		vv[f] = (v[0][f] - v[1][f]) * p + v[1][f];
		Dvv[f] = v[0][f] - v[1][f];
	}

	/* Delta L component */
	dl = vv[0] - vt[0];
	Ddlsq = 2.0 * dl * Dvv[0];
	Dlcomp = s->rev.lchw_sq[0] * Ddlsq;

	/* Delta CH component */
	da = vv[1] - vt[1];
	db = vv[2] - vt[2];
	Ddasq = 2.0 * da * Dvv[1];
	Ddbsq = 2.0 * db * Dvv[2];
	Ddchsq = Ddasq + Ddbsq;
	Dchcomp = s->rev.lchw_sq[2] * Ddchsq;

	/* Compute chromanance for the two colors */	
	ct = sqrt(vt[1] * vt[1] + vt[2] * vt[2]);	
	cv = sqrt(vv[1] * vv[1] + vv[2] * vv[2]);	
	dc = cv - ct;	
	Dvv1sq = 2.0 * vv[1] * Dvv[1];
	Dvv2sq = 2.0 * vv[2] * Dvv[2];
	Dcv = 0.5/cv * (Dvv1sq + Dvv2sq);
	Ddcsq = 2.0 * dc * Dcv;	
	Dccomp = s->rev.lchw_chsq * Ddcsq;

	Dde = Dlcomp + Dchcomp + Dccomp;

	return Dde;
}

/* return weighted 2nd derivative of delta squared for target to edge at param value p */
static double lchw_edge_DDp_sq(rspl *s, double *vt, double v[MXRI+1][MXRO+1], double p) {
	int f, fdi = s->fdi;
	double vv[MXRO];			/* Point at parameter location */
	double Dvv[MXRO];			/* Derivative wrt p of vv */
	double DDvvsq[MXRO];		/* 2nd Derivative wrt p of vv */
	double DDdchsq;
	double ct, cv, Dcv, DDcv, dc, Dvv1sq, Dvv2sq, DDdcsq;
	double DDlcomp, DDchcomp, DDccomp;
	double DDde;

	/* Compute point at parameter location */
	for (f = 0; f < fdi; f++) {
		vv[f] = (v[0][f] - v[1][f]) * p + v[1][f];
		Dvv[f] = v[0][f] - v[1][f];
		DDvvsq[f] = 2.0 * Dvv[f] * Dvv[f];
	}

	/* Delta L component */
	DDlcomp = s->rev.lchw_sq[0] * DDvvsq[0];

	/* Delta CH component */
	DDdchsq = DDvvsq[1] + DDvvsq[2];
	DDchcomp = s->rev.lchw_sq[2] * DDdchsq;

	/* Compute chromanance for the two colors */	
	ct = sqrt(vt[1] * vt[1] + vt[2] * vt[2]);	
	cv = sqrt(vv[1] * vv[1] + vv[2] * vv[2]);	
	dc = cv - ct;	
	Dvv1sq = 2.0 * vv[1] * Dvv[1];
	Dvv2sq = 2.0 * vv[2] * Dvv[2];

	Dcv = 0.5/cv * (Dvv1sq + Dvv2sq);
	DDcv = -0.5/(cv * cv) * Dcv * (Dvv1sq + Dvv2sq) + 0.5/cv * (DDvvsq[1] + DDvvsq[2]);

	DDdcsq = 2.0 * (Dcv * Dcv + dc * DDcv);	
	DDccomp = s->rev.lchw_chsq * DDdcsq;

	DDde = DDlcomp + DDchcomp + DDccomp;

	return DDde;
}

/* Solve for an edge. Return nz of solution. */
static int lchw_edge_solve(rspl *s, double *vv, double *p, double *vt, double v[MXRI+1][MXRO+1]) {
	int i, f, fdi = s->fdi;
	double pp, ee, dedp;
	double e0, e1; 

	/* Decide whether there is a solution on this edge. */
	/* This is reliable, and saves any itters in the loop. */
	e0 = lchw_edge_Dp_sq(s, vt, v, 0.0);
	e1 = lchw_edge_Dp_sq(s, vt, v, 1.0);

	if ((e0 < 0.0 && e1 < 0.0)
	 || (e0 > 0.0 && e1 > 0.0)) {
		return 0;
	}

	pp = 0.5;
	for (i = 0; i < 30; i++) {
		ee = lchw_edge_Dp_sq(s, vt, v, pp);
		dedp = lchw_edge_DDp_sq(s, vt, v, pp);
		pp -= ee/dedp;  

		if (fabs(ee) < 1e-6)
			break;
	}
	ee = lchw_edge_Dp_sq(s, vt, v, pp);

	if (fabs(ee) > 1e-6 || pp < -EPS || pp > (1.0 + EPS)) {
		return 0;
	}

	/* Return solution (output space) */
	for (f = 0; f < fdi; f++)
		vv[f] = (v[0][f] - v[1][f]) * pp + v[1][f];

	/* Return solution (simplex parameter space) */
	*p = pp;

	return 1;
}

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Triangle lchw Newton itteration code */ 

/* return weighted delta squared for target to triangle at param values p */
/* [ 0 <= p0 <= p1 <= 1 ] */
static double lchw_tri_sq(rspl *s, double *vt, double v[MXRI+1][MXRO+1], double *p) {
	int f, fdi = s->fdi;
	double vv[MXRO];			/* Point at parameter location */
	double dlsq;				/* Delta L squared */
	double da, db, dchsq;		/* Delta CH squared */
	double ct, cv, dc, dcsq;	/* Delta C squared */
	double lcomp, chcomp, ccomp;
	double de;

	/* Compute point at parameter location */
	for (f = 0; f < fdi; f++)
		vv[f] = (v[0][f] - v[1][f]) * p[0]
		      + (v[1][f] - v[2][f]) * p[1]
		      +  v[2][f];

	/* Delta L component */
	dlsq = vv[0] - vt[0];
	dlsq = dlsq * dlsq;
	lcomp = s->rev.lchw_sq[0] * dlsq;

	/* Delta CH component */
	da = vv[1] - vt[1];
	db = vv[2] - vt[2];
	dchsq = da * da + db * db;
	chcomp = s->rev.lchw_sq[2] * dchsq;

	/* Compute chromanance for the two colors */	
	ct = sqrt(vt[1] * vt[1] + vt[2] * vt[2]);	
	cv = sqrt(vv[1] * vv[1] + vv[2] * vv[2]);	
	dc = ct - cv;	
	dcsq = dc * dc;	
	
	ccomp = s->rev.lchw_chsq * dcsq;		/* w = cw - hw because dh = dch - dc */

	de = lcomp + chcomp + ccomp;

	return de;
}

/* return weighted two 1st derivativ of delta squared for target to edge at param value p */
static void lchw_tri_Dp_sq(rspl *s, double Dde[2], double *vt, double v[MXRI+1][MXRO+1], double *p) {
	int f, fdi = s->fdi;
	double vv[MXRO];			/* Point at parameter location */
	double Dvv[2][MXRO];			/* Derivative wrt p of vv */
	double dl, Ddl[2], Ddlsq[2];			/* Delta L squared */
	double da, Ddasq[2], db, Ddbsq[2], Ddchsq[2];		/* Delta CH squared */
	double ct, cv, Dcv[2], dc, Dvv1sq[2], Dvv2sq[2], Ddcsq[2];	/* Delta C squared */
	double Dlcomp[2], Dchcomp[2], Dccomp[2];

	/* Compute point at parameter location */
	for (f = 0; f < fdi; f++) {
		vv[f] = (v[0][f] - v[1][f]) * p[0]
		      + (v[1][f] - v[2][f]) * p[1]
		      +  v[2][f];
		Dvv[0][f] = v[0][f] - v[1][f];
		Dvv[1][f] = v[1][f] - v[2][f];
	}

	/* Delta L component */
	dl = vv[0] - vt[0];
	Ddlsq[0] = 2.0 * dl * Dvv[0][0];
	Ddlsq[1] = 2.0 * dl * Dvv[1][0];
	Dlcomp[0] = s->rev.lchw_sq[0] * Ddlsq[0];
	Dlcomp[1] = s->rev.lchw_sq[0] * Ddlsq[1];

	/* Delta CH component */
	da = vv[1] - vt[1];
	db = vv[2] - vt[2];
	Ddasq[0] = 2.0 * da * Dvv[0][1];
	Ddasq[1] = 2.0 * da * Dvv[1][1];
	Ddbsq[0] = 2.0 * db * Dvv[0][2];
	Ddbsq[1] = 2.0 * db * Dvv[1][2];
	Ddchsq[0] = Ddasq[0] + Ddbsq[0];
	Ddchsq[1] = Ddasq[1] + Ddbsq[1];
	Dchcomp[0] = s->rev.lchw_sq[2] * Ddchsq[0];
	Dchcomp[1] = s->rev.lchw_sq[2] * Ddchsq[1];

	/* Compute chromanance for the two colors */	
	ct = sqrt(vt[1] * vt[1] + vt[2] * vt[2]);	
	cv = sqrt(vv[1] * vv[1] + vv[2] * vv[2]);	
	dc = cv - ct;	
	Dvv1sq[0] = 2.0 * vv[1] * Dvv[0][1];
	Dvv1sq[1] = 2.0 * vv[1] * Dvv[1][1];
	Dvv2sq[0] = 2.0 * vv[2] * Dvv[0][2];
	Dvv2sq[1] = 2.0 * vv[2] * Dvv[1][2];
	Dcv[0] = 0.5/cv * (Dvv1sq[0] + Dvv2sq[0]);
	Dcv[1] = 0.5/cv * (Dvv1sq[1] + Dvv2sq[1]);
	Ddcsq[0] = 2.0 * dc * Dcv[0];	
	Ddcsq[1] = 2.0 * dc * Dcv[1];	
	Dccomp[0] = s->rev.lchw_chsq * Ddcsq[0];
	Dccomp[1] = s->rev.lchw_chsq * Ddcsq[1];

	Dde[0] = Dlcomp[0] + Dchcomp[0] + Dccomp[0];
	Dde[1] = Dlcomp[1] + Dchcomp[1] + Dccomp[1];
}

/* return weighted four 2nd derivatives of delta squared for target to edge at param value p */
/* ([first][second]) */
static void lchw_tri_DDp_sq(rspl *s, double DDde[2][2], double *vt, double v[MXRI+1][MXRO+1], double *p) {
	int f, fdi = s->fdi;
	double vv[MXRO];					/* Point at parameter location */
	double Dvv[2][MXRO];				/* Derivative wrt p of vv */
	double DDvvsq[2][2][MXRO];			/* 2nd Derivative wrt p of vv */
	double DDdchsq[2][2];				/* Delta CH squared */
	double ct, cv, Dcv[2], DDcv[2][2], dc, Dvv1sq[2], Dvv2sq[2], DDdcsq[2][2];
	double DDlcomp[2][2], DDchcomp[2][2], DDccomp[2][2];

	/* Due to comutivity, [0][1] == [1][0], so we omit */
	/* those redundant calculations. */

	/* Compute point at parameter location */
	for (f = 0; f < fdi; f++) {
		vv[f] = (v[0][f] - v[1][f]) * p[0]
		      + (v[1][f] - v[2][f]) * p[1]
		      +  v[2][f];
		Dvv[0][f] = v[0][f] - v[1][f];
		Dvv[1][f] = v[1][f] - v[2][f];

		DDvvsq[0][0][f] = 2.0 * Dvv[0][f] * Dvv[0][f];
		DDvvsq[1][0][f] = 2.0 * Dvv[1][f] * Dvv[0][f];
//		DDvvsq[0][1][f] = 2.0 * Dvv[0][f] * Dvv[1][f];
		DDvvsq[1][1][f] = 2.0 * Dvv[1][f] * Dvv[1][f];
	}

	/* Delta L component */
	DDlcomp[0][0] = s->rev.lchw_sq[0] * DDvvsq[0][0][0];
	DDlcomp[1][0] = s->rev.lchw_sq[0] * DDvvsq[1][0][0];
//	DDlcomp[0][1] = s->rev.lchw_sq[0] * DDvvsq[0][1][0];
	DDlcomp[1][1] = s->rev.lchw_sq[0] * DDvvsq[1][1][0];

	/* Delta CH component */
	DDdchsq[0][0] = DDvvsq[0][0][1] + DDvvsq[0][0][2];
	DDdchsq[1][0] = DDvvsq[1][0][1] + DDvvsq[1][0][2];
//	DDdchsq[0][1] = DDvvsq[0][1][1] + DDvvsq[0][1][2];
	DDdchsq[1][1] = DDvvsq[1][1][1] + DDvvsq[1][1][2];

	DDchcomp[0][0] = s->rev.lchw_sq[2] * DDdchsq[0][0];
	DDchcomp[1][0] = s->rev.lchw_sq[2] * DDdchsq[1][0];
//	DDchcomp[0][1] = s->rev.lchw_sq[2] * DDdchsq[0][1];
	DDchcomp[1][1] = s->rev.lchw_sq[2] * DDdchsq[1][1];

	/* Compute chromanance for the two colors */	
	ct = sqrt(vt[1] * vt[1] + vt[2] * vt[2]);	
	cv = sqrt(vv[1] * vv[1] + vv[2] * vv[2]);	
	dc = cv - ct;	

	Dvv1sq[0] = 2.0 * vv[1] * Dvv[0][1];
	Dvv1sq[1] = 2.0 * vv[1] * Dvv[1][1];

	Dvv2sq[0] = 2.0 * vv[2] * Dvv[0][2];
	Dvv2sq[1] = 2.0 * vv[2] * Dvv[1][2];

	Dcv[0] = 0.5/cv * (Dvv1sq[0] + Dvv2sq[0]);
	Dcv[1] = 0.5/cv * (Dvv1sq[1] + Dvv2sq[1]);


	DDcv[0][0] = -0.5/(cv * cv) * Dcv[0] * (Dvv1sq[0] + Dvv2sq[0])
	            + 0.5/cv * (DDvvsq[0][0][1] + DDvvsq[0][0][2]);

	DDcv[1][0] = -0.5/(cv * cv) * Dcv[0] * (Dvv1sq[1] + Dvv2sq[1])
		        + 0.5/cv * (DDvvsq[1][0][1] + DDvvsq[1][0][2]);

//	DDcv[0][1] = -0.5/(cv * cv) * Dcv[1] * (Dvv1sq[0] + Dvv2sq[0])
//	            + 0.5/cv * (DDvvsq[0][1][1] + DDvvsq[0][1][2]);

	DDcv[1][1] = -0.5/(cv * cv) * Dcv[1] * (Dvv1sq[1] + Dvv2sq[1])
		        + 0.5/cv * (DDvvsq[1][1][1] + DDvvsq[1][1][2]);

	DDdcsq[0][0] = 2.0 * (Dcv[0] * Dcv[0] + dc * DDcv[0][0]);	
	DDdcsq[1][0] = 2.0 * (Dcv[1] * Dcv[0] + dc * DDcv[1][0]);	
//	DDdcsq[0][1] = 2.0 * (Dcv[0] * Dcv[1] + dc * DDcv[0][1]);	
	DDdcsq[1][1] = 2.0 * (Dcv[1] * Dcv[1] + dc * DDcv[1][1]);	

	DDccomp[0][0] = s->rev.lchw_chsq * DDdcsq[0][0];
	DDccomp[1][0] = s->rev.lchw_chsq * DDdcsq[1][0];
//	DDccomp[0][1] = s->rev.lchw_chsq * DDdcsq[0][1];
	DDccomp[1][1] = s->rev.lchw_chsq * DDdcsq[1][1];

	DDde[0][0] = DDlcomp[0][0] + DDchcomp[0][0] + DDccomp[0][0];
	DDde[1][0] = DDlcomp[1][0] + DDchcomp[1][0] + DDccomp[1][0];
//	DDde[0][1] = DDlcomp[0][1] + DDchcomp[0][1] + DDccomp[0][1];
	DDde[0][1] = DDde[1][0];
	DDde[1][1] = DDlcomp[1][1] + DDchcomp[1][1] + DDccomp[1][1];
}


/* Solve for a triangle face. Return nz of solution. */
static int lchw_tri_solve(rspl *s, double *vv, double *p, double *vt, double v[MXRI+1][MXRO+1]) {
	int f, fdi = s->fdi;
	int i, j, k;
	double pp[2], ee[2], dedp[2][2];
	int ff1 = 0, ff2 = 0, fit = -1;

	/* Decide whether there is a solution in this triangle */
	j = k = 0;
	pp[0] = 0.0; pp[1] = 0.0;
	lchw_tri_Dp_sq(s, ee, vt, v, pp);
	if (ee[0] < 0.0) j++;
	if (ee[1] < 0.0) k++;

	pp[0] = 0.0; pp[1] = 1.0;
	lchw_tri_Dp_sq(s, ee, vt, v, pp);
	if (ee[0] < 0.0) j++;
	if (ee[1] < 0.0) k++;

	if (j != 1 || k != 1) {
		pp[0] = 1.0; pp[1] = 1.0;
		lchw_tri_Dp_sq(s, ee, vt, v, pp);
		if (ee[0] < 0.0) j++;
		if (ee[1] < 0.0) k++;
	
		/* Making this || filters out lots more for an avg itter of 0.74, */
		/* but has a failure rate of 1 in 50000. */
		/* This less stringent filter has an avg itter of 2.0 and 0 failure rate. */
		if ((j == 0 || j == 3) && (k == 0 || k == 3)) {
			return 0;
		}
	}

	pp[0] = 0.3333; pp[1] = 0.6667;

	for (i = 0; i < 30; i++) {
		double det;

		lchw_tri_Dp_sq(s, ee, vt, v, pp);
		lchw_tri_DDp_sq(s, dedp, vt, v, pp);

		/* Correct the point using inverse of dedp */
		det = (dedp[0][0] * dedp[1][1] - dedp[0][1] * dedp[1][0]);
		if (fabs(det) < 1e-20)
			break;			/* Hmm. */

		det = 1.0/det;
		pp[0] -= det * ( dedp[1][1] * ee[0] - dedp[0][1] * ee[1]);
		pp[1] -= det * (-dedp[1][0] * ee[0] + dedp[0][0] * ee[1]);

		/* If we're sufficiently close to zero point */
		if (fabs(ee[0]) < 1e-6 && fabs(ee[1]) < 1e-6)
			break;

#ifdef NEVER
#define THR 0.25
		/* If we're too far out of bounds, give up */
		/* (Speeds things up by about 40% at the cost of failing */
		/*  some that would suceed.) */
		if (i >= 2 && (pp[0] < -THR || pp[0] > (1.0 + THR) || pp[1] < -THR || pp[1] > (1.0 + THR) 
		  || pp[1] < (pp[0]-THR))) {
			return 0;
		}
#undef THR
#endif
	}

	lchw_tri_Dp_sq(s, ee, vt, v, pp);

	if (fabs(ee[0]) > 1e-6 || fabs(ee[1]) > 1e-6
	 || pp[0] < -EPS || pp[1] < (pp[0]-EPS) || pp[1] > (1.0 + EPS)) {
		return 0;
	}

	/* Return solution (output space) */
	for (f = 0; f < fdi; f++) {
		vv[f] = (v[0][f] - v[1][f]) * pp[0]
		      + (v[1][f] - v[2][f]) * pp[1]
		      +  v[2][f];
	}

	/* Return solution (simplex parameter space) */
	p[0] = pp[0];
	p[1] = pp[1];

	return 1;
}

/* -------------------------------------------------------- */
/* Cell/simplex object lower level code */

/* Utility to get or calculate a vertexes ink limit value */
static double get_limitv(
schbase *b,			/* Base search information */
int ix,				/* fwd index of cell */
float *fcb,			/* Pointer to base of vertex value array (ix is used if NULL) */
double *p			/* Array of input values (can be NULL to compute) */
) {
	rspl *s = b->s;
	float *base = fcb;
	double lv;
	if (base == NULL)
		base = s->g.a + ix * s->g.pss;
	lv = base[-1];					/* Fetch existing ink limit function value */
	if ((float)lv == L_UNINIT) {			/* Not been computed yet */
		if (p != NULL) {
			lv = INKSCALE * s->limitf(s->lcntx, p);	/* Do it */
			base[-1] = (float)lv;
		} else {
			int e, di = s->di;
			double pp[MXRI];			/* Copy from float to double */
			int tix;					/* Temp fwd cell index */

			for (tix = ix, e = 0; e < di; e++) {
				int dix;
				dix = tix % s->g.res[e];
				tix /= s->g.res[e];
				pp[e] = s->g.l[e] + (double)dix * s->g.w[e];	/* Base point */
			}
			lv = INKSCALE * s->limitf(s->lcntx, pp);	/* Do it */
			base[-1] = (float)lv;
		}
		s->g.limitv_cached = 1;			/* At least one limit value is cached */
	}
	return lv;
}

/* Utility to invalidate all the ink limit values */
/* cached in the main rspl array */
static void clear_limitv(
rspl *s
) {
	int i;
	float *gp;		/* Grid point pointer */

	if (s->g.limitv_cached != 0) {	/* If any have been set */
		/* Unset them all */
		for (i = 0, gp = s->g.a; i < s->g.no; i++, gp += s->g.pss) {
			gp[-1] = L_UNINIT;
		}
		s->g.limitv_cached = 0;
	}
}

/* Cell code */

static void free_cell_contents(fxcell *c);
static fxcell *cache_fxcell(revcache *r, int ix, int force);
static void uncache_fxcell(revcache *r, fxcell *cp);

/* Return a pointer to an appropriate fxcell cache structure. */
/* None of the sub simplex lists will be initialised. */
/* NOTE: must unget_cell() (== uncache_fxcell()) when fxcell */
/* is no longer needed */
/* Return NULL if we ran out of room in the cache. */
static fxcell *get_fxcell(
schbase *b,			/* Base search information */
int ix,				/* fwd index of cell */
int force			/* if nz, force memory allocation, so that we have at least one cell */
) {
	rspl *s = b->s;
	int ee, e, di = s->di;
	int p2di = (1<<di);
	int ff, f, fdi = s->fdi;
	fxcell *c;

	c = cache_fxcell(s->rev.cache, ix, force);	/* Fetch it from the cache and lock it */
	if (c == NULL)
		return NULL;

	if (!(c->flags & CELL_FLAG_1)) {			/* Have to (re)initialize cell & simplexes */
		int tix;								/* Temp fwd cell index */
		float *fcb = s->g.a + ix * s->g.pss;	/* Pointer to base float of fwd cell */

		/* Compute basic Cell info and vertex output values */
		for (ee = 0; ee < p2di; ee++) {
			float *vp = fcb + s->g.fhi[ee];
			for (f = 0; f < fdi; f++)		/* Transfer cell verticy values from grid */
				c->v[ee][f] = vp[f];

			/* ~~ reset any other cell info that will be stale */
		}

		/* Convert from cell index, to absolute fwd coord base values */
		c->limmin = INF_DIST;				/* and min/max values */
		c->limmax = -INF_DIST;
		for (tix = ix, e = 0; e < di; e++) {
			int dix;
			dix = tix % s->g.res[e];
			tix /= s->g.res[e];
			c->p[0][e] = s->g.l[e] + (double)dix * s->g.w[e];	/* Base point */
		}
		if (s->limitf != NULL) {			/* Compute ink limit values at base verticy */
			double lv = get_limitv(b, ix, fcb, c->p[0]); /* Fetch or generate limit value */
			c->v[0][fdi] = lv;
			if (lv < c->limmin)	/* And min/max for this cell */
				c->limmin = lv;
			if (lv > c->limmax)
				c->limmax = lv;
		}
			
		/* Setup cube verticy input position values, and ink limit values */
		for (ee = 1; ee < p2di; ee++) {
			for (e = 0; e < di; e++) {
				c->p[ee][e] = c->p[0][e];
				if (ee & (1 << e))
					c->p[ee][e] += s->g.w[e];		/* In input space offset */
			}
			if (s->limitf != NULL) {			/* Compute ink limit values at cell verticies */
				double lv = get_limitv(b, ix, fcb + s->g.fhi[ee], c->p[ee]);
				c->v[ee][fdi] = lv;
				if (lv < c->limmin)	/* And min/max for this cell */
					c->limmin = lv;
				if (lv > c->limmax)
					c->limmax = lv;
			}
		}
		
		/* Compute the output bounding group for fast rejection testing */
		{
			double *vp[POW2MXRI];

			/* Make array of pointers to double vectors */
			for (ee = 0; ee < p2di; ee++)
				vp[ee] = c->v[ee];

			nn_grpinit(s, &c->g, vp, p2di, NULL);
		}
		c->flags = CELL_FLAG_1;
	}

	return c;
}

void free_simplex_info(fxcell *c, int dof);

/* Free up any allocated simplexes in a cell, */
/* and set the pointers to NULL. */
/* Nothing else is changed (ie. it's NOT removed from */
/* the cache index or unthrheaded from the mru list). */
static void
free_cell_contents(
fxcell *c
) {
	int nsdi;
	
	/* Free up all the simplexes */
	if (c->s != NULL) {
		for (nsdi = 0; nsdi <= c->s->di; nsdi++) {
			if (c->sx[nsdi] != NULL) {
				free_simplex_info(c, nsdi);
				c->sx[nsdi] = NULL;
			}
		}
	}
	/* ~~ free any other cell information */
}

/* - - - - - -  */
/* Simplex code */

/* Simplex and Cell hash index size increments */
int primes[] = {
	367,
	853,
	1489,
	3373,
	6863,
	12919,
	23333,
	43721,
	97849,
	146221,
	254941,
	407843,
	756869,
	999983,
	-1
};

/* Compute a simplex hash index */
unsigned int simplex_hash(revcache *rc, int sdi, int efdi, int *vix) {
	unsigned int hash = 0;
	int i;

	for (i = 0; i <= sdi; i++)
		hash = hash * 17 + vix[i];
	hash = hash * 17 + sdi;
	hash = hash * 17 + efdi;

	hash %= rc->spx_hash_size;
	return hash;
}

/* Allocate and do the basic initialisation for a DOF list of simplexes */
void alloc_simplexes(
fxcell *c,
int nsdi			/* Non limited sub simplex dimensionality */
) {
	rspl *s = c->s;
	schbase *b = s->rev.sb;
	revcache *rc = s->rev.cache;
	int ee, e, di = s->di;
	int f, fdi = s->fdi;
	int lsdi;			/* Ink limited Sub-simplex sdi */
	int tsxno;			/* Total number of DOF simplexes */
	int nsxno;			/* Number of non-ink limited DOF simplexes */
	int si, so;			/* simplex index in and out */

	DBG(("Allocating level %d sub simplexes in cell %d\n",nsdi,c->ix));
	if (c->sx[nsdi] != NULL)
		error("rspl rev, internal, trying allocate already allocated simplexes\n");

	/* Figure out how many simplexes will be at this nsdi */
	lsdi = nsdi + 1;	/* Ink limit simplexes sdi */

	tsxno = nsxno = s->rev.sspxi[nsdi].nospx;

	if (s->limitf != NULL && lsdi <= di)
		tsxno += s->rev.sspxi[lsdi].nospx;		/* Second set with extra input dimension */

	/* Make sure there is enough space in temp simplex filter list */
	if (b->lsxfilt < tsxno) {	/* Allocate more space if needed */

		if (b->lsxfilt > 0) {	/* Free old space before allocating new */
			free(b->sxfilt);
			DECSZ(b->s, b->lsxfilt * sizeof(char));
		}
		b->lsxfilt = 0;
		/* Allocate enough space for all the candidate cells */
		if ((b->sxfilt = (char *)rev_malloc(s, tsxno * sizeof(char))) == NULL)
			error("rev: malloc failed - temp simplex filter list, count %d",tsxno);
		b->lsxfilt = tsxno;	/* Current allocated space */
		INCSZ(b->s, b->lsxfilt * sizeof(char));
	}
		
	/* Figure out the number of simplexes that will actually be needed */
	for (si = so = 0; si < tsxno; si++) {
		psxinfo *psxi = NULL;
		int *icomb, *offs;
		int sdi = nsdi;
		int efdi = fdi;
		int ssi = si;
		int isclip = 0;
		if (si >= nsxno) {				/* If limit boundary simplex */
			sdi++;						/* One more dimension */
			efdi++;						/* One more constraint */
			ssi -= nsxno;				/* In second half of list */
			isclip++;					/* Limit clipped simplex */
		}
		psxi = &s->rev.sspxi[sdi].spxi[ssi];
		icomb = psxi->icomb;
		offs  = psxi->offs;

		b->sxfilt[si] = 0;				/* Assume simplex won't be used */

		/* Check if simplex should be discared due to the ink limit */
		if (s->limitf != NULL) {
			double max = -INF_DIST;
			double min =  INF_DIST;

			/* Find the range of ink limit values covered by simplex */
			for (e = 0; e <= sdi; e++) {		/* For all the simplex verticies */
				int i = offs[e];
				double vv = c->v[i][fdi];		/* Ink limit value */
				if (vv < min)
					min = vv;
				if (vv > max)
					max = vv;
			}
			
//if ((max - min) > EPS) printf("~1 Found simplex sdi %d, efdi %d, min = %f, max = %f, limitv = %f\n", sdi, efdi, min,max,s->limitv);
			if (isclip) {	/* Limit clipped simplex */
				/* (Make sure it straddles the limit boundary) */
				if (max <= s->limitv || min > s->limitv)
					continue;		/* Discard this simplex - it can't straddle the ink limit */
//printf("~1 using sub simplex sdi %d, efdi %d, min = %f, max = %f, limitv = %f\n", sdi, efdi, min,max,s->limitv);
			} else {
				if (min > s->limitv)
					continue;		/* Discard this simplex - it is above the ink limit */
			}
		}

		b->sxfilt[si] |= 1;		/* This cell will be OK */
		so++;
	}

	DBG(("There are %d level %d sub simplexes\n",so, nsdi));
	/* Allocate space for all the DOF simplexes that will be used */
	if (so > 0) {
		if ((c->sx[nsdi] = (simplex **) rev_calloc(s, so, sizeof(simplex *))) == NULL)
			error("rspl malloc failed - fxcell simplexes - list of pointers");
		INCSZ(s, so * sizeof(simplex *));
	}

	/* Setup SPLX_FLAG_1 level information in the simplex */
	for (si = so = 0; si < tsxno; si++) {
		simplex *x;
		psxinfo *psxi = NULL;
		int *icomb;
		int sdi, efdi;
		int ssi;
		int vix[MXRI+1];            /* fwd cell vertex indexes of this simplex [sdi+1] */

		if (b->sxfilt[si] == 0)		/* Decided not to use this one */
			continue;

#ifdef STATS
		s->rev.st[b->op].sinited++;
#endif /* STATS */

		sdi = nsdi;
		efdi = fdi;
		ssi = si;
		if (si >= nsxno) {				/* If limit boundary simplex */
			sdi++;						/* One more dimension */
			efdi++;						/* One more constraint */
			ssi -= nsxno;				/* In second half of list */
		}

		psxi = &s->rev.sspxi[sdi].spxi[ssi];
		icomb = psxi->icomb;

		/* Compute simplex vertexes so we can match it in the cache */
		for (e = 0; e <= sdi; e++) 
			vix[e] = c->ix + s->g.hi[psxi->offs[e]];

		x = c->sx[nsdi][so];

		/* If this is a shared face simplex, see if we already have it in another fxcell */
		if (x == NULL && psxi->face) {
			unsigned int hash;
//printf("~1 looking for existing simplex nsdi = %d\n",nsdi);
			hash = simplex_hash(rc, sdi, efdi, vix);
			for (x = rc->spxhashtop[hash]; x != NULL; x = x->hlink) {
				if (x->sdi != sdi
				 || x->efdi != efdi)
					continue;			/* miss */
				for (e = 0; e <= sdi; e++) {
					if (x->vix[e] != vix[e])
						break;			/* miss */
				}
				if (e > sdi)
					break;				/* hit */
			}
			if (x != NULL) {
				x->refcount++;
//printf("~1 found hit in simplex face list hash %d, refcount = %d\n",hash,x->refcount);
			}
		}
		/* Doesn't already exist */
		if (x == NULL) {
			if ((x = (simplex *) rev_calloc(s, 1, sizeof(simplex))) == NULL)
				error("rspl malloc failed - fxcell simplexes - base simplex %d bytes",sizeof(simplex));
			INCSZ(s, sizeof(simplex));
			x->refcount = 1;
			x->touch = s->rev.stouch-1;
			x->flags = 0;

			if (si >= nsxno) {				/* If limit boundary simplex */
				x->flags |= SPLX_CLIPSX;	/* Limit clipped simplex */
			}

			/* Fill in the other simplex details */
			x->s    = s;					/* Parent rspl */
			x->ix   = c->ix;				/* Construction cube base index */
			for (e = 0; e <= sdi; e++)		/* Indexs of fwd verticies that make up this simplex */
				x->vix[e] = vix[e];
			x->psxi = psxi;					/* Pointer to constant per simplex info */
//printf("~1 set simplex 0x%x psxi = 0x%x\n",x,x->psxi);
			x->si   = so;					/* Diagnostic, simplex offset in list */
			x->sdi  = sdi;					/* Copy of simplex dimensionaity */
			x->efdi = efdi;					/* Copy of effective output dimensionality */

			/* Copy cell simplex vertex output and limit values */
			for (e = 0; e <= sdi; e++) {		/* For all the simplex verticies */
				int i = x->psxi->offs[e];

				for (f = 0; f <= fdi; f++)		/* Copy vertex value + ink sum */
					x->v[e][f] = c->v[i][f];

				/* Setup output bounding box values (the hard way) */
				if (e == 0) {						/* Init to first vertex of simplex */
					for (f = 0; f <= fdi; f++)		/* Output space */
						x->min[f] = x->max[f] = c->v[i][f];
				} else {
					for (f = 0; f <= fdi; f++) {	/* Output space + ink sum */
						double vv;
//						if (f == fdi && s->limit == NULL)
//							continue;			/* Skip ink */
						vv = c->v[i][f];
						if (vv < x->min[f])
							x->min[f] = vv;
						else if (vv > x->max[f])
							x->max[f] = vv;
					}
				}
			}
			/* Add a margin */
			for (f = 0; f <= fdi; f++) {	/* Output space + ink sum */
				x->min[f] -= EPS;
				x->max[f] += EPS;
			}

			/* Setup input bounding box value pointers (the easy way) */
			for (ee = 0; ee < di; ee++) {
				x->p0[ee]   = c->p[0][ee];		/* Construction base cube origin */
				x->pmin[ee] = c->p[x->psxi->pmino[ee]][ee] - EPS;
				x->pmax[ee] = c->p[x->psxi->pmaxo[ee]][ee] + EPS;
			}

			x->flags |= SPLX_FLAG_1;		/* vv & iv done, nothing else */

			x->aloc2 = x->aloc5 = NULL;		/* Matrix allocations not done yet */

			/* Add it to the shared face simplex hash index */
			if (x->psxi->face) {
				unsigned int hash;
				int i;
				/* See if we should re-size the simplex hash index */
				if (++rc->nspx > (HASH_FILL_RATIO * rc->spx_hash_size)) {
					for (i = 0; primes[i] > 0 && primes[i] <= rc->spx_hash_size; i++)
						;
					if (primes[i] > 0) {
						int spx_hash_size = rc->spx_hash_size;	/* Old */
						simplex **spxhashtop = rc->spxhashtop;

						rc->spx_hash_size = primes[i];

						DBG(("Increasing face simplex hash index to %d\n",spx_hash_size));
//printf("~1 increasing simplex hash index size to %d\n",spx_hash_size);
						/* Allocate a new index */
						if ((rc->spxhashtop = (simplex **) rev_calloc(s, rc->spx_hash_size,
						                                   sizeof(simplex *))) == NULL)
							error("rspl malloc failed - reverse simplex cache index");
						INCSZ(s, rc->spx_hash_size * sizeof(simplex *));

						/* Transfer all the simplexes to the new index */
						for (i = 0; i < spx_hash_size; i++) {
							simplex *x, *nx;
							for (x = spxhashtop[i]; x != NULL; x = nx) {
								nx = x->hlink;
								hash = simplex_hash(rc, x->sdi, x->efdi, x->vix);	/* New hash */
								x->hlink = rc->spxhashtop[hash];	/* Add to new hash index */
								rc->spxhashtop[hash] = x;
							}
						}
						free(spxhashtop); /* Done with old index */
						DECSZ(s, spx_hash_size * sizeof(simplex *));
					}
				}
				hash = simplex_hash(rc, sdi, efdi, vix);

				/* Add this to hash index */
				x->hlink = rc->spxhashtop[hash];
				rc->spxhashtop[hash] = x;
//printf("~1 Added simplex to hash %d, rc->nspx = %d\n",hash,rc->nspx);
			}

//if (rc->nunlocked == 0 && rc->s->rev.sz > rc->s->rev.max_sz)
//printf("~1 unable to decrease_revcache 1\n");

			/* keep memory in check */
			while (rc->nunlocked > 0 && rc->s->rev.sz > rc->s->rev.max_sz) {
				if (decrease_revcache(rc) == 0)
					break;
			}
		}
		c->sx[nsdi][so] = x;
		so++;
	}
	c->sxno[nsdi] = so;				/* Record actual number in list */
	c->flags |= CELL_FLAG_2;		/* Note that cell now has simplexes */
}

/* Free up any allocated for a list of sub-simplexes */
void
free_simplex_info(
fxcell *c,
int nsdi			/* non limit sub simplex dimensionaity */
) {
	int si, sxno = c->sxno[nsdi];	/* Number of simplexes */

	for (si = 0; si < sxno; si++) { /* For all the simplexes */
		simplex *x = c->sx[nsdi][si];
		int dof = x->sdi - x->efdi;

//printf("~1 freeing simplex, refcount = %d\n",x->refcount);
		if (--x->refcount <= 0) {		/* Last reference to this simplex */

//printf("~1 freeing simplex 0x%x psxi = 0x%x\n",x,x->psxi);
			if (x->psxi->face) {
				unsigned int hash;
				revcache *rc = c->s->rev.cache;
				
				hash = simplex_hash(rc, x->sdi, x->efdi, x->vix);

				/* Free it from the hash list */
				if (rc->spxhashtop[hash] == x) {
					rc->spxhashtop[hash] = x->hlink;
					rc->nspx--;
//printf("~1 removed simplex from hash %d, nspx now = %d\n",hash,rc->nspx);
				} else {
					simplex *xx;
					for (xx = rc->spxhashtop[hash]; xx != NULL && xx->hlink != x; xx = xx->hlink)
						;
					if (xx != NULL) {		/* Found it */
						xx->hlink = x->hlink;
						rc->nspx--;
//printf("~1 removed simplex from hash %d, nspx now = %d\n",hash,rc->nspx);
					}
//else
//printf("~1 warning, failed to find face simplex hash %d, sdi = %d in cache index (nspx = %d)!!\n",hash,x->sdi,rc->nspx);
				}
			}
			if (x->aloc2 != NULL) {
				int adof = dof >= 0 ? dof : 0;		/* Allocation dof */
				int asize;
				if (dof == 0)
					asize = sizeof(double) * (x->efdi * x->sdi)
				          + sizeof(double *) * x->efdi 
				          + sizeof(int) * x->sdi;
				else
					asize = sizeof(double) * (x->sdi * (x->efdi + x->sdi + adof + 2) + x->efdi)
				          + sizeof(double *) * (x->efdi + 2 * x->sdi);
				free(x->aloc2);
				DECSZ(x->s, asize);
			}

			if (x->aloc5 != NULL) {
				int asize;
				if (x->naux == dof)
					asize = sizeof(double *) * x->naux
				          + sizeof(double) * (x->naux * dof)
				          + sizeof(int) * dof;
				else
					asize = sizeof(double *) * (x->naux + dof) 
					      + sizeof(double) * (dof * (x->naux + dof + 1));
				free(x->aloc5);
				DECSZ(x->s, asize);
			}

			/* ~~ free any other simplex information */

			free(x);
			DECSZ(c->s, sizeof(simplex));
			c->sx[nsdi][si] = NULL;
		}
	}
	free(c->sx[nsdi]);
	DECSZ(c->s, c->sxno[nsdi] * sizeof(simplex *));
	c->sx[nsdi] = NULL;
	c->sxno[nsdi] = 0;

	/* ~~ free any other cell information */
}

/* - - - - - - - - - - - - */
/* Check that an input space vector is within a given simplex, */
/* and that it meets any ink limit. */
/* Return zero if outside the simplex, */
/* 1 normally if within the simplex, */
/* and 2 if it would be over the ink limit if limit was enabled. */
static int
within_simplex(
simplex *x,				/* Simplex */
double *p				/* Input coords in simplex space */
) {
	rspl *s = x->s;
	int    fdi = s->fdi;
	int e, sdi = x->sdi;		/* simplex dimensionality */
	double cp, lp;
	int rv = 1;
	/* EPS is allowance for numeric error */
	/* (Don't want solutions falling down */
	/* the numerical cracks between the simplexes) */

	/* Check we are within baricentric limits */
	for (lp = 0.0, e = 0; e < sdi; e++) {
		cp = p[e];
		if ((cp+EPS) < lp) 		/* Outside baricentric or not in correct */
			return 0;			/* order for this simplex  */
		lp = cp;
	}
	if ((1.0+EPS) < lp)  		/* outside baricentric range */
		return 0;

	/* Compute limit using interp. - assume simplex would have been trivially rejected */
	if (s->limitf != NULL) {
		double sum = 0.0;			/* Might be over the limit */
		for (e = 0; e < sdi; e++)
			sum += p[e] * (x->v[e][fdi] - x->v[e+1][fdi]);
		sum += x->v[sdi][fdi];
		if (sum > s->limitv) {
			if (s->limiten != 0)
	 			return 0;			/* Exceeds ink limit */
			else
				rv = 2;				/* would have exceeded limit */
		}
	}

#ifdef NEVER
	/* Constrain to legal values */
	/* (Is this needed ?????) */
	for (e = 0; e < sdi; e++) {
		cp = p[e];
		if (cp < 0.0)
			p[e] = 0.0;
		else if (cp > 1.0)
			p[e] = 1.0;
	}
#endif
	return rv;
}

/* Check that an input space vector of a simplex meets the ink limit. */
/* Return zero if outside the simplex, */
/* 1 normally if within the simplex, */
/* and 2 if it would be over the ink limit if limit was enabled. */
/* This is the same as within_simplex() but only checks the ink limit. */
static int
within_simplex_limit(
simplex *x,				/* Simplex */
double *p				/* Input coords in simplex space */
) {
	rspl *s = x->s;
	int    fdi = s->fdi;
	int e, sdi = x->sdi;		/* simplex dimensionality */
	int rv = 1;

	/* Compute limit using interp. - assume simplex would have been trivially rejected */
	if (s->limitf != NULL) {
		double sum = 0.0;			/* Might be over the limit */
		for (e = 0; e < sdi; e++)
			sum += p[e] * (x->v[e][fdi] - x->v[e+1][fdi]);
		sum += x->v[sdi][fdi];
		if (sum > s->limitv) {
			if (s->limiten != 0)
	 			return 0;			/* Exceeds ink limit */
			else
				rv = 2;				/* would have exceeded limit */
		}
	}
	return rv;
}

/* Similar check to within_simplex(), but with explicit simplex definition */
/* and no ink limit check. Returns 0 if outside, 1 if within */
static int
simple_within_simplex(
double v[MXRI+1][MXRO],	/* Vertex values */
double *p,				/* Input coords in simplex space */
int sdi					/* input dimensionality of simplex */
) {
	int e;
	double cp, lp;

	/* Check we are within baricentric limits */
	for (lp = 0.0, e = 0; e < sdi; e++) {
		cp = p[e];
		if ((cp+EPS) < lp) 		/* Outside baricentric or not in correct */
			return 0;			/* order for this simplex  */
		lp = cp;
	}
	if ((1.0+EPS) < lp)  		/* outside baricentric range */
		return 0;

	return 1;
}

/* Convert vector from simplex space to absolute cartesian space */
static void simplex_to_abs(
simplex *x,
double *out,	/* output in absolute space */
double *in		/* Input in simplex space */
) {
	rspl *s     = x->s;
	int e, di   = s->di;
	int *icomb  = x->psxi->icomb;	/* Coord combination order */

	for (e = 0; e < di; e++) {		/* For each absolute coord */
		double ov = x->p0[e];		/* Base value */
		int ee = icomb[e];			/* Simplex param index */
		if (ee >= 0)				/* Simplex param value */
			ov += s->g.w[e] * in[ee];
		else if (ee == -2)			/* 1 value */
			ov += s->g.w[e];
									/* Else 0 value */
		out[e] = ov;
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Given the parametric clip line equation, compute the */
/* implicit equation in terms of the absolute output space. */
/* Pad equation with target ink limit in case it is use */
/* with CLIPSX sub-simplexes. */
/* Note that no line equation values are returned if fdi = 1, */
/* since there is no such thing as an implicit line equation. */
/* (Re-usable version for lines in general) */
static void
init_line_eq_imp(
rspl *s,
schbase *b,			/* to set cdir, may be NULL if not needed. */
double ***pcla,		/* pointer to clip vector LHS implicit equation matrix */
double clb[MXRO+1],	/* Clip vector RHS implicit equation vector */
double st[MXRO],	/* Start point */
double de[MXRO],	/* Delta */
int inkeq			/* nz to add ink limit target equation if s->limitf != NULL */
) {
	int ff, f, fdi = s->fdi;
	int i, p;
	double lgst;
	double **cla = *pcla;

	DBG(("Computing clipping line implicit equation, dim = %d\n", fdi));
	
	/* Pick a pivot element */
	for (lgst = -1.0, p = -1, f = 0; f < fdi; f++) {
		double tt = de[f];
		if (b != NULL)
			b->cdir[f] = tt;		/* Stash this away */
		tt = fabs(tt);
		if (tt > lgst) {
			lgst = tt;
			p = f;
		}
	}
	if (p < 0)	/* Shouldn't happen */
		error("rspl rev, internal, trying to cope with zero length clip line\n");
	
	if (cla == NULL) {
		cla = dmatrix(0, fdi-1, 0, fdi);	/* Allow for ink limit supliment */
		*pcla = cla;
	}

	for (i = ff = 0;  ff < fdi; ff++) {	/* For the input rows */
		if (ff == p) {
			continue;					/* Skip pivot row */
		}
		for (f = 0; f < fdi; f++) {		/* For input & output columns */
			if (f == p) {
				cla[i][f] = -de[ff];	/* Last column is -ve delta value */
			} else if (f == ff) {
				cla[i][f] = de[p];	/* Diagonal is pivot value */
			} else {
				cla[i][f] = 0.0;		/* Else zero */
			}
		}
		clb[i] = de[p] * st[ff] - de[ff] * st[p];
		i++;
	}

	/* Add ink limit target equation - */
	/* interpolated ink value == target */
	if (inkeq && s->limitf != NULL) {
		for (i = 0;  i < (fdi-1); i++)
			cla[i][fdi] = 0.0;

		for (f = 0; f < fdi; f++) 
			cla[fdi-1][f] = 0.0;
		
		cla[fdi-1][fdi] = 1.0;
		clb[fdi-1] = s->limitv;
	}

#ifdef NEVER
/* Verify that the implicit equation is correct */
{
	double pnt[MXRO], v[MXRO];
	double pa;	/* Parameter */
	for (pa = 0.0; pa <= 1.0; pa += 0.125) {
		for (f = 0; f < fdi; f++) {
			pnt[f] = st[f] + pa * de[f];
		}

		/* Verify the implicit equation */
		for (ff = 0; ff < (fdi-1); ff++) {
			v[ff] = 0.0;
			for (f = 0; f < fdi; f++) {
				v[ff] += cla[ff][f] * pnt[f];
			}
			v[ff] -= clb[ff];
			if (v[ff] < 0.0)
				v[ff] = -v[ff];
			if (v[ff] > 0.000001) {
				printf("Point on clip line = %f %f %f\n",pnt[0],pnt[1],pnt[2]);
				printf("Implicit %d error of = %f\n",ff, v[ff]);
			}
		}
	}
}
#endif /* NEVER */

}

/* Version of above used to set vector clipping line up */
static void
init_line_eq(
schbase *b,
double st[MXRO],	/* Start point */
double de[MXRO]		/* Delta */
) {
	DBG(("Computing clipping line implicit equation, dim = %d\n", b->s->fdi));
	
	init_line_eq_imp(b->s, b, &b->cla, b->clb, st, de, 1); 
}

/* - - - - - -  */
/* Simpex solution info #2 */

/* Create the LU or SVD decomp needed to compute solution or locus. */
/* Return non-zero if it cannot be created */
static int
add_lu_svd(simplex *x) {

	if (x->flags & SPLX_FLAG_2F) {		/* Previously failed */
		return 1;
	}
	if (!(x->flags & SPLX_FLAG_2)) {
		int ee, e, sdi = x->sdi; 
		int f, efdi = x->efdi; 
		int dof = sdi-efdi;		/* Degree of freedom of locus, or -ve over specification */
		int adof = dof >= 0 ? dof : 0;		/* Allocation dof */
		int i;

		if (x->aloc2 == NULL) {	/* Allocate space for matricies and arrays */
			/* Do this in one hit to minimise malloc overhead */
			if (dof == 0) {
				int i;
				char *mem;
				int asize = sizeof(double) * (efdi * sdi)
				          + sizeof(double *) * efdi 
				          + sizeof(int) * sdi;

				if ((x->aloc2 = mem = (char *) rev_malloc(x->s, asize)) == NULL)
					error("rspl malloc failed - fxcell sub-simplex matricies");
				INCSZ(x->s, asize);

				/* Allocate biggest to smallest (double, pointers, ints) */
				/* to make sure that items lie on the natural boundaries. */

				/* Reserve matrix doubles */
				mem += efdi * sdi * sizeof(double);

				/* Allocate pointers */
				x->d_u = (double **)mem, mem += efdi * sizeof(double *);

				/* Allocate ints */
				x->d_w = (double *)mem, mem += sdi * sizeof(int);

#ifdef DEBUG
				if (mem != (x->aloc2 + asize))
					error("~1 aloc2a assert failed! Is %d, should be %d\n",mem - x->aloc2,asize);
#endif /* DEBUG */

				/* Reset and allocate matrix doubles */
				mem = x->aloc2; 
				for (i = 0; i < efdi; i++)
					x->d_u[i] = (double *)mem,	mem += sdi * sizeof(double);

			} else {
				int i;
				char *mem;
				int asize = sizeof(double) * (sdi * (efdi + sdi + adof + 2) + efdi)
				          + sizeof(double *) * (efdi + 2 * sdi);

				if ((x->aloc2 = mem = (char *) rev_malloc(x->s, asize)) == NULL)
					error("rspl malloc failed - fxcell sub-simplex matricies");
				INCSZ(x->s, asize);

				/* Allocate biggest to smallest (double, pointers, ints) */
				/* to make sure that items lie on the natural boundaries. */

				/* Reserve matrix doubles */
				mem += sdi * (efdi + sdi + adof) * sizeof(double);

				/* Allocate doubles */
				x->lo_xb = (double *)mem, mem += efdi * sizeof(double);
				x->lo_bd = (double *)mem; mem += sdi * sizeof(double);
				x->d_w = (double *)mem, mem += sdi * sizeof(double);

				/* Allocate pointers */
				x->d_u = (double **)mem, mem += efdi * sizeof(double *);
				x->d_v = (double **)mem, mem += sdi * sizeof(double *);
				x->lo_l = (double **)mem, mem += sdi * sizeof(double *);

#ifdef DEBUG
				if (mem != (x->aloc2 + asize))
					error("~1 aloc2b assert failed! Is %d, should be %d\n",mem - x->aloc2,asize);
#endif /* DEBUG */

				/* Reset and allocate matrix doubles */
				mem = x->aloc2;
				for (i = 0; i < efdi; i++)
					x->d_u[i] = (double *)mem,	mem += sdi * sizeof(double);
				for (i = 0; i < sdi; i++)
					x->d_v[i] = (double *)mem,	mem += sdi * sizeof(double);
				for (i = 0; i < sdi; i++)
					x->lo_l[i] = (double *)mem,	mem += adof * sizeof(double);

				/* Init any values that will be read before being written to. */
				for (f = 0; f < efdi; f++)
					x->lo_xb[f] = 1e100;		/* Silly value */
			}
		}

		/* Setup matrix from vertex values */
		for (f = 0; f < efdi; f++)
			for (e = 0; e < sdi; e++)
				x->d_u[f][e] = x->v[e][f] - x->v[e+1][f];

		if (dof == 0) {	/* compute LU */
			double rip;
#ifdef STATS
			x->s->rev.st[x->s->rev.sb->op].sinited2a++;
#endif /* STATS */
			if (lu_decomp(x->d_u, sdi, (int *)x->d_w, &rip)) {
				x->flags |= SPLX_FLAG_2F;	/* Failed */
				return 1;
			}
		} else {
//printf("~~ Creating SVD decomp, sdi = %d, efdi = %d\n", sdi, efdi);

#ifdef STATS
			x->s->rev.st[x->s->rev.sb->op].sinited2b++;
#endif /* STATS */
			if (svdecomp(x->d_u, x->d_w, x->d_v, efdi, sdi)) {
				x->flags |= SPLX_FLAG_2F;	/* Failed */
				return 1;
			}
	
			/* Threshold the singular values W[] */ 
			svdthresh(x->d_w, sdi);
	
			if (dof >= 0) {		/* If we expect a locus */
//printf("~~ got dif %d locus from SVD\n",dof);
				/* copy the locus direction coefficients out */
				for (i = e = 0; e < sdi; e++) {
					if (x->d_w[e] == 0.0) {		/* Found a zero W[] */
						if (i < dof) {
							for (ee = 0; ee < sdi; ee++) {	/* Copy column of V[][] */
								x->lo_l[ee][i] = x->d_v[ee][e];
							}
						}
						i++;
					}
				}
				if (i != dof) {
//printf("~~ got unexpected dof in svd\n");
					x->flags |= SPLX_FLAG_2F;	/* Failed */
					return 1;					/* Didn't get expected d.o.f. */
				}
			}
		}
		x->flags |= SPLX_FLAG_2;	/* Set flag so that it isn't attempted again */

//if (x->s->rev.cache->nunlocked == 0 && x->s->rev.sz > x->s->rev.max_sz)
//printf("~1 unable to decrease_revcache 2\n");

		/* keep memory in check */
		while (x->s->rev.cache->nunlocked > 0 && x->s->rev.sz > x->s->rev.max_sz) {
			if (decrease_revcache(x->s->rev.cache) == 0)
				break;
		}
	}
	return 0;
}

/* - - - - - -  */
/* Simplex solution info #4 */

/* Calculate the solution locus equation for this simplex and target */
/* (The direction was calculated by add_svd(), but now calculate */
/* the base solution point for this particular reverse lookup) */
/* Return non-zero if this point canot be calculated */
/* We are assuming that sdi > efdi */
static int
add_locus(
schbase *b,
simplex *x
) {
	int sdi = x->sdi; 
	int f, efdi = x->efdi; 
	int doback = 0;

#ifdef STATS
	x->s->rev.st[x->s->rev.sb->op].sinited4++;
#endif /* STATS */
	/* Use output of svdcmp() to solve overspecified and/or */
	/* singular equation A.x = b */

	/* Init the RHS B[] vector, and check if it doesn't match */
	/* that used to compute base value last time. */
	for (f = 0; f < efdi; f++) {
		double xb = b->v[f] - x->v[sdi][f];
		if (x->lo_xb[f] != xb) {
			x->lo_xb[f] = xb;
			doback = 1;			/* RHS differs, so re-compute */
		}
	}
	
#ifdef STATS
	if (doback && (x->flags & SPLX_FLAG_4))
		x->s->rev.st[x->s->rev.sb->op].sinited4i++;
#endif /* STATS */

	/* Compute locus */
	if (doback || !(x->flags & SPLX_FLAG_4))
		svdbacksub(x->d_u, x->d_w, x->d_v, x->lo_xb, x->lo_bd, efdi, sdi);
	
	x->flags |= SPLX_FLAG_4;

//if (x->s->rev.cache->nunlocked == 0 && x->s->rev.sz > x->s->rev.max_sz)
//printf("~1 unable to decrease_revcache 3\n");

	/* keep memory in check */
	while (x->s->rev.cache->nunlocked > 0 && x->s->rev.sz > x->s->rev.max_sz) {
		if (decrease_revcache(x->s->rev.cache) == 0)
			break;
	}

	return 0;
}

/* - - - - - -  */
/* Simplex solution info #5 */

/* Compute LU or SVD decomp of lo_l */
/* Allocates the memory for the various matricies */
/* Return non-zero if this canot be calculated. */
static int
add_auxil_lu_svd(
schbase *b,
simplex *x
) {
	int ee, sdi = x->sdi; 
	int f, efdi = x->efdi; 
	int dof = sdi-efdi;		/* Degree of freedom of locus */
	int naux = b->naux;		/* Number of auxiliaries actually available */

#ifdef STATS
	if (x->aaux != b->naux || x->auxbm != b->auxbm)
		x->s->rev.st[x->s->rev.sb->op].sinited5i++;
#endif /* STATS */

	if (x->aaux != b->naux) {	/* Number of auxiliaries has changed */
		if (x->aloc5 != NULL) {
			int asize;
			if (x->naux == dof)
				asize = sizeof(double *) * x->naux
			          + sizeof(double) * (x->naux * dof)
			          + sizeof(int) * dof;
			else
				asize = sizeof(double *) * (x->naux + dof) 
				      + sizeof(double) * (dof * (x->naux + dof + 1));
			free(x->aloc5);
			x->aloc5 = NULL;
			DECSZ(x->s, asize);
		}
		x->flags &= ~(SPLX_FLAG_5 | SPLX_FLAG_5F);	/* Force recompute */
	}
	
	if (x->auxbm != b->auxbm) {	/* Different selection of auxiliaries */
		x->flags &= ~(SPLX_FLAG_5 | SPLX_FLAG_5F);	/* Force recompute */
	}

	if (x->flags & SPLX_FLAG_5F) {		/* Previously failed */
		return 1;
	}
	if (!(x->flags & SPLX_FLAG_5)) {
		int *icomb = x->psxi->icomb; /* abs -> simplex coordinate translation */

		if (x->aloc5 == NULL) {	/* Allocate space for matricies and arrays */
			/* Do this in one hit to minimise malloc overhead */
			if (naux == dof) {
				int i;
				char *mem;
				int asize = sizeof(double *) * naux
				          + sizeof(double) * (naux * dof)
				          + sizeof(int) * dof;

				if ((x->aloc5 = mem = (char *) rev_malloc(x->s, asize)) == NULL)
					error("rspl malloc failed - fxcell sub-simplex matricies");
				INCSZ(x->s, asize);

				/* Allocate biggest to smallest (double, pointers, ints) */
				/* to make sure that items lie on the natural boundaries. */

				/* Reserve matrix doubles */
				mem += naux * dof * sizeof(double);

				/* Allocate pointers and ints */
				x->d_u = (double **)mem, mem += naux * sizeof(double *);
				x->d_w = (double *)mem, mem += dof * sizeof(int);

#ifdef DEBUG
				if (mem != (x->aloc5 + asize))
					error("aloc5a assert failed! Is %d, should be %d\n",mem - x->aloc5,asize);
#endif /* DEBUG */

				/* Reset and allocate matrix doubles */
				mem = x->aloc5;
				for (i = 0; i < naux; i++)
					x->d_u[i] = (double *)mem,	mem += dof * sizeof(double);
			} else {
				int i;
				char *mem;
				int asize = sizeof(double *) * (naux + dof) 
				          + sizeof(double) * (dof * (naux + dof + 1));

				if ((x->aloc5 = mem = (char *) rev_malloc(x->s, asize)) == NULL)
					error("rspl malloc failed - fxcell sub-simplex matricies");
				INCSZ(x->s, asize);

				/* Allocate biggest to smallest (double, pointers, ints) */
				/* to make sure that items lie on the natural boundaries. */

				/* Reserve matrix doubles */
				mem += dof * (naux + dof) * sizeof(double);

				/* Allocate doubles */
				x->ax_w = (double *)mem, mem += dof * sizeof(double);

				/* Allocate pointers, ints */
				x->ax_u = (double **)mem, mem += naux * sizeof(double *);
				x->ax_v = (double **)mem, mem += dof * sizeof(double *);

#ifdef DEBUG
				if (mem != (x->aloc5 + asize))
					error("aloc5b assert failed! Is %d, should be %d\n",mem - x->aloc5,asize);
#endif /* DEBUG */

				/* Reset and allocate matrix doubles */
				mem = x->aloc5;
				for (i = 0; i < naux; i++)
					x->ax_u[i] = (double *)mem,	mem += dof * sizeof(double);
				for (i = 0; i < dof; i++)
					x->ax_v[i] = (double *)mem, mem += dof * sizeof(double);
			}
			x->aaux = naux;				/* Number of auxiliaries allocated for */
		}
	
		/* Setup A[][] matrix to decompose, and figure number of auxiliaries actually needed */
		for (ee = naux = 0; ee < b->naux; ee++) {
			int ei = icomb[b->auxi[ee]];		/* Simplex relative auxiliary index */
			if (ei < 0)
				continue;		/* aux corresponds with fixed input value for this simplex */
			for (f = 0; f < dof; f++)
				x->ax_u[naux][f] = x->lo_l[ei][f];
			naux++;
		}
		x->naux = naux;					/* Number of auxiliaries actually available */
		x->auxbm = b->auxbm;			/* Mask of auxiliaries used */

		if (naux == dof) {				/* Use LU decomp to solve exactly */
			double rip;

#ifdef STATS
			x->s->rev.st[x->s->rev.sb->op].sinited5a++;
#endif /* STATS */
			if (lu_decomp(x->ax_u, dof, (int *)x->ax_w, &rip)) {
				x->flags |= SPLX_FLAG_5F;
				return 1;
			}

		} else if (naux > 0) {			/* Use SVD to solve least squares */

#ifdef STATS
			x->s->rev.st[x->s->rev.sb->op].sinited5b++;
#endif /* STATS */
			if (svdecomp(x->ax_u, x->ax_w, x->ax_v, naux, dof)) {
				x->flags |= SPLX_FLAG_5F;
				return 1;
			}
	
			/* Threshold the singular values W[] */ 
			svdthresh(x->ax_w, dof);
		} /* else naux == 0, don't setup anything */

		x->flags |= SPLX_FLAG_5;

//if (x->s->rev.cache->nunlocked == 0 && x->s->rev.sz > x->s->rev.max_sz)
//printf("~1 unable to decrease_revcache 4\n");

		/* keep memory in check */
		while (x->s->rev.cache->nunlocked > 0 && x->s->rev.sz > x->s->rev.max_sz) {
			if (decrease_revcache(x->s->rev.cache) == 0)
				break;
		}
	}
	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Initialise a static sub-simplex verticy information table */
void
rspl_init_ssimplex_info(
rspl *s,
ssxinfo *xip,				/* Pointer to sub-simplex info structure to init. */
int sdi						/* Sub-simplex dimensionality (range 0 - di) */
) {
	int e, di = s->di;		/* Dimensionality */
	int vi, nospx;			/* Number of sub-simplexes */
	XCOMBO(vcmb, MXDI, sdi+1, 1 << di);/* Simplex dimension sdi out of cube dimension di counter */

	DBG(("init_ssimplex_info called with sdi = %d\n",sdi));
	/* First count the number of sub-simplexes */
	nospx = 0;
	XCB_INIT(vcmb);
	while (!XCB_DONE(vcmb)) {
		nospx++;
		XCB_INC(vcmb);
	}

	xip->sdi = sdi;
	xip->nospx = nospx;
	if ((xip->spxi = (psxinfo *) rev_calloc(s, nospx, sizeof(psxinfo))) == NULL)
		error("rspl malloc failed - fxcell sub-simplex info array");
	INCSZ(s, nospx * sizeof(psxinfo));
	
	DBG(("Number of subsimplex = %d\n",nospx));
	/* For all sub-simplexes */
	XCB_INIT(vcmb);
	for (vi = 0; vi < nospx; vi++) {
		psxinfo *x = &xip->spxi[vi];
		int i;
		int andm, orm;

		/* XCOMB generates verticies in order from max to min offset */

		/* Compute Absolute -> Parameter mapping */
		for (e = 0; e < di; e++) {				/* For each absolute axis */

			if ((vcmb[sdi] & (1<<e)) != 0) {
				x->icomb[e] = -2;	/* This abs is always '1' */

			} else if ((vcmb[0] & (1<<e)) == 0) {
				x->icomb[e] = -1;	/* This abs is always '0' */

			} else {
				for (i = 0; i < sdi; i++) {	/* For each verticy in large to small order (!first) */
					if ((vcmb[i]   & (1<<e)) != 0 && 
					    (vcmb[i+1] & (1<<e)) == 0) {/* Transition from offset 1 to 0 */
						x->icomb[e] = i;	/* This is parameter */
						break;
					}
				}
			}
		}
		
		/* Compute fwd grid offsets for each simplex vertex in baricentric order */
		for (i = 0; i <= sdi; i++) {	/* For each verticy */
			int pmin[MXRI], pmax[MXRI];
			x->offs[i]  = vcmb[i];
			x->goffs[i] = s->g.hi[vcmb[i]];
			x->foffs[i] = s->g.fhi[vcmb[i]];

			/* Setup input coordinate bounding box value offsets */
			if (i == 0) {								/* Init to first vertex of simplex */
				for (e = 0; e < di; e++) {				/* Input space */
					x->pmino[e] = x->pmaxo[e] = vcmb[i];
					pmin[e] = pmax[e] = vcmb[i] & (1<<e);
				}
			} else {
				for (e = 0; e < di; e++) {			/* Input space */
					int vv = vcmb[i] & (1<<e);
					if (vv < pmin[e]) {				/* Adjust min/max offsets */
						x->pmino[e] = vcmb[i];
						pmin[e] = vv;
					} else if (vv > pmax[e]) {
						x->pmaxo[e] = vcmb[i];
						pmax[e] = vv;
					}
				}
			}
		}

		/* See if the sub-simplex lies on a cube face */
		andm = ~0;
		orm = 0;
		for (i = 0; i <= sdi; i++) {	/* For each verticy */
			andm &= vcmb[i];
			orm  |= vcmb[i];
		}
		/* If one coordinate is common (all 0 or all 1) to the verticies, */
		/* they must all be on the same cube face. */
		if (andm != 0 || orm != ((1 << di)-1))
			x->face = 1;
		else
			x->face = 0;

#ifdef DEBUG
		printf("Verticies   = ");
		for (i = 0; i <= sdi; i++)
			printf("%d ",vcmb[i]);
		printf("\n");
		
		printf("Face        = %s\n",x->face ? "True" : "False");
		
		printf("Abs -> Parm = ");
		for (e = 0; e < di; e++)
			printf("%d ",x->icomb[e]);
		printf("\n");
		
		printf("Grid Offset      = ");
		for (e = 0; e <= sdi; e++)
			printf("%d ",x->goffs[e]);
		printf("Float Offset      = ");
		for (e = 0; e <= sdi; e++)
			printf("%d ",x->foffs[e]);
		printf("\n");
		printf("\n");
#endif /* DEBUG */

		/* Increment the counter value */
		XCB_INC(vcmb);
	}
}

/* Free the given sub-simplex verticy information */
void
rspl_free_ssimplex_info(
rspl *s,
ssxinfo *xip		/* Pointer to sub-simplex info structure */
) {
	if (xip == NULL)	/* Assert */
		return;

	free(xip->spxi);
	DECSZ(s, xip->nospx * sizeof(psxinfo));
	xip->spxi = NULL;
}

/* ====================================================== */
/* Reverse cell cache code                                */

/* Allocate and initialise the fxcell cache */
static revcache *
alloc_revcache(
rspl *s
) {
	revcache *rc;

	DBG(("alloc_revcache called\n"));
	if ((rc = (revcache *) rev_calloc(s, 1, sizeof(revcache))) == NULL)
		error("rspl malloc failed - fxcell cache");
	INCSZ(s, sizeof(revcache));
	
	rc->s = s;		/* For stats */

	/* Allocate an initial cell hash index */
	rc->cell_hash_size = primes[0];

	if ((rc->hashtop = (fxcell **) rev_calloc(s, rc->cell_hash_size, sizeof(fxcell *))) == NULL)
		error("rspl malloc failed - fxcell cache index");
	INCSZ(s, rc->cell_hash_size * sizeof(fxcell *));

	/* Allocate an initial simplex face match hash index */
	rc->spx_hash_size = primes[0];

	if ((rc->spxhashtop = (simplex **) rev_calloc(s, rc->spx_hash_size, sizeof(simplex *))) == NULL)
		error("rspl malloc failed - reverse simplex cache index");
	INCSZ(s, rc->spx_hash_size * sizeof(simplex *));

	return rc;
}

/* Free the fxcell cache */
static void
free_revcache(revcache *rc) {
	int i;
	fxcell *cp, *ncp;

	/* Free any stuff allocated in the cell contents, and the cell itself. */
	for (cp = rc->mrubot; cp != NULL; cp = ncp) {
		ncp = cp->mruup;
		free_cell_contents(cp);
		free(cp);
		DECSZ(rc->s, sizeof(fxcell));
	}

	/* Free the hash indexes */
	free(rc->hashtop);
	DECSZ(rc->s, rc->cell_hash_size * sizeof(fxcell *));
	free(rc->spxhashtop);
	DECSZ(rc->s, rc->spx_hash_size * sizeof(simplex *));

	DECSZ(rc->s, sizeof(revcache));
	free(rc);
}

/* Invalidate the whole cache */
static void
invalidate_revcache(
revcache *rc)
{
	int i;
	fxcell *cp;

	rc->nunlocked = 0;

	/* Free any stuff allocated in the cell contents */
	for (cp = rc->mrubot; cp != NULL; cp = cp->mruup) {
		free_cell_contents(cp);
		cp->refcount = 0;		/* Make sure they can now be reused */
		cp->ix = 0;
		cp->flags = 0;			/* Contents needs re-initializing */
		rc->nunlocked++;
	}

	/* Clear the hash table so they can't be hit */
	for (i = 0; i < rc->cell_hash_size; i++) {
		rc->hashtop[i] = NULL;
	}

}

#define HASH(xx, yy) ((yy) % xx->cell_hash_size)

/* Allocate another cell, and add it to the cache. */
/* This may re-size the hash index too. */
/* Return the pointer to the new cell. */
/* (Note it's not our job here to honour the memory limit) */
static fxcell *
increase_revcache(
revcache *rc
) {
	fxcell *nxcell;		/* Newly allocated fxcell */
	int i;

//	DBG(("Adding another cell to cache\n"));

#ifdef NEVER	/* We may be called with force != 0 */
	if (rc->s->rev.sz >= rc->s->rev.max_sz)
		return NULL;
#endif

	if ((nxcell = (fxcell *) rev_calloc(rc->s, 1, sizeof(fxcell))) == NULL)
		error("rspl malloc failed - reverse fxcells");
	INCSZ(rc->s, sizeof(fxcell));

	nxcell->s = rc->s;

	/* Add cell to the bottom of the cache mru linked list */
	if (rc->mrutop == NULL)					/* List was empty */
		rc->mrutop = nxcell;
	else {
		rc->mrubot->mrudown = nxcell;	/* Splice into bottom */
		nxcell->mruup = rc->mrubot;
	}
	rc->mrubot = nxcell;
	rc->nacells++;
	rc->nunlocked++;

//	DBG(("cache is now %d cells\n",rc->nacells));

	/* See if the hash index should be re-sized */
	if (rc->nacells > (HASH_FILL_RATIO * rc->cell_hash_size)) {
		for (i = 0; primes[i] > 0 && primes[i] <= rc->cell_hash_size; i++)
			;
		if (primes[i] > 0) {
			int cell_hash_size = rc->cell_hash_size;	/* Old */
			fxcell **hashtop = rc->hashtop;

			rc->cell_hash_size = primes[i];

			DBG(("Increasing cell cache hash index to %d\n",cell_hash_size));
			/* Allocate a new index */
			if ((rc->hashtop = (fxcell **) rev_calloc(rc->s, rc->cell_hash_size, sizeof(fxcell *))) == NULL)
				error("rspl malloc failed - fxcell cache index");
			INCSZ(rc->s, rc->cell_hash_size * sizeof(fxcell *));

			/* Transfer all the cells to the new index */
			for (i = 0; i < cell_hash_size; i++) {
				fxcell *c, *nc;
				for (c = hashtop[i]; c != NULL; c = nc) {
					int hash;
					nc = c->hlink;
					hash = HASH(rc, c->ix); 		/* New hash */
					c->hlink = rc->hashtop[hash];	/* Add to new hash index */
					rc->hashtop[hash] = c;
				}
			}

			/* Done with old index */
			free(hashtop);
			DECSZ(rc->s, cell_hash_size * sizeof(fxcell *));
		}
	}
	
	return nxcell;
}

/* Reduce the cache memory usage by freeing the least recently used unlocked cell. */
/* Return nz if we suceeeded in freeing some memory. */
static int decrease_revcache(
revcache *rc		/* Reverse cache structure */
) {
	int hit = 0;
	int hash;
	fxcell *cp;
	
	DBG(("Decreasing cell cache memory allocation by freeing a cell\n"));

	/* Use the least recently used unlocked fxcell */
	for (cp = rc->mrubot; cp != NULL && cp->refcount > 0; cp = cp->mruup)
		;

	/* Run out of unlocked cells */
	if (cp == NULL) {
		DBG(("Failed to find unlocked cell to free\n"));
//printf("~1 failed to decrease memory\n");
		return 0;
	}
	
	/* If it has been used before, free up the simplexes */
	free_cell_contents(cp);

	/* Remove from current hash index (if it is in it) */
	hash = HASH(rc,cp->ix);			/* Old hash */
	if (rc->hashtop[hash] == cp) {
		rc->hashtop[hash] = cp->hlink;
	} else {
		fxcell *c;
		for (c = rc->hashtop[hash]; c != NULL && c->hlink != cp; c = c->hlink)
			;
		if (c != NULL)
			c->hlink = cp->hlink;
	}

	/* Free up this cell - Remove it from LRU list */
	if (rc->mrutop == cp)
		rc->mrutop = cp->mrudown;
	if (rc->mrubot == cp)
		rc->mrubot = cp->mruup;
	if (cp->mruup != NULL)
		cp->mruup->mrudown = cp->mrudown;
	if (cp->mrudown != NULL)
		cp->mrudown->mruup = cp->mruup;
	cp->mruup = cp->mrudown = NULL;
	free(cp);
	DECSZ(rc->s, sizeof(fxcell));
	rc->nacells--;
	rc->nunlocked--;

	DBG(("Freed a rev fxcell\n"));
	return 1;
}

/* Return a pointer to an appropriate fxcell */
/* cache structure. cell->flags will be 0 if the fxcell */
/* has been reallocated. cell contents will be 0 if */
/* never used before. */
/* The cell reference count is incremented, so that it */
/* can't be thrown out of the cache. The cell must be */
/* released with uncache_fxcell() when it's no longer needed. */
/* return NULL if we ran out of room in the cache */
static fxcell *cache_fxcell(
revcache *rc,		/* Reverse cache structure */
int ix,				/* fwd index of cell */
int force			/* if nz, force memory allocation, so that we have at least one fxcell */
) {
	int hit = 0;
	int hash;
	fxcell *cp;
	
	/* keep memory in check - fail if we're out of memory and can't free any */
	/* (Doesn't matter if it might be a hit, it will get picked up the next time) */
	if (!force && rc->s->rev.sz > rc->s->rev.max_sz && rc->nunlocked <= 0) {
		return NULL;
	}

//if (rc->nunlocked == 0 && rc->s->rev.sz > rc->s->rev.max_sz)
//printf("~1 unable to decrease_revcache 5\n");

	/* Free up memory to get below threshold */
	while (rc->nunlocked > 0 && rc->s->rev.sz > rc->s->rev.max_sz) {
		if (decrease_revcache(rc) == 0)
			break;
	}

	hash = HASH(rc,ix);		/* Compute hash of fwd cell index */

	/* See if we get a cache hit */
	for (cp = rc->hashtop[hash]; cp != NULL; cp = cp->hlink) {
		if (ix == cp->ix) {	/* Hit */
			hit = 1;
#ifdef STATS
			rc->s->rev.st[rc->s->rev.sb->op].chits++;
#endif /* STATS */
			break;
		}
	}
	if (!hit) {			/* No hit, use new cell or the least recently used fxcell */
		int ohash;

		/* If we haven't used all our memory, or if we are forced and have */
		/* no cell we can re-use, then allocate another fxcell */
		if (rc->s->rev.sz < rc->s->rev.max_sz
		 || (force && rc->nunlocked == 0)) {
			cp = increase_revcache(rc);
			hash = HASH(rc,ix);			/* Re-compute hash in case hash size changed */
//printf("~1 using new cell\n");
		} else {
//printf("~1 memory limit has been reached, using old cell\n");

			for (;;) {
				/* Use the least recently used unlocked fxcell */
				for (cp = rc->mrubot; cp != NULL && cp->refcount > 0; cp = cp->mruup)
					;
	
				/* Run out of unlocked cells */
				if (cp == NULL) {
//printf("~1 none available\n");
					return NULL;
				}
	
				/* If it has been used before, free up the simplexes */
				free_cell_contents(cp);

				/* Remove from current hash index (if it is in it) */
				ohash = HASH(rc,cp->ix);			/* Old hash */
				if (rc->hashtop[ohash] == cp) {
					rc->hashtop[ohash] = cp->hlink;
				} else {
					fxcell *c;
					for (c = rc->hashtop[ohash]; c != NULL && c->hlink != cp; c = c->hlink)
						;
					if (c != NULL)
						c->hlink = cp->hlink;
				}

				/* If we're now under the memory limit, use this fxcell */
				if (rc->s->rev.sz < rc->s->rev.max_sz) {
					break;
				}

//printf("~1 freeing a cell\n");
				/* Free up this cell and look for another one */
				/* Remove it from LRU list */
				if (rc->mrutop == cp)
					rc->mrutop = cp->mrudown;
				if (rc->mrubot == cp)
					rc->mrubot = cp->mruup;
				if (cp->mruup != NULL)
					cp->mruup->mrudown = cp->mrudown;
				if (cp->mrudown != NULL)
					cp->mrudown->mruup = cp->mruup;
				cp->mruup = cp->mrudown = NULL;
				free(cp);
				DECSZ(rc->s, sizeof(fxcell));
				rc->nacells--;
				rc->nunlocked--;
			}
		}

#ifdef STATS
		rc->s->rev.st[rc->s->rev.sb->op].cmiss++;
#endif /* STATS */

		/* Add this cell to hash index */
		cp->hlink = rc->hashtop[hash];
		rc->hashtop[hash] = cp;	/* Add to hash table and list */

		cp->ix = ix;
		cp->flags = 0;			/* Contents needs re-initializing */
//printf("~1 returning fresh cell\n");
	}
	
	/* Move slected cell to the top of the mru list */
	if (cp->mruup != NULL) {		/* This one wasn't already at top */
		cp->mruup->mrudown = cp->mrudown;
		if (cp->mrudown == NULL)	/* This was bottom */
			rc->mrubot = cp->mruup;	/* New bottom */
		else
			cp->mrudown->mruup = cp->mruup;
		/* Put this one at the top */
		rc->mrutop->mruup = cp;
		cp->mrudown = rc->mrutop;
		rc->mrutop = cp;
		cp->mruup = NULL;
	}
	if (cp->refcount == 0) {
		rc->nunlocked--;
	}

	cp->refcount++;

	return cp;
}

/* Tell the cache that we aren't using this cell anymore, */
/* but to keep it in case it is needed again. */
static void uncache_fxcell(
revcache *rc,		/* Reverse cache structure */
fxcell *cp
) {
	if (cp->refcount > 0) {
		cp->refcount--;
		if (cp->refcount == 0) {
			rc->nunlocked++;
		}
	} else
		warning("rspl cell cache assert: refcount overdecremented!");
}

/* ====================================================== */
/* Reverse rspl setup functions                           */

static void del_bxcell(rspl *s, bxcell *bx);
static void free_sharelist(rspl *s);
static void free_indexlist(rspl *s, int **rp);
static void free_surfhash(rspl *s, int del);
static void free_surflist(rspl *s);

/* Called by rspl initialisation */
/* Note that fxcell lookup tables are not */
/* allocated & created until the first call */
/* to a reverse interpolation function. */
void
init_rev(rspl *s) {

	/* First section */
	s->rev.inited = 0;
	s->rev.res = 0;
	s->rev.no = 0;
	s->rev.rev = NULL;

	/* Second section */
	s->rev.rev_valid = 0;
	s->rev.nnrev = NULL;

	/* Third section */
	s->rev.cache = NULL;

	/* Fourth section */
	s->rev.sb = NULL;

	/* Methods */
	s->rev_set_limit   = rev_set_limit_rspl;
	s->rev_get_limit   = rev_get_limit_rspl;
	s->rev_set_lchw    = rev_set_lchw;
	s->rev_interp      = rev_interp_rspl;
	s->rev_locus       = rev_locus_rspl;
	s->rev_locus_segs  = rev_locus_segs_rspl;
}

/* Free up all the reverse interpolation info */
void free_rev(
rspl *s		/* Pointer to rspl grid */
) {
	int e, di = s->di;
	int **rpp, *rp;
		
#ifdef STATS
	{
		int i, totcalls = 0;
		for (i = 0; i < 5; i++) {
			totcalls += s->rev.st[i].searchcalls;
		}

		printf("\n===============================\n");
		printf("di = %d, do = %d\n",s->di, s->fdi);
		for (i = 0; i < 5; i++) {
			int calls = s->rev.st[i].searchcalls;
			if (calls == 0) 
				continue;
			printf("\n- - - - - - - - - - - - - - - -\n");
			printf("Operation %s\n",opnames[i]);
			printf("Search calls = %d = %f%%\n",s->rev.st[i].searchcalls,
			100.0 * s->rev.st[i].searchcalls/totcalls);
			printf("Cells searched/call = %f\n",s->rev.st[i].csearched/(double)calls);
			printf("Simplexes searched/call = %f\n",s->rev.st[i].ssearched/(double)calls);
			printf("Simplexes inited level 1/call = %f\n",s->rev.st[i].sinited/(double)calls);
			printf("Simplexes inited level 2 (LU)/call = %f\n",s->rev.st[i].sinited2a/(double)calls);
			printf("Simplexes inited level 2 (SVD)/call = %f\n",s->rev.st[i].sinited2b/(double)calls);
			printf("Simplexes invalidated level 4/call = %f\n",s->rev.st[i].sinited4i/(double)calls);
			printf("Simplexes inited level 4/call = %f\n",s->rev.st[i].sinited4/(double)calls);
			printf("Simplexes invalidated level 5/call = %f\n",s->rev.st[i].sinited5i/(double)calls);
			printf("Simplexes inited level 5 (LU)/call = %f\n",s->rev.st[i].sinited5a/(double)calls);
			printf("Simplexes inited level 5 (SVD)/call = %f\n",s->rev.st[i].sinited5b/(double)calls);
			if ((s->rev.st[i].chits + s->rev.st[i].cmiss) == 0)
				printf("No cache calls\n");
			else
				printf("Cell hit rate = %f%%\n",
					100.0 * s->rev.st[i].chits/(double)(s->rev.st[i].chits + s->rev.st[i].cmiss));
		}
		printf("\n===============================\n");
	}
#endif /* STATS */

	/* Free up Fourth section */
	if (s->rev.sb != NULL) {
		free_search(s->rev.sb);
		s->rev.sb = NULL;
	}
	/* Free up the Third section */
	if (s->rev.cache != NULL) {
		free_revcache(s->rev.cache);	/* Reverse cell cache */
		s->rev.cache = NULL;
	}

	/* Free up the Second section */
	if (s->rev.nnrev != NULL) {

		/* Free up nn list sharelist records - this will free and set */
		/* any shared lists to NULL */
		free_sharelist(s);

		/* Free any remaining arrays at grid points */
		for (rpp = s->rev.nnrev; rpp < (s->rev.nnrev + s->rev.no); rpp++) {
			if (*rpp != NULL)
				free_indexlist(s, rpp);
		}
		free(s->rev.nnrev);
		DECSZ(s, s->rev.no * sizeof(int *));
		s->rev.nnrev = NULL;
	}

	if (di > 1 && s->rev.rev_valid) {
		rev_struct *rsi, **rsp;
		size_t ram_portion = g_avail_ram;

		/* Remove it from the linked list */
		for (rsp = &g_rev_instances; *rsp != NULL; rsp = &((*rsp)->next)) {
			if (*rsp == &s->rev) {
				*rsp = (*rsp)->next;
				break;
			}
		}

		/* Aportion the memory */
		g_no_rev_cache_instances--;

		if (g_no_rev_cache_instances > 0) {
			ram_portion /= g_no_rev_cache_instances; 
			for (rsi = g_rev_instances; rsi != NULL; rsi = rsi->next)
				rsi->max_sz = ram_portion;
			if (s->verbose)
				fprintf(stdout, "%cThere %s %d rev cache instance%s with %lu Mbytes limit\n",
				                cr_char,
								g_no_rev_cache_instances > 1 ? "are" : "is",
			                    g_no_rev_cache_instances,
								g_no_rev_cache_instances > 1 ? "s" : "",
			                    (unsigned long)ram_portion/1000000);
		}
	}

	s->rev.rev_valid = 0;

	if (s->rev.rev != NULL) {
		/* Free arrays at grid points */
		for (rpp = s->rev.rev; rpp < (s->rev.rev + s->rev.no); rpp++) {
			if (*rpp != NULL)
				free_indexlist(s, rpp);
		}
		free(s->rev.rev);
		DECSZ(s, s->rev.no * sizeof(int *));
		s->rev.rev = NULL;
	}

	/* If first section has been initialised */
	if (s->rev.inited != 0)	 {

		/* Sub-simplex information */
		for (e = 0; e <= di; e++) {
			rspl_free_ssimplex_info(s, &s->rev.sspxi[e]);
		}
		s->rev.res = 0;
		s->rev.no = 0;
		s->rev.inited = 0;
	}

	/* Free up surface linked list and the bxcells in it. */
	free_surflist(s);

	/* Free up surface bxcell hash index */
	free_surfhash(s, 0);

	DBG(("rev allocation left after free = %d bytes\n",s->rev.sz));

#ifdef CHECK_NNLU
	print_nnck(s);
#endif /* CHECK_NNLU */
}


/* ========================================================== */
/* reverse lookup acceleration structure initialisation code. */

/* The reverse lookup relies on a search of the fwd interpolation tables.
   To eliminate out of gamut points quickly, to provide a starting point for
   the search, and to guarantee that all possible reverse solutions are discovered,
   a spatial indexing structure is used to provide a list of starting candidate
   forward cell indexes for a given output value. (rev.rev[])
   The reverse structure contains two fdi dimensional bwd cell grids, each element of the
   cell grid holding the indexes of the forward interpolation grid.
   The rev[] grid holds fwd cell indexes which intersect that bwd cell's range of
   output values. A rev[] cell will be empty if there is no potential exact solution.
   The nnrev[] grid holds fwd cell indexes of those cells that may be the lch weighted
   closest to that bwd cell.
   The rev.nnrev[] array is almost a complement of the rev.rev[] array,
   with the exception of any overlap near the gamut surface.
   Since many of the nnrev[] bwd cells map to nearly the same surface region, many
   of the fwd cell lists are shared.

   When s->rev.fastsetup is set, then the rev.nnrev[] grid is left empty, and
   any call for nn lookup is satisfied by filling the requisite rev.nnrev[] on-demand,
   by an exaustive search of the surface bwd cells (rev.surflist)

   Note that unlike the forward grid which is composed of verticies,
   these rev lists are composed of fwd cells.

   The nnrev[] setup code identifies possible surface bwd revp[] cells
   by them being face neighbors of empty (out of gamut) bwd cells.
   It then converts the vertexes of the fwd cell list into a vertex list,
   and "thins" the list by deleting any vertex that is shaded by a triangle
   that other vertexes are part of. This is done on a backward cell basis,
   but includes vertexes of other possibly shadowed backward cells. 

   If ink limiting is being used, then over ink limit partners to
   the vertexes are added in, and then the list of vertexes is
   converted back into fwd cells in a way that ensures 2 dimensional
   connectivity of the cells, while minimizing the number of
   extra (non surface) vertexes implied by the fwd cells.

 */

/*
	The gamut hull fwcell finding code is not robust - it assumes visiblity
	of the surface from some center point(s).

	Perfect gamut hull finding approach would be something like this:
	(using vertex and triangle caching structures.)

	Add all triangles on device gamut surface with at least
	one vertex within ink limit.

	Add all triangles that are part of a full di simplex
	with at least one vertex within ink limit (and not on device gamut),
	where all the other verticies of the simplex are on one
	side of the triangle (non-monotonic surfaces).

	Add all triangles on the ink limit plane.
	(Will be 1 or more triangles per simplex that has
	1..di verticies that are over the ink limit.)

	Check all triangles for instersection with each other.
	Convert any such intersections into smaller, non-intersecting
	triangles that share verticies along intersection line.

	Delete triangles that have dangling edges (i.e. triangles that
	have edges with odd number of associated triangles).
	This is to eliminate "dangling" triangles. Should only be left
	with "bubles" in surface after this ?

	Bubbles join at edges where more than 2 triangles co-incide.
	Can internal bubles be "un-stitched" if we can decide which
	triangles are part of a bubble ????
	i.e. use even/odd inside rule for points between
	the triangles at the edge.
	
	Delete all vertexes and associated triangles that are
	inside the surface.
	Will odd/even test work ? - i.e. from vertex of triangle,
	is on surface if intersections in one direction are even, and
	other direction are odd. 

	Or "point within odd number of tetrahedrons formed with point on surface" ?	
	- seems to be the same as the odd/even rule. Can't detect connectivity.
		
	Or do this using a winding number algorithm
	with signed crossings optimization ?
	<Point in Polyhedron Testing Using Spherical Polygons, Graphics Gems V pp42>
    But do we have to order triangles in a consistent direction ?
  	How to do this when more than 2 triangles meet at an adge ???
	i.e. catch-22 - need to know which are inside triangles to
	set edge direction, but need edge direction to detect inside-outside.

*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if defined(REVTABLESTATS) || defined(DEBUG)
static int bxcount = 0;
static int maxbxcount = 0;
#endif

static void add2indexlist(rspl *s, int **rpp, int ix, int shrec);
static void comp_shadow_group(rspl *s, double *gcent, double *rgc, double *pcc,
	double *pdw, double *gc, double (*v)[MXRO], int nverta);

/* Allocate a new bx cell. */
/* (Doesn't add to hash or list) */
static bxcell *new_bxcell(
	rspl *s,
	int ix,				/* rev[] index of cell being created */
	int *gc,			/* Coord of rev[] cell being created */
	bxcell *ss,			/* search starting bxcell to commence with, this cell if NULL */
	double sdist,		/* Est. distance from this cell to six */
	char *vflag			/* If non-NULL, create a super-cell if far from seed */
) {
	int f, fdi = s->fdi;
	int i;
	bxcell *bx = NULL;
	DCOUNT(cc, MXRO, fdi, 0, 0, 2);		/* Vertex counter */

//printf("~1 creating new bxcell with index %d\n",ix);
	if ((bx = (bxcell *) rev_calloc(s, 1, sizeof(bxcell))) == NULL)
		error("rspl malloc failed - rev bxcell structs");
	INCSZ(s, sizeof(bxcell));

	bx->ix = ix;
	bx->tix = -1;
	for (f = 0; f < fdi; f++)
		bx->gc[f] = gc[f];
	bx->ss = (ss == NULL) ? bx : ss;
	bx->sdist = sdist;

//printf("~1 new_bxcell ix %d, co %s, base %s\n",ix,debPiv(s->fdi, bx->gc),debPdv(s->fdi, vp[0]));

	/* super-cell code (to speed filling) */
	if (vflag != NULL && (vflag[ix] & 2) == 0 && ss != NULL) {
		double codist = 0.0;

		/* Compute distance of seed from this cell */
		for (codist = 0.0, f = 0; f < fdi; f++) {
			int tt = bx->gc[f] - ss->gc[f];
			codist += tt * tt;
		}
		codist = sqrt(codist);

//printf("~1 codist %f, codist/s->rev.res = %f\n",codist,codist/s->rev.res);
		/* Create a super-cell if we are far enough from the seed. */
		/* (this determines what portion of filling uses super-cells) */
//		if (codist >= 1.0 && (codist/s->rev.res) > 0.05)
		if (codist >= 2.0)
		{
			int co[MXRO]; 
			DCOUNT(ss, MXRO, s->fdi, -1, -1, 2);
			double (*vp)[MXRO];
			double **vpp;
			int nverts;
//printf("~1 creating super-cell for bx %d\n",ix);

			/* Maximum number of verticies for all surrounders */
			for (nverts = (1 << fdi), f = 0; f < fdi; f++)
				nverts *= 3;

			if ((vp = (double(*)[MXRO]) rev_calloc(s, nverts, sizeof(double) * MXRO)) == NULL)
				error("rspl malloc failed - rev bxcell vertex list");
			INCSZ(s, nverts * sizeof(double) * MXRO);
	
			if ((vpp = (double **) rev_calloc(s, nverts, sizeof(double *))) == NULL)
				error("rspl malloc failed - rev bxcell vertex list");
			INCSZ(s, nverts * sizeof(double *));
	
			/* Search around this cell for other cells to be filled */
			i = 0;
			DC_INIT(ss);
			while (!DC_DONE(ss)) {
				int nix = ix;
				for (f = 0; f < fdi; f++) {
					nix += ss[f] * s->rev.coi[f];
					co[f] = bx->gc[f] + ss[f];
					if (co[f] < 0 || co[f] >= s->rev.res)
						break;
				}

				/* If within boundary and un-filled non-surface bxcell */
				if (f >= fdi && (vflag[nix] & 0xf) == 0) {
					add2indexlist(s, &bx->scell, nix, 0);
					vflag[nix] = (vflag[nix] & ~0xf) | 1;	/* Assume it's now on the seed list */

					/* Create vertex locations for this bxcell */
					DC_INIT(cc);
					while (!DC_DONE(cc)) {
						for (f = 0; f < fdi; f++)
							vp[i][f] = (co[f] + cc[f]) * s->rev.gw[f] + s->rev.gl[f];  
						vpp[i] = vp[i];
						DC_INC(cc);
						i++;
					}
				}
				DC_INC(ss);
			}

			/* Init the group boundary data */
			nn_grpinit(s, &bx->g, vpp, i, NULL);
	
			/* Compute the default shadowing test width and distance */
			/* (Not actually used for super-cell ?) */
			comp_shadow_group(s, s->rev.ocent, NULL, &bx->cc, &bx->dw, bx->g.bcent, vp, i);   

			free(vpp);
			DECSZ(s, nverts * sizeof(double *));

			free(vp);
			DECSZ(s, nverts * sizeof(double) * MXRO);
//printf(" - %d sub-cells\n",bx->scell[1]-3);
		}
	}

	if (bx->scell == NULL) {
		double vp[POW2MXRO][MXRO];
		double *vpp[POW2MXRO];

		/* Create vertex locations for this bxcell */
		i = 0;
		DC_INIT(cc);
		while (!DC_DONE(cc)) {
			for (f = 0; f < fdi; f++)
				vp[i][f] = (gc[f] + cc[f]) * s->rev.gw[f] + s->rev.gl[f];  
			vpp[i] = vp[i];
			DC_INC(cc);
			i++;
		}

		/* Init the group boundary data */
		nn_grpinit(s, &bx->g, vpp, i, NULL);

		/* Compute the default shadowing test width and distance */
		comp_shadow_group(s, s->rev.ocent, NULL, &bx->cc, &bx->dw, bx->g.bcent, vp, 1 << fdi);   
	}

//printf("~1 grp bcent %s, brad %f\n",debPdv(s->fdi, bx->g.bcent), bx->g.brad);
#if defined(REVTABLESTATS) || defined(DEBUG)
	bxcount++;
	if (bxcount > maxbxcount)
		maxbxcount = bxcount;
//printf("~1 now %d bxcells\n",bxcount);
#endif
	return bx;
}

/* Free a bxcell (up to caller to free bx->sl, remove from cache etc.) */
/* We free the super-cell info. */
static void del_bxcell(rspl *s, bxcell *bx) {
	if (bx->scell != NULL)	/* If this is a supercell */
		free_indexlist(s, &bx->scell);
	if (bx->dl != NULL)		/* We have a deleted fwd vertex list */
		free_indexlist(s, &bx->dl);
	free(bx);
	DECSZ(s, sizeof(bxcell));
#if defined(REVTABLESTATS) || defined(DEBUG)
	bxcount--;
//printf("~1 now %d bxcells\n",--bxcount);
#endif
}

/* Allocate the surflist hash index */
static void create_surfhash(rspl *s) { 

	s->rev.surf_hash_size = primes[2];		/* 1489 */
	if ((s->rev.surfhash = (bxcell **) rev_calloc(s, s->rev.surf_hash_size, sizeof(bxcell *))) == NULL)
		error("rspl malloc failed - reverse bxcell surface cache index");
	INCSZ(s, s->rev.surf_hash_size * sizeof(bxcell *));
}

/* Add a bxcell to the surface hash list */
static void add_bxcell_hash(rspl *s, bxcell *bx) {
	unsigned int hash = 0;

	hash = bx->ix % s->rev.surf_hash_size;
	bx->hlink = s->rev.surfhash[hash];
	s->rev.surfhash[hash] = bx;
}

/* Remove a bxcell from the surface hash list. */
/* Doesn't delete the bxcell though. */
static void rem_bxcell_hash(rspl *s, int ix) {
	unsigned int hash = 0;
	bxcell *bx = NULL, **pbx;

	hash = ix % s->rev.surf_hash_size;

	for (pbx = &s->rev.surfhash[hash], bx = *pbx; bx != NULL; pbx = &bx->hlink, bx = *pbx) {
		if (bx->ix == ix) {
			*pbx = bx->hlink;
			return;
		}
	}
}

/* Fetch a surface bxcell from the surface hash list, given its index, or */ 
/* Return NULL if none */
static bxcell *get_surface_bxcell(rspl *s, int ix) {
	unsigned int hash = 0;
	bxcell *bx = NULL;

	hash = ix % s->rev.surf_hash_size;

	for (bx = s->rev.surfhash[hash]; bx != NULL; bx = bx->hlink) {
		if (bx->ix == ix)
			return bx;
	}
	return NULL;
}

/* Free up surface linked list and delete the bxcells. */
/* (If we use this, don't use free_surfhash with del set.) */
static void free_surflist(rspl *s) {

	while (s->rev.surflist != NULL) {
		bxcell *this = s->rev.surflist;
		s->rev.surflist = s->rev.surflist->slist;
		if (this->sl != NULL)
			free_indexlist(s, &this->sl);
		del_bxcell(s, this);
	}
}


/* If del is set, free up all the bxcell cells in the hash index, */
/* then free the surfhash itself. */
/* (Use instead of surflist to manage allocation, */
/*  or to clean up hashlist after surflist has been freed.) */
static void free_surfhash(rspl *s, int del) {

	if (s->rev.surfhash != NULL) {
		if (del) {
			int i;
			for (i = 0; i < s->rev.surf_hash_size; i++) {
				bxcell *bx, *nbx;
				for (bx = s->rev.surfhash[i]; bx != NULL; bx = nbx) {
					nbx = bx->hlink;
					if (bx->sl != NULL)
						free_indexlist(s, &bx->sl);
					del_bxcell(s, bx);
				}
			}
		}
		free(s->rev.surfhash);
		DECSZ(s, s->rev.surf_hash_size * sizeof(bxcell *));
		s->rev.surfhash = NULL;
		s->rev.surf_hash_size = 0;
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Structure to cache prime vertex information when filtering surface */
/* cell vertex lists. */

/* vertex status */
typedef enum {
	vtx_norm = 0,		/* Normal vertex in primary bxcell - initial value */
	vtx_sha  = 1,		/* Vertex has been shadowed */
	vtx_del  = 2,		/* Vertex has been deleted because it's shadowed */
	vtx_oil  = 3		/* Vertex is over ink limit */
} vstat;

struct _vtxrec {
	int ix;						/* fwd index of vertex */
	int cix;					/* Fwd cell vertex is in index */
	double vv[MXRO];			/* Output value of vertex */
	double vl[MXRO];			/* Log compressed output value of vertex */
	double dist;				/* Distance from center point squared */
	int tcount;					/* Touch count for converting to fwd cells */
	int acount;					/* Actual count for converting to fwd cells */

	vstat status;					
	int tix;					/* Target vertex when being created */

	struct _vtxrec *hlink;		/* Linked list of vtxrecs with same ix hash */
	int rix;					/* nnrev[] index vertex falls into */
	int ival[MXRO];				/* nnrev[] coordinate rix */ 

	char prim;					/* nz when primary vertex of bx (not shadow bx) */
	char cross;					/* nz when part of suspected crossed triangle */
	char pres;					/* nz when preserved shadowed vertex from crossed triangle */

	char tflag;					/* nz when on tlist */
	struct _vtxrec *tlist;		/* Linked list of vertexes for nnrev[] cell/freelist */

#ifdef REVVRML
	int addvtx;					/* Vertex that caused a bxcell to be added */
	int vrmlix;					/* Index for plotting */
#endif

}; typedef struct _vtxrec vtxrec;

struct _vtxcache {
	vtxrec *vtxlist;		/* vertex list for soring/itterating selected nnrev cell. */
	int nilist;				/* Number of vertexes in the list */

	int hash_size;			/* Current size of vtxrec hash list */
	vtxrec **hash;			/* hash index list */

	vtxrec *freelist;		/* Unused vertex structures (to avoid memory allocs) */
}; typedef struct _vtxcache vtxcache;


/* Create the (empty) vertex list & hash */
static void create_vtxrec_list(rspl *s, vtxcache *vc) {
	vc->hash_size = primes[3];		/* 3373 */
	if ((vc->hash = (vtxrec **) rev_calloc(s, vc->hash_size, sizeof(vtxrec *))) == NULL)
		error("rspl malloc failed - vtxrec cache index");
	INCSZ(s, vc->hash_size * sizeof(vtxrec *));
	vc->vtxlist = NULL;
	vc->nilist = 0;
	vc->freelist = NULL;
}

/* Clear the vertex hash and list */
static void clear_vtxrec_lists(rspl *s, vtxcache *vc) {
	vtxrec *vp, *nvp;
	int i;

	/* Transfer all records in hash to freelist, */
	/* and clear hash. */
	for (i = 0; i < vc->hash_size; i++) {
		for (vp = vc->hash[i]; vp != NULL; vp = nvp) {
			nvp = vp->hlink;
			vp->tlist = vc->freelist;
			vc->freelist = vp;
		}
		vc->hash[i] = NULL;
	}

	vc->vtxlist = NULL;
	vc->nilist = 0;
}

/* Free the vertex list & hash */
static void free_vtxrec_list(rspl *s, vtxcache *vc) {
	clear_vtxrec_lists(s, vc);

	while (vc->freelist != NULL) {
		vtxrec *this = vc->freelist;
		vc->freelist = vc->freelist->tlist;
		free(this);
		DECSZ(s, sizeof(vtxrec));
	}
	free(vc->hash);
	DECSZ(s, vc->hash_size * sizeof(vtxrec *));
	vc->hash = NULL;
	vc->hash_size = 0;
}

/* Add a vtxrec to the vertex hash list */
static void add_vtxrec_hash(vtxcache *vc, vtxrec *vx) {
	unsigned int hash = 0;

	hash = vx->ix % vc->hash_size;
	vx->hlink = vc->hash[hash];
	vc->hash[hash] = vx;
}

/* Delete a vtxrec from the vertex hash list */
/* (Assume it's not part of vtxlist!) */
static void del_vtxrec_hash(vtxcache *vc, int ix) {
	unsigned int hash = 0;
	vtxrec *vx = NULL, **pvx;

	hash = ix % vc->hash_size;

	for (pvx = &vc->hash[hash], vx = *pvx; vx != NULL; pvx = &vx->hlink, vx = *pvx) {
		if (vx->ix == ix) {
			*pvx = vx->hlink;
			vx->tlist = vc->freelist;
			vc->freelist = vx;
			vx->hlink = NULL;
			return;
		}
	}
}

/* Fetch a surface vtxrec from the hash list, given its index */ 
/* Return NULL if none */
static vtxrec *get_vtxrec(vtxcache *vc, int ix) {
	unsigned int hash = 0;
	vtxrec *vx = NULL;

	hash = ix % vc->hash_size;

	for (vx = vc->hash[hash]; vx != NULL; vx = vx->hlink) { 
		if (vx->ix == ix)
			return vx;
	}
	return NULL;
}

/* Log compress an output value wrt to center point */
static void logcomp(
	rspl *s,
	double *out,
	double *in,
	double *cent
) {
	int f, fdi = s->fdi;
	double len;

	if (s->rev.surflin_en) {
#ifdef NEVER
		/* (This doesn't seem to improve things) */
		/* Calculate vector length */
		for (len = 0.0, f = 0; f < fdi; f++) {
			double tt= in[f] - cent[f]; 
			len += tt * tt;
		}
		len = sqrt(len);

		/* change length to log length */
		if (len > DBL_EPSILON) {
 
			len = 20.0 * pow(len, 0.25)/len;	/* Ratio */ 

			for (f = 0; f < fdi; f++) {
				double tt = in[f] - cent[f]; 
				out[f] = len * tt + cent[f];
			}
		}
#else
		if (s->rev.surflin != NULL) {
			co p;

			for (f = 0; f < fdi; f++)
				p.p[f] = in[f];
			s->rev.surflin->interp(s->rev.surflin, &p);
			for (f = 0; f < fdi; f++)
				out[f] = p.v[f] - s->rev.linoff[f];
		} else {
			for (f = 0; f < fdi; f++)
				out[f] = in[f];
		}
#endif
	} else {
		for (f = 0; f < fdi; f++)
			out[f] = in[f];
	}
}

/* Create a new vtxrec or return the current one. */
/* Allocates it, adds it to cache. */
/* DOESN"T add it to vtxlist. */
/* If new, sets status = vtx_norm */
static vtxrec *new_vtxrec(
	rspl *s,
	vtxcache *vc,
	int ix				/* fwd index of vertex */
) {
	int e, di = s->di;
	int f, fdi = s->fdi;
	vtxrec *vx = NULL;
	float *gp;
	int rix;
	int rgres_1 = s->rev.res -1;	/* rgres -1 == maximum base coord value */

	/* See if we've already got this vertex */
	if ((vx = get_vtxrec(vc, ix)) != NULL)
		return vx;

	/* Fetch or allocate a new structure */
	if (vc->freelist != NULL) {		/* Grab one from free list */
		vx = vc->freelist;
		vc->freelist = vx->tlist;
		memset((void *)vx, 0, sizeof(vtxrec));

	} else {
		if ((vx = (vtxrec *) rev_calloc(s, 1, sizeof(vtxrec))) == NULL)
			error("rspl malloc failed - rev vtxrec structs");
		INCSZ(s, sizeof(vtxrec));
	}

	/* Our fwd index */
	vx->ix = ix;

	/* Add it to the hash */
	add_vtxrec_hash(vc, vx);

	/* Fwd vertex array address */
	gp = s->g.a + ix * s->g.pss;

	/* Set cell index so that cell verticies don't exceed grid boundary */
	vx->cix = ix;
	for (e = 0; e < di; e++) {
		if (G_FL(gp, e) == 0)		/* At the top edge */
			vx->cix -= s->g.ci[e];	/* Move cell base down a row */
	}

	/* Get the output value */
	for (f = 0; f < fdi; f++)
		vx->vv[f] = gp[f];

	/* Set vl[] */
	logcomp(s, vx->vl, vx->vv, s->rev.ocent);

	/* Compute distance to overall center point squared */
	vx->dist = 0.0;
	for (f = 0; f < fdi; f++) {
		double tt = vx->vl[f] - s->rev.ocent[f]; 
		vx->dist += tt * tt;
	}

	/* Figure the actual nncell it lands in */
	for (rix = f = 0; f < fdi; f++) {
		double t;
		int mi;
		double gw = s->rev.gw[f];
		double gl = s->rev.gl[f];
		t = (vx->vv[f] - gl)/gw;
		mi = (int)floor(t);			/* Grid coordinate */
		if (mi < 0)					/* Limit to valid cube base index range */
			mi = 0;
		else if (mi > rgres_1)
			mi = rgres_1;
		vx->ival[f] = mi;
		rix += mi * s->rev.coi[f];
	}
	vx->rix = rix;

	return vx;
}

/* Add a vertex to the list. */
/* Don't add if already on list (if tflag set), or if shadowed) */
/* set prim flag to value */
static void add_vtxrec_list(vtxcache *vc, vtxrec *vx, int prim) {

	vx->prim = (char)prim;		/* Always set prim flag */

	if (vx->tflag || vx->status != vtx_norm)
		return;

	vx->tlist = vc->vtxlist;
	vc->vtxlist = vx;
	vx->tflag = 1;
	vc->nilist++;

}

int dumpvtxsort = 0;

/* Sort the vertex linked list by dist. */
/* Also reset the tflag */ 
static void sort_vtxrec_list(rspl *s, vtxcache *vc) {
	int i;
	vtxrec **sort, *vx;

	/* Create temporary array of pointers to vtxrec's in list */
	if ((sort = (vtxrec **) rev_calloc(s, vc->nilist, sizeof(vtxrec *))) == NULL)
		error("rspl malloc failed - rev vtxrec sort array");
	INCSZ(s, vc->nilist * sizeof(vtxrec *));

	for (i = 0, vx = vc->vtxlist; vx != NULL; vx = vx->tlist, i++)
		sort[i] = vx;

	/* Sort the list into ascending distance from center */
#define 	HEAP_COMPARE(A,B) (A->dist < B->dist)
	HEAPSORT(vtxrec *, sort, vc->nilist)
#undef 		HEAP_COMPARE

	/* Re-create the linked list in descending order */
	vc->vtxlist = NULL;
	for (i = 0; i < vc->nilist; i++) {
		vx = sort[i];
		vx->tlist = vc->vtxlist;
		vc->vtxlist = vx;
		vx->tflag = 0;
	}

	free(sort);
	DECSZ(s, vc->nilist * sizeof(vtxrec *));

#ifndef NEVER
	if (dumpvtxsort) {
		printf("sorted vertex list:\n");
		for (i = 0, vx = vc->vtxlist; vx != NULL; vx = vx->tlist, i++)
			printf("%d: ix %d, dist %f\n",i,vx->ix, sqrt(vx->dist));
	}
#endif
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Structure to cache surface triangle vertexes, to avoid repeated */
/* shadowing test */

struct _trirec{
	int ix[3];					/* vertex indexes of triangle in simplex order */
	struct _trirec *hlink;		/* Linked list of triangles in hash/freelist */
}; typedef struct _trirec trirec;

typedef struct {
	int hash_size;			/* Current size of trirec hash list */
	trirec **hash;			/* hash index list */
	trirec *freelist;		/* Unused trirec structures (to avoid memory allocs) */
} tricache;

/* Create the tricache list & hash. */
/* Set sm flag if we only want a small cache size */
static void create_trirec(rspl *s, tricache *tc, int sm) {
	if (sm)
		tc->hash_size = primes[1];		/* 853 */
	else
		tc->hash_size = primes[5];		/* 12919 */
	if ((tc->hash = (trirec **) rev_calloc(s, tc->hash_size, sizeof(trirec *))) == NULL)
		error("rspl malloc failed - trirec cache index");
	INCSZ(s, tc->hash_size * sizeof(trirec *));
	tc->freelist = NULL;
}

/* Clear the trirec list & hash */
static void clear_trirec(rspl *s, tricache *tc) {
	int i;
	trirec *tp, *ntp;

	/* Transfer all records in hash to freelist, */
	/* and clear hash. */
	for (i = 0; i < tc->hash_size; i++) {
		for (tp = tc->hash[i]; tp != NULL; tp = ntp) {
			ntp = tp->hlink;
			tp->hlink = tc->freelist;
			tc->freelist = tp;
		}
		tc->hash[i] = NULL;
	}
}

/* Free the triangle list & hash */
static void free_trirec(rspl *s, tricache *tc) {
	clear_trirec(s, tc);

	while (tc->freelist != NULL) {
		trirec *this = tc->freelist;
		tc->freelist = tc->freelist->hlink;
		free(this);
		DECSZ(s, sizeof(trirec));
	}
	free(tc->hash);
	DECSZ(s, tc->hash_size * sizeof(trirec *));
	tc->hash = NULL;
	tc->hash_size = 0;
}

/* Check if a triangle is in the cache. */
/* return nz if it is, and z if it isn't, and add it. */
static int check_trirec(rspl *s, tricache *tc, int *ix) {
	int i;
	unsigned int hash = 0;
	trirec *tp = NULL;

	hash = ix[0];
	hash = hash * 17 + ix[1];
	hash = hash * 17 + ix[2];
	hash %= tc->hash_size;

	for (tp = tc->hash[hash]; tp != NULL; tp = tp->hlink) {
		if (tp->ix[0] == ix[0]
		 && tp->ix[1] == ix[1]
		 && tp->ix[2] == ix[2]) {
//printf("check_trirec %d %d %d is in cache\n",ix[0], ix[1], ix[2]);
			return 1;
		}
	}
//printf("check_trirec %d %d %d NOT in cache\n",ix[0], ix[1], ix[2]);

	/* Allocate a new structure */
	if (tc->freelist != NULL) {		/* Grab one from free list */
		tp = tc->freelist;
		tc->freelist = tp->hlink;
		memset((void *)tp, 0, sizeof(trirec));

	} else {
		if ((tp = (trirec *) rev_calloc(s, 1, sizeof(trirec))) == NULL)
			error("rspl malloc failed - rev trirec structs");
		INCSZ(s, sizeof(trirec));
	}
	
	tp->ix[0] = ix[0];
	tp->ix[1] = ix[1];
	tp->ix[2] = ix[2];

	/* add it into the hash */
	tp->hlink = tc->hash[hash];
	tc->hash[hash] = tp;

	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Add another entry to an index/share list, taking care of any re-sizing */
/* Set shlist if this is a sharer record */
static void add2indexlist(rspl *s, int **rpp, int ix, int shrec) {
	int *rp = *rpp;

	if (rp == NULL) {
		if ((rp = (int *) rev_malloc(s, 6 * sizeof(int))) == NULL)
			error("rspl malloc failed - rev.grid list");
		INCSZ(s, 6 * sizeof(int));
		rp[0] = 6;		/* Allocation */
		rp[1] = 4;		/* Next free Cell */
		rp[2] = -1;		/* share list index - default none */
		rp[3] = ix;		/* Index added to list */
		rp[4] = -1;		/* End of list marker */
		*rpp = rp;		/* Update pointer */
	} else {
		int z = rp[1], ll = rp[0];
		if (z >= (ll-1)) {			/* Not enough space */
			if (!shrec && rp[2] != -1)
				error("Re-allocating shared fwd index list");
			INCSZ(s, ll * sizeof(int));
			ll *= 2;
			if ((rp = (int *) rev_realloc(s, rp, sizeof(int) * ll)) == NULL)
				error("rspl realloc failed - rev.grid list size %d",ll);
			rp[0] = ll;	/* Allocation */
			*rpp = rp;	/* Update pointer */
		}
		rp[z++] = ix;	/* Index added to list */
		rp[z] = -1;		/* End of list marker */
		rp[1] = z;		/* Next free Cell */
	}
}

/* Copy an index list (i.e. from nnrev[] to bxcell->sl) */
static void copy_indexlist(rspl *s, int **dp, int *sp) {
	if (sp == NULL)
		*dp = NULL;
	else {
		int i;
		if ((*dp = (int *) rev_malloc(s, sp[0] * sizeof(int))) == NULL)
			error("rspl malloc failed - rev.grid list");
		INCSZ(s, sp[0] * sizeof(int));
		for (i = 0; i <= sp[1]; i++)
			(*dp)[i] = sp[i];
		(*dp)[2] = -1;
	}
}

/* Free an index list, at set it to NULL */
static void free_indexlist(rspl *s, int **rp) {
	if (*rp != NULL) {
		DECSZ(s, (*rp)[0] * sizeof(int));
		free(*rp);
		*rp = NULL;
	}
}

/* Add a (fwd index list) sharer to share list. */
/* Record will be created if list[2] == -1, */
/* or incremented otherwise. */
/* sharerix is the index of the cell sharing the *list */
static void add2sharelist(rspl *s, int sharerix, int *list) {
	int *sharerec = NULL;

	/* Create a new record and add our (one) sharer to it */
	if (list[2] == -1) {
		if (s->rev.sharellen >= s->rev.sharelaloc) {
			/* Allocate another sharelist entry */
			INCSZ(s, (10 + s->rev.sharelaloc) * sizeof(int *));
			s->rev.sharelaloc = 10 + 2 * s->rev.sharelaloc;
			if ((s->rev.sharelist = (int **)rev_realloc(s, s->rev.sharelist,
				                           s->rev.sharelaloc * sizeof(int *))) == NULL)
				error("add2sharelist: realloc failed");
		}
		add2indexlist(s, &sharerec, sharerix, 1);
		s->rev.sharelist[s->rev.sharellen] = sharerec;
		list[2] = s->rev.sharellen;
		s->rev.sharellen++;

	/* Add the sharer to the existing sharer list */
	} else {
		if (list[2] >= s->rev.sharellen)
			error("add2sharelist got list with sharelist index out of range");
		sharerec = s->rev.sharelist[list[2]];
		add2indexlist(s, &sharerec, sharerix, 1);
		s->rev.sharelist[list[2]] = sharerec;
	}
}

/* Return the sharer list for the given (fwd cell) list */
/* Return NULL if not shared */
static int *getsharelist(rspl *s, int *list) {
	if (list[2] == -1)
		return NULL;
	if (list[2] >= s->rev.sharellen) {
		error("getsharelist got list with sharelist index out of range (%d > %d)",list[2],s->rev.sharellen);
	}
	return s->rev.sharelist[list[2]];
}

/* Free all the sharelist and the shared nnrev[] fwd cell lists as well */
static void free_sharelist(rspl *s) {
	if (s->rev.sharelist != NULL) {
		int i, j;
		for (i = 0; i < s->rev.sharellen; i++) {
			int *shrec = s->rev.sharelist[i];

			/* Free the shared fwd cell list */
			if (shrec[1] > 3) {
				int *clist = s->rev.nnrev[shrec[3]];
				DECSZ(s, clist[0] * sizeof (int));
				free(clist);
			}

			/* Make sure freeing of s->rev.nnrev[] doesn't free them twice */
			for (j = 3; shrec[j] != -1; j++) 
				s->rev.nnrev[shrec[j]] = NULL;

			DECSZ(s, s->rev.sharelist[i][0] * sizeof (int));
			free(s->rev.sharelist[i]);
		}
		DECSZ(s, s->rev.sharelaloc * sizeof(int *));
		free(s->rev.sharelist);
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* For shadow bxcell testing, compute the delta distance */
/* rev.ocent, and delta "width" between two vertex values. */ 

/* Compute shadow testing group values */ 
static void comp_shadow_group(
	rspl *s,
	double *gcent,			/* In gamut center point to compute from */			
	double *rgc,			/* Return group center if non-NULL */
	double *pcc,			/* Return distance of group center to gamut center */  
	double *pdw,			/* Return width of furthest point from group center */
	double *gc,				/* if not-NULL, the group center */
	double (*v)[MXRO],		/* Input verticies */
	int nvert				/* Number of verticies */
) {
	double _gc[MXRO];
	int i;
	int f, fdi = s->fdi;
	double cc;				/* gamut center to group center */
	double dw = -1.0;		/* Largest goup vertex width */

	/* if no group center given, compute one simply as an average */ 
	/* (Used for triangle) */
	if (gc == NULL) {
		gc = _gc;
		for (f = 0; f < fdi; f++)
			gc[f] = 0.0;
		for (i = 0; i < nvert; i++) {
			for (f = 0; f < fdi; f++) {
				gc[f] += v[i][f];
			}
		}
		for (f = 0; f < fdi; f++)
			gc[f] /= (double)nvert;
	}

	/* Return it if requested */
	if (rgc != NULL) {
		for (f = 0; f < fdi; f++)
			rgc[f] = gc[f];
	}

	/* Compute distance from gamut center to group center */ 
	for (cc = 0.0, f = 0; f < fdi; f++) {
		double tt = gcent[f] - gc[f];
		cc += tt * tt;
	}
	cc = sqrt(cc);

	if (pcc != NULL)
		*pcc = cc;

	/* Compute width for each vertex, and track maximum */
	for (i = 0; i < nvert; i++) {
		double vlen, scale;
		double sv[MXRO];		/* Vertex scaled to same distance as group center */
		double w;

		/* vertex length from gamut center */
		for (vlen= 0.0, f = 0; f < fdi; f++) {
			double tt = v[i][f] - gcent[f];
			vlen += tt * tt;
		}
		vlen = sqrt(vlen);
	
		if (vlen > 1e-6)
			scale = cc/vlen;
		else
			scale = 1.0;

		for (f = 0; f < fdi; f++)
			sv[f] = (scale * (v[i][f] - gcent[f])) + gcent[f];

		/* Distance from scaled vertex to group center */
		for (w = 0.0, f = 0; f < fdi; f++) {
			double tt = sv[f] - gc[f];
			w += tt * tt;
		}
		if (w > dw)
			dw = w;
		
	}
	dw = sqrt(dw);

	if (pdw != NULL)
		*pdw = dw;
}

/* Expand a bxcell's shadow testing group values based on it's vertex list */
static void extend_bxcell_shadow_group(
	rspl *s,
	vtxcache *vc,
	bxcell *bx
) {
	int *ip;
	int f, fdi = s->fdi;
	double dw;

	if (bx->sl == NULL)
		return;

	/* Current dw squared */
	dw = bx->dw * bx->dw;

	/* Compute width for each vertex, and track maximum */
	for (ip = bx->sl+3; *ip != -1; ip++) {
		vtxrec *vx;
		double vlen, scale;
		double sv[MXRO];		/* Vertex scaled to same distance as group center */
		double w;

		if ((vx = get_vtxrec(vc, *ip)) == NULL)
			continue;

		/* vertex length from gamut center */
		for (vlen= 0.0, f = 0; f < fdi; f++) {
			double tt = vx->vl[f] - s->rev.ocent[f];
			vlen += tt * tt;
		}
		vlen = sqrt(vlen);
	
		if (vlen > 1e-6)
			scale = bx->cc/vlen;
		else
			scale = 1.0;

		for (f = 0; f < fdi; f++)
			sv[f] = (scale * (vx->vl[f] - s->rev.ocent[f])) + s->rev.ocent[f];

		/* Distance from scaled vertex to group center */
		for (w = 0.0, f = 0; f < fdi; f++) {
			double tt = sv[f] - bx->g.bcent[f];
			w += tt * tt;
		}
		if (w > dw)
			dw = w;
		
	}
	if (dw > bx->dw)
		bx->dw = sqrt(dw);
}

/* Shadow group to group compare. Return nz if within range */
static int shadow_group_group(
	rspl *s,
	double *gcent,		/* Input gamut center point to compute from */			
	double *gc1,		/* Reference group center point */
	double cc1,			/* Reference point cc value */
	double dw1,			/* Reference point dw value */
	double *gc2,		/* Comparison group center point */
	double cc2,			/* Comparison point cc value */
	double dw2			/* Comparison point dw value */
) {
	int i;
	int f, fdi = s->fdi;
	double dot, scale;
	double sv[MXRO];		/* Comparison group center scaled to same distance as ref center */
	double w;

	/* Compute dot product of cc1 and cc2 */
	for (dot = 0.0, f = 0; f < fdi ; f++)
		dot += (gc1[f] - gcent[f]) * (gc2[f] - gcent[f]);

	/* If the groupls are not in the same direction, return false */
	if (dot < 0.0)
		return 0;
  
	if (cc2 > 1e-6)
		scale = cc1/cc2;
	else
		scale = 1.0;

	for (f = 0; f < fdi; f++)
		sv[f] = (scale * (gc2[f] - gcent[f])) + gcent[f];

	/* Distance from scaled group center to ref. group center */
	for (w = 0.0, f = 0; f < fdi; f++) {
		double tt = sv[f] - gc1[f];
		w += tt * tt;
	}
	w = sqrt(w);

	if (w <= (dw1 + (scale * dw2) + EPS))
		return 1;

	return 0;
}

/* Shadow group to vertex compare. Return nz if within range */
static int shadow_group_vertex(
	rspl *s,
	double *gcent,		/* Input gamut center point to compute from */			
	double *gc1,		/* Reference group center point */
	double cc1,			/* Reference point cc value */
	double dw1,			/* Reference point dw value */
	double *v			/* Comparison vertex location */
) {
	int i;
	int f, fdi = s->fdi;
	double vlen, dot, scale;
	double sv[MXRO];		/* Vertex scaled to same distance as group center */
	double w;

	/* Compute dot product of cc1 and cc2 */
	for (dot = 0.0, f = 0; f < fdi ; f++)

	/* vertex length from center */
	/* and dot product with group center vector */
	for (vlen= 0.0, f = 0; f < fdi; f++) {
		double tt = v[f] - gcent[f];
		vlen += tt * tt;
		dot += (gc1[f] - gcent[f]) * tt;
	}

	/* If the groupls are not in the same direction, return false */
	if (dot < 0.0)
		return 0;
  
	vlen = sqrt(vlen);
	
	if (vlen > 1e-6)
		scale = cc1/vlen;
	else
		scale = 1.0;

	for (f = 0; f < fdi; f++)
		sv[f] = (scale * (v[f] - gcent[f])) + gcent[f];

	/* Distance from scaled vertex to group center */
	for (w = 0.0, f = 0; f < fdi; f++) {
		double tt = sv[f] - gc1[f];
		w += tt * tt;
	}
	w = sqrt(w);

	if (w <= (dw1 + EPS))
		return 1;

	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Given a pointer to a bxcell, use the ->tlist to fill in the corresponding nnrev[]. */
/* Expand and share lists with nearby nnrev[] cells if they are similar. */
static void create_nnrev_list(
rspl *s,
bxcell *tx,			/* Target nnrev[] cell */
bxcell *ss,			/* Head of solution list of surface nnrev[] cells */
double emax			/* smallest emax in solution list */
) {
	int i, j;
	bxcell *bx;
	int *dp = NULL, *sp;
	double *eminlist;
	unsigned int hashk;

	DBG(("create_nnrev_list: di %d target cell ix %d co[] %s emax = %f\n",s->di, tx->ix,debPiv(s->fdi, tx->gc),emax));

	/* Update tx->ss and tx->sdist with best in tlist for future */
	/* searches from this surface cell. */
	tx->sdist = 1e200;
	for (bx = ss; bx != NULL; bx = bx->tlist) {
//printf("~1 checking ix %d\n",bx->ix);
		if (bx->sdist < tx->emin) {
			tx->ss = bx;
			tx->sdist = bx->emin;
			DBG((" Set target ss to ix %d, emin %f\n",tx->ss->ix, tx->sdist));
//printf(" Set target ss to ix %d, emin %f\n",tx->ss->ix, tx->sdist);
		}
	}
//printf("~1 set sdist\n");

#ifdef DEBUG2
	{
		int tot = 0;
		printf(" Initial fwd list from following surface cells:\n");
		for (bx = ss; bx != NULL; bx = bx->tlist) {
			if (bx->emin <= emax) {
				if (bx->sl == NULL)
					error("rev create_nnrev_list: found empty surface bxcell");
				printf("   ix %d co %s fwd count %d\n",bx->ix,debPiv(s->fdi, bx->gc),bx->sl[1]-3);
				tot += bx->sl[1]-3;
			}
		}
		printf(" Total fwd cells = %d\n",tot);
	}
#endif

	/* Create an initial list of fwd cells from all bxcells */
	/* on the solution list that have ->emin <= emax */
	for (bx = ss; bx != NULL; bx = bx->tlist) {
//printf("~1 checking solution cell ix %d co %s\n",bx->ix,debPiv(s->fdi, bx->gc));
		if (bx->emin <= emax) {
//printf("~1 solution cell has emin %f < emax %f\n",bx->emin, emax);
			sp = bx->sl;
			if (sp == NULL)
				error("rev create_nnrev_list: found empty surface bxcell %d",ss->ix);
			for (sp += 3; *sp != -1; sp++)
				add2indexlist(s, &dp, *sp, 0);
		}
	}

	if (dp == NULL)
		error("create_nnrev_list got NULL new list\n");

#ifdef DEBUG2
	printf(" Initial fwd list (length %d, alloc %d):\n",dp[1]-3,dp[0]);
	for (i = 3; dp[i] != -1; i++)
		printf("  %d: ix %d\n",i-3,dp[i]);
#endif

	/* Sort the list into ascending order */
#define 	HEAP_COMPARE(A,B) (A < B)
	HEAPSORT(int, dp + 3, dp[1]-3)
#undef 		HEAP_COMPARE

#ifdef DEBUG2
	printf(" After sorting:\n");
	for (i = 3; dp[i] != -1; i++)
		printf("  %d: ix %d\n",i-3,dp[i]);
#endif

	/* Delete any duplicates */
	for (i = 3, j = i+1; ; j++) {
		if (dp[i] != dp[j])
			dp[++i] = dp[j];
		if (dp[j] == -1)
			break;
	}
	dp[1] = i;

#ifdef DEBUG2
	printf(" After de-duplication (length %d, alloc %d):\n",dp[1],dp[0]);
	for (i = 3; dp[i] != -1; i++)
		printf("  %d: ix %d\n",i-3,dp[i]);
#endif

	/* Filter fwd cells against emin/emax. */
	/* (Don't bother for 1D, as there's no point in filling up the cache */
	/* at this point, since 1D ins't participating in RAM management ?) */
	if (s->fdi > 1) {
		/* Allocate a temporary array to hold fwd cell emin */
		if ((eminlist = (double *) rev_malloc(s, (dp[1]-3) * sizeof(double))) == NULL)
			error("rspl malloc failed - rev create_nnrev_list emin array");
		INCSZ(s, (dp[1]-3) * sizeof(double));
	
		for (i = 0; i < (dp[1]-3); i++)
			eminlist[i] = 1e200;
		
		/* Get an fxcell for each fwd index, and compute emin & emax for this target. */
		/* Tracl smallest maximum and record each fxcell emin */
		emax = 1e200;
		for (i = 3; dp[i] != -1; i++) {
			fxcell *fc;
			double em, ex;
	
			fc = get_fxcell(s->rev.sb, dp[i], 1);
	
			eminlist[i-3] = nn_grpgrp_est(s, &ex, &fc->g, &tx->g);
			if (ex < emax)
				emax = ex;
			
			unget_fxcell(s->rev.cache, fc);
		}
	
#ifdef DEBUG2
		printf(" Smallest emax = %f\n",emax);
		for (i = 3; dp[i] != -1; i++)
			printf("  %d: ix %d, emin %f\n",i-3,dp[i],eminlist[i-3]);
#endif
	
		/* Delete any fwd cells/indexes that have an emin > smallest emax */
		for (i = j = 3; dp[j] != -1; j++) {
			if (eminlist[j-3] <= emax)
				dp[i++] = dp[j];
		}
		dp[i] = -1;
		dp[1] = i;
	
		free(eminlist);
		DECSZ(s, sizeof(schbase));

#ifdef DEBUG2
		printf(" After removing too far cells (length %d, alloc %d):\n",i-3,dp[0]);
		for (i = 3; dp[i] != -1; i++)
			printf("  %d: ix %d\n",i-3,dp[i]);
#endif
	}

	/* If the size of the list has reduced substatially, reclaim some memory */
	if ((dp[1]+1) <= (dp[0]/2)) {
		int ll = dp[0];
		while (ll > (dp[1]+1))
			ll /= 2;
		ll *= 2;
		DBG((" Reducing list allocation from %d to %d entries\n",dp[0],ll));
		DECSZ(s, (dp[0] - ll) * sizeof(int));
		if ((dp = (int *) rev_realloc(s, dp, sizeof(int) * ll)) == NULL)
			error("rspl realloc failed - create_nnrev_list");
		dp[0] = ll;	/* New allocation */
	}

	/* Check if any neighbor lists are similar to the list we just created, */
	/* so that we can merge similar lists, greatly reducing memory usage */
	/* at the cost of slightly longer lists. */
	/* Don't do this if this is a super-cell. */
	/* [ This seems to increase nnrev fill time by about 5% ] */
	if (tx->scell == NULL) {
		DCOUNT(cc, MXRO, s->fdi, -1, -1, 2);	/* bwd neighborhood offset counter */
		int nn[MXRO];
		int shlim, lnlim;
		int sh, ln;
		int bnix = -1, *blist = NULL, bwhgt = 0x7ffffff, bsh, bln;
		int f, nix;

		/* Set limits of an acceptable match at 2% short, 15% long */
		/* This trades off list size against number of lists/memory */
		/* i.e. a 10% rise in average list length for a 100 x reduction in */
		/* number of lists. (vary lnlim for most effect) */
		shlim = (2 * (dp[1]-3) + 50)/100;
		lnlim = (15 * (dp[1]-3) + 50)/100;

		DC_INIT(cc);
		while (!DC_DONE(cc)) {

			nix = tx->ix;
			for (f = 0; f < s->fdi; f++) {
				nn[f] = tx->gc[f] + cc[f];
				if (nn[f] < 0 || nn[f] >= s->rev.res)
					break;					/* Out of bounds */
				nix += cc[f] * s->rev.coi[f];
			}
			if (nix == tx->ix)			/* Skip this cell */
				goto next_neighbor;

			/* If neighbor is in bounds and has a fwd cell list, */
			/* check what sort of match it is to this list */
			if (f >= s->fdi && s->rev.nnrev[nix] != NULL) {
				int *np = s->rev.nnrev[nix];
				int *shrecs = getsharelist(s, np);

				if (shrecs != NULL) {
					if (shrecs[2] == tx->ix)		/* Already looked at this list */
						goto next_neighbor;
					shrecs[2] = tx->ix;				/* Remember we've done this one */
				}

				/* See how much it is a super or sub-set */
//printf("~1 checking ix %d against nix %d\n",tx->ix, nix);
				if ((dp[1] - np[1]) > shlim
				 || (np[1] - dp[1]) > lnlim) {
					goto next_neighbor;				/* No possibility of being acceptable */
				}

				sh = ln = 0;
				for (j = i = 3; dp[i] != -1 || np[j] != -1;) {

//printf("1: dp[%d] %d - np[%d] %d\n",i,dp[i],j,np[j]);
					while (np[j] != -1 && (dp[i] == -1 || dp[i] > np[j])) {
						j++;
						ln++;
//printf("2: dp[%d] %d - np[%d] %d, ln %d\n",i,dp[i],j,np[j],ln);
						if (ln > lnlim)
							goto next_neighbor;		/* No possibility of being acceptable */
					}

					while (dp[i] != -1 && (np[j] == -1 || dp[i] < np[j])) {
						i++;
						sh++;
//printf("3: dp[%d] %d - np[%d] %d, sh %d\n",i,dp[i],j,np[j],sh);
						if (sh > shlim)
							goto next_neighbor;		/* No possibility of being acceptable */
					}

					while (dp[i] != -1 && np[j] != -1 && dp[i] == np[j]) {
						i++;
						j++;
//printf("4: dp[%d] %d - np[%d] %d\n",i,dp[i],j,np[j]);
					}
				}
//printf("~1 len %d, short %d, long %d\n",dp[1]-3,sh,ln);

				/* remember best similar list within our criteria */
				if (sh <= shlim && ln <= lnlim) {
					int whgt = 2 * sh + ln;
					if (whgt < bwhgt) {
						bnix = nix;
						blist = np;
						bwhgt = bwhgt;
						bsh = sh;
						bln = ln;
					}
				}
			}
		next_neighbor:;
			DC_INC(cc);
		}

		/* Got a list we want to share with */
		if (blist != NULL) {
			int *shrecs = NULL;
			int *exlist = NULL;

			DBG((" Found similar existing list (short %d, long %d)\n",bsh,bln));

#ifdef DEBUG2
			printf(" Similar list (length %d, alloc %d):\n",blist[1]-3,blist[0]);
//			for (i = 3; blist[i] != -1; i++)
//				printf("  %d: ix %d\n",i-3,blist[i]);
#endif
			/* If the neighbor list is not a super-set */
			if (bsh > 0) {

				/* But new list is superset of neighbor list */
				if (bln == 0) {
					DBG((" Using new list to share\n"));

					exlist = dp;			/* Use our new list */
					exlist[2] = blist[2];	/* Same sharers */
					dp = NULL;

					/* Free neighbor list */
					free_indexlist(s, &blist);

				/* Create superset list from new list and neighbor list */
				} else {
	
					DBG((" Creating superset list\n"));

					for (j = i = 3; dp[i] != -1 || blist[j] != -1;) {
		
						while (blist[j] != -1 && (dp[i] == -1 || dp[i] > blist[j])) {
							add2indexlist(s, &exlist, blist[j], 0);
							j++;
						}
		
						while (dp[i] != -1 && (blist[j] == -1 || dp[i] < blist[j])) {
							add2indexlist(s, &exlist, dp[i], 0);
							i++;
						}
		
						while (dp[i] != -1 && blist[j] != -1 && dp[i] == blist[j]) {
							add2indexlist(s, &exlist, dp[i], 0);
							i++;
							j++;
						}
					}
	
					exlist[2] = blist[2];		/* Same sharers */
	
					/* Free neighbor list */
					free_indexlist(s, &blist);

					/* Done with list we created for this nnrev[] */
					free_indexlist(s, &dp);
				}

			} else {
				DBG((" Using existing list to share\n"));

				exlist = blist;		/* blist is already a super-set */
				blist = NULL;		/* Done with neighbor list */

				/* Done with list we created for this nnrev[] */
				free_indexlist(s, &dp);
			}
	
#ifdef DEBUG2
			printf(" Superset list nnrev[%d] (length %d, alloc %d):\n",tx->ix,exlist[1]-3,exlist[0]);
//			for (i = 3; exlist[i] != -1; i++)
//				printf("  %d: ix %d\n",i-3,exlist[i]);
#endif
//if (s->fdi > 1 && (tx->ix == 19054 || tx->ix == 19055)) { 
//printf(" Superset list nnrev[%d] (length %d, alloc %d):\n",tx->ix,exlist[1]-3,exlist[0]);
//for (i = 3; exlist[i] != -1; i++)
//	printf("  %d: ix %d\n",i-3,exlist[i]);
//}

			/* If this list has not been shared before, create share record for it */
			if (getsharelist(s, exlist) == NULL) 
				add2sharelist(s, bnix, exlist);

			/* Add this cell as a sharer */
			add2sharelist(s, tx->ix, exlist);

			/* Update pointers for all sharers of this (possibly new) list */
			shrecs = getsharelist(s, exlist);
//printf("Number shared now %d\n", shrecs[1]-3);
			for (i = 3; shrecs[i] != -1; i++) {
				s->rev.nnrev[shrecs[i]] = exlist;
			}
	
		} else {
			DBG((" no matching existing list\n"));
//printf(" no matching existing list\n");
	
			/* Put list in place for target nnrev[]*/
			s->rev.nnrev[tx->ix] = dp;
		}
	} else {

		if (tx->scell != NULL) {

			/* Put list in place for all nnrev[]'s covered by super-cell */
			for (sp = tx->scell + 3; *sp != -1; sp++) {  

				/* Add this cell as a sharer */
				add2sharelist(s, *sp, dp);

				s->rev.nnrev[*sp] = dp;
			}

		} else {
			/* Put list in place for target nnrev[]*/
			s->rev.nnrev[tx->ix] = dp;
		}
	}

	DBG(("create_nnrev_list done, total fwd cells = %d\n",s->rev.nnrev[tx->ix][1]-3));
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* This is the routine used to fill nnrev[] cells on demand, */
/* because s->rev.fastsetup is set. */
/* This is similar to the code used in the normal case, except */
/* we search the rev[] bxcell list rather than use the bxcell surface list */
static void fill_nncell(
	rspl *s,
	int *co,	/* Integer coords of cell to be filled */
	int ix		/* Index of cell to be filled */
) {
	int f, fdi = s->fdi;
	DCOUNT(gg, MXRO, fdi, 0, 0, s->rev.res);	/* search seed coordinate */
	int i, six = -1, nn[MXRO];
	double bdist = 1e200;
	bxcell *tx, *ss;
	DCOUNT(cc, MXRO, fdi, -1, -1, 2);	/* bwd neighborhood offset counter */
	int nix;						/* Neighbor offset index */
	bxcell *xlist = NULL;			/* Linked list of cells being searched */
	bxcell *xlistend = NULL;		/* Last item on xlist */
	bxcell *tlist;					/* Linked list of cells being considered as soln. */
	double emax;					/* Current smallest estimated max weigted distance */

	DBG(("fill_nncell: (triggered on-demand)\n"));
	
	/* Allocate the bxcell hash index */
	create_surfhash(s);

	/* Locate a starting search cell. */
	/* We use a simple full search of rev[] for the cell */
	/* closest to our target. */
	DC_INIT(gg);
	for (i = 0; i < s->rev.no; i++) {
		if (s->rev.rev[i] != NULL) {
			double dist;
			for (dist = 0.0, f = 0; f < fdi; f++) {
				double tt = co[f] - gg[f];
				dist += tt * tt;
			}
			if (dist < bdist) {
				bdist = dist; 
				six = i;
				for (f = 0; f < fdi; f++)
					nn[f] = gg[f];
			}
		}
		DC_INC(gg);
	}
	if (six < 0)
		error("fill_nncell: rev[] is empty");

	/* Create search seed cell */
	ss = new_bxcell(s, six, nn, NULL, 0.0, NULL);
	add_bxcell_hash(s, ss);

	/* Create a target cell */
	tx = new_bxcell(s, ix, co, ss, 0.0, NULL);
	add_bxcell_hash(s, tx);

	DBG((" Target ix = %d, co[] %s\n",ix,debPiv(fdi, tx->gc)));
	DBG((" Search start ix = %d, co[] %s\n",six,debPiv(fdi, ss->gc)));
//printf(" Target ix = %d, co[] %s\n",ix,debPiv(fdi, tx->gc));
//printf(" Search start ix = %d, co[] %s\n",six,debPiv(fdi, ss->gc));

	emax = 1e200;				/* Smallest emax */
	ss->tix = tx->ix;			/* Mark this cell as being in search list */
	
	/* Make start cell the only entry in the search list */
	ss->xlist = NULL;
	xlist = ss;
	xlistend = ss;

	/* Clear the solution list */
	tlist = NULL;

	/* While there are cells to search for solutions */
	while (xlist != NULL) {
		double em, ex;

		ss = xlist;					/* Remove next search cell from linked list */
		xlist = xlist->xlist;

		/* Check if this cell could be in solution */
		em = nn_grpgrp_est(s, &ex, &tx->g, &ss->g);
		ss->emin = em;

		DBG(("Searching rev[%d] co %s, em %f, ex %f\n",ss->ix, debPiv(s->fdi, ss->gc), em, ex));
//printf("Searching rev[%d] co %s, em %f, ex %f\n",ss->ix, debPiv(s->fdi, ss->gc), em, ex);
	
		if (em < emax) {		/* Yes */

			/* Add it to the solution list */
			ss->tlist = tlist;
			tlist = ss;

			// copy rev[] list to ss->sl
			copy_indexlist(s, &ss->sl, s->rev.rev[ss->ix]);

			DBG(("Adding %d to solution list\n",ss->ix));
//printf("Adding %d to to solution list\n",ss->ix);

			/* Update smallest maximum */
			/* (Will cull existing bxcell solutions with emin > emax later) */
			if (ex < emax)
				emax = ex;

			/* Explore all neighbours, and add any surface cells that haven't been */
			/* searched for this target yet. */ 
			DC_INIT(cc);
			while (!DC_DONE(cc)) {
				bxcell *nbx;

				nix = ss->ix;
				for (f = 0; f < fdi; f++) {
					nn[f] = ss->gc[f] + cc[f];
					if (nn[f] < 0 || nn[f] >= s->rev.res)
						break;					/* Out of bounds */
					nix += cc[f] * s->rev.coi[f];
				}
				if (f < fdi || nix == ss->ix) {
//printf("Rejecting search neigbor co %s because out of bounds or current cell\n",debPiv(s->fdi,nn));
					goto next_neighbor;
				}

				/* Can only search within filled rev[] cells */
				if (s->rev.rev[nix] == NULL) {
					goto next_neighbor;
				}

				/* If neighbor is in bounds, and a surface bxcell*/
				{
					/* Make sure we have a bxcell for the neighbor */
					if ((nbx = get_surface_bxcell(s, nix)) == NULL) {
						nbx = new_bxcell(s, nix, nn, NULL, 0.0, NULL);
						add_bxcell_hash(s, nbx);
					}

					/* If not already in search list */
					if (nbx->tix != tx->ix) {
//						DBG(("Adding search neigbor nnrev[%d] co %s to search list\n",nbx->ix, debPiv(s->fdi, nbx->gc)));
//printf("Adding search neigbor nnrev[%d] co %s to search list\n",nbx->ix, debPiv(s->fdi, nbx->gc));
						/* Add neigbor to end of search list */
						nbx->tix = tx->ix;		/* Is now in search list */
						nbx->xlist = NULL;
						if (xlist == NULL)
							xlist = nbx;
						else
							xlistend->xlist = nbx;
						xlistend = nbx;
					}
//else 
//printf("Rejecting search neigbor nnrev[%d] co %s because already in list\n",nbx->ix, debPiv(s->fdi, nbx->gc));
				}
			next_neighbor:;
				DC_INC(cc);
			}
		}
//else
//printf("Rejected rev[%d] co %s, because em %f >= emax %f\n",ss->ix, debPiv(s->fdi, ss->gc), em, emax);
	}

	if (tlist == NULL)
		error("fill_nncell: search for rev[] cells failed");

//printf("Got solution list, filling in nnrev[] cell\n");
	/* Create the nnrev[] list from the candidate bxcell solutions */
	create_nnrev_list(s, tx, tlist, emax);  

//printf("nnrev[%d] list length = %d\n",tx->ix,s->rev.nnrev[tx->ix][1]-3);

	/* Free up bxcell hash index and all bxcell's we've created */
	free_surfhash(s, 1);

	DBG(("fill_nncell done\n"));
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Associated sub-simplex tables. For a given base vertex with a given fwd */
/* access flags FLV(), create is a pointer to a list of sub-simplex verticies */
/* offset from the base vertex that the base vertex is part of, in all possible */
/* directions. We do all possible directions to make the bxcell triangle */
/* search symetrical, and searching triangles using points outside */
/* the bxcell list seems to actually speed it up (by culling more effectively) */

typedef struct {
	int pos;			/* nz if this is a ssimplex that is only in the +ve direction */ 
	int ee;
	int goffs[MXDI+1];	/* Offsets to sub-simplex verticies within grid in simplex order. */
} assinfo;

typedef struct {
	int sdi;			/* Sub dimensionality */
	int no;				/* Number of sub-simplexes in list */
	assinfo *ti;		/* Per sub-simplex info array */
} assdire;

#if (FL_BITS != 3)
#error FL_BITS is not 3!
#endif

/* init the triangle/edge directory list and related assinfo tables */
/* sdi = 2 for triangles, 1 for edges */
static void init_assdir(rspl *s, assdire **passdir, int sdi) {
	assdire *assdir;
	int i, j, k;
	int e, ee, di = s->di;
	int dirsize;
	DCOUNT(cc, MXRI, di, -1, -1, 2);	/* Clipping values for each dim */
	ssxinfo *xip;		/* Pointer to sub-simplex info structure */

	DBG(("init_assdir called, di = %d\n",di));

	dirsize = (1 << (FL_BITS * di));

	if ((assdir = (assdire *) rev_calloc(s, dirsize, sizeof(assdire))) == NULL)
		error("rspl malloc failed - assdir");
	INCSZ(s, dirsize * sizeof(assdire));

	assdir->sdi = sdi;
	xip = &s->rev.sspxi[sdi];

#ifdef NEVER
	printf("simplex dim %d:\n",xip->sdi);
	for (i = 0; i < xip->nospx; i++) {
		printf("offs = %s\n", debPiv(sdi+1, xip->spxi[i].offs));
		printf("goffs = %s\n", debPiv(sdi+1, xip->spxi[i].goffs));
	}
#endif

	/* For each possible clip combination */
	/* (where < 0 == clipping lower edge, > 0 == clipping upper edge */
	DC_INIT(cc);
	while (!DC_DONE(cc)) {
		int trilaloc, trillen;
		assinfo *trilist;

		/* Start a new table, allocate the maximum possible number of entries. */
		trilaloc = (1 << di) * xip->nospx;
		if ((trilist = (assinfo *) rev_calloc(s, trilaloc, sizeof(assinfo))) == NULL)
			error("rspl malloc failed - trilist");
		INCSZ(s, trilaloc * sizeof(assinfo));
		trillen = 0;

		/* For all cube directions from base, 0 = +ve, 1 = -ve */
		for (ee = 0; ee < (1<<di); ee++) {

			/* For all the sub-simplexes in a cube */
			for (i = 0; i < xip->nospx; i++) {
				int gotbase = 0;

				/* Offset the sub-simplex by the direction, and check that the */
				/* base vertex is part of it. */
				trilist[trillen].ee = ee;
				trilist[trillen].pos = (ee == 0);
				for (j = 0; j < (sdi+1); j++) {
					trilist[trillen].goffs[j] = xip->spxi[i].goffs[j] - s->g.hi[ee];
					if (trilist[trillen].goffs[j] == 0)	/* Base vertex is present */
						gotbase = 1;
				}
				if (!gotbase) {
					continue;
				}

				/* See if the direction of each vertex of the sub-simplex is */
				/* compatible with the clipping. */
				for (j = 0; j < (sdi+1); j++) {
					for (e = 0; e < di; e++) {
						if (xip->spxi[i].offs[j] & (1<<e)) {

							if ((cc[e] < 0 && (ee & (1<<e)) != 0)	
							 || (cc[e] > 0 && (ee & (1<<e)) == 0)) {
								break;		/* not compatible */
							}
						}
					}
					if (e < di) { 		/* Not compatible */
						break;
					}
				}
				if (j < (sdi+1)) {
					continue;			/* Not compatible */
				}

				/* We end up with aliases due to the sspxi having all */
				/* sub-simplexes within a cube, so see if we already */
				/* created this one. */ 
				for (k = 0; k < trillen; k++) {
					for (j = 0; j < (sdi+1); j++) {
						if (trilist[k].goffs[j] != trilist[trillen].goffs[j])
							break;
					}
					if (j >= (sdi+1))
						break;		/* Redundant - don't add this point */
				}
				if (k < trillen) {
					continue;		/* Skip redundant combination */
				}

//printf(" Clip %s off %d tri %d goffs = %s\n", debPiv(di, cc), ee, trillen, debPiv(sdi+1, trilist[trillen].goffs));
				trillen++;
			}
		}

//printf("Got %d triangles for cc %s\n", trillen, debPiv(di, cc));

		/* Add table to all matching combination of FLV() */
		for (i = 0; i < dirsize; i++) {
			for (e = 0; e < di; e++) {
				int fl = (i >> (3 * e)) & 7;
				if (!								/* NOT: */
				    ((cc[e] > 0 && fl == 0) 			/* Top edge clip and on top edge */
				  || (cc[e] < 0 && fl == 4) 			/* Bottom edge clip and on bottom edge */
				  || (cc[e] == 0 && fl != 0 && fl != 4)))	/* No clipping and in middle */ 
					break;	/* Not a match */
			}
			if (e >= di) {	/* Table matches this FLV() */
				assdir[i].no = trillen;
				assdir[i].ti = trilist;
			}
		}
		DC_INC(cc);		/* Next clip combination */
	}

#ifdef NEVER
	/* Check that there is a list for every flag value */
	for (i = 0; i < dirsize; i++) {
		if (assdir[i].no == 0)
			error("init_assdir has fl %d entry with no sub-simplexes",i);
		else
			printf("fl %d has %d triangles\n",i,assdir[i].no);
	}
#endif

	*passdir = assdir;
}

static void free_assdir(rspl *s, assdire *assdir) {
	int i, j;
	int e, ee, di = s->di;
	int sdi = assdir->sdi;
	int dirsize = (1 << (FL_BITS * di));
	int trilaloc = (1 << di) * s->rev.sspxi[sdi].nospx;

	for (i = 0; i < dirsize; i++) {
		assinfo *trilist;

		if ((trilist = assdir[i].ti) == NULL)
			continue;

		/* Free all aliases of list */
		for (j = i; j < dirsize; j++) {
			if (trilist == assdir[j].ti) {
				assdir[j].ti = NULL;
			}
		}

		free(trilist);
		DECSZ(s, trilaloc * sizeof(assinfo));
	}
	free(assdir);
	DECSZ(s, dirsize * sizeof(assdire));
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Solve the 2x2 simultaneous linear equations A.X = B */
static int solve_se_2x2(double **ta, double *tb) {
	double b[2] = { tb[0], tb[1] };
	double det;
	int rv;

	det = (ta[0][0] * ta[1][1] - ta[0][1] * ta[1][0]);

	if (fabs(det) < 1e-20)
		return 1;

	det = 1.0/det;
	tb[0] = det * ( ta[1][1] * b[0] - ta[0][1] * b[1]);
	tb[1] = det * (-ta[1][0] * b[0] + ta[0][0] * b[1]);

	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef CHECK_NNLU

/* Debug code */

static int debug(int ix) {
	if (
		(ix == 619 || ix == 618 || ix == 329)
	 || (ix == 330 || ix == 329 || ix == 312)
	 || (ix == 619 || ix == 329 || ix == 312)
	 || (ix == 329 || ix == 312 || ix == 23)
	 || (ix == 329 || ix == 23 || ix == 22)
	 || (ix == 40 || ix == 23 || ix == 22)
	)
		return 1;
	return 0;
}

static int debug2(int *ix) {
	if (
		(ix[0] == 619 && ix[1] == 618 && ix[2] == 329)
	 || (ix[0] == 330 && ix[1] == 329 && ix[2] == 312)
	 || (ix[0] == 619 && ix[1] == 329 && ix[2] == 312)
	 || (ix[0] == 329 && ix[1] == 312 && ix[2] == 23)
	 || (ix[0] == 329 && ix[1] == 23 && ix[2] == 22)
	 || (ix[0] == 40 && ix[1] == 23 && ix[2] == 22)
	)
		return 1;
	return 0;
}

#endif /* CHECK_NNLU */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef REVVRML		/* Plotting routine declarations */

static void plot_bxfwcells(rspl *s, int dobxcells, int dofwcells, int dofwlabels);

static void plot_tri_check(rspl *s, int dobxcells, int dowait, bxcell *bx, int vtxix,
	int trii, int triix[3], int nvtxix, int sorv, int wsrv, int shdwd,
	double v[MXRI+1][MXRO], double de[MXRO], double pv[MXRO], double xv[MXRO]);

static void plot_vtx_surface(rspl *s, int dovtxlabels, int dodeleted, int doadded,
	int dopres, int dooil, int dobxcells, int dowait, vtxcache *vc, assdire *edgdir); 

static void plot_touched_bxcells(rspl *s, int bxix);

static void plot_fxcell_surface(rspl *s, int dofclabels, int dobxcells, int dowait);

#endif /* REVVRML */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*
	The basic strategy to thin the gamut surface as much as possible
	to reduce the nnrev[] list size for best memory consuption and
	rev lookup speed relies on being able to decide if a vertx
	is inside or on the surface of the gamut. A simple and definitive
	topological rule hasn't been forthcoming, so a simpler heursitic
	of visiblilty from a singe internal "focal" point is currently used.
	calc_ocent() attempts to choose a point with best visibility of
	all the gamut surfaces, since any self-shadowing results in
	gamut surface holes.

	Improvements would be to create per-axis mappings (and separate
	the shadow vertex locations from the real ones) to re-shape
	the gamut into a square as much as possible.
	Multiple external focal points could be used,
    a vertex being shadowed only when it can't be "seen" by any
	external focal point. It's hard to figure how to make the latter
    fast enough to be useful though, unless some breakthrough
    in the algorithm or spatial data structure can be developed.

	The code is still slower than desired. A possible avenue for
	improving the thinning would be to add an explicit triangle
	structure (similar to gamut ?), add a suitable spatial
	accelleration structure for shadow testing (BSP tree ??),
	and build the gamut surface incrementally from existing
	furthest points.

	(Have loop re-orderings been exausted ? i.e. can overlap
	 triangle processing "Do a first pass for each test vertex,
	 testing against just the triangles that are associated with
	 it's triangle" be used for main shadowing testing ?)
*/

/* Struct to slice locus points */
struct _slpoint {
	double v[MXRO];			/* Point location */
	double wrad;			/* Weighted radius */
	double rad;				/* Distance from ccent along slice */
	double minrad;			/* Minimum distance from ccent */

	double cvec[MXRO];		/* Vector from this point to ccent */
	double len;				/* Length of segment, -1 if no good */
	double trad;			/* Trial center point to this v[] radius */
}; typedef struct _slpoint slpoint;

/* Center finding context */
struct _ocenctx {
	rspl *s;
	int ares;					/* angle resolution */ 
	slpoint *p[MXRO];			/* Slice locus points */
	double ccent[MXRO];			/* Construction center point */
	double ret;					/* return value */
	int oog;					/* flag set if center is out of gamut */
	int debug;
}; typedef struct _ocenctx ocenctx;

/* Given a set of slice locus points and a proposed center point, */
/* compute the weighted average of the orthogonality of the point */
/* to each locus line segment. (smaller is better) */
/* (Used for optimizing the focal/center point.) */
static double aorthog(void *_ctx, double *cent) {
	ocenctx *ctx = (ocenctx *) _ctx;
	rspl *s = ctx->s;
	int f, ff, fdi = s->fdi;
	int aa, ares = ctx->ares;
	double tcent[MXRO];
	double ang, aang = 0.0;
	int naang = 0;
	
	ctx->oog = 0;

	if (ctx->debug) printf("aorthog called with cent %s\n",debPdv(fdi,cent));

	for (ff = 0; ff < fdi; ff++) {
		if (ctx->debug) printf(" Axis %d\n",ff);

		for (f = 0; f < fdi; f++)
			tcent[f] = cent[f];
		/* Flatten the points to lie on the notional center */
		tcent[ff] = ctx->ccent[ff]; 

		for (aa = 0; aa < ares; aa++) {
			double trad, nrad;
			double cvec[MXRO], dot;

			if (ctx->p[ff][aa].len < 0.0)
				continue;

			if (aa == 0) {
				/* Compute normalize vector from cent to locus to this point */ 
				trad = 0.0;
				for (f = 0; f < fdi; f++) {
					double tt = tcent[f] - ctx->p[ff][aa].v[f];
					trad += tt * tt;
				}
				trad = sqrt(trad);
			} else {		/* Was computed by previous */
				trad = ctx->p[ff][aa].trad;
			}

			/* Compute normalize vector from tcent to locus to next point */ 
			nrad = 0.0;
			for (f = 0; f < fdi; f++) {
				cvec[f] = tcent[f] - ctx->p[ff][aa+1].v[f];
				nrad += cvec[f] * cvec[f];
			}
			nrad = ctx->p[ff][aa+1].trad = sqrt(nrad);
			
			/* Normalized difference in distance over length */
			/* Compute dot product of cv and segment vector */
			/* (ang range 0.0 .. 1.0 */
			ang = fabs(trad - nrad)/ctx->p[ff][aa].len;
			if (ang > 1.0)
				ang = 1.0;

			if (ctx->debug) printf("  aa %d: trad %f nrad %f, diff %f, len %f, ang %f\n",aa,trad,nrad,fabs(trad - nrad),ctx->p[ff][aa].len,ang);

			/* Compute dot of next point vector from trial center */
			/* with vector from construction center, to detect */
			/* if the trial has wandered outside of gamut. */
			dot = 0.0;
			for (f = 0; f < fdi; f++)
				dot += cvec[f] * ctx->p[ff][aa+1].cvec[f];

			if (dot < 0.0) {
				if (ctx->debug) printf("  dot is %f\n",dot);
				ang = 50.0;				/* Big value */
				ctx->oog = 1;
			} else {
				ang = pow(ang, 50.0);	/* Weight high angles */
			}
			aang += ang;
			naang++;
		}
	}
	aang /= (double)naang;
	
	if (ctx->debug) printf(" returning %f\n",aang);

	ctx->ret = aang;

	return aang;
}

/* Determine a gamut center point, for surface triangle shadow testing. */
/* We assume that rev[] has been setup. */
/* The idea is to locate a point that best "sees" all internal */
/* surface of the gamut. */
static void calc_ocent(rspl *s) {
	int i, j, aa, mm;
	int e, ee, di = s->di;
	int f, ff, fdi = s->fdi;
	int rgres = s->rev.res;			/* number of bwd cells */
	double minmax[2][MXRO][MXRO];	/* Range min/max points for each axis */
	float *gp, *ep;
	int midix[MXDO];				/* Middle rev[] index */ 
	double mid[MXRO];				/* Middle of midix[] */
	double ss[MXRO];
	ocenctx ctx;					/* Context */
	double atanscale;
	int ici, nici;

	/* Scan the forward array for the min and max points of each axis */
	for (f = 0; f < fdi; f++) {
		minmax[0][f][f] = 1e200;
		minmax[1][f][f] = -1e200;
	}

	/* Scan the fwd Grid for min/max values */
	for (gp = s->g.a, ep = s->g.a + s->g.no * s->g.pss; gp < ep; gp += s->g.pss) {
		for (ff = 0; ff < fdi; ff++) {
			if (minmax[0][ff][ff] > gp[ff]) {
				for (f = 0; f < fdi; f++)
					 minmax[0][ff][f]= gp[f];
			}
			if (minmax[1][ff][ff] < gp[ff]) {
				for (f = 0; f < fdi; f++)
					 minmax[1][ff][f] = gp[f];
			}
		}
	}

	if (fdi == 1) {
		for (f = 0; f < fdi; f++)
			s->rev.ocent[f] = 0.5 * (minmax[0][f][f] + minmax[1][f][f]);
		DBG(("calc_ocent: got 1d ocent = %s\n",debPdv(fdi,s->rev.ocent)));
		return;
	}

	/* Aprox. mid point of gamut from average of min/max points */
	for (f = 0; f < fdi; f++)
		ctx.ccent[f] = 0.0;
	for (ff = 0; ff < fdi; ff++) {
		for (f = 0; f < fdi; f++) {
			if (f == ff)
				continue;
			for (mm = 0; mm < 2; mm++)
				ctx.ccent[f] += minmax[mm][ff][f];
		}
	}
	for (f = 0; f < fdi; f++)
		s->rev.ocent[f] = ctx.ccent[f] /= ((fdi-1) * 2.0);

	DBG(("calc_ocent: initial ccent = %s\n",debPdv(fdi,ctx.ccent)));
//printf("calc_ocent: initial ccent = %s\n",debPdv(fdi,ctx.ccent));

	/* If it's all to hard ... */
	if (fdi != 3) {
		return;
	}

	/* Index of data mid point in rev[] grid */
	for (f = 0; f < fdi; f++) {
		midix[f] = (int)((ctx.ccent[f] - s->rev.gl[f])/s->rev.gw[f] + 0.5);
		mid[f] = (midix[f]+0.5) * s->rev.gw[f] + s->rev.gl[f];
	}
//printf("calc_ocent: mid point = %s\n",debPdv(fdi,mid));

	/* Array for each slice values at angle (+ repeat at end) */
	ctx.debug = 0;
	ctx.s = s;
	ctx.ares = (rgres + 1) & ~1;	/* Make even so that there is an opposite angle */
	if (ctx.ares < 6)
		ctx.ares = 6;
	else if (ctx.ares > 20)
		ctx.ares = 20;
//printf("  ocent ares %d\n",ctx.ares);
	atanscale = ctx.ares/(2.0 * DBL_PI);
	for (ff = 0; ff < fdi; ff++) {
		if ((ctx.p[ff] = (slpoint *)rev_calloc(s, ctx.ares+1,sizeof(slpoint))) == NULL)
			error("rspl malloc failed - calc_ocent arrays");
		INCSZ(s, (ctx.ares+1) * sizeof(slpoint));
	}

	/* Use 5 passes to locate a more reliable initial center point */
	for (nici = 10, ici = 0; ici < nici; ici++) {
//printf("  locating center point iter %d\n",ici);

		/* Set initial radius values */
		for (ff = 0; ff < fdi; ff++) {
			for (aa = 0; aa < ctx.ares; aa++) {
				ctx.p[ff][aa].wrad = -1.0;
				ctx.p[ff][aa].rad = -1.0;
				ctx.p[ff][aa].minrad = 1e38;
			}
		}
				
		/* Take three slices through the rev[] array, plotting */
		/* the maximum circumference for the slice */

		/* For the axis we're slicing */
		for (ff = 0; ff < fdi; ff++) {
			FCOUNT(cc, MXRO, 3);	/* Counter through bwd cells */
			int start[3], endp1[3];
			double vv[MXRO];
			int aa;
			
//printf(" slice axis %d\n",ff);

			/* Setup "fat" slice range */
			for (f = 0; f < fdi; f++) {
				if (f == ff) {
					start[f] = midix[f]-1;
					if (start[f] < 0)
						start[f] = 0;
					endp1[f] = midix[f]+2;
					if (endp1[f] > rgres)
						endp1[f] = rgres;
				} else {
					start[f] = 0;
					endp1[f] = rgres;
				}
			}
			FRECONFA(cc, start, endp1);

//printf("  slice range %d - %d, %d - %d, %d - %d\n", start[0], endp1[0]-1, start[1], endp1[1]-1, start[2], endp1[2]-1);

			/* Scan this 3 thick, 2D slice of rev[] */
			FC_INIT(cc);
			while (!FC_DONE(cc)) {
				int ix;
				int slix[MXRO];		/* Indexes in slice direction + orthogonal */
				int *rp;

				/* Compute bx index */
				ix = 0;
				for (j = f = 0; f < fdi; f++) {
					ix += cc[f] * s->rev.coi[f];
					if (f != ff)
						slix[j++] = f;
				}
				slix[j++] = ff;
//printf("  bx %d, %d ix %d\n",cc[0],cc[1],cc[2],ix);

				if (s->rev.rev[ix] == NULL) {
//printf("  rev is empty\n");
					goto next_bx;
				}

				/* For all the cubes bx rev[] */
				for (rp = s->rev.rev[ix]+3; *rp != -1; rp++) {

					/* For each vertx of this cube */
					for (ee = 0; ee < (1<<di); ee++) {
						int vix = *rp + s->g.hi[ee];
						float *gp = s->g.a + vix * s->g.pss;	/* Pointer to float of fwd vertex */
						double fcb[MXRO];
						double x, y, z, wrad, rad, ang;

						/* Don't add over ink limit vertexes */
						if (s->limiten && gp[-1] > s->limitv) {
							continue;
						}

						for (f = 0; f < fdi; f++)
							fcb[f] = gp[f]; 

						/* (Don't) Convert output values to log values */
						/* logcomp(s, fcb, fcb, ctx.ccent); */

						/* Compute 2D radius and normalize */
						x = fcb[slix[0]] - ctx.ccent[slix[0]];
						y = fcb[slix[1]] - ctx.ccent[slix[1]];
						z = fcb[slix[2]] - ctx.ccent[slix[2]];
						/* wrad is "elipsoid" weighted radius in slice */
						wrad = x * x + y * y - 1.5 * z * z;
						wrad = sqrt(wrad < 0.0 ? 0.0 : wrad);
						rad = sqrt(x * x + y * y);
						if (rad < EPS || wrad < EPS)
							continue;

						/* Quantized angle this point is at */
						ang = atanscale * atan2(y, x);
						aa = (int)floor(ang);
						if (aa < 0)
							aa += ctx.ares;
						if (aa >= ctx.ares) 
							aa -= ctx.ares;

//printf("   slice %d vtx %f %f %f rad %f, ang %f aa %d\n", ff, fcb[0], fcb[1], fcb[2], rad, ang, aa);

						if (wrad > ctx.p[ff][aa].wrad) {
							ctx.p[ff][aa].wrad = wrad;
							ctx.p[ff][aa].rad = rad;

							/* Copy far point */
							for (f = 0; f < fdi; f++)
								ctx.p[ff][aa].v[f] = fcb[f];

							/* (don't) Flatten the points to lie on the notional center */
							/* ctx.p[ff][aa].v[ff] = ctx.ccent[ff]; */
						}

						/* Track min in case ccent is not within slice */
						if (rad < ctx.p[ff][aa].minrad) {
							ctx.p[ff][aa].minrad = rad;
						}
					}
				}
			  next_bx:;
				FC_INC(cc);
			}

			/* Repeat first in extra at end */
			ctx.p[ff][ctx.ares] = ctx.p[ff][0];		/* Structure copy */
		}

		/* Check if center point is within slice by looking for empty entries. */
		{
			double ccvec[MXRO];		/* Center correction vector */
			double ccount = 0.0;
			
			for (f = 0; f < fdi; f++)
				ccvec[f] = 0.0;

			for (ff = 0; ff < fdi; ff++) {
				for (aa = 0; aa < ctx.ares; aa++) {

//printf("   slice %d aa %d, vtx %s rad %f, irad %f\n", ff, aa, debPdv(fdi,ctx.p[ff][aa].v), ctx.p[ff][aa].rad, ctx.p[ff][aa].minrad);

					/* Either the grid is very sparse, or our center */
					/* is outside */
					if (ctx.p[ff][aa].rad < 0.0) {
						int oaa = aa + (ctx.ares/2);
						if (oaa >= ctx.ares) 
							oaa -= ctx.ares;

//printf("   oaa %d, vtx %s rad %f, irad %f\n", oaa, debPdv(fdi,ctx.p[ff][oaa].v), ctx.p[ff][oaa].rad, ctx.p[ff][oaa].minrad);

						/* If oposite side has an entry */
						if (ctx.p[ff][oaa].rad > 0.0) {
							double cor[MXRO];
							double prop = (3.0 * ctx.p[ff][oaa].minrad
							                   + ctx.p[ff][oaa].rad)/4.0;

							prop /= ctx.p[ff][oaa].rad;		/* Proportion of distance to v[] */

							for (f = 0; f < fdi; f++)
								cor[f] = prop * (ctx.p[ff][oaa].v[f] - ctx.ccent[f]);
//printf("   prop %f, corr %s\n",prop,debPdv(fdi,cor));
							for (f = 0; f < fdi; f++)
								ccvec[f] += cor[f];
							ccount++;
						}
					}
				}
			}

//printf("ccount %f\n",ccount);
			if (ccount > 0.0) {	/* Make adjustment */
				if (ici < (nici-1)) {
					for (f = 0; f < fdi; f++)
						ccvec[f] /= ccount;
//printf("Corecting center by %s\n",debPdv(fdi,ccvec));
					for (f = 0; f < fdi; f++)
						ctx.ccent[f] += ccvec[f];
//printf("cceny now %s\n",debPdv(fdi,ctx.ccent));
				} else {		/* Last round and correction needed */
					if (0.0 && s->verbose)
						fprintf(stdout, "%cFailed to locate aprox. gamut center\n",cr_char);
				}
			} else {
				break;		/* We're done */
			}
		}
	}

//printf("calc_ocent: refined ccent = %s\n",debPdv(fdi,ctx.ccent));

	/* Pre-compute point to point info to speed optimization */
	for (ff = 0; ff < fdi; ff++) {
		for (aa = 0; aa < ctx.ares; aa++) {
			double len = 0.0;

			for (f = 0; f < fdi; f++)
				ctx.p[ff][aa+1].cvec[f] = ctx.ccent[f] - ctx.p[ff][aa].v[f];

			for (f = 0; f < fdi; f++) {
				double tt = ctx.p[ff][aa+1].v[f] - ctx.p[ff][aa].v[f];
				len += tt * tt;
			}
			if (len < EPS)
				ctx.p[ff][aa].len = -1.0;
			else
				ctx.p[ff][aa].len = sqrt(len);
		}
	}

	/* Locate center point that maximised the orthogonallity to each */
	/* slice segment. This should maximize visibility of the inner of the */
	/* gamut surface, for shadow testing. */
	for (f = 0; f < fdi; f++)
		ss[f] = fabs(0.1 * (minmax[1][f][f] - minmax[0][f][f]));

//ctx.debug = 1;

	/* return 0 on sucess, 1 on failure due to excessive itterations */
	if (powell(NULL, fdi, s->rev.ocent, ss, 1e-3, 500, aorthog, (void *)&ctx, NULL, NULL)) {
		printf("calc_ocent powell failed\n");
		for (f = 0; f < fdi; f++)
			s->rev.ocent[f] = ctx.ccent[f];
	}

//ctx.debug = 1;

	/* Check result */
	aorthog(&ctx, ctx.ccent);

	/* Hmm. This isn't very reliable in detecting failure. */
	if (ctx.oog)
		printf("calc_ocent failed to return in-gamut focal point!\n");
			
//printf("Final angle = %f\n", ctx.ret);

#ifdef REVVRML		/* Plotting routine declarations */
	/* Diagnostic - dump the gamut slice locii */
	{
		vrml *wrl;
		double grey[3] = { 0.5, 0.5, 0.5 };
		double white[3] = { 1.0, 1.0, 1.0 };
		double red[3] = { 0.8, 0.1, 0.1 };
		double green[3] = { 0.1, 1.0, 0.1 };
		double blue[3] = { 0.1, 0.1, 0.8 };
		double *rgb[3] = { red, green, blue };

		wrl = new_vrml("section", 0, s->rev.probxyz ? vrml_xyz : vrml_lab);
		wrl->add_marker(wrl, s->rev.ocent, NULL, 1.0);

		/* Show vertex labels */
		for (ff = 0; ff < fdi; ff++) {
			char index[100];

			for (aa = 0; aa < ctx.ares; aa++) {
				if (ctx.p[ff][aa].rad > 0) {
					sprintf(index, "%d:%d",ff,aa);
					wrl->add_text(wrl, index, ctx.p[ff][aa].v, white, 1.0);
				}
			}
		}

		/* Axis we're slicing */
		for (ff = 0; ff < fdi; ff++) {
			int vix[100];

			for (aa = 0; aa < ctx.ares; aa++) {
				if (ctx.p[ff][aa].rad > 0)
					vix[aa] = wrl->add_vertex(wrl, 0, ctx.p[ff][aa].v);
			}
			vix[aa] = vix[0];

			for (aa = 0; aa < ctx.ares; aa++) {
				if (ctx.p[ff][aa].rad > 0
				 && ctx.p[ff][aa+1].rad > 0)
					wrl->add_col_line(wrl, 0, vix + aa, rgb[ff]);
			}
		}
		wrl->make_lines_vc(wrl, 0, 0.0);
		printf("Created %s\n",wrl->name);
		wrl->del(wrl);
	}
#endif /* REVVRML */

	/* Free up the context data */
	for (ff = 0; ff < fdi; ff++) {
		free(ctx.p[ff]);
		DECSZ(s, (ctx.ares+1) * sizeof(slpoint));
	}

	DBG(("calc_ocent: final ocent = %s\n",debPdv(fdi,s->rev.ocent)));
//printf("calc_ocent: final ocent = %s\n",debPdv(fdi,s->rev.ocent));
}

/* Create gamut surface linearization (surflin) transform. */
/* This is used by logcomp() to try and straighten out the */
/* device response so that the ocent is "visible" from */
/* any point on the surface. */ 
/* (We assume we are called at the correct point when bx->status == bx_uninit) */
static int calc_surflin(
rspl *s,
vtxcache *vc,		/* Vertexes */
assdire *edgdir		/* Edge lookup for vertex */
) {
	int i, j, g;
	int e, ee, di = s->di;
	int f, ff, fdi = s->fdi;
	vtxrec *vx, *nvx;
	int nitter, itter;

	double vxv[POW2MXRI][MXRO];		/* Overal fwd interp vertex values */
	int nvtx;
	cow *mpoints;
	int gres[MXDO];
	double min[MXDO], max[MXDO];
	double vmin[MXDO], vmax[MXDO];

	for (f = 0; f < fdi; f++)
		gres[f] = s->g.bres;		// ??

//printf("calc_surflin: rspl res %d\n",gres[0]);
	DBG(("calc_surflin: rspl res %d\n",gres[0]));

//printf("~1 gres = %d %d %d\n", s->g.res[0], s->g.res[1], s->g.res[2]);
//printf("~1 ci = %d %d %d\n", s->g.ci[0], s->g.ci[1], s->g.ci[2]);

	/* Lookup interpolation cube corners */
	for (ee = 0; ee < (1<<di); ee++) {
		int ix;
		float *fcb;
		for (ix = e = 0; e < di; e++) {
			if (ee & (1<<e))
				ix += s->g.ci[e] * (s->g.res[e] - 1);
		}
		fcb = s->g.a + ix * s->g.pss;
		for (f = 0; f < fdi; f++)
			vxv[ee][f] = fcb[f];
//printf("~1 cube corners %d = %f %f %f\n",ee, vxv[ee][0], vxv[ee][1], vxv[ee][2]);
	}

//printf("calc_surflin: counting number of vertexes:\n");
	/* Count the number of vertexes we may need */
	nvtx = 0;
	for (i = 0; i < vc->hash_size; i++) {
		for (vx = vc->hash[i]; vx != NULL; vx = vx->hlink) {
			nvtx++;
		}
	}
	DBG(("calc_surflin: %d mapping points\n",nvtx));
//printf("calc_surflin: %d mapping points\n",nvtx);

//printf("calc_surflin: computing goal values:\n");
	/* Currently the vertex vl = vv = output value of vertex. */
	/* Temporarily replace vv with the idealize (linear interp) output "goal" values. */
	for (i = 0; i < vc->hash_size; i++) {
		for (vx = vc->hash[i]; vx != NULL; vx = vx->hlink) {
			int tix;				/* Temp fwd cell index */
			double we[MXRI];		/* Vertex input position == 1.0 - Weight */
			double gw[POW2MXRI];	/* weight for each grid cube corner */
			double w;

			/* Compute this vertexes relative input position */
			for (tix = vx->ix, e = 0; e < di; e++) {
				int dix;
				dix = tix % s->g.res[e];
				tix /= s->g.res[e];
				we[e] = (double)dix/(s->g.res[e]-1.0);
			}
			
			/* Compute corner weights needed for interpolation */
			gw[0] = 1.0;
			for (e = 0, g = 1; e < di; g *= 2, e++) {
				for (j = 0; j < g; j++) {
					gw[g+j] = gw[j] * we[e];
					gw[j] *= (1.0 - we[e]);
				}
			}
		
			/* Linear interpolated output values */
			w = gw[0];
			for (f = 0; f < fdi; f++)			/* Base of cube */
				vx->vv[f] = w * vxv[0][f];
		
			for (g = 1; g < (1<<di); g++) {	/* For all other corners of cube */
				w = gw[g];
				for (f = 0; f < fdi; f++)
					vx->vv[f] += w * vxv[g][f];
			}
		}
	}

	/* Now itteratively adjust the vl values to better match the scaled */
	/* relative positions of the goal values. */

	/* Go through them again to get every line they are part of */
	/* (We're assuming we need exponentially more itters with finer point */
	/* spacing ?) */
	nitter = (int)(0.06 * s->g.bres * s->g.bres + 0.5);
	if (nitter < 1)
		nitter = 1;
	for (itter = 0; itter < nitter; itter++) {

		DBG(("calc_surflin: maping itter %d\n",itter));
//printf("calc_surflin: maping itter %d\n",itter);

		for (i = 0; i < vc->hash_size; i++) {
			for (vx = vc->hash[i]; vx != NULL; vx = vx->hlink) {
				assdire *edg;			/* Edge table */
				float *fp;
				int fl;
				double agrad, aorad;
				int nn;
				double scale;

				fp = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */
				fl = FLV(fp);		/* Edge flags for this vertex */
				edg = edgdir + fl;

				/* For vertexes at the end of all possible edges common with this vertex, */
				/* compute average radius */
				agrad = aorad = 0.0;
				for (j = 0; j < edgdir[fl].no; j++) {
					int eix;

					/* Index number of vertex other than the one we got it from */
					if (edg->ti[j].goffs[0] != 0)
						eix = vx->ix + edg->ti[j].goffs[0];
					else
						eix = vx->ix + edg->ti[j].goffs[1];

					if ((nvx = get_vtxrec(vc, eix)) != NULL) {
						double glen, olen;
						glen = olen = 0.0;
						for (f = 0; f < fdi; f++) {
							double tt;

							tt = vx->vv[f] - nvx->vv[f];
							glen += tt * tt;

							tt = vx->vl[f] - nvx->vl[f];
							olen += tt * tt;
						}
						glen = sqrt(glen);
						olen = sqrt(olen);
						agrad += glen;
						aorad += olen;
						nn++;
					}
				}

				if (nn == 0) {					/* Hmm. No neighors ? */
					vx->status = vtx_del;		/* Mark it as isolated */
					continue;
				}

				scale = aorad/agrad;	/* Local scale factor from goal to output */

				/* Reset the current vertex output value based on the relative */
				/* position of it in goal space */
				for (f = 0; f < fdi; f++)
					vx->vl[f] = 0.0;

				nn = 0;
				for (j = 0; j < edgdir[fl].no; j++) {
					int eix;

					/* Index number of vertex other than the one we got it from */
					if (edg->ti[j].goffs[0] != 0)
						eix = vx->ix + edg->ti[j].goffs[0];
					else
						eix = vx->ix + edg->ti[j].goffs[1];

					if ((nvx = get_vtxrec(vc, eix)) != NULL) {
						for (f = 0; f < fdi; f++)
							vx->vl[f] += nvx->vl[f] + scale * (vx->vv[f] - nvx->vv[f]);
						nn++;
					}
				}
				for (f = 0; f < fdi; f++)
					vx->vl[f] /= (double)nn;
			}
		}
	}

	DBG(("calc_surflin: creating rspl\n"));
//printf("calc_surflin: creating rspl\n");

	/* Now construct rspl setup mapping points from vertex normal output values */
	/* to adjusted vl values */

	if ((s->rev.surflin = new_rspl(RSPL_NOFLAGS, fdi, fdi)) == NULL)
		error("calc_surflin: new_rspl failed");

	nvtx++;		/* One for center point */

	/* Allocate rspl setup points */
	if ((mpoints = malloc(sizeof(cow) * 2 * nvtx)) == NULL)
//	if ((mpoints = malloc(sizeof(cow) * nvtx)) == NULL)
		error("calc_surflin: malloc of %d rspl setup points failed",nvtx);
	
	nvtx = 0;

#ifndef NEVER
	/* Center point */
	for (f = 0; f < fdi; f++) {
		mpoints[nvtx].p[f] = s->rev.ocent[f]; 
		mpoints[nvtx].v[f] = s->rev.ocent[f]; 
	}
	mpoints[nvtx].w = 10.0;
	nvtx++;
#endif

	/* Set the surface mapping points and restore vertexes values */
	for (i = 0; i < vc->hash_size; i++) {
		for (vx = vc->hash[i]; vx != NULL; vx = vx->hlink) {
			float *fcb;

			fcb = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */

			/* Actual output value as source of mapping */
			for (f = 0; f < fdi; f++)
				mpoints[nvtx].p[f] = vx->vv[f] = fcb[f];

			if (vx->status != vtx_norm) {	/* Skip isolated values */
				vx->status = vtx_norm;		/* Restore vtx contents */
				for (f = 0; f < fdi; f++)
					vx->vl[f] = vx->vv[f];
				continue;
			}

			for (f = 0; f < fdi; f++)
				mpoints[nvtx].v[f] = vx->vl[f];
			mpoints[nvtx].w = 1.0;
			nvtx++;

			vx->status = vtx_norm;		/* Restore vtx contents */
			for (f = 0; f < fdi; f++)
				vx->vl[f] = vx->vv[f];
			
//printf("Map[%d] %f %f %f -> %f %f %f\n", nvtx-1, mpoints[nvtx-1].p[0], mpoints[nvtx-1].p[1], mpoints[nvtx-1].p[2], mpoints[nvtx-1].v[0], mpoints[nvtx-1].v[1], mpoints[nvtx-1].v[2]);

#ifndef NEVER
			/* Add intermediate "fixed" point */
			for (f = 0; f < fdi; f++) {
				mpoints[nvtx].p[f] = mpoints[nvtx].v[f]
				                  = 0.5 * mpoints[nvtx-1].p[f] + 0.5 * s->rev.ocent[f];
			}
			mpoints[nvtx].w = 0.5;
			nvtx++;
#endif
		}
	}

	for (f = 0; f < fdi; f++) {
		min[f] = 1e38;
		max[f] = -1e38;
		vmin[f] = 1e38;
		vmax[f] = -1e38;
	}
	for (i = 0; i < nvtx; i++) {

#ifdef NEVER
		/* Blend with original values */
		for (f = 0; f < fdi; f++)
			mpoints[i].v[f] = 0.5 * mpoints[i].p[f] + 0.5 * mpoints[i].v[f];
#endif

		for (f = 0; f < fdi; f++) {
			if (mpoints[i].p[f] < min[f])
				min[f] = mpoints[i].p[f];
			if (mpoints[i].p[f] > max[f])
				max[f] = mpoints[i].p[f];

			if (mpoints[i].v[f] < vmin[f])
				vmin[f] = mpoints[i].v[f];
			if (mpoints[i].v[f] > vmax[f])
				vmax[f] = mpoints[i].v[f];
		}
	}

#ifdef REVVRML	/* Plot mapping vectors red->green */
	{
		vrml *wrl;
		double red[3] = { 1.0, 0.0, 0.0 };
		double green[3] = { 0.0, 1.0, 0.0 };

		wrl = new_vrml("suflinvecss", 0, vrml_lab);
		wrl->start_line_set(wrl, 0);

		for (i = 0; i < nvtx; i++) {
			wrl->add_col_vertex(wrl, 0, mpoints[i].p, red);
			wrl->add_col_vertex(wrl, 0, mpoints[i].v, green);
			
		}

		wrl->make_lines(wrl, 0, 2);
		wrl->del(wrl);
	}
#endif

	DBG(("calc_surflin: mapping points set, about to creat rspl:\n"));
//printf("calc_surflin: mapping points set, about to creat rspl:\n");

	/* Fit the rspl */
	s->rev.surflin->fit_rspl_w(s->rev.surflin, RSPL_NOFLAGS, mpoints, nvtx,  
	                           min, max, gres, vmin, vmax, 4.0, NULL, NULL);   

	DBG(("calc_surflin: mapping created\n"));
//printf("calc_surflin: mapping created\n");

#ifdef NEVER
	{
		double de;
		co p;
		extern double icmNorm33(double *, double *);

		/* Check fit */
		de = 0.0;

		for (i = 0; i < nvtx; i++) {
			for (f = 0; f < fdi; f++)
				p.p[f] = mpoints[i].p[f];

			s->rev.surflin->interp(s->rev.surflin, &p);
	
			de += icmNorm33(mpoints[i].v, p.v);
		}
		de = de/(double)nvtx;

		printf("Avg fit error = %f\n",de);
	}
#endif
	
	free(mpoints);

	/* Lookup ocent mapping offset */
	{
		co p;

		for (f = 0; f < fdi; f++)
			p.p[f] = s->rev.ocent[f];
		s->rev.surflin->interp(s->rev.surflin, &p);
//printf("opoint mapping %f %f %f -> %f %f %f\n", p.p[0], p.p[1], p.p[2], p.v[0], p.v[1], p.v[2]);
		for (f = 0; f < fdi; f++)
			s->rev.linoff[f] = p.v[f] - s->rev.ocent[f];
	}
	

#ifndef NEVER
	/* Put the transform into use */
	s->rev.surflin_en = 1;

	/* Transform all the vertexes */
	for (i = 0; i < vc->hash_size; i++) {
		for (vx = vc->hash[i]; vx != NULL; vx = vx->hlink) {
			logcomp(s, vx->vl, vx->vv, s->rev.ocent);

			/* Compute distance to overall center point squared */
			vx->dist = 0.0;
			for (f = 0; f < fdi; f++) {
				double tt = vx->vl[f] - s->rev.ocent[f]; 
				vx->dist += tt * tt;
			}
		}
	}
#endif

	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Initialise the rev Second section acceleration information. */
/* This is called when it is discovered on a call that s->rev.rev_valid == 0 */
static void init_revaccell(
rspl *s
) {
	int i, j;		/* Index of fwd grid point */
	int e, f, ee, ff;
	int di = s->di;
	int fdi = s->fdi;
	int gno = s->g.no;
	int rgno = s->rev.no;
	int rgres = s->rev.res;			/* number of bwd cells */
	int rgres_1 = rgres-1;			/* rgres -1 == maximum base coord value */

	schbase *b = s->rev.sb;		/* Base search information */
	char *vflag = NULL;			/* Per bwd vertex flag used during construction of nnrev */
								/* 0  nnrev[] cell empty,       not surface */
								/* 1  nnrev[] done/don't fill,  not surface */
								/* 2  nnrev[] cell empty,        on surface */
								/* 3  nnrev[] done,              on surface */
								/* 1X nnrev[] contains ink limited fwcells */
								/* Note that bit 1 can be set for cells that are not */
								/* to be explored because they are in the gamut interior, */
								/* and because they have already been added to the seedlist. */
	int pass;					/* Construction pass */
	float *gp;					/* Pointer to fwd grid points */

	DCOUNT(gg, MXRO, fdi, 0, 0, rgres);	/* Track the prime seed coordinate */
	int nn[MXRO];						/* bwd neighbor coordinate */

	vtxcache vc;				/* List + cache of vertexes being processed */
	tricache tc;				/* cache of surface triangles that have been processed */
	tricache stc;				/* small cache of surface triangles that have been processed */
	bxcell *bx, *nbx, **pbx;
	bxcell *xlist = NULL;		/* Linked list of added surface bxcells */
	assdire *tridir = NULL;		/* Triangle tables */
	assdire *edgdir = NULL;		/* Edge tables */
	double **cla = NULL;		/* Line LHS implicit equation matrix [fdi][fdi+1] */
	double *ta[MXRO], TA[MXRO][MXRO];	/* temp for intersection solving */

#if defined(REVTABLESTATS) || defined(DEBUG)
	/* Some statistics */
	unsigned long smsec;
	int nskcells = 0;					/* Number of skipped cells because over ink limit (debug) */
	int nascells = 0;					/* Number of added surface cells */
	int nrscells = 0;					/* Number of removed surface cells */
	int naoulvtxs = 0;					/* Number of added over ink limit vertexes */
	int revcells = 0;					/* Non-empty rev[] cells */
	int revcelldepth = 0;				/* Sum of rev[] list lengths */
	int ingamutcells = 0;				/* No of rev[] cells not on surface */
	int surfcells = 0;					/* No. surface cells */
	int emptycells = 0;					/* No. empty cells */
	int nnrevcells = 0;					/* Non-empty nnrev[] cells */
	int nnrevcellsearch = 0;			/* Sum of number of surface cells searched */
	int nnsinglefill = 0;				/* Number of nnrev[] cells seeded singly */
	int nnsuperfill = 0;				/* Number of nnrev[] cells seeded using supercell */
	int nnrevcelldepth = 0;				/* Sum of nnrev[] list lengths */
	int nnmxrevcelldepth = 0;			/* Maximum nnrev[] list lengths */
	int nnrevshare = 0;					/* Sum of nnrev[] list reference counts */
#endif
	datao rgmin, rgmax;

	DBG(("init_revaccell called, di = %d, fdi = %d, mgres = %d\n",di,fdi,(int)s->g.mres));

	/* To help VRML diagnostics, make a guess as to whether the output */
	/* space is XYZ like, or L*a*b* like */

	s->get_out_range(s, rgmin, rgmax);	/* overall output min/max */

	if (fdi >= 3
	 && rgmin[0] >= -1.0 && rgmax[0] < 3.0
	 && rgmin[1] >= -1.0 && rgmax[1] < 3.0
	 && rgmin[2] >= -1.0 && rgmax[2] < 3.0) {
		s->rev.probxyz = 1;
		if (s->verbose)
			fprintf(stdout, "%cLooks like an XYZ space\n",cr_char);
	}

	if (fdi > 1 && s->verbose)
		fprintf(stdout, "%cInitializing nnrev arrays...\n",cr_char);

	/* Add this instance into memory management */
	if (s->rev.rev_valid == 0 && di > 1) {
		rev_struct *rsi;
		size_t ram_portion = g_avail_ram;

		/* Add into linked list */
		s->rev.next = g_rev_instances;
		g_rev_instances = &s->rev;

		/* Aportion the memory, and reduce cache if it is over new limit. */
		g_no_rev_cache_instances++;
		ram_portion /= g_no_rev_cache_instances; 
		for (rsi = g_rev_instances; rsi != NULL; rsi = rsi->next) {
			revcache *rc = rsi->cache;

			rsi->max_sz = ram_portion;
			while (rc->nunlocked > 0 && rsi->sz > rsi->max_sz) {
				if (decrease_revcache(rc) == 0)
					break;
			}
//printf("~1 rev instance ram = %d MB\n",rsi->sz/1000000);
		}
		
		if (s->verbose)
			fprintf(stdout, "%cThere %s %d rev cache instance%s with %lu Mbytes limit\n",
			                    cr_char,
								g_no_rev_cache_instances > 1 ? "are" : "is",
			                    g_no_rev_cache_instances,
								g_no_rev_cache_instances > 1 ? "s" : "",
			                    (unsigned long)ram_portion/1000000);
	}

#if defined(REVTABLESTATS) || defined(DEBUG)
	smsec = msec_time();
#endif

	/* Temporary per bwd vertex/cell flag for nn setup */
	if ((vflag = (char *) rev_calloc(s, rgno, sizeof(char))) == NULL)
		error("rspl malloc failed - rev.vflag points");
	INCSZ(s, rgno * sizeof(char));

	/*
	 * The rev[] and nnrev[] grids contain pointers to lists of grid cube base indexes.
	 * If the pointer is NULL, then there are no base indexes in that list.
	 * A non NULL list uses element [0] to indicate the allocation size of the list,
	 * [1] contains the index of the next free location, [2] contains the reference
     * count (lists may be shared), the list starts at [3]. The last entry is marked with -1.
	 */

	/* We won't include any fwd cells that are over the ink limit, */
	/* so makes sure that the fwd cell nodes all have an ink limit value. */ 
	if (b != NULL && s->limiten) {
		ECOUNT(gc, MXDIDO, s->di, 0, s->g.res, 0);    /* coordinates */
		double iv[MXDI];				/* Input value corresponding to grid */

		DBG(("Looking up fwd vertex ink limit values\n"));
//printf("Looking up fwd vertex ink limit values\n");
//printf("s->limitv = %f\n",s->limitv);
		/* Calling the limit function for each fwd vertex could be bad */
		/* if the limit function is slow. Maybe an octree type algorithm */
		/* could be used if this is a problem ? */
		EC_INIT(gc);
		for (i = 0, gp = s->g.a; i < s->g.no; i++, gp += s->g.pss) {
			if (gp[-1] == L_UNINIT) {
				for (e = 0; e < di; e++)
					iv[e] = s->g.l[e] + gc[e] * s->g.w[e];  /* Input sample values */
				gp[-1] = (float)(INKSCALE * s->limitf(s->lcntx, iv));
//printf("~1 set ix %d limitv to %f\n",i,gp[-1]);
			}
//else printf("~1 ix %d limitv is %f\n",i,gp[-1]);
			EC_INC(gc);
		}
		s->g.limitv_cached = 1;
	}

	/* We then fill in the in-gamut reverse grid lookups, */
	/* and identify nnrev prime seed verticies to put in the surface bxcells. */

	DBG(("filling in rev.rev[] grid\n"));
	
	/* To create rev.rev[], for all fwd grid points, form the cube with that */
	/* point at its base, and determine the bounding box of the output values */
	/* that could intersect that fwd cube. Add that fwd index to the lists of */
	/* of all bwd cells that the bounding box intersects. */
	/* As a start for creating surface bxcell list, flag which bwd verticies */
	/* are covered by the fwd grid output range. */

	/* Pre-marking device edge rev cells creates many more initial cells, */
	/* but avoids having to discover them with multiple passes ? */
	for (gp = s->g.a, i = 0; i < gno; gp += s->g.pss, i++) {
		datao min, max;
		int imin[MXRO], imax[MXRO], gc[MXRO];
		int edge = 0;		/* This fwd cell contains a device edge */
		int uil;			/* One is under the ink limit */
		int oil;			/* One is over the ink limit */

//printf("~1 i = %d/%d\n",i,gno);
		/* Skip grid points on the upper edge of the grid, since there */
		/* is no further grid point to form a cube range with. */
		for (e = 0; e < di; e++) {
			int flags = G_FL(gp, e);

			if (flags == 0)		/* At the top edge */
				break;

			/* If we at the bottom edge, or one away from top edge */
			if (flags == 4 || flags == 1)
				edge = 1;				/* This fwd cell is on device gamut edge */
		}
		if (e < di) {	/* Top edge - skip this cube */
//printf("~1 skipping base vertex %d on top edge\n",i);
			continue;
		}

//printf("~1 adding to rev[]\n");

		/* Find the output value bounding box values for this grid cell */
		/* Start with base vertex */
		uil = oil = 0;
		for (f = 0; f < fdi; f++)	/* Init output min/max */
			min[f] = max[f] = gp[f];

		if (!s->limiten || gp[-1] <= s->limitv)
			uil = 1;
		else
			edge = oil = 1;		/* May be stradling ink limit edge */
	
		/* Then add all other fwd cube verticies */
		for (ee = 1; ee < (1 << di); ee++) {
			float *gt = gp + s->g.fhi[ee];	/* Pointer to cube vertex */
			
			if (!s->limiten || gt[-1] <= s->limitv)
				uil = 1;
			else
				edge = oil = 1;

			/* Update bounding box for this grid point */
			for (f = 0; f < fdi; f++) {
				if (min[f] > gt[f])	
					 min[f] = gt[f];
				if (max[f] < gt[f])
					 max[f] = gt[f];
			}
		}

		/* Skip any fwd cells that have every vertex over the ink limit */
		if (!uil) {
#if defined(REVTABLESTATS) || defined(DEBUG)
			nskcells++;
#endif
			continue;
		}

		/* Figure out intersection range in bwd cell grid */
		for (f = 0; f < fdi; f++) {
			double t;
			int mi;
			double gw = s->rev.gw[f];
			double gl = s->rev.gl[f];
			t = (min[f] - gl - EPS)/gw;
			mi = (int)floor(t);			/* Grid coordinate */
			if (mi < 0)					/* Limit to valid cube base index range */
				mi = 0;
			else if (mi > rgres_1)
				mi = rgres_1;
			imin[f] = mi;	
			t = (max[f] - gl + EPS)/gw;
			mi = (int)floor(t);			/* Grid coordinate */
			if (mi < 0)					/* Limit to valid cube base index range */
				mi = 0;
			else if (mi > rgres_1)
				mi = rgres_1;
			imax[f] = mi;	
		}

//printf(" Scanning over bwd cell range grid:\n");
//for (f = 0; f < fdi; f++)
//printf(" Min[%d] = %d -> Max[%d] = %d\n",f,imin[f],f,imax[f]);

		/* Now create forward index and vector with all the reverse grid cells */
		for (f = 0; f < fdi; f++)
			gc[f] = imin[f];	/* init coords */

		/* Until increment at bottom carries */
		for (f = 0; f < fdi;) {	/* For all of intersect bwd cube */
			int **rpp;
			char *vflagp;
			
			/* Compute pointer to bwd grid cell and vflag[] */
			for (rpp = s->rev.rev, vflagp = vflag, f = 0; f < fdi; f++) {
				int inc = gc[f] * s->rev.coi[f];
				rpp    += inc;
				vflagp += inc;
			}

#undef PRE_LOAD_SURFACE		/* [und] Makes it slower ? */
#ifdef PRE_LOAD_SURFACE				/* Pre-load device edge cells */
			if (edge) {
				*vflagp = 2;		/* This is definitely a gamut surface bwd cell */
									/* and so nnrev[] needs to be filled */
			} else
#endif
			if (*vflagp == 0) {
				*vflagp = 1;		/* This is possibly not a surface bwd cell, */
									/* and otherwise is an inside gamut bwd cell */
			}

			if (oil)
				*vflagp |= 0x10;	/* Contains over ink limit vertexes */

//printf("seting vflag[%d] to surface done (%x)\n",vflagp-vflag,*vflagp);

//printf("  Currently at grid ix %d, (vflag = %x) adding fwd %d:\n",vflagp-vflag,*vflagp,i);
//for (f = 0; f < fdi; f++)
//printf("  gc[%d] = %d\n",f,gc[f]);

			/* Add fwd cells to rev[] list */
			add2indexlist(s, rpp, i, 0);

			/* Increment index up to and including imax[] */
			for (f = 0; f < fdi; f++) {
				gc[f]++;
				if (gc[f] <= imax[f])
					break;	/* No carry */
				gc[f] = imin[f];
			}
		}	/* Next reverse grid point in intersecting cube */
	}	/* Next base grid point */

	DBG(("We skipped %d/%d cells that were over the limit\n",nskcells,gno));

#ifdef CHECK_NNLU
	if (fdi > 1) {
		/* Check that every flagged rev[] cell is filled */
		printf("Checking all %d flagged rev[] cells are filled\n",rgno);
		for (i = 0; i < rgno; i++) {
			if (   (vflag[i] & 1) != 0 
			     && (s->rev.rev[i] == NULL || s->rev.rev[i][1] == 3)) {
				printf("Found empty rev[%d] ?:\n",i);
				printf(" vflag %x\n",vflag[i]);
				if (s->rev.rev[i] == NULL)
					printf(" rev = NULL\n");
				else
					printf(" rev = length = %d\n",s->rev.rev[i][1]-3);
			}
		} 
	}
#endif	/* CHECK_NNLU */

	/* If doing fast setup, then this is all we need. */
	if (s->rev.fastsetup) {

		/* Free up flag array used for construction */
		if (vflag != NULL) {
			DECSZ(s, rgno * sizeof(char));
			free(vflag);
		}

		s->rev.rev_valid = 1;

		if (fdi > 1 && s->verbose)
			fprintf(stdout, "%cFast nnrev initialization done\n",cr_char);

		DBG(("init_revaccell fastsetup finished\n"));

#if defined(REVTABLESTATS) || defined(DEBUG)
		printf("Fastsetup took %f seconds\n",0.001 * (msec_time()-smsec));
#endif

		return;
	}

	/* Rough outline of overall nn setup process:

		Fill rev[] array by scanning fwd cells.

		Locating initial surface bwd cells.

		loop:
			Fill empty surface cells from rev[] list and convert to vertexes.

			(In two phases, first just against primary bx, second with all 
			 shadowed bx's:)

			Test all triangles against all vertexes and mark those that are shadowed.

			Remove vertexes from bx if they have been deleted, but leave
			them in the vertex cache for testing against.

			If any vertexes of a bx land outside it in a bx that is not
			part of the surface list, add that bx to the surface list and
			mark it for processing.

			Track which bx cells shadow newly added bx cells, 
			so that new bx cells get tested against all their shadowers,
			as well as being used to test against their shadowees.

		Locate and preserve all overlapping surface triangles.

		Delete any shadowed vertexes, and remove any empty bxcells.

		Add extra over ink limit vertexes.

		Convert vertexes back to minimum number of fwd cubes.

	*/

	calc_ocent(s);

	/* Locate and process the surface bxcells and fill the nnrev array if we */
	/* are not doing a fast setup. (fastsetup will instead fill the nnrev[] array */
	/* on demand, by searching the rev[] array.) */
	DBG(("Identifying surface rev cells\n"));

	/* Allocate the surflist hash index. */
	/* (Note that we track bxcells in the surface list rather */
	/* than the hash list, in this context.) */
	create_surfhash(s);

	/* Locate surface reverse cells */
	DC_INIT(gg);
	for (i = 0; i < rgno; i++) {

		if ((vflag[i] & 0xf) == 1) {	/* if filled rev[] cell but not surface */ 
			char *vflagp;

			/* Check face neighbors */
			int cc[MXDO];				/* Neigbor offset counter */

			/* Check if any of the face neigbors of this bwd cell are empty. */
			/* If so, mark it as a surface cell. */
			/* [This won't detect all surface nncells, but will hit most of them */
			/* without including too many false ones. The vertex filter code */
			/* should discover any surface nncells that are missed.] */
			for (f = 0; f < fdi; f++)
				cc[f] = gg[f];
			vflagp = vflag + i;

			for (ff = 0; ff < (fdi << 1); ff++) {
				f = ff >> 1;

				cc[f]  += (ff & 1) ? 1 : -1;
				vflagp += (ff & 1) ? s->rev.coi[f] : -s->rev.coi[f];

				/* Out of bounds or empty */
				if (cc[f] < 0 || cc[f] >= rgres || ((*vflagp & 0xf) == 0)) {
					vflag[i] = (vflag[i] & ~0xf) | 2;	/* Convert this one to empty surface cell */
//printf("seting vflag[%d] to surface cell (%x)\n",i,vflag[i]);

					/* Add a bxcell to surf hash. Initial status = bx_uninit */
					if ((bx = get_surface_bxcell(s, i)) == NULL) {
						/* Since it's a surface point, the seeding point is itself (NULL). */
						bx = new_bxcell(s, i, gg, NULL, 0.0, NULL);
						add_bxcell_hash(s, bx);
		
						/* Add to surface linked list */
						bx->slist = s->rev.surflist; 
						s->rev.surflist = bx;
//printf("~1 adding nnrev[%d] to surface list\n",bx->ix);
					}
					break;		
				}

				cc[f]  -= (ff & 1) ? 1 : -1;
				vflagp -= (ff & 1) ? s->rev.coi[f] : -s->rev.coi[f];
			}
		} 
#ifdef PRE_LOAD_SURFACE
		  else if ((vflag[i] & 0xf) == 2) {		/* Pre-marked surface rev cell */

			/* Add a bxcell to surf hash. Initial status = bx_uninit */
			if ((bx = get_surface_bxcell(s, i)) == NULL) {
				/* Since it's a surface point, the seeding point is itself (NULL). */
				bx = new_bxcell(s, i, gg, NULL, 0.0, NULL);
				add_bxcell_hash(s, bx);

				/* Add to surface linked list */
				bx->slist = s->rev.surflist; 
				s->rev.surflist = bx;
//printf("~1 adding pre-marked nnrev[%d] to surface list\n",bx->ix);
			}
		}
#endif /* PRE_LOAD_SURFACE */

#if defined(REVTABLESTATS) || defined(DEBUG)
		if (vflag[i] & 2) 
			surfcells++;
		else if ((vflag[i] & 0xf) != 0)
			ingamutcells++;
		else 
			emptycells++;

		if (s->rev.rev[i] != NULL) {
			revcells++;
			revcelldepth += s->rev.rev[i][1]-3;
		}
#endif
		DC_INC(gg);
	}

	if (di < 2)
	{
		/* Create surface fwd cell list */
		DBG(("create surface fwd cell lists\n"));

		/* For each rev[] containing fwd cells, */
		/* copy the cells to the corresponding surface bxcel cell */
		for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
			int *crp, *rp;
	
			if ((crp = s->rev.rev[bx->ix]) == NULL)
				error("Surface list bxcell ix %d has no vertexes",bx->ix);

			/* For each fwd cell in surface rev[] */ 
			for (rp = crp+3; *rp != -1; rp++) {
				add2indexlist(s, &bx->sl, *rp, 0);
			}
		}

	} else {		/* di >= 3 gamut surface finding */

#ifdef REVVRML
		/* Plot the initial surface bxcells & their fwd cells. */
		/* Rev cells? Fwd cells? Fwd cell base indexs? */
		if (0) plot_bxfwcells(s, 0, 1, 0);
#endif /* REVVRML */

		/* per reverse cell vertex cache */
		create_vtxrec_list(s, &vc);

		/* per reverse cell surface triangle cache */
		create_trirec(s, &tc, 0);

		/* small per reverse cell surface triangle cache */
		create_trirec(s, &stc, 1);
	
		/* create associated sub-simplex (triangle) lookup table */
		init_assdir(s, &tridir, 2);

		/* create associated sub-simplex (edge) lookup table */
		init_assdir(s, &edgdir, 1);
	
		/* Process surface bxcells */
		/* (Maintain current list of vtxrec's for all vertexes) */
			
		/* Setup temporary matrix */
		for (f = 0; f < 2; f++)
			ta[f] = TA[f];

		/* - - - - - - - - - - - - - - */
		/* fill, thin and add, until */
		/* there is no more work to do. */
		for (pass = 0;; pass++) {
			int phase;
			int morevtxadded = 0;
#if defined(REVTABLESTATS) || defined(DEBUG)
			unsigned long lmsec = msec_time();
			int thcount = 0, rethcount = 0;

//			printf("At top of gamut surface loop\n");
#endif

			/* For each surface bxcell, convert the corresponding */
			/* rev[] fwd cubes into vertices. */
			/* (Must keep bxcells even if none of their verticies */
			/*  are physically in them, so that those verticies get thinned. */
			/*  could only remove them if vertex was not in any surface cell ?) */
			for (pbx = &s->rev.surflist, bx = *pbx; bx != NULL; bx = nbx) {
				int *crp, *rp, *nrp;
				vtxrec *vx;

				nbx = bx->slist;
		
				if (bx->status != bx_uninit) {
					pbx = &bx->slist;
					continue;
				}

				if ((crp = s->rev.rev[bx->ix]) == NULL)
					error("Surface list bxcell ix %d has no vertexes",bx->ix);
		
//printf("Initializing bxcell %d with vertexes\n",bx->ix);
				/* For each fwd cell in surface rev[] */ 
				for (rp = crp+3; *rp != -1; rp++) {
		
//adding cube %d to bx %d\n",*rp, bx->ix);
					/* For each vertex of cube */
					for (ee = 0; ee < (1<<di); ee++) {
						int vix = *rp + s->g.hi[ee];
						float *fcb = s->g.a + vix * s->g.pss;	/* Pointer to base float of fwd cell */
						vtxrec *vx;
		
//printf("~1 adding cube %d vtx %d to bx %d\n",*rp, vix, bx->ix);

						/* Don't add over ink limit vertexes */
						/* (we'll re-add them in later) */
						if (s->limiten && fcb[-1] > s->limitv) {
//printf("Skipping vtx %d because over ink limit\n",vix);
							continue;
						}
		
						if ((vx = get_vtxrec(&vc, vix)) != NULL) {
							if (vx->rix == bx->ix) {
//printf("Already have vertex %d in bx %d\n",vx->ix,vx->rix);
							}

							/* Skip vertexes that we've already added to this bxcell */
							if (vx->tix == bx->ix) {
//printf("Skipping vtx %d because alread in bx %d\n",vix,bx->ix);
								continue;
							}
						} else {
							/* Create new vertex */
							vx = new_vtxrec(s, &vc, vix); 
							vx->tix = bx->ix;		/* Added to this bx */
//printf("Create vtx %d for bx %d (actually in bx %d)\n",vix,bx->ix,vx->rix);
						}
		
						/* Add vertex to bxcell sl list */
						add2indexlist(s, &bx->sl, vix, 0);

						if (vx->rix == bx->ix) {
//printf("Added vertex %d is in this bx %d\n",vx->ix,vx->rix);
						} else {
//printf("Added vertex %d is in different bx %d to this one %d\n",vx->ix,vx->rix,bx->ix);
						}
					}
				}
				/* Expand a bxcell's shadow testing group values based on it's vertex list */
				/* so that shadow testing works correctly for vertexes that don't */
				/* actually lie within the bxcell. (Note that in fact the triangle */
				/* testing creates triangles that are made of vertexes that may not */
				/* be in this bx's list, so the shadow size doesn't accuratly  reprsent */
				/* the possible shadow area. It's not clear what consequences this has, */
				/* if any. If we extanded the group to cover this, we would need to have */
				/* two groups, a shadower group including those vertexes, and a shadowee */
				/* goup for just those vertexes that are part of the bx. */
				extend_bxcell_shadow_group(s, &vc, bx);
				bx->status = bx_filled;
				pbx = &bx->slist;
				morevtxadded = 1;
			}
	
			/* Compute transform rspl that helps "unfold" any regions of the surface */
			/* that overlap from the perspective of ocent, to try and avoid gaps in */
			/* the final gamut surface. Existing vtxrec are converted to have vl */
			/* in the unfolded space. */
			if (pass == 0 /* && function flag set */) {
#ifndef EN_UNTWIST		/* Control using an environment variable */
				if (getenv("ARGYLL_UNTWIST_GAMUT_SURFACE") != NULL)
#endif
				{
					calc_surflin(s, &vc, edgdir);
				}
			}

			DBG(("thinning surface vertex lists and converting to cells\n"));

			/* (Sorting bxcells doesn't seem to make any performace difference.) */

			for (phase = 0; phase < 2; phase++) {

//printf("Phase %d\n",phase);

				/* For each surface bxcell, form triangles from vertexes */
				/* and mark as shadowed and other vertexes that are in the */
				/* triangles shadow. */
				/* rev[] fwd cubes into vertices. */
				for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
					int sdi = 2;						/* sub-simplexes are triangles */
					double clb[MXRO+1];	/* Line RHS implicit equation vector [fdi+1] */
					int *crp, *rp, *nrp;
					vtxrec *vx, *nvx;
					int aftercount;			/* vertex count after thinning */

					if (bx->status != bx_filled && bx->status != bx_rethinnd) {
//printf("~1 skipping bx %d because status = %d\n",bx->ix,bx->status);
						continue;
					}
//printf("~1 checking bx %d\n",bx->ix);

					/* Only do first pass through primary alone if never thinned before */
					if (phase == 0 && bx->status == bx_rethinnd) {
						continue;
					}

					/* If this bxcell is empty (because all it's vertexes are shadowed ?) */
					if (bx->sl == NULL || bx->sl[1] == 3) {
//printf("~1 skipping nnrev[%d] because it's empty\n",bx->ix);
						continue;
					}
//printf("Thinning bxcell %d\n",bx->ix);
					/* Create nnrev[] shadowing linked list. nnrev[] cells who's shadow in */
					/* the direction of rev.ocent[] touches another nnrev[], add that nnrev[] */
					/* to their shadow list. This allows us to filter vertexes in other */
					/* nnrev[] cells from triangles above them */
					bx->wlist = NULL;

					/* Only go through all shadowed bxcells once primary has been */
					/* thinned alone */
					if (phase == 1) {

						/* Use just extra list for re-thinning, for 10% speed advantage. */
						if (bx->status == bx_rethinnd && xlist != NULL) {
//printf("Adding shadows to bxcell %d from xlist\n",bx->ix);
							for (nbx = xlist; nbx != NULL; nbx = nbx->xlist) {

								if (nbx->status == bx_uninit)	/* Newly added cells (shouldn't happen) */
									break;

								if (nbx == bx)
									continue;
							
								/* If any of bx is further from nbx and their bounding */
								/* cylinders overlap in perspective from rev.ocenter, */
								/* assume nbx is a shadow */
								if (shadow_group_group(s, s->rev.ocent, bx->g.bcent, bx->cc,
									                 bx->dw, nbx->g.bcent, nbx->cc, nbx->dw)) {
									nbx->wlist = bx->wlist;
									bx->wlist = nbx;
//printf("~1 adding shadow nnrev[%d] from xlist\n",nbx->ix);
								}
							}
						} else {

//printf("Adding shadows to bxcell %d from surflist\n",bx->ix);
							for (nbx = s->rev.surflist; nbx != NULL; nbx = nbx->slist) {

//printf("Considering bx %d for shadow list\n",nbx->ix);
								if (nbx->status == bx_uninit) /* Newly added cells (shouldn't happen) */
									break;

								if (nbx == bx) 
									continue;
							
								/* If any of bx is further from nbx and their bounding */
								/* cylinders overlap in perspective from rev.ocenter, */
								/* assume nbx is a shadow */
								if (shadow_group_group(s, s->rev.ocent, bx->g.bcent, bx->cc, bx->dw,    
									                                    nbx->g.bcent, nbx->cc, nbx->dw))
								{
//printf("Added bx %d for shadow list, prim bx %d\n",nbx->ix,bx->ix);
									nbx->wlist = bx->wlist;
									bx->wlist = nbx;
								}
							}
						}
					}	/* if phase == 1 */

#if defined(REVTABLESTATS) || defined(DEBUG)
					if (bx->status == bx_rethinnd)
						rethcount++;
					else
						thcount++;
#endif

					/* Abort doing this cell until all its shadowees are filled */
					/* (Shouldn't happen ?) */
					if (nbx != NULL) {
//printf("Skipping thinning of bx %d because newly added bx %d is in surfce list\n",bx->ix,nbx->ix);
						continue;
					}
			
					/* Be able to detect triangles already tested */
					/* from this shadowing bxcell. */
					clear_trirec(s, &tc);

					/* Put just primary and shadows on vx->tlist */
					vc.vtxlist = NULL;
					vc.nilist = 0;

					/* Add all the secondary bxcell vertexes to the vtxlist */
					for (nbx = bx->wlist; nbx != NULL; nbx = nbx->wlist) {
//printf("Adding bx %d verticies\n",nbx->ix);
						for (rp = nbx->sl+3; *rp != -1; rp++) {

							if ((vx = get_vtxrec(&vc, *rp)) == NULL)
								error("Failed to find vertex %s in cache",*rp);

//printf("Checking ix %d from bx %d\n",vx->ix,nbx->ix);
							/* Check vertex falls within shadow of main bx */
							/* (just checking non-deleted vertexes (triangles) */
							/*  improves speed by 20%, but we end up with stray fwd cells */
							/*  and some holes, because crossed triangles vertexes get */
							/* marked deleted ??) */
							if (
//						    vx->status == vtx_norm &&
							    shadow_group_vertex(s, s->rev.ocent, bx->g.bcent, bx->cc,
								                                          bx->dw, vx->vl)) { 
									add_vtxrec_list(&vc, vx, 0);	/* Add if not deleted */
//printf(" Added ix %d from bx %d\n",vx->ix,nbx->ix);
							}
//else 
//printf(" Not added ix %d from bx %d because no within prim bx %d\n",vx->ix,nbx->ix,bx->ix);
						}
					}

					/* Add all the primary bxcell verticies to the list, and */
					/* mark them (override shadow mark) */
//printf("Adding bx %d verticies\n",bx->ix);
					for (rp = bx->sl+3; *rp != -1; rp++) {
						if ((vx = get_vtxrec(&vc, *rp)) == NULL)
							error("Failed to find vertex %s in cache",*rp);

						if (vx->status == vtx_norm &&
						    shadow_group_vertex(s, s->rev.ocent, bx->g.bcent, bx->cc,
							                                          bx->dw, vx->vl)) { 
								add_vtxrec_list(&vc, vx, 1);	/* Add if not hidden/deleted */
						}
					}

					aftercount = vc.nilist;

					/* sort vertexes by decending distance to center point */
					/* (and also reset list tflag) */
					sort_vtxrec_list(s, &vc);
		
					/* For vertexes of this bxcell and shadowers, */
					/* in order from largst to smallest distance from center. */
					for (vx = vc.vtxlist; vx != NULL; vx = vx->tlist) {
						float *fcb;				/* Vertex being tested */
						int fl;
						assdire *tri;			/* Triangle table */

//printf("~1 checking against vtx %d\n",vx->ix);

						/* Only check triangles using verticies of the primary bxcell, */
						/* not shadow bx's. */
						if (!vx->prim)
							continue;

//printf("~1 doing vertex %d at %s dist %f\n",vx->ix, debPdv(fdi,vx->v), sqrt(vx->dist));

						fcb = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */
						fl = FLV(fcb);		/* Edge flags for this vertex */

						tri = tridir + fl;
//printf("~1 fl %d = 0o%o, no triangles %d\n",fl, fl, tri->no);

						/* For all possible triangles that use this vertex */
						for (i = 0; i < tridir[fl].no; i++) {
							int triix[3];
							vtxrec *trivx[3];
							double v[MXRI+1][MXRO]; 	/* Triangle vertex values */
							double gc[MXRO], cc, dw;	/* Triangle shadow group info. */
							int ntvsh = 0;				/* Number of triangle verticies shadowed */
							double bdist = -1.0;

							/* Get triangle verticy values */
							for (e = 0; e <= sdi; e++) {
								triix[e] = vx->ix + tri->ti[i].goffs[e];
								
								if ((trivx[e] = get_vtxrec(&vc, triix[e])) == NULL)
									break;		/* Vertex doesn't exist in our set */

							    if (trivx[e]->status != vtx_norm)
									ntvsh++;

								if (trivx[e]->dist > bdist)
									bdist = trivx[e]->dist;
							}
//printf("~1 tri %d: vtxs %s goffs %s\n",i, debPiv(di,triix), debPiv(sdi+1, tri->ti[i].goffs));

							/* Don't test against triangle unless all vertexes */
							/* are in current surface, and whole triangle is visible. */
							if (e <= sdi || ntvsh >= 3)
								continue;

							/* If triangle has been done before for this bxcell, skip it. */
							if (check_trirec(s, &tc, triix))
								continue;

							for (e = 0; e <= sdi; e++) {
								for (f = 0; f < fdi; f++)
									v[e][f] = trivx[e]->vl[f];
							}

							/* Compute shadow group params of triangle for quick vertex test */
							comp_shadow_group(s, s->rev.ocent, gc, &cc, &dw, NULL, v, sdi+1);    

							/* For all vertexes */
							for (nvx = vc.vtxlist; nvx != NULL; nvx = nvx->tlist) {
								double pv[MXRO];	/* Vertex being tested */
								double de[MXRO];	/* Line delta */
								double tb[MXRI];	/* Solution point in input space */
								double xv[MXRO];	/* Solution point in output space */
								int g, sorv, wsrv;	/* Solved & within simplex return value */
								double dist;		/* distance to line origin */
								double dot;			/* dot product of solution to line */
								int shdwd;			/* whether vertex is shadowed */

								/* If vertex is above triangle, it can't be shadowed */
								if (nvx->dist > bdist)
									continue;

								/* If this other vertex has already been deleted, skip it */
								if (nvx->status != vtx_norm)
									continue;

								/* If this other vertex is part of the triangle, skip it */
								if (nvx->ix == triix[0]
								 || nvx->ix == triix[1]
								 || nvx->ix == triix[2]) {
									continue;
								}

//printf("~1 checking vertex %d against tri %s\n",nvx->ix,debPiv(3,triix));

								/* Do quick check against triangle */
								if (!shadow_group_vertex(s,
									     s->rev.ocent, gc, cc, dw, nvx->vl)) { 
//printf("~1 shadow group check shows no intersection\n");
									continue;
								}

//printf("~1 checking vertex %d at %s dist %f\n",nvx->ix, debPdv(fdi,nvx->v), sqrt(nvx->dist));

								/* Compute intersection: */
								shdwd = wsrv = 0;

								/* Compute line delta */
								fcb = s->g.a + nvx->ix * s->g.pss;
								for (f = 0; f < fdi; f++)
									pv[f] = fcb[f];

								logcomp(s, pv, pv, s->rev.ocent);

								for (f = 0; f < fdi; f++)
									de[f] = pv[f] - s->rev.ocent[f];

								/* Setup line cla and clb */
								init_line_eq_imp(s, NULL, &cla, clb, s->rev.ocent, de, 0);

								/* Solve line/triangle intersection using same */
								/* method as vnearest_clip_solve(). */
							
								/* LHS: ta[sdi][sdi] = cla[sdi][fdi] * vv[fdi][sdi] */
								/* RHS: tb[sdi] = clb[sdi] - cla[sdi][fdi] * vv_di[fdi] */
								for (f = 0; f < sdi; f++) {
									double tt;
									for (e = 0; e < sdi; e++) {
										for (tt = 0.0, g = 0; g < fdi; g++)
											tt += cla[f][g] * (v[e][g] - v[e+1][g]);
										ta[f][e] = tt;
									}
									for (tt = 0.0, g = 0; g < fdi; g++)
										tt += cla[f][g] * v[sdi][g];
									tb[f] = clb[f] - tt;
								}
							
								/* Compute the solution */
								/* (Solve the simultaneous linear equations A.X = B) */
//								sorv = !solve_se(ta, tb, sdi);
								sorv = !solve_se_2x2(ta, tb);		/* Saves a few % only */

								/* If it was solved */
								if (sorv) {
								
									/* Check that the solution is within the simplex & ink limit */
									if ((wsrv = simple_within_simplex(v, tb, sdi)) != 0) {

										/* Compute the output space solution point */
										for (f = 0; f < fdi; f++) {
											double tt = 0.0;
											for (e = 0; e < sdi; e++)
												tt += (v[e][f] - v[e+1][f]) * tb[e];
											xv[f] = tt + v[sdi][f];
										}
								
										/* Compute distance to gamut center squared, */
										/* as well as the dot product */
										for (dot = dist = 0.0, f = 0; f < fdi ; f++) {
											double tt = (xv[f] - s->rev.ocent[f]);
											dist += tt * tt;
											dot += de[f] * tt;
										}
//printf("~1 intersection at %s dist %f\n", debPdv(fdi,xv), sqrt(dist));
   
										/* If intersection distance is greater than vertex distance, */
										/* delete the vertex */
										if (dot > 0.0 && dist > (nvx->dist + EPS)) {
											shdwd = 1;
											nvx->status = vtx_sha;		/* Shadowed */
											aftercount--;
//printf("~1 deleting vx %d\n",nvx->ix);
										}
									}
								}
//if (!sorv) printf("~1 solve failed\n");
//if (sorv && !wsrv) printf("~1 %d not within simplex, tb = %s\n",nvx->ix, debPdv(sdi,tb));
//if (sorv && wsrv && shdwd) printf("~1 tri %s deleting vertex %d\n",debPiv(3,triix), nvx->ix);

#ifdef REVVRML
							/* Plot vertex & triangle check setup & solution */
							/* + the primary and shadow bxcells. */ 
							/* Plot prim & shadow bxcell cells ? Wait for user press ? */
							if (0 && phase && shdwd) plot_tri_check(s, 1, 1,
								bx, vx->ix, i, triix, nvx->ix, sorv, wsrv, shdwd, v, de, pv, xv); 
#endif /* REVVRML */

							}	/* Next other vertex */
						}		/* Next triangle */
					}			/* Next main vertex */

					if (phase == 1)
						bx->status = bx_thinned;
//printf("Thinned vertexes in bx %d from %d to %d (%d)\n",bx->ix, vc.nilist,aftercount, vc.nilist-aftercount);
				}			/* Next surface bx cell */
				/* Done with lists */
				vc.vtxlist = NULL;
				vc.nilist = 0;
				xlist = NULL;

				DBG(("deleting verticies in all bxcells\n"));

				/* The thinning may have deleted verticies from bxcell's that */
				/* were not involved in the thinning, so go though all bxcells */
				/* to do deletions. Look also for any needed additional surface bxcells. */
				for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
					int beforecount, aftercount;
					vtxrec *nvx;
					int *crp, *rp, *nrp;
			
					if (bx->status == bx_uninit)
						continue;

#ifdef REVVRML
					bx->debug = 0;		/* Not an addition */
#endif

					beforecount = bx->sl[1]-3;

#undef DELETE_SHAD	/* [und] try deleting shadowed vertexes with no un-shadowed neighbors. */
					/* Seems to actually slow things down though ? */

					/* Delete all the marked vertexes from bxcell list */
					for (nrp = rp = bx->sl+3; *rp != -1; rp++) {
						vtxrec *vx;
#ifdef DELETE_SHAD
						int nshad = 0, nnshad = 0;
#endif

						if ((vx = get_vtxrec(&vc, *rp)) == NULL)
							continue;			/* Already deleted */

#ifdef REVVRML
						vx->addvtx = 0;
#endif

#ifdef DELETE_SHAD
						/* Check all of its neighbor vertexes, to see if */
						/* it's safe to actually delete them. */
						if (vx->status >= vtx_sha) {	/* vertex to delete ? */
							float *fcb;
							int fl;
							assdire *edg;			/* Edge table */

//printf("Checking vx %d neighbors\n",vx->ix);
							fcb = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */
							fl = FLV(fcb);		/* Edge flags for this vertex */
							edg = edgdir + fl;

							/* For all possible edges that use this vertex */
							for (i = 0; i < edgdir[fl].no; i++) {
								int eix;

								/* Edge vertex index number of other vertex */
								if (edg->ti[i].goffs[0] != 0)
									eix = vx->ix + edg->ti[i].goffs[0];
								else
									eix = vx->ix + edg->ti[i].goffs[1];
				
								if ((nvx = get_vtxrec(&vc, eix)) != NULL) {
//printf("vx %d neighbor vx %d status %d\n",vx->ix,nvx->ix,nvx->status);
								 	if (nvx->status >= vtx_sha) {
										nshad++;
									} else {
										nnshad++;
									}
								}
							}
						}
//printf("vx %d nshad %d nnshad %d\n",vx->ix);
#endif	/* DELETE_SHAD */

						/* Keep un-shadowed vertexes and */
						/* shadowes ones that have non-shadows neigbors */
						if (vx->status == vtx_norm
#ifdef DELETE_SHAD
						 || vx->status >= vtx_sha && nnshad != 0
#endif
						) {
							*nrp++ = *rp;
//printf("~1 leaving vtx %d status %d in bxcell %d list\n",vx->ix,vx->status,bx->ix);

							if (phase == 1) {
#ifndef NEVER 					/* Do additions */
								/* If vertex doesn't land in a surface bxcell, */
								/* create a new surface bxcell for it. */
								if (vx->status == vtx_norm			// ????
								 && (vflag[vx->rix] & 2) == 0) {
									bxcell *nx;
				
//if (get_surface_bxcell(s, vx->rix) != NULL)
//error("new addition bx %d is already surface cell!\n",vx->rix);

#if defined(REVTABLESTATS) || defined(DEBUG)
									nascells++;
#endif
									/* Since it's a surface point, the seeding point is itself (NULL). */
									nx = new_bxcell(s, vx->rix, vx->ival, NULL, 0.0, NULL);
									add_bxcell_hash(s, nx);
									/* Convert to empty surface cell */
									vflag[nx->ix] = (vflag[nx->ix] & ~0xf) | 2;
			
									/* Add to surface linked list */
									nx->slist = s->rev.surflist;
									s->rev.surflist = nx;

									/* Add to additions list */
									nx->xlist = xlist;
									xlist = nx;
//printf("Added bxcell %d, status %d due to vx %d status %d\n",nx->ix, nx->status,vx->ix,vx->status);

#ifdef REVVRML
									vx->addvtx = 1;		/* Cause of added bxcell */
									nx->debug = 1;		/* Mark added bxcells */
#endif
								}
#if defined(REVTABLESTATS) || defined(DEBUG)
								/* Keep addvtx flag straight */
								else if (vx->status == vtx_norm) {
									bxcell *nx;
									if ((nx = get_surface_bxcell(s, vx->rix)) != NULL) {
										if (nx->status == bx_uninit)	/* Must be just added */
											vx->addvtx = 1;				/* Cause of added bxcell */
									}
								}
#endif
							}
#endif	/* Do additions */
						/* Omit vertex from bx list, and mark it as deleted, */
						/* and remove it if it has no un-shadowed neighbors */
						} else {
							vx->status = vtx_del;
//printf("~1 marking vtx %d status %d nnshad %d deleted bxcell %d list\n",vx->ix,vx->status,nnshad,bx->ix);
#ifdef DELETE_SHAD
							/* Remove it from cache if all its neighbors are */
							/* shadowed too. */
							if (nnshad == 0) {
//printf("~1 deleting vtx %d\n",vx->ix);
								del_vtxrec_hash(&vc, vx->ix);
								if (get_vtxrec(&vc, vx->ix) != NULL)
									error("get_vtxrec suceeded after del_vtxrec_hash!");
							}
#else /* !DELETE_SHAD */
							/* Keep track of deleted verticies that are in this bx, */
							/* so we can add back in crossing triangle vertexes */
							add2indexlist(s, &bx->dl, vx->ix, 0);
#endif /* !DELETE_SHAD */
						}
					}			/* Next vertex in bx's list */
					*nrp = -1;
					bx->sl[1] = nrp - bx->sl;

//aftercount = bx->sl[1]-3;
//if (beforecount != 0 && aftercount < beforecount) printf("Reduced bx from %d to %d verticies\n",beforecount,aftercount);
				}				/* Next bx */
			}				/* Next phase */

#ifdef REVVRML
			/* Main summary plot at each thinning round */
			/* Vtx ix tag ? Deleted vtxs ? Added vtxs ? Preserved vtxs ? oil ? bxcells ? Wait ? */
			if (0) plot_vtx_surface(s, 0, 0, 1, 0, 0, 0, 1, &vc, edgdir); 
#endif /* REVVRML */

			if (xlist == NULL) {
				break;			/* No added surface cells */
			}

			DBG(("reseting shadowers of new bxcells\n"));

			/* Locate all the bxcells that shadow the added bxcells, */
			/* and revert status to rethinned. */
			for (bx = xlist; bx != NULL; bx = bx->xlist) {
		
#ifdef REVVRML
				for (nbx = s->rev.surflist; nbx != NULL; nbx = nbx->slist)
					nbx->debug = 0;
#endif

				/* Locate the nnrev[] bxcells that shadow this added bxcell */
				bx->wlist = NULL;		/* For debug */
				for (nbx = s->rev.surflist; nbx != NULL; nbx = nbx->slist) {

					if (
#ifdef REVVRML
					    (nbx->status != bx_thinned && nbx->status != bx_filled)	 // Show all
#else
					    (nbx->status != bx_thinned)
#endif
					 || nbx == bx
					 || nbx->sl == NULL
					 || nbx->sl[1] == 3)
						continue;
				
					/* If any of nbx is further from bx and their bounding cylinders */
					/* overlap in perspective from rev.ocenter, assume nbx is a shadower. */
					if (shadow_group_group(s, s->rev.ocent, nbx->g.bcent, nbx->cc, nbx->dw,   
						                                    bx->g.bcent, bx->cc, bx->dw)) {
						nbx->status = bx_rethinnd;

#ifdef REVVRML
						bx->debug = 1;				/* rethinned bx */
						nbx->debug = 2;				/* added bx */
						nbx->wlist = bx->wlist;		/* For debug */
						bx->wlist = nbx;
#endif
//printf("~1 marking bxcell %d as un-thinned due to added bxcell %d\n",nbx->ix, bx->ix);
					}
				}

#ifdef REVVRML
				/* Plot bxcells touched by added cell */
				if (0) plot_touched_bxcells(s, bx->ix);
#endif	/* VRML */
			}

#if defined(REVTABLESTATS) || defined(DEBUG)
			printf(" %d bxcells thinned, %d re-thinned\n",thcount,rethcount);
			printf("Loop took %f seconds\n",0.001 * (msec_time()-lmsec));
#endif
		}		/* Loop until done */

#if defined(REVTABLESTATS) || defined(DEBUG)
		printf("Thinning took %f seconds\n",0.001 * (msec_time()-smsec));
#endif

		/* = = = = = = = = = = = = = = = = = = */
		DBG(("Preserving overlapping triangles\n"));
		{
#ifdef REVTABLESTATS
			int notverts = 0;			/* Number of possible crossed triangles/test verticies */
			int nopreserved = 0;		/* Number of verticies preseved for crossied triangles */
			unsigned long lmsec = msec_time();
#endif
			int sdi = 2;				/* sub-simplexes are triangles */
			int k, jj;
			vtxrec *vx;

			/* Struct to hold test vertex locations */
			struct _tvxrec {
				double v[MXRO];			/* Log output vertex value */
				double dist;			/* Distance from center point squared */
				int ix[MXRO+1];			/* Indexes of the triangle verticies */
				int shad;				/* Test result */
				struct _tvxrec *tlist;
			}; typedef struct _tvxrec tvxrec;
			tvxrec *tlist = NULL, *ftlist = NULL, *tvx, *ntvx;
			int nitlist = 0;

			/* For each surface bxcell, form triangles from vertexes */
			/* and detect possible crossed triangles */
			for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
				int sdi = 2;						/* sub-simplexes are triangles */
				double clb[MXRO+1];	/* Line RHS implicit equation vector [fdi+1] */
				int *crp, *rp, *nrp;
				vtxrec *vx, *nvx;
				int aftercount;			/* vertex count after thinning */

				/* Skip cell if empty */
				if (bx->sl == NULL || bx->sl[1] == 3)
					continue;

				/* Put the testing triangle verticies on the vtxlist */
				vc.vtxlist = NULL;
				vc.nilist = 0;

				/* Be able to detect triangles already tested */
				/* from this shadowing bxcell. */
				clear_trirec(s, &tc);

				/* See whether to add cell verticies to the list. */
				for (rp = bx->sl+3; *rp != -1; rp++) {
					assdire *tri;			/* Triangle table */
					float *fp;
					int fl;
					int added = 0;

					if ((vx = get_vtxrec(&vc, *rp)) == NULL)
						error("Failed to find vertex %s in cache",*rp);

					if (vx->status != vtx_norm)		// ???
						continue;

					fp = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */
					fl = FLV(fp);		/* Edge flags for this vertex */
					tri = tridir + fl;

					/* For all +ve triangles that use this vertex */
					for (k = 0; k < tridir[fl].no; k++) {
						int triix[MXRI+1];
						vtxrec *trivx[3];
						int ntvsh = 0;			/* Number of verticies shadowed */
						int nntvsh = 0;			/* Number of verticies not shadowed */

//printf("~1 tri %d: goffs = %s\n", k, debPiv(sdi+1, tri->ti[k].goffs));

						/* Triangle vertex index numbers */
						for (j = 0; j <= sdi; j++) {
							triix[j] = vx->ix + tri->ti[k].goffs[j];

							if ((trivx[j] = get_vtxrec(&vc, triix[j])) == NULL) {
								break;		/* Vertex doesn't exist */
							}
						    if (trivx[j]->status != vtx_norm)
								ntvsh++;
							else
								nntvsh++;
						}

						/* If a vertex isn't valid, or all vertexes are shadowed or not shadowed */
						if (j <= sdi
						 || ntvsh == (sdi+1)
						 || nntvsh == (sdi+1)) {
//printf("~1 vtx missing %d, ntvsh %d, nntvxsh %d\n",j <= sdi, ntvsh, nntvsh);
								continue;			/* Skip this triangle */
						}

						/* If triangle has been done before for this bxcell, skip it. */
						if (check_trirec(s, &tc, triix)) {
							continue;
						}

						/* We've decided to add triangle and test vertex */
						if (!added) {
							add_vtxrec_list(&vc, vx, 1);	/* Add vertex to list to test against */
							added = 1;
						}

						/* Create or re-use test vertex */
						if (ftlist != NULL) {		/* Grab one from free list */
							tvx = ftlist;
							ftlist = tvx->tlist;
							memset((void *)tvx, 0, sizeof(tvxrec));
					
						} else {
							if ((tvx = (tvxrec *) rev_calloc(s, 1, sizeof(tvxrec))) == NULL)
								error("rspl malloc failed - rev tvxrec structs");
							INCSZ(s, sizeof(tvxrec));
						}

						tvx->tlist = tlist;
						tlist = tvx;
						nitlist++;

						for (f = 0; f < fdi; f++)
							tvx->v[f] = 0.0;

						for (j = 0; j <= sdi;  j++) {
							if (trivx[j]->status == vtx_norm) {
								for (f = 0; f < fdi; f++)
									tvx->v[f] += 0.95/nntvsh * trivx[j]->vl[f];
							} else {
								for (f = 0; f < fdi; f++)
									tvx->v[f] += 0.05/ntvsh * trivx[j]->vl[f];
								trivx[j]->cross = 1;		/* For diagnostics */
							}
						}

						/* Compute distance of test vertex to overall center point squared */
						tvx->dist = 0.0;
						for (f = 0; f < fdi; f++) {
							double tt = tvx->v[f] - s->rev.ocent[f]; 
							tvx->dist += tt * tt;
						}

						/* Note the triangles vertexes indexes */
						for (j = 0; j <= sdi;  j++)
							tvx->ix[j] = trivx[j]->ix;
#ifdef REVTABLESTATS
						notverts++;
#endif
					}
				}

				/* Do a first pass for each test vertex, testing against */
				/* just the triangles that are associated with it's triangle. */
				/* (This quickly culls the test vertex list size, greatly */
				/* reducing the time taken in the second pass */

				/* For each test vertex */
				for (tvx = tlist; tvx != NULL; tvx = tvx->tlist) {
					double pv[MXRO];	/* Vertex being tested */
					double de[MXRO];	/* Line delta */

					clear_trirec(s, &stc);

					/* Compute line delta */
					for (f = 0; f < fdi; f++) {
						pv[f] = tvx->v[f];
						de[f] = pv[f] - s->rev.ocent[f];
					}

					/* Setup line cla and clb */
					init_line_eq_imp(s, NULL, &cla, clb, s->rev.ocent, de, 0);
	
					/* For each vertex of the test vertex triangle */
					for (jj = 0; jj <= sdi; jj++) {
						assdire *tri;			/* Triangle table */
						float *fp;
						int fl;

						if ((vx = get_vtxrec(&vc, tvx->ix[jj])) == NULL)
							error("rev crossing test - failed to get vertex");

						if (vx->status != vtx_norm)
							continue;

						fp = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */
						fl = FLV(fp);		/* Edge flags for this vertex */
						tri = tridir + fl;

						/* For all +ve triangles that use this vertex */
						for (k = 0; k < tridir[fl].no; k++) {
							int triix[MXRI+1];
							vtxrec *trivx[MXRI+1];
							double v[MXRI+1][MXRO]; 	/* Triangle vertex values */
							double gc[MXRO], cc, dw;	/* Triangle shadow group info. */
							int ntvsh = 0;			/* Number of verticies shadowed */
							double bdist = -1.0;
							double tb[MXRI];	/* Solution point in input space */
							double xv[MXRO];	/* Solution point in output space */
							int g, sorv, wsrv;	/* Solved & within simplex return value */
							double dist;		/* distance to line origin */
							double dot;			/* dot product of solution to line */

//printf("~1 tri %d: goffs = %s\n", i, debPiv(sdi+1, tri->ti[k].goffs));

							/* Triangle vertex index numbers */
							triix[0] = vx->ix + tri->ti[k].goffs[0];
							triix[1] = vx->ix + tri->ti[k].goffs[1];
							triix[2] = vx->ix + tri->ti[k].goffs[2];

							/* If triangle has been done before for this tvx, skip it. */
							if (check_trirec(s, &stc, triix)) {
								continue;
							}

							/* Triangle vertex index numbers */
							for (j = 0; j <= sdi; j++) {
//								triix[j] = vx->ix + tri->ti[k].goffs[j];

								if ((trivx[j] = get_vtxrec(&vc, triix[j])) == NULL) {
									break;		/* Vertex doesn't exist */
								}
							    if (trivx[j]->status != vtx_norm)
									ntvsh++;

								if (trivx[j]->dist > bdist)
									bdist = trivx[j]->dist;
							}

							/* If vertex is above triangle, it can't be shadowed */
							if (tvx->dist > bdist)
								continue;

							/* If a vertex isn't valid, or all vertexes are shadowed */
							if (j <= sdi
							 || ntvsh >= (sdi+1)) {
									continue;			/* Skip this triangle */
							}

							/* If this triangle is the test vertex triangle, skip it */
							if (tvx->ix[0] == triix[0]
							 && tvx->ix[1] == triix[1]
							 && tvx->ix[2] == triix[2]) {
								continue;
							}

							for (j = 0; j <= sdi; j++) {
								for (f = 0; f < fdi; f++)
									v[j][f] = trivx[j]->vl[f];
							}

							/* Compute shadow group params of triangle for quick vertex test */
							comp_shadow_group(s, s->rev.ocent, gc, &cc, &dw, NULL, v, sdi+1);    

							/* Do quick check against triangle */
							if (!shadow_group_vertex(s,
								     s->rev.ocent, gc, cc, dw, tvx->v)) { 
								continue;
							}

//printf("~1 checking vertex %d at %s dist %f\n",tvx->ix, debPdv(fdi,tvx->v), sqrt(tvx->dist));
							/* Compute intersection: */
							wsrv = 0;

							/* Solve line/triangle intersection using same */
							/* method as vnearest_clip_solve(). */
						
							/* LHS: ta[sdi][sdi] = cla[sdi][fdi] * vv[fdi][sdi] */
							/* RHS: tb[sdi] = clb[sdi] - cla[sdi][fdi] * vv_di[fdi] */
							for (f = 0; f < sdi; f++) {
								double tt;
								for (e = 0; e < sdi; e++) {
									for (tt = 0.0, g = 0; g < fdi; g++)
										tt += cla[f][g] * (v[e][g] - v[e+1][g]);
									ta[f][e] = tt;
								}
								for (tt = 0.0, g = 0; g < fdi; g++)
									tt += cla[f][g] * v[sdi][g];
								tb[f] = clb[f] - tt;
							}
						
							/* Compute the solution */
							/* (Solve the simultaneous linear equations A.X = B) */
//							sorv = !solve_se(ta, tb, sdi);
							sorv = !solve_se_2x2(ta, tb);		/* Saves a few % only */

							if (!sorv)
								continue;
							
							/* Check that the solution is within the simplex & meets ink limit */
							if ((wsrv = simple_within_simplex(v, tb, sdi)) != 0) {

								/* Compute the output space solution point */
								for (f = 0; f < fdi; f++) {
									double tt = 0.0;
									for (e = 0; e < sdi; e++)
										tt += (v[e][f] - v[e+1][f]) * tb[e];
									xv[f] = tt + v[sdi][f];
								}
						
								/* Compute distance to gamut center squared, */
								/* as well as the dot product */
								for (dot = dist = 0.0, f = 0; f < fdi ; f++) {
									double tt = (xv[f] - s->rev.ocent[f]);
									dist += tt * tt;
									dot += de[f] * tt;
								}
   
								/* If intersection distance is greater than vertex distance, */
								/* mark the test vertex as shadowed (== crossed triangle */
								/* is shadowed) */
								if (dot > 0.0 && dist > (tvx->dist + EPS)) {
									tvx->shad = 1;
									goto next_tvx;
								}
							}
						}	/* Next associated triangle */
					}		/* Next vertex of test triangle */
				  next_tvx:;
				}			/* Next test vertex */

				/* Delete shadowed tvx, and sort remaining tlist by distance so */
				/* that we have a better chance of shadowing it early ? */
				{
					int i;
					tvxrec **sort, *vx, *nvx;
				
					/* Create temporary array of pointers to tvxrec's in list */
					if ((sort = (tvxrec **) rev_calloc(s, nitlist, sizeof(tvxrec *))) == NULL)
						error("rspl malloc failed - rev tvxrec sort array");
					INCSZ(s, nitlist * sizeof(tvxrec *));
				
					for (i = 0, vx = tlist; vx != NULL; vx = nvx) {
						nvx = vx->tlist;
						if (!vx->shad) {
							sort[i++] = vx;
						} else {
							/* Put deleted tvxrec on the free list to re-use */
							vx->tlist = ftlist;
							ftlist = vx;
						}
					}
					nitlist = i;
				
					/* Sort the list into ascending distance from center */
#define 	HEAP_COMPARE(A,B) (A->dist < B->dist)
					HEAPSORT(tvxrec *, sort, nitlist)
#undef 		HEAP_COMPARE
				
					/* Re-create the linked list in descending order */
					tlist = NULL;
					for (i = 0; i < nitlist; i++) {
						vx = sort[i];
						vx->tlist = tlist;
						tlist = vx;
					}
				
					free(sort);
					DECSZ(s, nitlist * sizeof(tvxrec *));

#ifdef NEVER
					printf("sorted test vertex list:\n");
					for (i = 0, vx = tlist; vx != NULL; vx = vx->tlist, i++)
						printf("%d: ix %s dist %f\n",i,debPiv(3,vx->ix), sqrt(vx->dist));
#endif
				}

				/* Be able to detect triangles already tested */
				/* from this shadowing bxcell. */
				clear_trirec(s, &tc);

				/* sort vertexes by descending distance to center point */
				/* (and also reset list tflag), to detect shadowing early */
				sort_vtxrec_list(s, &vc);
	
				/* Check if the test points are shadowed by any triangle */
				for (vx = vc.vtxlist; vx != NULL; vx = vx->tlist) {
					assdire *tri;			/* Triangle table */
					float *fp;
					int fl;

					if (vx->status != vtx_norm)		// ???
						continue;

					fp = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */
					fl = FLV(fp);		/* Edge flags for this vertex */
					tri = tridir + fl;

					/* For all +ve triangles that use this vertex */
					for (k = 0; k < tridir[fl].no; k++) {
						int triix[MXRI+1];
						vtxrec *trivx[MXRI+1];
						double v[MXRI+1][MXRO]; 	/* Triangle vertex values */
						double gc[MXRO], cc, dw;	/* Triangle shadow group info. */
						int ntvsh = 0;			/* Number of verticies shadowed */
						double bdist = -1.0;

//printf("~1 tri %d: goffs = %s\n", i, debPiv(sdi+1, tri->ti[k].goffs));

						/* Triangle details */
						for (j = 0; j <= sdi; j++) {
							triix[j] = vx->ix + tri->ti[k].goffs[j];

							if ((trivx[j] = get_vtxrec(&vc, triix[j])) == NULL) {
								break;		/* Vertex doesn't exist */
							}
						    if (trivx[j]->status != vtx_norm)
								ntvsh++;

							if (trivx[j]->dist > bdist)
								bdist = trivx[j]->dist;
						}

						/* If a vertex isn't valid, or all vertexes are shadowed */
						if (j <= sdi
						 || ntvsh >= (sdi+1)) {
								continue;			/* Skip this triangle */
						}

						/* If triangle has been done before for this bxcell, skip it. */
						if (check_trirec(s, &tc, triix)) {
							continue;
						}

						for (j = 0; j <= sdi; j++) {
							for (f = 0; f < fdi; f++)
								v[j][f] = trivx[j]->vl[f];
						}

						/* Compute shadow group params of triangle for quick vertex test */
						comp_shadow_group(s, s->rev.ocent, gc, &cc, &dw, NULL, v, sdi+1);    

						/* For all test vertexes */
						for (tvx = tlist; tvx != NULL; tvx = tvx->tlist) {
							double pv[MXRO];	/* Vertex being tested */
							double de[MXRO];	/* Line delta */
							double tb[MXRI];	/* Solution point in input space */
							double xv[MXRO];	/* Solution point in output space */
							int g, sorv, wsrv;	/* Solved & within simplex return value */
							double dist;		/* distance to line origin */
							double dot;			/* dot product of solution to line */

							/* If vertex is above triangle, it can't be shadowed */
							if (tvx->dist > bdist)
								continue;

							/* If we have already determined this one is shadowed */
							if (tvx->shad)
								continue;

							/* If this vertex for this triangle, skip it */
							if (tvx->ix[0] == triix[0]
							 && tvx->ix[1] == triix[1]
							 && tvx->ix[2] == triix[2]) {
								continue;
							}

							/* Do quick check against triangle */
							if (!shadow_group_vertex(s, s->rev.ocent, gc, cc, dw, tvx->v))
								continue;
//printf("~1 checking vertex %d at %s dist %f\n",tvx->ix, debPdv(fdi,tvx->v), sqrt(tvx->dist));
							/* Compute intersection: */
							wsrv = 0;

							/* Compute line delta */
							for (f = 0; f < fdi; f++) {
								pv[f] = tvx->v[f];
								de[f] = pv[f] - s->rev.ocent[f];
							}

							/* Setup line cla and clb */
							init_line_eq_imp(s, NULL, &cla, clb, s->rev.ocent, de, 0);

							/* Solve line/triangle intersection using same */
							/* method as vnearest_clip_solve(). */
						
							/* LHS: ta[sdi][sdi] = cla[sdi][fdi] * vv[fdi][sdi] */
							/* RHS: tb[sdi] = clb[sdi] - cla[sdi][fdi] * vv_di[fdi] */
							for (f = 0; f < sdi; f++) {
								double tt;
								for (e = 0; e < sdi; e++) {
									for (tt = 0.0, g = 0; g < fdi; g++)
										tt += cla[f][g] * (v[e][g] - v[e+1][g]);
									ta[f][e] = tt;
								}
								for (tt = 0.0, g = 0; g < fdi; g++)
									tt += cla[f][g] * v[sdi][g];
								tb[f] = clb[f] - tt;
							}
						
							/* Compute the solution */
							/* (Solve the simultaneous linear equations A.X = B) */
//							sorv = !solve_se(ta, tb, sdi);
							sorv = !solve_se_2x2(ta, tb);		/* Saves a few % only */

							/* If it was solved */
							if (sorv) {
							
								/* Check that the solution is within the simplex & ink limit */
								if ((wsrv = simple_within_simplex(v, tb, sdi)) != 0) {

									/* Compute the output space solution point */
									for (f = 0; f < fdi; f++) {
										double tt = 0.0;
										for (e = 0; e < sdi; e++)
											tt += (v[e][f] - v[e+1][f]) * tb[e];
										xv[f] = tt + v[sdi][f];
									}
							
									/* Compute distance to gamut center squared, */
									/* as well as the dot product */
									for (dot = dist = 0.0, f = 0; f < fdi ; f++) {
										double tt = (xv[f] - s->rev.ocent[f]);
										dist += tt * tt;
										dot += de[f] * tt;
									}
//printf("~1 intersection at %s dist %f\n", debPdv(fdi,xv), sqrt(dist));
   
									/* If intersection distance is greater than vertex distance, */
									/* mark the test vertex as shadowed (== crossed triangle */
									/* is shadowed) */
									if (dot > 0.0 && dist > (tvx->dist + EPS)) {
										tvx->shad = 1;
									}
								}
							}
						}	/* Next test vertex */
					}		/* Next triangle from vertex */
				}			/* Next vertex */

				/* Go through test vertex results, and if it is un-shadowed, */
				/* mark all the corresponding triangle vertexes as un-shadowed. */
				/* For all test vertexes */
				for (tvx = tlist; tvx != NULL; tvx = ntvx) {
					ntvx = tvx->tlist;

					/* If the test point wasn't shadowed, assume it */
					/* is part of the gamut surface, and mark all its */
					/* vertexes as valid. */
					if (!tvx->shad) {
						for (j = 0; j <= sdi; j++) {
							if ((vx = get_vtxrec(&vc, tvx->ix[j])) == NULL)
								error("rev - failed to locate vertex %d\n",tvx->ix[j]);

							if (vx->status != vtx_norm) {
								vx->pres = 1;		/* Don't treat it as deleted */
							}
						}
					}
					/* Put all the tvxrec's on the free list to re-use */
					tvx->tlist = ftlist;
					ftlist = tvx;
				}
				tlist = NULL;
				nitlist = 0;

				/* If the preseved vertexes have been deleted from the bx list, */
				/* add them back in again */
				if (bx->dl != NULL) {
					for (nrp = rp = bx->dl+3; *rp != -1; rp++) {
						vtxrec *vx;
		
						if ((vx = get_vtxrec(&vc, *rp)) == NULL)
							continue;		/* Hmm. */
		
						/* If preserved, transfer it to the active bx list */
						if (vx->pres) {
							add2indexlist(s, &bx->sl, *rp, 0);
	
						/* Leave it in deleted list */
						} else {
							*nrp++ = *rp;
						}
					}
					*nrp = -1;
					bx->dl[1] = nrp - bx->dl;
	
					/* We don't need the deleted list now */
					free_indexlist(s, &bx->dl);
				}

			}		/* Next bxcell */

			/* Free up tvxrec's */
			while (ftlist != NULL) {
				tvxrec *this = ftlist;
				ftlist = ftlist->tlist;
				free(this);
				DECSZ(s, sizeof(tvxrec));
			}

#ifdef REVTABLESTATS
			/* Count the number of preserved vertexes */
			for (i = 0; i < vc.hash_size; i++) {
				for (vx = vc.hash[i]; vx != NULL; vx = vx->hlink) {
					if (vx->pres)
						nopreserved++;
				}
			}
			printf("%d crossed triangles tested\n",notverts);
			printf("%d hidden verticies retained for crossed triangles\n",nopreserved);
			printf("Took %f secs to preserving crossing triangless\n",0.001 * (msec_time()-lmsec));
#endif
		}	/* End of preserve shadowed triangles */

		/* = = = = = = = = = = = = = = = = = = */
		/* Delete any shadowed vertexes, and remove any empty bxcells. */
		for (pbx = &s->rev.surflist, bx = *pbx; bx != NULL; bx = nbx) {
			int *rp, *nrp;

			/* Delete all the shadowed or delted vertexes from bxcell list, */
			/* unless they are preserved because they are part of a crossed triangle. */
			for (nrp = rp = bx->sl+3; *rp != -1; rp++) {
				vtxrec *vx;

				if ((vx = get_vtxrec(&vc, *rp)) == NULL)
					continue;		/* Hmm. Delete it.*/

				/* Keep all the un-shadowed or preserved vertexes */
				if (vx->status == vtx_norm	
				  || vx->pres) {
					*nrp++ = *rp;
				} else {
					del_vtxrec_hash(&vc, vx->ix);
				}
			}
			*nrp = -1;
			bx->sl[1] = nrp - bx->sl;

			if (bx->sl == NULL	/* Missing or empty fwd index list */
			 || bx->sl[1] == 3
			) {
				/* Remove it from vflag array */
				if (s->rev.rev[bx->ix] != NULL) {
					vflag[bx->ix] = (vflag[bx->ix] & ~0xf) | 1;		/* Not surface and done */
				} else {
					vflag[bx->ix] = (vflag[bx->ix] & ~0xf) | 0;		/* Not surface and empty */
				}

				/* Remove it from hash */
				rem_bxcell_hash(s, bx->ix);

				/* Free fwd index list (none are shared at this point) */
				if (bx->sl != NULL)
					free_indexlist(s, &bx->sl);

				/* Remove it from surface list */
				*pbx = nbx = bx->slist;

				/* Free it */
				del_bxcell(s, bx);
#if defined(REVTABLESTATS) || defined(DEBUG)
				nrscells++;
#endif

			} else {						/* Move on to next */
				pbx = &bx->slist;
				nbx = bx->slist;
			}
		}

		/* Add extra over ink limit vertexes. */
		for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
			int sdi = 1;					/* sub-simplexes are edges */
			int *rp;
			vtxrec *vx, *nvx;
			int *exlist = NULL;
	
			/* Add over ink limit vertexes, so that fwd cells will straddle */
			/* the ink limit boundary. */
			/* Do this by checking all vertexes edge neighbors, */
			/* and adding any that are over the ink limit. */
			/* (Only do this for bx cells that are known to contain */
			/* over ink limit verticies.) */ 
			if (s->limiten && vflag[bx->ix] & 0x10) {
				int *rp;

//printf("~1 ink limitin is enabled bx %d\n", bx->ix);
				for (rp = bx->sl+3; *rp != -1; rp++) {
					float *vp, *evp;
					int fl;
					assdire *edg;			/* Edge table */

					if ((vx = get_vtxrec(&vc, *rp)) == NULL)
						continue;		/* Hmm. */

					/* Don't do this for preserved or oil vertexes */
					if (vx->status != vtx_norm)
						continue;
 
					vp = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */
					fl = FLV(vp);		/* Edge flags for this vertex */
					edg = edgdir + fl;

#ifdef CHECK_NNLU
					if (vp[-1] > s->limitv)
						error("Thinned vertex %d is over ink limit!",vx->ix);
#endif

//printf("~1 fl %d = 0o%o, no edges %d\n",fl, fl, edg->no);

					/* For all possible edges that use this vertex */
					for (i = 0; i < edgdir[fl].no; i++) {
						int eix;

//printf("~1 edg %d: goffs = %s\n", i, debPiv(sdi+1, edg->ti[i].goffs));

						/* Edge vertex index number of other vertex */
						if (edg->ti[i].goffs[0] != 0)
							eix = vx->ix + edg->ti[i].goffs[0];
						else
							eix = vx->ix + edg->ti[i].goffs[1];
		
						evp = s->g.a + eix * s->g.pss;	/* Other vertex in fwd grid */

//printf(" Checking edge %d (%f) -> %d (%f)\n", vx->ix, vp[-1], eix, evp[-1]);

						/* If over limit, add it to the expansion list */
						if (evp[-1] > s->limitv) {
//printf("~1 added over ink limit vertex %d\n",eix);

							if (get_vtxrec(&vc, eix) != NULL)
								continue;		/* Added by another bx */
							nvx = new_vtxrec(s, &vc, eix);
							nvx->status = vtx_oil; 
							add2indexlist(s, &exlist, eix, 0);
						}
					}
				}
		
				/* If we found over ink limit verticies, add them to our list */
				if (exlist != NULL) {
					for (rp = exlist+3; *rp != -1; rp++) {
						add2indexlist(s, &bx->sl, *rp, 0);
					}
					free_indexlist(s, &exlist);
				}
			}
		}

#ifdef REVTABLESTATS
		/* Count the number of over ink limit vertexes */
		for (i = 0; i < vc.hash_size; i++) {
			vtxrec *vx;
			for (vx = vc.hash[i]; vx != NULL; vx = vx->hlink) {
				if (vx->status == vtx_oil)
					naoulvtxs++;
			}
		}
#endif

#ifdef REVVRML
		/* Plot final vertex surface before converting to fwcells */
		/* Vtx ix tag ? Deleted vtxs ? Added vtxs ? Preserved vtxs ? oil vtxs ? bxcells ? Wait ? */
		if (1) plot_vtx_surface(s, 0, 0, 0, 1, 1, 0, 0, &vc, edgdir); 
#endif /* REVVRML */

		/* Convert vertexes to cube lists */
		for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
			int sdi = 1;					/* sub-simplexes are edges */
			int *crp, *rp, *nrp;
			int ttouch;
			vtxrec *vx, *nvx;
	
			/* If there are no vertexes left (i.e. they have all been deleted) */
			/* Don't try and convert to fwd cells. */
			if (bx->sl == NULL || bx->sl[1] == 3) {
				bx->status = bx_conv;
				continue;
			}
	
			/* Create cach list of vxrec's for just this nnrev[] */
			clear_vtxrec_lists(s, &vc);

			/* Add all this bxcell verticies to cache and list */
			for (rp = bx->sl+3; *rp != -1; rp++) {
				vx = new_vtxrec(s, &vc, *rp); 
				add_vtxrec_list(&vc, vx, 0);
			}

			/* Convert fwd index list into fwd cells list. Do this in */
			/* a way that minimizes the number of cells needed while still */
			/* ensuring that there is 2 dimensional connectivity for all the vertexes. */
	
			/* Count number of touches if we add a cube for each prime vertex */ 
//printf("~1 counting number of touches\n");
			crp = bx->sl;
			i = 0;
			for (rp = crp+3; *rp != -1; rp++) {
				vtxrec *vx;
	
				if ((vx = get_vtxrec(&vc, *rp)) == NULL)
					error("get_vtxrec() failed on surface vtx");
	
				i++;

				/* For each vertex of cube placed at vx->cix */
				for (ee = 0; ee < (1<<di); ee++) {
					int vix = vx->cix + s->g.hi[ee];
					vtxrec *nx;
	
					if ((nx = get_vtxrec(&vc, vix)) != NULL)
						vx->tcount++;
				}
			}
//printf("there were %d vertexes",i);
	
//printf("~1 adding cells in order of touch count\n");
			/* Add cells in order of touch count, i.e. from most necessary */
			/* to least necessary. Allow a maximum touch of 4, to ensure */
			/* 2 dimensional connectivity of the fwd cells */
			nrp = NULL;
			i = 0;
			for (ttouch = 1; ; ttouch++) {
				int more = 0;
//printf("~1 ttouch = %d\n",ttouch);
				for (rp = crp+3; *rp != -1; rp++) {
					vtxrec *vx = get_vtxrec(&vc, *rp); 
	
					if (vx->tcount == 0)
						continue;
	
					more = 1;
					if (vx->tcount > ttouch)
						continue;
	
					/* For each cube vertex placed at vx->cix */
					for (ee = 0; ee < (1<<di); ee++) {
						int vix = vx->cix + s->g.hi[ee];
						vtxrec *nx;
	
						/* Track touch count on creating cells, and */
						/* clear vertexes that have reached 4, */
						/* so that they don't get any more */
						if ((nx = get_vtxrec(&vc, vix)) != NULL) {
//printf("bx %d, adding fwcell vertex %d for vertex %d\n",bx->ix,vix,*rp);  
							vx->acount++;
							if (vx->acount >= 4)
								vx->tcount = 0;
						}
					}
					i++;
					add2indexlist(s, &nrp, vx->cix, 0);
				}
				if (!more)
					break;
			}
//printf(", now %d fwdcells\n",i);
//printf("~1 replacing vertex list with cell list\n");
	
			if (nrp == NULL)
				error("Surface list bxcell ix %d has no fwd cells",bx->ix);
	
			/* Replace vertex list with cell list */
			free_indexlist(s, &bx->sl);
			bx->sl = nrp;
	
			if (bx->sl == NULL)
				error("Surcface cell nnrev[%d] is empty!\n",bx->ix);
			bx->status = bx_conv;
		}

		if (s->rev.surflin != NULL) {	/* Don't need surflin anymore */
			s->rev.surflin->del(s->rev.surflin);
			s->rev.surflin = NULL;
			s->rev.surflin_en = 0;
		}
		if (cla != NULL)
			free_dmatrix(cla, 0, fdi-1, 0, fdi);
		free_trirec(s, &stc);
		free_trirec(s, &tc);
		free_vtxrec_list(s, &vc);
		free_assdir(s, edgdir);
		free_assdir(s, tridir);
	}

#if defined(REVTABLESTATS) || defined(DEBUG)
	if (fdi > 1) {
		bxcell *bx;
		int surfcelldepth = 0, surfcells = 0;
		for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
			if (bx->sl == NULL
			 || bx->sl[1] == 3)
				continue;
			surfcells++;
			surfcelldepth += bx->sl[1]-3;
		}

		printf("%d/%d surface cells\n",surfcells,rgno);
		printf("%d/%d non-surface cells\n",ingamutcells,rgno);
		printf("%d/%d empty cells\n",emptycells,rgno);
		printf("%d/%d used cells in rev[]\n",revcells,rgno);
		printf("%f average rev[] list length\n",(double)revcelldepth/(double)revcells);
		printf("%f average nnrev[] surface list length\n",(double)surfcelldepth/(double)surfcells);
		printf("%d added surface cells\n",nascells);
		printf("%d removed surface cells\n",nrscells);
		printf("%d added over ink limit vertexes\n",naoulvtxs);
	}
#endif

#ifdef REVVRML
	/* Plot the thinned surface fwd cells */
	/* fwd cell base ix's ? bxcells ? Wait ? */
	if (1 && fdi > 1) plot_fxcell_surface(s, 0, 0, 0);
#endif /* REVVRML */

	/* Fill the non-surface nnrev array from the surface list. */
	{
		bxcell *seedlist = NULL;		/* Linked list of active seeds */
		bxcell *seedlistend = NULL;		/* Last item on seedlist */
		bxcell *xlist = NULL;			/* Linked list of cells being searched */
		bxcell *xlistend = NULL;		/* Last item on xlist */
		bxcell *tlist;					/* Linked list of cells being considered as soln. */
		double emax;					/* Current smallest estimated max weigted distance */
#if defined(REVTABLESTATS) || defined(DEBUG)
		unsigned long smsec = msec_time();
#endif

		DBG(("Filling in rev.nnrev[] grid\n"));

		/* Start the seeding of the nnrev[] array with all the surface cells */
		{
			bxcell *ss;

			for (ss = s->rev.surflist; ss != NULL; ss = ss->slist) {
				/* Add to end of seedlist */
				ss->flist = NULL;
				if (seedlist == NULL)
					seedlist = ss;
				else
					seedlistend->flist = ss;
				seedlistend = ss;

				vflag[ss->ix] |= 1;		/* They are on seed list, so will be filled */
			}
		}

		/* While there are nnrev[] cells to fill */
		while (seedlist != NULL) {
			DCOUNT(cc, MXRO, fdi, -1, -1, 2);	/* bwd neighborhood offset counter */
			int nix;					/* Neighbor offset index */
			bxcell *ss, *tx;

			tx = seedlist;				/* Remove target cell from front of seed list */
			seedlist = tx->flist;

			if (s->rev.nnrev[tx->ix] != NULL)
				error("nncel[%d] in seed list is not empty\n",tx->ix);

#ifdef CHECK_NNLU
			if (tx->ss == NULL || (vflag[tx->ss->ix] & 2) == 0 ) {
				if (tx->ss == NULL)
					printf("nnrev[%d] has NULL seed\n",tx->ix);
				else
					printf("nnrev[%d] has seed %d with flag %x != 3\n",tx->ix,tx->ss->ix, vflag[tx->ss->ix]);
			}
#endif
				
			DBG(("Doing nnrev[%d] vflag %x co %s\n",tx->ix, vflag[tx->ix], debPiv(s->fdi, tx->gc)));
//printf("Doing nnrev[%d] vflag %x co %s\n",tx->ix, vflag[tx->ix], debPiv(s->fdi, tx->gc));

			emax = 1e200;				/* Smallest emax */
			ss = tx->ss;				/* Search start cell */ 
			ss->tix = tx->ix;			/* Mark this cell as being in search list */
			
			/* Make start cell the only entry in the search list */
			ss->xlist = NULL;
			xlist = ss;
			xlistend = ss;

			/* Clear the solution list */
			tlist = NULL;

			/* Note that filling an nnrev[] cell using a seeded search may miss fw cells */
			/* that should be in it, if they are in physically dis-continuous locations */
			/* due to gamut hull convexity. LCh weighting will reduce this somewhat, and */
			/* discontinuity is rarely a desired characteristic of a color conversion, so */
			/* we are ignoring this issue for now. */

			/* While there are cells to search for solutions */
			while (xlist != NULL) {
				double em, ex;

				ss = xlist;					/* Remove next search cell from linked list */
				xlist = xlist->xlist;

				/* Check if this cell could be in solution */
				em = nn_grpgrp_est(s, &ex, &tx->g, &ss->g);
				ss->emin = em;
#if defined(REVTABLESTATS) || defined(DEBUG)
				nnrevcellsearch++;
#endif

				DBG(("Searching rev[%d] co %s, em %f, ex %f\n",ss->ix, debPiv(s->fdi, ss->gc), em, ex));
//printf("Searching rev[%d] co %s, em %f, ex %f\n",ss->ix, debPiv(s->fdi, ss->gc), em, ex);
			
				if (em < emax) {		/* Yes */

					/* Add it to the solution list */
					ss->tlist = tlist;
					tlist = ss;

					DBG(("Adding it to solution list\n"));

					/* Update smallest maximum */
					/* (Will cull existing bxcell solutions with emin > emax later) */
					if (ex < emax)
						emax = ex;

					/* Explore all neighbours, and add any surface cells that haven't been */
					/* searched for this target yet. */ 
					DC_INIT(cc);
					while (!DC_DONE(cc)) {
						bxcell *nbx;

						nix = ss->ix;
						for (f = 0; f < fdi; f++) {
							nn[f] = ss->gc[f] + cc[f];
							if (nn[f] < 0 || nn[f] >= rgres)
								break;					/* Out of bounds */
							nix += cc[f] * s->rev.coi[f];
						}
						if (f < fdi || nix == ss->ix) {
//printf("Rejecting search neigbor co %s because out of bounds or current cell\n",debPiv(s->fdi,cc));
							goto next_neighbor;
						}

						/* We only search surface bxcells */
						if ((vflag[nix] & 2) == 0) {
//printf("Rejecting search neigbor nnrev[%d] co %s because flags = %x\n",nix, debPiv(s->fdi, cc),vflag[nix]);
							goto next_neighbor;
						}

						/* If neighbor is in bounds, and a surface bxcell*/
						{

							/* Expect all all surface bxcells to be in cache */
							if ((nbx = get_surface_bxcell(s, nix)) == NULL)
								error("rspl rev get_surface_bxcell %d failed",nix);

							/* If not already in search list */
							if (nbx->tix != tx->ix) {
//									DBG(("Adding search neigbor nnrev[%d] co %s to search list\n",nbx->ix, debPiv(s->fdi, nbx->gc)));
//printf("Adding search neigbor nnrev[%d] co %s to search list\n",nbx->ix, debPiv(s->fdi, nbx->gc));
								/* Add neigbor to end of search list */
								nbx->tix = tx->ix;		/* Is now in search list */
								nbx->xlist = NULL;
								if (xlist == NULL)
									xlist = nbx;
								else
									xlistend->xlist = nbx;
								xlistend = nbx;
							}
//else 
//printf("Rejecting search neigbor nnrev[%d] co %s because already in list\n",nbx->ix, debPiv(s->fdi, nbx->gc));
						}
					next_neighbor:;
						DC_INC(cc);
					}
				} 
//else
//printf("Rejected rev[%d] co %s, because em %f >= emax %f\n",ss->ix, debPiv(s->fdi, ss->gc), em, emax);
			}

			/* Create the nnrev[] list from the candidate bxcell solutions */
			if (tlist != NULL) {
				create_nnrev_list(s, tx, tlist, emax);  
			}
#if defined(REVTABLESTATS) || defined(DEBUG)
			nnrevcells++;
			nnrevcelldepth += s->rev.nnrev[tx->ix][1]-3;
			if (s->rev.nnrev[tx->ix][1]-3 > nnmxrevcelldepth)
				nnmxrevcelldepth = s->rev.nnrev[tx->ix][1]-3;
#endif

			/* If this was a super-cell, explore the 2nd row around this cell, */
			/* and locate any cells not on the seeding list */
			if (tx->scell != NULL) {
				DCOUNT(sc, MXRO, fdi, -3, -3, 4);
				DC_INIT(sc);
				while (!DC_DONE(sc)) {
					int co[MXRO]; 
					int ok = 0;
					int nix = tx->ix;

					for (f = 0; f < fdi; f++) {
						co[f] = tx->gc[f] + sc[f];
						if (co[f] < 0 || co[f] >= s->rev.res)
							break;
						nix += sc[f] * s->rev.coi[f];
						if (sc[f] == -3 || sc[f] == 3)
							ok = 1;				/* Just surface of +/- 2 */
					}
					if (!ok && sc[0] == -2)
						sc[0] = 2;			/* Skip center */

					/* Put this cell on list and stop searching. */
					if (f >= fdi && (vflag[nix] & 1) == 0) {

						if ((vflag[nix] & 2) != 0) {	/* If un-filled surface bxcell */
							/* Get surface bxcell from cache index for seed */
							if ((ss = get_surface_bxcell(s, nix)) == NULL)
								error("rspl rev get_surface_bxcell %d failed #2, vflag = %x",nix,vflag[nix]);
//printf("Fetched surface bxcell seed %d vflag %x\n",ss->ix, vflag[ss->ix]);
						} else {	/* If un-filled nnrev */
							if (get_surface_bxcell(s, nix) != NULL)
								error("vflag[%d] = %x, but cell is in surface list hash\n");

							/* Create new temporary (non-surface) bxcell seed. */
							/* If we are sufficiently far from the seed point, */
							/* a super-cell to improve seeding performance will be created. */
							ss = new_bxcell(s, nix, co, tx->ss, tx->sdist, vflag);
#if defined(REVTABLESTATS) || defined(DEBUG)
							if (tx->scell != NULL)
//								nnsuperfill += tx->scell[3]-3;
								nnsuperfill++;
							else
								nnsinglefill++;
#endif
//printf("Created temporary seed bxcell %d vflag %x\n",ss->ix, vflag[ss->ix]);
						}
						DBG(("Adding seed neighbor nnrev[%d] vflag %x co %s to seed list\n",ss->ix, vflag[ss->ix], debPiv(s->fdi, ss->gc)));
//printf("Adding seed neighbor nnrev[%d] vflag %x co %s to seed list\n",ss->ix, vflag[ss->ix], debPiv(s->fdi, ss->gc));

						/* Add to end of seedlist */
						ss->flist = NULL;
						if (seedlist == NULL)
							seedlist = ss;
						else
							seedlistend->flist = ss;
						seedlistend = ss;
						vflag[ss->ix] |= 1;		/* This is on seed list, so will be filled */
					}
					DC_INC(sc);
				}
			} else {
				/* Explore neighbours, and add any nnrev[] cells that haven't been */
				/* put on the seed list yet. */
				for (f = 0; f < fdi; f++)
					cc[f] = tx->gc[f];
				nix = tx->ix;

				for (ff = 0; ff < (fdi << 1); ff++) {
					f = ff >> 1;			/* Dimension being explored */

					cc[f] += (ff & 1) ? 1 : -1;
					nix   += (ff & 1) ? s->rev.coi[f] : -s->rev.coi[f];

					/* If found unfilled nnrev[] cell */
					if (cc[f] >= 0 && cc[f] < rgres && (vflag[nix] & 1) == 0) {

						if ((vflag[nix] & 2) != 0) {	/* If un-filled surface bxcell */
							/* Get surface bxcell from cache index for seed */
							if ((ss = get_surface_bxcell(s, nix)) == NULL)
								error("rspl rev get_surface_bxcell %d failed #2, vflag = %x",nix,vflag[nix]);
//printf("Fetched surface bxcell seed %d vflag %x\n",ss->ix, vflag[ss->ix]);
						} else {	/* If un-filled nnrev */
							if (get_surface_bxcell(s, nix) != NULL)
								error("vflag[%d] = %x, but cell is in surface list hash\n");

							/* Create new temporary (non-surface) bxcell seed. */
							/* If we are sufficiently far from the seed point, */
							/* a super-cell to improve seeding performance will be created. */
							ss = new_bxcell(s, nix, cc, tx->ss, tx->sdist, vflag);
#if defined(REVTABLESTATS) || defined(DEBUG)
							if (tx->scell != NULL)
//								nnsuperfill += tx->scell[3]-3;
								nnsuperfill++;
							else {
								nnsinglefill++;
							}
#endif
//printf("Created temporary seed bxcell %d vflag %x\n",ss->ix, vflag[ss->ix]);
						}
						DBG(("Adding seed neighbor nnrev[%d] vflag %x co %s to seed list\n",ss->ix, vflag[ss->ix], debPiv(s->fdi, ss->gc)));
//printf("Adding seed neighbor nnrev[%d] vflag %x co %s to seed list\n",ss->ix, vflag[ss->ix], debPiv(s->fdi, ss->gc));

						/* Add to end of seedlist */
						ss->flist = NULL;
						if (seedlist == NULL)
							seedlist = ss;
						else
							seedlistend->flist = ss;
						seedlistend = ss;
						vflag[ss->ix] |= 1;		/* This is on seed list, so will be filled */
					}

					cc[f] -= (ff & 1) ? 1 : -1;
					nix   -= (ff & 1) ? s->rev.coi[f] : -s->rev.coi[f];
				} 
			}

			/* if this is a temporary bxcell (i.e. not a surface bxcell), */
			/* we can now free it */
			if ((vflag[tx->ix] & 2) == 0) {
//printf("Done with non-surface bxcell %d vflag %x\n",tx->ix,vflag[tx->ix]);
				del_bxcell(s, tx);
			}
		}
		/* We've done the nnrev[] setup */
		DBG(("rev.nnrev[] grid done - cleaning up\n"));

#ifdef CHECK_NNLU
		if (fdi > 1) {
			/* Check that every nnrev[] cell is filled */
			printf("Checking all %d nnrev[] cells are filled\n",rgno);
			for (i = 0; i < rgno; i++) {
				if (   ((vflag[i] & 2) != 0 || s->rev.rev[i] == NULL || s->rev.rev[i][1] == 3)
				     && (s->rev.nnrev[i] == NULL || s->rev.nnrev[i][1] == 3)) {
					printf("Found empty nnrev[%d] ?:\n",i);
					printf(" vflag %x\n",vflag[i]);
					if (s->rev.nnrev[i] == NULL)
						printf(" nnrev = NULL\n");
					else
						printf(" nnrev length = %d\n",s->rev.nnrev[i][1]-3);
					if (s->rev.rev[i] == NULL)
						printf(" rev = NULL\n");
					else
						printf(" rev = length = %d\n",s->rev.rev[i][1]-3);
				}
			} 
		}
#endif	/* CHECK_NNLU */

		/* Free up flag array used for construction */
		if (vflag != NULL) {
			DECSZ(s, rgno * sizeof(char));
			free(vflag);
		}

#ifndef CHECK_NNLU
		/* Free up surface linked list and delete the bxcells. */
		free_surflist(s);
#endif

		/* Free up surface bxcell hash index */
		free_surfhash(s, 0);
	
#if defined(REVTABLESTATS) || defined(DEBUG)
		if (fdi > 1) {
			nnrevshare = nnrevcells;
			for (i = 0; i < s->rev.sharellen; i++)
				nnrevshare += (s->rev.sharelist[i][1]-4) * (s->rev.sharelist[i][1]-3);

			printf("%d/%d used cells in nnrev list\n",nnrevcells,rgno);
			printf("%f average cells searched\n",(double)nnrevcellsearch/(double)nnrevcells);
			printf("%d max bxcells used\n",maxbxcount);
			printf("%.1f%% super-cell filled\n",100.0 * nnsuperfill/(nnsuperfill+nnsinglefill));
			printf("%f average list length\n",(double)nnrevcelldepth/(double)nnrevcells);
			printf("%d max list length\n",nnmxrevcelldepth);
			printf("%f average shared lists\n",(double)nnrevshare/(double)nnrevcells);
			printf("Took %f seconds\n",0.001 * (msec_time()-smsec));
			printf("Overall took %f seconds\n",0.001 * (msec_time()-smsec));
		}
#endif
	}

	s->rev.rev_valid = 1;

	if (fdi > 1 && s->verbose)
		fprintf(stdout, "%cnnrev initialization done\n",cr_char);

	DBG(("init_revaccell finished\n"));
}

/* Invalidate the reverse acceleration structures (section Two) */
static void invalidate_revaccell(
rspl *s		/* Pointer to rspl grid */
) {
	int e, di = s->di;
	int **rpp, *rp;

	/* Invalidate the whole rev cache (Third section) */
	invalidate_revcache(s->rev.cache);

	/* Free up the contents of rev.rev[] and rev.nnrev[] */
	if (s->rev.rev != NULL) {
		for (rpp = s->rev.rev; rpp < (s->rev.rev + s->rev.no); rpp++) {
			if (*rpp != NULL)
				free_indexlist(s, rpp);
		}
	}
	if (s->rev.nnrev != NULL) {

		/* Free up nn list sharelist records - this will free and set */
		/* any shared lists to NULL */
		free_sharelist(s);

		for (rpp = s->rev.nnrev; rpp < (s->rev.nnrev + s->rev.no); rpp++) {
			if (*rpp != NULL)
				free_indexlist(s, rpp);
		}
	}

	if (di > 1 && s->rev.rev_valid) {
		rev_struct *rsi, **rsp;
		size_t ram_portion = g_avail_ram;

		/* Remove it from the linked list */
		for (rsp = &g_rev_instances; *rsp != NULL; rsp = &((*rsp)->next)) {
			if (*rsp == &s->rev) {
				*rsp = (*rsp)->next;
				break;
			}
		}

		/* Aportion the memory */
		g_no_rev_cache_instances--;

		if (g_no_rev_cache_instances > 0) {
			ram_portion /= g_no_rev_cache_instances; 
			for (rsi = g_rev_instances; rsi != NULL; rsi = rsi->next)
				rsi->max_sz = ram_portion;
			if (s->verbose)
				fprintf(stdout, "%cThere %s %d rev cache instance%s with %lu Mbytes limit\n",
				                cr_char,
								g_no_rev_cache_instances > 1 ? "are" : "is",
			                    g_no_rev_cache_instances,
								g_no_rev_cache_instances > 1 ? "s" : "",
			                    (unsigned long)ram_portion/1000000);
		}
	}
	s->rev.rev_valid = 0;
}

#ifdef CHECK_NNLU
/* ====================================================== */

/* Used exautive searches to check that nn lookup found a good solution */
static void check_nn(
rspl *s,
double *oval,	/* Un-clipped output target value */
co *cpp			/* Clipped output space value in cpp[0].v[] */
				/* nn solution in cpp[0].p[] */
) {
	int i, j;		/* Index of fwd grid point */
	int e, f, ee, ff;
	int di = s->di;
	int fdi = s->fdi;
	int gno = s->g.no;
	int good = 1;
	int found = 0;
	int printed = 0;

	ECOUNT(gc, MXRI, di, 0, s->g.res, 0);/* coordinates */
	float *gp;			/* Pointer to grid data */
	double iv[MXDI];
	double ov[MXDO];
	double chov[MXDO], de;
	
	int bix = -1;
	double bdist = 1e200;
	double biv[MXDI];
	double bov[MXDO];
	int six = -1;
	double sdist = 1e200;
	double siv[MXDI];
	double sov[MXDO];

	double odelta;
	double idelta;
	double fsdelta;
	double sodelta;
	double sidelta;

	s->rev.cknn_no++;

	/* Compute the given solutions de */
	de = sqrt(lchw_sq(s, oval, cpp[0].v));

	/* Go through every fwd vertex looking for closest and 2nd closest */
	EC_INIT(gc);
	for (gp = s->g.a, i = 0; i < gno; gp += s->g.pss, i++) {
		double dist;

		if (s->limiten && gp[-1] > s->limitv) {
			EC_INC(gc);
			continue;			/* Over the ink limit */
		}

		for (f = 0; f < fdi; f++)
			ov[f] = gp[f];

		dist = lchw_sq(s, oval, ov);

		if (dist < bdist) {
			six = bix;
			bix = i;
			for (e = 0; e < s->di; e++) {
				siv[e] = biv[e];
				biv[e] = s->g.l[e] + gc[e] * s->g.w[e];
			}
			for (f = 0; f < fdi; f++) {
				sov[f] = bov[f];
				bov[f] = ov[f];
			}
			sdist = bdist;
			bdist = dist;

		} else if (dist < sdist) {
			six = i;
			for (e = 0; e < s->di; e++)
				siv[e] = s->g.l[e] + gc[e] * s->g.w[e];
			for (f = 0; f < fdi; f++)
				sov[f] = ov[f];
			sdist = dist;

		}
		EC_INC(gc);
	}

	/* What is magnitude of target match ? */
	odelta = sqrt(lchw_sq(s, bov, oval));
	
	/* What is magnitude of solution match */
	idelta = 0.0;
	for (e = 0; e < s->di; e++) {
		double tt = biv[e] - cpp[0].p[e];
		idelta += tt * tt;
	}
	idelta = sqrt(idelta);
	
	/* What is scale of solution from closest to 2nd closest ? */
	fsdelta = 0.0;
	for (e = 0; e < s->di; e++) {
		double tt = biv[e] - siv[e];
		fsdelta += tt * tt;
	}
	fsdelta = sqrt(fsdelta);

	/* What is magnitude of target match to secondary ? */
	sodelta = sqrt(lchw_sq(s, sov, oval));
	
	/* What is magnitude of solution match to secondary ?*/
	sidelta = 0.0;
	for (e = 0; e < s->di; e++) {
		double tt = siv[e] - cpp[0].p[e];
		sidelta += tt * tt;
	}
	sidelta = sqrt(sidelta);
	
	/* If our exaustive search is better than the nn solution: */
	if (odelta < (de - 1e-6)) {
		double dde = de - odelta;
		if (dde > s->rev.cknn_we)
			s->rev.cknn_we = dde;
		s->rev.cknn_noerrs++;
		good = 0;
		printf("check_nn: target %s\n",debPdv(s->fdi,oval));
		printf("check_nn: cliped to %s, de %f\n",debPdv(s->di,cpp[0].v),de);
		printf("check_nn: solution %s\n",debPdv(s->di,cpp[0].p));
		printf("check_nn: check target %s, de %f\n",debPdv(s->fdi, bov),odelta);
		printf("check_nn: check solution %s, de %f @ix %d\n",debPdv(s->di, biv),idelta,bix);
		printf("check_nn: check 2nd target %s, de %f\n",debPdv(s->fdi, sov),sodelta);
		printf("check_nn: check 2nd solution %s, de %f @ ix %d\n",debPdv(s->di, siv),sidelta,six);
		printf("check_nn: excess delta %f\n",dde);
		printf("check_nn: first-second delta %f\n",fsdelta);
		if (six >= 0 && (de - sodelta) > 1e-6) {
			printf("check_nn: beyond 2nd best by %f!\n",de-sodelta);
			s->rev.cknn_nobsb++;
		}
		printed = 1;
	}

	/* Search surface nnrev cells, to make sure our best is in it somewhere */
	if (s->rev.surflist != NULL) {
		bxcell *ss;

		for (ss = s->rev.surflist; ss != NULL; ss = ss->slist) {
			int *flist = ss->sl;		/* List of fwd cells */

			if (flist == NULL)
				error("surflist nnrev[%d] is empty!",ss->ix);

			/* For each forward cell */
			for (flist += 3; *flist != -1; flist++) {
				/* For each cube vertex */
				for (ee = 0; ee < (1<<di); ee++) {
					int vix = *flist + s->g.hi[ee];
					if (vix == bix) {
						found = 1;
						if (!good)
							printf("check_nn: found best vertex in surf nnrev[%d] fwd %d \n",ss->ix,*flist);
						break;
					}
				}
			}
		}

		if (!found) {
			int rgno = s->rev.no;
			int **rpp;
			int revfound = 0;

			s->rev.cknn_nonis++;
			if (good) {
				printf("check_nn: target %s\n",debPdv(s->fdi,oval));
				printf("check_nn: cliped to %s, de %f\n",debPdv(s->di,cpp[0].v),de);
				printf("check_nn: solution %s\n",debPdv(s->di,cpp[0].p));
				printf("check_nn: check target %s, de %f\n",debPdv(s->fdi, bov),odelta);
				printf("check_nn: check solution %s, de %f\n",debPdv(s->di, biv),idelta);
				printf("check_nn: result is OK\n");
			}
			if (s->rev.surflist == NULL) {
				printf("check_nn: No surface list to check against\n");
			} else {
				printf("check_nn: DIDN'T find best vertex %d in nnrev[] surface list\n",bix);
			}
			printed = 1;

			/* See where it is in the rev[] list, and what the corresponding nnrev[] */
			/* looks like */
			for (rpp = s->rev.rev, i = 0; i < rgno; rpp++, i++) {
				int *flist = *rpp;

				if (flist == NULL)
					continue;

				/* For each forward cell */
				for (flist += 3; *flist != -1; flist++) {
					/* For each cube vertex */
					for (ee = 0; ee < (1<<di); ee++) {
						int vix = *flist + s->g.hi[ee];
						if (vix == bix) {
							revfound = 1;
							printf("check_nn: found best vertex in rev[%d] fwd %d",i,*flist);
							if (s->rev.nnrev[i] != NULL)
								printf(" - cspndg. nnrev has list\n");
							else
								printf(" - cspndg. nnrev is empty\n");
							break;
						}
					}
				}
			}
			if (!revfound) {
				printf("check_nn: DIDN'T find best vertex %d in rev list\n",bix);
			}
		}
	}

	/* Check if the nnrev[] cell for this target has the fwd cell */
	if (s->rev.surflist != NULL && (!good || !found)) {
		int mi[MXDO];
		int rgres_1 = s->rev.res - 1;
		int ix, *flist;
		int found2 = 0;

		for (ix = 0, f = 0; f < fdi; f++) {
			double t = (oval[f] - s->rev.gl[f])/s->rev.gw[f];
			mi[f] = (int)floor(t);			/* Grid coordinate */
			if (mi[f] < 0)			/* Clip to reverse range, so we always return a result  */
				mi[f] = 0;
			else if (mi[f] > rgres_1)
				mi[f] = rgres_1;
			ix += mi[f] * s->rev.coi[f];	/* Accumulate reverse grid index */
		}
		flist = s->rev.nnrev[ix];

		if (flist != NULL) {
			/* For each forward cell */
			for (flist += 3; *flist != -1; flist++) {
				/* For each cube vertex */
				for (ee = 0; ee < (1<<di); ee++) {
					int vix = *flist + s->g.hi[ee];
					if (vix == bix) {
						found2 = 1;
						printf("check_nn: found best vertex %d in expected nnrev[%d], fwd %d\n",bix,ix,*flist);
						printed = 1;
						break;
					}
				}
			}
		}
		if (!found2) {
			printf("check_nn: DIDN'T find best vertex %d in expected nnrev[%d] list\n",bix,ix);
			printed = 1;
		}
	}
	if (printed)
		printf("\n");
}

static void print_nnck(rspl *s) {
	printf("check_nn di %d fdi %d checked %d lookups:\n",s->di,s->fdi,s->rev.cknn_no);
	printf("check_nn got %d not as good as best vertex\n",s->rev.cknn_noerrs);
	printf("check_nn got %d not as good as 2nd best vertex\n",s->rev.cknn_nobsb);
	printf("check_nn got %d not in surface list\n",s->rev.cknn_nonis);
	printf("check_nn got %f worst excess de\n",s->rev.cknn_we);
	printf("\n");
}

#endif /* CHECK_NNLU */
/* ====================================================== */

/* Initialise the rev First section, basic information that doesn't change */
/* This is called on initial setup when s->rev.inited == 0 */
static void make_rev_one(
rspl *s
) {
	int i, j;		/* Index of fwd grid point */
	int e, f, ee, ff;
	int di = s->di;
	int fdi = s->fdi;
	int rgno, gno = s->g.no;
	int rgres;		/* bwd cell grid (rev[], nnrev[]) resolution */
	int rgres_1;	/* rgres -1 == maximum base coord value */
	datao rgmin, rgmax;

	DBG(("make_rev_one called, di = %d, fdi = %d, mgres = %d\n",di,fdi,(int)s->g.mres));

//printf("~1 nnb = %d\n",nnb);

	s->get_out_range(s, rgmin, rgmax);	/* overall output min/max */

	/* Expand out range to encompass declared range */
	/* The declared range is assumed to be the range over which */
	/* we may want an reasonably accurate nearest reverse lookup. */
	for (f = 0; f < fdi; f++) {
		if ((s->d.vl[f] + s->d.vw[f]) > rgmax[f])
				rgmax[f] = s->d.vl[f] + s->d.vw[f];
		if (s->d.vl[f] < rgmin[f])
				rgmin[f] = s->d.vl[f];
	}

	/* Expand out range slightly to allow for out of gamut points */
	for (f = 0; f < fdi; f++) {
		double del = (rgmax[f] - rgmin[f]) * 0.10;	/* Expand by +/- 10% */
		rgmax[f] += del;
		rgmin[f] -= del;
	}
//printf("~~got output range\n");

	/* Heuristic - reverse grid acceleration resolution ? */
	/* Should this really be adapted to be constant in output space ? */
	/* (ie. make the gw aprox equal ?) Would complicate code rev accell */
	/* indexing though. */
	{
		char *ev;
		double gresmul = REV_ACC_GRES_MUL;		/* Typically 2.0 */

		if ((gresmul * s->g.mres) > (double)REV_ACC_GRES_LIMIT) {
			gresmul = (double)REV_ACC_GRES_LIMIT/s->g.mres;		/* Limit target res to typ. 43. */
		}

		/* Allow the user to override if it causes memory consumption problems */
		/* or to speed things up if more memory is available */
		if ((ev = getenv("ARGYLL_REV_ACC_GRID_RES_MULT")) != NULL) {
			double mm;
			mm = atof(ev);
			if (mm > 0.1 && mm < 20.0)
				gresmul *= mm;
		}
		/* Less than 4 is not functional */
		if ((rgres = (int) gresmul * s->g.mres) < 4)
			rgres = 4;
	}
	s->rev.res = rgres;			/* == number of cells per side */
	rgres_1 = rgres-1;

	/* Number of elements in the rev.grid */
	for (rgno = 1, f = 0; f < fdi; f++, rgno *= rgres);
	s->rev.no = rgno;

//printf("~1 rgres = %d\n",rgres);
	/* Compute coordinate increments */
	s->rev.coi[0] = 1;
//printf("~1 coi[0] = %d\n",s->rev.coi[0]);
	for (f = 1; f < fdi; f++) {
		s->rev.coi[f] = s->rev.coi[f-1] * rgres;
//printf("~1 coi[%d] = %d\n",f,s->rev.coi[f]);
	}

	/* Compute index offsets from base of cube to other corners. */
	for (s->rev.hoi[0] = f = 0, j = 1; f < fdi; j *= 2, f++) {
		for (i = 0; i < j; i++)
			s->rev.hoi[j+i] = s->rev.hoi[i] + s->rev.coi[f];	/* In grid points */
	}
//for (ff = 0; ff < (1 << fdi); ff++)
//printf("~1 hoi[%d] = %d\n",ff,s->rev.hoi[ff]);

	/* Conversion from output value to cell indexes */
	for (f = 0; f < fdi; f++) {
		s->rev.gl[f] = rgmin[f];
		s->rev.gh[f] = rgmax[f];
		s->rev.gw[f] = (rgmax[f] - rgmin[f])/(double)rgres;
	}

	if ((s->rev.rev = (int **) rev_calloc(s, rgno, sizeof(int *))) == NULL)
		error("rspl malloc failed - rev.grid points");
	INCSZ(s, rgno * sizeof(int *));

	if ((s->rev.nnrev = (int **) rev_calloc(s, rgno, sizeof(int *))) == NULL)
		error("rspl malloc failed - rev.nngrid points");
	INCSZ(s, rgno * sizeof(int *));

	s->rev.inited = 1;
	s->rev.stouch = 1;

	DBG(("make_rev_one finished\n"));
}

/* ====================================================== */

/* First section of rev_struct init. */
/* Initialise the fxcell cache, sub simplex information */
/* and reverse lookup acceleration structures. */
/* This is called by a reverse interpolation call */
/* that discovers that the reverse index list haven't */
/* been initialised. */
static void make_rev(
rspl *s
) {
	int e, di = s->di;
	char *ev;
	size_t avail_ram = 256 * 1024 * 1024;	/* Default assumed RAM in the system */
	size_t ram1, ram2;						/* First Gig and rest */
	static int repsr = 0;					/* Have we reported system RAM size ? */
	size_t max_vmem = 0;

	DBG(("make_rev called, di = %d, fdi = %d, mgres = %d\n",di,s->fdi,(int)s->g.mres));

	/* Figure out how much RAM we can use for the rev cache. */
	/* (We compute this for each rev instance, to account for any VM */
	/* limit changes due to intervening allocations) */
	if (di > 1 || g_avail_ram == 0) {
	#ifdef NT 
		{
			BOOL (WINAPI* pGlobalMemoryStatusEx)(MEMORYSTATUSEX *) = NULL;
			MEMORYSTATUSEX mstat;
	
			pGlobalMemoryStatusEx = (BOOL (WINAPI*)(MEMORYSTATUSEX *))
			                        GetProcAddress(LoadLibrary("KERNEL32"), "GlobalMemoryStatusEx");
	
			if (pGlobalMemoryStatusEx == NULL)
				error("Unable to link to GlobalMemoryStatusEx()");
			mstat.dwLength = sizeof(MEMORYSTATUSEX);
			if ((*pGlobalMemoryStatusEx)(&mstat) != 0) {
				if (sizeof(avail_ram) < 8 && mstat.ullTotalPhys > 0xffffffffL)
					mstat.ullTotalPhys = 0xffffffffL;
				avail_ram = mstat.ullTotalPhys;
			} else {
				warning("%cWarning - Unable to get system memory size",cr_char);
			}
		}
	#else
	#ifdef __APPLE__
		{
			long long memsize;
			size_t memsize_sz = sizeof(long long);
			if (sysctlbyname("hw.memsize", &memsize, &memsize_sz, NULL, 0) == 0) {
				if (sizeof(avail_ram) < 8 && memsize > 0xffffffffL)
					memsize = 0xffffffff;
				avail_ram = memsize;
			} else {
				warning("%cWarning - Unable to get system memory size",cr_char);
			}
			
		}
	#else	/* Linux */
		{
			long long total;
			total = (long long)sysconf(_SC_PAGESIZE) * (long long)sysconf(_SC_PHYS_PAGES);
			if (sizeof(avail_ram) < 8 && total > 0xffffffffL)
				total = 0xffffffffL;
			avail_ram = total;
		}
	#endif
	#endif
		DBG(("System RAM = %d Mbytes\n",avail_ram/1000000));
	
		/* Make it sane */
		if (avail_ram < (256 * 1024 * 1024)) {
			warning("%cWarning - System RAM size seems very small (%d MBytes),"
			        " assuming 256Mb instead",cr_char,avail_ram/1000000);
			avail_ram = 256 * 1024 * 1024;
		}
		// avail_ram = -1;		/* Fake 4GB of RAM. This will swap! */
	
		ram1 = avail_ram;
		ram2 = 0;
		if (ram1 > (1024 * 1024 * 1024)) {
			ram1 = 1024 * 1024 * 1024;
			ram2 = avail_ram - ram1;
		}
	
		/* Default maximum reverse memory (typically 50% of the first Gig, 75% of the rest) */
		g_avail_ram = (size_t)(REV_MAX_MEM_RATIO * ram1
		            +          REV_MAX_MEM_RATIO2 * ram2);
	
		/* Many 32 bit systems have a virtual memory limit, so we'd better stay under it. */
		/* This is slightly dodgy though, since we don't know how much memory other */
		/* software will need to malloc. A more sophisticated approach would be to */
		/* replace all malloc/calloc/realloc calls in the exe with a version that on failure, */
		/* sets the current memory usage as the new limit, and then */
		/* frees up some rev cache space before re-trying. This is a non-trivial change */
		/* to the source code though, and really has to include all user mode */
		/* libraries we're linked to, making implementation problematic. */ 
		/* Instead we do a simple test to see what the maximum allocation is, and */
		/* then use 75% of that for cache, and free cache and retry if */
		/* malloc failes in rev.c. Too bad if 25% isn't enough, and a malloc fails */
		/* outside rev.c... */
		if (sizeof(avail_ram) < 8) {
			char *alocs[4 * 1024];
			size_t safe_max_vmem = 0;
			int i; 
	
#ifdef __APPLE__
			int old_stderr, new_stderr;

			/* OS X malloc() blabs about a malloc failure. This */
			/* will confuse users, so we temporarily redirect stdout */
			fflush(stderr);
			old_stderr = dup(fileno(stderr));
			new_stderr = open("/dev/null", O_WRONLY | O_APPEND);
			dup2(new_stderr, fileno(stderr));
#endif
			for (i = 0; (i < 4 * 1024);i++) {
				if ((alocs[i] = malloc(1024 * 1024)) == NULL) {
					break;
				}
				max_vmem = (i+1) * 1024 * 1024;
			}
			for (--i; i >= 0; i--) {
				free(alocs[i]);
			}
#ifdef __APPLE__
			fflush(stderr);
			dup2(old_stderr, fileno(stderr));	/* Restore stderr */
			close(new_stderr);
			close(old_stderr);
#endif
			/* To compute a true value, we need to allow for any VM already */
			/* used by any rev instances. */
			{
				rev_struct *rsi;

				for (rsi = g_rev_instances; rsi != NULL; rsi = rsi->next)
					max_vmem += rsi->sz;
			}
			
//fprintf(stdout,"~ Abs max VM = %d Mbytes\n",max_vmem/1000000);
			safe_max_vmem = (size_t)(0.85 * max_vmem);
			if (g_avail_ram > safe_max_vmem) {
				g_avail_ram = safe_max_vmem;
				if (s->verbose && repsr == 0)
					fprintf(stdout,"%cTrimmed maximum cache RAM to %lu Mbytes to allow for VM limit\n",cr_char,(unsigned long)g_avail_ram/1000000);
			}
		}
	
		/* Check for environment variable tweak  */
		if ((ev = getenv("ARGYLL_REV_CACHE_MULT")) != NULL) {
			double mm, gg;
			mm = atof(ev);
			if (mm < 0.01)			/* Make it sane */
				mm = 0.01;
			else if (mm > 100.0)
				mm = 100.0;
			gg = g_avail_ram * mm + 0.5;
			if (gg > (double)(((size_t)0)-1))
				gg  = (double)(((size_t)0)-1);
			g_avail_ram = (size_t)(gg);
		}
		if (max_vmem != 0 && g_avail_ram > max_vmem && repsr == 0) {
			g_avail_ram = (size_t)(0.95 * max_vmem);
			fprintf(stdout,"%cARGYLL_REV_CACHE_MULT * RAM trimmed to %lu Mbytes to allow for VM limit\n",cr_char,(unsigned long)g_avail_ram/1000000);
		}
	}

	/* Default - this will get aportioned as more instances appear */
	s->rev.max_sz = g_avail_ram;

	DBG(("reverse cache max memory = %d Mbytes\n",s->rev.max_sz/1000000));
	if (s->verbose && repsr == 0) {
		fprintf(stdout, "%cRev cache RAM = %lu Mbytes\n",cr_char,(unsigned long)g_avail_ram/1000000);
		repsr = 1;
	}

	/* Sub-simplex information for each sub dimension */
	for (e = 0; e <= di; e++) {
		if (s->rev.sspxi[e].spxi != NULL)	/* Assert */
			error("rspl rev, internal, init_ssimplex_info called on already init'd\n");

		rspl_init_ssimplex_info(s, &s->rev.sspxi[e], e);
	}

	make_rev_one(s);

	/* Reverse cell cache allocation */
	s->rev.cache = alloc_revcache(s);

	DBG(("make_rev finished\n"));
}

/* ====================================================== */

#if defined(DEBUG1) || defined(DEBUG2)

/* Utility - return a string containing a fwd cells output value range */
static char *pcellorange(fxcell *c) {
	static char buf[5][300];
	static ix = 0;
	char *bp;
	rspl *s = c->s;
	int di = s->di, fdi = s->fdi;
	int ee, e, f;
	
	datao min, max;

//	double p[POW2MXRI][MXRI]; /* Vertex input positions for this cube. */
//	double v[POW2MXRI][MXRO+1]; /* Vertex data for this cube. Copied to x->v[] */
//							/* v[][fdi] is the ink limit values, if relevant */

	for (f = 0; f < fdi; f++) {
		min[f] = 1e60;
		max[f] = -1e60; 
	}

	/* For all other grid points in the cube */
	for (ee = 0; ee < (1 << di); ee++) {
		
		/* Update bounding box for this grid point */
		for (f = 0; f < fdi; f++) {
			if (min[f] > c->v[ee][f])	
				 min[f] = c->v[ee][f];
			if (max[f] < c->v[ee][f])
				 max[f] = c->v[ee][f];
		}
	}
	if (++ix >= 5)
		ix = 0;
	bp = buf[ix];

	for (e = 0; e < fdi; e++) {
		if (e > 0)
			*bp++ = ' ';
		sprintf(bp, "%f:%f", min[e],max[e]); bp += strlen(bp);
	}
	return buf[ix];
}

#endif
/* ====================================================== */

#undef DEBUG
#undef DBGV
#undef DBG
#define DBGV(xxx)
#define DBG(xxx) 

#ifdef REVVRML
/* ====================================================== */
/* VRML diagnostic output functions */

/* Plot the initial surface rev cells */
static void plot_bxfwcells(
rspl *s,
int dobxcells,			/* Plot rev cells */
int dofwcells,			/* Plot fwd cells */
int dofwlabels			/* Plot fwd cell base indexs */ 
) {
	int i, j;		/* Index of fwd grid point */
	int e, f, ee, ff;
	int di = s->di;
	int fdi = s->fdi;
	bxcell *bx;
	vrml *wrl;
	double grey[3] = { 0.5, 0.5, 0.5 };
	double white[3] = { 1.0, 1.0, 1.0 };

	wrl = new_vrml("raw_bxfwcells", 0, s->rev.probxyz ? vrml_xyz : vrml_lab);
	wrl->add_marker(wrl, s->rev.ocent, NULL, 1.0);

	if (dofwlabels) {
		/* Put text for every base cube index */
		for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
			int vix[POW2MXRI];
			int *crp, *rp;
		
			crp = s->rev.rev[bx->ix];
		
			for (rp = crp+3; *rp != -1; rp++) {
				int ix = *rp;
				char index[100];
				double vv[MXRI];
				int off = 0;		// 0 .. 7, choose cube vertex
				float *fcb = s->g.a + (ix + s->g.hi[off]) * s->g.pss;
		
				for (e = 0; e < di; e++)
					vv[e] = fcb[e];
				sprintf(index, "%d",ix);
				wrl->add_text(wrl, index, vv, white, 0.3);
			}
		}
	}

	if (dobxcells) {
		for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
			int vix[POW2MXRO];
			DCOUNT(cc, MXRO, fdi, 0, 0, 2);		/* Vertex counter */
			int *crp, *rp;
		
			/* Plot bxcell's */
			i = 0;
			DC_INIT(cc);
			while (!DC_DONE(cc)) {
				double vv[MXRO];
				for (f = 0; f < fdi; f++)
					vv[f] = (bx->gc[f] + cc[f]) * s->rev.gw[f] + s->rev.gl[f];  
				vix[i] = wrl->add_vertex(wrl, 0, vv);
				DC_INC(cc);
				i++;
			}

			/* For each vertex */
			for (i = 0; i < (1 << fdi); i++) {
				int lix[2];

				lix[0] = vix[i];

				/* for each dimension */ 
				for (j = 0; j < fdi; j++) {
					if (i & (1<<j))
						continue;		/* Would go outside cube */

					lix[1] = vix[i | (1 << j)];
					if (dofwcells)
						wrl->add_col_line(wrl, 0, lix, grey);
					else
						wrl->add_line(wrl, 0, lix);
				}
			}
		}
	}

	if (dofwcells) {
		for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
			int vix[POW2MXRI];
			int *crp, *rp;
		
			/* Add fwd cells */
			crp = s->rev.rev[bx->ix];
			for (rp = crp+3; *rp != -1; rp++) {
				float *fcb = s->g.a + *rp * s->g.pss;

				/* Skip grid base points on the upper edge of the grid */
				for (e = 0; e < di; e++) {
					if (G_FL(fcb, e) == 0)		/* At the top edge */
						break;
				}
				if (e < di) {
					printf("Fwd cell base index %d is on upper edge!\n",*rp);
					continue;
				}

				/* For each vertex of cube */
				for (i = 0; i < (1<<di); i++) {
					double vv[MXRI];
					int ix = *rp + s->g.hi[i];
					fcb = s->g.a + ix * s->g.pss;

					if (!s->limiten || fcb[-1] <= s->limitv)
						break;
				}
				/* Skip any cubes that a completely over the ink limit */
				if (i >= (1<<di))
					continue;

				/* For each vertex of cube */
				for (i = 0; i < (1<<di); i++) {
					double vv[MXRI];
					int ix = *rp + s->g.hi[i];
					fcb = s->g.a + ix * s->g.pss;

					for (e = 0; e < di; e++)
						vv[e] = fcb[e];
					vix[i] = wrl->add_vertex(wrl, 1, vv);
				}
					
				/* For each vertex of cube */
				for (i = 0; i < (1<<di); i++) {
					int lix[2];
		
					lix[0] = vix[i];
		
					/* for each dimension */ 
					for (j = 0; j < di; j++) {
						if (i & (1<<j))
							continue;		/* Would go outside cube */
		
						lix[1] = vix[i | (1 << j)];
						wrl->add_line(wrl, 1, lix);
					}
				}
			}
		}
	}
	wrl->make_lines_vc(wrl, 0, 0.0);
	wrl->make_lines_vc(wrl, 1, 0.0);

	printf("Created %s\n",wrl->name);
	wrl->del(wrl);
}

/* Plot vertex & triangle check setup & solution */
/* + the primary and shadow bxcells. */ 
static void plot_tri_check(
rspl *s,
int dobxcells,			/* Plot prim & shadow bxcell cells */
int dowait,				/* Wait for the user to hit return */
bxcell *bx,				/* First bx cell (if dobxcells set) */
int vtxix,				/* triangle base vertex index (-1 if not applicable) */
int trii,				/* Triangle eneration */
int triix[3],			/* Triangle indexes */
int nvtxix,				/* test point vertex index number (may be -1 if not vtxrec) */
int sorv,				/* Intersection was solved ? */			
int wsrv,				/* Within simplex ? */
int shdwd,				/* Vertex is shadowed ? */
double v[MXRI+1][MXRO],	/* Triangle vertex values */
double de[MXRO],		/* Line delta */
double pv[MXRO],		/* Vertex being tested */
double xv[MXRO]			/* Intersection point */
) {
	int j;
	int e, f, ee, ff;
	int di = s->di;
	int fdi = s->fdi;
	vrml *wrl;
	bxcell *vbx;
	int first = 1;
	int ii, vix[POW2MXRO], lix[3];
	double vv[MXRO];
	double white[3] = { 1.0, 1.0, 1.0 };
	double grey[3] = { 0.5, 0.5, 0.5 };
	double green[3] = { 0.1, 1.0, 0.1 };
	double red[3] = { 0.8, 0.1, 0.1 };
	double blue[3] = { 0.1, 0.1, 0.8 };
	double yellow[3] = { 0.8, 0.8, 0.1 };

	wrl = new_vrml("tri_check", 0, s->rev.probxyz ? vrml_xyz : vrml_lab);

	/* Gamut center point marker */
	wrl->add_marker(wrl, s->rev.ocent, NULL, 1.0);

	/* point being tested marker */
	wrl->add_marker(wrl, pv, shdwd ? red : blue, 0.5);

	/* Intersection point */
	if (wsrv)
		wrl->add_marker(wrl, xv, blue, 0.2);

	/* Line from center through point being tested */
	lix[0] = wrl->add_vertex(wrl, 0, s->rev.ocent);
	for (ii = 0; ii < fdi; ii++)
		vv[ii] = s->rev.ocent[ii] + 10.0 * de[ii];
	lix[1] = wrl->add_vertex(wrl, 0, vv);
	wrl->add_col_line(wrl, 0, lix, grey);

	/* Triangle */
	lix[0] = wrl->add_vertex(wrl, 1, v[0]);
	lix[1] = wrl->add_vertex(wrl, 1, v[1]);
	lix[2] = wrl->add_vertex(wrl, 1, v[2]);
	wrl->add_col_triangle(wrl, 1, lix, green);
	/* And again to get both faces */
	lix[0] = wrl->add_vertex(wrl, 1, v[0]);
	lix[1] = wrl->add_vertex(wrl, 1, v[2]);
	lix[2] = wrl->add_vertex(wrl, 1, v[1]);
	wrl->add_col_triangle(wrl, 1, lix, green);

	if (dobxcells) {
//printf(" bx = %p\n",bx);
		for (vbx = bx; vbx != NULL; vbx = vbx->wlist) {
			DCOUNT(cc, MXRO, fdi, 0, 0, 2);		/* Vertex counter */
			int *crp, *rp;
		
//printf(" vrml adding bxcell %d\n",vbx->ix);
			/* Plot bxcell's */
			ii = 0;
			DC_INIT(cc);
			while (!DC_DONE(cc)) {
				for (f = 0; f < fdi; f++)
					vv[f] = (vbx->gc[f] + cc[f]) * s->rev.gw[f] + s->rev.gl[f];  
//printf("  vrml vtx %d from %s\n",vix[i], debPdv(3,vv));
				vix[ii] = wrl->add_vertex(wrl, 0, vv);
				DC_INC(cc);
				ii++;
			}

			/* For each vertex */
			for (ii = 0; ii < (1 << fdi); ii++) {

				lix[0] = vix[ii];

				/* for each dimension */ 
				for (j = 0; j < fdi; j++) {
					if (ii & (1<<j))
						continue;		/* Would go outside cube */

					lix[1] = vix[ii | (1 << j)];
//printf("  vrml line from vtx %d - %d\n",lix[0],lix[1]);
					wrl->add_col_line(wrl, 0, lix, first ? white : red);
				}
			}
			first = 0;
		}
	}

	wrl->make_lines_vc(wrl, 0, 0.0);
	wrl->make_triangles(wrl, 1, 0.0, NULL);
	printf("Created %s\n",wrl->name);
	wrl->del(wrl);

	printf(" Solved %s, Within triang %s, shadowed %s\n", sorv ? "true" : "false", wsrv ? "true" : "false", shdwd ? "true" : "false");
	printf("Testing against tri %d %d %d\n", triix[0], triix[1], triix[2]);

	printf(" bx %d vtx %d tri %d checking nvx %d, hit return key:\n",bx->ix, vtxix, trii, nvtxix);
	if (dowait) {
		printf(" hit return key to continue:\n");
		getchar();
	}
}

/* Main summary plot at each thinning round and at end. */
/* Show vertex surface & optional added or deleted vertexes, */
/* + optional bxcells. */
#define VV vv		/* Actual surface values */
//#define VV vl		/* Logf mapped surface values */
static void plot_vtx_surface(
rspl *s,
int dovtxlabels,		/* Show vertex index numbers */
int dodeleted,			/* Show deleted vertexes */
int doadded,			/* Show added vertexes */
int dopres,				/* Show preserved vertexes */
int dooil,				/* Show over ink limit vertexes */
int dobxcells,			/* Show bxcells */
int dowait,				/* Wait for a return key */
vtxcache *vc,			/* Vertexes */
assdire *edgdir		/* Edge lookup for vertex */
) {
	vtxrec *vx, *nvx;
	int i, j;
	int f, fdi = s->fdi;
	vrml *wrl;
	double grey[3] = { 0.5, 0.5, 0.5 };
	double red[3] = { 0.8, 0.1, 0.1 };
	double green[3] = { 0.2, 0.8, 0.2 };
	double blue[3] = { 0.2, 0.2, 0.8 };
	double white[3] = { 0.8, 0.8, 0.8 };
	double magenta[3] = { 0.8, 0.2, 0.8 };
	double cyan[3] = { 0.0, 1.0, 1.0 };
	double yellow[3] = { 1.0, 1.0, 0.0 };
	bxcell *vbx;

	if (dopres)
		wrl = new_vrml("last_surface", 0, s->rev.probxyz ? vrml_xyz : vrml_lab);
	else
		wrl = new_vrml("thinned_surface", 0, s->rev.probxyz ? vrml_xyz : vrml_lab);
	wrl->add_marker(wrl, s->rev.ocent, NULL, 1.0);

	if (dovtxlabels) {
		for (i = 0; i < vc->hash_size; i++) {
			for (vx = vc->hash[i]; vx != NULL; vx = vx->hlink) {
				char index[100];

				if (vx->status == vtx_norm
				 || (dodeleted && (vx->status == vtx_sha || vx->status == vtx_del))
				 || (doadded && vx->addvtx)
				 || (dopres && vx->pres)
				 || (dooil && vx->status == vtx_oil)) {
					sprintf(index, "%d",vx->ix);
					wrl->add_text(wrl, index, vx->VV, cyan, 0.3);
				}
			}
		}
	}

	/* Go through the vertex hash to set every vertex value */
	for (i = 0; i < vc->hash_size; i++) {
		for (vx = vc->hash[i]; vx != NULL; vx = vx->hlink) {
			double *col = NULL;

			if (vx->status != vtx_norm && vx->addvtx)
				error ("Found vertex that is both deleted and cause of added bxcell");

			if (doadded && vx->addvtx)				/* Cause of added bxcell */
				col = green;
			else if (dopres && vx->pres)			/* Preserved vertex */
				col = yellow;
			else if (dodeleted && (vx->status == vtx_sha || vx->status == vtx_del))
				col = red;
			else if (dooil && vx->status == vtx_oil)
				col = blue;
			else if (vx->status == vtx_norm)
				col = white;

			if (col != NULL) {
				vx->vrmlix = wrl->add_col_vertex(wrl, 0, vx->VV, col);
			}
		}
	}

	/* Go through them again to get every line they are part of */
	for (i = 0; i < vc->hash_size; i++) {
		for (vx = vc->hash[i]; vx != NULL; vx = vx->hlink) {
			assdire *edg;			/* Edge table */
			float *fp;
			int fl;
			int pline = 0;			/* Plotted at least 1 line */
			int lix[2];

			fp = s->g.a + vx->ix * s->g.pss;	/* This vertex in fwd grid */
			fl = FLV(fp);		/* Edge flags for this vertex */
			edg = edgdir + fl;

			/* For all possible edges that use this vertex */
			for (j = 0; j < edgdir[fl].no; j++) {
				int fix;
				int eix;

				/* Index of first vertex of the line */
				fix = vx->ix + edg->ti[j].goffs[0];

				/* Index number of vertex other than the one we got it from */
				if (edg->ti[j].goffs[0] != 0)
					eix = vx->ix + edg->ti[j].goffs[0];
				else
					eix = vx->ix + edg->ti[j].goffs[1];

				if ((nvx = get_vtxrec(vc, eix)) != NULL) {
					if ( (vx->status == vtx_norm
					  || (dodeleted && (vx->status == vtx_sha || vx->status == vtx_del))
					  || (doadded && vx->addvtx)
					  || (dopres && vx->pres)
					  || (dooil && vx->status == vtx_oil))
					 &&  (nvx->status == vtx_norm
					  || (dodeleted && (nvx->status == vtx_sha || nvx->status == vtx_del))
					  || (doadded && nvx->addvtx)
					  || (dopres && nvx->pres)
					  || (dooil && nvx->status == vtx_oil))) {

						pline = 1; /* Will/would plot this */

						/* Only plot the line once though */
						if (fix == vx->ix) {
							lix[0] = vx->vrmlix; 
							lix[1] = nvx->vrmlix;
							wrl->add_line(wrl, 0, lix);
						}
					}
				}
			}

			/* we have an orphan vertex */
			if (pline == 0
			 && (dodeleted || vx->status == vtx_norm)
			 && (doadded || !vx->addvtx)) {
				double vv[MXRO], off = 0.15, *col;

				if (doadded && vx->addvtx)				/* Cause of added bxcell */
					col = green;
				else if (dopres && vx->pres)			/* Preserved vertex */
					col = yellow;
				else if (dodeleted && vx->status != vtx_norm)
					col = red;
				else if (dooil && vx->status == vtx_oil)
					col = blue;
				else if (vx->status == vtx_norm)
					col = white;

				for (f = 0; f < fdi; f++)
					vv[f] = vx->VV[f] + off;
				lix[0] = wrl->add_vertex(wrl, 2, vv);
				for (f = 0; f < fdi; f++)
					vv[f] = vx->VV[f] - off;
				lix[1] = wrl->add_vertex(wrl, 2, vv);
				wrl->add_col_line(wrl, 2, lix, col);

				for (f = 0; f < fdi; f++)
					vv[f] = vx->VV[f] + ((f & 1) ? off : -off);
				lix[0] = wrl->add_vertex(wrl, 2, vv);
				for (f = 0; f < fdi; f++)
					vv[f] = vx->VV[f] - ((f & 1) ? off : -off);
				lix[1] = wrl->add_vertex(wrl, 2, vv);
				wrl->add_col_line(wrl, 2, lix, col);

				for (f = 0; f < fdi; f++)
					vv[f] = vx->VV[f] + ((f & 2) ? off : -off);
				lix[0] = wrl->add_vertex(wrl, 2, vv);
				for (f = 0; f < fdi; f++)
					vv[f] = vx->VV[f] - ((f & 2) ? off : -off);
				lix[1] = wrl->add_vertex(wrl, 2, vv);
				wrl->add_col_line(wrl, 2, lix, col);
			}
		}
	}
	wrl->make_lines_vc(wrl, 0, 0.0);
	wrl->make_lines_vc(wrl, 2, 0.0);

	/* Plot surface cells */
	if (dobxcells) {
		for (vbx = s->rev.surflist; vbx != NULL; vbx = vbx->slist) {
			DCOUNT(cc, MXRO, fdi, 0, 0, 2);		/* Vertex counter */
			int *crp, *rp, ii;
			double vv[MXRO];
			int vix[POW2MXRO], lix[2];
		
			ii = 0;
			DC_INIT(cc);
			while (!DC_DONE(cc)) {
				for (f = 0; f < fdi; f++) {
					vv[f] = (vbx->gc[f] + cc[f]) * s->rev.gw[f] + s->rev.gl[f];  
					vv[f] += d_rand(-0.05, 0.05);
				}
				vix[ii] = wrl->add_vertex(wrl, 1, vv);
				DC_INC(cc);
				ii++;
			}

			for (ii = 0; ii < (1 << fdi); ii++) {

				lix[0] = vix[ii];

				/* for each dimension */ 
				for (j = 0; j < fdi; j++) {
					if (ii & (1<<j))
						continue;		/* Would go outside cube */

					lix[1] = vix[ii | (1 << j)];
					if (vbx->debug) {		/* Added bxcell */
						wrl->add_col_line(wrl, 1, lix, magenta);
					} else {				/* Existing bxcell */
						wrl->add_col_line(wrl, 1, lix, grey);
					}
				}
			}
		}
		wrl->make_lines_vc(wrl, 1, 0.0);
	}

	printf("Created %s\n",wrl->name);
	wrl->del(wrl);
	if (dowait) {
		printf(" Thinned vertexes surface: Hit return to continue\n");
		getchar();
	}
}

/* Plot bxcells touched by added cell */
static void plot_touched_bxcells(
rspl *s,
int bxix			/* Index of bx cell causing touches */
) {
	int j, f, fdi = s->fdi;
	vrml *wrl;
	bxcell *vbx;
	int first = 1;
	int ii, vix[POW2MXRO], lix[3];
	double vv[MXRO];
	double green[3] = { 0.1, 0.6, 0.1 };
	double white[3] = { 1.0, 1.0, 1.0 };
	double red[3] = { 0.8, 0.1, 0.1 };

	wrl = new_vrml("add_touch_bxcells", 0, s->rev.probxyz ? vrml_xyz : vrml_lab);

	/* Gamut center point marker */
	wrl->add_marker(wrl, s->rev.ocent, NULL, 1.0);

	for (vbx = s->rev.surflist; vbx != NULL; vbx = vbx->slist) {
		DCOUNT(cc, MXRO, fdi, 0, 0, 2);		/* Vertex counter */
		int *crp, *rp;
	
		/* Plot bxcell's */
		ii = 0;
		DC_INIT(cc);
		while (!DC_DONE(cc)) {
			for (f = 0; f < fdi; f++) {
				vv[f] = (vbx->gc[f] + cc[f]) * s->rev.gw[f] + s->rev.gl[f];  
				if (vbx->debug == 2)
					vv[f] += 0.05;
				else if (vbx->debug == 1)
				vv[f] -= 0.05;
			}
			vix[ii] = wrl->add_vertex(wrl, 0, vv);
			DC_INC(cc);
			ii++;
		}

		/* For each vertex */
		for (ii = 0; ii < (1 << fdi); ii++) {

			lix[0] = vix[ii];

			/* for each dimension */ 
			for (j = 0; j < fdi; j++) {
				if (ii & (1<<j))
					continue;		/* Would go outside cube */

				lix[1] = vix[ii | (1 << j)];
				wrl->add_col_line(wrl, 0, lix,
				           vbx->debug == 2 ? white :  vbx->debug == 1 ? red : green);
			}
		}
	}

	wrl->make_lines_vc(wrl, 0, 0.0);
	printf("Created %s\n",wrl->name);
	wrl->del(wrl);

	printf(" Touched bx cells for bx %d: Hit return to continue\n",bxix);
	getchar();
}

/* Plot the thinned surface fwd cells */
static void plot_fxcell_surface(
rspl *s,
int dofclabels,			/* Show fwd cell base indexes */
int dobxcells,			/* Show bxcells */
int dowait				/* Wait for a return key */
) {
	bxcell *bx;
	int i, j;
	int e, di = s->di;
	int f, fdi = s->fdi;
	vrml *wrl;
	double grey[3] = { 0.5, 0.5, 0.5 };
	double white[3] = { 1.0, 1.0, 1.0 };

	wrl = new_vrml("thinned_fwcells", 0, s->rev.probxyz ? vrml_xyz : vrml_lab);
	wrl->add_marker(wrl, s->rev.ocent, NULL, 1.0);

	if (dofclabels) {
		/* Put text for every base cube index */
		for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
			int vix[POW2MXRI];
			int *crp, *rp;

			crp = bx->sl;

			for (rp = crp+3; *rp != -1; rp++) {
				int ix = *rp;
				char index[100];
				double vv[MXRI];
				int off = 0;		// 0 .. 7, choose cube vertex
				float *fcb = s->g.a + (ix + s->g.hi[off]) * s->g.pss;

				for (e = 0; e < di; e++)
					vv[e] = fcb[e];
				sprintf(index, "%d",ix + s->g.hi[off]);
				wrl->add_text(wrl, index, vv, white, 0.3);
			}
		}
	}

	for (bx = s->rev.surflist; bx != NULL; bx = bx->slist) {
		DCOUNT(cc, MXRO, fdi, 0, 0, 2);		/* Vertex counter */
		int vix[POW2MXRI];
		int *crp, *rp;

		if (dobxcells) {
			/* Plot bxcell's */
			i = 0;
			DC_INIT(cc);
			while (!DC_DONE(cc)) {
				double vv[MXRO];
				for (f = 0; f < fdi; f++)
					vv[f] = (bx->gc[f] + cc[f]) * s->rev.gw[f] + s->rev.gl[f];  
				vix[i] = wrl->add_vertex(wrl, 1, vv);
				DC_INC(cc);
				i++;
			}

			/* For each vertex */
			for (i = 0; i < (1 << fdi); i++) {
				int lix[2];

				lix[0] = vix[i];

				/* for each dimension */ 
				for (j = 0; j < fdi; j++) {
					if (i & (1<<j))
						continue;		/* Would go outside cube */

					lix[1] = vix[i | (1 << j)];
					wrl->add_col_line(wrl, 1, lix, white);
				}
			}
		}

		crp = bx->sl;

		for (rp = crp+3; *rp != -1; rp++) {
			float *fcb = s->g.a + *rp * s->g.pss;

			/* Skip grid base points on the upper edge of the grid */
			for (e = 0; e < di; e++) {
				if (G_FL(fcb, e) == 0)		/* At the top edge */
					break;
			}
			if (e < di) {
				printf("Fwd cell base index %d is on upper edge!\n",*rp);
				continue;
			}

			/* For each vertex of cube */
			for (i = 0; i < (1<<di); i++) {
				double vv[MXRI];
				int ix = *rp + s->g.hi[i];
				float *fcb = s->g.a + ix * s->g.pss;

				if (!s->limiten || fcb[-1] <= s->limitv)
					break;
			}
			/* Skip any cubes that a completely over the ink limit */
			if (i >= (1<<di))
				continue;

			/* For each vertex of cube */
			for (i = 0; i < (1<<di); i++) {
				double vv[MXRI];
				int ix = *rp + s->g.hi[i];
				float *fcb = s->g.a + ix * s->g.pss;

				for (e = 0; e < di; e++)
					vv[e] = fcb[e];
				vix[i] = wrl->add_vertex(wrl, 0, vv);
			}
				
			/* For each vertex of cube */
			for (i = 0; i < (1<<di); i++) {
				int lix[2];

				lix[0] = vix[i];

				/* for each dimension */ 
				for (j = 0; j < di; j++) {
					if (i & (1<<j))
						continue;		/* Would go outside cube */

					lix[1] = vix[i | (1<<j)];
					wrl->add_line(wrl, 0, lix);
				}
			}
		}
	}
	if (dobxcells)
		wrl->make_lines_vc(wrl, 1, 0.0);
	wrl->make_lines_vc(wrl, 0, 0.0);
	printf("Created %s\n",wrl->name);
	wrl->del(wrl);

	if (dowait) {
		printf(" Thinned fwd cell surface: Hit return to continue\n");
		getchar();
	}
}

/* ====================================================== */
#endif /* REVVRML */



