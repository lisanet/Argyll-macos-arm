#ifndef RSPL_REV_H
#define RSPL_REV_H

/* 
 * Argyll Color Management System
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

#undef STATS		/* Collect and print reverse cach stats */

/* Data structures used by reverse lookup code. */
/* Note that the reverse lookup code only supports a more limited */
/* dimension range than the general rspl code. */

/*
 * Note on simplex parameter space. 
 *
 * Simplex interpolation is normaly done in Baricentric space, where there
 * is one more baricentric coordinate than dimensions, and the sum of all
 * the baricentric coordinates must be 1. 
 *
 * To simplify things, we work in a "Simplex parameter" space, in which
 * there are only dimension parameters, and each directly corresponds
 * with a cartesian coordinate, but the parameter order corresponds with
 * the baricentric order.
 *
 * For example, given cartesian sub-coordinates D0, D1, D2
 * into a (3D) forward interpolation cube, these should be sorted
 * from smallest to largest, thereby choosing a particular
 * simplex within a cube, and allowing a correspondence to
 * the parameter coordinates, ie:
 *
 * D2 D0 D1		Smallest -> Largest cartesian sort
 * P0 P1 P2 	Corresponding Parameter coordinates
 * 
 * B0 = P0		Conversion to Baricentric weighting/coordinates
 * B1 = P1 - P0
 * B2 = P2 - P1
 * B3 = 1  - P2
 *
 * The 4 (tetrahedron) vertex values directly correspond to Baricentric
 * weighting/coordinates, giving the usual interpolation equation of:
 *
 *   VV0 * B0
 * + VV1 * B1
 * + VV2 * B2
 * + VV3 * B3
 *
 * Reversing the Parameter -> Baricentric equations gives the
 * following interpolation equation using Parameter coordinates:
 *
 *   (VV0 - VV1) * P0
 * + (VV1 - VV2) * P1
 * + (VV2 - VV3) * P2
 * +  VV3
 *
 * Note that withing the simplex, 0 <= P0 && P0 <= P1 && P1 <= P2 && P2 <= 1
 *
 * It is this which is used in rev.c for solving the reverse
 * interpolation problem.
 */

/* - - - - - - - - - - - - - - - - - - - - - */
/* Group size information for nn LCh weighted quick accept/reject testing */
typedef struct {
	double bcent[MXRO];	/* Group center location in output space */ 
	double brad;		/* Output value bounding shere radius */
	double bradsq;		/* Output value bounding shere radius squared */

	double maxDlc;		/* Maximum members weighted delta LC (+extra dims) squared */
	double maxDh;		/* Maximum members delta h squared */
	double maxDh_;		/* Maximum members delta h (not squared) */
	double sratio;		/* Minimum members C ratio to Gc squared */
	double bratio;		/* Maximum members C ratio to Gc squared */
	double Wsratio;		/* Minimum members C ratio to Gc squared - pre-weighted */
	double Wbratio;		/* Maximum members C ratio to Gc squared - pre-weighted */
	double Gc;			/* Group center Chrominance squared */
	double Gc_;			/* Group center Chrominance (not squared) */
} nn_grp;

/* - - - - - - - - - - - - - - - - - - - - - */
/* A structure to hold per simplex coordinate and vertex mapping for ssxinfo. */
/* This is relative to the construction cube. A face sub-simplex */
/* that is common between cubes, will have a different psxinfo */
/* depending on which cube created it. */
typedef struct {
	int face;			/* Flag, nz if simplex lies on cube surface */
	int icomb[MXDI];    /* icomb[] specifies the equation to convert simplex space */
						/* coordinates back into cartesian space. */
	                    /* Index by Absolute[di] -> Simplex Parameter[sdi], */
						/*                          -1 == value 0, -2 == value 1 */
						/* Note that many Abs can map to one Param to form a sum. */
	int offs[MXDI+1];	/* Offsets to simplex verticies within cube == bit per dim */
	int goffs[MXDI+1];	/* Offsets to simplex verticies within grid */
	int foffs[MXDI+1];	/* Fwd grid floating offsets to simplex verticies from cube base */
	int pmino[MXDI], pmaxo[MXDI]; /* Cube verticy offsets to setup simplex pmin[] and */
						/* pmax[] bounding box pointers. */
} psxinfo;

/* Sub simplexes of a cube information structure */
typedef struct {
	int sdi;			/* Sub-simplex dimensionality */
	int nospx;			/* Number of sub-simplexs per cube */
	psxinfo *spxi;		/* Per sub-simplex info array, NULL if not initialised */
} ssxinfo;

/* - - - - - - - - - - - - - - - - - - - - - */
/* NOTE :- This should really be re-arranged to be per-sub-simplex caching, */
/* rather than fxcell caching. Rather than stash the simplex info in the fxcells, */
/* create a separate cache or some other way of sharing the common simplexes. */
/* The code that ignores common face simplexes in fxcells could then be removed (?). */

/* Simplex definition. Each top level fwd interpolation cell, */
/* is decomposed into sub-simplexes. Sub-simplexes are of equal or */
/* lower dimensionality (ie. faces, edges, verticies) to the cube. */
struct _simplex {
	int refcount;				/* reference count */
	struct _rspl *s;			/* Pointer to parent rspl */
	int ix;						/* Construction Fwd cell index */
	int si;						/* Diagnostic - simplex number within level */
	int sdi;					/* Sub-simplex dimensionality */
	int efdi;					/* Effective fdi. This will be = fdi for a non clip */
								/* plane simplex, and fdi+1 for a clip plane simplex. */
								/* The DOF (Degress of freedom) within this simplex = sdi - efdi */

	psxinfo *psxi;				/* Generic per simplex info (construction cube relative) */
	int vix[MXRI+1];			/* fwd cell vertex indexes of this simplex [sdi+1] */
								/* This is a universal identification of this simplex. */
	struct _simplex *hlink;		/* Link to other cells with this hash */
	unsigned int touch;			/* Last touch count. */
	short flags;				/* Various flags */

#define SPLX_CLIPSX  0x01		/* This is a clip plane simplex */

#define SPLX_FLAG_1  0x04		/* v, linmin/max  initialised */
#define SPLX_FLAG_2  0x08		/* lu/svd initialised */
#define SPLX_FLAG_2F 0x10		/* lu/svd init. failed */
#define SPLX_FLAG_4  0x20		/* locus found */
#define SPLX_FLAG_5  0x40		/* auxiliary lu/svd initialised */
#define SPLX_FLAG_5F 0x80		/* auxiliary lu/svd init. failed */

#define SPLX_FLAGS  (SPLX_FLAG_1 | SPLX_FLAG_2 | SPLX_FLAG_2F \
                   | SPLX_FLAG_4 | SPLX_FLAG_5 | SPLX_FLAG_5F)

	double v[MXRI+1][MXRO+1]; 	/* Simplex Vertex values */
								/* v[0..sdi][0..fdi-1] are the output interpolation values */
								/* v[0..sdi][fdi] are the ink limit interpolation values */

								/* Baricentric vv[x][y] = (v[y][x] - v[y+1][x]) */
								/* and         vv[x][sdi] = v[sdi][x]           */

								/* Note that #num indicates appropriate flag number */
								/* and *num indicates a validator */

	double p0[MXRI]; 			/* Simplex base position = construction cube p[0] */
	double pmin[MXRI];			/* Simplex vertex input space min and */
	double pmax[MXRI];			/* max values [di] */

	double min[MXRO+1], max[MXRO+1]; /* Simplex vertex output space [fdi+1] and */
								/* ink limit bounding values at minmax[fdi] */

	/* If sdi == efdi, this holds the LU decomposition, */
	/* else this holds the SVD and solution locus info */

	char *aloc2;		/* Memory allocation for #2 & #4 */

	/* double **d_u;	   LU decomp of vv, U[0..efdi-1][0..sdi-1]		#2 */
	/* int *d_w;		   LU decomp of vv, W[0..sdi-1]					#2 */

	double **d_u;		/* SVD decomp of vv, U[0..efdi-1][0..sdi-1]		#2 */
	double *d_w;		/* SVD decomp of vv, W[0..sdi-1]				#2 */
	double **d_v;		/* SVD decomp of vv, V[0..sdi-1][0..sdi-1]		#2 */

						/* Degrees of freedom = dof = sdi - efdi */
	double **lo_l;		/* Locus coefficients,  [0..sdi-1][0..dof-1]	#2 */

	double *lo_xb;		/* RHS used to compute lo_bd [0..efdi-1]		*4 */
	double *lo_bd;		/* Locus base solution, [0..sdi-1]				#4 */

	unsigned auxbm;		/* aux bitmap mask for ax_lu and ax_svd 		*5 */
	int      aaux;		/* naux count for allocation					*5 */
	int      naux;		/* naux for calculation (may be < aaux ?)		*5 */

	/* if (sdi-efdi = dof) == naux this holds LU of lo_l */
	/* else this holds the SVD of lo_l */

	char *aloc5;		/* Memory allocation for #5 */

	/* double **ax_u;	   LU decomp of lo_l							#5 */
	/* int *ax_w;		   Pivot record for ax_lu decomp				#5 */

	double **ax_u;		/* SVD decomp of lo_l, U[0..naux-1][0..dof-1]	#5 */
	double *ax_w;		/* SVD decomp of lo_l, W[0..dof-1]				#5 */
	double **ax_v;		/* SVD decomp of lo_l, V[0..dof-1][0..dof-1]	#5 */

}; typedef struct _simplex simplex;

/* A candidate search (fwd) fxcell (cell cache entry structure) */
struct _fxcell {
	struct _rspl *s;		/* Pointer to parent rspl */

	/* Cache information */
	int ix;					/* Corresponding fwd cell index */
	struct _fxcell *hlink;	/* Link to other cells with this hash */
	struct _fxcell *mrudown;/* Links to next most recently used fxcell */
	struct _fxcell *mruup;
	int refcount;			/* Reference count */
	int flags;				/* Non-zero if the fxcell has been initialised */
#define CELL_FLAG_1  0x01	/* Basic initialisation, including nn_grp */
#define CELL_FLAG_2  0x02	/* Simplex information initialised */

	/* Use information */
	double sort;			/* Sort key */
	double limmin, limmax;	/* limit() min/max for this fxcell */

	/* Quick nn distance information */
	nn_grp g;	

//	double bcent[MXRO];		/* Output value bounding shere center */
//	double brad;			/* Output value bounding shere radius */
//	double bradsq;			/* Output value bounding shere radius squared */
//	double wbrad;			/* Output value weighted bounding shere radius */

	double p[POW2MXRI][MXRI]; /* Vertex input positions for this cube. */
							/* Copied to x->pmin/pmax[] & ink limit */

	double v[POW2MXRI][MXRO+1]; /* Vertex data for this cube. Copied to x->v[] */
							/* v[][fdi] is the ink limit values, if relevant */

	simplex **sx[MXRI+1];	/* Lists of simplexes that make up this fxcell. */
							/* Each list indexed by the non-limited simplex */
							/* dimensionality (similar to sspxi[]) */
							/* Each list will be NULL if it hasn't been created yet */
	int sxno[MXRI+1];		/* Corresponding count of each list */

}; typedef struct _fxcell fxcell;

/* surface bxcell sl status */
typedef enum {
	bx_uninit  = 0,		/* sl is not initialised */
	bx_filled  = 1,		/* sl has been filled with initial fwd cell vertexes */
	bx_rethinnd = 2,	/* sl vertexes need to be re-thinned */
	bx_thinned = 3,		/* sl vertexes have been thinned */
	bx_conv    = 4		/* sl vertexes have been converted to fwd cell indexes */
} bxstat;

/* Structure to hold bwd cell information for surface list, and also */
/* for seed fill bwd cell propogation. (Cells on surface will have two */
/* of these) */
struct _bxcell{
	int ix;						/* nnrev[] index of this bwd cell */
	int gc[MXRO];				/* coordinate of this bwd cell */
	nn_grp g;					/* Group nn calculation info */
	struct _bxcell *ss;			/* bwd surface cell to start search from */
	double sdist;				/* Est. wtd distance from this cell to ss */
	int tix;					/* target rev[] index being filled (visited check) */

	bxstat status;				/* State of sl list */					
	int *sl;					/* fwd vertex seed list for surface bxcells */
								/* or cell list after conversion to cells */
	int *dl;					/* deleted fwd vertex list for this bxcell */

	int *scell;					/* If non-NULL, this is a (non-surface) */
								/* seeding super bxcell, and scell contains */
								/* the list of bxcells covered */ 

	struct _bxcell *slist;		/* Linked list of all surface bxcells */
	struct _bxcell *hlink;		/* Linked list of surface bxcells with same ix hash */
	struct _bxcell *xlist;		/* Linked list of surface exploration search region */
	double emin;				/* estimated minimum wtd distance of this cell in current search */
	struct _bxcell *tlist;		/* Linked list of solution surface cells for current search. */

	struct _bxcell *flist;		/* Linked list for nnrev[] fill seeds */

	double cc;					/* Distance of group from gamut center */
	double dw;					/* Delta width from ocent of furthest point */
	struct _bxcell *wlist;		/* Linked list for shadow bxcells during thinning */

	int    debug;				/* debug flag - for VRML tagging */
}; typedef struct _bxcell bxcell;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


/* Enough space is needed to cache all the fxcells/simplexes */
/* for a full aux. locus, or the query will be processed in */
/* several "chunks" and will be slower. */
/* This sets the basic memory usage of the rev code. */

#define REV_ACC_GRES_MUL 2.0		/* 2.0 Reverse accelleration grid resolution */
									/* multiplier over fwd grid resolution */
#define REV_ACC_GRES_LIMIT 43		/* Reverse accelln. grid resolution limit before env. mult. */
#define REV_MAX_MEM_RATIO 0.3		/* 0.3 Proportion of first 1G of Ram to use */
#define REV_MAX_MEM_RATIO2 0.4		/* 0.4 Proportion of rest of Ram to use */
									/* rev as a fraction of the System RAM. */
#define HASH_FILL_RATIO 3			/* 3 Ratio of entries to hash size */

/* The structure where fxcells and simplexes are allocated and cached. */

/* Holds the fxcell and simplex match cache specific information */
typedef struct {
	struct _rspl *s;			/* Pointer to parent rspl */
	int nacells;				/* Number of allocated cells */
	int nunlocked;				/* Number of unlocked cells that could be freed */
	int cell_hash_size;			/* Current size of cell hash list */
	fxcell **hashtop;			/* Top of hash list [cell_hash_size] */
	fxcell *mrutop, *mrubot;	/* Top and bottom pointers of mru list */
								/* that tracks allocated fxcells */

	int spx_hash_size;			/* Current size of shared face simplex hash list */
	simplex **spxhashtop;		/* Shared face simplex hash index list */
	int nspx;					/* Number of simplexes in hash list */
} revcache;

/* common search information */
/* Type of (internal) reverse search */
enum ops {
	exact  = 0,		/* Search for all input values that exactly map to given output value */
	clipv  = 1,		/* Search for input values that map to outermost solution along a vector */
	clipn  = 2,		/* Search for input values that map to closest solution */
	auxil  = 3,		/* Search for input values that map to given output, and closest to auxiliary target */
	locus  = 4 };	/* Return range of auxiliary values that contains solution */

/* + possible exact with clip */

/* Structure to hold clip line state information */
typedef struct {
	struct _rspl *s;	/* Pointer to parent rspl */
	double st[MXRO];	/* start of line - reverse grid base value */
	double de[MXRO];	/* direction of line */
	int    di[MXRO];	/* incerement in line direction */
	int    ci[MXRO];	/* current rev grid coordinate */
	double t;			/* Parameter 0.0 - 1.0, line finished if t > 1.0 */
} line;


/* Structure to hold aux value of an intersection of a */
/* solution locus with a sub-simplex. Used when asegs flag is set */
typedef struct {
	double xval;		/* Auxiliary value */
	int nv;				/* Number of verticies valid */
	int vix[MXRI+1];	/* Verticy indexes of sub-simplex involved */
} axisec;

/* -------------------------------------------- */
/* Information needed/cached for reverse lookup */
struct _schbase {
	struct _rspl *s;	/* Pointer to parent rspl */

	int flags;			/* Hint flags */
	enum ops op;		/* Type of reverse search operation */
	int ixc;			/* Cube index of corner that holds maximum input values */

	int snsdi, ensdi;	/* Start and end extra sub-simplex dimensionality */
	int (*setsort)(struct _schbase *b, fxcell *c);	/* Function to triage & set cube sort index */
	int (*check)(struct _schbase *b, fxcell *c);		/* Function to recheck cube worth keeping */
	int (*compute)(struct _schbase *b, simplex *x);	/* Function to compute a simplex solution */

	double v[MXRO+1];	/* Target output value, + ink limit */
	double av[MXRI];	/* Target auxiliary values */
	int auxm[MXRI];		/* aux mask flags */
	unsigned auxbm;		/* aux bitmap mask */
	int naux;			/* Number of auxiliary target input values */
	int auxi[MXRI];		/* aux list of auxiliary target input values */
	double idist;		/* best input distance auxiliary target found (smaller is better) */
	int iabove;			/* Number of auxiliaries at or above zero */

	int canvecclip;		/* Non-zero if vector clip direction usable */
	double cdir[MXRO];	/* Clip vector direction and length wrt. v[] */
	double ncdir[MXRO];	/* Normalised clip vector */
	double **cla;		/* Clip vector LHS implicit equation matrix [fdi][fdi+1] (inc. ink tgt.) */
	double clb[MXRO+1];	/* Clip vector RHS implicit equation vector [fdi+1] (inc. ink tgt.) */
	double cdist;		/* Best clip locus distance found (aim is min +ve) :- weighted for nn */
	int    iclip;		/* NZ if result is above (disabled) ink limit */

	int mxsoln;			/* Maximum number of solutions that we want */
	int nsoln;			/* Current number of solutions found */
	co *cpp;			/* Store solutions here */

	int   lxi;			/* Locus search axiliary index */
	double min, max;	/* current extreme locus values for locus search */
	int    asegs;		/* flag - find all search segments */
	int    axisln;		/* Number of elements used in axisl[] */
	int    axislz;		/* Space allocated to axisl[] */
	axisec *axisl;		/* Auxiliary intersections */

	int lclistz;		/* Allocated space to fxcell sort list */
	fxcell **lclist;	/* Sorted list of pointers to candidate fxcells */

	int pauxcell;		/* Indexe of previous call solution fxcell, -1 if not relevant */
	int plmaxcell;		/* Indexe of previous call solution fxcell, -1 if not relevant */
	int plmincell;		/* Indexe of previous call solution fxcell, -1 if not relevant */

	int lsxfilt;		/* Allocated space of simplex filter list */
	char *sxfilt;		/* Flag for simplexes that should be in an fxcell */

	int rix;			/* Diagnostic - rev[] or nnrev[] index for this point */

}; typedef struct _schbase schbase;

/* ----------------------------------------- */

#ifdef STATS
struct _stats {
	int 	searchcalls;/* Number of top level searches */
	int 	csearched;	/* Cells searched */
	int 	ssearched;	/* Simplexes searched */
	int		sinited;	/* Simplexes initialised to base level */
	int		sinited2a;	/* Simplexes initialised to 2nd level with LU */
	int		sinited2b;	/* Simplexes initialised to 2nd level with SVD */
	int		sinited4i;	/* Simplexes invalidated at 4th level */
	int		sinited4;	/* Simplexes initialised to 4th level */
	int		sinited5i;	/* Simplexes invalidated at 5th level */
	int		sinited5a;	/* Simplexes initialised to 5th level with LU */
	int		sinited5b;	/* Simplexes initialised to 5th level with SVD */
	int		chits;		/* Cells hit in cache */
	int		cmiss;		/* Cells misses in cache */
}; typedef struct _stats stats;
#endif /* STATS */

/* ----------------------------------------- */
/* Reverse info stored in main rspl function */
struct _rev_struct {

	/* First section, basic information that doesn't change */
	/* Has been initialised if inited != 0 */

	int inited;			/* Non-zero if first section has been initialised */
						/* All other sections depend on this. */
	int fastsetup;		/* Flag - NZ if fast setup at cost of slow throughput */

	int probxyz;		/* Flag - NZ - guess if XYZ for VRML diagnostics use */

	int lchweighted;	/* Non-zero if nearest search is LCh weighted */
	double lchw[MXRO];	/* LCh nearest weighting */
	double lchw_sq[MXRO];	/* LCh nearest weighting squared */
	double lchw_chsq;		/* lchw_sq[1] - lchw_sq[2] */

	struct _rev_struct *next;	/* Linked list of global instances sharing memory */
	size_t max_sz;		/* Maximum size permitted */
	size_t sz;			/* Total memory current allocated by rev */

#ifdef NEVER
	int thissz, lastsz;	/* Debug reporting */
#endif

	/* Reverse grid lookup information */
	int res;			/* Reverse grid resolution == ncells on a side */
	int no;				/* Total number of points in reverse grid = rev.ares ^ fdi */
	int coi[MXRO];		/* Coordinate increments for each dimension */
	int hoi[1 << MXRO];	/* Coordinate increments for progress through cube */
	datao gl,gh,gw;		/* Reverse grid low, high, grid bwd cell width */

	/* Second section, accelleration information that changes with ink limit. */
	int rev_valid;		/* nz if this information in rev[] and nnrev[] is valid */
	int **rev;			/* Exact reverse lookup starting list. */
	int **nnrev;		/* Nearest reverse lookup starting list. */
						/* These lists are of fwd grid base indexes. */
						/* [0] is allocation length */
						/* [1] is the next free entry index (length + 3, not counting -1) */
						/* [2] is index into share lists, -1 if not shared. */
						/* Then follows cube indexes */
						/* Last entry is marked with -1 */

	double ocent[MXRO];	/* rev cell gamut "center" point for thinning and shadow testing. */

	int surflin_en;			/* Flag set when suflin is enabled */
	struct _rspl *surflin; /* gamut surface linearization transform used by logcomp() */
	double linoff[MXRO];	/* ocent offset after surflin mapping */

	bxcell *surflist;	/* Linked list of rev[] bwd cells that contain gamut surface fwd cells. */
						/* Used to speed up fill_nncell() when rev.fastsetup is set, else NULL */
	int surf_hash_size;	/* Current size of bxcell hash list */
	bxcell **surfhash;	/* bxcell hash index list */

	int **sharelist;	/* Array of pointers to shared (fwd grid list sharer) records. */ 
						/* Each record is same format as rev[]/nnrev[], except */
						/* [2] is used to detect scanning the same list. */
	int sharellen;		/* Size of sharelist */ 
	int sharelaloc;		/* Allocation size of sharelist */ 

	/* Third section */
	revcache *cache;	/* Where fxcells and simplexes are allocated and cached */
	/* Sub-dimension simplex information */
	ssxinfo sspxi[MXRI+1];/* One per sub dimenstionality at offset sdi */

	/* Fourth section */
	/* Has been initialise if sb != NULL */
	schbase *sb;		/* Structure holding calculated per-search call information */

	unsigned int stouch; /* Simplex touch count to avoid searching shared simplexs twice */
#ifdef STATS
	stats st[5];	/* Set of stats info indexed by enum ops */
#endif	/* STATS */

	int primsecwarn;	/* Not primary or secondary warning has been issued */

#ifdef CHECK_NNLU
int cknn_no;			/* Number checked */
double cknn_we;			/* Worst DE */
int cknn_noerrs;		/* Number not as good as closet vertex */
int cknn_nobsb;			/* Number not as good as second closest vertex */
int cknn_nonis;			/* Number not in surface */
#endif /* CHECK_NNLU */

}; typedef struct _rev_struct rev_struct;


/* ------------------------------------ */
/* Utility functions used by other parts of rspl implementation */

/* Initialise a static sub-simplex verticy information table */
void rspl_init_ssimplex_info(struct _rspl *s,	/* RSPL object */
ssxinfo *xip,				/* Pointer to sub-simplex info structure to init. */
int sdi);					/* Sub-simplex dimensionality (range 0 - di) */

/* Free the given sub-simplex verticy information */
void rspl_free_ssimplex_info(struct _rspl *s,
ssxinfo *xip);		/* Pointer to sub-simplex info structure */

#endif /* RSPL_REV_H */














