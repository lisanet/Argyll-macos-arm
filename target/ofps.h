
#ifndef OFPS_H

/* 
 * Argyll Color Correction System
 *
 * Optimised Farthest Point Sampling
 *
 * Author: Graeme W. Gill
 * Date:   6/9/2004
 *
 * Copyright 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifndef MXPD
#define MXPD 4			/* Maximum ofps dimentionality */
#endif

#define MXNIX (MXPD+3)	/* Maximum vertex node indexes + hash + ixm */ 

struct _acell;

/* A gamut surface plane equation. */
struct _pleq {
	double pe[MXPD+1];	/* Vertex plane equation. First di elements are normalized */
						/* outward pointing normal (from first sample point), last element */
						/* is constant. If point . pe > 0, then point is outside surface */
	int ix;				/* Index of fake node associated with gamut surface (-ve) */
}; typedef struct _pleq pleq;

/* A sub-surface set mask */
#define MXSMASKW 6		/* Maximum number of setmask words */
typedef struct {
	unsigned int m[MXSMASKW];
} setmask;

/* A vertex. This is a point that has the highest eserr within */
/* the local region of di+1 nodes. It is therefore a candidate for a new */
/* sample node during seeding, or a point at which the sampling */
/* error of the current node locations can be estimated. */
/* Vertexes are the vericies of the Voronoi polyhedra. */
/* Because these are based on the natural eserr neighborhood, */
/* they are not necessarily exactly poyhedra. */
/* (Non gamut boundary verticies are shared) */
struct _vtx {
	int    no;			/* Serial number for id */
	int nix[MXNIX];		/* di+1 Sample point node indexes involved in vertex */
						/* Index is -ve if it is a fake gamut boundary node, */
						/* sorted largest to smallest (so fake gamut nodes are last) */
						/* [MXPD+1] is hash of all the node indexes, */
						/* [MXPD+2] is the OR of all the node ixm's */
	double ce[MXPD+1];	/* Estimated curvature error from mid point to each real node */
						/* (ce's are compacted, skipping fake boundary nodes) */
	int nnv;			/* Number of neighbour verticies (vertex net) */ 
	int _nnv;			/* Number allocated */
	struct _vtx **nv;	/* List of neighbour verticies. This includes hidden verts. */
						/* (If we didn't have to support INDEP_SURFACE, then there would */
						/*  be exactly di neighbour verticies.) */

	double p[MXPD];		/* Vertex location */
	double v[MXPD];		/* Subjective value at vertex (Labj) ? */
	double eperr;		/* Estimated position error */
	double eserr;		/* Estimated sampling error */
	double p_eperr;		/* Previous estimated position error */
	char ghost;			/* Don't use for optimization or stats. */
	char ifake;			/* A fake inside node */
	char ofake;			/* A fake outside node */
	char used;			/* Set to nz if already used for seeding */
	char bch;			/* nz if added to batch update list */
	char del;			/* Marked for deletion (used by add_to_vsurf()) */
	char add;			/* 1 = add to node, 2 = update hm (used by add_to_vsurf()) */
	char par;			/* Marked for deletion because it's already a parent node */
	int sch;			/* Sanity check houskeeping */
	setmask buvm;		/* Batch update to sub-surface hidden setmask */
	setmask bdvm;		/* Batch delete change to sub-surface hidden setmask */
	struct _vtx *batch;	/* Batch update list */
	int nsp;			/* Number of gamut surface planes it touches */
	pleq *sp[MXPD+1];	/* List of gamut surface planes */
	unsigned int pmask;	/* Gamut surface plane mask, from its location */
	unsigned int cmask;	/* Gamut surface composition mask, from it's parent nodes */
	setmask vm;			/* Sub-surface visiblility setmask */
	struct _vtx *link;	/* Linked list of free/used vtx's */
	struct _vtx **plp;	/* Pointer to link pointer in used list */

	struct _vtx *n;		/* Next in acceleration list */
	struct _vtx **pn;	/* Pointer to link pointer in acceleration list */
	int pci;			/* Accelleration grid index */

	struct _vtx *chn;	/* Next in cache index list */
	struct _vtx **pchn;	/* Pointer to link pointer in cache index list */

	int fflag;			/* fchl set flag set from s->fflag */
	struct _vtx *fchl;	/* Next in fixup check list */
	struct _vtx **pfchl;/* Pointer to link pointer */
	int fupcount;		/* Number of times vertex has fixed in a round */
	double fuptol;		/* Tollerance for fixing this vertex */
	struct _node *hnode;/* Hit node */
	double hitmarg;		/* Hit margine to node */
	struct _vtx **psvtxs;	/* Pointer to entry in s->svtxs[] */

	int cflag;			/* Vertex checked for hit flag, set from s->flag */
	int sflag;			/* Vertex search before flag, set from s->flag */
	int hflag;			/* Vertex search after hit flag, set from s->flag */
	char opqsq;			/* flag, on post hit search queue */
	struct _vtx *slist;	/* Breadth first search list link */
	int disth;			/* Search distance from hit vertex, valid when hflag == s->flag */
	struct _vtx *nxh;	/* Vtxs hit by node (to be deleted) list set by ofps_check_vtx_vn() */
	double nba_eperr;	/* node being added eperr, valid if cflag == s->flag */

	struct _vtx *dell;	/* Deleted/Not Deleted list */
}; typedef struct _vtx vtx;

/* A mid point. This is a point that has the highest eserr directly */
/* between two neighboring nodes. It is used during optimization */
/* to try and encourage even spacing between nodes. */
struct _mid {
	int    no;			/* Serial number for id */
	double p[MXPD];		/* Midpoint location */
	int nix[2];			/* The two sample point node indexes involved in midpoint */
	double ce[2];		/* Estimated curvature error from mid point to two nodes */
	double np;			/* Interpolation point between nux[0] and nix[1] */
	double v[MXPD];		/* Subjective value at midpoint (Labj) ? */
	double eperr;		/* Estimated position error */
	double eserr;		/* Estimated sampling error */
	int refc;			/* Reference count */
	struct _mid *link;	/* Linked list of free/used mid's */
	struct _mid **plp;	/* Pointer to link pointer in used list */
}; typedef struct _mid mid;


/* A measurement sample point node. */
/* Sample points are the points around which the Voronoi polyhedra */
/* are constructed. If a network is constructed between nodes that */
/* form a vertex, the netork will be the Delaunay tesselation. */
struct _node {
	int    ix;			/* Index of node in s->n[] */
	int    ixm;			/* Hash mask of node ix */
	int    fx;			/* nz if point is fixed */
	double p[MXPD];		/* Device coordinate position */
	double v[MXPD];		/* Subjective value (Labk) ? */

	double np[MXPD];	/* Next device coordinates during opt */
	double nv[MXPD];	/* Next subjective coordinates during opt */

	double op[MXPD];	/* Previous device coordinates during opt */

	int nvv;			/* Number of Voronoi surface verticies */ 
	int _nvv;			/* Number allocated */
	vtx **vv;			/* List of Voronoi surface verticies */

	int nvn;			/* Number of Voronoi nodes & midpoints */ 
	int _nvn;			/* Number allocated */
	int *vn;			/* List of Voronoi nodes indexes. Doesn't include this node. */
						/* Index is -ve if it is a fake gamut boundary node. */
	mid **mm;			/* List of midpoints. Midpoints will be NULL if not created yet. */

	int nsp;			/* Number of touched gamut surface planes */
	pleq *sp[MXPD+1];	/* List of touched gamut surface planes */
	unsigned int pmask;	/* Gamut surface plane mask */

	struct _acell *cell;/* Pointer to cell vertex is in */
	struct _node *n;	/* Next in acceleration list */
	struct _node **pn;	/* Pointer to link pointer in acceleration list */
	int pci;			/* Accelleration grid index */

	int flag;			/* Node being added access flag, set from s->flag */
	int nvnflag;		/* node_recomp_nvn_dmxs access flag */
	struct _node *na;	/* Next in 'to be added' list */
	int upflag;			/* Set to s->flag if node has been added to ->nup list */
	struct _node *nup;	/* Next node in 'recomp_nvn' list */

}; typedef struct _node node;

#define BOUND_GFLAG ((unsigned int)-1)		/* Boundary cell gflag */

/* An acceleration structure cube */
struct _acell {
	unsigned int gflag;	/* Acceleration grid search touched & boundary flag */
	node *head;			/* List of nodes inside acceleration cell */
	vtx  *vhead;		/* List of verticies with all real nodes inside acceleration cell */
	int co[MXPD];		/* coordinate of cell */
	double p[MXPD];		/* Device position of base of cell */
	double v[MXPD];		/* Corresponfing perceptual value of base of cell */
	double cp[MXPD];	/* Device position of center of cell */
	double cv[MXPD];	/* Corresponfing perceptual value of center of cell */
	double eperr;		/* Worst case eperr from a corner to the center */
	struct _acell *slist;	/* Search list */
}; typedef struct _acell acell;

/* Storage for a node combination/new position */
struct _nodecomb {
	int nix[MXNIX];		/* di+1 Sample point node indexes involved in vertex */
	double ce[MXPD+1];	/* Estimated curvature error from mid point to two nodes */
	vtx **v1, **v2;		/* Deleted and non-deleted vertex involved */
	int _count;			/* v1/v2 allocation */
	int count;			/* Number of times this is generated, used of v1/v2 */
	double p[MXPD];		/* Position of resulting vertex or */
	double v[MXPD];		/* Value of resulting vertex or */
	double eperr;
	double eserr;
	double weserr;
	int pvalid;			/* Valid new position, else use deleted location & eserr */
	vtx *vv;			/* Existing vertex at this combination (if opt) */
	double ceperr;		/* Current eperr to be bettered */
	double oog;			/* Out of gamut value */
	setmask vm;			/* Sub-surface visibility setmask */
	int startex;		/* nz if existing p[] should be starting dnsqe point */
}; typedef struct _nodecomb nodecomb;

/* Vertex cache hash index/table size */
#define VTXCHSIZE 33037
//#define VTXCHSIZE 67493

/* Record of a set of gamut surface plane combination */
struct _surfcomb {
	unsigned int co;		/* Surface combination mask */
	int valid;				/* Valid flag */
	int nos;				/* Number surfaces */
	int smset;				/* i_sm has been set */
	setmask i_sm;			/* Indiviual set of this surface combination */
	setmask a_sm;			/* Accumulated set of this and higher dimensions */
	struct _surfcomb *ds;	/* Circular list of the disjoint set */
}; typedef struct _surfcomb surfcomb;

/* Main sample point object */
struct _ofps {
/* private: */
	int verb;		/* Verbose */

	int di;			/* Point dimensionality */
	double ilimit;	/* Ink limit - limit on sum of p[] */
	double imin[MXPD];	/* Ink limit - limit on min of p[], must be >= 0.0 */
	double imax[MXPD];	/* Ink limit - limit on min of p[], must be <= 1.0  */
	int good;			/* 0 = fast, 1 = good flag */
	double surftol;		/* Surface tollerance distance */
	int maxits;			/* Maximum itterative improvement passes */
	double lperterb;	/* level of random peturbation */
	double ssurfpref, esurfpref;	/* Start and end surface preference weightig */ 

	/* Error estimate model parameters */
	double devd_wght;	/* Device space weighting */
	double perc_wght;	/* Perceptual space weighting */
	double curv_wght;	/* Curvature weighting */

	int fxno;		/* Total number of fixed points provided, possibly non-unique */
	fxpos **ufx;	/* fnp randomized unique fixed points to add */

	int gnp;		/* Number of fake gamut nodes (-ve index) */
					/* -1 to -2di-1 are fake boundary node indexes, */
					/* with -2di-1 being the ink limit boundary. */
					/* -2di-2 is the fake inside node, */
					/* -2di-3 is the fake outside node, */
	int fnp;		/* Number of unique fixed points in list */
	int tinp;		/* Target number of total points in list, including fnp */
	int np;			/* Number of points currently in list */
	node *_n, **n;	/* tinp allocation of points, list of pointers to points */
	int nv;			/* Current number of verticies */
	int nxvno;		/* Next vertex serial number */
	int nxmno;		/* Next midpoint serial number */

	/* Gamut surface definition planes. */
	int nbp;			/* Number of boundary planes. Either 2di or 2di+1 */
	pleq gpeqs[2 + MXPD * 2 + 1];		/* Plane equations associated */
						/* with the fake gamut surface/boundary nodes */
						/* points. Index is 1-ix */
						/* (allow for other fakes, just in case) */
	int sminit;			/* Flag, nz if sc has been inited */
	surfcomb *sc;		/* Array of (1 << nbp) surface combinations */
	int smbits;			/* Total set mask bits */
	int bpsmw;			/* Bits per set mask word */
	int nsmw;			/* Number of setmask words */
	unsigned int lwmask;	/* Last word mask */

	/* Perceptual function handed in. All device values must have been */
	/* clipped before calling this. (On running, is replaced with rspl */
	/* cached version) */
	void (*percept)(void *od, double *out, double *in);
	void *od;			/* Opaque data for perceptual point */

	int pcache_res;		/* Grid resolution of pcache */ 
	rspl *pcache;		/* cache of perceptual lookup */
	
	/* Other info */
	int rix;			/* Next read index */
	double mn,mx,av;	/* serr stats */
	double smns;		/* Closest node spacing */
	double mxmvsq;		/* Maximum movement during optimisation */
	int optit;			/* Optimization itteration */
	vtx *nxh;			/* Vtxs hit by node (to be deleted) list set by ofps_check_vtx_vn() */
	vtx *nxp;			/* Vtxs that are already parents of added node list (fixups) */
	struct _vtx *batch;	/* Batch update list */
	int checklev;		/* check node recursion level */
	struct _vtx *fchl;	/* Next vertex in fixup check list */
	struct _vtx **svtxs;/* Vertex "to be fixed" list */ 
	int nsvtxs;			/* Number in vertex "to be fixed" list */
	int _nsvtxs;		/* Allocated size of vertex "to be fixed" list */
	struct _node *nup;	/* Next node in 'recomp_nvn' list */

	/* Used and free lists */
	struct _vtx *uvtx;	/* Linked list of used vtx's */
	struct _vtx *fvtx;	/* Linked list of free vtx's */
	struct _vtx *hvtx;	/* Linked list of hidden vtx's */
	struct _mid *umid;	/* Linked list of used mid's */
	struct _mid *fmid;	/* Linked list of free mid's */

	/* Unbounded perceptual model */
	double pmod[MXPD * (1 << MXPD)];
	int pmod_init;		/* It's been initialised */

	/* Acceleration structure */
	int agres;			/* Acceleration grid resolution (not including extra row) */
	double gw; 			/* Grid cell width */
	int gim[MXPD];		/* Grid index multiplier */
	double gcd;			/* Grid cell diagonal */
	int nig;			/* Number of cells in grid (including guard rows) */
	acell *_grid;		/* Pointer to allocated array of grid structures */
	acell *grid;		/* Pointer to base of array of grid structures */
	int agrid_init;		/* accell grid p[] and v[] have been inited */
	unsigned int gflag;	/* Acceleration grid search flag */
	int nacnl;			/* Number of bytes in Accelleration cell neighbour offset list */
	int *acnl;			/* Accelleration cell neighbour offset list */

	int flag;		/* Access flag associated with node being added */
	int nvnflag;	/* node_recomp_nvn_dmxs access flag */
	int fflag;		/* Fixup round flag */
	vtx *vch[VTXCHSIZE];	/* Vertex cache index */

	aat_atree_t *vtreep;	/* Binary tree of vertexes sorted by eperr */
	aat_atree_t *vtrees[MXPD+2];	/* Per nsp, binary tree of vertexes sorted by eserr */
							/* We get di+2 planes for fake initial nodes */

	/* Utility - avoid re-allocation/initialization */
	sobol *sob;
	nodecomb *combs;	/* New node combinations being created in add_to_vsurf() */
	int _ncombs;  /* Number of node combinations allocated. */

	/* Debug/stats */
	int nopstop;	/* Number of optimization passes before stopping with diagnostics */
	int ntostop;	/* Number of points before stopping with diagnostics */
	int positions;	/* Number of calls to locate vertex */
	int dnsqs;		/* Number of dnsq is called */
	int funccount;	/* Number of times dnsq callback function is called */
	int maxfunc;	/* Maximum function count per dnsq */
	int sucfunc;	/* Function count per sucessful dnsq */
	int sucdnsq;	/* Number of sucessful dnsqs */
	int maxretries;	/* Maximum retries used on sucessful dnsq */
	int posfails;	/* Number of position_vtx failures */
	int posfailstp;	/* Number of position_vtx failures this pass */
	int nvtxcreated;	/* Number of vertexes created */
	int nvtxdeleted;	/* Number of vertexes deleted */
	int add_hit;	/* Number of add_to_vsurf hits */
	int add_mis;	/* Number of add_to_vsurf misses */
	int fadd_hit;	/* Number of fixup add_to_vsurf hits */
	int fadd_mis;	/* Number of fixup add_to_vsurf misses */

	int nvcheckhits;	/* Number of vertexes hit durint ofps_check_node() */
	int vvchecks;		/* Number of vertexes checked for hit */
	int vvpchecks;		/* Number of vertexes that would be exaustively checked for hits */
	int naccsrch;		/* Number of accellerated searches */
	int ncellssch;		/* Number of accelleration cells searched */
	int nnschd;			/* Number of nodes nearest searched */
	int nnfschd;		/* Number of nodes for full nearest searched */
	int nvschd;			/* Number of vertexes nearest searched */
	int nvfschd;		/* Number of vertexes for full nearest searched */
	
	int nsurfadds;		/* Number of add_to_vsurf() calls */
	int nhitv;			/* Number of hit vertexes in add_to_vsurf() */
	int maxhitv;		/* Maximum number of hit vertexes in add_to_vsurf() */

	int nfseeds;		/* Number of times searched for vtx with largest eserr */
	int nfseedsvtx;		/* Number of vertexes searched for vtx with largest eserr */
	
	unsigned int l_mstime;	/* Last create surface pass timestamp */
	int l_positions;	/* Last Number of calls to locate vertex */
	int l_nvtxcreated;	/* Last Number of vertexes created */
	int l_nvtxdeleted;	/* Last Number of vertexes deleted */

	struct _vtx *i_uvtx;	/* Incrementallu created vertexes used by DEBUG_RESEED_AFTER_FIXUPS */

/* public: */
	/* return non-zero if the perceptual point is within the device gammut */
	int (*pig)(struct _ofps *s, double *p);

	/* Initialise, ready to read out all the points */
	void (*reset)(struct _ofps *s);

	/* Read the next set of non-fixed points values */
	/* return non-zero when no more points */
	/* p = position, v = value, either may be NULL */
	int (*read)(struct _ofps *s, double *p, double *v);

	/* Calculate and print stats */
	void (*stats)(struct _ofps *s);

	/* Destroy ourselves */
	void (*del)(struct _ofps *s);

	}; typedef struct _ofps ofps;


/* Constructor */
extern ofps *new_ofps(
	int verb,
	int di, double ilimit, int npoints,
	int good,
	double dadaptation,
	double devd_wght,
	double perc_wght,
	double curv_wght,
	fxpos *fxlist, int fxno, 		/* Existing, fixed point list */
	void (*percept)(void *od, double *out, double *in), void *od);

/* Extended constructor */
ofps *new_ofps_ex(
int verb,				/* Verbosity */
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
double *imin,			/* Ink limit - limit on min of p[], normally >= 0.0 */
double *imax,			/* Ink limit - limit on min of p[], normally <= 1.0  */
int tinp,				/* Total number of points to generate, including fixed */
int good,				/* 0 = fast, 1 = good */
double dadaptation,		/* Degree of adaptation to device characteristic 0.0 - 1.0, */
						/* use -ve value for explicit weightings: */
double devd_wght,		/* Device space weighting (if dad < 0) */
double perc_wght,		/* Perceptual space weighting (if dad < 0) */
double curv_wght,		/* Curvature weighting (if dad < 0) */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od,				/* context for Perceptual function */
int ntostop,			/* Debug - number of points until diagnostic stop */
int nopstop				/* Debug - number of optimizations until diagnostic stop, -1 = not */
);

#define OFPS_H
#endif /* OFPS_H */
