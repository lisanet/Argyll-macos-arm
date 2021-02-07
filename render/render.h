
#ifndef RENDER2D_H
#define RENDER2D_H

/* 
 * render2d
 *
 * Simple 2D raster rendering support.
 *
 * Author:  Graeme W. Gill
 * Date:    28/12/2005
 *
 * Copyright 2005, 2008, 2012, 2014 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This is basically a simple 2D ray tracing renderer, so it's not especially */
/* efficient, but it's simple and direct, easy to add new primitives or */
/* capabilities, is high quality, and has an accelleration algorithm that */
/* makes it fast enough for printed output. */

/* Mathematical coordinate in mm are used for primitives, ie. the origin is */
/* the bottom left corner. */
/* Device color values range from 0.0 to 1.0 */

#define MXCH2D  16			/* Maximum color channels */
#define TOTC2D  (MXCH2D+1)	/* Maximum total components */
#define PRIX2D  (MXCH2D)	/* Index of primitive kept with color value */

#define MXPATSIZE 4			/* Maximum color pattern size */

/* Color type */
/* Shouldn't this be an xcolorants mask ? */
typedef enum {
	w_2d,			/* Video style grey */
	k_2d,			/* Printing style grey */
	lab_2d,			/* Lab */
	rgb_2d,			/* RGB */
	cmyk_2d,		/* CMYK */
	ncol_2d,		/* N color */
	ncol_a_2d		/* N color with extra as alpha */
} colort2d;

/* Pixel depth */
typedef enum {
	bpc8_2d,		/* 8 bits per component */
	bpc16_2d		/* 16 bits per component */
} depth2d;

typedef double color2d[TOTC2D]; 

/* Font type */
typedef enum {
	rowman_s = 0,	/* Rownman, single stroke */
	rowman_d = 1,	/* Rownman, double stroke */
	rowman_t = 2,	/* Rownman, triple stroke */
	timesr   = 3,	/* Times Roman */
	timesr_b = 4,	/* Times Roman, Bold */
	futura_l = 5,	/* Futura, Light */
	futura_m = 6	/* Futura, Medium */
} font2d;

/* ------------------------------------ */

struct _render2d;

#define PRIM_STRUCT							\
/*	primt2d tag; */			/* Type of primitive */	\
	int    ix;				/* Index (order added) */ \
	int    ncc;				/* Number of color components */	\
	struct _prim2d *next;	/* Linked list to next primitive */ \
	struct _prim2d *yl0;	/* Previous lines Y list linked list */ \
	struct _prim2d *yl;		/* Active Y list linked list */ \
	struct _prim2d *xl;		/* Active X list linked list */ \
	double x0, y0, x1, y1;	/* Extent, top & left inclusive, bot & right non-inclusive */ \
	void (*del)(struct _prim2d *s);		/* Delete the object */ \
							/* Render the object at location. Return nz if in primitive */ \
	int (*rend)(struct _prim2d *s, color2d rv, double x, double y);

struct _prim2d {
	PRIM_STRUCT
}; typedef struct _prim2d prim2d;

/* ------------------------------------ */
/* Solid rectange primitive */
struct _rect2d {
	PRIM_STRUCT
	double rx0, ry0, rx1, ry1;	/* Rectangle verticies */
	color2d c;					/* Color of rectangle (if dpat == NULL) */
	double (*dpat)[MXPATSIZE][MXPATSIZE][TOTC2D];	/* Special for ChromeCast experiments */
	int dp_w, dp_h;				/* dpat dimensions */
}; typedef struct _rect2d rect2d;

prim2d *new_rect2d(struct _render2d *s, double x, double y, double w, double h, color2d c);

void set_rect2d_dpat(struct _rect2d *s, double (*pat)[MXPATSIZE][MXPATSIZE][TOTC2D], int w, int h);

/* ------------------------------------ */
/* Vertex shaded rectange */
struct _rectvs2d {
	PRIM_STRUCT
	double rx0, ry0, rx1, ry1;	/* Rectangle verticies */
	color2d c[4];	/* Bot left, bot right, top left, top right */
	int x_blend;	/* Blending rule flags, 0 = linear, 1 = spline, 2 = sine */
	int y_blend;
	int y_sine;
}; typedef struct _rectvs2d rectvs2d;

prim2d *new_rectvs2d(struct _render2d *s, double x, double y, double w, double h, color2d c[4]);

/* ------------------------------------ */
/* Vertex shaded triangle */
struct _trivs2d {
	PRIM_STRUCT
	double be[3][3];	/* baricentric equations */
	color2d c[3];		/* Color of each vertex */
}; typedef struct _trivs2d trivs2d;

prim2d *new_trivs2d(struct _render2d *s, double v[3][2], color2d c[3]);

/* ------------------------------------ */
/* Polygon */

/* Solid polygon primitive */
struct _poly2d {
	PRIM_STRUCT
	color2d c;			/* Color of polyangle (if dpat == NULL) */
	int n;				/* Number of verticies */
	double co[1][2];	/* Array of [n][2] verticies */
}; typedef struct _poly2d poly2d;

prim2d *new_poly2d(struct _render2d *s, int n, double v[][2], color2d c);

#ifdef NEVER
int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy) {
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

#endif // NEVER

/* ------------------------------------ */

/* Circular disk/line primitive */
struct _disk2d {
	PRIM_STRUCT
	double cx, cy;				/* Center */
	color2d c;					/* Color of disk/line */
	double orr, irr;			/* Outer radius squared, inner radius squared (0.0 if disk) */
}; typedef struct _disk2d disk2d;

/* Center and radius */
prim2d *new_disk2d(struct _render2d *s, double x, double y, double r, color2d c);

/* Center, radius and line width */
void add_circle2d(struct _render2d *s, double x, double y, double r, double w, color2d c);

/* ------------------------------------ */
/* A single line. */
	
struct _line2d {
	PRIM_STRUCT
	double lx0, ly0, lx1, ly1;	/* Line verticies */
	double ww;					/* half width of line squared */
	int cap;					/* 0 = butt, 1 = round, 2 = square */
	color2d c;					/* Color of the line */
	int t;						/* nz if line is degenerate */
	double vx, vy;				/* Vector relative to x0 y0 */
}; typedef struct _line2d line2d;

prim2d *new_line2d(struct _render2d *s, double x0, double y0, double x1, double y1, double w, int cap, color2d c);

/* ------------------------------------ */
/* Comound primitives */

/* add a dashed line */
void add_dashed_line2d(
struct _render2d *s,
double x0, double y0,
double x1, double y1,
double w,
double on, double off,
int cap,
color2d c);

/* Add a text character at the given location using lines */
void add_char2d(
struct _render2d *s,
double *xinc,		/* Return increment in position for next character */
double *yinc,
font2d fo,			/* Font to use */
char ch,			/* Character code to be printed */
double x, double y,	/* Location of bottom left of normal orientation text */
double h,			/* Height of text in normal orientation */
int or,				/* Orintation, 0 = right, 1 = down, 2 = left, 3 = up */
color2d c			/* Color of text */
);

/* Add a string from the given location using lines. */
void add_string2d(
struct _render2d *s,
double *xinc,		/* Return increment in position for next character */
double *yinc,
font2d fo,			/* Font to use */
char *string,		/* Character code to be printed */
double x, double y,	/* Location of bottom left of normal orientation text */
double h,			/* Height of text in normal orientation */
int or,				/* Orintation, 0 = right, 1 = down, 2 = left, 3 = up */
color2d c			/* Color of text */
);

/* Return the total width of the string without adding it */
void meas_string2d(
struct _render2d *s,
double *xinc,		/* Return increment in position for next character */
double *yinc,
font2d fo,			/* Font to use */
char *string,		/* Character code to be printed */
double h,			/* Height of text in normal orientation */
int or				/* Orintation, 0 = right, 1 = down, 2 = left, 3 = up */
);


/* Convert an angle in degrees (0 = right) */
/* into transform matrix */
void deg2mat(double mat[2][2], double a);

/* Convert a vector into a rotation matrix */
void vec2mat(double mat[2][2], double dx, double dy);

/* Add a text character at the given location using lines */
/* (matrix orientation version) */
void add_char2dmat(
struct _render2d *s,
double *xinc,		/* Return increment in position for next character */
double *yinc,
font2d fo,			/* Font to use */
char ch,			/* Character code to be printed */
double x, double y,	/* Location of bottom left of normal orientation text */
double h,			/* Height of text in normal orientation */
double mat[2][2],	/* Unity orientation matrix */
color2d c			/* Color of text */
);

/* Add a string from the given location using lines. */
/* (matrix orientation version) */
void add_string2dmat(
struct _render2d *s,
double *xinc,		/* Return increment in position for next character */
double *yinc,
font2d fo,			/* Font to use */
char *string,		/* Character code to be printed */
double x, double y,	/* Location of bottom left of normal orientation text */
double h,			/* Height of text in normal orientation */
double mat[2][2],	/* Unity orientation matrix */
color2d c			/* Color of text */
);

/* Return the total width of the string without adding it */
/* (matrix orientation version) */
void meas_string2dmat(
struct _render2d *s,
double *xinc,		/* Return increment in position for next character */
double *yinc,
font2d fo,			/* Font to use */
char *string,		/* Character code to be printed */
double h,			/* Height of text in normal orientation */
double mat[2][2]	/* Unity orientation matrix */
);

/* ------------------------------------ */

/* Type of output to save to. */
typedef enum {
	tiff_file,		/* Write a TIFF format file */
	png_file,		/* Write a PNG format file */
	png_mem			/* Write a PNG image to a memory buffer */
} rend_format;

/* ------------------------------------ */
/* Render object */

struct _render2d {

/* Private: */
	int ix;					/* Next primitive index */
	double fw, fh;			/* Page size in mm including margines */
	double lm, rm, tm, bm;	/* Page margines in mm */
	double w, h;			/* Page size in mm excluding margines */
	double hres, vres;		/* Page pixel resolution in pixels/mm */
	int pw, ph;				/* Page size in pixels */
	colort2d csp;			/* Color space */
	int      ncc;			/* Number of color components */
	depth2d  dpth;			/* Depth of the components */
	int    dither;     		/* Dither flag, 1 = ordered, 2 = error diffusion */
	int    noavg;     		/* Don't anti-alias or average 4 pixels together */
	int    dithfgo;			/* Dither F.G. only flag */
	void (*quant)(void *qcntx, double *out, double *in); /* optional quantization func. for edith */
	void *qcntx;
	double mxerr;			/* Maximum error diffusion error */

	color2d defc;			/* Default color value */

	void (*bgfunc)(void *cntx, color2d c, double x, double y);	/* BG color function */
	void *cntx;

	prim2d *head;			/* Start of list of primitives in rendering order */
	prim2d *yl;				/* Active Y list linked list head */
	prim2d *xl;				/* Active X list linked list head */

/* Public: */
	/* Methods */
	void (*del)(struct _render2d *s);					/* Free ourselves and all primitives */

	void (*set_defc)(struct _render2d *s, color2d c);	/* Set the default/background color */

	void (*set_bg_func)(struct _render2d *s,			/* Set background color function */
		void (*func)(void *cntx, color2d c, double x, double y),	/* Func can choose not to set */
		void *cntx
	);

	void (*add)(struct _render2d *s, prim2d *p);		/* Add a primitive */

	int (*write)(struct _render2d *s, char *filename, int comprn,
		unsigned char **obuf, size_t *olen,
	    rend_format fmt);
												/* Render and write to a TIFF or PNG file */
}; typedef struct _render2d render2d;

/* Constructor */
render2d *new_render2d(
	double w,		/* width in mm */
	double h,		/* height in mm */
	double ma[4],	/* Margines, left, right, top, bottom, NULL for zero in mm */
	double hres,	/* horizontal resolution in pixels/mm */
	double vres,	/* horizontal resolution in pixels/mm */
	colort2d csp,	/* Color type */
	int nd,			/* Number of channels if c = ncol */
	depth2d dpth,	/* Pixel depth */
	int dither,		/* Dither flag, 1 = ordered, 2 = error diffusion, | 0x8000 to dither FG only */
					/* | 0x4000 don't anti-alias by averaging pixels together. */
	void (*quant)(void *qcntx, double *out, double *in), /* optional quantization func. for edith */
	void *qcntx,
	double mxerr	/* Maximum error diffusion error */
);

#endif /* RENDER2D_H */

