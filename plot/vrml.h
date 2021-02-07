
#ifndef VRML_H

/*
 * Simple diagnostic VRML function library for debugging
 *
 * Copyright 2005 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* Set global format using environment variable:

	ARGYLL_3D_DISP_FORMAT to "VRML", "X3D" or "X3DOM".
*/

/* Problems with .x3d using FreeWRL:

	Latest version makes surface colors look emissive, and you can't see
	the poligon edges == no detail.

	Lighting is weird - some stuff is black unless -1-1-1 light intensity is turned up.

	Poor navigation.

	Text size is different to Cosmo.
*/

/* Output format */
typedef enum {
    fmt_uninit  = -1,
    fmt_vrml    = 0,
    fmt_x3d     = 1,
    fmt_x3dom   = 2
} vrml_fmt;

/* (Make sure all callers are using vrml_space enum. before changing any values) */
typedef enum {
    vrml_lab    = 0,	/* L*a*b* 0..100 */
    vrml_xyz    = 1,	/* XYZ scale 0..1 */
    vrml_rgb    = 2		/* RGB scale 0..1 */
} vrml_space;


struct vrml_point {
	double pp[3];			/* Vertex position */
	double cc[3];			/* Vertex color */
	int last;				/* Last vertex of line flag */
};

struct vrml_triquad {
	int ix[4];
	double cc[3];			/* Per polygon color if ppoly, natural if cc[0] < 0 */
};

struct _vrml {

/* Private: */
	char *name;			/* Name including extension */
	FILE *fp;

	int written;		/* Set to nz when file has been written */

	vrml_fmt fmt;		/* Format to output. Defaults to global format */
						/* (Have to add ext() and format() methods if we change this) */ 
	vrml_space ispace;	/* 0 if Lab, 1 if XYZ plot, (range 0..1), 2 if RGB (range 0..1) */
	double scale;		/* Scale factor to use (applied before off) */
	double off;			/* Gamut L center, usually 50, to put Z 0 at center of view */

	/* Expandable point and line/patch arrays */
	struct {
		int npoints;
		int paloc;
		struct vrml_point *pary;

		/* Expandable line/triangle/quad vertex index */
		/* (Line has ix[2] = -1, Triangle has ix[3] = -1) */
		int ntrqu;
		int taloc;
		struct vrml_triquad *tqary;

		int ppoly;		/* Use per poligon color rather than vertex */

	} set[10];		/* Up to ten sets */

/* Public: */

	/* Methods */

	/* Return this files format extension (i.e. ".wrl" */
	char *(*ext)(struct _vrml *s);

	/* Return this files format type name */
	char *(*format)(struct _vrml *s);

	/* Write the file out. Return nz on error */
	int (*flush)(struct _vrml *s);

	/* Finish writing the file and free ourselves */
	void (*del)(struct _vrml *s);


	/* Add a spherical marker point to the plot. col == NULL for natural color  */
	/* rad is in normalized delta E scale units */
	/* (Need to do this before or after start_line_set()/dd_vertex()/make_lines() !) */
	void (*add_marker)(struct _vrml *s, double pos[3], double col[3], double rad);

	/* Add a spherical marker with transparency */
	/* rad is in normalized delta E scale units */
	void (*add_marker_trans)(struct _vrml *s, double pos[3], double col[3], double trans, double rad);

	/* Add a cone marker to the plot. col == NULL for natural color  */
	/* rad is in normalized delta E scale units */
	/* (Need to do this before or after start_line_set()/dd_vertex()/make_lines() !) */
	void (*add_cone)(struct _vrml *s, double p0[3], double p1[3], double col[3], double rad);

	/* Add a text marker to the plot. col == NULL for natural color  */
	/* size is in normalized delta E scale units */
	/* (Need to do this before or after start_line_set()/dd_vertex()/make_lines() !) */
	void (*add_text)(struct _vrml *s, char *text, double p[3], double col[3], double size);


	/* Start building up verticies that will be converted to lines, triangles or quads. */
	/* Set can be from 0 - 9 (call this before add_vertx() etc) */
	void (*start_line_set)(struct _vrml *s, int set);

	/* Add a verticy (default natural color from pos.) */
	/* Return the index number */
	int (*add_vertex)(struct _vrml *s, int set, double pos[3]);

	/* Add a verticy with per vertex color */
	/* Return the index number */
	int (*add_col_vertex)(struct _vrml *s, int set, double pos[3], double col[3]);


	/* Turn all the vertexes into a set of points */
	void (*make_points)(struct _vrml *s, int set);

	//void (*make_col_points)(struct _vrml *s, int set, double col[3]);


	/* Turn the last added vertex into the last vertex of the line (incremental line creation) */
	/* (Sets .last flag on vertex) */
	void (*make_last_vertex)(struct _vrml *s, int set);

	/* Convert the verticies to lines, ppset verticies per line (or using .last flag) */
	/* and output all the lines, using per vertex color. */
	/* Use ppset > no verticies for just .last flag */
	void (*make_lines)(struct _vrml *s, int set, int ppset);


	/* Add a line defined by vertex indexes using per vertex color */
	void (*add_line)(struct _vrml *s, int set, int ix[2]);

	/* Add a line defined by vertex index, and set per line color. */
	/* col == NULL or col[0] < 0.0 not to set per line color */
	void (*add_col_line)(struct _vrml *s, int set, int ix[2], double col[3]);

	/* Output lines using per vertex or per line colors. */
	void (*make_lines_vc)(struct _vrml *s, int set, double trans);

	/* Output lines with overall per line color */
	/* col == NULL or col[0] < 0.0 not to set overall color */
	void (*make_lines_cc)(struct _vrml *s, int set, double trans, double col[3]);


	/* Add a triangles defined by vertex indexes using per vertex color */
	void (*add_triangle)(struct _vrml *s, int set, int ix[3]);

	/* Add a triangle defined by vertex indexes, and set per triangle color. */
	/* col == NULL or col[0] < 0.0 not to set per line color */
	void (*add_col_triangle)(struct _vrml *s, int set, int ix[3], double col[3]);

	/* Output triangles using per vertex colors. */
	void (*make_triangles_vc)(struct _vrml *s, int set, double trans);

	/* Output triangles with overall per line color */
	/* col == NULL or col[0] < 0.0 not to set overall color */
	void (*make_triangles)(struct _vrml *s, int set, double trans, double col[3]);


	/* Add a quad defined by vertex index using per vertex color */
	void (*add_quad)(struct _vrml *s, int set, int ix[4]);

	/* Add a quad defined by vertex indexes, and set per quad color. */
	/* col == NULL or col[0] < 0.0 not to set per quad color */
	void (*add_col_quad)(struct _vrml *s, int set, int ix[4], double col[3]);

	/* Output quads using per vertex colors. */
	void (*make_quads_vc)(struct _vrml *s, int set, double trans);

	/* Output quads with overall per line color. */
	/* col == NULL or col[0] < 0.0 not to set overall color */
	void (*make_quads)(struct _vrml *s, int set, double trans, double col[3]);


	/* Clear verticies and lines/triangles/quads */
	void (*clear)(struct _vrml *s);
	
	/* Helper :- convert a Lab value to RGB */
	void (*Lab2RGB)(struct _vrml *s, double *out, double *in);

	/* Helper :- convert a XYZ value to RGB */
	void (*XYZ2RGB)(struct _vrml *s, double *out, double *in);

	/* Always put these last */
#ifdef GAMUT_H		/* If gamut.h is #included ahead of us */
	/* Create a solid gamut surface from the given gamut */
	/* trans is trasparency, cc is surface color, cc[0] < 0.0 for natural (Uses set 9) */
	void (*make_gamut_surface)(struct _vrml *s, gamut *g, double trans, double cc[3]);

	/* Create a solid or wireframe gamut surface from the given gamut (Uses set 9) */
	/* trans is trasparency, cc is surface color, cc[0] < 0.0 for natural */
	void (*make_gamut_surface_2)(struct _vrml *s, gamut *g, double trans, int wire, double cc[3]);

	/* Add cusp markers from the gamut surface) */
	void (*add_cusps)(struct _vrml *s, gamut *g, double trans, double cc[3]);
#endif /* GAMUT_H */

}; typedef struct _vrml vrml;

/* Return the global file format extension (i.e. ".wrl" */
char *vrml_ext();

/* Return the global file format type name */
char *vrml_format();

/* Create a vrml/x3d plot object. */
/* Filename will have appropriate extension added automatically. */
vrml *new_vrml(char *name, int doaxes, vrml_space ispace);

/* Same as above but override default Z viewing distance */
vrml *new_vrml_vdist(char *name, int doaxes, vrml_space ispace, double vdist);


#define VRML_H
#endif /* VRML_H */
