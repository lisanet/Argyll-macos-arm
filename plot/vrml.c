
/*
 * Simple diagnostic VRML function library for debugging
 *
 * Copyright 2005 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD:
 * 
 *	It would be nice to be able to select a data format as output for use with
 *	other software, i.e. meshlab.
 *
 */

/* NOTES:

	X3DOM commands:

		Examine Mode (activate with key e):
			Left Button / Left Button + Shift 			Rotate
			Mid Button / Left Button + Ctl 				Pan
			Right Button / Wheel / Left Button + Alt 	Zoom
			Left double click 							Set center of rotation

		n		back to normal view
		e		examine mode
		a		show all
		u		upright
		<space>	stats pane
		<????>	log
		

		<ctrl>	move
		<alt>	zoom

		<scroller>	zoom
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "numlib.h"
#include "icc.h"
#include "gamut.h"
#include "vrml.h"

#ifdef NT		/* You'd think there might be some standards.... */
# ifndef __BORLANDC__
#  define stricmp _stricmp
# endif
#else
# define stricmp strcasecmp
#endif

/* Convert input values to x,y, z */
static void cs2xyz(vrml *s, double *out, double *in) {
	if (s->ispace == vrml_rgb) {			/* RGB */
		out[0] = s->scale * in[0];
		out[1] = s->scale * in[1];
		out[2] = s->scale * in[2];
	} else if (s->ispace == vrml_xyz) {		/* XYZ */
		out[0] = s->scale * in[1];
		out[1] = s->scale * in[2];
		out[2] = s->scale * in[0] - s->off;
	} else {								/* Lab */
		out[0] = s->scale * in[1];
		out[1] = s->scale * in[2];
		out[2] = s->scale * in[0] - s->off;
	}
} 

/* Add a sphere at the given location, with transparency. */
/* If col[] is NULL, use natural color. */
/* rad is in normalized delta E scale units */
/* Need to do this before or after start_line_set()/dd_vertex()/make_lines() ! */
static void add_marker_trans(vrml *s, double pos[3], double col[3], double trans, double rad) {
	double rgb[3], xyz[3];

	if (rad <= 0.0)
		rad = 1.0;

	if (col == NULL || col[0] < 0.0) {
		if (s->ispace == vrml_rgb)			/* RGB */
			icmCpy3(rgb, pos);
		else if (s->ispace == vrml_xyz)		/* XYZ */
			s->XYZ2RGB(s, rgb, pos);
		else								/* Lab */
			s->Lab2RGB(s, rgb, pos);	
	} else {
		rgb[0] = col[0];
		rgb[1] = col[1];
		rgb[2] = col[2];
	}

	cs2xyz(s, xyz, pos);
	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"    # Shere\n");
		fprintf(s->fp,"    Transform { translation %f %f %f\n", xyz[0], xyz[1], xyz[2]);
		fprintf(s->fp,"      children [\n");
		fprintf(s->fp,"        Shape{\n");
		fprintf(s->fp,"          geometry Sphere { radius %f }\n", rad);
		fprintf(s->fp,"          appearance Appearance { material Material { \n");
		if (trans > 0.0) {
			fprintf(s->fp,"              transparency %f, \n",trans);
		} 
		fprintf(s->fp,"                  diffuseColor %f %f %f } }\n", rgb[0], rgb[1], rgb[2]);
		fprintf(s->fp,"        }\n");
		fprintf(s->fp,"      ]\n");
		fprintf(s->fp,"    }\n");

	} else {
		fprintf(s->fp,"    <!-- Shere -->\n");
		fprintf(s->fp,"    <Transform translation='%f %f %f'>\n", xyz[0], xyz[1], xyz[2]);
		fprintf(s->fp,"      <Shape>\n");
		fprintf(s->fp,"        <Appearance>\n");
		if (trans > 0.0) {
			fprintf(s->fp,"          <Material diffuseColor='%f %f %f'\n", rgb[0], rgb[1], rgb[2]);
			fprintf(s->fp,"                    transparency='%f'></Material>\n", trans);
		} else {
			fprintf(s->fp,"          <Material diffuseColor='%f %f %f'></Material>\n", rgb[0], rgb[1], rgb[2]);
		}
		fprintf(s->fp,"        </Appearance>\n");
		fprintf(s->fp,"        <Sphere radius='%f'></Sphere>\n",rad);
		fprintf(s->fp,"      </Shape>\n");
		fprintf(s->fp,"    </Transform>\n");
	}
}

/* Add a sphere at the given location. */
/* if col[] is NULL, use natural color. */
/* rad is in normalized delta E scale units */
/* Need to do this before or after start_line_set()/dd_vertex()/make_lines() ! */
/* (Hasn't been fixed to work in RGB space) */
static void add_marker(vrml *s, double pos[3], double col[3], double rad) {
	add_marker_trans(s, pos, col, 0.0, rad);
}

/* Add a cone marker to the plot. col == NULL for natural color  */
/* rad is in normalized delta E scale units */
/* Need to do this before or after start_line_set()/dd_vertex()/make_lines() ! */
static void add_cone(vrml *s, double pp0[3], double pp1[3], double col[3], double rad) {
	double rgb[3];
	double p0[3], p1[3];

	icmScale3(p0, pp0, s->scale);
	icmScale3(p1, pp1, s->scale);

//printf("~1 cone %f %f %f -> %f %f %f rad %f\n", p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], rad);

	if (rad <= 0.0)
		rad = 1.0;

	if (col == NULL || col[0] < 0.0) {
		icmAdd3(rgb, p1, p0);
		icmScale3(rgb, rgb, 0.5);		/* Compute half way value */

		if (s->ispace == vrml_rgb)			/* RGB */
			icmCpy3(rgb, rgb);
		else if (s->ispace == vrml_xyz)		/* XYZ */
			s->XYZ2RGB(s, rgb, rgb);
		else								/* Lab */
			s->Lab2RGB(s, rgb, rgb);	
	} else {
		rgb[0] = col[0];
		rgb[1] = col[1];
		rgb[2] = col[2];
	}

	p0[0] -= s->off;
	p1[0] -= s->off;

	{
		double base[3] = { 0.0, 0.0, 1.0 };		/* Default orientation of cone is b axis */
		double len;
		double loc[3];
		double vec[3];
		double axis[3];		/* Axis to rotate around */
		double rot;			/* In radians */
		int j;

		//printf("~1 edge vert %d to %d\n",tp->v[0]->n, tp->v[1]->n);
		//printf("~1 edge %f %f %f to %f %f %f\n",
		//tp->v[0]->ch[0], tp->v[0]->ch[1], tp->v[0]->ch[2],
		//tp->v[1]->ch[0], tp->v[1]->ch[1], tp->v[1]->ch[2]);

		icmAdd3(loc, p1, p0);
		icmScale3(loc, loc, 0.5);		/* Compute half way value */
		icmSub3(vec, p1, p0);
		len = icmNorm3(vec);
		//printf("~1 loc = %f %f %f\n", loc[0], loc[1], loc[2]);
		//printf("~1 vec = %f %f %f\n", vec[0], vec[1], vec[2]);
		//printf("~1 len = %f\n", len);

		if (len < 0.1)
			len = 0.1;

		icmNormalize3(base, base, 1.0);
		icmNormalize3(vec, vec, 1.0);
		icmCross3(axis, base, vec);
		rot = icmDot3(base, vec);
		//printf("~1 base = %f %f %f\n", base[0], base[1], base[2]);
		//printf("~1 vec = %f %f %f\n", vec[0], vec[1], vec[2]);
		//printf("~1 axis = %f %f %f, rot = %f\n",axis[0],axis[1],axis[2],rot);
		if (icmNorm3sq(axis) < 1e-10) {		/* 0 or 180 degrees */
			double base2[3];
			int mxi = 0;
			//printf("~1 computing a different axis\n");
			base2[0] = vec[1];		/* Comute vector in a different direction */
			base2[1] = vec[2];
			base2[2] = vec[0];
			for (j = 1; j < 3; j++) {
				if (fabs(base2[j]) > fabs(base2[mxi])) 
					mxi = j;
			}
			base2[mxi] = -base2[mxi];
				
			icmCross3(axis, base2, vec);
			if (icmNorm3sq(axis) < 1e-10) {		/* 0 or 180 degrees */
				error("VRML rotate axis still too small");
			}
			if (rot < 0.0)
				rot = 3.1415926;
			else
				rot = 0.0;			
		} else {
			rot = acos(rot);
			//printf("~1 rotation %f\n",rot);
		}

		if (s->fmt == fmt_vrml) {
			fprintf(s->fp,"\n");
			fprintf(s->fp,"    # Cone\n");
			fprintf(s->fp,"    Transform {\n");
			fprintf(s->fp,"      rotation %f %f %f %f\n",axis[1], axis[2], axis[0], rot);
			fprintf(s->fp,"      translation %f %f %f\n",loc[1], loc[2], loc[0]);
			fprintf(s->fp,"      children [\n");
			fprintf(s->fp,"		Shape { \n");
			fprintf(s->fp,"		 geometry Cone { bottomRadius %f height %f }\n",rad,len);
			fprintf(s->fp,"        appearance Appearance { material Material { diffuseColor %f %f %f } }\n",rgb[0],rgb[1],rgb[2]);
			fprintf(s->fp,"		} \n");
			fprintf(s->fp,"      ]\n");
			fprintf(s->fp,"    }\n");
		} else {
			fprintf(s->fp,"\n");
			fprintf(s->fp,"    <!-- Cone -->\n");
			fprintf(s->fp,"    <Transform rotation='%f %f %f %f'\n", axis[1], axis[2], axis[0], rot);
			fprintf(s->fp,"               translation='%f %f %f'>\n", loc[1], loc[2], loc[0]);
			fprintf(s->fp,"      <Shape>\n");
			fprintf(s->fp,"        <Appearance>\n");
			fprintf(s->fp,"          <Material diffuseColor='%f %f %f'></Material>\n", rgb[0], rgb[1], rgb[2]);
			fprintf(s->fp,"        </Appearance>\n");
			fprintf(s->fp,"        <Cone bottomRadius='%f' height='%f'></Cone>\n",rad, len);
			fprintf(s->fp,"      </Shape>\n");
			fprintf(s->fp,"    </Transform>\n");
		}
	}
}

/* Add a text marker to the plot. col == NULL for natural color  */
/* size is in normalized delta E scale units */
/* (Need to do this before or after start_line_set()/dd_vertex()/make_lines() !) */
static void add_text(vrml *s, char *text, double p[3], double col[3], double size) {
	double rgb[3], xyz[3];

	if (size <= 0.0)
		size = 1.0;

	if (col == NULL || col[0] < 0.0) {
		if (s->ispace == vrml_rgb)			/* RGB */
			icmCpy3(rgb, p);
		else if (s->ispace == vrml_xyz)		/* XYZ */
			s->XYZ2RGB(s, rgb, p);
		else								/* Lab */
			s->Lab2RGB(s, rgb, p);	
	} else {
		rgb[0] = col[0];
		rgb[1] = col[1];
		rgb[2] = col[2];
	}
	cs2xyz(s, xyz, p);

	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"    # Text\n");
		fprintf(s->fp,"    Transform { translation %f %f %f\n", xyz[0], xyz[1], xyz[2]);
		fprintf(s->fp,"      children [\n");
		fprintf(s->fp,"        Shape{\n");
		fprintf(s->fp,"          geometry Text { string [\"%s\"]\n",text);
		fprintf(s->fp,"            fontStyle FontStyle { family \"SANS\" style \"BOLD\" size %f }\n",
		                                                                                      size);
		fprintf(s->fp,"                        }\n");
		fprintf(s->fp,"          appearance Appearance { material Material ");
		fprintf(s->fp,"{ diffuseColor %f %f %f } }\n", rgb[0], rgb[1], rgb[2]);
		fprintf(s->fp,"        }\n");
		fprintf(s->fp,"      ]\n");
		fprintf(s->fp,"    }\n");
	
	} else {
		fprintf(s->fp,"    <!-- Text -->\n");
		fprintf(s->fp,"    <Transform translation='%f %f %f'>\n", xyz[0], xyz[1], xyz[2]);
		fprintf(s->fp,"      <Shape>\n");
		fprintf(s->fp,"        <Appearance>\n");
		fprintf(s->fp,"          <Material diffuseColor='%f %f %f'></Material>\n", rgb[0], rgb[1], rgb[2]);
		fprintf(s->fp,"        </Appearance>\n");
		fprintf(s->fp,"        <Text string='\"%s\"'>\n",text);
		fprintf(s->fp,"          <FontStyle family='\"SANS\"' style='BOLD' size='%f'></FontStyle>\n", size);
		fprintf(s->fp,"        </Text>\n");
		fprintf(s->fp,"      </Shape>\n");
		fprintf(s->fp,"    </Transform>\n");
	}
}

/* Start building up verticies that will be converted to lines or patches. */
/* This clears the data from an existing set */
/* Set can be from 0 - 9 */
static void start_line_set(vrml *s, int set) {
	if (set < 0 || set > 9)
		error("vrml start_line_set set %d out of range",set);
	s->set[set].npoints = 0;
	s->set[set].ntrqu = 0;
	s->set[set].ppoly = 0;
}

/* Add a verticy with color. */
/* col == NULL or col[0] < 0.0 for natural color */
/* Return the index number */
static int add_col_vertex_l(vrml *s, int set, double pos[3], double col[3], int last) {

	if (set < 0 || set > 9)
		error("vrml add_col_vertex_l set %d out of range",set);

	if (s->set[set].npoints >= s->set[set].paloc) {
		s->set[set].paloc = (s->set[set].paloc + 10) * 2;
		if (s->set[set].pary == NULL)
			s->set[set].pary = malloc(s->set[set].paloc * sizeof(struct vrml_point));
		else
			s->set[set].pary = realloc(s->set[set].pary, s->set[set].paloc * sizeof(struct vrml_point));

		if (s->set[set].pary == NULL)
			error("VRML malloc failed at count %d\n",s->set[set].paloc);
	}
	s->set[set].pary[s->set[set].npoints].pp[0] = pos[0];
	s->set[set].pary[s->set[set].npoints].pp[1] = pos[1];
	s->set[set].pary[s->set[set].npoints].pp[2] = pos[2];

	/* Make lines/triangle etc. will convert to natural color if col[0] < 0.0 */
	if (col == NULL || col[0] < 0.0) {
		s->set[set].pary[s->set[set].npoints].cc[0] = -1.0;
	} else {
		s->set[set].pary[s->set[set].npoints].cc[0] = col[0];
		s->set[set].pary[s->set[set].npoints].cc[1] = col[1];
		s->set[set].pary[s->set[set].npoints].cc[2] = col[2];
	}

	s->set[set].pary[s->set[set].npoints].last = last;
	s->set[set].npoints++;

	return s->set[set].npoints-1;
}

/* Add a verticy (default natural color from pos.) */
/* Return the index number */
static int add_vertex(vrml *s, int set, double pos[3]) {
	return add_col_vertex_l(s, set, pos, NULL, 0);
}

/* Add a verticy with color */
/* Retun the index number */
static int add_col_vertex(vrml *s, int set, double pos[3], double col[3]) {
	return add_col_vertex_l(s, set, pos, col, 0);
}

/* Turn the last added vertex into the last vertex of the line */
static void make_last_vertex(vrml *s, int set) {

	if (set < 0 || set > 9)
		error("vrml make_last_vertex set %d out of range",set);

	if (s->set[set].npoints <= 0)
		warning("vrml plot: tried to set last point with no points added!\n");
	else
		s->set[set].pary[s->set[set].npoints-1].last = 1;
}

/* Turn all the vertexes into a set of points */
static void make_points(vrml *s, int set) {
	double xyz[3];
	int i, j;

	if (set < 0 || set > 9)
		error("vrml make_points set %d out of range",set);

	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"\n");
		fprintf(s->fp,"    # Points\n");
		fprintf(s->fp,"    Shape {\n");
		fprintf(s->fp,"      geometry PointSet { \n");
		fprintf(s->fp,"        coord Coordinate { \n");
		fprintf(s->fp,"          point [\n");
	} else {
		fprintf(s->fp,"\n");
		fprintf(s->fp,"    <!-- Points -->\n");
		fprintf(s->fp,"    <Shape>\n");
		fprintf(s->fp,"      <PointSet>\n");
		fprintf(s->fp,"        <Coordinate point ='\n");
	}

	for (i = 0; i < s->set[set].npoints; i++) {
		cs2xyz(s, xyz, s->set[set].pary[i].pp);
		if (s->fmt == fmt_vrml)
			fprintf(s->fp,"            %f %f %f,\n", xyz[0], xyz[1], xyz[2]); 
		else
			fprintf(s->fp,"          %f %f %f\n", xyz[0], xyz[1], xyz[2]); 
	}
	
	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"          ]\n");
		fprintf(s->fp,"        }\n");
	} else {
		fprintf(s->fp,"        '></Coordinate>\n");
	}

	/* Color */
	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"        color Color {\n");
		fprintf(s->fp,"          color [			# RGB colors of each vertex\n");
	} else {
		fprintf(s->fp,"        <Color color='\n");
	}

	for (i = 0; i < s->set[set].npoints; i++) {
		double rgb[3], Lab[3];

		if (s->set[set].pary[i].cc[0] < 0.0) {
			Lab[0] = s->set[set].pary[i].pp[0];
			Lab[1] = s->set[set].pary[i].pp[1];
			Lab[2] = s->set[set].pary[i].pp[2];
			if (s->ispace == vrml_rgb)			/* RGB */
				icmCpy3(rgb, Lab);
			else if (s->ispace == vrml_xyz)		/* XYZ */
				s->XYZ2RGB(s, rgb, Lab);
			else								/* Lab */
				s->Lab2RGB(s, rgb, Lab);	
		} else {
			icmCpy3(rgb, s->set[set].pary[i].cc);
		}
		if (s->fmt == fmt_vrml)
			fprintf(s->fp,"            %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
		else
			fprintf(s->fp,"          %f %f %f\n", rgb[0], rgb[1], rgb[2]);
	}
	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"          ] \n");
		fprintf(s->fp,"        }\n");
	} else {
		fprintf(s->fp,"        '></Color>\n");
	}
	/* End color */

	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"      }\n");
		fprintf(s->fp,"    } # end shape\n");
	} else {
		fprintf(s->fp,"      </PointSet>\n");
		fprintf(s->fp,"    </Shape>\n");
	}
}

/* Convert the verticies to lines, ppset verticies per line (or .last flag) */
static void make_lines(vrml *s, int set, int ppset) {
	double xyz[3];
	int i, j;

	if (set < 0 || set > 9)
		error("vrml make_lines set %d out of range",set);

	/* - - - - - - - - - - - - */
	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"\n");
		fprintf(s->fp,"    # Lines\n");
		fprintf(s->fp,"    Shape {\n");
		fprintf(s->fp,"      geometry IndexedLineSet { \n");
		fprintf(s->fp,"        coord Coordinate { \n");
		fprintf(s->fp,"          point [\n");

		for (i = 0; i < s->set[set].npoints; i++) {
			cs2xyz(s, xyz, s->set[set].pary[i].pp);
			fprintf(s->fp,"            %f %f %f,\n", xyz[0], xyz[1], xyz[2]); 
		}
	
		fprintf(s->fp,"          ]\n");
		fprintf(s->fp,"        }\n");
		fprintf(s->fp,"        coordIndex [\n");

		for (i = 0; i < s->set[set].npoints;) {
			fprintf(s->fp,"          ");
			for (j = 0; i < s->set[set].npoints  && j < ppset; j++) {
				fprintf(s->fp,"          %d, ", i++);
				if (s->set[set].pary[i-1].last != 0)
					break;
			}
			fprintf(s->fp,"          -1,\n");
		}
		fprintf(s->fp,"        ]\n");

		/* Color */
		fprintf(s->fp,"        colorPerVertex TRUE\n");
		fprintf(s->fp,"        color Color {\n");
		fprintf(s->fp,"          color [			# RGB colors of each vertex\n");

		for (i = 0; i < s->set[set].npoints; i++) {
			double rgb[3], Lab[3];

			if (s->set[set].pary[i].cc[0] < 0.0) {
				Lab[0] = s->set[set].pary[i].pp[0];
				Lab[1] = s->set[set].pary[i].pp[1];
				Lab[2] = s->set[set].pary[i].pp[2];
				if (s->ispace == vrml_rgb)			/* RGB */
					icmCpy3(rgb, Lab);
				else if (s->ispace == vrml_xyz)		/* XYZ */
					s->XYZ2RGB(s, rgb, Lab);
				else								/* Lab */
					s->Lab2RGB(s, rgb, Lab);	
			} else {
				icmCpy3(rgb, s->set[set].pary[i].cc);
			}
			fprintf(s->fp,"            %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
		}
		fprintf(s->fp,"          ] \n");
		fprintf(s->fp,"        }\n");
		/* End color */

		fprintf(s->fp,"      }\n");
		fprintf(s->fp,"    } # end shape\n");

	/* - - - - - - - - - - - - */
	} else {	/* x3d */
		fprintf(s->fp,"\n");
		fprintf(s->fp,"    <!-- Lines -->\n");
		fprintf(s->fp,"    <Shape>\n");
		fprintf(s->fp,"      <IndexedLineSet\n");

		fprintf(s->fp,"        colorPerVertex='true'\n");

		/* Indexes */
		fprintf(s->fp,"        coordIndex='\n");
		for (i = 0; i < s->set[set].npoints;) {
			fprintf(s->fp,"          ");
			for (j = 0; i < s->set[set].npoints  && j < ppset; j++) {
				fprintf(s->fp,"          %d ", i++);
				if (s->set[set].pary[i-1].last != 0)
					break;
			}
			fprintf(s->fp,"          -1\n");
		}
		fprintf(s->fp,"        '\n");
		fprintf(s->fp,"        >	<!-- CoordIndex -->\n");

		/* Coordinates */
		fprintf(s->fp,"        <Coordinate point='\n");
		for (i = 0; i < s->set[set].npoints; i++) {
			cs2xyz(s, xyz, s->set[set].pary[i].pp);
			fprintf(s->fp,"          %f %f %f\n", xyz[0], xyz[1], xyz[2]); 
		}
		fprintf(s->fp,"        '></Coordinate>\n");

		/* Color */
		fprintf(s->fp,"        <Color color='\n");
		for (i = 0; i < s->set[set].npoints; i++) {
			double rgb[3], Lab[3];
	
			if (s->set[set].pary[i].cc[0] < 0.0) {
				Lab[0] = s->set[set].pary[i].pp[0];
				Lab[1] = s->set[set].pary[i].pp[1];
				Lab[2] = s->set[set].pary[i].pp[2];
				if (s->ispace == vrml_rgb)			/* RGB */
					icmCpy3(rgb, Lab);
				else if (s->ispace == vrml_xyz)		/* XYZ */
					s->XYZ2RGB(s, rgb, Lab);
				else								/* Lab */
					s->Lab2RGB(s, rgb, Lab);	
			} else {
				icmCpy3(rgb, s->set[set].pary[i].cc);
			}
			fprintf(s->fp,"          %f %f %f\n", rgb[0], rgb[1], rgb[2]);
		}
		fprintf(s->fp,"        '></Color>\n");
		fprintf(s->fp,"      </IndexedLineSet>\n");
		fprintf(s->fp,"    </Shape>\n");
	}
}

/* Convert the verticies to lines, triangles or quads */
static void make_line_tri_quad(
vrml *s,
int set,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc == NULL or cc[0] < 0.0 */
				/* for previously set or vertex or natural color */
) {
	int i, j, nverts, ix;
	int v[3];
	int lines = 0;

	if (set < 0 || set > 9)
		error("vrml make_line_tri_quad set %d out of range",set);

	if (s->set[set].npoints > 0
	 && s->set[set].ntrqu > 0
	 && s->set[set].tqary[0].ix[2] < 0)	/* First is a line */
		lines = 1;		/* Assume all are lines */

	if (cc != NULL && cc[0] >= 0.0) {
		s->set[set].ppoly = 1;			/* Per poligon color */
	}

	/* - - - - - - - - - - - - */
	if (s->fmt == fmt_vrml) {
		if (lines) {
			fprintf(s->fp,"    # Lines\n");
		} else {
			fprintf(s->fp,"    # Triangles and Quads\n");
		}
	
		fprintf(s->fp,"      Shape { \n");
		if (lines)
			fprintf(s->fp,"        geometry IndexedLineSet {\n");
		else {
			fprintf(s->fp,"        geometry IndexedFaceSet {\n");
			fprintf(s->fp,"          ccw FALSE\n");
			fprintf(s->fp,"          convex TRUE\n");
			if (trans > 0.0)
				fprintf(s->fp,"          solid FALSE\n");
			else
				fprintf(s->fp,"          solid TRUE\n");
		}
		fprintf(s->fp,"\n");
		fprintf(s->fp,"          coord Coordinate { \n");
		fprintf(s->fp,"            point [			# Verticy coordinates\n");
	
		/* Spit out the point values, in order. */
		for (i = 0; i < s->set[set].npoints; i++) {
			double xyz[3];
			cs2xyz(s, xyz, s->set[set].pary[i].pp);
			fprintf(s->fp,"              %f %f %f,\n", xyz[0], xyz[1], xyz[2]);
		}
		fprintf(s->fp,"            ]\n");
		fprintf(s->fp,"          }\n");
		fprintf(s->fp,"\n");
		fprintf(s->fp,"          coordIndex [ 		# Indexes of %s Verticies \n",
		                                                  lines ? "line" : "polygon");
	
		/* Spit out the lines/triangles/quads */
		for (i = 0; i < s->set[set].ntrqu; i++) {
			if (s->set[set].tqary[i].ix[2] < 0)		/* Line */
				fprintf(s->fp,"            %d, %d, -1\n", s->set[set].tqary[i].ix[0],
				                                            s->set[set].tqary[i].ix[1]);
			else if (s->set[set].tqary[i].ix[3] < 0)		/* Triangle */
				fprintf(s->fp,"            %d, %d, %d, -1\n", s->set[set].tqary[i].ix[0],
				                                                s->set[set].tqary[i].ix[1],
				                                                s->set[set].tqary[i].ix[2]);
			else							/* Quad */
				fprintf(s->fp,"            %d, %d, %d, %d, -1\n", s->set[set].tqary[i].ix[0],
				                                                    s->set[set].tqary[i].ix[1],
				                                                    s->set[set].tqary[i].ix[2],
				                                                    s->set[set].tqary[i].ix[3]);
		}
	
		fprintf(s->fp,"          ]\n");
		fprintf(s->fp,"\n");
	
		if (s->set[set].ppoly) {
			fprintf(s->fp,"          colorPerVertex FALSE\n");
			fprintf(s->fp,"          color Color {\n");
			fprintf(s->fp,"          color [			# RGB colors of each line/tri/quad\n");
		
			/* Spit out the colors for each line/tri/quad */
			for (i = 0; i < s->set[set].ntrqu; i++) {
				double out[3];
				double rgb[3];
		
				/* Use line/patch overall supplied color */
				if (cc != NULL && cc[0] >= 0.0) {
					fprintf(s->fp,"            %f %f %f,\n", cc[0], cc[1], cc[2]);
		
				/* Use per line/tri/quad color */
				} else if (s->set[set].tqary[i].cc[0] >= 0.0) {
					fprintf(s->fp,"            %f %f %f,\n", s->set[set].tqary[i].cc[0],
					                                           s->set[set].tqary[i].cc[1],
					                                           s->set[set].tqary[i].cc[2]);
		
				/* Hmm. We can only have all per vertex or all per polygon - */
				/* use natural color of first vertex. */
				} else {
					int vx = s->set[set].tqary[i].ix[0];
					if (s->ispace == vrml_rgb)		/* RGB */
						icmCpy3(rgb, s->set[set].pary[vx].pp);
					else if (s->ispace == vrml_xyz)		/* XYZ */
						s->XYZ2RGB(s, rgb, s->set[set].pary[vx].pp);
					else								/* Lab */
						s->Lab2RGB(s, rgb, s->set[set].pary[vx].pp);	
					fprintf(s->fp,"            %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
				}
			}
			fprintf(s->fp,"            ] \n");
			fprintf(s->fp,"          }\n");
	
		/* Per vertex color */
		} else {
			fprintf(s->fp,"          colorPerVertex TRUE\n");
			fprintf(s->fp,"          color Color {\n");
			fprintf(s->fp,"          color [			# RGB colors of each vertex\n");
		
			/* Spit out the colors for each vertex */
			for (i = 0; i < s->set[set].npoints; i++) {
				double out[3];
				double rgb[3];
		
				/* Use vertex color */
				if (s->set[set].pary[i].cc[0] >= 0.0) {
					fprintf(s->fp,"            %f %f %f,\n",
						  s->set[set].pary[i].cc[0], s->set[set].pary[i].cc[1], s->set[set].pary[i].cc[2]);
		
				/* Use natural color of vertex */
				} else {
					if (s->ispace == vrml_rgb)		/* RGB */
						icmCpy3(rgb, s->set[set].pary[i].pp);
					else if (s->ispace == vrml_xyz)		/* XYZ */
						s->XYZ2RGB(s, rgb, s->set[set].pary[i].pp);
					else								/* Lab */
						s->Lab2RGB(s, rgb, s->set[set].pary[i].pp);	
		
					fprintf(s->fp,"            %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
				}
			}
			fprintf(s->fp,"            ] \n");
			fprintf(s->fp,"          }\n");
		}
	
		fprintf(s->fp,"        }\n");
		fprintf(s->fp,"        appearance Appearance { \n");
		fprintf(s->fp,"          material Material {\n");
		fprintf(s->fp,"            shininess 0.95\n");
		fprintf(s->fp,"            specularColor .6 .6 .6\n");
		if (trans > 0.0)
			fprintf(s->fp,"            transparency %f\n",trans);
		fprintf(s->fp,"          }\n");
		fprintf(s->fp,"        }\n");					 /* } */
		fprintf(s->fp,"      }	# end Shape\n");

	/* - - - - - - - - - - - - */
	/* x3d */
	} else {
		if (lines) {
			fprintf(s->fp,"    <!-- Lines -->\n");
		} else {
			fprintf(s->fp,"    <!-- Triangles and Quads -->\n");
		}
	
		fprintf(s->fp,"      <Shape>\n");
		if (lines)
			fprintf(s->fp,"        <IndexedLineSet\n");
		else {
			fprintf(s->fp,"        <IndexedFaceSet\n");
			fprintf(s->fp,"          convex='true'\n");
			fprintf(s->fp,"          ccw='false'\n");
			if (trans > 0.0)
				fprintf(s->fp,"          solid='false'\n");
			else
				fprintf(s->fp,"          solid='true'\n");
		}
		if (s->set[set].ppoly)
			fprintf(s->fp,"          colorPerVertex='false'\n");
		else
			fprintf(s->fp,"          colorPerVertex='true'\n");

		/* Indexes */
		fprintf(s->fp,"          coordIndex='\n");
		for (i = 0; i < s->set[set].ntrqu; i++) {
			fprintf(s->fp,"           ");
			for (j = 0; j < 4; j++) {
				if (s->set[set].tqary[i].ix[j] < 0)
					break;
				fprintf(s->fp," %d", s->set[set].tqary[i].ix[j]);
			}
			fprintf(s->fp," -1\n");
		}
		fprintf(s->fp,"          '>\n");

#ifdef NEVER
		if (s->set[set].ppoly) {
			/* colorIndex field is necessary for colorPerVertex=true ? */
			fprintf(s->fp,"          colorIndex='\n");
			for (i = 0; i < s->set[set].ntrqu; i++) {
				fprintf(s->fp,"            %d\n",i);
			}
			fprintf(s->fp,"            -1\n");
			fprintf(s->fp,"          '>\n");
		}
#endif

		/* Coordinates */
		fprintf(s->fp,"\n");
		fprintf(s->fp,"          <Coordinate point='\n");
		for (i = 0; i < s->set[set].npoints; i++) {
			double xyz[3];
			cs2xyz(s, xyz, s->set[set].pary[i].pp);
			fprintf(s->fp,"            %f %f %f\n", xyz[0], xyz[1], xyz[2]);
		}
		fprintf(s->fp,"          '></Coordinate>\n");

		/* Color */
		fprintf(s->fp,"\n");
		fprintf(s->fp,"          <Color color='\n");
		
		/* Per poligon color */
		if (s->set[set].ppoly) {
			/* Spit out the colors for each line/tri/quad */
			for (i = 0; i < s->set[set].ntrqu; i++) {
				double out[3];
				double rgb[3];
		
				/* Use line/patch overall supplied color */
				if (cc != NULL && cc[0] >= 0.0) {
					fprintf(s->fp,"            %f %f %f\n", cc[0], cc[1], cc[2]);
		
				/* Use per line/tri/quad color */
				} else if (s->set[set].tqary[i].cc[0] >= 0.0) {
					fprintf(s->fp,"            %f %f %f\n", s->set[set].tqary[i].cc[0],
					                                           s->set[set].tqary[i].cc[1],
					                                           s->set[set].tqary[i].cc[2]);
		
				/* Hmm. We can only have all per vertex or all per polygon - */
				/* use natural color of first vertex. */
				} else {
					int vx = s->set[set].tqary[i].ix[0];
					if (s->ispace == vrml_rgb)		/* RGB */
						icmCpy3(rgb, s->set[set].pary[vx].pp);
					else if (s->ispace == vrml_xyz)		/* XYZ */
						s->XYZ2RGB(s, rgb, s->set[set].pary[vx].pp);
					else								/* Lab */
						s->Lab2RGB(s, rgb, s->set[set].pary[vx].pp);	
					fprintf(s->fp,"            %f %f %f\n", rgb[0], rgb[1], rgb[2]);
				}
			}
	
		/* Per vertex color */
		} else {
			/* Spit out the colors for each vertex */
			for (i = 0; i < s->set[set].npoints; i++) {
				double out[3];
				double rgb[3];
		
				/* Use vertex color */
				if (s->set[set].pary[i].cc[0] >= 0.0) {
					fprintf(s->fp,"            %f %f %f\n",
						  s->set[set].pary[i].cc[0], s->set[set].pary[i].cc[1], s->set[set].pary[i].cc[2]);
		
				/* Use natural color of vertex */
				} else {
					if (s->ispace == vrml_rgb)		/* RGB */
						icmCpy3(rgb, s->set[set].pary[i].pp);
					else if (s->ispace == vrml_xyz)		/* XYZ */
						s->XYZ2RGB(s, rgb, s->set[set].pary[i].pp);
					else								/* Lab */
						s->Lab2RGB(s, rgb, s->set[set].pary[i].pp);	
		
					fprintf(s->fp,"            %f %f %f\n", rgb[0], rgb[1], rgb[2]);
				}
			}
		}

		fprintf(s->fp,"          '></Color>\n");
		if (lines)
			fprintf(s->fp,"        </IndexedLineSet>\n");
		else 
			fprintf(s->fp,"        </IndexedFaceSet>\n");
		fprintf(s->fp,"        <Appearance>\n");
		fprintf(s->fp,"          <Material shininess='0.95'\n");
		fprintf(s->fp,"                    specularColor='.6 .6 .6'\n");
		if (trans > 0.0)
			fprintf(s->fp,"                    transparency='%f'></Material>\n", trans);
		else
			fprintf(s->fp,"                    ></Material>\n");
		/* Hack to workaround bugs in x3dom trasparency */
		if (s->fmt == fmt_x3dom && trans > 0.0)	
			fprintf(s->fp,"          <DepthMode readOnly='true'></depthMode>\n");
		fprintf(s->fp,"        </Appearance>\n");
		fprintf(s->fp,"      </Shape>\n");
	}
}

/* Output lines using per vertex or per line colors. */
static void make_lines_vc(
vrml *s,
int set,
double trans	/* Transparency level */
) {
	make_line_tri_quad(s, set, trans, NULL);
}

/* Output lines with overall per line color */
/* col == NULL or col[0] < 0.0 not to set overall color */
static void make_lines_cc(
vrml *s,
int set,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc == NULL or cc[0] < 0.0 for per vertex color */
) {
	make_line_tri_quad(s, set, trans, cc);
}

/* Add a line defined by vertex index, and set per line color. */
/* col == NULL or col[0] < 0.0 not to set per line color */
static void add_col_line(vrml *s, int set, int ix[2], double col[3]) {

	if (set < 0 || set > 9)
		error("vrml add_col_line set %d out of range",set);

	if (s->set[set].ntrqu >= s->set[set].taloc) {
		s->set[set].taloc = (s->set[set].taloc + 10) * 2;
		if (s->set[set].tqary == NULL)
			s->set[set].tqary = malloc(s->set[set].taloc * sizeof(struct vrml_triquad));
		else
			s->set[set].tqary = realloc(s->set[set].tqary, s->set[set].taloc * sizeof(struct vrml_triquad));

		if (s->set[set].tqary == NULL)
			error("VRML malloc failed at count %d\n",s->set[set].taloc);
	}
	s->set[set].tqary[s->set[set].ntrqu].ix[0] = ix[0];
	s->set[set].tqary[s->set[set].ntrqu].ix[1] = ix[1];
	s->set[set].tqary[s->set[set].ntrqu].ix[2] = -1;
	s->set[set].tqary[s->set[set].ntrqu].ix[3] = -1;
	if (col != NULL && col[0] >= 0.0) {
		icmCpy3(s->set[set].tqary[s->set[set].ntrqu].cc, col);
		s->set[set].ppoly = 1;
	}
	s->set[set].ntrqu++;
}

/* Add a line defined by vertex indexes using per vertex color */
static void add_line(vrml *s, int set, int ix[2]) {
	add_col_line(s, set, ix, NULL);
}

/* Output triangles using per vertex colors. */
static void make_triangles_vc(
vrml *s,
int set,
double trans	/* Transparency level */
) {
	make_line_tri_quad(s, set, trans, NULL);
}

/* Output triangles with overall per line color */
/* col == NULL or col[0] < 0.0 not to set overall color */
static void make_triangles(
vrml *s,
int set,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc == NULL or cc[0] < 0.0 for natural color */
) {
	make_line_tri_quad(s, set, trans, cc);
}

/* Add a triangle defined by vertex indexes, and set per triangle color. */
/* col == NULL or col[0] < 0.0 not to set per line color */
static void add_col_triangle(vrml *s, int set, int ix[3], double col[3]) {

	if (set < 0 || set > 9)
		error("vrml add_col_triangle set %d out of range",set);

	if (s->set[set].ntrqu >= s->set[set].taloc) {
		s->set[set].taloc = (s->set[set].taloc + 10) * 2;
		if (s->set[set].tqary == NULL)
			s->set[set].tqary = malloc(s->set[set].taloc * sizeof(struct vrml_triquad));
		else
			s->set[set].tqary = realloc(s->set[set].tqary, s->set[set].taloc * sizeof(struct vrml_triquad));

		if (s->set[set].tqary == NULL)
			error("VRML malloc failed at count %d\n",s->set[set].taloc);
	}
	s->set[set].tqary[s->set[set].ntrqu].ix[0] = ix[0];
	s->set[set].tqary[s->set[set].ntrqu].ix[1] = ix[1];
	s->set[set].tqary[s->set[set].ntrqu].ix[2] = ix[2];
	s->set[set].tqary[s->set[set].ntrqu].ix[3] = -1;
	if (col != NULL && col[0] >= 0.0) {
		icmCpy3(s->set[set].tqary[s->set[set].ntrqu].cc, col);
		s->set[set].ppoly = 1;
	}
	s->set[set].ntrqu++;
}

/* Add a triangles defined by vertex indexes using per vertex color */
static void add_triangle(vrml *s, int set, int ix[3]) {
	add_col_triangle(s, set, ix, NULL);
}

/* Convert the verticies to quads with vertex color */
static void make_quads_vc(
vrml *s,
int set,
double trans	/* Transparency level */
) {
	make_line_tri_quad(s, set, trans, NULL);
}

/* Convert the verticies to quads with color */
static void make_quads(
vrml *s,
int set,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc == NULL or cc[0] < 0.0 for natural color */
) {
	make_line_tri_quad(s, set, trans, cc);
}

/* Add a quad defined by vertex indexes, and set per quad color. */
/* col == NULL or col[0] < 0.0 not to set per quad color */
static void add_col_quad(vrml *s, int set, int ix[4], double col[3]) {

	if (set < 0 || set > 9)
		error("vrml add_quad set %d out of range",set);

	if (s->set[set].ntrqu >= s->set[set].taloc) {
		s->set[set].taloc = (s->set[set].taloc + 10) * 2;
		if (s->set[set].tqary == NULL)
			s->set[set].tqary = malloc(s->set[set].taloc * sizeof(struct vrml_triquad));
		else
			s->set[set].tqary = realloc(s->set[set].tqary, s->set[set].taloc * sizeof(struct vrml_triquad));

		if (s->set[set].tqary == NULL)
			error("VRML malloc failed at count %d\n",s->set[set].taloc);
	}
	s->set[set].tqary[s->set[set].ntrqu].ix[0] = ix[0];
	s->set[set].tqary[s->set[set].ntrqu].ix[1] = ix[1];
	s->set[set].tqary[s->set[set].ntrqu].ix[2] = ix[2];
	s->set[set].tqary[s->set[set].ntrqu].ix[3] = ix[3];
	if (col != NULL && col[0] >= 0.0) {
		icmCpy3(s->set[set].tqary[s->set[set].ntrqu].cc, col);
		s->set[set].ppoly = 1;
	}
	s->set[set].ntrqu++;
}

/* Add a quad */
static void add_quad(vrml *s, int set, int ix[4]) {
	add_col_quad(s, set, ix, NULL);
}


/* Create a gamut surface solid or wireframe from the given gamut. */
/* Use the given transparency level. */
/* Display in natural colors if c[0] < 0.0, */
/* or the given color otherwise. Uses set 9 */
static void make_gamut_surface_2(
vrml *s,
gamut *g,
double trans,	/* Transparency level */
int wire,		/* Z for solid, NZ for wireframe */
double cc[3]	/* Surface color, cc == NULL or cc[0] < 0.0 for natural color */
) {
	int i, nverts, ix;
	int v[3];

	nverts = g->nverts(g);

	if (nverts == 0)
		return;

	s->start_line_set(s, 9);

	for (ix = i = 0; ix >= 0 && i < nverts; i++) {
		double out[3];

		ix = g->getvert(g, NULL, out, ix);
		s->add_vertex(s, 9, out);
	}

	g->startnexttri(g);
	while (g->getnexttri(g, v) == 0) {
		if (wire) {
			int ix[2];
			/* Only output 1 wire of two on an edge */
			if (v[0] < v[1])
				ix[0] = v[0], ix[1] = v[1];
			if (v[1] < v[2])
				ix[0] = v[1], ix[1] = v[2];
			if (v[2] < v[0])
				ix[0] = v[2], ix[1] = v[0];
			s->add_line(s, 9, ix);
		} else {
			s->add_triangle(s, 9, v);
		}
	}

	if (wire)
		s->make_lines_cc(s, 9, trans, cc);
	else
		s->make_triangles(s, 9, trans, cc);

	s->start_line_set(s, 9);
}

/* Create a gamut surface from the given gamut. */
/* Use the given transparency level. */
/* Display in natural colors if c[0] < 0.0, */
/* or the given color otherwise */
static void make_gamut_surface(
vrml *s,
gamut *g,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc == NULL or cc[0] < 0.0 for natural color */
) {
	s->make_gamut_surface_2(s, g, trans, 0, cc);
}

/* Add cusp markers from a gamut surface */
/* Use the given transparency level. */
/* Display in natural colors if c[0] < 0.0, */
/* or the given color otherwise. */
static void add_cusps(
vrml *s,
gamut *g,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc[0] < 0.0 for natural color, NULL for default */
) {
	double cusps[6][3];
	double ccolors[6][3] = {
		{ 1.0, 0.1, 0.1 },		/* Red */
		{ 1.0, 1.0, 0.1 },		/* Yellow */
		{ 0.1, 1.0, 0.1 },		/* Green */
		{ 0.1, 1.0, 1.0 },		/* Cyan */
		{ 0.1, 0.1, 1.0 },		/* Blue */
		{ 1.0, 0.1, 1.0 }		/* Magenta */
	};
	double rgb[3], xyz[3];
	double *cv = NULL;
	int i;
	int v[3];

	if (g->getcusps(g, cusps) != 0)
		return;

	for (i = 0; i < 6; i++) {
		if (cc == NULL) {
			cv = ccolors[i];
		} else if (cc[0] < 0.0) {
			if (s->ispace == vrml_rgb)			/* RGB */
				icmCpy3(rgb, cusps[i]);
			else if (s->ispace == vrml_xyz)		/* XYZ */
				s->XYZ2RGB(s, rgb, cusps[i]);
			else								/* Lab */
				s->Lab2RGB(s, rgb, cusps[i]);	

			cv = rgb;
		} else {
			cv = cc;
		}
		s->add_marker_trans(s, cusps[i], cc, trans, 2.0);
	}
}

/* Clear verticies and triangles */
static void clear(vrml *s) {
	int i;

	for (i = 0; i < 10; i++) {
		if (s->set[i].pary != NULL)
			free(s->set[i].pary);
		s->set[i].pary = NULL;
		s->set[i].npoints = s->set[i].paloc = 0;

		if (s->set[i].tqary != NULL)
			free(s->set[i].tqary);
		s->set[i].tqary = NULL;
		s->set[i].ntrqu = s->set[i].taloc = 0;
	}
}
	
/* Helper :- convert a Lab value to RGB for display purposes */
static void Lab2RGB(vrml *s, double *out, double *in) {
	double L = in[0], a = in[1], b = in[2];
	double x,y,z,fx,fy,fz;
	double R, G, B;

	/* Scale so that black is visible */
	L = L * (100 - 40.0)/100.0 + 40.0;

	/* First convert to XYZ using D50 white point */
	if (L > 8.0) {
		fy = (L + 16.0)/116.0;
		y = pow(fy,3.0);
	} else {
		y = L/903.2963058;
		fy = 7.787036979 * y + 16.0/116.0;
	}

	fx = a/500.0 + fy;
	if (fx > 24.0/116.0)
		x = pow(fx,3.0);
	else
		x = (fx - 16.0/116.0)/7.787036979;

	fz = fy - b/200.0;
	if (fz > 24.0/116.0)
		z = pow(fz,3.0);
	else
		z = (fz - 16.0/116.0)/7.787036979;

	x *= 0.9642;	/* Multiply by white point, D50 */
	y *= 1.0;
	z *= 0.8249;

	/* Now convert to sRGB values */
	R = x * 3.2410  + y * -1.5374 + z * -0.4986;
	G = x * -0.9692 + y * 1.8760  + z * 0.0416;
	B = x * 0.0556  + y * -0.2040 + z * 1.0570;

	if (R < 0.0)
		R = 0.0;
	else if (R > 1.0)
		R = 1.0;

	if (G < 0.0)
		G = 0.0;
	else if (G > 1.0)
		G = 1.0;

	if (B < 0.0)
		B = 0.0;
	else if (B > 1.0)
		B = 1.0;

	R = pow(R, 1.0/2.2);
	G = pow(G, 1.0/2.2);
	B = pow(B, 1.0/2.2);

	/* For a black background: */
//	R = R * 0.85 + 0.15;
//	G = G * 0.85 + 0.15;
//	B = B * 0.85 + 0.15;

	/* For a white background: */
//	R = R * 0.70 + 0.05;
//	G = G * 0.70 + 0.05;
//	B = B * 0.70 + 0.05;

	out[0] = R;
	out[1] = G;
	out[2] = B;
}

/* Helper :- convert an XYZ value to RGB for display purposes */
static void XYZ2RGB(vrml *s, double *out, double *in) {
	double x = in[0], y = in[1], z = in[2];
	double R, G, B;

	/* Now convert to sRGB values */
	R = x * 3.2410  + y * -1.5374 + z * -0.4986;
	G = x * -0.9692 + y * 1.8760  + z * 0.0416;
	B = x * 0.0556  + y * -0.2040 + z * 1.0570;

	if (R < 0.0)
		R = 0.0;
	else if (R > 1.0)
		R = 1.0;

	if (G < 0.0)
		G = 0.0;
	else if (G > 1.0)
		G = 1.0;

	if (B < 0.0)
		B = 0.0;
	else if (B > 1.0)
		B = 1.0;

	R = pow(R, 1.0/2.2);
	G = pow(G, 1.0/2.2);
	B = pow(B, 1.0/2.2);

	/* For a black background: */
//	R = R * 0.85 + 0.15;
//	G = G * 0.85 + 0.15;
//	B = B * 0.85 + 0.15;

	/* For a white background: */
	R = R * 0.70 + 0.05;
	G = G * 0.70 + 0.05;
	B = B * 0.70 + 0.05;

	out[0] = R;
	out[1] = G;
	out[2] = B;
}

static int do_flush(vrml *s);
static void del_vrml(vrml *s);

/* Global format */
static vrml_fmt g_fmt = fmt_uninit; 

/* Check and override the global format */
static void check_format() {
	if (g_fmt == fmt_uninit) {
		char *ev;
		g_fmt = fmt_x3dom;			/* Set global default */
		
		if ((ev = getenv("ARGYLL_3D_DISP_FORMAT")) != NULL) {
			if (stricmp(ev, "VRML") == 0
			 || stricmp(ev, "WRL") == 0)
				g_fmt = fmt_vrml;
			else if (stricmp(ev, "X3D") == 0)
				g_fmt = fmt_x3d;
			else if (stricmp(ev, "X3DOM") == 0)
				g_fmt = fmt_x3dom;
		}
	}
}

/* Return the global format file extension */
static char *ret_ext(vrml_fmt fmt) {
	if (fmt == fmt_uninit) {
		check_format();
		fmt = g_fmt;
	}

	if (fmt == fmt_x3dom)	
		return ".x3d.html";
	else if (fmt == fmt_x3d)	
		return ".x3d";
	else
		return ".wrl";
}

/* Return the global format type name */
static char *ret_format(vrml_fmt fmt) {
	if (fmt == fmt_uninit) {
		check_format();
		fmt = g_fmt;
	}

	if (fmt == fmt_x3dom)	
		return "X3DOM";
	else if (fmt == fmt_x3d)	
		return "X3D";
	else
		return "VRML";
}

/* Return this files format extension (i.e. ".wrl" */
static char *get_ext(vrml *s) {
	return ret_ext(s->fmt);
}

/* Return this files format type name */
static char *get_format(vrml *s) {
	return ret_format(s->fmt);
}

/* Constructor */
/* ispace = 0 = L*a*b* */
/* ispace = 1 = XYZ scale (range 0..1) */
/* ispace = 2 = RGB scale (range 0..1) */
vrml *new_vrml_vdist(
char *name,
int doaxes,
vrml_space ispace,
double vdist
) {
	vrml *s;

	int i, j;

	if ((s = (vrml *)calloc(1, sizeof(vrml))) == NULL) {
		warning("Malloc of vrml plot object failed");
		return NULL;
	}

	if ((s->name = (char *)malloc(strlen(name) + 10)) == NULL) {
		warning("Malloc of vrml name failed");
		free(s);
		return NULL;
	}

	s->ext                   = get_ext;
	s->format                = get_format;
	s->flush                 = do_flush;
	s->del                   = del_vrml;

	s->add_marker            = add_marker;
	s->add_marker_trans      = add_marker_trans;
	s->add_cone              = add_cone;
	s->add_text              = add_text;
	s->start_line_set        = start_line_set;
	s->add_vertex            = add_vertex;
	s->add_col_vertex        = add_col_vertex;
	s->make_last_vertex      = make_last_vertex;
	s->make_lines            = make_lines;
	s->make_points           = make_points;
	s->add_line              = add_line;
	s->add_col_line          = add_col_line;
	s->make_lines_vc         = make_lines_vc;
	s->make_lines_cc         = make_lines_cc;
	s->add_triangle          = add_triangle;
	s->add_col_triangle      = add_col_triangle;
	s->make_triangles_vc     = make_triangles_vc;
	s->make_triangles        = make_triangles;
	s->add_quad              = add_quad;
	s->add_col_quad          = add_col_quad;
	s->make_quads_vc         = make_quads_vc;
	s->make_quads            = make_quads;
	s->make_gamut_surface    = make_gamut_surface;
	s->make_gamut_surface_2  = make_gamut_surface_2;
	s->add_cusps             = add_cusps;
	s->clear                 = clear;
	s->Lab2RGB               = Lab2RGB;
	s->XYZ2RGB               = XYZ2RGB;

	if (g_fmt == fmt_uninit)
		check_format();
	s->fmt = g_fmt;			/* Use global format */

	s->ispace = ispace;

	if (s->ispace == vrml_rgb) {	/* RGB, scale 0..1 to 0..100 */
		s->scale = 100.0;
		s->off = 0.0;		

	} else if (s->ispace == vrml_xyz) {	/* XYZ, scale 0..1 to 0..100 */
		s->scale = 100.0;
		s->off = 50.0;			/* z axis offset */

	} else {					/* L*a*b*, leav 0..100 */
		s->scale = 1.0;		
		s->off = 50.0;			/* z axis offset */
	}

	/* Create filename with the right exension */
	{
		char *xl = NULL;
		strcpy(s->name, name);
		if ((xl = strrchr(s->name, '.')) != NULL) {		/* Found extension */
			/* Hmm. Could set format from extension ?? */
			if (stricmp(xl, ".wrl") != 0
			 && stricmp(xl, ".vrml") != 0
			 && stricmp(xl, ".x3d") != 0
			 && stricmp(xl, ".x3dom") != 0
			 && stricmp(xl, ".html") != 0)
				xl = NULL;		/* Don't override extension */
		}
		if (xl == NULL)
			xl = s->name + strlen(s->name);
		strcpy(xl, vrml_ext());
	}

	if ((s->fp = fopen(s->name,"w")) == NULL) {
		warning("Opening of vrml plot file '%s' for write failed",s->name);
		free(s);
		return NULL;
	}

	/* Output header and prolog */
	if (s->fmt == fmt_vrml) {
		fprintf(s->fp,"#VRML V2.0 utf8\n");
		fprintf(s->fp,"\n");
		fprintf(s->fp,"# Created by the Argyll CMS\n");
		fprintf(s->fp,"Transform {\n");
		fprintf(s->fp,"  children [\n");
		fprintf(s->fp,"\n");
		fprintf(s->fp,"    NavigationInfo {\n");
		fprintf(s->fp,"      type \"EXAMINE\"		# It's an object we examine\n");
		fprintf(s->fp,"      headlight FALSE\n");
		fprintf(s->fp,"    } # We'll add our own light\n");
		fprintf(s->fp,"\n");
		fprintf(s->fp,"    DirectionalLight {\n");
		fprintf(s->fp,"      intensity 0.4\n");
		fprintf(s->fp,"      ambientIntensity 0.05\n");
		fprintf(s->fp,"      direction 1 1 -1\n");
		fprintf(s->fp,"    }\n");
		fprintf(s->fp,"    DirectionalLight {\n");
		fprintf(s->fp,"      intensity 0.4\n");
		fprintf(s->fp,"      ambientIntensity 0.05\n");
		fprintf(s->fp,"      direction 0 -0.7 -1\n");
		fprintf(s->fp,"    }\n");
		fprintf(s->fp,"    DirectionalLight {\n");
		fprintf(s->fp,"      intensity 0.4\n");
		fprintf(s->fp,"      ambientIntensity 0.05\n");
		fprintf(s->fp,"      direction -0.7 0 -1\n");
		fprintf(s->fp,"    }\n");
		fprintf(s->fp,"    DirectionalLight {\n");
		fprintf(s->fp,"      intensity 0.4\n");
		fprintf(s->fp,"      ambientIntensity 0.05\n");
		fprintf(s->fp,"      direction -1 -1 1\n");
		fprintf(s->fp,"    }\n");
		fprintf(s->fp,"    DirectionalLight {\n");
		fprintf(s->fp,"      intensity 0.4\n");
		fprintf(s->fp,"      ambientIntensity 0.05\n");
		fprintf(s->fp,"      direction 0 0.7 1\n");
		fprintf(s->fp,"    }\n");
		fprintf(s->fp,"    DirectionalLight {\n");
		fprintf(s->fp,"      intensity 0.4\n");
		fprintf(s->fp,"      ambientIntensity 0.05\n");
		fprintf(s->fp,"      direction 0.7 0 1\n");
		fprintf(s->fp,"    }\n");
		fprintf(s->fp,"\n");
		fprintf(s->fp,"    Viewpoint {\n");
		fprintf(s->fp,"      position 0 0 %f      # Position we view from\n",vdist);
		fprintf(s->fp,"    }\n");
		fprintf(s->fp,"\n");

	} else {
		/* For some strange reason, x3dom can't handle some tags that don't have a */
		/* discrete closing tag. Is this XML or what ??? */

		if (s->fmt == fmt_x3dom) {
			fprintf(s->fp,"<!DOCTYPE html>\n");
			fprintf(s->fp,"<html>\n");
			fprintf(s->fp,"  <head>\n");
			fprintf(s->fp,"    <meta http-equiv='Content-Type' content='text/html;charset=utf-8'></meta>\n");
			fprintf(s->fp,"    <link rel='stylesheet' type='text/css' href='x3dom.css'></link> \n");
			fprintf(s->fp,"  </head>\n");
			fprintf(s->fp,"  <body>\n");
			fprintf(s->fp,"    <noscript><p>Please enable JavaScript</p></noscript>\n");
			fprintf(s->fp,"    <script type='text/javascript' src='x3dom.js'> </script> \n");
			fprintf(s->fp,"    <x3d style='width: 100%%; height: 70%%;'\n");
			fprintf(s->fp,"         x='0px' y='0px' width='100%%' height='70%%'\n");
			fprintf(s->fp,"         id='someUniqueId' showStat='false' showLog='false'>\n");
		} else {
			fprintf(s->fp,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
			fprintf(s->fp,"<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.0//EN\" \"http://www.web3d.org/specifications/x3d-3.0.dtd\">\n");
			fprintf(s->fp,"<X3D xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' profile='Immersive' version='3.0' xsd:noNamespaceSchemaLocation='http://www.web3d.org/specifications/x3d-3.0.xsd'>\n");
		}
		fprintf(s->fp,"  <Scene DEF='scene'>\n");
		if (s->fmt == fmt_x3dom) {
			/* To match other renderers, we turn Gamma encoded output off :-( :-( :-( */
			fprintf(s->fp,"    <Environment gammaCorrectionDefault='none'></Environment>\n");
		}
		fprintf(s->fp,"    <Background groundColor='0 0 0' skyColor='0 0 0'></Background>\n");
		fprintf(s->fp,"    <Transform>\n");
		fprintf(s->fp,"      <NavigationInfo type='\"EXAMINE\"' headlight='false'></NavigationInfo>\n");
		fprintf(s->fp,"      <DirectionalLight\n");
		fprintf(s->fp,"        intensity='0.4'\n");
		fprintf(s->fp,"        ambientIntensity='0.05'\n");
		fprintf(s->fp,"        direction='1 1 -1'\n");
		fprintf(s->fp,"      ></DirectionalLight>\n");
		fprintf(s->fp,"      <DirectionalLight\n");
		fprintf(s->fp,"        intensity='0.3'\n");
		fprintf(s->fp,"        ambientIntensity='0.05'\n");
		fprintf(s->fp,"        direction='0 -0.7 -1'\n");
		fprintf(s->fp,"      ></DirectionalLight>\n");
		fprintf(s->fp,"      <DirectionalLight\n");
		fprintf(s->fp,"        intensity='0.4'\n");
		fprintf(s->fp,"        ambientIntensity='0.05'\n");
		fprintf(s->fp,"        direction='-0.7 0 -1'\n");
		fprintf(s->fp,"      ></DirectionalLight>\n");
		fprintf(s->fp,"      <DirectionalLight\n");
		fprintf(s->fp,"        intensity='0.4'\n");
		fprintf(s->fp,"        ambientIntensity='0.05'\n");
		fprintf(s->fp,"        direction='-1 -1 1'\n");
		fprintf(s->fp,"      ></DirectionalLight>\n");
		fprintf(s->fp,"      <DirectionalLight\n");
		fprintf(s->fp,"        intensity='0.3'\n");
		fprintf(s->fp,"        ambientIntensity='0.05'\n");
		fprintf(s->fp,"        direction='0 0.7 1'\n");
		fprintf(s->fp,"      ></DirectionalLight>\n");
		fprintf(s->fp,"      <DirectionalLight\n");
		fprintf(s->fp,"        intensity='0.4'\n");
		fprintf(s->fp,"        ambientIntensity='0.05'\n");
		fprintf(s->fp,"        direction='0.7 0 1'\n");
		fprintf(s->fp,"      ></DirectionalLight>\n");
		fprintf(s->fp,"      <Viewpoint position='0 0 %f'></Viewpoint>\n",vdist);
	}

	if (doaxes != 0) {
		/* Axes definition */
		struct {
			char *label;
			double x, y, z;			/* == a,b,L or Y,Z,X */
			double wx, wy, wz;
			double r, g, b;
		} axes[3][6] = {
			{	/* Box coords are center and size: */
				{ "L",  0, 0,   50,  2, 2, 100,  .7, .7, .7 },	/* L axis */
				{ "+a", 50, 0,  0,  100, 2, 2,   1,  0,  0 },	/* +a (red) axis */
				{ "-b", 0, -50, 0,  2, 100, 2,   0,  0,  1 },	/* -b (blue) axis */
				{ "-a", -50, 0, 0,  100, 2, 2,   0,  1,  0 },	/* -a (green) axis */
				{ "+b", 0,  50, 0,  2, 100, 2,   1,  1,  0 },	/* +b (yellow) axis */
				{ NULL },
			}, {
				{ "X",  0,  0, 50,  2, 2, 100,  .7, .7, .7 },	/* X axis */
				{ "Y", 50,  0,  0,  100, 2, 2,   1,  0,  0 },	/* Y (red) axis */
				{ "Z",  0, 50,  0,  2, 100, 2,   0,  0,  1 },	/* Z (blue) axis */
				{ NULL },
			}, {
				{ "R", 50,  0,  0,  100, 2, 2,   1,  0,  0 },	/* R (red) */
				{ "G",  0, 50,  0,  2, 100, 2,   0,  1,  0 },	/* G (green) axis */
				{ "B",  0,  0, 50,  2, 2, 100,   0,  0,  1 },	/* B (blue) axis */
				{ NULL },
			}
		};

		if (s->ispace == 2) {
			j = 2;
			if (s->fmt == fmt_vrml)
				fprintf(s->fp,"    # RGB axes as boxes:\n");
			else
				fprintf(s->fp,"      <!--RGB axes as boxes -->\n");
		} else if (s->ispace == 1) {
			j = 1;
			if (s->fmt == fmt_vrml)
				fprintf(s->fp,"    # XYZ axes as boxes:\n");
			else
				fprintf(s->fp,"      <!--RGB axes as boxes -->\n");
		} else {
			j = 0;
			if (s->fmt == fmt_vrml)
				fprintf(s->fp,"    # Lab axes as boxes:\n");
			else
				fprintf(s->fp,"      <!--Lab axes as boxes -->\n");
		}
		for (i = 0; ; i++) {
			double toff[3] = { -3.0, -2.0, 0 };

			if (axes[j][i].label == NULL)
				break;

			if (s->fmt == fmt_vrml) {
				fprintf(s->fp,"\tTransform { translation %f %f %f\n",
				                              axes[j][i].x, axes[j][i].y, axes[j][i].z - s->off);
				fprintf(s->fp,"\t\tchildren [\n");
				fprintf(s->fp,"\t\t\tShape {\n");
				fprintf(s->fp,"\t\t\t\tgeometry Box { size %f %f %f }\n",
				                              axes[j][i].wx, axes[j][i].wy, axes[j][i].wz);
				fprintf(s->fp,"\t\t\t\tappearance Appearance {\n");
				fprintf(s->fp,"\t\t\t\t\tmaterial Material { diffuseColor %f %f %f }\n",
				                              axes[j][i].r, axes[j][i].g, axes[j][i].b);
				fprintf(s->fp,"\t\t\t\t}\n");
				fprintf(s->fp,"\t\t\t}\n");
				fprintf(s->fp,"\t\t]\n");
				fprintf(s->fp,"\t}\n");
			} else {
				fprintf(s->fp,"      <Transform translation='%f %f %f'>\n",
				                              axes[j][i].x, axes[j][i].y, axes[j][i].z - s->off);
				fprintf(s->fp,"        <Shape>\n");
				fprintf(s->fp,"          <Appearance>\n");
				fprintf(s->fp,"            <Material diffuseColor='%f %f %f'></Material>\n",
				                              axes[j][i].r, axes[j][i].g, axes[j][i].b);
				fprintf(s->fp,"          </Appearance>\n");
				fprintf(s->fp,"          <Box size='%f %f %f'></Box>\n",
				                              axes[j][i].wx, axes[j][i].wy, axes[j][i].wz);
				fprintf(s->fp,"        </Shape>\n");
				fprintf(s->fp,"      </Transform>\n");
			}

			if (fabs(axes[j][i].x) > fabs(axes[j][i].y) && fabs(axes[j][i].x) > fabs(axes[j][i].z)) {
				if (axes[j][i].x > 0.0)
					toff[0] += axes[j][i].x + 0.5 * axes[j][i].wx + 5.0;
				else
					toff[0] += axes[j][i].x - 0.5 * axes[j][i].wx - 5.0;
			} else if (fabs(axes[j][i].y) > fabs(axes[j][i].x) && fabs(axes[j][i].y) > fabs(axes[j][i].z)) {
				if (axes[j][i].y > 0.0)
					toff[1] += axes[j][i].y + 0.5 * axes[j][i].wy + 5.0;
				else
					toff[1] += axes[j][i].y - 0.5 * axes[j][i].wy - 5.0;
			} else { 
				if (axes[j][i].z > 0.0)
					toff[2] += axes[j][i].z + 0.5 * axes[j][i].wz + 5.0;
				else
					toff[2] += axes[j][i].z - 0.5 * axes[j][i].wz - 5.0;
			}

			if (s->fmt == fmt_vrml) {
				fprintf(s->fp,"\tTransform { translation %f %f %f\n",
				                                    toff[0], toff[1], toff[2] - s->off);
				fprintf(s->fp,"\t\tchildren [\n");
				fprintf(s->fp,"\t\t\tShape {\n");
				fprintf(s->fp,"\t\t\t\tgeometry Text { string [\"%s\"]\n",axes[j][i].label);
				fprintf(s->fp,"\t\t\t\t\tfontStyle FontStyle { family \"SANS\" style \"BOLD\" size %f }\n",
				                                            7.0);
				fprintf(s->fp,"\t\t\t\t\t}\n");
				fprintf(s->fp,"\t\t\t\tappearance Appearance { material Material ");
				fprintf(s->fp,"\t{ diffuseColor %f %f %f } }\n",
				                                    axes[j][i].r, axes[j][i].g, axes[j][i].b);
				fprintf(s->fp,"\t\t\t}\n");
				fprintf(s->fp,"\t\t]\n");
				fprintf(s->fp,"\t}\n");
			} else {
				fprintf(s->fp,"      <Transform translation='%f %f %f'>\n",
				                                            toff[0], toff[1], toff[2] - s->off);
				fprintf(s->fp,"        <Shape>\n");
				fprintf(s->fp,"          <Appearance>\n");
				fprintf(s->fp,"            <Material diffuseColor='%f %f %f'></Material>\n",
				                                            axes[j][i].r, axes[j][i].g, axes[j][i].b);
				fprintf(s->fp,"          </Appearance>\n");
				fprintf(s->fp,"          <Text string='\"%s\"'>\n",axes[j][i].label);
				fprintf(s->fp,"            <FontStyle family='\"SANS\"' style='BOLD' size='%f'></FontStyle>\n", 7.0);
				fprintf(s->fp,"          </Text>\n");
				fprintf(s->fp,"        </Shape>\n");
				fprintf(s->fp,"      </Transform>\n");
			}
		}
	}
	return s;
}

vrml *new_vrml(
char *name,
int doaxes,
vrml_space ispace
) {
	return new_vrml_vdist(name, doaxes, ispace, 340.0);
}

/* The X3DOM files */
unsigned char x3dom_css[] = {
#include "x3dom.css.h"
};
#define x3dom_css_len sizeof(x3dom_css)

unsigned char x3dom_js[] = {
#include "x3dom.js.h"
};
#define x3dom_js_len sizeof(x3dom_js)


/* Finish writing the file */
/* Return nz on error */
static int do_flush(vrml *s) {
	int rv = 0;

	if (!s->written) {
		
		if (s->fmt == fmt_vrml) {
			fprintf(s->fp,"\n");
			fprintf(s->fp,"  ] # end of children for world\n");
			fprintf(s->fp,"}\n");
		} else {
			fprintf(s->fp,"    </Transform>\n");
			fprintf(s->fp,"  </Scene>\n");
			if (s->fmt == fmt_x3dom) {
				fprintf(s->fp,"    </x3d>\n");
				fprintf(s->fp,"  </body>\n");
				fprintf(s->fp,"</html>\n");
			} else {
				fprintf(s->fp,"</X3D>\n");
			}
		}
	
		fflush(s->fp);
		rv = fclose(s->fp);

		/* Check that there are the x3dom files with the output file */
		if (s->fmt == fmt_x3dom) {
			char *xl, *x3name;
			struct sys_stat sbuf;
			FILE *fp;
			char *oflags = 
#if !defined(O_CREAT) && !defined(_O_CREAT)		// No O_BINARY possible
# error "Need to #include fcntl.h!"
#endif
#if defined(O_BINARY) || defined(_O_BINARY)
			"wb";
#else
			"w";
#endif

			if ((x3name = (char *)malloc(strlen(s->name) + 20)) == NULL) {
				warning("VRML: failed to malloc x3dom filename\n",rv);
				return -1;
			}

			strcpy(x3name, s->name);
			// Locate start of filename
			if ((xl = strrchr(x3name, '/')) == NULL
			 && (xl = strrchr(x3name, '\\')) == NULL
			 && (xl = strrchr(x3name, ':')) == NULL)
				xl = x3name;
			else
				xl++;

			strcpy(xl, "x3dom.css");
			if (sys_stat(x3name, &sbuf) != 0
			 || sbuf.st_size != x3dom_css_len) {
//				printf("Can't locate '%s' or wrong size\n",x3name);

				if ((fp = fopen(x3name, oflags)) == NULL) {
					warning("Opening '%s' for write failed",x3name);
					return -1;
				}
				
				if (fwrite((void *)x3dom_css, sizeof(char), x3dom_css_len, fp) != x3dom_css_len
				 || fclose(fp) != 0) {
					warning("Writing '%s'failed",x3name);
					return -1;
				}
//				printf("Written '%s' %d bytes\n",x3name,x3dom_css_len);
			}
			
			strcpy(xl, "x3dom.js");
			if (sys_stat(x3name, &sbuf) != 0
			 || sbuf.st_size != x3dom_js_len) {
//				printf("Can't locate '%s'\n",x3name);

				if ((fp = fopen(x3name, oflags)) == NULL) {
					warning("Opening '%s' for write failed",x3name);
					return -1;
				}
				
				if (fwrite((void *)x3dom_js, sizeof(char), x3dom_js_len, fp) != x3dom_js_len
				 || fclose(fp) != 0) {
					warning("Writing '%s'failed",x3name);
					return -1;
				}
//				printf("Written '%s' %d bytes\n",x3name,x3dom_js_len);
			}
			free(x3name);
		}

		s->written = 1;
	}
	return rv;
}

/* Return the global format file extension */
char *vrml_ext() {
	return ret_ext(g_fmt);
}

/* Return the global format type name */
char *vrml_format() {
	return ret_format(g_fmt);
}

/* Finish writing the file and free ourselves */
static void del_vrml(vrml *s) {
	int i, rv;

	if ((rv = do_flush(s)) != 0)
		error("VRML: Error %d closing VRML file\n",rv);

	for (i = 0; i < 10; i++) {
		if (s->set[i].pary)
			free(s->set[i].pary);
		if (s->set[i].tqary)
			free(s->set[i].tqary);
	}
	if (s->name != NULL)
		free(s->name);
    free(s);
}

