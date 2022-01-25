
/* 
 * Argyll Color Management System
 * ChromeCast Display target patch window
 *
 * Author: Graeme W. Gill
 * Date:   8/9/14
 *
 * Copyright 2013, 2014 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <string.h>
#include "copyright.h"
#include "aconfig.h"
#include "icc.h"
#include "numsup.h"
#include "cgats.h"
#include "conv.h"
#include "dispwin.h"
#include "conv.h"
#include "mongoose.h"
#include "ccast.h"
#include "ccwin.h"
#include "render.h"

#undef WRITE_PNG		/* [und] Write each test patch to "ccwin.png" */
#undef CCTEST_PATTERN	/* [und] Create dev. spatial test pattern */
#undef SEND_TEST_FILE 	/* [und] Send this file name to ChromeCast instead of pattern */
#undef DO_TIMING       	/* [und] Print rendering timing */
#undef DEBUG			/* [und] */

#define DDITHER 1		/*       0 = no dither - quantize to PNG 8 bit RGB value */
						/* [def] 1 = use error diffusion dithering with ccast quant model */
						/*       2 = use crafted 4x4 dither cell */

#define VWIDTH  1920.0	/* Video stream and display size ? */
#define VHEIGHT 1080.0

/* Definitions for cctyp_1 */
#define IWIDTH  1280.0	/* This is the native still image framebuffer size for ChromeCasts */
#define IHEIGHT 720.0

/* 2K Definitions for cctyp_Ultra ? */
#define IWIDTH_H  1920.0
#define IHEIGHT_H 1080.0

/* 4K Definitions for cctyp_Ultra ? */
#define IWIDTH_U  3840.0
#define IHEIGHT_U 2160.0

//#define STANDALONE_TEST

#ifdef DEBUG
# pragma message("######### ccwin DEBUG is defined! ########")
#define errout g_log,0
# define debug2(xx)	a1logd xx
# define debugr2(xx)	a1logd xx
# define debugrr2(xx)	a1logd xx
# define debugrr2l(lev, xx)	a1logd xx
#else
#define errout g_log,0
# define debug2(xx)													// Debug, args
# define debugr2(xx) if (p->ddebug) a1logd xx						// Run ->ddebug, args
# define debugrr2(xx) if (callback_ddebug) a1logd xx				// Run cback, args
# define debugrr2l(lev, xx) if (callback_ddebug >= lev) a1logd xx	// Run cback >= lev, args
#endif

/* ================================================================== */

/* Chromwin context and (possible) web server  */
typedef struct _chws {
	int verb;
	int ddebug;
	
	int direct;					/* End PNG directly, rather than using web server */

	struct mg_context *mg;		/* Mongoose context (if needed) */
	char *ws_url;				/* Web server URL for accessing server */

//	double hoff, voff;			/* Input position of test square */
	double x, y;				/* position of test square in pixels */
	double w, h;				/* size of test square in pixels */
	double bg[3];				/* Background color */
	int pno;					/* Index to generate a sequence of URLs */
	unsigned char *ibuf;		/* Memory image of .png file */
	size_t ilen;

	ccast *cc;					/* ChromeCast */

	/* Set a whole screen sized png image */
	int (*set)(struct _chws *p, unsigned char *ibuf, size_t ilen);

	/* Update the test patch png image */
	int (*update)(struct _chws *p, unsigned char *ibuf, size_t ilen, double *bg);

	/* Destroy ourselves */
	void (*del)(struct _chws *p);

} chws;

static void chws_del(chws *p) {

	/* delete mongoose, if we are using it */
	if (p->mg != NULL)
		mg_stop(p->mg);

	/* Delete ChromeCast */
	if (p->cc != NULL)
		p->cc->del(p->cc);

	if (p->ibuf != NULL)
		free(p->ibuf);

	if (p->ws_url != NULL)
		free(p->ws_url);

	free(p);
}

/* Set a whole screen .png (size is assumed to be large enough) */
/* Return nz on error */
static int chws_set(chws *p, unsigned char *ibuf, size_t ilen) {
	char url[200];
	double bg[3] = { 0.0, 0.0, 0.0 };

	debug2((errout,"\nUpdate png\n"));

	if (p->ibuf != NULL)
		free(p->ibuf);

	p->ibuf = ibuf;
	p->ilen = ilen;

	/* Send the PNG swatch direct */
	if (p->direct) {
		double x, y, w, h;
		/* Convert x,y,w,h to relative rather than pixel size */

		debugr2((errout,"Got x %f y %f w %f h %f\n", p->x, p->y, p->w, p->h));

		// Convert from quantized to direct loader parameters
		x = 0.0;
		y = 0.0;
		w = 10.0;		/* Scale */
		h = 10.0 * 9.0/16.0;

		debugr2((errout,"Sending direct x %f y %f w %f h %f\n", x, y, w, h));

		if (p->cc->load(p->cc, NULL, p->ibuf, p->ilen, bg, x, y, w, h)) {
			debugr2((errout,"ccwin_set direct load failed\n"));
			return 1;
		}

	/* Using web server */
	} else {

#ifdef SEND_TEST_FILE
		sprintf(url, "%s%s",p->ws_url, SEND_TEST_FILE);
#else
		sprintf(url, "%stpatch_%d.png",p->ws_url, ++p->pno); 
#endif
		if (p->cc->load(p->cc, url, NULL, 0, NULL,  0.0, 0.0, 0.0, 0.0)) {
			debugr2((errout,"ccwin_set server load failed\n"));
			return 1;
		}
	}
	return 0;
}

/* Change the test patch .png being served */
/* Return nz on error */
static int chws_update(chws *p, unsigned char *ibuf, size_t ilen, double *bg) {
	char url[200];

	debug2((errout,"\nUpdate png\n"));

	if (p->ibuf != NULL)
		free(p->ibuf);

	p->ibuf = ibuf;
	p->ilen = ilen;

	/* Send the PNG swatch direct */
	if (p->direct) {
		double x, y, w, h;
		/* Convert x,y,w,h to relative rather than pixel size */

		debugr2((errout,"Got x %f y %f w %f h %f\n", p->x, p->y, p->w, p->h));

		// Convert from quantized to direct loader parameters
		if (p->w < IWIDTH)
			x = p->x/(IWIDTH - p->w);
		else	
			x = 0.0;
		if (p->h < IHEIGHT)
			y = p->y/(IHEIGHT - p->h);
		else
			y = 0.0;
		w = p->w/(0.1 * IWIDTH);
		h = p->h/(0.1 * IWIDTH);

		debugr2((errout,"Sending direct x %f y %f w %f h %f\n", x, y, w, h));

		if (p->cc->load(p->cc, NULL, p->ibuf, p->ilen, bg, x, y, w, h)) {
			debugr2((errout,"ccwin_update direct load failed\n"));
			return 1;
		}

	/* Using web server */
	} else {

#ifdef SEND_TEST_FILE
		sprintf(url, "%s%s",p->ws_url, SEND_TEST_FILE);
#else
		sprintf(url, "%stpatch_%d.png",p->ws_url, ++p->pno); 
#endif
		if (p->cc->load(p->cc, url, NULL, 0, NULL,  0.0, 0.0, 0.0, 0.0)) {
			debugr2((errout,"ccwin_update server load failed\n"));
			return 1;
		}
	}
	return 0;
}

/* Web server event handler - return the current .png image */
static void *ccwin_ehandler(enum mg_event event,
                           struct mg_connection *conn) {
	const struct mg_request_info *request_info = mg_get_request_info(conn);
	chws *p = (chws *)mg_get_user_data(conn);
	char *cp;
	char sbuf[200];

	debugr2((errout,"ccwin_ehandler()\n"));

	if (event != MG_NEW_REQUEST) {
		return NULL;
	}

	debugr2((errout,"Event: uri = '%s'\n",request_info->uri));

#ifdef SEND_TEST_FILE
#pragma message("############################# ccwin.c SEND_TEST_FILE defined ! ##")
	return NULL;
#endif

	if (p->ibuf != NULL && p->ilen > 0
     && (cp = strrchr(request_info->uri, '.')) != NULL
	 && strcmp(cp, ".png") == 0) { 

		debugr2((errout,"Event: Loading %s\n",request_info->uri));

		debugr2((errout,"Returning current png size %d bytes\n",(int)p->ilen));
		sprintf(sbuf,
			"HTTP/1.1 200 OK\r\n"
			"Content-Type: image/png\r\n"
		    "Content-Length: %d\r\n"
		    "\r\n"
		    ,(int)p->ilen);
		
	    mg_write(conn, sbuf, strlen(sbuf));
	    mg_write(conn, p->ibuf, p->ilen);

	} else {
		debugr2((errout,"Bad request or png - returning 404\n"));
		sprintf(sbuf,
			"HTTP/1.0 404 Not Found\r\n"
		    "\r\n"
			"<html><body><h1>Error 404 - Not Found</h1></body></html>");
		
	    mg_write(conn, sbuf, strlen(sbuf));
	}

	return "yes";
}

/* Change the patch display parameters. */
/* Optional - may be NULL */
static int ccwin_set_patch_win(
dispwin *p, 
double hoff, double voff,		/* Offset from c. in fraction of screen, -1.0 .. 1.0 */
double area,					/* Patch area 0..1 */
dw_bg_type bge					/* Background */  
) {
	chws *ws = (chws *)p->pcntx;
	double width, height;

	/* Set background color handling */
	p->fullscreen = 1;
	p->bge = bge;

	if (area < 0.0)
		area = 0.0;
	else if (area > 1.0)
		area = 1.0;

	/* Can't do constant luma/power with larger than 50% area */
	if (bge == dw_bg_cvideo
	 || bge == dw_bg_clight) {
		if (area > 0.5)
			area = 0.5;
	}

	p->area = area;

	if (area < (IHEIGHT * IHEIGHT)/(IWIDTH * IHEIGHT)) {	// Not height limited
		width = height = sqrt(area * IWIDTH * IHEIGHT)/IWIDTH;

	} else {
		height = IHEIGHT/IWIDTH;
		width = area;
	}

	/* Setup window size and position */
	/* The default size is 10% of the width */
	ws->w = floor(width * IWIDTH + 0.5); 
	if (ws->w > IWIDTH)
		ws->w = IWIDTH;
	ws->h = floor(height * IWIDTH + 0.5); 
	if (ws->h > IHEIGHT)
		ws->h = IHEIGHT;

	ws->x = floor((hoff * 0.5 + 0.5) * (IWIDTH - ws->w) + 0.5);
	ws->y = floor((voff * 0.5 + 0.5) * (IHEIGHT - ws->h) + 0.5);

	// Make offset be on an even pixel boundary, so that we know
	// the up-filter phase.
	if (((int)ws->x) & 1)
		ws->x++;
	if (((int)ws->y) & 1)
		ws->y++;

	return 0;
}

chws *new_chws(
ccast_id *cc_id,				/* ChromeCast to open */
double width, double height,	/* Width and height as % */
double hoff, double voff,		/* Offset from center in fraction of screen, range -1.0 .. 1.0 */
int verb, int ddebug) {
	chws *p;
	const char *options[3];
	char port[50];
	int portno = 0;		/* Port number allocated */
	int forcedef = 0;	/* Force default reciever app. */

	if ((p = (chws *)calloc(sizeof(chws), 1)) == NULL) {
		error("new_chws: calloc failed");
		return NULL;
	}

	p->verb = verb;
	p->ddebug = ddebug;

	p->set = chws_set;
	p->update = chws_update;
	p->del = chws_del;

	/* We make sure we round the test patch size and */
	/* location to integer resolution so that we can know */
	/* it's exact relationship to the upsampled pixel locations. */

	/* Setup window size and position */
	/* The default size is 10% of the width */
	p->w = floor(width/100.0 * 0.1 * IWIDTH + 0.5); 
	if (p->w > IWIDTH)
		p->w = IWIDTH;
	p->h = floor(height/100.0 * 0.1 * IWIDTH + 0.5); 
	if (p->h > IHEIGHT)
		p->h = IHEIGHT;

	p->x = floor((hoff * 0.5 + 0.5) * (IWIDTH - p->w) + 0.5);
	p->y = floor((voff * 0.5 + 0.5) * (IHEIGHT - p->h) + 0.5);

	// Make offset be on an even pixel boundary, so that we know
	// the up-filter phase.
	if (((int)p->x) & 1)
		p->x++;
	if (((int)p->y) & 1)
		p->y++;

	if (verb) printf("Opening ChromeCast '%s'\n",cc_id->name);

#ifdef SEND_TEST_FILE
	forcedef = 1;
#endif

	/* Connect to the chrome cast */
	if ((p->cc = new_ccast(cc_id, forcedef)) == NULL) {
		debugr2((errout,"new_chws: new_ccast() failed\n"));
		return NULL;
	}

	p->direct = p->cc->get_direct_send(p->cc);

	if (!p->direct) {
		/* Create a web server */
		options[0] = "listening_ports";
//		sprintf(port,"%d", 0);			/* Use any available */
		sprintf(port,"%d", 8081);		/* Use fixed port for Linux firewall rule */
		options[1] = port;
		options[2] = NULL;

		p->mg = mg_start(&ccwin_ehandler, (void *)p, options);
	
		if ((p->ws_url = mg_get_url(p->mg)) == NULL) {
			debugr2((errout, "mg_get_url() failed\n"));
			chws_del(p);
			return NULL;
		}
		if (p->ddebug)
			printf("Created .png server at '%s'\n",p->ws_url);
	}

	return p;
}


/* ================================================================== */

/* Get RAMDAC values. ->del() when finished. */
/* Return NULL if not possible */
static ramdac *ccwin_get_ramdac(dispwin *p) {
	debugr2((errout,"webdisp doesn't have a RAMDAC\n")); 
	return NULL;
}

/* Set the RAMDAC values. */
/* Return nz if not possible */
static int ccwin_set_ramdac(dispwin *p, ramdac *r, int persist) {
	debugr2((errout,"webdisp doesn't have a RAMDAC\n")); 
	return 1;
}

/* ----------------------------------------------- */
/* Install a display profile and make */
/* it the default for this display. */
/* Return nz if failed */
int ccwin_install_profile(dispwin *p, char *fname, ramdac *r, p_scope scope) {
	debugr2((errout,"webdisp doesn't support installing profiles\n")); 
	return 1;
}

/* Un-Install a display profile */
/* Return nz if failed, */
int ccwin_uninstall_profile(dispwin *p, char *fname, p_scope scope) {
	debugr2((errout,"webdisp doesn't support uninstalling profiles\n")); 
	return 1;
}

/* Get the currently installed display profile. */
/* Return NULL if failed. */
icmFile *ccwin_get_profile(dispwin *p, char *name, int mxlen) {
	debugr2((errout,"webdisp doesn't support getting the current profile\n")); 
	return NULL;
}

/* ----------------------------------------------- */

/* Change the window color. */
/* Return 1 on error, 2 on window being closed */
/* inst_license, inst_licensenc, inst_tamper or inst_syscompat on licening problem */
static int ccwin_set_color(
dispwin *p,
double r, double g, double b	/* Color values 0.0 - 1.0 */
) {
	chws *ws = (chws *)p->pcntx;
	int j;
	double orgb[3];		/* Previous RGB value */
	double kr, kf;
	int update_delay = 0;

	debugr2((errout, "ccwin_set_color called with %f %f %f\n",r,g,b));

	if (p->nowin) {
		debugr2((errout,"ccwin_set_color: nowin - give up\n"));
		return 1;
	}

	orgb[0] = p->rgb[0]; p->rgb[0] = r;
	orgb[1] = p->rgb[1]; p->rgb[1] = g;
	orgb[2] = p->rgb[2]; p->rgb[2] = b;

	for (j = 0; j < 3; j++) {
		if (p->rgb[j] < 0.0)
			p->rgb[j] = 0.0;
		else if (p->rgb[j] > 1.0)
			p->rgb[j] = 1.0;
		p->r_rgb[j] = p->s_rgb[j] = p->rgb[j];
		if (p->out_tvenc) {
			p->r_rgb[j] = p->s_rgb[j] = ((235.0 - 16.0) * p->s_rgb[j] + 16.0)/255.0;

			/* For video encoding the extra bits of precision are created by bit shifting */
			/* rather than scaling, so we need to scale the fp value to account for this. */
			if (p->edepth > 8)
				p->r_rgb[j] = (p->s_rgb[j] * 255 * (1 << (p->edepth - 8)))
				            /((1 << p->edepth) - 1.0); 	
		}
	}

	/* This is probably not actually thread safe... */
	p->ncix++;

#if DDITHER != 1
# pragma message("############################# ccwin.c DDITHER != 1 ##")
#endif

	/* Turn the color into a png file */
	{
		/* We want a raster of IWIDTH x IHEIGHT pixels for web server, */
		/* or p->w x p->h for PNG direct. */
		render2d *rr;
		prim2d *rct;
		depth2d depth = bpc8_2d;
#if DDITHER == 1
		int dither = 0x8002;		/* 0x8002 = error diffuse FG only */
#elif DDITHER == 2
		int dither = 0x4000;		/* 0x4000 = no dither but don't average pixels */
									/* so as to allow pattern to come through. */
#else
		int dither = 0;				/* Don't dither in renderer */
#endif
		double hres = 1.0;					/* Resoltion in pix/mm */
		double vres = 1.0;					/* Resoltion in pix/mm */
		double iw, ih;						/* Size of page in mm (pixels) */
		color2d c;
		unsigned char *ibuf;		/* Memory image of .png file */
		size_t ilen;
		int rv;
#ifdef DO_TIMING
		int stime;
#endif

		if (ws->direct) {
			iw = ws->w;		/* Requested size */
			ih = ws->h;
		} else {
			iw = IWIDTH;
			ih = IHEIGHT;	/* Size of page in mm */
		}

		/* Full screen background: */
		if (p->fullscreen) {
			if (p->bge == dw_bg_grey) {
				ws->bg[0] = 0.2;
				ws->bg[1] = 0.2;
				ws->bg[2] = 0.2;
			} else if (p->bge == dw_bg_cvideo) {
				ws->bg[0] = p->area * (1.0 - r)/(1.0 - p->area); 
				ws->bg[1] = p->area * (1.0 - g)/(1.0 - p->area); 
				ws->bg[2] = p->area * (1.0 - b)/(1.0 - p->area); 
 
			} else if (p->bge == dw_bg_clight) {
				double gamma = 2.3;
				ws->bg[0] = pow(p->area * (1.0 - pow(r, gamma))/(1.0 - p->area), 1.0/gamma);
				ws->bg[1] = pow(p->area * (1.0 - pow(g, gamma))/(1.0 - p->area), 1.0/gamma);
				ws->bg[2] = pow(p->area * (1.0 - pow(b, gamma))/(1.0 - p->area), 1.0/gamma); 

			} else {		/* Assume dw_bg_black */
				ws->bg[0] = 0.0;
				ws->bg[1] = 0.0;
				ws->bg[2] = 0.0;
			}

		/* Use default dark gray background */ 
		} else {
			ws->bg[0] = 0.2;
			ws->bg[1] = 0.2;
			ws->bg[2] = 0.2;
		}

		debugr2((errout, "ccwin_set_color iw %f ih %f\n",iw,ih));

		if ((rr = new_render2d(iw, ih, NULL, hres, vres, rgb_2d,
		     0, depth, dither,
#if DDITHER == 1
			 ccastQuant, NULL, 3.0/255.0
#else
			 NULL, NULL, 0.0
#endif
			 )) == NULL) {
			a1loge(g_log, 1,"ccwin: new_render2d() failed\n");
			return 1;
		}
	
		/* Set the background color */
		c[0] = ws->bg[0];
		c[1] = ws->bg[1];
		c[2] = ws->bg[2];
		rr->set_defc(rr, c);
	
		c[0] = p->r_rgb[0];
		c[1] = p->r_rgb[1];
		c[2] = p->r_rgb[2];
		if (ws->direct)
			rr->add(rr, rct = new_rect2d(rr, 0.0, 0.0, ws->w, ws->h, c));
		else
			rr->add(rr, rct = new_rect2d(rr, ws->x, ws->y, ws->w, ws->h, c));

#if DDITHER == 2			/* Use dither pattern */
		{
			double rgb[3];
			double dpat[CCDITHSIZE][CCDITHSIZE][3];
			double (*cpat)[MXPATSIZE][MXPATSIZE][TOTC2D];
			int i, j;

			/* Get a chrome cast dither pattern to match target color */
			for (i = 0; i < 3; i++)
				rgb[i] = p->r_rgb[i] * 255.0;
			get_ccast_dith(dpat, rgb);

			if ((cpat = malloc(sizeof(double) * MXPATSIZE * MXPATSIZE * TOTC2D)) == NULL) {
				a1loge(g_log, 1, "ccwin: malloc of dither pattern failed\n");
				return 1;
			}
			
			for (i = 0; i < CCDITHSIZE; i++) {
				for (j = 0; j < CCDITHSIZE; j++) {
					int k = (((int)IHEIGHT-2) - j) % CCDITHSIZE;	/* Flip to origin bot left */
					(*cpat)[i][k][0] = dpat[i][j][0]/255.0;			/* (HEIGHT-2 is correct!) */
					(*cpat)[i][k][1] = dpat[i][j][1]/255.0;
					(*cpat)[i][k][2] = dpat[i][j][2]/255.0;
				}
			}
			
			set_rect2d_dpat((rect2d *)rct, cpat, CCDITHSIZE, CCDITHSIZE);
		}
#endif /* DDITHER == 2 */

#ifdef CCTEST_PATTERN
#pragma message("############################# ccwin.c TEST_PATTERN defined ! ##")
		if (getenv("ARGYLL_CCAST_TEST_PATTERN") != NULL) {
			verbose(0, "Writing test pattern to '%s'\n","testpattern.png");
			if (r->write(r, "testpattern.png", 1, NULL, NULL, png_file)) {
				a1loge(g_log, 1, "ccwin: render->write failed\n");
				return 1;
			}
		}
#else	/* !CCTEST_PATTERN */
# ifdef WRITE_PNG		/* Write it to a file so that we can look at it */
#  pragma message("############################# spectro/ccwin.c WRITE_PNG is enabled ######")
		if (r->write(rr, "ccwin.png", 1, NULL, NULL, png_file)) {
			a1loge(g_log, 1, "ccwin: render->write failed\n");
			return 1;
		}
# endif	/* WRITE_PNG */
#endif	/* !CCTEST_PATTERN */


#ifdef DO_TIMING
		stime = msec_time();
#endif

		rv = rr->write(rr, "MemoryBuf", 1, &ibuf, &ilen, png_mem);


		if (rv) {
			a1loge(g_log, 1, "ccwin: render->write failed\n");
			return 1;
		}
		rr->del(rr);
#ifdef DO_TIMING
		stime = msec_time() - stime;
		printf("render->write took %d msec\n",stime);
#endif
		if (ws->update(ws, ibuf, ilen, ws->bg)) {
			a1loge(g_log, 1, "ccwin: color update failed\n");
			return 1;
		}
		p->ccix = p->ncix;
	}


	/* If update is notified asyncronously ... */
	while(p->ncix != p->ccix) {
		msec_sleep(50);
	}
//printf("#################################################################\n");
//printf("#################     RGB update notified        ################\n");
//printf("#################################################################\n");

	/* Allow for display update & instrument delays */
	update_delay = dispwin_compute_delay(p, orgb);
	debugr2((errout, "ccwin_set_color delaying %d msec\n",update_delay));
	msec_sleep(update_delay);

	return 0;
}

/* Set/unset the full screen background color flag */
/* Return nz on error */
static int ccwin_set_fc(dispwin *p, int fullscreen) {
	p->fullscreen = fullscreen;

	return 0;
}


/* ----------------------------------------------- */
/* Set the shell set color callout */
void ccwin_set_callout(
dispwin *p,
char *callout
) {
	debugr2((errout,"ccwin_set_callout called with '%s'\n",callout));

	p->callout = strdup(callout);
}

/* ----------------------------------------------- */
/* Destroy ourselves */
static void ccwin_del(
dispwin *p
) {
	chws *ws;

	debugr2((errout,"ccwin_del called with %p\n",p));

	if (p == NULL)
		return;

	ws = (chws *)p->pcntx;

	if (ws != NULL)
		ws->del(ws);

	if (p->name != NULL)
		free(p->name);
	if (p->description != NULL)
		free(p->description);
	if (p->callout != NULL)
		free(p->callout);

	free(p);
}

/* ----------------------------------------------- */

/* Create a ChromeCast display test window, default grey */
dispwin *new_ccwin(
ccast_id *cc_id,				/* ChromeCast to open */
double width, double height,	/* Width and height in mm. (TV width assumed to b 1000mm) */
double hoff, double voff,		/* Offset from center in fraction of screen, range -1.0 .. 1.0 */
int nowin,						/* NZ if no window should be created - RAMDAC access only */
int native,						/* X0 = use current per channel calibration curve */
								/* X1 = set native linear output and use ramdac high precn. */
								/* 0X = use current color management cLut (MadVR) */
								/* 1X = disable color management cLUT (MadVR) */
int *noramdac,					/* Return nz if no ramdac access. native is set to X0 */
int *nocm,						/* Return nz if no CM cLUT access. native is set to 0X */
int out_tvenc,					/* 1 = use RGB Video Level encoding */
int fullscreen,					/* NZ if whole screen should be filled with black */
int noinitpatch,				/* NZ if no initial test patch should be shown */
int verb,						/* NZ for verbose prompts */
int ddebug						/* >0 to print debug statements to stderr */
) {
	dispwin *p = NULL;
	char *cp;
	chws *ws = NULL;
	const char *options[3];

	debug2((errout,"new_ccwin called\n"));

	if ((p = (dispwin *)calloc(sizeof(dispwin), 1)) == NULL) {
		debugr2((errout,"new_ccwin failed because malloc failed\n"));
		return NULL;
	}

	/* !!!! Make changes in dispwin.c & madvrwin.c as well !!!! */
	p->name = strdup("Web Window");
	p->width = width;
	p->height = height;
	p->nowin = nowin;
	p->native = native;
	p->out_tvenc = out_tvenc;
	p->fullscreen = fullscreen;
	p->ddebug = ddebug;
	p->get_ramdac          = ccwin_get_ramdac;
	p->set_ramdac          = ccwin_set_ramdac;
	p->install_profile     = ccwin_install_profile;
	p->uninstall_profile   = ccwin_uninstall_profile;
	p->get_profile         = ccwin_get_profile;
	p->set_color           = ccwin_set_color;
	p->set_fc              = ccwin_set_fc;
	p->set_patch_win       = ccwin_set_patch_win;
	p->set_update_delay    = dispwin_set_update_delay;
	p->set_settling_delay  = dispwin_set_settling_delay;
	p->enable_update_delay = dispwin_enable_update_delay;
	p->set_callout         = ccwin_set_callout;
	p->del                 = ccwin_del;

	if (noramdac != NULL)
		*noramdac = 1;
	p->native &= ~1;

	if (nocm != NULL)
		*nocm = 1;
	p->native &= ~2;

	p->rgb[0] = p->rgb[1] = p->rgb[2] = 0.5;	/* Set Grey as the initial test color */
	
	dispwin_set_default_delays(p);

	p->ncix = 1;

	p->fdepth = 8;				/* Assume this by API */
	p->rdepth = p->fdepth;		/* Assumed */
	p->ndepth = p->rdepth;		/* Assumed */
	p->nent = 0;				/* No ramdac */
	p->edepth = 8;				/* Assumed */

	/* Basic object is initialised, so create connection to ChromeCast */
	if ((ws = new_chws(cc_id, width, height, hoff, voff, verb, ddebug)) == NULL) {
		debugr2((errout,"new_ccwin failed - new_chws() failed\n"));
		p->del(p);
		return NULL;
	}

	/* Extra delay ccast adds after confirming load */
	p->extra_update_delay = ws->cc->get_load_delay(ws->cc) / 1000.0;

	p->pcntx = (void *)ws;

	/* Create a suitable description */
	{
		char buf[100];
		sprintf(buf,"ChromeCast '%s'",cc_id->name);
		p->description = strdup(buf);
	}

    // Set a default first color
	if (!noinitpatch && ccwin_set_color(p, 128.0, 128.0, 128.0)) {
		debugr2((errout,"new_ccwin failed because set_color() failed\n"));
		p->del(p);
		return NULL;
	}

	debugr2((errout,"new_ccwin: return sucessfully\n"));

	return p;
}

