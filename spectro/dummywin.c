

/* 
 * Argyll Color Management System
 * Dummy (non existant) target patch window.
 * Allows for shell callout function.
 *
 * Author: Graeme W. Gill
 * Date:   3/4/12
 *
 * Copyright 2013-2018 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <string.h>
#ifdef NT
# include <winsock2.h>
#endif
#include "copyright.h"
#include "aconfig.h"
#include "icc.h"
#include "numsup.h"
#include "cgats.h"
#include "conv.h"
#include "dispwin.h"
#include "dummywin.h"
#include "conv.h"
#include "mongoose.h"

#define ENABLE_RAMDAC

#undef DEBUG
//#define STANDALONE_TEST

#ifdef DEBUG
#define errout stderr
# define debug(xx)	fprintf(errout, xx )
# define debug2(xx)	fprintf xx
# define debugr(xx)	fprintf(errout, xx )
# define debugr2(xx)	fprintf xx
# define debugrr(xx)	fprintf(errout, xx )
# define debugrr2(xx)	fprintf xx
# define debugrr2l(lev, xx)	fprintf xx
#else
#define errout stderr
# define debug(xx) 
# define debug2(xx)
# define debugr(xx) if (p->ddebug) fprintf(errout, xx ) 
# define debugr2(xx) if (p->ddebug) fprintf xx
# define debugrr(xx) if (callback_ddebug) fprintf(errout, xx ) 
# define debugrr2(xx) if (callback_ddebug) fprintf xx
# define debugrr2l(lev, xx) if (callback_ddebug >= lev) fprintf xx
#endif

/* ----------------------------------------------- */

/* Get RAMDAC values. ->del() when finished. */
/* Return NULL if not possible */
static ramdac *dummywin_get_ramdac(dispwin *p) {
	debugr("dummydisp doesn't have a RAMDAC\n"); 
	return NULL;
}

/* Set the RAMDAC values. */
/* Return nz if not possible */
static int dummywin_set_ramdac(dispwin *p, ramdac *r, int persist) {
	debugr("dummydisp doesn't have a RAMDAC\n"); 
	return 1;
}

/* ----------------------------------------------- */
/* Install a display profile and make */
/* it the default for this display. */
/* Return nz if failed */
int dummywin_install_profile(dispwin *p, char *fname, ramdac *r, p_scope scope) {
	debugr("dummydisp doesn't support installing profiles\n"); 
	return 1;
}

/* Un-Install a display profile */
/* Return nz if failed, */
int dummywin_uninstall_profile(dispwin *p, char *fname, p_scope scope) {
	debugr("dummydisp doesn't support uninstalling profiles\n"); 
	return 1;
}

/* Get the currently installed display profile. */
/* Return NULL if failed. */
icmFile *dummywin_get_profile(dispwin *p, char *name, int mxlen) {
	debugr("dummydisp doesn't support getting the current profile\n"); 
	return NULL;
}

/* ----------------------------------------------- */

/* Change the window color. */
/* Return 1 on error, 2 on window being closed */
static int dummywin_set_color(
dispwin *p,
double r, double g, double b	/* Color values 0.0 - 1.0 */
) {
	int j;
	double orgb[3];		/* Previous RGB value */
	double kr, kf;
	int update_delay = 0;

	debugr("dummywin_set_color called\n");

	if (p->nowin)
		return 1;

	orgb[0] = p->rgb[0]; p->rgb[0] = r;
	orgb[1] = p->rgb[1]; p->rgb[1] = g;
	orgb[2] = p->rgb[2]; p->rgb[2] = b;

	if (p->callout != NULL) {
		int rv;
		char *cmd;

		if ((cmd = malloc(strlen(p->callout) + 200)) == NULL)
			error("Malloc of command string failed");

		sprintf(cmd, "%s %d %d %d %f %f %f",p->callout,
			        (int)(r * 255.0 + 0.5),(int)(g * 255.0 + 0.5),(int)(b * 255.0 + 0.5), r, g, b);
		if ((rv = system(cmd)) != 0)
			warning("System command '%s' failed with %d",cmd,rv); 
		free(cmd);
	}

	/* Allow for display update & instrument delays */
	update_delay = dispwin_compute_delay(p, orgb);
	debugr2((errout, "dummywin_set_color delaying %d msec\n",update_delay));
	msec_sleep(update_delay);

	return 0;
}

/* Set/unset the full screen black flag */
/* Return nz on error */
static int dummywin_set_fc(dispwin *p, int fullscreen) {
	int perc, bgperc, bgmode, border;

	p->fullscreen = fullscreen;

	return 0;
}

/* ----------------------------------------------- */
/* set patch info */
/* Return nz on error */
static int dummywin_set_pinfo(dispwin *p, int pno, int tno) {

	return 0;
}

/* ----------------------------------------------- */
/* Set the shell set color callout */
void dummywin_set_callout(
dispwin *p,
char *callout
) {
	debugr2((errout,"dummywin_set_callout called with '%s'\n",callout));

	p->callout = strdup(callout);
}

/* ----------------------------------------------- */
/* Destroy ourselves */
static void dummywin_del(
dispwin *p
) {

	debugr("dummywin_del called\n");

	if (p == NULL)
		return;

	if (p->name != NULL)
		free(p->name);
	if (p->description != NULL)
		free(p->description);
	if (p->callout != NULL)
		free(p->callout);

	/* Since we don't restore the display, delete these here */
	if (p->oor != NULL) {
		p->oor->del(p->oor);
		p->oor = NULL;
	}
	if (p->or != NULL) {
		p->or->del(p->or);
		p->or = NULL;
	}
	if (p->r != NULL) {
		p->r->del(p->r);
		p->r = NULL;
	}

	free(p);
}

/* ----------------------------------------------- */

/* Create a dummy display test window, default grey */
dispwin *new_dummywin(
double width, double height,	/* Width and height in mm - turned into % of screen */
double hoff, double voff,		/* Offset from center in fraction of screen - ignored */
int nowin,						/* NZ if no window should be created - RAMDAC access only */
int native,						/* X0 = use current per channel calibration curve */
								/* X1 = set native linear output and use ramdac high precn. */
								/* 0X = use current color management cLut (MadVR) */
								/* 1X = disable color management cLUT (MadVR) */
int *noramdac,					/* Return nz if no ramdac access. native is set to X0 */
int *nocm,						/* Return nz if no CM cLUT access. native is set to 0X */
int out_tvenc,					/* 1 = use RGB Video Level encoding */
int fullscreen,					/* NZ if whole screen should be filled with black */
int verb,						/* NZ for verbose prompts */
int ddebug						/* >0 to print debug statements to stderr */
) {
	dispwin *p = NULL;
	char *cp;
	const char *options[3];
	char port[50];

	debug("new_dummywin called with native = %d\n");

	if (out_tvenc) {
		if (ddebug) fprintf(stderr,"new_dummywin failed because out_tvenc set\n");
		return NULL;
	}

	if ((p = (dispwin *)calloc(sizeof(dispwin), 1)) == NULL) {
		if (ddebug) fprintf(stderr,"new_dummywin failed because malloc failed\n");
		return NULL;
	}

	/* !!!! Make changes in dispwin.c & webwin.c etc. as well !!!! */
	p->name = strdup("Web Window");
	p->width = width;
	p->height = height;
	p->nowin = nowin;
	p->native = native;
	p->out_tvenc = 0;
	p->fullscreen = fullscreen;
	p->ddebug = ddebug;
	p->get_ramdac          = dummywin_get_ramdac;
	p->set_ramdac          = dummywin_set_ramdac;
	p->install_profile     = dummywin_install_profile;
	p->uninstall_profile   = dummywin_uninstall_profile;
	p->get_profile         = dummywin_get_profile;
	p->set_color           = dummywin_set_color;
	p->set_fc              = dummywin_set_fc;
	p->set_pinfo           = dummywin_set_pinfo;
	p->set_update_delay    = dispwin_set_update_delay;
	p->set_settling_delay  = dispwin_set_settling_delay;
	p->enable_update_delay = dispwin_enable_update_delay;
	p->set_callout         = dummywin_set_callout;
	p->del                 = dummywin_del;

	debugr2((errout, "new_dummywin got native = %d\n",native));

#ifndef ENABLE_RAMDAC
	if (noramdac != NULL)
		*noramdac = 1;
	p->native &= ~1;
#endif

	p->rgb[0] = p->rgb[1] = p->rgb[2] = 0.5;	/* Set Grey as the initial test color */

	dispwin_set_default_delays(p);

	p->fdepth = 8;				/* Assume this */
	p->rdepth = p->fdepth;		/* Assumed */
	p->ndepth = p->rdepth;		/* Assumed */
#ifdef ENABLE_RAMDAC
	p->nent = (1 << p->ndepth);	
#else
	p->nent = 0;				/* No ramdac */
#endif
	p->edepth = 16;				/* Assumed */

	p->set_fc(p, fullscreen);

	/* Create a suitable description */
	{
		char buf[1000];

		sprintf(buf,"ArgyllCMS Patches");
		p->description = strdup(buf);

		if (verb)
			printf("Created dummy window\n");
	}

#ifdef ENABLE_RAMDAC

	/* Save the original ramdac, which gets restored on exit */
	if ((p->or = p->get_ramdac(p)) != NULL) {

		debugr("Saved original VideoLUT\n");

		if (noramdac != NULL)
			*noramdac = 0;

		/* Copy original ramdac that never gets altered */
		if ((p->oor = p->or->clone(p->or)) == NULL) {
			dummywin_del(p);
			debugr("ramdac clone failed - memory ?\n");
			return NULL;
		}

		/* Create a working ramdac for native or other use */
		if ((p->r = p->or->clone(p->or)) == NULL) {
			dummywin_del(p);
			debugr("ramdac clone failed - memory ?\n");
			return NULL;
		}

	} else {
		debugr2((errout,"Unable to access VideoLUT\n"));
		if (noramdac != NULL)
			*noramdac = 1;
		p->oor = p->or = p->r = NULL;
	}

	if (!p->nowin) {

		/* Make sure initial test color is displayed */
		dummywin_set_color(p, p->rgb[0], p->rgb[1], p->rgb[2]);
	}
#endif

	debugr("new_dummywin: return sucessfully\n");

	return p;
}

