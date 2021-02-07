
#ifndef CCWIN_H

/* 
 * Argyll Color Correction System
 * ChromeCast Display target patch window
 *
 * Author: Graeme W. Gill
 * Date:   8/9/14
 *
 * Copyright 2014 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/* See ccast/ccmdns.h for function to get list of available ChromeCasts */

/* Create a ChromeCast display test window, default grey */
dispwin *new_ccwin(
ccast_id *cc_id,                /* ChromeCast to open */
double width, double height,	/* Width and height as multiplier of 10% width default. */
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
);

#define CCWIN_H
#endif /* CCWIN_H */
