
#ifndef DUMMYWIN_H

/* 
 * Argyll Color Correction System
 * Virtual target patch window
 *
 * Author: Graeme W. Gill
 * Date:   26/6/13
 *
 * Copyright 2013-2018 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/* Create a dummy (i.e. not real) isplay test window */
dispwin *new_dummywin(
double width, double height,	/* Width and height in mm - turned into % of screen */
double hoff, double voff,		/* Offset from center in fraction of screen - ignored */
int nowin,						/* NZ if no window should be created - RAMDAC access only */
int native,						/* X0 = use current per channel calibration curve */
								/* X1 = set native linear output and use ramdac high precn. */
								/* 0X = use current color management cLut */
								/* 1X = disable color management cLUT */
int *noramdac,					/* Return nz if no ramdac access. native is set to X0 */
int *nocm,						/* Return nz if no CM cLUT access. native is set to 0X */
int out_tvenc,					/* 1 = use RGB Video Level encoding */
int fullscreen,					/* NZ if whole screen should be filled with black */
int verb,						/* NZ for verbose prompts */
int ddebug						/* >0 to print debug statements to stderr */
);

#define DUMMYWIN_H
#endif /* DUMMYWIN_H */
