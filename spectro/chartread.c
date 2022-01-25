
/* Spectrometer/Colorimeter target test chart reader */

/* 
 * Argyll Color Management System
 *
 * Author: Graeme W. Gill
 * Date:   4/10/96
 *
 * Copyright 1996 - 2013, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This program reads a reflective or transmissive print target chart */
/* using a spectrometer or colorimeter. */

/* TTBD
 *
 *  When there is a strip misread, should offer the option
 *	of skipping the strip.
 *
 *
 *	Should an a -X [xl] mode that reads a simple list of readings
 *  from a file.
 *
 *	Someone reported that XY mode (spectroscan) didn't work for one paper
 *  orientation ?
 *
 *  Should fix XY chart read to also allow interruption/save/resume,
 *  just like the strip reading code.
 *
 *	Should add verbose option to print average & max DE to expected value
 *  for each patch/strip read.
 *
 */
 
/*
 * Nomencalture:
 *
 *	Largely due to how the strip readers name things, the following terms
 *  are used for how patches are grouped:
 *
 *  Step:   One test patch in a pass, usually labelled with a number.
 *  Pass:	One row of patches in a strip. A pass is usually labeled
 *          with a unique alphabetic label. 
 *  Strip:  A group of passes that can be read by a strip reader.
 *          For an XY instrument, the strip is a complete sheet, and
 *          a each pass is one column. The rows of an XY chart are
 *          the step numbers within a pass.
 *  Sheet:  One sheet of paper, containing full and partial strips.
 *          For an XY instrument, there will be only one strip per sheet.
 *           
 */

#undef DEBUG

#define COMPORT 1		/* Default com port 1..4 */

#undef TEST_EVENT_CALLBACK      /* Report async event callbacks, and implement beep prompt there. */

#undef USESTRDELTA		/* [Und] Use patch delta's for correlation rather than match DE */
						/* Doesn't seem to work as well. Why ? */

#ifdef __MINGW32__
# define WINVER 0x0500
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#ifndef SALONEINSTLIB
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "xicc.h"
#include "conv.h"
#include "ui.h"
#include "icc.h"
#else /* SALONEINSTLIB */
#include "sa_config.h"
#include "numsup.h"
#include "cgats.h"
#include "rspl1.h"
#include "xspect.h"
#include "xcolorants.h"
#include "xcal.h"
#include "conv.h"
#include "sa_conv.h"
#endif /* SALONEINSTLIB */
#include "cgats.h"
#include "insttypes.h"
#include "icoms.h"
#include "inst.h"
#include "ccmx.h"
#include "ccss.h"
#ifndef SALONEINSTLIB
#include "dispwin.h"
#include "ccast.h"
#include "dispsup.h"
#endif /* !SALONEINSTLIB */
#include "alphix.h"
#include "sort.h"
#include "instappsup.h"

#include <stdarg.h>

#if defined (NT)
#include <conio.h>
#endif

#ifdef NEVER    /* Not currently used */

/* Convert control chars to ^[A-Z] notation in a string */
static char *
fix_asciiz(char *s) {
	static char buf [200];
	char *d;
	for(d = buf; ;) {
		if (*s < ' ' && *s > '\000') {
			*d++ = '^';
			*d++ = *s++ + '@';
		} else
		*d++ = *s++;
		if (s[-1] == '\000')
			break;
	}
	return buf;
}
#endif /* NEVER */

/* Return the normal Delta E, given two 100 scaled XYZ values */
static double xyzLabDE(double ynorm, double *pat, double *ref) {
	int i;
	double Lab1[3];
	double Lab2[3];

	for (i = 0; i < 3; i++) {
		Lab1[i] = ynorm * pat[i]/100.0;
		Lab2[i] = ref[i]/100.0;
	}

	icmXYZ2Lab(&icmD50, Lab1, Lab1);
	icmXYZ2Lab(&icmD50, Lab2, Lab2);

	return icmLabDE(Lab1, Lab2);
}

#ifdef USESTRDELTA
/* Return the -ve correlation of the delta E's between steps */
static double xyzLabcorr(double *pat0, double *ref0,
                         double *pat1, double *ref1) {
	double p0[3], p1[3], pd[3];
	double r0[3], r1[3], rd[3];

	icmXYZ2Lab(&icmD50, p0, pat0);
	icmXYZ2Lab(&icmD50, p1, pat1);
	ICMSUB3(pd, p0, p1);
	icmXYZ2Lab(&icmD50, r0, ref0);
	icmXYZ2Lab(&icmD50, r1, ref1);
	ICMSUB3(rd, r0, r1);

	return -icmDot3(rd, pd);
}
#endif

/* A chart read color structure */
/* This can hold all representations simultaniously */
typedef struct {
	char *id;					/* Id string (e.g. "1") */
	char *loc;					/* Location string (e.g. "A1") */
	int loci;					/* Location integer = pass * 256 + step */

	int n;						/* Number of colorants */
	double dev[ICX_MXINKS];		/* Value of colorants */
	double eXYZ[3];				/* Expected XYZ values (100.0 scale for ref.) */

	int rr;						/* nz if reading read (used for tracking unread patches) */

	inst_meas_type mtype;		/* Measurement type */
	inst_meas_cond mcond;		/* Measurement conditions */
	double XYZ[3];				/* Colorimeter readings (100.0 scale for ref.) */

	xspect sp;					/* Spectrum. sp.spec_n > 0 if valid, 100 scaled for ref. */
} chcol;

/* Convert a base 62 character into a number */
/* (This is used for converting the PASSES_IN_STRIPS string */
/* (Could convert this to using an alphix("0-9A-Za-Z")) */
static int b62_int(char *p) {
	int cv, rv;

	cv = *p;
	if (cv == '\000')
		rv = 0;
	else if (cv <= '9')
		rv = cv - '0';
	else if (cv <= 'Z')
		rv = cv - 'A' + 10;
	else
		rv = cv - 'a' + 36;
	return rv;
}

#ifdef TEST_EVENT_CALLBACK
void test_event_callback(void *cntx, inst_event_type event) {
	a1logd(g_log,0,"Got event_callback with 0x%x\n",event);

	if (event == inst_event_scan_ready)
		normal_beep();
}
#endif

/* Deal with an instrument error. */
/* Return 0 to retry, 1 to abort */
static int ierror(inst *it, inst_code ic) {
	int ch;
	empty_con_chars();
	printf("Got '%s' (%s) error.\nHit Esc or 'q' to give up, any other key to retry:",
	       it->inst_interp_error(it, ic), it->interp_error(it, ic));
	fflush(stdout);
	ch = next_con_char();
	printf("\n");
	if (ch == 0x03 || ch == 0x1b || ch == 'q' || ch == 'Q')	/* ^C, Escape or Q */
		return 1;
	return 0;
}

/* Read all the strips, and return nonzero on abort/error */
static int
read_strips(
instType itype,		/* Chart instrument type */
chcol **scols,		/* Location sorted pointers to cols (return values) */
instType *atype,	/* Return the instrument type used to read the chart */
int npat,			/* Total valid patches */
int totpa, 			/* Total passes (rows) */
int stipa,			/* Steps (patches) in pass (Excluding DTP51 Max/Min) */
int *pis,			/* Passes in each strip (rows in each sheet), 0 terminated */
alphix *paix,		/* Pass (row) index generators */
alphix *saix,		/* Step (patch) index generators */
int ixord,			/* Index order, 0 = pass then step */
int rstart,			/* Random start/chart id */
int rand,			/* Random order used  - can do auto strip ID & Bi-Dir */
int hex,			/* Hexagon test patches */
icompath *ipath,	/* Instrument path to open */
flow_control fc,	/* flow control */
double plen,		/* Patch length in mm (used by DTP20/41) */
double glen,		/* Gap length in mm  (used by DTP20/41) */
double tlen,		/* Trailer length in mm  (used by DTP41T) */
int trans,			/* Use transmission mode */
int emis,			/* Use emissive mode */
int displ,			/* 1 = Use display emissive mode, 2 = display bright rel. */
					/* 3 = display white rel. */
int ditype,			/* Display type selection charater */
inst_opt_filter fe,	/* Optional filter */
xcalstd scalstd,	/* X-Rite calibration standard to set */
xcalstd *ucalstd,	/* X-Rite calibration standard actually used */
int nocal,			/* Disable initial calibration */
int disbidi,		/* 1 = Disable automatic bi-directional strip recognition, 2 = force enable */
int highres,		/* Use high res spectral mode */
char *ccxxname,		/* Colorimeter Correction/Colorimeter Calibration name */
icxObserverType obType,	/* ccss observer */
xspect custObserver[3],	/* If obType = icxOT_custom */
double scan_tol,	/* Modify patch consistency tolerance */
int pbypatch,		/* Patch by patch measurement */
int xtern,			/* Use external (user supplied) values rather than instument read */
					/* 1 = Lab, 2 = XYZ */ 
int spectral,		/* Generate spectral info flag */
int uvmode,			/* ~~~ i1pro2 test mode ~~~ */
int accurate_expd,	/* Expected values can be assumed to be accurate */
int emit_warnings,	/* Emit warnings for wrong strip, unexpected value */
int doplot,			/* Plot each spectra in patch by patch mode */
a1log *log			/* verb, debug & error log */
) {
	inst *it = NULL;
	inst_mode  cap;		/* Mode capability */
	inst2_capability cap2;
	inst3_capability cap3;
	int n, i, j;
	int rmode = 0;		/* Read mode, 0 = spot, 1 = strip, 2 = xy, 3 = chart */
	int svdmode = 0;	/* Saved mode, 0 = no, 1 = use saved mode */
	inst_code rv;
	baud_rate br = baud_38400;	/* Target baud rate */
	int skipp = 0;		/* Initial strip readings to skip */
	int	nextrap = 0;	/* Number of extra patches for max and min */
	int ch;

	if (xtern == 0) {		/* Use instrument values */

		/* Instrument that the chart is set up for */
		if (itype == instDTP51) {
			skipp = 1;		/* First reading is the Max density patch */
			nextrap = 2;
		}

		if ((it = new_inst(ipath, 0, log, DUIH_FUNC_AND_CONTEXT)) == NULL) {
			printf("Unknown, inappropriate or no instrument detected\n");
			return -1;
		}

#ifdef TEST_EVENT_CALLBACK
		it->set_event_callback(it, test_event_callback, (void *)it);
#endif

		/* Establish communications */
		if ((rv = it->init_coms(it, br, fc, 15.0)) != inst_ok) {
			printf("Establishing communications with instrument failed with message '%s' (%s)\n",
			       it->inst_interp_error(it, rv), it->interp_error(it, rv));
			it->del(it);
			return -1;
		}

		/* Set it up the way we want */
		if ((rv = it->init_inst(it)) != inst_ok) {
			printf("Initialising instrument failed with message '%s' (%s)\n",
			       it->inst_interp_error(it, rv), it->interp_error(it, rv));
			it->del(it);
			return -1;
		}

		/* For reflective */
		if (emis == 0 && trans == 0) {

			/* set XRGA conversion */
			if (scalstd != xcalstd_none) {
				if ((rv = it->get_set_opt(it, inst_opt_set_xcalstd, scalstd)) != inst_ok) {
					printf("Warning: Setting calibration standard not supported by instrument\n");
				}
			}
	
			/* Get actual XRGA conversion */
			if (ucalstd != NULL) {
				if ((rv = it->get_set_opt(it, inst_opt_get_xcalstd, ucalstd)) != inst_ok) {
					*ucalstd = xcalstd_none;
				}
			}
		}

		*atype = it->get_itype(it);		/* Actual instrument type */
		if (*atype != itype)
			a1logv(log, 1, "Warning: chart is for %s, using instrument %s\n",inst_name(itype),inst_name(*atype));

		{
			int ccssset = 0;
			inst_mode mode = 0;

			it->capabilities(it, &cap, &cap2, &cap3);

			if (trans) {
				if (!IMODETST(cap, inst_mode_transmission)) {
					printf("Need transmission reading capability,\n");
					printf("and instrument doesn't support it\n");
					it->del(it);
					return -1;
				}

			} else if (emis || displ) {

				if (emis) {
					if (it->check_mode(it, inst_mode_emis_spot) != inst_ok
					 && it->check_mode(it, inst_mode_emis_strip) != inst_ok) {
						printf("Need emissive spot or strip reading capability\n");
						printf("and instrument doesn't support it\n");
						it->del(it);
						return -1;
					}
				} else {
					/* Should we allow for non-adaptive mode ? */
					if (it->check_mode(it, inst_mode_emis_spot) != inst_ok) {
						printf("Need emissive reading capability\n");
						printf("and instrument doesn't support it\n");
						it->del(it);
						return -1;
					}
				}

			} else {	/* reflectance */
				if (!IMODETST(cap, inst_mode_reflection)) {
					printf("Need reflection spot, strip, xy or chart reading capability,\n");
					printf("and instrument doesn't support it\n");
					it->del(it);
					return -1;
				}
			}

			/* Set display type or calibration mode */
			if (ditype != 0) {

				if (cap2 & inst2_disptype) {
					int ix;
					if ((ix = inst_get_disptype_index(it, ditype, 0)) < 0) {
						printf("Setting display type ix %d failed\n",ix);
						it->del(it);
						return -1;
					}
		
					if ((rv = it->set_disptype(it, ix)) != inst_ok) {
						printf("Setting display type ix %d not supported by instrument\n",ix);
						it->del(it);
						return -1;
					}
				} else
					printf("Display type ignored - instrument doesn't support display type selection\n");
			}

			if (spectral && !IMODETST(cap, inst_mode_spectral)) {
				printf("Warning: Instrument isn't capable of spectral measurement\n");
				spectral = 0;
			}

			/* Colorimeter Correction Matrix */
			if (ccxxname[0] != '\000') {
				ccss *cs = NULL;
				ccmx *cx = NULL;

				if ((cx = new_ccmx()) == NULL) {
					printf("\nnew_ccmx failed\n");
					it->del(it);
					return -1;
				}
				if (cx->read_ccmx(cx,ccxxname) == 0) {
					if ((cap2 & inst2_ccmx) == 0) {
						printf("\nInstrument doesn't have Colorimeter Correction Matrix capability\n");
						it->del(it);
						return -1;
					}
					if ((rv = it->col_cor_mat(it, cx->dtech, cx->cc_cbid, cx->matrix)) != inst_ok) {
						printf("\nSetting Colorimeter Correction Matrix failed with error :'%s' (%s)\n",
					     	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
						cx->del(cx);
						it->del(it);
						return -1;
					}
					cx->del(cx);
				} else {
					cx->del(cx);
					cx = NULL;
					
					/* CCMX failed, try CCSS */
					if ((cs = new_ccss()) == NULL) {
						printf("\nnew_ccss failed\n");
						it->del(it);
						return -1;
					}
					if (cs->read_ccss(cs,ccxxname)) {
						printf("\nReading CCMX/CCSS File '%s' failed with error %d:'%s'\n",
					     	       ccxxname, cs->errc, cs->err);
						cs->del(cs);
						it->del(it);
						return -1;
					}
					if ((cap2 & inst2_ccss) == 0) {
						printf("\nInstrument doesn't have Colorimeter Calibration Spectral Sample capability\n");
						cs->del(cs);
						it->del(it);
						return -1;
					}
					if ((rv = it->get_set_opt(it, inst_opt_set_ccss_obs, obType, custObserver)) != inst_ok) {
						printf("\nSetting CCSS observer failed with error :'%s' (%s)\n",
					     	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
						cs->del(cs);
						it->del(it);
						return -1;
					}
					if ((rv = it->col_cal_spec_set(it, cs->dtech, cs->samples, cs->no_samp)) != inst_ok) {
						printf("\nSetting Colorimeter Calibration Spectral Samples failed with error :'%s' (%s)\n",
					     	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
						cs->del(cs);
						it->del(it);
						return -1;
					}
					ccssset = 1;
					cs->del(cs);
				}
			}

			/* If non-standard observer wasn't set by a CCSS file above */
			if (obType != icxOT_default && (cap2 & inst2_ccss) && ccssset == 0) {
				if ((rv = it->get_set_opt(it, inst_opt_set_ccss_obs, obType, custObserver)) != inst_ok) {
					printf("\nSetting CCSS observer failed with error :'%s' (%s)\n",
				     	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
					it->del(it);
					return -1;
				}
			}

			/* set reflective mode filters  */
			if (fe != inst_opt_filter_unknown) {
				if ((rv = it->get_set_opt(it, inst_opt_set_filter, fe)) != inst_ok) {
					printf("Setting requested filter not supported by instrument\n");
					it->del(it);
					return -1;
				}
			}

			/* Disable initial calibration of machine if selected */
			if (nocal != 0){
				if ((rv = it->get_set_opt(it,inst_opt_noinitcalib, 0)) != inst_ok) {
					printf("Setting no-initial calibrate failed with '%s' (%s)\n",
				       it->inst_interp_error(it, rv), it->interp_error(it, rv));
					printf("Disable inital-calibrate not supported\n");
				}
			}

			/* If it battery powered, show the status of the battery */
			if ((cap2 & inst2_has_battery)) {
				double batstat = 0.0;
				if ((rv = it->get_set_opt(it, inst_stat_battery, &batstat)) != inst_ok) {
					printf("\nGetting instrument battery status failed with error :'%s' (%s)\n",
			       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
					it->del(it);
					return -1;
				}
				printf("The battery charged level is %.0f%%\n",batstat * 100.0);
			}

			if (scan_tol != 1.0) {
				if (cap2 & inst2_has_scan_toll) {
					inst_code ev;
					if ((ev = it->get_set_opt(it, inst_opt_scan_toll, scan_tol)) != inst_ok) {
						printf("\nSetting patch consistency tolerance to %f failed with error :'%s' (%s)\n",
				       	       scan_tol, it->inst_interp_error(it, ev), it->interp_error(it, ev));
						it->del(it);
						return -1;
					}
				} else {
					a1logv(log, 1, "Modified patch consistency tolerance ignored - instrument doesn't support it\n");
				}
			}

			/* Should look at instrument type & user spec ??? */
			if (trans) {
				if (pbypatch
				 && it->check_mode(it, inst_mode_trans_spot) == inst_ok) {
					mode = inst_mode_trans_spot;
					rmode = 0;
				} else if (it->check_mode(it, inst_mode_trans_chart) == inst_ok) {
					mode = inst_mode_trans_chart;
					rmode = 3;
				} else if (it->check_mode(it, inst_mode_trans_xy) == inst_ok) {
					mode = inst_mode_trans_xy;
					rmode = 2;
				} else if (it->check_mode(it, inst_mode_trans_strip) == inst_ok) {
					mode = inst_mode_trans_strip;
					rmode = 1;
				} else {
					mode = inst_mode_trans_spot;
					rmode = 0;
				}
			} else if (displ) {
				/* We assume a display mode will always be spot by spot */
				mode = inst_mode_emis_spot;
				rmode = 0;
			} else if (emis) {
				if (pbypatch
				 && it->check_mode(it, inst_mode_emis_spot) == inst_ok) {
					mode = inst_mode_emis_spot;
					rmode = 0;
				} else if (it->check_mode(it, inst_mode_emis_strip) == inst_ok) {
					mode = inst_mode_emis_strip;
					rmode = 1;
				} else {
					mode = inst_mode_emis_spot;
					rmode = 0;
				}
			} else {	/* Reflectance */
				inst_stat_savdrd sv = inst_stat_savdrd_none;

				/* See if instrument has a saved mode, and if it has data that */
				/* could match this chart */
				if (IMODETST(cap, inst_mode_s_reflection)) {

					a1logv(log, 2, "Instrument has a svaed chart mode\n");

					if ((rv = it->get_set_opt(it, inst_stat_saved_readings, &sv)) != inst_ok) {
						printf("Getting saved reading status failed with error :'%s' (%s)\n",
				       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
						it->del(it);
						return -1;
					}

					if (sv & inst_stat_savdrd_chart) {
						int no_patches, no_rows, pat_per_row, chart_id, missing_row;

						a1logv(log, 2, "There is a saved chart\n");

						if ((rv = it->get_set_opt(it, inst_stat_s_chart,
						    &no_patches, &no_rows, &pat_per_row, &chart_id, &missing_row))
							                                                      != inst_ok) {
							printf("Getting saved chart details failed with error :'%s' (%s)\n",
					       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
							it->del(it);
							return -1;
						}

						if (npat != no_patches || totpa != no_rows
						 || stipa != pat_per_row || rstart != chart_id) {
							printf("Can't use saved chart because it doesn't match\n");
							sv &= ~inst_stat_savdrd_chart;
						}

						if (missing_row >= 0) {
							printf("Can't use saved chart because row %d hasn't been read\n",missing_row);
							sv &= ~inst_stat_savdrd_chart;
						}
					}

					if (sv & inst_stat_savdrd_xy) {
						int nstr;
						int no_sheets, no_patches, no_rows, pat_per_row;

						/* Count the number of strips (sheets) */
						for (nstr = 0; pis[nstr] != 0; nstr++)
							;

						if ((rv = it->get_set_opt(it, inst_stat_s_xy,
						    &no_sheets, &no_patches, &no_rows, &pat_per_row)) != inst_ok) {
							printf("Getting saved sheet details failed with error :'%s' (%s)\n",
					       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
							it->del(it);
							return -1;
						}

						if (nstr != no_sheets || npat != no_patches
						 || totpa != no_rows || stipa != pat_per_row) {
							a1logv(log, 1, "Can't use saved sheets because they don't match chart\n");
							a1logv(log, 1, "Got %d sheets, expect %d. Got %d patches, expect %d.\n",
							       no_sheets,nstr,no_patches,npat);
							a1logv(log, 1, "Got %d rows, expect %d. Got %d patches per row, expect %d.\n",
							       no_rows,totpa,pat_per_row,stipa);
							sv &= ~inst_stat_savdrd_xy;
						}
					}

					if (sv & inst_stat_savdrd_strip) {
						int no_patches, no_rows, pat_per_row;

						if ((rv = it->get_set_opt(it, inst_stat_s_strip,
						    &no_patches, &no_rows, &pat_per_row)) != inst_ok) {
							printf("Getting saved strip details failed with error :'%s' (%s)\n",
					       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
							it->del(it);
							return -1;
						}

						if (npat != no_patches || totpa != no_rows || stipa != pat_per_row) {
							a1logv(log, 1, "Can't use saved strips because they don't match chart\n");
							sv &= ~inst_stat_savdrd_strip;
						}
					}

					if (sv & inst_stat_savdrd_spot) {
						int no_patches;

						if ((rv = it->get_set_opt(it, inst_stat_s_spot, &no_patches)) != inst_ok) {
							printf("Getting saved spot details failed with error :'%s' (%s)\n",
					       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
							it->del(it);
							return -1;
						}

						if (npat != no_patches) {
							a1logv(log, 1, "Can't use saved spots because they don't match chart - got %d patches, expect %d\n",no_patches,npat);
							sv &= ~inst_stat_savdrd_spot;
						}
					}
				}

				if (pbypatch
				        && it->check_mode(it, inst_mode_s_ref_spot) == inst_ok
				        && (sv & inst_stat_savdrd_spot)) {
					mode = inst_mode_s_ref_spot;
					svdmode = 1;
					rmode = 0;

				} else if (it->check_mode(it, inst_mode_s_ref_chart) == inst_ok
				        && (sv & inst_stat_savdrd_chart)) {
					mode = inst_mode_s_ref_chart;
					svdmode = 1;
					rmode = 3;

				} else if (it->check_mode(it, inst_mode_s_ref_xy) == inst_ok
				        && (sv & inst_stat_savdrd_xy)) {
					mode = inst_mode_s_ref_xy;
					svdmode = 1;
					rmode = 2;

				} else if (it->check_mode(it, inst_mode_s_ref_strip) == inst_ok
				        && (sv & inst_stat_savdrd_strip)) {
					mode = inst_mode_s_ref_strip;
					svdmode = 1;
					rmode = 1;

				} else if (it->check_mode(it, inst_mode_s_ref_spot) == inst_ok
				        && (sv & inst_stat_savdrd_spot)) {
					mode = inst_mode_s_ref_spot;
					svdmode = 1;
					rmode = 0;

				} else if (pbypatch 
				        && it->check_mode(it, inst_mode_ref_spot) == inst_ok) {
					mode = inst_mode_ref_spot;
					rmode = 0;

				} else if (it->check_mode(it, inst_mode_ref_chart) == inst_ok) {
					mode = inst_mode_ref_chart;
					rmode = 3;

				} else if (it->check_mode(it, inst_mode_ref_xy) == inst_ok) {
					mode = inst_mode_ref_xy;
					rmode = 2;

				} else if (it->check_mode(it, inst_mode_ref_strip) == inst_ok) {
					mode = inst_mode_ref_strip;
					rmode = 1;

				} else {
					mode = inst_mode_ref_spot;
					rmode = 0;
				}
			}
			if (spectral)
				mode |= inst_mode_spectral;

			if (highres) {
				if (IMODETST(cap, inst_mode_highres)) {
					mode |= inst_mode_highres;
				} else {
					a1logv(log, 1, "high resolution ignored - instrument doesn't support high res. mode\n");
				}
			}

			// ~~~ i1pro2 test code ~~~ */
			if (uvmode) {
				if (!IMODETST(cap, inst_mode_ref_uv)) {
					warning("UV measurement mode requested, but instrument doesn't support this mode");
					uvmode = 0;
				} else {
					mode |= inst_mode_ref_uv;
				}
			}

			if ((rv = it->set_mode(it, mode)) != inst_ok) {
				printf("\nSetting instrument mode failed with error :'%s' (%s)\n",
			     	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
				it->del(it);
				return -1;
			}
			it->capabilities(it, &cap, &cap2, &cap3);
		}
	}

	/* -------------------------------------------------- */
	if (rmode == 3) {		/* For chart mode, read all at once */
		int chid;			/* Chart ID number */
		ipatch *vals;		/* Array of values */

		if (svdmode)
			printf("Reading chart from values saved in instrument\n");
		else  {
			/* ~~999 ??? Need to setup trigger and wait for it appropriately ???  */
			printf("Reading the whole chart in one go\n");
		}

		/* Allocate space for patches */
		if ((vals = (ipatch *)calloc(sizeof(ipatch), npat)) == NULL)
			error("Malloc failed!");

		/* Initialise return values */
		for (i = 0; i < npat; i++) {
			strncpy(vals[i].loc, scols[i]->loc, ICOM_MAX_LOC_LEN-1);
			vals[i].loc[ICOM_MAX_LOC_LEN-1] = '\000';
			vals[i].XYZ_v = 0;
		}

		for (;;) {		/* retry loop */
			if ((rv = it->read_chart(it, npat, totpa, stipa, pis, rstart, vals)) == inst_ok)
				break;
			if (ierror(it, rv) == 0)	/* Ignore */
				continue;
			break;						/* Abort */
		}
		if (rv != inst_ok) {
			free(vals);
			it->del(it);
			return -1;
		}

		printf("Chart read OK\n");

		/* Transfer the values */
		/* We assume they are all in the right order */
		for (i = 0; i < npat; i++) {
			/* Copy XYZ */
			if (vals[i].XYZ_v == 0)
				error("Instrument didn't return XYZ value for patch %d, loc %s",i,scols[i]->loc);
			for (j = 0; j < 3; j++)
				scols[i]->XYZ[j] = vals[i].XYZ[j];

			/* Copy spectral */
			if (vals[i].sp.spec_n > 0) {
				scols[i]->sp = vals[i].sp;
			}
			scols[i]->mtype = vals[i].mtype;
			scols[i]->mcond = vals[i].mcond;
			scols[i]->rr = 1;		/* Has been read */
		}
		free(vals);

	/* -------------------------------------------------- */
	/* !!! Hmm. Should really allow user to navigate amongst the sheets, */
	/* !!! and skip any sheets already read. */
	} else if (rmode == 2) { /* For xy mode, read each sheet */
		ipatch *vals;
		int nsheets, sheet;		/* Total sheets/current sheet (sheet == pass) */
		int rpat = npat;		/* Remaining total patches */
		int pai;				/* Overall pass index */
		int sti;				/* Overall step index */
		char *pn[3] = { NULL, NULL, NULL} ;	/* Location 1/2/3 Pass name (letters) */
		char *sn[3] = { NULL, NULL, NULL} ;	/* Location 1/2/3 Step name (numbers) */
		int k;
		
		{	/* Figure the maximum number sheets and of patches in a sheet, for allocation */
			int lpaist = 0;			/* Largest number of passes in strip/sheet */
			for (nsheets = 0; pis[nsheets] != 0; nsheets++) {
				if (pis[nsheets] > lpaist)
					lpaist = pis[nsheets];
			}

			if ((vals = (ipatch *)calloc(sizeof(ipatch), (lpaist * stipa))) == NULL)
				error("Malloc failed!");
		}

		/* Make sure we can access the instrument table */
		if (cap2 & inst2_xy_holdrel) {
			it->xy_clear(it);
		}

		/* XY mode doesn't use the trigger mode */

		/* For each pass (==sheet) */
		for (sheet = 1, pai = sti = 0; pis[sheet-1] != 0; sheet++) {
			int paist;		/* Passes in current Strip (== columns in current sheet) */ 
			int rnpatch;	/* Rounded up (inc. padding) Patches in current pass (sheet) */
			int npatch;		/* Patches in pass (sheet), excluding padding */
			int fspipa;		/* First pass steps in pass */
			int nloc;		/* Number of fiducial locations needed */
			double ox = 0.0, oy = 0.0;		/* Origin */
			double ax = 1.0, ay = 0.0;		/* pass increment */
			double aax = 0.0, aay = 0.0;	/* pass offset for hex odd steps */
			double px = 0.0, py = 1.0;		/* step (==patch) increment */

			fspipa = stipa;
			paist = pis[sheet-1];			/* columns (letters) in sheet (strip) */
			npatch = rnpatch = paist * stipa;	/* Total patches including padding */
			if (npatch > rpat) {				/* This is a non-full pass */
				if (paist == 1) {
					fspipa -= (npatch - rpat);/* Last patch in first strip */
					if (fspipa < 1)
						error ("Assert in read_strips, fspipa = %d",fspipa);
				}
				npatch = rpat;				/* Total patches excluding padding */
			}

			nloc = 3;
			if (paist == 1) {
				nloc = 2;					/* Only one strip, so only 2 locations needed */
				if (fspipa == 1)
					nloc = 1;				/* Only one strip, one patch, so one location */
			}
			for (k = 0; k < 3; k++) {
				if (pn[k] != NULL)
					free(pn[k]);
				if (sn[k] != NULL)
					free(sn[k]);
			}
			pn[0] = paix->aix(paix, pai);			/* First pass (letter) */
			sn[0] = saix->aix(saix, 0);				/* First step (patch) (number) */
			pn[1] = paix->aix(paix, pai);			/* First pass (letter) */
			sn[1] = saix->aix(saix, 0 + fspipa-1);	/* Last step (patch) (number) */
			pn[2] = paix->aix(paix, pai + paist-1);	/* Last pass (letter) */
			sn[2] = saix->aix(saix, 0);				/* First step (patch) (number) */

			empty_con_chars();
			if (sheet == 1) {
				printf("Please place sheet %d of %d on table, then\n",sheet, nsheets);
			} else
				printf("\nPlease remove previous sheet, then place sheet %d of %d on table, then\n",sheet, nsheets);
				printf("hit return to continue, Esc or 'q' to give up"); fflush(stdout);
				if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
					printf("\n");
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->del(it);
					return -1;
				}
			printf("\n");

			if (cap2 & inst2_xy_holdrel) {

				/* Hold table */
				for (;;) {		/* retry loop */
					if ((rv = it->xy_sheet_hold(it)) == inst_ok)
						break;
					if (ierror(it, rv) == 0)	/* Ignore */
						continue;
					break;						/* Abort */
				}
				if (rv != inst_ok) {
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->del(it);
					return -1;
				}
			}

			if (cap2 & inst2_xy_locate) {
				int ll;
				double x[3], y[3];

				/* Allow user location of points */
				for (;;) {		/* retry loop */
					if ((rv = it->xy_locate_start(it)) == inst_ok)
						break;
					if (ierror(it, rv) == 0)	/* Ignore */
						continue;
					break;						/* Abort */
				}
				if (rv != inst_ok) {
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->del(it);
					return -1;
				}

				/* For each location point */
				for (ll = 0; ll < nloc; ll++) {
					empty_con_chars();
					printf("\nUsing the XY table controls, locate patch %s%s with the sight,\n",
					        pn[ll], sn[ll]);
					printf("then hit return to continue, Esc or 'q' to give up"); fflush(stdout);
						if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
							printf("\n");
							for (k = 0; k < 3; k++) {
								if (pn[k] != NULL)
									free(pn[k]);
								if (sn[k] != NULL)
									free(sn[k]);
							}
							it->del(it);
							return -1;
						}
					printf("\n");

					for (;;) {		/* retry loop */
						if ((rv = it->xy_get_location(it, &x[ll], &y[ll])) == inst_ok)
							break;
						if (ierror(it, rv) == 0)	/* Ignore */
							continue;
						break;						/* Abort */
					}
					if (rv != inst_ok) {
						for (k = 0; k < 3; k++) {
							if (pn[k] != NULL)
								free(pn[k]);
							if (sn[k] != NULL)
								free(sn[k]);
						}
						it->del(it);
						return -1;
					}
				}

				/* We're done with user control */
				for (;;) {		/* retry loop */
					if ((rv = it->xy_locate_end(it)) == inst_ok)
						break;
					if (ierror(it, rv) == 0)	/* Ignore */
						continue;
					break;						/* Abort */
				}
				if (rv != inst_ok) {
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->del(it);
					return -1;
				}

				/* Convert point locations into navigation values */
				ox = x[0];
				oy = y[0];
				if (hex) {
					double kk = sqrt(1.0/3.0);
					double nn = fspipa - 1.0;
					if (((fspipa-1) & 1) == 0) {	/* [0] & [1] are lined up */
						if (nloc == 3) {
							px = (x[1] - x[0])/nn;
							py = (y[1] - y[0])/nn;
							ax = (x[2] - x[0])/(paist-1);
							ay = (y[2] - y[0])/(paist-1);
							aax = 0.5 * ax;
							aay = 0.5 * ay;
						} else if (nloc == 2) {
							px = (x[1] - x[0])/nn;
							py = (y[1] - y[0])/nn;
							aax = kk * py;			/* Scale and rotate */
							aay = kk * -px;
						}
					} else {						/* [0] & [1] are offset by aa[xy] */
						if (nloc == 3) {
							ax = (x[2] - x[0])/(paist-1);
							ay = (y[2] - y[0])/(paist-1);
							aax = 0.5 * ax;
							aay = 0.5 * ay;
							px = (x[1] - x[0] - aax)/nn;
							py = (y[1] - y[0] - aay)/nn;
						} else if (nloc == 2) {
							px =  (nn * (x[1] - x[0]) - kk * (y[1] - y[0]))/(kk * kk + nn * nn);
							py =  (nn * (y[1] - y[0]) + kk * (x[1] - x[0]))/(kk * kk + nn * nn);
							aax = kk * py;			/* Scale and rotate */
							aay = kk * -px;
						}
					}

				} else {	/* Rectangular patches */
					if (paist > 1) {
						ax = (x[2] - x[0])/(paist-1);
						ay = (y[2] - y[0])/(paist-1);
					}
					if (fspipa > 1) {
						px = (x[1] - x[0])/(fspipa-1);
						py = (y[1] - y[0])/(fspipa-1);
					}
					aax = aay = 0.0;
				}
			}

			/* Read the sheets patches */
			for (;;) {		/* retry loop */
				if ((rv = it->read_xy(it, paist, stipa, npatch, pn[0], sn[0],
					                  ox, oy, ax, ay, aax, aay, px, py, vals)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok) {
				for (k = 0; k < 3; k++) {
					if (pn[k] != NULL)
						free(pn[k]);
					if (sn[k] != NULL)
						free(sn[k]);
				}
				it->del(it);
				return -1;
			}

			printf("Sheet %d of %d read OK\n",sheet, nsheets);

			/* Transfer the values */
			/* We assume they are all in the right order */
			for (i = 0; i < npatch; i++, sti++) {
				/* Copy XYZ */
				if (vals[i].XYZ_v == 0)
					error("Instrument didn't return XYZ value for patch %d, loc %s",i,scols[sti]->loc);
				for (j = 0; j < 3; j++)
					scols[sti]->XYZ[j] = vals[i].XYZ[j];

				/* Copy spectral */
				if (vals[i].sp.spec_n > 0) {
					scols[sti]->sp = vals[i].sp;
				}
				scols[sti]->rr = 1;		/* Has been read */
			}

			if (cap2 & inst2_xy_holdrel) {

				/* Release table and reset head */
				for (;;) {		/* retry loop */
					if ((rv = it->xy_clear(it)) == inst_ok)
						break;
					if (ierror(it, rv) == 0)	/* Ignore */
						continue;
					break;						/* Abort */
				}
				if (rv != inst_ok) {
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->del(it);
					return -1;
				}
			}

			pai += paist;		/* Tracj next first pass in strip */
			rpat -= npatch;		/* Track remaining patches */
		}	/* Next pass */

		for (k = 0; k < 3; k++) {
			if (pn[k] != NULL)
				free(pn[k]);
			if (sn[k] != NULL)
				free(sn[k]);
		}
		free(vals);

		printf("\nPlease remove last sheet from table\n"); fflush(stdout);

	/* -------------------------------------------------- */
	} else if (rmode == 1) {	/* For strip mode, simply read each strip */
		int uswitch = 0;		/* 0 if switch can be used, 1 if switch or keyboard */
		ipatch *vals;			/* Values read for a strip pass */
		int incflag = 0;		/* 0 = no change, 1 = increment, 2 = inc unread, */
								/* -1 = decrement, -2 = done */
		int stix;				/* Strip index */
		int pai;				/* Current pass in current strip */
		int oroi;				/* Overall row index */

		if (
		    itype != instDTP20 &&
	        !rand && disbidi == 2) {
			warning("Bi-directional strip recognition may not work without randomized patch locations");
		}

		/* Do any needed calibration before the user places the instrument on a desired spot */
		if (it->needs_calibration(it) & inst_calt_n_dfrble_mask) {
			if ((rv = inst_handle_calibrate(it, inst_calt_needed, inst_calc_none, NULL, NULL, 0))
			                                                                    != inst_ok) {
				printf("\nCalibration failed with error :'%s' (%s)\n",
	       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
				it->del(it);
				return -1;
			}
		}

		/* Enable switch or user via uicallback trigger if possible */
		if (cap2 & inst2_user_switch_trig) {
			rv = it->get_set_opt(it, inst_opt_trig_user_switch);
			uswitch = 2;

		/* Or use just switch trigger */
		} else if (cap2 & inst2_switch_trig) {
			rv = it->get_set_opt(it, inst_opt_trig_switch);
			uswitch = 1;

		/* Or go for user vi uicallback trigger */
		} else if (cap2 & inst2_user_trig) {
			rv = it->get_set_opt(it, inst_opt_trig_user);

		/* Or something is wrong with instrument capabilities */
		} else {
			printf("\nNo reasonable trigger mode avilable for this instrument\n");
			it->del(it);
			return -1;
		}
		if (rv != inst_ok) {
			printf("\nSetting trigger mode failed with error :'%s' (%s)\n",
	       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
			it->del(it);
			return -1;
		}

		/* Set so that return or any other key triggers, */
		/* but retain our abort keys */
		inst_set_uih(0x00, 0xff, DUIH_TRIG);
		inst_set_uih('f', 'f',   DUIH_CMND);
		inst_set_uih('F', 'F',   DUIH_CMND);
		inst_set_uih('b', 'b',   DUIH_CMND);
		inst_set_uih('B', 'B',   DUIH_CMND);
		inst_set_uih('n', 'n',   DUIH_CMND);
		inst_set_uih('N', 'N',   DUIH_CMND);
		inst_set_uih('g', 'g',   DUIH_CMND);
		inst_set_uih('G', 'G',   DUIH_CMND);
		inst_set_uih('d', 'd',   DUIH_CMND);
		inst_set_uih('D', 'D',   DUIH_CMND);
		inst_set_uih('q', 'q',   DUIH_ABORT);
		inst_set_uih('Q', 'Q',   DUIH_ABORT);
		inst_set_uih(0x03, 0x03, DUIH_ABORT);		/* ^c */
		inst_set_uih(0x1b, 0x1b, DUIH_ABORT);		/* Esc */

		/* Allocate space for values from a pass/strip */
		if ((vals = (ipatch *)calloc(sizeof(ipatch), (stipa+nextrap))) == NULL)
			error("Malloc failed!");

		/* Until we're done reading rows */
		incflag = 0;

		/* Skip to next unread if first has been read */
		if (pis[0] != 0 && scols[0]->rr != 0)
			incflag = 2;

		for (oroi = stix = pai = 0;pis[0] != 0;) {
			char *nn = NULL;	/* Pass name */
			int guide;
			chcol **scb;
			int boff = 0;			/* Best offset */
			int bdir = 0;			/* Best overall direction */
			int done = 0;			/* nz if there are no unread rows */

//printf("\n~1 incflag = %d, oroi %d, pai %d, stix %d\n", incflag, oroi, pai, stix);

			/* Increment or decrement to the next row */
			if (incflag > 0) {
				int s_oroi = oroi;

				/* Until we get to an unread pass */
				for (;;) {
					oroi++;
					if (++pai >= pis[stix]) {		/* Carry */
						if (pis[++stix] == 0) {		/* Carry */
							stix = 0;
							oroi = 0;
						}
						pai = 0;
					}
//printf("~1 stix = %d, pis[stix] = %d, oroi = %d, rr %d\n",stix, pis[stix],oroi,scols[oroi * stipa]->rr);
					// Note we aren't protecting agains a bodgy pis[] value
					if (incflag == 1 || scols[oroi * stipa]->rr == 0 || oroi == s_oroi)
						break;
				}

			/* Decrement the row */
			} else if (incflag < 0) {
				oroi--;
				if (--pai < 0) {			/* Carry */
					if (--stix < 0) {		/* Carry */
						for (oroi = stix = 0; pis[stix] != 0; stix++) {
							oroi += pis[stix];
						}
						stix--;
						oroi--;
					}
					pai = pis[stix]-1;
				}
			}
			incflag = 0;

			/* See if there are any unread patches */
			for (done = i = 0; i < npat; i += stipa) {
				if (scols[i]->rr == 0)
					break;					/* At least one patch read */
			}
			if (i >= npat)
				done = 1;
			
//printf("~1 oroi %d, pai %d, stix %d pis[stix] %d, rr = %d\n", oroi, pai, stix, pis[stix],scols[oroi * stipa]->rr);
			/* Convert overall pass number index into alpha label */
			if (nn != NULL)
				free(nn);
			nn = paix->aix(paix, oroi);

			guide = (pis[stix] - pai) * 5;		/* Mechanical guide offset */

			for (;;) {	/* Until we give up reading this row */

				/* Read a strip pass */
				printf("\nReady to read strip pass %s%s\n",nn, done ? " (!! ALL ROWS READ !!)" : scols[oroi *stipa]->rr ? " (This row has been read)" : "" );
				printf("Press 'f' to move forward, 'b' to move back, 'n' for next unread,\n");
				printf(" 'd' when done, Esc or 'q' to quit without saving.\n");
				
				if (uswitch == 1) {
					printf("Trigger instrument switch to start reading,");
				} else if (uswitch == 2) {
					printf("Trigger instrument switch or any other key to start:");
				} else {
					printf("Press any other key to start:");
				}
				fflush(stdout);
				if ((rv = it->read_strip(it, "STRIP", stipa+nextrap, nn, guide, plen, glen, tlen, vals)) != inst_ok
				 && (rv & inst_mask) != inst_user_trig) {

#ifdef DEBUG
					printf("read_strip returned '%s' (%s)\n",
				       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
					/* Deal with user abort or command */
					if ((rv & inst_mask) == inst_user_abort) {
						int keyc = inst_get_uih_char();

						/* Deal with a command */
						if (keyc & DUIH_CMND) {
							ch = keyc & 0xff;

							printf("\n");
							if (ch == 'f' || ch == 'F') {
								incflag = 1;
								break;
							} else if (ch == 'b' || ch == 'B') {
								incflag = -1;
								break;
							} else if (ch == 'n' || ch == 'N') {
								incflag = 2;
								break;
							} else {	/* Assume 'd' or 'D' */

								/* See if there are any unread patches */
								for (done = i = 0; i < npat; i += stipa) {
									if (scols[i]->rr == 0)
										break;					/* At least one patch read */
								}
								if (i >= npat)
									done = 1;
				
								if (done) {
									incflag = -2;
									break;
								}

								/* Not all rows have been read */
								empty_con_chars();
								printf("\nDone ? - At least one unread patch (%s, %s), Are you sure [y/n]: ",
								       scols[i]->id, scols[i]->loc);
								fflush(stdout);
								ch = next_con_char();
								printf("\n");
								if (ch == 'y' || ch == 'Y') {
									incflag = -2;
									break;
								}
								continue;
							}

						/* Deal with a user abort */
						} else if (keyc & DUIH_ABORT) {
							empty_con_chars();
							printf("\n\nStrip read stopped at user request!\n");
							printf("Hit Esc or 'q' to give up, any other key to retry:"); fflush(stdout);
							if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
								printf("\n");
								if (nn != NULL) free(nn);
								free(vals);
								it->del(it);
								return -1;
							}
							printf("\n");
							continue;
						}

					/* Deal with needs calibration */
					} else if ((rv & inst_mask) == inst_needs_cal) {
						inst_code ev;

						if (cap2 & inst2_no_feedback)
							bad_beep();
						printf("\nStrip read failed because instruments needs calibration\n");
						ev = inst_handle_calibrate(it, inst_calt_needed, inst_calc_none, NULL, NULL, 0);
						if (ev != inst_ok) {	/* Abort or fatal error */
							if (nn != NULL) free(nn);
							free(vals);
							it->del(it);
							return -1;
						}
						continue;

					/* Deal with a bad sensor position */
					} else if ((rv & inst_mask) == inst_wrong_config) {
						printf("\n\nSpot read failed due to the sensor being in the wrong position\n(%s)\n",it->interp_error(it, rv));
						continue;

					/* Deal with a misread */
					} else if ((rv & inst_mask) == inst_misread) {
						if (cap2 & inst2_no_feedback)
							bad_beep();
						empty_con_chars();
						printf("\nStrip read failed due to misread (%s)\n",it->interp_error(it, rv));
						printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
						if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
							printf("\n");
							if (nn != NULL) free(nn);
							free(vals);
							it->del(it);
							return -1;
						}
						printf("\n");
						continue;
					/* Deal with a communications error */
					} else if ((rv & inst_mask) == inst_coms_fail) {
						if (cap2 & inst2_no_feedback)
							bad_beep();
						empty_con_chars();
						printf("\nStrip read failed due to communication problem.\n");
						printf("Hit Esc or 'q' to give up, any other key to retry:"); fflush(stdout);
						if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
							printf("\n");
							if (nn != NULL) free(nn);
							free(vals);
							it->del(it);
							return -1;
						}
						printf("\n");
						if ((it->icom->port_type(it->icom) & icomt_serial)
						 && !(it->icom->port_attr(it->icom) & icomt_fastserial)) {
							/* Allow retrying at a lower baud rate */
							int tt = it->last_scomerr(it);
							if (tt & (ICOM_BRK | ICOM_FER | ICOM_PER | ICOM_OER)) {
								if (br == baud_57600) br = baud_38400;
								else if (br == baud_38400) br = baud_9600;
								else if (br == baud_9600) br = baud_4800;
								else if (br == baud_9600) br = baud_4800;
								else if (br == baud_2400) br = baud_1200;
								else br = baud_1200;
							}
							if ((rv = it->init_coms(it, br, fc, 15.0)) != inst_ok) {
#ifdef DEBUG
								printf("init_coms returned '%s' (%s)\n",
							       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
								if (nn != NULL) free(nn);
								free(vals);
								it->del(it);
								return -1;
							}
						}
						continue;
					} else {
						/* Some other error. Treat it as fatal */
						if (cap2 & inst2_no_feedback)
							bad_beep();
						printf("\nStrip read failed due unexpected error :'%s' (%s)\n",
				       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
						printf("Hit Esc or 'q' to give up, any other key to retry:"); fflush(stdout);
						if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
							printf("\n");
							if (nn != NULL) free(nn);
							free(vals);
							it->del(it);
							return -1;
						}
						printf("\n");
						continue;
					}

				/* Successfully read the strip */
				/* See which expected row best correlates with the one we've read. */
				/* Figure out if there is an "off by one" error for a DTP51 */
				} else {
					int choroi;				/* Check overall row index */
					double corr;			/* Correlation with expected value */
					int loff = 0, hoff = 0;	/* DTP51 offset test range */
					int toff;				/* Test offset */
					int dir, dirrg = 1;		/* Direction range, 1 = forward, 2 = fwd & bwd */

					int boroi = -1;			/* Best overall row index */

					double bcorr = 1e9;		/* Best correlation value */
					double werror = 0.0;	/* Worst case error in best correlation strip */

					double xbcorr = 1e9;	/* Expected pass correlation value */
					int xboff;				/* Expected pass best offset */
					int xbdir;				/* Expected pass overall pass direction */
					double xwerror = 0.0;	/* Expected pass worst patcch error */

					if (((rand && disbidi == 0) || disbidi == 2)
					 && (cap2 & inst2_bidi_scan))
						dirrg = 2;			/* Enable bi-directional strip recognition */

					/* DTP51 has a nasty habit of misaligning test squares by +/- 1 */
					/* See if this might have happened */
					if (it->dtype == instDTP51) {
						loff = -1;
						hoff = 1;
					}

					/* Explor the rows */
					for (choroi = 0; choroi < totpa; choroi++) {

						/* Explore strip direction */
						for (dir = 0; dir < dirrg; dir++) {
							double pwerr;	/* This rows worst error */
							scb = &scols[choroi * stipa];

							/* Explore off by +/-1 error for DTP51 */
							for (toff = loff; toff <= hoff; toff++) {
								double ynorm = 1.0;

								/* Compute a Y scaling value to give correlation */
								/* a chance for absolute readings */
								if (vals[skipp+toff].XYZ_v != 0) {
									double refnorm = 0.0;
									ynorm = 0.0;
									for (i = 0; i < stipa; i++) {
										int ix = i+skipp+toff;
										if (dir != 0)
											ix = stipa - 1 - ix;
										refnorm += scb[i]->eXYZ[1];
										ynorm += vals[ix].XYZ[1];
									}
									if (fabs(ynorm) > 1e-6)
										ynorm = refnorm/ynorm;
									else
										ynorm = 1.0;
								}

								/* Compare just sample patches (not padding Max/Min) */
#ifdef USESTRDELTA
								for (pwerr = corr = 0.0, n = 0, i = 0; i < (stipa-1); i++, n++)
#else
								for (pwerr = corr = 0.0, n = 0, i = 0; i < stipa; i++, n++)
#endif
								{
									double vcorr;
									int ix = i+skipp+toff, ix1 = ix+1;
									if (dir != 0) {
										ix = stipa - 1 - ix;
										ix1 = stipa - 1 - ix1;
									}
									if (vals[ix].XYZ_v == 0)
										error("Instrument didn't return XYZ value");
#ifdef USESTRDELTA
									vcorr = xyzLabcorr(vals[ix].XYZ, scb[i]->eXYZ,
									                   vals[ix1].XYZ, scb[i+1]->eXYZ);
#else
									vcorr = xyzLabDE(ynorm, vals[ix].XYZ, scb[i]->eXYZ);
#endif
//printf("DE %f from vals[%d] %f %f %f and scols[%d] %f %f %f\n", vcorr, ix, vals[ix].XYZ[0], vals[ix].XYZ[1], vals[ix].XYZ[2], i + choroi * stipa, scb[i]->eXYZ[0], scb[i]->eXYZ[1], scb[i]->eXYZ[2]);
									corr += vcorr;
									if (vcorr > pwerr)	/* Worsed patch error */
										pwerr = vcorr;
								}
								corr /= (double)n;
#ifdef DEBUG
								printf("  Strip %d dir %d offset %d correlation = %f\n",choroi,dir,toff,corr);

#endif
								/* If this is the expected strip, */
								/* note correlation and best fit to off by 1 and direction */
								if (choroi == oroi && corr < xbcorr) { 
									xbcorr = corr;		/* Expected strip correlation */
									xboff = toff;
									xbdir = dir;
									xwerror = pwerr;	/* Expected passes worst error */
								}

								/* Best overall matched strip correlation */
								if (corr < bcorr) {
									boroi = choroi;
									bcorr = corr;
									boff = toff;
									bdir = dir;
									werror = pwerr;
								}
							}
						}
					}
					if (emit_warnings != 0 && rand && boroi != oroi) {	/* Looks like the wrong strip */
						char *mm = NULL;
						mm = paix->aix(paix, boroi);
#ifdef DEBUG
						printf("Strip pass %s (%d) seems to have a better correlation that strip %s (%d)\n",
						mm, boroi, nn, oroi);
#endif
						if (cap2 & inst2_no_feedback)
							bad_beep();
						empty_con_chars();
						printf("\n(Warning) Seem to have read strip pass %s rather than %s!\n",mm,nn);
						printf("Hit Return to use it anyway, any other key to retry, Esc or 'q' to give up:"); fflush(stdout);
						free(mm);
						if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
							printf("\n");
							if (nn != NULL) free(nn);
							free(vals);
							it->del(it);
							return -1;
						}
						if (ch != 0x0d && ch != 0x0a) { 			/* !(CR or LF) */
							printf("\n");
							continue;		/* Try again */
						}
						printf("\n");

						/* Switch to state for expected strip */
						bcorr = xbcorr;
						boff = xboff;
						bdir = xbdir;
						werror = xwerror;
					}
					/* Arbitrary threshold. Good seems about 15-35, bad 95-130 */
					if (emit_warnings != 0 && accurate_expd != 0 && werror >= 30.0) {
#ifdef DEBUG
						printf("(Warning) Patch error %f (>35 not good, >95 bad)\n",werror);
#endif
						if (cap2 & inst2_no_feedback)
							bad_beep();
						empty_con_chars();
						printf("\nThere is at least one patch with an very unexpected response! (DeltaE %f)\n",werror);
						printf("Hit Return to use it anyway, any other key to retry, Esc or  'q' to give up:"); fflush(stdout);
						if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
							printf("\n");
							if (nn != NULL) free(nn);
							free(vals);
							it->del(it);
							return -1;
						}
						if (ch != 0x0d && ch != 0x0a) { 			/* !Cr */
							printf("\n");
							continue;
						}
						printf("\n");
						break;
					}

					/* Must be OK - save the readings */
					if (cap2 & inst2_no_feedback)
						good_beep();
					printf(" Strip read OK");
					if (boff != 0)
						printf(" (DTP51 offset fix of %d applied)",boff);
					if (bdir != 0)
						printf(" (Strip read in reverse direction)");
					printf("\n");
					break;		/* Break out of retry loop */
				}
			}

			if (nn != NULL)		/* Finished with strip alpha index */
				free(nn);
			nn = NULL;

			/* If we're done */
			if (incflag == -2)
				break;

			/* If we are moving row, rather than having read one. */
			if (incflag != 0)
				continue;

			/* Transfer the values (including DTP51 offset) */
			scb = &scols[oroi * stipa];
			for (n = 0, i = 0; i < stipa; i++, n++) {
				int ix = i+skipp+boff;
				if (bdir != 0)
					ix = stipa - 1 - ix;

				/* Copy XYZ */
				if (vals[ix].XYZ_v == 0)
					error("Instrument didn't return XYZ value");
				for (j = 0; j < 3; j++)
					scb[i]->XYZ[j] = vals[ix].XYZ[j];

				/* Copy spectral */
				if (vals[ix].sp.spec_n > 0) {
					scb[i]->sp = vals[ix].sp;
				}
				scb[i]->rr = 1;		/* Has been read */
			}
			incflag = 2;		/* Skip to next unread */
		}		/* Go around to read another row */
		free(vals);

	/* -------------------------------------------------- */
	/* Spot mode. This will be used if xtern != 0 */
	} else {
		int pix = 0;
		int uswitch = 0;		/* nz if switch can be used */
		int incflag = 0;		/* 0 = no change, 1 = increment, 2 = inc by 10, */
								/* 3 = inc next unread, -1 = decrement, -2 = dec by 10 */
								/* 4 = goto specific patch */
		inst_opt_type omode;	/* The option mode used */
		ipatch val;

		if (xtern == 0) {	/* Instrument patch by patch */

			/* Do any needed calibration before the user places the instrument on a desired spot */
			if (it->needs_calibration(it) & inst_calt_n_dfrble_mask) {
				if ((rv = inst_handle_calibrate(it, inst_calt_needed, inst_calc_none, NULL, NULL, 0))
				                                                                    != inst_ok) {
					printf("\nCalibration failed with error :'%s' (%s)\n",
		       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
					it->del(it);
					return -1;
				}
			}
	
			/* Enable switch or user via uicallback trigger if possible */
			if (cap2 & inst2_user_switch_trig) {
				omode = inst_opt_trig_user_switch;
				rv = it->get_set_opt(it, omode);
				uswitch = 1;
	
			/* Or go for user via uicallback trigger */
			} else if (cap2 & inst2_user_trig) {
				omode = inst_opt_trig_user;
				rv = it->get_set_opt(it, omode);
	
			/* Or something is wrong with instrument capabilities */
			} else {
				printf("\nNo reasonable trigger mode avilable for this instrument\n");
				it->del(it);
				return -1;
			}
			if (rv != inst_ok) {
				printf("\nSetting trigger mode failed with error :'%s' (%s)\n",
		       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
				it->del(it);
				return -1;
			}

			/* Setup the keyboard trigger to return our commands */
			inst_set_uih('f', 'f', DUIH_CMND);
			inst_set_uih('F', 'F', DUIH_CMND);
			inst_set_uih('b', 'b', DUIH_CMND);
			inst_set_uih('B', 'B', DUIH_CMND);
			inst_set_uih('n', 'n', DUIH_CMND);
			inst_set_uih('N', 'N', DUIH_CMND);
			inst_set_uih('g', 'g', DUIH_CMND);
			inst_set_uih('G', 'G', DUIH_CMND);
			inst_set_uih('d', 'd', DUIH_CMND);
			inst_set_uih('D', 'D', DUIH_CMND);
			inst_set_uih('k', 'k',   DUIH_CMND);
			inst_set_uih('q', 'q', DUIH_ABORT);
			inst_set_uih('Q', 'Q', DUIH_ABORT);
			inst_set_uih(0xd, 0xd, DUIH_TRIG);	/* Return */
			inst_set_uih(' ', ' ', DUIH_TRIG);
		}

		/* Skip to next unread if first has been read */
		/* !!! would be nice to skip padding patches !!! */
		incflag = 0;
		if (npat > 0 && scols[0]->rr != 0)
			incflag = 3;

		/* Until we're  done */
		for (;pix < npat;) {
			char buf[200], *bp = NULL, *ep = NULL;
			char ch = 0;

			/* Adjust the location */
			if (incflag > 0 && incflag <= 2) {	/* Incremente by 1 or 10 */

				if (incflag == 2)
					pix += 10;
				else
					pix++;
				pix = pix % npat;

			} else if (incflag < 0 && incflag >= -2) {	/* Decrement by 1 or 10 */

				if (incflag == -2)
					pix -= 10;
				else
					pix--;
				pix = pix % npat;
				if (pix < 0)
					pix += npat;

			} else if (incflag == 3) {		/* Increment to next unread */
				int opix = pix;

				if (pix >= npat)
					pix = 0;
				for (;;) {
					if (scols[pix]->rr == 0 && strcmp(scols[pix]->id, "0") != 0)
						break;
					pix++;
					if (pix >= npat)
						pix = 0;
					if (pix == opix)
						break;
				}
			} else if (incflag == 4) {		/* Goto specific patch */
				printf("\nEnter patch to go to: "); fflush(stdout);

				/* Read in the next line from stdin. */
				if (fgets(buf, 200, stdin) == NULL) {
					printf("Error - unrecognised input\n");
				} else {
					int opix = pix;

					/* Skip whitespace */
					for (bp = buf; *bp != '\000' && isspace(*bp); bp++)
						;

					/* Skip non-whitespace */
					for (ep = bp; *ep != '\000' && !isspace(*ep); ep++)
						;
					*ep = '\000';
				
					if (pix >= npat)
						pix = 0;
					for (;;) {
						if (stricmp(scols[pix]->loc, bp) == 0)
							break;
						pix++;
						if (pix >= npat)
							pix = 0;
						if (pix == opix) {
							printf("Patch '%s' not found\n",bp);
							break;
						}
					}
				}
			}
			incflag = 0;

			/* See if there are any unread patches */
			for (i = 0; i < npat; i++) {
				if (scols[i]->rr == 0 && strcmp(scols[i]->id, "0") != 0)
					break;					/* At least one patch read */
			}

			if (xtern != 0) {	/* User entered values */
				printf("\nReady to read patch '%s' at '%s'%s\n",scols[pix]->id, scols[pix]->loc,
				       i >= npat ? "(All patches read!)" :
				       strcmp(scols[pix]->id, "0") == 0 ? " (Padding Patch)" :
				       scols[pix]->rr ? " (Already read)" : "");
				printf("Enter %s value (separated by spaces), or\n",
				       xtern == 1 ? "L*a*b*" : "XYZ");
				printf("  'f' to move forward, 'F' move forward 10,\n");
				printf("  'b' to move back, 'B; to move back 10,\n");
				printf("  'n' for next unread, 'g' to goto patch,\n");
				printf("  'd' when done, 'q' to abort, then press <return>: ");
				fflush(stdout);
	
				/* Read in the next line from stdin. */
				if (fgets(buf, 200, stdin) == NULL) {
					printf("Error - unrecognised input\n");
					continue;
				}
				/* Skip whitespace */
				for (bp = buf; *bp != '\000' && isspace(*bp); bp++)
					;

				ch = *bp;
				if (ch == '\000') {
					printf("Error - unrecognised input\n");
					continue;
				}

			} else {	/* Using instrument */

				empty_con_chars();

				printf("\nReady to read patch '%s' at '%s'%s\n",scols[pix]->id, scols[pix]->loc,
				       i >= npat ? " (All patches read!)" :
				       strcmp(scols[pix]->id, "0") == 0 ? " (Padding Patch)" :
				       scols[pix]->rr ? " (Already read)" : "");

				printf("hit 'f' to move forward, 'F' move forward 10,\n");
				printf("    'b' to move back,    'B; to move back 10,\n");
				printf("    'n' for next unread, 'g' to goto patch,\n");
				printf("    'd' when done,       'k' to calibrate, <esc> to abort,\n");

				if (uswitch)
					printf("    Instrument switch,   <return> or <space> to read:");
				else
					printf("    <return> or <space> to read:");
				fflush(stdout);

				rv = it->read_sample(it, "SPOT", &val, 1);

				/* Deal with reading */
				if (rv == inst_ok) {
					/* Read OK */
					if (cap2 & inst2_no_feedback)
						good_beep();
					ch = '0';

				/* Deal with user trigger via user interface callback function */
				} else if ((rv & inst_mask) == inst_user_trig) {
					if (cap2 & inst2_no_feedback)
						good_beep();
					ch = inst_get_uih_char();

				/* Deal with a abort or command */
				} else if ((rv & inst_mask) == inst_user_abort) {
					int keyc = inst_get_uih_char();

					/* User issued a command */
					if (keyc & DUIH_CMND) {
						ch = keyc & 0xff;
						printf("\n");

					/* User aborted */
					} else if (keyc & DUIH_ABORT) {
						empty_con_chars();
						printf("\n\nSpot read stopped at user request!\n");
						printf("Hit Esc or 'q' to give up, any other key to retry:"); fflush(stdout);
						if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
							printf("\n");
							it->del(it);
							return -1;
						}
						printf("\n");
						continue;
					}

				/* Deal with needs calibration */
				} else if ((rv & inst_mask) == inst_needs_cal) {
					inst_code ev;

					if (cap2 & inst2_no_feedback)
						bad_beep();
					printf("\nSpot read failed because instruments needs calibration\n");
					ev = inst_handle_calibrate(it, inst_calt_needed, inst_calc_none, NULL, NULL, 0);
					if (ev != inst_ok) {	/* Abort or fatal error */
						it->del(it);
						return -1;
					}
					continue;
				/* Deal with a misread */
				} else if ((rv & inst_mask) == inst_misread) {
					if (cap2 & inst2_no_feedback)
						bad_beep();
					empty_con_chars();
					printf("\nStrip read failed due to misread (%s)\n",it->interp_error(it, rv));
					printf("Hit Esc or 'q' to give up, any other key to retry:"); fflush(stdout);
					if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
						printf("\n");
						it->del(it);
						return -1;
					}
					printf("\n");
					continue;
				/* Deal with a communications error */
				} else if ((rv & inst_mask) == inst_coms_fail) {
					if (cap2 & inst2_no_feedback)
						bad_beep();
					empty_con_chars();
					printf("\nStrip read failed due to communication problem.\n");
					printf("Hit Esc or 'q' to give up, any other key to retry:"); fflush(stdout);
					if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
						printf("\n");
						it->del(it);
						return -1;
					}
					printf("\n");
					if ((it->icom->port_type(it->icom) & icomt_serial)
					 && !(it->icom->port_attr(it->icom) & icomt_fastserial)) {
						/* Allow retrying at a lower baud rate */
						int tt = it->last_scomerr(it);
						if (tt & (ICOM_BRK | ICOM_FER | ICOM_PER | ICOM_OER)) {
							if (br == baud_57600) br = baud_38400;
							else if (br == baud_38400) br = baud_9600;
							else if (br == baud_9600) br = baud_4800;
							else if (br == baud_9600) br = baud_4800;
							else if (br == baud_2400) br = baud_1200;
							else br = baud_1200;
						}
						if ((rv = it->init_coms(it, br, fc, 15.0)) != inst_ok) {
#ifdef DEBUG
							printf("init_coms returned '%s' (%s)\n",
						       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
							it->del(it);
							return -1;
						}
					}
					continue;


				} else {
					/* Some other error. Treat it as fatal */
					if (cap2 & inst2_no_feedback)
						bad_beep();
					printf("\nPatch read failed due unexpected error :'%s' (%s)\n",
			       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
					printf("Hit Esc or 'q' to give up, any other key to retry:"); fflush(stdout);
					if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
						printf("\n");
						it->del(it);
						return -1;
					}
					printf("\n");
					continue;
				}
			}
	
			if (ch == 'q' || ch == 0x1b || ch == 0x03) {	/* q or Esc or ^C */
				empty_con_chars();
				printf("\nAbort ? - Are you sure ? [y/n]:"); fflush(stdout);
				if ((ch = next_con_char()) == 'y' || ch == 'Y') {
					printf("\n");
					it->del(it);
					return -1;
				}
				printf("\n");
				continue;
			} else if (ch == 'k') {
				inst_code ev;

				ev = inst_handle_calibrate(it, inst_calt_available, inst_calc_none, NULL, NULL, 0);
				if (ev != inst_ok) {	/* Abort or fatal error */
					it->del(it);
					return -1;
				}
				continue;
			} else if (ch == 'f') {
				incflag = 1;
				continue;
			} else if (ch == 'F') {
				incflag = 2;
				continue;
			} else if (ch == 'b') {
				incflag = -1;
				continue;
			} else if (ch == 'B') {
				incflag = -2;
				continue;
			} else if (ch == 'n' || ch == 'N') {
				incflag = 3;
				continue;
			} else if (ch == 'g' || ch == 'G') {
				incflag = 4;
				continue;
			} else if (ch == 'd' || ch == 'D') {
				int i;
				for (i = 0; i < npat; i++) {
					if (scols[i]->rr == 0 && strcmp(scols[i]->id, "0") != 0)
						break;
				}
				if (i >= npat)
					break;			/* None unread, so done */

				/* Not all patches have been read */
				empty_con_chars();
				printf("\nDone ? - At least one unread patch (%s, %s), Are you sure [y/n]: ",
				       scols[i]->id, scols[i]->loc);
				fflush(stdout);
				if ((ch = next_con_char()) == 0x1b) {
					printf("\n");
					it->del(it);
					return -1;
				}
				printf("\n");
				if (ch == 'y' || ch == 'Y')
					break;
				continue;

			/* Read the external sample */
	        } else if (xtern != 0 && (isdigit(*bp) || ch == '-' || ch == '+' || ch == '.')) {
				int i;

				/* For each input number */
				for (i = 0; *bp != '\000' && i < 3; i++) {
					char *tp, *nbp;

					/* Find the start of the number */
					while(*bp != '\000' && !isdigit(*bp)
					      && *bp != '-' && *bp != '+' && *bp != '.')
						bp++;
					if (!isdigit(*bp) && *bp != '-' && *bp != '+' && *bp != '.')
						break;

					/* Find the end of the number */
					for (tp = bp+1; isdigit(*tp) || *tp == 'e' || *tp == 'E'
					                || *tp == '-' || *tp == '+' || *tp == '.'; tp++)
						;
					if (*tp != '\000')
						nbp = tp+1;
					else
						nbp = tp;
					*tp = '\000';

					/* Read the number */
					scols[pix]->XYZ[i] = atof(bp);

					bp = nbp;
				}
				if (i < 3) {	/* Didn't find 3 numbers */
					printf("Error - unrecognised input\n");
					continue;
				}
				if (xtern == 1) {
					icmLab2XYZ(&icmD50, scols[pix]->XYZ,scols[pix]->XYZ);
					scols[pix]->XYZ[0] *= 100.0;
					scols[pix]->XYZ[1] *= 100.0;
					scols[pix]->XYZ[2] *= 100.0;
				}

				scols[pix]->rr = 1;		/* Has been read */
				printf(" Got XYZ value %f %f %f\n",scols[pix]->XYZ[0], scols[pix]->XYZ[1], scols[pix]->XYZ[2]);

				/* Advance to next patch. */
				incflag = 1;
				
			/* We've read the spot sample */
			} else if (xtern == 0 && (ch == '0' || ch == ' ' || ch == '\r')) {	

				/* Save the reading */
				if (val.XYZ_v == 0)
					error("Instrument didn't return XYZ value");

				for (j = 0; j < 3; j++)
					scols[pix]->XYZ[j] = val.XYZ[j];

				/* Copy spectral */
				if (val.sp.spec_n > 0) {
					scols[pix]->sp = val.sp;
				}
				scols[pix]->rr = 1;		/* Has been read */
				printf(" Patch read OK\n");

#ifndef SALONEINSTLIB
				if (doplot && val.sp.spec_n > 0)
					xspect_plot_w(&val.sp, NULL, NULL, 0);
#endif

				/* Advance to next patch. */
				incflag = 1;
			} else {	/* Unrecognised response */
				continue;
			}
		}
	}
	/* clean up */
	if (it != NULL)
		it->del(it);
	return 0;
}

void
usage() {
	icompaths *icmps;
	inst2_capability cap2 = 0;
	fprintf(stderr,"Read Target Test Chart, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: chartread [-options] outfile\n");
	fprintf(stderr," -v             Verbose mode\n");
	fprintf(stderr," -c listno      Set communication port from the following list (default %d)\n",COMPORT);
	if ((icmps = new_icompaths(g_log)) != NULL) {
		icompath **paths;
		if ((paths = icmps->paths) != NULL) {
			int i;
			for (i = 0; ; i++) {
				if (paths[i] == NULL)
					break;
				fprintf(stderr,"    %d = '%s'\n",i+1,paths[i]->name);
			}
		} else
			fprintf(stderr,"    ** No ports found **\n");
	}
	fprintf(stderr," -t              Use transmission measurement mode\n");
	fprintf(stderr," -d              Use display measurement mode (white Y relative results)\n");
	cap2 = inst_show_disptype_options(stderr, " -y              ", icmps, 0, 0);
	fprintf(stderr," -e              Emissive for transparency on a light box\n");
	fprintf(stderr," -p              Measure patch by patch rather than strip\n");
	fprintf(stderr," -x [lx]         Take external values, either L*a*b* (-xl) or XYZ (-xx).\n");
	fprintf(stderr," -n              Don't save spectral information (default saves spectral)\n");
	fprintf(stderr," -l              Save CIE as D50 L*a*b* rather than XYZ\n");
	fprintf(stderr," -L              Save CIE as D50 L*a*b* as well as XYZ\n");
	fprintf(stderr," -r              Resume reading partly read chart\n");
	fprintf(stderr," -I file.cal     Override calibration info from .ti2 in resulting .ti3\n");
	fprintf(stderr," -F filter            Set filter configuration (if aplicable):\n");
	fprintf(stderr,"    n                  None (M0)\n");
	fprintf(stderr,"    5                  D50 (M1)\n");
	fprintf(stderr,"    6                  D65\n");
	fprintf(stderr,"    u                  U.V. Cut (M2)\n");
	fprintf(stderr,"    p                  Polarising filter (M3)\n");
	fprintf(stderr," -A N|A|X|G      XRGA conversion (default N)\n");
	fprintf(stderr," -N              Disable initial calibration of instrument if possible\n");
	fprintf(stderr," -B              Disable auto bi-directional strip recognition\n");
	fprintf(stderr," -b              Force auto bi-directional strip recognition\n");
	fprintf(stderr," -H              Use high resolution spectrum mode (if available)\n");
	if (cap2 & inst2_ccmx)
		fprintf(stderr," -X file.ccmx    Apply Colorimeter Correction Matrix\n");
	if (cap2 & inst2_ccss) {
		int i;
		fprintf(stderr," -X file.ccss    Use Colorimeter Calibration Spectral Samples for calibration\n");
		fprintf(stderr," -Q observ       Choose CIE Observer for CCSS instrument:\n");
#ifdef SALONEINSTLIB
		fprintf(stderr,"                  1931_2 (def), 1964_10, 2012_2, 2012_10\n");
#else /* !SALONEINSTLIB */
		fprintf(stderr,"                  1931_2 (def), 1964_10, 2012_2, 2012_10, S&B 1955_2, shaw, J&V 1978_2 or file.cmf\n");
#endif /* !SALONEINSTLIB */
	}
	fprintf(stderr," -T ratio        Modify strip patch consistency tolerance by ratio\n");
	fprintf(stderr," -S              Suppress wrong strip & unexpected value warnings\n");
//	fprintf(stderr," -Y U                 Test i1pro2 UV measurement mode\n");
	fprintf(stderr," -W n|h|x        Override serial port flow control: n = none, h = HW, x = Xon/Xoff\n");
#ifndef SALONEINSTLIB
	fprintf(stderr," -P              Plot spectral if patch by patch\n");
#endif
	fprintf(stderr," -D [level]      Print debug diagnostics to stderr\n");
	fprintf(stderr," outfile         Base name for input[ti2]/output[ti3] file\n");
	if (icmps != NULL)
		icmps->del(icmps);
	exit(1);
}

int main(int argc, char *argv[]) {
	int i, j;
	int fa, nfa, mfa;				/* current argument we're looking at */
	int verb = 0;
	int debug = 0;
	int comport = COMPORT;			/* COM port used */
	icompaths *icmps = NULL;
	icompath *ipath = NULL;
	flow_control fc = fc_nc;		/* Default flow control */
	instType itype = instUnknown;	/* Instrument chart is targeted to */
	instType atype = instUnknown;	/* Instrument used to read the chart */
	int trans = 0;					/* Use transmission mode */
	int emis = 0;					/* Use emissive mode */
	int displ = 0;					/* 1 = Use display emissive mode, 2 = display bright rel. */
	                                /* 3 = display white rel. */
	int ditype = 0;					/* Display type selection charater(s) */
	inst_opt_filter fe = inst_opt_filter_unknown;
	int pbypatch = 0;				/* Read patch by patch */
	int disbidi = 0;				/* 1 = Disable bi-directional strip recognition, 2 = force enable */
	int highres = 0;				/* Use high res mode if available */
	double scan_tol = 1.0;			/* Patch consistency tolerance modification */
	int xtern = 0;					/* Take external values, 1 = Lab, 2 = XYZ */
	int spectral = 1;				/* Save spectral information */
	int uvmode = 0;					/* ~~~ i1pro2 test mode ~~~ */
	xcalstd scalstd = xcalstd_none;	/* X-Rite calibration standard to set */
	xcalstd ucalstd = xcalstd_none;	/* X-Rite calibration standard actually used */
	int accurate_expd = 0;			/* Expected value assumed to be accurate */
	int emit_warnings = 1;			/* Emit warnings for wrong strip, unexpected value */
	int dolab = 0;					/* 1 = Save CIE as Lab, 2 = Save CIE as XYZ and Lab */
	int doresume = 0;				/* Resume reading a chart */
	int nocal = 0;					/* Disable initial calibration */
	int doplot = 0;					/* Plot spectral of patch by patch */
	char ccxxname[MAXNAMEL+1] = "\000";  /* Colorimeter Correction/Colorimeter Calibration name */
	icxObserverType obType = icxOT_default;		/* ccss observer */
	xspect custObserver[3];			/* If obType = icxOT_custom */
	static char inname[MAXNAMEL+1] = { 0 };	/* Input cgats file base name */
	static char outname[MAXNAMEL+1] = { 0 };	/* Output cgats file base name */
	cgats *icg;					/* input cgats structure */
	cgats *ocg;					/* output cgats structure */
	static char calname[MAXNAMEL+1] = { 0 };	/* User supplied calibration filename */
	xcal *cal = NULL;			/* Any calibration to be output as well */
	int nmask = 0;				/* Device colorant mask */
	char *pixpat = "A-Z, A-Z";			/* Pass index pattern */		
	char *sixpat = "0-9,@-9,@-9;1-999";	/* Step index pattern */		
	alphix *paix, *saix;		/* Pass and Step index generators */
	int ixord = 0;				/* Index order, 0 = pass then step */
	int rstart = 0;				/* Random start/chart id */
	int rand = 0;				/* Random patch order, - can use auto strip ID and Bi-Di */
	int hex = 0;				/* Hexagon pattern layout */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	chcol *cols;				/* Internal storage of all the patch colors */
	chcol **scols;				/* Location sorted pointers to cols */
	int nchan = 0;				/* Number of device chanels */
	int npat;					/* Number of overall patches */
	int *pis;					/* Passes in eachstrip, zero terminated */
	int stipa;					/* Steps in each Pass */
	int totpa;					/* Total Passes Needed */
	int runpat;					/* Rounded Up to (totpa * stipa) Number of patches */
	int wpat;					/* Set to index of white patch for display */
	int si;						/* Sample id index */
	int li;						/* Location id index */
	int ti;						/* Temp index */
	int fi;						/* Colorspace index */
	double plen = 7.366, glen = 2.032, tlen = 18.8;	/* Patch, gap and trailer length in mm */

	set_exe_path(argv[0]);		/* Set global exe_path and error_program */
	check_if_not_interactive();

	if (argc <= 1)
		usage();

	/* Process the arguments */
	mfa = 1;        /* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {		/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1+mfa) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage();

			/* Verbose */
			else if (argv[fa][1] == 'v') {
				verb = 1;
				g_log->verb = verb;

			/* No auto calibration */
			} else if (argv[fa][1] == 'N')
				nocal = 1;

			/* Disable bi-directional strip recognition */
			else if (argv[fa][1] == 'B')
				disbidi = 1;

			/* Force enable bi-directional strip recognition */
			else if (argv[fa][1] == 'b')
				disbidi = 2;

			/* High res mode */
			else if (argv[fa][1] == 'H')
				highres = 1;

			/* Colorimeter Correction Matrix or */
			/* or Colorimeter Calibration Spectral Samples */
			else if (argv[fa][1] == 'X') {
				int ix;
				fa = nfa;
				if (na == NULL) usage();
				strncpy(ccxxname,na,MAXNAMEL-1); ccxxname[MAXNAMEL-1] = '\000';

			/* CCSS Spectral Observer type */
			} else if (argv[fa][1] == 'Q') {
				fa = nfa;
				if (na == NULL) usage();
				if (strcmp(na, "1931_2") == 0) {			/* Classic 2 degree */
					obType = icxOT_CIE_1931_2;
				} else if (strcmp(na, "1964_10") == 0) {	/* Classic 10 degree */
					obType = icxOT_CIE_1964_10;
				} else if (strcmp(na, "2012_2") == 0) {		/* Latest 2 degree */
					obType = icxOT_CIE_2012_2;
				} else if (strcmp(na, "2012_10") == 0) {	/* Latest 10 degree */
					obType = icxOT_CIE_2012_10;
#ifndef SALONEINSTLIB
				} else if (strcmp(na, "1955_2") == 0) {		/* Stiles and Burch 1955 2 degree */
					obType = icxOT_Stiles_Burch_2;
				} else if (strcmp(na, "1978_2") == 0) {		/* Judd and Voss 1978 2 degree */
					obType = icxOT_Judd_Voss_2;
				} else if (strcmp(na, "shaw") == 0) {		/* Shaw and Fairchilds 1997 2 degree */
					obType = icxOT_Shaw_Fairchild_2;
#endif /* !SALONEINSTLIB */
				} else {	/* Assume it's a filename */
					obType = icxOT_custom;
					if (read_cmf(custObserver, na) != 0)
						usage();
				}
			}

			/* Scan tolerance ratio */
			else if (argv[fa][1] == 'T') {
				fa = nfa;
				if (na == NULL)
					usage();
				scan_tol = atof(na);

			/* Suppress warnings */
			} else if (argv[fa][1] == 'S') {
				emit_warnings = 0;

			/* Serial port flow control */
			} else if (argv[fa][1] == 'W') {
				fa = nfa;
				if (na == NULL) usage();
				if (na[0] == 'n' || na[0] == 'N')
					fc = fc_None;
				else if (na[0] == 'h' || na[0] == 'H')
					fc = fc_Hardware;
				else if (na[0] == 'x' || na[0] == 'X')
					fc = fc_XonXOff;
				else
					usage();


			/* Debug coms */
			} else if (argv[fa][1] == 'D') {
				debug = 1;
				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					debug = atoi(na);
					fa = nfa;
				}
				g_log->debug = debug;

			/* COM port  */
			} else if (argv[fa][1] == 'c') {
				fa = nfa;
				if (na == NULL) usage();
				comport = atoi(na);
				if (comport < 1 || comport > 99) usage();

			/* Request transmission measurement */
			} else if (argv[fa][1] == 't') {
				emis = 0;
				trans = 1;
				displ = 0;

			/* Request display measurement */
			} else if (argv[fa][1] == 'd') {

				emis = 0;
				trans = 0;
				displ = 2;

			/* Request emissive measurement */
			} else if (argv[fa][1] == 'e') {
				emis = 1;
				trans = 0;
				displ = 0;

			/* Display type */
			} else if (argv[fa][1] == 'y') {
				fa = nfa;
				if (na == NULL) usage();
				ditype = na[0];
				if (ditype == '_' && na[1] != '\000')
					ditype = ditype << 8 | na[1];


			/* Request patch by patch measurement */
			} else if (argv[fa][1] == 'p') {
				pbypatch = 1;

			/* Request external values */
			} else if (argv[fa][1] == 'x') {
				fa = nfa;
				if (na == NULL) usage();

				if (na[0] == 'l' || na[0] == 'L')
					xtern = 1;
				else if (na[0] == 'x' || na[0] == 'X')
					xtern = 2;
				else
					usage();

			/* Turn off spectral measurement */
			} else if (argv[fa][1] == 'n')
				spectral = 0;

			/* Save as Lab */
			else if (argv[fa][1] == 'l')
				dolab = 1;

			/* Save as Lab */
			else if (argv[fa][1] == 'L')
				dolab = 2;

			/* Resume reading a chart */
			else if (argv[fa][1] == 'r')
				doresume = 1;

			/* Printer calibration info */
			else if (argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage();
				strncpy(calname,na,MAXNAMEL); calname[MAXNAMEL] = '\000';

			/* Filter configuration */
			} else if (argv[fa][1] == 'F') {
				fa = nfa;
				if (na == NULL) usage();
				if (na[0] == 'n' || na[0] == 'N')
					fe = inst_opt_filter_none;
				else if (na[0] == '5')
					fe = inst_opt_filter_D50;
				else if (na[0] == '6')
					fe = inst_opt_filter_D65;
				else if (na[0] == 'u' || na[0] == 'U')
					fe = inst_opt_filter_UVCut;
				else if (na[0] == 'p' || na[0] == 'P')
					fe = inst_opt_filter_pol;
				else
					usage();

			/* XRGA conversion */
			} else if (argv[fa][1] == 'A') {
				fa = nfa;
				if (na == NULL) usage();
				if (na[0] == 'N')
					scalstd = xcalstd_none;
				else if (na[0] == 'A')
					scalstd = xcalstd_xrga;
				else if (na[0] == 'X')
					scalstd = xcalstd_xrdi;
				else if (na[0] == 'G')
					scalstd = xcalstd_gmdi;
				else
					usage();

#ifndef SALONEINSTLIB
			/* Plot spectral patch by patch */
			} else if (argv[fa][1] == 'P') {
				doplot = 1;
#endif

			/* Extra flags */
			} else if (argv[fa][1] == 'Y') {
				if (na == NULL)
					usage();
			
				/* ~~~ i1pro2 test code ~~~ */
				if (na[0] == 'U') {
					uvmode = 1;
				} else {
					usage();
				}
				fa = nfa;

			} else 
				usage();
		} else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(inname,argv[fa]);
	strcat(inname,".ti2");
	strcpy(outname,argv[fa]);
	strcat(outname,".ti3");

	/* See if there is an environment variable ccxx */
	if (ccxxname[0] == '\000') {
		char *na;
		if ((na = getenv("ARGYLL_COLMTER_CAL_SPEC_SET")) != NULL) {
			strncpy(ccxxname,na,MAXNAMEL-1); ccxxname[MAXNAMEL-1] = '\000';

		} else if ((na = getenv("ARGYLL_COLMTER_COR_MATRIX")) != NULL) {
			strncpy(ccxxname,na,MAXNAMEL-1); ccxxname[MAXNAMEL-1] = '\000';
		}
	}

	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI2"); 	/* our special input type is Calibration Target Information 2 */
	icg->add_other(icg, "CAL");		/* There may be a calibration too */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI2 format file");
	if (icg->ntables < 1)
		error ("Input file doesn't contain at least one table");

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data in first table");

	/* Setup output cgats file */
	ocg = new_cgats();				/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI3"); 		/* our special type is Calibration Target Information 3 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll chartread", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	if (displ != 0)
		ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);	/* What sort of device this is */
	else
		ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */

	if (itype == instUnknown) {
		if ((ti = icg->find_kword(icg, 0, "TARGET_INSTRUMENT")) >= 0) {

				if ((itype = inst_enum(icg->t[0].kdata[ti])) == instUnknown)
					error ("Unrecognised chart target instrument '%s'", icg->t[0].kdata[ti]);
		} else {
			itype = instI1Pro;		/* Default chart target instrument */
		}
	}

	if ((ti = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
		ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[ti], NULL);
	
	if ((ti = icg->find_kword(icg, 0, "ACCURATE_EXPECTED_VALUES")) >= 0
	 && strcmp(icg->t[0].kdata[ti], "true") == 0)
		accurate_expd = 1;

	if ((ti = icg->find_kword(icg, 0, "STEPS_IN_PASS")) < 0)
		error ("Input file doesn't contain keyword STEPS_IN_PASS");
	stipa = atoi(icg->t[0].kdata[ti]);

	/* Old style */
	if ((ti = icg->find_kword(icg, 0, "PASSES_IN_STRIPS")) >= 0) {
		char *paists = icg->t[0].kdata[ti];
		int nstr;

		/* Count the number of strips (sheets) */
		for (nstr = 0; paists[nstr] != '\000'; nstr++)
			;

		/* Allocate space for passes in strips */
		if ((pis = (int *)calloc(sizeof(int), nstr+1)) == NULL)
			error("Malloc failed!");

		/* Set the number or passes per strip */
		for (i = 0; i < nstr; i++) {
			pis[i] = b62_int(&paists[i]);
		}
		pis[i] = 0;

	/* New style */
	} else if ((ti = icg->find_kword(icg, 0, "PASSES_IN_STRIPS2")) >= 0) {
		char *cp, *paists = icg->t[0].kdata[ti];
		int nstr;


		/* Count the number of strips (sheets) */
		for (nstr = 1, cp = paists; *cp != '\000'; cp++) {
			if (*cp == ',') {
				nstr++;
				*cp = '\000';
			}
		}


		/* Allocate space for passes in strips */
		if ((pis = (int *)calloc(sizeof(int), nstr+1)) == NULL)
			error("Malloc failed!");

		/* Set the number or passes per strip */
		for (i = 0, cp = paists; i < nstr; i++) {
			pis[i] = atoi(cp);
			cp += strlen(cp) + 1;
		}
		pis[i] = 0;
		
	} else
		error ("Input file doesn't contain keyword PASSES_IN_STRIPS");

	/* Get specified location indexing patterns */
	if ((ti = icg->find_kword(icg, 0, "STRIP_INDEX_PATTERN")) >= 0)
		pixpat = icg->t[0].kdata[ti];
	if ((ti = icg->find_kword(icg, 0, "PATCH_INDEX_PATTERN")) >= 0)
		sixpat = icg->t[0].kdata[ti];
	if ((ti = icg->find_kword(icg, 0, "INDEX_ORDER")) >= 0) {
		if (strcmp(icg->t[0].kdata[ti], "PATCH_THEN_STRIP") == 0)
			ixord = 1;
	}

	if ((ti = icg->find_kword(icg, 0, "RANDOM_START")) >= 0) {
		rstart = atoi(icg->t[0].kdata[ti]);
		rand = 1;
	} else if ((ti = icg->find_kword(icg, 0, "CHART_ID")) >= 0) {
		rstart = atoi(icg->t[0].kdata[ti]);
		rand = 0;
	}

	if ((ti = icg->find_kword(icg, 0, "HEXAGON_PATCHES")) >= 0)
		hex = 1;

	if ((paix = new_alphix(pixpat)) == NULL)
		error("Strip indexing pattern '%s' doesn't parse",pixpat);

	if ((saix = new_alphix(sixpat)) == NULL)
		error("Patch in strip indexing pattern '%s' doesn't parse",sixpat);

	if ((ti = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0)
		ocg->add_kword(ocg, 0, "TOTAL_INK_LIMIT",icg->t[0].kdata[ti], NULL);

	if (itype == instDTP20
	 || itype == instDTP41) {	/* DTP20/41 specific info */
		if ((ti = icg->find_kword(icg, 0, "PATCH_LENGTH")) < 0)
			error ("Input file doesn't contain keyword PATCH_LENGTH");
		plen = atof(icg->t[0].kdata[ti]);
		if ((ti = icg->find_kword(icg, 0, "GAP_LENGTH")) < 0)
			error ("Input file doesn't contain keyword GAP_LENGTH");
		glen = atof(icg->t[0].kdata[ti]);
		if ((ti = icg->find_kword(icg, 0, "TRAILER_LENGTH")) < 0) {
			tlen = 18.0;	/* Default for backwards compatibility */
		} else 
			tlen = atof(icg->t[0].kdata[ti]);
	}

	if (verb) {
		printf("Steps in each Pass = %d\n",stipa);
		printf("Passes in each Strip = ");
		for (i = 0; pis[i] != 0; i++) {
			printf("%s%d",i > 0 ? ", " : "", pis[i]);
		}
		printf("\n");
	}

	/* Fields we want */

	if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
		error ("Input file doesn't contain field SAMPLE_ID");
	if (icg->t[0].ftype[si] != nqcs_t)
		error ("Field SAMPLE_ID is wrong type");
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);

	if ((li = icg->find_field(icg, 0, "SAMPLE_LOC")) < 0)
		error ("Input file doesn't contain field SAMPLE_LOC");
	if (icg->t[0].ftype[li] != cs_t)
		error ("Field SAMPLE_LOC is wrong type");
	ocg->add_field(ocg, 0, "SAMPLE_LOC", cs_t);

	totpa = (npat + stipa -1)/stipa;	/* Total passes for all strips */
	runpat = stipa * totpa;				/* Rounded up totao number of patches */
	if (runpat < npat)		/* (If pattern doesn't match patches) */
		runpat = npat;
	if ((cols = (chcol *)malloc(sizeof(chcol) * runpat)) == NULL)
		error("Malloc failed!");
	if ((scols = (chcol **)calloc(sizeof(chcol *), runpat)) == NULL)
		error("Malloc failed!");

	/* Figure out the color space */
	if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file doesn't contain keyword COLOR_REP");

	if ((nmask = icx_char2inkmask(icg->t[0].kdata[fi])) != 0) {
		int i, j, ii;
		int chix[ICX_MXINKS];	/* Device chanel indexes */
		int xyzix[3];			/* XYZ/Lab chanel indexes */
		char *ident;			/* Full ident */
		char *bident;			/* Base ident */
		char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };
		char *labfname[3] = { "LAB_L", "LAB_A", "LAB_B" };
		int gotexyz = 1;		/* Flag set if input file already has approx XYZ */ 
		icxColorantLu *clu;	/* Xcolorants model based Device -> CIE */

		nchan = icx_noofinks(nmask);
		ident = icx_inkmask2char(nmask, 1); 
		bident = icx_inkmask2char(nmask, 0); 

		/* Device channels */
		for (j = 0; j < nchan; j++) {
			int imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : bident,
			                      icx_ink2char(imask));

			if ((ii = icg->find_field(icg, 0, fname)) < 0)
				error ("Input file doesn't contain field %s",fname);
			if (icg->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type - expect float",fname);
	
			ocg->add_field(ocg, 0, fname, r_t);
			chix[j] = ii;
		}

		/* Approximate XYZ */
		for (j = 0; j < 3; j++) {
			if ((ii = icg->find_field(icg, 0, xyzfname[j])) >= 0) {
				if (icg->t[0].ftype[ii] != r_t)
					error ("Field %s is wrong type - expect float",xyzfname[j]);
				xyzix[j] = ii;
			} else {
				gotexyz = 0;
			}
		}

		/* Measured XYZ and/or Lab */
		if (dolab == 0 || dolab == 2) {
			for (j = 0; j < 3; j++)
				ocg->add_field(ocg, 0, xyzfname[j], r_t);
		} 
		if (dolab == 1 || dolab == 2) {
			for (j = 0; j < 3; j++)
				ocg->add_field(ocg, 0, labfname[j], r_t);
		}

		if ((clu = new_icxColorantLu(nmask)) == NULL)
			error ("Creation of xcolorant lu object failed");

		{
			char fname[100];
			if (dolab)
				sprintf(fname, "%s_LAB", ident);
			else
				sprintf(fname, "%s_XYZ", ident);
			ocg->add_kword(ocg, 0, "COLOR_REP", fname, NULL);
		}

		/* Read all the test patches in */
		for (i = 0; i < npat; i++) {
			cols[i].id = ((char *)icg->t[0].fdata[i][si]);
			cols[i].loc = ((char *)icg->t[0].fdata[i][li]);
			cols[i].n  = nchan;
			for (j = 0; j < nchan; j++)
				cols[i].dev[j] = *((double *)icg->t[0].fdata[i][chix[j]]) / 100.0;
			if (gotexyz) {
				for (j = 0; j < 3; j++)
					cols[i].eXYZ[j] = *((double *)icg->t[0].fdata[i][xyzix[j]]);
			} else {
				clu->dev_to_XYZ(clu, cols[i].eXYZ, cols[i].dev);
				for (j = 0; j < 3; j++)
					cols[i].eXYZ[j] *= 100.0;
			}
			cols[i].XYZ[0] = cols[i].XYZ[1] = cols[i].XYZ[2] = -1.0;
		}

		for (; i < runpat; i++) {
			cols[i].id = cols[i].loc = "-1";
			for (j = 0; j < nchan; j++)
				cols[i].dev[j] = 0.0;
			clu->dev_to_XYZ(clu, cols[i].eXYZ, cols[i].dev);
			for (j = 0; j < 3; j++)
				cols[i].eXYZ[j] *= 100.0;
			cols[i].XYZ[0] = cols[i].XYZ[1] = cols[i].XYZ[2] = -1.0;
		}

		clu->del(clu);
		free(ident);
		free(bident);
	} else
		error ("Input file keyword COLOR_REP has unknown value");

	/* Read any user supplied calibration information */
	if (calname[0] != '\000') {
		if ((cal = new_xcal()) == NULL)
			error("new_xcal failed");
		if ((cal->read(cal, calname)) != 0)
			error("%s",cal->err);
	}

	/* If the user hasn't overridden it, get any calibration in the .ti2 */
	if (cal == NULL) {			/* No user supplied calibration info */
		int oi, tab;

		oi = icg->get_oi(icg, "CAL");

		for (tab = 0; tab < icg->ntables; tab++) {
			if (icg->t[tab].tt == tt_other && icg->t[tab].oi == oi) {
				break;
			}
		}
		if (tab < icg->ntables) {
			if ((cal = new_xcal()) == NULL) {
				error("new_xcal failed");
			}
			if (cal->read_cgats(cal, icg, tab, inname) != 0)  {
				error("%s",cal->err);
			}
		}
	}

	/* If there is calibration information, write it to the .ti3 */
	if (cal != NULL) {			/* No user supplied calibration info */
		if (cal->write_cgats(cal, ocg) != 0) {
			error("%s",cal->err);
		}
		cal->del(cal);
		cal = NULL;
	}

	/* Set up the location sorted array of pointers. */
	/* If the order is not randomized, we don't care what form */
	/* the location identifiers take - i.e. they can be arbitrary. */
	{
		int badloc = 0;

		for (i = 0; i < npat; i++) {
			scols[i] = &cols[i];
			if ((cols[i].loci = patch_location_order(paix, saix, ixord, cols[i].loc)) < 0)
				badloc = 1;
		}
		for (; i < runpat; i++) {	/* Extra on end */
			scols[i] = &cols[i];
			cols[i].loci = (totpa-1) * (256 - stipa) + i; 
/* printf("~~extra = %d, %d\n",cols[i].loci >> 8, cols[i].loci & 255); */
		}

		/* Reset 'read' flag and all data */
		for (i = 0; i < runpat; i++) {
			cols[i].rr = 0;
			cols[i].XYZ[0] = -1.0;
			cols[i].XYZ[1] = -1.0;
			cols[i].XYZ[2] = -1.0;
			cols[i].sp.spec_n = 0;
		}

		if (rand && badloc)
			error ("Bad location field value '%s' on patch %d", cols[i].loc, i);

		if (!badloc) {
#define HEAP_COMPARE(A,B) (A->loci < B->loci)
			HEAPSORT(chcol *, scols, npat);
		}
	}

	/* If we're resuming a chartread, fill in all the patches that */
	/* have been read. */
	if (doresume) {
		cgats *rcg;                 /* output cgats structure */
		int nrpat;					/* Number of resumed patches */
		int lix;					/* Patch location index */
		int islab = 0;				/* nz if Lab, z if XYZ */
		int cieix[3];				/* CIE value indexes */
		int hasspec = 0;			/* nz if has spectral */
		xspect sp;					/* Parameters of spectrum */
		int spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
		char *fname[2][3] = { { "XYZ_X", "XYZ_Y", "XYZ_Z" },
		                      { "LAB_L", "LAB_A", "LAB_B" } };
		int k, ii;
		char buf[200];

		/* Open and look at the .ti3 profile patches file */
		rcg = new_cgats();			/* Create a CGATS structure */
		rcg->add_other(rcg, "CTI3"); 	/* our special input type is Calibration Target Information 3 */
		rcg->add_other(rcg, "CAL"); 	/* our special device Calibration state */
	
		if (rcg->read_name(rcg, outname))
			error("Unable to read chart being resumed '%s' : %s",outname, rcg->err);
	
		if (rcg->ntables == 0 || rcg->t[0].tt != tt_other || rcg->t[0].oi != 0)
			error ("Resumed file '%s' isn't a CTI3 format file",outname);
		if (rcg->ntables < 1)
			error ("Resumed file '%s' doesn't contain at least one table",outname);
	
		if ((lix = rcg->find_field(rcg, 0, "SAMPLE_LOC")) < 0)
			error ("Resumed file '%s' doesn't contain SAMPLE_LOC field",outname);
		if (rcg->t[0].ftype[lix] != cs_t)
			error("Field SAMPLE_LOC is wrong type - corrupted file ?");

		/* Get the CIE field indexes */
		if (rcg->find_field(rcg, 0, "LAB_L") >= 0)
			islab = 1;
			
		for (j = 0; j < 3; j++) {
			if ((cieix[j] = rcg->find_field(rcg, 0, fname[islab][j])) < 0)
				error("Input file doesn't contain field %s",fname[islab][j]);
			if (rcg->t[0].ftype[cieix[j]] != r_t)
				error("Field %s is wrong type - corrupted file ?",fname[islab][j]);
		} 
	
		if ((ii = rcg->find_kword(rcg, 0, "SPECTRAL_BANDS")) >= 0) {
			hasspec = 1;
			sp.spec_n = atoi(rcg->t[0].kdata[ii]);
			if ((ii = rcg->find_kword(rcg, 0, "SPECTRAL_START_NM")) < 0)
				error ("Resumed file '%s' doesn't contain keyword SPECTRAL_START_NM",outname);
			sp.spec_wl_short = atof(rcg->t[0].kdata[ii]);
			if ((ii = rcg->find_kword(rcg, 0, "SPECTRAL_END_NM")) < 0)
				error ("Resumed file '%s' doesn't contain keyword SPECTRAL_END_NM",outname);
			sp.spec_wl_long = atof(rcg->t[0].kdata[ii]);
	
			/* Find the fields for spectral values */
			for (j = 0; j < sp.spec_n; j++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
				            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);
	
				if ((spi[j] = rcg->find_field(rcg, 0, buf)) < 0)
					error("Resumed file '%s' doesn't contain field %s",outname,buf);

				if (rcg->t[0].ftype[spi[j]] != r_t)
					error("Field %s is wrong type - expect float",buf);
			}
		}

		nrpat = rcg->t[0].nsets;

		/* Now see if we can match the previously read patches. */
		/* We'll use the patch location to do this. */
		for (i = 0; i < runpat; i++) {
			int k;
			for (k = 0; k < nrpat; k++) {
				if (strcmp(cols[i].loc, ((char *)rcg->t[0].fdata[k][lix])) == 0)
					break;
			}
			if (k >= nrpat)
				continue;
		
#ifdef DEBUG
			printf("Recovering patch '%s' value from .ti3 file\n",cols[i].loc);
#endif
			cols[i].XYZ[0] = *((double *)rcg->t[0].fdata[k][cieix[0]]);
			cols[i].XYZ[1] = *((double *)rcg->t[0].fdata[k][cieix[1]]);
			cols[i].XYZ[2] = *((double *)rcg->t[0].fdata[k][cieix[2]]);
			if (islab) {
				icmLab2XYZ(&icmD50, cols[i].XYZ, cols[i].XYZ);
				cols[i].XYZ[0] *= 100.0;
				cols[i].XYZ[1] *= 100.0;
				cols[i].XYZ[2] *= 100.0;
			}
			if (hasspec) {
				cols[i].sp.spec_n = sp.spec_n;
				cols[i].sp.spec_wl_short = sp.spec_wl_short;
				cols[i].sp.spec_wl_long = sp.spec_wl_long;
				for (j = 0; j < sp.spec_n; j++)
					cols[i].sp.spec[j] = *((double *)rcg->t[0].fdata[k][spi[j]]);
			}
			cols[i].rr = 1;
		}
		rcg->del(rcg);
	}

	/* We can't fiddle white point with spectral data, */
	/* so turn spectral off for display with white point relative. */
	if (displ == 2 || displ == 3) {
		spectral = 0;

		/* Check that there is a white patch, so that we can compute Y relative */
		/* Read all the test patches in */
		if (nmask != ICX_RGB)
			error("Don't know how to handle non-RGB display space");

		for (wpat = 0; wpat < npat; wpat++) {
			if (cols[wpat].dev[0] > 0.9999999 && 
			    cols[wpat].dev[1] > 0.9999999 && 
			    cols[wpat].dev[2] > 0.9999999) {
				break;
			}
		}
		if (wpat >= npat) {	/* Create a white patch */
			error("Can't compute white Y relative display values without a white test patch");
		}
	}

	if (!xtern) {
		if ((icmps = new_icompaths(g_log)) == NULL)
			error("Finding instrument paths failed");
		if ((ipath = icmps->get_path(icmps, comport)) == NULL)
			error("No instrument at port %d",comport);
	}

	/* Read all of the strips in */
	if (read_strips(itype, scols, &atype, npat, totpa, stipa, pis, paix,
	                saix, ixord, rstart, rand, hex, ipath, fc, plen, glen, tlen,
	                trans, emis, displ, ditype, fe, scalstd, &ucalstd, nocal, disbidi, highres,
		            ccxxname, obType, custObserver,
	                scan_tol, pbypatch, xtern, spectral, uvmode, accurate_expd,
	                emit_warnings, doplot, g_log) == 0) {
		/* And save the result */

		int nrpat;				/* Number of read patches */
		int vpix = 0;			/* Valid patch index, if nrpatch > 0 */
		int nsetel = 0;
		cgats_set_elem *setel;	/* Array of set value elements */

		/* Note what instrument the chart was read with */
		ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", inst_name(atype) , NULL);

		/* X-Rite calibration standard (If reflective mode) */
		if (displ == 0 && trans == 0 && ucalstd != xcalstd_none) 
			ocg->add_kword(ocg, 0, "DEVCALSTD",xcalstd2str(ucalstd), NULL);

		if (fe == inst_opt_filter_pol)
			ocg->add_kword(ocg, 0, "INSTRUMENT_FILTER", "POLARIZED", NULL);
		else if (fe == inst_opt_filter_D65)
			ocg->add_kword(ocg, 0, "INSTRUMENT_FILTER", "D65", NULL);
		else if (fe == inst_opt_filter_UVCut)
			ocg->add_kword(ocg, 0, "INSTRUMENT_FILTER", "UVCUT", NULL);

		/* Count patches actually read */
		for (nrpat = i = 0; i < npat; i++) {
			if (cols[i].rr) {
				vpix = i;
				nrpat++;
			}
		}

		/* If we've used a display white relative mode, record the absolute white */
		if (displ == 2 || displ == 3) {
			double nn[3];
			char buf[100];

			if (cols[wpat].rr == 0) {
				error("Can't compute white Y relative display values without reading a white test patch");
			}
			sprintf(buf,"%f %f %f", cols[wpat].XYZ[0], cols[wpat].XYZ[1], cols[wpat].XYZ[2]);
			ocg->add_kword(ocg, 0, "LUMINANCE_XYZ_CDM2",buf, NULL);

			/* Normalise to white Y 100 */
			if (displ == 2) {
				nn[0] = 100.0 / cols[wpat].XYZ[1];
				nn[1] = 100.0 / cols[wpat].XYZ[1];
				nn[2] = 100.0 / cols[wpat].XYZ[1];
			/* Normalise to the white point */
			} else {
				nn[0] = 100.0 * icmD50.X / cols[wpat].XYZ[0];
				nn[1] = 100.0 * icmD50.Y / cols[wpat].XYZ[1];
				nn[2] = 100.0 * icmD50.Z / cols[wpat].XYZ[2];
			}

			for (i = 0; i < npat; i++) {
				if (cols[i].rr) {
					cols[i].XYZ[0] *= nn[0];
					cols[i].XYZ[1] *= nn[1];
					cols[i].XYZ[2] *= nn[2];
				}
			}
		}

		nsetel += 1;		/* For id */
		nsetel += 1;		/* For loc */
		nsetel += nchan;	/* For device values */
		nsetel += 3;		/* For XYZ or Lab */
		if (dolab == 2)
			nsetel += 3;	/* For XYZ and Lab */

		/* If we have spectral information, output it too */
		if (nrpat > 0 && cols[vpix].sp.spec_n > 0) {
			char buf[100];

			nsetel += cols[vpix].sp.spec_n;		/* Spectral values */
			sprintf(buf,"%d", cols[vpix].sp.spec_n);
			ocg->add_kword(ocg, 0, "SPECTRAL_BANDS",buf, NULL);
			sprintf(buf,"%f", cols[vpix].sp.spec_wl_short);
			ocg->add_kword(ocg, 0, "SPECTRAL_START_NM",buf, NULL);
			sprintf(buf,"%f", cols[vpix].sp.spec_wl_long);
			ocg->add_kword(ocg, 0, "SPECTRAL_END_NM",buf, NULL);

			/* Generate fields for spectral values */
			for (i = 0; i < cols[vpix].sp.spec_n; i++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(cols[vpix].sp.spec_wl_short + ((double)i/(cols[vpix].sp.spec_n-1.0))
				            * (cols[vpix].sp.spec_wl_long - cols[vpix].sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);
				ocg->add_field(ocg, 0, buf, r_t);
			}
		}

		if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * nsetel)) == NULL)
			error("Malloc failed!");

		/* Write out the patch info to the output CGATS file */
		for (i = 0; i < npat; i++) {
			int k = 0;

			if (cols[i].rr == 0					/* If this patch wasn't read */
			 || strcmp(cols[i].id, "0") == 0)	/* or it is a padding patch. */
				continue;			/* Skip it */

			setel[k++].c = cols[i].id;
			setel[k++].c = cols[i].loc;
			
			for (j = 0; j < nchan; j++)
				setel[k++].d = 100.0 * cols[i].dev[j];

			if (dolab == 0 || dolab == 2) {
				setel[k++].d = cols[i].XYZ[0];
				setel[k++].d = cols[i].XYZ[1];
				setel[k++].d = cols[i].XYZ[2];
			}
			if (dolab == 1 || dolab == 2) {
				double lab[3];
				double xyz[3];

				xyz[0] = cols[i].XYZ[0]/100.0;
				xyz[1] = cols[i].XYZ[1]/100.0;
				xyz[2] = cols[i].XYZ[2]/100.0;
				icmXYZ2Lab(&icmD50, lab, xyz);
				setel[k++].d = lab[0];
				setel[k++].d = lab[1];
				setel[k++].d = lab[2];
			}

			/* Check that the spectral matches, in case we're resuming */
			if (   cols[i].sp.spec_n != cols[vpix].sp.spec_n
				|| fabs(cols[i].sp.spec_wl_short - cols[vpix].sp.spec_wl_short) > 0.01
				|| fabs(cols[i].sp.spec_wl_long - cols[vpix].sp.spec_wl_long) > 0.01 ) {
				error("The resumed spectral type seems to have changed!");
			}

			for (j = 0; j < cols[i].sp.spec_n; j++) {
				setel[k++].d = cols[i].sp.spec[j];
			}

			ocg->add_setarr(ocg, 0, setel);
		}

		free(setel);

		if (ocg->write_name(ocg, outname))
			error("Write error : %s",ocg->err);
	}

	free(scols);
	icmps->del(icmps);
	free(pis);
	saix->del(saix);
	paix->del(paix);
	free(cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */

	return 0;
}


