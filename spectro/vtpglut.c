
/* Abstract Video Test Patch Generator and/or 3dLUT class implemenation */

/* 
 * Argyll Color Correction System
 *
 * Author: Graeme W. Gill
 * Date:   10/3/2001
 *
 * Copyright 2001 - 2013 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 2 or later :-
 * see the License2.txt file for licencing details.
 */

/*
	TTBD:

	Note that there is some support for HW 3dLUTs in X11 xrandr
	in a patch for xf86-videeo-amdgpu drmmode_display.c
	submitted on 3 May 2018. 
	See <https://lists.freedesktop.org/archives/amd-gfx/2018-May/022007.html>

*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#ifndef SALONEINSTLIB
#include "copyright.h"
#include "aconfig.h"
#else
#include "sa_config.h"
#endif /* !SALONEINSTLIB */
#include "numsup.h"
#include "cgats.h"
#include "xspect.h"
#include "conv.h"

#include "insttypes.h"
#include "icoms.h"
#include "vtpglut.h"
#include "rspec.h"
#include "vtpgluttypes.h"
#include "sort.h"

icom_type dev_category(instType itype);

/* --------------------------------------------------------- */

#if defined(ENABLE_FAST_SERIAL)
extern devType fast_ser_dev_type(icoms *p, int tryhard, 
       inst_code (*uicallback)(void *cntx, inst_ui_purp purp), void *cntx);
# if defined(ENABLE_SERIAL)
static devType ser_vtpglut_type(icoms *p, 
       vtpglut_code (*uicallback)(void *cntx, vtpglut_ui_purp purp), void *cntx);

/* Translate vtpglut uicallback into inst callback */
typedef struct {
    void *cntx;
    vtpglut_code (*uicallback)(void *cntx, vtpglut_ui_purp purp);
} vtpglut_callback_cntx;

static inst_code ser_uicallback(void *cntx, inst_ui_purp purp) {
	vtpglut_callback_cntx *p = (vtpglut_callback_cntx *)cntx;
	vtpglut_ui_purp ser_purp;

	if (purp == vtpglut_negcoms)
		ser_purp = inst_negcoms;
	else
		return inst_internal_error; 
	return p->uicallback(p->cntx, ser_purp);
}

# endif /* ENABLE_SERIAL */
#endif /* ENABLE_FAST_SERIAL */


/* ------------------------------------ */
/* Default methods for instrument class */

/* Establish communications at the indicated baud rate. */
/* Timout in to seconds, and return non-zero error code */
static vtpglut_code init_coms(
vtpglut *p) {		/* Timeout */
	return vtpglut_unsupported;
}

/* Initialise or re-initialise the INST */
/* return non-zero on an error, with inst error code */
static vtpglut_code init_vtpglut(  
vtpglut *p) {
	return vtpglut_unsupported;
}

/* Return the device type */
static devType get_dtype(vtpglut *p) {
	if (p != NULL)
		return p->dtype;
	return devUnknown;
}

/* Return the device serial number. */
/* (This will be an empty string if there is no serial no) */
static char *get_serial_no(vtpglut *p) {
	return "";
}

/* Return the device mode/capabilities */
static void capabilities(vtpglut *p, vtpglut_mode *cap1,
                         vtpglut_capability *cap2) {
	if (cap1 != NULL)
		*cap1 = vtpglut_mode_none;
	if (cap2 != NULL)
		*cap2 = vtpglut2_none;
}

/* Check the device mode */                                       
static vtpglut_code check_mode(
vtpglut *p,
vtpglut_mode m) {		/* Requested mode */
	return vtpglut_unsupported;
}

/* Set the device mode */                                       
static vtpglut_code set_mode(
vtpglut *p,
vtpglut_mode m) {		/* Requested mode */
	return vtpglut_unsupported;
}

/* Get a status or set or get an option (default implementation) */
vtpglut_code vtpglut_get_set_opt_def(
vtpglut *p,
vtpglut_opt_type m,	/* Requested option type */
va_list args) {		/* Option parameters */                             
	return vtpglut_unsupported;
}

/* Get a status or set or get an option */
static vtpglut_code get_set_opt(
vtpglut *p,
vtpglut_opt_type m,	/* Requested status type */
...) {				/* Status parameters */                             
	vtpglut_code rv;
	va_list args;

	va_start(args, m);
	rv = vtpglut_get_set_opt_def(p, m, args);	/* Call the above */
	va_end(args);

	return rv;
}

/* Supply a user interaction callback function. */
static void set_uicallback(
vtpglut *pp,
vtpglut_code (*uicallback)(void *cntx, vtpglut_ui_purp purp),
void *cntx)	 {
	pp->uicallback = uicallback;
	pp->uic_cntx = cntx;
}

/* Supply an asynchronous event callback function. */
static void set_event_callback(
vtpglut *pp,
void (*eventcallback)(void *cntx, vtpglut_event_type event),
void *cntx)	 {
	pp->eventcallback = eventcallback;
	pp->event_cntx = cntx;
}

/* Generic inst error codes interpretation */
static char *vtpglut_interp_error(vtpglut *p, vtpglut_code ec) {
	switch(ec & vtpglut_mask) {
		case vtpglut_ok:
			return "No error";
		case vtpglut_notify:
			return "Notification";
		case vtpglut_warning:
			return "Warning";
		case vtpglut_no_coms:
			return "Internal error - communications needed but not established";
		case vtpglut_no_init:
			return "Internal error - initialisation needed but not done";
		case vtpglut_unsupported:
			return "Unsupported function";
		case vtpglut_internal_error:
			return "Internal software error";
		case vtpglut_coms_fail:
			return "Communications failure";
		case vtpglut_unknown_model:
			return "Not expected device model";
		case vtpglut_protocol_error:
			return "Communication protocol breakdown";
		case vtpglut_user_abort:
			return "User hit Abort Key";
		case vtpglut_user_trig:
			return "User hit Trigger Key";
		case vtpglut_unexpected_reply:
			return "Unexpected Reply";
		case vtpglut_wrong_setup:
			return "Wrong or conflicting setup";
		case vtpglut_hardware_fail:
			return "Hardware Failure";
		case vtpglut_system_error:
			return "Operating System Error";
		case vtpglut_bad_parameter:
			return "Bad Parameter Value";
		case vtpglut_other_error:
			return "Non-specific error";
	}
	return "Unknown inst error code";
}

/* Instrument specific error codes interpretation */
static char *interp_error(vtpglut *p, int ec) {
	return "interp_error is uninplemented in base class!";
}

/* Return the last serial communication error code */
/* (This is used for deciding fallback/retry strategies) */
static int last_scomerr(vtpglut *p) {
	return p->icom->lserr;
}

/* ---------------------------------------------- */

/* Delete things set/done by new_vtpglut() */
static vtpglut_code virtual_del(vtpglut *p) {

#if defined(UNIX_APPLE)
	osx_latencycritical_end();
#endif

	return vtpglut_ok;
}


/* Virtual constructor. */
/* Return NULL for unknown instrument, */
/* or serial instrument if nocoms == 0. */
extern vtpglut *new_vtpglut(
icompath *path,		/* Device path to instrument */
int nocoms,			/* Don't open if communications are needed to establish inst type */
a1log *log,			/* log to use */
vtpglut_code (*uicallback)(void *cntx, vtpglut_ui_purp purp),	/* optional uicallback for abort */
void *cntx			/* Context for callback */
) {
	devType itype = devUnknown;		/* Actual type */
	icoms *icom;
	vtpglut *p = NULL;

	if (path == NULL) {
		a1logd(log, 2, "new_vtpglut: got NULL path\n");
		return NULL;
	}

	a1logd(log, 2, "new_vtpglut: called with path '%s' type '%s'\n",path->name,inst_sname(path->itype));

	if ((icom = new_icoms(path, log)) == NULL) {
		a1logd(log, 2, "new_vtpglut: new_icoms failed to open it\n");
		return NULL;
	}



	/* Set device type from USB port, if not specified */
	itype = icom->itype;		/* Instrument type if its known from usb/hid */

#if defined(ENABLE_FAST_SERIAL)
	if (itype == instUnknown && !nocoms && (icom->dctype & icomt_fastserial)) {
		vtpglut_callback_cntx wcntx;

		/* Wrap vtpglut callback in inst callback */
		wcntx.cntx = cntx;
		wcntx.uicallback = uicallback;

		itype = fast_ser_dev_type(icom, 1, ser_uicallback, &wcntx);		/* Else type from serial */
		icom->dctype = (icom->dctype & ~icomt_cat_mask) | dev_category(itype);
		a1logd(log, 8, "new_vtpglut: fast set '%s' dctype 0x%x\n",icom->name,icom->dctype);
	}
#endif /* ENABLE_FAST_SERIAL */

#if defined(ENABLE_SERIAL)
	if (itype == instUnknown && !nocoms) {
		itype = ser_vtpglut_type(icom, uicallback, cntx);		/* Else type from serial */
		icom->dctype = (icom->dctype & ~icomt_cat_mask) | dev_category(itype);
		a1logd(log, 8, "new_vtpglut: set '%s' dctype 0x%x\n",icom->name,icom->dctype);
	}
#endif /* ENABLE_SERIAL */


#ifdef ENABLE_SERIAL
//	if (itype == devRadiance)
//		p = (vtpglut *)new_radiance(icom, itype);
#endif /* ENABLE_SERIAL */

#ifdef ENABLE_USB
//	if (itype == instXXXX)
//		p = (vtpglut *)new_XXXX(icom, itype);
#endif /* ENABLE_USB */


	/* Nothing matched */
	if (p == NULL) {
		a1logd(log, 2, "new_vtpglut: device type not recognised\n");
		icom->del(icom);
		return NULL;
	}

	p->vdel = virtual_del;

	/* Add default methods if constructor did not supply them */
	if (p->init_coms == NULL)
		p->init_coms = init_coms;
	if (p->init_vtpglut == NULL)
		p->init_vtpglut = init_vtpglut;
	if (p->get_dtype == NULL)
		p->get_dtype = get_dtype;
	if (p->get_serial_no == NULL)
		p->get_serial_no = get_serial_no;
	if (p->capabilities == NULL)
		p->capabilities = capabilities;
	if (p->check_mode == NULL)
		p->check_mode = check_mode;
	if (p->set_mode == NULL)
		p->set_mode = set_mode;
	if (p->get_set_opt == NULL)
		p->get_set_opt = get_set_opt;
	if (p->set_uicallback == NULL)
		p->set_uicallback = set_uicallback;
	if (p->set_event_callback == NULL)
		p->set_event_callback = set_event_callback;
	if (p->vtpglut_interp_error == NULL)
		p->vtpglut_interp_error = vtpglut_interp_error;
	if (p->interp_error == NULL)
		p->interp_error = interp_error;
	if (p->last_scomerr == NULL)
		p->last_scomerr = last_scomerr;

	/* Set the provided user interaction callback */
	p->set_uicallback(p, uicallback, cntx);

#if defined(UNIX_APPLE)
	osx_latencycritical_start();
#endif

	return p;
}

/* ============================================================= */
/* Detect serial device */

#ifdef ENABLE_SERIAL
static void hex2bin(char *buf, int len);

/* Heuristicly determine the device type for */
/* a serial connection, and devUnknown if not serial. */
/* Set it in icoms and also return it. */
static devType ser_vtpglut_type(
	icoms *p,
	vtpglut_code (*uicallback)(void *cntx, vtpglut_ui_purp purp),		/* optional uicallback */
	void *cntx			/* Context for callback */
) {
	devType rv = instUnknown;
#define BUFSZ (128 + 10)
	char buf[BUFSZ];
	baud_rate brt[] = { baud_9600, baud_57600, baud_115200, baud_230400, baud_nc };
	unsigned int etime;
	unsigned int bi, i;
	int se, len, bread;
	int lumagen = 0, lgchsum = 0;
	
#ifdef ENABLE_USB
	if (p->usbd != NULL || p->hidd != NULL)
		return p->itype;
#endif /* ENABLE_USB */

	bi = 0;

	/* The tick to give up on */
	etime = msec_time() + (long)(20.0 * 1000.0 + 0.5);

	a1logd(p->log, 1, "ser_vtpglut_type: Trying different baud rates (%u msec to go)\n",etime - msec_time());

	/* Until we time out, find the correct baud rate */
	for (i = bi; msec_time() < etime; i++) {
		lgchsum = 0;

		if (brt[i] == baud_nc)
			i = 0;
		if ((se = p->set_ser_port(p, fc_None, brt[i], parity_none,
			                         stop_1, length_8)) != ICOM_OK) { 
			a1logd(p->log, 5, "ser_vtpglut_type: set_ser_port failed with 0x%x\n",se);
			return instUnknown;		/* Give up */
		}

		a1logd(p->log, 5, "Trying %s baud\n",baud_rate_to_str(brt[i]));
		bread = 0;

		/* If Lumagen baud rate */
		if (brt[i] == baud_9600
		 || brt[i] == baud_57600
		 || brt[i] == baud_115200
		 || brt[i] == baud_230400) {

			/* Try a spectrolino/spectroscan command first, so as not to upset it, */
			/* but do first character only and see if there is an echo */
			p->write_read_ex(p, ";", 1, buf, BUFSZ-1, &bread, "\r", 1, 0.2, 1);
			if (bread == 1 && buf[0] == ';')
				goto check_lumagen;		/* It may be a Lumagen, so skip it. */
	
			/* Send the rest of the spectrolino command */
			p->write_read_ex(p, "D024\r\n", 0, buf, BUFSZ-1, &bread, "\r", 1, 0.5, 1);


			if (bread == 0) {
				char *bp;
				a1logd(p->log, 5, "ser_vtpglut_type: Spectrolino command returned nothing\n");

				/* It could be a Lumagen Radiance with echo off, */
				/* so poke it and see if it responds. */
				/* (Unfortunately the Lumagen delimeters modes aren't */
				/*  backwards compatible, so we may have to poke it twice...) */

				/* Send "X" first, to get it out of menu mode ? */
				p->write_read_ex(p, "X", 1, buf, BUFSZ, NULL, "\n", 1, 0.1, 1);		// Menu off

				a1logd(p->log, 5, "ser_vtpglut_type: Checking for Lumagen Radiance\n");
				p->write_read_ex(p, "#ZQS00\r", 0, buf, BUFSZ, &bread, "\n", 1, 0.5, 1);
				if (bread == 0) {
					p->write_read_ex(p, "#0ZQS008E\r", 0, buf, BUFSZ, &bread, "\n", 1, 0.5, 1);
					lgchsum = 1;
				}

				if (bread > 0
				 && (bp = strrchr(buf, '!')) != NULL
				 && strlen(bp) >= 4
				 && strncmp(bp,"!S00",4) == 0) {
					goto check_lumagen;
				}

				/* Nope - look for something at a different baud rate */
				goto continue_looking;
			}
			buf[bread] = '\000';
			len = strlen(buf);

			/* Is this a Lumagen Radiance with echo on, it responds with ";D024..." */
			if ((len >= 4 && strncmp(buf, ";D024", 4) == 0)
			 || (len >= 4 && strncmp(buf, "!N\n\r", 4) == 0)) {
				char *bp;

			  check_lumagen:;

				/* Get the Lumagen device information */
				p->write_read_ex(p, lgchsum ? "#0ZQS018F\r" : "#ZQS01\r",
					                    0, buf, BUFSZ, &bread, "\n", 1, 0.5, 1);

				/* Might have echo with checksum, so lgchsum not set correctly */
				if (!lgchsum && bread > 0 && strstr(buf, "!N") != NULL) {
					p->write_read_ex(p, "#0ZQS018F\r", 0, buf, BUFSZ, &bread, "\n", 1, 0.5, 1);
					if (bread >= 11 && strncmp(buf, "#0ZQS018F!Y", 11) == 0)
						lgchsum = 1;
				}
			
				/* returns something like "ZQS01!S01,Radiance2020,030115,1016,001309\r\m" */
				if ((bp = strrchr(buf, '!')) != NULL && strlen(bp) >= 13) {
				    if (strncmp(bp,"!S01,Radiance",13) == 0) {
						a1logd(p->log, 5, "ser_vtpglut_type: Found Lumagen Radiance\n");
						lumagen = 1;
						break;
					}
				}
				a1logd(p->log, 5, "ser_vtpglut_type: Not Lumagen Radiance\n");
				goto continue_looking;
			}
		}	/* Possibly Lumagen */

	  continue_looking:;

		/* Check for user abort */
		if (uicallback != NULL) {
			vtpglut_code ev;
			if ((ev = uicallback(cntx, vtpglut_negcoms)) == vtpglut_user_abort) {
				a1logd(p->log, 5, "ser_vtpglut_type: User aborted\n");
				return instUnknown;
			}
		}
	}	/* next baud */

	if (rv == instUnknown
	 && lumagen == 0
	 && msec_time() >= etime) {		/* We haven't established comms */
		a1logd(p->log, 5, "ser_vtpglut_type: Failed to establish coms\n");
		return instUnknown;
	}

	a1logd(p->log, 5, "ser_vtpglut_type: Got coms with device\n");

	if (lumagen) {
		char *bp;

		/* Get the Lumagen device information */
		if ((se = p->write_read_ex(p, "ZQS01", 0, buf, BUFSZ, NULL, "\n", 1, 2.5, 1)) != 0)
			return instUnknown;
	
		/* returns something like "ZQS01!S01,Radiance2020,030115,1016,001309\r\m" */
		if ((bp = strstr(buf, "!")) != NULL && strlen(bp) >= 13) {
			if (strncmp(bp,"!S01,Radiance",13) == 0) {
				rv = devRadiance;
			}
		}
	}

	a1logd(p->log, 5, "ser_vtpglut_type: Device type is '%s'\n", inst_name(rv));

	p->close_port(p);	/* Or should we leave it open ?? */

	p->itype = rv;

	return rv;
}
#undef BUFSZ

#endif /* ENABLE_SERIAL */

#if defined(ENABLE_SERIAL) || defined(ENABLE_FAST_SERIAL)

/* Convert an ASCII Hex character to an integer. */
static int h2b(char c) {
	if (c >= '0' && c <= '9')
		return (c-(int)'0');
	if (c >= 'A' && c <= 'F')
		return (10 + c-(int)'A');
	if (c >= 'a' && c <= 'f')
		return (10 + c-(int)'a');
	return 0;
}

/* Convert a Hex encoded buffer into binary. */
/* len is number of bytes out */
static void hex2bin(char *buf, int len) {
	int i;

	for (i = 0; i < len; i++) {
		buf[i] = (char)((h2b(buf[2 * i + 0]) << 4)
		              | (h2b(buf[2 * i + 1]) << 0));
	}
}

#endif /* ENABLE_SERIAL */























