
#ifndef VTPGLUT_H

/*
 * v3dlut/vtpg API definition.
 *
 * Abstract base class for common color Video Test Pattern Generator interface
 * and Video 3DLut devices.
 */

/* 
 * Argyll Color Management System
 *
 * Author: Graeme W. Gill
 * Date:   15/3/2001
 *
 * Copyright 2001 - 2013 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 2 or later :-
 * see the License2.txt file for licencing details.
 *
 */

#include "dev.h"			/* Device base class */
#include "insttypes.h"		/* libinst Includes this functionality */
#include "icoms.h"			/* libinst Includes this functionality */
#include "conv.h"

#ifdef __cplusplus
	extern "C" {
#endif

/* ------------------------------------------------- */
/* aprox. debug level guide:

	1,2  Applications, internal errors
	2,3  High level instrument drivers
	4,5  High level instrument communications
	6,7  High level serial/USB communications
	8,9  Low level serial/USB communications

*/

/* ---------------------------------------- */
/* Device interface abstract base class. */

/* This is used for all device types except color instruments */

/* Abstract return codes in ms 8bits. */
/* Device dependant codes in ls 16 bits. */
/* Note :- update vtpglut_interp_error() in other.c if anything here is changed. */
/* and also check all the device specific XXX_interp_code() routines too. */
typedef enum {
	vtpglut_ok                = 0x000000,
	vtpglut_notify            = 0x010000,	/* A Notification */
	vtpglut_warning           = 0x020000,	/* A Warning */
	vtpglut_no_coms           = 0x030000,	/* init_coms() hasn't been called yet */
	vtpglut_no_init           = 0x040000,	/* init_dev() hasn't been called yet */
	vtpglut_unsupported       = 0x050000,	/* Unsupported function */
	vtpglut_internal_error    = 0x060000,	/* Internal software error */
	vtpglut_coms_fail         = 0x070000,	/* Communication failure */
	vtpglut_unknown_model     = 0x080000,	/* Not the expected device */
	vtpglut_protocol_error    = 0x090000, 	/* Read or Write protocol error */
	vtpglut_user_abort        = 0x0A0000,	/* User hit escape */
	vtpglut_user_trig         = 0x0C0000,	/* User hit trigger key */
	vtpglut_unexpected_reply  = 0x140000,	/* Unexpected Reply */
	vtpglut_wrong_setup       = 0x150000,	/* Setup is wrong or conflicting */
	vtpglut_hardware_fail     = 0x160000,	/* Hardware failure */
	vtpglut_system_error      = 0x170000,	/* System call (ie malloc) fail */
	vtpglut_bad_parameter     = 0x180000,	/* Bad parameter value */
	vtpglut_other_error       = 0x190000,	/* Some other error */
	vtpglut_mask              = 0xff0000,	/* vtpglut_code mask value */
	vtpglut_dmask             = 0x00ffff	/* device specific mask value */
} vtpglut_code;

/* Device capabilities & modes */
/* Note that due to the binary combinations, capabilities is not definititive */
/* as to valid modes. check_mode() is definitive. */
/* #defines are for saving modes in a version independent way. */
/* Note :- update vtpglut_mode_sym[] table in vtpglut.c if anything here is changed. */
typedef enum {
	vtpglut_mode_none                = 0x00000000, /* No capability/mode */
} vtpglut_mode;

typedef enum {
	vtpglut2_none                    = 0x00000000, /* No capability */
} vtpglut_capability;


/* Device options for get_set_opt() */
typedef enum {
	vtpglut_opt_unknown            = 0x0000,	/* Option not specified */
} vtpglut_opt_type;

/* User interaction callback (uicallback()) function purpose */
/* (Make sure vtpglut.c inst_code ser_uicallbac() is updated if this changes. */
typedef enum {
    vtpglut_negcoms			/* Negotiating communications - can abort */
} vtpglut_ui_purp;

/* Asynchronous event callback type */
typedef enum {
    inst_event_none		/* (Placeholder) */
} vtpglut_event_type;

/* ---------------------------------------- */
/* Device interface abstract base class */

# define EXTRA_VTPGLUT_OBJ

/* vtpg/v3dlut interface base object */
/* Note that some methods work after creation, while many */
/* will return an error if communications hasn't been established and */
/* the device initialised. Some may change their response before and */
/* after initialisation. */
#define VTPGLUT_OBJ_BASE														\
																				\
	DEV_OBJ_BASE																\
																				\
	EXTRA_VTPGLUT_OBJ															\
																				\
	vtpglut_code (*uicallback)(void *cntx, vtpglut_ui_purp purp);				\
	void *uic_cntx;		/* User interaction callback function */				\
	void (*eventcallback)(void *cntx, vtpglut_event_type event);				\
	void *event_cntx;	/* Asynchronous event callback function */				\
																				\
	/* Virtual delete. Cleans up things done by new_vtpglut(). */				\
	vtpglut_code (*vdel)(														\
        struct _vtpglut *p);													\
																				\
	/* Establish communications at the indicated baud rate. */					\
	/* (Serial parameters are ignored for USB/Network device) */				\
	/* Timout in to seconds, and return non-zero error code */					\
	vtpglut_code (*init_coms)(													\
        struct _vtpglut *p);													\
																				\
	/* Initialise or re-initialise the vtpglut. */								\
	/* return non-zero on an error, with dev error code */						\
	vtpglut_code (*init_vtpglut)(  												\
        struct _vtpglut *p);													\
																				\
	/* Return the device type */												\
	/* (this could concievably change after init_vtpglut()) */					\
	/* Can be called before init */												\
	devType (*get_dtype)(  														\
        struct _vtpglut *p);													\
																				\
	/* Return the device serial number. */										\
	/* (This will be an empty string if there is no serial no) */               \
	char *(*get_serial_no)(  													\
        struct _vtpglut *p);													\
																				\
	/* Return the avilable devices modes and capabilities. */					\
	/* Can be called before init, but may be different to */					\
	/* what's returned after initilisation. */									\
	/* Note that these may change with the mode. */								\
	/* Arguments may be NULL */													\
	void (*capabilities)(struct _vtpglut *p,									\
	        vtpglut_mode *cap1,													\
	        vtpglut_capability *cap2);											\
																				\
    /* Check that the particular device measurement mode is valid, */           \
	/* since it's not possible to be 100% sure from capabilities */				\
    vtpglut_code (*check_mode)(													\
        struct _vtpglut *p,														\
        vtpglut_mode m);		/* Requested mode */							\
																				\
    /* Set the device mode */													\
	/* Note that this may change the capabilities. */							\
    vtpglut_code (*set_mode)(													\
        struct _vtpglut *p,														\
        vtpglut_mode m);		/* Requested mode */							\
																				\
    /* Get a status or get or set an option */                                  \
	/* option state. */															\
	/* Some options can be set before init */									\
	/* See vtpglut_opt_type typedef for list of mode types */					\
    vtpglut_code (*get_set_opt)(												\
        struct _vtpglut *p,														\
        vtpglut_opt_type m,	/* Requested option mode */							\
		...);				/* Option parameters */                             \
																				\
	/* Supply a user interaction callback function.								\
	 * (Nothing currentlt defined) 												\
	 *																			\
	 * NULL can be set to disable the callback.									\
	 */																			\
	void (*set_uicallback)(struct _vtpglut *p,										\
		vtpglut_code (*uicallback)(void *cntx, vtpglut_ui_purp purp), 			\
		void *cntx);															\
																				\
	/* Supply an aynchronous event callback function.							\
	 * This is called from a different thread with the following possible events:	\
	 *																			\
     * inst_event_XXX: (placeholder)											\
	 *																			\
	 * NULL can be set to disable the callback.									\
	 */																			\
	void (*set_event_callback)(struct _vtpglut *p,									\
		void (*eventcallback)(void *cntx, vtpglut_event_type event),			\
		void *cntx);															\
																				\
	/* Generic device error codes interpretation */								\
	char * (*vtpglut_interp_error)(struct _vtpglut *p, vtpglut_code ec);		\
																				\
	/* Instrument specific error codes interpretation */						\
	char * (*interp_error)(struct _vtpglut *p, int ec);							\
																				\
	/* Return the last serial communication error code */						\
	/* (This is used for deciding fallback/retry strategies) */					\
	int (*last_scomerr)(struct _vtpglut *p);									\
																				\
	/* Destroy ourselves */														\
	void (*del)(struct _vtpglut *p);											\


/* The base object type */
struct _vtpglut {
	VTPGLUT_OBJ_BASE
	}; typedef struct _vtpglut vtpglut;

/* Virtual constructor. */
/* Return NULL for unknown device, */
/* or serial device if nocoms == 0. */
/* (Doesn't copy icompaths log!) */
extern vtpglut *new_vtpglut(
	icompath *path,		/* Device path for this device */
	int nocoms,			/* Don't open if communications are needed to establish device type */
	a1log *log,			/* Log to use */
	vtpglut_code (*uicallback)(void *cntx, vtpglut_ui_purp purp),
	                    /* optional uicallback for abort */
	void *cntx			/* Context for callback */
);

/* - - - - - - - - - - - - - - - - - - -- */

#ifdef __cplusplus
	}
#endif

#define VTPGLUT_H
#endif /* VTPGLUT_H */
