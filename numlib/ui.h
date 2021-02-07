
#ifndef __UI_H__

/*
 * Do OS specific setup for using UI windows.
 * Include this only in file with main().
 *
 * Copyright 2014 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Typically we need to set things up and then call the
 * "normal" main, called "uimain" in ArgyllCMS utils.
 */

/* Because linkers are dumb (they pass over a library containing main() */
/* and then complain that there is no entry point!), we force linking */
/* of libui if ui.h gets #included. */
extern int ui_initialized;		
static int *pui_initialized = &ui_initialized;

/* Call this if we decide we are actually going to display */
/* something in the GUI */
void ui_UsingGUI();

#ifdef UNIX
# ifdef __APPLE__

extern pthread_t ui_thid;		/* Thread ID of main thread running io run loop */
extern pthread_t ui_main_thid;	/* Thread ID of thread running application main() */

/* Run a function in the main thread and return when it is complete */
void ui_runInMainThreadAndWait(void *cntx, void (*function)(void *context));

/* We are about to change the UI */
void ui_aboutToWait();

/* Wait until we are sure our UI change is complete */
void ui_waitForEvents();

#ifndef __UI_C__
# define main uimain
#endif

#else	/* Linux etc. */

#endif	/* Linux etc. */
#endif	/* UNIX */

#ifdef NT

#ifdef NEVER 	/* Not practical - see ui.c for full explanation */
#ifndef __UI_C__
# define main uimain
#endif
#endif

#endif

#define __UI_H__
#endif /* __UI_H__ */
