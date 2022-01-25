

#define __UI_C__

/*
 * Do OS specific setup for using UI windows.
 *
 * Copyright 2014 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Typically we need to set things up and then call the
 * "normal" main, called "uimain" in ArgyllCMS utils,
 * created by ui.h #defining main to uimain.
 */

#ifdef UNIX

#if __APPLE__

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* OS X */

/*

	OS X is dumb in relying on an event loop that
	_has_ to be in the main thread to operate.
	Since we can invoke UI operations from any thread,
	the way we work around this is to intercept the
	main thread, spawn a secondary thread to run the
	application main(), and then do nothing but
	service the events in the main thread.
	We can also pass functions back to be run in the main thread,
	and wait for events in the main thread to be processed.

	Note though that Cocoa has poor thread safety :-
	ie. NSRunLoop can't be used to access events - use CFRunLoop

	Should bracket all drawing code between the lockFocusIfCanDraw and
	unlockFocus methods of NSView ? Or is this only a problem if
	different threads try and draw to the same window ?

 */

# include <stdio.h>
# include <stdlib.h>
# include <pthread.h>
# include "ui.h"

# include <Foundation/Foundation.h>
# include <AppKit/AppKit.h>

/* (Duplicate declaration from numsup.h) */

/* Tell App Nap that this is user initiated */
void osx_userinitiated_start();

/* Done with user initiated */
void osx_userinitiated_end();

/* This is a mechanism to force libui to link */
int ui_initialized = 0;

static int g_argc;
static char **g_argv;

pthread_t ui_thid = 0;		/* Thread ID of main thread running io run loop */
pthread_t ui_main_thid = 0;	/* Thread ID of thread running application main() */

static pthread_mutex_t ui_lock1, ui_lock2;	/* Protect wait code */
static pthread_cond_t ui_cond1, ui_cond2;	/* Signal to waiting thread */
static int ui_event1 = 0, ui_event2 = 0;	/* Sync event was received */

extern int uimain(int argc, char *argv[]);

/* Thread that calls the real application main() */
static void *callMain(void *p) {
	int rv;

	ui_main_thid = pthread_self();

	NSAutoreleasePool *tpool = [NSAutoreleasePool new];

	/* Turn App Nap off */
	osx_userinitiated_start();

	/* Should we catch and report exceptions ? */

	rv = uimain(g_argc, g_argv);

	/* Restore App Nap state */
	osx_userinitiated_end();

	[tpool release];

	/* Cleanup, since main thread won't return */
	pthread_cond_destroy(&ui_cond1);
	pthread_cond_destroy(&ui_cond2);
	pthread_mutex_destroy(&ui_lock1);
	pthread_mutex_destroy(&ui_lock2);

	exit(rv);
}

/* Dumy method for NSThread to start */
@interface MainClass : NSObject
+(void)dummyfunc:(id)param;
@end

@implementation MainClass
+(void)dummyfunc:(id)param{
}

@end

int main(int argc, char ** argv) {

	pthread_mutex_init(&ui_lock1, NULL);
	pthread_mutex_init(&ui_lock2, NULL);
	pthread_cond_init(&ui_cond1, NULL);
	pthread_cond_init(&ui_cond2, NULL);

	ui_thid = pthread_self();

	/* Create an NSApp */
	static NSAutoreleasePool *pool = nil;

	pool = [NSAutoreleasePool new];

	[NSApplication sharedApplication];	/* Creates NSApp */

	[NSApp finishLaunching];

	/* We need to create at least one NSThread to tell Cocoa that we are using */
	/* threads, and to protect Cococa objects. (We don't actually have to start the thread.) */
    [NSThread detachNewThreadSelector:@selector(dummyfunc:) toTarget:[MainClass class] withObject:nil];

	/* Call the real main() in another thread */
	int rv;
	pthread_attr_t stackSzAtrbt;
	pthread_t thid;
	g_argc = argc;
	g_argv = argv;

	ui_initialized = 1;

	/* Default stack size is 512K - this is a bit small - raise it */
	if ((rv = pthread_attr_init(&stackSzAtrbt)) != 0
	 || (rv = pthread_attr_setstacksize(&stackSzAtrbt, 8 * 1024 * 1024)) != 0) {
		fprintf(stderr,"ui: thread_attr_setstacksize failed with %d\n",rv);
		return -1;
	}

	if ((rv = pthread_create(&thid, &stackSzAtrbt, callMain, (void *)NULL)) != 0) {
		fprintf(stderr,"ui: pthread_create failed with %d\n",rv);
		return -1;
	}

	/* Service the run queue */
	{
		NSEvent *event;
		NSDate *to;

		/* Process events, looking for application events. */
		for (;;) {
			/* Hmm. Assume to is autorelease */
			to = [NSDate dateWithTimeIntervalSinceNow:1.0];
			/* Hmm. Assume event is autorelease */
			if ((event = [NSApp nextEventMatchingMask:NSAnyEventMask
			             untilDate:to inMode:NSDefaultRunLoopMode dequeue:YES]) != nil) {

				/* call function message */
				if ([event type] == NSApplicationDefined
				 && [event subtype] == 1) {
					void *cntx = (void *)[event data1];
					void (*function)(void *cntx) = (void (*)(void *)) [event data2];

					function(cntx);

					pthread_mutex_lock(&ui_lock1);
					ui_event1 = 1;
					pthread_cond_signal(&ui_cond1);
					pthread_mutex_unlock(&ui_lock1);
					[event release];

				/* event flush message */
				} else if ([event type] == NSApplicationDefined
				 && [event subtype] == 2) {
					pthread_mutex_lock(&ui_lock2);
					ui_event2 = 1;
					pthread_cond_signal(&ui_cond2);
					pthread_mutex_unlock(&ui_lock2);
					[event release];

				/* Everything else */
				} else {
					[NSApp sendEvent:event];
				}
			}
		}
   	 }

	/* Note that we don't actually clean this up on exit - */
	/* possibly we can't. */
//	[NSApp terminate: nil];
}

/* Call this if we decide we are actually going to display something in the GUI. */
/* We switch to "interact with the Dock" mode. */
void ui_UsingGUI() {
	static int attached = 0;
	ProcessSerialNumber psn = { 0, 0 };

	if (!attached) {
#if MAC_OS_X_VERSION_MAX_ALLOWED >= 1060
		/* Make the application appear in the Dock, and  interact with the desktop properly. */
		/* (Unbundled applications default to NSApplicationActivationPolicyProhibited) */
		[NSApp setActivationPolicy:NSApplicationActivationPolicyRegular];
	
#else
# if MAC_OS_X_VERSION_MAX_ALLOWED >= 1030
		/* Make the application appear in the Dock, and  interact with the desktop properly. */
		/* We don't need resources or a bundle if we do this. */
		if (GetCurrentProcess(&psn) == noErr) {
			OSStatus stat;
			if (psn.lowLongOfPSN != 0 && (stat = TransformProcessType(&psn,
				               kProcessTransformToForegroundApplication)) != noErr) {
				fprintf(stderr,"TransformProcess failed with code %d\n",stat);
	
				/* An older trick uses an undocumented API:
				CPSEnableForegroundOperation(&processSerialNum, 0, 0, 0, 0);
				*/
			} else {
	//			fprintf(stderr,"TransformProcess suceeded\n");
			}
		}
# endif /* OS X 10.3 */
#endif /* !OS X 10.6 */
	
		/* We seem to need this, because otherwise we don't get focus automatically */
		[NSApp activateIgnoringOtherApps: YES];

		attached = 1;
	}
}

/* Run a function in the main thread and return when it is complete. */
/* (It's up to the function to record it's result status in its context) */
void ui_runInMainThreadAndWait(void *cntx, void (*function)(void *cntx)) {

	NSEvent *event;
	NSPoint point = { 0.0, 0.0 };
	int rv;

	pthread_mutex_lock(&ui_lock1);
	ui_event1 = 0;

	event = [NSEvent otherEventWithType:NSApplicationDefined
                               location:point
                          modifierFlags:0
                              timestamp:0.0
                           windowNumber:0
                                context:nil
                                subtype:1
                                  data1:(long)cntx				/* long same size as * */
                                  data2:(long)function];
	[NSApp postEvent:event atStart:NO];

	// unlock and wait for signal
	for (;;) {	/* Ignore spurious wakeups */
		if ((rv = pthread_cond_wait(&ui_cond1, &ui_lock1)) != 0) {
			break;		// Hmm.
		}
		if (ui_event1)	/* Got what we were waiting for */
			break;
	}
	pthread_mutex_unlock(&ui_lock1);
}


/* We are about to change the UI */
void ui_aboutToWait() {

	pthread_mutex_lock(&ui_lock2);
	ui_event2 = 0;
}

/* Wait until we are sure our UI change is complete, */
/* because our event has trickled through. */
void ui_waitForEvents() {
	NSEvent *event;
	NSPoint point = { 0.0, 0.0 };
	int rv;

	event = [NSEvent otherEventWithType:NSApplicationDefined
                               location:point
                          modifierFlags:0
                              timestamp:0.0
                           windowNumber:0
                                context:nil
                                subtype:2
                                  data1:0
                                  data2:0];
	[NSApp postEvent:event atStart:NO];

	// unlock and wait for signal
	for (;;) {	/* Ignore spurious wakeups */
		if ((rv = pthread_cond_wait(&ui_cond2, &ui_lock2)) != 0) {
			break;		// Hmm.
		}
		if (ui_event2)	/* Got what we were waiting for */
			break;
	}
	pthread_mutex_unlock(&ui_lock2);
}

#else /* !APPLE */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* UNIX */

/* This is a mechanism to force libui to link */
int ui_initialized = 1;			/* Nothing needs initializing */

/* Call this if we decide we are actually going to display */
/* something in the GUI */
void ui_UsingGUI() {
}

#endif /* !APPLE */

#endif /* UNIX */

#ifdef NT
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* NT */

#ifdef NEVER

/*
	This code isn't generally usable for console programs, because
	Microsoft changed the default behaviour of cmd.exe in Vista to
	not wait for a gui mode .exe to terminate. A non-gui .exe always
	pops a new console if run from explorer.

	So this approach is still viable for primarily gui programs
	to output debug to stdout, but a different shell or "cmd.exe /E:OFF"
	is needed to interact with it via the shell.

	Alternatives are messy - mark the exe as console and
	have it shut down the popup console in gui mode (looks
	ugly), or have a stub app.bat that invokes the
	real .exe, forcing cmd.exe to operate in the mode where
	it waits for execution to finish.
*/

/* This is a mechanism to force libui to link */
int ui_initialized = 0;

/*

	On MSWin we can rely on WinMain to be called instead of
	main(), so that we can re-attache the stdio so that the
	resulting exe works the same when involked either from
	a shell, or directly from explorer.

 */

/* May have to add link flag -Wl,-subsystem,windows */
/* since MingW is stupid about noticing WinMain or pragma */

#include <windows.h>
#include <stdio.h>

# pragma comment( linker, "/subsystem:windows" )
//# pragma comment( linker, "/subsystem:console /ENTRY:WinMainCRTStartup" )
//# pragma comment( linker, "/subsystem:windows /ENTRY:mainCRTStartup" )
//# pragma comment( linker, "/subsystem:windows /ENTRY:WinMainCRTStartup" )

extern int uimain(int argc, char *argv[]);

APIENTRY WinMain(
    HINSTANCE hInstance,
    HINSTANCE hPrevInstance,
    LPSTR lpCmdLine,
    int nCmdShow
) {
	{	/* Only works on >= XP though */
		BOOL (WINAPI *AttachConsole)(DWORD dwProcessId);

		*(FARPROC *)&AttachConsole = 
          GetProcAddress(LoadLibraryA("kernel32.dll"), "AttachConsole");

		if (AttachConsole != NULL && AttachConsole(((DWORD)-1)))
		{
			if (_fileno(stdout) < 0)
				freopen("CONOUT$","wb",stdout);
			if (_fileno(stderr) < 0)
				freopen("CONOUT$","wb",stderr);
			if (_fileno(stdin) < 0)
				freopen("CONIN$","rb",stdin);
#ifdef __cplusplus 
			// make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog point to console as well
			std::ios::sync_with_stdio();
#endif
		}
	}

	ui_initialized = 1;

	return uimain(__argc, __argv);
}

#else /* !NEVER */

/* This is a mechanism to force libui to link */
int ui_initialized = 1;			/* Nothing needs initializing */

#endif /* !NEVER */

/* Call this if we decide we are actually going to display */
/* something in the GUI */
void ui_UsingGUI() {
}

#endif	/* NT */
