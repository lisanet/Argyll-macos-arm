
#ifndef UCMM_H
#define UCMM_H

/* 
 * Unix micro-cmm to manage X11 display
 * calibration and profile loading.
 */

/*************************************************************************
 Copyright 2008 Graeme W. Gill

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.

 *************************************************************************/

typedef enum {
	ucmm_ok = 0,
	ucmm_resource,				/* Resource failure (e.g. out of memory) */
	ucmm_invalid_profile,		/* Profile is not a valid display ICC profile */
	ucmm_no_profile,			/* There is no associated profile */
	ucmm_no_home,				/* There is no HOME environment variable defined */
	ucmm_no_edid_or_display,	/* There is no edid or display name */
	ucmm_profile_copy,			/* There was an error copying the profile */
	ucmm_open_config,			/* There was an error opening the config file */
	ucmm_access_config,			/* There was an error accessing the config information */
	ucmm_set_config,			/* There was an error setting the config file */
	ucmm_save_config,			/* There was an error saving the config file */
	ucmm_monitor_not_found,		/* The EDID or display wasn't matched */
	ucmm_delete_key,			/* Delete_key failed */
	ucmm_delete_profile,		/* Delete_key failed */
} ucmm_error;

/* Install scope */
typedef enum {
	ucmm_user,
	ucmm_local_system
} ucmm_scope;

/* Install a profile for a given monitor */
/* Either EDID or display_name may be NULL, but not both. */
/* Any existing association is overwritten. Installed profiles */
/* are not deleted. */
ucmm_error ucmm_install_monitor_profile(
	ucmm_scope scope,				/* Scope of instalation */
	unsigned char *edid,	/* Primary device identifier, NULL if none. */
	int edid_len,			/* Length of edid data */
	char *display_name,	/* Fall back device association, */
	                            /* the X11 display name */
	char *profile			/* Path to profile to be installed. */
);

/* Un-install a profile for a given monitor. */
/* Either EDID or display_name may be NULL, but not both. */
/* The monitor is left with no profile association. */
/* Return an error code */
ucmm_error ucmm_uninstall_monitor_profile(
	ucmm_scope scope,				/* Scope of instalation */
	unsigned char *edid,	/* Primary device identifier, NULL if none. */
	int edid_len,			/* Length of edid data */
	char *display_name		/* Fall back device association, */
);

/* Get an associated monitor profile. */
/* Return ucmm_no_profile if there is no installed profile for this */
/* monitor. */
/* Return an error code */
ucmm_error ucmm_get_monitor_profile(
	unsigned char *edid,	/* Primary device identifier, NULL if none. */
	int edid_len,			/* Length of edid data */
	char *display_name,	/* Fall back device association, */
	char **profile		/* Return path to profile. free() afterwards. */
);

/* Return an ASCII error message string interpretation of scope enum */
char *ucmm_scope_string(ucmm_scope scope);

/* Return an ASCII error message string interpretation of an error number */
char *ucmm_error_string(ucmm_error erno);

#endif /* UCMM_H */
