
#ifndef JCONF_H
#define JCONF_H

/*
 * JSON based configuration format class.
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


/* General description: 

	Key names are UNIX style '/' separated paths, stored
	in UTF8 format. The paths should not start with a '/'.

	Duplicate key names are supported, but typically
	will be avoided by the library user.

	Five types are supported, NULL, 32 bit boolean, 64 bit real,
	64 bit integer and string. Strings are always nul terminated.

 */

/* jcnf error codes */
typedef enum {
	jc_ok		    = 0,		/* No error */
	jc_malloc,					/* malloc, calloc or realloc failed */
	jc_lock_error,				/* Error opening lock file */
	jc_locked,					/* File is locked by someone else */
	jc_unlock,					/* Unlock failed */
	jc_noexisting,				/* No existing file to read */
	jc_stat,					/* Unable to stat file */
	jc_changed,					/* File has changed since it was read */
	jc_read_fail,				/* Failure to read from the file */
	jc_parse_fail,				/* Failure to parse from the file */
	jc_write_open,				/* Failed to open file for writing */
	jc_write_fail,				/* Failed to write to file */
	jc_write_close,				/* Failed to close file after writing */
	jc_update_nomod,			/* Attempt to update file that wasn't opened for modifications */
	jc_bad_addkey_params,		/* Bad add_key() parameters */
	jc_unknown_key_type,		/* An unrecognisd key type was encountered */
	jc_ix_oorange,				/* Key index is out of range */
	jc_no_keyname,				/* No key name provided when it is expected */
	jc_string_not_terminated	/* String doesn't include nul */
} jc_error;

/* Argument types */
typedef enum {
	jc_read		= 0,		/* Just read the config file */
	jc_modify	= 1			/* Read the config file in preparation to write it */
} jc_mod;

typedef enum {
	jc_no_create	= 0,	/* Don't create the config if it doesn't exist */
	jc_create		= 1		/* Create the config if it doesn't exist */
} jc_crte;

/* Internal jcnf structure */

/* The different type of values supported */
typedef enum {
	jc_null		= 0,		/* Null value */
	jc_boolean	= 1,		/* Boolean */
	jc_real     = 2,		/* double floating point */
	jc_integer	= 3,		/* 64 bit integer */
	jc_string	= 4			/* UTF8 string, nul terminated */
} jc_type;

/* A value */
struct _jc_key {
	char *key;						/* Key path */
	jc_type type;					/* Type of value */
	char *c_comment;				/* C Comment */
	char *cpp_comment;				/* C++ Comment */
	unsigned char *data;			/* Pointer to data */
	size_t dataSize;				/* Size of data */
}; typedef struct _jc_key jc_key;

/* A recursion depth record used during parsing */
struct _jc_recd {
	char *key;						/* Key name, or */
	int aix;						/* Array index, -2 = no array */
}; typedef struct _jc_recd jc_recd;

/* jcnf Object, representing the keys in a jcnf file */
struct _jcnf {
	jc_key **keys;		/* Array of pointers to keys */
	int nkeys;			/* Number of valid key pointers */
	int akeys;			/* Number of allocated key poiters */
	jc_key *lk;			/* Last key created */

	/* Parsing support, key recursion depth */
	jc_recd *recds;	
	int nrecd;			/* Number of ised recd */
	int arecd;			/* Allocate rec depth */

	/* Config & status */
	char *fname;		/* filename */
	FILE *fp;			/* opened, locked file */
	off_t rsize;		/* Size of file when read */
	time_t rtime;		/* Time the file was read */
	int modify;			/* Opened for modifications */
	int create;			/* Create if it doesn't exist */
	int locked;			/* nz if file is locked */
	int modified;		/* nz if keys have been modified */

	/* Locate the index of the next key matching the key name, starting */
	/* at the given index. Update the index to the matching key. */ 
	/* Look for an exact match if exact != 0, or leading match if exact = 0 */
	/* Search backwards if bwd != 0 or forwards if bwd = 0 */
	/* Set *ix = -1 to begin search from the end. */
	/* Return jc_ix_oorange if no more matchs. */
	jc_error (*locate_key)(struct _jcnf *p, int *ix, char *key, int exact, int bwd);

	/* Retrieve a keys information. Return pointers may be NULL. */
	/* If ix >= 0, return the key of the given index. */
	/* jc_ix_oorange is returned when past end. */
	/* If ix == -1, return the first from the beginning matching key name. */
	/* (Returned data is internal to jcnf object, so call must copy it). */
	jc_error (*get_key)(struct _jcnf *p, int ix, char **key, jc_type *type, unsigned char **data,
                   size_t *dataSize, char **comment);

	/* Set a keys information. */
	/* If ix >= 0, set the key of the given index. */
	/* jc_ix_oorange is returned when past end. */
	/* If ix == -1, overwrite an existing key with the same name, */
	/* or add a new key with that name at the end if there is no existing key. */
	jc_error (*set_key)(struct _jcnf *p, int ix, char *key, jc_type type, unsigned char *data,
                   size_t dataSize, char *comment);

	/* Add a key value to the jcnf at the end, irrespective of whether there is */
	/* an existing key with that name. */
	jc_error (*add_key)(struct _jcnf *p, char *key, jc_type type, unsigned char *data,
	               size_t dataSize, char *comment);

	/* Delete a key. */
	/* If ix >= 0, delete the key of the given index. */
	/* jc_ix_oorange is returned when past end. */
	/* If ix == -1, delete the key with the given name. */ 
	jc_error (*delete_key)(struct _jcnf *p, int ix, char *key);

	/* Diagnostic - Print the value of a key */
	jc_error (*print_key)(struct _jcnf *p, int ix, FILE *fp);

	/* Switch from read only to update of the config file. */
	/* (This re-opens the file and checks that it hasn't been */
	/*  modified since it was read) */
	jc_error (*enable_modify)(struct _jcnf *p);

	/* Save and changes out to the file, unlock it and and close it. */
	/* It can't be udated again after this. */
	jc_error (*update)(struct _jcnf *p);

	/* Delete this object */
	void (*del)(struct _jcnf *p);

}; typedef struct _jcnf jcnf;

/* Create a new jcnf and read it's keys from the file. */
/* Return NULL on error */
jcnf *new_jcnf(
	jc_error *pev,		/* return error code on error */
	char *fname,		/* Corresponding filename */
	jc_mod modify,		/* Flag, nz to open for modification */
	jc_crte create		/* Flag, nz to create if it doesn't exist (modify must be set) */
);

/* Utilities */

/* Return a pointer to the nth element of the key name. */
/* Return null if it is out of range or malloc failed. */
/* Free the returned value when done. */
char *jc_get_nth_elem(char *path, int n);

#endif /* JCNF_H */


