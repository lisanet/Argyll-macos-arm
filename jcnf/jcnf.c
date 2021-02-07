
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

#define LOCK_RETRIES 5

/* We convert JSON arrays to numerical sub paths on reading, */
/* but they get converted back to maps on writing. */
/* Paths are grouped under common paths on writing, */
/* but aren't sorted. */

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifndef NT
# include <sys/file.h>
# include <unistd.h>
#endif
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "yajl_common.h"
#include "yajl_gen.h"
#include "yajl_parse.h"

#include "jcnf.h"

#define BUF_SIZE 2048

/* - - - - - - - - - - - - - - - - */
/* Utilities */

/* Construct a string that represents the current key path */
/* Return NULL on error (malloc fail) */
/* Free the returned string when done */
static char *cur_keypath(jcnf *p) {
	char *ckey = NULL;
	int i, len;

	if (p->nrecd <= 0) {
		return strdup("");
	}
	for (len = i = 0; i < p->nrecd; i++) {
		if (p->recds[i].aix > -2) {
			len += 3 + (int)log10((double)p->recds[p->nrecd-1].aix);
		} else {
			len += strlen(p->recds[i].key) + 1;
		}
	}

	if ((ckey = malloc(len)) == NULL)
		return NULL;
	
	for (len = i = 0; i < p->nrecd; i++) {
		int sl;

		if (p->recds[i].aix > -2) {
			char num[13];
			sprintf(num, "%d",p->recds[i].aix);
			sl = strlen(num);
			memmove(ckey + len, num, sl);
		} else {
			sl = strlen(p->recds[i].key);
			memmove(ckey + len, p->recds[i].key, sl);
		}
		len += sl;

		if ((i+1) >= p->nrecd)
			ckey[len] = '\000';
		else
			ckey[len] = '/';
		len++;
	}
	return ckey;
}

/* Diagnostic - Print the value of a key */
static jc_error jcnf_print_key(jcnf *p, int ix, FILE *fp) {
	jc_key *kp;

	if (ix < 0 || ix >= p->nkeys)
		return jc_ix_oorange;

	kp = p->keys[ix];

	fprintf(fp,"Key '%s' has value",kp->key);

	switch(kp->type) {
		case jc_null: {
			fprintf(fp," null");
			break;
		}
		case jc_boolean: {
			fprintf(fp," %s",*((int *)kp->data) ? "true" : "false");
			break;
		}
		case jc_real: {
			fprintf(fp," %lf",*((double *)kp->data));
			break;
		}
		case jc_integer: {
			fprintf(fp," %ld",*((long *)kp->data));
			break;
		}
		case jc_string: {
			fprintf(fp," '%s'",kp->data);
			break;
		}
		default: {
			fprintf(fp," unknown type %d",kp->type);
			break;
		}
	}
	if (kp->c_comment != NULL) {
		fprintf(fp," C comment = '%s'",kp->c_comment);
	}
	if (kp->cpp_comment != NULL) {
		fprintf(fp," C++ comment = '%s'",kp->cpp_comment);
	}
	fprintf(fp,"\n");
	return jc_ok;
}

/* - - - - - - - - - - - - - - - - */
/* Free keys contents */
static void free_key_contents(jc_key *p) {
	if (p->key != NULL) {
		free(p->key);
		p->key = NULL;
	}
	if (p->c_comment != NULL) {
		free(p->c_comment);
		p->c_comment = NULL;
	}
	if (p->cpp_comment != NULL) {
		free(p->cpp_comment);
		p->cpp_comment = NULL;
	}
	if (p->data != NULL) {
		free(p->data);
		p->data = NULL;
	}
	p->type = -1;
}

/* Free key and its contents */
static void free_key(jc_key *p) {
	free_key_contents(p);
	free(p);
}

/* Set a keyvalues (internal) */
/* All parameters are copied */
/* return NZ on error */
static jc_error jcnf_set_key_internal(
	jcnf *p,
	jc_key *kp,			/* Pointer to key */
	char *key,			/* Key path name */
	jc_type type,
	void *data,
	size_t dataSize,	/* Data size (including string nul) */
	char *comment		/* C++ style comment */
) {
	free_key_contents(kp);

	if ((kp->key = strdup(key)) == NULL)
		return jc_malloc;
	kp->type = type;

	if (type != jc_null) {
		if (type == jc_string && ((char *)data)[dataSize-1] != '\000')
			return jc_string_not_terminated;
		if ((kp->data = malloc(dataSize)) == NULL)
			return jc_malloc;
		kp->dataSize = dataSize;
		memmove(kp->data, data, dataSize);
	}

	if (comment != NULL) {
		if ((kp->cpp_comment = strdup(comment)) == NULL)
			return jc_malloc;
	}
	return jc_ok;
}

/* Add a key value to the jcnf (internal) */
/* All parameters are copied */
/* return NZ on error */
static jc_error jcnf_add_key_internal(
	jcnf *p,
	char *key,			/* Key path name */
	jc_type type,
	void *data,
	size_t dataSize,	/* Data size (including string nul) */
	char *comment
) {
	jc_error ev;
	jc_key *kp;

	if (key == NULL || (type != jc_null && (data == NULL || dataSize == 0)))
		return jc_bad_addkey_params;

	if (p->nkeys >= p->akeys) {	 /* Need more pointer space */
		p->akeys = p->akeys ? 2 * p->akeys : 10;
		if ((p->keys = realloc(p->keys, p->akeys * sizeof(jc_key*))) == NULL) {
			return jc_malloc;
		}
	}
	if ((kp = p->keys[p->nkeys] = calloc(1, sizeof(jc_key))) == NULL) {
		return jc_malloc;
	}
	p->nkeys++;

	if ((ev = jcnf_set_key_internal(p, kp, key, type, data, dataSize, comment)) != jc_ok)
		return ev;
	p->lk = kp;
	return jc_ok;
}

/* Locate the index of the next key matching the key name, starting */
/* at the given index. Update the index to the matching key. */ 
/* Look for an exact match if exact != 0, or leading match if exact = 0 */
/* Search backwards if bwd != 0 or forwards if bwd = 0 */
/* Set *ix = -1 to begin search from the start/end for fwd/bwd. */
/* Return jc_ix_oorange if no more match. */
static jc_error jcnf_locate_key(
	jcnf *p,
	int *pix,
	char *key,
	int exact,
	int bwd
) {
	int ix = *pix;
	int sl;

	if (ix == -1) {
		if (bwd)
			ix = p->nkeys-1;
		else
			ix = 0;
	}

	if (ix < 0 || ix >= p->nkeys)
		return jc_ix_oorange;

	if (key == NULL)
		return jc_no_keyname;

	sl = strlen(key);

	/* We're doing a simple linear search */
	for (; ix >= 0 && ix < p->nkeys; ix += bwd ? -1 : 1) {
		if (exact) {
			if (strcmp(key, p->keys[ix]->key) == 0)
				break;
		} else {
			if (strncmp(key, p->keys[ix]->key, sl) == 0)
				break;
		}
	}

	if (ix < 0 || ix >= p->nkeys)
		return jc_ix_oorange;

	*pix = ix;
	return jc_ok;
}

/* Retrieve a keys information. Return pointers may be NULL. */
/* If ix >= 0, return the key of the given index. */
/* If ix == -1, return the first from the beginning exactly matching key name. */
/* jc_ix_oorange is returned when past end. */
/* (Do not modify anything returned) */
static jc_error jcnf_get_key(
	jcnf *p,
	int ix,	
	char **key,	
	jc_type *type,	
	unsigned char **data,
	size_t *dataSize,		/* Data size (including string nul) */
	char **comment
) {
	jc_error ev;

	if (ix == -1) {
		if (key == NULL || *key == NULL)
			return jc_no_keyname;
		if ((ev = jcnf_locate_key(p, &ix, *key, 1, 0)) != jc_ok)
			return ev;
	} else if (ix < 0 || ix >= p->nkeys)
		return jc_ix_oorange;

	if (key != NULL)
		*key = p->keys[ix]->key;
	if (type != NULL)
		*type = p->keys[ix]->type;
	if (data != NULL)
		*data = p->keys[ix]->data;
	if (dataSize != NULL)
		*dataSize = p->keys[ix]->dataSize;
	if (comment != NULL)
		*comment = p->keys[ix]->cpp_comment;

	return jc_ok;
}

/* Set a keys information. */
/* If ix >= 0, set the key of the given index. */
/* jc_ix_oorange is returned when past end. */
/* If ix == -1, overwrite an existing key with the same name, */
/* or add a new key with that name at the end if there is no existing key. */
static jc_error jcnf_set_key(
	jcnf *p,
	int ix,
	char *key,
	jc_type type,
	unsigned char *data,
	size_t dataSize,	/* Data size (including string nul) */
	char *comment
) {
	jc_error ev;

	if (ix == -1) {
		if (key == NULL)
			return jc_no_keyname;
		if ((ev = jcnf_locate_key(p, &ix, key, 1, 0)) != jc_ok) {
			if (ev != jc_ix_oorange)
				return ev;
			if (p->modify == 0)
				return jc_update_nomod;
			p->modified = 1;
			return jcnf_add_key_internal(p, key, type, data, dataSize, comment);
		}
	} else if (ix < 0 || ix >= p->nkeys)
		return jc_ix_oorange;

	if (p->modify == 0)
		return jc_update_nomod;
	p->modified = 1;
	return jcnf_set_key_internal(p, p->keys[ix], key, type, data, dataSize, comment);
}

/* Add a key value to the jcnf at the end, irrespective of whether there is */
/* an existing key with that name. */
static jc_error jcnf_add_key(
	jcnf *p,
	char *key,
	jc_type type,
	unsigned char *data,
	size_t dataSize,		/* Data size (including string nul) */
	char *comment
) {
	jc_error ev;
	if (p->modify == 0)
		return jc_update_nomod;
	if ((ev = jcnf_add_key_internal(p,key,type,data,dataSize,comment)) == jc_ok) {
		p->modified = 1;
	}
	return ev;
}

/* Delete a key. */
/* If ix >= 0, delete the key of the given index. */
/* jc_ix_oorange is returned when past end. */
/* If ix == -1, delete the  key with the given name. */ 
static jc_error jcnf_delete_key(
	jcnf *p,
	int ix,
	char *key
) {
	jc_error ev;

	if (ix == -1) {
		if (key == NULL)
			return jc_no_keyname;
		if ((ev = jcnf_locate_key(p, &ix, key, 1, 0)) != jc_ok)
			return ev;
	} else if (ix < 0 || ix >= p->nkeys)
		return jc_ix_oorange;

	if (p->modify == 0)
		return jc_update_nomod;
	free_key_contents(p->keys[ix]);
	
	if ((ix+1) < p->nkeys) {
		memmove(p->keys+ix, p->keys+ix+1,sizeof(jc_key *) * p->nkeys-ix-1);
	}
	free(p->keys[p->nkeys-1]);
	p->nkeys--;
	p->modified = 1;

	return jc_ok;
}

/* - - - - - - - - - - - - - - - - */
/* yajl parser callbacks */

static int jcnf_yajl_null(void *ctx) {
	jcnf *p = (jcnf *)ctx;
	char *t1;

//    printf("null\n");

	/* If current level is an array, bump the index */
	if (p->nrecd > 0 && p->recds[p->nrecd-1].aix > -2)
		p->recds[p->nrecd-1].aix++;
	if ((t1 = cur_keypath(p)) == NULL)
		return 0;

	if (jcnf_add_key_internal(p, t1, jc_null, NULL, 0, NULL) != jc_ok) {
		free(t1);
		return 0;
	}
	free(t1);

    return 1;
}

static int jcnf_yajl_boolean(void * ctx, int boolVal) {
	jcnf *p = (jcnf *)ctx;
	char *t1;

//    printf("bool: %s\n", boolVal ? "true" : "false");

	/* If current level is an array, bump the index */
	if (p->nrecd > 0 && p->recds[p->nrecd-1].aix > -2)
		p->recds[p->nrecd-1].aix++;
	if ((t1 = cur_keypath(p)) == NULL)
		return 0;

	if (jcnf_add_key_internal(p, t1, jc_boolean, (void *)&boolVal, sizeof(int), NULL) != jc_ok) {
		free(t1);
		return 0;
	}
	free(t1);

    return 1;
}

static int jcnf_yajl_integer(void *ctx, longlong integerVal) {
	jcnf *p = (jcnf *)ctx;
	char *t1;

//    printf("integer: %lld\n", integerVal);

	/* If current level is an array, bump the index */
	if (p->nrecd > 0 && p->recds[p->nrecd-1].aix > -2)
		p->recds[p->nrecd-1].aix++;
	if ((t1 = cur_keypath(p)) == NULL)
		return 0;

	if (jcnf_add_key_internal(p, t1, jc_integer, (void *)&integerVal, sizeof(long), NULL)  != jc_ok) {
		free(t1);
		return 0;
	}
	free(t1);

    return 1;
}

static int jcnf_yajl_double(void *ctx, double doubleVal) {
	jcnf *p = (jcnf *)ctx;
	char *t1;

//    printf("double: %lf\n", doubleVal);

	/* If current level is an array, bump the index */
	if (p->nrecd > 0 && p->recds[p->nrecd-1].aix > -2)
		p->recds[p->nrecd-1].aix++;
	if ((t1 = cur_keypath(p)) == NULL)
		return 0;

	if (jcnf_add_key_internal(p, t1, jc_real, (void *)&doubleVal, sizeof(double), NULL)  != jc_ok) {
		free(t1);
		return 0;
	}
	free(t1);

    return 1;
}

static int jcnf_yajl_string(void *ctx, const unsigned char *stringVal,
                     size_t stringLen) {
	jcnf *p = (jcnf *)ctx;
	char *t1, *t2;

//	printf("string: '");
//	fwrite(stringVal, 1, stringLen, stdout);
//	printf("'\n");    

	/* If current level is an array, bump the index */
	if (p->nrecd > 0 && p->recds[p->nrecd-1].aix > -2)
		p->recds[p->nrecd-1].aix++;
	if ((t1 = cur_keypath(p)) == NULL)
		return 0;

    if ((t2 = malloc(stringLen + 1)) == NULL)		/* Room to add implicit nul */
		return 0;
    memmove(t2, stringVal, stringLen);
    t2[stringLen] = 0;

	if (jcnf_add_key_internal(p, t1, jc_string, (void *)t2, stringLen+1, NULL)  != jc_ok) {
		free(t2);
		free(t1);
		return 0;
	}
	free(t2);
	free(t1);

    return 1;
}

static int jcnf_yajl_c_comment(void *ctx, const unsigned char * stringVal,
                     unsigned int stringLen) {
	jcnf *p = (jcnf *)ctx;

//	printf("c_comment: '");
//	fwrite(stringVal, 1, stringLen, stdout);
//	printf("'\n");    

	if (p->lk != NULL && p->lk->c_comment == NULL) {
		char *t1;
		if ((t1 = malloc(stringLen + 1)) == NULL)
			return 0;
		memmove(t1, stringVal, stringLen);
		t1[stringLen] = 0;

		p->lk->c_comment = t1;
	}

    return 1;
}

static int jcnf_yajl_cpp_comment(void *ctx, const unsigned char * stringVal,
                     unsigned int stringLen) {
	jcnf *p = (jcnf *)ctx;

//	printf("cpp_comment: '");
//	fwrite(stringVal, 1, stringLen, stdout);
//	printf("'\n");    

	if (p->lk != NULL && p->lk->cpp_comment == NULL) {
		char *t1;
		if ((t1 = malloc(stringLen + 1)) == NULL)
			return 0;
		memmove(t1, stringVal, stringLen);
		t1[stringLen] = 0;

		p->lk->cpp_comment = t1;
	}

    return 1;
}

static int jcnf_yajl_start_map(void *ctx) {
	jcnf *p = (jcnf *)ctx;

	/* Start another recursion level */
	if (p->nrecd >= p->arecd) {
		p->arecd *= 2;
		if ((p->recds = (jc_recd *)realloc(p->recds, p->arecd * sizeof(jc_recd))) == NULL) {
			/* realloc failed */
			return 0;
		}
	}

	/* If current level is an array, bump the index */
	if (p->nrecd > 0 && p->recds[p->nrecd-1].aix > -2)
		p->recds[p->nrecd-1].aix++;

	/* Move to new level */
    p->recds[p->nrecd].key = NULL;
    p->recds[p->nrecd].aix = -2;
	p->nrecd++;

#ifdef NEVER
    printf("map open '{'\n");
#endif

    return 1;
}

/* Callback from yajl */
static int jcnf_yajl_map_key(void *ctx, const unsigned char * stringVal,
                     size_t stringLen) {
	jcnf *p = (jcnf *)ctx;
	int i;

	if (stringLen == 0) {
		/* Zero length key */
		return 0;
	}

	if (p->recds[p->nrecd-1].key != NULL) {
		free(p->recds[p->nrecd-1].key);
		p->recds[p->nrecd-1].key = NULL;
		p->recds[p->nrecd-1].aix = -2;
	}
    if ((p->recds[p->nrecd-1].key = malloc(stringLen + 1)) == NULL)
		return 0;

    p->recds[p->nrecd-1].key[stringLen] = 0;
    memmove(p->recds[p->nrecd-1].key, stringVal, stringLen);

#ifdef NEVER
	{
		char *tt;

	   	printf("Added path = '%s'\n", p->recds[p->nrecd-1].key);
		if ((tt = cur_keypath(p)) != NULL) {
			printf("Current path = '%s'\n",tt);
			free(tt);
		}
	}
#endif

    return 1;
}

static int jcnf_yajl_end_map(void *ctx) {
	jcnf *p = (jcnf *)ctx;

	if (p->nrecd == 0) {
		/* End of map without start of map */
		return 0;
	}
	if (p->recds[p->nrecd-1].key != NULL)
		free(p->recds[p->nrecd-1].key);
	p->nrecd--;

#ifdef NEVER
	{
		char *tt;
   	 printf("map close '}'\n");
		if ((tt = cur_keypath(p)) != NULL)
			printf("Current path = '%s'\n",tt);
		free(tt);
	}
#endif

    return 1;
}

static int jcnf_yajl_start_array(void *ctx) {
	jcnf *p = (jcnf *)ctx;

#ifdef NEVER
    printf("array open '['\n");
#endif

	/* Start another recursion level */
	if (p->nrecd >= p->arecd) {
		p->arecd *= 2;
		if ((p->recds = (jc_recd *)realloc(p->recds, p->arecd * sizeof(jc_recd))) == NULL) {
			/* realloc failed */
			return 0;
		}
	}

	/* If current level is an array, bump the index */
	if (p->nrecd > 0 && p->recds[p->nrecd-1].aix > -2)
		p->recds[p->nrecd-1].aix++;

	/* Move to new level */
    p->recds[p->nrecd].key = NULL;
    p->recds[p->nrecd].aix = -1;
	p->nrecd++;

    return 1;
}

static int jcnf_yajl_end_array(void *ctx) {
	jcnf *p = (jcnf *)ctx;

#ifdef NEVER
    printf("array close ']'\n");
#endif
	if (p->nrecd == 0) {
		/* End of map without start of map */
		return 0;
	}
	p->nrecd--;

    return 1;
}

static yajl_callbacks callbacks = {
    jcnf_yajl_null,
    jcnf_yajl_boolean,
    jcnf_yajl_integer,
    jcnf_yajl_double,
	NULL,				/* number */
    jcnf_yajl_string,
    jcnf_yajl_c_comment,
    jcnf_yajl_cpp_comment,
    jcnf_yajl_start_map,
    jcnf_yajl_map_key,
    jcnf_yajl_end_map,
    jcnf_yajl_start_array,
    jcnf_yajl_end_array
};

/* - - - - - - - - - - - - - - - - */
/* Lock the open fp */
static jc_error jcnf_lock_file(jcnf *p) {
#ifndef NT
	int i, fh;
	int lop;

	fh = fileno(p->fp);

	if (p->modify) {
		lop = LOCK_EX | LOCK_NB;
	} else {
		lop = LOCK_SH | LOCK_NB;
	}

	for (i = 0; i < LOCK_RETRIES; i++) {
		if (flock(fh, lop) == 0)
			break;
		sleep(1);
	}
	if (i >= LOCK_RETRIES)
		return jc_locked;
#endif
	p->locked = 1;

	return jc_ok;
}

/* Read a file into the object. */
/* (Doesn't do locking) */
/* Return NZ on error */
static jc_error jcnf_read(
	jcnf *p
) {
	jc_error ev;
	yajl_handle hand;
	unsigned char buf[BUF_SIZE];	
	struct stat sbuf;

	if ((p->fp = fopen(p->fname, p->modify ? "r+" : "r")) == NULL) {
		if (!p->modify)
			return jc_noexisting;
	
		if ((p->fp = fopen(p->fname, "w")) == NULL)
			return jc_write_open;
		if ((ev = jcnf_lock_file(p)) != jc_ok)
			return ev;
		return jc_ok;
	}

	if ((ev = jcnf_lock_file(p)) != jc_ok)
		return ev;

	hand = yajl_alloc(&callbacks, NULL, (void *)p);

	yajl_config(hand, yajl_allow_comments, 1);
	yajl_config(hand, yajl_dont_validate_strings, 1);

	/* Parse the file */
	for(;;) {
		size_t rd;

		rd = fread(buf, 1, BUF_SIZE, p->fp);
		
		if (rd == 0) {
			if (!feof(p->fp)) {
				fprintf(stderr, "error reading from '%s'\n", p->fname);
				return jc_read_fail;
			}
			break;
		} else {
			yajl_status stat;

			/* read file data, pass to parser */
			stat = yajl_parse(hand, buf, rd);
			if (stat != yajl_status_ok) {
				unsigned char * str = yajl_get_error(hand, 1, buf, rd);
				fflush(stdout);
				fprintf(stderr, "%s", (char *) str);
				yajl_free_error(hand, str);
				return jc_parse_fail;
			}
		}
	} 

	/* Record some things about the file so that we can */
	/* recognized if it's been modified */
	if (stat(p->fname, &sbuf) != 0) {
		return jc_stat;
	}
	p->rsize = sbuf.st_size;
	p->rtime = time(NULL);

	yajl_free(hand);

	/* We're not modifying it, so close it */
	if (!p->modify) {
		fclose(p->fp);
		p->fp = NULL;
		p->locked = 0;
	}

	return jc_ok;
}

/* Write an object into a file */
/* Unlock it & close it afterwards. */
/* Return NZ on error */
static jc_error jcnf_write(
	jcnf *p
) {
	FILE *fp;		/* For temporary file */
	char *tname = NULL;
	yajl_gen g;
    yajl_status stat;
	const unsigned char * buf;
	size_t len;
	int clevel = 0;			/* Current level */
	char *pkey = "";		/* Previous key */
	char *ckey;				/* Current key */
	int i;

	if (!p->locked && p->fp != NULL) {
		return jc_update_nomod;
	}

#ifndef NT
	{
		int fh;

		if ((tname = malloc(strlen(p->fname) + 8)) == NULL)
			return jc_malloc;

		/* Create temporary file, open it and lock it LOCK_EX */
		strcpy(tname, p->fname);
		strcat(tname,"-XXXXXX");
		if ((fh = mkstemp(tname)) == -1) {
			free(tname);
			return jc_write_open;
		}
		if (fchmod(fh, 0644) != 0) {
			free(tname);
			return jc_write_open;
		}
		if ((fp = fdopen(fh, "w")) == NULL) {
			free(tname);
			return jc_write_open;
		}
	}
#else
	/* Open a temporary file in the same directory to write to */
	if ((tname = malloc(strlen(p->fname) + 8)) == NULL)
		return jc_malloc;

	if (tmpnam(tname) == NULL) { 
		free(tname);
		return jc_write_open;
	}
	if ((fp = fopen(tname, "w")) == NULL) {
		free(tname);
		return jc_write_open;
	}
#endif

	g = yajl_gen_alloc(NULL);

	yajl_gen_config(g, yajl_gen_beautify, 1);
	yajl_gen_config(g, yajl_gen_indent_string, "\t");
	yajl_gen_config(g, yajl_gen_validate_utf8, 1);

	/* Generate the file */
	for (i = 0; i < p->nkeys; i++, pkey = ckey) { 
		char *pc, *cc, *dc;
		int nplev;				/* Number of previous different levels */
		int ndlev;				/* Number of new different levels */
		int same;

		ckey = p->keys[i]->key;
		
		/* See how many keys are not in common */
		nplev = ndlev = 0;
		same = 1;
		for(pc = pkey, dc = cc = ckey; *pc != '\000' || *cc != '\000';) {

//printf("~1 pc = '%c', cc = '%c', same = %d\n",*pc,*cc,same);
			if (same == 0 || *pc != *cc) {
				same = 0;
				if (*cc == '/') {
					ndlev++;
//printf("~1 ndlev now %d\n",ndlev);
				}
				if (*pc == '/') {
					nplev++;
//printf("~1 ndlev now %d\n",ndlev);
				}
			
			} else {
				if (*cc == '/') {
					dc = cc+1;
				}
			}
			if (*pc != '\000') {
				pc++;
				if (same == 0 && *pc == '\000') {
					nplev++;
//printf("~1 nplev now %d\n",nplev);
				}
			}
			if (*cc != '\000') {
				cc++;
				if (same == 0 && *cc == '\000') {
					ndlev++;
//printf("~1 ndlev now %d\n",ndlev);
				}
			}
		}
//printf("~1 Prev = '%s'\n",pkey);
//printf("~1 Curr = '%s'\n",ckey);
//printf("~1 New = '%s'\n",dc);
//printf("~1 Old different = %d, new different = %d\n\n",nplev,ndlev);

		while(nplev > 0) {
			if (nplev > 1)
				yajl_gen_map_close(g);
			nplev--;
			clevel--;
		}
		while(ndlev > 0) {
			yajl_gen_map_open(g);
			for (cc = dc; *cc != '\000' && *cc != '/'; cc++)
				;
			yajl_gen_string(g, dc, cc-dc);
			if (*cc != '\000')
				dc = cc + 1;
			ndlev--;
			clevel++;
		}

		switch(p->keys[i]->type) {
//printf("~1 key %d = type %d\n",i,p->keys[i]->type);
			case jc_null:
				yajl_gen_null(g);
				break;
			
			case jc_boolean:
				yajl_gen_bool(g, *((int *)p->keys[i]->data));
				break;
			
			case jc_real:
				yajl_gen_double(g, *((double *)p->keys[i]->data));
				break;
			
			case jc_integer:
				yajl_gen_integer(g, *((long *)p->keys[i]->data));
				break;
			
			case jc_string:
				yajl_gen_string(g, (char *)p->keys[i]->data, p->keys[i]->dataSize-1);
				break;

			default: {
				free(tname);
				return jc_unknown_key_type; 
			}
		}

		if (p->keys[i]->cpp_comment != NULL) {
			yajl_gen_cpp_comment(g, p->keys[i]->cpp_comment, strlen(p->keys[i]->cpp_comment));
		}
		if (p->keys[i]->c_comment != NULL) {
			yajl_gen_c_comment(g, p->keys[i]->c_comment, strlen(p->keys[i]->c_comment), 1);
		}

#ifdef NEVER
		yajl_gen_map_open(g);
		yajl_gen_string(g, "test", strlen("test"));
		yajl_gen_string(g, "test value", strlen("test value"));
		yajl_gen_c_comment(g, " A comment ", strlen(" A comment "));
		yajl_gen_map_close(g);
#endif

		/* Do some writing */
		yajl_gen_get_buf(g, &buf, &len);
		if (len >=  BUF_SIZE) {
			if (fwrite(buf, 1, len, fp) != len)
				return jc_write_fail;
			yajl_gen_clear(g);
		}
	} 

	while(clevel > 0) {
		yajl_gen_map_close(g);
		clevel--;
	}
	
	yajl_gen_get_buf(g, &buf, &len);
	if (len > 0) {
		if (fwrite(buf, 1, len, fp) != len)
			return jc_write_fail;
		yajl_gen_clear(g);
	}
	yajl_gen_free(g);
	if (fflush(fp) != 0) {
		free(tname);
		return jc_write_close;
	}

#ifdef NT
	/* MSWindows rename won't replace existing or open files. */
	/* Lucky this is just for testing, as this leaves a window */
	if (fp != NULL) {
		fclose(fp);
		fp = NULL;
	}
	if (p->fp != NULL) {
		fclose(p->fp);
		p->fp = NULL;
	}
	unlink(p->fname);
#endif

//printf("~1 about to rename '%s' to '%s'\n",tname,p->fname);
	/* Now atomicaly rename the file to replace the file we read */
	if (rename(tname, p->fname) != 0) {
		free(tname);
		return jc_write_close;
	}

//printf("~1 closing files\n");
	/* Close our files and release the locks */
	if (fp != NULL && fclose(fp) != 0) {
		free(tname);
		return jc_write_close;
	}
	if (p->fp != NULL && fclose(p->fp) != 0) {
		free(tname);
		return jc_write_close;
	}
	p->fp = NULL;
	p->locked = 0;
	p->modify = 0;

	free(tname);
	return jc_ok;
}

/* Switch from read only to update of the config. */
static jc_error jcnf_enable_modify(jcnf *p) {
	jc_error ev;
	struct stat sbuf;

	if (p->modify)
		return jc_ok;		/* Nothing to do */

	/* We need to re-open the file and lock it for modification */
	if ((p->fp = fopen(p->fname, "r+")) == NULL) {
		return jc_changed;
	}
	p->modify = 1;

	if ((ev = jcnf_lock_file(p)) != jc_ok) {
		p->modify = 0;
		fclose(p->fp);
		p->fp = NULL;
		return ev;
	}

	/* Check that it hasn't been modified since it was first read */
	if (stat(p->fname, &sbuf) != 0) {
		p->modify = 0;
		fclose(p->fp);
		p->fp = NULL;
		return jc_stat;
	}

	if (sbuf.st_size != p->rsize
	 || sbuf.st_mtime > p->rtime) {
		p->modify = 0;
		fclose(p->fp);
		p->fp = NULL;
		return jc_changed;
	}

	return jc_ok;
}


/* Update the file that was opened for modification, and unlock it. */
static jc_error jcnf_update(jcnf *p) {
	jc_error ev;
	struct stat sbuf;

	if (p->modified) {
		if (p->modify == 0 || p->locked == 0) {

			return jc_update_nomod;				/* File wasn't opened for modification */
		}
		if ((ev = jcnf_write(p)) != jc_ok) {
			return ev;
		}
		p->modified = 0;
	}
	p->modify = 0;

	return jc_ok;
}
	
/* free the object */
static void jcnf_del(jcnf *p) {
	
	if (p->fp)
		fclose(p->fp);

	if (p->keys != NULL) {
		int i;
		for (i = 0; i < p->nkeys; i++) { 
			free_key(p->keys[i]);
		}
		free(p->keys);
	}

	if (p->recds != NULL) {
		int i;
		for (i = 0; i < p->nrecd; i++) { 
			free(p->recds[i].key);
		}
		free(p->recds);
	}

	if (p->fname)
		free(p->fname);

	free(p);
}

/* Create a new jconf. */
/* Return NULL on error */
jcnf *new_jcnf(
	jc_error *pev,		/* return error code on error */
	char *fname,	/* Corresponding filename */
	jc_mod modify,		/* Flag, nz to open for modification */
	jc_crte create		/* Flag, nz to create if it doesn't exist (modify must be set) */
) {
	jcnf *p;
	jc_error ev;

	if ((p = (jcnf *) calloc(1, sizeof(jcnf))) == NULL) {
		if (pev != NULL) *pev = jc_malloc;
		return NULL;
	}

	p->arecd = 10;
	if ((p->recds = (jc_recd *) calloc(p->arecd, sizeof(jc_recd))) == NULL) {
		if (pev != NULL) *pev = jc_malloc;
		p->del(p);
		return NULL;
	}

	if ((p->fname = strdup(fname)) == NULL) {
		if (pev != NULL) *pev = jc_malloc;
		p->del(p);
		return NULL;
	}

	p->modify = modify == jc_modify ? 1 : 0;
	p->create = create == jc_create ? 1 : 0;

	p->locate_key = jcnf_locate_key;
	p->get_key    = jcnf_get_key;
	p->set_key    = jcnf_set_key;
	p->add_key    = jcnf_add_key;
	p->delete_key = jcnf_delete_key;
	p->print_key  = jcnf_print_key;

	p->enable_modify = jcnf_enable_modify;
	p->update = jcnf_update;
	p->del = jcnf_del;

	if ((ev = jcnf_read(p)) != jc_ok) {
		if (ev != jc_noexisting) {
			if (pev != NULL) *pev = ev;
			p->del(p);
			return NULL;
		}
	}

	if (pev != NULL) *pev = jc_ok;

	return p;
}

/* ------------------------------- */
/* Return a pointer to the nth element of the key name. */
/* Return null if it is out of range or malloc failed. */
/* Free the returned value when done. */
char *jc_get_nth_elem(char *path, int n) {
	int i;
	char *p1, *p2, *rv;

	if (path == NULL)
		return NULL;

	p1 = path;
	if (*p1 == '/')
		p1++;

	for (i = 0; *p1 != '\000'; p1 = p2 + 1, i++) {
		if ((p2 = strchr(p1, '/')) == NULL)
			p2 = p1 + strlen(p1);
		if (i >= n) {
			if ((rv = malloc(p2 - p1 + 1)) == NULL)
				return NULL;
			strncpy(rv, p1, p2 - p1);
			rv[p2-p1] = '\000';
			return rv;
		}
		if (*p2 == '\000')
			break;
	}
	return NULL;
}


