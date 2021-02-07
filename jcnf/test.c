

/* Test out the jcnf object */

/* Reads the file test.jcnf, and copies it to testo.jcnf */

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


#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "jcnf.h"

#undef TEST_PATH_PARSING
#define TEST_ENABLE_MODIFY

void error(char *fmt, ...);

int
main(int argc, char *argv[]) {
    int fa,nfa;             /* argument we're looking at */

	jcnf *jci, *jco;
	jc_error ev;
	int ix;

	printf("Hi there\n");

#ifdef TEST_PATH_PARSING
	{
		char *str;
		int i, j;

		for (i = 0; i < 4; i++) {
			if (i == 0) 
				str = "this/is/a/test";
			else if (i == 1) 
				str = "/this/is/another//test/";
			else if (i == 2) 
				str = "/";
			else 
				str = "";

			printf("String = '%s'\n",str);
			for (j = 0; j < 10; j++) {
				char *s;

				s = jc_get_nth_elem(str, j);
				if (s == NULL)
					break;
				printf("%d th element = '%s'\n",j,s);
				free(s);
			}
		}
		return 0;
	}
#endif /* TEST_PATH_PARSING */

	if ((jci = new_jcnf(&ev, "test.jcnf",jc_read,jc_no_create)) == NULL) {
		error("new_jcnf '%s' failed with error %d\n","test.jcnf",ev);
	}

#ifdef TEST_ENABLE_MODIFY
	if ((jco = new_jcnf(&ev, "testo.jcnf",jc_read,jc_no_create)) == NULL) {
		jci->del(jci);
		error("new_jcnf '%s' failed with error %d\n","testo.jcnf",ev);
	}
	printf("Read output file, hit return to continue\n");
	getchar();
	if ((ev = jco->enable_modify(jco)) != jc_ok) {
		jci->del(jci);
		error("enable_modify '%s' failed with error %d\n","testo.jcnf",ev);
	}
#else
	if ((jco = new_jcnf(&ev, "testo.jcnf",jc_modify,jc_create)) == NULL) {
		jci->del(jci);
		error("new_jcnf '%s' failed with error %d\n","testo.jcnf",ev);
	}
#endif

	/* Delete everything from the second file */
	for (ix = jco->nkeys-1; ix >= 0; ix--) {
		if ((ev = jco->delete_key(jco, ix, NULL)) != jc_ok) {
			jco->del(jco);
			jci->del(jci);
			error("deleting keys from '%s' failed with error %d\n","testo.jcnf",ev);
		}
	}

	for (ix = 0; ; ix++) {
		char *key;
		jc_type type;
		unsigned char *data;
		size_t dataSize;
		char *comment;

		if ((ev = jci->get_key(jci, ix, &key, &type, &data, &dataSize, &comment)) != jc_ok) {
			if (ev != jc_ix_oorange) {
				jco->del(jco);
				jci->del(jci);
				error("get_key failed with %d\n",ev);
			}
			break;
		}
		jci->print_key(jci,ix, stderr);
		if ((ev = jci->add_key(jco, key, type, data, dataSize, comment)) != jc_ok) {
			jco->del(jco);
			jci->del(jci);
			error("add_key failed with %d\n",ev);
		}
	}
	if ((ev = jco->update(jco)) != 0) {
		jco->del(jco);
		jci->del(jci);
		error("jcnf write failed with error %d\n",ev);
	}

	jci->del(jci);
	jco->del(jco);

	printf("We're done\n");
	return 0;
}

void
error(char *fmt, ...)
{
    va_list args;

    fprintf(stderr,"test: Error - ");
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    exit (-1);
}
