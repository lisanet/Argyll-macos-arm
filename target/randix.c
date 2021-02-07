
/* 
 * Argyll Color Correction System
 *
 * Random index routines.
 *
 * Author: Graeme W. Gill
 * Date:   5/10/96
 *
 * Copyright 1996, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include "randix.h"

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

/* A table entry structure */
typedef struct {
	int bits;		/* Number of bits */
	int length;		/* Length of sequence */
	int xorm;		/* Xor mask */
} tabe;

		
static tabe table[] = {
	{ 1,	(1<<1)-1,	0x1},
	{ 2,	(1<<2)-1,	0x3},
	{ 3,	(1<<3)-1,	0x3},
	{ 4,	(1<<4)-1,	0x13},
	{ 5,	(1<<5)-1,	0x1b},
	{ 6,	(1<<6)-1,	0x1b},
	{ 7,	(1<<7)-1,	0x65},
	{ 8,	(1<<8)-1,	0xc3},
	{ 9,	(1<<9)-1,	0x1b5},
	{ 10,	(1<<10)-1,	0x1c7},
	{ 11,	(1<<11)-1,	0x6bb},
	{ 12,	(1<<12)-1,	0x5c5},
	{ 13,	(1<<13)-1,	0x15b9},
	{ 14,	(1<<14)-1,	0x36d1},
	{ 15,	(1<<15)-1,	0x376b},
	{ 16,	(1<<16)-1,	0x5bab},
	{ 17,	(1<<17)-1,	0x1b5c5},
	{ 18,	(1<<18)-1,	0x15561},
	{ 19,	(1<<19)-1,	0x5b4db},
	{ 20,	(1<<20)-1,	0xf6e01},
	{ 21,	(1<<21)-1,	0x186517},
	{ 22,	(1<<22)-1,	0x6bf4f},
	{ 23,	(1<<23)-1,	0x376623},
	{ 24,	(1<<24)-1,	0xf54e35},
	{ 0,	0,	0}
};

#define PSRAND(S, XOV, TBIT, MASK) ((((S) & TBIT) ? (((S) << 1) ^ (XOV)) : ((S) << 1)) & MASK)

static int
randix_next(struct _randix *p) {
	int rv = p->ss-1;	/* Return start value first */
	do {
		p->ss = PSRAND(p->ss, p->xorm, p->tbit, p->mask);
	} while(p->ss >= p->length);
	return rv;
}


/* Destroy ourselves */
static void
randix_del(randix *p) {
	free (p);
}

/* Creator */
randix *new_randix(int length, int start)
{
	int i;
	randix *p;
	if ((p = (randix *)malloc(sizeof(randix))) == NULL)
		error ("randix: malloc failed");

	p->next = randix_next;
	p->del  = randix_del;

	if (length == 0)
		error ("randix: Can't handle length %d",length);

	start %= length;

	p->length = length+1;
	for (i = 0; table[i].bits != 0; i++) {
		if (length <= table[i].length) {
			p->tbit = 1 << (table[i].bits-1);
			p->mask = (p->tbit << 1)-1;
			p->xorm = table[i].xorm;
			p->ss = start + 1;
			return p;
			break;
		}
	}
	error ("randix: Can't handle length %d",length);
	return NULL;
}





#ifdef NEVER
main(argc,argv)
int argc;
char *argv[];
{
	int length = 11;
	int i, iv, cv;
	randix *r;

	if (argc > 1)
		length = atoi(argv[1]);

	r = new_randix(length, 0);

	iv = r->next(r);
	printf("First Val = %d\n",iv);
	i = 0;
	do {
		cv = r->next(r);
		i++;
	} while (cv != iv);
	printf("Last Val = %d\n",cv);
	printf("Count for length %d = %d\n",length,i);
	return 0;
}

/* Basic printf type error() and warning() routines */
void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"targen: Error - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	fflush(stdout);
	exit (-1);
}
#endif /* NEVER */













