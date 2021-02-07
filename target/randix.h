
#ifndef RANDIX_H
/* 
 * Argyll Color Correction System
 *
 * Random array indexing class
 *
 * Author: Graeme W. Gill
 * Date:   16/10/96
 *
 * Copyright 1996, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

struct _randix {
/* private: */
	int tbit;		/* Top bit mask */
	int mask;		/* Overall mask */
	int xorm;		/* Xor value */
	int length;	/* Length needed */
	int ss;		/* Current value */

/* public: */
	/* return the next in the sequence */
	int (*next)(struct _randix *p);

	/* Destroy ourselves */
	void (*del)(struct _randix *p);

	}; typedef struct _randix randix;

/* Creator */
/* Counts withing range 0 to length-1 */
extern randix *new_randix(int length, int start);


#define RANDIX_H
#endif /* RANDIX_H */
