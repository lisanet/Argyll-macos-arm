/* Test code that generates optimised C sorting */
/* code, suitable for classifying an input vector */
/* into a particular simplex. */

/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* Intended for inclusion in IMDI for 16 bit support */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/* Sorting data */
#define MAXSV 50		/* Maximum sort variables */
typedef struct {
	int bot;			/* 0 if we are doing ==, 1 if >, 2 if else */			
	int dval;			/* Difference value to apply */
	int vals[MAXSV];	/* Variable values */
	int sort[MAXSV];	/* Sorted index */
} sortd;

/* Context */
typedef struct {
	FILE *of;			/* Output file */
	int indt;			/* Indent */

	/* Sort data */
	int nv;				/* Number of variables to sort */
	int sl;				/* Non-zero to generate code for smallest to largest */
	int eq;				/* Non-zero to generate code for equal case */
	int la;				/* Non-zero to generate labels */
	int rl;				/* Recursion level */
	sortd sd[MAXSV];	/* Data at each recursion level */
	int ln;				/* Label number */

} fileo;

void line(fileo *f, char *fmt, ...);	/* Output one line */
void niline(fileo *f, char *fmt, ...);	/* Output one line, no indent */
void sline(fileo *f, char *fmt, ...);	/* Output start of line line */
void mline(fileo *f, char *fmt, ...);	/* Output middle of line */
void eline(fileo *f, char *fmt, ...);	/* Output end of line */
	
void cr(fileo *f) { line(f,""); }		/* Output a blank line */

void inc(fileo *f) { f->indt++; }		/* Increment the indent level */

void dec(fileo *f) { f->indt--; }		/* Decrement the indent level */

int
main(void) {
	fileo f[1];
	
	f->of = fopen("ttt.c","w");
	f->indt = 0;

	line(f,"/* This is the top */");
	cr(f);
	line(f,"#include <stdio.h>");
	cr(f);
	line(f,"void main(void) {");
	inc(f);
	line(f,"int v0, v1, v2, v3, v4;");
	line(f,"printf(\"Hi there!\\n\");");

	line(f,"v0 = 1;");
	line(f,"v1 = 9;");
	line(f,"v2 = 3;");
	line(f,"v3 = 2;");
	line(f,"v4 = 0;");

	/* ================ */
	/* Sort development */
	{
int xx;
		int i, j;

		f->nv = 6;			/* No variables */
		f->sl = 1;			/* Smallest to largest */
		f->eq = 0;			/* No equals case */
		f->la = 0;			/* No labels */
		f->rl = 0;			/* Top recursion level */
		f->ln = 0;			/* Start label number */
		f->sd[f->rl].bot = f->eq ? 0 : 1;
		f->sd[f->rl].dval = 256;
		for (i = 0; i < f->nv; i++) {
			f->sd[f->rl].vals[i] = -99999;	/* Initial sort value  = unknown */
			f->sd[f->rl].sort[i] = i;	/* Initial sort order */
		}
		f->sd[f->rl].vals[0] = 0;	/* Make sure one is known */

		/* First label is for no known sort */
		if (f->la)
			niline(f,"	label%d:;",f->ln++);

		/* Until we have done every sort condition */
		for (;f->rl >= 0;) {
			int ix, nix;		/* Variable index that needs comparison */
			int six, snix;		/* Sorted references to ix, nix */

//printf("\n~~recursion level %d\n",f->rl);
//for (i = 0; i < f->nv; i++) {
//xx = f->sd[f->rl].sort[i];
//printf("Current = %d (index %d)\n",f->sd[f->rl].vals[xx],xx);
//}
			/* see if we have resolved a sort */
			for (i = 1; i < f->nv; i++) {
				six = i-1;
				ix = f->sd[f->rl].sort[i-1];
				snix = i;
				nix = f->sd[f->rl].sort[i];
				if (f->sd[f->rl].vals[ix] == f->sd[f->rl].vals[nix]
				 || f->sd[f->rl].vals[nix] == -99999) {
					break;		/* Not distinguishable or not known */
				}
			}
//if (i < f->nv)
//printf("~~Unresolved sort indexes %d %d\n",ix,nix);
//else
//printf("~~Resolved sort\n");

			if (i < f->nv) {
				/* We haven't fully sorted the variables yet. */
				/* Compare ix to nix */
				if (f->sd[f->rl].bot == 0) {
//printf("~~ Unresolved at level %d, doing comparison for equals\n",f->rl);
					line(f,"if (v%d == v%d) {",ix,nix);		/* } */
					inc(f);
					f->sd[f->rl+1] = f->sd[f->rl];	/* Structure copy */
					f->sd[f->rl+1].bot = 0;
					f->sd[f->rl].bot = 1;
					f->rl++;
// ~~~~~9
					/* Find next lowest value from nix */
					for (i = snix+1; i < f->nv; i++) {
						int si = f->sd[f->rl].sort[i];
						if (f->sd[f->rl].vals[si] == -99999)
							continue;
						if (f->sd[f->rl].vals[si] < f->sd[f->rl].vals[nix]) {
							f->sd[f->rl].vals[nix] =
							              (f->sd[f->rl].vals[si] + f->sd[f->rl].vals[nix])/2;
							break;
						}
					}
					if (i >= f->nv) {		/* Not between next lowest */
						f->sd[f->rl].vals[nix] = f->sd[f->rl].vals[ix] -256;
					}
				} else if (f->sd[f->rl].bot == 1) {		/* { */
//printf("~~ Unresolved at level %d, doing comparison for greater\n",f->rl);
					if (f->eq) {
						dec(f);
						line(f,"} else if (v%d %c v%d) {",ix, f->sl ? '<':'>',nix);		/* } */
					} else
						line(f,"if (v%d %c v%d) {",ix,f->sl ? '<':'>', nix);		/* } */
					inc(f);
					if (f->la)
						niline(f,"	label%d:;",f->ln++);
					f->sd[f->rl+1] = f->sd[f->rl];	/* Structure copy */
					f->sd[f->rl+1].bot = f->eq ? 0 : 1;
					f->sd[f->rl].bot = 2;
					f->rl++;
					/* Find next lowest value from nix */
					for (i = snix+1; i < f->nv; i++) {
						int si = f->sd[f->rl].sort[i];
						if (f->sd[f->rl].vals[si] == -99999)
							continue;
						if (f->sd[f->rl].vals[si] < f->sd[f->rl].vals[nix]) {
							f->sd[f->rl].vals[nix] =
							              (f->sd[f->rl].vals[si] + f->sd[f->rl].vals[nix])/2;
							break;
						}
					}
					if (i >= f->nv) {		/* Not between next lowest */
						f->sd[f->rl].vals[nix] = f->sd[f->rl].vals[ix] -256;
					}
				} else if (f->sd[f->rl].bot == 2) {		/* { */
//printf("~~ Unresolved at level %d, doing else\n",f->rl);
					dec(f);
					line(f,"} else {	/* v%d %c v%d  */",ix,f->sl ? '>':'<', nix);	/* } */
					inc(f);
					if (f->la)
						niline(f,"	label%d:;",f->ln++);
					f->sd[f->rl+1] = f->sd[f->rl];	/* Structure copy */
					f->sd[f->rl+1].bot = f->eq ? 0 : 1;
					f->sd[f->rl].bot = 3;
					f->rl++;
					if (six-1 >= 0) {	/* If there is a higher value */
						int si = f->sd[f->rl].sort[six-1];
						/* Make equal to next highest, since never been compared */
						f->sd[f->rl].vals[nix] = f->sd[f->rl].vals[si];
					} else {	/* Nothing above */
						f->sd[f->rl].vals[nix] = f->sd[f->rl].vals[ix] + 256;
					}
					/* sort nix above ix */
					j = f->sd[f->rl].sort[six];
					f->sd[f->rl].sort[six] = f->sd[f->rl].sort[snix];
					f->sd[f->rl].sort[snix] = j;
				} else {	/* We've done both cases - return from recursion */
//printf("~~ Unresolved at level %d, going up a level\n",f->rl);
					dec(f);				/* { */
					line(f,"}");
					f->rl--;			/* Back a recursion level */
				}
			} else {
//printf("~~ Resolved at level %d, outputing sort code, up a level\n",f->rl);
				/* We have a resolved sort */
				/* So output the code for this sort combination */
				sline(f,"/* Sorted ");
				for (i = 0; i < f->nv; i++) {
					mline(f," %d",f->sd[f->rl].sort[i]);
				}
				eline(f," */");

#ifdef NEVER
				sline(f,"printf(\"Sorted = %%d %%d %%d %%d %%d\\n\"");
				for (i = 0; i < f->nv; i++) {
					mline(f,",v%d",f->sd[f->rl].sort[i]);
				}
				eline(f,");");
#else
				for (i = 0; i < f->nv; i++) {
					sline(f,"");
					mline(f,"op[%d] = v%d",i,f->sd[f->rl].sort[i]);
					eline(f,";");
				}
				line(f,"ip += %d;",f->nv);
				line(f,"op += %d;",f->nv);
				line(f,"continue;");
#endif
				f->rl--;			/* Back a recursion level */
			}
		}
	}
	/* ================ */
	dec(f);
	line(f,"}");

	fclose(f->of);

	return 0;
}


/* Output a line to the file (including trailing \n) */
void
line(fileo *f, char *fmt, ...)
{
	int i;
	va_list args;

	/* Indent to the correct level */
	for (i = 0; i < f->indt; i++)
		fprintf(f->of,"	");
	
	va_start(args, fmt);
	vfprintf(f->of, fmt, args);
	va_end(args);
	fprintf(f->of, "\n");
}

/* Output a line to the file (including trailing \n) */
/* No indent */
void
niline(fileo *f, char *fmt, ...)
{
	int i;
	va_list args;

	va_start(args, fmt);
	vfprintf(f->of, fmt, args);
	va_end(args);
	fprintf(f->of, "\n");
}

/* Output the start of a line to the file) */
void
sline(fileo *f, char *fmt, ...)
{
	int i;
	va_list args;

	/* Indent to the correct level */
	for (i = 0; i < f->indt; i++)
		fprintf(f->of,"	");
	
	va_start(args, fmt);
	vfprintf(f->of, fmt, args);
	va_end(args);
}

/* Output the middle of a line to the file) */
void
mline(fileo *f, char *fmt, ...)
{
	int i;
	va_list args;

	va_start(args, fmt);
	vfprintf(f->of, fmt, args);
	va_end(args);
}

/* Output the end of a line to the file (including trailing \n) */
void
eline(fileo *f, char *fmt, ...)
{
	int i;
	va_list args;

	va_start(args, fmt);
	vfprintf(f->of, fmt, args);
	va_end(args);
	fprintf(f->of, "\n");
}

