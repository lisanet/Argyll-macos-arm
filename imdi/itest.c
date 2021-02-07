
/* Verify and benchmark the imdi code */
/*
 * Copyright 2000 - 2006 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "imdi.h"
#include "refi.h"

/* Test parameters */
#undef TEST1					/* Test just one combination */
#define FULL					/* Test full range */
#undef VERBOSE					/* Report every conversion */
#undef REPORT_ERRORS			/* Report conversion with large errors */
#define TRAND					/* Test random in/out table contents */
#define TCRAND					/* Test random clut table contents */
#undef QUANTIZE					/* Quantize the target table values */
#define TBUFSIZE (2 * 1024 * 1024)	/* Default number of input bytes to test */
#define ITERS 10					/* Itterations */

double trans1(double in, double t);
double trans2(double in, double t);

/* Context for refi setup callbacks */
typedef struct {
	int id;
	int od;
	double icurve[MXDI];		/* Input curve factors */
	double clut[MXDO][MXDI];	/* clut matrix factors */
	double ocurve[MXDO];		/* Output curve factors */
} rcntx;

/* Input curve function */
static void input_curves(
	void *cntx,
	double *out_vals,
	double *in_vals
) {
	int i;
	rcntx *rx = (rcntx *)cntx;

	for (i = 0; i < rx->id; i++) {
		double val = in_vals[i];
#ifdef TRAND
		val = trans1(val, rx->icurve[i]);
#ifdef QUANTIZE
		val = ((int)(val * 255.0 + 0.5))/255.0;	/* Quantize to 8 bits */
#endif
#endif
		out_vals[i] = val;
	}
}

/* Multi-dim table function */
static void md_table(
	void *cntx,
	double *out_vals,
	double *in_vals
) {
	rcntx *rx = (rcntx *)cntx;
	int i, j;

#ifdef TCRAND
	for (j = 0; j < rx->od; j++) {
		double val = 0.0;
		for (i = 0; i < rx->id; i++) {
			val += rx->clut[j][i] * in_vals[i];
		}
#ifdef QUANTIZE
		val = ((int)(val * 255.0 + 0.5))/255.0;	/* Quantize to 8 bits */
#endif
		out_vals[j] = val;
	}
#else
	for (j = 0; j < rx->od; j++)
		out_vals[j] = in_vals[j % rx->id];
#endif
//printf("~1 out %f %f %f %f from %f\n", out_vals[0], out_vals[1], out_vals[2], out_vals[3], in_vals[0]);
}

/* Output curve function */
static void output_curves(
	void *cntx,
	double *out_vals,
	double *in_vals
) {
	int i;
	rcntx *rx = (rcntx *)cntx;

	for (i = 0; i < rx->od; i++) {
		double val = in_vals[i];
#ifdef TRAND
		val = trans1(val, rx->ocurve[i]);
#ifdef QUANTIZE
		val = ((int)(val * 255.0 + 0.5))/255.0;	/* Quantize to 8 bits */
#endif
#endif
		out_vals[i] = val;
	}
}

/* Complete reference interpolation */
void refi_interp(refi *r, double *out_vals, double *in_vals);

void usage(void) {
	fprintf(stderr,"Regression test imdi code Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"usage: itest [-q] [-s]\n");
	fprintf(stderr," -q            Quick test\n");
	fprintf(stderr," -s            Stop on error\n");
	fprintf(stderr," -r bits       Specify bits of randomness in input data\n");
	exit(1);
}


int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	int quick = 0;
	int stop = 0;
	int rbits = 16;
	int n, i, j, e;
	int pix, iix, oix;
	int ip, op, id, od;
	int ires, cres, ores;
	clock_t stime, ttime;	/* Start and total times */
	double xtime;			/* Total execution time in seconds */
	double npixels;			/* Number of pixels processed */
	rcntx rx;

	refi *r;
	double ribuf[MXDI];
	double robuf[MXDO];

	imdi *s;
	int iters;
	unsigned long tbufsize;
	unsigned char *ibuf;
	unsigned char *obuf;
	unsigned char *inp[MXDI];
	unsigned char *outp[MXDO];
	unsigned short *ibuf2;			/* 16 bit references */
	unsigned short *obuf2;

	double omxerr = 0;		/* Overall max error /1.0 */

	/* Define combinations to test */
#ifdef TEST1
#pragma message("!!!!!!!!!!!!!!!!! TEST1 is defined !!!!!!!!!!!!!!!1")
	int ids[] = { 3, 0 };			/* Input dimensions */
	int ods[] = { 3, 0 };			/* Output dimensions */
	int iprs[] = { 16, 0 };
	int oprs[] = { 16, 0 };
#else
#ifndef FULL
	int ids[] = { 1, 3, 4, 8, 0 };
	int ods[] = { 1, 3, 4, 8, 0 };
	int iprs[] = { 8, 8,  16, 0};
	int oprs[] = { 8, 16, 16, 0};
#else
	int ids[] = { 1, 3, 4, 5, 6, 7, /* 8, */ 0 };
	int ods[] = { 1, 3, 4, 5, 6, 7, /* 8, */ 0 };
	int iprs[] = { 8, 8,  16, 0};
	int oprs[] = { 8, 16, 16, 0};
#endif
#endif

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')	{	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage();

			/* Quick test */
			else if (argv[fa][1] == 'q' || argv[fa][1] == 'Q') {
				quick = 1;
			}

			/* Stop on error */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				stop = 1;
			}

			/* Degree of data randomness */
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				fa = nfa;
				if (na == NULL) usage();
				rbits = atoi(na);
				if (rbits < 0)
					rbits = 1;
				else if (rbits > 16)
					rbits = 16;
			}
			else 
				usage();
		} else
			break;
	}

	for (iix = 0; ids[iix] > 0; iix++) {			/* All input dimensions */
		for (oix = 0; ods[oix] > 0; oix++) {		/* All output dimensions */
			for (pix = 0; iprs[pix] > 0; pix++) {	/* All precisions */
				int rmask = 0;
				ip = iprs[pix];		/* Input precision */
				op = oprs[pix];		/* Output precision */
				id = ids[iix];		/* Input dimensions */
				od = ods[oix];		/* Outpu dtimensions */

				if (ip != 8 && ip != 16)
					error("Can't handle input precision of %d in testing\n",ip);
				if (op != 8 && op != 16)
					error("Can't handle output precision of %d in testing\n",op);
				if (id > MXDI)
					error("Can't handle id greater than %d in testing\n",MXDI);
				if (od > MXDO)
					error("Can't handle od greater than %d in testing\n",MXDO);

				printf("Testing id = %d, od = %d, ip = %d, op = %d\n",id,od,ip,op);

				ires = 256;

				switch (id) {
					case 1:
						cres = 257;
						break;
					case 2:
					case 3:
						cres = 33;
						break;
					case 4:
						cres = 17;
						break;
					case 5:
						cres = 7;
						break;
					case 6:
						cres = 7;
						break;
					case 7:
						cres = 5;
						break;
					case 8:
						cres = 5;
						break;
					default:
						cres = 3;
				}
				ores = 256;

				/* Setup error messages */
				error_program = "itest";

				/* Create the reference first */
				printf("About to create refi\n");

				/* Create test curves */
				rx.id = id;
				rx.od = od;

				rand32(0x1234);
				for (i = 0; i < id; i++) {
					rx.icurve[i] = d_rand(-1.0, 1.0);
				}

//printf("Matrix params =\n");
				rand32(0x2345);
				for (j = 0; j < od; j++) {
					double s = 0.0;
					for (i = 0; i < id; i++) {
						double v = d_rand(1.0, 0.0);
						rx.clut[j][i] = v;			/* Make them sum to 1.0 */
						s += v;
					}
					for (i = 0; i < id; i++) {
							rx.clut[j][i] /= s;
//printf(" %f",rx.clut[j][i]);
					}
//printf("\n");
				}

				rand32(0x3456);
				for (j = 0; j < od; j++) {
					rx.ocurve[j] = d_rand(-1.0, 1.0);
				}

				r = new_refi(
					id,				/* Number of input dimensions */
					od,				/* Number of output dimensions */
					ires,			/* Input table resolution */
					cres,			/* clut resolution */
					ores,			/* Output table resolution */
					input_curves,	/* Callback functions */
					md_table,
					output_curves,
					(void *)&rx		/* Context to callbacks */
				);

				if (r == NULL) {
					error("new_refi failed");
				}

				/* Now crate the imdi we want to test, using the */
				/* reference as a template */
				printf("About to create imdi\n");

				s = new_imdi(
					id,				/* Number of input dimensions */
					od,				/* Number of output dimensions */
					ip == 8 ? pixint8 : pixint16,	/* Input pixel representation */
					0x0,			/* Treat every channel as unsigned */
					NULL,			/* No raster to callback mapping */
					prec_min,		/* Minimum of input and output precision */
					op == 8 ? pixint8 : pixint16,	/* Output pixel representation */
					0x0,			/* Treat every channel as unsigned */
					NULL,			/* No raster to callback mapping */
					cres,			/* Desired table resolution */
					oopts_none,		/* Desired per channel output options */
					NULL,			/* Output channel check values */
					opts_none,		/* Desired processing direction and stride support */
					refi_input,		/* Callback functions */
					refi_clut,
					refi_output,
					(void *)r		/* Context to callbacks */
				);

				if (s == NULL) {
					error("new_imdi failed");
				}

				if (quick) {
					iters = 1;
					tbufsize = 4096/id;
				} else {
					iters = ITERS;
					tbufsize = TBUFSIZE/id;
				}

				/* Allocate the test buffers */
				if ((ibuf = malloc(sizeof(unsigned char) * ip/8 * id * tbufsize)) == NULL)
					error("Malloc of input buffer failed");
				ibuf2 = (unsigned short *)ibuf;

				if ((obuf = malloc(sizeof(unsigned char) * op/8 * od * tbufsize)) == NULL)
					error("Malloc of output buffer failed");
				obuf2 = (unsigned short *)obuf;

				/* Initialise the input buffer contents */
				rand32(0x12345678);
				if (ip == 8) {
					unsigned ui;
					int rr = rbits;
					if (rr > 8)
						rr = 8;
					rmask = ((1 << rr) -1) << (8-rr);
					for (ui = 0; ui < tbufsize; ui += id) {
						for (e = 0; e < id; e++) {
							unsigned long ran = rand32(0);
							ibuf[ui + e] = (unsigned char)(ran & rmask);
						}
					}
				} else {
					unsigned ui;
					rmask = ((1 << rbits) -1) << (16-rbits);
					for (ui = 0; ui < tbufsize; ui += id) {
						for (e = 0; e < id; e++) {
							unsigned long ran = rand32(0);
							ibuf2[ui + e] = (unsigned short)(ran & rmask);
						}
					}
				}

				/* We are assuming packed pixel interleaved */
				inp[0]  = ibuf;
				outp[0] = obuf;

				/* Benchmark it */
				stime = clock();
				for (n = 0; n < iters; n++) {
					s->interp(s, (void **)outp, 0, (void **)inp, 0, tbufsize);
				}
				ttime = clock() - stime;
				xtime = (double)ttime/(double)CLOCKS_PER_SEC;
				npixels = (double)iters * (double)tbufsize;
				
				if (xtime > 0.0)
					printf("Speed = rate = %f Mpix/sec\n",1e-6 * npixels / xtime);
				else
					printf("Speed - too fast!\n");

				{
					unsigned ui;
					double mxerr = 0.0;
					double avgerr = 0.0;

					/* Verify the accuracy against refi of each sample */
					for (ui = j = 0; ui < tbufsize; ui += id, j += od) {
						int mxserr;

						if (ip == 8) {
							for (e = 0; e < id; e++) {
								ribuf[e] = ibuf[ui + e]/255.0;
							}
						} else {
							for (e = 0; e < id; e++) {
								ribuf[e] = ibuf2[ui + e]/65535.0;
							}
						}
						refi_interp(r, robuf, ribuf);
				
						mxserr = 0;

						if (op == 8) {
							for (e = 0; e < od; e++) {
								double err = robuf[e] * 255.0 - obuf[j + e];
								if (err < 0)
									err = -err;
								if (err > mxerr)
									mxerr = err;
								if (err > mxserr)
									mxserr = err;
								avgerr += err;
							}
						} else {
							for (e = 0; e < od; e++) {
								double err = robuf[e] * 65535.0 - obuf2[j + e];
								if (err < 0)
									err = -err;
								if (err > mxerr)
									mxerr = err;
								if (err > mxserr)
									mxserr = err;
								avgerr += err;
							}
						}
#if defined(REPORT_ERRORS) || defined(VERBOSE)
						if (mxserr >= 37) {
							if (ip == 8) {
								if (id == 1)
									printf("in %d, ", ibuf[ui+0]);
								if (id == 2)
									printf("in %d %d, ", ibuf[ui+0], ibuf[ui+1]);
								if (id == 3)
									printf("in %d %d %d, ", ibuf[ui+0], ibuf[ui+1], ibuf[ui+2]);
								if (id == 4)
									printf("in %d %d %d, ", 
											ibuf[ui+0], ibuf[ui+1], ibuf[ui+2], ibuf[ui+3]);
							} else {
								if (id == 1)
									printf("in %d, ", ibuf2[ui+0]);
								if (id == 2)
									printf("in %d %d, ", ibuf2[ui+0], ibuf2[ui+1]);
								if (id == 3)
									printf("~in %d %d %d, ", ibuf2[ui+0], ibuf2[ui+1], ibuf2[ui+2]);
								if (id == 4)
									printf("in %d %d %d, ",
											ibuf2[ui+0], ibuf2[ui+1], ibuf2[ui+2], ibuf2[ui+3]);
							}
							if (ip == 8) {
								if (od == 1)
									printf("is %d, should be %d\n",
											obuf[j+0],
											(int)(robuf[0] * 255.0 + 0.5));
								if (od == 2)
									printf("is %d %d, should be %d %d\n",
											obuf[j+0], obuf[j+1],
											(int)(robuf[0] * 255.0 + 0.5),
											(int)(robuf[1] * 255.0 + 0.5));
								if (od == 3)
									printf("is %d %d %d, should be %d %d %d\n",
											obuf[j+0], obuf[j+1], obuf[j+2],
											(int)(robuf[0] * 255.0 + 0.5),
											(int)(robuf[1] * 255.0 + 0.5),
											(int)(robuf[2] * 255.0 + 0.5));
								if (od == 4)
									printf("is %d %d %d %d, should be %d %d %d %d\n",
											obuf[j+0], obuf[j+1], obuf[j+2], obuf[j+3],
											(int)(robuf[0] * 255.0 + 0.5),
											(int)(robuf[1] * 255.0 + 0.5),
											(int)(robuf[2] * 255.0 + 0.5),
											(int)(robuf[3] * 255.0 + 0.5));
							} else {
								if (od == 1)
									printf("is %d, should be %d\n",
											obuf2[j+0],
											(int)(robuf[0] * 65535.0 + 0.5));
								if (od == 2)
									printf("is %d %d, should be %d %d\n",
											obuf2[j+0], obuf2[j+1],
											(int)(robuf[0] * 65535.0 + 0.5),
											(int)(robuf[1] * 65535.0 + 0.5));
								if (od == 3)
									printf("is %d %d %d, should be %d %d %d\n",
											obuf2[j+0], obuf2[j+1], obuf2[j+2],
											(int)(robuf[0] * 65535.0 + 0.5),
											(int)(robuf[1] * 65535.0 + 0.5),
											(int)(robuf[2] * 65535.0 + 0.5));
								if (od == 4)
									printf("is %d %d %d %d, should be %d %d %d %d\n",
											obuf2[j+0], obuf2[j+1], obuf2[j+2], obuf2[j+3],
											(int)(robuf[0] * 65535.0 + 0.5),
											(int)(robuf[1] * 65535.0 + 0.5),
											(int)(robuf[2] * 65535.0 + 0.5),
											(int)(robuf[3] * 65535.0 + 0.5));
							}
						}
						
#endif /* VERBOSE || REPORT+ERRORS */
					}
					avgerr /= ((double)tbufsize * od);

					{
						double fmxerr;	/* Relative maximum error */
						double favgerr;	/* Relative average error */
					
						/* Always relative to basic precision */
						fmxerr = mxerr / ((1 << op) - 1.0);
						favgerr = avgerr / ((1 << op) - 1.0);

						printf("Worst error = %f = %f%%, average error = %f%%\n",
						       mxerr, 100.0 * fmxerr, 100.0 * favgerr);
						printf("\n");
						
						if (fmxerr > omxerr)
							omxerr = fmxerr;

					}
				}
				/* Free everything up */
				free(ibuf);
				free(obuf);
				refi_free(r);
				s->del(s);

				if (stop &&	(100.0 * omxerr) > 1.2) {
					goto quit;
				} 
			}
		}
	}

 quit:;
	printf("Overall worst error = %f%%\n", 100.0 * omxerr);

	return 0;
}




/* ------------------------------------------------- */
/* Support */

/* Run an interpolation though the whole refi */
void refi_interp(
refi *r,
double *out_vals,
double *in_vals
) {
#ifdef QUANTIZE
	int e;
#endif
	double ivals[MXDI];
	double ovals[MXDO];

	refi_input((void *)r, ivals, in_vals);
#ifdef QUANTIZE
	for (e = 0; e < r->id; e++)
		ivals[e] = ((int)(ivals[e] * 255.0 + 0.5))/255.0;	/* Quantize to 8 bits */
#endif

	refi_clut((void *)r, ovals, ivals);

#ifdef QUANTIZE
	for (e = 0; e < r->od; e++)
		ovals[e] = ((int)(ovals[e] * 255.0 + 0.5))/255.0;	/* Quantize to 8 bits */
#endif
	refi_output((void *)r, out_vals, ovals);
}


/* ------------------------------------------------- */
/* Some test functions */

/* Single bump version */
double
trans1(
double in,
double t			/* Non-linearity factor. 0.0 to +/-1.0 */
) {
	double out;

	out = in + t * in * (1.0 - in);
	return out;
}

/* Double bump version */
double
trans2(
double in,
double t			/* Non-linearity factor. 0.0 to +/-1.0 */
) {
	double nf, out;

	if (in >= 0.5) {
		nf = -2.0 * (in - 0.5) * (1.0 - in);
	} else {
		nf = 2.0 * in * (0.5 - in);
	}
	out = in + t * nf;
	return out;
}


















