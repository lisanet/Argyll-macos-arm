
			

/* 
 * Argyll Color Correction System
 * Print Device calibration curve generator.
 *
 * Author: Graeme W. Gill
 * Date:   2008/3/3
 *
 * Copyright 1996-2008 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in the colorent wedge test chart
 * points, and creates a set of per channel correction curves.
 */

/*
 * TTBD:
 *		Allow auto max threshold to be scaled on command line ?
 *		ie. -m# set % to go below the default optimal maximum.
 */


/*
	Additive spaces are handled by inverting the device values internally.
	(Such a space should probably have ICX_INVERTED set as well, indicating
	 that the underlying device is actually subtractive.) 
 */

#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <time.h>
#include "copyright.h"
#include "aconfig.h"
#include "cgats.h"
#include "numlib.h"
#include "sort.h"
#include "rspl.h"
#include "xicc.h"
#include "plot.h"
#include "ui.h"


#define RSPLFLAGS (0 /* | RSPL_2PASSSMTH | RSPL_EXTRAFIT2 */)

#define RSPLSMOOTH 2.0	/* RSPL Smoothness factor use on measured device points */

#define TCURVESMOOTH 1.0 /* RSPL smoothness factor for target aim points */

#define GRES 256		/* Rspl grid resolution */
#define SLOPE_NORM 70.0	/* Normalized delta E for below thresholds */
#define MIN_SLOPE_A 8.0	/* Criteria for Auto max, DE/dDev at max */
#define MIN_SLOPE_O 3.0	/* Criteria for Auto max, min DE/dDev below max */

#define CAL_RES 256		/* Resolution saved to .cal file */

#define PRES 256		/* Plotting resolution */

void usage(char *diag, ...) {
	int i;
	fprintf(stderr,"Create printer calibration, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: %s [-options] [prevcal] inoutfile\n",error_program);
	fprintf(stderr," -v verbosity    Verbose mode\n");
	fprintf(stderr," -p              Plot graphs.\n");
	fprintf(stderr," -i              Initial calibration, set targets, create .cal\n");
	fprintf(stderr," -r              Re-calibrate against previous .cal and create new .cal\n");
	fprintf(stderr," -e              Verify against previous .cal\n");
	fprintf(stderr," -I              Create imitation target from .ti3 and null calibration\n");
	fprintf(stderr," -d              Go through the motions but don't write any files\n");
	fprintf(stderr," -s smoothing    Extra curve smoothing (default 1.0)\n");
	fprintf(stderr," -A manufacturer Set the manufacturer description string\n");
	fprintf(stderr," -M model        Set the model description string\n");
	fprintf(stderr," -D description  Set the profile Description string\n");
	fprintf(stderr," -C copyright    Set the copyright string\n");
	fprintf(stderr," -x# percent     Set initial maximum device %% target (override auto)\n");
	fprintf(stderr," -m# percent     Set initial dev target to %% of auto maximum\n");
	fprintf(stderr," -n# deltaE      Set initial white minimum deltaE target\n");
	fprintf(stderr," -t# percent     Set initial 50%% transfer curve percentage target\n");
	fprintf(stderr,"   # = c, r, 0	 First channel\n");
	fprintf(stderr,"       m, g, 1	 Second channel\n");
	fprintf(stderr,"       y, b, 2	 Third channel\n");
	fprintf(stderr,"       k,    3	 Fourth channel, etc.\n");
	fprintf(stderr," -a              Create an Adobe Photoshop .AMP file as well as a .cal\n");
	fprintf(stderr," prevcal         Base name of previous .cal file for recal or verify.\n");
	fprintf(stderr," inoutname       Base name of input .ti3 file, output .cal file\n");
	exit(1);
}

/* - - - - - - - - - - - - - - - - - - - - - - - */
typedef struct {
	double loc;					/* Location up the curve 0.0 - 1.0 */
	double val[MAX_CHAN];		/* Value at that location 0.0 - 1.0 */
} trans_point;

/* Class to hold a print calibration target */
struct _pcaltarg {
	inkmask devmask;            	/* ICX ink mask of device space */
	
	/* Note that with all of these, a value < 0.0 */
	/* indicates no value set. */
	int devmaxset;						/* Flag - nz if the devmax is set */
	double devmax[MAX_CHAN];			/* Device value maximum 0.0 - 1.0 */

	int ademaxset;						/* Flag - nz if the ademax is set */
	double ademax[MAX_CHAN];			/* abs DE maximum for each channel */

	int ademinset;						/* Flag - nz if the ademin is set */
	double ademin[MAX_CHAN];			/* abs DE minimum for each channel */

	int no_tpoints;						/* Number of transfer curve points */
	trans_point *tpoints;				/* Array of transfer curve points */

	char err[500];						/* Error message from diagnostics */

	/* Methods */
	void (*del)(struct _pcaltarg *p);

	/* Save/restore to a CGATS file */
	int (*write)(struct _pcaltarg *p, cgats *cg, int tab);	/* return nz on error */
	int (*read)(struct _pcaltarg *p, cgats *cg, int tab);	/* return nz on error */

	/* Set values in the target */
	void (*update_devmax)(struct _pcaltarg *p, int chan, double val);
	void (*update_ademax)(struct _pcaltarg *p, int chan, double val);
	void (*update_ademin)(struct _pcaltarg *p, int chan, double val);
	void (*update_tcurve)(struct _pcaltarg *p, int chan, double loc, double val);

	/* Reurn nz if the target has been set */
	int (*is_set)(struct _pcaltarg *p);

	/* Update settings or from one from another */
	void (*update)(struct _pcaltarg *p, struct _pcaltarg *s);

}; typedef struct _pcaltarg pcaltarg; 

static void pcaltarg_del(pcaltarg *p) {
	if (p != NULL) {
		free(p);
	}
}

/* Write the cal target to a givent cgats table */
static int pcaltarg_write(pcaltarg *p, cgats *cg, int tab) {
	int i, j;
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	char *ident = icx_inkmask2char(p->devmask, 1); 
	char *bident = icx_inkmask2char(p->devmask, 0); 
	int devchan = icx_noofinks(p->devmask);
	int nsetel = 0;
	cgats_set_elem *setel;	/* Array of set value elements */
	char buf[100];

	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */

	/* Setup output cgats file */
	cg->add_table(cg, tt_other, 0);   /* Add a table for Calibration TarGet values */
	cg->add_kword(cg, tab, "DESCRIPTOR", "Argyll Calibration Target Definition File",NULL);
	cg->add_kword(cg, tab, "ORIGINATOR", "Argyll printcal", NULL);
	cg->add_kword(cg, tab, "CREATED",atm, NULL);
	cg->add_kword(cg, tab, "COLOR_REP", ident, NULL);

	/* Setup the table, which holds all the model parameters. */
	/* There is always a parameter per X Y Z or spectral band */
	cg->add_field(cg, tab, "PARAMTYPE", nqcs_t);
	nsetel++;
	sprintf(buf, "%s_I",bident);
	cg->add_field(cg, tab, buf, r_t);
	nsetel++;
	for (i = 0; i < devchan; i++) {
		inkmask imask = icx_index2ink(p->devmask, i);
		sprintf(buf, "%s_%s",bident,icx_ink2char(imask));
		cg->add_field(cg, tab, buf, r_t);
		nsetel++;
	}

	if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * nsetel)) == NULL) {
		free(ident);
		free(bident);
		sprintf(p->err,"ctg_write: malloc of setel failed");
		return 1;
	}

	/* Write out the values */
	if (p->devmaxset) {
		/* This is informational only */
		setel[0].c = "DEVMAX_USED";
		setel[1].d = 0.0;		/* Not used */
			
		if (p->devmask & ICX_ADDITIVE) {
			for (i = 0; i < devchan; i++)
				setel[2+i].d = 1.0 - p->devmax[i];
		} else {
			for (i = 0; i < devchan; i++)
				setel[2+i].d = p->devmax[i];
		}
		cg->add_setarr(cg, tab, setel);
	}
	if (p->ademaxset) {
		setel[0].c = "DELMAX_AIM";
		setel[1].d = 0.0;		/* Not used */
			
		for (i = 0; i < devchan; i++)
			setel[2+i].d = p->ademax[i];
		cg->add_setarr(cg, tab, setel);
	}
	if (p->ademinset) {
		setel[0].c = "DELMIN_AIM";
		setel[1].d = 0.0;		/* Not used */
			
		for (i = 0; i < devchan; i++)
			setel[2+i].d = p->ademin[i];
		cg->add_setarr(cg, tab, setel);
	}
	for (j = 0; j < p->no_tpoints; j++) {
		setel[0].c = "TRANS_PNT";
		setel[1].d = p->tpoints[j].loc;

		for (i = 0; i < devchan; i++)
			setel[2+i].d = p->tpoints[j].val[i];
		cg->add_setarr(cg, tab, setel);
	}
	free(setel);
	free(ident);
	free(bident);

	return 0;
}

/* Read the cal target from a given cgats table */
static int pcaltarg_read(pcaltarg *p, cgats *cg, int tab) {
	char *bident;
	int devchan;
	int i, j, ix;
	int ti;					/* Temporary CGATs index */
	int spi[2+MAX_CHAN];	/* CGATS indexes for each field */
	char buf[100];

	if ((ti = cg->find_kword(cg, tab, "COLOR_REP")) < 0) {
		sprintf(p->err, "ctg_read: Can't fint COLOR_REP");
		return 1;
	}

	if ((p->devmask = icx_char2inkmask(cg->t[tab].kdata[ti]) ) == 0) {
		sprintf(p->err, "ctg_read: unrecognized COLOR_REP '%s'",cg->t[tab].kdata[ti]);
		return 1;
	}
	devchan = icx_noofinks(p->devmask);
	bident = icx_inkmask2char(p->devmask, 0); 

	/* Figure out the indexes of all the fields */
	if ((spi[0] = cg->find_field(cg, tab, "PARAMTYPE")) < 0) {
		sprintf(p->err, "ctg_read: Can't find field PARAMTYPE");
		free(bident);
		return 1;
	}
	sprintf(buf, "%s_I",bident);
	if ((spi[1] = cg->find_field(cg, tab, buf)) < 0) {
		sprintf(p->err, "ctg_read: Can't find field %s",buf);
		free(bident);
		return 1;
	}
	for (i = 0; i < devchan; i++) {
		inkmask imask = icx_index2ink(p->devmask, i);
		sprintf(buf, "%s_%s",bident,icx_ink2char(imask));
		if ((spi[2+i] = cg->find_field(cg, tab, buf)) < 0) {
			sprintf(p->err, "ctg_read: Can't find field %s",buf);
			free(bident);
			return 1;
		}
	}

	/* Go through all the entries in the table, putting them in the right place */
	for (ix = 0; ix < cg->t[tab].nsets; ix++) {

		if (strcmp((char *)cg->t[tab].fdata[ix][spi[0]], "DELMAX_AIM") == 0) { 
			for (i = 0; i < devchan; i++)
				p->ademax[i] = *((double *)cg->t[tab].fdata[ix][spi[2+i]]);
			p->ademaxset = 1;

		} else if (strcmp((char *)cg->t[tab].fdata[ix][spi[0]], "DELMIN_AIM") == 0) { 
			for (i = 0; i < devchan; i++)
				p->ademin[i] = *((double *)cg->t[tab].fdata[ix][spi[2+i]]);
			p->ademinset = 1;

		} else if (strcmp((char *)cg->t[tab].fdata[ix][spi[0]], "TRANS_PNT") == 0) { 
			if ((p->tpoints = (trans_point *)realloc(p->tpoints, sizeof(trans_point)
			                                                   * (p->no_tpoints+1))) == NULL)
				error("Realloc of tpoints");
			p->tpoints[p->no_tpoints].loc = *((double *)cg->t[tab].fdata[ix][spi[1]]);
			for (i = 0; i < devchan; i++)
				p->tpoints[p->no_tpoints].val[i] = *((double *)cg->t[tab].fdata[ix][spi[2+i]]);
			p->no_tpoints++;
		}
	}
	free(bident);

	return 0;
}

/* Update an individual setting. Use chan < 0 to set all to default */
void pcaltarg_update_devmax(struct _pcaltarg *p, int chan, double val) {
	int i;
	if (p->devmaxset == 0) {
		for (i = 0; i < MAX_CHAN; i++)
			p->devmax[i] = -1.0;
		p->devmaxset = 1;
	}
	if (chan >= 0)
		p->devmax[chan] = val;
}
void pcaltarg_update_ademax(struct _pcaltarg *p, int chan, double val) {
	int i;
	if (p->ademaxset == 0) {
		for (i = 0; i < MAX_CHAN; i++)
			p->ademax[i] = -1.0;
		p->ademaxset = 1;
	}
	if (chan >= 0)
		p->ademax[chan] = val;
}
void pcaltarg_update_ademin(struct _pcaltarg *p, int chan, double val) {
	int i;
	if (p->ademinset == 0) {
		for (i = 0; i < MAX_CHAN; i++)
			p->ademin[i] = -1.0;
		p->ademinset = 1;
	}
	if (chan >= 0)
		p->ademin[chan] = val;
}
void pcaltarg_update_tcurve(struct _pcaltarg *p, int chan, double loc, double val) {
	int i, j;

	/* See if a transfer curve point already exists */
	for (j = 0; j < p->no_tpoints; j++) {
		if (p->tpoints[j].loc == loc)
			break;
	}
	/* If not, allocate a new one */
	if (j >= p->no_tpoints) {
		p->no_tpoints++;
		if ((p->tpoints = (trans_point *)realloc(p->tpoints, sizeof(trans_point) * p->no_tpoints)) == NULL)
			error("Realloc of tpoints");
		p->tpoints[j].loc = loc;
		for (i = 0; i < MAX_CHAN; i++)
			p->tpoints[j].val[i] = -1.0;
	}
	p->tpoints[j].val[chan] = val;

	if (p->no_tpoints > 0) {
		/* Sort the transfer points into loc order */
#define HEAP_COMPARE(A,B) ((A).loc < (B).loc)
		HEAPSORT(trans_point, p->tpoints, p->no_tpoints);
#undef HEAP_COMPARE
	}
}

/* Reurn nz if the target has been set */
static int pcaltarg_is_set(pcaltarg *p) {
	if (p->devmaxset != 0
	 || p->ademaxset != 0
	 || p->ademinset != 0
	 || p->no_tpoints > 0)
		return 1;
	return 0;
}

/* Update one from another */
static void pcaltarg_update(pcaltarg *p, pcaltarg *s) {
	int i, j, k;

	if (s->devmaxset) {
		if (p->devmaxset == 0) {
			for (i = 0; i < MAX_CHAN; i++)
				p->devmax[i] = -1.0;
			p->devmaxset = 1;
		}
		for (i = 0; i < MAX_CHAN; i++) {
			if (s->devmax[i] >= 0.0)
				p->devmax[i] = s->devmax[i];
		}
	}

	if (s->ademaxset) {
		if (p->ademaxset == 0) {
			for (i = 0; i < MAX_CHAN; i++)
				p->ademax[i] = -1.0;
			p->ademaxset = 1;
		}
		for (i = 0; i < MAX_CHAN; i++) {
			if (s->ademax[i] >= 0.0)
				p->ademax[i] = s->ademax[i];
		}
	}

	if (s->ademinset) {
		if (p->ademinset == 0) {
			for (i = 0; i < MAX_CHAN; i++)
				p->ademin[i] = -1.0;
			p->ademinset = 1;
		}
		for (i = 0; i < MAX_CHAN; i++) {
			if (s->ademin[i] >= 0.0)
				p->ademin[i] = s->ademin[i];
		}
	}

	/* For each source transfer curve point */
	for (k = 0; k < s->no_tpoints; k++) {

		/* See if a transfer curve point already exists */
		for (j = 0; j < p->no_tpoints; j++) {
			if (p->tpoints[j].loc == s->tpoints[k].loc)
				break;
		}
		/* If not, allocate a new one */
		if (j >= p->no_tpoints) {
			p->no_tpoints++;
			if ((p->tpoints = (trans_point *)realloc(p->tpoints, sizeof(trans_point) * p->no_tpoints)) == NULL)
				error("Realloc of tpoints");
			p->tpoints[j].loc = s->tpoints[k].loc;
			for (i = 0; i < MAX_CHAN; i++)
				p->tpoints[j].val[i] = -1.0;
		}
		for (i = 0; i < MAX_CHAN; i++) {
			if (s->tpoints[k].val[i] >= 0.0)
				p->tpoints[j].val[i] = s->tpoints[k].val[i]; 
		}
	}
	if (s->no_tpoints > 0) {
		/* Sort the transfer points into loc order */
#define HEAP_COMPARE(A,B) ((A).loc < (B).loc)
		HEAPSORT(trans_point, p->tpoints, p->no_tpoints);
#undef HEAP_COMPARE
	}
}

/* Create a new, empty pcaltarget */
/* Return NULL on error */
pcaltarg *new_pcaltarg() {
	pcaltarg *p;

	if ((p = (pcaltarg *)calloc(1, sizeof(pcaltarg))) == NULL) {
		return NULL;
	}

	/* Set method pointers */
	p->del = pcaltarg_del;
	p->write = pcaltarg_write;
	p->read = pcaltarg_read;
	p->update_devmax = pcaltarg_update_devmax;
	p->update_ademax = pcaltarg_update_ademax;
	p->update_ademin = pcaltarg_update_ademin;
	p->update_tcurve = pcaltarg_update_tcurve;
	p->is_set = pcaltarg_is_set;
	p->update = pcaltarg_update;

	return p;
}

/* - - - - - - - - - - - - - - - - - - - - - - - */

/* A wedge sample value */
typedef struct {
	double inv;			/* Input value (cal table) */
	double dev;			/* Device value */
	double XYZ[3];		/* XYZ value */
	double Lab[3];		/* Lab value */
	double del;			/* Absolute delta (to white) */
} wval;

#define MAX_INVSOLN	10	/* Rspl maximum reverse solutions */

/* rspl setting functions */
static void rsplset1(void *cbntx, double *out, double *in) {
	co *dpoints = (co *)cbntx;
	int ix;

	ix = *((int*)&in[-0-1]);	/* Get grid index being looked up */
	out[0] = dpoints[ix].v[0];
}

/* Do an inverse lookup of an rspl. Return -1.0 on error. */
/* dir is value to favour if there are multiple solutions. */
static double rspl_ilookup(rspl *r, double dir, double in) {
	int nsoln;				/* Number of solutions found */
	co pp[MAX_INVSOLN];		/* Room for all the solutions found */
	int k;					/* Chosen solution */

	pp[0].v[0] = in;

	nsoln = r->rev_interp (
		r,				 	/* this */
		RSPL_NEARCLIP,		/* Clip to nearest (faster than vector) */
		MAX_INVSOLN,		/* Maximum number of solutions allowed for */
		NULL, 				/* No auxiliary input targets */
		NULL,				/* Clip vector direction and length */
		pp);				/* Input and output values */

	nsoln &= RSPL_NOSOLNS;		/* Get number of solutions */

	if (nsoln == 1) { /* Exactly one solution */
		k = 0;
	} else if (nsoln == 0) {	/* Zero solutions. This is unexpected. */
		return -1.0;
	} else {		/* Multiple solutions */
		double bdist = 1e300;
		int bsoln = 0;
//		warning("Multiple solutions for curve %d for DE %f",j,pp[0].v[0]); 
		for (k = 0; k < nsoln; k++) {
			double tt;
			tt = pp[k].p[0] - dir;
			tt *= tt;
			if (tt < bdist) {	/* Better solution */
				bdist = tt;
				bsoln = k;
			}
		}
		k = bsoln;
	}
	return pp[k].p[0];
}

int main(int argc, char *argv[]) {
	int fa,nfa,mfa;				/* current argument we're looking at */
	int verb = 0;
	int doplot = 0;
	int initial = 0;			/* Do initial creation of cal target and calibration */
	int recal = 0;				/* Do recalibrate/use cal target. */
	int verify = 0;				/* Do verification */
	int imitate = 0;			/* Do target directly from input */
	int dowrite = 1;			/* Write to files */
	int doamp = 0;				/* Write Adobe Photoshop .AMP file */
	profxinf xpi;				/* Extra profile/calibration information */
	pcaltarg *upct = NULL;		/* User settings of print calibration target */
	pcaltarg *pct = NULL;		/* Settings of print calibration target */
	double smooth = RSPLSMOOTH;	/* RSPL Smoothness factor */
	double xsmooth = 1.0;		/* Smoothing multiplier */
	double ver_maxde = 2.0;		/* Verify maximum Delta E (1.0 for smooth == 1.0) */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType observ = icxOT_CIE_1931_2;	/* The classic observer */

	char baname[MAXNAMEL+1] = "";	/* Input & Output base name */
	char inname[MAXNAMEL+1] = "";	/* new .ti3 input file name */
	char calname[MAXNAMEL+1] = "";	/* previous .cal input file name */
	char outname[MAXNAMEL+1] = "";	/* new .cal output file name */
	char ampname[MAXNAMEL+1] = "";	/* new .amp output file name */
	double maxscale[MAX_CHAN];	/* Scale auto device maximum to % */
	cgats *icg = NULL;			/* .ti3 input cgats structure */
	int ti;						/* Temporary CGATs index */
	inkmask devmask;			/* ICX ink mask of device space */
	int devchan;				/* Number of chanels in device space */
	int isLab = 0;				/* Flag indicating whether PCS is XYZ or Lab */
	int n_pvals[MAX_CHAN];		/* Number of measurement values */
	wval *pvals[MAX_CHAN];		/* Patch measurement values */
	wval white;					/* Average white value */
	int n_white = 0;			/* Number of values to average */
	icmXYZNumber wht;			/* White value */
	rspl *raw[MAX_CHAN];		/* Raw Lab values fitted to rspl */
	rspl *ade[MAX_CHAN];		/* Absolute delta E */
	rspl *rde[MAX_CHAN];		/* Relative delta E */
	rspl *pcade[MAX_CHAN];		/* Previous calibrated absolute delta E */
	double mxade[MAX_CHAN];		/* Maximum ade value */
	double idpow[MAX_CHAN] = { -1.0 };		/* Ideal power-like of targen values */
	int n_cvals;				/* Number of calibration curve values */
	wval *cvals[MAX_CHAN];		/* Calibration curve tables */
	rspl *tcurves[MAX_CHAN];	/* Tweak target curves */

	int i, j;

	/* Init pointers to NULL */
	for (j = 0; j < MAX_CHAN; j++) {
		maxscale[j] = -1.0;
		pvals[j] = NULL;
		raw[j] = NULL;
		ade[j] = NULL;
		rde[j] = NULL;
		pcade[j] = NULL;
		cvals[j] = NULL;
		tcurves[j] = NULL;
	}

	error_program = argv[0];
	memset((void *)&xpi, 0, sizeof(profxinf));	/* Init extra profile info to defaults */
	if ((upct = new_pcaltarg()) == NULL || (pct = new_pcaltarg()) == NULL)
		error("new_caltarg failed");

	if (argc < 3)
		usage("Too few arguments, got %d expect at least %d",argc-1,2);

#ifdef NEVER
	{
		double src, dst, pp;

		src = 0.5;
		dst = 0.25;
		pp = icx_powlike_needed(src, dst); 
		printf("%f -> %f needs %f, check %f\n",src,dst,pp,icx_powlike(src,pp));
		
		src = 0.25;
		dst = 0.5;
		pp = icx_powlike_needed(src, dst); 
		printf("%f -> %f needs %f, check %f\n",src,dst,pp,icx_powlike(src,pp));

		src = 0.5;
		dst = 0.707106;
		pp = icx_powlike_needed(src, dst); 
		printf("%f -> %f needs %f, check %f\n",src,dst,pp,icx_powlike(src,pp));

		src = 0.5;
		dst = 0.5;
		pp = icx_powlike_needed(src, dst); 
		printf("%f -> %f needs %f, check %f\n",src,dst,pp,icx_powlike(src,pp));
	}
#endif // NEVER

	/* Process the arguments */
	mfa = 1;		/* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1+mfa) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage("Usage requested");

			else if (argv[fa][1] == 'v') {
				if (na != NULL) {
					fa = nfa;
					verb = atoi(na);
				} else
					verb = 1;
			}

			else if (argv[fa][1] == 'p') {
				if (na != NULL) {
					fa = nfa;
					doplot = atoi(na);
				} else
					doplot = 1;				/* Plot various graphs */
			}

			else if (argv[fa][1] == 'i') {
				initial = 1;			/* Initial calibration */
				recal = 0;
				verify = 0;
				imitate = 0;
				mfa = 1;
			}

			else if (argv[fa][1] == 'r') {
				initial = 0;
				recal = 1;			/* Recalibrate */
				verify = 0;
				imitate = 0;
				mfa = 2;
			}

			else if (argv[fa][1] == 'e') {
				initial = 0;
				recal = 0;
				verify = 1;			/* Verify */
				imitate = 0;
				mfa = 2;
			}

			else if (argv[fa][1] == 'I') {
				initial = 0;
				recal = 0;
				verify = 0;
				imitate = 1;			/* Imitation target */
				mfa = 1;
			}

			else if (argv[fa][1] == 'd')
				dowrite = 0;			/* Don't write to files */

			else if (argv[fa][1] == 'a')
				doamp = 1;			/* write AMP file */

			/* Smoothing modfider */
			else if (argv[fa][1] == 's') {
				if (na == NULL) usage("Expect argument to smoothing flag -s");
				fa = nfa;
				xsmooth = atof(na);
			}

			/* Manufacturer description string */
			else if (argv[fa][1] == 'A') {
				if (na == NULL) usage("Expect argument to manufacturer description flag -A");
				fa = nfa;
				xpi.deviceMfgDesc = na;
			}

			/* Model description string */
			else if (argv[fa][1] == 'M') {
				if (na == NULL) usage("Expect argument to model description flag -M");
				fa = nfa;
				xpi.modelDesc = na;
			}

			/* Profile Description */
			else if (argv[fa][1] == 'D') {
				if (na == NULL) usage("Expect argument to profile description flag -D");
				fa = nfa;
				xpi.profDesc = na;
			}

			/* Copyright string */
			else if (argv[fa][1] == 'C') {
				if (na == NULL) usage("Expect argument to copyright flag -C");
				fa = nfa;
				xpi.copyright = na;
			}

			/* Per channel target modifiers */
			else if (argv[fa][1] == 'x'
			      || argv[fa][1] == 'm'
			      || argv[fa][1] == 'n' 
			      || argv[fa][1] == 't') {
				char fch = argv[fa][1];
				int chan = -1;
				double val = -1.0;
				if (na == NULL)
					usage("Expect channel flag after flag -%c",argv[fa][1]);
				fa = nfa;
    			switch (na[0]) {
					case 'c': case 'r': case '0':
						chan = 0;
						break;
					case 'm': case 'g': case '1':
						chan = 1;
						break;
					case 'y': case 'b': case '2':
						chan = 2;
						break;
					case 'k': case '3':
						chan = 3;
						break;
					case '4':
						chan = 4;
						break;
					case '5':
						chan = 5;
						break;
					case '6':
						chan = 6;
						break;
					case '7':
						chan = 7;
						break;
					case '8':
						chan = 8;
						break;
					case '9':
						chan = 9;
						break;
					case 'A':
						chan = 10;
						break;
					case 'B':
						chan = 11;
						break;
					case 'C':
						chan = 12;
						break;
					case 'D':
						chan = 13;
						break;
					case 'E':
						chan = 14;
						break;
					case 'F':
						chan = 15;
						break;
					default:
						usage("Unknown channel flag '%s' after flag -%c",argv[fa][2],argv[fa][1]);
				}
				++fa;
				if (fa >= argc || argv[fa][0] == '-') usage("Expect argument after flag -%c%c",fch,na[0]);
				val = atof(argv[fa]);
			
				if (fch == 'x') {
					if (val < 0.0 || val > 100.0)
						usage("Argument to -%c%c %f from '%s' is out of range",fch,na[0],val,argv[fa]);
					val /= 100.0;
					upct->update_devmax(upct, chan, val);

				} else if (fch == 'm') {
					if (val < 0.0 || val > 100.0)
						usage("Argument to -%c%c %f from '%s' is out of range",fch,na[0],val,argv[fa]);
					val /= 100.0;
					maxscale[chan] = val;

				} else if (fch == 'n') {
					upct->update_ademin(upct, chan, val);

				} else if (fch == 't') {
					if (val < 0.0 || val > 100.0)
						usage("Argument to -%c%c %f from '%s' is out of range",fch,na[0],val,argv[fa]);
					val /= 100.0;
					upct->update_tcurve(upct, chan, 0.5, val);
				}
			}
			else 
				usage("Unknown flag '%c'",argv[fa][1]);
		} else
			break;
	}

	smooth *= xsmooth;

	if (!(   (initial && !recal && !verify && !imitate)
	      || (!initial && recal && !verify && !imitate)
	      || (!initial && !recal && verify && !imitate)
	      || (!initial && !recal && !verify && imitate)))
		error("One of -i, -r -e or -I must be set");

	/* Get the file name arguments */
	if (verify || recal) {
		if (fa >= argc || argv[fa][0] == '-') usage("Missing prevoius .cal basename");
		strncpy(calname,argv[fa++],MAXNAMEL-4); calname[MAXNAMEL-4] = '\000';
		strcat(calname,".cal");
	}
	if (fa >= argc || argv[fa][0] == '-') usage("Missing .ti3 and new .cal basename");
	strncpy(baname,argv[fa++],MAXNAMEL-4); baname[MAXNAMEL-4] = '\000';
	strcpy(inname,baname);		/* new .ti3 file */
	strcat(inname,".ti3");
	strcpy(outname,baname);		/* New .cal file */
	strcat(outname,".cal");
	strcpy(ampname,baname);		/* New .amp file */
	strcat(ampname,".amp");

	if (fa < argc) usage("Too many arguments ('%s')",argv[fa]);

	/* Open and look at the .ti3 profile patches file */
	icg = new_cgats();				/* Create a CGATS structure */
	icg->add_other(icg, "CTI3"); 	/* our special input type is Calibration Target Information 3 */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI3 format file");
	if (icg->ntables < 1)
		error ("Input file doesn't contain at least one table");

	/* See if CIE is actually available - some sources of .TI3 don't provide it */
	if (!spec
	 && icg->find_field(icg, 0, "LAB_L") < 0
	 && icg->find_field(icg, 0, "XYZ_X") < 0) {

		if (icg->find_kword(icg, 0, "SPECTRAL_BANDS") < 0)
			error ("Neither CIE nor spectral data found in file '%s'",inname);

		/* Switch to using spectral information */
		if (verb)
			printf("No CIE data found, switching to spectral with standard observer & D50\n");
		spec = 1;
		illum = icxIT_D50;
		observ = icxOT_CIE_1931_2;
	}
	
	/* If we requested spectral, check that it is available */
	if (spec) {
		if (icg->find_kword(icg, 0, "SPECTRAL_BANDS") < 0)
			error ("Requested spectral interpretation when data not available");
	}

	/* Get colorspace information from input CGATS file */
	{
		char *buf;
		char *inc, *outc;

		if ((ti = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
			error("Input file doesn't contain keyword COLOR_REPS");

		if ((buf = strdup(icg->t[0].kdata[ti])) == NULL)
			error("Malloc failed - color rep");

		/* Split COLOR_REP into device and PCS space */
		inc = buf;
		if ((outc = strchr(buf, '_')) == NULL)
			error("COLOR_REP '%s' invalid", icg->t[0].kdata[ti]);
		*outc++ = '\000';

		if (strcmp(outc, "XYZ") == 0)
			isLab = 0;
		else if (strcmp(outc, "LAB") == 0)
			isLab = 1;
		else
			error("COLOR_REP '%s' invalid (Neither XYZ nor LAB)", icg->t[0].kdata[ti]);

		devmask = icx_char2inkmask(inc); 
		devchan = icx_noofinks(devmask);
		
		if (devchan == 0)
			error("COLOR_REP '%s' invalid (No matching devmask)", icg->t[0].kdata[ti]);
		
		if ((devmask & ICX_ADDITIVE) && !(devmask & ICX_INVERTED))
			warning("COLOR_REP '%s' is probably not suitable for print calibration!", icg->t[0].kdata[ti]);

		free(buf);
	}

	if (verify || recal || imitate) {
		if (upct->is_set(upct)) {
			warning("Command line calibration target paramers ignored on re-calibrate, verify and imitate!");
		}
	}

	/* For recalibrate or verify, load the previous calibration file */
	if (verify || recal) {
		cgats *tcg;						/* Previous .cal file */

		tcg = new_cgats();				/* Create a CGATS structure */
		tcg->add_other(tcg, "CAL"); 	/* our special input type is Calibration Target */

		if (tcg->read_name(tcg, calname))
			error("No cal target '%s' found for re-calibrate (%s)\n",calname,tcg->err);

		/* Check that this is an output cal file */
		if ((ti = tcg->find_kword(tcg, 0, "DEVICE_CLASS")) < 0)
			error ("Calibration file '%s'doesn't contain keyword DEVICE_CLASS",calname);
		if (strcmp(tcg->t[0].kdata[ti],"OUTPUT") != 0)
			error ("Calibration file '%s' doesn't has DEVICE_CLASS that is not OUTPUT",calname);

		if (pct->read(pct, tcg, 1) != 0)
			error("Reading cal target '%s' failed",calname);

		if (pct->devmask != devmask)
			error("Target '%s' colorspace '%s' doesn't match '%s' colorspace '%s'",
				   calname,icx_inkmask2char(pct->devmask, 1),inname,icx_inkmask2char(devmask, 1)); 

		/* Load the previous expected absolute DE response */
		/* It will be in the third table with other type "CAL" */
		if (tcg->ntables >= 3 && tcg->t[2].tt == tt_other && tcg->t[0].oi == 0) {
			int ti;
			char *bident;
			int spi[1+MAX_CHAN];	/* CGATS indexes for each field */
			char buf[100];

			bident = icx_inkmask2char(pct->devmask, 0); 

			if (tcg->t[2].nsets <= 0)
				error ("No Calibration Expected DE Response in '%s'",calname);

			/* Figure out the indexes of all the fields */
			sprintf(buf, "%s_I_DE",bident);
			if ((spi[0] = tcg->find_field(tcg, 2, buf)) < 0)
				error("Can't find field %s in '%s'",buf,calname);

			for (i = 0; i < devchan; i++) {
				inkmask imask = icx_index2ink(pct->devmask, i);
				sprintf(buf, "%s_%s_DE",bident,icx_ink2char(imask));
				if ((spi[1+i] = tcg->find_field(tcg, 2, buf)) < 0)
					error("Can't find field %s in '%s'",buf,calname);
			}

			/* Read in each channels values and put them in a rspl */
			for (j = 0; j < devchan; j++) {
				datai low,high;
				int gres[MXDI];
				co *dpoints;

				low[0] = 0.0;
				high[0] = 1.0;
				gres[0] = tcg->t[2].nsets;

				if ((pcade[j] = new_rspl(RSPL_NOFLAGS,1, 1)) == NULL)
					error("new_rspl() failed");

				if ((dpoints = malloc(sizeof(co) * gres[0])) == NULL)
					error("malloc dpoints[%d] failed",gres[0]);

				/* Copy the points to our array */
				if (devmask & ICX_ADDITIVE) {
					for (i = 0; i < gres[0]; i++) {
						dpoints[i].p[0] = 1.0 - i/(double)(gres[0]-1);
						dpoints[i].v[0] = *((double *)tcg->t[2].fdata[gres[0]-1-i][spi[1+j]]);
					}
				} else {
					for (i = 0; i < gres[0]; i++) {
						dpoints[i].p[0] = i/(double)(gres[0]-1);
						dpoints[i].v[0] = *((double *)tcg->t[2].fdata[i][spi[1+j]]);
					}
				}

				pcade[j]->set_rspl(pcade[j],
						   0, 
						   (void *)dpoints,		/* Read points */
						   rsplset1,			/* Setting function */
						   low, high, gres,		/* Low, high, resolution of grid */
						   NULL, NULL			/* Default data scale */
						   );
				free(dpoints);
			}
			free(bident);
		}
		tcg->del(tcg);

	} else {	/* Must be an initial or Imitation calibration */

		pct->devmask = devmask;

		/* Set the cal target from any user supplied parameters */
		pct->update(pct, upct);

		/* No previous absolute de reference */
		for (j = 0; j < devchan; j++)
			pcade[j] = NULL;
	}

	/* Common processing: */

	/* Read in the patch data */
	{
		char buf[100];
		char *pcsfname[2][3] = { { "XYZ_X", "XYZ_Y", "XYZ_Z" },
		                         { "LAB_L", "LAB_A", "LAB_B" } };
		int dvi[MAX_CHAN];	/* CGATS indexes for each device field */
		int pcsix[3];		/* XYZ/Lab chanel indexes */
		xsp2cie *sp2cie = NULL;		/* Spectral conversion object */
		xspect sp;
		int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
		char *bident = icx_inkmask2char(devmask, 0); 
	
		/* Figure out the indexes of all the device fields */
		for (j = 0; j < devchan; j++) {
			inkmask imask = icx_index2ink(devmask, j);
			sprintf(buf, "%s_%s",bident,icx_ink2char(imask));
			if ((dvi[j] = icg->find_field(icg, 0, buf)) < 0)
				error("Can't find field %s in '%s'",buf,inname);
#ifdef DEBUG
			printf("devn chan %d field %s = %d\n",j,buf,dvi[j]);
#endif
		}
		free(bident);

		if (spec) {
			int ii;
			char buf[100];

			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_BANDS")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_BANDS");
			sp.spec_n = atoi(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_START_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_START_NM");
			sp.spec_wl_short = atof(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_END_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_END_NM");
			sp.spec_wl_long = atof(icg->t[0].kdata[ii]);
			sp.norm = 100.0;

			/* Find the fields for spectral values */
			for (j = 0; j < sp.spec_n; j++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
				            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);

				if ((spi[j] = icg->find_field(icg, 0, buf)) < 0)
					error("Input file doesn't contain field %s",buf);
			}

			/* Create a spectral conversion object to XYZ */
			if ((sp2cie = new_xsp2cie(illum, 0.0, &cust_illum, observ, NULL, icSigXYZData, icxClamp)) == NULL)
				error("Creation of spectral conversion object failed");

			/* To add FWA comp. would have to locate/create spectral white here, */
			/* then set the FWA comp. on. */
			/* See profout.c */

		} else {
			/* Figure out the indexes of the PCS fields */
			for (j = 0; j < 3; j++) {
				if ((i = icg->find_field(icg, 0, pcsfname[isLab][j])) >= 0) {
					if (icg->t[0].ftype[i] != r_t)
						error ("Field %s is wrong type",pcsfname[isLab][j]);
					pcsix[j] = i;
#ifdef DEBUG
					printf("PCS chan %d field %s = %d\n",j,pcsfname[isLab][j],pcsix[j]);
#endif
				} else {
					error ("Failed to find field %s",pcsfname[isLab][j]);
				}
			}
		}
	
		n_cvals = 0;
		for (j = 0; j < devchan; j++) {
			pvals[j] = NULL;
			n_pvals[j] = 0;
			cvals[j] = NULL;
		}

		/* Read all the test patches in */
		for (i = 0; i < icg->t[0].nsets; i++) {
			double maxv = -1.0;
			int maxch;

#ifdef DEBUG
			printf("Reading patch %d\n",i);
#endif

			/* Locate the maximum device value of any channel */
			for (j = 0; j < devchan; j++) {
				double val = *((double *)icg->t[0].fdata[i][dvi[j]]) / 100.0;
				if (devmask & ICX_ADDITIVE)
					val = 1.0 - val;
				if (val > maxv) {
					maxv = val;
					maxch = j;
				}
			}
#ifdef DEBUG
			printf("max %f at chan %d\n",maxv,maxch);
#endif
			/* Treat white specially, and take it out of the list */
			if (maxv < 1e-6) {
				double wxyz[3];
				if (n_white == 0) {
					white.dev = 0.0;
					for (j = 0; j < 3; j++)
						white.XYZ[j] = 0.0;
				}
				if (spec) {
					/* Read and convert the spectral value */
					for (j = 0; j < sp.spec_n; j++)
						sp.spec[j] = *((double *)icg->t[0].fdata[i][spi[j]]);

					sp2cie->convert(sp2cie, wxyz, &sp);

				} else {
					/* Read the CIE value */
					for (j = 0; j < 3; j++)
						wxyz[j] = *((double *)icg->t[0].fdata[i][pcsix[j]]);

					/* And convert to XYZ 0..1 */
					if (isLab) {
						icmLab2XYZ(&icmD50, wxyz, wxyz);
					} else {
						for (j = 0; j < 3; j++)
							wxyz[j] /= 100.0;
					}
				}
			
				white.dev += maxv;
				for (j = 0; j < 3; j++)
					white.XYZ[j] += wxyz[j];
				n_white++;
#ifdef DEBUG
				printf("  white: dev %f,XYZ %f %f %f\n",
				white.dev,white.XYZ[0]/n_white, white.XYZ[1]/n_white, white.XYZ[2]/n_white);
#endif
			} else {
				wval *vp;
				/* Check that all the non-max value channels are zero */
				for (j = 0; j < devchan; j++) {
					double val = *((double *)icg->t[0].fdata[i][dvi[j]]) / 100.0;
					if (devmask & ICX_ADDITIVE)
						val = 1.0 - val;
					if (j == maxch)
						continue;
					if (val > 0.001)
						break;
				}
				if (j < devchan) {
#ifdef DEBUG
				printf("Skipping patch\n");
#endif
					continue;		/* Ignore this patch */
				}

				if ((pvals[maxch] = (wval *)realloc(pvals[maxch],
				                       sizeof(wval) * (n_pvals[maxch]+1))) == NULL)
					error("Realloc of pvals failed");
				vp = &pvals[maxch][n_pvals[maxch]];
				vp->dev = maxv;

				if (spec) {
					/* Read and convert the spectral value */
					for (j = 0; j < sp.spec_n; j++)
						sp.spec[j] = *((double *)icg->t[0].fdata[i][spi[j]]);

					sp2cie->convert(sp2cie, vp->XYZ, &sp);

				} else {
					/* Read the CIE value */
					for (j = 0; j < 3; j++)
						vp->XYZ[j] = *((double *)icg->t[0].fdata[i][pcsix[j]]);

					/* And convert to XYZ 0..1 */
					if (isLab) {
						icmLab2XYZ(&icmD50, vp->XYZ, vp->XYZ);
					} else {
						for (j = 0; j < 3; j++)
							vp->XYZ[j] /= 100.0;
					}
				}
				/* Temporary D50 Lab */
				icmXYZ2Lab(&icmD50, vp->Lab, vp->XYZ);
#ifdef DEBUG
				printf("  patch %d: dev %f,XYZ %f %f %f, D50 Lab %f %f %f\n",
				n_pvals[maxch], vp->dev,vp->XYZ[0], vp->XYZ[1], vp->XYZ[2],
				vp->Lab[0], vp->Lab[1], vp->Lab[2]);
#endif
				n_pvals[maxch]++;
			}
		}

		/* Average the white */
		if (n_white == 0)
			error("Can't find even one white patch in '%s'",inname);
#ifdef DEBUG
			printf("% white patches\n",n_white);
#endif
		for (j = 0; j < 3; j++) {
			white.dev    /= (double)n_white;
			white.XYZ[j] /= (double)n_white;
		}
		icmAry2XYZ(wht, white.XYZ);
		icmXYZ2Lab(&icmD50, white.Lab, white.XYZ);

		/* Convert the Lab white reference to absolute */
		wht.X /= wht.Y;
		wht.Z /= wht.Y;
		wht.Y /= wht.Y;

		if (verb) {
			icmXYZ2Lab(&icmD50, white.Lab, white.XYZ);
			printf("Average white = XYZ %f %f %f, D50 Lab %f %f %f\n",
			white.XYZ[0], white.XYZ[1], white.XYZ[2], white.Lab[0], white.Lab[1], white.Lab[2]);
		}

		for (j = 0; j < devchan; j++) {
			wval *wp;

			/* Add averaged white back into each channel */
			if ((pvals[j] = (wval *)realloc(pvals[j],
			                       sizeof(wval) * (n_pvals[j]+1))) == NULL)
				error("Realloc (%d) of pvals failed",n_pvals[j]+1);
			wp = &pvals[j][n_pvals[j]];
			wp->dev = white.dev;
			wp->XYZ[0] = white.XYZ[0];
			wp->XYZ[1] = white.XYZ[1];
			wp->XYZ[2] = white.XYZ[2];
			n_pvals[j]++;

			/* Convert all the XYZ values to Lab paper relative */
			for (i = 0; i < n_pvals[j]; i++) {
				wp = &pvals[j][i];
				icmXYZ2Lab(&wht, wp->Lab, wp->XYZ);
			}

			/* Sort the channel acording to device value */ 
			/* For a consistent result for identical device values, */
			/* secondary sort by inverse CIE value */
//#define HEAP_COMPARE(A,B) ((A).dev < (B).dev)
#define HEAP_COMPARE(A,B) ((A).dev != (B).dev ? ((A).dev < (B).dev) : ((A).Lab[0] > (B).Lab[0])) 
			HEAPSORT(wval, pvals[j], n_pvals[j]);
#undef HEAP_COMPARE

			/* Check the maximum value looks OK */
			if (n_pvals[j] < 5)
				warning("Channel %d has only %d test patches",n_pvals[j]);
			if (pvals[j][n_pvals[j]-1].dev < 0.99)
				warning("Channel %d has max test patch value of %f",pvals[j][n_pvals[j]-1]); 

			if (verb > 1) {
				printf("Chan %d has %d raw values:\n",j,n_pvals[j]);
				for (i = 0; i < n_pvals[j]; i++) {
					wp = &pvals[j][i];
					printf("  %d: dev %f,XYZ %f %f %f, Lab %f %f %f\n",
					       i,wp->dev,wp->XYZ[0], wp->XYZ[1], wp->XYZ[2],
					         wp->Lab[0], wp->Lab[1], wp->Lab[2]);
				}
			}
		}

		if (sp2cie != NULL)
			sp2cie->del(sp2cie);
	}
	icg->del(icg);		/* Clean up */

	/* Interpolate Lab using rspl */
	for (j = 0; j < devchan; j++) {
		datai low,high;
		datao olow,ohigh;
		int gres[MXDI];
		double avgdev[MXDO];
		cow *dpoints;

		low[0] = 0.0;
		high[0] = 1.0;
		gres[0] = GRES;
		olow[0] = 0.0;
		ohigh[0] = 100.0;
		olow[1] = olow[2] = -128.0; 
		ohigh[1] = ohigh[2] = 128.0;
		avgdev[0] = 0.0025;
		avgdev[1] = 0.005;
		avgdev[2] = 0.005;

		if ((raw[j] = new_rspl(RSPL_NOFLAGS,1, 3)) == NULL)
			error("new_rspl() failed");

		if ((dpoints = (cow *)malloc(sizeof(cow) * n_pvals[j])) == NULL)
			error("malloc dpoints[%d] failed",n_pvals[j]);

		for (i = 0; i < n_pvals[j]; i++) {
			dpoints[i].p[0] = pvals[j][i].dev;
			dpoints[i].v[0] = pvals[j][i].Lab[0];
			dpoints[i].v[1] = pvals[j][i].Lab[1];
			dpoints[i].v[2] = pvals[j][i].Lab[2];
			if (i == 0)
				dpoints[i].w = (double)n_white;
			else
				dpoints[i].w = 1.0;
		}

		raw[j]->fit_rspl_w(raw[j],
		           RSPLFLAGS,
		           dpoints,			/* Test points */
		           n_pvals[j],			/* Number of test points */
		           low, high, gres,		/* Low, high, resolution of grid */
		           olow, ohigh,			/* Default data scale */
		           smooth,				/* Smoothing */
		           avgdev,				/* Average deviation */
		           NULL);				/* iwidth */


		/* Compute & show fit quality */
		if (verb > 0) {
			double avgde = 0.0, maxde = 0.0;
			for (i = 0; i < n_pvals[j]; i++) {
				co tp;	/* Test point */
				double de;
				tp.p[0] = pvals[j][i].dev;
				raw[j]->interp(raw[j], &tp);
				de = icmLabDE(pvals[j][i].Lab, tp.v);

				avgde += de;
				if (de > maxde)
					maxde = de;
			}
			avgde /= (double)n_pvals[j];
			printf("Chan %d raw fit avg DE %f, max %f\n",j,avgde,maxde);
		}

		free(dpoints);
	}

	/* Plot the raw curves */
	if (doplot > 1) {
		double xx[PRES];
		double yy[3][PRES];

		for (j = 0; j < devchan; j++) {
			printf("Chan %d raw L*a*b*:\n",j);
			for (i = 0; i < PRES; i++) {
				co tp;	/* Test point */
				xx[i] = i/(double)(PRES-1);
				tp.p[0] = xx[i];
				raw[j]->interp(raw[j], &tp);
				yy[0][i] = tp.v[0];
				yy[1][i] = tp.v[1];
				yy[2][i] = tp.v[2];
			}
			do_plot(xx, yy[0], yy[1], yy[2], PRES);
		}
	}

	/* Create a RSPL of absolute deltaE and relative deltaE '94 */ 
	for (j = 0; j < devchan; j++) {
		datai low,high;
		int gres[MXDI];
		double avgdev[MXDO];
		co *dpoints_a;
		co *dpoints_r;
		double wh[3], prev[3], tot;

		low[0] = 0.0;
		high[0] = 1.0;
		gres[0] = GRES;
		avgdev[0] = 0.0;

		if ((ade[j] = new_rspl(RSPL_NOFLAGS,1, 1)) == NULL)
			error("new_rspl() failed");
		if (imitate) {
			if ((pcade[j] = new_rspl(RSPL_NOFLAGS,1, 1)) == NULL)
				error("new_rspl() failed");
		}
		if ((rde[j] = new_rspl(RSPL_NOFLAGS,1, 1)) == NULL)
			error("new_rspl() failed");

		if ((dpoints_a = malloc(sizeof(co) * GRES)) == NULL)
			error("malloc dpoints[%d] failed",GRES);
		if ((dpoints_r = malloc(sizeof(co) * GRES)) == NULL)
			error("malloc dpoints[%d] failed",GRES);

//printf("~1 Chan %d:\n",j);
		for (i = 0; i < GRES; i++) {
			co tp;	/* Test point */

			tp.p[0] = i/(double)(GRES-1);
			raw[j]->interp(raw[j], &tp);

			dpoints_a[i].p[0] = tp.p[0];
			dpoints_r[i].p[0] = tp.p[0];
			if (i == 0) {
//printf("~1 wht = %f %f %f\n",tp.v[0],tp.v[1],tp.v[2]);
				tot = 0.0;
				prev[0] = wh[0] = tp.v[0];
				prev[1] = wh[1] = tp.v[1];
				prev[2] = wh[2] = tp.v[2];
				dpoints_a[i].v[0] = 0.0;
				dpoints_r[i].v[0] = 0.0;
			} else {
//printf("~1 samp %d = %f %f %f\n",i,tp.v[0],tp.v[1],tp.v[2]);
				/* Use Euclidean for large DE: (CIE94 stuffs up here) */
				dpoints_a[i].v[0] = icmLabDE(tp.v, wh);
				/* And CIE94 for small: */
				tot += icmCIE94(tp.v, prev);
				prev[0] = tp.v[0];
				prev[1] = tp.v[1];
				prev[2] = tp.v[2];
				dpoints_r[i].v[0] = tot;
			}
//printf("~1   %d: dev %f, ade %f, rde %f\n",i,tp.p[0],dpoints_a[i].v[0],dpoints_r[i].v[0]);
		}

		ade[j]->set_rspl(ade[j],
		           0, 
		           (void *)dpoints_a,	/* Test points */
		           rsplset1,			/* Setting function */
		           low, high, gres,		/* Low, high, resolution of grid */
		           NULL, NULL			/* Default data scale */
		           );
		if (imitate) {
			pcade[j]->set_rspl(pcade[j],
			           0, 
			           (void *)dpoints_a,	/* Test points */
			           rsplset1,			/* Setting function */
			           low, high, gres,		/* Low, high, resolution of grid */
			           NULL, NULL			/* Default data scale */
			           );
		}
		rde[j]->set_rspl(rde[j],
		           0, 
		           (void *)dpoints_r,		/* Test points */
		           rsplset1,			/* Setting function */
		           low, high, gres,		/* Low, high, resolution of grid */
		           NULL, NULL			/* Default data scale */
		           );
		free(dpoints_a);
		free(dpoints_r);
	}

	if (initial) {
		/* Establish the ademax values */
		pct->update_devmax(pct, -1, -1.0);		/* Make sure there is a value for each */
		pct->update_ademax(pct, -1, -1.0);
		pct->update_ademin(pct, -1, -1.0);
		for (j = 0; j < devchan; j++) {
			co tp;	/* Test point */

			if (pct->devmax[j] < 0.0) {		/* Auto */
				double maxd, maxde, maxix;
				/* Locate the point of maximum aDE */
				for (maxde = -1.0, i = 0; i < GRES; i++) {

					tp.p[0] = i/(GRES-1.0);
					ade[j]->interp(ade[j], &tp);
					if (tp.v[0] > maxde) {
						maxd = tp.p[0];
						maxde = tp.v[0];
						maxix = i;
					}
				}
				pct->devmax[j] = maxd;
				pct->ademax[j] = maxde;		/* Temporary */
//printf("Chan %d, dev %f, max de = %f\n", j, maxd, maxde);

				if (maxd < 0.2) {
					warning("Chan %d, max DE point %f is below < 0.2 - ignored\n", j, maxd);
					maxix = GRES-1;
				}
				/* Then locate the point below that where the slope */
				/* becomes reasonable. */ 
				for (i = maxix; i >= 40; i--) {
					double aslope, minslope = 1e6;
					double naslope, nminslope;
					int k;

					/* Compute the minimum over a span of 20/GRES */
					for (k = 0; k < 40; k++) {
						double dp, dv, slope;

						tp.p[0] = (i-k)/(GRES-1.0);
						dp = tp.p[0];
						ade[j]->interp(ade[j], &tp);
						dv = tp.v[0];

						tp.p[0] = (i-k-1)/(GRES-1.0);
						ade[j]->interp(ade[j], &tp);
						slope = (dv - tp.v[0])/(dp - (i-k-1)/(GRES-1.0));
						if (k == 0)
							aslope = slope;
//printf("  Chan %d, dev %f, dv = %f, slope = %f\n", j, (i-k)/(GRES-1.0),dv - tp.v[0],slope);
						if (slope < minslope)
							minslope = slope;
					}
//printf("Chan %d, dev %f, aslope = %f, min slope = %f\n", j, i/(GRES-1.0),aslope,minslope);
					/* Normalize the slopes */
					naslope = aslope * SLOPE_NORM/pct->ademax[j];
					nminslope = minslope * SLOPE_NORM/pct->ademax[j];
//printf("Chan %d, dev %f, norm aslope = %f, min slope = %f\n", j, i/(GRES-1.0),naslope,nminslope);

					if (naslope > MIN_SLOPE_A && nminslope >= MIN_SLOPE_O)
						break;
				}
				pct->devmax[j] = i/(GRES-1.0);

				/* Scale auto max device value */
				if (maxscale[j] >= 0.0)
					pct->devmax[j] *= maxscale[j];

			/* Manually set initial dev max */
			} else {
				if (maxscale[j] >= 0.0)
					warning("Chan %d, scale %.1f%% of auto max ignored since max override used\n", j, maxscale[j] * 100.0);
			}

			/* Lookup devmax to set ademax */
			tp.p[0] = pct->devmax[j];
			ade[j]->interp(ade[j], &tp);
			pct->ademax[j] = tp.v[0];

	
			/* Establish a default ademin value */
			if (pct->ademin[j] < 0.0)
				pct->ademin[j] = 0.0;
		}

	} else if (recal) {

		/* Since the plot markers use devmax, look it up */
		for (j = 0; j < devchan; j++) {
			if ((pct->devmax[j] = rspl_ilookup(ade[j], 0.5, pct->ademax[j])) < 0.0)
				error("Unexpected failure to invert curve %d for ADE %f",j,pct->ademax[j]); 
		}
	}

	/* Find the maximum aDE value for each curve */
	for (j = 0; j < devchan; j++) {
		co tp;	/* Test point */
		mxade[j] = -1e6;
		for (i = 0; i < PRES; i++) {
			tp.p[0] = i/(double)(PRES-1);
			ade[j]->interp(ade[j], &tp);
			if (tp.v[0] > mxade[j])
				mxade[j] = tp.v[0];
		}
	}

	if (initial || recal) {
		/* Compute an ideal power-like value for test target */
		for (j = 0; j < devchan; j++) {
			double hdv;			/* Half device value */
			double hdvrde;		/* Half device rDE value */
			double thdvrde;		/* Target half device rDE value */
			double thdv;		/* Target half device value */
			double fdvrde;		/* Full device rDE value */
			co tp;	/* Test point */
	
			/* full rDE */
			tp.p[0] = pct->devmax[j];
			rde[j]->interp(rde[j], &tp);
			fdvrde = tp.v[0];

			/* Half device value of maximum */
			hdv = pct->devmax[j] * 0.5;
	
			/* rDE value half the device value */
			tp.p[0] = hdv;
			rde[j]->interp(rde[j], &tp);
			hdvrde = tp.v[0];
	
			/* rDE value we'd like at half the device value */
			thdvrde = 0.5 * fdvrde;
			
			/* Device value to get the rDE value we'd like at half */
			if ((thdv = rspl_ilookup(rde[j], 0.5, thdvrde)) < 0.0)
				error("Unexpected failure to invert curve %d for ADE %f",j,thdvrde); 

//printf("hdv %f, hdvrde %f, thdvrde %f, fdvrde %f, thdv %f\n",hdv,hdvrde,thdvrde, fdvrde,thdv);
	
			/* Power like value needed to get rDE value we'd like at hald device */
			idpow[j] = icx_powlike_needed(hdv, thdv); 
			
		}
	}

	if (verb > 1) {
		printf("Abs DE values:\n");
		for (i = 0; i < PRES; i++) {
			co tp;	/* Test point */
			tp.p[0]= i/(double)(PRES-1);

			printf("  dev %f, aDE",tp.p[0]);
			for (j = 0; j < 6 && j < devchan; j++) {
				ade[j]->interp(ade[j], &tp);
				printf(" %f",tp.v[0]);
			}
			printf("\n");
		}
		printf("Rel DE values:\n");
		for (i = 0; i < PRES; i++) {
			co tp;	/* Test point */
			tp.p[0]= i/(double)(PRES-1);

			printf("  dev %f, rdev",tp.p[0]);
			for (j = 0; j < 6 && j < devchan; j++) {
				rde[j]->interp(rde[j], &tp);
				printf(" %f",tp.v[0]);
			}
			printf("\n");
		}
	}

	if (initial || recal) {
		if (verb && pct->is_set(pct))  {
			for (j = 0; j < devchan; j++) {
				printf("Chan %d Dev max %f, aDE Max %f, aDE Min %f\n",j,pct->devmax[j],pct->ademax[j],pct->ademin[j]);
			}
		}

		if (verb)  {
			double avgpow = 0.0;
			for (j = 0; j < devchan; j++) {
				printf("Chan %d ideal targen power = %f\n",j,idpow[j]);
				avgpow += idpow[j];
			}
			avgpow /= (double)devchan;
			printf("Average ideal targen power = %f\n",avgpow);
		}
	}

	/* Plot both the delta E curves, and markers */
	if (doplot) {
		co tp;	/* Test point */
		double xx[PRES];
		double yy[10][PRES];
		double cx[10], cy[10];
		int nmark;

		printf("Absolute DE plot:\n");

		for (i = 0; i < PRES; i++) {
			xx[i] = i/(double)(PRES-1);

			for (j = 0; j < 10 && j < devchan; j++) {
				tp.p[0] = xx[i];
				ade[j]->interp(ade[j], &tp);
				yy[j][i] = tp.v[0];
			}
		}
		nmark = 0;
		if (pct->is_set(pct)) {
			/* Add markers for deMax */
			for (j = 0; j < 10 && j < devchan; j++) {
				cx[j] = pct->devmax[j];
				cy[j] = pct->ademax[j];
				nmark++;
			}
		}
		do_plot10p(xx, devchan > 3 ? yy[3] : NULL,
		               devchan > 1 ? yy[1] : NULL,
		               devchan > 4 ? yy[4] : NULL,
		               devchan > 0 ? yy[0] : NULL,
		               devchan > 2 ? yy[2] : NULL,
		               devchan > 5 ? yy[5] : NULL,
		               devchan > 6 ? yy[6] : NULL,
		               devchan > 7 ? yy[7] : NULL,
		               devchan > 8 ? yy[8] : NULL,
		               devchan > 9 ? yy[9] : NULL,
		               PRES,
					   cx, cy, verify ? 0 : nmark);

		printf("Relative DE plot:\n");

		for (i = 0; i < PRES; i++) {
			xx[i] = i/(double)(PRES-1.0);

			for (j = 0; j < 10 && j < devchan; j++) {
				tp.p[0] = xx[i];
				rde[j]->interp(rde[j], &tp);
				yy[j][i] = tp.v[0];
			}
		}
		nmark = 0;
		if (pct->is_set(pct)) {
			/* Add markers for deMax */
			for (j = 0; j < 10 && j < devchan; j++) {
				cx[j] = pct->devmax[j];
				tp.p[0] = cx[j];
				rde[j]->interp(rde[j], &tp);
				cy[j] = tp.v[0];
				nmark++;
			}
		}
		do_plot10p(xx, devchan > 3 ? yy[3] : NULL,
		             devchan > 1 ? yy[1] : NULL,
		             devchan > 4 ? yy[4] : NULL,
		             devchan > 0 ? yy[0] : NULL,
		             devchan > 2 ? yy[2] : NULL,
		             devchan > 5 ? yy[5] : NULL,
		             devchan > 6 ? yy[6] : NULL,
		             devchan > 7 ? yy[7] : NULL,
		             devchan > 8 ? yy[8] : NULL,
		             devchan > 9 ? yy[9] : NULL,
		             PRES,
					 cx, cy, verify ? 0 : nmark);

		if (idpow[0] > 0.0) {
			printf("Relative DE plot with ideal targen power applied:\n");

			for (i = 0; i < PRES; i++) {
				xx[i] = i/(double)(PRES-1.0);
	
				for (j = 0; j < 10 && j < devchan; j++) {
					tp.p[0] = icx_powlike(xx[i],idpow[j]);
					rde[j]->interp(rde[j], &tp);
					yy[j][i] = tp.v[0];
				}
			}
			nmark = 0;
			if (pct->is_set(pct)) {
				/* Add markers for deMax */
				for (j = 0; j < 10 && j < devchan; j++) {
					cx[j] = pct->devmax[j];
					tp.p[0] = icx_powlike(cx[j],idpow[j]);
					rde[j]->interp(rde[j], &tp);
					cy[j] = tp.v[0];
					nmark++;
				}
			}
			do_plot10p(xx, devchan > 3 ? yy[3] : NULL,
			             devchan > 1 ? yy[1] : NULL,
			             devchan > 4 ? yy[4] : NULL,
			             devchan > 0 ? yy[0] : NULL,
			             devchan > 2 ? yy[2] : NULL,
			             devchan > 5 ? yy[5] : NULL,
			             devchan > 6 ? yy[6] : NULL,
			             devchan > 7 ? yy[7] : NULL,
			             devchan > 8 ? yy[8] : NULL,
			             devchan > 9 ? yy[9] : NULL,
			             PRES,
						 cx, cy, verify ? 0 : nmark);
		}
	}

	/* Compare the previous expected aDE against the current one */
	if (verify) {
		co tp;	/* Test point */
		double avg[MAX_CHAN];
		double max[MAX_CHAN];
		double rms[MAX_CHAN];
		int verified = 1;

		/* Verify each channel */
		for (j = 0; j < devchan; j++) {
			co tp;

			avg[j] = 0.0;
			max[j] = -1.0;
			rms[j] = 0.0;

			/* Sample it at GRES */
			for (i = 0; i < GRES; i++) {
				double iv, targ, val, tt;

				iv = i/(GRES-1.0);

				/* Lookup the ade that we expect */ 
				tp.p[0] = iv;
				pcade[j]->interp(pcade[j], &tp);
				targ = tp.v[0];

				/* Lookup the ade that we have */ 
				tp.p[0] = iv;
				ade[j]->interp(ade[j], &tp);
				val = tp.v[0];

				/* Compute the stats */
				tt = fabs(targ - val);
				avg[j] += tt;
				rms[j] += tt * tt;
				if (tt > max[j])
					max[j] = tt;
//printf("~1 chan %d, ix %d, iv %f, targ %f, actual %f, err %f\n",j,i,iv,targ,val,tt);
			}
			avg[j] /= (double)GRES;
			rms[j] /= (double)GRES;
			rms[j] = sqrt(rms[j]);
			if (max[j] > ver_maxde)
				verified = 0;
		}
		if (verb) {
			for (j = 0; j < devchan; j++) {
				printf("Verify results:\n");
				printf("Channel %d has DE avg %.1f, rms %.1f, max %.1f\n",j,avg[j],rms[j],max[j]);
			}
			if (verified)
				printf("Verified OK\n");
			else
				printf("Verification FAILED\n");
		}

		/* Plot the verification curves */
		if (doplot) {
			double xx[PRES];
			double yy[10][PRES];
	
			printf("Verification match plot:\n");

			for (j = 0; j < 6 && j < devchan; j++) {
				co tp;	/* Test point */
				double max;

				/* Establish the scale */
				tp.p[0] = 1.0;
				pcade[j]->interp(pcade[j], &tp);
				max = tp.v[0];

				for (i = 0; i < PRES; i++) {
					xx[i] = i/(double)(PRES-1);

					/* Convert ade target to device */
					tp.v[0] = max * xx[i];
					if ((tp.p[0] = rspl_ilookup(pcade[j], 1.0, tp.v[0])) < 0.0)
						error("Unexpected failure to invert curve %d for pcADE %f",j,tp.v[0]); 
					/* Convert device to actual ade */
					ade[j]->interp(ade[j], &tp);
					if (fabs(max) > 0.1)
						yy[j][i] = tp.v[0]/max;
					else
						yy[j][i] = 0.0;
				}
			}
			do_plot10(xx, devchan > 3 ? yy[3] : NULL,
						  devchan > 1 ? yy[1] : NULL,
						  devchan > 4 ? yy[4] : NULL,
						  devchan > 0 ? yy[0] : NULL,
						  devchan > 2 ? yy[2] : NULL,
						  devchan > 5 ? yy[5] : NULL,
		                  devchan > 6 ? yy[6] : NULL,
		                  devchan > 7 ? yy[7] : NULL,
		                  devchan > 8 ? yy[8] : NULL,
		                  devchan > 9 ? yy[9] : NULL,
			              PRES, 0);
			if (!verified)
				exit(1);
		}
	} else if (initial) {

		/* Convert any transfer curve target points into a smooth curve */ 
		if (pct->no_tpoints > 0) {
			int gres[MXDI] = { GRES };
			co *pnts;
		
			if ((pnts = (co *)calloc(pct->no_tpoints + 2, sizeof(co))) == NULL)
				error ("Malloc of rspl points failed");
	
			for (j = 0; j < devchan; j++) {
				int npts;
				int gotmin, gotmax;
	
				/* Count the number of valid points */
				for (npts = i = 0; i < pct->no_tpoints; i++) {
					if (pct->tpoints[i].val[j] >= 0.0)
						npts++;
				}
				if (npts == 0)
					continue;		/* No target curve for this channel */
	
				if ((tcurves[j] = new_rspl(RSPL_NOFLAGS, 1, 1)) == NULL)
					error("new_rspl(1,1) failed");
				
				gotmin = gotmax = 0;
				for (npts = i = 0; i < pct->no_tpoints; i++) {
					if (pct->tpoints[i].val[j] < 0.0)
						continue;
	
					pnts[npts].p[0] = pct->tpoints[i].loc;
					pnts[npts].v[0] = pct->tpoints[i].val[j];
					if (pnts[npts].p[0] < 0.0)
						pnts[npts].p[0] = 0.0;
					else if (pnts[npts].p[0] > 1.0)
						pnts[npts].p[0] = 1.0;
					if (pnts[npts].v[0] < 0.0)
						pnts[npts].v[0] = 0.0;
					else if (pnts[npts].v[0] > 1.0)
						pnts[npts].v[0] = 1.0;
	
					if (pnts[npts].p[0] <= 0.05)
						gotmin = 1;
					if (pnts[npts].p[0] >= 0.95)
						gotmax = 1;
					npts++;
				}
				/* Add default anchors if there are none supplied */
				if (gotmin == 0) {
					pnts[npts].p[0] = 0.0;
					pnts[npts++].v[0] = 0.0;
				}
				if (gotmax == 0) {
					pnts[npts].p[0] = 1.0;
					pnts[npts++].v[0] = 1.0;
				}
	
				/* Fit the curve to the given points */
				tcurves[j]->fit_rspl(tcurves[j], RSPLFLAGS, pnts, npts,
				                NULL, NULL, gres, NULL, NULL, TCURVESMOOTH, NULL, NULL); 
			}
			free(pnts);
	
			/* Plot the target curves */
			if (doplot) {
				double xx[PRES];
				double yy[10][PRES];
	
				printf("Target curves plot:\n");
	
				for (i = 0; i < (PRES-1); i++) {
					co tp;	/* Test point */
					double pp = i/(PRES-1.0);
	
					xx[i] = pp;
					for (j = 0; j < 10 && j < devchan; j++) {
						if (tcurves[j] != NULL) {
							tp.p[0] = pp;
							tcurves[j]->interp(tcurves[j], &tp);
							yy[j][i] = tp.v[0];
						} else
							yy[j][i] = pp;
					}
				}
				do_plot10(xx, devchan > 3 ? yy[3] : NULL,
							  devchan > 1 ? yy[1] : NULL,
							  devchan > 4 ? yy[4] : NULL,
							  devchan > 0 ? yy[0] : NULL,
							  devchan > 2 ? yy[2] : NULL,
							  devchan > 5 ? yy[5] : NULL,
		                      devchan > 6 ? yy[6] : NULL,
		                      devchan > 7 ? yy[7] : NULL,
		                      devchan > 8 ? yy[8] : NULL,
		                      devchan > 9 ? yy[9] : NULL,
				              PRES, 0);
			}
		}

		/* Do inverse lookup to create relative linearization curves */
		n_cvals = CAL_RES;
		for (j = 0; j < devchan; j++) {
			co tp;
			double rdemin, rdemax;	/* Relative DE min and max targets */

			/* Convert absolute de aims to relative */
			if ((rdemax = rspl_ilookup(ade[j], 0.0, pct->ademax[j])) < 0.0)
				error("Unexpected failure to invert curve %d for DE %f",j,pct->ademax[j]); 

			tp.p[0] = rdemax;
			rde[j]->interp(rde[j], &tp);
			rdemax = tp.v[0];

			if ((rdemin = rspl_ilookup(ade[j], 1.0, pct->ademin[j])) < 0.0)
				error("Unexpected failure to invert curve %d for DE %f",j,pct->ademax[j]); 

			tp.p[0] = rdemin;
			rde[j]->interp(rde[j], &tp);
			rdemin = tp.v[0];

			if (verb > 0)
				printf("Chan %d: rDE Max = %f, rDE Min = %f\n",j,rdemax,rdemin);

			if ((cvals[j] = (wval *)malloc(sizeof(wval) * n_cvals)) == NULL)
				error("Malloc of %d cvals failed",n_cvals);

			/* Convert relative delta E aim to device value */
			for (i = 0; i < n_cvals; i++) {
				double x = i/(n_cvals-1.0);
				double inv;

				cvals[j][i].inv = x;

				/* Apply any aim tweak curve */
				if (tcurves[j] != NULL) {
					tp.p[0] = x;
					tcurves[j]->interp(tcurves[j], &tp);
					x = tp.v[0];
				}

				inv = x * (rdemax - rdemin) + rdemin;

				if ((cvals[j][i].dev = rspl_ilookup(rde[j], 0.5, inv)) < 0.0)
					error("Unexpected failure to invert curve %d for DE %f",j,inv); 
//printf("~1 chan %d, step %d, inv %f, detarg %f, got dev %f\n",j,i,x,pp[0].v[0],pp[k].p[0]);
			}
		}

	} else if (recal || imitate) {

		n_cvals = CAL_RES;
		for (j = 0; j < devchan; j++) {
			co tp;

			if ((cvals[j] = (wval *)malloc(sizeof(wval) * n_cvals)) == NULL)
				error("Malloc of %d cvals failed",n_cvals);

			/* Lookup the expected ade for each input device value, and */
			/* then translate it into the required output device value */
			for (i = 0; i < n_cvals; i++) {
				double x = i/(n_cvals-1.0);

				cvals[j][i].inv = tp.p[0] = x;
				pcade[j]->interp(pcade[j], &tp);

				if ((cvals[j][i].dev = rspl_ilookup(ade[j], 0.5, tp.v[0])) < 0.0)
					error("Unexpected failure to invert curve %d for DE %f",j,tp.v[0]); 
//printf("~1 chan %d, ix %d, inv %f, pcade %f, iade %f\n",j,i,x,tp.v[0],cvals[j][i].dev);
			}
		}
	}

	if (initial || recal || imitate) {

		if (verb > 1) {
			printf("Calibration curve values:\n");
			for (i = 0; i < n_cvals; i++) {
				printf("  inv %f, dev",cvals[0][i].inv);
				for (j = 0; j < devchan; j++) {
					printf(" %f",cvals[j][i].dev);
				}
				printf("\n");
			}
		}

		/* Plot the calibration curves */
		if (doplot) {
			double xx[PRES];
			double yy[10][PRES];

			printf("Calibration curve plot:\n");

			for (i = 0; i < n_cvals; i++) {
				xx[i] = cvals[0][i].inv;
				for (j = 0; j < 10 && j < devchan; j++) {
					yy[j][i] = cvals[j][i].dev;
				}
			}
			do_plot10(xx, devchan > 3 ? yy[3] : NULL,	/* Black */
			              devchan > 1 ? yy[1] : NULL,	/* Red */
			              devchan > 4 ? yy[4] : NULL,	/* Green */
			              devchan > 0 ? yy[0] : NULL,	/* Blue */
			              devchan > 2 ? yy[2] : NULL,	/* Yellow */
			              devchan > 5 ? yy[5] : NULL,	/* Purple */
		                  devchan > 6 ? yy[6] : NULL,	/* Brown */
		                  devchan > 7 ? yy[7] : NULL,	/* Orange */
		                  devchan > 8 ? yy[8] : NULL,	/* Grey */
		                  devchan > 9 ? yy[9] : NULL,	/* White */
			              PRES, 0);
		}

		/* Write out an Argyll .CAL file */
		if (dowrite) {
			cgats *ocg;						/* output cgats structure */
			time_t clk = time(0);
			struct tm *tsp = localtime(&clk);
			char *atm = asctime(tsp);		/* Ascii time */
			char *ident = icx_inkmask2char(devmask, 1); 
			char *bident = icx_inkmask2char(devmask, 0); 
			cgats_set_elem *setel;			/* Array of set value elements */
			int nsetel = 0;
			int ncps;						/* Number of curve parameters */
			double *cps[3];					/* Arrays of curve parameters */
			char *bp = NULL, buf[100];		/* Buffer to sprintf into */
			co tp;

			ocg = new_cgats();				/* Create a CGATS structure */
			ocg->add_other(ocg, "CAL"); 	/* our special type is Calibration file */

			ocg->add_table(ocg, tt_other, 0);	/* Add a table for RAMDAC values */
			ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Device Calibration Curves",NULL);
			ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll printcal", NULL);
			atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
			ocg->add_kword(ocg, 0, "CREATED",atm, NULL);

			ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);
			ocg->add_kword(ocg, 0, "COLOR_REP", ident, NULL);

			if (xpi.deviceMfgDesc != NULL)
				ocg->add_kword(ocg, 0, "MANUFACTURER", xpi.deviceMfgDesc, NULL);
			if (xpi.modelDesc != NULL)
				ocg->add_kword(ocg, 0, "MODEL", xpi.modelDesc, NULL);
			if (xpi.profDesc != NULL)
				ocg->add_kword(ocg, 0, "DESCRIPTION", xpi.profDesc, NULL);
			if (xpi.copyright != NULL)
				ocg->add_kword(ocg, 0, "COPYRIGHT", xpi.copyright, NULL);

			/* Setup the table which holds the translation from calibrated */
			/* device value "I" to the raw device channel value */ 
			sprintf(buf, "%s_I",bident);
			ocg->add_field(ocg, 0, buf, r_t);
			nsetel++;
			for (j = 0; j < devchan; j++) {
				inkmask imask = icx_index2ink(devmask, j);
				sprintf(buf, "%s_%s",bident,icx_ink2char(imask));
				ocg->add_field(ocg, 0, buf, r_t);
				nsetel++;
			}
			if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * nsetel)) == NULL)
				error("Malloc failed!");

			/* Write the per channel device to device loolup curve values */
			if (devmask & ICX_ADDITIVE) {
				for (i = n_cvals-1; i >= 0; i--) {
		
					setel[0].d = 1.0 - cvals[0][i].inv;
					for (j = 0; j < devchan; j++)
						setel[1+j].d = 1.0 - cvals[j][i].dev;
					ocg->add_setarr(ocg, 0, setel);
				}
			} else {
				for (i = 0; i < n_cvals; i++) {
		
					setel[0].d = cvals[0][i].inv;
					for (j = 0; j < devchan; j++)
						setel[1+j].d = cvals[j][i].dev;
					ocg->add_setarr(ocg, 0, setel);
				}
			}

			free(setel);

			/* Write the calibration target information to a second table */
			if (pct->write(pct, ocg, 1) != 0)
				error("Writing cal target info to cal file '%s'",outname);

			/* Add a third table which is the expected absolute DE response */
			/* of the calibrated device. */
			ocg->add_table(ocg, tt_other, 0);	/* Add a table for RAMDAC values */
			ocg->add_kword(ocg, 2, "DESCRIPTOR", "Argyll Output Calibration Expected DE Response",NULL);
			ocg->add_kword(ocg, 2, "ORIGINATOR", "Argyll printcal", NULL);
			atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
			ocg->add_kword(ocg, 2, "CREATED",atm, NULL);

			ocg->add_kword(ocg, 2, "DEVICE_CLASS","OUTPUT", NULL);
			ocg->add_kword(ocg, 2, "COLOR_REP", ident, NULL);

			/* Setup the table which holds the translation from calibrated */
			/* device value "I" to the expected absolute deltaE value for */
			/* each colorant. */
			sprintf(buf, "%s_I_DE",bident);
			ocg->add_field(ocg, 2, buf, r_t);
			nsetel++;
			for (j = 0; j < devchan; j++) {
				inkmask imask = icx_index2ink(devmask, j);
				sprintf(buf, "%s_%s_DE",bident,icx_ink2char(imask));
				ocg->add_field(ocg, 2, buf, r_t);
				nsetel++;
			}
			if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * nsetel)) == NULL)
				error("Malloc failed!");

			for (i = 0; i < n_cvals; i++) {
				int ix = i;

				/* Calibrated device value */
				if (devmask & ICX_ADDITIVE) {
					ix = n_cvals -1 - i;
					setel[0].d = 1.0 - cvals[0][ix].inv;
				} else
					setel[0].d = cvals[0][ix].inv;

				for (j = 0; j < devchan; j++) {
					tp.p[0] = cvals[j][ix].dev;		/* Raw device value */
					ade[j]->interp(ade[j], &tp);	/* Corresponding ade value */
					setel[1+j].d = tp.v[0];
				}
				ocg->add_setarr(ocg, 2, setel);
			}

			free(setel);

			if (ocg->write_name(ocg, outname))
				error("Write error to file '%s': %s",outname,ocg->err);

			if (verb)
				printf("Written calibration file '%s'\n",outname);

			ocg->del(ocg);		/* Clean up */
			free(ident);
			free(bident);
		}

		/* 
			The structure of *.AMP is very simple.
	
				It has 5 tables that have
					256 entries that are 8-bit (byte), even for 16 bit mode.
			
						The first table 0000..00FFh is the CMYK channel.
						The second table 0100..01FFh is the C channel.
						The third table 0200..02FFh is the M channel.
						The fourth table 0300..03FFh is the Y channel.
						The fifth table 0300..03FFh is the K channel.
			
						Table position nn00h is the black end and table position nnFFh
						is the white end.
			
		*/

		/* Write an Adobe map format (.AMP) file */
		/* (It's not clear if more than 4 channels is allowed) */
		if (dowrite && doamp) {
			FILE *fp;
			cgatsFile *p;
			char nmode[50] = { '\000' };

			strcpy(nmode, "w");
#if defined(O_BINARY) || defined(_O_BINARY)
			strcat(nmode, "b");
#endif
			if ((fp = fopen(ampname, nmode)) == NULL)
				error("Couldn't open '%s' for writing",ampname);
			
			/* CMYK table is unity */
			for (i = 0; i < 256; i++) {
				if (putc(i,fp) == EOF)
					error("Error writing to fle '%s'",ampname);
			}
			for (j = 0; j < devchan; j++) {
				for (i = 0; i < 256; i++) {
					int x;
					if (devmask & ICX_ADDITIVE)
						x = (int)(cvals[j][i].dev * 255.0 + 0.5);		/* ??? */
					else
						x = 255 - (int)(cvals[j][255 - i].dev * 255.0 + 0.5);
					if (putc(x,fp) == EOF)
						error("Error writing to fle '%s'",ampname);
//printf("~1 chan %d, inv %d, dev %d\n",j,i,x);
				}
			}
			/* Extra 1:1 table */
			for (i = 0; i < 256; i++) {
				if (putc(i,fp) == EOF)
					error("Error writing to fle '%s'",ampname);
			}

			if (fclose(fp) != 0)
				error("Closing '%s' failed",ampname);

			if (verb)
				printf("Written calibration curves to '%s'\n",ampname);
		}
	}

	/* Free up various possible allocations */
	for (j = 0; j < devchan; j++) {
		if (pvals[j] != NULL)
			free(pvals[j]);
		if (raw[j] != NULL)
			raw[j]->del(raw[j]);
		if (ade[j] != NULL)
			ade[j]->del(ade[j]);
		if (rde[j] != NULL)
			rde[j]->del(rde[j]);
		if (pcade[j] != NULL)
			pcade[j]->del(pcade[j]);
		if (cvals[j] != NULL)
			free(cvals[j]);
		if (tcurves[j] != NULL)
			tcurves[j]->del(tcurves[j]);
	}

	return 0;
}










