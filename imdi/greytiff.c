
/* 
 * Convert a TIFF to monochrome in a colorimetrically correct way.
 *
 * Author:  Graeme W. Gill
 * Date:    01/8/29
 * Version: 1.00
 *
 * Copyright 2000, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * Thanks to Neil Okamoto for the 16 bit TIFF mods.
 */

/* TTBD:
 *
 */

/*

	This program is a framework that exercises the
	IMDI code.  It can also do the conversion using the
    floating point code in ICCLIB as a reference.

 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "aconfig.h"
#include "tiffio.h"
#include "icc.h"
#include "numlib.h"
#include "xicc.h"
#include "imdi.h"

#undef DO_CHECK		/* Do floating point check */

void usage(void) {
	fprintf(stderr,"Convert a TIFF file to monochrome using an ICC device profile, V%s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: greytiff [-v level] profile.icm infile.tif outfile.tif\n");
	fprintf(stderr," -v            Verbose\n");
	fprintf(stderr," -p            Use slow precise correction\n");
	fprintf(stderr," -j            Use CIECAM02\n");
	exit(1);
}

/* Convert an ICC colorspace to the corresponding TIFF Photometric tag */
/* return 0xffff if not possible. */

int
ColorSpaceSignature2TiffPhotometric(
icColorSpaceSignature cspace
) {
	switch(cspace) {
		case icSigGrayData:
			return PHOTOMETRIC_MINISBLACK;
		case icSigRgbData:
			return PHOTOMETRIC_RGB;
		case icSigCmykData:
			return PHOTOMETRIC_SEPARATED;
		case icSigYCbCrData:
			return PHOTOMETRIC_YCBCR;
		case icSigLabData:
			return PHOTOMETRIC_CIELAB;

		case icSigXYZData:
		case icSigLuvData:
		case icSigYxyData:
		case icSigHsvData:
		case icSigHlsData:
		case icSigCmyData:
		case icSig2colorData:
		case icSig3colorData:
		case icSig4colorData:
		case icSig5colorData:
		case icSigMch5Data:
		case icSig6colorData:
		case icSigMch6Data:
		case icSig7colorData:
		case icSigMch7Data:
		case icSig8colorData:
		case icSigMch8Data:
		case icSig9colorData:
		case icSig10colorData:
		case icSig11colorData:
		case icSig12colorData:
		case icSig13colorData:
		case icSig14colorData:
		case icSig15colorData:
		default:
			return 0xffff;
	}
	return 0xffff;
}

char *
Photometric2str(
int pmtc
) {
	static char buf[80];
	switch (pmtc) {
		case PHOTOMETRIC_MINISWHITE:
			return "Subtractive Gray";
		case PHOTOMETRIC_MINISBLACK:
			return "Additive Gray";
		case PHOTOMETRIC_RGB:
			return "RGB";
		case PHOTOMETRIC_PALETTE:
			return "Indexed";
		case PHOTOMETRIC_MASK:
			return "Transparency Mask";
		case PHOTOMETRIC_SEPARATED:
			return "CMYK";
		case PHOTOMETRIC_YCBCR:
			return "YCbCr";
		case PHOTOMETRIC_CIELAB:
			return "CIELab";
		case PHOTOMETRIC_LOGL:
			return "CIELog2L";
		case PHOTOMETRIC_LOGLUV:
			return "CIELog2Luv";
	}
	sprintf(buf,"Unknonw Tag %d",pmtc);
	return buf;
}

/* Callbacks used to initialise imdi */

/* Context for imdi setup callbacks */
typedef struct {
	int id, od;
	icxLuBase *flu;		/* Device -> Jab/Lab */
	icxLuBase *blu;		/* Jab/Lab -> Device */
} sucntx;

/* Input curve function */
void input_curve(
	void *cntx,
	double *out_vals,
	double *in_vals
) {
	sucntx *rx = (sucntx *)cntx;
	int e;

	for (e = 0; e < rx->id; e++) 
		out_vals[e] = in_vals[e];
}

/* Multi-dim table function */
void md_table(
void *cntx,
double *out_vals,
double *in_vals
) {
	sucntx *rx = (sucntx *)cntx;
	double Lab[3];
	
	rx->flu->lookup(rx->flu, Lab, in_vals);
	Lab[1] = Lab[2] = 0.0;
	rx->blu->lookup(rx->blu, out_vals, Lab);
}


/* Output curve function */
void output_curve(
	void *cntx,
	double *out_vals,
	double *in_vals
) {
	sucntx *rx = (sucntx *)cntx;
	int e;

	for (e = 0; e < rx->od; e++) 
		out_vals[e] = in_vals[e];
}

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char prof_name[100];
	char in_name[100];
	char out_name[100];

	icmFile *p_fp;
	icc *icco;
	xicc *xicco;
	int verb = 0;
	int slow = 0;
	int rv = 0;
	icColorSpaceSignature pcsor = icSigLabData;

	TIFF *rh = NULL, *wh = NULL;
	int x, y, width, height;					/* Size of image */
	uint16 samplesperpixel, bitspersample;
	uint16 pconfig, photometric, pmtc;
	uint16 resunits;
	float resx, resy;
	tdata_t *inbuf, *outbuf, *checkbuf;

	icColorSpaceSignature ins;		/* Type of input spaces */
	int inn;						/* Number of device components */

	/* IMDI */
	imdi *s = NULL;
	sucntx su;		/* Setup context */
	unsigned char *inp[MAX_CHAN];
	unsigned char *outp[MAX_CHAN];

#ifdef DO_CHECK
	/* Error check */
	int mxerr = 0;
	double avgerr = 0.0;
	double avgcount = 0.0;
#endif /* DO_CHECK */

	if (argc < 2)
		usage();

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

			/* Slow, Precise */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				slow = 1;
			}

			/* Use CIECAM02 */
			else if (argv[fa][1] == 'j' || argv[fa][1] == 'J') {
				pcsor = icxSigJabData;
			}

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			}

			else 
				usage();
		} else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(prof_name,argv[fa++]);

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(in_name,argv[fa++]);

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(out_name,argv[fa++]);

	/* - - - - - - - - - - - - - - - - */
	/* Open up the profile for reading */
	if ((p_fp = new_icmFileStd_name(prof_name,"r")) == NULL)
		error ("Can't open file '%s'",prof_name);

	if ((icco = new_icc()) == NULL)
		error ("Creation of ICC object failed");

	/* Wrap with an expanded icc */
	if ((xicco = new_xicc(icco)) == NULL)
		error ("Creation of xicc failed");

	if ((rv = icco->read(icco,p_fp,0)) != 0)
		error ("%d, %s",rv,icco->err);

	if (verb) {
		icmFile *op;
		if ((op = new_icmFileStd_fp(stdout)) == NULL)
			error ("Can't open stdout");
		icco->header->dump(icco->header, op, 1);
		op->del(op);
	}

	/* Check that the profile is appropriate */
	if (icco->header->deviceClass != icSigInputClass
	 && icco->header->deviceClass != icSigDisplayClass
	 && icco->header->deviceClass != icSigOutputClass
	 && icco->header->deviceClass != icSigColorSpaceClass)	/* For sRGB etc. */
		error("Profile isn't a device profile");

	/* Get a expanded color conversion object */
	if ((su.flu = xicco->get_luobj(xicco, ICX_CLIP_NEAREST, icmFwd, icRelativeColorimetric, pcsor, icmLuOrdNorm, NULL, NULL)) == NULL)
		error ("%d, %s",xicco->errc, xicco->err);

	/* Get details of conversion (Arguments may be NULL if info not needed) */
	su.flu->spaces(su.flu, &ins, &inn, NULL, NULL, NULL, NULL, NULL, NULL);

	su.id = inn;
	su.od = inn;

	/* Get a bwd conversion object */
	if ((su.blu = xicco->get_luobj(xicco, ICX_CLIP_NEAREST, icmBwd, icRelativeColorimetric, pcsor, icmLuOrdNorm, NULL, NULL)) == NULL)
		error ("%d, %s",xicco->errc, xicco->err);

	/* - - - - - - - - - - - - - - - */
	/* Open up input tiff file ready for reading */
	/* Got arguments, so setup to process the file */
	if ((rh = TIFFOpen(in_name, "r")) == NULL)
		error("error opening read file '%s'",in_name);

	TIFFGetField(rh, TIFFTAG_IMAGEWIDTH,  &width);
	TIFFGetField(rh, TIFFTAG_IMAGELENGTH, &height);

	TIFFGetField(rh, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	if (bitspersample != 8 && bitspersample != 16) {
		error("TIFF Input file must be 8 or 16 bit/channel");
	}

	TIFFGetField(rh, TIFFTAG_PHOTOMETRIC, &photometric);
	if  ((pmtc = ColorSpaceSignature2TiffPhotometric(ins)) == 0xffff)
		error("ICC  input colorspace '%s' can't be handled by a TIFF file!",
		      icm2str(icmColorSpaceSignature, ins));
	if (pmtc != photometric)
		error("ICC  input colorspace '%s' doesn't match TIFF photometric '%s'!",
		      icm2str(icmColorSpaceSignature, ins), Photometric2str(photometric));

	TIFFGetField(rh, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
	if (inn != samplesperpixel)
		error ("TIFF Input file has %d input channels mismatched to colorspace '%s'",
		       samplesperpixel, icm2str(icmColorSpaceSignature, ins));

	TIFFGetField(rh, TIFFTAG_PLANARCONFIG, &pconfig);
	if (pconfig != PLANARCONFIG_CONTIG)
		error ("TIFF Input file must be planar");

	TIFFGetField(rh, TIFFTAG_RESOLUTIONUNIT, &resunits);
	TIFFGetField(rh, TIFFTAG_XRESOLUTION, &resx);
	TIFFGetField(rh, TIFFTAG_YRESOLUTION, &resy);

	/* - - - - - - - - - - - - - - - */
	if ((wh = TIFFOpen(out_name, "w")) == NULL)
		error("Can\'t create TIFF file '%s'!",out_name);
	
	TIFFSetField(wh, TIFFTAG_IMAGEWIDTH,  width);
	TIFFSetField(wh, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(wh, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(wh, TIFFTAG_SAMPLESPERPIXEL, inn);
	TIFFSetField(wh, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(wh, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	if  ((pmtc = ColorSpaceSignature2TiffPhotometric(ins)) == 0xffff)
		error("TIFF file can't handle output colorspace '%s'!",
		      icm2str(icmColorSpaceSignature, ins));
	TIFFSetField(wh, TIFFTAG_PHOTOMETRIC, pmtc);
	TIFFSetField(wh, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	if (resunits) {
		TIFFSetField(wh, TIFFTAG_RESOLUTIONUNIT, resunits);
		TIFFSetField(wh, TIFFTAG_XRESOLUTION, resx);
		TIFFSetField(wh, TIFFTAG_YRESOLUTION, resy);
	}
	TIFFSetField(wh, TIFFTAG_IMAGEDESCRIPTION, "Color corrected by Argyll");

	/* - - - - - - - - - - - - - - - */
	/* Setup the imdi */

	if (!slow) {
		s = new_imdi(
			inn,			/* Number of input dimensions */
			inn,			/* Number of output dimensions */
			bitspersample == 8 ? pixint8 : pixint16,
							/* Output pixel representation */
			0x0,			/* Treat every channel as unsigned */
			NULL,			/* No raster to callback channel mapping */
			prec_min,		/* Minimum of input and output precision */
			bitspersample == 8 ? pixint8 : pixint16,
			0x0,			/* Treat every channel as unsigned */
			NULL,			/* No raster to callback channel mapping */
			17,				/* Desired table resolution. 33 is also a good number */
			oopts_none,		/* Desired per channel output options */
			NULL,			/* Output channel check values */
			opts_none,		/* Desired processing direction and stride support */
			input_curve,	/* Callback functions */
			md_table,
			output_curve,
			(void *)&su		/* Context to callbacks */
		);
	
		if (s == NULL)
			error("new_imdi failed");
	}

	/* - - - - - - - - - - - - - - - */
	/* Process colors to translate */
	/* (Should fix this to process a group of lines at a time ?) */

	inbuf  = _TIFFmalloc(TIFFScanlineSize(rh));
	outbuf = _TIFFmalloc(TIFFScanlineSize(wh));
	checkbuf = _TIFFmalloc(TIFFScanlineSize(wh));

	inp[0] = (unsigned char *)inbuf;
	outp[0] = (unsigned char *)outbuf;

	if (!slow) {		/* Fast */
		for (y = 0; y < height; y++) {

			/* Read in the next line */
			if (TIFFReadScanline(rh, inbuf, y, 0) < 0)
				error ("Failed to read TIFF line %d",y);

			/* Do fast conversion */
			s->interp(s, (void **)outp, 0, (void **)inp, 0, width);
			
#ifdef DO_CHECK
			/* Do floating point conversion */
			for (x = 0; x < width; x++) {
				int i;
				double in[MAX_CHAN], out[MAX_CHAN];
				double Lab[3];
				
				if (bitspersample == 8)
					for (i = 0; i < inn; i++)
						in[i] = ((unsigned char *)inbuf)[x * inn + i]/255.0;
				else
					for (i = 0; i < inn; i++)
						in[i] = ((unsigned short *)inbuf)[x * inn + i]/65535.0;
					
				if ((rv = su.flu->lookup(su.flu, Lab, in)) > 1)
					error ("%d, %s",icco->errc,icco->err);

				Lab[1] = Lab[2] = 0.0;

				if ((rv = su.blu->lookup(su.blu, out, Lab)) > 1)
					error ("%d, %s",icco->errc,icco->err);

				if (bitspersample == 8)
					for (i = 0; i < inn; i++)
						((unsigned char *)checkbuf)[x * inn + i] = (int)(out[i] * 255.0 + 0.5);
				else
					for (i = 0; i < inn; i++)
						((unsigned short *)checkbuf)[x * inn + i] = (int)(out[i] * 65535.0 + 0.5);
			}
			/* Compute the errors */
			for (x = 0; x < (width * inn); x++) {
				int err;

				if (bitspersample == 8)
					err = ((unsigned char *)outbuf)[x] - ((unsigned char *)checkbuf)[x];
				else
					err = ((unsigned short *)outbuf)[x] - ((unsigned short *)checkbuf)[x];
				if (err < 0)
					err = -err;
				if (err > mxerr)
					mxerr = err;
				avgerr += (double)err;
				avgcount++;
			}
#endif /* DO_CHECK */
				
			if (TIFFWriteScanline(wh, outbuf, y, 0) < 0)
				error ("Failed to write TIFF line %d",y);

		}

	} else {	/* Slow but precise */
		if (bitspersample == 8) {
			for (y = 0; y < height; y++) {

				/* Read in the next line */
				if (TIFFReadScanline(rh, inbuf, y, 0) < 0)
					error ("Failed to read TIFF line %d",y);

				/* Do floating point conversion */
				for (x = 0; x < width; x++) {
					int i;
					double in[MAX_CHAN], out[MAX_CHAN];
					double Lab[3];
					
					for (i = 0; i < inn; i++) {
						in[i] = ((unsigned char *)inbuf)[x * inn + i]/255.0;
					}
					
					if ((rv = su.flu->lookup(su.flu, Lab, in)) > 1)
						error ("%d, %s",icco->errc,icco->err);
	
					Lab[1] = Lab[2] = 0.0;
	
					if ((rv = su.blu->lookup(su.blu, out, Lab)) > 1)
						error ("%d, %s",icco->errc,icco->err);

					for (i = 0; i < inn; i++) {
						((unsigned char *)outbuf)[x * inn + i] = (int)(out[i] * 255.0 + 0.5);
					}
				}
				if (TIFFWriteScanline(wh, outbuf, y, 0) < 0)
					error ("Failed to write TIFF line %d",y);
			}
		} else if (bitspersample == 16) {
			for (y = 0; y < height; y++) {

				/* Read in the next line */
				if (TIFFReadScanline(rh, inbuf, y, 0) < 0)
					error ("Failed to read TIFF line %d",y);

				/* Do floating point conversion */
				for (x = 0; x < width; x++) {
					int i;
					double in[MAX_CHAN], out[MAX_CHAN];
					double Lab[3];
					
					for (i = 0; i < inn; i++) {
						in[i] = ((unsigned short *)inbuf)[x * inn + i]/65535.0;
					}
					
					if ((rv = su.flu->lookup(su.flu, Lab, in)) > 1)
						error ("%d, %s",icco->errc,icco->err);
	
					Lab[1] = Lab[2] = 0.0;
	
					if ((rv = su.blu->lookup(su.blu, out, Lab)) > 1)
						error ("%d, %s",icco->errc,icco->err);

					for (i = 0; i < inn; i++) {
						((unsigned short *)outbuf)[x * inn + i] = (int)(out[i] * 65535.0 + 0.5);
					}
				}
				if (TIFFWriteScanline(wh, outbuf, y, 0) < 0)
					error ("Failed to write TIFF line %d",y);
			}
		}
	}


#ifdef DO_CHECK
	printf("Worst error = %d bits, average error = %f bits\n", mxerr, avgerr/avgcount);
	if (bitspersample == 8)
		printf("Worst error = %f%%, average error = %f%%\n",
		       mxerr/2.55, avgerr/(2.55 * avgcount));
	else
		printf("Worst error = %f%%, average error = %f%%\n",
		       mxerr/655.35, avgerr/(655.35 * avgcount));
#endif /* DO_CHECK */

	/* Done with lookup object */
	if (s != NULL)
		s->del(s);
	su.flu->del(su.flu);
	su.blu->del(su.blu);
	xicco->del(xicco);
	icco->del(icco);
	p_fp->del(p_fp);

	TIFFClose(rh);		/* Close Input file */
	TIFFClose(wh);		/* Close Output file */

	return 0;
}

