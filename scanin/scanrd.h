
/* 
 * Argyll Color Correction System
 *
 * Scanrd: Scan chart reader
 * This is the core chart recognition code.
 *
 * Author: Graeme W. Gill
 * Date:   28/9/96
 *
 * Copyright 1995, 1996, 2008, Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * Public Interface and object file for scanrd.c
 */

/* Operation flags */
#define SI_BUILD_REF 	      0x10000	/* Build the reference file */
#define SI_PERSPECTIVE 	      0x20000	/* Allow perspective correction */
#define SI_GENERAL_ROT 	      0x40000	/* Allow general rotation, else assume zero degrees */
#define SI_ASISIFFAIL 	      0x80000	/* Read patch values "as is" if everything else failes */

/* Scanrd diagnostic flags */
#define SI_SHOW_FLAGS         0xffff	/* Mask for all SHOW flags */
#define SI_SHOW_IMAGE	      0x0001	/* Show B&W version of input image in output */
#define SI_SHOW_DIFFSH        0x0002	/* Show the horizontal edges detected */
#define SI_SHOW_DIFFSV	      0x0004	/* Show the vertical edges detected  */
#define SI_SHOW_GROUPS	      0x0008	/* Show the groups detected */
#define SI_SHOW_LINES         0x0010	/* Show the lines detected */
#define SI_SHOW_PERS          0x0020	/* Show the lines un-perspective */
#define SI_SHOW_ROT           0x0040	/* Show the lines rotated */
#define SI_SHOW_IMPL          0x0080	/* Show the lines used for improvements */
#define SI_SHOW_ALL_LINES     0x0100	/* Show all lines, valid and invalid */
#define SI_SHOW_SBOX          0x0200	/* Show aligned sample box info in diagnostic raster */
#define SI_SHOW_SBOX_OUTLINES 0x0400	/* Show the sample box outlines */
#define SI_SHOW_SBOX_NAMES    0x0800	/* Show sample boxes names */
#define SI_SHOW_SBOX_AREAS    0x1000	/* Show sample boxes sample areas */
#define SI_SHOW_SAMPLED_AREA  0x2000	/* Show pixels areas sampled */

/* Error flags */
#define SI_QUAL_ERR(flag) ((flag & 0xf0000000) == 0)
#define SI_FIND_PERSPECTIVE_FAILED   0x00000001
#define SI_FIND_ROTATION_FAILED      0x00000002
#define SI_POOR_MATCH                0x00000003

#define SI_FILE_ERR(flag) ((flag & 0xf0000000) == 0x10000000)
#define SI_RAST_READ_ERR             0x10000001
#define SI_DIAG_WRITE_ERR            0x10000002
#define SI_REF_WRITE_ERR             0x10000003
#define SI_REF_READ_ERR              0x10000004
#define SI_REF_FORMAT_ERR            0x10000005
#define SI_PIX_DEPTH_ERR             0x10000006
#define SI_BIT_DEPTH_ERR             0x10000007
#define SI_NO_FIDUCIALS_ERR          0x10000008
#define SI_BAD_FIDUCIALS_ERR         0x10000009		/* Not really file error */

#define SI_MALLOC_ERR(flag) ((flag & 0xf0000000) == 0x80000000)
#define SI_MALLOC_DIAG_RAST          0x80000001
#define SI_MALLOC_INPUT_BUF          0x80000002
#define SI_MALLOC_POINT2LINE         0x80000003
#define SI_MALLOC_ELIST              0x80000004
#define SI_MALLOC_REFREAD            0x80000005
#define SI_MALLOC_SETUP_BOXES        0x80000006
#define SI_MALLOC_VALUE_SCAN         0x80000007
#define SI_MALLOC_VREGION            0x80000008
#define SI_MALLOC_POINTS             0x80000009
#define SI_REALLOC_POINTS            0x8000000A
#define SI_MALLOC_AAINIT             0x8000000B

#define SI_INTERNAL_ERR(flag) ((flag & 0xf0000000) == 0xA0000000)
#define SI_INTERNAL                  0xA0000001

/* The public scanrd object */
struct _scanrd {
	/*** Public methods ***/
	/* Initialise, ready to read out all the values */
	/* return the total number of values */
	int (*reset)(struct _scanrd *s);

	/* Read the next samples values */
	/* return non-zero when no more points */
	int (*read)(struct _scanrd *s,
		char *id,									/* patch id copied to here */
		double *P,		/* Robust mean values */
		double *mP,		/* Raw Mean values */
		double *sdP,	/* Standard deviation */
		int *cnt);									/* Pixel count */

	/* Return the error flag, and set the message pointer */
	unsigned int (*error)(struct _scanrd *s, char **errm);

	/* Free up the structure */
	void (*free)(struct _scanrd *s);

	}; typedef struct _scanrd scanrd;


/* Read in a chart */
/* Then use reset() and read() to get values read */
scanrd *do_scanrd(
	int flags,			/* option flags */
	int verb,			/* verbosity level */

	double gamma,		/* Approximate gamma encoding of image (0.0 = default 1.7) */

	double *sfid,		/* Specified fiducuals x1,y1, x2,y2, x3,y3, x4,y4,  */
						/* Typically clockwise from top left, NULL if auto recognition */

	int w, int h, 		/* Width and Height of input raster in pixels */
	int d, int td, int p,	/* Useful plane depth, Total depth, Bit presision of input pixels */

	int (*read_line)(void *fdata, int y, char *dst),	/* Read pixel interleaved line of source */
	void *fdata,		/* Opaque data for read_line */

	char *refname,		/* reference file name */

	int (*write_line)(void *ddata, int y, char *src),	/* Write 8bpp RGB line of diag file */
	void *ddata			/* Opaque data for write_line */
);


