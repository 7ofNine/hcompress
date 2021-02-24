/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* hdecomp.c	Decompress image file that was compressed using hcomp.
 *
 * This version of the program can produce output in 4 formats:
 *
 * (1) raw: I*2 pixel values with bytes in machine-dependent order,
 *			i.e. no byte swapping is done on output.  Image size may
 *			be specified with the -r and -c parameters.
 *
 * (2) net: I*2 pixel values with bytes with bytes in "network" order:
 *			high byte first, then low byte for each I*2 pixel.  Byte-swapping
 *			is done on output if needed.  Note that this is the same as
 *			raw format on some machines (e.g. Suns) but is different on
 *			others (e.g. VAXes).  Files in net format can be transferred
 *			from one machine to another without modification, but files
 *			in raw format cannot.
 *
 * (3) fits: FITS format image.  Header gives image size.  Bytes are
 *			in network order.
 *
 * (4) hhh:	This is identical to raw format, but is included for consistency
 *			with the input formats used by hcomp.
 *
 * The compressed input file may have two different formats, which are
 * recognized automatically by the program.  If the compressed file was
 * produced from a FITS input file, then the compressed file includes the
 * (uncompressed) FITS header at the start of the file, followed by the
 * compressed image data.  Otherwise the compressed input file contains
 * only the compressed image data.
 *
 * If the compressed file has a FITS header, then the default output format
 * is fits.  Otherwise the default output format is raw.
 *
 * Programmer: R. White		Date: 30 June 1994
 */
#include "hcompress.h"

extern int  optind;
extern char *optarg;

/*
 * a[nx][ny] is the image array (gets allocated in decode)
 * Note that ny is the fast-varying dimension
 *
 * scale is the scale factor for digitization
 */
int  verbose;
static int  *a, nx, ny;
static int  scale;
static int  smooth;
static char *format;

#ifdef VMS
/*
 * buffers to speed VMS IO
 */
#define BUFSIZE 64000
static char inbuf[BUFSIZE+1];
static char outbuf[BUFSIZE+1];

#endif

#ifdef __STDC__
static void usage(char *argv[]);
static void get_args(int argc, char *argv[]);
#else
static void usage();
static void get_args();
#endif

#ifdef __STDC__
int main (int argc, char *argv[])
#else
main (argc,argv)
int argc;
char *argv[];
#endif
{
#ifdef VMS
	/*
	 * set buffer sizes to improve IO speeds on VMS
	 */
	if (setvbuf(stdin, inbuf, _IOFBF, BUFSIZE) != 0) {
		perror("setvbuf(stdin) failed");
	}
	if (setvbuf(stdout, outbuf, _IOFBF, BUFSIZE) != 0) {
		perror("setvbuf(stdout) failed");
	}
#endif
	/*
	 * Get command line arguments
	 */
	get_args(argc, argv);
	/*
	 * Read from stdin, passing header through to stdout for FITS format,
	 * and decode.  Returns address, size, scale,
	 * and (possibly) format
	 */
	decode(stdin,stdout,&a,&nx,&ny,&scale,&format);
	/*
	 * Inverse H-transform
	 */
	hinv(a,nx,ny,smooth,scale);
	/*
	 * Write data
	 */
	put_data(stdout,a,nx,ny,format);
	free(a);
	if (verbose) {
		if (smooth) {
			fprintf(stderr,
		"Image size (%d,%d)  Scale factor %d  Smoothed  Output in %s format\n",
				ny,nx,scale,format);
		} else {
			fprintf(stderr,
		"Image size (%d,%d)  Scale factor %d  Output in %s format\n",
				ny,nx,scale,format);
		}
	}
	return(0);
}

#ifdef __STDC__
static void usage(char *argv[])
#else
static void usage(argv)
char *argv[];
#endif
{
	fprintf(stderr, "%s version %s\n", PACKAGE, VERSION);
	fprintf(stderr,
		"Usage: %s [-v] [-s] [-o raw|net|fits|hhh]\n",
		argv[0]);
	exit(-1);
}

/* GET COMMAND LINE ARGUMENTS */

#ifdef __STDC__
static void get_args(int argc, char *argv[])
#else
static void get_args(argc,argv)
int argc;
char *argv[];
#endif
{
int c;

	/*
	 * default values
	 */
	verbose = 0;
	smooth = 0;
	format = "";
	/*
	 * get options
	 */
	while ((c = getopt(argc,argv,"vso:h")) != -1) {
		switch (c) {
		case 'v':
			/*
			 * verbose flag -v
			 */
			verbose = 1;
			break;
		case 'h':
			/*
			 * help flag -h
			 */
			usage(argv);
			break;
		case 's':
			/*
			 * smoothing flag -s
			 */
			smooth = 1;
			break;
		case 'o':
			/*
			 * -o <format> = raw, net, fits, or hhh
			 */
			format = optarg;
			if (strcmp(format,"raw")  != 0 &&
				strcmp(format,"net")  != 0 &&
				strcmp(format,"fits") != 0 &&
				strcmp(format,"hhh")  != 0) {
				fprintf(stderr, "Illegal input format %s\n", format);
				usage(argv);
			}
			break;
		case '?':
			usage(argv);
		}
	}
	/*
	 * make sure there aren't any trailing parameters being ignored
	 */
	if (optind < argc) {
		fprintf(stderr, "Too many parameters: %s ...\n", argv[optind]);
		usage(argv);
	}
}
