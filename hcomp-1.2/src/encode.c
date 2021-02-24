/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* encode.c		encode H-transform and write to outfile
 *
 * Programmer: R. White		Date: 15 June 1994
 */
#include "hcompress.h"

int verbose;
int bitcount;

static char code_magic[2] = { (char)0xDD, (char)0x99 };

#ifdef __STDC__
static int arraymax(int a[], int nx, int ny, int ndim);
#else
static int arraymax();
#endif

#ifdef __STDC__
extern void encode(MYFILE *outfile, int a[], int nx, int ny, int scale)
#else
extern void encode(outfile,a,nx,ny,scale)
MYFILE *outfile;						/* Output file						*/
int a[];								/* input H-transform array (nx,ny)	*/
int nx,ny;								/* size of H-transform array		*/
int scale;								/* scale factor for digitization	*/
#endif
{
int nel, nx2, ny2, i, q, vmax[3], vtmp, nsign, bits_to_go;
unsigned char nbitplanes[3];
unsigned char *signbits, buffer;

	nel = nx*ny;
	/*
	 * write magic value
	 */
	qwrite(outfile, code_magic, sizeof(code_magic));
	writeint(outfile, nx);				/* size of image					*/
	writeint(outfile, ny);
	writeint(outfile, scale);			/* scale factor for digitization	*/
	/*
	 * write first value of A (sum of all pixels -- the only value
	 * which does not compress well)
	 */
	writeint(outfile, a[0]);
	a[0] = 0;
	/*
	 * allocate array for sign bits and save values, 8 per byte
	 */
	signbits = (unsigned char *) malloc((nel+7)/8);
	if (signbits == (unsigned char *) NULL) {
		fprintf(stderr, "encode: insufficient memory\n");
		exit(-1);
	}
	nsign = 0;
	bits_to_go = 8;
	buffer = 0;
	for (i = 0; i < nel; i++) {
		if (a[i] > 0) {
			/*
			 * positive element, put zero at end of buffer
			 */
			buffer <<= 1;
			if (--bits_to_go == 0) {
				/*
				 * filled up this byte, go to the next one
				 */
				signbits[nsign++] = buffer;
				bits_to_go = 8;
				buffer = 0;
			}
		} else if (a[i] < 0) {
			/*
			 * negative element, shift in a one
			 * replace a by absolute value
			 */
			a[i] = - a[i];
			buffer = (buffer<<1) | 1;
			if (--bits_to_go == 0) {
				signbits[nsign++] = buffer;
				bits_to_go = 8;
				buffer = 0;
			}
		}
	}
	if (bits_to_go != 8) {
		/*
		 * some bits in last element
		 * move bits in last byte to bottom and increment nsign
		 */
		signbits[nsign++] = buffer<<bits_to_go;
	}
	/*
	 * calculate number of bit planes for 3 quadrants
	 *
	 * quadrant 0=bottom left, 1=bottom right or top left, 2=top right, 
	 */
	for (q=0; q<3; q++) vmax[q] = 0;
	/*
	 * get maximum value in each quadrant
	 */
	nx2 = (nx+1)/2;
	ny2 = (ny+1)/2;

	vmax[0] = arraymax(  a,                  nx2,      ny2, ny);
	vmax[1] = arraymax( &a[ny2],             nx2, ny - ny2, ny);
	vtmp    = arraymax( &a[ny*nx2],     nx - nx2,      ny2, ny);
	vmax[2] = arraymax( &a[ny*nx2+ny2], nx - nx2, ny - ny2, ny);
	if (vtmp > vmax[1]) vmax[1] = vtmp;

	/*
	 * now calculate number of bits for each quadrant
	 */
	for (q = 0; q < 3; q++) {
		nbitplanes[q] = log((float) (vmax[q]+1))/log(2.0)+0.5;
		if ( (vmax[q]+1) > (int) (1 << nbitplanes[q]) ) {
			nbitplanes[q] += 1;
		}
	}
	/*
	 * write nbitplanes
	 */
	qwrite(outfile, (char *) nbitplanes, sizeof(nbitplanes));
	/*
	 * write coded array
	 */
	doencode(outfile, a, nx, ny, nbitplanes);
	/*
	 * write sign bits
	 */
	if (nsign > 0) qwrite(outfile, (char *) signbits, nsign);
	if (verbose) {
		/*
		 * total number of bits written to file
		 */
		i=bitcount+
			8*(nsign+sizeof(code_magic)+sizeof(nbitplanes)+4*sizeof(int));
		fprintf(stderr, "%6.3f bits/pixel, compression factor %5.1f\n",
			((float) i)/nel,
			16.0*nel/((float) i));
	}
	free(signbits);
}

#ifdef __STDC__
static int arraymax(int a[], int nx, int ny, int ndim)
#else
static int arraymax(a, nx, ny, ndim)
int a[];
int nx;
int ny;
int ndim;
#endif
{
int i, *p, amax;
	amax = 0;
	for (i=0; i<nx; i++) {
		for (p = &a[i*ndim]; p < &a[i*ndim+ny]; p++) {
			if (*p > amax) amax = *p;
		}
	}
	return (amax);
}
