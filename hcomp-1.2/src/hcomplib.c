/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 *
 * hcompress: compress image using hcompress algorithm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "config.h"

int errno;

#define FAILURE (1)

/* MYFILE structure for compression to a file
 * Defined this way so same code can be used for file or buffer.
 */
typedef FILE MYFILE;
#define mygetc getc
#define myputc putc

/* function prototypes */

/* hcompress routines */
static void htrans(int a[], int nx, int ny, int *status);
static void digitize(int a[], int nx, int ny, int scale, int *status);
static void encode(MYFILE *outfile, int a[], int nx, int ny, int scale, int verbose, int *status);
static int doencode(MYFILE *outfile, int a[], int nx, int ny, unsigned char nbitplanes[3], int *status);
static void qtree_encode(MYFILE *outfile, int a[], int n, int nqx, int nqy, int nbitplanes, int *status);

/* data I/O routines */
static void writeint(MYFILE *outfile, int a, int *status);
static void qwrite(MYFILE *outfile, char *a, int n, int *status);
static int mywrite(MYFILE *file, char buffer[], int n);

/* bit I/O */
static void start_outputing_bits(void);
static void output_bit(MYFILE *outfile, int bit);
static void output_nbits(MYFILE *outfile, int bits, int n);
static int done_outputing_bits(MYFILE *outfile);

/* secondary routines */
static void xshuffle(int a[], int nx, int ny, int nydim, int *status);
static void yshuffle(int a[], int nx, int ny, int nydim, int *status);
static int arraymax(int a[], int nx, int ny, int ndim);
static void qtree_onebit(int a[], int n, int nx, int ny, unsigned char b[], int bit);
static void qtree_reduce(unsigned char a[], int n, int nx, int ny, unsigned char b[]);
static int  bufcopy(unsigned char a[], int n, unsigned char **bptr, unsigned char *bufend);
static void write_bdirect(MYFILE *outfile, unsigned char a[], int nqx, int nqy);
static void printerr(const char *format, ...);

extern int hcompress(int *a, int nx, int ny, int scale, int verbose, MYFILE *outfile, int *status)
{
	if (*status != 0) return(*status);
	/*
	 * H-transform
	 */
	htrans(a,nx,ny,status);
	/*
	 * Digitize
	 */
	digitize(a,nx,ny,scale,status);
	/*
	 * Encode and write to stdout
	 */
	encode(outfile,a,nx,ny,scale,verbose,status);
	return(*status);
}

/* htrans.c   H-transform of NX x NY integer image
 *
 * Programmer: R. White		Date: 13 June 1994
 */

static void htrans(int a[], int nx, int ny, int *status)
{
int nmax, log2n, nxtop, nytop, i, k;
int shift, mask, mask2, prnd, prnd2, nrnd2;
int h0, hx, hy, hc, sum1, sum2, sum3, sum4;
int *p00, *p10, *pend;

	if (*status != 0) return;
	/*
	 * log2n is log2 of max(nx,ny) rounded up to next power of 2
	 */
	nmax = (nx>ny) ? nx : ny;
	log2n = log((float) nmax)/log(2.0)+0.5;
	if ( nmax > (1<<log2n) ) {
		log2n += 1;
	}
	/*
	 * set up rounding and shifting masks
	 */
	shift = 0;
	mask  = -2;
	mask2 = mask << 1;
	prnd  = 1;
	prnd2 = prnd << 1;
	nrnd2 = prnd2 - 1;
	/*
	 * do log2n reductions
	 *
	 * We're indexing a as a 2-D array with dimensions (nx,ny).
	 */
	nxtop = nx;
	nytop = ny;
	for (k = 0; k<log2n; k++) {
		for (i = 0; i<nxtop-1; i += 2) {
			pend = &a[i*ny+nytop-1];
			for (p00 = &a[i*ny], p10 = p00+ny; p00<pend; p00 += 2, p10 += 2) {
				/*
				 * Divide h0,hx,hy,hc by 2 (1 the first time through).
				 */
				sum1 = *(p10+1) + *p10;
				sum2 = *(p00+1) + *p00;
				sum3 = *(p10+1) - *p10;
				sum4 = *(p00+1) - *p00;
				if (shift==0) {
					h0 = (sum1 + sum2);
					hx = (sum1 - sum2);
					hy = (sum3 + sum4);
					hc = (sum3 - sum4);
				} else {
					h0 = (sum1 + sum2) >> 1;
					hx = (sum1 - sum2) >> 1;
					hy = (sum3 + sum4) >> 1;
					hc = (sum3 - sum4) >> 1;
				}
				/*
				 * Throw away the 2 bottom bits of h0, bottom bit of hx,hy.
				 * To get rounding to be same for positive and negative
				 * numbers, nrnd2 = prnd2 - 1.
				 */
				*(p10+1) = hc;
				*(p10  ) = ( (hx>=0) ? (hx+prnd)  :  hx        ) & mask ;
				*(p00+1) = ( (hy>=0) ? (hy+prnd)  :  hy        ) & mask ;
				*(p00  ) = ( (h0>=0) ? (h0+prnd2) : (h0+nrnd2) ) & mask2;
			}
			if (p00 == pend) {
				/*
				 * do last element in row if row length is odd
				 * (p00+1), (p10+1) are off edge
				 */
				h0 = (*p10 + *p00) << (1-shift);
				hx = (*p10 - *p00) << (1-shift);
				*p10 = ( (hx>=0) ? (hx+prnd)  :  hx        ) & mask ;
				*p00 = ( (h0>=0) ? (h0+prnd2) : (h0+nrnd2) ) & mask2;
			}
		}
		if (i<nxtop) {
			/*
			 * do last row if column length is odd
			 * p10, (p10+1) are off edge
			 */
			pend = &a[i*ny+nytop-1];
			for (p00 = &a[i*ny]; p00<pend; p00 += 2) {
				h0 = (*(p00+1) + *p00) << (1-shift);
				hy = (*(p00+1) - *p00) << (1-shift);
				*(p00+1) = ( (hy>=0) ? (hy+prnd)  :  hy        ) & mask ;
				*(p00  ) = ( (h0>=0) ? (h0+prnd2) : (h0+nrnd2) ) & mask2;
			}
			if (p00==pend) {
				/*
				 * do corner element if both row and column lengths are odd
				 * (p00+1), p10, (p10+1) are off edge
				 */
				h0 = *p00 << (2-shift);
				*p00 = ( (h0>=0) ? (h0+prnd2) : (h0+nrnd2) ) & mask2;
			}
		}
		/*
		 * now shuffle in each dimension to group coefficients by order
		 */
		xshuffle(a,nxtop,nytop,ny,status);
		if (*status != 0) return;
		yshuffle(a,nxtop,nytop,ny,status);
		if (*status != 0) return;
		/*
		 * image size reduced by 2 (round up if odd)
		 */
		nxtop = (nxtop+1)>>1;
		nytop = (nytop+1)>>1;
		/*
		 * divisor doubles after first reduction
		 */
		shift = 1;
		/*
		 * masks, rounding values double after each iteration
		 */
		mask  = mask2;
		prnd  = prnd2;
		mask2 = mask2 << 1;
		prnd2 = prnd2 << 1;
		nrnd2 = prnd2 - 1;
	}
}

/*
 * shuffle in x direction to get even elements at beginning of
 * array, odd elements at end.  procedure is to copy odd elements
 * to temporary array, compress even elements to beginning of row,
 * and this copy odd elements to end of array.
 */

static void xshuffle(int a[], int nx, int ny, int nydim, int *status)
{
int j, *p1, *p2, *pt, *pend, *tmp;

	if (*status != 0) return;
	/*
	 * get temporary storage for shuffling elements
	 */
	tmp  = (int *) malloc(((ny+1)/2)*sizeof(int));
	if(tmp == (int *) NULL) {
		printerr("htrans: insufficient memory\n");
		*status = FAILURE;
		return;
	}
	for (j = 0; j<nx; j++) {
		/*
		 * copy odd elements to tmp
		 */
		pend = &a[nydim*j+ny-1];
		for (pt = tmp, p1 = &a[nydim*j+1]; p1 <= pend; p1 += 2, pt += 1)
			*pt = *p1;
		/*
		 * compress even elements into first half of A
		 */
		for (p1 = &a[nydim*j+1], p2 = &a[nydim*j+2]; p2 <= pend; p1 += 1, p2 += 2)
			*p1 = *p2;
		/*
		 * put odd elements into 2nd half
		 */
		(void) memcpy(p1, tmp, (ny/2)*sizeof(int));
	}
	free(tmp);
}

/*
 * shuffle in y direction to get even elements at beginning of
 * array, odd elements at end.  This is done using a somewhat
 * complicated method for efficiency.  The straightforward
 * approach is slow because the scattered memory accesses don't
 * take advantage of the cache (I think.)  This version does
 * operations a row at a time so that consecutive memory locations
 * are accessed.
 */

static void yshuffle(int a[], int nx, int ny, int nydim, int *status)
{
int j, k, tt, oddoffset, *tmp;
int *p, *pt;
unsigned char *flag;

	if (*status != 0) return;
	/*
	 * get temporary storage for shuffling elements
	 */
	tmp  =           (int *) malloc(ny*sizeof(int));
	flag = (unsigned char *) malloc(nx*sizeof(unsigned char));
	if(tmp == (int *) NULL || flag == (unsigned char *) NULL) {
		printerr("htrans: insufficient memory\n");
		*status = FAILURE;
		return;
	}
	/*
	 * initialize flag array telling whether row is done
	 */
	for (j=0; j<nx; j++) flag[j] = 1;
	oddoffset = (nx+1)/2;
	/*
	 * shuffle each row to appropriate location
	 * row 0 is already in right location
	 */
	for (j=1; j<nx; j++) {
		if (flag[j]) {
			flag[j] = 0;
			/*
			 * where does this row belong?
			 * second factor = 0 for even rows, oddoffset for odd rows
			 */
			k = (j>>1) + ((j & 1) ? oddoffset : 0);
			if (j != k) {
				/*
				 * copy the row
				 */
				(void) memcpy(tmp,&a[nydim*j],ny*sizeof(int));
				/*
				 * keep shuffling until we reach a row that is done
				 */
				while (flag[k]) {
					flag[k] = 0;
					/*
					 * do the exchange
					 */
					for (p = &a[nydim*k], pt=tmp; p < &a[nydim*k+ny]; p++, pt++) {
						tt = *p;
						*p = *pt;
						*pt = tt;
					}
					k = (k>>1) + ((k & 1) ? oddoffset : 0);
				}
				/*
				 * copy the last row into place
				 * this should always end up with j=k
				 */
				(void) memcpy(&a[nydim*k],tmp,ny*sizeof(int));
				if (j != k) {
					printerr("error: yshuffle failed!\nj=%d k=%d\n", j, k);
					*status = FAILURE;
					return;
				}
			}
		}
	}
	free(tmp);
	free(flag);
}

/* digitize.c	digitize H-transform
 *
 * Programmer: R. White		Date: 15 June 1994
 */

static void digitize(int a[], int nx, int ny, int scale, int *status)
{
int d, *p;

	if (*status != 0) return;
	/*
	 * round to multiple of scale
	 */
	if (scale <= 1) return;
	d=(scale+1)/2-1;
	if (d == 0) {
		for (p=a; p <= &a[nx*ny-1]; p++) *p = *p / scale;
	} else {
		for (p=a; p <= &a[nx*ny-1]; p++) {
			if (*p > 0) {
				*p = (*p+d) / scale;
			} else {
				*p = (*p-d) / scale;
			}
		}
	}
}

/* encode.c		encode H-transform and write to outfile
 *
 * Programmer: R. White		Date: 15 June 1994
 */

static char code_magic[2] = { (char)0xDD, (char)0x99 };

static void encode(MYFILE *outfile, int a[], int nx, int ny, int scale, int verbose, int *status)
{
int nel, nx2, ny2, i, q, vmax[3], vtmp, nsign, lbits_to_go, lbitcount;
unsigned char nbitplanes[3];
unsigned char *signbits, lbitbuffer;

	if (*status != 0) return;
	nel = nx*ny;
	/*
	 * write magic value
	 */
	qwrite(outfile, code_magic, sizeof(code_magic), status);
	writeint(outfile, nx, status);				/* size of image					*/
	writeint(outfile, ny, status);
	writeint(outfile, scale, status);			/* scale factor for digitization	*/
	/*
	 * write first value of A (sum of all pixels -- the only value
	 * which does not compress well)
	 */
	writeint(outfile, a[0], status);
	if (*status != 0) return;
	a[0] = 0;

	/*
	 * allocate array for sign bits and save values, 8 per byte
	 */
	signbits = (unsigned char *) malloc((nel+7)/8);
	if (signbits == (unsigned char *) NULL) {
		printerr("encode: insufficient memory\n");
		*status = FAILURE;
		return;
	}
	nsign = 0;
	lbits_to_go = 8;
	lbitbuffer = 0;
	for (i = 0; i < nel; i++) {
		if (a[i] > 0) {
			/*
			 * positive element, put zero at end of buffer
			 */
			lbitbuffer <<= 1;
			if (--lbits_to_go == 0) {
				/*
				 * filled up this byte, go to the next one
				 */
				signbits[nsign++] = lbitbuffer;
				lbits_to_go = 8;
				lbitbuffer = 0;
			}
		} else if (a[i] < 0) {
			/*
			 * negative element, shift in a one
			 * replace a by absolute value
			 */
			a[i] = - a[i];
			lbitbuffer = (lbitbuffer<<1) | 1;
			if (--lbits_to_go == 0) {
				signbits[nsign++] = lbitbuffer;
				lbits_to_go = 8;
				lbitbuffer = 0;
			}
		}
	}
	if (lbits_to_go != 8) {
		/*
		 * some bits in last element
		 * move bits in last byte to bottom and increment nsign
		 */
		signbits[nsign++] = lbitbuffer<<lbits_to_go;
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
	qwrite(outfile, (char *) nbitplanes, sizeof(nbitplanes), status);
	if (*status != 0) return;
	/*
	 * write coded array
	 */
	lbitcount = doencode(outfile, a, nx, ny, nbitplanes, status);
	if (*status != 0) return;
	/*
	 * write sign bits
	 */
	if (nsign > 0) {
		qwrite(outfile, (char *) signbits, nsign, status);
		if (*status != 0) return;
	}
	if (verbose) {
		/*
		 * total number of bits written to file
		 */
		lbitcount += 8*(nsign+sizeof(code_magic)+sizeof(nbitplanes)+4*sizeof(int));
		printerr("%6.3f bits/pixel, compression factor %5.1f\n",
			((float) lbitcount)/nel,
			16.0*nel/((float) lbitcount));
	}
	free(signbits);
}

static int arraymax(int a[], int nx, int ny, int ndim)
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

/* doencode.c	Encode 2-D array and write stream of characters on outfile
 *
 * This version assumes that A is positive.
 *
 * Returns total number of bits written
 *
 * Programmer: R. White		Date: 7 May 1991
 */


static int doencode(MYFILE *outfile, int a[], int nx, int ny, unsigned char nbitplanes[3], int *status)
{
int nx2, ny2;

	if (*status != 0) return 0;
	nx2 = (nx+1)/2;
	ny2 = (ny+1)/2;
	/*
	 * Initialize bit output
	 */
	start_outputing_bits();
	/*
	 * write out the bit planes for each quadrant
	 */
	qtree_encode(outfile, &a[0],          ny, nx2,  ny2,  nbitplanes[0], status);
	qtree_encode(outfile, &a[ny2],        ny, nx2,  ny/2, nbitplanes[1], status);
	qtree_encode(outfile, &a[ny*nx2],     ny, nx/2, ny2,  nbitplanes[1], status);
	qtree_encode(outfile, &a[ny*nx2+ny2], ny, nx/2, ny/2, nbitplanes[2], status);
	if (*status != 0) return 0;
	/*
	 * Add zero nybble as an EOF symbol
	 */
	output_nbits(outfile, 0, 4);
	return done_outputing_bits(outfile);
}

/* qtree_encode.c	Encode values in quadrant of 2-D array using binary
 *					quadtree coding for each bit plane.  Assumes array is
 *					positive.
 *
 * Programmer: R. White		Date: 14 June 1994
 */

/*
 * Huffman code values and number of bits in each code
 */
static int code[16] =
	{
	0x3e, 0x00, 0x01, 0x08, 0x02, 0x09, 0x1a, 0x1b,
	0x03, 0x1c, 0x0a, 0x1d, 0x0b, 0x1e, 0x3f, 0x0c
	};
static int ncode[16] =
	{
	6,    3,    3,    4,    3,    4,    5,    5,
	3,    5,    4,    5,    4,    5,    6,    4
	};

/*
 * variables for bit output to buffer when Huffman coding
 */
static int bitbuffer2;
static int bits_to_go2;

static void qtree_encode(MYFILE *outfile, int a[], int n, int nqx, int nqy, int nbitplanes, int *status)
{
int log2n, k, bit, bmax, nqmax, nqx2, nqy2, nx, ny;
unsigned char *scratch, *scr1, *sinput, *buffer, *bptr, *bufend;

	if (*status != 0) return;
	/*
	 * log2n is log2 of max(nqx,nqy) rounded up to next power of 2
	 */
	nqmax = (nqx>nqy) ? nqx : nqy;
	log2n = log((float) nqmax)/log(2.0)+0.5;
	if (nqmax > (1<<log2n)) {
		log2n += 1;
	}
	/*
	 * initialize buffer point, max buffer size
	 */
	nqx2 = (nqx+1)/2;
	nqy2 = (nqy+1)/2;
	bmax = (nqx2*nqy2+1)/2;
	/*
	 * We're indexing A as a 2-D array with dimensions (nqx,nqy).
	 * Scratch is 2-D with dimensions (nqx/2,nqy/2) rounded up.
	 * Scr1 is used to store first level of quadtree in case direct
	 * coding is needed.
	 * Buffer is used to store string of codes for output.
	 */
	scratch = (unsigned char *) malloc(2*bmax);
	scr1    = (unsigned char *) malloc(2*bmax);
	buffer  = (unsigned char *) malloc(bmax);
	if ((scratch == (unsigned char *) NULL) ||
		(buffer  == (unsigned char *) NULL) ||
		(scr1    == (unsigned char *) NULL)) {
		printerr("qtree_encode: insufficient memory\n");
		*status = FAILURE;
		return;
	}
	bufend = &buffer[bmax];
	/*
	 * now encode each bit plane, starting with the top
	 */
	for (bit=nbitplanes-1; bit >= 0; bit--) {
		/*
		 * initialize bit buffer
		 */
		bptr = buffer;
		bitbuffer2 = 0;
		bits_to_go2 = 0;
		/*
		 * on first pass copy A to scr1 array
		 */
		qtree_onebit(a,n,nqx,nqy,scr1,bit);
		nx = (nqx+1)>>1;
		ny = (nqy+1)>>1;
		/*
		 * copy non-zero values to output buffer, which will be written
		 * in reverse order
		 */
		if (bufcopy(scr1,nx*ny,&bptr,bufend)) {
			/*
			 * quadtree is expanding data,
			 * change warning code and just fill buffer with bit-map
			 */
			write_bdirect(outfile,scr1,nqx,nqy);
			goto bitplane_done;
		}
		/*
		 * do log2n reductions
		 * sinput is pointer to input scratch array; it changes to
		 * scratch after the first pass.
		 */
		sinput = scr1;
		for (k = 1; k<log2n; k++) {
			qtree_reduce(sinput,ny,nx,ny,scratch);
			sinput = scratch;
			nx = (nx+1)>>1;
			ny = (ny+1)>>1;
			if (bufcopy(scratch,nx*ny,&bptr,bufend)) {
				write_bdirect(outfile,scr1,nqx,nqy);
				goto bitplane_done;
			}
		}
		/*
		 * OK, we've got the code in buffer
		 * Write quadtree warning code, then write buffer in reverse order
		 */
		output_nbits(outfile,0xF,4);
		if (bptr == buffer) {
			if (bits_to_go2>0) {
				/*
				 * put out the last few bits
				 */
				output_nbits(outfile, bitbuffer2 & ((1<<bits_to_go2)-1),
					bits_to_go2);
			} else {
				/*
				 * have to write a zero nybble if there are no 1's in array
				 */
				output_nbits(outfile,code[0],ncode[0]);
			}
		} else {
			if (bits_to_go2>0) {
				/*
				 * put out the last few bits
				 */
				output_nbits(outfile, bitbuffer2 & ((1<<bits_to_go2)-1),
					bits_to_go2);
			}
			/*
			 * write in blocks of 24 bits to speed things up
			 */
			for (bptr--; bptr >= buffer+2; bptr -= 3)
				output_nbits(outfile,
					(*bptr<<16) | (*(bptr-1)<<8) | *(bptr-2),24);
			for ( ; bptr >= buffer; bptr--)
				output_nbits(outfile,*bptr,8);
		}
		bitplane_done: ;
	}
	free(buffer);
	free(scratch);
	free(scr1);
}

/*
 * copy non-zero codes from array to buffer
 */

static int bufcopy(unsigned char a[], int n, unsigned char **bptr, unsigned char *bufend)
{
unsigned char *p, *buffer;
/* local copies of global variables */
int lbitbuffer;
int lbits_to_go;

	buffer = *bptr;
	lbits_to_go = bits_to_go2;
	lbitbuffer = bitbuffer2;
	for (p = a; p < &a[n]; p++) {
		if (*p != 0) {
			/*
			 * add Huffman code for a[i] to buffer
			 */
			lbitbuffer |= code[*p] << lbits_to_go;
			lbits_to_go += ncode[*p];
			if (lbits_to_go >= 24) {
				/*
				 * return warning code if we're going to fill buffer
				 */
				if (buffer+3 >= bufend) return(1);
				/*
				 * move 3 bytes to buffer
				 */
				*buffer++ =  lbitbuffer      & 0xFF;
				*buffer++ = (lbitbuffer>> 8) & 0xFF;
				*buffer++ = (lbitbuffer>>16) & 0xFF;
				lbitbuffer >>= 24;
				lbits_to_go -= 24;
			}
		}
	}
	*bptr = buffer;
	bitbuffer2 = lbitbuffer;
	bits_to_go2 = lbits_to_go;
	return(0);
}

/*
 * Do first quadtree reduction step on bit BIT of array A.
 * Results put into B.
 * 
 */

static void qtree_onebit(int a[], int n, int nx, int ny, unsigned char b[], int bit)
{
int i, *p, *pend;
unsigned char *pb, *pb0;
int mbit, bitm2;

	/*
	 * mask to get selected bit
	 */
	mbit = 1<<bit;
	pb = b;
	bitm2 = bit - 2;
	for (i = 0; i<nx; i += 2) {
		pb0 = pb;
		pend = &a[n*i+ny-1];
		switch (bit) {
			case 0:
				for (p = &a[n*i]; p < pend; p += 2, pb += 1)
					*pb = ((*p & mbit)<<3) | ((*(p+1) & mbit)<<2);
				break;
			case 1:
				for (p = &a[n*i]; p < pend; p += 2, pb += 1)
					*pb = ((*p & mbit)<<2) | ((*(p+1) & mbit)<<1);
				break;
			case 2:
				for (p = &a[n*i]; p < pend; p += 2, pb += 1)
					*pb = ((*p & mbit)<<1) | ((*(p+1) & mbit)   );
				break;
			default:
				for (p = &a[n*i]; p < pend; p += 2, pb += 1)
					*pb = ( ((*p & mbit)<<1) | (*(p+1) & mbit) ) >> bitm2;
		}
		if (p == pend) {
			/*
			 * row size is odd, do last element in row
			 * *(p+1) is off edge
			 */
			*pb = ((*p & mbit)<<3) >> bit;
			pb += 1;
		}
		if (i < nx-1) {
			/*
			 * not on last row, add in next row
			 */
			pb = pb0;
			pend = &a[n*(i+1)+ny-1];
			switch (bit) {
				case 0:
					for (p = &a[n*(i+1)]; p < pend; p += 2, pb += 1)
						*pb = ( ((*p & mbit)<<1) |  (*(p+1) & mbit)     ) | *pb;
					break;
				case 1:
					for (p = &a[n*(i+1)]; p < pend; p += 2, pb += 1)
						*pb = (  (*p & mbit)     | ((*(p+1) & mbit)>>1) ) | *pb;
					break;
				default:
					for (p = &a[n*(i+1)]; p < pend; p += 2, pb += 1)
						*pb = ( ( ((*p & mbit)<<1) | (*(p+1) & mbit) ) >> bit) | *pb ;
					break;
			}
			if (p == pend) {
				/* odd row size */
				*pb = ( ((*p & mbit)<<1) >> bit) | *pb ;
				pb += 1;
			}
		}
	}
}

/*
 * do one quadtree reduction step on array a
 * results put into b (which may be the same as a)
 */

static void qtree_reduce(unsigned char a[], int n, int nx, int ny, unsigned char b[])
{
int i;
unsigned char *p, *pend, *pb, *pb0;

	pb = b;
	for (i = 0; i<nx; i += 2) {
		pb0 = pb;
		pend = &a[n*i+ny-1];
		for (p = &a[n*i]; p < pend; p += 2, pb += 1)
			*pb = ( ((*p != 0) << 3) | ((*(p+1) != 0) << 2) );
		if (p == pend) {
			/*
			 * row size is odd, do last element in row
			 * *(p+1) is off edge
			 */
			*pb = (*p != 0) << 3;
			pb += 1;
		}
		if (i < nx-1) {
			/*
			 * not on last row, add in next row
			 */
			pb = pb0;
			pend = &a[n*(i+1)+ny-1];
			for (p = &a[n*(i+1)]; p < pend; p += 2, pb += 1)
				*pb = ((*p != 0) << 1) | (*(p+1) != 0) | *pb;
			if (p == pend) {
				/* odd row size */
				*pb = ((*p != 0) << 1) | *pb;
				pb += 1;
			}
		}
	}
}

static void write_bdirect(MYFILE *outfile, unsigned char a[], int nqx, int nqy)
{
int i;
int ilast;

	/*
	 * Write the direct bitmap warning code
	 */
	output_nbits(outfile,0x0,4);
	/*
	 * write quads for original bitplane to outfile
	 * write in blocks of 24 bits to speed things up
	 */
	ilast = ((nqx+1)/2) * ((nqy+1)/2);
	for (i = 0; i < ilast-5; i += 6) {
		output_nbits(outfile,
			(a[i  ]<<20) |
			(a[i+1]<<16) |
			(a[i+2]<<12) |
			(a[i+3]<< 8) |
			(a[i+4]<< 4) |
			 a[i+5],
			24);
	}
	for (; i < ilast; i++) output_nbits(outfile,a[i],4);
}

/* qwrite.c Write binary data
 *
 * Programmer: R. White     Date: 11 March 1991
 */

static void writeint(MYFILE *outfile, int a, int *status)
{
int i;
unsigned char b[4];

	if (*status != 0) return;
	/* Write integer A one byte at a time to outfile.
	 *
	 * This is portable from Vax to Sun since it eliminates the
	 * need for byte-swapping.
	 */
	for (i=3; i>=0; i--) {
		b[i] = a & 0x000000ff;
		a >>= 8;
	}
	for (i=0; i<4; i++) {
		qwrite(outfile,(char *) &b[i],1,status);
	}
	return;
}

static void qwrite(MYFILE *outfile, char *a, int n, int *status)
{
	if (*status != 0) return;
	if(mywrite(outfile, a, n) != n) {
		/* perror("qwrite"); */
		printerr("qwrite: %s\n",strerror(errno));
		*status = FAILURE;
	}
}

/*
 * write n bytes from buffer into file
 * returns number of bytes written (=n) if successful, <=0 if not
 */

static int mywrite(MYFILE *file, char buffer[], int n)
{
#ifdef TO_A_BUFFER
int nmax;
	/*
	 * this version used when doing I/O to a buffer
	 */
	if (file->current+n > file->end) {
		nmax = file->end - file->current;
		if (nmax > 0) (void) memcpy(file->current, buffer, nmax);
		file->current = file->end;
		printerr("mywrite: buffer overflow\n");
		return (-1);
	}
	(void) memcpy(file->current, buffer, n);
	file->current += n;
	return (n);
#else
	/*
	 * this version used when doing I/O to a file
	 */
	return ( fwrite(buffer, 1, n, file) );
#endif
}

/* Bit output routines */

/* Bit buffer */

static int bitbuffer1;				/* Bits buffered for output			*/
static int bits_to_go1;				/* Number of bits free in bitbuffer	*/
static int bitcount;				/* Count of total bits written		*/

/* Initialize for bit output */

static void start_outputing_bits(void)
{
	bitbuffer1 = 0;					/* Buffer is empty to start	*/
	bits_to_go1 = 8;				/* with						*/
	bitcount = 0;
}


/* Output a bit */

static void output_bit(MYFILE *outfile, int bit)
{
	bitbuffer1 <<= 1;								/* Put bit at end of buffer */
	if (bit) bitbuffer1 |= 1;
	bits_to_go1 -= 1;
	if (bits_to_go1 == 0) {							/* Output buffer if it is	*/
		(void) myputc(bitbuffer1 & 0xff,outfile);	/* now full					*/
		bits_to_go1 = 8;
		bitbuffer1 = 0;
		bitcount += 8;
	}
}


/* Output N bits (N must be <= 24) */

static void output_nbits(MYFILE *outfile, int bits, int n)
{
/* local copies */
int lbitbuffer;
int lbits_to_go;

	/*
	 * insert bits at end of buffer
	 */
	lbitbuffer = bitbuffer1;
	lbits_to_go = bits_to_go1;
	lbitbuffer <<= n;
	lbitbuffer |= ( bits & ((1<<n)-1) );
	lbits_to_go -= n;
	while (lbits_to_go <= 0) {
		/*
		 * buffer full, put out top 8 bits
		 */
		(void) myputc((lbitbuffer>>(-lbits_to_go)) & 0xff,outfile);
		lbits_to_go += 8;
		bitcount += 8;
	}
	bitbuffer1 = lbitbuffer;
	bits_to_go1 = lbits_to_go;
}


/*
 * Flush out the last bits
 * Returns total number of bits written
 */

static int done_outputing_bits(MYFILE *outfile)
{
	if(bits_to_go1 < 8) {
		(void) myputc(bitbuffer1<<bits_to_go1,outfile);
		/* count the garbage bits too */
		bitcount += 8;
	}
	return bitcount;
}

/* Print error messages */

static void printerr(const char *format, ...)
{
va_list args;

	va_start (args, format);
	vfprintf(stderr, format, args);
	va_end(args);
}
