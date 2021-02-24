/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 *
 * hdecompress: decompress image using hcompress algorithm
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

/*
 * address buffer type (could be unsigned char if image size < 512 pixels)
 */
typedef short ADDRTYPE;

/* function prototypes */

/* hcompress routines */
static void hinv(int a[], int nx, int ny, int smooth, int scale, int *status);
static void hsmooth(int a[], int nxtop, int nytop, int ny, int scale);
static void decode(MYFILE *infile, int  **a, int  *nx, int *ny, int  *scale, int *status);
static void dodecode(MYFILE *infile, int a[], int nx, int ny, unsigned char nbitplanes[3], int scale, int *status);
static void qtree_decode(MYFILE *infile, int a[], int n, int nqx, int nqy, int nbitplanes, int *status);

/* data I/O routines */
static int readint(MYFILE *infile, int *status);
static void qread(MYFILE *infile, char *a, int n, int *status);
static int myread(MYFILE *file, char buffer[], int n);

/* bit I/O */
static void start_inputing_bits(void);
static int input_bit(MYFILE *infile, int *status);
static int input_nbits(MYFILE *infile, int n, int *status);

/* secondary routines */

static void xunshuffle(int a[], int nx, int ny, int nydim, int *status);
static void yunshuffle(int a[], int nx, int ny, int nydim, int *status);
static int qaddr_expand(MYFILE *infile, ADDRTYPE *xthis, ADDRTYPE *ythis,
	int nt, ADDRTYPE *xnext, ADDRTYPE *ynext, ADDRTYPE *xscr, ADDRTYPE *yscr);
static void qaddr_bitins(MYFILE *infile, ADDRTYPE *xthis, ADDRTYPE *ythis,
	int nt, int a[], int n, int bit);
static void read_bdirect(MYFILE *infile, int b[], int n, int nx,
	int ny, int bit);
static int input_huffman(MYFILE *infile);
static void erreof(const char *str, int *status);
static void printerr(const char *format, ...);

extern int hdecompress(int **a, int *nx, int *ny, int *scale, int smooth, MYFILE *infile, int *status)
{
	if (*status != 0) return(*status);
	/*
	 * Read from stdin, passing header through to stdout for FITS format,
	 * and decode.  Returns address, size, scale,
	 * and (possibly) format
	 */
	decode(infile,a,nx,ny,scale,status);
	/*
	 * Inverse H-transform
	 */
	hinv(*a,*nx,*ny,smooth,*scale,status);

	return(*status);
}

/* hinv.c   Inverse H-transform of NX x NY integer image
 *
 * Programmer: R. White		Date: 16 June 1994
 */

/* smooth = 0 for no smoothing, else smooth during inversion
 * scale is used if smoothing is specified
 */
static void hinv(int a[], int nx, int ny, int smooth, int scale, int *status)
{
int nmax, log2n, i, k;
int nxtop,nytop,nxf,nyf,c;
int bit0, bit1, bit2, mask0, mask1, mask2,
	prnd0, prnd1, prnd2, nrnd0, nrnd1, nrnd2, lowbit0, lowbit1;
int h0, hx, hy, hc, sum1, sum2;
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
	 * set up masks, rounding parameters
	 */
	bit0   = 1 << (log2n - 1);
	bit1   = bit0 << 1;
	bit2   = bit0 << 2;
	mask0  = -bit0;
	mask1  = mask0 << 1;
	mask2  = mask0 << 2;
	prnd0  = bit0 >> 1;
	prnd1  = bit1 >> 1;
	prnd2  = bit2 >> 1;
	nrnd0  = prnd0 - 1;
	nrnd1  = prnd1 - 1;
	nrnd2  = prnd2 - 1;
	/*
	 * round h0 to multiple of bit2
	 */
	a[0] = (a[0] + ((a[0] >= 0) ? prnd2 : nrnd2)) & mask2;
	/*
	 * do log2n expansions
	 *
	 * We're indexing a as a 2-D array with dimensions (nx,ny).
	 */
	nxtop = 1;
	nytop = 1;
	nxf = nx;
	nyf = ny;
	c = 1<<log2n;
	for (k = log2n-1; k>0; k--) {
		/*
		 * this somewhat cryptic code generates the sequence
		 * ntop[k-1] = (ntop[k]+1)/2, where ntop[log2n] = n
		 */
		c = c>>1;
		nxtop = nxtop<<1;
		nytop = nytop<<1;
		if (nxf <= c) { nxtop -= 1; } else { nxf -= c; }
		if (nyf <= c) { nytop -= 1; } else { nyf -= c; }
		/*
		 * unshuffle in each dimension to interleave coefficients
		 */
		xunshuffle(a,nxtop,nytop,ny,status);
		yunshuffle(a,nxtop,nytop,ny,status);
		/*
		 * smooth by interpolating coefficients if SMOOTH != 0
		 */
		if (smooth) hsmooth(a,nxtop,nytop,ny,scale);
		for (i = 0; i<nxtop-1; i += 2) {
			pend = &a[ny*i+nytop-1];
			for (p00 = &a[ny*i], p10 = p00+ny; p00 < pend; p00 += 2, p10 += 2) {
				h0 = *(p00  );
				hx = *(p10  );
				hy = *(p00+1);
				hc = *(p10+1);
				/*
				 * round hc to multiple of bit0
				 */
				hc = (hc + ((hc >= 0) ? prnd0 : nrnd0)) & mask0;
				/*
				 * round hx and hy to multiple of bit1
				 * h0 is already a multiple of bit2
				 * propagate bit0 of hc to hx,hy
				 */
				lowbit0 = hc & bit0;
				if (lowbit0==0) {
					/*
					 * lowbit0 = 0
					 */
					if (hx >= 0) {
						hx = ((hx + prnd1) & mask1);
					} else {
						hx = ((hx + nrnd1) & mask1);
					}
					if (hy >= 0) {
						hy = ((hy + prnd1) & mask1);
					} else {
						hy = ((hy + nrnd1) & mask1);
					}
					/*
					 * Propagate bits 0 and 1 of hc,hx,hy to h0.
					 * This could be simplified if we assume h0>0, but then
					 * the inversion would not be lossless for images with
					 * negative pixels.
					 */
					lowbit1 = (hc ^ hx ^ hy) & bit1;
					h0 = (h0 >= 0) ? (h0 - lowbit1) : (h0 + lowbit1);
				} else {
					/*
					 * lowbit0 != 0
					 */
					if (hx >= 0) {
						hx = ((hx + prnd1) & mask1) - lowbit0;
					} else {
						hx = ((hx + nrnd1) & mask1) + lowbit0;
					}
					if (hy >= 0) {
						hy = ((hy + prnd1) & mask1) - lowbit0;
					} else {
						hy = ((hy + nrnd1) & mask1) + lowbit0;
					}
					lowbit1 = (hc ^ hx ^ hy) & bit1;
					h0 = h0 + lowbit0 - lowbit1;
				}
				/*
				 * Divide sums by 2
				 */
				sum1 = h0+hx;
				sum2 = hy+hc;
				*(p10+1) = (sum1 + sum2) >> 1;
				*(p10  ) = (sum1 - sum2) >> 1;
				sum1 = h0-hx;
				sum2 = hy-hc;
				*(p00+1) = (sum1 + sum2) >> 1;
				*(p00  ) = (sum1 - sum2) >> 1;
			}
			if (p00==pend) {
				/*
				 * do last element in row if row length is odd
				 * p00+1, p10+1 are off edge
				 */
				h0 = *(p00  );
				hx = *(p10  );
				if (hx >= 0) {
					hx = ((hx + prnd1) & mask1);
				} else {
					hx = ((hx + nrnd1) & mask1);
				}
				lowbit1 = hx & bit1;
				h0 = (h0 >= 0) ? (h0 - lowbit1) : (h0 + lowbit1);
				*p10 = (h0 + hx) >> 1;
				*p00 = (h0 - hx) >> 1;
			}
		}
		if (i<nxtop) {
			/*
			 * do last row if column length is odd
			 * p10, p10+1 are off edge
			 */
			pend = &a[ny*i+nytop-1];
			for (p00 = &a[ny*i]; p00 < pend; p00 += 2) {
				h0 = *(p00  );
				hy = *(p00+1);
				if (hy >= 0) {
					hy = ((hy + prnd1) & mask1);
				} else {
					hy = ((hy + nrnd1) & mask1);
				}
				lowbit1 = hy & bit1;
				h0 = (h0 >= 0) ? (h0 - lowbit1) : (h0 + lowbit1);
				*(p00+1) = (h0 + hy) >> 1;
				*(p00  ) = (h0 - hy) >> 1;
			}
			if (p00==pend) {
				/*
				 * do corner element if both row and column lengths are odd
				 * p00+1, p10, p10+1 are off edge
				 */
				*p00 = *p00 >> 1;
			}
		}
		/*
		 * divide all the masks and rounding values by 2
		 */
		bit2 = bit1;
		bit1 = bit0;
		bit0 = bit0 >> 1;
		mask1 = mask0;
		mask0 = mask0 >> 1;
		prnd1 = prnd0;
		prnd0 = prnd0 >> 1;
		nrnd1 = nrnd0;
		nrnd0 = prnd0 - 1;
	}
	/*
	 * Last pass (k=0) has some differences:
	 *
	 * Shift by 2 instead of 1
	 *
	 * Use explicit values for all variables to avoid unnecessary shifts etc:
	 *
	 *   N    bitN maskN prndN nrndN
	 *   0     1    -1     0     0  (note nrnd0 != prnd0-1)
	 *   1     2    -2     1     0
	 *   2     4    -4     2     1
	 */

	/*
	 * Check nxtop=nx, nytop=ny
	 */
	c = c>>1;
	nxtop = nxtop<<1;
	nytop = nytop<<1;
	if (nxf <= c) { nxtop -= 1; } else { nxf -= c; }
	if (nyf <= c) { nytop -= 1; } else { nyf -= c; }
	if (nxtop != nx || nytop != ny) {
		printerr("hinv: error, final image size is %dx%d, not %dx%d\n",
			nxtop, nytop, nx, ny);
		*status = FAILURE;
		return;
	}
	/*
	 * unshuffle in each dimension to interleave coefficients
	 */
	xunshuffle(a,nx,ny,ny,status);
	yunshuffle(a,nx,ny,ny,status);
	/*
	 * smooth by interpolating coefficients if SMOOTH != 0
	 */
	if (smooth) hsmooth(a,nx,ny,ny,scale);
	for (i = 0; i<nx-1; i += 2) {
		pend = &a[ny*i+ny-1];
		for (p00 = &a[ny*i], p10 = p00+ny; p00 < pend; p00 += 2, p10 += 2) {
			h0 = *(p00  );
			hx = *(p10  );
			hy = *(p00+1);
			hc = *(p10+1);
			lowbit0 = hc & 1;
			if (lowbit0==0) {
				/*
				 * lowbit0 = 0
				 */
				if (hx>=0) { hx = (hx+1) & -2; } else { hx = hx & -2; }
				if (hy>=0) { hy = (hy+1) & -2; } else { hy = hy & -2; }
				lowbit1 = (hc ^ hx ^ hy) & 2;
				h0 = (h0 >= 0) ? (h0 - lowbit1) : (h0 + lowbit1);
			} else {
				/*
				 * lowbit0 = 1
				 */
				if (hx>=0) { hx = ((hx+1) & -2)-1; } else { hx = (hx & -2)+1; }
				if (hy>=0) { hy = ((hy+1) & -2)-1; } else { hy = (hy & -2)+1; }
				lowbit1 = (hc ^ hx ^ hy) & 2;
				h0 = h0 + 1 - lowbit1;
			}
			/*
			 * Divide sums by 4
			 */
			sum1 = h0+hx;
			sum2 = hy+hc;
			*(p10+1) = (sum1 + sum2) >> 2;
			*(p10  ) = (sum1 - sum2) >> 2;
			sum1 = h0-hx;
			sum2 = hy-hc;
			*(p00+1) = (sum1 + sum2) >> 2;
			*(p00  ) = (sum1 - sum2) >> 2;
		}
		if (p00==pend) {
			/*
			 * do last element in row if row length is odd
			 * p00+1, p10+1 are off edge
			 */
			h0 = *p00;
			hx = *p10;
			if (hx >= 0) { hx = (hx + 1) & -2; } else { hx = hx & -2; }
			lowbit1 = hx & 2;
			h0 = (h0 >= 0) ? (h0 - lowbit1) : (h0 + lowbit1);
			*p10 = (h0 + hx) >> 2;
			*p00 = (h0 - hx) >> 2;
		}
	}
	if (i<nx) {
		/*
		 * do last row if column length is odd
		 * p10, p10+1 are off edge
		 */
		pend = &a[ny*i+ny-1];
		for (p00 = &a[ny*i]; p00 < pend; p00 += 2) {
			h0 = *(p00  );
			hy = *(p00+1);
			if (hy >= 0) { hy = (hy + 1) & -2; } else { hy = hy & -2; }
			lowbit1 = hy & 2;
			h0 = (h0 >= 0) ? (h0 - lowbit1) : (h0 + lowbit1);
			*(p00+1) = (h0 + hy) >> 2;
			*(p00  ) = (h0 - hy) >> 2;
		}
		if (p00==pend) {
			/*
			 * do corner element if both row and column lengths are odd
			 * p00+1, p10, p10+1 are off edge
			 */
			*p00 = *p00 >> 2;
		}
	}
}
 
static void xunshuffle(int a[], int nx, int ny, int nydim, int *status)
{
int j;
int nhalf;
int *p1, *p2, *pt, *pend, *tmp;
 
	if (*status != 0) return;
	/*
	 * get temporary storage for shuffling elements
	 */  
	tmp = (int *) malloc(((ny+1)/2)*sizeof(int));
	if (tmp == (int *) NULL) {
		printerr("hinv: insufficient memory\n");
		*status = FAILURE;
		return;
	}
	nhalf = (ny+1)>>1;
	for (j = 0; j<nx; j++) {
		/* unshuffle(&a[nydim*j],ny,1,tmp); */

		/*
		 * copy 2nd half of array to tmp
		 */
		(void) memcpy(tmp, &a[j*nydim+nhalf], (ny-nhalf)*sizeof(int));
		/*
		 * distribute 1st half of array to even elements
		 */
		pend = &a[j*nydim];
		for (p2 = &a[j*nydim+nhalf-1], p1 = &a[j*nydim+((nhalf-1)<<1)];
				p2 >= pend; p1 -= 2, p2 -= 1) {
			*p1 = *p2;
		}
		/*
		 * now distribute 2nd half of array (in tmp) to odd elements
		 */
		pend = &a[j*nydim+ny];
		for (pt = tmp, p1 = &a[j*nydim+1]; p1<pend; p1 += 2, pt += 1) {
			*p1 = *pt;
		}
	}
	free(tmp);
}

/*
 * unshuffle in y direction: take even elements from beginning of
 * array, odd elements from end and interleave them.  This is done
 * using a somewhat complicated method for efficiency.  The straightforward
 * approach is slow because the scattered memory accesses don't
 * take advantage of the cache (I think.)  This version does
 * operations a row at a time so that consecutive memory locations
 * are accessed.
 */

static void yunshuffle(int a[], int nx, int ny, int nydim, int *status)
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
		printerr("hinv: insufficient memory\n");
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
			 */
			if (j >= oddoffset) {
				/* odd row */
				k = ((j-oddoffset)<<1) + 1;
			} else {
				/* even row */
				k = j<<1;
			}
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
					if (k >= oddoffset) {
						k = ((k-oddoffset)<<1) + 1;
					} else {
						k = k<<1;
					}
				}
				/*
				 * copy the last row into place
				 * this should always end up with j=k
				 */
				(void) memcpy(&a[nydim*k],tmp,ny*sizeof(int));
				if (j != k) {
					printerr("error: yunshuffle failed!\nj=%d k=%d\n", j, k);
					*status = FAILURE;
					return;
				}
			}
		}
	}
	free(tmp);
	free(flag);
}

/* hsmooth.c	Smooth H-transform image by adjusting coefficients toward
 *				interpolated values
 *
 * Programmer: R. White		Date: 13 April 1992
 */

#define min(a,b) (((a)<(b)) ? (a) : (b))
#define max(a,b) (((a)>(b)) ? (a) : (b))

static void hsmooth(int a[], int nxtop, int nytop, int ny, int scale)
{
int i, j;
int ny2, s10, s00, diff, dmax, dmin, s, smax;
int hm, h0, hp, hmm, hpm, hmp, hpp, hx2, hy2;
int m1,m2;

	/*
	 * Maximum change in coefficients is determined by scale factor.
	 * Since we rounded during division (see digitize.c), the biggest
	 * permitted change is scale/2.
	 */
	smax = (scale >> 1);
	if (smax <= 0) return;
	ny2 = ny << 1;
	/*
	 * We're indexing a as a 2-D array with dimensions (nxtop,ny) of which
	 * only (nxtop,nytop) are used.  The coefficients on the edge of the
	 * array are not adjusted (which is why the loops below start at 2
	 * instead of 0 and end at nxtop-2 instead of nxtop.)
	 */
	/*
	 * Adjust x difference hx
	 */
	for (i = 2; i<nxtop-2; i += 2) {
		s00 = ny*i;				/* s00 is index of a[i,j]	*/
		s10 = s00+ny;			/* s10 is index of a[i+1,j]	*/
		for (j = 0; j<nytop; j += 2) {
			/*
			 * hp is h0 (mean value) in next x zone, hm is h0 in previous x zone
			 */
			hm = a[s00-ny2];
			h0 = a[s00];
			hp = a[s00+ny2];
			/*
			 * diff = 8 * hx slope that would match h0 in neighboring zones
			 */
			diff = hp-hm;
			/*
			 * monotonicity constraints on diff
			 */
			dmax = max( min( (hp-h0), (h0-hm) ), 0 ) << 2;
			dmin = min( max( (hp-h0), (h0-hm) ), 0 ) << 2;
			/*
			 * if monotonicity would set slope = 0 then don't change hx.
			 * note dmax>=0, dmin<=0.
			 */
			if (dmin < dmax) {
				diff = max( min(diff, dmax), dmin);
				/*
				 * Compute change in slope limited to range +/- smax.
				 * Careful with rounding negative numbers when using
				 * shift for divide by 8.
				 */
				s = diff-(a[s10]<<3);
				s = (s>=0) ? (s>>3) : ((s+7)>>3) ;
				s = max( min(s, smax), -smax);
				a[s10] = a[s10]+s;
			}
			s00 += 2;
			s10 += 2;
		}
	}
	/*
	 * Adjust y difference hy
	 */
	for (i = 0; i<nxtop; i += 2) {
		s00 = ny*i+2;
		s10 = s00+ny;
		for (j = 2; j<nytop-2; j += 2) {
			hm = a[s00-2];
			h0 = a[s00];
			hp = a[s00+2];
			diff = hp-hm;
			dmax = max( min( (hp-h0), (h0-hm) ), 0 ) << 2;
			dmin = min( max( (hp-h0), (h0-hm) ), 0 ) << 2;
			if (dmin < dmax) {
				diff = max( min(diff, dmax), dmin);
				s = diff-(a[s00+1]<<3);
				s = (s>=0) ? (s>>3) : ((s+7)>>3) ;
				s = max( min(s, smax), -smax);
				a[s00+1] = a[s00+1]+s;
			}
			s00 += 2;
			s10 += 2;
		}
	}
	/*
	 * Adjust curvature difference hc
	 */
	for (i = 2; i<nxtop-2; i += 2) {
		s00 = ny*i+2;
		s10 = s00+ny;
		for (j = 2; j<nytop-2; j += 2) {
			/*
			 * ------------------    y
			 * | hmp |    | hpp |    |
			 * ------------------    |
			 * |     | h0 |     |    |
			 * ------------------    -------x
			 * | hmm |    | hpm |
			 * ------------------
			 */
			hmm = a[s00-ny2-2];
			hpm = a[s00+ny2-2];
			hmp = a[s00-ny2+2];
			hpp = a[s00+ny2+2];
			h0  = a[s00];
			/*
			 * diff = 64 * hc value that would match h0 in neighboring zones
			 */
			diff = hpp + hmm - hmp - hpm;
			/*
			 * 2 times x,y slopes in this zone
			 */
			hx2 = a[s10  ]<<1;
			hy2 = a[s00+1]<<1;
			/*
			 * monotonicity constraints on diff
			 */
			m1 = min(max(hpp-h0,0)-hx2-hy2, max(h0-hpm,0)+hx2-hy2);
			m2 = min(max(h0-hmp,0)-hx2+hy2, max(hmm-h0,0)+hx2+hy2);
			dmax = min(m1,m2) << 4;
			m1 = max(min(hpp-h0,0)-hx2-hy2, min(h0-hpm,0)+hx2-hy2);
			m2 = max(min(h0-hmp,0)-hx2+hy2, min(hmm-h0,0)+hx2+hy2);
			dmin = max(m1,m2) << 4;
			/*
			 * if monotonicity would set slope = 0 then don't change hc.
			 * note dmax>=0, dmin<=0.
			 */
			if (dmin < dmax) {
				diff = max( min(diff, dmax), dmin);
				/*
				 * Compute change in slope limited to range +/- smax.
				 * Careful with rounding negative numbers when using
				 * shift for divide by 64.
				 */
				s = diff-(a[s10+1]<<6);
				s = (s>=0) ? (s>>6) : ((s+63)>>6) ;
				s = max( min(s, smax), -smax);
				a[s10+1] = a[s10+1]+s;
			}
			s00 += 2;
			s10 += 2;
		}
	}
}

/* decode.c		read codes from infile and construct array
 *
 * Programmer: R. White		Date: 16 June 1994
 */

static char code_magic[2] = { (char)0xDD, (char)0x99 };

static void decode(MYFILE *infile, int  **a, int  *nx, int *ny, int  *scale, int *status)
{
int nel, sumall;
unsigned char nbitplanes[3];
char tmagic[2];

	if (*status != 0) return;
	/*
	 * File must start with special 2-byte magic code
	 */
	qread(infile, tmagic, sizeof(tmagic), status);
	/*
	 * check for correct magic code value
	 */
	if (memcmp(tmagic,code_magic,sizeof(code_magic)) != 0) {
		printerr("bad file format\n");
		*status = FAILURE;
		return;
	}
	*nx =readint(infile,status);				/* x size of image					*/
	*ny =readint(infile,status);				/* y size of image					*/
	*scale=readint(infile,status);				/* scale factor for digitization	*/
	/*
	 * allocate memory for array
	 * use calloc so it gets initialized to zero
	 */
	nel = (*nx) * (*ny);
	*a = (int *) calloc(nel,sizeof(int));
	if (*a == (int *) NULL) {
		printerr("decode: insufficient memory\n");
		*status = FAILURE;
		return;
	}
	/* sum of all pixels	*/
	sumall=readint(infile,status);
	/* # bits in quadrants	*/
	qread(infile, (char *) nbitplanes, sizeof(nbitplanes), status);
	dodecode(infile, *a, *nx, *ny, nbitplanes, *scale, status);
	/*
	 * put sum of all pixels back into pixel 0
	 */
	if (*scale > 1) {
		(*a)[0] = sumall * (*scale);
	} else {
		(*a)[0] = sumall;
	}
}

/* dodecode.c	Decode stream of characters on infile and return array
 *
 * This version encodes the different quadrants separately
 *
 * Programmer: R. White		Date: 16 June 1994
 */


static void dodecode(MYFILE *infile, int a[], int nx, int ny,
	unsigned char nbitplanes[3], int scale, int *status)
{
int i, nel, nx2, ny2;

	if (*status != 0) return;
	nel = nx*ny;
	nx2 = (nx+1)/2;
	ny2 = (ny+1)/2;
	/*
	 * Initialize bit input
	 */
	start_inputing_bits();
	/*
	 * read bit planes for each quadrant
	 */
	qtree_decode(infile, &a[0],          ny, nx2,  ny2,  nbitplanes[0], status);
	qtree_decode(infile, &a[ny2],        ny, nx2,  ny/2, nbitplanes[1], status);
	qtree_decode(infile, &a[ny*nx2],     ny, nx/2, ny2,  nbitplanes[1], status);
	qtree_decode(infile, &a[ny*nx2+ny2], ny, nx/2, ny/2, nbitplanes[2], status);
	/*
	 * make sure there is an EOF symbol (nybble=0) at end
	 */
	if (input_nbits(infile,4,status) != 0) {
		printerr("dodecode: bad bit plane values\n");
		*status = FAILURE;
		return;
	}
	/*
	 * now get the sign bits and do scaling at same time
	 * Re-initialize bit input
	 */
	start_inputing_bits();
	if (scale>1) {
		for (i=0; i<nel; i++) {
			if (a[i] != 0) {
				if (input_bit(infile,status) != 0) {
					a[i] = -scale*a[i];
				} else {
					a[i] =  scale*a[i];
				}
			}
		}
	} else {
		/* no scaling if scale <= 1 */
		for (i=0; i<nel; i++) {
			if (a[i] != 0) {
				if (input_bit(infile,status) != 0) a[i] = -a[i];
			}
		}
	}
}

/* Bit input routines */

/* Bit buffer */

static int buffer;					/* Bits waiting to be input				*/
static int bits_to_go;				/* Number of bits still in buffer		*/

/* Initialize bit input */

static void start_inputing_bits(void)
{
	/*
	 * Buffer starts out with no bits in it
	 */
	bits_to_go = 0;
}


/* Input a bit */

static int input_bit(MYFILE *infile, int *status)
{
	if (*status != 0) return(0);
	if (bits_to_go == 0) {
		/* Read next byte if no bits left in buffer */
		buffer = mygetc(infile);
		if (buffer == EOF) {
			erreof("input_bit", status);
			return(0);
		}
		bits_to_go = 8;
	}
	/*
	 * Return the next bit
	 */
	bits_to_go -= 1;
	return((buffer>>bits_to_go) & 1);
}


/* Input N bits (N must be <= 24) */

static int input_nbits(MYFILE *infile, int n, int *status)
{
int c;
/* local copies */
int lbuffer;
int lbits_to_go;

	if (*status != 0) return(0);
	lbuffer = buffer;
	lbits_to_go = bits_to_go;
	while (lbits_to_go < n) {
		/*
		 * need another byte's worth of bits
		 */
		lbuffer <<= 8;
		c = getc(infile);
		if (c == EOF) {
			erreof("input_nbits", status);
			return(0);
		}
		lbuffer |= c;
		lbits_to_go += 8;
	}
	/*
	 * now pick off the first n bits
	 */
	lbits_to_go -= n;
	c = (lbuffer>>lbits_to_go) & ((1<<n)-1);

	bits_to_go = lbits_to_go;
	buffer = lbuffer;
	return( c );
}

/* qtree_decode.c	Read stream of codes from infile and construct bit planes
 *					in quadrant of 2-D array using binary quadtree coding.
 *					New version using address approach.
 *
 * Programmer: R. White		Date: 12 September 1994
 */

/*
 * input 1 bit without EOF error checking
 * (EOFs get picked up in dodecode after calls to qtree_decode are done)
 */
#define input_bit_noerr(infile) (bits_to_go ? \
			((buffer>>(--bits_to_go)) & 1) : \
			(((buffer = mygetc(infile))>>(bits_to_go = 7)) & 1))
/*
 * input N bits without EOF error checking
 */
#define input_nbits_noerr(infile,n) ((bits_to_go<n) ? \
			(((buffer = (buffer<<8) | getc(infile)) \
				>>  (bits_to_go += 8-n)) & ((1<<n)-1)) : \
			((buffer>>(bits_to_go -= n)) & ((1<<n)-1)) )

static void qtree_decode(MYFILE *infile, int a[], int n, int nqx, int nqy, int nbitplanes, int *status)
{
int log2n, k, bit, b, nqmax;
int nqx2, nqy2;
/*
 * address buffers
 */
ADDRTYPE *xaddr1, *xaddr2, *yaddr1, *yaddr2, *xscr, *yscr;
ADDRTYPE *xthis, *ythis, *xnext, *ynext, *ptmp;
int nt;

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
	 * allocate buffers for addresses
	 */
	nqx2=(nqx+1)/2;
	nqy2=(nqy+1)/2;
	xaddr1 = (ADDRTYPE *) malloc(nqx2*nqy2*sizeof(ADDRTYPE));
	yaddr1 = (ADDRTYPE *) malloc(nqx2*nqy2*sizeof(ADDRTYPE));
	xaddr2 = (ADDRTYPE *) malloc(((nqx2+1)/2)*((nqy2+1)/2)*sizeof(ADDRTYPE));
	yaddr2 = (ADDRTYPE *) malloc(((nqx2+1)/2)*((nqy2+1)/2)*sizeof(ADDRTYPE));
	xscr   = (ADDRTYPE *) malloc(nqmax*sizeof(ADDRTYPE));
	yscr   = (ADDRTYPE *) malloc(nqmax*sizeof(ADDRTYPE));
	if ((xaddr1 == (ADDRTYPE *) NULL) || (yaddr1 == (ADDRTYPE *) NULL) ||
		(xaddr2 == (ADDRTYPE *) NULL) || (yaddr2 == (ADDRTYPE *) NULL) ||
		(xscr   == (ADDRTYPE *) NULL) || (yscr   == (ADDRTYPE *) NULL)) {
		printerr("qtree_decode: insufficient memory\n");
		*status = FAILURE;
		return;
	}
	/*
	 * now decode each bit plane, starting at the top
	 * A is assumed to be initialized to zero
	 */
	for (bit = nbitplanes-1; bit >= 0; bit--) {
		/*
		 * Was bitplane was quadtree-coded or written directly?
		 */
		b = input_nbits_noerr(infile,4);
		if(b == 0) {
			/*
			 * bit map was written directly
			 */
			read_bdirect(infile,a,n,nqx,nqy,bit);
		} else if (b != 0xf) {
			printerr("qtree_decode: bad format code %x\n",b);
			*status = FAILURE;
			return;
		} else {
			/*
			 * bitmap was quadtree-coded, do log2n expansions
			 * read first code
			 */
			b = input_huffman(infile);
			/*
			 * if code is zero, implying all bits are zero, just
			 * skip the rest for this plane
			 */
			if (b != 0) {
				/*
				 * initialize pointers to buffer pairs
				 * want to end up writing to buffer 1
				 */
				if ((log2n & 1) == 0) {
					xthis = xaddr1;
					ythis = yaddr1;
					xnext = xaddr2;
					ynext = yaddr2;
				} else {
					xthis = xaddr2;
					ythis = yaddr2;
					xnext = xaddr1;
					ynext = yaddr1;
				}
				nt = 0;
				if ((b & 1) != 0) { xthis[nt] = 1; ythis[nt++] = 1; }
				if ((b & 2) != 0) { xthis[nt] = 0; ythis[nt++] = 1; }
				if ((b & 4) != 0) { xthis[nt] = 1; ythis[nt++] = 0; }
				if ((b & 8) != 0) { xthis[nt] = 0; ythis[nt++] = 0; }
				/*
				 * now do log2n expansions, reading codes from file as necessary
				 */
				for (k = 1; k<log2n-1; k++) {
					nt = qaddr_expand(infile, xthis, ythis, nt,
						xnext, ynext, xscr, yscr);
					/* swap buffers */
					ptmp  = xthis; xthis = xnext; xnext = ptmp;
					ptmp  = ythis; ythis = ynext; ynext = ptmp;
				}
				/*
				 * now copy last set of 4-bit codes to bitplane bit of array a
				 */
				qaddr_bitins(infile, xthis, ythis, nt, a, n, bit);
			}
		}
	}
	free(xaddr1);
	free(yaddr1);
	free(xaddr2);
	free(yaddr2);
	free(xscr);
	free(yscr);
}


/*
 * do one quadtree address expansion step
 * current address list in xthis, ythis is expanded and put into xnext, ynext
 * return value is new number of elements
 */

static int qaddr_expand(MYFILE *infile, ADDRTYPE *xthis, ADDRTYPE *ythis, int nt, ADDRTYPE *xnext, ADDRTYPE *ynext, ADDRTYPE *xscr, ADDRTYPE *yscr)
{
int b, i, j, k, m;
ADDRTYPE ylast, xoff, yoff;

	/*
	 * read 1 quad for each element of xthis,ythis
	 * keep second row expansions in xscr,yscr and copy them to next buffer
	 * after each row is finished.
	 */
	ylast = ythis[0];
	k = 0;
	for (i = 0, j = 0; i < nt; i++) {
		if (ylast != ythis[i]) {
			for (m = 0; m < k; m++, j++) {
				xnext[j] = xscr[m];
				ynext[j] = yscr[m];
			}
			k = 0;
			ylast = ythis[i];
		}
		b = input_huffman(infile);
		xoff = xthis[i] << 1;
		yoff = ythis[i] << 1;
		if ((b & 1) != 0) { xnext[j] = xoff | 1; ynext[j++] = yoff | 1; }
		if ((b & 2) != 0) { xnext[j] = xoff    ; ynext[j++] = yoff | 1; }
		if ((b & 4) != 0) { xscr[k]  = xoff | 1; yscr[k++]  = yoff    ; }
		if ((b & 8) != 0) { xscr[k]  = xoff    ; yscr[k++]  = yoff    ; }
	}
	for (m = 0; m < k; m++, j++) {
		xnext[j] = xscr[m];
		ynext[j] = yscr[m];
	}
	return (j);
}

/*
 * Read quads and insert in address locations in bitplane BIT of array A.
 */

static void qaddr_bitins(MYFILE *infile, ADDRTYPE *xthis, ADDRTYPE *ythis, int nt, int a[], int n, int bit)
{
int b, i, *p, bitval;

	bitval = 1 << bit;
	for (i = 0; i < nt; i++) {
		b = input_huffman(infile);

		p = &a[ (n * ythis[i] + xthis[i]) << 1 ];

		if ((b & 8) != 0) *p       |= bitval;
		if ((b & 4) != 0) *(p+1)   |= bitval;
		if ((b & 2) != 0) *(p+n)   |= bitval;
		if ((b & 1) != 0) *(p+n+1) |= bitval;
	}
}

/*
 * Read 4-bit codes from infile, expanding each value to 2x2 pixels and
 * inserting into bitplane BIT of B.
 */

static void read_bdirect(MYFILE *infile, int b[], int n, int nx, int ny, int bit)
{
int i, tmp, bitval;
int *p00, *pend;

	bitval = 1 << bit;
	for (i = 0; i<nx; i += 2) {
		p00 = &b[n*i];
		pend = p00 + ny;
		for ( ; p00 < pend; p00 += 2) {
			tmp = input_nbits_noerr(infile,4);
			if ((tmp & 8) != 0) *(p00    ) |= bitval;
			if ((tmp & 4) != 0) *(p00+1  ) |= bitval;
			if ((tmp & 2) != 0) *(p00+n  ) |= bitval;
			if ((tmp & 1) != 0) *(p00+n+1) |= bitval;
		}
	}
}

/*
 * Huffman decoding for fixed codes
 *
 * Coded values range from 0-15
 *
 * Huffman code values (hex):
 *
 *	3e, 00, 01, 08, 02, 09, 1a, 1b,
 *	03, 1c, 0a, 1d, 0b, 1e, 3f, 0c
 *
 * and number of bits in each code:
 *
 *	6,  3,  3,  4,  3,  4,  5,  5,
 *	3,  5,  4,  5,  4,  5,  6,  4
 */

/*
 * table of Huffman code translated values
 * -1 means no such code
 */
static int tabhuff[31] =
	/* 00  01  02  03  04  05  06  07  08  09  0a  0b  0c  0d  0e  0f */
	 {  1,  2,  4,  8, -1, -1, -1, -1,  3,  5, 10, 12, 15, -1, -1, -1,
	   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  6,  7,  9, 11, 13     };
	/* 10  11  12  13  14  15  16  17  18  19  1a  1b  1c  1d  1e     */

static int input_huffman(MYFILE *infile)
{
int c;

	/*
	 * get first 3 bits to start
	 */
	c = input_nbits_noerr(infile,3);
	if (c < 4) return(tabhuff[c]);
	/*
	 * get the next bit
	 */
	c = input_bit_noerr(infile) | (c<<1);
	if (c < 13) return(tabhuff[c]);
	/*
	 * get yet another bit
	 */
	c = input_bit_noerr(infile) | (c<<1);
	if (c < 31) return(tabhuff[c]);
	/*
	 * the 6th bit decides
	 */
	if (input_bit_noerr(infile)) {
		/*
		 * c = 63
		 */
		return(14);
	} else {
		/*
		 * c = 62
		 */
		return(0);
	}
}

/* qread.c	Read binary data
 *
 * Programmer: R. White		Date: 11 March 1991
 */

static int readint(MYFILE *infile, int *status)
{
int a,i;
unsigned char b[4];

	if (*status != 0) return(0);
    /* Read integer A one byte at a time from infile.
     *
     * This is portable from Vax to Sun since it eliminates the
     * need for byte-swapping.
     */
    for (i=0; i<4; i++) qread(infile,(char *) &b[i],1,status);
	if (*status != 0) return(0);
    a = b[0];
    for (i=1; i<4; i++) a = (a<<8) + b[i];
    return(a);
}

static void qread(MYFILE *infile, char *a, int n, int *status)
{
	if (*status != 0) return;
    if(myread(infile, a, n) != n) {
        printerr("qread: %s\n", strerror(errno));
		*status = FAILURE;
    }
}

/*
 * read n bytes from file into buffer
 * returns number of bytes read (=n) if successful, <=0 if not
 */

static int myread(MYFILE *file, char buffer[], int n)
{
#ifdef TO_A_BUFFER
    /*
     * this version used when doing I/O to a buffer
     */
    if (file->current + n > file->end) return(-1);
    (void) memcpy(buffer, file->current, n);
    file->current += n;
    return(n);
#else
    /*
     * this version used when doing I/O to a file
     *
     * this version is for VMS C: each read may return
     * fewer than n bytes, so multiple reads may be needed
     * to fill the buffer.
     *
     * I think this is unnecessary for Sun Unix, but it won't hurt
     * either, so I'll leave it in.
     */
    int nread, total;
    /* keep reading until we've read n bytes */
    total = 0;
    while ( (nread = fread(&buffer[total], 1, n-total, file)) > 0) {
        total += nread;
        if (total==n) return(total);
    }
    /* end-of-file or error occurred if we got to here */
    return(nread);
#endif
}

/* Print error messages */

static void erreof(const char *str, int *status)
{
	/*
	 * end of file is an error for this application
	 */
	printerr("%s: unexpected end-of-file\n", str);
	*status = FAILURE;
}

static void printerr(const char *format, ...)
{
va_list args;

	va_start (args, format);
	vfprintf(stderr, format, args);
	va_end(args);
}
