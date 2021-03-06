/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* hinv.c   Inverse H-transform of NX x NY integer image
 *
 * Programmer: R. White		Date: 16 June 1994
 */
#include "hcompress.h"

#ifdef __STDC__
static void xunshuffle(int a[], int nx, int ny, int nydim);
static void yunshuffle(int a[], int nx, int ny, int nydim);
#else
static void xunshuffle();
static void yunshuffle();
#endif

#ifdef __STDC__
extern void hinv(int a[], int nx, int ny, int smooth, int scale)
#else
extern void hinv(a,nx,ny,smooth,scale)
int a[];
int nx,ny;
int smooth;   /* 0 for no smoothing, else smooth during inversion */
int scale;    /* used if smoothing is specified */
#endif
{
int nmax, log2n, i, k;
int nxtop,nytop,nxf,nyf,c;
int bit0, bit1, bit2, mask0, mask1, mask2,
	prnd0, prnd1, prnd2, nrnd0, nrnd1, nrnd2, lowbit0, lowbit1;
int h0, hx, hy, hc, sum1, sum2;
int *p00, *p10, *pend;

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
		xunshuffle(a,nxtop,nytop,ny);
		yunshuffle(a,nxtop,nytop,ny);
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
		fprintf(stderr,
			"hinv: error, final image size is %dx%d, not %dx%d\n",
			nxtop, nytop, nx, ny);
		exit(-1);
	}
	/*
	 * unshuffle in each dimension to interleave coefficients
	 */
	xunshuffle(a,nx,ny,ny);
	yunshuffle(a,nx,ny,ny);
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
 
#ifdef __STDC__
static void xunshuffle(int a[], int nx, int ny, int nydim)
#else
static void xunshuffle(a,nx,ny,nydim)
int a[];				/* array to unshuffle				*/
int nx;					/* number of elements in column		*/
int ny;					/* number of elements in row		*/
int nydim;				/* physical length of row in array	*/
#endif
{
int j;
int nhalf;
int *p1, *p2, *pt, *pend, *tmp;
 
	/*
	 * get temporary storage for shuffling elements
	 */  
	tmp = (int *) malloc(((ny+1)/2)*sizeof(int));
	if (tmp == (int *) NULL) {
		fprintf(stderr, "hinv: insufficient memory\n");
		exit(-1);
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

#ifdef __STDC__
static void yunshuffle(int a[], int nx, int ny, int nydim)
#else
static void yunshuffle(a,nx,ny,nydim)
int a[];				/* array to unshuffle				*/
int nx;					/* number of elements in column		*/
int ny;					/* number of elements in row		*/
int nydim;				/* actual length of row in array	*/
#endif
{
int j, k, tt, oddoffset, *tmp;
int *p, *pt;
unsigned char *flag;

	/*
	 * get temporary storage for shuffling elements
	 */
	tmp  =           (int *) malloc(ny*sizeof(int));
	flag = (unsigned char *) malloc(nx*sizeof(unsigned char));
	if(tmp == (int *) NULL || flag == (unsigned char *) NULL) {
		fprintf(stderr, "hinv: insufficient memory\n");
		exit(-1);
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
					fprintf(stderr, "error: yunshuffle failed!\nj=%d k=%d\n",
						j, k);
					exit(-1);
				}
			}
		}
	}
	free(tmp);
	free(flag);
}
