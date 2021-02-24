/* Copyright (c) 1994 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* qtree_decode.c	Read stream of codes from infile and construct bit planes
 *					in quadrant of 2-D array using binary quadtree coding.
 *					New version using address approach.
 *
 * Programmer: R. White		Date: 12 September 1994
 */
#include "hcompress.h"

/* bit buffer from bit_input.c */

int buffer;					/* Bits waiting to be input				*/
int bits_to_go;				/* Number of bits still in buffer		*/

/*
 * input 1 bit without EOF error checking
 * (EOFs get picked up in dodecode after calls to qtree_decode are done)
 */
#define input_bit(infile) (bits_to_go ? \
			((buffer>>(--bits_to_go)) & 1) : \
			(((buffer = mygetc(infile))>>(bits_to_go = 7)) & 1))
/*
 * input N bits without EOF error checking
 */
#define input_nbits(infile,n) ((bits_to_go<n) ? \
			(((buffer = (buffer<<8) | getc(infile)) \
				>>  (bits_to_go += 8-n)) & ((1<<n)-1)) : \
			((buffer>>(bits_to_go -= n)) & ((1<<n)-1)) )
/*
 * input 4 bits
 */
#define input_nybble(infile)    input_nbits(infile,4)

/*
 * address buffer type (could be unsigned char if image size < 512 pixels)
 */
typedef short ADDRTYPE;

#ifdef __STDC__

static int qaddr_expand(MYFILE *infile, ADDRTYPE *xthis, ADDRTYPE *ythis,
	int nt, ADDRTYPE *xnext, ADDRTYPE *ynext, ADDRTYPE *xscr, ADDRTYPE *yscr);
static void qaddr_bitins(MYFILE *infile, ADDRTYPE *xthis, ADDRTYPE *ythis,
	int nt, int a[], int n, int bit);
static void read_bdirect(MYFILE *infile, int b[], int n, int nx,
	int ny, int bit);
static int input_huffman(MYFILE *infile);

#else

static int qaddr_expand();
static void qaddr_bitins();
static void read_bdirect();
static int input_huffman();

#endif

#ifdef __STDC__
extern void qtree_decode(MYFILE *infile, int a[], int n, int nqx, int nqy, int nbitplanes)
#else
extern void qtree_decode(infile,a,n,nqx,nqy,nbitplanes)
MYFILE *infile;
int a[];							/* a is 2-D array with dimensions (n,n)	*/
int n;								/* length of full row in a				*/
int nqx;							/* partial length of row to decode		*/
int nqy;							/* partial length of column (<=n)		*/
int nbitplanes;						/* number of bitplanes to decode		*/
#endif
{
int log2n, k, bit, b, nqmax;
int nqx2, nqy2;
/*
 * address buffers
 */
ADDRTYPE *xaddr1, *xaddr2, *yaddr1, *yaddr2, *xscr, *yscr;
ADDRTYPE *xthis, *ythis, *xnext, *ynext, *ptmp;
int nt;

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
		fprintf(stderr, "qtree_decode: insufficient memory\n");
		exit(-1);
	}
	/*
	 * now decode each bit plane, starting at the top
	 * A is assumed to be initialized to zero
	 */
	for (bit = nbitplanes-1; bit >= 0; bit--) {
		/*
		 * Was bitplane was quadtree-coded or written directly?
		 */
		b = input_nybble(infile);
		if(b == 0) {
			/*
			 * bit map was written directly
			 */
			read_bdirect(infile,a,n,nqx,nqy,bit);
		} else if (b != 0xf) {
			fprintf(stderr, "qtree_decode: bad format code %x\n",b);
			exit(-1);
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

#ifdef __STDC__
static int qaddr_expand(MYFILE *infile, ADDRTYPE *xthis, ADDRTYPE *ythis, int nt, ADDRTYPE *xnext, ADDRTYPE *ynext, ADDRTYPE *xscr, ADDRTYPE *yscr)
#else
static int qaddr_expand(infile, xthis, ythis, nt, xnext, ynext, xscr, yscr)
MYFILE *infile;
ADDRTYPE *xthis;
ADDRTYPE *ythis;
int    nt;
ADDRTYPE *xnext;
ADDRTYPE *ynext;
ADDRTYPE *xscr; /* needs to have nx elements */
ADDRTYPE *yscr;
#endif
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

#ifdef __STDC__
static void qaddr_bitins(MYFILE *infile, ADDRTYPE *xthis, ADDRTYPE *ythis, int nt, int a[], int n, int bit)
#else
static void qaddr_bitins(infile, xthis, ythis, nt, a, n, bit)
MYFILE *infile;
ADDRTYPE *xthis;
ADDRTYPE *ythis;
int nt;
int a[];
int n;		/* declared y dimension of a */
int bit;
#endif
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

#ifdef __STDC__
static void read_bdirect(MYFILE *infile, int b[], int n, int nx, int ny, int bit)
#else
static void read_bdirect(infile,b,n,nx,ny,bit)
MYFILE *infile;
int b[];
int n;		/* declared y dimension of b */
int nx;
int ny;
int bit;
#endif
{
int i, tmp, bitval;
int *p00, *pend;

	bitval = 1 << bit;
	for (i = 0; i<nx; i += 2) {
		p00 = &b[n*i];
		pend = p00 + ny;
		for ( ; p00 < pend; p00 += 2) {
			tmp = input_nybble(infile);
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

#ifdef __STDC__
static int input_huffman(MYFILE *infile)
#else
static int input_huffman(infile)
MYFILE *infile;
#endif
{
int c;

	/*
	 * get first 3 bits to start
	 */
	c = input_nbits(infile,3);
	if (c < 4) return(tabhuff[c]);
	/*
	 * get the next bit
	 */
	c = input_bit(infile) | (c<<1);
	if (c < 13) return(tabhuff[c]);
	/*
	 * get yet another bit
	 */
	c = input_bit(infile) | (c<<1);
	if (c < 31) return(tabhuff[c]);
	/*
	 * the 6th bit decides
	 */
	if (input_bit(infile)) {
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
