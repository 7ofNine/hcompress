/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* dodecode.c	Decode stream of characters on infile and return array
 *
 * This version encodes the different quadrants separately
 *
 * Programmer: R. White		Date: 16 June 1994
 */

#include "hcompress.h"

#define input_nybble(infile)	input_nbits(infile,4)

#ifdef __STDC__
extern void dodecode(MYFILE *infile, int a[], int nx, int ny,
	unsigned char nbitplanes[3], int scale)
#else
extern void dodecode(infile,a,nx,ny,nbitplanes,scale)
MYFILE *infile;
int a[];							/* Array of values to decode (zero)		*/
int nx,ny;							/* Array dimensions are [nx][ny]		*/
unsigned char nbitplanes[3];		/* Number of bit planes in quadrants	*/
int scale;							/* Scale factor used in digitize		*/
#endif
{
int i, nel, nx2, ny2;

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
	qtree_decode(infile, &a[0],          ny, nx2,  ny2,  nbitplanes[0]);
	qtree_decode(infile, &a[ny2],        ny, nx2,  ny/2, nbitplanes[1]);
	qtree_decode(infile, &a[ny*nx2],     ny, nx/2, ny2,  nbitplanes[1]);
	qtree_decode(infile, &a[ny*nx2+ny2], ny, nx/2, ny/2, nbitplanes[2]);
	/*
	 * make sure there is an EOF symbol (nybble=0) at end
	 */
	if (input_nybble(infile) != 0) {
		fprintf(stderr, "dodecode: bad bit plane values\n");
		exit(-1);
	}
	/*
	 * now get the sign bits and do scaling at same time
	 * Re-initialize bit input
	 */
	start_inputing_bits();
	if (scale>1) {
		for (i=0; i<nel; i++) {
			if (a[i] != 0) {
				if (input_bit(infile) != 0) {
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
				if (input_bit(infile) != 0) a[i] = -a[i];
			}
		}
	}
}
