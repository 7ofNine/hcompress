/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* digitize.c	digitize H-transform
 *
 * Programmer: R. White		Date: 15 June 1994
 */
#include "hcompress.h"

#ifdef __STDC__
extern void digitize(int a[], int nx, int ny, int scale)
#else
extern void digitize(a,nx,ny,scale)
int a[];
int nx,ny;
int scale;
#endif
{
int d, *p;

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
