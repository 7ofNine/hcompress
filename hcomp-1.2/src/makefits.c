/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* makefits.c	Write a simple FITS header for a 2-D image to outfile.
 *
 * Programmer: R. White		Date: 24 April 1992
 *
 * example of header
 * lines must be exactly 80 characters, with no newline at the end
 *
0123456789 123456789 123456789 123456789
SIMPLE  =                    T /Standard FITS format
BITPIX  =                   16 /
NAXIS   =                    2 /Number of axes
NAXIS1  =                  256 /
NAXIS2  =                  256 /
DATATYPE= 'INTEGER*2'          /
END
 *
 */
#include "hcompress.h"

#ifdef __STDC__
extern void makefits(FILE *outfile, int nx, int ny, int bitpix, char *datatype)
#else
extern void makefits(outfile,nx,ny,bitpix,datatype)
FILE *outfile;
int nx,ny;
int bitpix;
char *datatype;
#endif
{
char line[81];
int i;

	fprintf(outfile, "%-80.80s",
		"SIMPLE  =                    T /Standard FITS format");

	(void) strcpy(line,"BITPIX  =                      /");
	(void) sprintf(&line[10], "%20d", bitpix);
	line[30] = ' ';
	fprintf(outfile, "%-80.80s", line);

	fprintf(outfile, "%-80.80s",
		"NAXIS   =                    2 /Number of axes");

	(void) strcpy(line,"NAXIS1  =                      /");
	(void) sprintf(&line[10], "%20d", ny);
	line[30] = ' ';
	fprintf(outfile, "%-80.80s", line);

	(void) strcpy(line,"NAXIS2  =                      /");
	(void) sprintf(&line[10], "%20d", nx);
	line[30] = ' ';
	fprintf(outfile, "%-80.80s", line);

	(void) strcpy(line,"DATATYPE=                      /");
	(void) sprintf(&line[10], "'%*.*s'", (int) strlen(datatype),
		(int) strlen(datatype), datatype);
	line[12+strlen(datatype)] = ' ';
	fprintf(outfile, "%-80.80s", line);

	fprintf(outfile, "%-80.80s", "END");

	/*
	 * pad with blank lines to get multiple of 36 (2880 bytes)
	 */
	for (i=0; i<80; i++) line[i] = ' ';
	line[80] = '\0';
	for (i=7; i<36; i++) (void) fputs(line,outfile);
}
