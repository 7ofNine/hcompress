/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* decode.c		read codes from infile and construct array
 *
 * Programmer: R. White		Date: 16 June 1994
 */
#include "hcompress.h"

static char code_magic[2] = { (char)0xDD, (char)0x99 };

#ifdef __STDC__
extern void decode(MYFILE *infile, FILE *outfile, int  **a,
	int  *nx, int *ny, int  *scale, char **format)
#else
extern void decode(infile,outfile,a,nx,ny,scale,format)
MYFILE *infile;				/* input file							*/
FILE *outfile;				/* output file (NULL for none)			*/
int  **a;					/* address of output array [nx][ny]		*/
int  *nx,*ny;				/* size of output array					*/
int  *scale;				/* scale factor for digitization		*/
char **format;				/* output file format					*/
							/* (may be changed if empty)			*/
#endif

{
int nel, sumall, newfits = 0;
unsigned char nbitplanes[3];
char tmagic[2], line[81];

	/*
	 * File starts either with special 2-byte magic code or with
	 * FITS keyword "SIMPLE  ="
	 */
	qread(infile, tmagic, sizeof(tmagic));
	/*
	 * Check for FITS
	 */
	if (strncmp(tmagic,"SI",2)==0) {
		/*
		 * read rest of line and make sure the whole keyword is correct
		 */
		(void) strncpy(line,"SI",2);
		if ((myread(infile, &line[2], 78) != 78) || 
			(strncmp(line, "SIMPLE  =", 9) != 0) ) {
			fprintf(stderr, "bad file format\n");
			exit(-1);
		}
		/*
		 * set output format to default "fits" if it is empty
		 */
		if ((*format)[0] == '\0') *format = "fits";
		/*
		 * if fits output format and outfile != NULL, write this line to
		 * outfile and then copy rest of FITS header; else just read past
		 * FITS header.
		 */
		if (strcmp(*format, "fits") == 0 && outfile != (FILE *) NULL) {
            if (fwrite(line,1,80,outfile) != 80) {
                perror("stdout");
                fprintf(stderr, "decode: error writing output file\n");
                exit(-1);
            }
			fitspass(infile,1,outfile);
		} else {
			fitspass(infile,0,outfile);
		}
		/*
		 * now read the first two characters again -- this time they
		 * MUST be the magic code!
		 */
		qread(infile, tmagic, sizeof(tmagic));
	} else {
		/*
		 * set default format to raw if it is not specified
		 */
		if((*format)[0] == '\0') *format = "raw";
		/*
		 * if input format is not FITS but output format is FITS, set
		 * a flag so we generate a FITS header once we know how big
		 * the image must be.
		 */
		if (strcmp(*format, "fits") == 0) newfits = 1;
	}
	/*
	 * check for correct magic code value
	 */
	if (memcmp(tmagic,code_magic,sizeof(code_magic)) != 0) {
		fprintf(stderr, "bad file format\n");
		exit(-1);
	}
	*nx =readint(infile);				/* x size of image					*/
	*ny =readint(infile);				/* y size of image					*/
	*scale=readint(infile);				/* scale factor for digitization	*/
	/*
	 * write the new fits header to outfile if needed
	 */
	if (newfits && (outfile != (FILE *) NULL))
		makefits(outfile,*nx,*ny,16,"INTEGER*2");
	/*
	 * allocate memory for array
	 * use calloc so it gets initialized to zero
	 */
	nel = (*nx) * (*ny);
	*a = (int *) calloc(nel,sizeof(int));
	if (*a == (int *) NULL) {
		fprintf(stderr, "decode: insufficient memory\n");
		exit(-1);
	}
	/* sum of all pixels	*/
	sumall=readint(infile);
	/* # bits in quadrants	*/
	qread(infile, (char *) nbitplanes, sizeof(nbitplanes));
	dodecode(infile, *a, *nx, *ny, nbitplanes, *scale);
	/*
	 * put sum of all pixels back into pixel 0
	 */
	if (*scale > 1) {
		(*a)[0] = sumall * (*scale);
	} else {
		(*a)[0] = sumall;
	}
}
