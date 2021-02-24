/* Copyright (c) 1993 Association of Universities for Research
 * in Astronomy. All rights reserved. Produced under National
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */

/* BIT INPUT ROUTINES */

#include "hcompress.h"

/* THE BIT BUFFER */

int buffer;					/* Bits waiting to be input				*/
int bits_to_go;				/* Number of bits still in buffer		*/

#ifdef __STDC__
static void erreof();
#else
static void erreof();
#endif

/* INITIALIZE BIT INPUT */

#ifdef __STDC__
extern void start_inputing_bits(void)
#else
extern void start_inputing_bits()
#endif
{
	/*
	 * Buffer starts out with no bits in it
	 */
	bits_to_go = 0;
}


/* INPUT A BIT */

#ifdef __STDC__
extern int input_bit(MYFILE *infile)
#else
extern int input_bit(infile)
MYFILE *infile;
#endif
{
	if (bits_to_go == 0) {
		/* Read next byte if no bits left in buffer */
		buffer = mygetc(infile);
		if (buffer == EOF) erreof("input_bit");
		bits_to_go = 8;
	}
	/*
	 * Return the next bit
	 */
	bits_to_go -= 1;
	return((buffer>>bits_to_go) & 1);
}


/* INPUT N BITS (N must be <= 24) */

#ifdef __STDC__
extern int input_nbits(MYFILE *infile, int n)
#else
extern int input_nbits(infile,n)
MYFILE *infile;
int n;
#endif
{
int c;
/* local copies */
int lbuffer;
int lbits_to_go;

	lbuffer = buffer;
	lbits_to_go = bits_to_go;
	while (lbits_to_go < n) {
		/*
		 * need another byte's worth of bits
		 */
		lbuffer <<= 8;
		c = getc(infile);
		if (c == EOF) erreof("input_nbits");
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

#ifdef __STDC__
static void erreof(char *str)
#else
static void erreof(str)
char *str;
#endif
{
	/*
	 * end of file is an error for this application
	 */
	fprintf(stderr, "%s: unexpected end-of-file\n", str);
	exit(-1);
}
