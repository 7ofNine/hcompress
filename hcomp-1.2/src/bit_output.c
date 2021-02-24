/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* BIT OUTPUT ROUTINES */

#include "hcompress.h"

int bitcount;

/* THE BIT BUFFER */

static int buffer;					/* Bits buffered for output			*/
static int bits_to_go;				/* Number of bits free in buffer	*/


/* INITIALIZE FOR BIT OUTPUT */

#ifdef __STDC__
extern void start_outputing_bits(void)
#else
extern void start_outputing_bits()
#endif
{
	buffer = 0;						/* Buffer is empty to start	*/
	bits_to_go = 8;					/* with						*/
	bitcount = 0;
}


/* OUTPUT A BIT */

#ifdef __STDC__
extern void output_bit(MYFILE *outfile, int bit)
#else
extern void output_bit(outfile,bit)
MYFILE *outfile;
int bit;
#endif
{
	buffer <<= 1;								/* Put bit at end of buffer */
	if (bit) buffer |= 1;
	bits_to_go -= 1;
	if (bits_to_go == 0) {						/* Output buffer if it is	*/
		(void) myputc(buffer & 0xff,outfile);	/* now full					*/
		bits_to_go = 8;
		buffer = 0;
		bitcount += 8;
	}
}


/* OUTPUT N BITS (N must be <= 24) */

#ifdef __STDC__
extern void output_nbits(MYFILE *outfile, int bits, int n)
#else
extern void output_nbits(outfile,bits,n)
MYFILE *outfile;
int bits;
int n;
#endif
{
/* local copies */
int lbuffer;
int lbits_to_go;

	/*
	 * insert bits at end of buffer
	 */
	lbuffer = buffer;
	lbits_to_go = bits_to_go;
	lbuffer <<= n;
	lbuffer |= ( bits & ((1<<n)-1) );
	lbits_to_go -= n;
	while (lbits_to_go <= 0) {
		/*
		 * buffer full, put out top 8 bits
		 */
		(void) myputc((lbuffer>>(-lbits_to_go)) & 0xff,outfile);
		lbits_to_go += 8;
		bitcount += 8;
	}
	buffer = lbuffer;
	bits_to_go = lbits_to_go;
}


/* FLUSH OUT THE LAST BITS */

#ifdef __STDC__
extern void done_outputing_bits(MYFILE *outfile)
#else
extern void done_outputing_bits(outfile)
MYFILE *outfile;
#endif
{
	if(bits_to_go < 8) {
		(void) myputc(buffer<<bits_to_go,outfile);
		/* count the garbage bits too */
		bitcount += 8;
	}
}
