/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* qwrite.c Write binary data
 *
 * Programmer: R. White     Date: 11 March 1991
 */
#include "hcompress.h"

#ifdef __STDC__
extern void writeint(MYFILE *outfile, int a)
#else
extern void writeint(outfile,a)
MYFILE *outfile;
int a;
#endif
{
int i;
unsigned char b[4];

    /* Write integer A one byte at a time to outfile.
     *
     * This is portable from Vax to Sun since it eliminates the
     * need for byte-swapping.
     */
    for (i=3; i>=0; i--) {
        b[i] = a & 0x000000ff;
        a >>= 8;
    }
    for (i=0; i<4; i++) qwrite(outfile,(char *) &b[i],1);
}

#ifdef __STDC__
extern void qwrite(MYFILE *outfile, char *a, int n)
#else
extern void qwrite(outfile,a,n)
MYFILE *outfile;
char *a;
int n;
#endif
{
    if(mywrite(outfile, a, n) != n) {
        perror("qwrite");
        exit(-1);
    }
}

/*
 * write n bytes from buffer into file
 * returns number of bytes written (=n) if successful, <=0 if not
 */

#ifdef __STDC__
extern int mywrite(MYFILE *file, char buffer[], int n)
#else
extern int mywrite(file, buffer, n)
MYFILE *file;
char buffer[];
int n;
#endif
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
		fprintf(stderr, "mywrite: buffer overflow\n");
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
