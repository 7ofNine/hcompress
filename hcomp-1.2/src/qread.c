/* Copyright (c) 1993 Association of Universities for Research 
 * in Astronomy. All rights reserved. Produced under National   
 * Aeronautics and Space Administration Contract No. NAS5-26555.
 */
/* qread.c	Read binary data
 *
 * Programmer: R. White		Date: 11 March 1991
 */
#include "hcompress.h"

#ifdef __STDC__
extern int readint(MYFILE *infile)
#else
extern int readint(infile)
MYFILE *infile;
#endif
{
int a,i;
unsigned char b[4];

    /* Read integer A one byte at a time from infile.
     *
     * This is portable from Vax to Sun since it eliminates the
     * need for byte-swapping.
     */
    for (i=0; i<4; i++) qread(infile,(char *) &b[i],1);
    a = b[0];
    for (i=1; i<4; i++) a = (a<<8) + b[i];
    return(a);
}

#ifdef __STDC__
extern void qread(MYFILE *infile, char *a, int n)
#else
extern void qread(infile,a,n)
MYFILE *infile;
char *a;
int n;
#endif
{
    if(myread(infile, a, n) != n) {
        perror("qread");
        exit(-1);
    }
}

/*
 * read n bytes from file into buffer
 * returns number of bytes read (=n) if successful, <=0 if not
 */

#ifdef __STDC__
extern int myread(MYFILE *file, char buffer[], int n)
#else
extern int myread(file, buffer, n)
MYFILE *file;
char buffer[];
int n;
#endif
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
