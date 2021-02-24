#ifndef _HCOMPRESS_H
#define _HCOMPRESS_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "config.h"

/* MYFILE structure for compression to a file
 * Defined this way so same code can be used for file or buffer.
 */
typedef FILE MYFILE;
#define mygetc getc
#define myputc putc

/* function prototypes */

#ifdef __STDC__

/* ANSI C versions */

/* bit I/O */
extern void start_inputing_bits(void);
extern int input_bit(MYFILE *infile);
extern int input_nbits(MYFILE *infile, int n);
extern void start_outputing_bits(void);
extern void output_bit(MYFILE *outfile, int bit);
extern void output_nbits(MYFILE *outfile, int bits, int n);
extern void done_outputing_bits(MYFILE *outfile);

/* hcompress routines */
extern void decode(MYFILE *infile, FILE *outfile,
	int  **a, int  *nx, int *ny, int  *scale, char **format);
extern void digitize(int a[], int nx, int ny, int scale);
extern void dodecode(MYFILE *infile, int a[], int nx, int ny,
	unsigned char nbitplanes[3], int scale);
extern void doencode(MYFILE *outfile, int a[], int nx, int ny,
	unsigned char nbitplanes[3]);
extern void encode(MYFILE *outfile, int a[], int nx, int ny,
	int scale);
extern void hinv(int a[], int nx, int ny, int smooth, int scale);
extern void hsmooth(int a[], int nxtop, int nytop, int ny, int scale);
extern void htrans(int a[], int nx, int ny);
extern void qtree_decode(MYFILE *infile, int a[], int n, int nqx, int nqy,
	int nbitplanes);
extern void qtree_encode(MYFILE *outfile, int a[], int n, int nqx, int nqy,
	int nbitplanes);

/* data I/O routines */
extern void fitspass(MYFILE *infile, int  passthru, FILE *outfile);
extern void fitsread(FILE *infile, char *inname, FILE *outfile,
	int  *nx, int *ny, int  passthru, int  padded, int  nlterm);
extern void get_data(FILE *infile[2], char *inname[2], FILE *outfile,
	int  **a, int  *nx, int *ny, char *format);
extern void makefits(FILE *outfile, int nx, int ny, int bitpix, char *datatype);
extern void put_data(FILE *outfile, int a[], int nx, int ny, char *format);
extern int readint(MYFILE *infile);
extern void qread(MYFILE *infile, char *a, int n);
extern int myread(MYFILE *file, char buffer[], int n);
extern void writeint(MYFILE *outfile, int a);
extern void qwrite(MYFILE *outfile, char *a, int n);
extern int mywrite(MYFILE *file, char buffer[], int n);

/* general utility routines */
extern int test_swap(void);
extern void swap_bytes(unsigned char a[], int n);

#else

/* non-ANSI C versions */

/* bit I/O */
extern void start_inputing_bits();
extern int input_bit();
extern int input_nbits();
extern void start_outputing_bits();
extern void output_bit();
extern void output_nbits();
extern void done_outputing_bits();

/* hcompress routines */
extern void decode();
extern void digitize();
extern void dodecode();
extern void doencode();
extern void encode();
extern void hinv();
extern void hsmooth();
extern void htrans();
extern void qtree_decode();
extern void qtree_encode();

/* data I/O routines */
extern void fitspass();
extern void fitsread();
extern void get_data();
extern void makefits();
extern void put_data();
extern int readint();
extern void qread();
extern int myread();
extern void writeint();
extern void qwrite();
extern int mywrite();

/* general utility routines */
extern int test_swap();
extern void swap_bytes();

#endif

#endif
