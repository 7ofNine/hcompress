# makefile for hcompress programs
# R. White, 11 May 2009
#
#-----------------------------------------------------------

AUTOMAKE_OPTIONS = gnu
LDADD =
bin_PROGRAMS = hcomp hdecomp
bin_SCRIPTS = hcompress hdecompress fcompress fdecompress
EXTRA_DIST = $(bin_SCRIPTS)

hcomp_SOURCES = hcomp.c hcompress.h \
	bit_output.c \
	digitize.c \
	doencode.c \
	encode.c \
	fitspass.c \
	fitsread.c \
	get_data.c \
	htrans.c \
	qtree_encode.c \
	qread.c \
	qwrite.c \
	swap_bytes.c

hdecomp_SOURCES = hdecomp.c hcompress.h \
	bit_input.c \
	decode.c \
	dodecode.c \
	fitspass.c \
	hinv.c \
	hsmooth.c \
	makefits.c \
	put_data.c \
	qread.c \
	qwrite.c \
	qtree_decode.c \
	swap_bytes.c

noinst_HEADERS = hcompress.h
