# Makefile for HMMER
# 
##########
# HMMER - Biological sequence analysis with HMMs
# Copyright (C) 1992-1995 Sean R. Eddy
#
#   This source code is distributed under the terms of the 
#   GNU General Public License. See the files COPYING and 
#   LICENSE for details.
#    
###########
# Changed by wym6912 
# For anaysising the function of HMM Training
###########

## where you want things installed.
#  In general BINDIR and MANDIR are all you need to think about 
#  changing in the Makefile.
#
prefix = /usr/local
BINDIR = $(prefix)/bin
MANDIR = $(prefix)/man

#######
## You only need to change stuff below this if
## you experience compilation difficulties, or if you
## want to change the default compiler or compiler options.
#######

RELEASE     = "1.8.5"
RELEASEDATE = "February 2006"

## your compiler. This must be an ANSI-compliant compiler.
## In particular, SunOS cc is not ANSI. On SunOS, you should use gcc.
##

CC = gcc

## any special compiler flags you want
##
CFLAGS = -O3 -g -Wall -fopenmp
#CFLAGS = -O -uniproc -w   # Alliant FX/2800

## how to install the man pages 
## either cp -- to copy unformatted man page
## or a script with identical syntax to cp, to format & install formatted page
## Alternatively, you can set INSTMAN to something harmless like
## "echo", and hand-install the man pages.
##
INSTMAN   = cp
# this is fine on most machines 
MANSUFFIX = 1

#######
## You should not need to modify below this line.
#######
SHELL  = /bin/sh
LIBS   = -lm
# -ltcmalloc -lprofiler -lunwind

MANS  = hmmt

PROGS = hmmt

OBJ = sqerror.o sqio.o states.o sre_ctype.o sre_math.o sre_string.o misc.o prior.o \
      types.o iupac.o dayhoff.o alignio.o selex.o maxmodelmaker.o viterbi.o hmmio.o \
	  trace.o msf.o forback.o weeviterbi.o align.o saviterbi.o weight.o cluster.o

all: $(PROGS)
	@echo done.

hmmt: $(OBJ) hmmt.o
	$(CC) $(CFLAGS) -o hmmt hmmt.o $(OBJ) $(LIBS)

install: all 
	test -d $(BINDIR) || mkdir -p $(BINDIR)
	test -d $(MANDIR)/man$(MANSUFFIX) || mkdir -p $(MANDIR)/man$(MANSUFFIX)
	cp $(PROGS) $(BINDIR)/
	for manfile in $(MANS); do\
		$(INSTMAN) $$manfile.man $(MANDIR)/man$(MANSUFFIX)/$$manfile.$(MANSUFFIX);\
	done

uninstall:
	for prog in $(PROGS); do\
		rm $(BINDIR)/$$prog;\
	done
	for manfile in $(MANS); do\
		rm $(MANDIR)/man$(MANSUFFIX)/$$manfile.$(MANSUFFIX); \
	done

clean:
	-rm -f *.o *~ $(PROGS)

.o:
	$(CC) $(CFLAGS) -c $<		
