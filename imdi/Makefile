
# Boilerplate Makefile for compiling imdi

# Copyright 2000 - 2007 Graeme W. Gill
# This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
# see the License.txt file for licencing details.

# "include" the right environment for your system,
# by uncommenting the appropriate line:

# Microsoft VC++, WinNT setup
include Makefile.WNT

# Generic UNIX setup
#include Makefile.UNIX

# Apple OS X
#include Makefile.OSX


###############################

#Set optimisation on
CCFLAGS = $(CCFLAGSDEF) $(CCOPTFLAG) $(CCDEFINES) $(BCONFIG)

#Set debugging on
#CCFLAGS = $(CCFLAGSDEF) $(CCDEBUGFLAG) $(CCDEFINES) $(BCONFIG)
# debugging #define
#CCFLAGS = $(CCFLAGSDEF) $(CCDEBUGFLAG) $(CCDEFINES) $(DEFFLAG)DEBUG
LINKFLAGS = $(LINKFLAGSDEF) $(LINKDEBUGFLAG)

STDHDRS = $(INCFLAG)$(STDHDRSDEF)

all:: libimdi$(SUFLIB)

# Used by both code generator and runtime
imdi_make$(SUFEXE): imdi_make$(SUFOBJ) imdi_gen$(SUFOBJ) cgen$(SUFOBJ)
	$(LINK) $(LINKOF)imdi_make$(SUFEXE) imdi_make$(SUFOBJ) imdi_gen$(SUFOBJ) cgen$(SUFOBJ)


# The code generator program

cgen$(SUFOBJ): cgen.c imdi_utl.h imdi_arch.h imdi_gen.h imdi_tab.h
	$(CC) cgen.c

imdi_gen$(SUFOBJ): imdi_gen.c imdi_utl.h imdi_arch.h imdi_gen.h
	$(CC) imdi_gen.c

# Generate the kernel files
imdi_k.h imdi_k.c : imdi_make$(SUFEXE)
	.$(SLASH)imdi_make$(SUFEXE) 


# imdi runtime library

imdi$(SUFOBJ): imdi.c imdi.h imdi_tab.h imdi_k.h imdi_k.c
	$(CC) imdi.c

libimdi$(SUFLIB): imdi$(SUFOBJ) imdi_tab$(SUFOBJ)
	$(LIBU) $(LIBOF)$@ imdi$(SUFOBJ) imdi_tab$(SUFOBJ)
	$(RANLIB) libimdi$(SUFLIB)





