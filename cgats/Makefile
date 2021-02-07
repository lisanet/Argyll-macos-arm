# UNIX style makefile, for icclib and friends.
# "include" the right environment for your system,
# by uncommenting the appropriate line:

# Microsoft C++, WinNT setup
#include Makefile.WNT

# IBM C++, WinNT setup
include Makefile.IBMNT

# Generic UNIX setup
#include Makefile.UNIX

# Apple OSX
#include Makefile.OSX

###############################

#Compile with separate stndard malloc & file io
#CCDEFINES = $(DEFFLAG)SEPARATE_STD

#Set optimisation on
CCFLAGS = $(CCFLAGSDEF) $(CCOPTFLAG) $(CCDEFINES)

#Set debugging on
#CCFLAGS = $(CCFLAGSDEF) $(CCDEBUGFLAG) $(CCDEFINES)

STDHDRS = $(STDHDRSDEF)
LINKFLAGS = $(LINKFLAGSDEF) $(LINKDEBUGFLAG)

all:: libcgats$(SUFLIB) pars$(SUFEXE) cgats$(SUFEXE)


# Separate for executables
parsstd$(SUFOBJ): parsstd.c pars.h
	$(CC) $(CCOF)parsstd$(SUFOBJ) parsstd.c

cgatsstd$(SUFOBJ): cgatsstd.c cgats.h pars.h
	$(CC) $(CCOF)cgatsstd$(SUFOBJ) cgatsstd.c


pars$(SUFOBJ): pars.c pars.h
	$(CC) $(CCOF)pars$(SUFOBJ) pars.c

cgats$(SUFOBJ): cgats.c cgats.h pars.h
	$(CC) $(CCOF)cgats$(SUFOBJ) cgats.c


libcgats$(SUFLIB): pars$(SUFOBJ) cgats$(SUFOBJ)
	$(LIBU) $(LIBOF)libcgats$(SUFLIB) cgats$(SUFOBJ) pars$(SUFOBJ)
	$(RANLIB) libcgats$(SUFLIB)


sa_pars$(SUFOBJ): pars.c pars.h
	$(CC) $(CCOF)sa_pars$(SUFOBJ) $(DEFFLAG)STANDALONE_TEST pars.c

pars$(SUFEXE): sa_pars$(SUFOBJ) parsstd$(SUFOBJ)
	$(LINK) $(LINKOF)pars$(SUFEXE) sa_pars$(SUFOBJ) parsstd$(SUFOBJ) $(LINKLIBS)


sa_cgats$(SUFOBJ): cgats.c cgats.h pars.h
	$(CC) $(CCOF)sa_cgats$(SUFOBJ) $(DEFFLAG)STANDALONE_TEST cgats.c

cgats$(SUFEXE): sa_cgats$(SUFOBJ) pars$(SUFOBJ) parsstd$(SUFOBJ) cgatsstd$(SUFOBJ)
	$(LINK) $(LINKOF)cgats$(SUFEXE) sa_cgats$(SUFOBJ) pars$(SUFOBJ) parsstd$(SUFOBJ) \
                                    cgatsstd$(SUFOBJ) $(LINKLIBS)


