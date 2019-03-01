# ============================================================================
#
#  Makefile:  MPI-based Black Box Multigrid 
#
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This is a makefile for the extras that are needed for the MPI-based
#  implementation of Joel Dendy's (Symmetric with Standard Coarsening) 
#  Black Box multigrid.   Specifically, we have included (pieces of)
#  
#     MSG (Message-Passing Tools for Structured Grid Communications)
#     
#     http://www.cerca.umontreal.ca/~malevsky/MSG/MSG.html
#
#  There have been minor modifications, if you're concerned run
#  a recursive diff against the original package.  
#
#  In addition, pieces of the well known packages LAPACK 
#  (with self-contained BLAS) and a simple portable random number 
#  generator from SLATEC, are included.  These codes are all freely 
#  available from netlib. Finally, a serial version of the BoxMG 
#  library is included to handle problems on a large number of 
#  processors (i.e., a problem that generates are large system on 
#  the coarsest parallel grid).
#
#  These routines are provided only for your convenience, optimized 
#  versions may be available on your system, updates may be available
#  from the web.  
#
#  To build and install the extras library make sure you
#  have done the following,
# 
#    - set the appropriate environment variables 
#      (see below and/or README)
#    - customized the compiler flags as necessary
#      (see README)
#
#  Then at the prompt type the commands,
#
#    make
#    make install
#
#  At the moment there are know tests for this portion
#  of the library.  If you have concerns there may be appropriate
#  tests available from netlib.
#
#  Environment Variables:
#  ----------------------
#
#  EXTRASdist   - location of boxmg-extras distribution
#  EXTRASlibp   - location to build the boxmg-extras library
#  EXTRASlibi   - location to install the library
#            
#  OS          - Operating system (e.g., Linux, SunOS )
#  CPU         - CPU (e.g., intel, mips, alpha, sparc)
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

#
#  BoxMG MPI distribution
#
#  Setting the full path is necessary for packages that
#  use make recursively with multi-level source trees 
#  ( e.g., the serial boxmg code ).
#
ifndef BOXMGdist
  BOXMGdist = $(shell cd ../; pwd )
endif

ifndef EXTRASdist
  EXTRASdist = $(BOXMGdist)/extras
endif

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# -------------------------
#  Compiler, Linker:
# -------------------------

include $(BOXMG_ARCHmake)/ARCH.$(BOXMG_ARCH)

# -------------------------
#  Path Definitions:
# -------------------------

MSG   = $(EXTRASdist)/msg/src
BMG   = $(EXTRASdist)/boxmg-serial
LAPK  = $(EXTRASdist)/lapack
RAND  = $(EXTRASdist)/random
MACH  = $(EXTRASdist)/machine/$(BOXMG_ARCH)

# -------------------------
#  Select Components
# -------------------------

PIECES = BUILD_MACHINE 

ifndef BUILD_MSG 
  MSGLIB=yes
endif
ifeq ($(MSGLIB),yes)
  PIECES += BUILD_MSG
endif

ifndef BOXMG
  BOXMG=yes
endif
ifeq ($(BOXMG),yes)
  PIECES += BUILD_BOXMG
  #
  # Set library target
  #
  EXTRAS_BOXMGlibp = BOXMGlibp=$(BOXMGlibp)
  EXTRAS_BOXMGlibn = BOXMGlibn=$(EXTRASlibn)
  LIB_BOXMG = $(EXTRAS_BOXMGlibp) $(EXTRAS_BOXMGlibn)
  #
  DEF_COMPILERS = BOXMG_COMPILER=$(BOXMG_COMPILER)
  DEF_CLEVEL = BOXMG_CLEVEL=$(BOXMG_CLEVEL)
  COMP_OPTS = $(DEF_COMPILERS) $(DEF_CLEVEL) WRAPPER=no
  #
  BUILD_BOXMG_OPTS  = BOXMG_ARCHmake=$(BOXMG_ARCHmake)
  BUILD_BOXMG_OPTS += $(COMP_BOXMG) $(LIB_BOXMG)
  CLEAN_BOXMG_OPTS  = BOXMG_ARCHmake=$(BOXMG_ARCHmake)
  CLEAN_BOXMG_OPTS += $(LIB_BOXMG)
  #
endif

ifndef LAPACK
  LAPACK=yes
endif
ifeq ($(LAPACK),yes) 
  PIECES += BUILD_LAPACK
endif

ifndef RANDOM
  RANDOM=yes
endif
ifeq ($(RANDOM),yes)
  PIECES += BUILD_RANDOM
endif

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: $(PIECES)

BUILD_MSG:
	@( $(CD) $(MSG); $(MAKE) $(MFLAGS) )

BUILD_BOXMG:
	@( $(CD) $(BMG); $(MAKE) $(MFLAGS) $(BUILD_BOXMG_OPTS) )

BUILD_LAPACK:
	@( $(CD) $(LAPK); $(MAKE) $(MFLAGS) )

BUILD_RANDOM:
	@( $(CD) $(RAND); $(MAKE) $(MFLAGS) )

BUILD_MACHINE:
	@( $(CD) $(MACH); $(MAKE) $(MFLAGS) )

install: 
	@echo "Installing $(EXTRASlibp)/$(NEWLIB) with ... "
	@echo "$(MV) $(EXTRASlibp)/$(NEWLIB) $(EXTRASlibi)"
	@$(MV) $(EXTRASlibp)/$(NEWLIB) $(EXTRASlibi)

uninstall: 
	@echo "Uninstalling $(NEWLIB) from $(EXTRASlibi) with ... "
	@echo "$(MV) $(EXTRASlibi)/$(NEWLIB) $(EXTRASlibp)"
	@$(MV) $(EXTRASlibi)/$(NEWLIB) $(EXTRASlibp)

clean:
	@echo ""
	@echo "Cleaning ... "
	@echo ""
	@echo "*Cleaning MSG ... "
	@echo ""
	@( $(CD) $(MSG); $(MAKE) $(MFLAGS) clean )
	@echo ""
	@echo "*Cleaning SERIAL BoxMG ... "
	@echo ""
	@( $(CD) $(BMG); $(MAKE) $(MFLAGS) $(CLEAN_BOXMG_OPTS) clean )
	@echo ""
	@echo "*Cleaning LAPACK ... "
	@echo ""
	@( $(CD) $(LAPK); $(MAKE) $(MFLAGS) clean )
	@echo ""
	@echo "*Cleaning random ... "
	@echo ""
	@( $(CD) $(RAND); $(MAKE) $(MFLAGS) clean )
	@echo ""
	@echo "*Cleaning machine ... "
	@echo ""
	@( $(CD) $(MACH); $(MAKE) $(MFLAGS) clean )

distclean: clean
	@echo ""
	@echo "Removing library: " $(EXTRASlibn)
	@/bin/rm -f $(EXTRASlibp)/$(EXTRASlibn)
	@echo ""

# ===========================================================================




