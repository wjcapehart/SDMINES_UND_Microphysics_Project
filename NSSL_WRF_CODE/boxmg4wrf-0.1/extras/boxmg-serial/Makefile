# ============================================================================
#
#  Makefile:  Dendy's Symmetric Black Box Multigrid (2D and 3D)
#
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This is a makefile for Joel Dendy's (Symmetric, Standard Coarsening)
#  Black Box multigrid. To build and install the library you may need to
# 
#    - set the appropriate environment variables 
#      (see below and/or INSTALL)
#    - customize the compiler flags as necessary
#      (see INSTALL)
#
#  However, the build system provides defaults that should work on most
#  systems (particularly, Linux/GNU systems). The most common customization
#  is to specify a compiler suite that you'd like to use (see INSTALL).
#  Now, at the prompt simply type
#
#    make
#
#  The libraries ( libboxmg_$(CLEVEL).a, libboxmg-extras_$(CLEVEL).a )
#  will be built and by default put in ./lib.
#
#  If you want to run the tests, do the following
#
#    make check
#
#  or change directories to the appropriate test
# 
#    cd tests/boxmg-sym-std-3D/ex_direct_1
#    make 
#
#  This will check that the library is installed, up to date, and then run 
#  the tests.  Examine the output of the tests to confirm that everything
#  is working.
#
#  -----------------------
#  Environment Variables:
#  -----------------------
#
#  BOXMGdist   - location of BOXMG distribution
#  BOXMGlibp   - location to build the BOXMG library
#            
#  OS          - Operating system (e.g., Linux, OSF1 )
#  CPU         - CPU (e.g., intel, AMD64 )
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

ifndef BOXMGdist
  BOXMGdist = .
endif

#
# Only use for building serial tests, so we need to point to
# the extras on level up.
#
ifndef EXTRASdist
  EXTRASdist = ..
endif

# ============================================================================

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

#
# Library build path
#
BOXMGlibp = $(BOXMGdist)/lib

#
# BoxMG tests 
#
CHECK = $(BOXMGdist)/tests

# -------------------------
#  Extras build options
# -------------------------

EXTRAS_OPTS = BOXMG=no MSGLIB=no EXTRASlibp=$(BOXMGdist)/lib

# -------------------------
#  Select Components
# -------------------------

PIECES = BUILD_2D 

ifndef SER_3D
  SER_3D=yes
endif
ifeq ($(SER_3D),yes)
  PIECES += BUILD_3D
endif

PIECES += BUILD_UTILS

ifndef SER_EXTRAS
  SER_EXTRAS=no
endif
ifeq ($(SER_EXTRAS),yes)
  PIECES += lib-extras
endif

# ========================================================================= 
# -------------------------
#   TARGETS:
# -------------------------

makewithflags:
	@echo ""
	@echo " * Building the SERIAL boxmg library .... "
	@( $(MAKE) $(MFLAGS) library )
	@echo ""

environment:
	@echo ""
	@echo "  Environment Variables "
	@echo ""
	@echo "  BOXMG_OS        = $(BOXMG_OS)"
	@echo "  BOXMG_CPU       = $(BOXMG_CPU)"
	@echo "  BOXMG_ARCH      = $(BOXMG_ARCH)"
	@echo ""
	@echo "  BOXMG_COMPILER  = $(BOXMG_COMPILER)"
	@echo "  BOXMG_GCC_VER   = $(BOXMG_GCC_VER)"
	@echo "  BOXMG_CLEVEL    = $(BOXMG_CLEVEL)"
	@echo ""
	@echo "  FC              = $(FC)"
	@echo "  FFLAGS          = $(FFLAGS)"
	@echo "  F90             = $(F90)"
	@echo "  F90FLAGS        = $(FFLAGS)"
	@echo ""
	@echo "  CC              = $(CC)"
	@echo "  CFLAGS          = $(CFLAGS)"
	@echo "  CPPFLAGS        = $(CPPFLAGS)"
	@echo ""
	@echo "  SRC2D           = $(SRC2D)"
	@echo "  SRC3D           = $(SRC3D)"
	@echo "  SRC_U           = $(SRC_U)"


# -------------------------
#  BoxMG library
# -------------------------

library: $(PIECES)

BUILD_2D:
	@( $(CD) $(SRC2D); $(MAKE) $(MFLAGS) )

BUILD_3D:
	@( $(CD) $(SRC3D); $(MAKE) $(MFLAGS) )

BUILD_UTILS:
	@( $(CD) $(SRC_U); $(MAKE) $(MFLAGS) )


clean: check-clean library-banner
	@( $(CD) $(SRC2D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(SRC3D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(SRC_U); $(MAKE) $(MFLAGS) clean )
	@echo ""

distclean: check-distclean library-banner
	@( $(CD) $(SRC2D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(SRC3D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(SRC_U); $(MAKE) $(MFLAGS) clean )
	@echo ""
	@echo "Removing the libraries in $(BOXMGlibp) "
	@( $(CD) $(BOXMGlibp); $(RM) -f *.a )
	@echo ""

library-banner:
	@echo ""
	@echo "SERIAL LIBRARY SOURCE:"
	@echo ""
	@echo "Removing objects from: "

# -------------------------
#  Extras Library
# -------------------------

lib-extras:
	@( $(CD) $(EXTRASdist); $(MAKE) $(MFLAGS) $(EXTRAS_OPTS) )

lib-extras-clean: lib-extras-banner
	@( $(CD) $(EXTRASdist); $(MAKE) $(MFLAGS) $(EXTRAS_OPTS) clean )

lib-extras-distclean: lib-extras-banner
	@( $(CD) $(EXTRASdist); $(MAKE) $(MFLAGS) $(EXTRAS_OPTS) distclean )

lib-extras-banner:
	@echo ""
	@echo "EXTRAS SUBDIRECTORY:"
	@echo ""

# -------------------------
#  Tests
# -------------------------

check:
	@( $(CD) $(CHECK); $(MAKE) $(MFLAGS) )

check-clean: check-banner
	@( $(CD) $(CHECK); $(MAKE) $(MFLAGS) clean )

check-distclean: check-banner
	@( $(CD) $(CHECK); $(MAKE) $(MFLAGS) distclean )

check-banner:
	@echo ""
	@echo "SERIAL TESTS:"
	@echo ""

# ===========================================================================

