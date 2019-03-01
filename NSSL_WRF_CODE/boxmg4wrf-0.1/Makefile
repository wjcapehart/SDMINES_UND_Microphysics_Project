# ============================================================================
#
#  Makefile:  MPI-based Symmetric Black Box Multigrid (2D and 3D)
#
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This is a makefile for a MPI-based parallel implementation of Joel Dendy's
#  (Symmetric, Standard Coarsening) Black Box multigrid. To build and install 
#  the library you may need to
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
#  will be in ./lib.
#
#  If you want to run the tests, you first need to build them with
#
#    make check
#
#  or change directories to the desired test/example
# 
#    cd tests/boxmg-sym-std-3D/ex_direct_1
#    make 
#
#  This will check that the library is installed, up to date, and build
#  the executable.  Since, running on parallel systems depends on the 
#  local environment and batch system, you will need to execute the
#  test manually.  Notes are included in each tests subdirectory regarding
#  the input and output of the test.
#
#  Environment Variables:
#  ----------------------
#
#  BOXMGdist   - location of BOXMG distribution
#  BOXMGlibp   - location to build the BOXMG library
#            
#  OS          - Operating system (e.g., Linux, SunOS )
#  CPU         - CPU (e.g., intel, sun4u)
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

ifndef BOXMGdist
  BOXMGdist = .
endif

ifndef EXTRASdist
  EXTRASdist = $(BOXMGdist)/extras
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

# --------------------
#  Path Definitions:
# --------------------

#
# Library build path
#
BOXMGlibp = $(BOXMGdist)/lib

#
# BoxMG tests 
#
CHECK = $(BOXMGdist)/tests

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

makewithflags:
	@echo ""
	@echo "Building the boxmg library .... "
	@( $(MAKE) $(MFLAGS) library )
	@echo ""
	@echo "Building the boxmg-extras library .... "
	@( $(MAKE) $(MFLAGS) lib-extras )
	@echo ""

# ------------------------------------
#   BoxMG Environment Debugging:
# ------------------------------------

include $(BOXMGdist)/make/environment

# -------------------------
#  BoxMG library
# -------------------------

library:
	@( $(CD) $(SRC2D); $(MAKE) $(MFLAGS) MPI=yes )
	@( $(CD) $(SRC3D); $(MAKE) $(MFLAGS) MPI=yes )

clean: check-clean library-banner
	@( $(CD) $(SRC2D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(SRC3D); $(MAKE) $(MFLAGS) clean )
	@echo ""
	@echo "========================================================"
	@echo ""

distclean: check-distclean lib-extras-distclean library-banner
	@( $(CD) $(SRC2D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(SRC3D); $(MAKE) $(MFLAGS) clean )
	@echo ""
	@echo "Removing the libraries in $(BOXMGlibp) "
	@( $(CD) $(BOXMGlibp); $(RM) -f *.a )
	@echo ""
	@echo "========================================================"
	@echo ""

library-banner:
	@echo "========================================================"
	@echo "---------------------"
	@echo "LIBRARY SOURCE:"
	@echo "---------------------"
	@echo ""
	@echo "Removing objects from: "

# -------------------------
#  Extras Library
# -------------------------

lib-extras:
	@( $(CD) $(EXTRASdist); $(MAKE) $(MFLAGS) ) 

lib-extras-clean: lib-extras-banner
	@( $(CD) $(EXTRASdist); $(MAKE) $(MFLAGS) clean )

lib-extras-distclean: lib-extras-banner
	@( $(CD) $(EXTRASdist); $(MAKE) $(MFLAGS) distclean )

lib-extras-banner:
	@echo "========================================================"
	@echo "---------------------"
	@echo "EXTRAS SUBDIRECTORY: "
	@echo "---------------------"

# -------------------------
#  Tests
# -------------------------

check: 
	@( $(CD) $(CHECK); $(MAKE) $(MFLAGS) )

check-clean: check-banner
	@( $(CD) $(CHECK); $(MAKE) $(MFLAGS) clean )
	@echo ""

check-distclean: check-banner
	@( $(CD) $(CHECK); $(MAKE) $(MFLAGS) distclean )
	@echo ""

check-banner:
	@echo ""
	@echo "========================================================"
	@echo "---------------------"
	@echo "TESTS:               "
	@echo "---------------------"
	@echo ""

# ===========================================================================
