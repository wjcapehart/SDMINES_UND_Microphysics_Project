# ============================================================================
#
#  Makefile:  Tests for BOXMG
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes the tests for the two and three
#  dimensional black box codes.  Assuming you have built and installed 
#  the library (see README) at the command line enter
#
#  make 
#
#  the output from each test is stored in a file in the tests subdirectory.
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE = make
endif

ifndef BOXMGdist
  BOXMGdist = ..
endif

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# --------------------
#  Architecture:
# --------------------

include $(BOXMG_ARCHmake)/ARCH.$(BOXMG_ARCH)

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

tests:  test-2D test-3D

test-2D:
	@( $(CD) $(TESTS_2D); $(MAKE) $(MFLAGS) )

test-3D:	
	@( $(CD) $(TESTS_3D); $(MAKE) $(MFLAGS) )

clean:
	@( $(CD) $(TESTS_2D); $(MAKE) $(MFLAGS) clean )
	@( $(CD) $(TESTS_3D); $(MAKE) $(MFLAGS) clean )

distclean:
	@( $(CD) $(TESTS_2D); $(MAKE) $(MFLAGS) distclean )
	@( $(CD) $(TESTS_3D); $(MAKE) $(MFLAGS) distclean )
