# ============================================================================
#
#  Makefile:  Tests for BoxMG
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds the tests for the two- and three-dimensional 
#  BoxMG codes.  Assuming you have built and installed the library 
#  (see README) at the command line enter
#
#  make 
#
#  Running the code depends on your parallel environment.  However, some
#  common scenarios are given in the next level down, as well as in each
#  tests own Makefile.
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

test:  
	@( $(MAKE) test-2D $(MFLAGS) )
	@( $(MAKE) test-3D $(MFLAGS) )

test-2D:
	@( $(CD) ./boxmg-sym-std-2D; $(MAKE) $(MFLAGS) )

test-3D:	
	@( $(CD) ./boxmg-sym-std-3D; $(MAKE) $(MFLAGS) )

clean:
	@( $(CD) ./boxmg-sym-std-2D; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./boxmg-sym-std-3D; $(MAKE) $(MFLAGS) clean )

distclean:
	@( $(CD) ./boxmg-sym-std-2D; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./boxmg-sym-std-3D; $(MAKE) $(MFLAGS) distclean )


