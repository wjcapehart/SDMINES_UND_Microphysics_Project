# ============================================================================
#
#  Makefile:  Tests for the three-dimensional BoxMG
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes the tests for the three-dimensional
#  black box multigrid code.  On the command line you can just type
#
#  make 
#
#  the output from each test is stored in a file in the tests subdirectory.
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

ifndef BOXMGdist
  BOXMGdist = ../..
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

tests:
	@( $(CD) ./ex_direct_1; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_direct_1_Q; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_dir_pcg; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_sym_1; $(MAKE) $(MFLAGS) )

clean:  banner
	@( $(CD) ./ex_direct_1; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_direct_1_Q; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_dir_pcg; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_sym_1; $(MAKE) $(MFLAGS) clean )

distclean: banner
	@( $(CD) ./ex_direct_1; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_direct_1_Q; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_dir_pcg; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_sym_1; $(MAKE) $(MFLAGS) distclean )

banner:
	@echo "THREE-DIMENSIONAL:"
