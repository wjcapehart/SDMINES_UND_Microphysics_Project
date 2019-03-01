# ============================================================================
#
#  Makefile:  Tests for the two-dimensional BOXMG
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes the tests for the two-dimensional
#  black box multigrid code.  Assuming you have built the library 
#  (see README) at the command line enter
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

tests: tests-np tests-p

tests-np:  
#	@( $(CD) ./ex_bbmg_1; $(MAKE) $(MFLAGS) )
#	@( $(CD) ./ex_bbmg_2; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_direct_1; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_direct_1_Q; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_direct_pcg; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_sym_1; $(MAKE) $(MFLAGS) )

tests-p:  
#	@( $(CD) ./ex_p_bbmgp_1; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_p_direct_1; $(MAKE) $(MFLAGS) )

clean: banner clean-np clean-p

clean-np:
#	@( $(CD) ./ex_bbmg_1; $(MAKE) $(MFLAGS) clean )
#	@( $(CD) ./ex_bbmg_2; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_direct_1; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_direct_1_Q; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_direct_pcg; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_sym_1; $(MAKE) $(MFLAGS) clean )

clean-p:
#	@( $(CD) ./ex_p_bbmgp_1; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_p_direct_1; $(MAKE) $(MFLAGS) clean )

distclean: banner distclean-np distclean-p

distclean-np:
#	@( $(CD) ./ex_bbmg_1; $(MAKE) $(MFLAGS) distclean )
#	@( $(CD) ./ex_bbmg_2; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_direct_1; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_direct_1_Q; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_direct_pcg; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_sym_1; $(MAKE) $(MFLAGS) distclean )

distclean-p: 
#	@( $(CD) ./ex_p_bbmgp_1; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_p_direct_1; $(MAKE) $(MFLAGS) distclean )


banner:
	@echo "TWO-DIMENSIONAL:"

