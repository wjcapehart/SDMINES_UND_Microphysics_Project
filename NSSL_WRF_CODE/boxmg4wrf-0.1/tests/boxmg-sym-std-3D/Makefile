# ============================================================================
#
#  Makefile:  Tests for the Three-Dimensional BoxMG solver
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds the tests for the three-dimensional Black Box 
#  multigrid code (BoxMG).  The target in each subdirectory checks that 
#  the library is current, so just enter
#
#  make 
#
#  To run an example you need to do the following:
#
#    - change to the directory you're interested, e.g., 
#
#      cd boxmg-sym-std-2D/ex_direct_1_Q
#
#    - look a the list of files in the "input" directory
#
#    - use whatever mechanism your MPI-based system uses to start
#      parallel jobs and redirect standard in from input/FileName.
#      For example, 
#
#      mpirun -np 4 < input/FileList-2x2-small ex_direct_1_Q
#
#  There are many ways to specify the list of machines or nodes
#  that the job will run on, as well as the interconnect that the
#  job will use for MPI communication.  This information you need to
#  obtain for your system, and translate into additional options here.
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
	@( $(CD) ./ex_direct_1_f90; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_direct_pcg; $(MAKE) $(MFLAGS) )
	@( $(CD) ./ex_direct_pcg_f90; $(MAKE) $(MFLAGS) )
#	@( $(CD) ex_sym_1; $(MAKE) $(MFLAGS) )

clean: banner
	@( $(CD) ./ex_direct_1; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_direct_1_f90; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_direct_pcg; $(MAKE) $(MFLAGS) clean )
	@( $(CD) ./ex_direct_pcg_f90; $(MAKE) $(MFLAGS) clean )
#	@( $(CD) ex_sym_1; $(MAKE) $(MFLAGS) clean )

distclean: banner
	@( $(CD) ./ex_direct_1; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_direct_1_f90; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_direct_pcg; $(MAKE) $(MFLAGS) distclean )
	@( $(CD) ./ex_direct_pcg_f90; $(MAKE) $(MFLAGS) distclean )
#	@( $(CD) ex_sym_1; $(MAKE) $(MFLAGS) distclean )

banner:
	@echo "THREE-DIMENSIONAL:"
