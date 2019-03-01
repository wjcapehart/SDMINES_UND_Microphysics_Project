# ============================================================================
#
#  Makefile:  Test 1 - calling BoxMG ( BMG2_SymStd_SOLVE_boxmg )
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes an example that uses a direct
#  call to the BoxMG solver.  The problem is Poisson, discretized
#  on vertices, and with the source based on the discrete solution.
#  So this is not a real example, just handy for our debugging.  If 
#  you want to build this test, just enter
#
#  make
#
#  To run the test you will need to know how to start parallel jobs
#  on your machine, then redirect input from files in the input
#  directory, e.g.,
#
#  mpirun -np 4 < input/FileList-2x2-small
#
#  will run using the list of input files in "FileList-2x2-small".
#  To change the behavior, size of the problem etc., simply edit
#  the contents of these input files, or create a new set of input
#  files using the same format.
#
#  ----------------------
#  Environment Variables:
#  ----------------------
#
#  BOXMGdist   - location of BOXMG distribution
#  BOXMGlibp   - location of the BOXMG and BOXMG-EXTRAS library
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
  BOXMGdist = ../../..
endif

ifndef EXTRASdist
  EXTRASdist = $(BOXMGdist)/extras
endif

# -------------------------
#  Parallel
# -------------------------

BOXMG_MPI=yes

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# --------------------
#  Architecture:
# --------------------

include $(BOXMG_ARCHmake)/ARCH.$(BOXMG_ARCH)

# --------------------
#  Includes:
# --------------------

ifndef BOXMG_FINC
  BOXMG_FINC = -I./include -I$(BOXMGdist)/include -I$(EXTRASdist)/msg/include
endif

# ----------------
#  LIBRARIES:
# ----------------

LIBBMG   = boxmg_$(BOXMG_CLEVEL)
LIBBMGx  = boxmg-extras_$(BOXMG_CLEVEL)

BMGlink  = $(LPOPT)$(BOXMGlibp) $(LLOPT)$(LIBBMG)
BMGxlink = $(LPOPT)$(EXTRASlibp) $(LLOPT)$(LIBBMGx)

# ======================================================================
# ----------------------
#   Object Modules:
# ----------------------

MAIN = main.o diff.o 
PUTF =  ExactSol.o putf.o 

# ----------------------
#   Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS)

EXE = ex_direct_1

IN_DIR = input

OUT_DIR = output

OUT_FILES = $(OUT_DIR)/out* $(OUT_DIR)/Stencil.* $(OUT_DIR)/Restrict.*

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

test: BEGIN $(EXE)

$(EXE): $(OBJS) library $(LIBCHEK)
	@echo "Linking objects .... " ;\
	$(F77_LD) $(F77_LFLAGS) $(OBJS) $(LIBS) -o $(EXE); \
	echo "The executable $(EXE) has been created."

BEGIN:
	@echo "Building/checking the executable $(EXE) ... "

library:
	@( $(CD) $(BOXMGdist); make $(MFLAGS) library )

clean:
	@echo "Removing objects .... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE) and output files ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUT_FILES)

# ============================================================================








