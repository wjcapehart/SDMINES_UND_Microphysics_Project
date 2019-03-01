# ============================================================================
#
#  Makefile:  Test 3 - direct call to BoxMG ( BMG2_SER_SymStd_SOLVE_boxmg )
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes an example that uses a direct
#  call to the BoxMG solver.  The actual problem is the same as 
#  test 1 (ex_bbmg_1).  This will check that the libraries have the 
#  necessary objects, so just enter
#
#  make
#
#  the required input is redirected from $(INPUT) and the output is 
#  redirected into the file $(OUTPUT).  This program demonstrates the
#  use of the BMG2_SER_SymStd_PtrWORK to compute the required workspace
#  dimensions, and hence, run with variable problem sizes. 
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
  BOXMG_FINC = -I$(BOXMGdist)/include
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
PUTF =  putf.o 

# ----------------------
#   Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS)

EXE = ex_direct_1

IN_DIR = input
INPUTa = $(IN_DIR)/input-30x60.txt

OUT_DIR = output
OUTPUTa = $(OUT_DIR)/output-30x60.txt

OUT_FILES = $(OUT_DIR)/out* $(OUT_DIR)/Stencil.* $(OUT_DIR)/Restrict.*

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

test: BEGIN $(EXE)
	@echo "Executing the command ./$(EXE) < $(INPUTa) > $(OUTPUTa)" ; \
	./$(EXE) < $(INPUTa) > $(OUTPUTa)

$(EXE): $(OBJS) library
	@echo " - linking objects ... " ;\
	$(F77_LD) $(F77_LFLAGS) $(OBJS) $(LIBS) -o $(EXE); \
	echo " - the executable $(EXE) has been created."

BEGIN:
	@echo "Building/checking the executable $(EXE) ... "

library:
	@( $(CD) $(SRC2D); make $(MFLAGS) )

clean:
	@echo "Removing objects for $(EXE) ... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE) and the output files ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUT_FILES)

# ============================================================================

