# ============================================================================
#
#  Makefile:  Test 4 - direct call to BoxMG solver( BMG2_SER_SymStd_SOLVE_boxmg )
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes an example that uses a direct
#  call to the BoxMG solver.  The problem is provides a mechanism
#  to define regions of a domain with different coefficients, thereby
#  allowing testing with discontinuous coefficents.  A similar test
#  is provided for three dimensions.  Just type 
#
#  make
#
#  the required input is redirected from $(INPUT) and the output is 
#  redirected into the file $(OUTPUT).  
#
#  ----------------------
#  Environment Variables:
#  ----------------------
#
#  BOXMGdist   - location of BOXMG distribution
#  BOXMGlibp   - location of the BOXMG and BOXMG-EXTRAS library
#            
#  OS          - Operating system (e.g., Linux, OSF1 )
#  CPU         - CPU (e.g., intel, AMD64)
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
  BOXMG_FINC = -I./include -I$(BOXMGdist)/include
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

INCS = ./include/BMG_SER_dim_parms.h ./include/common2.h

MAIN =  main.o               \
	EX_SETUP_BMG_SER_parms.o \
	EX_SETUP_PDE_parms.o \
	EX_SETUP_GRID_size.o

PUTF =  putf.o rhs.o da.o dr.o

# ----------------------
#   Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS)

EXE = ex_direct_1_Q

IN_DIR = input
INPUTa = $(IN_DIR)/FileList-small

OUT_DIR = output
OUTPUTa = $(OUT_DIR)/output-small.txt

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

$(OBJS): $(INCS)

BEGIN:
	@echo "Building/checking the executable $(EXE) ... "

library:
	@( $(CD) $(BOXMGdist); make $(MFLAGS) library )

clean:
	@echo "Removing objects for $(EXE) ... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE) and the output files ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUT_FILES)

# ============================================================================

