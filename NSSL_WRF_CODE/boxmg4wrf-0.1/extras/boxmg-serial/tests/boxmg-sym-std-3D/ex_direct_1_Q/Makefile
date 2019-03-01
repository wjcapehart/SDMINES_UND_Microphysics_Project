# ============================================================================
#
#  Makefile:  BoxMG Solver Test 2 ( BMG3_SER_SymStd_SOLVE_boxmg )
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes an example that uses a direct
#  call to the BoxMG solver.  The problem is provides a mechanism
#  to define regions of a domain with different coefficients, thereby
#  allowing testing with discontinuous coefficents.  A similar test
#  is provided for two dimensions.  Just type 
#
#  make
#
#  the required input is redirected from $(INPUT) and the output is 
#  redirected into the file $(OUTPUT).   
#
#  This program demonstrates a modular setup of PDE and solver
#  parameters.  In addition, it demonstrates the use of the
#  BMG3_SER_SymStd_PtrWORK to compute the required workspace dimensions,
#  and hence, run with variable problem sizes.
#
#  -----------------------
#  Files:
#  -----------------------
#
#    main.f                 - driver
#    EX_SETUP_BMG_SER_parms     - setup BMG cycle parameters
#    EX_SETUP_PDE_parms     - setup PDE parameters for discretization
#    EX_SETUP_GRID_size     - setup grid size and parameters.
#    putf.f                 - discretization of the PDE
#    da.f                   - diffusion tensor (1,1) entry
#    dr.f                   - diffusion tensor (2,2) entry
#    dz.f                   - diffusion tensor (3,3) entry
#    rhs.f                  - the source and removal
#
#  -----------------------
#  Environment Variables:
#  -----------------------
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
  BOXMG_FINC = -I$(BOXMGdist)/include -I./include
endif

# ----------------
#  LIBRARIES:
# ----------------

LIBBMG   = boxmg_$(BOXMG_CLEVEL)
LIBBMGx  = boxmg-extras_$(BOXMG_CLEVEL)

BMGlink  = $(LPOPT)$(BOXMGlibp) $(LLOPT)$(LIBBMG)
BMGxlink = $(LPOPT)$(EXTRASlibp) $(LLOPT)$(LIBBMGx)

# ----------------------------------------------------------------------
#  OBJECT MODULES:
# ----------------------

INCS = ./include/BMG_SER_dim_parms.h ./include/common2.h

MAIN =  main.o               \
	EX_SETUP_BMG_SER_parms.o \
	EX_SETUP_PDE_parms.o \
	EX_SETUP_GRID_size.o

PUTF =  putf.o rhs.o da.o dr.o dz.o

# ----------------------
#  Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS)

EXE  = ex_direct_1_Q

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


