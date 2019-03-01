# ============================================================================
#
#  Makefile:  BoxMG Solver Test 1 ( BMG3_SER_SymStd_SOLVE_boxmg )
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
#  redirected into the file $(OUTPUT).  This program demonstrates the
#  use of the BMG3_SER_SymStd_PtrWORK to compute the required workspace
#  dimensions, and hence, run with variable problem sizes. 
#
#  -----------------------
#  Files:
#  -----------------------
#
#    main.f       - driver
#    putf.f       - discretization of the PDE
#    da.f         - diffusion tensor (1,1) entry
#    dr.f         - diffusion tensor (2,2) entry
#    dz.f         - diffusion tensor (3,3) entry
#    rhs.f        - the source and removal
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
  BOXMG_FINC = -I$(BOXMGdist)/include -I../../../include -I./include
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

MAIN = main.o 
PUTF = da.o dr.o dz.o putf.o rhs.o

# ----------------------
#  Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS) -llapack -lblas

EXE  = ex_direct_1

# ----------------------
#  Input Files:
# ----------------------

IN_DIR  = input

INPUTa  = $(IN_DIR)/id.9
INPUTb  = $(IN_DIR)/id.16

INPUTc  = $(IN_DIR)/id.17-FMG

INPUTcx = $(IN_DIR)/id.17-Xeps
INPUTcy = $(IN_DIR)/id.17-Yeps
INPUTcz = $(IN_DIR)/id.17-Zeps

INPUTd  = $(IN_DIR)/id.35x35x35
INPUTdx = $(IN_DIR)/id.17x35x35
INPUTdy = $(IN_DIR)/id.35x17x35
INPUTdz = $(IN_DIR)/id.35x35x17

# ---------------------- 
#  Output Files:
# ----------------------

OUT_DIR  = output

OUTPUTa  = $(OUT_DIR)/output.09.txt
OUTPUTb  = $(OUT_DIR)/output.16.txt

OUTPUTc  = $(OUT_DIR)/output.17-FMG.txt

OUTPUTcx = $(OUT_DIR)/output.17-Xeps.txt
OUTPUTcy = $(OUT_DIR)/output.17-Yeps.txt
OUTPUTcz = $(OUT_DIR)/output.17-Zeps.txt

OUTPUTd  = $(OUT_DIR)/output.35x35x35.txt
OUTPUTdx = $(OUT_DIR)/output.17x35x35.txt
OUTPUTdy = $(OUT_DIR)/output.35x17x35.txt
OUTPUTdz = $(OUT_DIR)/output.35x35x17.txt

OUT_FILES = $(OUT_DIR)/out* $(OUT_DIR)/Stencil.* $(OUT_DIR)/Restrict.*

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

test: BEGIN $(EXE)
#	@echo "Executing the command ./$(EXE) < $(INPUTa) > $(OUTPUTa)" ; \
#	./$(EXE) < $(INPUTa) > $(OUTPUTa)
#	@echo "Executing the command ./$(EXE) < $(INPUTb) > $(OUTPUTb)" ; \
#	./$(EXE) < $(INPUTb) > $(OUTPUTb)
#	@echo "Executing the command ./$(EXE) < $(INPUTc) > $(OUTPUTc)" ; \
#	./$(EXE) < $(INPUTc) > $(OUTPUTc)
#	@echo "Executing the command ./$(EXE) < $(INPUTcx) > $(OUTPUTcx)" ; \
#	./$(EXE) < $(INPUTcx) > $(OUTPUTcx)
#	@echo "Executing the command ./$(EXE) < $(INPUTcy) > $(OUTPUTcy)" ; \
#	./$(EXE) < $(INPUTcy) > $(OUTPUTcy)
#	@echo "Executing the command ./$(EXE) < $(INPUTcz) > $(OUTPUTcz)" ; \
#	./$(EXE) < $(INPUTcz) > $(OUTPUTcz)
#	@echo "Executing the command ./$(EXE) < $(INPUTd) > $(OUTPUTd)" ; \
#	./$(EXE) < $(INPUTd) > $(OUTPUTd)
#	@echo "Executing the command ./$(EXE) < $(INPUTdx) > $(OUTPUTdx)" ; \
#	./$(EXE) < $(INPUTdx) > $(OUTPUTdx)
#	@echo "Executing the command ./$(EXE) < $(INPUTdy) > $(OUTPUTdy)" ; \
#	./$(EXE) < $(INPUTdy) > $(OUTPUTdy)
#	@echo "Executing the command ./$(EXE) < $(INPUTdz) > $(OUTPUTdz)" ; \
#	./$(EXE) < $(INPUTdz) > $(OUTPUTdz)

$(EXE): $(OBJS) library
	@echo " - linking objects .... " ;\
	$(F77_LD) $(F77_LFLAGS) $(OBJS) $(LIBS) -o $(EXE); \
	echo " - the executable $(EXE) has been created."

BEGIN:
	@echo "Building/checking the executable $(EXE) ... "

library:
	@( $(CD) $(SRC2D); make $(MFLAGS) )
	@( $(CD) $(SRC3D); make $(MFLAGS) )

clean:
	@echo "Removing objects for $(EXE) ... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE), and the output files ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUT_FILES)

# ============================================================================


