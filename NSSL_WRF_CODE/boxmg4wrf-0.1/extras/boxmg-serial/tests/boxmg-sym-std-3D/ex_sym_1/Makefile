# ============================================================================
#
#  Makefile:  Test V-Cycle Symmetry 
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes an example that tests the ability
#  of the BoxMG solver to perform V-cycles in a symmetric fashion.
#  This is necessary for using BoxMG as a preconditioner for a symmetric
#  Krylov method such as Conjugate Gradients (CG).  The actual problem 
#  is the same as test 1 (ex_bbmg_1).  This will check that the libraries 
#  have the necessary objects, so just enter
#
#  make
#
#  the required input is redirected from $(INPUT).
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

MAIN = main.o initvecs.o vecprod.o
PUTF = da.o dr.o dz.o putf.o rhs.o

# ----------------------
#   Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS)

EXE  = ex_sym_1

INPUTa  = id.9

OUTPUTa  = output.09.txt

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

test: BEGIN $(EXE)
	@echo "Executing the command ./$(EXE) < $(INPUTa) > $(OUTPUTa)" ; \
	./$(EXE) < $(INPUTa) > $(OUTPUTa)

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
	$(RM) -f $(OBJS) $(EXE) $(OUTPUTa) $(OUTPUTb) \
                 $(OUTPUTcx) $(OUTPUTcy) $(OUTPUTcz)

# ============================================================================


