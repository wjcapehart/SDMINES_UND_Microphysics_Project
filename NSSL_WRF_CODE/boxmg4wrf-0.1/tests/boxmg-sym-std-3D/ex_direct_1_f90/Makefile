# ============================================================================
#
#  Makefile:  Test 2 - Fortran90/95 setup and call to BoxMG solver
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes an example that uses a Fortran90/95
#  driver and dynamic memory allocation to call to the BoxMG solver.
#  The problem is provides a mechanism to define regions of a domain with 
#  different coefficients, thereby allowing testing with discontinuous 
#  coefficents.  A similar test is provided for three dimensions. If you 
#  want to build this test, just enter
#
#  make
#
#  To run the test you will need to know how to start parallel jobs
#  on your machine, then redirect input from files in the input
#  directory, e.g.,
#
#  mpirun -np 8 < input/FileList-2x2x2-small ./ex_direct_1_f90
#
#  will run using the list of input files in "FileList-2x2x2-small".
#  To change the behavior, size of the problem etc., simply edit
#  the contents of these input files, or create a new set of input
#  files using the same format.
#
#  ** NOTE:
#
#  One important difference between this example and ex_direct_1
#  is that the source is set to zero here, and a random initial guess
#  is used. This facilitates accurate estimates of the asymptotic 
#  convergence rate and more accurate timing estimates (i.e., with
#  larger number of iterations).
#
#  ----------------------
#  Environment Variables:
#  ----------------------
#
#  BOXMGdist   - location of BoxMG distribution
#  BOXMGlibp   - location of the BoxMG and BoxMG-EXTRAS library
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

# --------------------
#  Libraries:
# --------------------

LIBBMG   = boxmg_$(BOXMG_CLEVEL)
LIBBMGx  = boxmg-extras_$(BOXMG_CLEVEL)

BMGlink  = $(LPOPT)$(BOXMGlibp) $(LLOPT)$(LIBBMG)
BMGxlink = $(LPOPT)$(EXTRASlibp) $(LLOPT)$(LIBBMGx)

# ----------------------------------------------------------------------
#  OBJECT MODULES:
# ----------------------

MAIN = main.o               \
       EX_SETUP_BMG_parms.o \
       EX_SETUP_PDE_parms.o 

PUTF = da.o dr.o dz.o putf.o rhs.o

# ----------------------
#   Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
#LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS)  -L/opt/intel/Compiler/11.0/054/Frameworks/mkl/lib/universal -lmkl_intel_lp64  -lmkl_core
#LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS) -llapack -lblas
LIBS = $(BMGlink) $(MATHLIBS) $(BMGxlink) $(SYSLIBS)

EXE  = ex_direct_1_f90

OUTPUT = output/Stencil* output/Vector* output/out*

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

test: BEGIN $(EXE)

$(EXE): $(OBJS) library
	@echo " - linking objects .... " ;\
	$(F90_LD) $(F90_LFLAGS) $(OBJS) $(LIBS) -o $(EXE); \
	echo " - the executable $(EXE) has been created."

BEGIN:
	@echo "Building/checking the executable $(EXE) ... "

library:
	@( $(CD) $(BOXMGdist); make $(MFLAGS) library )

clean:
	@echo "Removing objects for $(EXE) ... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE), and the output files ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUTPUT)


# ============================================================================


