# ============================================================================
#
#  Example:  Conjugate Gradient Solution of 2D Laplace's Equation
#
#  Author:  J. David Moulton
#  Date:    January 13, 2000.
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
MAKE=make
endif

ifndef EXTRASdist
EXTRASdist = ../../..
endif

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(EXTRASdist)/make/global

# --------------------
#  Architecture:
# --------------------

include $(ARCHmake)/ARCH.$(MACHINE)

# ---------------------
#  MSG specific:
# ---------------------

include $(EXTRASdist)/msg/ARCH/ARCH.$(MACHINE)

FINC = -I../../include
CINC = -I../../include

# --------------------
#  LIBRARIES:
# --------------------

LIBBMGx  = boxmg-extras-mpi_$(CLEVEL)
LIBCHEK  = $(EXTRASdist)/lib/lib$(LIBBMGx).a

LIBPATH  = $(LPOPT)$(EXTRASdist)/lib
LIBLIST  = $(LLOPT)$(LIBBMGx) 

# ============================================================================

# ----------------------
#   Object Modules:
# ----------------------

MAIN = cg2.o diststencilhalo2D.o

# ----------------------
#   Collect the Macros:
# ----------------------

OBJS = $(MAIN)
LIBS = $(LIBPATH) $(LIBLIST) $(MPI_SYSLIBS) $(SYSLIBS) 

EXE = cg2D

INPUT = input.txt
OUTPUT = output.txt

# ========================================================================= 

# -----------------
#   TARGETS:
# -----------------

test: $(EXE)
#	@echo "Executing the command:  $(MPI_RUN) $(MPI_RFLAGS) $(EXE)" 
#	$(MPI_RUN) $(MPIRUN_FLAGS) $(EXE) 

$(EXE): $(OBJS) cg.h $(LIBCHEK)
	@echo "Linking objects .... " ;\
	$(MPI_LD) $(MPI_LFLAGS) $(MPI_FFLAGS) $(OBJS) $(LIBS) -o $(EXE);\
	echo "The executable $(EXE) has been created."

clean:
	@echo "Removing objects .... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE) and $(OUTPUT) ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUTPUT)

# ============================================================================


