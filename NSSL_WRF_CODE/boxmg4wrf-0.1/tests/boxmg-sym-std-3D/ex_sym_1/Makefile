# ============================================================================
#
#  Makefile:  3D Black Box Example
# 
#  Author: J. David Moulton
#  Date:   May 3, 1999.
# 
#
#  ****************************************
#
#  FILE LIST:
#  ~~~~~~~~~~ 
#
#    main.f       - driver
#    putf.f       - discretization of the PDE
#    da.f         - discretization
#    dr.f         - discretization
#    dz.f         - discretization
#    rhs.f        - discretization of the source
#
#  This makefile builds and executes the "original" example.
#  It assumes that you have built and installed the library 
#  (see README) at the command line enter
#
#  make 
#
#  the output is redirected into the file $(OUTPUT).
#
#  ----------------------
#  Environment Variables:
#  ----------------------
#
#  BOXMGdist   - location of BOXMG distribution
#  BOXMGlibp   - location of the BOXMG and BOXMG-EXTRAS library
#
#            
#  OS          - Operating system (e.g., Linux, SunOS )
#  VERSION     - OS version (e.g., 2.0.x, 5.5.x)
#  CPU         - CPU (e.g., i686, sun4u)
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

ifndef ARCHmake
ARCHmake = $(BOXMGdist)/ARCH
endif

# --------------------
#  Includes:
# --------------------

ifndef FINC
FINC = -I$(BOXMGdist)/include -I./include
endif

# --------------------
#  Path Definitions:
# --------------------

# ----------------
#  LIBRARIES:
# ----------------

LIBBMG   = boxmg_$(CLEVEL).a
LIBBMGx  = boxmg-extras_$(CLEVEL).a
LIBCHEK  = $(BOXMGlibp)/lib$(LIBBMG) $(BOXMGlibp)/lib$(LIBBMGx)

LIBPATH  = $(LPOPT)$(BOXMGlibp)
LIBLIST  = $(LLOPT)boxmg_$(CLEVEL) $(LLOPT)boxmg-extras_$(CLEVEL)

# --------------------
#  Compiler, Linker:
# --------------------

include $(ARCHmake)/ARCH.$(MACHINE)

# ----------------------------------------------------------------------
#  OBJECT MODULES:
# ----------------------

MAIN = main.o initvecs.o vecprod.o
PUTF = da.o dr.o dz.o putf.o rhs.o

# ----------------------
#   Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
LIBS = $(LIBPATH) $(LIBLIST) $(SYSLIBS)

EXE  = ex_sym_1

INPUTa  = id.9

OUTPUTa  = output.09.txt

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

test: $(EXE)
	@echo "Executing the command ./$(EXE) < $(INPUTa) > $(OUTPUTa)" ; \
	./$(EXE) < $(INPUTa) > $(OUTPUTa)

$(EXE): $(OBJS) $(LIBCHEK)
	@echo "Linking objects .... " ;\
	$(F77_LD) $(F77_LFLAGS) $(OBJS) $(LIBS) -o $(EXE); \
	echo "The executable $(EXE) has been created."

clean:
	@echo "Removing objects .... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE), and the output files ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUTPUTa) $(OUTPUTb) \
                 $(OUTPUTcx) $(OUTPUTcy) $(OUTPUTcz)

# ============================================================================


