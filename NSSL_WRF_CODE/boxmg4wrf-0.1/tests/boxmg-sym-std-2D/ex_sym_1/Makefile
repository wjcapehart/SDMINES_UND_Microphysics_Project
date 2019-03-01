# ============================================================================
#
#  Makefile:  Example 3 
# 
#  Author: J. David Moulton
#  Date:   November 13, 1993.
# 
#  21/10/96: updated machine entries
#  01/03/98: updated to use includes and conditionals
#  09/10/98: updated environment variables
#  12/14/98: updated to link to extras library
#  01/22/98  Explicitly define SHELL and set MAKE
#            - fixed specification of library link option 
#  11/08/99: use single ARCH file and global macro file 
#  ****************************************
#
#  This makefile builds and executes an example that uses a direct
#  call to boxmg.  The actual problem is the same as example 1.
#  Assuming you have built and installed the library (see README) 
#  at the command line enter
#
#  make
#
#  the required input is redirected from $(INPUT) and the output is 
#  redirected into the file $(OUTPUT).  This program demonstrates the
#  use of the BOXspace subroutine with variable problem sizes. 
#
#  ----------------------
#  Environment Variables:
#  ----------------------
#
#  BOXMGdist   - location of BOXMG distribution
#  BOXMGlibp   - location of the BOXMG and BOXMG-EXTRAS library
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

include $(ARCHmake)/ARCH.$(MACHINE)

# --------------------
#  Includes:
# --------------------

ifndef FINC
FINC = -I$(BOXMGdist)/include
endif

# ----------------
#  LIBRARIES:
# ----------------

LIBBMG   = boxmg_$(CLEVEL).a
LIBBMGx  = boxmg-extras_$(CLEVEL).a
LIBCHEK  = $(BOXMGlibp)/lib$(LIBBMG) $(BOXMGlibp)/lib$(LIBBMGx)

LIBPATH  = $(LPOPT)$(BOXMGlibp)
LIBLIST  = $(LLOPT)boxmg_$(CLEVEL) $(LLOPT)boxmg-extras_$(CLEVEL)

# ======================================================================
# ----------------------
#   Object Modules:
# ----------------------

MAIN = main.o diff.o initvecs.o vecprod.o
PUTF =  putf.o 

# ----------------------
#   Collect the Macros:
# ----------------------

OBJS = $(MAIN) $(PUTF)
LIBS = $(LIBPATH) $(LIBLIST) $(SYSLIBS)

EXE = ex_sym_1

INPUT = input.txt
OUTPUT = output.txt

# ============================================================================
# -----------------
#   TARGETS:
# -----------------

test: $(EXE)
	@echo "Executing the command ./$(EXE) > $(OUTPUT)" ; \
	./$(EXE) < $(INPUT) > $(OUTPUT)

$(EXE): $(OBJS) $(LIBCHEK)
	@echo "Linking objects .... " ;\
	$(F77_LD) $(F77_LFLAGS) $(OBJS) $(LIBS) -o $(EXE); \
	echo "The executable $(EXE) has been created."

clean:
	@echo "Removing objects .... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE) and $(OUTPUT) ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUTPUT)

# ============================================================================








