# ============================================================================
#
#  Makefile:  Test 1 - bbmg driver (Example 1)
# 
# =======================================================================
# $license_flag$
# ============================================================================
#
#  This makefile builds and executes the fist example in Victor Bandy's 
#  bbmg driver/boxmg release from mgnet.  This will check that the 
#  libraries have the necessary objects, so just enter
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

LIBBMG   = boxmg_$(CLEVEL)
LIBBMGx  = boxmg-extras_$(CLEVEL)

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

EXE = ex_bbmg_1

OUTPUT = output.txt

# ============================================================================
# -----------------
#   Targets:
# -----------------

test: BEGIN $(EXE)
	@echo "Executing the command ./$(EXE) > $(OUTPUT)" ; \
	./$(EXE) > $(OUTPUT)

$(EXE): $(OBJS) library
	@echo " - linking objects .... "; \
	$(F77_LD) $(F77_LFLAGS) $(OBJS) $(LIBS) -o $(EXE); \
	echo " - the executable $(EXE) has been created."

BEGIN:
	@echo "Building/checking the executable $(EXE) ... "

library:
	@( $(CD) $(SRC2D); make $(MFLAGS) )
	@( $(CD) $(EXTRASdist); make $(MFLAGS) BBMG=yes ) 

clean:
	@echo "Removing objects for $(EXE) ... " ;\
	$(RM) -f $(OBJS)

distclean:
	@echo "Removing objects, $(EXE) and $(OUTPUT) ... ";\
	$(RM) -f $(OBJS) $(EXE) $(OUTPUT)

# ============================================================================



