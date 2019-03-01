# ============================================================================
#
#  Makefile:  Psuedo-Random Number Generator
#
# =======================================================================
# $license_flag$
# ============================================================================
#
#  Please refer to the README in this directory for the license
#  information relating to BMG_rand, which is from the SLATEC library.
#
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

ifndef BOXMGdist
  BOXMGdist = ../..
endif

ifndef EXTRASdist
  EXTRASdist = $(BOXMGdist)/extras
endif

# -------------------------
#   Serial 
# -------------------------

MPI=no

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# --------------------
#  Architecture:
# --------------------

include $(BOXMG_ARCHmake)/ARCH.$(BOXMG_ARCH)

# ============================================================================

OBJS = BMG_rand.o

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(EXTRASlibp)/$(EXTRASlibn)($(obj)))

BEGIN:
	@echo " - building/checking random objects in $(EXTRASlibn) ... "

clean:
	@echo " - the source subdirectory: random"
	@$(RM) -f $(OBJS)


