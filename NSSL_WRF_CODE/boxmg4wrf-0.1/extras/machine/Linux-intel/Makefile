# ============================================================================
#
#  Makefile:  Linux-Intel Machine Extras
#
# =======================================================================
# $license_flag$
# ============================================================================
#
#  Please refer to the README in the extras/machine directory for the 
#  license information relating to i1mach, and d1mach, which are from
#  the SLATEC library.
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
  EXTRASdist = $(BOXMGdist)
endif

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# --------------------
#  Architecture:
# --------------------

include $(BOXMG_ARCHmake)/ARCH.$(BOXMG_ARCH)

# ============================================================================

OBJS = i1mach.o d1mach.o

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(EXTRASlibp)/$(EXTRASlibn)($(obj)))

BEGIN:
	@echo " - building machine objects in $(EXTRASlibn) ... "

clean:
	@echo " - the source subdirectory: machine/$(BOXMG_ARCH)"
	@$(RM) -f $(OBJS)


