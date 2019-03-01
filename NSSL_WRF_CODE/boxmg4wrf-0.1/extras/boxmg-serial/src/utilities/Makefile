# ============================================================================
#  
#  Makefile: Utilities for Symmetric Black Box Multigrid
#
# =======================================================================
# $license_flag$
# ============================================================================

SHELL=/bin/sh

ifndef MAKE
  MAKE=make
endif

ifndef BOXMGdist
  BOXMGdist = ../..
endif

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# --------------------
#  Architecture:
# --------------------

include $(BOXMG_ARCHmake)/ARCH.$(BOXMG_ARCH)

# ---------------------
#  Includes:
# ---------------------

ifndef BOXMG_FINC
  BOXMG_FINC = -I$(BOXMGdist)/include
endif

ifndef BOXMG_CINC
  BOXMG_CINC = -I$(BOXMGdist)/include -I.
endif

# ============================================================================

# ------------------
#  OBJECT MODULES:
# ------------------

O_UTILS = BMG_timer.o 

OBJS = $(O_UTILS)

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(BOXMGlibp)/$(BOXMGlibn)($(obj)))

BEGIN:
	@echo " - building boxmg-utility in $(BOXMGlibn) ... "

clean:
	@echo " - the source subdirectory: utilities"
	@$(RM) -f $(OBJS)







