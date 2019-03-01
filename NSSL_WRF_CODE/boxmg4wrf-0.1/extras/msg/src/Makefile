# ============================================================================
#
#  Makefile:  MSG Extras 
#
#  Author:  J. David Moulton
#  Date:    January 13, 2000.
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

BOXMG_MPI=yes

# -------------------------
#  Standard Global Macros:
# -------------------------

include $(BOXMGdist)/make/global

# --------------------
#  Architecture:
# --------------------

include $(BOXMG_ARCHmake)/ARCH.$(BOXMG_ARCH)

# ---------------------
#  MSG specific
# ---------------------

include $(EXTRASdist)/msg/ARCH/ARCH.$(BOXMG_ARCH)

ifeq ($(BOXMG_FC_UNDERSCORES),double)
  CFLAGS += -DG77
else
ifeq ($(BOXMG_FC_UNDERSCORES),single)
  CFLAGS += -DNO_SECOND_UNDERSCORE
else
ifeq ($(BOXMG_FC_UNDERSCORES),none)
  CFLAGS += -DNO_UNDERSCORE
endif
endif
endif

BOXMG_FINC = -I../include
BOXMG_CINC = -I../include

# ============================================================================

OBJS = mpi_msg.o tp_setup.o wc_timer.o

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(EXTRASlibp)/$(EXTRASlibn)($(obj)))

BEGIN:
	@echo " - building/checking MSG objects in $(EXTRASlibn) ... "

clean:
	@echo " - the source subdirectory: msg/src"
	@$(RM) -f $(OBJS)
