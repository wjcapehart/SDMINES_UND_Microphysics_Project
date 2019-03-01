# ============================================================================
#  
#  Makefile: Two-Dimensional Symmetric Black Box Multigrid
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

# ============================================================================

# ------------------
#  OBJECT MODULES:
# ------------------

O_BMG   = BMG2_SER_SymStd_SOLVE_boxmg.o        \
	  BMG2_SER_SymStd_SOLVE_cg.o           \
          BMG2_SER_SymStd_ncycle.o             \
          BMG2_SER_SymStd_updown.o             \
          BMG2_SER_SymStd_GET_pointers.o       \
	  BMG2_SER_SymStd_SETUP_parms.o        \
	  BMG2_SER_SymStd_SETUP_PtrWork.o      \
	  BMG2_SER_SymStd_SETUP_PtrGrid.o      \
	  BMG2_SER_SymStd_SETUP_space.o        \
	  BMG2_SER_SymStd_SETUP_recip.o        \
          BMG2_SER_SymStd_SETUP_parts.o        \
          BMG2_SER_SymStd_SETUP_interp_OI.o    \
          BMG2_SER_SymStd_SETUP_ITLI_ex.o      \
          BMG2_SER_SymStd_SETUP_ITLI_bl.o      \
	  BMG2_SER_SymStd_SETUP_cg_LU.o        \
          BMG2_SER_SymStd_SETUP_relax.o        \
	  BMG2_SER_SymStd_SETUP_lines_x.o      \
          BMG2_SER_SymStd_SETUP_lines_y.o      \
          BMG2_SER_SymStd_restrict.o           \
          BMG2_SER_SymStd_interp_add.o         \
          BMG2_SER_SymStd_relax.o              \
          BMG2_SER_SymStd_relax_GS.o           \
          BMG2_SER_SymStd_relax_lines_x.o      \
          BMG2_SER_SymStd_relax_lines_y.o      \
          BMG2_SER_SymStd_residual.o           \
          BMG2_SER_SymStd_ErrTrap.o

# No time to support periodic stuff (2000/07/24 - JDM)
#
O_BMGp  = BMG2_SER_PerSymStd_SETUP_cofp.o      \
          BMG2_SER_PerSymStd_SETUP_cg_LU.o     \
          BMG2_SER_PerSymStd_interp_add.o      \
          BMG2_SER_PerSymStd_residual.o        \
          BMG2_SER_PerSymStd_relax.o           \
          BMG2_SER_PerSymStd_SOLVE_cg.o

O_PCG   = BMG2_SER_SymStd_SOLVE_pcg.o          \
	  BMG2_SER_SymStd_PRECON_diag.o        \
	  BMG2_SER_SymStd_PRECON_boxmg.o


O_UTILS = BMG2_SER_SymStd_UTILS_rV_zero.o      \
	  BMG2_SER_SymStd_UTILS_norm_l2.o      \
	  BMG2_SER_SymStd_UTILS_dot_l2.o       \
	  BMG2_SER_SymStd_UTILS_matvec.o       \
	  BMG2_SER_SymStd_UTILS_daxpy.o        \
	  BMG2_SER_SymStd_UTILS_dscal.o        \
	  BMG2_SER_SymStd_UTILS_dcopy.o        \
	  BMG2_SER_SymStd_UTILS_dxpby.o        \
          BMG2_SER_SymStd_UTILS_reswt.o

O_DUMP  = BMG2_SER_SymStd_DUMP_parms.o         \
          BMG2_SER_SymStd_DUMP_stencil.o       \
          BMG2_SER_SymStd_DUMP_restrict.o      \
	  BMG2_SER_SymStd_DUMP_vector.o


OBJS = $(O_BMG) $(O_BMGp) $(O_PCG) $(O_UTILS) $(O_DUMP)

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(BOXMGlibp)/$(BOXMGlibn)($(obj)))

BEGIN:
	@echo " - building boxmg-sym-std-2D objects in $(BOXMGlibn) ... "

clean:
	@echo " - the source subdirectory: boxmg-sym-std-2D"
	@$(RM) -f $(OBJS)







