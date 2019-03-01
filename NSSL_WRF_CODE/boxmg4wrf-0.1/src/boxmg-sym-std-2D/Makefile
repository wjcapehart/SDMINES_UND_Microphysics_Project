# ============================================================================
#  
#    Makefile: MPI-based Two-Dimensional Black Box Multigrid
#
# ===========================================================================

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

# ---------------------
#  Includes:
# ---------------------

ifndef BOXMG_FINC
  BOXMG_FINC =  -I$(BOXMGdist)/include
  BOXMG_FINC += -I$(EXTRASdist)/msg/include
  BOXMG_FINC += -I$(EXTRASdist)/boxmg-serial/include
endif

# ============================================================================

# ------------------
#  OBJECT MODULES:
# ------------------

O_BMG   = BMG2_SymStd_SOLVE_boxmg.o          \
	  BMG2_SymStd_GET_pointers.o         \
          BMG2_SymStd_ncycle.o               \
          BMG2_SymStd_updown.o               \
	  BMG2_SymStd_SETUP_parms.o          \
	  BMG2_SymStd_SETUP_PtrWork.o        \
	  BMG2_SymStd_SETUP_PtrGrid.o        \
	  BMG2_SymStd_SETUP_space.o          \
	  BMG2_SymStd_SETUP_nog.o            \
          BMG2_SymStd_SETUP_parts.o          \
          BMG2_SymStd_SETUP_interp_OI.o      \
          BMG2_SymStd_SETUP_ITLI_ex.o        \
          BMG2_SymStd_SETUP_ITLI_bl.o        \
          BMG2_SymStd_SETUP_relax.o          \
	  BMG2_SymStd_SETUP_recip.o          \
	  BMG2_SymStd_SETUP_lines_x.o        \
          BMG2_SymStd_SETUP_lines_y.o        \
          BMG2_SymStd_restrict.o             \
          BMG2_SymStd_interp_add.o           \
          BMG2_SymStd_relax.o                \
          BMG2_SymStd_relax_GS.o             \
          BMG2_SymStd_relax_lines_x.o        \
          BMG2_SymStd_relax_lines_y.o        \
	  BMG2_SymStd_LineSolve_A.o          \
	  BMG2_SymStd_LineSolve_B.o          \
	  BMG2_SymStd_LineSolve_C.o          \
          BMG2_SymStd_residual.o             

O_CGS   = BMG2_SymStd_SETUP_cg_LU.o          \
	  BMG2_SymStd_SETUP_cg_parms.o       \
	  BMG2_SymStd_SETUP_cg_boxmg.o       \
          BMG2_SymStd_SOLVE_cg_LU.o          \
	  BMG2_SymStd_SOLVE_cg_boxmg.o       \
	  BMG2_SymStd_COPY_cg_WS_SO.o        \
	  BMG2_SymStd_COPY_cg_WS_RHS.o       \
	  BMG2_SymStd_COPY_cg_rV_G_L.o

O_PCG   = BMG2_SymStd_SOLVE_pcg.o            \
	  BMG2_SymStd_PRECON_diag.o          \
	  BMG2_SymStd_PRECON_boxmg.o         

O_MSG   = BMG2_SymStd_SETUP_ProcGrid_file.o  \
	  BMG2_SymStd_SETUP_fine_stencil.o   \
          BMG2_SymStd_SETUP_MSG.o            \
          BMG2_SymStd_SETUP_PtrMSG.o         \
          BMG2_SymStd_SETUP_PtrMSGSO.o       \
          BMG2_SymStd_SETUP_MSGGridSO.o      \
          BMG2_SymStd_SETUP_MSGGrid.o        \
          BMG2_SymStd_SETUP_PtrLS.o          \
          BMG2_SymStd_SETUP_LS.o             \
	  BMG2_SymStd_SETUP_LSGrid.o         

O_UTILS = BMG2_SymStd_UTILS_matvec.o         \
	  BMG2_SymStd_UTILS_dot_l2.o         \
	  BMG2_SymStd_UTILS_norm_l2.o        \
	  BMG2_SymStd_UTILS_daxpy.o          \
	  BMG2_SymStd_UTILS_dxpby.o          \
	  BMG2_SymStd_UTILS_dscal.o          \
	  BMG2_SymStd_UTILS_dcopy.o          \
	  BMG2_SymStd_UTILS_rV_zero.o        \
	  BMG2_SymStd_UTILS_update_ghosts.o  \
	  BMG2_SymStd_UTILS_zero_times.o     \
	  BMG2_SymStd_UTILS_bbody5.o         \
	  BMG2_SymStd_UTILS_bbody9.o         \
	  BMG2_SymStd_UTILS_reswt.o          \
	  BMG2_SymStd_OUTPUT_times.o         \
	  BMG2_SymStd_ErrTrap.o

O_DUMP  = BMG2_SymStd_DUMP_vector.o          \
	  BMG2_SymStd_DUMP_stencil.o         \
	  BMG2_SymStd_DUMP_restrict.o        \
	  BMG2_SymStd_DUMP_parms.o           \
	  BMG2_SymStd_DUMP_interp.o

OBJS = $(O_BMG) $(O_CGS) $(O_PCG) $(O_MSG) $(O_UTILS) $(O_DUMP)

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



