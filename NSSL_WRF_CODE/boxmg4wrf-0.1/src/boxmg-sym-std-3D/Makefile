# ============================================================================
#  
#    Makefile:  MPI-based Three-Dimensional Black Box Multigrid 
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

O_BMG   = BMG3_SymStd_SOLVE_boxmg.o         \
          BMG3_SymStd_GET_pointers.o        \
          BMG3_SymStd_ncycle.o              \
          BMG3_SymStd_updown.o              \
	  BMG3_SymStd_SETUP_parms.o         \
	  BMG3_SymStd_SETUP_PtrWork.o       \
	  BMG3_SymStd_SETUP_PtrGrid.o       \
	  BMG3_SymStd_SETUP_space.o         \
	  BMG3_SymStd_SETUP_nog.o           \
	  BMG3_SymStd_SETUP_parts.o         \
          BMG3_SymStd_SETUP_interp_OI.o     \
	  BMG3_SymStd_SETUP_ITLI07_ex.o     \
	  BMG3_SymStd_SETUP_ITLI27_ex.o     \
	  BMG3_SymStd_SETUP_ITLI_bl.o       \
	  BMG3_SymStd_SETUP_ITLI_Izyx.o     \
	  BMG3_SymStd_SETUP_relax.o         \
	  BMG3_SymStd_SETUP_recip.o         \
	  BMG3_SymStd_SETUP_planes_xy.o     \
	  BMG3_SymStd_SETUP_planes_yz.o     \
	  BMG3_SymStd_SETUP_planes_xz.o     \
          BMG3_SymStd_restrict.o            \
          BMG3_SymStd_interp_add.o          \
          BMG3_SymStd_relax.o               \
          BMG3_SymStd_relax_GS.o            \
	  BMG3_SymStd_relax_planes_xy.o     \
          BMG3_SymStd_relax_planes_yz.o     \
	  BMG3_SymStd_relax_planes_xz.o     \
          BMG3_SymStd_residual.o            \
	  BMG3_SymStd_update_yo.o           \
	  BMG3_SymStd_update_zo.o           \
	  BMG3_SymStd_update_sor2.o  

O_CGS   = BMG3_SymStd_SETUP_cg_LU.o         \
	  BMG3_SymStd_SETUP_cg_parms.o      \
	  BMG3_SymStd_SETUP_cg_boxmg.o      \
          BMG3_SymStd_SOLVE_cg_LU.o         \
	  BMG3_SymStd_SOLVE_cg_boxmg.o      \
	  BMG3_SymStd_COPY_cg_WS_SO.o       \
	  BMG3_SymStd_COPY_cg_WS_RHS.o      \
	  BMG3_SymStd_COPY_cg_rV_G_L.o

O_PCG   = BMG3_SymStd_SOLVE_pcg.o	    \
	  BMG3_SymStd_PRECON_diag.o	    \
	  BMG3_SymStd_PRECON_boxmg.o	    \

O_MSG   = BMG3_SymStd_SETUP_ProcGrid_file.o \
	  BMG3_SymStd_SETUP_fine_stencil.o  \
          BMG3_SymStd_SETUP_MSG.o           \
          BMG3_SymStd_SETUP_PtrMSG.o        \
          BMG3_SymStd_SETUP_PtrMSGSO.o      \
          BMG3_SymStd_SETUP_MSGGrid.o       \
	  BMG3_SymStd_SETUP_nMSGGrid.o

O_UTILS = BMG3_SymStd_UTILS_matvec.o        \
	  BMG3_SymStd_UTILS_daxpy.o	    \
	  BMG3_SymStd_UTILS_dcopy.o	    \
	  BMG3_SymStd_UTILS_dot_l2.o        \
	  BMG3_SymStd_UTILS_dscal.o	    \
	  BMG3_SymStd_UTILS_dxpby.o	    \
          BMG3_SymStd_UTILS_rV_zero.o       \
          BMG3_SymStd_UTILS_norm_l2.o       \
	  BMG3_SymStd_UTILS_bbody07.o        \
	  BMG3_SymStd_UTILS_bbody27.o        \
	  BMG3_SymStd_UTILS_reswt.o         \
          BMG3_SymStd_UTILS_update_ghosts.o \
	  BMG3_SymStd_ErrTrap.o             \
          BMG3_SymStd_UTILS_zero_times.o

O_COPY  = BMG3_SymStd_COPY_SO_xy.o          \
	  BMG3_SymStd_COPY_SO_yz.o          \
	  BMG3_SymStd_COPY_SO_xz.o          \
          BMG3_SymStd_COPY_RHS_xy.o         \
          BMG3_SymStd_COPY_RHS_yz.o         \
          BMG3_SymStd_COPY_RHS_xz.o         \
          BMG3_SymStd_COPY_rV_23_xy.o       \
          BMG3_SymStd_COPY_rV_23_yz.o       \
          BMG3_SymStd_COPY_rV_23_xz.o       \
          BMG3_SymStd_COPY_rV_32_xy.o       \
          BMG3_SymStd_COPY_rV_32_yz.o       \
          BMG3_SymStd_COPY_rV_32_xz.o

O_DUMP  = BMG3_SymStd_DUMP_parms.o         \
          BMG3_SymStd_DUMP_stencil.o       \
          BMG3_SymStd_DUMP_restrict.o      \
          BMG3_SymStd_DUMP_vector.o        \
	  BMG3_SymStd_OUTPUT_times.o


OBJS = $(O_BMG) $(O_CGS) $(O_PCG) $(O_MSG) $(O_UTILS) $(O_COPY) $(O_DUMP)

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(BOXMGlibp)/$(BOXMGlibn)($(obj)))

BEGIN:
	@echo " - building boxmg-sym-std-3D objects in $(BOXMGlibn) ... "

clean:
	@echo " - the source subdirectory: boxmg-sym-std-3D"
	@$(RM) -f $(OBJS)




