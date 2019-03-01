# ============================================================================
#
#  Makefile:  Lapack Extras
#
# =======================================================================
# $license_flag$
# ============================================================================
#
#  Please refer to the README in this directory for the open source
#  license information relating to these LAPACK routines.
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
#   Serial LAPACK
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

# --------------------
#  Source PATH
# --------------------

BLAS_SRC=$(EXTRASdist)/lapack/blas
SNGL_SRC=$(EXTRASdist)/lapack/lapack/single
DOUB_SRC=$(EXTRASdist)/lapack/lapack/double
UTIL_SRC=$(EXTRASdist)/lapack/lapack/util

# ============================================================================

BLAS_OBJ = $(BLAS_SRC)/ddot.o  $(BLAS_SRC)/dgemm.o \
           $(BLAS_SRC)/dgemv.o $(BLAS_SRC)/dscal.o \
           $(BLAS_SRC)/dsyr.o  $(BLAS_SRC)/dsyrk.o \
           $(BLAS_SRC)/dtbsv.o $(BLAS_SRC)/dtrsm.o \
           $(BLAS_SRC)/sdot.o  $(BLAS_SRC)/sgemm.o \
           $(BLAS_SRC)/sgemv.o $(BLAS_SRC)/sscal.o \
           $(BLAS_SRC)/ssyr.o  $(BLAS_SRC)/ssyrk.o \
           $(BLAS_SRC)/stbsv.o $(BLAS_SRC)/strsm.o \
           $(BLAS_SRC)/xerbla.o

SNGL_OBJ = $(SNGL_SRC)/spbsv.o  $(SNGL_SRC)/spbtf2.o \
           $(SNGL_SRC)/spbtrf.o $(SNGL_SRC)/spbtrs.o \
           $(SNGL_SRC)/spotf2.o $(SNGL_SRC)/spotrf.o \
           $(SNGL_SRC)/spotrs.o $(SNGL_SRC)/spttrf.o \
           $(SNGL_SRC)/spttrs.o $(SNGL_SRC)/sptts2.o

DOUB_OBJ = $(DOUB_SRC)/dpbsv.o  $(DOUB_SRC)/dpbtf2.o \
           $(DOUB_SRC)/dpbtrf.o $(DOUB_SRC)/dpbtrs.o \
           $(DOUB_SRC)/dpotf2.o $(DOUB_SRC)/dpotrf.o \
           $(DOUB_SRC)/dpotrs.o $(DOUB_SRC)/dpttrf.o \
           $(DOUB_SRC)/dpttrs.o $(DOUB_SRC)/dptts2.o

UTIL_OBJ = $(UTIL_SRC)/ieeeck.o $(UTIL_SRC)/ilaenv.o \
           $(UTIL_SRC)/lsame.o  $(UTIL_SRC)/sisnan.o \
           $(UTIL_SRC)/slaisnan.o

OBJS = $(BLAS_OBJ) $(DOUB_OBJ) $(SNGL_OBJ) $(UTIL_OBJ)

# ========================================================================= 
# -----------------
#   TARGETS:
# -----------------

library: BEGIN $(foreach obj, $(OBJS), $(EXTRASlibp)/$(EXTRASlibn)($(obj)))

BEGIN:
	@echo " - building/checking LAPACK objects in $(EXTRASlibn) ... "

clean:
	@echo " - the source subdirectory: lapack"
	@$(RM) -f $(OBJS)

