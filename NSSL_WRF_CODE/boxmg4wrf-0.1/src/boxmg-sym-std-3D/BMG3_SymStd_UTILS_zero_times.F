      subroutine BMG3_SymStd_UTILS_zero_times (
     &                BMG_rPARMS )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Zero all of the internal timers.
C     
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'
#include    "BMG_constants.h"
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      REAL*RKIND  BMG_rPARMS(NBMG_rPARMS)

C ----------------------------
C     Local Declarations
C
 
C =========================================================================

      BMG_rPARMS(id_BMG3_TIME_SETUP_FINE_STENCIL) = rZERO
      BMG_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI)      = rZERO
      BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI)    = rZERO
      BMG_rPARMS(id_BMG3_TIME_SETUP_RELAX)        = rZERO
      BMG_rPARMS(id_BMG3_TIME_SETUP_CG_LU)        = rZERO
      BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)           = rZERO
      BMG_rPARMS(id_BMG3_TIME_relax)              = rZERO
      BMG_rPARMS(id_BMG3_TIME_restrict)           = rZERO
      BMG_rPARMS(id_BMG3_TIME_interp_add)         = rZERO
      BMG_rPARMS(id_BMG3_TIME_SETUP_MSG)          = rZERO
      BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID)     = rZERO
      BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS)        = rZERO
      BMG_rPARMS(id_BMG3_TIME_SOLVE_TOTAL)        = rZERO
      BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL)        = rZERO
      BMG_rPARMS(id_BMG3_TIME_PCG_TOTAL)          = rZERO
      BMG_rPARMS(id_BMG3_TIME_PCG_PRECON)         = rZERO

C =========================================================================

      RETURN
      END


