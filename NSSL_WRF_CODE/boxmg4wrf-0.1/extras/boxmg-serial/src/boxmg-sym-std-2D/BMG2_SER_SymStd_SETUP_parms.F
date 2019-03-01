      SUBROUTINE BMG2_SER_SymStd_SETUP_parms( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAGS
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_SETUP_parms.f sets up default cycle parameters
C     which should work for most problems.
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
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_SER_rPARMS)
      LOGICAL  BMG_IOFLAGS(NBMG_SER_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER i

C ==========================================================================

C -------------------------------------
C    Memory allocation:
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_POINTERS) = BMG_SER_NO_pointers

C -------------------------------------
C    Stencil and BCs:
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_STENCIL)  = BMG_SER_STENCIL_5pt
      BMG_iPARMS(id_BMG2_SER_BC)       = BMG_SER_BCs_definite

C -------------------------------------
C     Setup:
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_SETUP) = BMG_SER_SETUP_Ptrs_Opers

C -------------------------------------
C     Relaxation:
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_RELAX)      = BMG_SER_GS_RB_x_y_lines
      BMG_iPARMS(id_BMG2_SER_RELAX_SYM ) = BMG_SER_RELAX_SYM

      BMG_iPARMS(id_BMG2_SER_NRELAX_DOWN ) = 1
      BMG_iPARMS(id_BMG2_SER_NRELAX_UP )   = 1
      BMG_iPARMS(id_BMG2_SER_NRELAX_FG )   = 0   ! not used at the moment

C -------------------------------------
C     Cycle Class and Type
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS) = BMG_SER_N_CYCLE
      BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE) = BMG_SER_W_CYCLE
      BMG_iPARMS(id_BMG2_SER_FMG_NNCYCLE) = 1

C -------------------------------------
C     Stopping Criteria
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_MAX_ITERS) = 10
      BMG_iPARMS(id_BMG2_SER_STOP_TEST) = BMG_SER_STOP_REL_RES_L2
      BMG_rPARMS(id_BMG2_SER_STOP_TOL)  = 1D-6

C -------------------------------------
C     Coarsening limit
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_CG_MIN_DIM) = 3

C -------------------------------------
C     Coarse-Grid Operator TYPE
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_CG_TYPE) = BMG_SER_CG_ITLI

C -------------------------------------
C     Coarse-Grid Operator CONSTRUCTION
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_CG_CONSTRUCT) = BMG_SER_CG_CONS_explicit

C -------------------------------------
C     Minimum number of grids
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_MIN_NOG) = 2

C -------------------------------------
C     Start with error code = 0
C -------------------------------------

      BMG_iPARMS(id_BMG2_SER_Err_Code) = 0
      BMG_iPARMS(id_BMG2_SER_Ext_Err_Code) = 0

C -------------------------------------
C     I/O Parameters
C -------------------------------------

      DO i=1, NBMG_SER_IOFLAG
         BMG_IOFLAGS(i)=.FALSE.
      ENDDO

C ==========================================================================

      RETURN
      END
