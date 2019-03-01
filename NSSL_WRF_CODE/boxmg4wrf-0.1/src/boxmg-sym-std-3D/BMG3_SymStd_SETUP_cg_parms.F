      SUBROUTINE BMG3_SymStd_SETUP_cg_parms(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_SER_iPARMS, BMG_SER_rPARMS, BMG_SER_IOFLAG,
     &                NGx, NGy, NGz, NOG
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG3_SymStd_SETUP_cg_parms sets up the cycling parameters
C   for the call to the serial BoxMG on the coarsest parallel grid.
C   It first sets the defaults from the serial, and then overrides
C   with current values. 
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

C ---------------------------
C     Includes
C

      INCLUDE 'mpif.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

      INCLUDE 'BMG_SER_parameters.h'
      INCLUDE 'BMG_SER_workspace-partial.h'

C ---------------------------
C     Argument Declarations:
C
      !
      ! BoxMG MPI Global Coarse-Grid Dimensions with ghosts.
      !
      INTEGER  NGx, NGy, NGz
      INTEGER  NOG

      !
      ! BoxMG MPI Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! BoxMG Serial Parameters
      !
      INTEGER  BMG_SER_iPARMS(NBMG_SER_iPARMS)
      REAL*RKIND   BMG_SER_rPARMS(NBMG_SER_rPARMS)
      LOGICAL  BMG_SER_IOFLAG(NBMG_SER_IOFLAG)

C ==========================================================================

      !
      ! Setup the default parameters
      !
      CALL BMG3_SER_SymStd_SETUP_parms( 
     &          BMG_SER_iPARMS, BMG_SER_rPARMS, BMG_SER_IOFLAG
     &          )

      !
      ! Match 3D parameters
      !
      BMG_SER_iPARMS(id_BMG3_SER_RELAX)
     &    = BMG_iPARMS(id_BMG3_RELAX)
      BMG_SER_iPARMS(id_BMG3_SER_RELAX_SYM )
     &    = BMG_iPARMS(id_BMG3_RELAX_SYM )

      BMG_SER_iPARMS(id_BMG3_SER_NRELAX_DOWN )
     &    = BMG_iPARMS(id_BMG3_NRELAX_DOWN ) 
      BMG_SER_iPARMS(id_BMG3_SER_NRELAX_UP )
     &    = BMG_iPARMS(id_BMG3_NRELAX_UP )

      BMG_SER_iPARMS(id_BMG3_SER_CYCLE_CLASS)
     &    = BMG_iPARMS(id_BMG3_CYCLE_CLASS)
      BMG_SER_iPARMS(id_BMG3_SER_NCYCLE_TYPE)
     &    = BMG_iPARMS(id_BMG3_NCYCLE_TYPE)
      BMG_SER_iPARMS(id_BMG3_SER_FMG_NNCYCLE)
     &    = BMG_iPARMS(id_BMG3_FMG_NNCYCLE)

      ! Coarse Grid Operator Construction (explicit, block)
      BMG_SER_iPARMS(id_BMG3_SER_CG_CONSTRUCT)
     &    = BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
      
      ! Coarse-grid operator type (ITLI, ITLI_IzIyIx)
      BMG_SER_iPARMS(id_BMG3_SER_CG_TYPE) 
     &    = BMG_iPARMS(id_BMG3_CG_TYPE)

      !
      ! Match 2D parameters
      !
      BMG_SER_iPARMS(id_BMG2_SER_RELAX)
     &    = BMG_iPARMS(id_BMG2_RELAX)
      BMG_SER_iPARMS(id_BMG2_SER_RELAX_SYM )
     &    = BMG_iPARMS(id_BMG2_RELAX_SYM )

      BMG_SER_iPARMS(id_BMG2_SER_NRELAX_DOWN )
     &    = BMG_iPARMS(id_BMG2_NRELAX_DOWN ) 
      BMG_SER_iPARMS(id_BMG2_SER_NRELAX_UP )
     &    = BMG_iPARMS(id_BMG2_NRELAX_UP )

      BMG_SER_iPARMS(id_BMG2_SER_CYCLE_CLASS)
     &    = BMG_iPARMS(id_BMG2_CYCLE_CLASS)
      BMG_SER_iPARMS(id_BMG2_SER_NCYCLE_TYPE)
     &    = BMG_iPARMS(id_BMG2_NCYCLE_TYPE)
      BMG_SER_iPARMS(id_BMG2_SER_FMG_NNCYCLE)
     &    = BMG_iPARMS(id_BMG2_FMG_NNCYCLE)

      !
      ! Just the bottom of the global cycle ...
      !
      IF ( NOG.GT.1 ) THEN
         BMG_SER_iPARMS(id_BMG3_SER_MAX_ITERS) = 1
      ELSEIF ( NOG.EQ.1 ) THEN
         BMG_SER_iPARMS(id_BMG3_SER_MAX_ITERS) 
     &       = BMG_iPARMS(id_BMG3_MAX_ITERS) 
         BMG_SER_rPARMS(id_BMG3_SER_STOP_TOL) 
     &       = BMG_rPARMS(id_BMG3_STOP_TOL)
      ENDIF

         
      BMG_SER_iPARMS(id_BMG3_SER_MIN_NOG)   = 1

      !
      ! Coarse-grid stencil is always 27 point
      !
      BMG_SER_iPARMS(id_BMG3_SER_STENCIL) = BMG_SER_STENCIL_27pt

      !
      ! To test just one grid in the serial BoxMG
      !
C     BMG_SER_iPARMS(id_BMG3_SER_CG_MIN_DIM) = NGx-2

      BMG_SER_IOFLAG(iBMG3_SER_BUG_RES_RELAX)
     &    = BMG_IOFLAG(iBMG3_BUG_RES_RELAX)
      
      IF ( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) ) THEN
         BMG_SER_IOFLAG(iBMG3_SER_BUG_RES_CG_SOLVE) = .TRUE.
      ELSE
         BMG_SER_IOFLAG(iBMG3_SER_BUG_RES_CG_SOLVE) = .FALSE.
      ENDIF

C ==========================================================================

      RETURN
      END

