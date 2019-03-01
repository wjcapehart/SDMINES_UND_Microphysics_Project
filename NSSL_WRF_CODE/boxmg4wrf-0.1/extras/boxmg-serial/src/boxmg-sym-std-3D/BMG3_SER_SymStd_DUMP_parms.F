      SUBROUTINE BMG3_SER_SymStd_DUMP_parms( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_DUMP_parms outputs the parameter arrays.
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
      INCLUDE 'BMG_SER_workspace.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER BMG_iPARMS(NBMG_SER_iPARMS)
      REAL*RKIND  BMG_rPARMS(NBMG_SER_rPARMS)
      LOGICAL BMG_IOFLAG(NBMG_SER_IOFLAG)

C --------------------------
C     Local Declarations:
C

C ==========================================================================

      WRITE(*,*) '****  Parmeters: 3D '
      WRITE(*,*) 
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_POINTERS)    = ',
     &           BMG_iPARMS(id_BMG3_SER_POINTERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_STENCIL)     = ', 
     &           BMG_iPARMS(id_BMG3_SER_STENCIL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_SETUP)       = ',
     &           BMG_iPARMS(id_BMG3_SER_SETUP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_RELAX)       = ',
     &           BMG_iPARMS(id_BMG3_SER_RELAX)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_RELAX_SYM)   = ',
     &           BMG_iPARMS(id_BMG3_SER_RELAX_SYM )
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_NRELAX_DOWN) = ',
     &           BMG_iPARMS(id_BMG3_SER_NRELAX_DOWN)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_NRELAX_UP)   = ',
     &           BMG_iPARMS(id_BMG3_SER_NRELAX_UP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_NRELAX_FG)   = ',
     &           BMG_iPARMS(id_BMG3_SER_NRELAX_FG)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_CYCLE_CLASS) = ',
     &           BMG_iPARMS(id_BMG3_SER_CYCLE_CLASS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_NCYCLE_TYPE) = ',
     &           BMG_iPARMS(id_BMG3_SER_NCYCLE_TYPE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_FMG_NNCYCLE) = ',
     &           BMG_iPARMS(id_BMG3_SER_FMG_NNCYCLE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_MAX_ITERS)   = ', 
     &           BMG_iPARMS(id_BMG3_SER_MAX_ITERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_STOP_TEST)   = ',
     &           BMG_iPARMS(id_BMG3_SER_STOP_TEST)
      WRITE(*,*) 'BMG_rPARMS(id_BMG3_SER_STOP_TOL)    = ',
     &           BMG_rPARMS(id_BMG3_SER_STOP_TOL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SER_CG_MIN_DIM)  = ', 
     &           BMG_iPARMS(id_BMG3_SER_CG_MIN_DIM)


      WRITE(*,*) '****  Parameters: 2D '
      WRITE(*,*) 
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_POINTERS)    = ',
     &           BMG_iPARMS(id_BMG2_SER_POINTERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_STENCIL)     = ', 
     &           BMG_iPARMS(id_BMG2_SER_STENCIL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_SETUP)       = ',
     &           BMG_iPARMS(id_BMG2_SER_SETUP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_RELAX)       = ',
     &           BMG_iPARMS(id_BMG2_SER_RELAX)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_RELAX_SYM)   = ',
     &           BMG_iPARMS(id_BMG2_SER_RELAX_SYM )
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_NRELAX_DOWN) = ',
     &           BMG_iPARMS(id_BMG2_SER_NRELAX_DOWN)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_NRELAX_UP)   = ',
     &           BMG_iPARMS(id_BMG2_SER_NRELAX_UP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_NRELAX_FG)   = ',
     &           BMG_iPARMS(id_BMG2_SER_NRELAX_FG)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS) = ',
     &           BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE) = ',
     &           BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_FMG_NNCYCLE) = ',
     &           BMG_iPARMS(id_BMG2_SER_FMG_NNCYCLE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_MAX_ITERS)   = ', 
     &           BMG_iPARMS(id_BMG2_SER_MAX_ITERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_STOP_TEST)   = ',
     &           BMG_iPARMS(id_BMG2_SER_STOP_TEST)
      WRITE(*,*) 'BMG_rPARMS(id_BMG2_SER_STOP_TOL)    = ',
     &           BMG_rPARMS(id_BMG2_SER_STOP_TOL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SER_CG_MIN_DIM)  = ', 
     &           BMG_iPARMS(id_BMG2_SER_CG_MIN_DIM)

C ==========================================================================

      RETURN
      END
