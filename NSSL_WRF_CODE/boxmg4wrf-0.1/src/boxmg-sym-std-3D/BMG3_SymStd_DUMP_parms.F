      SUBROUTINE BMG3_SymStd_DUMP_parms( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG3_SymStd_DUMP_parms outputs the parameter arrays.
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
C   OUTPUT:
C  --------------------
C
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
#include    "BMG_constants.h"
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND  BMG_rPARMS(NBMG_rPARMS)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C --------------------------
C     Local Declarations:
C

C ==========================================================================


      WRITE(*,*) 
      WRITE(*,*) '****  Parmeters: 3D '
      WRITE(*,*) 
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_POINTERS)    = ',
     &           BMG_iPARMS(id_BMG3_POINTERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_STENCIL)     = ', 
     &           BMG_iPARMS(id_BMG3_STENCIL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_SETUP)       = ',
     &           BMG_iPARMS(id_BMG3_SETUP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_RELAX)       = ',
     &           BMG_iPARMS(id_BMG3_RELAX)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_RELAX_SYM)   = ',
     &           BMG_iPARMS(id_BMG3_RELAX_SYM )
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_NRELAX_DOWN) = ',
     &           BMG_iPARMS(id_BMG3_NRELAX_DOWN)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_NRELAX_UP)   = ',
     &           BMG_iPARMS(id_BMG3_NRELAX_UP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_NRELAX_FG)   = ',
     &           BMG_iPARMS(id_BMG3_NRELAX_FG)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_CYCLE_CLASS) = ',
     &           BMG_iPARMS(id_BMG3_CYCLE_CLASS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = ',
     &           BMG_iPARMS(id_BMG3_NCYCLE_TYPE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_FMG_NNCYCLE) = ',
     &           BMG_iPARMS(id_BMG3_FMG_NNCYCLE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_MAX_ITERS)   = ', 
     &           BMG_iPARMS(id_BMG3_MAX_ITERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_STOP_TEST)   = ',
     &           BMG_iPARMS(id_BMG3_STOP_TEST)
      WRITE(*,*) 'BMG_rPARMS(id_BMG3_STOP_TOL)    = ',
     &           BMG_rPARMS(id_BMG3_STOP_TOL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_CG_MIN_DIM)  = ', 
     &           BMG_iPARMS(id_BMG3_CG_MIN_DIM)
      WRITE(*,*) 'BMG_iPARMS(id_BMG3_CG_SOLVER)  = ', 
     &           BMG_iPARMS(id_BMG3_CG_SOLVER)

      
      WRITE(*,*)
      WRITE(*,*) '****  Parameters: 2D '
      WRITE(*,*) 
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_POINTERS)    = ',
     &           BMG_iPARMS(id_BMG2_POINTERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_STENCIL)     = ', 
     &           BMG_iPARMS(id_BMG2_STENCIL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SETUP)       = ',
     &           BMG_iPARMS(id_BMG2_SETUP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_RELAX)       = ',
     &           BMG_iPARMS(id_BMG2_RELAX)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_RELAX_SYM)   = ',
     &           BMG_iPARMS(id_BMG2_RELAX_SYM )
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_NRELAX_DOWN) = ',
     &           BMG_iPARMS(id_BMG2_NRELAX_DOWN)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_NRELAX_UP)   = ',
     &           BMG_iPARMS(id_BMG2_NRELAX_UP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_NRELAX_FG)   = ',
     &           BMG_iPARMS(id_BMG2_NRELAX_FG)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_CYCLE_CLASS) = ',
     &           BMG_iPARMS(id_BMG2_CYCLE_CLASS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = ',
     &           BMG_iPARMS(id_BMG2_NCYCLE_TYPE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_FMG_NNCYCLE) = ',
     &           BMG_iPARMS(id_BMG2_FMG_NNCYCLE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_MAX_ITERS)   = ', 
     &           BMG_iPARMS(id_BMG2_MAX_ITERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_STOP_TEST)   = ',
     &           BMG_iPARMS(id_BMG2_STOP_TEST)
      WRITE(*,*) 'BMG_rPARMS(id_BMG2_STOP_TOL)    = ',
     &           BMG_rPARMS(id_BMG2_STOP_TOL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_CG_MIN_DIM)  = ', 
     &           BMG_iPARMS(id_BMG2_CG_MIN_DIM)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_CG_SOLVER)  = ', 
     &           BMG_iPARMS(id_BMG2_CG_SOLVER)


      WRITE(*,*)
      WRITE(*,*) '****  I/O Parameters: 3D '
      WRITE(*,*) 
      WRITE(*,*) 'BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE)  = ',
     &           BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE)
      WRITE(*,*) 'BMG_IOFLAG(iBMG3_OUT_WSPACE_POINT) = ',
     &           BMG_IOFLAG(iBMG3_OUT_WSPACE_POINT)
      WRITE(*,*) 'BMG_IOFLAG(iBMG3_OUT_TIME_SETUP)   = ',
     &           BMG_IOFLAG(iBMG3_OUT_TIME_SETUP)
      WRITE(*,*) 'BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) = ',
     &           BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING)
      WRITE(*,*) 'BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL)   = ',
     &           BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL)
      WRITE(*,*) 'BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = ',
     &           BMG_IOFLAG(iBMG3_OUT_ITERATIONS)


      WRITE(*,*)
      WRITE(*,*) '****  I/O Parameters: 2D '
      WRITE(*,*) 
      WRITE(*,*) 'BMG_IOFLAG(iBMG2_OUT_WSPACE_SIZE)  = ',
     &           BMG_IOFLAG(iBMG2_OUT_WSPACE_SIZE)
      WRITE(*,*) 'BMG_IOFLAG(iBMG2_OUT_WSPACE_POINT) = ',
     &           BMG_IOFLAG(iBMG2_OUT_WSPACE_POINT)
      WRITE(*,*) 'BMG_IOFLAG(iBMG2_OUT_TIME_SETUP)   = ',
     &           BMG_IOFLAG(iBMG2_OUT_TIME_SETUP)
      WRITE(*,*) 'BMG_IOFLAG(iBMG2_OUT_TIME_CYCLING) = ',
     &           BMG_IOFLAG(iBMG2_OUT_TIME_CYCLING)
      WRITE(*,*) 'BMG_IOFLAG(iBMG2_OUT_TIME_TOTAL)   = ',
     &           BMG_IOFLAG(iBMG2_OUT_TIME_TOTAL)
      WRITE(*,*) 'BMG_IOFLAG(iBMG2_OUT_ITERATIONS)   = ',
     &           BMG_IOFLAG(iBMG2_OUT_ITERATIONS)

C ==========================================================================

      RETURN
      END
