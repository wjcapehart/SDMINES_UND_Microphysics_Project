      SUBROUTINE BMG2_SymStd_UTILS_dot_l2( 
     &               u, v, Nx, Ny, l2dot, MPICOMM 
     &               )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_dot_l2.f computes the dot product of a grid function u
C     with grid function v.  It assumes that ghost points should be
C     neglected.
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
#if RKIND == 4
#define FLOAT_PRECISION MPI_REAL
#else
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION
#endif

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'
      
#include      "BMG_constants.h"

C ----------------------------
C     Argument Declarations
C 
      INTEGER  Nx, Ny, MPICOMM
      REAL*RKIND   l2dot, u(Nx,Ny), v(Nx,Ny)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, ierror
      REAL*RKIND   buffer

C ==========================================================================

      l2dot=rZERO

      DO j = 2, Ny-1
         DO i = 2, Nx-1
            l2dot=l2dot+u(i,j)*v(i,j)
         END DO
      END DO

      call MPI_AllReduce(l2dot, buffer,1,    
     &     FLOAT_PRECISION,MPI_SUM,MPICOMM,ierror)

      l2dot=buffer

C ==========================================================================

      RETURN
      END
