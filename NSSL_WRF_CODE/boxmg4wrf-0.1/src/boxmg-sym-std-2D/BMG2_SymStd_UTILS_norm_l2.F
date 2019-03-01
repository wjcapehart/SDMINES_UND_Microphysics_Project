      SUBROUTINE BMG2_SymStd_UTILS_norm_l2( u, Nx, Ny, l2norm, MPICOMM )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_norm_l2.f computes the l2-norm of a grid function
C     (the vector u).  It assumes that ghost points should be neglected.
C
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
      REAL*RKIND   l2norm, u(Nx,Ny)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, ierror
      REAL*RKIND   buffer

      INTEGER  MyProc, MPI_IERROR

C ==========================================================================

      l2norm=rZERO

      DO j = 2, Ny-1
         DO i = 2, Nx-1
            l2norm=l2norm+u(i,j)**2D0
         END DO
      END DO
      
      CALL MPI_AllReduce( l2norm, buffer, 1, FLOAT_PRECISION,
     &                    MPI_SUM, MPICOMM, ierror )

      l2norm=SQRT(buffer)

C ==========================================================================

      RETURN
      END
