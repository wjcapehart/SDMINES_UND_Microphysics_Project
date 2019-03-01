      SUBROUTINE BMG3_SER_SymStd_UTILS_dot_l2( u, v, Nx, Ny, Nz, l2dot )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_UTILS_dot_l2.f computes the dot product of a grid
C     function u with grid function v.  It assumes that ghost points
C     should be neglected.
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

C ----------------------------
C     Argument Declarations
C 
      INTEGER  Nx, Ny, Nz
      REAL*RKIND   l2dot, u(Nx,Ny,Nz), v(Nx,Ny,Nz)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, k

C =========================================================================

      l2dot=rZERO

      DO k = 2, Nz-1
         DO j = 2, Ny-1
            DO i = 2, Nx-1
               l2dot=l2dot+u(i,j,k)*v(i,j,k)
            END DO
         END DO
      END DO


C =========================================================================

      RETURN
      END
