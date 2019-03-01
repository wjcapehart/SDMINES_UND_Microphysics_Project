      SUBROUTINE BMG3_SER_SymStd_UTILS_norm_l2( 
     &                u, Nx, Ny, Nz, l2norm
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_UTILS_norm_l2 computes the l2-norm of a grid function
C     (the vector u).  It assumes that ghost points should be neglected.
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
      REAL*RKIND  u(Nx,Ny,Nz)
      REAL*8  l2norm

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, k

C =========================================================================

      l2norm=dZERO

      DO k = 2, Nz-1
         DO j = 2, Ny-1
            DO i = 2, Nx-1
               l2norm=l2norm+u(i,j,k)*u(i,j,k)
            END DO
         END DO
      END DO

C     This seems very wrong.
 
      l2norm=SQRT(l2norm)

C =========================================================================

      RETURN
      END
