      SUBROUTINE BMG2_SymStd_SETUP_lines_y( SO, SOR, Nx, Ny, NStncl)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG2_SymStd_SETUP_lines_y performs a factorization of the tridiagonal
C   matrix that arises in y-line relaxation.  It assumes that the system
C   is diagonally dominant and it works directly with the stencil.
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
#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER   Nx, Ny, NStncl
      REAL*RKIND    SO(Nx+1,Ny+1,NStncl), SOR(Ny,Nx,2)

C ----------------------------
C     Local Declarations
C
      INTEGER i, j

C ==========================================================================

      DO i=2, Nx-1
         DO j=1, Ny
            SOR(j,i,2) = -SO(i,j,KS)  ! off diagonal
            SOR(j,i,1) =  SO(i,j,KO)  ! diagonal
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END

