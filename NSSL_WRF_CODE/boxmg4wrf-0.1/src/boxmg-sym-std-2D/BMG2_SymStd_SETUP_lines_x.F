      SUBROUTINE BMG2_SymStd_SETUP_lines_x( SO, SOR, Nx, Ny, NStncl)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG2_SymStd_SETUP_lines_x performs a factorization of the tridiagonal
C   matrix that arises in x-line relaxation.  It assumes that the system
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
      REAL*RKIND    SO(Nx+1,Ny+1,NStncl), SOR(Nx,Ny,2)

C ----------------------------
C     Local Declarations
C
      INTEGER i, j

C ==========================================================================

      DO j=2, Ny-1
         DO i=1, Nx
            SOR(i,j,2) = -SO(i,j,KW)  ! off diagonal
            SOR(i,j,1) =  SO(i,j,KO)  ! diagonal
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END

