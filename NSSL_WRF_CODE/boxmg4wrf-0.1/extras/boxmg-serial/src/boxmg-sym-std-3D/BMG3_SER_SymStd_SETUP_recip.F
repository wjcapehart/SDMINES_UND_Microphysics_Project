      SUBROUTINE BMG3_SER_SymStd_SETUP_recip( 
     &                SO, SOR, Nx, Ny, Nz, NStncl, NSORv
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_SETUP_recip computes the reciprocal of the central
C     stencil coefficient for use in GS relaxation.
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
      INCLUDE 'BMG_SER_stencils.h'

C ---------------------------
C    Argument Declarations:
C
      INTEGER   NSORv, NStncl, Nx, Ny, Nz
      REAL*RKIND    SO(Nx,Ny,Nz,NStncl), SOR(Nx,Ny,Nz,NSORv)

C --------------------------
C     Local Declarations:
C
      INTEGER   i, j, k

C ==========================================================================

      DO k=2, Nz-1
         DO j = 2, Ny-1
            DO i = 2, Nx-1
               SOR(i,j,k,msor)=rONE/SO(i,j,k,kp)
            ENDDO
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END
