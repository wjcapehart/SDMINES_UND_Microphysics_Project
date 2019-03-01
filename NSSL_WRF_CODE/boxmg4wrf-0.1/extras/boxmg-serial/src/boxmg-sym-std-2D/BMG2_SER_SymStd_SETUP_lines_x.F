      SUBROUTINE BMG2_SER_SymStd_SETUP_lines_x(
     &                SO, SOR, Nx, Ny, NStncl 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_SETUP_lines_x performs a factorization of the tridiagonal
C     matrix that arises in x-line relaxation.  It assumes that the system
C     is diagonally dominant and it works directly with the stencil.
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

C ----------------------------
C     Argument Declarations
C 
      INTEGER   Nx, Ny, NStncl
      REAL*RKIND    SO(Nx,Ny,NStncl), SOR(Nx,Ny,2)

C ----------------------------
C     Local Declarations
C
      INTEGER i, j, INFO

C ==========================================================================

      DO j=2, Ny-1
         DO i=2, Nx-1
            SOR(i,j,2) = -SO(i,j,KW)  ! off diagonal
            SOR(i,j,1) =  SO(i,j,KO)  ! diagonal
         ENDDO
      ENDDO

C     calculate the L*D*L' factorizations for the lines in x-direction
      DO j=2,Ny-1
#if RKIND == 4
         CALL SPTTRF (Nx-2, SOR(2,j,1), SOR(3,j,2), INFO)
#else
         CALL DPTTRF (Nx-2, SOR(2,j,1), SOR(3,j,2), INFO)
#endif
         IF (INFO .NE. 0) THEN
            WRITE(*,*) 'SETUP lines - x INFO = ',INFO
         ENDIF
      ENDDO

c      DO j=2, Ny-1
c         SOR(2,j,MSOR)=rONE/SO(2,j,KO)
c         DO i=3, Nx-1
c            SOR(i,j,MSOR)=rONE
c     &                   /( SO(i,j,KO) - SOR(i-1,j,MSOR)*SO(i,j,KW)**2 )
c         ENDDO
c      ENDDO

C ==========================================================================

      RETURN
      END

