      SUBROUTINE BMG2_SER_SymStd_SETUP_lines_y(
     &                SO, SOR, Nx, Ny, NStncl
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_SETUP_lines_y.f performs a factorization of the
C     tridiagonal matrix that arises in y-line relaxation.  It assumes
C     that the system is diagonally dominant and it works directly with
C     the stencil.
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
      REAL*RKIND    SO(Nx,Ny,NStncl), SOR(Ny,Nx,2)

C ----------------------------
C     Local Declarations
C
      INTEGER i, j, INFO

C ==========================================================================

      DO i=2, Nx-1
         DO j=2, Ny-1
            SOR(j,i,2) = -SO(i,j,KS)  ! off diagonal
            SOR(j,i,1) =  SO(i,j,KO)  ! diagonal
         ENDDO
      ENDDO

C     calculate the L*D*L' factorizations for the lines in y-direction
      DO i=2, Nx-1
#if RKIND == 4
         CALL SPTTRF (Ny-2, SOR(2,i,1), SOR(3,i,2), INFO)
#else
         CALL DPTTRF (Ny-2, SOR(2,i,1), SOR(3,i,2), INFO)
#endif
         IF (INFO .NE. 0) THEN
            write(*,*) 'SETUP lines - y INFO = ',INFO
         ENDIF
      ENDDO

c      DO i=2, Nx-1
c         SOR(i,2,MSOS)=rONE/SO(i,2,KO)
c         DO j=3, Ny-1
c            SOR(i,j,MSOS)=rONE
c     &                   /( SO(i,j,KO) - SOR(i,j-1,MSOS)*SO(i,j,KS)**2 )
c         ENDDO
c      ENDDO

C ==========================================================================

      RETURN
      END

