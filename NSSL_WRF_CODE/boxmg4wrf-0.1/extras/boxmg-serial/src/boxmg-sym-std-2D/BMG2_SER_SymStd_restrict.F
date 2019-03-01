      SUBROUTINE BMG2_SER_SymStd_restrict( 
     &                       KF, KC, Q, QC, CI, 
     &                       Nx, Ny, Nxc, Nyc, IBC 
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_restrict computes the restriction of a vector on the
C     fine grid, Q, to a vector on the coarse grid, QC.  The weights
C     involve the transpose of the interpolation operator from the
C     coarse grid to the fine grid.
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
      INTEGER  Nx, Nxc, Ny, Nyc

      INTEGER  IBC, KC, KF
      REAL*RKIND   CI(Nxc,Nyc,8), Q(Nx,Ny), QC(Nxc,Nyc)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, ic, j, jc

C ==========================================================================

      
C ---------------------------------------------
C     Restrict the vector Q -> QC:
C ---------------------------------------------

      j=0
      DO jc=2, Nyc-1
         j=j+2
         i=0
         DO ic=2, Nxc-1
            i=i+2
            QC(ic,jc) = CI(ic,jc,LNE)*Q(i-1,j-1)
     &                + CI(ic,jc,LA)*Q(i,j-1)
     &                + CI(ic+1,jc,LNW)*Q(i+1,j-1)
     &                + CI(ic,jc,LR)*Q(i-1,j)
     &                + Q(i,j)
     &                + CI(ic+1,jc,LL)*Q(i+1,j)
     &                + CI(ic,jc+1,LSE)*Q(i-1,j+1)
     &                + CI(ic,jc+1,LB)*Q(i,j+1)
     &                + CI(ic+1,jc+1,LSW)*Q(i+1,j+1)
          ENDDO 
       ENDDO

C ===========================================================================

       RETURN
       END
