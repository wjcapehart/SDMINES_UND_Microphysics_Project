      SUBROUTINE BMG3_SER_SymStd_UTILS_matvec(
     &                Q, QF, SO, Nx, Ny, Nz, NStncl
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_UTILS_matvec performs a matrix multiplication that is
C     used within the pcg routine called BMG3_SER_SymStd_SOLVE_pcg.
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
      INTEGER IFD, kg, Nx, Ny, Nz, NOG, NStncl
      REAL*RKIND  Q(Nx,Ny,Nz), QF(Nx,Ny,Nz), SO(Nx,Ny,Nz,NStncl)

C ----------------------------
C     Local Declarations
C
      INTEGER i, j, k

C ==========================================================================

      IF( NStncl.EQ.14 ) THEN

         DO k=2, Nz-1
            DO j=2, Ny-1
               DO i=2, Nx-1

                  QF(i,j,k) =  SO(i,j,k,kp)*Q(i,j,k)
     &                       - SO(i,j,k,kpw)*Q(i-1,j,k)
     &                       - SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                       - SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                       - SO(i+1,j+1,k,kpsw) *Q(i+1,j+1,k)
     &                       - SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                       - SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                       - SO(i,j,k,kps)*Q(i,j-1,k)
     &                       - SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                       - SO(i,j,k,kb)*Q(i,j,k-1)
     &                       - SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                       - SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                       - SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                       - SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                       - SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                       - SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                       - SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                       - SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                       - SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                       - SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                       - SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                       - SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                       - SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                       - SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                       - SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                       - SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                       - SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
               ENDDO
            ENDDO
         ENDDO

      ELSE

         DO k=2, Nz-1
            DO j=2,Ny-1
               DO i=2, Nx-1

                  QF(i,j,k) =  SO(i,j,k,kp)*Q(i,j,k)
     &                       - SO(i,j,k,kpw)*Q(i-1,j,k)
     &                       - SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                       - SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                       - SO(i,j,k,kps)*Q(i,j-1,k)
     &                       - SO(i,j,k,kb)*Q(i,j,k-1)
     &                       - SO(i,j,k+1,kb)*Q(i,j,k+1)

               ENDDO
            ENDDO
         ENDDO
         
      ENDIF

C ==========================================================================

      RETURN
      END
