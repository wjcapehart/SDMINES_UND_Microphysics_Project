      SUBROUTINE BMG3_SER_SymStd_residual(
     &                kg, NOG, ifd, q, qf, so, RES, ii, jj, kk, NStncl
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_residual computes the residual on level kg.
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
      INTEGER  ifd, ii, jj, kg, kk, NOG, NStncl
      REAL*RKIND   q(ii,jj,kk), qf(ii,jj,kk), RES(ii,jj,kk),
     &         SO(ii,jj,kk,NStncl)

C ----------------------------
C     Local Declarations
C
      integer  i, i1, j, j1, k, k1

C ==========================================================================

      i1=ii-1
      j1=jj-1
      k1=kk-1

      IF( kg.lt.NOG .or. IFD.ne.1 ) THEN

         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  RES(i,j,k) = qf(i,j,k)
     &                       + so(i,j,k,kpw)*q(i-1,j,k)
     &                       + so(i,j+1,k,kpnw)*q(i-1,j+1,k)
     &                       + so(i,j+1,k,kps)*q(i,j+1,k)
     &                       + so(i+1,j+1,k,kpsw) *q(i+1,j+1,k)
     &                       + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                       + so(i+1,j,k,kpnw)*q(i+1,j-1,k)
     &                       + so(i,j,k,kps)*q(i,j-1,k)
     &                       + so(i,j,k,kpsw)*q(i-1,j-1,k)
     &                       + so(i,j,k,kb)*q(i,j,k-1)
     &                       + so(i,j,k,kbw)*q(i-1,j,k-1)
     &                       + so(i,j+1,k,kbnw)*q(i-1,j+1,k-1)
     &                       + so(i,j+1,k,kbn)*q(i,j+1,k-1)
     &                       + so(i+1,j+1,k,kbne)*q(i+1,j+1,k-1)
     &                       + so(i+1,j,k,kbe)*q(i+1,j,k-1)
     &                       + so(i+1,j,k,kbse)*q(i+1,j-1,k-1)
     &                       + so(i,j,k,kbs)*q(i,j-1,k-1)
     &                       + so(i,j,k,kbsw)*q(i-1,j-1,k-1)
     &                       + so(i,j,k+1,kb)*q(i,j,k+1)
     &                       + so(i,j,k+1,kbe)*q(i-1,j,k+1)
     &                       + so(i,j+1,k+1,kbse)*q(i-1,j+1,k+1)
     &                       + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)
     &                       + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)
     &                       + so(i+1,j,k+1,kbw)*q(i+1,j,k+1)
     &                       + so(i+1,j,k+1,kbnw)*q(i+1,j-1,k+1)
     &                       + so(i,j,k+1,kbn)*q(i,j-1,k+1)
     &                       + so(i,j,k+1,kbne)*q(i-1,j-1,k+1)
     &                       - so(i,j,k,kp)*q(i,j,k)
               ENDDO
            ENDDO
         ENDDO

      ELSE

         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  RES(i,j,k) = qf(i,j,k)
     &                       + so(i,j,k,kpw)*q(i-1,j,k)
     &                       + so(i,j+1,k,kps)*q(i,j+1,k)
     &                       + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                       + so(i,j,k,kps)*q(i,j-1,k)
     &                       + so(i,j,k,kb)*q(i,j,k-1)
     &                       + so(i,j,k+1,kb)*q(i,j,k+1)
     &                       - so(i,j,k,kp)*q(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         
      ENDIF

C ==========================================================================

      RETURN
      END
