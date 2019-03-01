      SUBROUTINE BMG3_SER_SymStd_relax_GS( 
     &                kg, so, qf, q, RES, sor, ii, jj, kk,
     &                RES_L2, BMG_IOFLAG, NOG, ifd, NStncl, NSORv, 
     &                iRELAX, iRELAX_SYM, UPDOWN 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_relax_GS performs one sweep of Gauss Seidel (with the
C     correct ordering depending on whether we are on the way down, or
C     up in the symmetric cycling case)
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

C ------------------------------------------------
C     Includes
C
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER ii, jj, kk, NSORv, NStncl

      INTEGER ifd, iRELAX, iRELAX_SYM, kg, NOG, UPDOWN
      REAL*RKIND  q(ii,jj,kk), qf(ii,jj,kk), RES(ii,jj,kk),
     &        so(ii,jj,kk,NStncl), sor(ii,jj,kk,NSORv)
      REAL*8 RES_L2
      LOGICAL BMG_IOFLAG(NBMG_SER_IOFLAG)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, i1, ibeg, iend, 
     &        j, j1,   k, k1,
     &        ptstart, ptend, ptstride,
     &        pts

C ==========================================================================

      j1 = jj-1
      i1 = ii-1
      k1 = kk-1

      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN 
         !
         !   27-point relaxation (8-color)
         !

         IF ( UPDOWN.EQ.BMG_SER_UP 
     &       .OR. IRELAX_SYM.EQ.BMG_SER_RELAX_NONSYM ) THEN
            ptstart = 1
            ptend   = 8
            ptstride = 1
         ELSE
            ptstart = 8
            ptend = 1
            ptstride = -1
         ENDIF         
        
         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

            DO k=2+mod((pts-1)/4,2),K1,2
               !
               DO j=2+mod(mod((pts-1)/2,2),2),J1,2
                  !
                  DO i=2+mod(pts-1,2),I1,2

                     q(i,j,k) = ( qf(i,j,k) 
     &                    + so(i,j,k,kpw)*q(i-1,j,k)
     &                    + so(i,j+1,k,kpnw)*q(i-1,j+1,k)
     &                    + so(i,j+1,k,kps)*q(i,j+1,k)
     &                    + so(i+1,j+1,k,kpsw)*q(i+1,j+1,k)
     &                    + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                    + so(i+1,j,k,kpnw)*q(i+1,j-1,k)
     &                    + so(i,j,k,kps)*q(i,j-1,k)
     &                    + so(i,j,k,kpsw)*q(i-1,j-1,k)
     &                    + so(i,j,k,kb)*q(i,j,k-1)
     &                    + so(i,j,k,kbw)*q(i-1,j,k-1)
     &                    + so(i,j+1,k,kbnw)*q(i-1,j+1,k-1)
     &                    + so(i,j+1,k,kbn)*q(i,j+1,k-1)
     &                    + so(i+1,j+1,k,kbne)*q(i+1,j+1,k-1)
     &                    + so(i+1,j,k,kbe)*q(i+1,j,k-1)
     &                    + so(i+1,j,k,kbse)*q(i+1,j-1,k-1)
     &                    + so(i,j,k,kbs)*q(i,j-1,k-1)
     &                    + so(i,j,k,kbsw)*q(i-1,j-1,k-1)
     &                    + so(i,j,k+1,kb)*q(i,j,k+1)
     &                    + so(i,j,k+1,kbe)*q(i-1,j,k+1)
     &                    + so(i,j+1,k+1,kbse)*q(i-1,j+1,k+1)
     &                    + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)
     &                    + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)
     &                    + so(i+1,j,k+1,kbw)*q(i+1,j,k+1)
     &                    + so(i+1,j,k+1,kbnw)*q(i+1,j-1,k+1)
     &                    + so(i,j,k+1,kbn)*q(i,j-1,k+1)
     &                    + so(i,j,k+1,kbne)*q(i-1,j-1,k+1)
     &                    )*sor(i,j,k,msor)
                  END DO
      
               END DO
      
            END DO

         ENDDO ! loop over colors
      
      ELSE
         !  
         !  7 point relaxation ( four colors )
         !
         IF ( (UPDOWN.eq.BMG_SER_UP)
     &        .OR. (IRELAX_SYM.EQ.BMG_SER_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF         

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod(j+k+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2

                     q(i,j,k) = ( qf(i,j,k)
     &                    + so(i,j,k,kpw)*q(i-1,j,k)
     &                    + so(i,j+1,k,kps)*q(i,j+1,k)
     &                    + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                    + so(i,j,k,kps)*q(i,j-1,k)
     &                    + so(i,j,k,kb)*q(i,j,k-1)
     &                    + so(i,j,k+1,kb)*q(i,j,k+1)
     &                    )*sor(i,j,k,msor)

                  END DO

               END DO
               !
            END DO
            !
         END DO ! loop over colors 
c$$$            CALL BMG3_SER_SymStd_DUMP_vector( 
c$$$     &                BMG_IOFLAG, Q, II, JJ, KK, kg, NOG
c$$$     &           )
c$$$
c$$$            STOP
         !
      ENDIF

      IF( BMG_IOFLAG(iBMG3_SER_BUG_RES_RELAX) ) THEN
         CALL BMG3_SER_SymStd_residual( 
     &             kg, NOG, IFD, Q, QF, SO, RES, II, JJ, KK, NStncl
     &             )
         CALL BMG3_SER_SymStd_UTILS_norm_l2( 
     &             RES, II, JJ, KK, RES_L2
     &             )
      ENDIF

C ==========================================================================

      RETURN
      END

