      SUBROUTINE BMG2_SER_SymStd_UTILS_matvec( 
     &                K, SO, QF, Q, II, JJ,
     &                KF, IFD, NStncl, IBC,
     &                IRELAX, IRELAX_SYM, UPDOWN 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG_SER_SymStd_UTILS_matvec computes the matrix-vector product
C     for the CG routine using the BoxMG stencil.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     K         K is the grid number.
C     SO        Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     QF        Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     RES       Refer to BMG2_SER_SymStd_SOLVE_boxmg      
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C     KF        index of the finest grid
C     IFD       Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     IRELAX    Refer to BMG2_SER_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     Q         Refer to BOXMG.
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
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
      INTEGER II, JJ, NStncl

      INTEGER IBC, IFD, IRELAX, IRELAX_SYM, K, KF, UPDOWN
      REAL*RKIND  Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1

C =========================================================================

      J1=JJ-1
      I1=II-1

C     -------------------------------------------------------------------

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         ! 
         DO J=2,J1
            DO I=2,I1
               QF(I,J) =   SO(I  ,J  ,KO )*Q(I  ,J)
     &                   - SO(I  ,J  ,KW )*Q(I-1,J)
     &                   - SO(I+1,J  ,KW )*Q(I+1,J)
     &                   - SO(I  ,J  ,KS )*Q(I  ,J-1)
     &                   - SO(I  ,J+1,KS )*Q(I  ,J+1)
     &                   - SO(I  ,J  ,KSW)*Q(I-1,J-1)
     &                   - SO(I+1,J  ,KNW)*Q(I+1,J-1)
     &                   - SO(I  ,J+1,KNW)*Q(I-1,J+1)
     &                   - SO(I+1,J+1,KSW)*Q(I+1,J+1)
             ENDDO
         ENDDO
         !
      ELSE
         !
         !  5-point stencil
         !        
         DO J=2,J1
            DO I=2,I1
               QF(I,J) =  SO(I  ,J  ,KO)*Q(I  ,J)
     &                  - SO(I  ,J  ,KW)*Q(I-1,J)
     &                  - SO(I+1,J  ,KW)*Q(I+1,J)
     &                  - SO(I  ,J  ,KS)*Q(I  ,J-1)
     &                  - SO(I  ,J+1,KS)*Q(I  ,J+1)

            ENDDO
         ENDDO
         !
      ENDIF

C =========================================================================

      RETURN
      END

