      SUBROUTINE BMG2_SER_SymStd_relax_lines_y( 
     &                K, SO, QF, Q, SOR, B, II, JJ, 
     &                KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Perform zebra-line relaxation in y.
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
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER II, JJ, NStncl

      INTEGER IFD, IRELAX_SYM, K, KF, UPDOWN
      REAL*RKIND  B(JJ), Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl), 
     &            SOR(JJ,II,2)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, IBEG, IEND, J, J1
      INTEGER IBEG_START, IBEG_END, IBEG_STRIDE
      INTEGER INFO

C =========================================================================

      J1=JJ-1
      I1=II-1
      
      IF (IRELAX_SYM.EQ.BMG_SER_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_SER_DOWN 
     &     .AND.IRELAX_SYM.EQ.BMG_SER_RELAX_SYM)) THEN
         !
         ! Relax red lines, then black lines.
         !     
         IBEG_START = 3
         IBEG_END   = 2
         IBEG_STRIDE= -1
      ELSE
         !
         ! Relax black lines, then red lines.
         !     
         IBEG_START = 2
         IBEG_END   = 3
         IBEG_STRIDE= 1
      ENDIF

      
      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !     
         !  9 pt. operator
         !  Relax red lines, then black lines.
         !
         DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE
            IEND=2*((I1-IBEG)/2)+IBEG
            DO  I=IBEG,IEND,2
               DO  J=2,J1
                  B(J)= QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)
     &                 *Q(I+1,J)+SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)
     &                 *Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)
     &                 +SO(I+1,J+1,KSW)*Q(I+1,J+1)
               ENDDO

#if RKIND == 4
               CALL SPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), 
     &              B(2), J1-1, INFO)
#else
               CALL DPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), 
     &              B(2), J1-1, INFO)
#endif
               DO j=2,J1
                  Q(I,J) = B(J)
               ENDDO
            ENDDO

c            DO  J=2,J1
c               DO  I=IBEG,IEND,2
c                  Q(I,J)=Q(I,J)+SO(I,J,KS)*SOR(I,J-1,MSOS)*Q(I,J-1)
c               ENDDO
c            ENDDO
c            DO  J=2,J1
c               DO  I=IBEG,IEND,2
c                  Q(I,J1-J+2)=SOR(I,J1-J+2,MSOS)*(Q(I,J1-J+2)
c     &                 +SO(I,J1-J+3,KS)*Q(I,J1-J+3))
c               ENDDO
c            ENDDO
         ENDDO
         !
      ELSE
         !     
         ! 5 pt. operator
         ! Relax red lines, then black lines.
         !
         DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE
            IEND=2*((I1-IBEG)/2)+IBEG
            DO J=2,J1
               DO I=IBEG,IEND,2
                  Q(I,J)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)
     &                 *Q(I+1,J)
               ENDDO
            ENDDO

            DO I=IBEG,IEND,2
               DO J=2,J1
                  B(J) = Q(I,J)
               ENDDO
#if RKIND == 4
               CALL SPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), 
     &              B(2), J1-1, INFO)
#else
               CALL DPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), 
     &              B(2), J1-1, INFO)
#endif
               DO j=2,J1
                  Q(I,J) = B(J)
               ENDDO
            ENDDO
            
c            DO J=2,J1
c               DO I=IBEG,IEND,2
c                  Q(I,J)=Q(I,J)+SO(I,J,KS)*SOR(I,J-1,MSOS)*Q(I,J-1)
c               ENDDO
c            ENDDO
c            DO J=2,J1
c               DO I=IBEG,IEND,2
c                  Q(I,J1-J+2)=SOR(I,J1-J+2,MSOS)*(Q(I,J1-J+2)
c     &                 +SO(I,J1-J+3,KS)*Q(I,J1-J+3))
c               ENDDO
c            ENDDO
         ENDDO
      ENDIF
      
      RETURN
      END
