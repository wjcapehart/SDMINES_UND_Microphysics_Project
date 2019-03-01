      SUBROUTINE BMG2_SER_SymStd_relax_GS ( 
     &                K, SO, QF, Q, SOR, II, JJ, 
     &                KF, IFD, NStncl, NSORv, IRELAX_SYM, UPDOWN 
     &                )


C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Perform red-black relaxation on 5-point stencils, and 
C     4-color relaxation on 9-point stencils.
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
      INTEGER II, JJ, NStncl, NSORv

      INTEGER IFD, IRELAX_SYM, K, KF, UPDOWN
      REAL*RKIND Q(II,JJ),QF(II,JJ),SO(II,JJ,NStncl),SOR(II,JJ,NSORv)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, I2, IBEG, IEND, J, J1, J2, JBEG, JEND, JO
      INTEGER LSTART, LEND, LSTRIDE

C =========================================================================

      J1=JJ-1
      I1=II-1
      J2=JJ-2
      I2=II-2

      
      IF (IRELAX_SYM.EQ.BMG_SER_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_SER_DOWN 
     &     .AND.IRELAX_SYM.EQ.BMG_SER_RELAX_SYM)) THEN
         LSTART = 2
         LEND   = 3
         LSTRIDE= 1
      ELSE
         LSTART = 3
         LEND   = 2
         LSTRIDE=-1
      ENDIF
      
      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         ! 9-point stencil
         !
         DO JBEG=LSTART,LEND,LSTRIDE
            JEND=2*((J1-JBEG)/2)+JBEG
            DO  J=JBEG,JEND,2
               DO  IBEG=LSTART,LEND,LSTRIDE
                  IEND=2*((I1-IBEG)/2)+IBEG
                  DO  I=IBEG,IEND,2

                     Q(I,J) = ( QF(I,J) 
     &                         + SO(I,J,KW)*Q(I-1,J)
     &                         + SO(I+1,J,KW)*Q(I+1,J)
     &                         + SO(I,J,KS)*Q(I,J-1)
     &                         + SO(I,J+1,KS)*Q(I,J+1)
     &                         + SO(I,J,KSW)*Q(I-1,J-1)
     &                         + SO(I+1,J,KNW)*Q(I+1,J-1)
     &                         + SO(I,J+1,KNW)*Q(I-1,J+1)
     &                         + SO(I+1,J+1,KSW)*Q(I+1,J+1)
     &                        )*SOR(I,J,MSOR)

                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !
      ELSE
         !
         ! 5-point stencil
         !
         DO JO=LSTART,LEND,LSTRIDE
            DO J=2,J1
               IBEG=MOD(J+JO,2)+2
               IEND=2*((I1-IBEG)/2)+IBEG
               DO I=IBEG,IEND,2

                  Q(I,J) = ( QF(I,J) 
     &                      + SO(I,J,KW)*Q(I-1,J)
     &                      + SO(I+1,J,KW)*Q(I+1,J)
     &                      + SO(I,J,KS)*Q(I,J-1)
     &                      + SO(I,J+1,KS)*Q(I,J+1)
     &                     )*SOR(I,J,MSOR)

               ENDDO
            ENDDO
         ENDDO
         !
      ENDIF

      RETURN
      END
