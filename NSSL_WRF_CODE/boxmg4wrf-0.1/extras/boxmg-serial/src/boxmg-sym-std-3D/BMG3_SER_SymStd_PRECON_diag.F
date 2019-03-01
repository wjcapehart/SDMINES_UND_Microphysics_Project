      SUBROUTINE BMG3_SER_SymStd_PRECON_diag( 
     &                SO, QF, Q, II, JJ, KK, NStncl
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_PRECON_diag solves the system D*Q = QF where D is the
C     diagonal of the matrix in SO.
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
      INTEGER  II, JJ, KK, NStncl
      REAL*RKIND   Q(II,JJ,KK), QF(II,JJ,KK), SO(II,JJ,KK,NStncl)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, I1, J, J1, K, K1, IERR
      REAL*RKIND   tiny

      PARAMETER ( tiny = 1e-10 )

C =========================================================================

      K1=KK-1
      J1=JJ-1
      I1=II-1

C -------------------------------------------------------------------

C     ----------------------------------------
C     PRE-CHECK DIAGONAL FOR ZEROS
C     ----------------------------------------

      DO K=2,K1
         DO J=2,J1
            DO I=2,I1
               IF ( ABS(SO(I,J,K,KO)).LT.tiny ) THEN  
                  WRITE(*,*)'Error in BMG3_SER_SymStd_PRECON_diag:  '
                  WRITE(*,*)'divide by zero.  Set D = I'
                  IERR = 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(IERR.EQ.0) THEN

C       ----------------------------------------
C       SOLVE D*Q = QF where D is diagonal of A.
C       ----------------------------------------

         DO K=2,K1
            DO J=2,J1
               DO I=2,I1
                  Q(I,J,K) = QF(I,J,K)/SO(I,J,K,KO)  
               ENDDO
            ENDDO
         ENDDO

      ELSE
         
C       ----------------------------------------
C       SET Q = QF by calling dcopy routine.
C       ----------------------------------------
        
         CALL BMG3_SER_SymStd_UTILS_dcopy(
     &             QF, Q, II, JJ, KK 
     &             )

      ENDIF


C =========================================================================

      RETURN
      END

