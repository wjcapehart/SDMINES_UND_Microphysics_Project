      SUBROUTINE BMG3_SymStd_PRECON_diag( 
     &                       KG, SO, QF, Q, II, JJ, KK, NStncl,
     &                       iWork, pMSG, buffer, MPICOMM
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_PRECON_diag solves the system D*Q = QF where D is the
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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER KG, II, JJ, KK, NStncl, MPICOMM,
     &        iWork(*), pMSG(NBMG_pMSG,*)

      REAL*RKIND  Q(II,JJ,KK), QF(II,JJ,KK), SO(II+1,JJ+1,KK+1,NStncl), 
     &        buffer(*)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1, K, K1, ERR, ptrn, ierror
      REAL*RKIND  tiny

      PARAMETER ( tiny = 1e-10 )

C ==========================================================================

      J1=JJ-1
      I1=II-1
      K1=KK-1
      ERR=0

C -------------------------------------------------------------------

C     ----------------------------------------
C     PRE-CHECK DIAGONAL FOR ZEROS
C     ----------------------------------------
      DO K=2,K1
         DO J=2,J1
            DO I=2,I1
#if RKIND == 4
               IF ( abs(SO(I,J,K,KO)).LT.tiny ) THEN  
#else
               IF ( dabs(SO(I,J,K,KO)).LT.tiny ) THEN  
#endif
                  WRITE(*,*)'Error in BMG3_SymStd_PRECON_diag:  '
                  WRITE(*,*)'divide by zero.  Set D = I'
                  ERR = 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(ERR.EQ.0) THEN

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

         ! Note: we need to update the ghost bdry of y here,
         ! since we eventually use in matrix multiplication.
      
         ptrn = 1

         call MSG_tbdx_send(q, buffer, 
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &        iWork(pMSG(ipL_MSG_Proc,KG)),
     &        iWork(pMSG(ipL_MSG_Ipr,KG)),
     &        iWork(pMSG(ipL_MSG_Index,KG)),
     &        ptrn, ierror)

         call MSG_tbdx_receive(q, buffer,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &        iWork(pMSG(ipL_MSG_Proc,KG)),
     &        iWork(pMSG(ipL_MSG_Ipr,KG)),
     &        iWork(pMSG(ipL_MSG_Index,KG)),
     &        ptrn, ierror)

         call MSG_tbdx_close(q, buffer,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &        iWork(pMSG(ipL_MSG_Proc,KG)),
     &        iWork(pMSG(ipL_MSG_Ipr,KG)),
     &        iWork(pMSG(ipL_MSG_Index,KG)),
     &        ptrn, ierror)

      ELSE
         
C       ----------------------------------------
C       SET Q = QF by calling dcopy routine.
C       ----------------------------------------
        
         CALL BMG3_SymStd_UTILS_dcopy( QF, Q, II, JJ, KK )

      ENDIF

C ==========================================================================

      RETURN
      END

