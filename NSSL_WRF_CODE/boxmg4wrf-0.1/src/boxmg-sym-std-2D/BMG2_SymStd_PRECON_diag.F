      SUBROUTINE BMG2_SymStd_PRECON_diag( 
     &                       K, SO, QF, Q, II, JJ, NStncl,
     &                       iWork, pMSG, buffer, MPICOMM
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_PRECON_diag solves the system D*Q = QF where D is the
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

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER K, II, JJ, NStncl, MPICOMM,
     &        iWork(*), pMSG(NBMG_pMSG,*)

      REAL*RKIND  Q(II,JJ), QF(II,JJ), SO(II+1,JJ+1,NStncl), buffer(*)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1, ERR, ptrn, ierror
      REAL*RKIND  tiny

      PARAMETER ( tiny = 1e-10 )

C ==========================================================================

      J1=JJ-1
      I1=II-1
      ERR=0

C -------------------------------------------------------------------

C     ----------------------------------------
C     PRE-CHECK DIAGONAL FOR ZEROS
C     ----------------------------------------
      DO 10 J=2,J1
      DO 10 I=2,I1
#if RKIND == 4
         IF ( abs(SO(I,J,KO)).LT.tiny ) THEN  
#else
         IF ( dabs(SO(I,J,KO)).LT.tiny ) THEN  
#endif
            WRITE(*,*)'Error in BMG2_SymStd_PRECON_diag:  '
            WRITE(*,*)'divide by zero.  Set D = I'
            ERR = 1
         ENDIF
 10   CONTINUE
      

      IF(ERR.EQ.0) THEN

C       ----------------------------------------
C       SOLVE D*Q = QF where D is diagonal of A.
C       ----------------------------------------

         DO 20 J=2,J1
         DO 20 I=2,I1
            Q(I,J) = QF(I,J)/SO(I,J,KO) 
 20      CONTINUE

         ! Note: we need to update the ghost bdry of y here,
         ! since we eventually use in matrix multiplication.
      
         ptrn = 1

         call MSG_tbdx_send(q, buffer, 
     &        iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &        iWork(pMSG(ipL_MSG_Proc,K)),
     &        iWork(pMSG(ipL_MSG_Ipr,K)),
     &        iWork(pMSG(ipL_MSG_Index,K)),
     &        ptrn, ierror)

         call MSG_tbdx_receive(q, buffer,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &        iWork(pMSG(ipL_MSG_Proc,K)),
     &        iWork(pMSG(ipL_MSG_Ipr,K)),
     &        iWork(pMSG(ipL_MSG_Index,K)),
     &        ptrn, ierror)

         call MSG_tbdx_close(q, buffer,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &        iWork(pMSG(ipL_MSG_Proc,K)),
     &        iWork(pMSG(ipL_MSG_Ipr,K)),
     &        iWork(pMSG(ipL_MSG_Index,K)),
     &        ptrn, ierror)

      ELSE
         
C       ----------------------------------------
C       SET Q = QF by calling dcopy routine.
C       ----------------------------------------
        
         CALL BMG2_SymStd_UTILS_dcopy( QF, Q, II, JJ )

      ENDIF

C ==========================================================================

      RETURN
      END

