      SUBROUTINE BMG3_SymStd_SETUP_fine_stencil( 
     &                KF, SO, 
     &                NLx, NLy, NLz, NStncl,
     &                iWork, NMSGi, pMSGSO, BUFFER, NMSGr,
     &                MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Perform communication setup (ghosts) for the fine-grid stencil.
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
      
C ---------------------------
C    Argument Declarations:
C
      INTEGER   NLx, NLy, NLz, NMSGi, NMSGr, NStncl

      INTEGER   KF, iWork(NMSGi), pMSGSO(NBMG_pMSG,KF), MPICOMM 
      REAL*RKIND    BUFFER(NMSGr), SO(NLx+1,NLy+1,NLz+1,NStncl)

C --------------------------
C     Local Declarations:
C
      INTEGER   kst, ptrn, ierror

C ==========================================================================

C ------------------------------
C     Update halo region:
C ------------------------------

      DO kst = 1, NStncl

         ptrn = 1
         
         CALL MSG_tbdx_send(SO(1,1,1,kst), buffer, 
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),
     &        ptrn, ierror)
         
         CALL MSG_tbdx_receive(SO(1,1,1,kst), buffer,
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),
     &        ptrn, ierror)
         
         CALL MSG_tbdx_close(SO(1,1,1,kst), buffer,
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),
     &        ptrn, ierror)

      ENDDO

C ==========================================================================

      RETURN
      END
