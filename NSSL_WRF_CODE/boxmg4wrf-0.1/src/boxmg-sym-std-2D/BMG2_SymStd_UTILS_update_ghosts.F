      SUBROUTINE BMG2_SymStd_UTILS_update_ghosts(
     &                       K, x, Nx, Ny, iWork, pMSG,
     &                       buffer, MPICOMM
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG2_SymStd_UTILS_update_ghosts.
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
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  K, Nx, Ny, MPICOMM, 
     &         iWork(*), pMSG(NBMG_pMSG,*)
      REAL*RKIND   x(Nx,Ny), buffer(*)

C ----------------------------
C     Local Declarations
C
      INTEGER  ptrn, ierror

C ==========================================================================

      ! Note: we need to update the ghost bdry of y here,
      ! since we eventually use in matrix multiplication.
      
      ptrn = 1

      call MSG_tbdx_send(x, buffer, 
     &     iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &     iWork(pMSG(ipL_MSG_Proc,K)),
     &     iWork(pMSG(ipL_MSG_Ipr,K)),
     &     iWork(pMSG(ipL_MSG_Index,K)),
     &     ptrn, ierror)

      call MSG_tbdx_receive(x, buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &     iWork(pMSG(ipL_MSG_Proc,K)),
     &     iWork(pMSG(ipL_MSG_Ipr,K)),
     &     iWork(pMSG(ipL_MSG_Index,K)),
     &     ptrn, ierror)
      
      call MSG_tbdx_close(x, buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &     iWork(pMSG(ipL_MSG_Proc,K)),
     &     iWork(pMSG(ipL_MSG_Ipr,K)),
     &     iWork(pMSG(ipL_MSG_Index,K)),
     &     ptrn, ierror)
      
C ==========================================================================

      RETURN
      END
