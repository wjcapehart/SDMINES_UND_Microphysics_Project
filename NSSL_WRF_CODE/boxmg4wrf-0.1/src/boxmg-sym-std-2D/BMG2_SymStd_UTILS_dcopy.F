      SUBROUTINE BMG2_SymStd_UTILS_dcopy( 
     &                       K, x, y, Nx, Ny, PASS, iWork, pMSG, 
     &                       buffer, MPICOMM
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_dcopy.f replaces y with x.
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
      INTEGER  K, Nx, Ny, MPICOMM, PASS,
     &         iWork(*), pMSG(NBMG_pMSG,*)
      REAL*RKIND   x(Nx,Ny), y(Nx,Ny), buffer(*)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, ptrn, ierror

C ==========================================================================

      DO j = 2, Ny-1
         DO i = 2, Nx-1
            y(i,j) = x(i,j)
         END DO
      END DO

      ! Note: we need to update the ghost bdry of y here,
      ! since we eventually use in matrix multiplication.
      
      IF ( PASS.EQ.1 ) THEN

         ptrn = 1

         call MSG_tbdx_send(y, buffer, 
     &        iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &        iWork(pMSG(ipL_MSG_Proc,K)),
     &        iWork(pMSG(ipL_MSG_Ipr,K)),
     &        iWork(pMSG(ipL_MSG_Index,K)),
     &        ptrn, ierror)

         call MSG_tbdx_receive(y, buffer,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &        iWork(pMSG(ipL_MSG_Proc,K)),
     &        iWork(pMSG(ipL_MSG_Ipr,K)),
     &        iWork(pMSG(ipL_MSG_Index,K)),
     &        ptrn, ierror)

         call MSG_tbdx_close(y, buffer,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &        iWork(pMSG(ipL_MSG_Proc,K)),
     &        iWork(pMSG(ipL_MSG_Ipr,K)),
     &        iWork(pMSG(ipL_MSG_Index,K)),
     &        ptrn, ierror)

      ENDIF

C ==========================================================================

      RETURN
      END
