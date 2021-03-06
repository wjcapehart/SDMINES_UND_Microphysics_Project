      SUBROUTINE BMG3_SymStd_UTILS_dxpby( 
     &                KG, beta, x, y, Nx, Ny, Nz,
     &                PASS, iWork, pMSG, BUFFER, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_UTILS_dxpby replaces y with x + beta*y.
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
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  KG, Nx, Ny, Nz, MPICOMM, PASS,
     &         iWork(*), pMSG(NBMG_pMSG,*)
      REAL*RKIND   beta, x(Nx,Ny,Nz), y(Nx,Ny,Nz), buffer(*)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, k, ptrn, ierror

C ==========================================================================

      DO k = 2, Nz-1
         DO j = 2, Ny-1
            DO i = 2, Nx-1
               y(i,j,k)=beta*y(i,j,k) + x(i,j,k)
            ENDDO
         ENDDO
      ENDDO

      ! Note: we need to update the ghost bdry of y here,
      ! since we eventually use in matrix multiplication.
      
      IF ( PASS.EQ.1 ) THEN

         ptrn = 1
         
         call MSG_tbdx_send(y, buffer, 
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &        iWork(pMSG(ipL_MSG_Proc,KG)),
     &        iWork(pMSG(ipL_MSG_Ipr,KG)),
     &        iWork(pMSG(ipL_MSG_Index,KG)),
     &        ptrn, ierror)
         
         call MSG_tbdx_receive(y, buffer,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &        iWork(pMSG(ipL_MSG_Proc,KG)),
     &        iWork(pMSG(ipL_MSG_Ipr,KG)),
     &        iWork(pMSG(ipL_MSG_Index,KG)),
     &        ptrn, ierror)

         call MSG_tbdx_close(y, buffer,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &        iWork(pMSG(ipL_MSG_Proc,KG)),
     &        iWork(pMSG(ipL_MSG_Ipr,KG)),
     &        iWork(pMSG(ipL_MSG_Index,KG)),
     &        ptrn, ierror)

      ENDIF

C ==========================================================================

      RETURN
      END
