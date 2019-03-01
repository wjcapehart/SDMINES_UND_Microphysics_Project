      SUBROUTINE BMG2_SymStd_residual( 
     &                K, KF, NOG, RES_L2, 
     &                SO, QF, Q, RES,
     &                II, JJ, IFD, NStncl, IRELAX, IRELAX_SYM,
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG_SymStd_resl2 calculates the l2 residual on grid K.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     K         index of the current grid
C     KF        index of the finest grid
C
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C
C     SO        Refer to BMG2_SymStd_SOLVE_boxmg.
C     QF        Refer to BMG2_SymStd_SOLVE_boxmg.
C
C     RES       Refer to BMG2_SymStd_SOLVE_boxmg.
C     IFD       Refer to BMG2_SymStd_SOLVE_boxmg.
C     IRELAX    Refer to BMG2_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     Q         Refer to BMG2_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
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
      INTEGER II, JJ, NMSGi, NMSGr, NOG, NStncl, MPICOMM

      INTEGER iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOG)

      INTEGER IFD, IRELAX, IRELAX_SYM, K, KF

      REAL*RKIND  Q(II,JJ), QF(II,JJ), RES_L2,
     &        SO(II+1,JJ+1,NStncl), RES(II,JJ), 
     &        MSG_Buffer(NMSGr)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1, ptrn, ierror

C ==========================================================================

      J1=JJ-1
      I1=II-1

C     -------------------------------------------------------------------
     
      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         ! 

         DO J=2,J1
            DO I=2,I1
               RES(I,J) = QF(I,J)
     &                  + SO(I  ,J  ,KW )*Q(I-1,J)
     &                  + SO(I+1,J  ,KW )*Q(I+1,J)
     &                  + SO(I  ,J  ,KS )*Q(I  ,J-1)
     &                  + SO(I  ,J+1,KS )*Q(I  ,J+1)
     &                  + SO(I  ,J  ,KSW)*Q(I-1,J-1)
     &                  + SO(I+1,J  ,KNW)*Q(I+1,J-1)
     &                  + SO(I  ,J+1,KNW)*Q(I-1,J+1)
     &                  + SO(I+1,J+1,KSW)*Q(I+1,J+1)
     &                  - SO(I  ,J  ,KO )*Q(I  ,J)
            ENDDO
         ENDDO      

         !
      ELSE
         !
         !  5-point stencil
         !      
         DO J=2,J1
            DO I=2,I1
               RES(I,J) = QF(I,J)
     &                  + SO(I  ,J  ,KW)*Q(I-1,J)
     &                  + SO(I+1,J  ,KW)*Q(I+1,J)
     &                  + SO(I  ,J  ,KS)*Q(I  ,J-1)
     &                  + SO(I  ,J+1,KS)*Q(I  ,J+1)
     &                  - SO(I  ,J  ,KO)*Q(I  ,J)
            ENDDO
         ENDDO
         !
      ENDIF


      ! Note: we need to update the ghost bdry of RES here,
      ! since in the next step we restrict RES to the next
      ! coarser mesh
      
      ptrn = 1

      call MSG_tbdx_send(RES, MSG_Buffer, 
     &     iWorkMSG(pMSG(ipL_MSG_NumAdjProc,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Proc,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Ipr,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Index,k)),
     &     ptrn, ierror)

      call MSG_tbdx_receive(RES, MSG_Buffer,
     &     iWorkMSG(pMSG(ipL_MSG_NumAdjProc,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Proc,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Ipr,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Index,k)),
     &     ptrn, ierror)

      call MSG_tbdx_close(RES, MSG_Buffer,
     &     iWorkMSG(pMSG(ipL_MSG_NumAdjProc,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Proc,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Ipr,k)),
     &     iWorkMSG(pMSG(ipL_MSG_Index,k)),
     &     ptrn, ierror)

C ==========================================================================
      
      RETURN
      END

