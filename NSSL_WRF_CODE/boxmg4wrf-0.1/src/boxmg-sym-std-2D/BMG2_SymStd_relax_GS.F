      SUBROUTINE BMG2_SymStd_relax_GS ( 
     &                K, SO, QF, Q, SOR, II, JJ, 
     &                KF, IFD, NStncl, IRELAX_SYM, UPDOWN,
     &                iGs, jGs, iWork, NMSGi, pMSG, 
     &                buffer, NMSGr, MPICOMM
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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER II, JJ, NMSGi, NMSGr, NStncl, MPICOMM,
     &        IFD, IRELAX_SYM, K, KF, UPDOWN, iGs, jGs
      INTEGER iWork(NMSGi), pMSG(NBMG_pMSG,KF)

      REAL*RKIND Q(II,JJ), QF(II,JJ), SO(II+1,JJ+1,NStncl),SOR(II,JJ,*)
      REAL*RKIND  buffer(NMSGr)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, I2, IBEG, J, J1, J2
      INTEGER IB(4), JB(4), pts

      INTEGER ptstart, ptend, ptstride, ILSTART
      INTEGER ptrn, ierror

C ==========================================================================

      J1=JJ-1
      I1=II-1
      J2=JJ-2
      I2=II-2


      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN

         !
         ! 9-point stencil
         !
         IB(mod(iGs,2)+2*mod(jGs,2)+1) = 2
         JB(mod(iGs,2)+2*mod(jGs,2)+1) = 2

         IB(mod(iGs+1,2)+2*mod(jGs,2)+1) = 3
         JB(mod(iGs+1,2)+2*mod(jGs,2)+1) = 2

         IB(mod(iGs,2)+2*mod(jGs+1,2)+1) = 2
         JB(mod(iGs,2)+2*mod(jGs+1,2)+1) = 3

         IB(mod(iGs+1,2)+2*mod(jGs+1,2)+1) = 3
         JB(mod(iGs+1,2)+2*mod(jGs+1,2)+1) = 3
         
         IF (IRELAX_SYM.EQ.BMG_RELAX_SYM .AND. UPDOWN.EQ.BMG_UP) THEN
            ptstart = 1
            ptend   = 4
            ptstride = 1
         ELSE
            ptstart = 4
            ptend = 1
            ptstride = -1
         ENDIF         

         
         DO pts=ptstart,ptend,ptstride
            DO  J=JB(pts),J1,2
               DO  I=IB(pts),I1,2
                  Q(I,J)=(QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)
     &                 *Q(I+1,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)
     &                 *Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)
     &                 +SO(I+1,J,KNW)*Q(I+1,J-1)+SO(I,J+1,KNW)
     &                 *Q(I-1,J+1)+SO(I+1,J+1,KSW)*Q(I+1,J+1))
     &                 *SOR(I,J,MSOR)

               ENDDO
            ENDDO

            ptrn = 1
            
            call MSG_tbdx_send(Q, buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &           iWork(pMSG(ipL_MSG_Proc,K)),
     &           iWork(pMSG(ipL_MSG_Ipr,K)),
     &           iWork(pMSG(ipL_MSG_Index,K)),
     &           ptrn, ierror)
            
            call MSG_tbdx_receive(Q, buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &           iWork(pMSG(ipL_MSG_Proc,K)),
     &           iWork(pMSG(ipL_MSG_Ipr,K)),
     &           iWork(pMSG(ipL_MSG_Index,K)),
     &           ptrn, ierror)
            
         ENDDO

         !
      ELSE
         !
         ! 5-point stencil
         !

         IF  (    ((mod(iGs,2).eq.1).and.(mod(jGs,2).eq.1))
     &        .or.((mod(iGs,2).eq.0).and.(mod(jGs,2).eq.0)) ) THEN
            ILSTART = 2
         ELSE
            ILSTART = 3
         ENDIF

         IF (IRELAX_SYM.EQ.BMG_RELAX_SYM .AND. UPDOWN.EQ.BMG_UP) THEN
            ptstart = 1
            ptend = 0            
            ptstride = -1
         ELSE
            ptstart = 0
            ptend = 1
            ptstride = 1
         ENDIF
         
         
         DO pts=ptstart, ptend, ptstride
c         DO pts=0,1,1
            DO j = 2, j1
               
               IF (mod(pts+j,2).eq.0) THEN
                  IBEG = ILSTART
               ELSE
                  IBEG = mod(ILSTART+1,2)+2
               ENDIF
               
               DO  I=IBEG,I1,2
                  Q(I,J)=(QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)
     &                 *Q(I+1,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)
     &                 *Q(I,J+1))*SOR(I,J,MSOR)
               ENDDO
            ENDDO      


            ptrn = 1
            
            call MSG_tbdx_send(Q, buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &           iWork(pMSG(ipL_MSG_Proc,K)),
     &           iWork(pMSG(ipL_MSG_Ipr,K)),
     &           iWork(pMSG(ipL_MSG_Index,K)),
     &           ptrn, ierror)
            
            call MSG_tbdx_receive(Q, buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &           iWork(pMSG(ipL_MSG_Proc,K)),
     &           iWork(pMSG(ipL_MSG_Ipr,K)),
     &           iWork(pMSG(ipL_MSG_Index,K)),
     &           ptrn, ierror)

         ENDDO
      ENDIF


      ptrn = 1

      call MSG_tbdx_close(Q, buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,K)),
     &     iWork(pMSG(ipL_MSG_Proc,K)),
     &     iWork(pMSG(ipL_MSG_Ipr,K)),
     &     iWork(pMSG(ipL_MSG_Index,K)),
     &     ptrn, ierror)

C ==========================================================================

      RETURN
      END
