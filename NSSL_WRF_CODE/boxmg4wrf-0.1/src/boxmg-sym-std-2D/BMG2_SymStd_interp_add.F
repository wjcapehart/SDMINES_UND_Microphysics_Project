      SUBROUTINE BMG2_SymStd_interp_add(
     &                KC, KF, NOG, 
     &                Q ,QC, RES, SO, NStncl, CI,
     &                IIC, JJC, IIF, JJF, iGs, jGs, 
     &                iWorkMSG, NMSGi, pMSG,
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_interp_add.f interpolates Q from the coarse mesh, KC,
C     to the fine mesh, KF, and adds the result to Q on fine mesh.
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

      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER iGs, IIC, IIF, jGs, JJC, JJF, KC, KF, NOG,
     &        NMSGi, NMSGr, NStncl, MPICOMM

      INTEGER iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOG) 
      REAL*RKIND  CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC), 
     &        SO(IIF+1,JJF+1,NStncl), RES(IIF,JJF),
     &        MSG_Buffer(NMSGr)

C ----------------------------
C     Local Declarations
C
      INTEGER IC, I, IICF, IICF1, IIC1, IIF1,
     &        JC, J, JJCF, JJCF1, JJC1, JJF1,
     &        ptrn, ierror, ICSTART, JCSTART,
     &        ISTART, JSTART
      REAL*RKIND  A, AQ

C ==========================================================================

C -------------------------------------------------
C     Useful index bounds:
C -------------------------------------------------


      IIF1=IIF-1
      JJF1=JJF-1

      IIC1=IIC-1
      JJC1=JJC-1

      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      

      IF ((mod(IIF,2).eq.0).and.(mod(iGs,2).eq.0)) THEN
         IICF = IICF - 1
      ENDIF

      IF ((mod(JJF,2).eq.0).and.(mod(jGs,2).eq.0)) THEN
         JJCF = JJCF - 1
      ENDIF
      

      IICF1=IICF-1
      JJCF1=JJCF-1



C --------------------------------------------------
C     NB: division is only possible in the interior
C --------------------------------------------------

      DO j = 2, JJF1
         DO i = 2, IIF1
            RES(i,j)=RES(i,j)/SO(i,j,ko)
         ENDDO
      ENDDO

C --------------------------------------------------
C   interpolate answers from coarse to fine mesh 
C   and add to answers on fine mesh.
C --------------------------------------------------


C     odd values of iGs or jGs indicate the
C     the case similar to the serial version

      IF (mod(iGs,2).eq.1) THEN
         ISTART = 2
         ICSTART = 2
      ELSE
         ISTART = 1
         ICSTART = 1
      ENDIF

      IF (mod(jGs,2).eq.1) THEN
         JSTART = 2
         JCSTART = 2
      ELSE
         JSTART = 1
         JCSTART = 1
      ENDIF
      
      J=JSTART

      I=ISTART
      
      IF ( (I.eq.2) .and. (J.eq.2) )  THEN 
         Q(I,J)=Q(I,J)+QC(ICSTART,JCSTART)
      ENDIF

      IF (J.eq.2) THEN
         DO IC=ICSTART+1,IICF1
            I=I+2
            Q(I,J)=Q(I,J)+QC(IC,JCSTART)
            A=CI(IC,JCSTART,LR)*QC(IC,JCSTART)
     &           +CI(IC,JCSTART,LL)*QC(IC-1,JCSTART)
            Q(I-1,J)=Q(I-1,J)+A +RES(I-1,J)
         ENDDO
      ENDIF
      
      DO JC=JCSTART+1,JJCF1
         J=J+2
         I=ISTART
         IF (I.gt.1) THEN
            Q(I,J)=Q(I,J)+QC(ICSTART,JC)
            AQ=CI(ICSTART,JC,LA)*QC(ICSTART,JC)
     &           +CI(ICSTART,JC,LB)*QC(ICSTART,JC-1)
            Q(ISTART,J-1)=Q(ISTART,J-1)+AQ +RES(ISTART,J-1)
         ENDIF
         DO IC=ICSTART+1,IICF1
            I=I+2
            Q(I,J)=Q(I,J)+QC(IC,JC)
            A=CI(IC,JC,LR)*QC(IC,JC)+CI(IC,JC,LL)*QC(IC-1,JC)
            Q(I-1,J)=Q(I-1,J)+A +RES(I-1,J)
            AQ=CI(IC,JC,LA)*QC(IC,JC)+CI(IC,JC,LB)*QC(IC,JC-1)
            Q(I,J-1)=Q(I,J-1)+AQ +RES(I,J-1)
            A=CI(IC,JC,LSW)*QC(IC-1,JC-1)+CI(IC,JC,LNW)
     &           *QC(IC-1,JC)+CI(IC,JC,LNE)*QC(IC,JC)
     &           +CI(IC,JC,LSE)*QC(IC,JC-1)
            Q(I-1,J-1)=Q(I-1,J-1)+A+RES(I-1,J-1)
         ENDDO
      ENDDO


      ptrn = 1
      
      call MSG_tbdx_send(Q, MSG_Buffer, 
     &     iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Proc,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Ipr,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Index,KF)),
     &     ptrn, ierror)
      
      call MSG_tbdx_receive(Q, MSG_Buffer,
     &     iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Proc,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Ipr,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Index,KF)),
     &     ptrn, ierror)
      
      call MSG_tbdx_close(Q, MSG_Buffer,
     &     iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Proc,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Ipr,KF)),
     &     iWorkMSG(pMSG(ipL_MSG_Index,KF)),
     &     ptrn, ierror)

C ==========================================================================

      RETURN
      END
