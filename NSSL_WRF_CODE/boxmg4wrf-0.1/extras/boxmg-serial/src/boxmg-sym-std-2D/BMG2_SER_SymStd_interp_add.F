      SUBROUTINE BMG2_SER_SymStd_interp_add(
     &                       KC, KF, Q ,QC, RES, SO, CI,
     &                       IIC, JJC, IIF, JJF, NStncl, IBC 
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_interp_add.f interpolates Q from the coarse mesh, KC,
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
      INCLUDE 'BMG_SER_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER IBC, IIC, IIF, JJC, JJF, KC, KF, NStncl
      REAL*RKIND  CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC), 
     &        SO(IIF,JJF,NStncl), RES(IIF,JJF)

C ----------------------------
C     Local Declarations
C
      INTEGER IC, I, IICF, IICF1, IIC1, IIF1,
     &        JC, J, JJCF, JJCF1, JJC1, JJF1
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

      J=2
      I=2
      Q(2,J)=Q(2,J)+QC(2,2)
      DO IC=3,IICF1
         I=I+2
         Q(I,J)=Q(I,J)+QC(IC,2)
         A=CI(IC,2,LR)*QC(IC,2)+CI(IC,2,LL)*QC(IC-1,2)
         Q(I-1,J)=Q(I-1,J)+A+RES(I-1,J)
      ENDDO

      DO JC=3,JJCF1
         J=J+2
         I=2
         Q(2,J)=Q(2,J)+QC(2,JC)
         AQ=CI(2,JC,LA)*QC(2,JC)+CI(2,JC,LB)*QC(2,JC-1)
         Q(2,J-1)=Q(2,J-1)+AQ+RES(2,J-1)
         DO IC=3,IICF1
            I=I+2
            Q(I,J)=Q(I,J)+QC(IC,JC)
            A=CI(IC,JC,LR)*QC(IC,JC)+CI(IC,JC,LL)*QC(IC-1,JC)
            Q(I-1,J)=Q(I-1,J)+A+RES(I-1,J)
            AQ=CI(IC,JC,LA)*QC(IC,JC)+CI(IC,JC,LB)*QC(IC,JC-1)
            Q(I,J-1)=Q(I,J-1)+AQ+RES(I,J-1)
            A=CI(IC,JC,LSW)*QC(IC-1,JC-1)+CI(IC,JC,LNW)
     &           *QC(IC-1,JC)+CI(IC,JC,LNE)*QC(IC,JC)
     &           +CI(IC,JC,LSE)*QC(IC,JC-1)
            Q(I-1,J-1)=Q(I-1,J-1)+A+RES(I-1,J-1)
         ENDDO
      ENDDO

C ---------------------------------------
C     Periodicity copying:
C ---------------------------------------

      IF ( ABS(IBC).EQ.1 .OR. ABS(IBC).EQ.3 ) THEN 
         DO i=1, IIF
            Q(i,1)   = Q(i,JJF1)
            Q(i,JJF) = Q(i,2)
         ENDDO
      ENDIF

      IF( ABS(IBC).EQ.2 .OR. ABS(IBC).EQ.3 ) THEN
         DO j=1, JJF
            Q(1,j)   = Q(IIF1,j)
            Q(IIF,j) = Q(2,j)
         ENDDO
      ENDIF

C ==========================================================================

C =======================

      RETURN
      END
