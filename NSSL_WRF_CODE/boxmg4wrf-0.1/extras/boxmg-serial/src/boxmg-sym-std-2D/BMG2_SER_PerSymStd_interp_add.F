      SUBROUTINE BMG2_SER_PerSymStd_interp_add(
     &                KC, KF, Q, QC, SOR, CI,
     &                IIC, JJC, IIF, JJF, NSORv, JPN 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_PerSymStd_interp_add interpolates Q from the coarse mesh KC
C     to the fine mesh KF and adds result to Q on fine mesh.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     KC        KC is grid number for coarse grid.
C     KF        KF is grid number for fine grid.
C
C     IIC       Number of grid points in x direction on coarse grid,
C               including two fictitious points.
C     JJC       Number of grid points in y direction on coarse grid,
C               including two fictitious points.
C
C     IIF       Number of grid points in x direction on fine grid,
C               including two fictitious points.
C     JJF       Number of grid points in y direction on fine grid,
C               including two fictitious points.
C
C     QC        Q for coarse grid.
C
C     SOR       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     CI        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C     IPN       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     Q         Refer to BMG2_SER_SymStd_SOLVE_boxmg.
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
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  IIC, IIF, JJC, JJF, JPN, KC, KF, NSORv
      REAL*RKIND   CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC),
     &         SOR(IIF,JJF,NSORv) 

C ----------------------------
C     Local Declarations
C
      INTEGER  IC, IF, IICF, IICF1, IIC1, IIF1, JC, JF, JJCF, JJCF1, 
     &         JJC1, JJF1, IPN
      REAL*RKIND   A, AQ

C =========================================================================

C
C   interpolate answers from coarse to fine mesh and add
C   to answers on fine mesh.
C

      IPN=IABS(JPN)
      JJF1=JJF-1
      IIF1=IIF-1
      IIC1=IIC-1
      JJC1=JJC-1
      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      IICF1=IICF-1
      JJCF1=JJCF-1
      JF=2
      IF=2
      Q(2,JF)=Q(2,JF)+QC(2,2)
      DO 10 IC=3,IICF1
      IF=IF+2
      Q(IF,JF)=Q(IF,JF)+QC(IC,2)
      A=CI(IC,2,LR)*QC(IC,2)+CI(IC,2,LL)*QC(IC-1,2)
      Q(IF-1,JF)=Q(IF-1,JF)+A+SOR(IF-1,JF,MTOT)
   10 CONTINUE
      DO 20 JC=3,JJCF1
      JF=JF+2
      IF=2
      Q(2,JF)=Q(2,JF)+QC(2,JC)
      AQ=CI(2,JC,LA)*QC(2,JC)+CI(2,JC,LB)*QC(2,JC-1)
      Q(2,JF-1)=Q(2,JF-1)+AQ+SOR(2,JF-1,MTOT)
      DO 20 IC=3,IICF1
      IF=IF+2
      Q(IF,JF)=Q(IF,JF)+QC(IC,JC)
      A=CI(IC,JC,LR)*QC(IC,JC)+CI(IC,JC,LL)*QC(IC-1,JC)
      Q(IF-1,JF)=Q(IF-1,JF)+A+SOR(IF-1,JF,MTOT)
      AQ=CI(IC,JC,LA)*QC(IC,JC)+CI(IC,JC,LB)*QC(IC,JC-1)
      Q(IF,JF-1)=Q(IF,JF-1)+AQ+SOR(IF,JF-1,MTOT)
      A=CI(IC-1,JC-1,LSW)*QC(IC-1,JC-1)+CI(IC-1,JC-1,LNW)*QC(IC-1,JC)
     1 +CI(IC-1,JC-1,LNE)*QC(IC,JC)+CI(IC-1,JC-1,LSE)*QC(IC,JC-1)
      Q(IF-1,JF-1)=Q(IF-1,JF-1)+A+SOR(IF-1,JF-1,MTOT)
   20 CONTINUE
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 30
      DO 21 IF=1,IIF
      Q(IF,1)=Q(IF,JJF1)
      Q(IF,JJF)=Q(IF,2)
   21 CONTINUE
   30 IF(IPN.NE.2.AND.IPN.NE.3)RETURN
      DO 31 JF=1,JJF
      Q(1,JF)=Q(IIF1,JF)
      Q(IIF,JF)=Q(2,JF)
   31 CONTINUE

C =========================================================================

      RETURN
      END
