      SUBROUTINE BMG2_SER_PerSymStd_residual(
     &                KF, KC, IMULT, SO, SOC, QF, QFC,
     &                Q, QC, CI, SOR, IIC, JJC, IIF, JJF,
     &                TAUERR, NStncl, NSORv, JPN 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C 
C     BMG2_SER_PerSymStd_residual computes weighted averages of residuals on
C     the fine mesh to lay down on the coarse mesh.  (The static
C     residuals are computed in BMG2_SER_PerSymStd_relax.)  The weights
C     involve the transpose of the interpolation operator from the
C     coarse grid to the fine grid. In one mode BMG2_SER_PerSymStd_residual
C     computes a quantity to be used in estimating truncation
C     error. BMG2_SER_PerSymStd_residual also does some prliminary setup for
C     the interpolation routine BMG2_SER_PerSymStd_interp_add.
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
C     IMULT     If IMULT.ne.0, MGRCAP computes a weighted average of
C               residuals on the fine grid to lay down on the
C               coarse grid. If IMULT.eq.0 MGRCAP computes
C               a quantity TAUERR used in estimating truncation error.
C
C     SO        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     SOC       SO on grid KC.
C
C     QF        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     Q         Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     QC        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     CI        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     IPN       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     SOR       Refer to BOXMG.
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     QFC       QF on coarse grid.
C     TAUERR    TAUERR is a quantity used in estimating truncation
C               error.
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C ---------------------------
C     Includes
C
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'

C ----------------------------
C     Argument Declarations
C
      INTEGER IIC, IIF, IMULT, JJC, JJF, JPN, KC, KF, NSORv, NStncl
      REAL*RKIND  CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC), 
     &        QF(IIF,JJF), QFC(IIC,JJC), SO(IIF,JJF,NStncl), 
     &        SOC(IIC,JJC,5), SOR(IIF,JJF,NSORv), TAUERR

C ----------------------------
C     Local Declarations
C
      INTEGER I, IC, IIC1, IIC2, IIF1, IIF2, IPN, 
     &        J, JC, JJC1, JJC2, JJF1, JJF2
      REAL*RKIND  XITAU, TEMP

C ==========================================================================

      IPN=IABS(JPN)
      XITAU=IMULT
      IIC1=IIC-1
      JJC1=JJC-1
      IIF1=IIF-1
      JJF1=JJF-1
      IIC2=IIC-2
      JJC2=JJC-2
      IIF2=IIF-2
      JJF2=JJF-2

      TAUERR=RZERO

      IF (IMULT.EQ.0) THEN
         DO J=2,JJF1
            DO I=2,IIF1
               SOR(I,J,MTOT)=QF(I,J)
            ENDDO
         ENDDO
      ENDIF
C
C   compute weighted averages of residuals on fine mesh to lay
C   down on coarse mesh. the weights are i(transpose) where i is the
C   the interpolation operator from grid kc to grid kf.
C
      IF(IPN.EQ.1.OR.IPN.EQ.3) THEN
         DO I=1,IIF
            SOR(I,1,MTOT)=SOR(I,JJF1,MTOT)
            SOR(I,JJF,MTOT)=SOR(I,2,MTOT)
         ENDDO
      ENDIF
      IF(IPN.EQ.2.OR.IPN.EQ.3) THEN
         DO J=1,JJF
            SOR(1,J,MTOT)=SOR(IIF1,J,MTOT)
            SOR(IIF,J,MTOT)=SOR(2,J,MTOT)
         ENDDO
      ENDIF
      J=0
      DO JC=2,JJC1
         J=J+2
         I=0
         DO IC=2,IIC1
            I=I+2
            QFC(IC,JC)=CI(IC-1,JC-1,LNE)*SOR(I-1,J-1,MTOT)+CI(IC,JC,LA)
     &           *SOR(I,J-1,MTOT)+CI(IC,JC-1,LNW)*SOR(I+1,J-1,MTOT)
     &           +CI(IC,JC,LR)*SOR(I-1,J,MTOT)+SOR(I,J,MTOT)
     &           +CI(IC+1,JC,LL)*SOR(I+1,J,MTOT)+CI(IC-1,JC,LSE)
     &           *SOR(I-1,J+1,MTOT)+CI(IC,JC+1,LB)*SOR(I,J+1,MTOT)
     &           +CI(IC,JC,LSW)*SOR(I+1,J+1,MTOT)
            
         ENDDO
      ENDDO
C
C   if imult=0, estimate truncation errors. otherwise, set up
C   quantities for mginadp.
C
      IF (IMULT.NE.0) THEN
         DO J=2,JJF1
            DO I=2,IIF1
               SOR(I,J,MTOT)=SOR(I,J,MTOT)/SO(I,J,KO)
            ENDDO
         ENDDO
      ELSE
         J=0
         DO JC=2,JJC1
            J=J+2
            I=0
            DO IC=2,IIC1
               I=I+2
               QC(IC,JC)=Q(I,J)
            ENDDO
         ENDDO
         IF(IPN.EQ.1.OR.IPN.EQ.3) THEN
            DO IC=1,IIC
               QC(IC,1)=QC(IC,JJC1)
               QC(IC,JJC)=QC(IC,2)
            ENDDO
         ENDIF
         IF(IPN.EQ.2.OR.IPN.EQ.3) THEN
            DO JC=1,JJC
               QC(IIC,JC)=QC(2,JC)
               QC(1,JC)=QC(IIC1,JC)
            ENDDO
         ENDIF
         DO JC=2,JJC1
            DO IC=2,IIC1
               TEMP=SOC(IC,JC,KW)*QC(IC-1,JC)+SOC(IC+1,JC,KW)
     &              *QC(IC+1,JC)+SOC(IC,JC,KS)*QC(IC,JC-1)
     &              +SOC(IC,JC+1,KS)*QC(IC,JC+1)+SOC(IC,JC,KSW)
     &              *QC(IC-1,JC-1)+SOC(IC+1,JC,KNW)*QC(IC+1,JC-1)
     &              +SOC(IC,JC+1,KNW)*QC(IC-1,JC+1)+SOC(IC+1,JC+1,KSW)
     &              *QC(IC+1,JC+1)-SOC(IC,JC,KO)*QC(IC,JC)
               QFC(IC,JC)=QFC(IC,JC)+TEMP
               TAUERR=TAUERR+QFC(IC,JC)**2
            ENDDO
         ENDDO
         TAUERR=SQRT(TAUERR)
      ENDIF
     
C ==========================================================================

      RETURN
      END
