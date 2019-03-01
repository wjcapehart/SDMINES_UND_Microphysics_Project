      SUBROUTINE BMG2_SER_PerSymStd_relax(
     &                K, SO, QF, Q, SOR, II, JJ, ERR_L2,
     &                ICOEF, IFD, NStncl, NSORv, IRELAX, JPN 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_PerSymStd_relax performs point or line relaxation on the fine
C     grid depending on the value of irelax.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     K         the index of the current grid
C
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C
C     SO        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     QF        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     SOR       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     ICOEF     Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     IFD       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     IRELAX    Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     IPN       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     Q         Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     ERR_L2    ERR_L2 is the l2 error of the resiuals.
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
      INTEGER  ICOEF, IFD, II, IRELAX, JJ, JPN, K, NSORv, NStncl
      REAL*RKIND   ERR_L2, Q(II,JJ), QF(II,JJ),
     &         SO(II,JJ,NStncl), SOR(II,JJ,NSORv)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, IBEG, IEND, I1, I2, IPN, J, JO, J1, J2, JBEG, JEND

C ==========================================================================

      IPN=IABS(JPN)
      J1=JJ-1
      I1=II-1
      J2=JJ-2
      I2=II-2
      ERR_L2=RZERO
      GO TO (10,60,130,60), IRELAX
   10 IF (K.GE.ICOEF.AND.IFD.EQ.1) GO TO 35
C
C   9 pt. operator
C   Relax red, then black, then green, then yellow points.
C
      IF(IPN.EQ.1.OR.IPN.EQ.2.OR.IPN.EQ.3) GO TO 11

      DO 20 JBEG=2,3

         JEND=2*((J1-JBEG)/2)+JBEG

         DO 19 J=JBEG,JEND,2
            DO 15 IBEG=2,3

               IEND=2*((I1-IBEG)/2)+IBEG

               DO 15 I=IBEG,IEND,2

                  Q(I,J) = ( QF(I,J) 
     &                 + SO(I,J,KW)*Q(I-1,J)
     &                 + SO(I+1,J,KW)*Q(I+1,J)
     &                 + SO(I,J,KS)*Q(I,J-1)
     &                 + SO(I,J+1,KS)*Q(I,J+1)
     &                 + SO(I,J,KSW)*Q(I-1,J-1)
     &                 + SO(I+1,J,KNW)*Q(I+1,J-1)
     &                 + SO(I,J+1,KNW)*Q(I-1,J+1)
     &                 + SO(I+1,J+1,KSW)*Q(I+1,J+1)
     &                 )*SOR(I,J,MSOR)

 15         CONTINUE
 19      CONTINUE
 20   CONTINUE

      IF(IPN.NE.1.AND.IPN.NE.2.AND.IPN.NE.3)GO TO 28
C
C   9 pt. operator, periodic boundary conditions
C   Relax red, then black, then green, then yellow points.
C
   11 DO 23 JBEG=2,3
      JEND=2*((J1-JBEG)/2)+JBEG
      DO 24 J=JBEG,JEND,2
      DO 26 IBEG=2,3
      IEND=2*((I1-IBEG)/2)+IBEG
Cdir$ ivdep
      DO 26 I=IBEG,IEND,2
      Q(I,J)=(QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+
     & SO(I+1,J,KNW)*Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)+
     & SO(I+1,J+1,KSW)*Q(I+1,J+1))*SOR(I,J,MSOR)
   26 CONTINUE
      Q(1,J)=Q(I1,J)
      Q(II,J)=Q(2,J)
      IF(J.NE.2.AND.J.NE.J1)GO TO 22
      DO 21 I=1,II
      Q(I,1)=Q(I,J1)
      Q(I,JJ)=Q(I,2)
   21 CONTINUE
   22 CONTINUE
   24 CONTINUE
   23 CONTINUE
C
C   Compute residuals.
C
   28 DO 30 J=2,J1
      DO 25 I=2,I1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+
     & SO(I+1,J,KNW)*Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)+
     & SO(I+1,J+1,KSW)*Q(I+1,J+1)-SO(I,J,KO)*Q(I,J)
   25 CONTINUE
      DO 27 I=2,I1
      ERR_L2=ERR_L2+SOR(I,J,MTOT)**2
   27 CONTINUE
   30 CONTINUE
      ERR_L2=SQRT(ERR_L2) 
      RETURN
C
C   5 pt. operator
C   Relax red points, then black points.
C
   35 IF(IPN.EQ.1.OR.IPN.EQ.2.OR.IPN.EQ.3)GO TO 39
      DO 55 JO=2,3
      DO 54 J=2,J1
      IBEG=MOD(J+JO,2)+2
      IEND=2*((I1-IBEG)/2)+IBEG
Cdir$ ivdep
      DO 50 I=IBEG,IEND,2
      Q(I,J)=(QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1))*SOR(I,J,MSOR)
      SOR(I,J,MTOT)=RZERO
   50 CONTINUE
   54 CONTINUE
   55 CONTINUE
      IF(IPN.NE.1.AND.IPN.NE.2.AND.IPN.NE.3)GO TO 56
C
C   5 pt. operator, periodic boundary conditions
C   Relax red points, then black points.
C
   39 DO 36 JO=2,3
      DO 37 J=2,J1
      IBEG=MOD(J+JO,2)+2
      IEND=2*((I1-IBEG)/2)+IBEG
Cdir$ ivdep
      DO 38 I=IBEG,IEND,2

      Q(I,J)=(QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1))*SOR(I,J,MSOR)
      SOR(I,J,MTOT)=RZERO
   38 CONTINUE
      Q(1,J)=Q(I1,J)
      Q(II,J)=Q(2,J)
      IF(J.NE.2.AND.J.NE.J1)GO TO 52
      DO 53 I=1,II
      Q(I,1)=Q(I,J1)
      Q(I,JJ)=Q(I,2)
   53 CONTINUE
   52 CONTINUE
   37 CONTINUE
   36 CONTINUE
C
C   Compute residuals for red points.
C
   56 DO 51 J=2,J1
      IBEG=MOD(J+2,2)+2
      IEND=2*((I1-IBEG)/2)+IBEG
      DO 57 I=IBEG,IEND,2
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)-SO(I,J,KO)*Q(I,J)
   57 CONTINUE
      DO 58 I=IBEG,IEND,2
      ERR_L2=ERR_L2+SOR(I,J,MTOT)**2
   58 CONTINUE
   51 CONTINUE
      ERR_L2=SQRT(ERR_L2)
      RETURN
   60 IF(IRELAX.EQ.4)GO TO 62
      DO 63 J=1,JJ
      DO 63 I=2,I1
      SOR(I,J,MTOT)=Q(I,J)
   63 CONTINUE
   62 CONTINUE
      IF(K.GE.ICOEF.AND.IFD.EQ.1)GO TO 90
C
C   9 pt. operator
C   Relax red lines, then black lines.
C
      DO 80 JBEG=3,2,-1
      JEND=2*((J1-JBEG)/2)+JBEG
      DO 65 I=2,I1
Cdir$ ivdep
      DO 65 J=JBEG,JEND,2
      Q(I,J)=QF(I,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)+
     & SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)*Q(I+1,J-1)+
     & SO(I,J+1,KNW)*Q(I-1,J+1)+SO(I+1,J+1,KSW)*Q(I+1,J+1)
   65 CONTINUE
      DO 70 I=2,I1
      DO 70 J=JBEG,JEND,2
      Q(I,J)=Q(I,J)+SO(I,J,KW)*SOR(I-1,J,MSOR)*Q(I-1,J)
   70 CONTINUE
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 74
      DO 71 J=JBEG,JEND,2
      QF(1,J)=SO(II,J,KW)
      QF(II,J)=SO(2,J,KW)*SOR(2,J,MSOR)
      Q(I1,J)=Q(I1,J)+QF(1,J)*SOR(2,J,MSOR)*Q(2,J)
   71 CONTINUE
      DO 72 I=3,I2
      DO 72 J=JBEG,JEND,2
      QF(1,J)=QF(1,J)*SO(I,J,KW)*SOR(I-1,J,MSOR)
      Q(I1,J)=Q(I1,J)+QF(1,J)*SOR(I,J,MSOR)*Q(I,J)
      QF(II,J)=QF(II,J)*SO(I,J,KW)*SOR(I,J,MSOR)
   72 CONTINUE
      GO TO 76
   74 DO 75 I=2,I1
      DO 75 J=JBEG,JEND,2
      Q(I1-I+2,J)=SOR(I1-I+2,J,MSOR)*(Q(I1-I+2,J)+
     & SO(I1-I+3,J,KW)*Q(I1-I+3,J))
   75 CONTINUE
      GO TO 82
   76 DO 77 J=JBEG,JEND,2
      Q(I1,J)=SOR(I1,J,MSOR)*Q(I1,J)
      Q(I2,J)=SOR(I2,J,MSOR)*(Q(I2,J)+SO(I1,J,KW)*Q(I1,J))+QF(II,J)
     &*Q(I1,J)
   77 CONTINUE
      DO 78 I=4,I1
      DO 78 J=JBEG,JEND,2
      QF(II,J)=QF(II,J)/(SO(I1-I+3,J,KW)*SOR(I1-I+3,J,MSOR))
      Q(I1-I+2,J)=SOR(I1-I+2,J,MSOR)*(Q(I1-I+2,J)+
     & SO(I1-I+3,J,KW)*Q(I1-I+3,J))+QF(II,J)*Q(I1,J)
   78 CONTINUE
   82 IF(IPN.NE.1.AND.IPN.NE.3)GO TO 81
      DO 79 I=1,II
      Q(I,1)=Q(I,J1)
      Q(I,JJ)=Q(I,2)
   79 CONTINUE
   81 IF(IPN.NE.2.AND.IPN.NE.3)GO TO 80
      DO 83 J=1,JJ
      Q(1,J)=Q(I1,J)
      Q(II,J)=Q(2,J)
      QF(1,J)=RZERO
      QF(II,J)=RZERO
   83 CONTINUE
   80 CONTINUE
      IF(IRELAX.EQ.4)GO TO 130
      JBEG=3
      JEND=2*((J1-3)/2)+3
      DO 85 I=2,I1
Cdir$ ivdep
      DO 85 J=JBEG,JEND,2
      SOR(I,J,MTOT)=SO(I,J,KS)*(Q(I,J-1)-SOR(I,J-1,MTOT))
     &+SO(I,J+1,KS)*(Q(I,J+1)-SOR(I,J+1,MTOT))+SO(I,J,KSW)*
     &(Q(I-1,J-1)-SOR(I-1,J-1,MTOT))+SO(I+1,J,KNW)*
     &(Q(I+1,J-1)-SOR(I+1,J-1,MTOT))+SO(I,J+1,KNW)*
     &(Q(I-1,J+1)-SOR(I-1,J+1,MTOT))+SO(I+1,J+1,KSW)*
     &(Q(I+1,J+1)-SOR(I+1,J+1,MTOT))
   85 CONTINUE
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 120
      J=2
      DO 86 I=2,I1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+
     & SO(I+1,J,KNW)*Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)+
     & SO(I+1,J+1,KSW)*Q(I+1,J+1)-SO(I,J,KO)*Q(I,J)
   86 CONTINUE
      J=J1
      DO 87 I=2,I1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+
     & SO(I+1,J,KNW)*Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)+
     & SO(I+1,J+1,KSW)*Q(I+1,J+1)-SO(I,J,KO)*Q(I,J)
   87 CONTINUE
      GO TO 120
   90 DO 110 JBEG=3,2,-1
      JEND=2*((J1-JBEG)/2)+JBEG
      DO 95 I=2,I1
Cdir$ ivdep
      DO 95 J=JBEG,JEND,2
      Q(I,J)=QF(I,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)
   95 CONTINUE
      DO 100 I=2,I1
      DO 100 J=JBEG,JEND,2
      Q(I,J)=Q(I,J)+SO(I,J,KW)*SOR(I-1,J,MSOR)*Q(I-1,J)
  100 CONTINUE
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 104
      DO 106 J=JBEG,JEND,2
      QF(1,J)=SO(II,J,KW)
      QF(II,J)=SO(2,J,KW)*SOR(2,J,MSOR)
      Q(I1,J)=Q(I1,J)+QF(1,J)*SOR(2,J,MSOR)*Q(2,J)
  106 CONTINUE
      DO 107 I=3,I2
      DO 107 J=JBEG,JEND,2
      QF(1,J)=QF(1,J)*SO(I,J,KW)*SOR(I-1,J,MSOR)
      Q(I1,J)=Q(I1,J)+QF(1,J)*SOR(I,J,MSOR)*Q(I,J)
      QF(II,J)=QF(II,J)*SO(I,J,KW)*SOR(I,J,MSOR)
  107 CONTINUE
      GO TO 101
  104 DO 105 I=2,I1
      DO 105 J=JBEG,JEND,2
      Q(I1-I+2,J)=SOR(I1-I+2,J,MSOR)*(Q(I1-I+2,J)+
     & SO(I1-I+3,J,KW)*Q(I1-I+3,J))
  105 CONTINUE
      GO TO 112
  101 DO 108 J=JBEG,JEND,2
      Q(I1,J)=SOR(I1,J,MSOR)*Q(I1,J)
      Q(I2,J)=SOR(I2,J,MSOR)*(Q(I2,J)+SO(I1,J,KW)*Q(I1,J))+
     & QF(II,J)*Q(I1,J)
  108 CONTINUE
      DO 111 I=4,I1
      DO 111 J=JBEG,JEND,2
      QF(II,J)=QF(II,J)/(SO(I1-I+3,J,KW)*SOR(I1-I+3,J,MSOR))
      Q(I1-I+2,J)=SOR(I1-I+2,J,MSOR)*(Q(I1-I+2,J)+
     & SO(I1-I+3,J,KW)*Q(I1-I+3,J))+QF(II,J)*Q(I1,J)
  111 CONTINUE
  112 IF(IPN.NE.1.AND.IPN.NE.3)GO TO 113
      DO 109 I=1,II
      Q(I,1)=Q(I,J1)
      Q(I,JJ)=Q(I,2)
  109 CONTINUE
  113 IF(IPN.NE.2.AND.IPN.NE.3)GO TO 110
      DO 114 J=1,JJ
      Q(1,J)=Q(I1,J)
      Q(II,J)=Q(2,J)
      QF(1,J)=RZERO
      QF(II,J)=RZERO
  114 CONTINUE
  110 CONTINUE
      IF(IRELAX.EQ.4)GO TO 130
      JBEG=3
      JEND=2*((J1-3)/2)+3
      DO 115 I=2,I1
Cdir$ ivdep
      DO 115 J=JBEG,JEND,2
      SOR(I,J,MTOT)=SO(I,J,KS)*(Q(I,J-1)-SOR(I,J-1,MTOT))
     &+SO(I,J+1,KS)*(Q(I,J+1)-SOR(I,J+1,MTOT))
  115 CONTINUE
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 120
      J=2
      DO 116 I=2,I1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)-SO(I,J,KO)*Q(I,J)
  116 CONTINUE
      J=J1
      DO 117 I=2,I1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)-SO(I,J,KO)*Q(I,J)
  117 CONTINUE
  120 DO 121 I=2,I1
      DO 121 J=JBEG,JEND,2
      ERR_L2=ERR_L2+SOR(I,J,MTOT)**2
  121 CONTINUE
      JBEG=2
      JEND=2*((J1-2)/2)+2
      IF(IPN.EQ.1.OR.IPN.EQ.3)JBEG=4
      IF((IPN.EQ.1.OR.IPN.EQ.3).AND.JEND.EQ.J1)JEND=JEND-2
      DO 122 I=2,I1
      DO 122 J=JBEG,JEND,2
      SOR(I,J,MTOT)=RZERO
  122 CONTINUE
      ERR_L2=SQRT(ERR_L2)
C
C   If irelax.eq.2, put zeroes into residuals for black lines
C   (5 pt. or 9 pt. operator.).
C
      IF (IRELAX.NE.4) RETURN
  130 DO 132 J=1,JJ
      DO 132 I=2,I1
      SOR(I,J,MTOT)=Q(I,J)
  132 CONTINUE
      IF(K.GE.ICOEF.AND.IFD.EQ.1)GO TO 160
C
C   9 pt. operator
C   Relax red lines, then black lines.
C
      DO 150 IBEG=3,2,-1
      IEND=2*((I1-IBEG)/2)+IBEG
      DO 135 J=2,J1
Cdir$ ivdep
      DO 135 I=IBEG,IEND,2
      Q(I,J)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)*Q(I+1,J-1)+SO(I,J+1,KNW)*
     & Q(I-1,J+1)+SO(I+1,J+1,KSW)*Q(I+1,J+1)
  135 CONTINUE
      DO 140 J=2,J1
      DO 140 I=IBEG,IEND,2
      Q(I,J)=Q(I,J)+SO(I,J,KS)*SOR(I,J-1,MSOS)*Q(I,J-1)
  140 CONTINUE
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 144
      DO 141 I=IBEG,IEND,2
      QF(I,1)=SO(I,JJ,KS)
      QF(I,JJ)=SO(I,2,KS)*SOR(I,2,MSOS)
      Q(I,J1)=Q(I,J1)+QF(I,1)*SOR(I,2,MSOS)*Q(I,2)
  141 CONTINUE
      DO 142 J=3,J2
      DO 142 I=IBEG,IEND,2
      QF(I,1)=QF(I,1)*SO(I,J,KS)*SOR(I,J-1,MSOS)
      Q(I,J1)=Q(I,J1)+QF(I,1)*SOR(I,J,MSOS)*Q(I,J)
      QF(I,JJ)=QF(I,JJ)*SO(I,J,KS)*SOR(I,J,MSOS)
  142 CONTINUE
      GO TO 148
  144 DO 145 J=2,J1
      DO 145 I=IBEG,IEND,2
      Q(I,J1-J+2)=SOR(I,J1-J+2,MSOS)*(Q(I,J1-J+2)+
     & SO(I,J1-J+3,KS)*Q(I,J1-J+3))
  145 CONTINUE
      GO TO 149
  148 DO 146 I=IBEG,IEND,2
      Q(I,J1)=SOR(I,J1,MSOS)*Q(I,J1)
      Q(I,J2)=SOR(I,J2,MSOS)*(Q(I,J2)+SO(I,J1,KS)*Q(I,J1))+
     & QF(I,JJ)*Q(I,J1)
  146 CONTINUE
      DO 147 J=4,J1
      DO 147 I=IBEG,IEND,2
      QF(I,JJ)=QF(I,JJ)/(SO(I,J1-J+3,KS)*SOR(I,J1-J+3,MSOS))
      Q(I,J1-J+2)=SOR(I,J1-J+2,MSOS)*(Q(I,J1-J+2)+
     & SO(I,J1-J+3,KS)*Q(I,J1-J+3))+QF(I,JJ)*Q(I,J1)
  147 CONTINUE
  149 IF(IPN.NE.1.AND.IPN.NE.3)GO TO 154
      DO 151 I=1,II
      Q(I,1)=Q(I,J1)
      Q(I,JJ)=Q(I,2)
      QF(I,1)=RZERO
      QF(I,JJ)=RZERO
  151 CONTINUE
  154 IF(IPN.NE.2.AND.IPN.NE.3)GO TO 150
      DO 153 J=1,JJ
      Q(1,J)=Q(I1,J)
      Q(II,J)=Q(2,J)
  153 CONTINUE
  150 CONTINUE
C
C   Compute residuals for red lines.
C
      IBEG=3
      IEND=2*((I1-3)/2)+3
      DO 155 J=2,J1
Cdir$ ivdep
      DO 155 I=IBEG,IEND,2
      SOR(I,J,MTOT)=SO(I,J,KW)*(Q(I-1,J)-SOR(I-1,J,MTOT))
     &+SO(I+1,J,KW)*(Q(I+1,J)-SOR(I+1,J,MTOT))+SO(I,J,KSW)*
     &(Q(I-1,J-1)-SOR(I-1,J-1,MTOT))+SO(I+1,J,KNW)*
     &(Q(I+1,J-1)-SOR(I+1,J-1,MTOT))+SO(I,J+1,KNW)*
     &(Q(I-1,J+1)-SOR(I-1,J+1,MTOT))+SO(I+1,J+1,KSW)*
     &(Q(I+1,J+1)-SOR(I+1,J+1,MTOT))
  155 CONTINUE
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 185
      I=2
      DO 156 J=2,J1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+
     & SO(I+1,J,KNW)*Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)+
     & SO(I+1,J+1,KSW)*Q(I+1,J+1)-SO(I,J,KO)*Q(I,J)
  156 CONTINUE
      I=I1
      DO 157 J=2,J1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+
     & SO(I+1,J,KNW)*Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)+
     & SO(I+1,J+1,KSW)*Q(I+1,J+1)-SO(I,J,KO)*Q(I,J)
  157 CONTINUE
      GO TO 185
C
C   5 pt. operator
C   Relax red lines, then black lines.
C
  160 DO 175 IBEG=3,2,-1
      IEND=2*((I1-IBEG)/2)+IBEG
      DO 165 J=2,J1
Cdir$ ivdep
      DO 165 I=IBEG,IEND,2
      Q(I,J)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)
  165 CONTINUE
      DO 168 J=2,J1
      DO 168 I=IBEG,IEND,2
      Q(I,J)=Q(I,J)+SO(I,J,KS)*SOR(I,J-1,MSOS)*Q(I,J-1)
  168 CONTINUE
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 164
      DO 166 I=IBEG,IEND,2
      QF(I,1)=SO(I,JJ,KS)
      QF(I,JJ)=SO(I,2,KS)*SOR(I,2,MSOS)
      Q(I,J1)=Q(I,J1)+QF(I,1)*SOR(I,2,MSOS)*Q(I,2)
  166 CONTINUE
      DO 167 J=3,J2
      DO 167 I=IBEG,IEND,2
      QF(I,1)=QF(I,1)*SO(I,J,KS)*SOR(I,J-1,MSOS)
      Q(I,J1)=Q(I,J1)+QF(I,1)*SOR(I,J,MSOS)*Q(I,J)
      QF(I,JJ)=QF(I,JJ)*SO(I,J,KS)*SOR(I,J,MSOS)
  167 CONTINUE
      GO TO 174
  164 DO 170 J=2,J1
      DO 170 I=IBEG,IEND,2
      Q(I,J1-J+2)=SOR(I,J1-J+2,MSOS)*(Q(I,J1-J+2)+SO(I,J1-J+3,KS)*
     & Q(I,J1-J+3))
  170 CONTINUE
      GO TO 173
  174 DO 171 I=IBEG,IEND,2
      Q(I,J1)=SOR(I,J1,MSOS)*Q(I,J1)
      Q(I,J2)=SOR(I,J2,MSOS)*(Q(I,J2)+SO(I,J1,KS)*Q(I,J1))+QF(I,JJ)
     &*Q(I,J1)
  171 CONTINUE
      DO 172 J=4,J1
      DO 172 I=IBEG,IEND,2
      QF(I,JJ)=QF(I,JJ)/(SO(I,J1-J+3,KS)*SOR(I,J1-J+3,MSOS))
      Q(I,J1-J+2)=SOR(I,J1-J+2,MSOS)*(Q(I,J1-J+2)+SO(I,J1-J+3,KS)
     &*Q(I,J1-J+3))+QF(I,JJ)*Q(I,J1)
  172 CONTINUE
  173 IF(IPN.NE.1.AND.IPN.NE.3)GO TO 178
      DO 176 I=1,II
      Q(I,1)=Q(I,J1)
      Q(I,JJ)=Q(I,2)
      QF(I,1)=RZERO
      QF(I,JJ)=RZERO
  176 CONTINUE
  178 IF(IPN.NE.2.AND.IPN.NE.3)GO TO 175
      DO 179 J=1,JJ
      Q(1,J)=Q(I1,J)
      Q(II,J)=Q(2,J)
  179 CONTINUE
  175 CONTINUE
C   Compute residuals for red lines.
      IBEG=3
      IEND=2*((I1-3)/2)+3
      DO 180 J=2,J1
Cdir$ ivdep
      DO 180 I=IBEG,IEND,2
      SOR(I,J,MTOT)=SO(I,J,KW)*(Q(I-1,J)-SOR(I-1,J,MTOT))
     &+SO(I+1,J,KW)*(Q(I+1,J)-SOR(I+1,J,MTOT))
  180 CONTINUE
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 185
      I=2
      DO 181 J=2,J1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)-SO(I,J,KO)*Q(I,J)
  181 CONTINUE
      I=I1
      DO 182 J=2,J1
      SOR(I,J,MTOT)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)*Q(I+1,J)+
     & SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)*Q(I,J+1)-SO(I,J,KO)*Q(I,J)
  182 CONTINUE
  185 DO 190 J=2,J1
      DO 190 I=IBEG,IEND,2
      ERR_L2=ERR_L2+SOR(I,J,MTOT)**2
  190 CONTINUE
C
C   Put zeroes into residuals for black lines (5 pt. and 9 pt.
C
      IBEG=2
      IEND=2*((I1-2)/2)+2
      IF(IPN.EQ.2.OR.IPN.EQ.3)IBEG=4
      IF((IPN.EQ.2.OR.IPN.EQ.3).AND.IEND.EQ.I1)IEND=IEND-2
      DO 195 J=2,J1
      DO 195 I=IBEG,IEND,2
      SOR(I,J,MTOT)=RZERO
  195 CONTINUE
      ERR_L2=SQRT(ERR_L2)

C ==========================================================================

      RETURN
      END
