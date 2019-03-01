      SUBROUTINE BMG2_SER_PerSymStd_SETUP_cofp(
     &                KF, KC, SO, SOC, QF, QFC, SOR, SORC, CI,
     &                IIF, JJF, IIC, JJC, ICOEF, IFD, NStncl,
     &                NSORv, JPN, IRELAX, ISTRT, ISKIP 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_PerSymStd_SETUP_cofp computes the interpolation operator from
C     the coarse grid to the fine grid. It also calculates the
C     difference operator on the coarse grid and computes the right hand
C     side for the coarse grid if required. Coef will also compute lu
C     decompositions if line relaxation is specified.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     KF        fine grid number
C     KC        coarse grid number
C
C     IIF       Number of grid points in x direction on fine grid,
C               including two fictitious points.
C     JJF       Number of grid points in y direction on fine grid,
C               including two fictitious points.
C
C     IIC       Number of grid points in x direction on coarse grid,
C               including two fictitious points.
C     JJC       Number of grid points in y direction on coarse grid,
C               including two fictitious points.
C
C     ICOEF     Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     IFD       Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     IRELAX    Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     ISTRT     Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     SO        Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     SOR       Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     CI        Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     ISKIP     Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     IPN       Refer to BMG2_SER_SymStd_SOLVE_boxmg
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
C     SOC       SO for coarse grid
C     QFC       QF for coarse grid
C     SOR       SOR for coarse grid
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================
     
      IMPLICIT NONE

C ----------------------------
C     Includes
C ----------------------------

#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'

C ---------------------------
C    Argument Declarations:
C ---------------------------

      INTEGER   ICOEF, IIC, IIF, IFD, IRELAX, ISTRT, ISKIP, 
     &          JJC, JJF, KC, KF, NSORv, NStncl

      REAL*RKIND    CI(IIC,JJC,8), QF(IIF,JJF), QFC(IIC,JJC),
     &          SO(IIF,JJF,NStncl), SOC(IIC,JJC,5), 
     &          SOR(IIF,JJF,NSORv), SORC(IIC,JJC,NSORv)

C --------------------------
C     Local Declarations:
C --------------------------

      INTEGER   IBEGC, IC, IF, IIC1, IICF, IICF1, IICF2, IIF1, IIF2, 
     &          INDX, IPN, JBEGC, JC, JF, JJC1, JJCF, JJCF1, JJCF2, 
     &          JJF1, JJF2, JPN
      REAL*RKIND A, B, CE, CEA, CENW, CFNW, CN, CNE, CNW, CO, COA,CONW,
     &          COSW, CS, CSA, CSE, CSEA, CSENW, CSNW, CSWA, CSSW, CSW,
     &          CSWSW, CW, CWA, CWSW, D1MACH, EP, ONE, SUM, S, ZEPS,
     &          ZERO

C
C***FIRST EXECUTABLE STATEMENT  MGCOFP
C
      ZERO = 0
      ONE  = 1
      ZEPS = D1MACH(3)
C
      IPN=IABS(JPN)
      IIC1=IIC-1
      JJC1=JJC-1
      IIF1=IIF-1
      JJF1=JJF-1
      IIF2=IIF-2
      JJF2=JJF-2
      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      IICF1=IICF-1
      JJCF1=JJCF-1
      IICF2=IICF-2
      JJCF2=JJCF-2
C
C   for finest grid compute reciprocal of so and store in
C   sor.
C
      DO 10 JF=2,JJF1
      DO 10 IF=2,IIF1
      SOR(IF,JF,MSOR)=ONE/SO(IF,JF,KO)
      SOR(IF,JF,MTOT)=QF(IF,JF)*SOR(IF,JF,MSOR)
   10 CONTINUE
      IF(ISKIP.EQ.2)GO TO 230
C
C   if kf difference operator is five point and kf.ge.icoef, go to 1600
C
      IF (IFD.EQ.1.AND.KF.GE.ICOEF) GO TO 120
C******************************
C   begin computation of i when kf difference operator is nine point
C
      JF=0
      DO 50 JC=2,JJC1
      JF=JF+2
      IF=2
      IBEGC=3
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 39
      IF=0
      IBEGC=2
   39 DO 40 IC=IBEGC,IICF1
      IF=IF+2
      A=SO(IF,JF,KW)+SO(IF,JF,KNW)+SO(IF,JF+1,KSW)
      B=SO(IF-1,JF,KW)+SO(IF-1,JF,KSW)+SO(IF-1,JF+1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF-1,JF,KS)+SO(IF-1,JF+1,KS)
      SUM=A+B+(SO(IF-1,JF,KO)-SUM)*MAX(SO(IF-1,JF,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF-1,JF,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
   40 CONTINUE
   50 CONTINUE
      IF(IIC.NE.IICF.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 45
      IF=3
      IC=IIC
      JF=0
      DO 43 JC=2,JJC1
      JF=JF+2
      A=SO(IF,JF,KW)+SO(IF,JF,KNW)+SO(IF,JF+1,KSW)
      B=SO(IF-1,JF,KW)+SO(IF-1,JF,KSW)+SO(IF-1,JF+1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF-1,JF,KS)+SO(IF-1,JF+1,KS)
      SUM=A+B+(SO(IF-1,JF,KO)-SUM)*MAX(SO(IF-1,JF,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF-1,JF,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
   43 CONTINUE
   45 IF(IPN.NE.1.AND.IPN.NE.3)GO TO 70
      IF(JJCF.EQ.JJC)GO TO 53
      DO 51 IC=2,IIC
      CI(IC,1,LL)=CI(IC,JJC1,LL)
      CI(IC,1,LR)=CI(IC,JJC1,LR)
      CI(IC,JJC,LL)=CI(IC,2,LL)
      CI(IC,JJC,LR)=CI(IC,2,LR)
   51 CONTINUE
      GO TO 70
   53 JF=JJF2
      JC=1
      INDX=0
      GO TO 55
   54 JF=3
      JC=JJC
   55 IBEGC=3
      IF=2
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 57
      IBEGC=2
      IF=0
   57 DO 56 IC=IBEGC,IICF1
      IF=IF+2
      A=SO(IF,JF,KW)+SO(IF,JF,KNW)+SO(IF,JF+1,KSW)
      B=SO(IF-1,JF,KW)+SO(IF-1,JF,KSW)+SO(IF-1,JF+1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF-1,JF,KS)+SO(IF-1,JF+1,KS)
      SUM=A+B+(SO(IF-1,JF,KO)-SUM)*MAX(SO(IF-1,JF,KO)-(ONE+EP)*SUM,ZERO)
     1/(ABS(SO(IF-1,JF,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
   56 CONTINUE
      IF(IICF.NE.IIC.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 58
      IF=3
      IC=IIC
      A=SO(IF,JF,KW)+SO(IF,JF,KNW)+SO(IF,JF+1,KSW)
      B=SO(IF-1,JF,KW)+SO(IF-1,JF,KSW)+SO(IF-1,JF+1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF-1,JF,KS)+SO(IF-1,JF+1,KS)
      SUM=A+B+(SO(IF-1,JF,KO)-SUM)*MAX(SO(IF-1,JF,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF-1,JF,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
   58 IF(INDX.NE.0)GO TO 70
      INDX=1
      GO TO 54
   70 JF=2
      JBEGC=3
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 75
      JF=0
      JBEGC=2
   75 DO 80 JC=JBEGC,JJCF1
      JF=JF+2
      IF=0
      DO 79 IC=2,IIC1
      IF=IF+2
      A=SO(IF,JF,KS)+SO(IF,JF,KNW)+SO(IF+1,JF,KSW)
      B=SO(IF,JF-1,KS)+SO(IF,JF-1,KSW)+SO(IF+1,JF-1,KNW)
      SUM=A+B+SO(IF,JF-1,KW)+SO(IF+1,JF-1,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+(SO(IF,JF-1,KO)-SUM)*MAX(SO(IF,JF-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF,JF-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LA)=A*SUM
      CI(IC,JC,LB)=B*SUM
   79 CONTINUE
      IF(IICF.NE.IIC.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 59
      INDX=0
      IF=IIF2
      IC=1
      GO TO 78
   77 IF=3
      IC=IIC
   78 A=SO(IF,JF,KS)+SO(IF,JF,KNW)+SO(IF+1,JF,KSW)
      B=SO(IF,JF-1,KS)+SO(IF,JF-1,KSW)+SO(IF+1,JF-1,KNW)
      SUM=A+B+SO(IF,JF-1,KW)+SO(IF+1,JF-1,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+(SO(IF,JF-1,KO)-SUM)*MAX(SO(IF,JF-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF,JF-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LA)=A*SUM
      CI(IC,JC,LB)=B*SUM
      IF(INDX.NE.0)GO TO 80
      INDX=1
      GO TO 77
   59 IF(IICF.EQ.IIC.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 80
      CI(1,JC,LA)=CI(IIC1,JC,LA)
      CI(1,JC,LB)=CI(IIC1,JC,LB)
      CI(IIC,JC,LA)=CI(2,JC,LA)
      CI(IIC,JC,LB)=CI(2,JC,LB)
   80 CONTINUE
      IF(JJC.NE.JJCF.OR.(IPN.NE.1.AND.IPN.NE.3))GO TO 84
      JF=3
      IF=0
      DO 81 IC=2,IIC1
      IF=IF+2
      A=SO(IF,JF,KS)+SO(IF,JF,KNW)+SO(IF+1,JF,KSW)
      B=SO(IF,JF-1,KS)+SO(IF,JF-1,KSW)+SO(IF+1,JF-1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF,JF-1,KW)+SO(IF+1,JF-1,KW)
      SUM=A+B+(SO(IF,JF-1,KO)-SUM)*MAX(SO(IF,JF-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF,JF-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JJC,LA)=A*SUM
      CI(IC,JJC,LB)=B*SUM
   81 CONTINUE
      IF(IICF.NE.IIC.OR.IPN.NE.3)GO TO 835
      INDX=0
      IF=IIF2
      IC=1
      GO TO 83
   82 IF=3
      IC=IIC
   83 A=SO(IF,JF,KS)+SO(IF,JF,KNW)+SO(IF+1,JF,KSW)
      B=SO(IF,JF-1,KS)+SO(IF,JF-1,KSW)+SO(IF+1,JF-1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF,JF-1,KW)+SO(IF+1,JF-1,KW)
      SUM=A+B+(SO(IF,JF-1,KO)-SUM)*MAX(SO(IF,JF-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF,JF-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JJC,LA)=A*SUM
      CI(IC,JJC,LB)=B*SUM
      IF(INDX.NE.0)GO TO 84
      INDX=1
      GO TO 82
 835  IF (IICF.EQ.IIC .OR. IPN.NE.3) GOTO 84
      CI(1,JJC,LA)=CI(IIC1,JJC,LA)
      CI(1,JJC,LB)=CI(IIC1,JJC,LB)
      CI(IIC,JJC,LA)=CI(2,JJC,LA)
      CI(IIC,JJC,LB)=CI(2,JJC,LB)      
   84 JF=0
      JBEGC=2
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 85
      JBEGC=1
      JF=-2
   85 DO 100 JC=JBEGC,JJCF2
      JF=JF+2
      IF=0
      IBEGC=2
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 89
      IBEGC=1
      IF=-2
   89 DO 90 IC=IBEGC,IICF2
      IF=IF+2
      SUM=SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KNW)+SO(IF+1,JF+2,KS)
     1+SO(IF+2,JF+2,KSW)+SO(IF+2,JF+1,KW)+SO(IF+2,JF+1,KNW)
     2+SO(IF+1,JF+1,KS)+SO(IF+1,JF+1,KSW)
      EP=MIN(ABS(SO(IF+1,JF+1,KSW)+SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KNW)),
     &       ABS(SO(IF+1,JF+2,KNW)+SO(IF+1,JF+2,KS)+SO(IF+2,JF+2,KSW)),
     &       ABS(SO(IF+2,JF+2,KSW)+SO(IF+2,JF+1,KW)+SO(IF+2,JF+1,KNW)),
     &       ABS(SO(IF+2,JF+1,KNW)+SO(IF+1,JF+1,KS)+SO(IF+1,JF+1,KSW)),
     &       ONE)
      SUM=SUM+(SO(IF+1,JF+1,KO)-SUM)*MAX(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*
     & CI(IC,JC+1,LB)+SO(IF+1,JF+1,KSW))*S
      CI(IC,JC,LSE)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LB)+SO(IF+2,JF+1,KNW))*S
      CI(IC,JC,LNW)=(SO(IF+1,JF+1,KW)*CI(IC,JC+1,LA)+SO(IF+1,JF+2,KS)*
     &CI(IC+1,JC+1,LL) +SO(IF+1,JF+2,KNW))*S
      CI(IC,JC,LNE)=(SO(IF+1,JF+2,KS)*CI(IC+1,JC+1,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LA)+SO(IF+2,JF+2,KSW))*S
   90 CONTINUE
      IF(IIC.NE.IICF.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 100
      IF=1
      IC=IIC1
      SUM=SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KNW)+SO(IF+1,JF+2,KS)
     &+SO(IF+2,JF+2,KSW)+SO(IF+2,JF+1,KW)+SO(IF+2,JF+1,KNW)
     &+SO(IF+1,JF+1,KS)+SO(IF+1,JF+1,KSW)
      EP=MIN(ABS(SO(IF+1,JF+1,KSW)+SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KNW)),
     &       ABS(SO(IF+1,JF+2,KNW)+SO(IF+1,JF+2,KS)+SO(IF+2,JF+2,KSW)),
     &       ABS(SO(IF+2,JF+2,KSW)+SO(IF+2,JF+1,KW)+SO(IF+2,JF+1,KNW)),
     &       ABS(SO(IF+2,JF+1,KNW)+SO(IF+1,JF+1,KS)+SO(IF+1,JF+1,KSW)),
     &       ONE)
      SUM=SUM+(SO(IF+1,JF+1,KO)-SUM)*MAX(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*
     & CI(IC,JC+1,LB)+SO(IF+1,JF+1,KSW))*S
      CI(IC,JC,LSE)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LB)+SO(IF+2,JF+1,KNW))*S
      CI(IC,JC,LNW)=(SO(IF+1,JF+1,KW)*CI(IC,JC+1,LA)+SO(IF+1,JF+2,KS)*
     &CI(IC+1,JC+1,LL) +SO(IF+1,JF+2,KNW))*S
      CI(IC,JC,LNE)=(SO(IF+1,JF+2,KS)*CI(IC+1,JC+1,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LA)+SO(IF+2,JF+2,KSW))*S
  100 CONTINUE
      IF(JJC.NE.JJCF .OR. (IPN.NE.1 .AND. IPN.NE.3))GO TO 105
      JF=1
      JC=JJC1
      IF=0
      IBEGC=2
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 103
      IBEGC=1
      IF=-2
  103 DO 101 IC=IBEGC,IICF2
      IF=IF+2
      SUM=SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KNW)+SO(IF+1,JF+2,KS)
     &+SO(IF+2,JF+2,KSW)+SO(IF+2,JF+1,KW)+SO(IF+2,JF+1,KNW)
     &+SO(IF+1,JF+1,KS)+SO(IF+1,JF+1,KSW)
      EP=MIN(ABS(SO(IF+1,JF+1,KSW)+SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KNW)),
     &       ABS(SO(IF+1,JF+2,KNW)+SO(IF+1,JF+2,KS)+SO(IF+2,JF+2,KSW)),
     &       ABS(SO(IF+2,JF+2,KSW)+SO(IF+2,JF+1,KW)+SO(IF+2,JF+1,KNW)),
     &       ABS(SO(IF+2,JF+1,KNW)+SO(IF+1,JF+1,KS)+SO(IF+1,JF+1,KSW)),
     &       ONE)
      SUM=SUM+(SO(IF+1,JF+1,KO)-SUM)*MAX(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*
     & CI(IC,JC+1,LB)+SO(IF+1,JF+1,KSW))*S
      CI(IC,JC,LSE)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LB)+SO(IF+2,JF+1,KNW))*S
      CI(IC,JC,LNW)=(SO(IF+1,JF+1,KW)*CI(IC,JC+1,LA)+SO(IF+1,JF+2,KS)*
     &CI(IC+1,JC+1,LL) +SO(IF+1,JF+2,KNW))*S
      CI(IC,JC,LNE)=(SO(IF+1,JF+2,KS)*CI(IC+1,JC+1,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LA)+SO(IF+2,JF+2,KSW))*S
  101 CONTINUE
      IF(IIC.NE.IICF.OR.IPN.NE.3)GO TO 105
      IF=1
      IC=IIC1
      SUM=SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KNW)+SO(IF+1,JF+2,KS)
     &+SO(IF+2,JF+2,KSW)+SO(IF+2,JF+1,KW)+SO(IF+2,JF+1,KNW)
     &+SO(IF+1,JF+1,KS)+SO(IF+1,JF+1,KSW)
      EP=MIN(ABS(SO(IF+1,JF+1,KSW)+SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KNW)),
     &       ABS(SO(IF+1,JF+2,KNW)+SO(IF+1,JF+2,KS)+SO(IF+2,JF+2,KSW)),
     &       ABS(SO(IF+2,JF+2,KSW)+SO(IF+2,JF+1,KW)+SO(IF+2,JF+1,KNW)),
     &       ABS(SO(IF+2,JF+1,KNW)+SO(IF+1,JF+1,KS)+SO(IF+1,JF+1,KSW)),
     &       ONE)
      SUM=SUM+(SO(IF+1,JF+1,KO)-SUM)*MAX(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*
     & CI(IC,JC+1,LB)+SO(IF+1,JF+1,KSW))*S
      CI(IC,JC,LSE)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LB)+SO(IF+2,JF+1,KNW))*S
      CI(IC,JC,LNW)=(SO(IF+1,JF+1,KW)*CI(IC,JC+1,LA)+SO(IF+1,JF+2,KS)*
     & CI(IC+1,JC+1,LL) +SO(IF+1,JF+2,KNW))*S
      CI(IC,JC,LNE)=(SO(IF+1,JF+2,KS)*CI(IC+1,JC+1,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LA)+SO(IF+2,JF+2,KSW))*S
  105 CONTINUE
C   end of computation of i when kf difference operator is nine point
C******************************
C   begin computation of grid kc difference operator when kf difference
C   operator is nine point unless kc. ge. icoef
C
      IF (KC.GE.ICOEF) GO TO 230
      JF=0
      DO 110 JC=2,JJC1
      JF=JF+2
      IF=0
      DO 109 IC=2,IIC1
      IF=IF+2
      CO=SO(IF,JF+1,KNW)*CI(IC-1,JC,LSW)+SO(IF,JF,KW)*CI(IC,JC,LL)+
     & SO(IF,JF,KSW)*CI(IC-1,JC-1,LNW)
      CS=SO(IF,JF-1,KW)*CI(IC-1,JC-1,LNW)+SO(IF,JF,KNW)*CI(IC,JC,LL)
      CSW=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LNW)+SO(IF-1,JF-1,KW)*
     &CI(IC-1,JC,LA)+SO(IF-1,JF,KNW)+SO(IF-1,JF,KS)*CI(IC,JC,LL)
      CW=-SO(IF-1,JF,KO)*CI(IC,JC,LL)+SO(IF-1,JF,KS)*CI(IC-1,JC-1,LNW)+
     & SO(IF-1,JF,KSW)*CI(IC-1,JC,LA)+SO(IF-1,JF,KW)+SO(IF-1,JF+1,KNW)*
     & CI(IC-1,JC+1,LB)+SO(IF-1,JF+1,KS)*CI(IC-1,JC,LSW)
      CNW=-SO(IF-1,JF+1,KO)*CI(IC-1,JC,LSW)+
     & SO(IF-1,JF+1,KS)*CI(IC,JC,LL)+SO(IF-1,JF+1,KSW)+
     & SO(IF-1,JF+1,KW)*CI(IC-1,JC+1,LB)
      CN=SO(IF,JF+1,KSW)*CI(IC,JC,LL)+SO(IF,JF+1,KW)*CI(IC-1,JC,LSW)
      SOC(IC,JC,KW)=CO+CI(IC,JC,LA)*CS+CI(IC-1,JC-1,LNE)*CSW+
     & CI(IC,JC,LR)*CW+CI(IC-1,JC,LSE)*CNW+CI(IC,JC+1,LB)*CN
      COSW=SO(IF,JF,KSW)*CI(IC-1,JC-1,LSW)
      CSSW=SO(IF,JF-1,KSW)*CI(IC,JC-1,LL)+
     & SO(IF,JF-1,KW)*CI(IC-1,JC-1,LSW)
      CSWSW=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LSW)+SO(IF-1,JF-1,KS)*
     & CI(IC,JC-1,LL)+SO(IF-1,JF-1,KSW)+SO(IF-1,JF-1,KW)*CI(IC-1,JC,LB)
      CWSW=SO(IF-1,JF,KS)*CI(IC-1,JC-1,LSW)+
     & SO(IF-1,JF,KSW)*CI(IC-1,JC,LB)
      SOC(IC,JC,KSW)=COSW+CI(IC,JC,LA)*CSSW+CI(IC-1,JC-1,LNE)*CSWSW+
     & CI(IC,JC,LR)*CWSW
      COA=SO(IF,JF,KSW)*CI(IC-1,JC-1,LSE)+SO(IF,JF,KS)*CI(IC,JC,LB)+
     &SO(IF+1,JF,KNW)*CI(IC,JC-1,LSW)
      CEA=SO(IF+1,JF,KS)*CI(IC,JC-1,LSW)+SO(IF+1,JF,KSW)*CI(IC,JC,LB)
      CSEA=-SO(IF+1,JF-1,KO)*CI(IC,JC-1,LSW)+SO(IF+1,JF-1,KS)*
     & CI(IC+1,JC-1,LL)+SO(IF+1,JF-1,KSW)+SO(IF+1,JF-1,KW)*CI(IC,JC,LB)
      CSA=-SO(IF,JF-1,KO)*CI(IC,JC,LB)+SO(IF+1,JF-1,KW)*CI(IC,JC-1,LSW)
     &+SO(IF+1,JF-1,LSW)*CI(IC+1,JC-1,LL)+SO(IF,JF-1,KS)+SO(IF,JF-1,KSW)
     &*CI(IC,JC-1,LR)+SO(IF,JF-1,KW)*CI(IC-1,JC-1,LSE)
      CSWA=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LSE)+SO(IF,JF-1,KW)*
     & CI(IC,JC,LB)+SO(IF,JF-1,KNW)+SO(IF-1,JF-1,KS)*CI(IC,JC-1,LR)
      CWA=SO(IF,JF,KNW)*CI(IC,JC,LB)+SO(IF-1,JF,KS)*CI(IC-1,JC-1,LSE)
      SOC(IC,JC,KS)=COA+CI(IC+1,JC,LL)*CEA+CI(IC,JC-1,LNW)*CSEA+
     & CI(IC,JC,LA)*CSA+CI(IC-1,JC-1,LNE)*CSWA+CI(IC,JC,LR)*CWA
      CONW=SO(IF-1,JF,KNW)*CI(IC-1,JC-1,LSE)
      CENW=SO(IF,JF,KNW)*CI(IC,JC,LB)+SO(IF-1,JF,KS)*CI(IC-1,JC-1,LSE)
      CSENW=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LSE)+SO(IF,JF-1,KW)*
     & CI(IC,JC,LB)+SO(IF,JF-1,KNW)+SO(IF-1,JF-1,KS)*CI(IC,JC-1,LR)
      CFNW=SO(IF-1,JF-1,KW)*CI(IC-1,JC-1,LSE)+SO(IF-1,JF-1,KNW)*
     & CI(IC,JC-1,LR)
      SOC(IC,JC,KNW)=CONW+CI(IC,JC,LL)*CENW+CI(IC-1,JC-1,LNW)*CSENW+
     & CI(IC-1,JC,LA)*CFNW
  109 CONTINUE
  110 CONTINUE

  118 JF=0
      DO 112 JC=2,JJC1
      JF=JF+2
      IF=0
      DO 111 IC=2,IIC1
      IF=IF+2
      CO=SO(IF,JF,KW)*CI(IC,JC,LR)+SO(IF,JF+1,KNW)*CI(IC-1,JC,LSE)+
     &SO(IF,JF+1,KS)*CI(IC,JC+1,LB)+SO(IF+1,JF+1,KSW)*CI(IC,JC,LSW)
     &+SO(IF+1,JF,KW)*CI(IC+1,JC,LL)+SO(IF+1,JF,KNW)*CI(IC,JC-1,LNW)
     &+SO(IF,JF,KS)*CI(IC,JC,LA)+SO(IF,JF,KSW)*CI(IC-1,JC-1,LNE)-
     & SO(IF,JF,KO)
      CW=-SO(IF-1,JF,KO)*CI(IC,JC,LR)+SO(IF-1,JF+1,KS)*CI(IC-1,JC,LSE)
     &+SO(IF,JF+1,KSW)*CI(IC,JC+1,LB)+SO(IF,JF,KW)+SO(IF,JF,KNW)*
     &CI(IC,JC,LA)+SO(IF-1,JF,KS)*CI(IC-1,JC-1,LNE)
      CNW=-SO(IF-1,JF+1,KO)*CI(IC-1,JC,LSE)+
     & SO(IF,JF+1,KW)*CI(IC,JC+1,LB)+SO(IF,JF+1,KNW)+
     & SO(IF-1,JF+1,KS)*CI(IC,JC,LR)
      CN=-SO(IF,JF+1,KO)*CI(IC,JC+1,LB)+SO(IF+1,JF+1,KNW)*CI(IC+1,JC,LL)
     & +SO(IF,JF+1,KS)+SO(IF,JF+1,KSW)*CI(IC,JC,LR)+SO(IF,JF+1,KW)*
     & CI(IC-1,JC,LSE)+SO(IF+1,JF+1,KW)*CI(IC,JC,LSW)
      CNE=-SO(IF+1,JF+1,KO)*CI(IC,JC,LSW)+
     & SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KSW)+
     & SO(IF+1,JF+1,KW)*CI(IC,JC+1,LB)
      CE=-SO(IF+1,JF,KO)*CI(IC+1,JC,LL)+SO(IF+1,JF,KS)*CI(IC,JC-1,LNW)+
     & SO(IF+1,JF,KSW)*CI(IC,JC,LA)+SO(IF+1,JF,KW)+SO(IF+1,JF+1,KNW)*
     & CI(IC,JC+1,LB)+SO(IF+1,JF+1,KS)*CI(IC,JC,LSW)
      CSE=-SO(IF+1,JF-1,KO)*CI(IC,JC-1,LNW)+
     & SO(IF+1,JF-1,KW)*CI(IC,JC,LA)+SO(IF+1,JF,KNW)+
     & SO(IF+1,JF,KS)*CI(IC+1,JC,LL)
      CS=-SO(IF,JF-1,KO)*CI(IC,JC,LA)+SO(IF,JF-1,KW)*CI(IC-1,JC-1,LNE)+
     &SO(IF,JF,KNW)*CI(IC,JC,LR)+SO(IF,JF,KS)+SO(IF+1,JF,KSW)*
     &CI(IC+1,JC,LL)+SO(IF+1,JF-1,KW)*CI(IC,JC-1,LNW)
      CSW=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LNE)+
     & SO(IF-1,JF,KS)*CI(IC,JC,LR)+SO(IF,JF,KSW)+
     & SO(IF,JF-1,KW)*CI(IC,JC,LA)
      SOC(IC,JC,KO)=-CI(IC-1,JC,LSE)*CNW-CI(IC,JC+1,LB)*CN-CI(IC,JC,LSW)
     & *CNE-CI(IC,JC,LR)*CW-CO-CI(IC+1,JC,LL)*CE-CI(IC-1,JC-1,LNE)*CSW-
     & CI(IC,JC,LA)*CS-CI(IC,JC-1,LNW)*CSE
      SORC(IC,JC,MSOR)=ONE/SOC(IC,JC,KO)
  111 CONTINUE
  112 CONTINUE
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 114
      DO 113 IC=1,IIC
      SOC(IC,JJC,KS)=SOC(IC,2,KS)
      SOC(IC,JJC,KW)=SOC(IC,2,KW)
      SOC(IC,JJC,KNW)=SOC(IC,2,KNW)
      SOC(IC,JJC,KSW)=SOC(IC,2,KSW)
      SOC(IC,JJC,KO)=SOC(IC,2,KO)
      SOC(IC,1,KO)=SOC(IC,JJC1,KO)
      SOC(IC,1,KW)=SOC(IC,JJC1,KW)
      SOC(IC,1,KS)=SOC(IC,JJC1,KS)
      SOC(IC,1,KSW)=SOC(IC,JJC1,KSW)
      SOC(IC,1,KNW)=SOC(IC,JJC1,KNW)
  113 CONTINUE
  114 IF(IPN.NE.2.AND.IPN.NE.3)GO TO 116
      DO 115 JC=1,JJC
      SOC(IIC,JC,KW)=SOC(2,JC,KW)
      SOC(IIC,JC,KS)=SOC(2,JC,KS)
      SOC(IIC,JC,KNW)=SOC(2,JC,KNW)
      SOC(IIC,JC,KSW)=SOC(2,JC,KSW)
      SOC(IIC,JC,KO)=SOC(2,JC,KO)
      SOC(1,JC,KW)=SOC(IIC1,JC,KW)
      SOC(1,JC,KS)=SOC(IIC1,JC,KS)
      SOC(1,JC,KSW)=SOC(IIC1,JC,KSW)
      SOC(1,JC,KNW)=SOC(IIC1,JC,KNW)
      SOC(1,JC,KO)=SOC(IIC1,JC,KO)
  115 CONTINUE
  116 CONTINUE
      GO TO 230
C   end of computation of kc difference operator when kf difference
C   operator is nine point
C******************************
C   begin computation of i when kf difference operator is five point
C
  120 JF=0
      DO 160 JC=2,JJC1
      JF=JF+2
      IF=2
      IBEGC=3
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 149
      IF=0
      IBEGC=2
  149 DO 150 IC=IBEGC,IICF1
      IF=IF+2
      A=SO(IF,JF,KW)
      B=SO(IF-1,JF,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF-1,JF,KS)+SO(IF-1,JF+1,KS)
      SUM=A+B+(SO(IF-1,JF,KO)-SUM)*MAX(SO(IF-1,JF,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF-1,JF,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
  150 CONTINUE
  160 CONTINUE
      IF(IIC.NE.IICF.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 155
      IF=3
      IC=IIC
      JF=0
      DO 153 JC=2,JJC1
      JF=JF+2
      A=SO(IF,JF,KW)
      B=SO(IF-1,JF,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF-1,JF,KS)+SO(IF-1,JF+1,KS)
      SUM=A+B+(SO(IF-1,JF,KO)-SUM)*MAX(SO(IF-1,JF,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF-1,JF,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
  153 CONTINUE
  155 IF(IPN.NE.1 .AND. IPN.NE.3)GO TO 180
      IF(JJCF.EQ.JJC)GO TO 163
      DO 161 IC=2,IIC
      CI(IC,1,LL)=CI(IC,JJC1,LL)
      CI(IC,1,LR)=CI(IC,JJC1,LR)
      CI(IC,JJC,LL)=CI(IC,2,LL)
      CI(IC,JJC,LR)=CI(IC,2,LR)
  161 CONTINUE
      GO TO 180
  163 JF=JJF2
      JC=1
      INDX=0
      GO TO 165
  164 JF=3
      JC=JJC
  165 IBEGC=3
      IF=2
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 167
      IBEGC=2
      IF=0
  167 DO 166 IC=IBEGC,IICF1
      IF=IF+2
      A=SO(IF,JF,KW)
      B=SO(IF-1,JF,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF-1,JF,KS)+SO(IF-1,JF+1,KS)
      SUM=A+B+(SO(IF-1,JF,KO)-SUM)*MAX(SO(IF-1,JF,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF-1,JF,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
  166 CONTINUE
      IF(IICF.NE.IIC.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 168
      IF=3
      IC=IIC
      A=SO(IF,JF,KW)
      B=SO(IF-1,JF,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF-1,JF,KS)+SO(IF-1,JF+1,KS)
      SUM=A+B+(SO(IF-1,JF,KO)-SUM)*MAX(SO(IF-1,JF,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF-1,JF,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
  168 IF(INDX.NE.0)GO TO 180
      INDX=1
      GO TO 164
  180 JF=2
      JBEGC=3
      IF(IPN.NE.1 .AND. IPN.NE.3)GO TO 175
      JF=0
      JBEGC=2
  175 DO 190 JC=JBEGC,JJCF1
      JF=JF+2
      IF=0
      DO 189 IC=2,IIC1
      IF=IF+2
      A=SO(IF,JF,KS)
      B=SO(IF,JF-1,KS)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF,JF-1,KW)+SO(IF+1,JF-1,KW)
      SUM=A+B+(SO(IF,JF-1,KO)-SUM)*MAX(SO(IF,JF-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF,JF-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LA)=A*SUM
      CI(IC,JC,LB)=B*SUM
  189 CONTINUE
      IF(IICF.NE.IIC.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 169
      INDX=0
      IF=IIF2
      IC=1
      GO TO 188
  187 IF=3
      IC=IIC
  188 A=SO(IF,JF,KS)
      B=SO(IF,JF-1,KS)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF,JF-1,KW)+SO(IF+1,JF-1,KW)
      SUM=A+B+(SO(IF,JF-1,KO)-SUM)*MAX(SO(IF,JF-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF,JF-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LA)=A*SUM
      CI(IC,JC,LB)=B*SUM
      IF(INDX.NE.0)GO TO 190
      INDX=1
      GO TO 187
  169 IF(IICF.EQ.IIC.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 190
      CI(1,JC,LA)=CI(IIC1,JC,LA)
      CI(1,JC,LB)=CI(IIC1,JC,LB)
      CI(IIC,JC,LA)=CI(2,JC,LA)
      CI(IIC,JC,LB)=CI(2,JC,LB)
  190 CONTINUE
      IF(JJC.NE.JJCF.OR.(IPN.NE.1.AND.IPN.NE.3))GO TO 194
      JF=3
      IF=0
      DO 191 IC=2,IIC1
      IF=IF+2
      A=SO(IF,JF,KS)
      B=SO(IF,JF-1,KS)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF,JF-1,KW)+SO(IF+1,JF-1,KW)
      SUM=A+B+(SO(IF,JF-1,KO)-SUM)*MAX(SO(IF,JF-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF,JF-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JJC,LA)=A*SUM
      CI(IC,JJC,LB)=B*SUM
  191 CONTINUE
      IF(IICF.NE.IIC.OR.IPN.NE.3)GO TO 1935
      INDX=0
      IF=IIF2
      IC=1
      GO TO 193
  192 IF=3
      IC=IIC
  193 A=SO(IF,JF,KS)
      B=SO(IF,JF-1,KS)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(IF,JF-1,KW)+SO(IF+1,JF-1,KW)
      SUM=A+B+(SO(IF,JF-1,KO)-SUM)*MAX(SO(IF,JF-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(IF,JF-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JJC,LA)=A*SUM
      CI(IC,JJC,LB)=B*SUM
      IF(INDX.NE.0)GO TO 194
      INDX=1
      GO TO 192
 1935 IF (IICF.EQ.IIC .OR. IPN.NE.3) GOTO 194
      CI(1,JJC,LA)=CI(IIC1,JJC,LA)
      CI(1,JJC,LB)=CI(IIC1,JJC,LB)
      CI(IIC,JJC,LA)=CI(2,JJC,LA)
      CI(IIC,JJC,LB)=CI(2,JJC,LB)      

  194 JF=0
      JBEGC=2
      IF(IPN.NE.1.AND.IPN.NE.3)GO TO 195
      JBEGC=1
      JF=-2
  195 DO 210 JC=JBEGC,JJCF2
      JF=JF+2
      IF=0
      IBEGC=2
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 209
      IBEGC=1
      IF=-2
  209 DO 200 IC=IBEGC,IICF2
      IF=IF+2
      SUM=SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KS)+SO(IF+2,JF+1,KW)
     &+SO(IF+1,JF+1,KS)
      EP=MIN(ABS(SO(IF+1,JF+1,KW)),ABS(SO(IF+1,JF+2,KS)),
     &       ABS(SO(IF+2,JF+1,KW)),ABS(SO(IF+1,JF+1,KS)),ONE)
      SUM=SUM+(SO(IF+1,JF+1,KO)-SUM)*MAX(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*
     & CI(IC,JC+1,LB))*S
      CI(IC,JC,LSE)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LB))*S
      CI(IC,JC,LNW)=(SO(IF+1,JF+1,KW)*CI(IC,JC+1,LA)+SO(IF+1,JF+2,KS)*
     & CI(IC+1,JC+1,LL))*S
      CI(IC,JC,LNE)=(SO(IF+1,JF+2,KS)*CI(IC+1,JC+1,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LA))*S
  200 CONTINUE
      IF(IIC.NE.IICF.OR.(IPN.NE.2.AND.IPN.NE.3))GO TO 210
      IF=1
      IC=IIC1
      SUM=SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KS)+SO(IF+2,JF+1,KW)
     &+SO(IF+1,JF+1,KS)
      EP=MIN(ABS(SO(IF+1,JF+1,KW)),ABS(SO(IF+1,JF+2,KS)),
     &       ABS(SO(IF+2,JF+1,KW)),ABS(SO(IF+1,JF+1,KS)),ONE)
      SUM=SUM+(SO(IF+1,JF+1,KO)-SUM)*MAX(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*
     & CI(IC,JC+1,LB))*S
      CI(IC,JC,LSE)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LB))*S
      CI(IC,JC,LNW)=(SO(IF+1,JF+1,KW)*CI(IC,JC+1,LA)+SO(IF+1,JF+2,KS)*
     & CI(IC+1,JC+1,LL))*S
      CI(IC,JC,LNE)=(SO(IF+1,JF+2,KS)*CI(IC+1,JC+1,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LA))*S
  210 CONTINUE
      IF(JJC.NE.JJCF.OR.(IPN.NE.1.AND.IPN.NE.3))GO TO 215
      JF=1
      JC=JJC1
      IF=0
      IBEGC=2
      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 205
      IBEGC=1
      IF=-2
  205 DO 211 IC=IBEGC,IICF2
      IF=IF+2
      SUM=SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KS)+SO(IF+2,JF+1,KW)
     &+SO(IF+1,JF+1,KS)
      EP=MIN(ABS(SO(IF+1,JF+1,KW)),ABS(SO(IF+1,JF+2,KS)),
     &       ABS(SO(IF+2,JF+1,KW)),ABS(SO(IF+1,JF+1,KS)),ONE)
      SUM=SUM+(SO(IF+1,JF+1,KO)-SUM)*MAX(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*
     & CI(IC,JC+1,LB))*S
      CI(IC,JC,LSE)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LB))*S
      CI(IC,JC,LNW)=(SO(IF+1,JF+1,KW)*CI(IC,JC+1,LA)+SO(IF+1,JF+2,KS)*
     & CI(IC+1,JC+1,LL))*S
      CI(IC,JC,LNE)=(SO(IF+1,JF+2,KS)*CI(IC+1,JC+1,LR)+SO(IF+2,JF+1,KW)*
     & CI(IC+1,JC+1,LA))*S
  211 CONTINUE
      IF(IIC.NE.IICF.OR.IPN.NE.3)GO TO 215
      IF=1
      IC=IIC1
      SUM=SO(IF+1,JF+1,KW)+SO(IF+1,JF+2,KS)+SO(IF+2,JF+1,KW)
     &+SO(IF+1,JF+1,KS)
      EP=MIN(ABS(SO(IF+1,JF+1,KW)),ABS(SO(IF+1,JF+2,KS)),
     &       ABS(SO(IF+2,JF+1,KW)),ABS(SO(IF+1,JF+1,KS)),ONE)
      SUM=SUM+(SO(IF+1,JF+1,KO)-SUM)*MAX(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(IF+1,JF+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*
     1 CI(IC,JC+1,LB))*S
      CI(IC,JC,LSE)=(SO(IF+1,JF+1,KS)*CI(IC+1,JC,LR)+SO(IF+2,JF+1,KW)*
     1 CI(IC+1,JC+1,LB))*S
      CI(IC,JC,LNW)=(SO(IF+1,JF+1,KW)*CI(IC,JC+1,LA)+SO(IF+1,JF+2,KS)*
     1 CI(IC+1,JC+1,LL))*S
      CI(IC,JC,LNE)=(SO(IF+1,JF+2,KS)*CI(IC+1,JC+1,LR)+SO(IF+2,JF+1,KW)*
     1 CI(IC+1,JC+1,LA))*S
  215 CONTINUE
C   end computation of i when kf difference operator is five point
C******************************
C   begin computation of kc difference operator when kf difference
C   operator is five point unless kc.ge.icoef
C
      IF (KC.GE.ICOEF) GO TO 230
      JF=0
      DO 220 JC=2,JJC1
      JF=JF+2
      IF=0
      DO 219 IC=2,IIC1
      IF=IF+2
      CO=SO(IF,JF,KW)*CI(IC,JC,LL)
      CS=SO(IF,JF-1,KW)*CI(IC-1,JC-1,LNW)
      CSW=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LNW)+SO(IF-1,JF-1,KW)*
     & CI(IC-1,JC,LA)+SO(IF-1,JF,KS)*CI(IC,JC,LL)
      CW=-SO(IF-1,JF,KO)*CI(IC,JC,LL)+SO(IF-1,JF,KS)*CI(IC-1,JC-1,LNW)+
     & SO(IF-1,JF,KW)+SO(IF-1,JF+1,KS)*CI(IC-1,JC,LSW)
      CNW=-SO(IF-1,JF+1,KO)*CI(IC-1,JC,LSW)+
     & SO(IF-1,JF+1,KS)*CI(IC,JC,LL)+SO(IF-1,JF+1,KW)*CI(IC-1,JC+1,LB)
      CN=SO(IF,JF+1,KW)*CI(IC-1,JC,LSW)
      SOC(IC,JC,KW)=CO+CI(IC,JC,LA)*CS+CI(IC-1,JC-1,LNE)*CSW+
     & CI(IC,JC,LR)*CW+CI(IC-1,JC,LSE)*CNW+CI(IC,JC+1,LB)*CN
      CSSW=SO(IF,JF-1,KW)*CI(IC-1,JC-1,LSW)
      CSWSW=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LSW)+SO(IF-1,JF-1,KS)*
     & CI(IC,JC-1,LL)+SO(IF-1,JF-1,KW)*CI(IC-1,JC,LB)
      CWSW=SO(IF-1,JF,KS)*CI(IC-1,JC-1,LSW)
      SOC(IC,JC,KSW)=CI(IC,JC,LA)*CSSW+CI(IC-1,JC-1,LNE)*CSWSW+
     & CI(IC,JC,LR)*CWSW
      COA=SO(IF,JF,KS)*CI(IC,JC,LB)
      CEA=SO(IF+1,JF,KS)*CI(IC,JC-1,LSW)
      CSEA=-SO(IF+1,JF-1,KO)*CI(IC,JC-1,LSW)+SO(IF+1,JF-1,KS)*
     & CI(IC+1,JC-1,LL)+SO(IF+1,JF-1,KW)*CI(IC,JC,LB)
      CSA=-SO(IF,JF-1,KO)*CI(IC,JC,LB)+SO(IF+1,JF-1,KW)*CI(IC,JC-1,LSW)+
     & SO(IF,JF-1,KS)+SO(IF,JF-1,KW)*CI(IC-1,JC-1,LSE)
      CSWA=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LSE)+
     & SO(IF,JF-1,KW)*CI(IC,JC,LB)+SO(IF-1,JF-1,KS)*CI(IC,JC-1,LR)
      CWA=SO(IF-1,JF,KS)*CI(IC-1,JC-1,LSE)
      SOC(IC,JC,KS)=COA+CI(IC+1,JC,LL)*CEA+CI(IC,JC-1,LNW)*CSEA+
     & CI(IC,JC,LA)*CSA+CI(IC-1,JC-1,LNE)*CSWA+CI(IC,JC,LR)*CWA
      CO=SO(IF,JF,KW)*CI(IC,JC,LR)+SO(IF,JF+1,KS)*CI(IC,JC+1,LB)+
     & SO(IF+1,JF,KW)*CI(IC+1,JC,LL)+SO(IF,JF,KS)*CI(IC,JC,LA)-
     & SO(IF,JF,KO)
      CW=-SO(IF-1,JF,KO)*CI(IC,JC,LR)+SO(IF-1,JF+1,KS)*CI(IC-1,JC,LSE)+
     & SO(IF,JF,KW)+SO(IF-1,JF,KS)*CI(IC-1,JC-1,LNE)
      CNW=-SO(IF-1,JF+1,KO)*CI(IC-1,JC,LSE)+SO(IF,JF+1,KW)*
     & CI(IC,JC+1,LB)+SO(IF-1,JF+1,KS)*CI(IC,JC,LR)
      CN=-SO(IF,JF+1,KO)*CI(IC,JC+1,LB)+SO(IF,JF+1,KS)+SO(IF,JF+1,KW)*
     & CI(IC-1,JC,LSE)+SO(IF+1,JF+1,KW)*CI(IC,JC,LSW)
      CNE=-SO(IF+1,JF+1,KO)*CI(IC,JC,LSW)+
     & SO(IF+1,JF+1,KS)*CI(IC+1,JC,LL)+SO(IF+1,JF+1,KW)*CI(IC,JC+1,LB)
      CE=-SO(IF+1,JF,KO)*CI(IC+1,JC,LL)+SO(IF+1,JF,KS)*CI(IC,JC-1,LNW)+
     & SO(IF+1,JF,KW)+SO(IF+1,JF+1,KS)*CI(IC,JC,LSW)
      CSE=-SO(IF+1,JF-1,KO)*CI(IC,JC-1,LNW)+
     & SO(IF+1,JF-1,KW)*CI(IC,JC,LA)+SO(IF+1,JF,KS)*CI(IC+1,JC,LL)
      CS=-SO(IF,JF-1,KO)*CI(IC,JC,LA)+SO(IF,JF-1,KW)*CI(IC-1,JC-1,LNE)+
     & SO(IF,JF,KS)+SO(IF+1,JF-1,KW)*CI(IC,JC-1,LNW)
      CSW=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LNE)+
     & SO(IF-1,JF,KS)*CI(IC,JC,LR)+SO(IF,JF-1,KW)*CI(IC,JC,LA)
      SOC(IC,JC,KO)=-CI(IC-1,JC,LSE)*CNW-CI(IC,JC+1,LB)*CN-CI(IC,JC,LSW)
     & *CNE-CI(IC,JC,LR)*CW-CO-CI(IC+1,JC,LL)*CE-CI(IC-1,JC-1,LNE)*CSW-
     & CI(IC,JC,LA)*CS-CI(IC,JC-1,LNW)*CSE
      SORC(IC,JC,MSOR)=ONE/SOC(IC,JC,KO)
      CENW=SO(IF-1,JF,KS)*CI(IC-1,JC-1,LSE)
      CSENW=-SO(IF-1,JF-1,KO)*CI(IC-1,JC-1,LSE)+SO(IF,JF-1,KW)*
     & CI(IC,JC,LB)+SO(IF-1,JF-1,KS)*CI(IC,JC-1,LR)
      CSNW=SO(IF-1,JF-1,KW)*CI(IC-1,JC-1,LSE)
      SOC(IC,JC,KNW)=CI(IC,JC,LL)*CENW+CI(IC-1,JC-1,LNW)*CSENW+
     & CI(IC-1,JC,LA)*CSNW
  219 CONTINUE
  220 CONTINUE
  226 IF(IPN.NE.1.AND.IPN.NE.3)GO TO 222
      DO 221 IC=1,IIC
      SOC(IC,JJC,KS)=SOC(IC,2,KS)
      SOC(IC,JJC,KW)=SOC(IC,2,KW)
      SOC(IC,JJC,KNW)=SOC(IC,2,KNW)
      SOC(IC,JJC,KSW)=SOC(IC,2,KSW)
      SOC(IC,JJC,KO)=SOC(IC,2,KO)
      SOC(IC,1,KO)=SOC(IC,JJC1,KO)
      SOC(IC,1,KW)=SOC(IC,JJC1,KW)
      SOC(IC,1,KS)=SOC(IC,JJC1,KS)
      SOC(IC,1,KSW)=SOC(IC,JJC1,KSW)
      SOC(IC,1,KNW)=SOC(IC,JJC1,KNW)
  221 CONTINUE
  222 IF(IPN.NE.2.AND.IPN.NE.3)GO TO 224
      DO 223 JC=1,JJC
      SOC(IIC,JC,KW)=SOC(2,JC,KW)
      SOC(IIC,JC,KS)=SOC(2,JC,KS)
      SOC(IIC,JC,KNW)=SOC(2,JC,KNW)
      SOC(IIC,JC,KSW)=SOC(2,JC,KSW)
      SOC(IIC,JC,KO)=SOC(2,JC,KO)
      SOC(1,JC,KW)=SOC(IIC1,JC,KW)
      SOC(1,JC,KS)=SOC(IIC1,JC,KS)
      SOC(1,JC,KSW)=SOC(IIC1,JC,KSW)
      SOC(1,JC,KNW)=SOC(IIC1,JC,KNW)
      SOC(1,JC,KO)=SOC(IIC1,JC,KO)
  223 CONTINUE
  224 CONTINUE
C   end of computation of grid kc difference operator, when kf
C   difference operator is five point
C******************************

  230 IF (ISTRT.LT.0) GO TO 250  ! IF FMG cycles then restrict QF

C
C   unless istrt.lt.0 form right hand side for grid kc
C   by weighting right hand side of grid kf with i(transpose).
C

C     Ensure periodicity prior to restriction 
C     (consistent with mgrcap.f)

      IF ( IPN.EQ.2 .OR. IPN.EQ.3 ) THEN
         DO JF=2,JJF1
            QF(1,JF)  = QF(IIF1,JF)
            QF(IIF,JF)= QF(2,JF)
         ENDDO
      ENDIF
         
      IF ( IPN.EQ.1 .OR. IPN.EQ.3 ) THEN
         DO IF=2,IIF1
            QF(IF,1)  = QF(IF,JJF1)
            QF(IF,JJF)= QF(IF,2)
         ENDDO
      ENDIF

      JF=0
      DO 241 JC=2,JJC1
         !
         JF=JF+2
         IF=0
         !
         DO 240 IC=2,IIC1
            !
            IF=IF+2
            !
            QFC(IC,JC) = CI(IC-1,JC-1,LNE)*QF(IF-1,JF-1)
     &                 + CI(IC,JC,LA)*QF(IF,JF-1)
     &                 + CI(IC,JC-1,LNW)*QF(IF+1,JF-1)
     &                 + CI(IC,JC,LR)*QF(IF-1,JF)
     &                 + QF(IF,JF)
     &                 + CI(IC+1,JC,LL)*QF(IF+1,JF)
     &                 + CI(IC-1,JC,LSE)*QF(IF-1,JF+1)
     &                 + CI(IC,JC+1,LB)*QF(IF,JF+1)
     &                 + CI(IC,JC,LSW)*QF(IF+1,JF+1)
            !
            SORC(IC,JC,MTOT) = QFC(IC,JC)
            !
 240     CONTINUE
 241  CONTINUE

      IF( ISKIP.EQ.2 ) RETURN   ! Skip setup of relaxation 
C
C   if irelax=1, point relaxation. hence, return.
C   if irelax=2, relaxation by lines in x. form lu decomposition
C   of tridiagonal matrices along y=const. lines. if irelax=3,
C   relaxation by lines in y. form lu decomposition of tridiagonal
C   matrices along x=const. lines. if irelax=4, relaxation by lines in x
C   and y. hence need lu decomposition of both sets of tridiagonal
C   matrices.
C
 250  GO TO (260,270,290,270), IRELAX

 260  RETURN                    ! point relaxation so return

      ! Factorization for lines in x
 270  DO 281 JF=2,JJF1
         SOR(2,JF,MSOR)=ONE/SO(2,JF,KO)
         SOR(2,JF,MTOT)=QF(2,JF)*SOR(2,JF,MSOR)
         DO 280 IF=3,IIF1
            SOR(IF,JF,MSOR)=ONE
     &           /(SO(IF,JF,KO)-SOR(IF-1,JF,MSOR)*SO(IF,JF,KW)**2)
            SOR(IF,JF,MTOT)=QF(IF,JF)/SO(IF,JF,KO)
 280     CONTINUE
 281  CONTINUE

      IF(IPN.NE.2.AND.IPN.NE.3)GO TO 296  ! not periodic in x jump to y

      ! Periodic fixup of factorization (use ghost part of QF for workspace)
      DO 294 JF=2,JJF1
         QF(1,JF)=SO(IIF,JF,KW)
         QF(IIF,JF)=SO(2,JF,KW)*SOR(2,JF,MSOR)
         SOR(IIF1,JF,MSOR)=SO(IIF1,JF,KO)-QF(1,JF)*QF(IIF,JF)
         DO 293 IF=3,IIF2
            QF(IIF,JF)=QF(IIF,JF)*SO(IF,JF,KW)*SOR(IF,JF,MSOR)
            QF(1,JF)=QF(1,JF)*SO(IF,JF,KW)*SOR(IF-1,JF,MSOR)
            SOR(IIF1,JF,MSOR)=SOR(IIF1,JF,MSOR)-QF(IIF,JF)*QF(1,JF)
 293     CONTINUE
         SOR(IIF1,JF,MSOR)=ONE/( SOR(IIF1,JF,MSOR) 
     &        -(SO(IIF1,JF,KW)+QF(1,JF))
     &        *SOR(IIF2,JF,MSOR)*SO(IIF1,JF,KW)
     &        -SO(IIF1,JF,KW)*QF(IIF,JF) )
 294  CONTINUE
      ! Restore periodicity in QF workspace
      DO 295 JF=2,JJF1
         QF(1,JF)  = QF(IIF1,JF)
         QF(IIF,JF)= QF(2,JF)
 295  CONTINUE

 296  IF (IRELAX.NE.4) RETURN  ! just line in x so return

      ! Factorization for lines in y
 290  DO 301 IF=2,IIF1
         SOR(IF,2,MSOS)=ONE/SO(IF,2,KO)
         SOR(IF,2,MTOT)=QF(IF,2)*SOR(IF,2,MSOS)
         DO 300 JF=3,JJF1
            SOR(IF,JF,MSOS)=ONE
     &           /(SO(IF,JF,KO)-SOR(IF,JF-1,MSOS)*SO(IF,JF,KS)**2)
            SOR(IF,JF,MTOT)=QF(IF,JF)/SO(IF,JF,KO)
 300     CONTINUE
 301  CONTINUE

      IF (IPN.NE.1 .AND. IPN.NE.3) RETURN ! not periodic in y

      ! Periodic fixup of factorization (use ghost part of QF for workspace)
      DO 310 IF=2,IIF1
         QF(IF,1)=SO(IF,JJF,KS)
         QF(IF,JJF)=SO(IF,2,KS)*SOR(IF,2,MSOS)
         SOR(IF,JJF1,MSOS)=SO(IF,JJF1,KO)-QF(IF,1)*QF(IF,JJF)
         DO 305 JF=3,JJF2
            QF(IF,JJF)=QF(IF,JJF)*SO(IF,JF,KS)*SOR(IF,JF,MSOS)
            QF(IF,1)=QF(IF,1)*SO(IF,JF,KS)*SOR(IF,JF-1,MSOS)
            SOR(IF,JJF1,MSOS)=SOR(IF,JJF1,MSOS)-QF(IF,JJF)*QF(IF,1)
 305     CONTINUE
         SOR(IF,JJF1,MSOS)=ONE/( SOR(IF,JJF1,MSOS)
     &        -(SO(IF,JJF1,KS)+QF(IF,1))
     &        *SOR(IF,JJF2,MSOS)*SO(IF,JJF1,KS)
     &        -SO(IF,JJF1,KS)*QF(IF,JJF) )
 310  CONTINUE
      ! Restore periodicity in QF workspace
      DO 315 IF=2,IIF1
         QF(IF,1)  = QF(IF,JJF1)
         QF(IF,JJF)= QF(IF,2)
 315  CONTINUE
C     
      RETURN
      END
