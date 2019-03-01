      SUBROUTINE BMG2_SER_PerSymStd_SOLVE_cg(
     &                Q, QF, II, JJ, ABD, BBD, NABD1, NABD2, JPN 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_PerSymStd_SOLVE_cg does a direct solve on the coarsest grid. 
C     It uses the LAPACK routine DPBTRS or DPOTRS. 
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     QF        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C     ABD       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C     NABD1     Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     IPN       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     ABD       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     BBD       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     NABD1     Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     Q         Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C     CALLING ARGUMENTS

      integer II, JJ, JPN, NABD1, NABD2
      REAL*RKIND ABD(NABD1,NABD2), BBD(NABD2), Q(II,JJ), QF(II,JJ)

C     LOCAL VARIABLES

      integer I, I1, I2, IPN, J, J1, KK, N, INFO
      REAL*RKIND C, CINT, QINT, RZERO
C
C   direct solve on coarsest grid
C
C***FIRST EXECUTABLE STATEMENT  MGSADP
C
      RZERO = 0

      IPN=IABS(JPN)
      I1=II-1
      J1=JJ-1
      I2=I1-1
      N=I2*(J1-1)
      KK=0

      DO 11 J=2,J1
         DO 10 I=2,I1
            KK=KK+1
            BBD(KK)=QF(I,J)
 10      CONTINUE
 11   CONTINUE

C
C   spbsl and sposl are linpack routines.
C
      IF(IPN.NE.1.AND.IPN.NE.2.AND.IPN.NE.3) THEN
#if RKIND == 4
         CALL SPBTRS('U',KK,I1,1,ABD,NABD1,BBD,NABD2,INFO)
#else
         CALL DPBTRS('U',KK,I1,1,ABD,NABD1,BBD,NABD2,INFO)
C         CALL DPBSL(ABD,NABD1,N,I1,BBD)
#endif
      ELSE
#if RKIND == 4
         CALL SPOTRS('U',KK,1,ABD,NABD1,BBD,NABD2,INFO)
#else
         CALL DPOTRS('U',KK,1,ABD,NABD1,BBD,NABD2,INFO)
C         CALL DPOSL(ABD,NABD1,N,BBD)
#endif
      ENDIF
      
      KK=0
      DO 21 J=2,J1
         DO 20 I=2,I1
            KK=KK+1
            Q(I,J)=BBD(KK)
 20      CONTINUE
 21   CONTINUE

      IF ( JPN.GE.0 ) GO TO 60

      CINT=RZERO
      QINT=RZERO
      DO 71 J=2,J1
         DO 70 I=2,I1
            QINT=QINT+Q(I,J)
            CINT=CINT+1
 70      CONTINUE
 71   CONTINUE
      C=-QINT/CINT
      DO 81 J=2,J1
         DO 80 I=2,I1
            Q(I,J)=Q(I,J)+C
 80      CONTINUE
 81   CONTINUE

 60   CONTINUE

      IF(IPN.NE.1.AND.IPN.NE.3) GO TO 40

      DO 30 I=2,I1
         Q(I,JJ)=Q(I,2)
         Q(I,1)=Q(I,J1)
 30   CONTINUE

 40   CONTINUE

      IF ( IPN.NE.2 .AND. IPN.NE.3 ) RETURN

      DO 50 J=2,J1
         Q(II,J)=Q(2,J)
         Q(1,J)=Q(I1,J)
 50   CONTINUE

      IF ( IPN.NE.3 ) RETURN

      Q(1,1)   = Q(I1,J1)
      Q(II,1)  = Q(2,J1)
      Q(1,JJ)  = Q(I1,2)
      Q(II,JJ) = Q(2,2)

C ==========================================================================

      RETURN
      END
